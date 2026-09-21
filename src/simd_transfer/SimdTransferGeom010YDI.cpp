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


#include "SimdTransferGeom010YDI.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_geom_010y_di_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                const size_t target, const size_t pi_1, const size_t pi_0,
                                const size_t pk_1, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);

        const auto *pi_1_0 = buffer.data(pi_1 + 0 * ncomps + c);
        const auto *pi_1_1 = buffer.data(pi_1 + 1 * ncomps + c);
        const auto *pi_1_2 = buffer.data(pi_1 + 2 * ncomps + c);
        const auto *pi_1_3 = buffer.data(pi_1 + 3 * ncomps + c);
        const auto *pi_1_4 = buffer.data(pi_1 + 4 * ncomps + c);
        const auto *pi_1_5 = buffer.data(pi_1 + 5 * ncomps + c);
        const auto *pi_1_6 = buffer.data(pi_1 + 6 * ncomps + c);
        const auto *pi_1_7 = buffer.data(pi_1 + 7 * ncomps + c);
        const auto *pi_1_8 = buffer.data(pi_1 + 8 * ncomps + c);
        const auto *pi_1_9 = buffer.data(pi_1 + 9 * ncomps + c);
        const auto *pi_1_10 = buffer.data(pi_1 + 10 * ncomps + c);
        const auto *pi_1_11 = buffer.data(pi_1 + 11 * ncomps + c);
        const auto *pi_1_12 = buffer.data(pi_1 + 12 * ncomps + c);
        const auto *pi_1_13 = buffer.data(pi_1 + 13 * ncomps + c);
        const auto *pi_1_14 = buffer.data(pi_1 + 14 * ncomps + c);
        const auto *pi_1_15 = buffer.data(pi_1 + 15 * ncomps + c);
        const auto *pi_1_16 = buffer.data(pi_1 + 16 * ncomps + c);
        const auto *pi_1_17 = buffer.data(pi_1 + 17 * ncomps + c);
        const auto *pi_1_18 = buffer.data(pi_1 + 18 * ncomps + c);
        const auto *pi_1_19 = buffer.data(pi_1 + 19 * ncomps + c);
        const auto *pi_1_20 = buffer.data(pi_1 + 20 * ncomps + c);
        const auto *pi_1_21 = buffer.data(pi_1 + 21 * ncomps + c);
        const auto *pi_1_22 = buffer.data(pi_1 + 22 * ncomps + c);
        const auto *pi_1_23 = buffer.data(pi_1 + 23 * ncomps + c);
        const auto *pi_1_24 = buffer.data(pi_1 + 24 * ncomps + c);
        const auto *pi_1_25 = buffer.data(pi_1 + 25 * ncomps + c);
        const auto *pi_1_26 = buffer.data(pi_1 + 26 * ncomps + c);
        const auto *pi_1_27 = buffer.data(pi_1 + 27 * ncomps + c);
        const auto *pi_1_28 = buffer.data(pi_1 + 28 * ncomps + c);
        const auto *pi_1_29 = buffer.data(pi_1 + 29 * ncomps + c);
        const auto *pi_1_30 = buffer.data(pi_1 + 30 * ncomps + c);
        const auto *pi_1_31 = buffer.data(pi_1 + 31 * ncomps + c);
        const auto *pi_1_32 = buffer.data(pi_1 + 32 * ncomps + c);
        const auto *pi_1_33 = buffer.data(pi_1 + 33 * ncomps + c);
        const auto *pi_1_34 = buffer.data(pi_1 + 34 * ncomps + c);
        const auto *pi_1_35 = buffer.data(pi_1 + 35 * ncomps + c);
        const auto *pi_1_36 = buffer.data(pi_1 + 36 * ncomps + c);
        const auto *pi_1_37 = buffer.data(pi_1 + 37 * ncomps + c);
        const auto *pi_1_38 = buffer.data(pi_1 + 38 * ncomps + c);
        const auto *pi_1_39 = buffer.data(pi_1 + 39 * ncomps + c);
        const auto *pi_1_40 = buffer.data(pi_1 + 40 * ncomps + c);
        const auto *pi_1_41 = buffer.data(pi_1 + 41 * ncomps + c);
        const auto *pi_1_42 = buffer.data(pi_1 + 42 * ncomps + c);
        const auto *pi_1_43 = buffer.data(pi_1 + 43 * ncomps + c);
        const auto *pi_1_44 = buffer.data(pi_1 + 44 * ncomps + c);
        const auto *pi_1_45 = buffer.data(pi_1 + 45 * ncomps + c);
        const auto *pi_1_46 = buffer.data(pi_1 + 46 * ncomps + c);
        const auto *pi_1_47 = buffer.data(pi_1 + 47 * ncomps + c);
        const auto *pi_1_48 = buffer.data(pi_1 + 48 * ncomps + c);
        const auto *pi_1_49 = buffer.data(pi_1 + 49 * ncomps + c);
        const auto *pi_1_50 = buffer.data(pi_1 + 50 * ncomps + c);
        const auto *pi_1_51 = buffer.data(pi_1 + 51 * ncomps + c);
        const auto *pi_1_52 = buffer.data(pi_1 + 52 * ncomps + c);
        const auto *pi_1_53 = buffer.data(pi_1 + 53 * ncomps + c);
        const auto *pi_1_54 = buffer.data(pi_1 + 54 * ncomps + c);
        const auto *pi_1_55 = buffer.data(pi_1 + 55 * ncomps + c);
        const auto *pi_1_56 = buffer.data(pi_1 + 56 * ncomps + c);
        const auto *pi_1_57 = buffer.data(pi_1 + 57 * ncomps + c);
        const auto *pi_1_58 = buffer.data(pi_1 + 58 * ncomps + c);
        const auto *pi_1_59 = buffer.data(pi_1 + 59 * ncomps + c);
        const auto *pi_1_60 = buffer.data(pi_1 + 60 * ncomps + c);
        const auto *pi_1_61 = buffer.data(pi_1 + 61 * ncomps + c);
        const auto *pi_1_62 = buffer.data(pi_1 + 62 * ncomps + c);
        const auto *pi_1_63 = buffer.data(pi_1 + 63 * ncomps + c);
        const auto *pi_1_64 = buffer.data(pi_1 + 64 * ncomps + c);
        const auto *pi_1_65 = buffer.data(pi_1 + 65 * ncomps + c);
        const auto *pi_1_66 = buffer.data(pi_1 + 66 * ncomps + c);
        const auto *pi_1_67 = buffer.data(pi_1 + 67 * ncomps + c);
        const auto *pi_1_68 = buffer.data(pi_1 + 68 * ncomps + c);
        const auto *pi_1_69 = buffer.data(pi_1 + 69 * ncomps + c);
        const auto *pi_1_70 = buffer.data(pi_1 + 70 * ncomps + c);
        const auto *pi_1_71 = buffer.data(pi_1 + 71 * ncomps + c);
        const auto *pi_1_72 = buffer.data(pi_1 + 72 * ncomps + c);
        const auto *pi_1_73 = buffer.data(pi_1 + 73 * ncomps + c);
        const auto *pi_1_74 = buffer.data(pi_1 + 74 * ncomps + c);
        const auto *pi_1_75 = buffer.data(pi_1 + 75 * ncomps + c);
        const auto *pi_1_76 = buffer.data(pi_1 + 76 * ncomps + c);
        const auto *pi_1_77 = buffer.data(pi_1 + 77 * ncomps + c);
        const auto *pi_1_78 = buffer.data(pi_1 + 78 * ncomps + c);
        const auto *pi_1_79 = buffer.data(pi_1 + 79 * ncomps + c);
        const auto *pi_1_80 = buffer.data(pi_1 + 80 * ncomps + c);
        const auto *pi_1_81 = buffer.data(pi_1 + 81 * ncomps + c);
        const auto *pi_1_82 = buffer.data(pi_1 + 82 * ncomps + c);
        const auto *pi_1_83 = buffer.data(pi_1 + 83 * ncomps + c);

        const auto *pi_0_28 = buffer.data(pi_0 + 28 * ncomps + c);
        const auto *pi_0_29 = buffer.data(pi_0 + 29 * ncomps + c);
        const auto *pi_0_30 = buffer.data(pi_0 + 30 * ncomps + c);
        const auto *pi_0_31 = buffer.data(pi_0 + 31 * ncomps + c);
        const auto *pi_0_32 = buffer.data(pi_0 + 32 * ncomps + c);
        const auto *pi_0_33 = buffer.data(pi_0 + 33 * ncomps + c);
        const auto *pi_0_34 = buffer.data(pi_0 + 34 * ncomps + c);
        const auto *pi_0_35 = buffer.data(pi_0 + 35 * ncomps + c);
        const auto *pi_0_36 = buffer.data(pi_0 + 36 * ncomps + c);
        const auto *pi_0_37 = buffer.data(pi_0 + 37 * ncomps + c);
        const auto *pi_0_38 = buffer.data(pi_0 + 38 * ncomps + c);
        const auto *pi_0_39 = buffer.data(pi_0 + 39 * ncomps + c);
        const auto *pi_0_40 = buffer.data(pi_0 + 40 * ncomps + c);
        const auto *pi_0_41 = buffer.data(pi_0 + 41 * ncomps + c);
        const auto *pi_0_42 = buffer.data(pi_0 + 42 * ncomps + c);
        const auto *pi_0_43 = buffer.data(pi_0 + 43 * ncomps + c);
        const auto *pi_0_44 = buffer.data(pi_0 + 44 * ncomps + c);
        const auto *pi_0_45 = buffer.data(pi_0 + 45 * ncomps + c);
        const auto *pi_0_46 = buffer.data(pi_0 + 46 * ncomps + c);
        const auto *pi_0_47 = buffer.data(pi_0 + 47 * ncomps + c);
        const auto *pi_0_48 = buffer.data(pi_0 + 48 * ncomps + c);
        const auto *pi_0_49 = buffer.data(pi_0 + 49 * ncomps + c);
        const auto *pi_0_50 = buffer.data(pi_0 + 50 * ncomps + c);
        const auto *pi_0_51 = buffer.data(pi_0 + 51 * ncomps + c);
        const auto *pi_0_52 = buffer.data(pi_0 + 52 * ncomps + c);
        const auto *pi_0_53 = buffer.data(pi_0 + 53 * ncomps + c);
        const auto *pi_0_54 = buffer.data(pi_0 + 54 * ncomps + c);
        const auto *pi_0_55 = buffer.data(pi_0 + 55 * ncomps + c);
        const auto *pi_0_56 = buffer.data(pi_0 + 56 * ncomps + c);
        const auto *pi_0_57 = buffer.data(pi_0 + 57 * ncomps + c);
        const auto *pi_0_58 = buffer.data(pi_0 + 58 * ncomps + c);
        const auto *pi_0_59 = buffer.data(pi_0 + 59 * ncomps + c);
        const auto *pi_0_60 = buffer.data(pi_0 + 60 * ncomps + c);
        const auto *pi_0_61 = buffer.data(pi_0 + 61 * ncomps + c);
        const auto *pi_0_62 = buffer.data(pi_0 + 62 * ncomps + c);
        const auto *pi_0_63 = buffer.data(pi_0 + 63 * ncomps + c);
        const auto *pi_0_64 = buffer.data(pi_0 + 64 * ncomps + c);
        const auto *pi_0_65 = buffer.data(pi_0 + 65 * ncomps + c);
        const auto *pi_0_66 = buffer.data(pi_0 + 66 * ncomps + c);
        const auto *pi_0_67 = buffer.data(pi_0 + 67 * ncomps + c);
        const auto *pi_0_68 = buffer.data(pi_0 + 68 * ncomps + c);
        const auto *pi_0_69 = buffer.data(pi_0 + 69 * ncomps + c);
        const auto *pi_0_70 = buffer.data(pi_0 + 70 * ncomps + c);
        const auto *pi_0_71 = buffer.data(pi_0 + 71 * ncomps + c);
        const auto *pi_0_72 = buffer.data(pi_0 + 72 * ncomps + c);

        const auto *pk_1_0 = buffer.data(pk_1 + 0 * ncomps + c);
        const auto *pk_1_1 = buffer.data(pk_1 + 1 * ncomps + c);
        const auto *pk_1_2 = buffer.data(pk_1 + 2 * ncomps + c);
        const auto *pk_1_3 = buffer.data(pk_1 + 3 * ncomps + c);
        const auto *pk_1_4 = buffer.data(pk_1 + 4 * ncomps + c);
        const auto *pk_1_5 = buffer.data(pk_1 + 5 * ncomps + c);
        const auto *pk_1_6 = buffer.data(pk_1 + 6 * ncomps + c);
        const auto *pk_1_7 = buffer.data(pk_1 + 7 * ncomps + c);
        const auto *pk_1_8 = buffer.data(pk_1 + 8 * ncomps + c);
        const auto *pk_1_9 = buffer.data(pk_1 + 9 * ncomps + c);
        const auto *pk_1_10 = buffer.data(pk_1 + 10 * ncomps + c);
        const auto *pk_1_11 = buffer.data(pk_1 + 11 * ncomps + c);
        const auto *pk_1_12 = buffer.data(pk_1 + 12 * ncomps + c);
        const auto *pk_1_13 = buffer.data(pk_1 + 13 * ncomps + c);
        const auto *pk_1_14 = buffer.data(pk_1 + 14 * ncomps + c);
        const auto *pk_1_15 = buffer.data(pk_1 + 15 * ncomps + c);
        const auto *pk_1_16 = buffer.data(pk_1 + 16 * ncomps + c);
        const auto *pk_1_17 = buffer.data(pk_1 + 17 * ncomps + c);
        const auto *pk_1_18 = buffer.data(pk_1 + 18 * ncomps + c);
        const auto *pk_1_19 = buffer.data(pk_1 + 19 * ncomps + c);
        const auto *pk_1_20 = buffer.data(pk_1 + 20 * ncomps + c);
        const auto *pk_1_21 = buffer.data(pk_1 + 21 * ncomps + c);
        const auto *pk_1_22 = buffer.data(pk_1 + 22 * ncomps + c);
        const auto *pk_1_23 = buffer.data(pk_1 + 23 * ncomps + c);
        const auto *pk_1_24 = buffer.data(pk_1 + 24 * ncomps + c);
        const auto *pk_1_25 = buffer.data(pk_1 + 25 * ncomps + c);
        const auto *pk_1_26 = buffer.data(pk_1 + 26 * ncomps + c);
        const auto *pk_1_27 = buffer.data(pk_1 + 27 * ncomps + c);
        const auto *pk_1_36 = buffer.data(pk_1 + 36 * ncomps + c);
        const auto *pk_1_37 = buffer.data(pk_1 + 37 * ncomps + c);
        const auto *pk_1_38 = buffer.data(pk_1 + 38 * ncomps + c);
        const auto *pk_1_39 = buffer.data(pk_1 + 39 * ncomps + c);
        const auto *pk_1_40 = buffer.data(pk_1 + 40 * ncomps + c);
        const auto *pk_1_41 = buffer.data(pk_1 + 41 * ncomps + c);
        const auto *pk_1_42 = buffer.data(pk_1 + 42 * ncomps + c);
        const auto *pk_1_43 = buffer.data(pk_1 + 43 * ncomps + c);
        const auto *pk_1_44 = buffer.data(pk_1 + 44 * ncomps + c);
        const auto *pk_1_45 = buffer.data(pk_1 + 45 * ncomps + c);
        const auto *pk_1_46 = buffer.data(pk_1 + 46 * ncomps + c);
        const auto *pk_1_47 = buffer.data(pk_1 + 47 * ncomps + c);
        const auto *pk_1_48 = buffer.data(pk_1 + 48 * ncomps + c);
        const auto *pk_1_49 = buffer.data(pk_1 + 49 * ncomps + c);
        const auto *pk_1_50 = buffer.data(pk_1 + 50 * ncomps + c);
        const auto *pk_1_51 = buffer.data(pk_1 + 51 * ncomps + c);
        const auto *pk_1_52 = buffer.data(pk_1 + 52 * ncomps + c);
        const auto *pk_1_53 = buffer.data(pk_1 + 53 * ncomps + c);
        const auto *pk_1_54 = buffer.data(pk_1 + 54 * ncomps + c);
        const auto *pk_1_55 = buffer.data(pk_1 + 55 * ncomps + c);
        const auto *pk_1_56 = buffer.data(pk_1 + 56 * ncomps + c);
        const auto *pk_1_57 = buffer.data(pk_1 + 57 * ncomps + c);
        const auto *pk_1_58 = buffer.data(pk_1 + 58 * ncomps + c);
        const auto *pk_1_59 = buffer.data(pk_1 + 59 * ncomps + c);
        const auto *pk_1_60 = buffer.data(pk_1 + 60 * ncomps + c);
        const auto *pk_1_61 = buffer.data(pk_1 + 61 * ncomps + c);
        const auto *pk_1_62 = buffer.data(pk_1 + 62 * ncomps + c);
        const auto *pk_1_63 = buffer.data(pk_1 + 63 * ncomps + c);
        const auto *pk_1_64 = buffer.data(pk_1 + 64 * ncomps + c);
        const auto *pk_1_65 = buffer.data(pk_1 + 65 * ncomps + c);
        const auto *pk_1_66 = buffer.data(pk_1 + 66 * ncomps + c);
        const auto *pk_1_67 = buffer.data(pk_1 + 67 * ncomps + c);
        const auto *pk_1_68 = buffer.data(pk_1 + 68 * ncomps + c);
        const auto *pk_1_69 = buffer.data(pk_1 + 69 * ncomps + c);
        const auto *pk_1_70 = buffer.data(pk_1 + 70 * ncomps + c);
        const auto *pk_1_72 = buffer.data(pk_1 + 72 * ncomps + c);
        const auto *pk_1_73 = buffer.data(pk_1 + 73 * ncomps + c);
        const auto *pk_1_74 = buffer.data(pk_1 + 74 * ncomps + c);
        const auto *pk_1_75 = buffer.data(pk_1 + 75 * ncomps + c);
        const auto *pk_1_76 = buffer.data(pk_1 + 76 * ncomps + c);
        const auto *pk_1_77 = buffer.data(pk_1 + 77 * ncomps + c);
        const auto *pk_1_78 = buffer.data(pk_1 + 78 * ncomps + c);
        const auto *pk_1_79 = buffer.data(pk_1 + 79 * ncomps + c);
        const auto *pk_1_80 = buffer.data(pk_1 + 80 * ncomps + c);
        const auto *pk_1_81 = buffer.data(pk_1 + 81 * ncomps + c);
        const auto *pk_1_82 = buffer.data(pk_1 + 82 * ncomps + c);
        const auto *pk_1_83 = buffer.data(pk_1 + 83 * ncomps + c);
        const auto *pk_1_84 = buffer.data(pk_1 + 84 * ncomps + c);
        const auto *pk_1_85 = buffer.data(pk_1 + 85 * ncomps + c);
        const auto *pk_1_86 = buffer.data(pk_1 + 86 * ncomps + c);
        const auto *pk_1_87 = buffer.data(pk_1 + 87 * ncomps + c);
        const auto *pk_1_88 = buffer.data(pk_1 + 88 * ncomps + c);
        const auto *pk_1_89 = buffer.data(pk_1 + 89 * ncomps + c);
        const auto *pk_1_90 = buffer.data(pk_1 + 90 * ncomps + c);
        const auto *pk_1_91 = buffer.data(pk_1 + 91 * ncomps + c);
        const auto *pk_1_92 = buffer.data(pk_1 + 92 * ncomps + c);
        const auto *pk_1_93 = buffer.data(pk_1 + 93 * ncomps + c);
        const auto *pk_1_94 = buffer.data(pk_1 + 94 * ncomps + c);
        const auto *pk_1_95 = buffer.data(pk_1 + 95 * ncomps + c);
        const auto *pk_1_96 = buffer.data(pk_1 + 96 * ncomps + c);
        const auto *pk_1_97 = buffer.data(pk_1 + 97 * ncomps + c);
        const auto *pk_1_98 = buffer.data(pk_1 + 98 * ncomps + c);
        const auto *pk_1_99 = buffer.data(pk_1 + 99 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pi_1_0, pi_1_1, pi_1_2, pi_1_3, \
                         pi_1_4, pk_1_0, pk_1_1, pk_1_2, pk_1_3, \
                         pk_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * pi_1_0[k]
                     + pk_1_0[k];

            t_1[k] = -ab_x[k] * pi_1_1[k]
                     + pk_1_1[k];

            t_2[k] = -ab_x[k] * pi_1_2[k]
                     + pk_1_2[k];

            t_3[k] = -ab_x[k] * pi_1_3[k]
                     + pk_1_3[k];

            t_4[k] = -ab_x[k] * pi_1_4[k]
                     + pk_1_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, pi_1_5, pi_1_6, pi_1_7, pi_1_8, \
                         pi_1_9, pk_1_5, pk_1_6, pk_1_7, pk_1_8, \
                         pk_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * pi_1_5[k]
                     + pk_1_5[k];

            t_6[k] = -ab_x[k] * pi_1_6[k]
                     + pk_1_6[k];

            t_7[k] = -ab_x[k] * pi_1_7[k]
                     + pk_1_7[k];

            t_8[k] = -ab_x[k] * pi_1_8[k]
                     + pk_1_8[k];

            t_9[k] = -ab_x[k] * pi_1_9[k]
                     + pk_1_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, pi_1_10, pi_1_11, pi_1_12, \
                         pi_1_13, pi_1_14, pk_1_10, pk_1_11, pk_1_12, pk_1_13, \
                         pk_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * pi_1_10[k]
                      + pk_1_10[k];

            t_11[k] = -ab_x[k] * pi_1_11[k]
                      + pk_1_11[k];

            t_12[k] = -ab_x[k] * pi_1_12[k]
                      + pk_1_12[k];

            t_13[k] = -ab_x[k] * pi_1_13[k]
                      + pk_1_13[k];

            t_14[k] = -ab_x[k] * pi_1_14[k]
                      + pk_1_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, pi_1_15, pi_1_16, pi_1_17, \
                         pi_1_18, pi_1_19, pk_1_15, pk_1_16, pk_1_17, pk_1_18, \
                         pk_1_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * pi_1_15[k]
                      + pk_1_15[k];

            t_16[k] = -ab_x[k] * pi_1_16[k]
                      + pk_1_16[k];

            t_17[k] = -ab_x[k] * pi_1_17[k]
                      + pk_1_17[k];

            t_18[k] = -ab_x[k] * pi_1_18[k]
                      + pk_1_18[k];

            t_19[k] = -ab_x[k] * pi_1_19[k]
                      + pk_1_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, pi_1_20, pi_1_21, pi_1_22, \
                         pi_1_23, pi_1_24, pk_1_20, pk_1_21, pk_1_22, pk_1_23, \
                         pk_1_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * pi_1_20[k]
                      + pk_1_20[k];

            t_21[k] = -ab_x[k] * pi_1_21[k]
                      + pk_1_21[k];

            t_22[k] = -ab_x[k] * pi_1_22[k]
                      + pk_1_22[k];

            t_23[k] = -ab_x[k] * pi_1_23[k]
                      + pk_1_23[k];

            t_24[k] = -ab_x[k] * pi_1_24[k]
                      + pk_1_24[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, pi_1_25, pi_1_26, pi_1_27, \
                         pi_1_28, pi_1_29, pk_1_25, pk_1_26, pk_1_27, pk_1_36, \
                         pk_1_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * pi_1_25[k]
                      + pk_1_25[k];

            t_26[k] = -ab_x[k] * pi_1_26[k]
                      + pk_1_26[k];

            t_27[k] = -ab_x[k] * pi_1_27[k]
                      + pk_1_27[k];

            t_28[k] = -ab_x[k] * pi_1_28[k]
                      + pk_1_36[k];

            t_29[k] = -ab_x[k] * pi_1_29[k]
                      + pk_1_37[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, pi_1_30, pi_1_31, pi_1_32, \
                         pi_1_33, pi_1_34, pk_1_38, pk_1_39, pk_1_40, pk_1_41, \
                         pk_1_42 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * pi_1_30[k]
                      + pk_1_38[k];

            t_31[k] = -ab_x[k] * pi_1_31[k]
                      + pk_1_39[k];

            t_32[k] = -ab_x[k] * pi_1_32[k]
                      + pk_1_40[k];

            t_33[k] = -ab_x[k] * pi_1_33[k]
                      + pk_1_41[k];

            t_34[k] = -ab_x[k] * pi_1_34[k]
                      + pk_1_42[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, pi_1_35, pi_1_36, pi_1_37, \
                         pi_1_38, pi_1_39, pk_1_43, pk_1_44, pk_1_45, pk_1_46, \
                         pk_1_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * pi_1_35[k]
                      + pk_1_43[k];

            t_36[k] = -ab_x[k] * pi_1_36[k]
                      + pk_1_44[k];

            t_37[k] = -ab_x[k] * pi_1_37[k]
                      + pk_1_45[k];

            t_38[k] = -ab_x[k] * pi_1_38[k]
                      + pk_1_46[k];

            t_39[k] = -ab_x[k] * pi_1_39[k]
                      + pk_1_47[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, pi_1_40, pi_1_41, pi_1_42, \
                         pi_1_43, pi_1_44, pk_1_48, pk_1_49, pk_1_50, pk_1_51, \
                         pk_1_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * pi_1_40[k]
                      + pk_1_48[k];

            t_41[k] = -ab_x[k] * pi_1_41[k]
                      + pk_1_49[k];

            t_42[k] = -ab_x[k] * pi_1_42[k]
                      + pk_1_50[k];

            t_43[k] = -ab_x[k] * pi_1_43[k]
                      + pk_1_51[k];

            t_44[k] = -ab_x[k] * pi_1_44[k]
                      + pk_1_52[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, pi_1_45, pi_1_46, pi_1_47, \
                         pi_1_48, pi_1_49, pk_1_53, pk_1_54, pk_1_55, pk_1_56, \
                         pk_1_57 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * pi_1_45[k]
                      + pk_1_53[k];

            t_46[k] = -ab_x[k] * pi_1_46[k]
                      + pk_1_54[k];

            t_47[k] = -ab_x[k] * pi_1_47[k]
                      + pk_1_55[k];

            t_48[k] = -ab_x[k] * pi_1_48[k]
                      + pk_1_56[k];

            t_49[k] = -ab_x[k] * pi_1_49[k]
                      + pk_1_57[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, pi_1_50, pi_1_51, pi_1_52, \
                         pi_1_53, pi_1_54, pk_1_58, pk_1_59, pk_1_60, pk_1_61, \
                         pk_1_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * pi_1_50[k]
                      + pk_1_58[k];

            t_51[k] = -ab_x[k] * pi_1_51[k]
                      + pk_1_59[k];

            t_52[k] = -ab_x[k] * pi_1_52[k]
                      + pk_1_60[k];

            t_53[k] = -ab_x[k] * pi_1_53[k]
                      + pk_1_61[k];

            t_54[k] = -ab_x[k] * pi_1_54[k]
                      + pk_1_62[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, pi_1_55, pi_1_56, pi_1_57, \
                         pi_1_58, pi_1_59, pk_1_63, pk_1_72, pk_1_73, pk_1_74, \
                         pk_1_75 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * pi_1_55[k]
                      + pk_1_63[k];

            t_56[k] = -ab_x[k] * pi_1_56[k]
                      + pk_1_72[k];

            t_57[k] = -ab_x[k] * pi_1_57[k]
                      + pk_1_73[k];

            t_58[k] = -ab_x[k] * pi_1_58[k]
                      + pk_1_74[k];

            t_59[k] = -ab_x[k] * pi_1_59[k]
                      + pk_1_75[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, pi_1_60, pi_1_61, pi_1_62, \
                         pi_1_63, pi_1_64, pk_1_76, pk_1_77, pk_1_78, pk_1_79, \
                         pk_1_80 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * pi_1_60[k]
                      + pk_1_76[k];

            t_61[k] = -ab_x[k] * pi_1_61[k]
                      + pk_1_77[k];

            t_62[k] = -ab_x[k] * pi_1_62[k]
                      + pk_1_78[k];

            t_63[k] = -ab_x[k] * pi_1_63[k]
                      + pk_1_79[k];

            t_64[k] = -ab_x[k] * pi_1_64[k]
                      + pk_1_80[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, pi_1_65, pi_1_66, pi_1_67, \
                         pi_1_68, pi_1_69, pk_1_81, pk_1_82, pk_1_83, pk_1_84, \
                         pk_1_85 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * pi_1_65[k]
                      + pk_1_81[k];

            t_66[k] = -ab_x[k] * pi_1_66[k]
                      + pk_1_82[k];

            t_67[k] = -ab_x[k] * pi_1_67[k]
                      + pk_1_83[k];

            t_68[k] = -ab_x[k] * pi_1_68[k]
                      + pk_1_84[k];

            t_69[k] = -ab_x[k] * pi_1_69[k]
                      + pk_1_85[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, pi_1_70, pi_1_71, pi_1_72, \
                         pi_1_73, pi_1_74, pk_1_86, pk_1_87, pk_1_88, pk_1_89, \
                         pk_1_90 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * pi_1_70[k]
                      + pk_1_86[k];

            t_71[k] = -ab_x[k] * pi_1_71[k]
                      + pk_1_87[k];

            t_72[k] = -ab_x[k] * pi_1_72[k]
                      + pk_1_88[k];

            t_73[k] = -ab_x[k] * pi_1_73[k]
                      + pk_1_89[k];

            t_74[k] = -ab_x[k] * pi_1_74[k]
                      + pk_1_90[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, pi_1_75, pi_1_76, pi_1_77, \
                         pi_1_78, pi_1_79, pk_1_91, pk_1_92, pk_1_93, pk_1_94, \
                         pk_1_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * pi_1_75[k]
                      + pk_1_91[k];

            t_76[k] = -ab_x[k] * pi_1_76[k]
                      + pk_1_92[k];

            t_77[k] = -ab_x[k] * pi_1_77[k]
                      + pk_1_93[k];

            t_78[k] = -ab_x[k] * pi_1_78[k]
                      + pk_1_94[k];

            t_79[k] = -ab_x[k] * pi_1_79[k]
                      + pk_1_95[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, ab_x, pi_1_80, pi_1_81, pi_1_82, pi_1_83, \
                         pk_1_96, pk_1_97, pk_1_98, pk_1_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * pi_1_80[k]
                      + pk_1_96[k];

            t_81[k] = -ab_x[k] * pi_1_81[k]
                      + pk_1_97[k];

            t_82[k] = -ab_x[k] * pi_1_82[k]
                      + pk_1_98[k];

            t_83[k] = -ab_x[k] * pi_1_83[k]
                      + pk_1_99[k];
        }

#pragma omp simd aligned(t_84, t_85, t_86, ab_y, pi_1_28, pi_1_29, pi_1_30, pi_0_28, pi_0_29, \
                         pi_0_30, pk_1_37, pk_1_39, pk_1_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_84[k] = -ab_y[k] * pi_1_28[k]
                      + pi_0_28[k]
                      + pk_1_37[k];

            t_85[k] = -ab_y[k] * pi_1_29[k]
                      + pi_0_29[k]
                      + pk_1_39[k];

            t_86[k] = -ab_y[k] * pi_1_30[k]
                      + pi_0_30[k]
                      + pk_1_40[k];
        }

#pragma omp simd aligned(t_87, t_88, t_89, ab_y, pi_1_31, pi_1_32, pi_1_33, pi_0_31, pi_0_32, \
                         pi_0_33, pk_1_42, pk_1_43, pk_1_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_87[k] = -ab_y[k] * pi_1_31[k]
                      + pi_0_31[k]
                      + pk_1_42[k];

            t_88[k] = -ab_y[k] * pi_1_32[k]
                      + pi_0_32[k]
                      + pk_1_43[k];

            t_89[k] = -ab_y[k] * pi_1_33[k]
                      + pi_0_33[k]
                      + pk_1_44[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, ab_y, pi_1_34, pi_1_35, pi_1_36, pi_0_34, pi_0_35, \
                         pi_0_36, pk_1_46, pk_1_47, pk_1_48 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_y[k] * pi_1_34[k]
                      + pi_0_34[k]
                      + pk_1_46[k];

            t_91[k] = -ab_y[k] * pi_1_35[k]
                      + pi_0_35[k]
                      + pk_1_47[k];

            t_92[k] = -ab_y[k] * pi_1_36[k]
                      + pi_0_36[k]
                      + pk_1_48[k];
        }

#pragma omp simd aligned(t_93, t_94, t_95, ab_y, pi_1_37, pi_1_38, pi_1_39, pi_0_37, pi_0_38, \
                         pi_0_39, pk_1_49, pk_1_51, pk_1_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_93[k] = -ab_y[k] * pi_1_37[k]
                      + pi_0_37[k]
                      + pk_1_49[k];

            t_94[k] = -ab_y[k] * pi_1_38[k]
                      + pi_0_38[k]
                      + pk_1_51[k];

            t_95[k] = -ab_y[k] * pi_1_39[k]
                      + pi_0_39[k]
                      + pk_1_52[k];
        }

#pragma omp simd aligned(t_96, t_97, t_98, ab_y, pi_1_40, pi_1_41, pi_1_42, pi_0_40, pi_0_41, \
                         pi_0_42, pk_1_53, pk_1_54, pk_1_55 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_96[k] = -ab_y[k] * pi_1_40[k]
                      + pi_0_40[k]
                      + pk_1_53[k];

            t_97[k] = -ab_y[k] * pi_1_41[k]
                      + pi_0_41[k]
                      + pk_1_54[k];

            t_98[k] = -ab_y[k] * pi_1_42[k]
                      + pi_0_42[k]
                      + pk_1_55[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, ab_y, pi_1_43, pi_1_44, pi_1_45, pi_0_43, \
                         pi_0_44, pi_0_45, pk_1_57, pk_1_58, pk_1_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = -ab_y[k] * pi_1_43[k]
                      + pi_0_43[k]
                      + pk_1_57[k];

            t_100[k] = -ab_y[k] * pi_1_44[k]
                       + pi_0_44[k]
                       + pk_1_58[k];

            t_101[k] = -ab_y[k] * pi_1_45[k]
                       + pi_0_45[k]
                       + pk_1_59[k];
        }

#pragma omp simd aligned(t_102, t_103, t_104, ab_y, pi_1_46, pi_1_47, pi_1_48, pi_0_46, \
                         pi_0_47, pi_0_48, pk_1_60, pk_1_61, pk_1_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_102[k] = -ab_y[k] * pi_1_46[k]
                       + pi_0_46[k]
                       + pk_1_60[k];

            t_103[k] = -ab_y[k] * pi_1_47[k]
                       + pi_0_47[k]
                       + pk_1_61[k];

            t_104[k] = -ab_y[k] * pi_1_48[k]
                       + pi_0_48[k]
                       + pk_1_62[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, ab_y, pi_1_49, pi_1_50, pi_1_51, pi_0_49, \
                         pi_0_50, pi_0_51, pk_1_64, pk_1_65, pk_1_66 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_y[k] * pi_1_49[k]
                       + pi_0_49[k]
                       + pk_1_64[k];

            t_106[k] = -ab_y[k] * pi_1_50[k]
                       + pi_0_50[k]
                       + pk_1_65[k];

            t_107[k] = -ab_y[k] * pi_1_51[k]
                       + pi_0_51[k]
                       + pk_1_66[k];
        }

#pragma omp simd aligned(t_108, t_109, t_110, ab_y, pi_1_52, pi_1_53, pi_1_54, pi_0_52, \
                         pi_0_53, pi_0_54, pk_1_67, pk_1_68, pk_1_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_108[k] = -ab_y[k] * pi_1_52[k]
                       + pi_0_52[k]
                       + pk_1_67[k];

            t_109[k] = -ab_y[k] * pi_1_53[k]
                       + pi_0_53[k]
                       + pk_1_68[k];

            t_110[k] = -ab_y[k] * pi_1_54[k]
                       + pi_0_54[k]
                       + pk_1_69[k];
        }

#pragma omp simd aligned(t_111, t_112, t_113, ab_y, pi_1_55, pi_1_56, pi_1_57, pi_0_55, \
                         pi_0_56, pi_0_57, pk_1_70, pk_1_73, pk_1_75 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_111[k] = -ab_y[k] * pi_1_55[k]
                       + pi_0_55[k]
                       + pk_1_70[k];

            t_112[k] = -ab_y[k] * pi_1_56[k]
                       + pi_0_56[k]
                       + pk_1_73[k];

            t_113[k] = -ab_y[k] * pi_1_57[k]
                       + pi_0_57[k]
                       + pk_1_75[k];
        }

#pragma omp simd aligned(t_114, t_115, t_116, ab_y, pi_1_58, pi_1_59, pi_1_60, pi_0_58, \
                         pi_0_59, pi_0_60, pk_1_76, pk_1_78, pk_1_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_114[k] = -ab_y[k] * pi_1_58[k]
                       + pi_0_58[k]
                       + pk_1_76[k];

            t_115[k] = -ab_y[k] * pi_1_59[k]
                       + pi_0_59[k]
                       + pk_1_78[k];

            t_116[k] = -ab_y[k] * pi_1_60[k]
                       + pi_0_60[k]
                       + pk_1_79[k];
        }

#pragma omp simd aligned(t_117, t_118, t_119, ab_y, pi_1_61, pi_1_62, pi_1_63, pi_0_61, \
                         pi_0_62, pi_0_63, pk_1_80, pk_1_82, pk_1_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_117[k] = -ab_y[k] * pi_1_61[k]
                       + pi_0_61[k]
                       + pk_1_80[k];

            t_118[k] = -ab_y[k] * pi_1_62[k]
                       + pi_0_62[k]
                       + pk_1_82[k];

            t_119[k] = -ab_y[k] * pi_1_63[k]
                       + pi_0_63[k]
                       + pk_1_83[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, ab_y, pi_1_64, pi_1_65, pi_1_66, pi_0_64, \
                         pi_0_65, pi_0_66, pk_1_84, pk_1_85, pk_1_87 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_y[k] * pi_1_64[k]
                       + pi_0_64[k]
                       + pk_1_84[k];

            t_121[k] = -ab_y[k] * pi_1_65[k]
                       + pi_0_65[k]
                       + pk_1_85[k];

            t_122[k] = -ab_y[k] * pi_1_66[k]
                       + pi_0_66[k]
                       + pk_1_87[k];
        }

#pragma omp simd aligned(t_123, t_124, t_125, ab_y, pi_1_67, pi_1_68, pi_1_69, pi_0_67, \
                         pi_0_68, pi_0_69, pk_1_88, pk_1_89, pk_1_90 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_123[k] = -ab_y[k] * pi_1_67[k]
                       + pi_0_67[k]
                       + pk_1_88[k];

            t_124[k] = -ab_y[k] * pi_1_68[k]
                       + pi_0_68[k]
                       + pk_1_89[k];

            t_125[k] = -ab_y[k] * pi_1_69[k]
                       + pi_0_69[k]
                       + pk_1_90[k];
        }

#pragma omp simd aligned(t_126, t_127, t_128, ab_y, pi_1_70, pi_1_71, pi_1_72, pi_0_70, \
                         pi_0_71, pi_0_72, pk_1_91, pk_1_93, pk_1_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_126[k] = -ab_y[k] * pi_1_70[k]
                       + pi_0_70[k]
                       + pk_1_91[k];

            t_127[k] = -ab_y[k] * pi_1_71[k]
                       + pi_0_71[k]
                       + pk_1_93[k];

            t_128[k] = -ab_y[k] * pi_1_72[k]
                       + pi_0_72[k]
                       + pk_1_94[k];
        }
    }
}

static auto
compute_hrr_geom_010y_di_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                const size_t target, const size_t pi_1, const size_t pi_0,
                                const size_t pk_1, const size_t ncomps,
                                const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
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

        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *pi_1_56 = buffer.data(pi_1 + 56 * ncomps + c);
        const auto *pi_1_57 = buffer.data(pi_1 + 57 * ncomps + c);
        const auto *pi_1_58 = buffer.data(pi_1 + 58 * ncomps + c);
        const auto *pi_1_59 = buffer.data(pi_1 + 59 * ncomps + c);
        const auto *pi_1_60 = buffer.data(pi_1 + 60 * ncomps + c);
        const auto *pi_1_61 = buffer.data(pi_1 + 61 * ncomps + c);
        const auto *pi_1_62 = buffer.data(pi_1 + 62 * ncomps + c);
        const auto *pi_1_63 = buffer.data(pi_1 + 63 * ncomps + c);
        const auto *pi_1_64 = buffer.data(pi_1 + 64 * ncomps + c);
        const auto *pi_1_65 = buffer.data(pi_1 + 65 * ncomps + c);
        const auto *pi_1_66 = buffer.data(pi_1 + 66 * ncomps + c);
        const auto *pi_1_67 = buffer.data(pi_1 + 67 * ncomps + c);
        const auto *pi_1_68 = buffer.data(pi_1 + 68 * ncomps + c);
        const auto *pi_1_69 = buffer.data(pi_1 + 69 * ncomps + c);
        const auto *pi_1_70 = buffer.data(pi_1 + 70 * ncomps + c);
        const auto *pi_1_71 = buffer.data(pi_1 + 71 * ncomps + c);
        const auto *pi_1_72 = buffer.data(pi_1 + 72 * ncomps + c);
        const auto *pi_1_73 = buffer.data(pi_1 + 73 * ncomps + c);
        const auto *pi_1_74 = buffer.data(pi_1 + 74 * ncomps + c);
        const auto *pi_1_75 = buffer.data(pi_1 + 75 * ncomps + c);
        const auto *pi_1_76 = buffer.data(pi_1 + 76 * ncomps + c);
        const auto *pi_1_77 = buffer.data(pi_1 + 77 * ncomps + c);
        const auto *pi_1_78 = buffer.data(pi_1 + 78 * ncomps + c);
        const auto *pi_1_79 = buffer.data(pi_1 + 79 * ncomps + c);
        const auto *pi_1_80 = buffer.data(pi_1 + 80 * ncomps + c);
        const auto *pi_1_81 = buffer.data(pi_1 + 81 * ncomps + c);
        const auto *pi_1_82 = buffer.data(pi_1 + 82 * ncomps + c);
        const auto *pi_1_83 = buffer.data(pi_1 + 83 * ncomps + c);

        const auto *pi_0_73 = buffer.data(pi_0 + 73 * ncomps + c);
        const auto *pi_0_74 = buffer.data(pi_0 + 74 * ncomps + c);
        const auto *pi_0_75 = buffer.data(pi_0 + 75 * ncomps + c);
        const auto *pi_0_76 = buffer.data(pi_0 + 76 * ncomps + c);
        const auto *pi_0_77 = buffer.data(pi_0 + 77 * ncomps + c);
        const auto *pi_0_78 = buffer.data(pi_0 + 78 * ncomps + c);
        const auto *pi_0_79 = buffer.data(pi_0 + 79 * ncomps + c);
        const auto *pi_0_80 = buffer.data(pi_0 + 80 * ncomps + c);
        const auto *pi_0_81 = buffer.data(pi_0 + 81 * ncomps + c);
        const auto *pi_0_82 = buffer.data(pi_0 + 82 * ncomps + c);
        const auto *pi_0_83 = buffer.data(pi_0 + 83 * ncomps + c);

        const auto *pk_1_74 = buffer.data(pk_1 + 74 * ncomps + c);
        const auto *pk_1_76 = buffer.data(pk_1 + 76 * ncomps + c);
        const auto *pk_1_77 = buffer.data(pk_1 + 77 * ncomps + c);
        const auto *pk_1_79 = buffer.data(pk_1 + 79 * ncomps + c);
        const auto *pk_1_80 = buffer.data(pk_1 + 80 * ncomps + c);
        const auto *pk_1_81 = buffer.data(pk_1 + 81 * ncomps + c);
        const auto *pk_1_83 = buffer.data(pk_1 + 83 * ncomps + c);
        const auto *pk_1_84 = buffer.data(pk_1 + 84 * ncomps + c);
        const auto *pk_1_85 = buffer.data(pk_1 + 85 * ncomps + c);
        const auto *pk_1_86 = buffer.data(pk_1 + 86 * ncomps + c);
        const auto *pk_1_88 = buffer.data(pk_1 + 88 * ncomps + c);
        const auto *pk_1_89 = buffer.data(pk_1 + 89 * ncomps + c);
        const auto *pk_1_90 = buffer.data(pk_1 + 90 * ncomps + c);
        const auto *pk_1_91 = buffer.data(pk_1 + 91 * ncomps + c);
        const auto *pk_1_92 = buffer.data(pk_1 + 92 * ncomps + c);
        const auto *pk_1_94 = buffer.data(pk_1 + 94 * ncomps + c);
        const auto *pk_1_95 = buffer.data(pk_1 + 95 * ncomps + c);
        const auto *pk_1_96 = buffer.data(pk_1 + 96 * ncomps + c);
        const auto *pk_1_97 = buffer.data(pk_1 + 97 * ncomps + c);
        const auto *pk_1_98 = buffer.data(pk_1 + 98 * ncomps + c);
        const auto *pk_1_99 = buffer.data(pk_1 + 99 * ncomps + c);
        const auto *pk_1_100 = buffer.data(pk_1 + 100 * ncomps + c);
        const auto *pk_1_101 = buffer.data(pk_1 + 101 * ncomps + c);
        const auto *pk_1_102 = buffer.data(pk_1 + 102 * ncomps + c);
        const auto *pk_1_103 = buffer.data(pk_1 + 103 * ncomps + c);
        const auto *pk_1_104 = buffer.data(pk_1 + 104 * ncomps + c);
        const auto *pk_1_105 = buffer.data(pk_1 + 105 * ncomps + c);
        const auto *pk_1_106 = buffer.data(pk_1 + 106 * ncomps + c);
        const auto *pk_1_107 = buffer.data(pk_1 + 107 * ncomps + c);

#pragma omp simd aligned(t_129, t_130, t_131, ab_y, pi_1_73, pi_1_74, pi_1_75, pi_0_73, \
                         pi_0_74, pi_0_75, pk_1_95, pk_1_96, pk_1_97 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_129[k] = -ab_y[k] * pi_1_73[k]
                       + pi_0_73[k]
                       + pk_1_95[k];

            t_130[k] = -ab_y[k] * pi_1_74[k]
                       + pi_0_74[k]
                       + pk_1_96[k];

            t_131[k] = -ab_y[k] * pi_1_75[k]
                       + pi_0_75[k]
                       + pk_1_97[k];
        }

#pragma omp simd aligned(t_132, t_133, t_134, ab_y, pi_1_76, pi_1_77, pi_1_78, pi_0_76, \
                         pi_0_77, pi_0_78, pk_1_98, pk_1_100, \
                         pk_1_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_132[k] = -ab_y[k] * pi_1_76[k]
                       + pi_0_76[k]
                       + pk_1_98[k];

            t_133[k] = -ab_y[k] * pi_1_77[k]
                       + pi_0_77[k]
                       + pk_1_100[k];

            t_134[k] = -ab_y[k] * pi_1_78[k]
                       + pi_0_78[k]
                       + pk_1_101[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, ab_y, pi_1_79, pi_1_80, pi_1_81, pi_0_79, \
                         pi_0_80, pi_0_81, pk_1_102, pk_1_103, \
                         pk_1_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_y[k] * pi_1_79[k]
                       + pi_0_79[k]
                       + pk_1_102[k];

            t_136[k] = -ab_y[k] * pi_1_80[k]
                       + pi_0_80[k]
                       + pk_1_103[k];

            t_137[k] = -ab_y[k] * pi_1_81[k]
                       + pi_0_81[k]
                       + pk_1_104[k];
        }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, ab_y, ab_z, pi_1_56, pi_1_57, pi_1_82, \
                         pi_1_83, pi_0_82, pi_0_83, pk_1_74, pk_1_76, pk_1_105, \
                         pk_1_106 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_138[k] = -ab_y[k] * pi_1_82[k]
                       + pi_0_82[k]
                       + pk_1_105[k];

            t_139[k] = -ab_y[k] * pi_1_83[k]
                       + pi_0_83[k]
                       + pk_1_106[k];

            t_140[k] = -ab_z[k] * pi_1_56[k]
                       + pk_1_74[k];

            t_141[k] = -ab_z[k] * pi_1_57[k]
                       + pk_1_76[k];
        }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, ab_z, pi_1_58, pi_1_59, pi_1_60, \
                         pi_1_61, pi_1_62, pk_1_77, pk_1_79, pk_1_80, pk_1_81, \
                         pk_1_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_142[k] = -ab_z[k] * pi_1_58[k]
                       + pk_1_77[k];

            t_143[k] = -ab_z[k] * pi_1_59[k]
                       + pk_1_79[k];

            t_144[k] = -ab_z[k] * pi_1_60[k]
                       + pk_1_80[k];

            t_145[k] = -ab_z[k] * pi_1_61[k]
                       + pk_1_81[k];

            t_146[k] = -ab_z[k] * pi_1_62[k]
                       + pk_1_83[k];
        }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, ab_z, pi_1_63, pi_1_64, pi_1_65, \
                         pi_1_66, pi_1_67, pk_1_84, pk_1_85, pk_1_86, pk_1_88, \
                         pk_1_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_147[k] = -ab_z[k] * pi_1_63[k]
                       + pk_1_84[k];

            t_148[k] = -ab_z[k] * pi_1_64[k]
                       + pk_1_85[k];

            t_149[k] = -ab_z[k] * pi_1_65[k]
                       + pk_1_86[k];

            t_150[k] = -ab_z[k] * pi_1_66[k]
                       + pk_1_88[k];

            t_151[k] = -ab_z[k] * pi_1_67[k]
                       + pk_1_89[k];
        }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, ab_z, pi_1_68, pi_1_69, pi_1_70, \
                         pi_1_71, pi_1_72, pk_1_90, pk_1_91, pk_1_92, pk_1_94, \
                         pk_1_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_152[k] = -ab_z[k] * pi_1_68[k]
                       + pk_1_90[k];

            t_153[k] = -ab_z[k] * pi_1_69[k]
                       + pk_1_91[k];

            t_154[k] = -ab_z[k] * pi_1_70[k]
                       + pk_1_92[k];

            t_155[k] = -ab_z[k] * pi_1_71[k]
                       + pk_1_94[k];

            t_156[k] = -ab_z[k] * pi_1_72[k]
                       + pk_1_95[k];
        }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, ab_z, pi_1_73, pi_1_74, pi_1_75, \
                         pi_1_76, pi_1_77, pk_1_96, pk_1_97, pk_1_98, pk_1_99, \
                         pk_1_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_157[k] = -ab_z[k] * pi_1_73[k]
                       + pk_1_96[k];

            t_158[k] = -ab_z[k] * pi_1_74[k]
                       + pk_1_97[k];

            t_159[k] = -ab_z[k] * pi_1_75[k]
                       + pk_1_98[k];

            t_160[k] = -ab_z[k] * pi_1_76[k]
                       + pk_1_99[k];

            t_161[k] = -ab_z[k] * pi_1_77[k]
                       + pk_1_101[k];
        }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, ab_z, pi_1_78, pi_1_79, pi_1_80, \
                         pi_1_81, pi_1_82, pk_1_102, pk_1_103, pk_1_104, pk_1_105, \
                         pk_1_106 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_162[k] = -ab_z[k] * pi_1_78[k]
                       + pk_1_102[k];

            t_163[k] = -ab_z[k] * pi_1_79[k]
                       + pk_1_103[k];

            t_164[k] = -ab_z[k] * pi_1_80[k]
                       + pk_1_104[k];

            t_165[k] = -ab_z[k] * pi_1_81[k]
                       + pk_1_105[k];

            t_166[k] = -ab_z[k] * pi_1_82[k]
                       + pk_1_106[k];
        }

#pragma omp simd aligned(t_167, ab_z, pi_1_83, pk_1_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_167[k] = -ab_z[k] * pi_1_83[k]
                       + pk_1_107[k];
        }
    }
}

auto
compute_hrr_geom_010y_di(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                         const size_t target, const size_t pi_1, const size_t pi_0,
                         const size_t pk_1, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_geom_010y_di_piece0(buffer, coordinates, target, pi_1, pi_0, pk_1, ncomps,
                                    nmax);

    compute_hrr_geom_010y_di_piece1(buffer, coordinates, target, pi_1, pi_0, pk_1, ncomps,
                                    nmax);
}

}  // namespace simdtrf
