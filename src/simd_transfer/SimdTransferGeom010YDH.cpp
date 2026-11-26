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


#include "SimdTransferGeom010YDH.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_010y_dh(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                         const size_t target, const size_t ph_1, const size_t ph_0,
                         const size_t pi_1, const size_t ncomps, const size_t nmax) -> void
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ph_1_0 = buffer.data(ph_1 + 0 * ncomps + c);
        const auto *ph_1_1 = buffer.data(ph_1 + 1 * ncomps + c);
        const auto *ph_1_2 = buffer.data(ph_1 + 2 * ncomps + c);
        const auto *ph_1_3 = buffer.data(ph_1 + 3 * ncomps + c);
        const auto *ph_1_4 = buffer.data(ph_1 + 4 * ncomps + c);
        const auto *ph_1_5 = buffer.data(ph_1 + 5 * ncomps + c);
        const auto *ph_1_6 = buffer.data(ph_1 + 6 * ncomps + c);
        const auto *ph_1_7 = buffer.data(ph_1 + 7 * ncomps + c);
        const auto *ph_1_8 = buffer.data(ph_1 + 8 * ncomps + c);
        const auto *ph_1_9 = buffer.data(ph_1 + 9 * ncomps + c);
        const auto *ph_1_10 = buffer.data(ph_1 + 10 * ncomps + c);
        const auto *ph_1_11 = buffer.data(ph_1 + 11 * ncomps + c);
        const auto *ph_1_12 = buffer.data(ph_1 + 12 * ncomps + c);
        const auto *ph_1_13 = buffer.data(ph_1 + 13 * ncomps + c);
        const auto *ph_1_14 = buffer.data(ph_1 + 14 * ncomps + c);
        const auto *ph_1_15 = buffer.data(ph_1 + 15 * ncomps + c);
        const auto *ph_1_16 = buffer.data(ph_1 + 16 * ncomps + c);
        const auto *ph_1_17 = buffer.data(ph_1 + 17 * ncomps + c);
        const auto *ph_1_18 = buffer.data(ph_1 + 18 * ncomps + c);
        const auto *ph_1_19 = buffer.data(ph_1 + 19 * ncomps + c);
        const auto *ph_1_20 = buffer.data(ph_1 + 20 * ncomps + c);
        const auto *ph_1_21 = buffer.data(ph_1 + 21 * ncomps + c);
        const auto *ph_1_22 = buffer.data(ph_1 + 22 * ncomps + c);
        const auto *ph_1_23 = buffer.data(ph_1 + 23 * ncomps + c);
        const auto *ph_1_24 = buffer.data(ph_1 + 24 * ncomps + c);
        const auto *ph_1_25 = buffer.data(ph_1 + 25 * ncomps + c);
        const auto *ph_1_26 = buffer.data(ph_1 + 26 * ncomps + c);
        const auto *ph_1_27 = buffer.data(ph_1 + 27 * ncomps + c);
        const auto *ph_1_28 = buffer.data(ph_1 + 28 * ncomps + c);
        const auto *ph_1_29 = buffer.data(ph_1 + 29 * ncomps + c);
        const auto *ph_1_30 = buffer.data(ph_1 + 30 * ncomps + c);
        const auto *ph_1_31 = buffer.data(ph_1 + 31 * ncomps + c);
        const auto *ph_1_32 = buffer.data(ph_1 + 32 * ncomps + c);
        const auto *ph_1_33 = buffer.data(ph_1 + 33 * ncomps + c);
        const auto *ph_1_34 = buffer.data(ph_1 + 34 * ncomps + c);
        const auto *ph_1_35 = buffer.data(ph_1 + 35 * ncomps + c);
        const auto *ph_1_36 = buffer.data(ph_1 + 36 * ncomps + c);
        const auto *ph_1_37 = buffer.data(ph_1 + 37 * ncomps + c);
        const auto *ph_1_38 = buffer.data(ph_1 + 38 * ncomps + c);
        const auto *ph_1_39 = buffer.data(ph_1 + 39 * ncomps + c);
        const auto *ph_1_40 = buffer.data(ph_1 + 40 * ncomps + c);
        const auto *ph_1_41 = buffer.data(ph_1 + 41 * ncomps + c);
        const auto *ph_1_42 = buffer.data(ph_1 + 42 * ncomps + c);
        const auto *ph_1_43 = buffer.data(ph_1 + 43 * ncomps + c);
        const auto *ph_1_44 = buffer.data(ph_1 + 44 * ncomps + c);
        const auto *ph_1_45 = buffer.data(ph_1 + 45 * ncomps + c);
        const auto *ph_1_46 = buffer.data(ph_1 + 46 * ncomps + c);
        const auto *ph_1_47 = buffer.data(ph_1 + 47 * ncomps + c);
        const auto *ph_1_48 = buffer.data(ph_1 + 48 * ncomps + c);
        const auto *ph_1_49 = buffer.data(ph_1 + 49 * ncomps + c);
        const auto *ph_1_50 = buffer.data(ph_1 + 50 * ncomps + c);
        const auto *ph_1_51 = buffer.data(ph_1 + 51 * ncomps + c);
        const auto *ph_1_52 = buffer.data(ph_1 + 52 * ncomps + c);
        const auto *ph_1_53 = buffer.data(ph_1 + 53 * ncomps + c);
        const auto *ph_1_54 = buffer.data(ph_1 + 54 * ncomps + c);
        const auto *ph_1_55 = buffer.data(ph_1 + 55 * ncomps + c);
        const auto *ph_1_56 = buffer.data(ph_1 + 56 * ncomps + c);
        const auto *ph_1_57 = buffer.data(ph_1 + 57 * ncomps + c);
        const auto *ph_1_58 = buffer.data(ph_1 + 58 * ncomps + c);
        const auto *ph_1_59 = buffer.data(ph_1 + 59 * ncomps + c);
        const auto *ph_1_60 = buffer.data(ph_1 + 60 * ncomps + c);
        const auto *ph_1_61 = buffer.data(ph_1 + 61 * ncomps + c);
        const auto *ph_1_62 = buffer.data(ph_1 + 62 * ncomps + c);

        const auto *ph_0_21 = buffer.data(ph_0 + 21 * ncomps + c);
        const auto *ph_0_22 = buffer.data(ph_0 + 22 * ncomps + c);
        const auto *ph_0_23 = buffer.data(ph_0 + 23 * ncomps + c);
        const auto *ph_0_24 = buffer.data(ph_0 + 24 * ncomps + c);
        const auto *ph_0_25 = buffer.data(ph_0 + 25 * ncomps + c);
        const auto *ph_0_26 = buffer.data(ph_0 + 26 * ncomps + c);
        const auto *ph_0_27 = buffer.data(ph_0 + 27 * ncomps + c);
        const auto *ph_0_28 = buffer.data(ph_0 + 28 * ncomps + c);
        const auto *ph_0_29 = buffer.data(ph_0 + 29 * ncomps + c);
        const auto *ph_0_30 = buffer.data(ph_0 + 30 * ncomps + c);
        const auto *ph_0_31 = buffer.data(ph_0 + 31 * ncomps + c);
        const auto *ph_0_32 = buffer.data(ph_0 + 32 * ncomps + c);
        const auto *ph_0_33 = buffer.data(ph_0 + 33 * ncomps + c);
        const auto *ph_0_34 = buffer.data(ph_0 + 34 * ncomps + c);
        const auto *ph_0_35 = buffer.data(ph_0 + 35 * ncomps + c);
        const auto *ph_0_36 = buffer.data(ph_0 + 36 * ncomps + c);
        const auto *ph_0_37 = buffer.data(ph_0 + 37 * ncomps + c);
        const auto *ph_0_38 = buffer.data(ph_0 + 38 * ncomps + c);
        const auto *ph_0_39 = buffer.data(ph_0 + 39 * ncomps + c);
        const auto *ph_0_40 = buffer.data(ph_0 + 40 * ncomps + c);
        const auto *ph_0_41 = buffer.data(ph_0 + 41 * ncomps + c);
        const auto *ph_0_42 = buffer.data(ph_0 + 42 * ncomps + c);
        const auto *ph_0_43 = buffer.data(ph_0 + 43 * ncomps + c);
        const auto *ph_0_44 = buffer.data(ph_0 + 44 * ncomps + c);
        const auto *ph_0_45 = buffer.data(ph_0 + 45 * ncomps + c);
        const auto *ph_0_46 = buffer.data(ph_0 + 46 * ncomps + c);
        const auto *ph_0_47 = buffer.data(ph_0 + 47 * ncomps + c);
        const auto *ph_0_48 = buffer.data(ph_0 + 48 * ncomps + c);
        const auto *ph_0_49 = buffer.data(ph_0 + 49 * ncomps + c);
        const auto *ph_0_50 = buffer.data(ph_0 + 50 * ncomps + c);
        const auto *ph_0_51 = buffer.data(ph_0 + 51 * ncomps + c);
        const auto *ph_0_52 = buffer.data(ph_0 + 52 * ncomps + c);
        const auto *ph_0_53 = buffer.data(ph_0 + 53 * ncomps + c);
        const auto *ph_0_54 = buffer.data(ph_0 + 54 * ncomps + c);
        const auto *ph_0_55 = buffer.data(ph_0 + 55 * ncomps + c);
        const auto *ph_0_56 = buffer.data(ph_0 + 56 * ncomps + c);
        const auto *ph_0_57 = buffer.data(ph_0 + 57 * ncomps + c);
        const auto *ph_0_58 = buffer.data(ph_0 + 58 * ncomps + c);
        const auto *ph_0_59 = buffer.data(ph_0 + 59 * ncomps + c);
        const auto *ph_0_60 = buffer.data(ph_0 + 60 * ncomps + c);
        const auto *ph_0_61 = buffer.data(ph_0 + 61 * ncomps + c);
        const auto *ph_0_62 = buffer.data(ph_0 + 62 * ncomps + c);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ph_1_0, ph_1_1, ph_1_2, ph_1_3, \
                         ph_1_4, pi_1_0, pi_1_1, pi_1_2, pi_1_3, \
                         pi_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * ph_1_0[k]
                     + pi_1_0[k];

            t_1[k] = -ab_x[k] * ph_1_1[k]
                     + pi_1_1[k];

            t_2[k] = -ab_x[k] * ph_1_2[k]
                     + pi_1_2[k];

            t_3[k] = -ab_x[k] * ph_1_3[k]
                     + pi_1_3[k];

            t_4[k] = -ab_x[k] * ph_1_4[k]
                     + pi_1_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ph_1_5, ph_1_6, ph_1_7, ph_1_8, \
                         ph_1_9, pi_1_5, pi_1_6, pi_1_7, pi_1_8, \
                         pi_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * ph_1_5[k]
                     + pi_1_5[k];

            t_6[k] = -ab_x[k] * ph_1_6[k]
                     + pi_1_6[k];

            t_7[k] = -ab_x[k] * ph_1_7[k]
                     + pi_1_7[k];

            t_8[k] = -ab_x[k] * ph_1_8[k]
                     + pi_1_8[k];

            t_9[k] = -ab_x[k] * ph_1_9[k]
                     + pi_1_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, ph_1_10, ph_1_11, ph_1_12, \
                         ph_1_13, ph_1_14, pi_1_10, pi_1_11, pi_1_12, pi_1_13, \
                         pi_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * ph_1_10[k]
                      + pi_1_10[k];

            t_11[k] = -ab_x[k] * ph_1_11[k]
                      + pi_1_11[k];

            t_12[k] = -ab_x[k] * ph_1_12[k]
                      + pi_1_12[k];

            t_13[k] = -ab_x[k] * ph_1_13[k]
                      + pi_1_13[k];

            t_14[k] = -ab_x[k] * ph_1_14[k]
                      + pi_1_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ph_1_15, ph_1_16, ph_1_17, \
                         ph_1_18, ph_1_19, pi_1_15, pi_1_16, pi_1_17, pi_1_18, \
                         pi_1_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * ph_1_15[k]
                      + pi_1_15[k];

            t_16[k] = -ab_x[k] * ph_1_16[k]
                      + pi_1_16[k];

            t_17[k] = -ab_x[k] * ph_1_17[k]
                      + pi_1_17[k];

            t_18[k] = -ab_x[k] * ph_1_18[k]
                      + pi_1_18[k];

            t_19[k] = -ab_x[k] * ph_1_19[k]
                      + pi_1_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, ph_1_20, ph_1_21, ph_1_22, \
                         ph_1_23, ph_1_24, pi_1_20, pi_1_28, pi_1_29, pi_1_30, \
                         pi_1_31 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * ph_1_20[k]
                      + pi_1_20[k];

            t_21[k] = -ab_x[k] * ph_1_21[k]
                      + pi_1_28[k];

            t_22[k] = -ab_x[k] * ph_1_22[k]
                      + pi_1_29[k];

            t_23[k] = -ab_x[k] * ph_1_23[k]
                      + pi_1_30[k];

            t_24[k] = -ab_x[k] * ph_1_24[k]
                      + pi_1_31[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ph_1_25, ph_1_26, ph_1_27, \
                         ph_1_28, ph_1_29, pi_1_32, pi_1_33, pi_1_34, pi_1_35, \
                         pi_1_36 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * ph_1_25[k]
                      + pi_1_32[k];

            t_26[k] = -ab_x[k] * ph_1_26[k]
                      + pi_1_33[k];

            t_27[k] = -ab_x[k] * ph_1_27[k]
                      + pi_1_34[k];

            t_28[k] = -ab_x[k] * ph_1_28[k]
                      + pi_1_35[k];

            t_29[k] = -ab_x[k] * ph_1_29[k]
                      + pi_1_36[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, ph_1_30, ph_1_31, ph_1_32, \
                         ph_1_33, ph_1_34, pi_1_37, pi_1_38, pi_1_39, pi_1_40, \
                         pi_1_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * ph_1_30[k]
                      + pi_1_37[k];

            t_31[k] = -ab_x[k] * ph_1_31[k]
                      + pi_1_38[k];

            t_32[k] = -ab_x[k] * ph_1_32[k]
                      + pi_1_39[k];

            t_33[k] = -ab_x[k] * ph_1_33[k]
                      + pi_1_40[k];

            t_34[k] = -ab_x[k] * ph_1_34[k]
                      + pi_1_41[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, ph_1_35, ph_1_36, ph_1_37, \
                         ph_1_38, ph_1_39, pi_1_42, pi_1_43, pi_1_44, pi_1_45, \
                         pi_1_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * ph_1_35[k]
                      + pi_1_42[k];

            t_36[k] = -ab_x[k] * ph_1_36[k]
                      + pi_1_43[k];

            t_37[k] = -ab_x[k] * ph_1_37[k]
                      + pi_1_44[k];

            t_38[k] = -ab_x[k] * ph_1_38[k]
                      + pi_1_45[k];

            t_39[k] = -ab_x[k] * ph_1_39[k]
                      + pi_1_46[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, ph_1_40, ph_1_41, ph_1_42, \
                         ph_1_43, ph_1_44, pi_1_47, pi_1_48, pi_1_56, pi_1_57, \
                         pi_1_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * ph_1_40[k]
                      + pi_1_47[k];

            t_41[k] = -ab_x[k] * ph_1_41[k]
                      + pi_1_48[k];

            t_42[k] = -ab_x[k] * ph_1_42[k]
                      + pi_1_56[k];

            t_43[k] = -ab_x[k] * ph_1_43[k]
                      + pi_1_57[k];

            t_44[k] = -ab_x[k] * ph_1_44[k]
                      + pi_1_58[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ph_1_45, ph_1_46, ph_1_47, \
                         ph_1_48, ph_1_49, pi_1_59, pi_1_60, pi_1_61, pi_1_62, \
                         pi_1_63 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * ph_1_45[k]
                      + pi_1_59[k];

            t_46[k] = -ab_x[k] * ph_1_46[k]
                      + pi_1_60[k];

            t_47[k] = -ab_x[k] * ph_1_47[k]
                      + pi_1_61[k];

            t_48[k] = -ab_x[k] * ph_1_48[k]
                      + pi_1_62[k];

            t_49[k] = -ab_x[k] * ph_1_49[k]
                      + pi_1_63[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, ph_1_50, ph_1_51, ph_1_52, \
                         ph_1_53, ph_1_54, pi_1_64, pi_1_65, pi_1_66, pi_1_67, \
                         pi_1_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * ph_1_50[k]
                      + pi_1_64[k];

            t_51[k] = -ab_x[k] * ph_1_51[k]
                      + pi_1_65[k];

            t_52[k] = -ab_x[k] * ph_1_52[k]
                      + pi_1_66[k];

            t_53[k] = -ab_x[k] * ph_1_53[k]
                      + pi_1_67[k];

            t_54[k] = -ab_x[k] * ph_1_54[k]
                      + pi_1_68[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, ph_1_55, ph_1_56, ph_1_57, \
                         ph_1_58, ph_1_59, pi_1_69, pi_1_70, pi_1_71, pi_1_72, \
                         pi_1_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * ph_1_55[k]
                      + pi_1_69[k];

            t_56[k] = -ab_x[k] * ph_1_56[k]
                      + pi_1_70[k];

            t_57[k] = -ab_x[k] * ph_1_57[k]
                      + pi_1_71[k];

            t_58[k] = -ab_x[k] * ph_1_58[k]
                      + pi_1_72[k];

            t_59[k] = -ab_x[k] * ph_1_59[k]
                      + pi_1_73[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, ab_x, ab_y, ph_1_21, ph_1_60, ph_1_61, \
                         ph_1_62, ph_0_21, pi_1_29, pi_1_74, pi_1_75, \
                         pi_1_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * ph_1_60[k]
                      + pi_1_74[k];

            t_61[k] = -ab_x[k] * ph_1_61[k]
                      + pi_1_75[k];

            t_62[k] = -ab_x[k] * ph_1_62[k]
                      + pi_1_76[k];

            t_63[k] = -ab_y[k] * ph_1_21[k]
                      + ph_0_21[k]
                      + pi_1_29[k];
        }

#pragma omp simd aligned(t_64, t_65, t_66, ab_y, ph_1_22, ph_1_23, ph_1_24, ph_0_22, ph_0_23, \
                         ph_0_24, pi_1_31, pi_1_32, pi_1_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_64[k] = -ab_y[k] * ph_1_22[k]
                      + ph_0_22[k]
                      + pi_1_31[k];

            t_65[k] = -ab_y[k] * ph_1_23[k]
                      + ph_0_23[k]
                      + pi_1_32[k];

            t_66[k] = -ab_y[k] * ph_1_24[k]
                      + ph_0_24[k]
                      + pi_1_34[k];
        }

#pragma omp simd aligned(t_67, t_68, t_69, ab_y, ph_1_25, ph_1_26, ph_1_27, ph_0_25, ph_0_26, \
                         ph_0_27, pi_1_35, pi_1_36, pi_1_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_67[k] = -ab_y[k] * ph_1_25[k]
                      + ph_0_25[k]
                      + pi_1_35[k];

            t_68[k] = -ab_y[k] * ph_1_26[k]
                      + ph_0_26[k]
                      + pi_1_36[k];

            t_69[k] = -ab_y[k] * ph_1_27[k]
                      + ph_0_27[k]
                      + pi_1_38[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, ab_y, ph_1_28, ph_1_29, ph_1_30, ph_0_28, ph_0_29, \
                         ph_0_30, pi_1_39, pi_1_40, pi_1_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_y[k] * ph_1_28[k]
                      + ph_0_28[k]
                      + pi_1_39[k];

            t_71[k] = -ab_y[k] * ph_1_29[k]
                      + ph_0_29[k]
                      + pi_1_40[k];

            t_72[k] = -ab_y[k] * ph_1_30[k]
                      + ph_0_30[k]
                      + pi_1_41[k];
        }

#pragma omp simd aligned(t_73, t_74, t_75, ab_y, ph_1_31, ph_1_32, ph_1_33, ph_0_31, ph_0_32, \
                         ph_0_33, pi_1_43, pi_1_44, pi_1_45 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_73[k] = -ab_y[k] * ph_1_31[k]
                      + ph_0_31[k]
                      + pi_1_43[k];

            t_74[k] = -ab_y[k] * ph_1_32[k]
                      + ph_0_32[k]
                      + pi_1_44[k];

            t_75[k] = -ab_y[k] * ph_1_33[k]
                      + ph_0_33[k]
                      + pi_1_45[k];
        }

#pragma omp simd aligned(t_76, t_77, t_78, ab_y, ph_1_34, ph_1_35, ph_1_36, ph_0_34, ph_0_35, \
                         ph_0_36, pi_1_46, pi_1_47, pi_1_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_76[k] = -ab_y[k] * ph_1_34[k]
                      + ph_0_34[k]
                      + pi_1_46[k];

            t_77[k] = -ab_y[k] * ph_1_35[k]
                      + ph_0_35[k]
                      + pi_1_47[k];

            t_78[k] = -ab_y[k] * ph_1_36[k]
                      + ph_0_36[k]
                      + pi_1_49[k];
        }

#pragma omp simd aligned(t_79, t_80, t_81, ab_y, ph_1_37, ph_1_38, ph_1_39, ph_0_37, ph_0_38, \
                         ph_0_39, pi_1_50, pi_1_51, pi_1_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_79[k] = -ab_y[k] * ph_1_37[k]
                      + ph_0_37[k]
                      + pi_1_50[k];

            t_80[k] = -ab_y[k] * ph_1_38[k]
                      + ph_0_38[k]
                      + pi_1_51[k];

            t_81[k] = -ab_y[k] * ph_1_39[k]
                      + ph_0_39[k]
                      + pi_1_52[k];
        }

#pragma omp simd aligned(t_82, t_83, t_84, ab_y, ph_1_40, ph_1_41, ph_1_42, ph_0_40, ph_0_41, \
                         ph_0_42, pi_1_53, pi_1_54, pi_1_57 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_82[k] = -ab_y[k] * ph_1_40[k]
                      + ph_0_40[k]
                      + pi_1_53[k];

            t_83[k] = -ab_y[k] * ph_1_41[k]
                      + ph_0_41[k]
                      + pi_1_54[k];

            t_84[k] = -ab_y[k] * ph_1_42[k]
                      + ph_0_42[k]
                      + pi_1_57[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, ab_y, ph_1_43, ph_1_44, ph_1_45, ph_0_43, ph_0_44, \
                         ph_0_45, pi_1_59, pi_1_60, pi_1_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_y[k] * ph_1_43[k]
                      + ph_0_43[k]
                      + pi_1_59[k];

            t_86[k] = -ab_y[k] * ph_1_44[k]
                      + ph_0_44[k]
                      + pi_1_60[k];

            t_87[k] = -ab_y[k] * ph_1_45[k]
                      + ph_0_45[k]
                      + pi_1_62[k];
        }

#pragma omp simd aligned(t_88, t_89, t_90, ab_y, ph_1_46, ph_1_47, ph_1_48, ph_0_46, ph_0_47, \
                         ph_0_48, pi_1_63, pi_1_64, pi_1_66 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_88[k] = -ab_y[k] * ph_1_46[k]
                      + ph_0_46[k]
                      + pi_1_63[k];

            t_89[k] = -ab_y[k] * ph_1_47[k]
                      + ph_0_47[k]
                      + pi_1_64[k];

            t_90[k] = -ab_y[k] * ph_1_48[k]
                      + ph_0_48[k]
                      + pi_1_66[k];
        }

#pragma omp simd aligned(t_91, t_92, t_93, ab_y, ph_1_49, ph_1_50, ph_1_51, ph_0_49, ph_0_50, \
                         ph_0_51, pi_1_67, pi_1_68, pi_1_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_91[k] = -ab_y[k] * ph_1_49[k]
                      + ph_0_49[k]
                      + pi_1_67[k];

            t_92[k] = -ab_y[k] * ph_1_50[k]
                      + ph_0_50[k]
                      + pi_1_68[k];

            t_93[k] = -ab_y[k] * ph_1_51[k]
                      + ph_0_51[k]
                      + pi_1_69[k];
        }

#pragma omp simd aligned(t_94, t_95, t_96, ab_y, ph_1_52, ph_1_53, ph_1_54, ph_0_52, ph_0_53, \
                         ph_0_54, pi_1_71, pi_1_72, pi_1_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_94[k] = -ab_y[k] * ph_1_52[k]
                      + ph_0_52[k]
                      + pi_1_71[k];

            t_95[k] = -ab_y[k] * ph_1_53[k]
                      + ph_0_53[k]
                      + pi_1_72[k];

            t_96[k] = -ab_y[k] * ph_1_54[k]
                      + ph_0_54[k]
                      + pi_1_73[k];
        }

#pragma omp simd aligned(t_97, t_98, t_99, ab_y, ph_1_55, ph_1_56, ph_1_57, ph_0_55, ph_0_56, \
                         ph_0_57, pi_1_74, pi_1_75, pi_1_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_97[k] = -ab_y[k] * ph_1_55[k]
                      + ph_0_55[k]
                      + pi_1_74[k];

            t_98[k] = -ab_y[k] * ph_1_56[k]
                      + ph_0_56[k]
                      + pi_1_75[k];

            t_99[k] = -ab_y[k] * ph_1_57[k]
                      + ph_0_57[k]
                      + pi_1_77[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, ab_y, ph_1_58, ph_1_59, ph_1_60, ph_0_58, \
                         ph_0_59, ph_0_60, pi_1_78, pi_1_79, pi_1_80 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_y[k] * ph_1_58[k]
                       + ph_0_58[k]
                       + pi_1_78[k];

            t_101[k] = -ab_y[k] * ph_1_59[k]
                       + ph_0_59[k]
                       + pi_1_79[k];

            t_102[k] = -ab_y[k] * ph_1_60[k]
                       + ph_0_60[k]
                       + pi_1_80[k];
        }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, ab_y, ab_z, ph_1_42, ph_1_43, ph_1_61, \
                         ph_1_62, ph_0_61, ph_0_62, pi_1_58, pi_1_60, pi_1_81, \
                         pi_1_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_103[k] = -ab_y[k] * ph_1_61[k]
                       + ph_0_61[k]
                       + pi_1_81[k];

            t_104[k] = -ab_y[k] * ph_1_62[k]
                       + ph_0_62[k]
                       + pi_1_82[k];

            t_105[k] = -ab_z[k] * ph_1_42[k]
                       + pi_1_58[k];

            t_106[k] = -ab_z[k] * ph_1_43[k]
                       + pi_1_60[k];
        }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, ab_z, ph_1_44, ph_1_45, ph_1_46, \
                         ph_1_47, ph_1_48, pi_1_61, pi_1_63, pi_1_64, pi_1_65, \
                         pi_1_67 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_107[k] = -ab_z[k] * ph_1_44[k]
                       + pi_1_61[k];

            t_108[k] = -ab_z[k] * ph_1_45[k]
                       + pi_1_63[k];

            t_109[k] = -ab_z[k] * ph_1_46[k]
                       + pi_1_64[k];

            t_110[k] = -ab_z[k] * ph_1_47[k]
                       + pi_1_65[k];

            t_111[k] = -ab_z[k] * ph_1_48[k]
                       + pi_1_67[k];
        }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, ab_z, ph_1_49, ph_1_50, ph_1_51, \
                         ph_1_52, ph_1_53, pi_1_68, pi_1_69, pi_1_70, pi_1_72, \
                         pi_1_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_112[k] = -ab_z[k] * ph_1_49[k]
                       + pi_1_68[k];

            t_113[k] = -ab_z[k] * ph_1_50[k]
                       + pi_1_69[k];

            t_114[k] = -ab_z[k] * ph_1_51[k]
                       + pi_1_70[k];

            t_115[k] = -ab_z[k] * ph_1_52[k]
                       + pi_1_72[k];

            t_116[k] = -ab_z[k] * ph_1_53[k]
                       + pi_1_73[k];
        }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, ab_z, ph_1_54, ph_1_55, ph_1_56, \
                         ph_1_57, ph_1_58, pi_1_74, pi_1_75, pi_1_76, pi_1_78, \
                         pi_1_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_117[k] = -ab_z[k] * ph_1_54[k]
                       + pi_1_74[k];

            t_118[k] = -ab_z[k] * ph_1_55[k]
                       + pi_1_75[k];

            t_119[k] = -ab_z[k] * ph_1_56[k]
                       + pi_1_76[k];

            t_120[k] = -ab_z[k] * ph_1_57[k]
                       + pi_1_78[k];

            t_121[k] = -ab_z[k] * ph_1_58[k]
                       + pi_1_79[k];
        }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, ab_z, ph_1_59, ph_1_60, ph_1_61, ph_1_62, \
                         pi_1_80, pi_1_81, pi_1_82, pi_1_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_122[k] = -ab_z[k] * ph_1_59[k]
                       + pi_1_80[k];

            t_123[k] = -ab_z[k] * ph_1_60[k]
                       + pi_1_81[k];

            t_124[k] = -ab_z[k] * ph_1_61[k]
                       + pi_1_82[k];

            t_125[k] = -ab_z[k] * ph_1_62[k]
                       + pi_1_83[k];
        }
    }
}

}  // namespace simdtrf
