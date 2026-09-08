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


#include "SimdTransferFG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_fg_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t dg, const size_t dh,
                   const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.875 * std::sqrt(14.0);
    const auto f_1 = 0.625 * std::sqrt(14.0);
    const auto f_2 = 5.625 * std::sqrt(7.0);
    const auto f_3 = 1.875 * std::sqrt(7.0);
    const auto f_4 = 0.625 * std::sqrt(7.0);
    const auto f_5 = 1.875 * std::sqrt(2.0);
    const auto f_6 = 11.25 * std::sqrt(2.0);
    const auto f_7 = 0.625 * std::sqrt(2.0);
    const auto f_8 = 3.75 * std::sqrt(2.0);
    const auto f_9 = 0.28125 * std::sqrt(10.0);
    const auto f_10 = 0.5625 * std::sqrt(10.0);
    const auto f_11 = 2.25 * std::sqrt(10.0);
    const auto f_12 = 0.75 * std::sqrt(10.0);
    const auto f_13 = 0.09375 * std::sqrt(10.0);
    const auto f_14 = 0.1875 * std::sqrt(10.0);
    const auto f_15 = 0.25 * std::sqrt(10.0);
    const auto f_16 = 0.9375 * std::sqrt(2.0);
    const auto f_17 = 5.625 * std::sqrt(2.0);
    const auto f_18 = 0.3125 * std::sqrt(2.0);
    const auto f_19 = 0.46875 * std::sqrt(14.0);
    const auto f_20 = 2.8125 * std::sqrt(14.0);
    const auto f_21 = 0.15625 * std::sqrt(14.0);
    const auto f_22 = 0.9375 * std::sqrt(14.0);
    const auto f_23 = 2.5 * std::sqrt(21.0);
    const auto f_24 = 3.75 * std::sqrt(42.0);
    const auto f_25 = 1.25 * std::sqrt(42.0);
    const auto f_26 = 2.5 * std::sqrt(3.0);
    const auto f_27 = 15.0 * std::sqrt(3.0);
    const auto f_28 = 3.75 * std::sqrt(6.0);
    const auto f_29 = 5.0 * std::sqrt(6.0);
    const auto f_30 = 0.375 * std::sqrt(15.0);
    const auto f_31 = 0.75 * std::sqrt(15.0);
    const auto f_32 = 3.0 * std::sqrt(15.0);
    const auto f_33 = std::sqrt(15.0);
    const auto f_34 = 1.25 * std::sqrt(3.0);
    const auto f_35 = 7.5 * std::sqrt(3.0);
    const auto f_36 = 0.625 * std::sqrt(21.0);
    const auto f_37 = 3.75 * std::sqrt(21.0);
    const auto f_38 = 0.125 * std::sqrt(210.0);
    const auto f_39 = 0.5 * std::sqrt(210.0);
    const auto f_40 = 0.375 * std::sqrt(105.0);
    const auto f_41 = 0.125 * std::sqrt(105.0);
    const auto f_42 = 1.5 * std::sqrt(105.0);
    const auto f_43 = 0.5 * std::sqrt(105.0);
    const auto f_44 = 0.125 * std::sqrt(30.0);
    const auto f_45 = 0.75 * std::sqrt(30.0);
    const auto f_46 = 0.5 * std::sqrt(30.0);
    const auto f_47 = 3.0 * std::sqrt(30.0);
    const auto f_48 = 0.5 * std::sqrt(15.0);
    const auto f_49 = 1.5 * std::sqrt(15.0);
    const auto f_50 = 2.0 * std::sqrt(15.0);
    const auto f_51 = 0.09375 * std::sqrt(6.0);
    const auto f_52 = 0.1875 * std::sqrt(6.0);
    const auto f_53 = 0.75 * std::sqrt(6.0);
    const auto f_54 = 0.25 * std::sqrt(6.0);
    const auto f_55 = 0.375 * std::sqrt(6.0);
    const auto f_56 = 3.0 * std::sqrt(6.0);
    const auto f_57 = std::sqrt(6.0);
    const auto f_58 = 0.0625 * std::sqrt(30.0);
    const auto f_59 = 0.375 * std::sqrt(30.0);
    const auto f_60 = 0.25 * std::sqrt(30.0);
    const auto f_61 = 1.5 * std::sqrt(30.0);
    const auto f_62 = 0.03125 * std::sqrt(210.0);
    const auto f_63 = 0.1875 * std::sqrt(210.0);
    const auto f_64 = 0.75 * std::sqrt(210.0);
    const auto f_65 = 0.75 * std::sqrt(35.0);
    const auto f_66 = 0.5 * std::sqrt(35.0);
    const auto f_67 = 1.125 * std::sqrt(70.0);
    const auto f_68 = 0.375 * std::sqrt(70.0);
    const auto f_69 = 0.75 * std::sqrt(70.0);
    const auto f_70 = 0.25 * std::sqrt(70.0);
    const auto f_71 = 0.75 * std::sqrt(5.0);
    const auto f_72 = 4.5 * std::sqrt(5.0);
    const auto f_73 = 0.5 * std::sqrt(5.0);
    const auto f_74 = 3.0 * std::sqrt(5.0);
    const auto f_75 = 1.125 * std::sqrt(10.0);
    const auto f_76 = 1.5 * std::sqrt(10.0);
    const auto f_77 = std::sqrt(10.0);
    const auto f_78 = 0.375 * std::sqrt(5.0);
    const auto f_79 = 2.25 * std::sqrt(5.0);
    const auto f_80 = 0.25 * std::sqrt(5.0);
    const auto f_81 = 1.5 * std::sqrt(5.0);
    const auto f_82 = 0.1875 * std::sqrt(35.0);
    const auto f_83 = 1.125 * std::sqrt(35.0);
    const auto f_84 = 0.125 * std::sqrt(35.0);
    const auto f_85 = 1.25 * std::sqrt(21.0);
    const auto f_86 = 1.875 * std::sqrt(42.0);
    const auto f_87 = 0.625 * std::sqrt(42.0);
    const auto f_88 = 1.875 * std::sqrt(6.0);
    const auto f_89 = 2.5 * std::sqrt(6.0);
    const auto f_90 = 0.1875 * std::sqrt(15.0);
    const auto f_91 = 0.625 * std::sqrt(3.0);
    const auto f_92 = 3.75 * std::sqrt(3.0);
    const auto f_93 = 0.3125 * std::sqrt(21.0);
    const auto f_94 = 1.875 * std::sqrt(21.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;
    auto *g_9 = values + 9 * nvalues;
    auto *g_10 = values + 10 * nvalues;
    auto *g_11 = values + 11 * nvalues;
    auto *g_12 = values + 12 * nvalues;
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;
    auto *g_25 = values + 25 * nvalues;
    auto *g_26 = values + 26 * nvalues;
    auto *g_27 = values + 27 * nvalues;
    auto *g_28 = values + 28 * nvalues;
    auto *g_29 = values + 29 * nvalues;
    auto *g_30 = values + 30 * nvalues;
    auto *g_31 = values + 31 * nvalues;
    auto *g_32 = values + 32 * nvalues;
    auto *g_33 = values + 33 * nvalues;
    auto *g_34 = values + 34 * nvalues;
    auto *g_35 = values + 35 * nvalues;
    auto *g_36 = values + 36 * nvalues;
    auto *g_37 = values + 37 * nvalues;
    auto *g_38 = values + 38 * nvalues;
    auto *g_39 = values + 39 * nvalues;
    auto *g_40 = values + 40 * nvalues;
    auto *g_41 = values + 41 * nvalues;
    auto *g_42 = values + 42 * nvalues;
    auto *g_43 = values + 43 * nvalues;
    auto *g_44 = values + 44 * nvalues;
    auto *g_45 = values + 45 * nvalues;
    auto *g_46 = values + 46 * nvalues;
    auto *g_47 = values + 47 * nvalues;
    auto *g_48 = values + 48 * nvalues;
    auto *g_49 = values + 49 * nvalues;
    auto *g_50 = values + 50 * nvalues;
    auto *g_51 = values + 51 * nvalues;
    auto *g_52 = values + 52 * nvalues;
    auto *g_53 = values + 53 * nvalues;
    auto *g_54 = values + 54 * nvalues;
    auto *g_55 = values + 55 * nvalues;
    auto *g_56 = values + 56 * nvalues;
    auto *g_57 = values + 57 * nvalues;
    auto *g_58 = values + 58 * nvalues;
    auto *g_59 = values + 59 * nvalues;
    auto *g_60 = values + 60 * nvalues;
    auto *g_61 = values + 61 * nvalues;
    auto *g_62 = values + 62 * nvalues;

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_43 = buffer.data(dg + 43);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_47 = buffer.data(dg + 47);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_52 = buffer.data(dg + 52);
    const auto *dg_53 = buffer.data(dg + 53);
    const auto *dg_54 = buffer.data(dg + 54);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_56 = buffer.data(dg + 56);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_58 = buffer.data(dg + 58);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_60 = buffer.data(dg + 60);
    const auto *dg_61 = buffer.data(dg + 61);
    const auto *dg_62 = buffer.data(dg + 62);
    const auto *dg_63 = buffer.data(dg + 63);
    const auto *dg_64 = buffer.data(dg + 64);
    const auto *dg_65 = buffer.data(dg + 65);
    const auto *dg_66 = buffer.data(dg + 66);
    const auto *dg_67 = buffer.data(dg + 67);
    const auto *dg_68 = buffer.data(dg + 68);
    const auto *dg_69 = buffer.data(dg + 69);
    const auto *dg_70 = buffer.data(dg + 70);
    const auto *dg_71 = buffer.data(dg + 71);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_73 = buffer.data(dg + 73);
    const auto *dg_74 = buffer.data(dg + 74);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_76 = buffer.data(dg + 76);
    const auto *dg_77 = buffer.data(dg + 77);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_79 = buffer.data(dg + 79);
    const auto *dg_80 = buffer.data(dg + 80);
    const auto *dg_81 = buffer.data(dg + 81);
    const auto *dg_82 = buffer.data(dg + 82);
    const auto *dg_83 = buffer.data(dg + 83);
    const auto *dg_84 = buffer.data(dg + 84);
    const auto *dg_85 = buffer.data(dg + 85);
    const auto *dg_86 = buffer.data(dg + 86);
    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_88 = buffer.data(dg + 88);
    const auto *dg_89 = buffer.data(dg + 89);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_55 = buffer.data(dh + 55);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_64 = buffer.data(dh + 64);
    const auto *dh_65 = buffer.data(dh + 65);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_67 = buffer.data(dh + 67);
    const auto *dh_68 = buffer.data(dh + 68);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_70 = buffer.data(dh + 70);
    const auto *dh_71 = buffer.data(dh + 71);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_73 = buffer.data(dh + 73);
    const auto *dh_74 = buffer.data(dh + 74);
    const auto *dh_75 = buffer.data(dh + 75);
    const auto *dh_76 = buffer.data(dh + 76);
    const auto *dh_77 = buffer.data(dh + 77);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_79 = buffer.data(dh + 79);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_82 = buffer.data(dh + 82);
    const auto *dh_84 = buffer.data(dh + 84);
    const auto *dh_85 = buffer.data(dh + 85);
    const auto *dh_86 = buffer.data(dh + 86);
    const auto *dh_87 = buffer.data(dh + 87);
    const auto *dh_88 = buffer.data(dh + 88);
    const auto *dh_89 = buffer.data(dh + 89);
    const auto *dh_90 = buffer.data(dh + 90);
    const auto *dh_91 = buffer.data(dh + 91);
    const auto *dh_92 = buffer.data(dh + 92);
    const auto *dh_93 = buffer.data(dh + 93);
    const auto *dh_94 = buffer.data(dh + 94);
    const auto *dh_95 = buffer.data(dh + 95);
    const auto *dh_96 = buffer.data(dh + 96);
    const auto *dh_97 = buffer.data(dh + 97);
    const auto *dh_98 = buffer.data(dh + 98);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_106 = buffer.data(dh + 106);
    const auto *dh_107 = buffer.data(dh + 107);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_109 = buffer.data(dh + 109);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_111 = buffer.data(dh + 111);
    const auto *dh_112 = buffer.data(dh + 112);
    const auto *dh_113 = buffer.data(dh + 113);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_115 = buffer.data(dh + 115);
    const auto *dh_116 = buffer.data(dh + 116);
    const auto *dh_117 = buffer.data(dh + 117);
    const auto *dh_118 = buffer.data(dh + 118);
    const auto *dh_119 = buffer.data(dh + 119);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_121 = buffer.data(dh + 121);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_124 = buffer.data(dh + 124);
    const auto *dh_125 = buffer.data(dh + 125);

#pragma omp simd aligned(ab_x, ab_y, dg_16, dg_21, dg_46, dg_51, dh_22, dh_27, dh_66, \
                         dh_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -f_0 * ab_x[k] * dg_16[k]
                 + f_0 * ab_x[k] * dg_21[k]
                 + f_1 * ab_y[k] * dg_46[k]
                 - f_1 * ab_y[k] * dg_51[k]
                 + f_0 * dh_22[k]
                 - f_0 * dh_27[k]
                 - f_1 * dh_66[k]
                 + f_1 * dh_73[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_19, dg_26, dg_49, dg_56, dh_25, dh_32, dh_70, \
                         dh_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = -f_2 * ab_x[k] * dg_19[k]
                 + f_3 * ab_x[k] * dg_26[k]
                 + f_3 * ab_y[k] * dg_49[k]
                 - f_4 * ab_y[k] * dg_56[k]
                 + f_2 * dh_25[k]
                 - f_3 * dh_32[k]
                 - f_3 * dh_70[k]
                 + f_4 * dh_79[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_16, dg_21, dg_23, dg_46, dg_51, dg_53, dh_22, dh_27, \
                         dh_29, dh_66, dh_73, dh_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = f_5 * ab_x[k] * dg_16[k]
                 + f_5 * ab_x[k] * dg_21[k]
                 - f_6 * ab_x[k] * dg_23[k]
                 - f_7 * ab_y[k] * dg_46[k]
                 - f_7 * ab_y[k] * dg_51[k]
                 + f_8 * ab_y[k] * dg_53[k]
                 - f_5 * dh_22[k]
                 - f_5 * dh_27[k]
                 + f_6 * dh_29[k]
                 + f_7 * dh_66[k]
                 + f_7 * dh_73[k]
                 - f_8 * dh_75[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_19, dg_26, dg_28, dg_49, dg_56, dg_58, dh_25, dh_32, \
                         dh_34, dh_70, dh_79, dh_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = 5.625 * ab_x[k] * dg_19[k]
                 + 5.625 * ab_x[k] * dg_26[k]
                 - 7.5 * ab_x[k] * dg_28[k]
                 - 1.875 * ab_y[k] * dg_49[k]
                 - 1.875 * ab_y[k] * dg_56[k]
                 + 2.5 * ab_y[k] * dg_58[k]
                 - 5.625 * dh_25[k]
                 - 5.625 * dh_32[k]
                 + 7.5 * dh_34[k]
                 + 1.875 * dh_70[k]
                 + 1.875 * dh_79[k]
                 - 2.5 * dh_81[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_15, dg_18, dg_20, dg_25, dg_27, dg_29, dg_45, dg_48, \
                         dg_50, dg_55, dg_57, dg_59, dh_21, dh_24, dh_26, dh_31, dh_33, dh_35, \
                         dh_64, dh_69, dh_71, dh_78, dh_80, dh_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_9 * ab_x[k] * dg_15[k]
                 - f_10 * ab_x[k] * dg_18[k]
                 + f_11 * ab_x[k] * dg_20[k]
                 - f_9 * ab_x[k] * dg_25[k]
                 + f_11 * ab_x[k] * dg_27[k]
                 - f_12 * ab_x[k] * dg_29[k]
                 + f_13 * ab_y[k] * dg_45[k]
                 + f_14 * ab_y[k] * dg_48[k]
                 - f_12 * ab_y[k] * dg_50[k]
                 + f_13 * ab_y[k] * dg_55[k]
                 - f_12 * ab_y[k] * dg_57[k]
                 + f_15 * ab_y[k] * dg_59[k]
                 + f_9 * dh_21[k]
                 + f_10 * dh_24[k]
                 - f_11 * dh_26[k]
                 + f_9 * dh_31[k]
                 - f_11 * dh_33[k]
                 + f_12 * dh_35[k]
                 - f_13 * dh_64[k]
                 - f_14 * dh_69[k]
                 + f_12 * dh_71[k]
                 - f_13 * dh_78[k]
                 + f_12 * dh_80[k]
                 - f_15 * dh_82[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_17, dg_22, dg_24, dg_47, dg_52, dg_54, dh_23, dh_28, \
                         dh_30, dh_67, dh_74, dh_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = 5.625 * ab_x[k] * dg_17[k]
                 + 5.625 * ab_x[k] * dg_22[k]
                 - 7.5 * ab_x[k] * dg_24[k]
                 - 1.875 * ab_y[k] * dg_47[k]
                 - 1.875 * ab_y[k] * dg_52[k]
                 + 2.5 * ab_y[k] * dg_54[k]
                 - 5.625 * dh_23[k]
                 - 5.625 * dh_28[k]
                 + 7.5 * dh_30[k]
                 + 1.875 * dh_67[k]
                 + 1.875 * dh_74[k]
                 - 2.5 * dh_76[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_15, dg_20, dg_25, dg_27, dg_45, dg_50, dg_55, dg_57, \
                         dh_21, dh_26, dh_31, dh_33, dh_64, dh_71, dh_78, \
                         dh_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_16 * ab_x[k] * dg_15[k]
                 - f_17 * ab_x[k] * dg_20[k]
                 - f_16 * ab_x[k] * dg_25[k]
                 + f_17 * ab_x[k] * dg_27[k]
                 - f_18 * ab_y[k] * dg_45[k]
                 + f_5 * ab_y[k] * dg_50[k]
                 + f_18 * ab_y[k] * dg_55[k]
                 - f_5 * ab_y[k] * dg_57[k]
                 - f_16 * dh_21[k]
                 + f_17 * dh_26[k]
                 + f_16 * dh_31[k]
                 - f_17 * dh_33[k]
                 + f_18 * dh_64[k]
                 - f_5 * dh_71[k]
                 - f_18 * dh_78[k]
                 + f_5 * dh_80[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_17, dg_22, dg_47, dg_52, dh_23, dh_28, dh_67, \
                         dh_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_3 * ab_x[k] * dg_17[k]
                 + f_2 * ab_x[k] * dg_22[k]
                 + f_4 * ab_y[k] * dg_47[k]
                 - f_3 * ab_y[k] * dg_52[k]
                 + f_3 * dh_23[k]
                 - f_2 * dh_28[k]
                 - f_4 * dh_67[k]
                 + f_3 * dh_74[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_15, dg_18, dg_25, dg_45, dg_48, dg_55, dh_21, dh_24, \
                         dh_31, dh_64, dh_69, dh_78 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_19 * ab_x[k] * dg_15[k]
                 + f_20 * ab_x[k] * dg_18[k]
                 - f_19 * ab_x[k] * dg_25[k]
                 + f_21 * ab_y[k] * dg_45[k]
                 - f_22 * ab_y[k] * dg_48[k]
                 + f_21 * ab_y[k] * dg_55[k]
                 + f_19 * dh_21[k]
                 - f_20 * dh_24[k]
                 + f_19 * dh_31[k]
                 - f_21 * dh_64[k]
                 + f_22 * dh_69[k]
                 - f_21 * dh_78[k];
    }

#pragma omp simd aligned(ab_x, dg_61, dg_64, dg_66, dg_68, dg_71, dh_85, dh_88, dh_90, dh_92, \
                         dh_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_23 * ab_x[k] * dg_61[k]
                 + f_23 * ab_x[k] * dg_66[k]
                 + f_23 * dh_85[k]
                 - f_23 * dh_90[k];

        g_10[k] = -f_24 * ab_x[k] * dg_64[k]
                  + f_25 * ab_x[k] * dg_71[k]
                  + f_24 * dh_88[k]
                  - f_25 * dh_95[k];

        g_11[k] = f_26 * ab_x[k] * dg_61[k]
                  + f_26 * ab_x[k] * dg_66[k]
                  - f_27 * ab_x[k] * dg_68[k]
                  - f_26 * dh_85[k]
                  - f_26 * dh_90[k]
                  + f_27 * dh_92[k];
    }

#pragma omp simd aligned(ab_x, dg_64, dg_71, dg_73, dh_88, dh_95, \
                         dh_97 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_28 * ab_x[k] * dg_64[k]
                  + f_28 * ab_x[k] * dg_71[k]
                  - f_29 * ab_x[k] * dg_73[k]
                  - f_28 * dh_88[k]
                  - f_28 * dh_95[k]
                  + f_29 * dh_97[k];
    }

#pragma omp simd aligned(ab_x, dg_60, dg_63, dg_65, dg_70, dg_72, dg_74, dh_84, dh_87, dh_89, \
                         dh_94, dh_96, dh_98 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_30 * ab_x[k] * dg_60[k]
                  - f_31 * ab_x[k] * dg_63[k]
                  + f_32 * ab_x[k] * dg_65[k]
                  - f_30 * ab_x[k] * dg_70[k]
                  + f_32 * ab_x[k] * dg_72[k]
                  - f_33 * ab_x[k] * dg_74[k]
                  + f_30 * dh_84[k]
                  + f_31 * dh_87[k]
                  - f_32 * dh_89[k]
                  + f_30 * dh_94[k]
                  - f_32 * dh_96[k]
                  + f_33 * dh_98[k];
    }

#pragma omp simd aligned(ab_x, dg_62, dg_67, dg_69, dh_86, dh_91, \
                         dh_93 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_28 * ab_x[k] * dg_62[k]
                  + f_28 * ab_x[k] * dg_67[k]
                  - f_29 * ab_x[k] * dg_69[k]
                  - f_28 * dh_86[k]
                  - f_28 * dh_91[k]
                  + f_29 * dh_93[k];
    }

#pragma omp simd aligned(ab_x, dg_60, dg_62, dg_65, dg_67, dg_70, dg_72, dh_84, dh_86, dh_89, \
                         dh_91, dh_94, dh_96 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_34 * ab_x[k] * dg_60[k]
                  - f_35 * ab_x[k] * dg_65[k]
                  - f_34 * ab_x[k] * dg_70[k]
                  + f_35 * ab_x[k] * dg_72[k]
                  - f_34 * dh_84[k]
                  + f_35 * dh_89[k]
                  + f_34 * dh_94[k]
                  - f_35 * dh_96[k];

        g_16[k] = -f_25 * ab_x[k] * dg_62[k]
                  + f_24 * ab_x[k] * dg_67[k]
                  + f_25 * dh_86[k]
                  - f_24 * dh_91[k];
    }

#pragma omp simd aligned(ab_x, dg_60, dg_63, dg_70, dh_84, dh_87, \
                         dh_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_36 * ab_x[k] * dg_60[k]
                  + f_37 * ab_x[k] * dg_63[k]
                  - f_36 * ab_x[k] * dg_70[k]
                  + f_36 * dh_84[k]
                  - f_37 * dh_87[k]
                  + f_36 * dh_94[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_16, dg_21, dg_46, dg_51, dg_76, dg_81, dh_22, dh_27, \
                         dh_66, dh_73, dh_108, dh_115 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_38 * ab_x[k] * dg_16[k]
                  - f_38 * ab_x[k] * dg_21[k]
                  + f_38 * ab_y[k] * dg_46[k]
                  - f_38 * ab_y[k] * dg_51[k]
                  - f_39 * ab_y[k] * dg_76[k]
                  + f_39 * ab_y[k] * dg_81[k]
                  - f_38 * dh_22[k]
                  + f_38 * dh_27[k]
                  - f_38 * dh_66[k]
                  + f_38 * dh_73[k]
                  + f_39 * dh_108[k]
                  - f_39 * dh_115[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_19, dg_26, dg_49, dg_56, dg_79, dg_86, dh_25, dh_32, \
                         dh_70, dh_79, dh_112, dh_121 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_40 * ab_x[k] * dg_19[k]
                  - f_41 * ab_x[k] * dg_26[k]
                  + f_40 * ab_y[k] * dg_49[k]
                  - f_41 * ab_y[k] * dg_56[k]
                  - f_42 * ab_y[k] * dg_79[k]
                  + f_43 * ab_y[k] * dg_86[k]
                  - f_40 * dh_25[k]
                  + f_41 * dh_32[k]
                  - f_40 * dh_70[k]
                  + f_41 * dh_79[k]
                  + f_42 * dh_112[k]
                  - f_43 * dh_121[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_16, dg_21, dg_23, dg_46, dg_51, dg_53, dg_76, dg_81, \
                         dg_83, dh_22, dh_27, dh_29, dh_66, dh_73, dh_75, dh_108, dh_115, \
                         dh_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_44 * ab_x[k] * dg_16[k]
                  - f_44 * ab_x[k] * dg_21[k]
                  + f_45 * ab_x[k] * dg_23[k]
                  - f_44 * ab_y[k] * dg_46[k]
                  - f_44 * ab_y[k] * dg_51[k]
                  + f_45 * ab_y[k] * dg_53[k]
                  + f_46 * ab_y[k] * dg_76[k]
                  + f_46 * ab_y[k] * dg_81[k]
                  - f_47 * ab_y[k] * dg_83[k]
                  + f_44 * dh_22[k]
                  + f_44 * dh_27[k]
                  - f_45 * dh_29[k]
                  + f_44 * dh_66[k]
                  + f_44 * dh_73[k]
                  - f_45 * dh_75[k]
                  - f_46 * dh_108[k]
                  - f_46 * dh_115[k]
                  + f_47 * dh_117[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_19, dg_26, dg_28, dg_49, dg_56, dg_58, dg_79, dg_86, \
                         dg_88, dh_25, dh_32, dh_34, dh_70, dh_79, dh_81, dh_112, dh_121, \
                         dh_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_30 * ab_x[k] * dg_19[k]
                  - f_30 * ab_x[k] * dg_26[k]
                  + f_48 * ab_x[k] * dg_28[k]
                  - f_30 * ab_y[k] * dg_49[k]
                  - f_30 * ab_y[k] * dg_56[k]
                  + f_48 * ab_y[k] * dg_58[k]
                  + f_49 * ab_y[k] * dg_79[k]
                  + f_49 * ab_y[k] * dg_86[k]
                  - f_50 * ab_y[k] * dg_88[k]
                  + f_30 * dh_25[k]
                  + f_30 * dh_32[k]
                  - f_48 * dh_34[k]
                  + f_30 * dh_70[k]
                  + f_30 * dh_79[k]
                  - f_48 * dh_81[k]
                  - f_49 * dh_112[k]
                  - f_49 * dh_121[k]
                  + f_50 * dh_123[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_15, dg_18, dg_20, dg_25, dg_27, dg_29, dg_45, dg_48, \
                         dg_50, dg_55, dg_57, dg_59, dg_75, dg_78, dg_80, dg_85, dg_87, dg_89, \
                         dh_21, dh_24, dh_26, dh_31, dh_33, dh_35, dh_64, dh_69, dh_71, dh_78, \
                         dh_80, dh_82, dh_106, dh_111, dh_113, dh_120, dh_122, \
                         dh_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_51 * ab_x[k] * dg_15[k]
                  + f_52 * ab_x[k] * dg_18[k]
                  - f_53 * ab_x[k] * dg_20[k]
                  + f_51 * ab_x[k] * dg_25[k]
                  - f_53 * ab_x[k] * dg_27[k]
                  + f_54 * ab_x[k] * dg_29[k]
                  + f_51 * ab_y[k] * dg_45[k]
                  + f_52 * ab_y[k] * dg_48[k]
                  - f_53 * ab_y[k] * dg_50[k]
                  + f_51 * ab_y[k] * dg_55[k]
                  - f_53 * ab_y[k] * dg_57[k]
                  + f_54 * ab_y[k] * dg_59[k]
                  - f_55 * ab_y[k] * dg_75[k]
                  - f_53 * ab_y[k] * dg_78[k]
                  + f_56 * ab_y[k] * dg_80[k]
                  - f_55 * ab_y[k] * dg_85[k]
                  + f_56 * ab_y[k] * dg_87[k]
                  - f_57 * ab_y[k] * dg_89[k]
                  - f_51 * dh_21[k]
                  - f_52 * dh_24[k]
                  + f_53 * dh_26[k]
                  - f_51 * dh_31[k]
                  + f_53 * dh_33[k]
                  - f_54 * dh_35[k]
                  - f_51 * dh_64[k]
                  - f_52 * dh_69[k]
                  + f_53 * dh_71[k]
                  - f_51 * dh_78[k]
                  + f_53 * dh_80[k]
                  - f_54 * dh_82[k]
                  + f_55 * dh_106[k]
                  + f_53 * dh_111[k]
                  - f_56 * dh_113[k]
                  + f_55 * dh_120[k]
                  - f_56 * dh_122[k]
                  + f_57 * dh_124[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_17, dg_22, dg_24, dg_47, dg_52, dg_54, dg_77, dg_82, \
                         dg_84, dh_23, dh_28, dh_30, dh_67, dh_74, dh_76, dh_109, dh_116, \
                         dh_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_30 * ab_x[k] * dg_17[k]
                  - f_30 * ab_x[k] * dg_22[k]
                  + f_48 * ab_x[k] * dg_24[k]
                  - f_30 * ab_y[k] * dg_47[k]
                  - f_30 * ab_y[k] * dg_52[k]
                  + f_48 * ab_y[k] * dg_54[k]
                  + f_49 * ab_y[k] * dg_77[k]
                  + f_49 * ab_y[k] * dg_82[k]
                  - f_50 * ab_y[k] * dg_84[k]
                  + f_30 * dh_23[k]
                  + f_30 * dh_28[k]
                  - f_48 * dh_30[k]
                  + f_30 * dh_67[k]
                  + f_30 * dh_74[k]
                  - f_48 * dh_76[k]
                  - f_49 * dh_109[k]
                  - f_49 * dh_116[k]
                  + f_50 * dh_118[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_15, dg_20, dg_25, dg_27, dg_45, dg_50, dg_55, dg_57, \
                         dg_75, dg_80, dg_85, dg_87, dh_21, dh_26, dh_31, dh_33, dh_64, dh_71, \
                         dh_78, dh_80, dh_106, dh_113, dh_120, dh_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -f_58 * ab_x[k] * dg_15[k]
                  + f_59 * ab_x[k] * dg_20[k]
                  + f_58 * ab_x[k] * dg_25[k]
                  - f_59 * ab_x[k] * dg_27[k]
                  - f_58 * ab_y[k] * dg_45[k]
                  + f_59 * ab_y[k] * dg_50[k]
                  + f_58 * ab_y[k] * dg_55[k]
                  - f_59 * ab_y[k] * dg_57[k]
                  + f_60 * ab_y[k] * dg_75[k]
                  - f_61 * ab_y[k] * dg_80[k]
                  - f_60 * ab_y[k] * dg_85[k]
                  + f_61 * ab_y[k] * dg_87[k]
                  + f_58 * dh_21[k]
                  - f_59 * dh_26[k]
                  - f_58 * dh_31[k]
                  + f_59 * dh_33[k]
                  + f_58 * dh_64[k]
                  - f_59 * dh_71[k]
                  - f_58 * dh_78[k]
                  + f_59 * dh_80[k]
                  - f_60 * dh_106[k]
                  + f_61 * dh_113[k]
                  + f_60 * dh_120[k]
                  - f_61 * dh_122[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_17, dg_22, dg_47, dg_52, dg_77, dg_82, dh_23, dh_28, \
                         dh_67, dh_74, dh_109, dh_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_41 * ab_x[k] * dg_17[k]
                  - f_40 * ab_x[k] * dg_22[k]
                  + f_41 * ab_y[k] * dg_47[k]
                  - f_40 * ab_y[k] * dg_52[k]
                  - f_43 * ab_y[k] * dg_77[k]
                  + f_42 * ab_y[k] * dg_82[k]
                  - f_41 * dh_23[k]
                  + f_40 * dh_28[k]
                  - f_41 * dh_67[k]
                  + f_40 * dh_74[k]
                  + f_43 * dh_109[k]
                  - f_42 * dh_116[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_15, dg_18, dg_25, dg_45, dg_48, dg_55, dg_75, dg_78, \
                         dg_85, dh_21, dh_24, dh_31, dh_64, dh_69, dh_78, dh_106, dh_111, \
                         dh_120 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_62 * ab_x[k] * dg_15[k]
                  - f_63 * ab_x[k] * dg_18[k]
                  + f_62 * ab_x[k] * dg_25[k]
                  + f_62 * ab_y[k] * dg_45[k]
                  - f_63 * ab_y[k] * dg_48[k]
                  + f_62 * ab_y[k] * dg_55[k]
                  - f_38 * ab_y[k] * dg_75[k]
                  + f_64 * ab_y[k] * dg_78[k]
                  - f_38 * ab_y[k] * dg_85[k]
                  - f_62 * dh_21[k]
                  + f_63 * dh_24[k]
                  - f_62 * dh_31[k]
                  - f_62 * dh_64[k]
                  + f_63 * dh_69[k]
                  - f_62 * dh_78[k]
                  + f_38 * dh_106[k]
                  - f_64 * dh_111[k]
                  + f_38 * dh_120[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dg_31, dg_36, dg_61, dg_66, dg_76, dg_81, dh_43, \
                         dh_48, dh_87, dh_94, dh_109, dh_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = f_65 * ab_x[k] * dg_31[k]
                  - f_65 * ab_x[k] * dg_36[k]
                  + f_65 * ab_y[k] * dg_61[k]
                  - f_65 * ab_y[k] * dg_66[k]
                  - f_66 * ab_z[k] * dg_76[k]
                  + f_66 * ab_z[k] * dg_81[k]
                  - f_65 * dh_43[k]
                  + f_65 * dh_48[k]
                  - f_65 * dh_87[k]
                  + f_65 * dh_94[k]
                  + f_66 * dh_109[k]
                  - f_66 * dh_116[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dg_34, dg_41, dg_64, dg_71, dg_79, dg_86, dh_46, \
                         dh_53, dh_91, dh_100, dh_113, dh_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_67 * ab_x[k] * dg_34[k]
                  - f_68 * ab_x[k] * dg_41[k]
                  + f_67 * ab_y[k] * dg_64[k]
                  - f_68 * ab_y[k] * dg_71[k]
                  - f_69 * ab_z[k] * dg_79[k]
                  + f_70 * ab_z[k] * dg_86[k]
                  - f_67 * dh_46[k]
                  + f_68 * dh_53[k]
                  - f_67 * dh_91[k]
                  + f_68 * dh_100[k]
                  + f_69 * dh_113[k]
                  - f_70 * dh_122[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dg_31, dg_36, dg_38, dg_61, dg_66, dg_68, dg_76, \
                         dg_81, dg_83, dh_43, dh_48, dh_50, dh_87, dh_94, dh_96, dh_109, \
                         dh_116, dh_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = -f_71 * ab_x[k] * dg_31[k]
                  - f_71 * ab_x[k] * dg_36[k]
                  + f_72 * ab_x[k] * dg_38[k]
                  - f_71 * ab_y[k] * dg_61[k]
                  - f_71 * ab_y[k] * dg_66[k]
                  + f_72 * ab_y[k] * dg_68[k]
                  + f_73 * ab_z[k] * dg_76[k]
                  + f_73 * ab_z[k] * dg_81[k]
                  - f_74 * ab_z[k] * dg_83[k]
                  + f_71 * dh_43[k]
                  + f_71 * dh_48[k]
                  - f_72 * dh_50[k]
                  + f_71 * dh_87[k]
                  + f_71 * dh_94[k]
                  - f_72 * dh_96[k]
                  - f_73 * dh_109[k]
                  - f_73 * dh_116[k]
                  + f_74 * dh_118[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dg_34, dg_41, dg_43, dg_64, dg_71, dg_73, dg_79, \
                         dg_86, dg_88, dh_46, dh_53, dh_55, dh_91, dh_100, dh_102, dh_113, \
                         dh_122, dh_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_75 * ab_x[k] * dg_34[k]
                  - f_75 * ab_x[k] * dg_41[k]
                  + f_76 * ab_x[k] * dg_43[k]
                  - f_75 * ab_y[k] * dg_64[k]
                  - f_75 * ab_y[k] * dg_71[k]
                  + f_76 * ab_y[k] * dg_73[k]
                  + f_12 * ab_z[k] * dg_79[k]
                  + f_12 * ab_z[k] * dg_86[k]
                  - f_77 * ab_z[k] * dg_88[k]
                  + f_75 * dh_46[k]
                  + f_75 * dh_53[k]
                  - f_76 * dh_55[k]
                  + f_75 * dh_91[k]
                  + f_75 * dh_100[k]
                  - f_76 * dh_102[k]
                  - f_12 * dh_113[k]
                  - f_12 * dh_122[k]
                  + f_77 * dh_124[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dg_30, dg_33, dg_35, dg_40, dg_42, dg_44, dg_60, \
                         dg_63, dg_65, dg_70, dg_72, dg_74, dg_75, dg_78, dg_80, dg_85, dg_87, \
                         dg_89, dh_42, dh_45, dh_47, dh_52, dh_54, dh_56, dh_85, dh_90, dh_92, \
                         dh_99, dh_101, dh_103, dh_107, dh_112, dh_114, dh_121, dh_123, \
                         dh_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = 0.5625 * ab_x[k] * dg_30[k]
                  + 1.125 * ab_x[k] * dg_33[k]
                  - 4.5 * ab_x[k] * dg_35[k]
                  + 0.5625 * ab_x[k] * dg_40[k]
                  - 4.5 * ab_x[k] * dg_42[k]
                  + 1.5 * ab_x[k] * dg_44[k]
                  + 0.5625 * ab_y[k] * dg_60[k]
                  + 1.125 * ab_y[k] * dg_63[k]
                  - 4.5 * ab_y[k] * dg_65[k]
                  + 0.5625 * ab_y[k] * dg_70[k]
                  - 4.5 * ab_y[k] * dg_72[k]
                  + 1.5 * ab_y[k] * dg_74[k]
                  - 0.375 * ab_z[k] * dg_75[k]
                  - 0.75 * ab_z[k] * dg_78[k]
                  + 3.0 * ab_z[k] * dg_80[k]
                  - 0.375 * ab_z[k] * dg_85[k]
                  + 3.0 * ab_z[k] * dg_87[k]
                  - ab_z[k] * dg_89[k]
                  - 0.5625 * dh_42[k]
                  - 1.125 * dh_45[k]
                  + 4.5 * dh_47[k]
                  - 0.5625 * dh_52[k]
                  + 4.5 * dh_54[k]
                  - 1.5 * dh_56[k]
                  - 0.5625 * dh_85[k]
                  - 1.125 * dh_90[k]
                  + 4.5 * dh_92[k]
                  - 0.5625 * dh_99[k]
                  + 4.5 * dh_101[k]
                  - 1.5 * dh_103[k]
                  + 0.375 * dh_107[k]
                  + 0.75 * dh_112[k]
                  - 3.0 * dh_114[k]
                  + 0.375 * dh_121[k]
                  - 3.0 * dh_123[k]
                  + dh_125[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dg_32, dg_37, dg_39, dg_62, dg_67, dg_69, dg_77, \
                         dg_82, dg_84, dh_44, dh_49, dh_51, dh_88, dh_95, dh_97, dh_110, \
                         dh_117, dh_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_75 * ab_x[k] * dg_32[k]
                  - f_75 * ab_x[k] * dg_37[k]
                  + f_76 * ab_x[k] * dg_39[k]
                  - f_75 * ab_y[k] * dg_62[k]
                  - f_75 * ab_y[k] * dg_67[k]
                  + f_76 * ab_y[k] * dg_69[k]
                  + f_12 * ab_z[k] * dg_77[k]
                  + f_12 * ab_z[k] * dg_82[k]
                  - f_77 * ab_z[k] * dg_84[k]
                  + f_75 * dh_44[k]
                  + f_75 * dh_49[k]
                  - f_76 * dh_51[k]
                  + f_75 * dh_88[k]
                  + f_75 * dh_95[k]
                  - f_76 * dh_97[k]
                  - f_12 * dh_110[k]
                  - f_12 * dh_117[k]
                  + f_77 * dh_119[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dg_30, dg_35, dg_40, dg_42, dg_60, dg_65, dg_70, \
                         dg_72, dg_75, dg_80, dg_85, dg_87, dh_42, dh_47, dh_52, dh_54, dh_85, \
                         dh_92, dh_99, dh_101, dh_107, dh_114, dh_121, \
                         dh_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_78 * ab_x[k] * dg_30[k]
                  + f_79 * ab_x[k] * dg_35[k]
                  + f_78 * ab_x[k] * dg_40[k]
                  - f_79 * ab_x[k] * dg_42[k]
                  - f_78 * ab_y[k] * dg_60[k]
                  + f_79 * ab_y[k] * dg_65[k]
                  + f_78 * ab_y[k] * dg_70[k]
                  - f_79 * ab_y[k] * dg_72[k]
                  + f_80 * ab_z[k] * dg_75[k]
                  - f_81 * ab_z[k] * dg_80[k]
                  - f_80 * ab_z[k] * dg_85[k]
                  + f_81 * ab_z[k] * dg_87[k]
                  + f_78 * dh_42[k]
                  - f_79 * dh_47[k]
                  - f_78 * dh_52[k]
                  + f_79 * dh_54[k]
                  + f_78 * dh_85[k]
                  - f_79 * dh_92[k]
                  - f_78 * dh_99[k]
                  + f_79 * dh_101[k]
                  - f_80 * dh_107[k]
                  + f_81 * dh_114[k]
                  + f_80 * dh_121[k]
                  - f_81 * dh_123[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dg_32, dg_37, dg_62, dg_67, dg_77, dg_82, dh_44, \
                         dh_49, dh_88, dh_95, dh_110, dh_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_68 * ab_x[k] * dg_32[k]
                  - f_67 * ab_x[k] * dg_37[k]
                  + f_68 * ab_y[k] * dg_62[k]
                  - f_67 * ab_y[k] * dg_67[k]
                  - f_70 * ab_z[k] * dg_77[k]
                  + f_69 * ab_z[k] * dg_82[k]
                  - f_68 * dh_44[k]
                  + f_67 * dh_49[k]
                  - f_68 * dh_88[k]
                  + f_67 * dh_95[k]
                  + f_70 * dh_110[k]
                  - f_69 * dh_117[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dg_30, dg_33, dg_40, dg_60, dg_63, dg_70, dg_75, \
                         dg_78, dg_85, dh_42, dh_45, dh_52, dh_85, dh_90, dh_99, dh_107, \
                         dh_112, dh_121 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_82 * ab_x[k] * dg_30[k]
                  - f_83 * ab_x[k] * dg_33[k]
                  + f_82 * ab_x[k] * dg_40[k]
                  + f_82 * ab_y[k] * dg_60[k]
                  - f_83 * ab_y[k] * dg_63[k]
                  + f_82 * ab_y[k] * dg_70[k]
                  - f_84 * ab_z[k] * dg_75[k]
                  + f_65 * ab_z[k] * dg_78[k]
                  - f_84 * ab_z[k] * dg_85[k]
                  - f_82 * dh_42[k]
                  + f_83 * dh_45[k]
                  - f_82 * dh_52[k]
                  - f_82 * dh_85[k]
                  + f_83 * dh_90[k]
                  - f_82 * dh_99[k]
                  + f_84 * dh_107[k]
                  - f_65 * dh_112[k]
                  + f_84 * dh_121[k];
    }

#pragma omp simd aligned(ab_x, dg_1, dg_6, dg_46, dg_51, dg_76, dg_81, dh_1, dh_6, dh_64, \
                         dh_69, dh_106, dh_111 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_38 * ab_x[k] * dg_1[k]
                  - f_38 * ab_x[k] * dg_6[k]
                  + f_38 * ab_x[k] * dg_46[k]
                  - f_38 * ab_x[k] * dg_51[k]
                  - f_39 * ab_x[k] * dg_76[k]
                  + f_39 * ab_x[k] * dg_81[k]
                  - f_38 * dh_1[k]
                  + f_38 * dh_6[k]
                  - f_38 * dh_64[k]
                  + f_38 * dh_69[k]
                  + f_39 * dh_106[k]
                  - f_39 * dh_111[k];
    }

#pragma omp simd aligned(ab_x, dg_4, dg_11, dg_49, dg_56, dg_79, dg_86, dh_4, dh_11, dh_67, \
                         dh_74, dh_109, dh_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_40 * ab_x[k] * dg_4[k]
                  - f_41 * ab_x[k] * dg_11[k]
                  + f_40 * ab_x[k] * dg_49[k]
                  - f_41 * ab_x[k] * dg_56[k]
                  - f_42 * ab_x[k] * dg_79[k]
                  + f_43 * ab_x[k] * dg_86[k]
                  - f_40 * dh_4[k]
                  + f_41 * dh_11[k]
                  - f_40 * dh_67[k]
                  + f_41 * dh_74[k]
                  + f_42 * dh_109[k]
                  - f_43 * dh_116[k];
    }

#pragma omp simd aligned(ab_x, dg_1, dg_6, dg_8, dg_46, dg_51, dg_53, dg_76, dg_81, dg_83, \
                         dh_1, dh_6, dh_8, dh_64, dh_69, dh_71, dh_106, dh_111, \
                         dh_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_44 * ab_x[k] * dg_1[k]
                  - f_44 * ab_x[k] * dg_6[k]
                  + f_45 * ab_x[k] * dg_8[k]
                  - f_44 * ab_x[k] * dg_46[k]
                  - f_44 * ab_x[k] * dg_51[k]
                  + f_45 * ab_x[k] * dg_53[k]
                  + f_46 * ab_x[k] * dg_76[k]
                  + f_46 * ab_x[k] * dg_81[k]
                  - f_47 * ab_x[k] * dg_83[k]
                  + f_44 * dh_1[k]
                  + f_44 * dh_6[k]
                  - f_45 * dh_8[k]
                  + f_44 * dh_64[k]
                  + f_44 * dh_69[k]
                  - f_45 * dh_71[k]
                  - f_46 * dh_106[k]
                  - f_46 * dh_111[k]
                  + f_47 * dh_113[k];
    }

#pragma omp simd aligned(ab_x, dg_4, dg_11, dg_13, dg_49, dg_56, dg_58, dg_79, dg_86, dg_88, \
                         dh_4, dh_11, dh_13, dh_67, dh_74, dh_76, dh_109, dh_116, \
                         dh_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_30 * ab_x[k] * dg_4[k]
                  - f_30 * ab_x[k] * dg_11[k]
                  + f_48 * ab_x[k] * dg_13[k]
                  - f_30 * ab_x[k] * dg_49[k]
                  - f_30 * ab_x[k] * dg_56[k]
                  + f_48 * ab_x[k] * dg_58[k]
                  + f_49 * ab_x[k] * dg_79[k]
                  + f_49 * ab_x[k] * dg_86[k]
                  - f_50 * ab_x[k] * dg_88[k]
                  + f_30 * dh_4[k]
                  + f_30 * dh_11[k]
                  - f_48 * dh_13[k]
                  + f_30 * dh_67[k]
                  + f_30 * dh_74[k]
                  - f_48 * dh_76[k]
                  - f_49 * dh_109[k]
                  - f_49 * dh_116[k]
                  + f_50 * dh_118[k];
    }

#pragma omp simd aligned(ab_x, dg_0, dg_3, dg_5, dg_10, dg_12, dg_14, dg_45, dg_48, dg_50, \
                         dg_55, dg_57, dg_59, dg_75, dg_78, dg_80, dg_85, dg_87, dg_89, dh_0, \
                         dh_3, dh_5, dh_10, dh_12, dh_14, dh_63, dh_66, dh_68, dh_73, dh_75, \
                         dh_77, dh_105, dh_108, dh_110, dh_115, dh_117, \
                         dh_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_51 * ab_x[k] * dg_0[k]
                  + f_52 * ab_x[k] * dg_3[k]
                  - f_53 * ab_x[k] * dg_5[k]
                  + f_51 * ab_x[k] * dg_10[k]
                  - f_53 * ab_x[k] * dg_12[k]
                  + f_54 * ab_x[k] * dg_14[k]
                  + f_51 * ab_x[k] * dg_45[k]
                  + f_52 * ab_x[k] * dg_48[k]
                  - f_53 * ab_x[k] * dg_50[k]
                  + f_51 * ab_x[k] * dg_55[k]
                  - f_53 * ab_x[k] * dg_57[k]
                  + f_54 * ab_x[k] * dg_59[k]
                  - f_55 * ab_x[k] * dg_75[k]
                  - f_53 * ab_x[k] * dg_78[k]
                  + f_56 * ab_x[k] * dg_80[k]
                  - f_55 * ab_x[k] * dg_85[k]
                  + f_56 * ab_x[k] * dg_87[k]
                  - f_57 * ab_x[k] * dg_89[k]
                  - f_51 * dh_0[k]
                  - f_52 * dh_3[k]
                  + f_53 * dh_5[k]
                  - f_51 * dh_10[k]
                  + f_53 * dh_12[k]
                  - f_54 * dh_14[k]
                  - f_51 * dh_63[k]
                  - f_52 * dh_66[k]
                  + f_53 * dh_68[k]
                  - f_51 * dh_73[k]
                  + f_53 * dh_75[k]
                  - f_54 * dh_77[k]
                  + f_55 * dh_105[k]
                  + f_53 * dh_108[k]
                  - f_56 * dh_110[k]
                  + f_55 * dh_115[k]
                  - f_56 * dh_117[k]
                  + f_57 * dh_119[k];
    }

#pragma omp simd aligned(ab_x, dg_2, dg_7, dg_9, dg_47, dg_52, dg_54, dg_77, dg_82, dg_84, \
                         dh_2, dh_7, dh_9, dh_65, dh_70, dh_72, dh_107, dh_112, \
                         dh_114 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_30 * ab_x[k] * dg_2[k]
                  - f_30 * ab_x[k] * dg_7[k]
                  + f_48 * ab_x[k] * dg_9[k]
                  - f_30 * ab_x[k] * dg_47[k]
                  - f_30 * ab_x[k] * dg_52[k]
                  + f_48 * ab_x[k] * dg_54[k]
                  + f_49 * ab_x[k] * dg_77[k]
                  + f_49 * ab_x[k] * dg_82[k]
                  - f_50 * ab_x[k] * dg_84[k]
                  + f_30 * dh_2[k]
                  + f_30 * dh_7[k]
                  - f_48 * dh_9[k]
                  + f_30 * dh_65[k]
                  + f_30 * dh_70[k]
                  - f_48 * dh_72[k]
                  - f_49 * dh_107[k]
                  - f_49 * dh_112[k]
                  + f_50 * dh_114[k];
    }

#pragma omp simd aligned(ab_x, dg_0, dg_5, dg_10, dg_12, dg_45, dg_50, dg_55, dg_57, dg_75, \
                         dg_80, dg_85, dg_87, dh_0, dh_5, dh_10, dh_12, dh_63, dh_68, dh_73, \
                         dh_75, dh_105, dh_110, dh_115, dh_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_58 * ab_x[k] * dg_0[k]
                  + f_59 * ab_x[k] * dg_5[k]
                  + f_58 * ab_x[k] * dg_10[k]
                  - f_59 * ab_x[k] * dg_12[k]
                  - f_58 * ab_x[k] * dg_45[k]
                  + f_59 * ab_x[k] * dg_50[k]
                  + f_58 * ab_x[k] * dg_55[k]
                  - f_59 * ab_x[k] * dg_57[k]
                  + f_60 * ab_x[k] * dg_75[k]
                  - f_61 * ab_x[k] * dg_80[k]
                  - f_60 * ab_x[k] * dg_85[k]
                  + f_61 * ab_x[k] * dg_87[k]
                  + f_58 * dh_0[k]
                  - f_59 * dh_5[k]
                  - f_58 * dh_10[k]
                  + f_59 * dh_12[k]
                  + f_58 * dh_63[k]
                  - f_59 * dh_68[k]
                  - f_58 * dh_73[k]
                  + f_59 * dh_75[k]
                  - f_60 * dh_105[k]
                  + f_61 * dh_110[k]
                  + f_60 * dh_115[k]
                  - f_61 * dh_117[k];
    }

#pragma omp simd aligned(ab_x, dg_2, dg_7, dg_47, dg_52, dg_77, dg_82, dh_2, dh_7, dh_65, \
                         dh_70, dh_107, dh_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_41 * ab_x[k] * dg_2[k]
                  - f_40 * ab_x[k] * dg_7[k]
                  + f_41 * ab_x[k] * dg_47[k]
                  - f_40 * ab_x[k] * dg_52[k]
                  - f_43 * ab_x[k] * dg_77[k]
                  + f_42 * ab_x[k] * dg_82[k]
                  - f_41 * dh_2[k]
                  + f_40 * dh_7[k]
                  - f_41 * dh_65[k]
                  + f_40 * dh_70[k]
                  + f_43 * dh_107[k]
                  - f_42 * dh_112[k];
    }

#pragma omp simd aligned(ab_x, dg_0, dg_3, dg_10, dg_45, dg_48, dg_55, dg_75, dg_78, dg_85, \
                         dh_0, dh_3, dh_10, dh_63, dh_66, dh_73, dh_105, dh_108, \
                         dh_115 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_62 * ab_x[k] * dg_0[k]
                  - f_63 * ab_x[k] * dg_3[k]
                  + f_62 * ab_x[k] * dg_10[k]
                  + f_62 * ab_x[k] * dg_45[k]
                  - f_63 * ab_x[k] * dg_48[k]
                  + f_62 * ab_x[k] * dg_55[k]
                  - f_38 * ab_x[k] * dg_75[k]
                  + f_64 * ab_x[k] * dg_78[k]
                  - f_38 * ab_x[k] * dg_85[k]
                  - f_62 * dh_0[k]
                  + f_63 * dh_3[k]
                  - f_62 * dh_10[k]
                  - f_62 * dh_63[k]
                  + f_63 * dh_66[k]
                  - f_62 * dh_73[k]
                  + f_38 * dh_105[k]
                  - f_64 * dh_108[k]
                  + f_38 * dh_115[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_31, dg_36, dg_61, dg_66, dh_43, dh_48, dh_87, \
                         dh_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_85 * ab_x[k] * dg_31[k]
                  + f_85 * ab_x[k] * dg_36[k]
                  + f_85 * ab_y[k] * dg_61[k]
                  - f_85 * ab_y[k] * dg_66[k]
                  + f_85 * dh_43[k]
                  - f_85 * dh_48[k]
                  - f_85 * dh_87[k]
                  + f_85 * dh_94[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_34, dg_41, dg_64, dg_71, dh_46, dh_53, dh_91, \
                         dh_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_86 * ab_x[k] * dg_34[k]
                  + f_87 * ab_x[k] * dg_41[k]
                  + f_86 * ab_y[k] * dg_64[k]
                  - f_87 * ab_y[k] * dg_71[k]
                  + f_86 * dh_46[k]
                  - f_87 * dh_53[k]
                  - f_86 * dh_91[k]
                  + f_87 * dh_100[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_31, dg_36, dg_38, dg_61, dg_66, dg_68, dh_43, dh_48, \
                         dh_50, dh_87, dh_94, dh_96 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = f_34 * ab_x[k] * dg_31[k]
                  + f_34 * ab_x[k] * dg_36[k]
                  - f_35 * ab_x[k] * dg_38[k]
                  - f_34 * ab_y[k] * dg_61[k]
                  - f_34 * ab_y[k] * dg_66[k]
                  + f_35 * ab_y[k] * dg_68[k]
                  - f_34 * dh_43[k]
                  - f_34 * dh_48[k]
                  + f_35 * dh_50[k]
                  + f_34 * dh_87[k]
                  + f_34 * dh_94[k]
                  - f_35 * dh_96[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_34, dg_41, dg_43, dg_64, dg_71, dg_73, dh_46, dh_53, \
                         dh_55, dh_91, dh_100, dh_102 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_88 * ab_x[k] * dg_34[k]
                  + f_88 * ab_x[k] * dg_41[k]
                  - f_89 * ab_x[k] * dg_43[k]
                  - f_88 * ab_y[k] * dg_64[k]
                  - f_88 * ab_y[k] * dg_71[k]
                  + f_89 * ab_y[k] * dg_73[k]
                  - f_88 * dh_46[k]
                  - f_88 * dh_53[k]
                  + f_89 * dh_55[k]
                  + f_88 * dh_91[k]
                  + f_88 * dh_100[k]
                  - f_89 * dh_102[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_30, dg_33, dg_35, dg_40, dg_42, dg_44, dg_60, dg_63, \
                         dg_65, dg_70, dg_72, dg_74, dh_42, dh_45, dh_47, dh_52, dh_54, dh_56, \
                         dh_85, dh_90, dh_92, dh_99, dh_101, dh_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_90 * ab_x[k] * dg_30[k]
                  - f_30 * ab_x[k] * dg_33[k]
                  + f_49 * ab_x[k] * dg_35[k]
                  - f_90 * ab_x[k] * dg_40[k]
                  + f_49 * ab_x[k] * dg_42[k]
                  - f_48 * ab_x[k] * dg_44[k]
                  + f_90 * ab_y[k] * dg_60[k]
                  + f_30 * ab_y[k] * dg_63[k]
                  - f_49 * ab_y[k] * dg_65[k]
                  + f_90 * ab_y[k] * dg_70[k]
                  - f_49 * ab_y[k] * dg_72[k]
                  + f_48 * ab_y[k] * dg_74[k]
                  + f_90 * dh_42[k]
                  + f_30 * dh_45[k]
                  - f_49 * dh_47[k]
                  + f_90 * dh_52[k]
                  - f_49 * dh_54[k]
                  + f_48 * dh_56[k]
                  - f_90 * dh_85[k]
                  - f_30 * dh_90[k]
                  + f_49 * dh_92[k]
                  - f_90 * dh_99[k]
                  + f_49 * dh_101[k]
                  - f_48 * dh_103[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_32, dg_37, dg_39, dg_62, dg_67, dg_69, dh_44, dh_49, \
                         dh_51, dh_88, dh_95, dh_97 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = f_88 * ab_x[k] * dg_32[k]
                  + f_88 * ab_x[k] * dg_37[k]
                  - f_89 * ab_x[k] * dg_39[k]
                  - f_88 * ab_y[k] * dg_62[k]
                  - f_88 * ab_y[k] * dg_67[k]
                  + f_89 * ab_y[k] * dg_69[k]
                  - f_88 * dh_44[k]
                  - f_88 * dh_49[k]
                  + f_89 * dh_51[k]
                  + f_88 * dh_88[k]
                  + f_88 * dh_95[k]
                  - f_89 * dh_97[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_30, dg_35, dg_40, dg_42, dg_60, dg_65, dg_70, dg_72, \
                         dh_42, dh_47, dh_52, dh_54, dh_85, dh_92, dh_99, \
                         dh_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_91 * ab_x[k] * dg_30[k]
                  - f_92 * ab_x[k] * dg_35[k]
                  - f_91 * ab_x[k] * dg_40[k]
                  + f_92 * ab_x[k] * dg_42[k]
                  - f_91 * ab_y[k] * dg_60[k]
                  + f_92 * ab_y[k] * dg_65[k]
                  + f_91 * ab_y[k] * dg_70[k]
                  - f_92 * ab_y[k] * dg_72[k]
                  - f_91 * dh_42[k]
                  + f_92 * dh_47[k]
                  + f_91 * dh_52[k]
                  - f_92 * dh_54[k]
                  + f_91 * dh_85[k]
                  - f_92 * dh_92[k]
                  - f_91 * dh_99[k]
                  + f_92 * dh_101[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_32, dg_37, dg_62, dg_67, dh_44, dh_49, dh_88, \
                         dh_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_87 * ab_x[k] * dg_32[k]
                  + f_86 * ab_x[k] * dg_37[k]
                  + f_87 * ab_y[k] * dg_62[k]
                  - f_86 * ab_y[k] * dg_67[k]
                  + f_87 * dh_44[k]
                  - f_86 * dh_49[k]
                  - f_87 * dh_88[k]
                  + f_86 * dh_95[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dg_30, dg_33, dg_40, dg_60, dg_63, dg_70, dh_42, dh_45, \
                         dh_52, dh_85, dh_90, dh_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_93 * ab_x[k] * dg_30[k]
                  + f_94 * ab_x[k] * dg_33[k]
                  - f_93 * ab_x[k] * dg_40[k]
                  + f_93 * ab_y[k] * dg_60[k]
                  - f_94 * ab_y[k] * dg_63[k]
                  + f_93 * ab_y[k] * dg_70[k]
                  + f_93 * dh_42[k]
                  - f_94 * dh_45[k]
                  + f_93 * dh_52[k]
                  - f_93 * dh_85[k]
                  + f_94 * dh_90[k]
                  - f_93 * dh_99[k];
    }

#pragma omp simd aligned(ab_x, dg_1, dg_6, dg_46, dg_51, dh_1, dh_6, dh_64, \
                         dh_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_1 * ab_x[k] * dg_1[k]
                  + f_1 * ab_x[k] * dg_6[k]
                  + f_0 * ab_x[k] * dg_46[k]
                  - f_0 * ab_x[k] * dg_51[k]
                  + f_1 * dh_1[k]
                  - f_1 * dh_6[k]
                  - f_0 * dh_64[k]
                  + f_0 * dh_69[k];
    }

#pragma omp simd aligned(ab_x, dg_4, dg_11, dg_49, dg_56, dh_4, dh_11, dh_67, \
                         dh_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_3 * ab_x[k] * dg_4[k]
                  + f_4 * ab_x[k] * dg_11[k]
                  + f_2 * ab_x[k] * dg_49[k]
                  - f_3 * ab_x[k] * dg_56[k]
                  + f_3 * dh_4[k]
                  - f_4 * dh_11[k]
                  - f_2 * dh_67[k]
                  + f_3 * dh_74[k];
    }

#pragma omp simd aligned(ab_x, dg_1, dg_6, dg_8, dg_46, dg_51, dg_53, dh_1, dh_6, dh_8, dh_64, \
                         dh_69, dh_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_7 * ab_x[k] * dg_1[k]
                  + f_7 * ab_x[k] * dg_6[k]
                  - f_8 * ab_x[k] * dg_8[k]
                  - f_5 * ab_x[k] * dg_46[k]
                  - f_5 * ab_x[k] * dg_51[k]
                  + f_6 * ab_x[k] * dg_53[k]
                  - f_7 * dh_1[k]
                  - f_7 * dh_6[k]
                  + f_8 * dh_8[k]
                  + f_5 * dh_64[k]
                  + f_5 * dh_69[k]
                  - f_6 * dh_71[k];
    }

#pragma omp simd aligned(ab_x, dg_4, dg_11, dg_13, dg_49, dg_56, dg_58, dh_4, dh_11, dh_13, \
                         dh_67, dh_74, dh_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = 1.875 * ab_x[k] * dg_4[k]
                  + 1.875 * ab_x[k] * dg_11[k]
                  - 2.5 * ab_x[k] * dg_13[k]
                  - 5.625 * ab_x[k] * dg_49[k]
                  - 5.625 * ab_x[k] * dg_56[k]
                  + 7.5 * ab_x[k] * dg_58[k]
                  - 1.875 * dh_4[k]
                  - 1.875 * dh_11[k]
                  + 2.5 * dh_13[k]
                  + 5.625 * dh_67[k]
                  + 5.625 * dh_74[k]
                  - 7.5 * dh_76[k];
    }

#pragma omp simd aligned(ab_x, dg_0, dg_3, dg_5, dg_10, dg_12, dg_14, dg_45, dg_48, dg_50, \
                         dg_55, dg_57, dg_59, dh_0, dh_3, dh_5, dh_10, dh_12, dh_14, dh_63, \
                         dh_66, dh_68, dh_73, dh_75, dh_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_13 * ab_x[k] * dg_0[k]
                  - f_14 * ab_x[k] * dg_3[k]
                  + f_12 * ab_x[k] * dg_5[k]
                  - f_13 * ab_x[k] * dg_10[k]
                  + f_12 * ab_x[k] * dg_12[k]
                  - f_15 * ab_x[k] * dg_14[k]
                  + f_9 * ab_x[k] * dg_45[k]
                  + f_10 * ab_x[k] * dg_48[k]
                  - f_11 * ab_x[k] * dg_50[k]
                  + f_9 * ab_x[k] * dg_55[k]
                  - f_11 * ab_x[k] * dg_57[k]
                  + f_12 * ab_x[k] * dg_59[k]
                  + f_13 * dh_0[k]
                  + f_14 * dh_3[k]
                  - f_12 * dh_5[k]
                  + f_13 * dh_10[k]
                  - f_12 * dh_12[k]
                  + f_15 * dh_14[k]
                  - f_9 * dh_63[k]
                  - f_10 * dh_66[k]
                  + f_11 * dh_68[k]
                  - f_9 * dh_73[k]
                  + f_11 * dh_75[k]
                  - f_12 * dh_77[k];
    }

#pragma omp simd aligned(ab_x, dg_2, dg_7, dg_9, dg_47, dg_52, dg_54, dh_2, dh_7, dh_9, dh_65, \
                         dh_70, dh_72 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = 1.875 * ab_x[k] * dg_2[k]
                  + 1.875 * ab_x[k] * dg_7[k]
                  - 2.5 * ab_x[k] * dg_9[k]
                  - 5.625 * ab_x[k] * dg_47[k]
                  - 5.625 * ab_x[k] * dg_52[k]
                  + 7.5 * ab_x[k] * dg_54[k]
                  - 1.875 * dh_2[k]
                  - 1.875 * dh_7[k]
                  + 2.5 * dh_9[k]
                  + 5.625 * dh_65[k]
                  + 5.625 * dh_70[k]
                  - 7.5 * dh_72[k];
    }

#pragma omp simd aligned(ab_x, dg_0, dg_5, dg_10, dg_12, dg_45, dg_50, dg_55, dg_57, dh_0, \
                         dh_5, dh_10, dh_12, dh_63, dh_68, dh_73, \
                         dh_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_18 * ab_x[k] * dg_0[k]
                  - f_5 * ab_x[k] * dg_5[k]
                  - f_18 * ab_x[k] * dg_10[k]
                  + f_5 * ab_x[k] * dg_12[k]
                  - f_16 * ab_x[k] * dg_45[k]
                  + f_17 * ab_x[k] * dg_50[k]
                  + f_16 * ab_x[k] * dg_55[k]
                  - f_17 * ab_x[k] * dg_57[k]
                  - f_18 * dh_0[k]
                  + f_5 * dh_5[k]
                  + f_18 * dh_10[k]
                  - f_5 * dh_12[k]
                  + f_16 * dh_63[k]
                  - f_17 * dh_68[k]
                  - f_16 * dh_73[k]
                  + f_17 * dh_75[k];
    }

#pragma omp simd aligned(ab_x, dg_2, dg_7, dg_47, dg_52, dh_2, dh_7, dh_65, \
                         dh_70 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_4 * ab_x[k] * dg_2[k]
                  + f_3 * ab_x[k] * dg_7[k]
                  + f_3 * ab_x[k] * dg_47[k]
                  - f_2 * ab_x[k] * dg_52[k]
                  + f_4 * dh_2[k]
                  - f_3 * dh_7[k]
                  - f_3 * dh_65[k]
                  + f_2 * dh_70[k];
    }

#pragma omp simd aligned(ab_x, dg_0, dg_3, dg_10, dg_45, dg_48, dg_55, dh_0, dh_3, dh_10, \
                         dh_63, dh_66, dh_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_21 * ab_x[k] * dg_0[k]
                  + f_22 * ab_x[k] * dg_3[k]
                  - f_21 * ab_x[k] * dg_10[k]
                  + f_19 * ab_x[k] * dg_45[k]
                  - f_20 * ab_x[k] * dg_48[k]
                  + f_19 * ab_x[k] * dg_55[k]
                  + f_21 * dh_0[k]
                  - f_22 * dh_3[k]
                  + f_21 * dh_10[k]
                  - f_19 * dh_63[k]
                  + f_20 * dh_66[k]
                  - f_19 * dh_73[k];
    }
}

auto
compute_hrr_fg(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t dg, const size_t dh, const size_t nmax) -> void
{
    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);
    auto *t_9 = buffer.data(target + 9);
    auto *t_10 = buffer.data(target + 10);
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);
    auto *t_15 = buffer.data(target + 15);
    auto *t_16 = buffer.data(target + 16);
    auto *t_17 = buffer.data(target + 17);
    auto *t_18 = buffer.data(target + 18);
    auto *t_19 = buffer.data(target + 19);
    auto *t_20 = buffer.data(target + 20);
    auto *t_21 = buffer.data(target + 21);
    auto *t_22 = buffer.data(target + 22);
    auto *t_23 = buffer.data(target + 23);
    auto *t_24 = buffer.data(target + 24);
    auto *t_25 = buffer.data(target + 25);
    auto *t_26 = buffer.data(target + 26);
    auto *t_27 = buffer.data(target + 27);
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);
    auto *t_30 = buffer.data(target + 30);
    auto *t_31 = buffer.data(target + 31);
    auto *t_32 = buffer.data(target + 32);
    auto *t_33 = buffer.data(target + 33);
    auto *t_34 = buffer.data(target + 34);
    auto *t_35 = buffer.data(target + 35);
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);
    auto *t_45 = buffer.data(target + 45);
    auto *t_46 = buffer.data(target + 46);
    auto *t_47 = buffer.data(target + 47);
    auto *t_48 = buffer.data(target + 48);
    auto *t_49 = buffer.data(target + 49);
    auto *t_50 = buffer.data(target + 50);
    auto *t_51 = buffer.data(target + 51);
    auto *t_52 = buffer.data(target + 52);
    auto *t_53 = buffer.data(target + 53);
    auto *t_54 = buffer.data(target + 54);
    auto *t_55 = buffer.data(target + 55);
    auto *t_56 = buffer.data(target + 56);
    auto *t_57 = buffer.data(target + 57);
    auto *t_58 = buffer.data(target + 58);
    auto *t_59 = buffer.data(target + 59);
    auto *t_60 = buffer.data(target + 60);
    auto *t_61 = buffer.data(target + 61);
    auto *t_62 = buffer.data(target + 62);
    auto *t_63 = buffer.data(target + 63);
    auto *t_64 = buffer.data(target + 64);
    auto *t_65 = buffer.data(target + 65);
    auto *t_66 = buffer.data(target + 66);
    auto *t_67 = buffer.data(target + 67);
    auto *t_68 = buffer.data(target + 68);
    auto *t_69 = buffer.data(target + 69);
    auto *t_70 = buffer.data(target + 70);
    auto *t_71 = buffer.data(target + 71);
    auto *t_72 = buffer.data(target + 72);
    auto *t_73 = buffer.data(target + 73);
    auto *t_74 = buffer.data(target + 74);
    auto *t_75 = buffer.data(target + 75);
    auto *t_76 = buffer.data(target + 76);
    auto *t_77 = buffer.data(target + 77);
    auto *t_78 = buffer.data(target + 78);
    auto *t_79 = buffer.data(target + 79);
    auto *t_80 = buffer.data(target + 80);
    auto *t_81 = buffer.data(target + 81);
    auto *t_82 = buffer.data(target + 82);
    auto *t_83 = buffer.data(target + 83);
    auto *t_84 = buffer.data(target + 84);
    auto *t_85 = buffer.data(target + 85);
    auto *t_86 = buffer.data(target + 86);
    auto *t_87 = buffer.data(target + 87);
    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);
    auto *t_90 = buffer.data(target + 90);
    auto *t_91 = buffer.data(target + 91);
    auto *t_92 = buffer.data(target + 92);
    auto *t_93 = buffer.data(target + 93);
    auto *t_94 = buffer.data(target + 94);
    auto *t_95 = buffer.data(target + 95);
    auto *t_96 = buffer.data(target + 96);
    auto *t_97 = buffer.data(target + 97);
    auto *t_98 = buffer.data(target + 98);
    auto *t_99 = buffer.data(target + 99);
    auto *t_100 = buffer.data(target + 100);
    auto *t_101 = buffer.data(target + 101);
    auto *t_102 = buffer.data(target + 102);
    auto *t_103 = buffer.data(target + 103);
    auto *t_104 = buffer.data(target + 104);
    auto *t_105 = buffer.data(target + 105);
    auto *t_106 = buffer.data(target + 106);
    auto *t_107 = buffer.data(target + 107);
    auto *t_108 = buffer.data(target + 108);
    auto *t_109 = buffer.data(target + 109);
    auto *t_110 = buffer.data(target + 110);
    auto *t_111 = buffer.data(target + 111);
    auto *t_112 = buffer.data(target + 112);
    auto *t_113 = buffer.data(target + 113);
    auto *t_114 = buffer.data(target + 114);
    auto *t_115 = buffer.data(target + 115);
    auto *t_116 = buffer.data(target + 116);
    auto *t_117 = buffer.data(target + 117);
    auto *t_118 = buffer.data(target + 118);
    auto *t_119 = buffer.data(target + 119);
    auto *t_120 = buffer.data(target + 120);
    auto *t_121 = buffer.data(target + 121);
    auto *t_122 = buffer.data(target + 122);
    auto *t_123 = buffer.data(target + 123);
    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);
    auto *t_135 = buffer.data(target + 135);
    auto *t_136 = buffer.data(target + 136);
    auto *t_137 = buffer.data(target + 137);
    auto *t_138 = buffer.data(target + 138);
    auto *t_139 = buffer.data(target + 139);
    auto *t_140 = buffer.data(target + 140);
    auto *t_141 = buffer.data(target + 141);
    auto *t_142 = buffer.data(target + 142);
    auto *t_143 = buffer.data(target + 143);
    auto *t_144 = buffer.data(target + 144);
    auto *t_145 = buffer.data(target + 145);
    auto *t_146 = buffer.data(target + 146);
    auto *t_147 = buffer.data(target + 147);
    auto *t_148 = buffer.data(target + 148);
    auto *t_149 = buffer.data(target + 149);

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_43 = buffer.data(dg + 43);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_47 = buffer.data(dg + 47);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_52 = buffer.data(dg + 52);
    const auto *dg_53 = buffer.data(dg + 53);
    const auto *dg_54 = buffer.data(dg + 54);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_56 = buffer.data(dg + 56);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_58 = buffer.data(dg + 58);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_60 = buffer.data(dg + 60);
    const auto *dg_61 = buffer.data(dg + 61);
    const auto *dg_62 = buffer.data(dg + 62);
    const auto *dg_63 = buffer.data(dg + 63);
    const auto *dg_64 = buffer.data(dg + 64);
    const auto *dg_65 = buffer.data(dg + 65);
    const auto *dg_66 = buffer.data(dg + 66);
    const auto *dg_67 = buffer.data(dg + 67);
    const auto *dg_68 = buffer.data(dg + 68);
    const auto *dg_69 = buffer.data(dg + 69);
    const auto *dg_70 = buffer.data(dg + 70);
    const auto *dg_71 = buffer.data(dg + 71);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_73 = buffer.data(dg + 73);
    const auto *dg_74 = buffer.data(dg + 74);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_76 = buffer.data(dg + 76);
    const auto *dg_77 = buffer.data(dg + 77);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_79 = buffer.data(dg + 79);
    const auto *dg_80 = buffer.data(dg + 80);
    const auto *dg_81 = buffer.data(dg + 81);
    const auto *dg_82 = buffer.data(dg + 82);
    const auto *dg_83 = buffer.data(dg + 83);
    const auto *dg_84 = buffer.data(dg + 84);
    const auto *dg_85 = buffer.data(dg + 85);
    const auto *dg_86 = buffer.data(dg + 86);
    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_88 = buffer.data(dg + 88);
    const auto *dg_89 = buffer.data(dg + 89);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_55 = buffer.data(dh + 55);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_64 = buffer.data(dh + 64);
    const auto *dh_65 = buffer.data(dh + 65);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_67 = buffer.data(dh + 67);
    const auto *dh_68 = buffer.data(dh + 68);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_70 = buffer.data(dh + 70);
    const auto *dh_71 = buffer.data(dh + 71);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_73 = buffer.data(dh + 73);
    const auto *dh_74 = buffer.data(dh + 74);
    const auto *dh_75 = buffer.data(dh + 75);
    const auto *dh_76 = buffer.data(dh + 76);
    const auto *dh_77 = buffer.data(dh + 77);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_79 = buffer.data(dh + 79);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_82 = buffer.data(dh + 82);
    const auto *dh_84 = buffer.data(dh + 84);
    const auto *dh_85 = buffer.data(dh + 85);
    const auto *dh_86 = buffer.data(dh + 86);
    const auto *dh_87 = buffer.data(dh + 87);
    const auto *dh_88 = buffer.data(dh + 88);
    const auto *dh_89 = buffer.data(dh + 89);
    const auto *dh_90 = buffer.data(dh + 90);
    const auto *dh_91 = buffer.data(dh + 91);
    const auto *dh_92 = buffer.data(dh + 92);
    const auto *dh_93 = buffer.data(dh + 93);
    const auto *dh_94 = buffer.data(dh + 94);
    const auto *dh_95 = buffer.data(dh + 95);
    const auto *dh_96 = buffer.data(dh + 96);
    const auto *dh_97 = buffer.data(dh + 97);
    const auto *dh_98 = buffer.data(dh + 98);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_106 = buffer.data(dh + 106);
    const auto *dh_107 = buffer.data(dh + 107);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_109 = buffer.data(dh + 109);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_111 = buffer.data(dh + 111);
    const auto *dh_112 = buffer.data(dh + 112);
    const auto *dh_113 = buffer.data(dh + 113);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_115 = buffer.data(dh + 115);
    const auto *dh_116 = buffer.data(dh + 116);
    const auto *dh_117 = buffer.data(dh + 117);
    const auto *dh_118 = buffer.data(dh + 118);
    const auto *dh_119 = buffer.data(dh + 119);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_121 = buffer.data(dh + 121);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_124 = buffer.data(dh + 124);
    const auto *dh_125 = buffer.data(dh + 125);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, dg_0, dg_1, dg_2, dg_3, dg_4, dh_0, \
                         dh_1, dh_2, dh_3, dh_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * dg_0[k]
                 + dh_0[k];

        t_1[k] = -ab_x[k] * dg_1[k]
                 + dh_1[k];

        t_2[k] = -ab_x[k] * dg_2[k]
                 + dh_2[k];

        t_3[k] = -ab_x[k] * dg_3[k]
                 + dh_3[k];

        t_4[k] = -ab_x[k] * dg_4[k]
                 + dh_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, dg_5, dg_6, dg_7, dg_8, dg_9, dh_5, \
                         dh_6, dh_7, dh_8, dh_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * dg_5[k]
                 + dh_5[k];

        t_6[k] = -ab_x[k] * dg_6[k]
                 + dh_6[k];

        t_7[k] = -ab_x[k] * dg_7[k]
                 + dh_7[k];

        t_8[k] = -ab_x[k] * dg_8[k]
                 + dh_8[k];

        t_9[k] = -ab_x[k] * dg_9[k]
                 + dh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, dg_10, dg_11, dg_12, dg_13, \
                         dg_14, dh_10, dh_11, dh_12, dh_13, dh_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * dg_10[k]
                  + dh_10[k];

        t_11[k] = -ab_x[k] * dg_11[k]
                  + dh_11[k];

        t_12[k] = -ab_x[k] * dg_12[k]
                  + dh_12[k];

        t_13[k] = -ab_x[k] * dg_13[k]
                  + dh_13[k];

        t_14[k] = -ab_x[k] * dg_14[k]
                  + dh_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, dg_15, dg_16, dg_17, dg_18, \
                         dg_19, dh_21, dh_22, dh_23, dh_24, dh_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * dg_15[k]
                  + dh_21[k];

        t_16[k] = -ab_x[k] * dg_16[k]
                  + dh_22[k];

        t_17[k] = -ab_x[k] * dg_17[k]
                  + dh_23[k];

        t_18[k] = -ab_x[k] * dg_18[k]
                  + dh_24[k];

        t_19[k] = -ab_x[k] * dg_19[k]
                  + dh_25[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, dg_20, dg_21, dg_22, dg_23, \
                         dg_24, dh_26, dh_27, dh_28, dh_29, dh_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * dg_20[k]
                  + dh_26[k];

        t_21[k] = -ab_x[k] * dg_21[k]
                  + dh_27[k];

        t_22[k] = -ab_x[k] * dg_22[k]
                  + dh_28[k];

        t_23[k] = -ab_x[k] * dg_23[k]
                  + dh_29[k];

        t_24[k] = -ab_x[k] * dg_24[k]
                  + dh_30[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, dg_25, dg_26, dg_27, dg_28, \
                         dg_29, dh_31, dh_32, dh_33, dh_34, dh_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * dg_25[k]
                  + dh_31[k];

        t_26[k] = -ab_x[k] * dg_26[k]
                  + dh_32[k];

        t_27[k] = -ab_x[k] * dg_27[k]
                  + dh_33[k];

        t_28[k] = -ab_x[k] * dg_28[k]
                  + dh_34[k];

        t_29[k] = -ab_x[k] * dg_29[k]
                  + dh_35[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, dg_30, dg_31, dg_32, dg_33, \
                         dg_34, dh_42, dh_43, dh_44, dh_45, dh_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * dg_30[k]
                  + dh_42[k];

        t_31[k] = -ab_x[k] * dg_31[k]
                  + dh_43[k];

        t_32[k] = -ab_x[k] * dg_32[k]
                  + dh_44[k];

        t_33[k] = -ab_x[k] * dg_33[k]
                  + dh_45[k];

        t_34[k] = -ab_x[k] * dg_34[k]
                  + dh_46[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, dg_35, dg_36, dg_37, dg_38, \
                         dg_39, dh_47, dh_48, dh_49, dh_50, dh_51 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * dg_35[k]
                  + dh_47[k];

        t_36[k] = -ab_x[k] * dg_36[k]
                  + dh_48[k];

        t_37[k] = -ab_x[k] * dg_37[k]
                  + dh_49[k];

        t_38[k] = -ab_x[k] * dg_38[k]
                  + dh_50[k];

        t_39[k] = -ab_x[k] * dg_39[k]
                  + dh_51[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, dg_40, dg_41, dg_42, dg_43, \
                         dg_44, dh_52, dh_53, dh_54, dh_55, dh_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * dg_40[k]
                  + dh_52[k];

        t_41[k] = -ab_x[k] * dg_41[k]
                  + dh_53[k];

        t_42[k] = -ab_x[k] * dg_42[k]
                  + dh_54[k];

        t_43[k] = -ab_x[k] * dg_43[k]
                  + dh_55[k];

        t_44[k] = -ab_x[k] * dg_44[k]
                  + dh_56[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, dg_45, dg_46, dg_47, dg_48, \
                         dg_49, dh_63, dh_64, dh_65, dh_66, dh_67 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_x[k] * dg_45[k]
                  + dh_63[k];

        t_46[k] = -ab_x[k] * dg_46[k]
                  + dh_64[k];

        t_47[k] = -ab_x[k] * dg_47[k]
                  + dh_65[k];

        t_48[k] = -ab_x[k] * dg_48[k]
                  + dh_66[k];

        t_49[k] = -ab_x[k] * dg_49[k]
                  + dh_67[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, dg_50, dg_51, dg_52, dg_53, \
                         dg_54, dh_68, dh_69, dh_70, dh_71, dh_72 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_x[k] * dg_50[k]
                  + dh_68[k];

        t_51[k] = -ab_x[k] * dg_51[k]
                  + dh_69[k];

        t_52[k] = -ab_x[k] * dg_52[k]
                  + dh_70[k];

        t_53[k] = -ab_x[k] * dg_53[k]
                  + dh_71[k];

        t_54[k] = -ab_x[k] * dg_54[k]
                  + dh_72[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, dg_55, dg_56, dg_57, dg_58, \
                         dg_59, dh_73, dh_74, dh_75, dh_76, dh_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_x[k] * dg_55[k]
                  + dh_73[k];

        t_56[k] = -ab_x[k] * dg_56[k]
                  + dh_74[k];

        t_57[k] = -ab_x[k] * dg_57[k]
                  + dh_75[k];

        t_58[k] = -ab_x[k] * dg_58[k]
                  + dh_76[k];

        t_59[k] = -ab_x[k] * dg_59[k]
                  + dh_77[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, dg_60, dg_61, dg_62, dg_63, \
                         dg_64, dh_84, dh_85, dh_86, dh_87, dh_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_x[k] * dg_60[k]
                  + dh_84[k];

        t_61[k] = -ab_x[k] * dg_61[k]
                  + dh_85[k];

        t_62[k] = -ab_x[k] * dg_62[k]
                  + dh_86[k];

        t_63[k] = -ab_x[k] * dg_63[k]
                  + dh_87[k];

        t_64[k] = -ab_x[k] * dg_64[k]
                  + dh_88[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, dg_65, dg_66, dg_67, dg_68, \
                         dg_69, dh_89, dh_90, dh_91, dh_92, dh_93 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = -ab_x[k] * dg_65[k]
                  + dh_89[k];

        t_66[k] = -ab_x[k] * dg_66[k]
                  + dh_90[k];

        t_67[k] = -ab_x[k] * dg_67[k]
                  + dh_91[k];

        t_68[k] = -ab_x[k] * dg_68[k]
                  + dh_92[k];

        t_69[k] = -ab_x[k] * dg_69[k]
                  + dh_93[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, dg_70, dg_71, dg_72, dg_73, \
                         dg_74, dh_94, dh_95, dh_96, dh_97, dh_98 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = -ab_x[k] * dg_70[k]
                  + dh_94[k];

        t_71[k] = -ab_x[k] * dg_71[k]
                  + dh_95[k];

        t_72[k] = -ab_x[k] * dg_72[k]
                  + dh_96[k];

        t_73[k] = -ab_x[k] * dg_73[k]
                  + dh_97[k];

        t_74[k] = -ab_x[k] * dg_74[k]
                  + dh_98[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, dg_75, dg_76, dg_77, dg_78, \
                         dg_79, dh_105, dh_106, dh_107, dh_108, \
                         dh_109 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = -ab_x[k] * dg_75[k]
                  + dh_105[k];

        t_76[k] = -ab_x[k] * dg_76[k]
                  + dh_106[k];

        t_77[k] = -ab_x[k] * dg_77[k]
                  + dh_107[k];

        t_78[k] = -ab_x[k] * dg_78[k]
                  + dh_108[k];

        t_79[k] = -ab_x[k] * dg_79[k]
                  + dh_109[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, dg_80, dg_81, dg_82, dg_83, \
                         dg_84, dh_110, dh_111, dh_112, dh_113, \
                         dh_114 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = -ab_x[k] * dg_80[k]
                  + dh_110[k];

        t_81[k] = -ab_x[k] * dg_81[k]
                  + dh_111[k];

        t_82[k] = -ab_x[k] * dg_82[k]
                  + dh_112[k];

        t_83[k] = -ab_x[k] * dg_83[k]
                  + dh_113[k];

        t_84[k] = -ab_x[k] * dg_84[k]
                  + dh_114[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, dg_85, dg_86, dg_87, dg_88, \
                         dg_89, dh_115, dh_116, dh_117, dh_118, \
                         dh_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = -ab_x[k] * dg_85[k]
                  + dh_115[k];

        t_86[k] = -ab_x[k] * dg_86[k]
                  + dh_116[k];

        t_87[k] = -ab_x[k] * dg_87[k]
                  + dh_117[k];

        t_88[k] = -ab_x[k] * dg_88[k]
                  + dh_118[k];

        t_89[k] = -ab_x[k] * dg_89[k]
                  + dh_119[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_y, dg_45, dg_46, dg_47, dg_48, \
                         dg_49, dh_64, dh_66, dh_67, dh_69, dh_70 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = -ab_y[k] * dg_45[k]
                  + dh_64[k];

        t_91[k] = -ab_y[k] * dg_46[k]
                  + dh_66[k];

        t_92[k] = -ab_y[k] * dg_47[k]
                  + dh_67[k];

        t_93[k] = -ab_y[k] * dg_48[k]
                  + dh_69[k];

        t_94[k] = -ab_y[k] * dg_49[k]
                  + dh_70[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_y, dg_50, dg_51, dg_52, dg_53, \
                         dg_54, dh_71, dh_73, dh_74, dh_75, dh_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = -ab_y[k] * dg_50[k]
                  + dh_71[k];

        t_96[k] = -ab_y[k] * dg_51[k]
                  + dh_73[k];

        t_97[k] = -ab_y[k] * dg_52[k]
                  + dh_74[k];

        t_98[k] = -ab_y[k] * dg_53[k]
                  + dh_75[k];

        t_99[k] = -ab_y[k] * dg_54[k]
                  + dh_76[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_y, dg_55, dg_56, dg_57, dg_58, \
                         dg_59, dh_78, dh_79, dh_80, dh_81, dh_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = -ab_y[k] * dg_55[k]
                   + dh_78[k];

        t_101[k] = -ab_y[k] * dg_56[k]
                   + dh_79[k];

        t_102[k] = -ab_y[k] * dg_57[k]
                   + dh_80[k];

        t_103[k] = -ab_y[k] * dg_58[k]
                   + dh_81[k];

        t_104[k] = -ab_y[k] * dg_59[k]
                   + dh_82[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_y, dg_60, dg_61, dg_62, dg_63, \
                         dg_64, dh_85, dh_87, dh_88, dh_90, dh_91 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = -ab_y[k] * dg_60[k]
                   + dh_85[k];

        t_106[k] = -ab_y[k] * dg_61[k]
                   + dh_87[k];

        t_107[k] = -ab_y[k] * dg_62[k]
                   + dh_88[k];

        t_108[k] = -ab_y[k] * dg_63[k]
                   + dh_90[k];

        t_109[k] = -ab_y[k] * dg_64[k]
                   + dh_91[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_y, dg_65, dg_66, dg_67, dg_68, \
                         dg_69, dh_92, dh_94, dh_95, dh_96, dh_97 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = -ab_y[k] * dg_65[k]
                   + dh_92[k];

        t_111[k] = -ab_y[k] * dg_66[k]
                   + dh_94[k];

        t_112[k] = -ab_y[k] * dg_67[k]
                   + dh_95[k];

        t_113[k] = -ab_y[k] * dg_68[k]
                   + dh_96[k];

        t_114[k] = -ab_y[k] * dg_69[k]
                   + dh_97[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_y, dg_70, dg_71, dg_72, dg_73, \
                         dg_74, dh_99, dh_100, dh_101, dh_102, dh_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = -ab_y[k] * dg_70[k]
                   + dh_99[k];

        t_116[k] = -ab_y[k] * dg_71[k]
                   + dh_100[k];

        t_117[k] = -ab_y[k] * dg_72[k]
                   + dh_101[k];

        t_118[k] = -ab_y[k] * dg_73[k]
                   + dh_102[k];

        t_119[k] = -ab_y[k] * dg_74[k]
                   + dh_103[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_y, dg_75, dg_76, dg_77, dg_78, \
                         dg_79, dh_106, dh_108, dh_109, dh_111, \
                         dh_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = -ab_y[k] * dg_75[k]
                   + dh_106[k];

        t_121[k] = -ab_y[k] * dg_76[k]
                   + dh_108[k];

        t_122[k] = -ab_y[k] * dg_77[k]
                   + dh_109[k];

        t_123[k] = -ab_y[k] * dg_78[k]
                   + dh_111[k];

        t_124[k] = -ab_y[k] * dg_79[k]
                   + dh_112[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_y, dg_80, dg_81, dg_82, dg_83, \
                         dg_84, dh_113, dh_115, dh_116, dh_117, \
                         dh_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = -ab_y[k] * dg_80[k]
                   + dh_113[k];

        t_126[k] = -ab_y[k] * dg_81[k]
                   + dh_115[k];

        t_127[k] = -ab_y[k] * dg_82[k]
                   + dh_116[k];

        t_128[k] = -ab_y[k] * dg_83[k]
                   + dh_117[k];

        t_129[k] = -ab_y[k] * dg_84[k]
                   + dh_118[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_y, dg_85, dg_86, dg_87, dg_88, \
                         dg_89, dh_120, dh_121, dh_122, dh_123, \
                         dh_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_130[k] = -ab_y[k] * dg_85[k]
                   + dh_120[k];

        t_131[k] = -ab_y[k] * dg_86[k]
                   + dh_121[k];

        t_132[k] = -ab_y[k] * dg_87[k]
                   + dh_122[k];

        t_133[k] = -ab_y[k] * dg_88[k]
                   + dh_123[k];

        t_134[k] = -ab_y[k] * dg_89[k]
                   + dh_124[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_z, dg_75, dg_76, dg_77, dg_78, \
                         dg_79, dh_107, dh_109, dh_110, dh_112, \
                         dh_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = -ab_z[k] * dg_75[k]
                   + dh_107[k];

        t_136[k] = -ab_z[k] * dg_76[k]
                   + dh_109[k];

        t_137[k] = -ab_z[k] * dg_77[k]
                   + dh_110[k];

        t_138[k] = -ab_z[k] * dg_78[k]
                   + dh_112[k];

        t_139[k] = -ab_z[k] * dg_79[k]
                   + dh_113[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_z, dg_80, dg_81, dg_82, dg_83, \
                         dg_84, dh_114, dh_116, dh_117, dh_118, \
                         dh_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_140[k] = -ab_z[k] * dg_80[k]
                   + dh_114[k];

        t_141[k] = -ab_z[k] * dg_81[k]
                   + dh_116[k];

        t_142[k] = -ab_z[k] * dg_82[k]
                   + dh_117[k];

        t_143[k] = -ab_z[k] * dg_83[k]
                   + dh_118[k];

        t_144[k] = -ab_z[k] * dg_84[k]
                   + dh_119[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_z, dg_85, dg_86, dg_87, dg_88, \
                         dg_89, dh_121, dh_122, dh_123, dh_124, \
                         dh_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_145[k] = -ab_z[k] * dg_85[k]
                   + dh_121[k];

        t_146[k] = -ab_z[k] * dg_86[k]
                   + dh_122[k];

        t_147[k] = -ab_z[k] * dg_87[k]
                   + dh_123[k];

        t_148[k] = -ab_z[k] * dg_88[k]
                   + dh_124[k];

        t_149[k] = -ab_z[k] * dg_89[k]
                   + dh_125[k];
    }
}

}  // namespace simdtrf
