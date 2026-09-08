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


#include "SimdTransferGF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_gf_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t gd, const size_t hd,
                   const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.875 * std::sqrt(14.0);
    const auto f_1 = 0.625 * std::sqrt(14.0);
    const auto f_2 = 2.5 * std::sqrt(21.0);
    const auto f_3 = 0.125 * std::sqrt(210.0);
    const auto f_4 = 0.5 * std::sqrt(210.0);
    const auto f_5 = 0.75 * std::sqrt(35.0);
    const auto f_6 = 0.5 * std::sqrt(35.0);
    const auto f_7 = 1.25 * std::sqrt(21.0);
    const auto f_8 = 5.625 * std::sqrt(7.0);
    const auto f_9 = 1.875 * std::sqrt(7.0);
    const auto f_10 = 0.625 * std::sqrt(7.0);
    const auto f_11 = 3.75 * std::sqrt(42.0);
    const auto f_12 = 1.25 * std::sqrt(42.0);
    const auto f_13 = 0.375 * std::sqrt(105.0);
    const auto f_14 = 1.5 * std::sqrt(105.0);
    const auto f_15 = 0.125 * std::sqrt(105.0);
    const auto f_16 = 0.5 * std::sqrt(105.0);
    const auto f_17 = 1.125 * std::sqrt(70.0);
    const auto f_18 = 0.75 * std::sqrt(70.0);
    const auto f_19 = 0.375 * std::sqrt(70.0);
    const auto f_20 = 0.25 * std::sqrt(70.0);
    const auto f_21 = 1.875 * std::sqrt(42.0);
    const auto f_22 = 0.625 * std::sqrt(42.0);
    const auto f_23 = 1.875 * std::sqrt(2.0);
    const auto f_24 = 0.625 * std::sqrt(2.0);
    const auto f_25 = 11.25 * std::sqrt(2.0);
    const auto f_26 = 3.75 * std::sqrt(2.0);
    const auto f_27 = 2.5 * std::sqrt(3.0);
    const auto f_28 = 15.0 * std::sqrt(3.0);
    const auto f_29 = 0.125 * std::sqrt(30.0);
    const auto f_30 = 0.5 * std::sqrt(30.0);
    const auto f_31 = 0.75 * std::sqrt(30.0);
    const auto f_32 = 3.0 * std::sqrt(30.0);
    const auto f_33 = 0.75 * std::sqrt(5.0);
    const auto f_34 = 0.5 * std::sqrt(5.0);
    const auto f_35 = 4.5 * std::sqrt(5.0);
    const auto f_36 = 3.0 * std::sqrt(5.0);
    const auto f_37 = 1.25 * std::sqrt(3.0);
    const auto f_38 = 7.5 * std::sqrt(3.0);
    const auto f_39 = 3.75 * std::sqrt(6.0);
    const auto f_40 = 5.0 * std::sqrt(6.0);
    const auto f_41 = 0.375 * std::sqrt(15.0);
    const auto f_42 = 1.5 * std::sqrt(15.0);
    const auto f_43 = 0.5 * std::sqrt(15.0);
    const auto f_44 = 2.0 * std::sqrt(15.0);
    const auto f_45 = 1.125 * std::sqrt(10.0);
    const auto f_46 = 0.75 * std::sqrt(10.0);
    const auto f_47 = 1.5 * std::sqrt(10.0);
    const auto f_48 = std::sqrt(10.0);
    const auto f_49 = 1.875 * std::sqrt(6.0);
    const auto f_50 = 2.5 * std::sqrt(6.0);
    const auto f_51 = 0.28125 * std::sqrt(10.0);
    const auto f_52 = 0.09375 * std::sqrt(10.0);
    const auto f_53 = 0.5625 * std::sqrt(10.0);
    const auto f_54 = 0.1875 * std::sqrt(10.0);
    const auto f_55 = 2.25 * std::sqrt(10.0);
    const auto f_56 = 0.25 * std::sqrt(10.0);
    const auto f_57 = 0.75 * std::sqrt(15.0);
    const auto f_58 = 3.0 * std::sqrt(15.0);
    const auto f_59 = std::sqrt(15.0);
    const auto f_60 = 0.09375 * std::sqrt(6.0);
    const auto f_61 = 0.375 * std::sqrt(6.0);
    const auto f_62 = 0.1875 * std::sqrt(6.0);
    const auto f_63 = 0.75 * std::sqrt(6.0);
    const auto f_64 = 3.0 * std::sqrt(6.0);
    const auto f_65 = 0.25 * std::sqrt(6.0);
    const auto f_66 = std::sqrt(6.0);
    const auto f_67 = 0.1875 * std::sqrt(15.0);
    const auto f_68 = 0.9375 * std::sqrt(2.0);
    const auto f_69 = 0.3125 * std::sqrt(2.0);
    const auto f_70 = 5.625 * std::sqrt(2.0);
    const auto f_71 = 0.0625 * std::sqrt(30.0);
    const auto f_72 = 0.25 * std::sqrt(30.0);
    const auto f_73 = 0.375 * std::sqrt(30.0);
    const auto f_74 = 1.5 * std::sqrt(30.0);
    const auto f_75 = 0.375 * std::sqrt(5.0);
    const auto f_76 = 0.25 * std::sqrt(5.0);
    const auto f_77 = 2.25 * std::sqrt(5.0);
    const auto f_78 = 1.5 * std::sqrt(5.0);
    const auto f_79 = 0.625 * std::sqrt(3.0);
    const auto f_80 = 3.75 * std::sqrt(3.0);
    const auto f_81 = 0.46875 * std::sqrt(14.0);
    const auto f_82 = 0.15625 * std::sqrt(14.0);
    const auto f_83 = 2.8125 * std::sqrt(14.0);
    const auto f_84 = 0.9375 * std::sqrt(14.0);
    const auto f_85 = 0.625 * std::sqrt(21.0);
    const auto f_86 = 3.75 * std::sqrt(21.0);
    const auto f_87 = 0.03125 * std::sqrt(210.0);
    const auto f_88 = 0.1875 * std::sqrt(210.0);
    const auto f_89 = 0.75 * std::sqrt(210.0);
    const auto f_90 = 0.1875 * std::sqrt(35.0);
    const auto f_91 = 0.125 * std::sqrt(35.0);
    const auto f_92 = 1.125 * std::sqrt(35.0);
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

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_44 = buffer.data(gd + 44);
    const auto *gd_45 = buffer.data(gd + 45);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);
    const auto *gd_48 = buffer.data(gd + 48);
    const auto *gd_49 = buffer.data(gd + 49);
    const auto *gd_50 = buffer.data(gd + 50);
    const auto *gd_51 = buffer.data(gd + 51);
    const auto *gd_52 = buffer.data(gd + 52);
    const auto *gd_53 = buffer.data(gd + 53);
    const auto *gd_54 = buffer.data(gd + 54);
    const auto *gd_55 = buffer.data(gd + 55);
    const auto *gd_56 = buffer.data(gd + 56);
    const auto *gd_57 = buffer.data(gd + 57);
    const auto *gd_58 = buffer.data(gd + 58);
    const auto *gd_59 = buffer.data(gd + 59);
    const auto *gd_60 = buffer.data(gd + 60);
    const auto *gd_61 = buffer.data(gd + 61);
    const auto *gd_62 = buffer.data(gd + 62);
    const auto *gd_63 = buffer.data(gd + 63);
    const auto *gd_64 = buffer.data(gd + 64);
    const auto *gd_65 = buffer.data(gd + 65);
    const auto *gd_66 = buffer.data(gd + 66);
    const auto *gd_67 = buffer.data(gd + 67);
    const auto *gd_68 = buffer.data(gd + 68);
    const auto *gd_69 = buffer.data(gd + 69);
    const auto *gd_70 = buffer.data(gd + 70);
    const auto *gd_71 = buffer.data(gd + 71);
    const auto *gd_72 = buffer.data(gd + 72);
    const auto *gd_73 = buffer.data(gd + 73);
    const auto *gd_74 = buffer.data(gd + 74);
    const auto *gd_75 = buffer.data(gd + 75);
    const auto *gd_76 = buffer.data(gd + 76);
    const auto *gd_77 = buffer.data(gd + 77);
    const auto *gd_78 = buffer.data(gd + 78);
    const auto *gd_79 = buffer.data(gd + 79);
    const auto *gd_80 = buffer.data(gd + 80);
    const auto *gd_81 = buffer.data(gd + 81);
    const auto *gd_82 = buffer.data(gd + 82);
    const auto *gd_83 = buffer.data(gd + 83);
    const auto *gd_84 = buffer.data(gd + 84);
    const auto *gd_85 = buffer.data(gd + 85);
    const auto *gd_86 = buffer.data(gd + 86);
    const auto *gd_87 = buffer.data(gd + 87);
    const auto *gd_88 = buffer.data(gd + 88);
    const auto *gd_89 = buffer.data(gd + 89);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_62 = buffer.data(hd + 62);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_64 = buffer.data(hd + 64);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_67 = buffer.data(hd + 67);
    const auto *hd_68 = buffer.data(hd + 68);
    const auto *hd_69 = buffer.data(hd + 69);
    const auto *hd_70 = buffer.data(hd + 70);
    const auto *hd_71 = buffer.data(hd + 71);
    const auto *hd_72 = buffer.data(hd + 72);
    const auto *hd_73 = buffer.data(hd + 73);
    const auto *hd_74 = buffer.data(hd + 74);
    const auto *hd_75 = buffer.data(hd + 75);
    const auto *hd_76 = buffer.data(hd + 76);
    const auto *hd_77 = buffer.data(hd + 77);
    const auto *hd_78 = buffer.data(hd + 78);
    const auto *hd_79 = buffer.data(hd + 79);
    const auto *hd_80 = buffer.data(hd + 80);
    const auto *hd_81 = buffer.data(hd + 81);
    const auto *hd_82 = buffer.data(hd + 82);
    const auto *hd_83 = buffer.data(hd + 83);
    const auto *hd_84 = buffer.data(hd + 84);
    const auto *hd_85 = buffer.data(hd + 85);
    const auto *hd_86 = buffer.data(hd + 86);
    const auto *hd_87 = buffer.data(hd + 87);
    const auto *hd_88 = buffer.data(hd + 88);
    const auto *hd_89 = buffer.data(hd + 89);
    const auto *hd_93 = buffer.data(hd + 93);
    const auto *hd_94 = buffer.data(hd + 94);
    const auto *hd_95 = buffer.data(hd + 95);
    const auto *hd_99 = buffer.data(hd + 99);
    const auto *hd_100 = buffer.data(hd + 100);
    const auto *hd_101 = buffer.data(hd + 101);
    const auto *hd_105 = buffer.data(hd + 105);
    const auto *hd_106 = buffer.data(hd + 106);
    const auto *hd_107 = buffer.data(hd + 107);
    const auto *hd_111 = buffer.data(hd + 111);
    const auto *hd_112 = buffer.data(hd + 112);
    const auto *hd_113 = buffer.data(hd + 113);
    const auto *hd_117 = buffer.data(hd + 117);
    const auto *hd_118 = buffer.data(hd + 118);
    const auto *hd_119 = buffer.data(hd + 119);
    const auto *hd_125 = buffer.data(hd + 125);

#pragma omp simd aligned(ab_x, ab_y, gd_7, gd_9, gd_10, gd_37, gd_39, gd_40, hd_7, hd_10, \
                         hd_21, hd_37, hd_40, hd_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ab_x[k] * gd_7[k]
                 - f_1 * ab_y[k] * gd_9[k]
                 - f_0 * ab_x[k] * gd_37[k]
                 + f_1 * ab_y[k] * gd_39[k]
                 + f_0 * hd_7[k]
                 - f_1 * hd_21[k]
                 - f_0 * hd_37[k]
                 + f_1 * hd_63[k];

        g_1[k] = f_2 * ab_x[k] * gd_10[k]
                 - f_2 * ab_x[k] * gd_40[k]
                 + f_2 * hd_10[k]
                 - f_2 * hd_40[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_7, gd_9, gd_11, gd_37, gd_39, gd_41, hd_7, hd_21, \
                         hd_23, hd_37, hd_63, hd_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_3 * ab_x[k] * gd_7[k]
                 - f_3 * ab_y[k] * gd_9[k]
                 + f_4 * ab_y[k] * gd_11[k]
                 + f_3 * ab_x[k] * gd_37[k]
                 + f_3 * ab_y[k] * gd_39[k]
                 - f_4 * ab_y[k] * gd_41[k]
                 - f_3 * hd_7[k]
                 - f_3 * hd_21[k]
                 + f_4 * hd_23[k]
                 + f_3 * hd_37[k]
                 + f_3 * hd_63[k]
                 - f_4 * hd_65[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gd_8, gd_10, gd_11, gd_38, gd_40, gd_41, hd_8, \
                         hd_22, hd_29, hd_38, hd_64, hd_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_5 * ab_x[k] * gd_8[k]
                 - f_5 * ab_y[k] * gd_10[k]
                 + f_6 * ab_z[k] * gd_11[k]
                 + f_5 * ab_x[k] * gd_38[k]
                 + f_5 * ab_y[k] * gd_40[k]
                 - f_6 * ab_z[k] * gd_41[k]
                 - f_5 * hd_8[k]
                 - f_5 * hd_22[k]
                 + f_6 * hd_29[k]
                 + f_5 * hd_38[k]
                 + f_5 * hd_64[k]
                 - f_6 * hd_71[k];
    }

#pragma omp simd aligned(ab_x, gd_6, gd_9, gd_11, gd_36, gd_39, gd_41, hd_6, hd_9, hd_11, \
                         hd_36, hd_39, hd_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_3 * ab_x[k] * gd_6[k]
                 - f_3 * ab_x[k] * gd_9[k]
                 + f_4 * ab_x[k] * gd_11[k]
                 + f_3 * ab_x[k] * gd_36[k]
                 + f_3 * ab_x[k] * gd_39[k]
                 - f_4 * ab_x[k] * gd_41[k]
                 - f_3 * hd_6[k]
                 - f_3 * hd_9[k]
                 + f_4 * hd_11[k]
                 + f_3 * hd_36[k]
                 + f_3 * hd_39[k]
                 - f_4 * hd_41[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_8, gd_10, gd_38, gd_40, hd_8, hd_22, hd_38, \
                         hd_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_7 * ab_x[k] * gd_8[k]
                 - f_7 * ab_y[k] * gd_10[k]
                 - f_7 * ab_x[k] * gd_38[k]
                 + f_7 * ab_y[k] * gd_40[k]
                 + f_7 * hd_8[k]
                 - f_7 * hd_22[k]
                 - f_7 * hd_38[k]
                 + f_7 * hd_64[k];
    }

#pragma omp simd aligned(ab_x, gd_6, gd_9, gd_36, gd_39, hd_6, hd_9, hd_36, \
                         hd_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_1 * ab_x[k] * gd_6[k]
                 - f_0 * ab_x[k] * gd_9[k]
                 - f_1 * ab_x[k] * gd_36[k]
                 + f_0 * ab_x[k] * gd_39[k]
                 + f_1 * hd_6[k]
                 - f_0 * hd_9[k]
                 - f_1 * hd_36[k]
                 + f_0 * hd_39[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_25, gd_27, gd_28, gd_67, gd_69, gd_70, hd_25, hd_28, \
                         hd_45, hd_67, hd_70, hd_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_8 * ab_x[k] * gd_25[k]
                 - f_9 * ab_y[k] * gd_27[k]
                 - f_9 * ab_x[k] * gd_67[k]
                 + f_10 * ab_y[k] * gd_69[k]
                 + f_8 * hd_25[k]
                 - f_9 * hd_45[k]
                 - f_9 * hd_67[k]
                 + f_10 * hd_99[k];

        g_8[k] = f_11 * ab_x[k] * gd_28[k]
                 - f_12 * ab_x[k] * gd_70[k]
                 + f_11 * hd_28[k]
                 - f_12 * hd_70[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_25, gd_27, gd_29, gd_67, gd_69, gd_71, hd_25, hd_45, \
                         hd_47, hd_67, hd_99, hd_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_13 * ab_x[k] * gd_25[k]
                 - f_13 * ab_y[k] * gd_27[k]
                 + f_14 * ab_y[k] * gd_29[k]
                 + f_15 * ab_x[k] * gd_67[k]
                 + f_15 * ab_y[k] * gd_69[k]
                 - f_16 * ab_y[k] * gd_71[k]
                 - f_13 * hd_25[k]
                 - f_13 * hd_45[k]
                 + f_14 * hd_47[k]
                 + f_15 * hd_67[k]
                 + f_15 * hd_99[k]
                 - f_16 * hd_101[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gd_26, gd_28, gd_29, gd_68, gd_70, gd_71, hd_26, \
                         hd_46, hd_53, hd_68, hd_100, hd_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_17 * ab_x[k] * gd_26[k]
                  - f_17 * ab_y[k] * gd_28[k]
                  + f_18 * ab_z[k] * gd_29[k]
                  + f_19 * ab_x[k] * gd_68[k]
                  + f_19 * ab_y[k] * gd_70[k]
                  - f_20 * ab_z[k] * gd_71[k]
                  - f_17 * hd_26[k]
                  - f_17 * hd_46[k]
                  + f_18 * hd_53[k]
                  + f_19 * hd_68[k]
                  + f_19 * hd_100[k]
                  - f_20 * hd_107[k];
    }

#pragma omp simd aligned(ab_x, gd_24, gd_27, gd_29, gd_66, gd_69, gd_71, hd_24, hd_27, hd_29, \
                         hd_66, hd_69, hd_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_13 * ab_x[k] * gd_24[k]
                  - f_13 * ab_x[k] * gd_27[k]
                  + f_14 * ab_x[k] * gd_29[k]
                  + f_15 * ab_x[k] * gd_66[k]
                  + f_15 * ab_x[k] * gd_69[k]
                  - f_16 * ab_x[k] * gd_71[k]
                  - f_13 * hd_24[k]
                  - f_13 * hd_27[k]
                  + f_14 * hd_29[k]
                  + f_15 * hd_66[k]
                  + f_15 * hd_69[k]
                  - f_16 * hd_71[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_26, gd_28, gd_68, gd_70, hd_26, hd_46, hd_68, \
                         hd_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_21 * ab_x[k] * gd_26[k]
                  - f_21 * ab_y[k] * gd_28[k]
                  - f_22 * ab_x[k] * gd_68[k]
                  + f_22 * ab_y[k] * gd_70[k]
                  + f_21 * hd_26[k]
                  - f_21 * hd_46[k]
                  - f_22 * hd_68[k]
                  + f_22 * hd_100[k];
    }

#pragma omp simd aligned(ab_x, gd_24, gd_27, gd_66, gd_69, hd_24, hd_27, hd_66, \
                         hd_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_9 * ab_x[k] * gd_24[k]
                  - f_8 * ab_x[k] * gd_27[k]
                  - f_10 * ab_x[k] * gd_66[k]
                  + f_9 * ab_x[k] * gd_69[k]
                  + f_9 * hd_24[k]
                  - f_8 * hd_27[k]
                  - f_10 * hd_66[k]
                  + f_9 * hd_69[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_7, gd_9, gd_37, gd_39, gd_49, gd_51, hd_7, hd_21, \
                         hd_37, hd_49, hd_63, hd_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_23 * ab_x[k] * gd_7[k]
                  + f_24 * ab_y[k] * gd_9[k]
                  - f_23 * ab_x[k] * gd_37[k]
                  + f_24 * ab_y[k] * gd_39[k]
                  + f_25 * ab_x[k] * gd_49[k]
                  - f_26 * ab_y[k] * gd_51[k]
                  - f_23 * hd_7[k]
                  + f_24 * hd_21[k]
                  - f_23 * hd_37[k]
                  + f_25 * hd_49[k]
                  + f_24 * hd_63[k]
                  - f_26 * hd_75[k];
    }

#pragma omp simd aligned(ab_x, gd_10, gd_40, gd_52, hd_10, hd_40, \
                         hd_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_27 * ab_x[k] * gd_10[k]
                  - f_27 * ab_x[k] * gd_40[k]
                  + f_28 * ab_x[k] * gd_52[k]
                  - f_27 * hd_10[k]
                  - f_27 * hd_40[k]
                  + f_28 * hd_52[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_7, gd_9, gd_11, gd_37, gd_39, gd_41, gd_49, gd_51, \
                         gd_53, hd_7, hd_21, hd_23, hd_37, hd_49, hd_63, hd_65, hd_75, \
                         hd_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_29 * ab_x[k] * gd_7[k]
                  + f_29 * ab_y[k] * gd_9[k]
                  - f_30 * ab_y[k] * gd_11[k]
                  + f_29 * ab_x[k] * gd_37[k]
                  + f_29 * ab_y[k] * gd_39[k]
                  - f_30 * ab_y[k] * gd_41[k]
                  - f_31 * ab_x[k] * gd_49[k]
                  - f_31 * ab_y[k] * gd_51[k]
                  + f_32 * ab_y[k] * gd_53[k]
                  + f_29 * hd_7[k]
                  + f_29 * hd_21[k]
                  - f_30 * hd_23[k]
                  + f_29 * hd_37[k]
                  - f_31 * hd_49[k]
                  + f_29 * hd_63[k]
                  - f_30 * hd_65[k]
                  - f_31 * hd_75[k]
                  + f_32 * hd_77[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gd_8, gd_10, gd_11, gd_38, gd_40, gd_41, gd_50, \
                         gd_52, gd_53, hd_8, hd_22, hd_29, hd_38, hd_50, hd_64, hd_71, hd_76, \
                         hd_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_33 * ab_x[k] * gd_8[k]
                  + f_33 * ab_y[k] * gd_10[k]
                  - f_34 * ab_z[k] * gd_11[k]
                  + f_33 * ab_x[k] * gd_38[k]
                  + f_33 * ab_y[k] * gd_40[k]
                  - f_34 * ab_z[k] * gd_41[k]
                  - f_35 * ab_x[k] * gd_50[k]
                  - f_35 * ab_y[k] * gd_52[k]
                  + f_36 * ab_z[k] * gd_53[k]
                  + f_33 * hd_8[k]
                  + f_33 * hd_22[k]
                  - f_34 * hd_29[k]
                  + f_33 * hd_38[k]
                  - f_35 * hd_50[k]
                  + f_33 * hd_64[k]
                  - f_34 * hd_71[k]
                  - f_35 * hd_76[k]
                  + f_36 * hd_83[k];
    }

#pragma omp simd aligned(ab_x, gd_6, gd_9, gd_11, gd_36, gd_39, gd_41, gd_48, gd_51, gd_53, \
                         hd_6, hd_9, hd_11, hd_36, hd_39, hd_41, hd_48, hd_51, \
                         hd_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_29 * ab_x[k] * gd_6[k]
                  + f_29 * ab_x[k] * gd_9[k]
                  - f_30 * ab_x[k] * gd_11[k]
                  + f_29 * ab_x[k] * gd_36[k]
                  + f_29 * ab_x[k] * gd_39[k]
                  - f_30 * ab_x[k] * gd_41[k]
                  - f_31 * ab_x[k] * gd_48[k]
                  - f_31 * ab_x[k] * gd_51[k]
                  + f_32 * ab_x[k] * gd_53[k]
                  + f_29 * hd_6[k]
                  + f_29 * hd_9[k]
                  - f_30 * hd_11[k]
                  + f_29 * hd_36[k]
                  + f_29 * hd_39[k]
                  - f_30 * hd_41[k]
                  - f_31 * hd_48[k]
                  - f_31 * hd_51[k]
                  + f_32 * hd_53[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_8, gd_10, gd_38, gd_40, gd_50, gd_52, hd_8, hd_22, \
                         hd_38, hd_50, hd_64, hd_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_37 * ab_x[k] * gd_8[k]
                  + f_37 * ab_y[k] * gd_10[k]
                  - f_37 * ab_x[k] * gd_38[k]
                  + f_37 * ab_y[k] * gd_40[k]
                  + f_38 * ab_x[k] * gd_50[k]
                  - f_38 * ab_y[k] * gd_52[k]
                  - f_37 * hd_8[k]
                  + f_37 * hd_22[k]
                  - f_37 * hd_38[k]
                  + f_38 * hd_50[k]
                  + f_37 * hd_64[k]
                  - f_38 * hd_76[k];
    }

#pragma omp simd aligned(ab_x, gd_6, gd_9, gd_36, gd_39, gd_48, gd_51, hd_6, hd_9, hd_36, \
                         hd_39, hd_48, hd_51 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_24 * ab_x[k] * gd_6[k]
                  + f_23 * ab_x[k] * gd_9[k]
                  - f_24 * ab_x[k] * gd_36[k]
                  + f_23 * ab_x[k] * gd_39[k]
                  + f_26 * ab_x[k] * gd_48[k]
                  - f_25 * ab_x[k] * gd_51[k]
                  - f_24 * hd_6[k]
                  + f_23 * hd_9[k]
                  - f_24 * hd_36[k]
                  + f_23 * hd_39[k]
                  + f_26 * hd_48[k]
                  - f_25 * hd_51[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_25, gd_27, gd_67, gd_69, gd_79, gd_81, hd_25, hd_45, \
                         hd_67, hd_79, hd_99, hd_111 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -5.625 * ab_x[k] * gd_25[k]
                  + 1.875 * ab_y[k] * gd_27[k]
                  - 5.625 * ab_x[k] * gd_67[k]
                  + 1.875 * ab_y[k] * gd_69[k]
                  + 7.5 * ab_x[k] * gd_79[k]
                  - 2.5 * ab_y[k] * gd_81[k]
                  - 5.625 * hd_25[k]
                  + 1.875 * hd_45[k]
                  - 5.625 * hd_67[k]
                  + 7.5 * hd_79[k]
                  + 1.875 * hd_99[k]
                  - 2.5 * hd_111[k];
    }

#pragma omp simd aligned(ab_x, gd_28, gd_70, gd_82, hd_28, hd_70, \
                         hd_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_39 * ab_x[k] * gd_28[k]
                  - f_39 * ab_x[k] * gd_70[k]
                  + f_40 * ab_x[k] * gd_82[k]
                  - f_39 * hd_28[k]
                  - f_39 * hd_70[k]
                  + f_40 * hd_82[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_25, gd_27, gd_29, gd_67, gd_69, gd_71, gd_79, gd_81, \
                         gd_83, hd_25, hd_45, hd_47, hd_67, hd_79, hd_99, hd_101, hd_111, \
                         hd_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_41 * ab_x[k] * gd_25[k]
                  + f_41 * ab_y[k] * gd_27[k]
                  - f_42 * ab_y[k] * gd_29[k]
                  + f_41 * ab_x[k] * gd_67[k]
                  + f_41 * ab_y[k] * gd_69[k]
                  - f_42 * ab_y[k] * gd_71[k]
                  - f_43 * ab_x[k] * gd_79[k]
                  - f_43 * ab_y[k] * gd_81[k]
                  + f_44 * ab_y[k] * gd_83[k]
                  + f_41 * hd_25[k]
                  + f_41 * hd_45[k]
                  - f_42 * hd_47[k]
                  + f_41 * hd_67[k]
                  - f_43 * hd_79[k]
                  + f_41 * hd_99[k]
                  - f_42 * hd_101[k]
                  - f_43 * hd_111[k]
                  + f_44 * hd_113[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gd_26, gd_28, gd_29, gd_68, gd_70, gd_71, gd_80, \
                         gd_82, gd_83, hd_26, hd_46, hd_53, hd_68, hd_80, hd_100, hd_107, \
                         hd_112, hd_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_45 * ab_x[k] * gd_26[k]
                  + f_45 * ab_y[k] * gd_28[k]
                  - f_46 * ab_z[k] * gd_29[k]
                  + f_45 * ab_x[k] * gd_68[k]
                  + f_45 * ab_y[k] * gd_70[k]
                  - f_46 * ab_z[k] * gd_71[k]
                  - f_47 * ab_x[k] * gd_80[k]
                  - f_47 * ab_y[k] * gd_82[k]
                  + f_48 * ab_z[k] * gd_83[k]
                  + f_45 * hd_26[k]
                  + f_45 * hd_46[k]
                  - f_46 * hd_53[k]
                  + f_45 * hd_68[k]
                  - f_47 * hd_80[k]
                  + f_45 * hd_100[k]
                  - f_46 * hd_107[k]
                  - f_47 * hd_112[k]
                  + f_48 * hd_119[k];
    }

#pragma omp simd aligned(ab_x, gd_24, gd_27, gd_29, gd_66, gd_69, gd_71, gd_78, gd_81, gd_83, \
                         hd_24, hd_27, hd_29, hd_66, hd_69, hd_71, hd_78, hd_81, \
                         hd_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_41 * ab_x[k] * gd_24[k]
                  + f_41 * ab_x[k] * gd_27[k]
                  - f_42 * ab_x[k] * gd_29[k]
                  + f_41 * ab_x[k] * gd_66[k]
                  + f_41 * ab_x[k] * gd_69[k]
                  - f_42 * ab_x[k] * gd_71[k]
                  - f_43 * ab_x[k] * gd_78[k]
                  - f_43 * ab_x[k] * gd_81[k]
                  + f_44 * ab_x[k] * gd_83[k]
                  + f_41 * hd_24[k]
                  + f_41 * hd_27[k]
                  - f_42 * hd_29[k]
                  + f_41 * hd_66[k]
                  + f_41 * hd_69[k]
                  - f_42 * hd_71[k]
                  - f_43 * hd_78[k]
                  - f_43 * hd_81[k]
                  + f_44 * hd_83[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_26, gd_28, gd_68, gd_70, gd_80, gd_82, hd_26, hd_46, \
                         hd_68, hd_80, hd_100, hd_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_49 * ab_x[k] * gd_26[k]
                  + f_49 * ab_y[k] * gd_28[k]
                  - f_49 * ab_x[k] * gd_68[k]
                  + f_49 * ab_y[k] * gd_70[k]
                  + f_50 * ab_x[k] * gd_80[k]
                  - f_50 * ab_y[k] * gd_82[k]
                  - f_49 * hd_26[k]
                  + f_49 * hd_46[k]
                  - f_49 * hd_68[k]
                  + f_50 * hd_80[k]
                  + f_49 * hd_100[k]
                  - f_50 * hd_112[k];
    }

#pragma omp simd aligned(ab_x, gd_24, gd_27, gd_66, gd_69, gd_78, gd_81, hd_24, hd_27, hd_66, \
                         hd_69, hd_78, hd_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -1.875 * ab_x[k] * gd_24[k]
                  + 5.625 * ab_x[k] * gd_27[k]
                  - 1.875 * ab_x[k] * gd_66[k]
                  + 5.625 * ab_x[k] * gd_69[k]
                  + 2.5 * ab_x[k] * gd_78[k]
                  - 7.5 * ab_x[k] * gd_81[k]
                  - 1.875 * hd_24[k]
                  + 5.625 * hd_27[k]
                  - 1.875 * hd_66[k]
                  + 5.625 * hd_69[k]
                  + 2.5 * hd_78[k]
                  - 7.5 * hd_81[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_1, gd_3, gd_19, gd_21, gd_31, gd_33, gd_61, gd_63, \
                         gd_73, gd_75, gd_85, gd_87, hd_1, hd_9, hd_19, hd_31, hd_39, hd_51, \
                         hd_61, hd_73, hd_85, hd_93, hd_105, hd_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_51 * ab_x[k] * gd_1[k]
                  - f_52 * ab_y[k] * gd_3[k]
                  + f_53 * ab_x[k] * gd_19[k]
                  - f_54 * ab_y[k] * gd_21[k]
                  - f_55 * ab_x[k] * gd_31[k]
                  + f_46 * ab_y[k] * gd_33[k]
                  + f_51 * ab_x[k] * gd_61[k]
                  - f_52 * ab_y[k] * gd_63[k]
                  - f_55 * ab_x[k] * gd_73[k]
                  + f_46 * ab_y[k] * gd_75[k]
                  + f_46 * ab_x[k] * gd_85[k]
                  - f_56 * ab_y[k] * gd_87[k]
                  + f_51 * hd_1[k]
                  - f_52 * hd_9[k]
                  + f_53 * hd_19[k]
                  - f_55 * hd_31[k]
                  - f_54 * hd_39[k]
                  + f_46 * hd_51[k]
                  + f_51 * hd_61[k]
                  - f_55 * hd_73[k]
                  + f_46 * hd_85[k]
                  - f_52 * hd_93[k]
                  + f_46 * hd_105[k]
                  - f_56 * hd_117[k];
    }

#pragma omp simd aligned(ab_x, gd_4, gd_22, gd_34, gd_64, gd_76, gd_88, hd_4, hd_22, hd_34, \
                         hd_64, hd_76, hd_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_41 * ab_x[k] * gd_4[k]
                  + f_57 * ab_x[k] * gd_22[k]
                  - f_58 * ab_x[k] * gd_34[k]
                  + f_41 * ab_x[k] * gd_64[k]
                  - f_58 * ab_x[k] * gd_76[k]
                  + f_59 * ab_x[k] * gd_88[k]
                  + f_41 * hd_4[k]
                  + f_57 * hd_22[k]
                  - f_58 * hd_34[k]
                  + f_41 * hd_64[k]
                  - f_58 * hd_76[k]
                  + f_59 * hd_88[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_1, gd_3, gd_5, gd_19, gd_21, gd_23, gd_31, gd_33, \
                         gd_35, gd_61, gd_63, gd_65, gd_73, gd_75, gd_77, gd_85, gd_87, gd_89, \
                         hd_1, hd_9, hd_11, hd_19, hd_31, hd_39, hd_41, hd_51, hd_53, hd_61, \
                         hd_73, hd_85, hd_93, hd_95, hd_105, hd_107, hd_117, \
                         hd_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_60 * ab_x[k] * gd_1[k]
                  - f_60 * ab_y[k] * gd_3[k]
                  + f_61 * ab_y[k] * gd_5[k]
                  - f_62 * ab_x[k] * gd_19[k]
                  - f_62 * ab_y[k] * gd_21[k]
                  + f_63 * ab_y[k] * gd_23[k]
                  + f_63 * ab_x[k] * gd_31[k]
                  + f_63 * ab_y[k] * gd_33[k]
                  - f_64 * ab_y[k] * gd_35[k]
                  - f_60 * ab_x[k] * gd_61[k]
                  - f_60 * ab_y[k] * gd_63[k]
                  + f_61 * ab_y[k] * gd_65[k]
                  + f_63 * ab_x[k] * gd_73[k]
                  + f_63 * ab_y[k] * gd_75[k]
                  - f_64 * ab_y[k] * gd_77[k]
                  - f_65 * ab_x[k] * gd_85[k]
                  - f_65 * ab_y[k] * gd_87[k]
                  + f_66 * ab_y[k] * gd_89[k]
                  - f_60 * hd_1[k]
                  - f_60 * hd_9[k]
                  + f_61 * hd_11[k]
                  - f_62 * hd_19[k]
                  + f_63 * hd_31[k]
                  - f_62 * hd_39[k]
                  + f_63 * hd_41[k]
                  + f_63 * hd_51[k]
                  - f_64 * hd_53[k]
                  - f_60 * hd_61[k]
                  + f_63 * hd_73[k]
                  - f_65 * hd_85[k]
                  - f_60 * hd_93[k]
                  + f_61 * hd_95[k]
                  + f_63 * hd_105[k]
                  - f_64 * hd_107[k]
                  - f_65 * hd_117[k]
                  + f_66 * hd_119[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gd_2, gd_4, gd_5, gd_20, gd_22, gd_23, gd_32, \
                         gd_34, gd_35, gd_62, gd_64, gd_65, gd_74, gd_76, gd_77, gd_86, gd_88, \
                         gd_89, hd_2, hd_10, hd_17, hd_20, hd_32, hd_40, hd_47, hd_52, hd_59, \
                         hd_62, hd_74, hd_86, hd_94, hd_101, hd_106, hd_113, hd_118, \
                         hd_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -0.5625 * ab_x[k] * gd_2[k]
                  - 0.5625 * ab_y[k] * gd_4[k]
                  + 0.375 * ab_z[k] * gd_5[k]
                  - 1.125 * ab_x[k] * gd_20[k]
                  - 1.125 * ab_y[k] * gd_22[k]
                  + 0.75 * ab_z[k] * gd_23[k]
                  + 4.5 * ab_x[k] * gd_32[k]
                  + 4.5 * ab_y[k] * gd_34[k]
                  - 3.0 * ab_z[k] * gd_35[k]
                  - 0.5625 * ab_x[k] * gd_62[k]
                  - 0.5625 * ab_y[k] * gd_64[k]
                  + 0.375 * ab_z[k] * gd_65[k]
                  + 4.5 * ab_x[k] * gd_74[k]
                  + 4.5 * ab_y[k] * gd_76[k]
                  - 3.0 * ab_z[k] * gd_77[k]
                  - 1.5 * ab_x[k] * gd_86[k]
                  - 1.5 * ab_y[k] * gd_88[k]
                  + ab_z[k] * gd_89[k]
                  - 0.5625 * hd_2[k]
                  - 0.5625 * hd_10[k]
                  + 0.375 * hd_17[k]
                  - 1.125 * hd_20[k]
                  + 4.5 * hd_32[k]
                  - 1.125 * hd_40[k]
                  + 0.75 * hd_47[k]
                  + 4.5 * hd_52[k]
                  - 3.0 * hd_59[k]
                  - 0.5625 * hd_62[k]
                  + 4.5 * hd_74[k]
                  - 1.5 * hd_86[k]
                  - 0.5625 * hd_94[k]
                  + 0.375 * hd_101[k]
                  + 4.5 * hd_106[k]
                  - 3.0 * hd_113[k]
                  - 1.5 * hd_118[k]
                  + hd_125[k];
    }

#pragma omp simd aligned(ab_x, gd_0, gd_3, gd_5, gd_18, gd_21, gd_23, gd_30, gd_33, gd_35, \
                         gd_60, gd_63, gd_65, gd_72, gd_75, gd_77, gd_84, gd_87, gd_89, hd_0, \
                         hd_3, hd_5, hd_18, hd_21, hd_23, hd_30, hd_33, hd_35, hd_60, hd_63, \
                         hd_65, hd_72, hd_75, hd_77, hd_84, hd_87, \
                         hd_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_60 * ab_x[k] * gd_0[k]
                  - f_60 * ab_x[k] * gd_3[k]
                  + f_61 * ab_x[k] * gd_5[k]
                  - f_62 * ab_x[k] * gd_18[k]
                  - f_62 * ab_x[k] * gd_21[k]
                  + f_63 * ab_x[k] * gd_23[k]
                  + f_63 * ab_x[k] * gd_30[k]
                  + f_63 * ab_x[k] * gd_33[k]
                  - f_64 * ab_x[k] * gd_35[k]
                  - f_60 * ab_x[k] * gd_60[k]
                  - f_60 * ab_x[k] * gd_63[k]
                  + f_61 * ab_x[k] * gd_65[k]
                  + f_63 * ab_x[k] * gd_72[k]
                  + f_63 * ab_x[k] * gd_75[k]
                  - f_64 * ab_x[k] * gd_77[k]
                  - f_65 * ab_x[k] * gd_84[k]
                  - f_65 * ab_x[k] * gd_87[k]
                  + f_66 * ab_x[k] * gd_89[k]
                  - f_60 * hd_0[k]
                  - f_60 * hd_3[k]
                  + f_61 * hd_5[k]
                  - f_62 * hd_18[k]
                  - f_62 * hd_21[k]
                  + f_63 * hd_23[k]
                  + f_63 * hd_30[k]
                  + f_63 * hd_33[k]
                  - f_64 * hd_35[k]
                  - f_60 * hd_60[k]
                  - f_60 * hd_63[k]
                  + f_61 * hd_65[k]
                  + f_63 * hd_72[k]
                  + f_63 * hd_75[k]
                  - f_64 * hd_77[k]
                  - f_65 * hd_84[k]
                  - f_65 * hd_87[k]
                  + f_66 * hd_89[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_2, gd_4, gd_20, gd_22, gd_32, gd_34, gd_62, gd_64, \
                         gd_74, gd_76, gd_86, gd_88, hd_2, hd_10, hd_20, hd_32, hd_40, hd_52, \
                         hd_62, hd_74, hd_86, hd_94, hd_106, hd_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_67 * ab_x[k] * gd_2[k]
                  - f_67 * ab_y[k] * gd_4[k]
                  + f_41 * ab_x[k] * gd_20[k]
                  - f_41 * ab_y[k] * gd_22[k]
                  - f_42 * ab_x[k] * gd_32[k]
                  + f_42 * ab_y[k] * gd_34[k]
                  + f_67 * ab_x[k] * gd_62[k]
                  - f_67 * ab_y[k] * gd_64[k]
                  - f_42 * ab_x[k] * gd_74[k]
                  + f_42 * ab_y[k] * gd_76[k]
                  + f_43 * ab_x[k] * gd_86[k]
                  - f_43 * ab_y[k] * gd_88[k]
                  + f_67 * hd_2[k]
                  - f_67 * hd_10[k]
                  + f_41 * hd_20[k]
                  - f_42 * hd_32[k]
                  - f_41 * hd_40[k]
                  + f_42 * hd_52[k]
                  + f_67 * hd_62[k]
                  - f_42 * hd_74[k]
                  + f_43 * hd_86[k]
                  - f_67 * hd_94[k]
                  + f_42 * hd_106[k]
                  - f_43 * hd_118[k];
    }

#pragma omp simd aligned(ab_x, gd_0, gd_3, gd_18, gd_21, gd_30, gd_33, gd_60, gd_63, gd_72, \
                         gd_75, gd_84, gd_87, hd_0, hd_3, hd_18, hd_21, hd_30, hd_33, hd_60, \
                         hd_63, hd_72, hd_75, hd_84, hd_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_52 * ab_x[k] * gd_0[k]
                  - f_51 * ab_x[k] * gd_3[k]
                  + f_54 * ab_x[k] * gd_18[k]
                  - f_53 * ab_x[k] * gd_21[k]
                  - f_46 * ab_x[k] * gd_30[k]
                  + f_55 * ab_x[k] * gd_33[k]
                  + f_52 * ab_x[k] * gd_60[k]
                  - f_51 * ab_x[k] * gd_63[k]
                  - f_46 * ab_x[k] * gd_72[k]
                  + f_55 * ab_x[k] * gd_75[k]
                  + f_56 * ab_x[k] * gd_84[k]
                  - f_46 * ab_x[k] * gd_87[k]
                  + f_52 * hd_0[k]
                  - f_51 * hd_3[k]
                  + f_54 * hd_18[k]
                  - f_53 * hd_21[k]
                  - f_46 * hd_30[k]
                  + f_55 * hd_33[k]
                  + f_52 * hd_60[k]
                  - f_51 * hd_63[k]
                  - f_46 * hd_72[k]
                  + f_55 * hd_75[k]
                  + f_56 * hd_84[k]
                  - f_46 * hd_87[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_13, gd_15, gd_43, gd_45, gd_55, gd_57, hd_13, hd_27, \
                         hd_43, hd_55, hd_69, hd_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -5.625 * ab_x[k] * gd_13[k]
                  + 1.875 * ab_y[k] * gd_15[k]
                  - 5.625 * ab_x[k] * gd_43[k]
                  + 1.875 * ab_y[k] * gd_45[k]
                  + 7.5 * ab_x[k] * gd_55[k]
                  - 2.5 * ab_y[k] * gd_57[k]
                  - 5.625 * hd_13[k]
                  + 1.875 * hd_27[k]
                  - 5.625 * hd_43[k]
                  + 7.5 * hd_55[k]
                  + 1.875 * hd_69[k]
                  - 2.5 * hd_81[k];
    }

#pragma omp simd aligned(ab_x, gd_16, gd_46, gd_58, hd_16, hd_46, \
                         hd_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = -f_39 * ab_x[k] * gd_16[k]
                  - f_39 * ab_x[k] * gd_46[k]
                  + f_40 * ab_x[k] * gd_58[k]
                  - f_39 * hd_16[k]
                  - f_39 * hd_46[k]
                  + f_40 * hd_58[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_13, gd_15, gd_17, gd_43, gd_45, gd_47, gd_55, gd_57, \
                         gd_59, hd_13, hd_27, hd_29, hd_43, hd_55, hd_69, hd_71, hd_81, \
                         hd_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_41 * ab_x[k] * gd_13[k]
                  + f_41 * ab_y[k] * gd_15[k]
                  - f_42 * ab_y[k] * gd_17[k]
                  + f_41 * ab_x[k] * gd_43[k]
                  + f_41 * ab_y[k] * gd_45[k]
                  - f_42 * ab_y[k] * gd_47[k]
                  - f_43 * ab_x[k] * gd_55[k]
                  - f_43 * ab_y[k] * gd_57[k]
                  + f_44 * ab_y[k] * gd_59[k]
                  + f_41 * hd_13[k]
                  + f_41 * hd_27[k]
                  - f_42 * hd_29[k]
                  + f_41 * hd_43[k]
                  - f_43 * hd_55[k]
                  + f_41 * hd_69[k]
                  - f_42 * hd_71[k]
                  - f_43 * hd_81[k]
                  + f_44 * hd_83[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gd_14, gd_16, gd_17, gd_44, gd_46, gd_47, gd_56, \
                         gd_58, gd_59, hd_14, hd_28, hd_35, hd_44, hd_56, hd_70, hd_77, hd_82, \
                         hd_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_45 * ab_x[k] * gd_14[k]
                  + f_45 * ab_y[k] * gd_16[k]
                  - f_46 * ab_z[k] * gd_17[k]
                  + f_45 * ab_x[k] * gd_44[k]
                  + f_45 * ab_y[k] * gd_46[k]
                  - f_46 * ab_z[k] * gd_47[k]
                  - f_47 * ab_x[k] * gd_56[k]
                  - f_47 * ab_y[k] * gd_58[k]
                  + f_48 * ab_z[k] * gd_59[k]
                  + f_45 * hd_14[k]
                  + f_45 * hd_28[k]
                  - f_46 * hd_35[k]
                  + f_45 * hd_44[k]
                  - f_47 * hd_56[k]
                  + f_45 * hd_70[k]
                  - f_46 * hd_77[k]
                  - f_47 * hd_82[k]
                  + f_48 * hd_89[k];
    }

#pragma omp simd aligned(ab_x, gd_12, gd_15, gd_17, gd_42, gd_45, gd_47, gd_54, gd_57, gd_59, \
                         hd_12, hd_15, hd_17, hd_42, hd_45, hd_47, hd_54, hd_57, \
                         hd_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_41 * ab_x[k] * gd_12[k]
                  + f_41 * ab_x[k] * gd_15[k]
                  - f_42 * ab_x[k] * gd_17[k]
                  + f_41 * ab_x[k] * gd_42[k]
                  + f_41 * ab_x[k] * gd_45[k]
                  - f_42 * ab_x[k] * gd_47[k]
                  - f_43 * ab_x[k] * gd_54[k]
                  - f_43 * ab_x[k] * gd_57[k]
                  + f_44 * ab_x[k] * gd_59[k]
                  + f_41 * hd_12[k]
                  + f_41 * hd_15[k]
                  - f_42 * hd_17[k]
                  + f_41 * hd_42[k]
                  + f_41 * hd_45[k]
                  - f_42 * hd_47[k]
                  - f_43 * hd_54[k]
                  - f_43 * hd_57[k]
                  + f_44 * hd_59[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_14, gd_16, gd_44, gd_46, gd_56, gd_58, hd_14, hd_28, \
                         hd_44, hd_56, hd_70, hd_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_49 * ab_x[k] * gd_14[k]
                  + f_49 * ab_y[k] * gd_16[k]
                  - f_49 * ab_x[k] * gd_44[k]
                  + f_49 * ab_y[k] * gd_46[k]
                  + f_50 * ab_x[k] * gd_56[k]
                  - f_50 * ab_y[k] * gd_58[k]
                  - f_49 * hd_14[k]
                  + f_49 * hd_28[k]
                  - f_49 * hd_44[k]
                  + f_50 * hd_56[k]
                  + f_49 * hd_70[k]
                  - f_50 * hd_82[k];
    }

#pragma omp simd aligned(ab_x, gd_12, gd_15, gd_42, gd_45, gd_54, gd_57, hd_12, hd_15, hd_42, \
                         hd_45, hd_54, hd_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -1.875 * ab_x[k] * gd_12[k]
                  + 5.625 * ab_x[k] * gd_15[k]
                  - 1.875 * ab_x[k] * gd_42[k]
                  + 5.625 * ab_x[k] * gd_45[k]
                  + 2.5 * ab_x[k] * gd_54[k]
                  - 7.5 * ab_x[k] * gd_57[k]
                  - 1.875 * hd_12[k]
                  + 5.625 * hd_15[k]
                  - 1.875 * hd_42[k]
                  + 5.625 * hd_45[k]
                  + 2.5 * hd_54[k]
                  - 7.5 * hd_57[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_1, gd_3, gd_31, gd_33, gd_61, gd_63, gd_73, gd_75, \
                         hd_1, hd_9, hd_31, hd_51, hd_61, hd_73, hd_93, \
                         hd_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_68 * ab_x[k] * gd_1[k]
                  + f_69 * ab_y[k] * gd_3[k]
                  + f_70 * ab_x[k] * gd_31[k]
                  - f_23 * ab_y[k] * gd_33[k]
                  + f_68 * ab_x[k] * gd_61[k]
                  - f_69 * ab_y[k] * gd_63[k]
                  - f_70 * ab_x[k] * gd_73[k]
                  + f_23 * ab_y[k] * gd_75[k]
                  - f_68 * hd_1[k]
                  + f_69 * hd_9[k]
                  + f_70 * hd_31[k]
                  - f_23 * hd_51[k]
                  + f_68 * hd_61[k]
                  - f_70 * hd_73[k]
                  - f_69 * hd_93[k]
                  + f_23 * hd_105[k];
    }

#pragma omp simd aligned(ab_x, gd_4, gd_34, gd_64, gd_76, hd_4, hd_34, hd_64, \
                         hd_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_37 * ab_x[k] * gd_4[k]
                  + f_38 * ab_x[k] * gd_34[k]
                  + f_37 * ab_x[k] * gd_64[k]
                  - f_38 * ab_x[k] * gd_76[k]
                  - f_37 * hd_4[k]
                  + f_38 * hd_34[k]
                  + f_37 * hd_64[k]
                  - f_38 * hd_76[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_1, gd_3, gd_5, gd_31, gd_33, gd_35, gd_61, gd_63, \
                         gd_65, gd_73, gd_75, gd_77, hd_1, hd_9, hd_11, hd_31, hd_51, hd_53, \
                         hd_61, hd_73, hd_93, hd_95, hd_105, hd_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_71 * ab_x[k] * gd_1[k]
                  + f_71 * ab_y[k] * gd_3[k]
                  - f_72 * ab_y[k] * gd_5[k]
                  - f_73 * ab_x[k] * gd_31[k]
                  - f_73 * ab_y[k] * gd_33[k]
                  + f_74 * ab_y[k] * gd_35[k]
                  - f_71 * ab_x[k] * gd_61[k]
                  - f_71 * ab_y[k] * gd_63[k]
                  + f_72 * ab_y[k] * gd_65[k]
                  + f_73 * ab_x[k] * gd_73[k]
                  + f_73 * ab_y[k] * gd_75[k]
                  - f_74 * ab_y[k] * gd_77[k]
                  + f_71 * hd_1[k]
                  + f_71 * hd_9[k]
                  - f_72 * hd_11[k]
                  - f_73 * hd_31[k]
                  - f_73 * hd_51[k]
                  + f_74 * hd_53[k]
                  - f_71 * hd_61[k]
                  + f_73 * hd_73[k]
                  - f_71 * hd_93[k]
                  + f_72 * hd_95[k]
                  + f_73 * hd_105[k]
                  - f_74 * hd_107[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gd_2, gd_4, gd_5, gd_32, gd_34, gd_35, gd_62, \
                         gd_64, gd_65, gd_74, gd_76, gd_77, hd_2, hd_10, hd_17, hd_32, hd_52, \
                         hd_59, hd_62, hd_74, hd_94, hd_101, hd_106, \
                         hd_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_75 * ab_x[k] * gd_2[k]
                  + f_75 * ab_y[k] * gd_4[k]
                  - f_76 * ab_z[k] * gd_5[k]
                  - f_77 * ab_x[k] * gd_32[k]
                  - f_77 * ab_y[k] * gd_34[k]
                  + f_78 * ab_z[k] * gd_35[k]
                  - f_75 * ab_x[k] * gd_62[k]
                  - f_75 * ab_y[k] * gd_64[k]
                  + f_76 * ab_z[k] * gd_65[k]
                  + f_77 * ab_x[k] * gd_74[k]
                  + f_77 * ab_y[k] * gd_76[k]
                  - f_78 * ab_z[k] * gd_77[k]
                  + f_75 * hd_2[k]
                  + f_75 * hd_10[k]
                  - f_76 * hd_17[k]
                  - f_77 * hd_32[k]
                  - f_77 * hd_52[k]
                  + f_78 * hd_59[k]
                  - f_75 * hd_62[k]
                  + f_77 * hd_74[k]
                  - f_75 * hd_94[k]
                  + f_76 * hd_101[k]
                  + f_77 * hd_106[k]
                  - f_78 * hd_113[k];
    }

#pragma omp simd aligned(ab_x, gd_0, gd_3, gd_5, gd_30, gd_33, gd_35, gd_60, gd_63, gd_65, \
                         gd_72, gd_75, gd_77, hd_0, hd_3, hd_5, hd_30, hd_33, hd_35, hd_60, \
                         hd_63, hd_65, hd_72, hd_75, hd_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_71 * ab_x[k] * gd_0[k]
                  + f_71 * ab_x[k] * gd_3[k]
                  - f_72 * ab_x[k] * gd_5[k]
                  - f_73 * ab_x[k] * gd_30[k]
                  - f_73 * ab_x[k] * gd_33[k]
                  + f_74 * ab_x[k] * gd_35[k]
                  - f_71 * ab_x[k] * gd_60[k]
                  - f_71 * ab_x[k] * gd_63[k]
                  + f_72 * ab_x[k] * gd_65[k]
                  + f_73 * ab_x[k] * gd_72[k]
                  + f_73 * ab_x[k] * gd_75[k]
                  - f_74 * ab_x[k] * gd_77[k]
                  + f_71 * hd_0[k]
                  + f_71 * hd_3[k]
                  - f_72 * hd_5[k]
                  - f_73 * hd_30[k]
                  - f_73 * hd_33[k]
                  + f_74 * hd_35[k]
                  - f_71 * hd_60[k]
                  - f_71 * hd_63[k]
                  + f_72 * hd_65[k]
                  + f_73 * hd_72[k]
                  + f_73 * hd_75[k]
                  - f_74 * hd_77[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_2, gd_4, gd_32, gd_34, gd_62, gd_64, gd_74, gd_76, \
                         hd_2, hd_10, hd_32, hd_52, hd_62, hd_74, hd_94, \
                         hd_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_79 * ab_x[k] * gd_2[k]
                  + f_79 * ab_y[k] * gd_4[k]
                  + f_80 * ab_x[k] * gd_32[k]
                  - f_80 * ab_y[k] * gd_34[k]
                  + f_79 * ab_x[k] * gd_62[k]
                  - f_79 * ab_y[k] * gd_64[k]
                  - f_80 * ab_x[k] * gd_74[k]
                  + f_80 * ab_y[k] * gd_76[k]
                  - f_79 * hd_2[k]
                  + f_79 * hd_10[k]
                  + f_80 * hd_32[k]
                  - f_80 * hd_52[k]
                  + f_79 * hd_62[k]
                  - f_80 * hd_74[k]
                  - f_79 * hd_94[k]
                  + f_80 * hd_106[k];
    }

#pragma omp simd aligned(ab_x, gd_0, gd_3, gd_30, gd_33, gd_60, gd_63, gd_72, gd_75, hd_0, \
                         hd_3, hd_30, hd_33, hd_60, hd_63, hd_72, \
                         hd_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_69 * ab_x[k] * gd_0[k]
                  + f_68 * ab_x[k] * gd_3[k]
                  + f_23 * ab_x[k] * gd_30[k]
                  - f_70 * ab_x[k] * gd_33[k]
                  + f_69 * ab_x[k] * gd_60[k]
                  - f_68 * ab_x[k] * gd_63[k]
                  - f_23 * ab_x[k] * gd_72[k]
                  + f_70 * ab_x[k] * gd_75[k]
                  - f_69 * hd_0[k]
                  + f_68 * hd_3[k]
                  + f_23 * hd_30[k]
                  - f_70 * hd_33[k]
                  + f_69 * hd_60[k]
                  - f_68 * hd_63[k]
                  - f_23 * hd_72[k]
                  + f_70 * hd_75[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_13, gd_15, gd_16, gd_43, gd_45, gd_46, hd_13, hd_16, \
                         hd_27, hd_43, hd_46, hd_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_9 * ab_x[k] * gd_13[k]
                  - f_10 * ab_y[k] * gd_15[k]
                  - f_8 * ab_x[k] * gd_43[k]
                  + f_9 * ab_y[k] * gd_45[k]
                  + f_9 * hd_13[k]
                  - f_10 * hd_27[k]
                  - f_8 * hd_43[k]
                  + f_9 * hd_69[k];

        g_50[k] = f_12 * ab_x[k] * gd_16[k]
                  - f_11 * ab_x[k] * gd_46[k]
                  + f_12 * hd_16[k]
                  - f_11 * hd_46[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_13, gd_15, gd_17, gd_43, gd_45, gd_47, hd_13, hd_27, \
                         hd_29, hd_43, hd_69, hd_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_15 * ab_x[k] * gd_13[k]
                  - f_15 * ab_y[k] * gd_15[k]
                  + f_16 * ab_y[k] * gd_17[k]
                  + f_13 * ab_x[k] * gd_43[k]
                  + f_13 * ab_y[k] * gd_45[k]
                  - f_14 * ab_y[k] * gd_47[k]
                  - f_15 * hd_13[k]
                  - f_15 * hd_27[k]
                  + f_16 * hd_29[k]
                  + f_13 * hd_43[k]
                  + f_13 * hd_69[k]
                  - f_14 * hd_71[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gd_14, gd_16, gd_17, gd_44, gd_46, gd_47, hd_14, \
                         hd_28, hd_35, hd_44, hd_70, hd_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_19 * ab_x[k] * gd_14[k]
                  - f_19 * ab_y[k] * gd_16[k]
                  + f_20 * ab_z[k] * gd_17[k]
                  + f_17 * ab_x[k] * gd_44[k]
                  + f_17 * ab_y[k] * gd_46[k]
                  - f_18 * ab_z[k] * gd_47[k]
                  - f_19 * hd_14[k]
                  - f_19 * hd_28[k]
                  + f_20 * hd_35[k]
                  + f_17 * hd_44[k]
                  + f_17 * hd_70[k]
                  - f_18 * hd_77[k];
    }

#pragma omp simd aligned(ab_x, gd_12, gd_15, gd_17, gd_42, gd_45, gd_47, hd_12, hd_15, hd_17, \
                         hd_42, hd_45, hd_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_15 * ab_x[k] * gd_12[k]
                  - f_15 * ab_x[k] * gd_15[k]
                  + f_16 * ab_x[k] * gd_17[k]
                  + f_13 * ab_x[k] * gd_42[k]
                  + f_13 * ab_x[k] * gd_45[k]
                  - f_14 * ab_x[k] * gd_47[k]
                  - f_15 * hd_12[k]
                  - f_15 * hd_15[k]
                  + f_16 * hd_17[k]
                  + f_13 * hd_42[k]
                  + f_13 * hd_45[k]
                  - f_14 * hd_47[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_14, gd_16, gd_44, gd_46, hd_14, hd_28, hd_44, \
                         hd_70 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_22 * ab_x[k] * gd_14[k]
                  - f_22 * ab_y[k] * gd_16[k]
                  - f_21 * ab_x[k] * gd_44[k]
                  + f_21 * ab_y[k] * gd_46[k]
                  + f_22 * hd_14[k]
                  - f_22 * hd_28[k]
                  - f_21 * hd_44[k]
                  + f_21 * hd_70[k];
    }

#pragma omp simd aligned(ab_x, gd_12, gd_15, gd_42, gd_45, hd_12, hd_15, hd_42, \
                         hd_45 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_10 * ab_x[k] * gd_12[k]
                  - f_9 * ab_x[k] * gd_15[k]
                  - f_9 * ab_x[k] * gd_42[k]
                  + f_8 * ab_x[k] * gd_45[k]
                  + f_10 * hd_12[k]
                  - f_9 * hd_15[k]
                  - f_9 * hd_42[k]
                  + f_8 * hd_45[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_1, gd_3, gd_19, gd_21, gd_61, gd_63, hd_1, hd_9, \
                         hd_19, hd_39, hd_61, hd_93 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_81 * ab_x[k] * gd_1[k]
                  - f_82 * ab_y[k] * gd_3[k]
                  - f_83 * ab_x[k] * gd_19[k]
                  + f_84 * ab_y[k] * gd_21[k]
                  + f_81 * ab_x[k] * gd_61[k]
                  - f_82 * ab_y[k] * gd_63[k]
                  + f_81 * hd_1[k]
                  - f_82 * hd_9[k]
                  - f_83 * hd_19[k]
                  + f_84 * hd_39[k]
                  + f_81 * hd_61[k]
                  - f_82 * hd_93[k];
    }

#pragma omp simd aligned(ab_x, gd_4, gd_22, gd_64, hd_4, hd_22, hd_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_85 * ab_x[k] * gd_4[k]
                  - f_86 * ab_x[k] * gd_22[k]
                  + f_85 * ab_x[k] * gd_64[k]
                  + f_85 * hd_4[k]
                  - f_86 * hd_22[k]
                  + f_85 * hd_64[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_1, gd_3, gd_5, gd_19, gd_21, gd_23, gd_61, gd_63, \
                         gd_65, hd_1, hd_9, hd_11, hd_19, hd_39, hd_41, hd_61, hd_93, \
                         hd_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_87 * ab_x[k] * gd_1[k]
                  - f_87 * ab_y[k] * gd_3[k]
                  + f_3 * ab_y[k] * gd_5[k]
                  + f_88 * ab_x[k] * gd_19[k]
                  + f_88 * ab_y[k] * gd_21[k]
                  - f_89 * ab_y[k] * gd_23[k]
                  - f_87 * ab_x[k] * gd_61[k]
                  - f_87 * ab_y[k] * gd_63[k]
                  + f_3 * ab_y[k] * gd_65[k]
                  - f_87 * hd_1[k]
                  - f_87 * hd_9[k]
                  + f_3 * hd_11[k]
                  + f_88 * hd_19[k]
                  + f_88 * hd_39[k]
                  - f_89 * hd_41[k]
                  - f_87 * hd_61[k]
                  - f_87 * hd_93[k]
                  + f_3 * hd_95[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gd_2, gd_4, gd_5, gd_20, gd_22, gd_23, gd_62, \
                         gd_64, gd_65, hd_2, hd_10, hd_17, hd_20, hd_40, hd_47, hd_62, hd_94, \
                         hd_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_90 * ab_x[k] * gd_2[k]
                  - f_90 * ab_y[k] * gd_4[k]
                  + f_91 * ab_z[k] * gd_5[k]
                  + f_92 * ab_x[k] * gd_20[k]
                  + f_92 * ab_y[k] * gd_22[k]
                  - f_5 * ab_z[k] * gd_23[k]
                  - f_90 * ab_x[k] * gd_62[k]
                  - f_90 * ab_y[k] * gd_64[k]
                  + f_91 * ab_z[k] * gd_65[k]
                  - f_90 * hd_2[k]
                  - f_90 * hd_10[k]
                  + f_91 * hd_17[k]
                  + f_92 * hd_20[k]
                  + f_92 * hd_40[k]
                  - f_5 * hd_47[k]
                  - f_90 * hd_62[k]
                  - f_90 * hd_94[k]
                  + f_91 * hd_101[k];
    }

#pragma omp simd aligned(ab_x, gd_0, gd_3, gd_5, gd_18, gd_21, gd_23, gd_60, gd_63, gd_65, \
                         hd_0, hd_3, hd_5, hd_18, hd_21, hd_23, hd_60, hd_63, \
                         hd_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_87 * ab_x[k] * gd_0[k]
                  - f_87 * ab_x[k] * gd_3[k]
                  + f_3 * ab_x[k] * gd_5[k]
                  + f_88 * ab_x[k] * gd_18[k]
                  + f_88 * ab_x[k] * gd_21[k]
                  - f_89 * ab_x[k] * gd_23[k]
                  - f_87 * ab_x[k] * gd_60[k]
                  - f_87 * ab_x[k] * gd_63[k]
                  + f_3 * ab_x[k] * gd_65[k]
                  - f_87 * hd_0[k]
                  - f_87 * hd_3[k]
                  + f_3 * hd_5[k]
                  + f_88 * hd_18[k]
                  + f_88 * hd_21[k]
                  - f_89 * hd_23[k]
                  - f_87 * hd_60[k]
                  - f_87 * hd_63[k]
                  + f_3 * hd_65[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gd_2, gd_4, gd_20, gd_22, gd_62, gd_64, hd_2, hd_10, \
                         hd_20, hd_40, hd_62, hd_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_93 * ab_x[k] * gd_2[k]
                  - f_93 * ab_y[k] * gd_4[k]
                  - f_94 * ab_x[k] * gd_20[k]
                  + f_94 * ab_y[k] * gd_22[k]
                  + f_93 * ab_x[k] * gd_62[k]
                  - f_93 * ab_y[k] * gd_64[k]
                  + f_93 * hd_2[k]
                  - f_93 * hd_10[k]
                  - f_94 * hd_20[k]
                  + f_94 * hd_40[k]
                  + f_93 * hd_62[k]
                  - f_93 * hd_94[k];
    }

#pragma omp simd aligned(ab_x, gd_0, gd_3, gd_18, gd_21, gd_60, gd_63, hd_0, hd_3, hd_18, \
                         hd_21, hd_60, hd_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_82 * ab_x[k] * gd_0[k]
                  - f_81 * ab_x[k] * gd_3[k]
                  - f_84 * ab_x[k] * gd_18[k]
                  + f_83 * ab_x[k] * gd_21[k]
                  + f_82 * ab_x[k] * gd_60[k]
                  - f_81 * ab_x[k] * gd_63[k]
                  + f_82 * hd_0[k]
                  - f_81 * hd_3[k]
                  - f_84 * hd_18[k]
                  + f_83 * hd_21[k]
                  + f_82 * hd_60[k]
                  - f_81 * hd_63[k];
    }
}

}  // namespace simdtrf
