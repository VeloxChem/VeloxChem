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


#include "SimdTransferFH.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hrr_fh_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t dh, const size_t di,
                   const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.40625 * std::sqrt(35.0);
    const auto f_1 = 2.8125 * std::sqrt(35.0);
    const auto f_2 = 0.28125 * std::sqrt(35.0);
    const auto f_3 = 0.46875 * std::sqrt(35.0);
    const auto f_4 = 0.9375 * std::sqrt(35.0);
    const auto f_5 = 0.09375 * std::sqrt(35.0);
    const auto f_6 = 5.625 * std::sqrt(14.0);
    const auto f_7 = 1.875 * std::sqrt(14.0);
    const auto f_8 = 1.40625 * std::sqrt(7.0);
    const auto f_9 = 0.9375 * std::sqrt(7.0);
    const auto f_10 = 11.25 * std::sqrt(7.0);
    const auto f_11 = 0.46875 * std::sqrt(7.0);
    const auto f_12 = 3.75 * std::sqrt(7.0);
    const auto f_13 = 0.3125 * std::sqrt(7.0);
    const auto f_14 = 0.15625 * std::sqrt(7.0);
    const auto f_15 = 1.25 * std::sqrt(7.0);
    const auto f_16 = 1.875 * std::sqrt(42.0);
    const auto f_17 = 3.75 * std::sqrt(42.0);
    const auto f_18 = 0.625 * std::sqrt(42.0);
    const auto f_19 = 1.25 * std::sqrt(42.0);
    const auto f_20 = 0.46875 * std::sqrt(6.0);
    const auto f_21 = 0.9375 * std::sqrt(6.0);
    const auto f_22 = 5.625 * std::sqrt(6.0);
    const auto f_23 = 3.75 * std::sqrt(6.0);
    const auto f_24 = 0.15625 * std::sqrt(6.0);
    const auto f_25 = 0.3125 * std::sqrt(6.0);
    const auto f_26 = 1.875 * std::sqrt(6.0);
    const auto f_27 = 1.25 * std::sqrt(6.0);
    const auto f_28 = 1.40625 * std::sqrt(10.0);
    const auto f_29 = 2.8125 * std::sqrt(10.0);
    const auto f_30 = 3.75 * std::sqrt(10.0);
    const auto f_31 = 0.75 * std::sqrt(10.0);
    const auto f_32 = 0.46875 * std::sqrt(10.0);
    const auto f_33 = 0.9375 * std::sqrt(10.0);
    const auto f_34 = 1.25 * std::sqrt(10.0);
    const auto f_35 = 0.25 * std::sqrt(10.0);
    const auto f_36 = 0.9375 * std::sqrt(42.0);
    const auto f_37 = 0.3125 * std::sqrt(42.0);
    const auto f_38 = 1.40625 * std::sqrt(14.0);
    const auto f_39 = 8.4375 * std::sqrt(14.0);
    const auto f_40 = 0.46875 * std::sqrt(14.0);
    const auto f_41 = 2.8125 * std::sqrt(14.0);
    const auto f_42 = 0.9375 * std::sqrt(210.0);
    const auto f_43 = 1.875 * std::sqrt(210.0);
    const auto f_44 = 0.1875 * std::sqrt(210.0);
    const auto f_45 = 7.5 * std::sqrt(21.0);
    const auto f_46 = 7.5 * std::sqrt(42.0);
    const auto f_47 = 2.5 * std::sqrt(42.0);
    const auto f_48 = 7.5 * std::sqrt(7.0);
    const auto f_49 = 15.0 * std::sqrt(7.0);
    const auto f_50 = 1.875 * std::sqrt(15.0);
    const auto f_51 = 3.75 * std::sqrt(15.0);
    const auto f_52 = 5.0 * std::sqrt(15.0);
    const auto f_53 = std::sqrt(15.0);
    const auto f_54 = 1.875 * std::sqrt(21.0);
    const auto f_55 = 11.25 * std::sqrt(21.0);
    const auto f_56 = 0.46875 * std::sqrt(21.0);
    const auto f_57 = 0.9375 * std::sqrt(21.0);
    const auto f_58 = 0.09375 * std::sqrt(21.0);
    const auto f_59 = 3.75 * std::sqrt(21.0);
    const auto f_60 = 0.375 * std::sqrt(21.0);
    const auto f_61 = 0.375 * std::sqrt(210.0);
    const auto f_62 = 1.5 * std::sqrt(210.0);
    const auto f_63 = 0.09375 * std::sqrt(105.0);
    const auto f_64 = 0.0625 * std::sqrt(105.0);
    const auto f_65 = 0.75 * std::sqrt(105.0);
    const auto f_66 = 0.03125 * std::sqrt(105.0);
    const auto f_67 = 0.25 * std::sqrt(105.0);
    const auto f_68 = 0.375 * std::sqrt(105.0);
    const auto f_69 = 3.0 * std::sqrt(105.0);
    const auto f_70 = 0.125 * std::sqrt(105.0);
    const auto f_71 = std::sqrt(105.0);
    const auto f_72 = 0.375 * std::sqrt(70.0);
    const auto f_73 = 0.75 * std::sqrt(70.0);
    const auto f_74 = 1.5 * std::sqrt(70.0);
    const auto f_75 = 3.0 * std::sqrt(70.0);
    const auto f_76 = 0.09375 * std::sqrt(10.0);
    const auto f_77 = 0.1875 * std::sqrt(10.0);
    const auto f_78 = 1.125 * std::sqrt(10.0);
    const auto f_79 = 0.375 * std::sqrt(10.0);
    const auto f_80 = 4.5 * std::sqrt(10.0);
    const auto f_81 = 3.0 * std::sqrt(10.0);
    const auto f_82 = 0.25 * std::sqrt(6.0);
    const auto f_83 = 5.0 * std::sqrt(6.0);
    const auto f_84 = std::sqrt(6.0);
    const auto f_85 = 0.1875 * std::sqrt(70.0);
    const auto f_86 = 0.09375 * std::sqrt(210.0);
    const auto f_87 = 0.5625 * std::sqrt(210.0);
    const auto f_88 = 2.25 * std::sqrt(210.0);
    const auto f_89 = 0.28125 * std::sqrt(14.0);
    const auto f_90 = 0.9375 * std::sqrt(14.0);
    const auto f_91 = 0.1875 * std::sqrt(14.0);
    const auto f_92 = 2.25 * std::sqrt(35.0);
    const auto f_93 = 1.5 * std::sqrt(35.0);
    const auto f_94 = 0.28125 * std::sqrt(70.0);
    const auto f_95 = 2.25 * std::sqrt(70.0);
    const auto f_96 = 0.09375 * std::sqrt(70.0);
    const auto f_97 = 0.125 * std::sqrt(70.0);
    const auto f_98 = 0.0625 * std::sqrt(70.0);
    const auto f_99 = 0.5 * std::sqrt(70.0);
    const auto f_100 = 1.5 * std::sqrt(105.0);
    const auto f_101 = 0.5 * std::sqrt(105.0);
    const auto f_102 = 0.1875 * std::sqrt(15.0);
    const auto f_103 = 0.375 * std::sqrt(15.0);
    const auto f_104 = 2.25 * std::sqrt(15.0);
    const auto f_105 = 1.5 * std::sqrt(15.0);
    const auto f_106 = 0.125 * std::sqrt(15.0);
    const auto f_107 = 0.25 * std::sqrt(15.0);
    const auto f_108 = 0.5625 * std::sqrt(35.0);
    const auto f_109 = 3.375 * std::sqrt(35.0);
    const auto f_110 = 0.375 * std::sqrt(35.0);
    const auto f_111 = 0.46875 * std::sqrt(210.0);
    const auto f_112 = 0.46875 * std::sqrt(42.0);
    const auto f_113 = 0.15625 * std::sqrt(42.0);
    const auto f_114 = 0.9375 * std::sqrt(15.0);
    const auto f_115 = 2.5 * std::sqrt(15.0);
    const auto f_116 = 0.5 * std::sqrt(15.0);
    const auto f_117 = 1.875 * std::sqrt(7.0);
    const auto f_118 = 5.625 * std::sqrt(21.0);

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
    auto *g_63 = values + 63 * nvalues;
    auto *g_64 = values + 64 * nvalues;
    auto *g_65 = values + 65 * nvalues;
    auto *g_66 = values + 66 * nvalues;
    auto *g_67 = values + 67 * nvalues;
    auto *g_68 = values + 68 * nvalues;
    auto *g_69 = values + 69 * nvalues;
    auto *g_70 = values + 70 * nvalues;
    auto *g_71 = values + 71 * nvalues;
    auto *g_72 = values + 72 * nvalues;
    auto *g_73 = values + 73 * nvalues;
    auto *g_74 = values + 74 * nvalues;
    auto *g_75 = values + 75 * nvalues;
    auto *g_76 = values + 76 * nvalues;

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
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
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_41 = buffer.data(dh + 41);
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
    const auto *dh_57 = buffer.data(dh + 57);
    const auto *dh_58 = buffer.data(dh + 58);
    const auto *dh_59 = buffer.data(dh + 59);
    const auto *dh_60 = buffer.data(dh + 60);
    const auto *dh_61 = buffer.data(dh + 61);
    const auto *dh_62 = buffer.data(dh + 62);
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
    const auto *dh_83 = buffer.data(dh + 83);
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
    const auto *dh_104 = buffer.data(dh + 104);
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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_57 = buffer.data(di + 57);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_63 = buffer.data(di + 63);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_66 = buffer.data(di + 66);
    const auto *di_67 = buffer.data(di + 67);
    const auto *di_68 = buffer.data(di + 68);
    const auto *di_69 = buffer.data(di + 69);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_71 = buffer.data(di + 71);
    const auto *di_72 = buffer.data(di + 72);
    const auto *di_73 = buffer.data(di + 73);
    const auto *di_74 = buffer.data(di + 74);
    const auto *di_75 = buffer.data(di + 75);
    const auto *di_76 = buffer.data(di + 76);
    const auto *di_84 = buffer.data(di + 84);
    const auto *di_85 = buffer.data(di + 85);
    const auto *di_86 = buffer.data(di + 86);
    const auto *di_87 = buffer.data(di + 87);
    const auto *di_88 = buffer.data(di + 88);
    const auto *di_89 = buffer.data(di + 89);
    const auto *di_90 = buffer.data(di + 90);
    const auto *di_91 = buffer.data(di + 91);
    const auto *di_92 = buffer.data(di + 92);
    const auto *di_93 = buffer.data(di + 93);
    const auto *di_94 = buffer.data(di + 94);
    const auto *di_95 = buffer.data(di + 95);
    const auto *di_96 = buffer.data(di + 96);
    const auto *di_97 = buffer.data(di + 97);
    const auto *di_98 = buffer.data(di + 98);
    const auto *di_99 = buffer.data(di + 99);
    const auto *di_100 = buffer.data(di + 100);
    const auto *di_101 = buffer.data(di + 101);
    const auto *di_102 = buffer.data(di + 102);
    const auto *di_103 = buffer.data(di + 103);
    const auto *di_104 = buffer.data(di + 104);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_106 = buffer.data(di + 106);
    const auto *di_107 = buffer.data(di + 107);
    const auto *di_108 = buffer.data(di + 108);
    const auto *di_109 = buffer.data(di + 109);
    const auto *di_110 = buffer.data(di + 110);
    const auto *di_112 = buffer.data(di + 112);
    const auto *di_113 = buffer.data(di + 113);
    const auto *di_114 = buffer.data(di + 114);
    const auto *di_115 = buffer.data(di + 115);
    const auto *di_116 = buffer.data(di + 116);
    const auto *di_117 = buffer.data(di + 117);
    const auto *di_118 = buffer.data(di + 118);
    const auto *di_119 = buffer.data(di + 119);
    const auto *di_120 = buffer.data(di + 120);
    const auto *di_121 = buffer.data(di + 121);
    const auto *di_122 = buffer.data(di + 122);
    const auto *di_123 = buffer.data(di + 123);
    const auto *di_124 = buffer.data(di + 124);
    const auto *di_125 = buffer.data(di + 125);
    const auto *di_126 = buffer.data(di + 126);
    const auto *di_127 = buffer.data(di + 127);
    const auto *di_128 = buffer.data(di + 128);
    const auto *di_129 = buffer.data(di + 129);
    const auto *di_130 = buffer.data(di + 130);
    const auto *di_131 = buffer.data(di + 131);
    const auto *di_132 = buffer.data(di + 132);
    const auto *di_133 = buffer.data(di + 133);
    const auto *di_134 = buffer.data(di + 134);
    const auto *di_135 = buffer.data(di + 135);
    const auto *di_136 = buffer.data(di + 136);
    const auto *di_137 = buffer.data(di + 137);
    const auto *di_138 = buffer.data(di + 138);
    const auto *di_140 = buffer.data(di + 140);
    const auto *di_141 = buffer.data(di + 141);
    const auto *di_142 = buffer.data(di + 142);
    const auto *di_143 = buffer.data(di + 143);
    const auto *di_144 = buffer.data(di + 144);
    const auto *di_145 = buffer.data(di + 145);
    const auto *di_146 = buffer.data(di + 146);
    const auto *di_147 = buffer.data(di + 147);
    const auto *di_148 = buffer.data(di + 148);
    const auto *di_149 = buffer.data(di + 149);
    const auto *di_150 = buffer.data(di + 150);
    const auto *di_151 = buffer.data(di + 151);
    const auto *di_152 = buffer.data(di + 152);
    const auto *di_153 = buffer.data(di + 153);
    const auto *di_154 = buffer.data(di + 154);
    const auto *di_155 = buffer.data(di + 155);
    const auto *di_156 = buffer.data(di + 156);
    const auto *di_157 = buffer.data(di + 157);
    const auto *di_158 = buffer.data(di + 158);
    const auto *di_159 = buffer.data(di + 159);
    const auto *di_160 = buffer.data(di + 160);
    const auto *di_161 = buffer.data(di + 161);
    const auto *di_162 = buffer.data(di + 162);
    const auto *di_163 = buffer.data(di + 163);
    const auto *di_164 = buffer.data(di + 164);
    const auto *di_165 = buffer.data(di + 165);
    const auto *di_166 = buffer.data(di + 166);
    const auto *di_167 = buffer.data(di + 167);

#pragma omp simd aligned(ab_x, ab_y, dh_22, dh_27, dh_36, dh_64, dh_69, dh_78, di_29, di_34, \
                         di_43, di_87, di_94, di_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -f_0 * ab_x[k] * dh_22[k]
                 + f_1 * ab_x[k] * dh_27[k]
                 - f_2 * ab_x[k] * dh_36[k]
                 + f_3 * ab_y[k] * dh_64[k]
                 - f_4 * ab_y[k] * dh_69[k]
                 + f_5 * ab_y[k] * dh_78[k]
                 + f_0 * di_29[k]
                 - f_1 * di_34[k]
                 + f_2 * di_43[k]
                 - f_3 * di_87[k]
                 + f_4 * di_94[k]
                 - f_5 * di_105[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_25, dh_32, dh_67, dh_74, di_32, di_39, di_91, \
                         di_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = -f_6 * ab_x[k] * dh_25[k]
                 + f_6 * ab_x[k] * dh_32[k]
                 + f_7 * ab_y[k] * dh_67[k]
                 - f_7 * ab_y[k] * dh_74[k]
                 + f_6 * di_32[k]
                 - f_6 * di_39[k]
                 - f_7 * di_91[k]
                 + f_7 * di_100[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_22, dh_27, dh_29, dh_36, dh_38, dh_64, dh_69, dh_71, \
                         dh_78, dh_80, di_29, di_34, di_36, di_43, di_45, di_87, di_94, di_96, \
                         di_105, di_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = f_8 * ab_x[k] * dh_22[k]
                 + f_9 * ab_x[k] * dh_27[k]
                 - f_10 * ab_x[k] * dh_29[k]
                 - f_11 * ab_x[k] * dh_36[k]
                 + f_12 * ab_x[k] * dh_38[k]
                 - f_11 * ab_y[k] * dh_64[k]
                 - f_13 * ab_y[k] * dh_69[k]
                 + f_12 * ab_y[k] * dh_71[k]
                 + f_14 * ab_y[k] * dh_78[k]
                 - f_15 * ab_y[k] * dh_80[k]
                 - f_8 * di_29[k]
                 - f_9 * di_34[k]
                 + f_10 * di_36[k]
                 + f_11 * di_43[k]
                 - f_12 * di_45[k]
                 + f_11 * di_87[k]
                 + f_13 * di_94[k]
                 - f_12 * di_96[k]
                 - f_14 * di_105[k]
                 + f_15 * di_107[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_25, dh_32, dh_34, dh_67, dh_74, dh_76, di_32, di_39, \
                         di_41, di_91, di_100, di_102 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_16 * ab_x[k] * dh_25[k]
                 + f_16 * ab_x[k] * dh_32[k]
                 - f_17 * ab_x[k] * dh_34[k]
                 - f_18 * ab_y[k] * dh_67[k]
                 - f_18 * ab_y[k] * dh_74[k]
                 + f_19 * ab_y[k] * dh_76[k]
                 - f_16 * di_32[k]
                 - f_16 * di_39[k]
                 + f_17 * di_41[k]
                 + f_18 * di_91[k]
                 + f_18 * di_100[k]
                 - f_19 * di_102[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_22, dh_27, dh_29, dh_36, dh_38, dh_40, dh_64, dh_69, \
                         dh_71, dh_78, dh_80, dh_82, di_29, di_34, di_36, di_43, di_45, di_47, \
                         di_87, di_94, di_96, di_105, di_107, di_109 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_20 * ab_x[k] * dh_22[k]
                 - f_21 * ab_x[k] * dh_27[k]
                 + f_22 * ab_x[k] * dh_29[k]
                 - f_20 * ab_x[k] * dh_36[k]
                 + f_22 * ab_x[k] * dh_38[k]
                 - f_23 * ab_x[k] * dh_40[k]
                 + f_24 * ab_y[k] * dh_64[k]
                 + f_25 * ab_y[k] * dh_69[k]
                 - f_26 * ab_y[k] * dh_71[k]
                 + f_24 * ab_y[k] * dh_78[k]
                 - f_26 * ab_y[k] * dh_80[k]
                 + f_27 * ab_y[k] * dh_82[k]
                 + f_20 * di_29[k]
                 + f_21 * di_34[k]
                 - f_22 * di_36[k]
                 + f_20 * di_43[k]
                 - f_22 * di_45[k]
                 + f_23 * di_47[k]
                 - f_24 * di_87[k]
                 - f_25 * di_94[k]
                 + f_26 * di_96[k]
                 - f_24 * di_105[k]
                 + f_26 * di_107[k]
                 - f_27 * di_109[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_23, dh_28, dh_30, dh_37, dh_39, dh_41, dh_65, dh_70, \
                         dh_72, dh_79, dh_81, dh_83, di_30, di_35, di_37, di_44, di_46, di_48, \
                         di_88, di_95, di_97, di_106, di_108, di_110 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_28 * ab_x[k] * dh_23[k]
                 - f_29 * ab_x[k] * dh_28[k]
                 + f_30 * ab_x[k] * dh_30[k]
                 - f_28 * ab_x[k] * dh_37[k]
                 + f_30 * ab_x[k] * dh_39[k]
                 - f_31 * ab_x[k] * dh_41[k]
                 + f_32 * ab_y[k] * dh_65[k]
                 + f_33 * ab_y[k] * dh_70[k]
                 - f_34 * ab_y[k] * dh_72[k]
                 + f_32 * ab_y[k] * dh_79[k]
                 - f_34 * ab_y[k] * dh_81[k]
                 + f_35 * ab_y[k] * dh_83[k]
                 + f_28 * di_30[k]
                 + f_29 * di_35[k]
                 - f_30 * di_37[k]
                 + f_28 * di_44[k]
                 - f_30 * di_46[k]
                 + f_31 * di_48[k]
                 - f_32 * di_88[k]
                 - f_33 * di_95[k]
                 + f_34 * di_97[k]
                 - f_32 * di_106[k]
                 + f_34 * di_108[k]
                 - f_35 * di_110[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_21, dh_24, dh_26, dh_31, dh_33, dh_35, dh_63, dh_66, \
                         dh_68, dh_73, dh_75, dh_77, di_28, di_31, di_33, di_38, di_40, di_42, \
                         di_85, di_90, di_92, di_99, di_101, di_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_20 * ab_x[k] * dh_21[k]
                 - f_21 * ab_x[k] * dh_24[k]
                 + f_22 * ab_x[k] * dh_26[k]
                 - f_20 * ab_x[k] * dh_31[k]
                 + f_22 * ab_x[k] * dh_33[k]
                 - f_23 * ab_x[k] * dh_35[k]
                 + f_24 * ab_y[k] * dh_63[k]
                 + f_25 * ab_y[k] * dh_66[k]
                 - f_26 * ab_y[k] * dh_68[k]
                 + f_24 * ab_y[k] * dh_73[k]
                 - f_26 * ab_y[k] * dh_75[k]
                 + f_27 * ab_y[k] * dh_77[k]
                 + f_20 * di_28[k]
                 + f_21 * di_31[k]
                 - f_22 * di_33[k]
                 + f_20 * di_38[k]
                 - f_22 * di_40[k]
                 + f_23 * di_42[k]
                 - f_24 * di_85[k]
                 - f_25 * di_90[k]
                 + f_26 * di_92[k]
                 - f_24 * di_99[k]
                 + f_26 * di_101[k]
                 - f_27 * di_103[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_23, dh_30, dh_37, dh_39, dh_65, dh_72, dh_79, dh_81, \
                         di_30, di_37, di_44, di_46, di_88, di_97, di_106, \
                         di_108 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_36 * ab_x[k] * dh_23[k]
                 - f_16 * ab_x[k] * dh_30[k]
                 - f_36 * ab_x[k] * dh_37[k]
                 + f_16 * ab_x[k] * dh_39[k]
                 - f_37 * ab_y[k] * dh_65[k]
                 + f_18 * ab_y[k] * dh_72[k]
                 + f_37 * ab_y[k] * dh_79[k]
                 - f_18 * ab_y[k] * dh_81[k]
                 - f_36 * di_30[k]
                 + f_16 * di_37[k]
                 + f_36 * di_44[k]
                 - f_16 * di_46[k]
                 + f_37 * di_88[k]
                 - f_18 * di_97[k]
                 - f_37 * di_106[k]
                 + f_18 * di_108[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_21, dh_24, dh_26, dh_31, dh_33, dh_63, dh_66, dh_68, \
                         dh_73, dh_75, di_28, di_31, di_33, di_38, di_40, di_85, di_90, di_92, \
                         di_99, di_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_11 * ab_x[k] * dh_21[k]
                 - f_9 * ab_x[k] * dh_24[k]
                 - f_12 * ab_x[k] * dh_26[k]
                 - f_8 * ab_x[k] * dh_31[k]
                 + f_10 * ab_x[k] * dh_33[k]
                 - f_14 * ab_y[k] * dh_63[k]
                 + f_13 * ab_y[k] * dh_66[k]
                 + f_15 * ab_y[k] * dh_68[k]
                 + f_11 * ab_y[k] * dh_73[k]
                 - f_12 * ab_y[k] * dh_75[k]
                 - f_11 * di_28[k]
                 + f_9 * di_31[k]
                 + f_12 * di_33[k]
                 + f_8 * di_38[k]
                 - f_10 * di_40[k]
                 + f_14 * di_85[k]
                 - f_13 * di_90[k]
                 - f_15 * di_92[k]
                 - f_11 * di_99[k]
                 + f_12 * di_101[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_23, dh_28, dh_37, dh_65, dh_70, dh_79, di_30, di_35, \
                         di_44, di_88, di_95, di_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_38 * ab_x[k] * dh_23[k]
                 + f_39 * ab_x[k] * dh_28[k]
                 - f_38 * ab_x[k] * dh_37[k]
                 + f_40 * ab_y[k] * dh_65[k]
                 - f_41 * ab_y[k] * dh_70[k]
                 + f_40 * ab_y[k] * dh_79[k]
                 + f_38 * di_30[k]
                 - f_39 * di_35[k]
                 + f_38 * di_44[k]
                 - f_40 * di_88[k]
                 + f_41 * di_95[k]
                 - f_40 * di_106[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_21, dh_24, dh_31, dh_63, dh_66, dh_73, di_28, di_31, \
                         di_38, di_85, di_90, di_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_2 * ab_x[k] * dh_21[k]
                  + f_1 * ab_x[k] * dh_24[k]
                  - f_0 * ab_x[k] * dh_31[k]
                  + f_5 * ab_y[k] * dh_63[k]
                  - f_4 * ab_y[k] * dh_66[k]
                  + f_3 * ab_y[k] * dh_73[k]
                  + f_2 * di_28[k]
                  - f_1 * di_31[k]
                  + f_0 * di_38[k]
                  - f_5 * di_85[k]
                  + f_4 * di_90[k]
                  - f_3 * di_99[k];
    }

#pragma omp simd aligned(ab_x, dh_85, dh_88, dh_90, dh_95, dh_99, di_113, di_116, di_118, \
                         di_123, di_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_42 * ab_x[k] * dh_85[k]
                  + f_43 * ab_x[k] * dh_90[k]
                  - f_44 * ab_x[k] * dh_99[k]
                  + f_42 * di_113[k]
                  - f_43 * di_118[k]
                  + f_44 * di_127[k];

        g_12[k] = -f_45 * ab_x[k] * dh_88[k]
                  + f_45 * ab_x[k] * dh_95[k]
                  + f_45 * di_116[k]
                  - f_45 * di_123[k];
    }

#pragma omp simd aligned(ab_x, dh_85, dh_90, dh_92, dh_99, dh_101, di_113, di_118, di_120, \
                         di_127, di_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_36 * ab_x[k] * dh_85[k]
                  + f_18 * ab_x[k] * dh_90[k]
                  - f_46 * ab_x[k] * dh_92[k]
                  - f_37 * ab_x[k] * dh_99[k]
                  + f_47 * ab_x[k] * dh_101[k]
                  - f_36 * di_113[k]
                  - f_18 * di_118[k]
                  + f_46 * di_120[k]
                  + f_37 * di_127[k]
                  - f_47 * di_129[k];
    }

#pragma omp simd aligned(ab_x, dh_88, dh_95, dh_97, di_116, di_123, \
                         di_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_48 * ab_x[k] * dh_88[k]
                  + f_48 * ab_x[k] * dh_95[k]
                  - f_49 * ab_x[k] * dh_97[k]
                  - f_48 * di_116[k]
                  - f_48 * di_123[k]
                  + f_49 * di_125[k];
    }

#pragma omp simd aligned(ab_x, dh_85, dh_90, dh_92, dh_99, dh_101, dh_103, di_113, di_118, \
                         di_120, di_127, di_129, di_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -1.875 * ab_x[k] * dh_85[k]
                  - 3.75 * ab_x[k] * dh_90[k]
                  + 22.5 * ab_x[k] * dh_92[k]
                  - 1.875 * ab_x[k] * dh_99[k]
                  + 22.5 * ab_x[k] * dh_101[k]
                  - 15.0 * ab_x[k] * dh_103[k]
                  + 1.875 * di_113[k]
                  + 3.75 * di_118[k]
                  - 22.5 * di_120[k]
                  + 1.875 * di_127[k]
                  - 22.5 * di_129[k]
                  + 15.0 * di_131[k];
    }

#pragma omp simd aligned(ab_x, dh_86, dh_91, dh_93, dh_100, dh_102, dh_104, di_114, di_119, \
                         di_121, di_128, di_130, di_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -f_50 * ab_x[k] * dh_86[k]
                  - f_51 * ab_x[k] * dh_91[k]
                  + f_52 * ab_x[k] * dh_93[k]
                  - f_50 * ab_x[k] * dh_100[k]
                  + f_52 * ab_x[k] * dh_102[k]
                  - f_53 * ab_x[k] * dh_104[k]
                  + f_50 * di_114[k]
                  + f_51 * di_119[k]
                  - f_52 * di_121[k]
                  + f_50 * di_128[k]
                  - f_52 * di_130[k]
                  + f_53 * di_132[k];
    }

#pragma omp simd aligned(ab_x, dh_84, dh_87, dh_89, dh_94, dh_96, dh_98, di_112, di_115, \
                         di_117, di_122, di_124, di_126 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -1.875 * ab_x[k] * dh_84[k]
                  - 3.75 * ab_x[k] * dh_87[k]
                  + 22.5 * ab_x[k] * dh_89[k]
                  - 1.875 * ab_x[k] * dh_94[k]
                  + 22.5 * ab_x[k] * dh_96[k]
                  - 15.0 * ab_x[k] * dh_98[k]
                  + 1.875 * di_112[k]
                  + 3.75 * di_115[k]
                  - 22.5 * di_117[k]
                  + 1.875 * di_122[k]
                  - 22.5 * di_124[k]
                  + 15.0 * di_126[k];
    }

#pragma omp simd aligned(ab_x, dh_86, dh_93, dh_100, dh_102, di_114, di_121, di_128, \
                         di_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_12 * ab_x[k] * dh_86[k]
                  - f_48 * ab_x[k] * dh_93[k]
                  - f_12 * ab_x[k] * dh_100[k]
                  + f_48 * ab_x[k] * dh_102[k]
                  - f_12 * di_114[k]
                  + f_48 * di_121[k]
                  + f_12 * di_128[k]
                  - f_48 * di_130[k];
    }

#pragma omp simd aligned(ab_x, dh_84, dh_87, dh_89, dh_94, dh_96, di_112, di_115, di_117, \
                         di_122, di_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_37 * ab_x[k] * dh_84[k]
                  - f_18 * ab_x[k] * dh_87[k]
                  - f_47 * ab_x[k] * dh_89[k]
                  - f_36 * ab_x[k] * dh_94[k]
                  + f_46 * ab_x[k] * dh_96[k]
                  - f_37 * di_112[k]
                  + f_18 * di_115[k]
                  + f_47 * di_117[k]
                  + f_36 * di_122[k]
                  - f_46 * di_124[k];
    }

#pragma omp simd aligned(ab_x, dh_84, dh_86, dh_87, dh_91, dh_94, dh_100, di_112, di_114, \
                         di_115, di_119, di_122, di_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_54 * ab_x[k] * dh_86[k]
                  + f_55 * ab_x[k] * dh_91[k]
                  - f_54 * ab_x[k] * dh_100[k]
                  + f_54 * di_114[k]
                  - f_55 * di_119[k]
                  + f_54 * di_128[k];

        g_21[k] = -f_44 * ab_x[k] * dh_84[k]
                  + f_43 * ab_x[k] * dh_87[k]
                  - f_42 * ab_x[k] * dh_94[k]
                  + f_44 * di_112[k]
                  - f_43 * di_115[k]
                  + f_42 * di_122[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_22, dh_27, dh_36, dh_64, dh_69, dh_78, dh_106, dh_111, \
                         dh_120, di_29, di_34, di_43, di_87, di_94, di_105, di_143, di_150, \
                         di_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_56 * ab_x[k] * dh_22[k]
                  - f_57 * ab_x[k] * dh_27[k]
                  + f_58 * ab_x[k] * dh_36[k]
                  + f_56 * ab_y[k] * dh_64[k]
                  - f_57 * ab_y[k] * dh_69[k]
                  + f_58 * ab_y[k] * dh_78[k]
                  - f_54 * ab_y[k] * dh_106[k]
                  + f_59 * ab_y[k] * dh_111[k]
                  - f_60 * ab_y[k] * dh_120[k]
                  - f_56 * di_29[k]
                  + f_57 * di_34[k]
                  - f_58 * di_43[k]
                  - f_56 * di_87[k]
                  + f_57 * di_94[k]
                  - f_58 * di_105[k]
                  + f_54 * di_143[k]
                  - f_59 * di_150[k]
                  + f_60 * di_161[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_25, dh_32, dh_67, dh_74, dh_109, dh_116, di_32, di_39, \
                         di_91, di_100, di_147, di_156 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_61 * ab_x[k] * dh_25[k]
                  - f_61 * ab_x[k] * dh_32[k]
                  + f_61 * ab_y[k] * dh_67[k]
                  - f_61 * ab_y[k] * dh_74[k]
                  - f_62 * ab_y[k] * dh_109[k]
                  + f_62 * ab_y[k] * dh_116[k]
                  - f_61 * di_32[k]
                  + f_61 * di_39[k]
                  - f_61 * di_91[k]
                  + f_61 * di_100[k]
                  + f_62 * di_147[k]
                  - f_62 * di_156[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_22, dh_27, dh_29, dh_36, dh_38, dh_64, dh_69, dh_71, \
                         dh_78, dh_80, dh_106, dh_111, dh_113, dh_120, dh_122, di_29, di_34, \
                         di_36, di_43, di_45, di_87, di_94, di_96, di_105, di_107, di_143, \
                         di_150, di_152, di_161, di_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -f_63 * ab_x[k] * dh_22[k]
                  - f_64 * ab_x[k] * dh_27[k]
                  + f_65 * ab_x[k] * dh_29[k]
                  + f_66 * ab_x[k] * dh_36[k]
                  - f_67 * ab_x[k] * dh_38[k]
                  - f_63 * ab_y[k] * dh_64[k]
                  - f_64 * ab_y[k] * dh_69[k]
                  + f_65 * ab_y[k] * dh_71[k]
                  + f_66 * ab_y[k] * dh_78[k]
                  - f_67 * ab_y[k] * dh_80[k]
                  + f_68 * ab_y[k] * dh_106[k]
                  + f_67 * ab_y[k] * dh_111[k]
                  - f_69 * ab_y[k] * dh_113[k]
                  - f_70 * ab_y[k] * dh_120[k]
                  + f_71 * ab_y[k] * dh_122[k]
                  + f_63 * di_29[k]
                  + f_64 * di_34[k]
                  - f_65 * di_36[k]
                  - f_66 * di_43[k]
                  + f_67 * di_45[k]
                  + f_63 * di_87[k]
                  + f_64 * di_94[k]
                  - f_65 * di_96[k]
                  - f_66 * di_105[k]
                  + f_67 * di_107[k]
                  - f_68 * di_143[k]
                  - f_67 * di_150[k]
                  + f_69 * di_152[k]
                  + f_70 * di_161[k]
                  - f_71 * di_163[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_25, dh_32, dh_34, dh_67, dh_74, dh_76, dh_109, dh_116, \
                         dh_118, di_32, di_39, di_41, di_91, di_100, di_102, di_147, di_156, \
                         di_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_72 * ab_x[k] * dh_25[k]
                  - f_72 * ab_x[k] * dh_32[k]
                  + f_73 * ab_x[k] * dh_34[k]
                  - f_72 * ab_y[k] * dh_67[k]
                  - f_72 * ab_y[k] * dh_74[k]
                  + f_73 * ab_y[k] * dh_76[k]
                  + f_74 * ab_y[k] * dh_109[k]
                  + f_74 * ab_y[k] * dh_116[k]
                  - f_75 * ab_y[k] * dh_118[k]
                  + f_72 * di_32[k]
                  + f_72 * di_39[k]
                  - f_73 * di_41[k]
                  + f_72 * di_91[k]
                  + f_72 * di_100[k]
                  - f_73 * di_102[k]
                  - f_74 * di_147[k]
                  - f_74 * di_156[k]
                  + f_75 * di_158[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_22, dh_27, dh_29, dh_36, dh_38, dh_40, dh_64, dh_69, \
                         dh_71, dh_78, dh_80, dh_82, dh_106, dh_111, dh_113, dh_120, dh_122, \
                         dh_124, di_29, di_34, di_36, di_43, di_45, di_47, di_87, di_94, \
                         di_96, di_105, di_107, di_109, di_143, di_150, di_152, di_161, \
                         di_163, di_165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_76 * ab_x[k] * dh_22[k]
                  + f_77 * ab_x[k] * dh_27[k]
                  - f_78 * ab_x[k] * dh_29[k]
                  + f_76 * ab_x[k] * dh_36[k]
                  - f_78 * ab_x[k] * dh_38[k]
                  + f_31 * ab_x[k] * dh_40[k]
                  + f_76 * ab_y[k] * dh_64[k]
                  + f_77 * ab_y[k] * dh_69[k]
                  - f_78 * ab_y[k] * dh_71[k]
                  + f_76 * ab_y[k] * dh_78[k]
                  - f_78 * ab_y[k] * dh_80[k]
                  + f_31 * ab_y[k] * dh_82[k]
                  - f_79 * ab_y[k] * dh_106[k]
                  - f_31 * ab_y[k] * dh_111[k]
                  + f_80 * ab_y[k] * dh_113[k]
                  - f_79 * ab_y[k] * dh_120[k]
                  + f_80 * ab_y[k] * dh_122[k]
                  - f_81 * ab_y[k] * dh_124[k]
                  - f_76 * di_29[k]
                  - f_77 * di_34[k]
                  + f_78 * di_36[k]
                  - f_76 * di_43[k]
                  + f_78 * di_45[k]
                  - f_31 * di_47[k]
                  - f_76 * di_87[k]
                  - f_77 * di_94[k]
                  + f_78 * di_96[k]
                  - f_76 * di_105[k]
                  + f_78 * di_107[k]
                  - f_31 * di_109[k]
                  + f_79 * di_143[k]
                  + f_31 * di_150[k]
                  - f_80 * di_152[k]
                  + f_79 * di_161[k]
                  - f_80 * di_163[k]
                  + f_81 * di_165[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_23, dh_28, dh_30, dh_37, dh_39, dh_41, dh_65, dh_70, \
                         dh_72, dh_79, dh_81, dh_83, dh_107, dh_112, dh_114, dh_121, dh_123, \
                         dh_125, di_30, di_35, di_37, di_44, di_46, di_48, di_88, di_95, \
                         di_97, di_106, di_108, di_110, di_144, di_151, di_153, di_162, \
                         di_164, di_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = f_20 * ab_x[k] * dh_23[k]
                  + f_21 * ab_x[k] * dh_28[k]
                  - f_27 * ab_x[k] * dh_30[k]
                  + f_20 * ab_x[k] * dh_37[k]
                  - f_27 * ab_x[k] * dh_39[k]
                  + f_82 * ab_x[k] * dh_41[k]
                  + f_20 * ab_y[k] * dh_65[k]
                  + f_21 * ab_y[k] * dh_70[k]
                  - f_27 * ab_y[k] * dh_72[k]
                  + f_20 * ab_y[k] * dh_79[k]
                  - f_27 * ab_y[k] * dh_81[k]
                  + f_82 * ab_y[k] * dh_83[k]
                  - f_26 * ab_y[k] * dh_107[k]
                  - f_23 * ab_y[k] * dh_112[k]
                  + f_83 * ab_y[k] * dh_114[k]
                  - f_26 * ab_y[k] * dh_121[k]
                  + f_83 * ab_y[k] * dh_123[k]
                  - f_84 * ab_y[k] * dh_125[k]
                  - f_20 * di_30[k]
                  - f_21 * di_35[k]
                  + f_27 * di_37[k]
                  - f_20 * di_44[k]
                  + f_27 * di_46[k]
                  - f_82 * di_48[k]
                  - f_20 * di_88[k]
                  - f_21 * di_95[k]
                  + f_27 * di_97[k]
                  - f_20 * di_106[k]
                  + f_27 * di_108[k]
                  - f_82 * di_110[k]
                  + f_26 * di_144[k]
                  + f_23 * di_151[k]
                  - f_83 * di_153[k]
                  + f_26 * di_162[k]
                  - f_83 * di_164[k]
                  + f_84 * di_166[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_21, dh_24, dh_26, dh_31, dh_33, dh_35, dh_63, dh_66, \
                         dh_68, dh_73, dh_75, dh_77, dh_105, dh_108, dh_110, dh_115, dh_117, \
                         dh_119, di_28, di_31, di_33, di_38, di_40, di_42, di_85, di_90, \
                         di_92, di_99, di_101, di_103, di_141, di_146, di_148, di_155, di_157, \
                         di_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_76 * ab_x[k] * dh_21[k]
                  + f_77 * ab_x[k] * dh_24[k]
                  - f_78 * ab_x[k] * dh_26[k]
                  + f_76 * ab_x[k] * dh_31[k]
                  - f_78 * ab_x[k] * dh_33[k]
                  + f_31 * ab_x[k] * dh_35[k]
                  + f_76 * ab_y[k] * dh_63[k]
                  + f_77 * ab_y[k] * dh_66[k]
                  - f_78 * ab_y[k] * dh_68[k]
                  + f_76 * ab_y[k] * dh_73[k]
                  - f_78 * ab_y[k] * dh_75[k]
                  + f_31 * ab_y[k] * dh_77[k]
                  - f_79 * ab_y[k] * dh_105[k]
                  - f_31 * ab_y[k] * dh_108[k]
                  + f_80 * ab_y[k] * dh_110[k]
                  - f_79 * ab_y[k] * dh_115[k]
                  + f_80 * ab_y[k] * dh_117[k]
                  - f_81 * ab_y[k] * dh_119[k]
                  - f_76 * di_28[k]
                  - f_77 * di_31[k]
                  + f_78 * di_33[k]
                  - f_76 * di_38[k]
                  + f_78 * di_40[k]
                  - f_31 * di_42[k]
                  - f_76 * di_85[k]
                  - f_77 * di_90[k]
                  + f_78 * di_92[k]
                  - f_76 * di_99[k]
                  + f_78 * di_101[k]
                  - f_31 * di_103[k]
                  + f_79 * di_141[k]
                  + f_31 * di_146[k]
                  - f_80 * di_148[k]
                  + f_79 * di_155[k]
                  - f_80 * di_157[k]
                  + f_81 * di_159[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_23, dh_30, dh_37, dh_39, dh_65, dh_72, dh_79, dh_81, \
                         dh_107, dh_114, dh_121, dh_123, di_30, di_37, di_44, di_46, di_88, \
                         di_97, di_106, di_108, di_144, di_153, di_162, \
                         di_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = -f_85 * ab_x[k] * dh_23[k]
                  + f_72 * ab_x[k] * dh_30[k]
                  + f_85 * ab_x[k] * dh_37[k]
                  - f_72 * ab_x[k] * dh_39[k]
                  - f_85 * ab_y[k] * dh_65[k]
                  + f_72 * ab_y[k] * dh_72[k]
                  + f_85 * ab_y[k] * dh_79[k]
                  - f_72 * ab_y[k] * dh_81[k]
                  + f_73 * ab_y[k] * dh_107[k]
                  - f_74 * ab_y[k] * dh_114[k]
                  - f_73 * ab_y[k] * dh_121[k]
                  + f_74 * ab_y[k] * dh_123[k]
                  + f_85 * di_30[k]
                  - f_72 * di_37[k]
                  - f_85 * di_44[k]
                  + f_72 * di_46[k]
                  + f_85 * di_88[k]
                  - f_72 * di_97[k]
                  - f_85 * di_106[k]
                  + f_72 * di_108[k]
                  - f_73 * di_144[k]
                  + f_74 * di_153[k]
                  + f_73 * di_162[k]
                  - f_74 * di_164[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_21, dh_24, dh_26, dh_31, dh_33, dh_63, dh_66, dh_68, \
                         dh_73, dh_75, dh_105, dh_108, dh_110, dh_115, dh_117, di_28, di_31, \
                         di_33, di_38, di_40, di_85, di_90, di_92, di_99, di_101, di_141, \
                         di_146, di_148, di_155, di_157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_66 * ab_x[k] * dh_21[k]
                  + f_64 * ab_x[k] * dh_24[k]
                  + f_67 * ab_x[k] * dh_26[k]
                  + f_63 * ab_x[k] * dh_31[k]
                  - f_65 * ab_x[k] * dh_33[k]
                  - f_66 * ab_y[k] * dh_63[k]
                  + f_64 * ab_y[k] * dh_66[k]
                  + f_67 * ab_y[k] * dh_68[k]
                  + f_63 * ab_y[k] * dh_73[k]
                  - f_65 * ab_y[k] * dh_75[k]
                  + f_70 * ab_y[k] * dh_105[k]
                  - f_67 * ab_y[k] * dh_108[k]
                  - f_71 * ab_y[k] * dh_110[k]
                  - f_68 * ab_y[k] * dh_115[k]
                  + f_69 * ab_y[k] * dh_117[k]
                  + f_66 * di_28[k]
                  - f_64 * di_31[k]
                  - f_67 * di_33[k]
                  - f_63 * di_38[k]
                  + f_65 * di_40[k]
                  + f_66 * di_85[k]
                  - f_64 * di_90[k]
                  - f_67 * di_92[k]
                  - f_63 * di_99[k]
                  + f_65 * di_101[k]
                  - f_70 * di_141[k]
                  + f_67 * di_146[k]
                  + f_71 * di_148[k]
                  + f_68 * di_155[k]
                  - f_69 * di_157[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_23, dh_28, dh_37, dh_65, dh_70, dh_79, dh_107, dh_112, \
                         dh_121, di_30, di_35, di_44, di_88, di_95, di_106, di_144, di_151, \
                         di_162 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = f_86 * ab_x[k] * dh_23[k]
                  - f_87 * ab_x[k] * dh_28[k]
                  + f_86 * ab_x[k] * dh_37[k]
                  + f_86 * ab_y[k] * dh_65[k]
                  - f_87 * ab_y[k] * dh_70[k]
                  + f_86 * ab_y[k] * dh_79[k]
                  - f_61 * ab_y[k] * dh_107[k]
                  + f_88 * ab_y[k] * dh_112[k]
                  - f_61 * ab_y[k] * dh_121[k]
                  - f_86 * di_30[k]
                  + f_87 * di_35[k]
                  - f_86 * di_44[k]
                  - f_86 * di_88[k]
                  + f_87 * di_95[k]
                  - f_86 * di_106[k]
                  + f_61 * di_144[k]
                  - f_88 * di_151[k]
                  + f_61 * di_162[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_21, dh_24, dh_31, dh_63, dh_66, dh_73, dh_105, dh_108, \
                         dh_115, di_28, di_31, di_38, di_85, di_90, di_99, di_141, di_146, \
                         di_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_58 * ab_x[k] * dh_21[k]
                  - f_57 * ab_x[k] * dh_24[k]
                  + f_56 * ab_x[k] * dh_31[k]
                  + f_58 * ab_y[k] * dh_63[k]
                  - f_57 * ab_y[k] * dh_66[k]
                  + f_56 * ab_y[k] * dh_73[k]
                  - f_60 * ab_y[k] * dh_105[k]
                  + f_59 * ab_y[k] * dh_108[k]
                  - f_54 * ab_y[k] * dh_115[k]
                  - f_58 * di_28[k]
                  + f_57 * di_31[k]
                  - f_56 * di_38[k]
                  - f_58 * di_85[k]
                  + f_57 * di_90[k]
                  - f_56 * di_99[k]
                  + f_60 * di_141[k]
                  - f_59 * di_146[k]
                  + f_54 * di_155[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dh_43, dh_48, dh_57, dh_85, dh_90, dh_99, dh_106, \
                         dh_111, dh_120, di_57, di_62, di_71, di_115, di_122, di_133, di_144, \
                         di_151, di_162 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_38 * ab_x[k] * dh_43[k]
                  - f_41 * ab_x[k] * dh_48[k]
                  + f_89 * ab_x[k] * dh_57[k]
                  + f_38 * ab_y[k] * dh_85[k]
                  - f_41 * ab_y[k] * dh_90[k]
                  + f_89 * ab_y[k] * dh_99[k]
                  - f_90 * ab_z[k] * dh_106[k]
                  + f_7 * ab_z[k] * dh_111[k]
                  - f_91 * ab_z[k] * dh_120[k]
                  - f_38 * di_57[k]
                  + f_41 * di_62[k]
                  - f_89 * di_71[k]
                  - f_38 * di_115[k]
                  + f_41 * di_122[k]
                  - f_89 * di_133[k]
                  + f_90 * di_144[k]
                  - f_7 * di_151[k]
                  + f_91 * di_162[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dh_46, dh_53, dh_88, dh_95, dh_109, dh_116, di_60, \
                         di_67, di_119, di_128, di_148, di_157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_92 * ab_x[k] * dh_46[k]
                  - f_92 * ab_x[k] * dh_53[k]
                  + f_92 * ab_y[k] * dh_88[k]
                  - f_92 * ab_y[k] * dh_95[k]
                  - f_93 * ab_z[k] * dh_109[k]
                  + f_93 * ab_z[k] * dh_116[k]
                  - f_92 * di_60[k]
                  + f_92 * di_67[k]
                  - f_92 * di_119[k]
                  + f_92 * di_128[k]
                  + f_93 * di_148[k]
                  - f_93 * di_157[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dh_43, dh_48, dh_50, dh_57, dh_59, dh_85, dh_90, \
                         dh_92, dh_99, dh_101, dh_106, dh_111, dh_113, dh_120, dh_122, di_57, \
                         di_62, di_64, di_71, di_73, di_115, di_122, di_124, di_133, di_135, \
                         di_144, di_151, di_153, di_162, di_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_94 * ab_x[k] * dh_43[k]
                  - f_85 * ab_x[k] * dh_48[k]
                  + f_95 * ab_x[k] * dh_50[k]
                  + f_96 * ab_x[k] * dh_57[k]
                  - f_73 * ab_x[k] * dh_59[k]
                  - f_94 * ab_y[k] * dh_85[k]
                  - f_85 * ab_y[k] * dh_90[k]
                  + f_95 * ab_y[k] * dh_92[k]
                  + f_96 * ab_y[k] * dh_99[k]
                  - f_73 * ab_y[k] * dh_101[k]
                  + f_85 * ab_z[k] * dh_106[k]
                  + f_97 * ab_z[k] * dh_111[k]
                  - f_74 * ab_z[k] * dh_113[k]
                  - f_98 * ab_z[k] * dh_120[k]
                  + f_99 * ab_z[k] * dh_122[k]
                  + f_94 * di_57[k]
                  + f_85 * di_62[k]
                  - f_95 * di_64[k]
                  - f_96 * di_71[k]
                  + f_73 * di_73[k]
                  + f_94 * di_115[k]
                  + f_85 * di_122[k]
                  - f_95 * di_124[k]
                  - f_96 * di_133[k]
                  + f_73 * di_135[k]
                  - f_85 * di_144[k]
                  - f_97 * di_151[k]
                  + f_74 * di_153[k]
                  + f_98 * di_162[k]
                  - f_99 * di_164[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dh_46, dh_53, dh_55, dh_88, dh_95, dh_97, dh_109, \
                         dh_116, dh_118, di_60, di_67, di_69, di_119, di_128, di_130, di_148, \
                         di_157, di_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = -f_65 * ab_x[k] * dh_46[k]
                  - f_65 * ab_x[k] * dh_53[k]
                  + f_100 * ab_x[k] * dh_55[k]
                  - f_65 * ab_y[k] * dh_88[k]
                  - f_65 * ab_y[k] * dh_95[k]
                  + f_100 * ab_y[k] * dh_97[k]
                  + f_101 * ab_z[k] * dh_109[k]
                  + f_101 * ab_z[k] * dh_116[k]
                  - f_71 * ab_z[k] * dh_118[k]
                  + f_65 * di_60[k]
                  + f_65 * di_67[k]
                  - f_100 * di_69[k]
                  + f_65 * di_119[k]
                  + f_65 * di_128[k]
                  - f_100 * di_130[k]
                  - f_101 * di_148[k]
                  - f_101 * di_157[k]
                  + f_71 * di_159[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dh_43, dh_48, dh_50, dh_57, dh_59, dh_61, dh_85, \
                         dh_90, dh_92, dh_99, dh_101, dh_103, dh_106, dh_111, dh_113, dh_120, \
                         dh_122, dh_124, di_57, di_62, di_64, di_71, di_73, di_75, di_115, \
                         di_122, di_124, di_133, di_135, di_137, di_144, di_151, di_153, \
                         di_162, di_164, di_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_102 * ab_x[k] * dh_43[k]
                  + f_103 * ab_x[k] * dh_48[k]
                  - f_104 * ab_x[k] * dh_50[k]
                  + f_102 * ab_x[k] * dh_57[k]
                  - f_104 * ab_x[k] * dh_59[k]
                  + f_105 * ab_x[k] * dh_61[k]
                  + f_102 * ab_y[k] * dh_85[k]
                  + f_103 * ab_y[k] * dh_90[k]
                  - f_104 * ab_y[k] * dh_92[k]
                  + f_102 * ab_y[k] * dh_99[k]
                  - f_104 * ab_y[k] * dh_101[k]
                  + f_105 * ab_y[k] * dh_103[k]
                  - f_106 * ab_z[k] * dh_106[k]
                  - f_107 * ab_z[k] * dh_111[k]
                  + f_105 * ab_z[k] * dh_113[k]
                  - f_106 * ab_z[k] * dh_120[k]
                  + f_105 * ab_z[k] * dh_122[k]
                  - f_53 * ab_z[k] * dh_124[k]
                  - f_102 * di_57[k]
                  - f_103 * di_62[k]
                  + f_104 * di_64[k]
                  - f_102 * di_71[k]
                  + f_104 * di_73[k]
                  - f_105 * di_75[k]
                  - f_102 * di_115[k]
                  - f_103 * di_122[k]
                  + f_104 * di_124[k]
                  - f_102 * di_133[k]
                  + f_104 * di_135[k]
                  - f_105 * di_137[k]
                  + f_106 * di_144[k]
                  + f_107 * di_151[k]
                  - f_105 * di_153[k]
                  + f_106 * di_162[k]
                  - f_105 * di_164[k]
                  + f_53 * di_166[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dh_44, dh_49, dh_51, dh_58, dh_60, dh_62, dh_86, \
                         dh_91, dh_93, dh_100, dh_102, dh_104, dh_107, dh_112, dh_114, dh_121, \
                         dh_123, dh_125, di_58, di_63, di_65, di_72, di_74, di_76, di_116, \
                         di_123, di_125, di_134, di_136, di_138, di_145, di_152, di_154, \
                         di_163, di_165, di_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = 2.8125 * ab_x[k] * dh_44[k]
                  + 5.625 * ab_x[k] * dh_49[k]
                  - 7.5 * ab_x[k] * dh_51[k]
                  + 2.8125 * ab_x[k] * dh_58[k]
                  - 7.5 * ab_x[k] * dh_60[k]
                  + 1.5 * ab_x[k] * dh_62[k]
                  + 2.8125 * ab_y[k] * dh_86[k]
                  + 5.625 * ab_y[k] * dh_91[k]
                  - 7.5 * ab_y[k] * dh_93[k]
                  + 2.8125 * ab_y[k] * dh_100[k]
                  - 7.5 * ab_y[k] * dh_102[k]
                  + 1.5 * ab_y[k] * dh_104[k]
                  - 1.875 * ab_z[k] * dh_107[k]
                  - 3.75 * ab_z[k] * dh_112[k]
                  + 5.0 * ab_z[k] * dh_114[k]
                  - 1.875 * ab_z[k] * dh_121[k]
                  + 5.0 * ab_z[k] * dh_123[k]
                  - ab_z[k] * dh_125[k]
                  - 2.8125 * di_58[k]
                  - 5.625 * di_63[k]
                  + 7.5 * di_65[k]
                  - 2.8125 * di_72[k]
                  + 7.5 * di_74[k]
                  - 1.5 * di_76[k]
                  - 2.8125 * di_116[k]
                  - 5.625 * di_123[k]
                  + 7.5 * di_125[k]
                  - 2.8125 * di_134[k]
                  + 7.5 * di_136[k]
                  - 1.5 * di_138[k]
                  + 1.875 * di_145[k]
                  + 3.75 * di_152[k]
                  - 5.0 * di_154[k]
                  + 1.875 * di_163[k]
                  - 5.0 * di_165[k]
                  + di_167[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dh_42, dh_45, dh_47, dh_52, dh_54, dh_56, dh_84, \
                         dh_87, dh_89, dh_94, dh_96, dh_98, dh_105, dh_108, dh_110, dh_115, \
                         dh_117, dh_119, di_56, di_59, di_61, di_66, di_68, di_70, di_113, \
                         di_118, di_120, di_127, di_129, di_131, di_142, di_147, di_149, \
                         di_156, di_158, di_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_102 * ab_x[k] * dh_42[k]
                  + f_103 * ab_x[k] * dh_45[k]
                  - f_104 * ab_x[k] * dh_47[k]
                  + f_102 * ab_x[k] * dh_52[k]
                  - f_104 * ab_x[k] * dh_54[k]
                  + f_105 * ab_x[k] * dh_56[k]
                  + f_102 * ab_y[k] * dh_84[k]
                  + f_103 * ab_y[k] * dh_87[k]
                  - f_104 * ab_y[k] * dh_89[k]
                  + f_102 * ab_y[k] * dh_94[k]
                  - f_104 * ab_y[k] * dh_96[k]
                  + f_105 * ab_y[k] * dh_98[k]
                  - f_106 * ab_z[k] * dh_105[k]
                  - f_107 * ab_z[k] * dh_108[k]
                  + f_105 * ab_z[k] * dh_110[k]
                  - f_106 * ab_z[k] * dh_115[k]
                  + f_105 * ab_z[k] * dh_117[k]
                  - f_53 * ab_z[k] * dh_119[k]
                  - f_102 * di_56[k]
                  - f_103 * di_59[k]
                  + f_104 * di_61[k]
                  - f_102 * di_66[k]
                  + f_104 * di_68[k]
                  - f_105 * di_70[k]
                  - f_102 * di_113[k]
                  - f_103 * di_118[k]
                  + f_104 * di_120[k]
                  - f_102 * di_127[k]
                  + f_104 * di_129[k]
                  - f_105 * di_131[k]
                  + f_106 * di_142[k]
                  + f_107 * di_147[k]
                  - f_105 * di_149[k]
                  + f_106 * di_156[k]
                  - f_105 * di_158[k]
                  + f_53 * di_160[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dh_44, dh_51, dh_58, dh_60, dh_86, dh_93, dh_100, \
                         dh_102, dh_107, dh_114, dh_121, dh_123, di_58, di_65, di_72, di_74, \
                         di_116, di_125, di_134, di_136, di_145, di_154, di_163, \
                         di_165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_68 * ab_x[k] * dh_44[k]
                  + f_65 * ab_x[k] * dh_51[k]
                  + f_68 * ab_x[k] * dh_58[k]
                  - f_65 * ab_x[k] * dh_60[k]
                  - f_68 * ab_y[k] * dh_86[k]
                  + f_65 * ab_y[k] * dh_93[k]
                  + f_68 * ab_y[k] * dh_100[k]
                  - f_65 * ab_y[k] * dh_102[k]
                  + f_67 * ab_z[k] * dh_107[k]
                  - f_101 * ab_z[k] * dh_114[k]
                  - f_67 * ab_z[k] * dh_121[k]
                  + f_101 * ab_z[k] * dh_123[k]
                  + f_68 * di_58[k]
                  - f_65 * di_65[k]
                  - f_68 * di_72[k]
                  + f_65 * di_74[k]
                  + f_68 * di_116[k]
                  - f_65 * di_125[k]
                  - f_68 * di_134[k]
                  + f_65 * di_136[k]
                  - f_67 * di_145[k]
                  + f_101 * di_154[k]
                  + f_67 * di_163[k]
                  - f_101 * di_165[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dh_42, dh_45, dh_47, dh_52, dh_54, dh_84, dh_87, \
                         dh_89, dh_94, dh_96, dh_105, dh_108, dh_110, dh_115, dh_117, di_56, \
                         di_59, di_61, di_66, di_68, di_113, di_118, di_120, di_127, di_129, \
                         di_142, di_147, di_149, di_156, di_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_96 * ab_x[k] * dh_42[k]
                  + f_85 * ab_x[k] * dh_45[k]
                  + f_73 * ab_x[k] * dh_47[k]
                  + f_94 * ab_x[k] * dh_52[k]
                  - f_95 * ab_x[k] * dh_54[k]
                  - f_96 * ab_y[k] * dh_84[k]
                  + f_85 * ab_y[k] * dh_87[k]
                  + f_73 * ab_y[k] * dh_89[k]
                  + f_94 * ab_y[k] * dh_94[k]
                  - f_95 * ab_y[k] * dh_96[k]
                  + f_98 * ab_z[k] * dh_105[k]
                  - f_97 * ab_z[k] * dh_108[k]
                  - f_99 * ab_z[k] * dh_110[k]
                  - f_85 * ab_z[k] * dh_115[k]
                  + f_74 * ab_z[k] * dh_117[k]
                  + f_96 * di_56[k]
                  - f_85 * di_59[k]
                  - f_73 * di_61[k]
                  - f_94 * di_66[k]
                  + f_95 * di_68[k]
                  + f_96 * di_113[k]
                  - f_85 * di_118[k]
                  - f_73 * di_120[k]
                  - f_94 * di_127[k]
                  + f_95 * di_129[k]
                  - f_98 * di_142[k]
                  + f_97 * di_147[k]
                  + f_99 * di_149[k]
                  + f_85 * di_156[k]
                  - f_74 * di_158[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dh_44, dh_49, dh_58, dh_86, dh_91, dh_100, dh_107, \
                         dh_112, dh_121, di_58, di_63, di_72, di_116, di_123, di_134, di_145, \
                         di_152, di_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_108 * ab_x[k] * dh_44[k]
                  - f_109 * ab_x[k] * dh_49[k]
                  + f_108 * ab_x[k] * dh_58[k]
                  + f_108 * ab_y[k] * dh_86[k]
                  - f_109 * ab_y[k] * dh_91[k]
                  + f_108 * ab_y[k] * dh_100[k]
                  - f_110 * ab_z[k] * dh_107[k]
                  + f_92 * ab_z[k] * dh_112[k]
                  - f_110 * ab_z[k] * dh_121[k]
                  - f_108 * di_58[k]
                  + f_109 * di_63[k]
                  - f_108 * di_72[k]
                  - f_108 * di_116[k]
                  + f_109 * di_123[k]
                  - f_108 * di_134[k]
                  + f_110 * di_145[k]
                  - f_92 * di_152[k]
                  + f_110 * di_163[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, dh_42, dh_45, dh_52, dh_84, dh_87, dh_94, dh_105, \
                         dh_108, dh_115, di_56, di_59, di_66, di_113, di_118, di_127, di_142, \
                         di_147, di_156 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_89 * ab_x[k] * dh_42[k]
                  - f_41 * ab_x[k] * dh_45[k]
                  + f_38 * ab_x[k] * dh_52[k]
                  + f_89 * ab_y[k] * dh_84[k]
                  - f_41 * ab_y[k] * dh_87[k]
                  + f_38 * ab_y[k] * dh_94[k]
                  - f_91 * ab_z[k] * dh_105[k]
                  + f_7 * ab_z[k] * dh_108[k]
                  - f_90 * ab_z[k] * dh_115[k]
                  - f_89 * di_56[k]
                  + f_41 * di_59[k]
                  - f_38 * di_66[k]
                  - f_89 * di_113[k]
                  + f_41 * di_118[k]
                  - f_38 * di_127[k]
                  + f_91 * di_142[k]
                  - f_7 * di_147[k]
                  + f_90 * di_156[k];
    }

#pragma omp simd aligned(ab_x, dh_1, dh_6, dh_15, dh_64, dh_69, dh_78, dh_106, dh_111, dh_120, \
                         di_1, di_6, di_15, di_85, di_90, di_99, di_141, di_146, \
                         di_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_56 * ab_x[k] * dh_1[k]
                  - f_57 * ab_x[k] * dh_6[k]
                  + f_58 * ab_x[k] * dh_15[k]
                  + f_56 * ab_x[k] * dh_64[k]
                  - f_57 * ab_x[k] * dh_69[k]
                  + f_58 * ab_x[k] * dh_78[k]
                  - f_54 * ab_x[k] * dh_106[k]
                  + f_59 * ab_x[k] * dh_111[k]
                  - f_60 * ab_x[k] * dh_120[k]
                  - f_56 * di_1[k]
                  + f_57 * di_6[k]
                  - f_58 * di_15[k]
                  - f_56 * di_85[k]
                  + f_57 * di_90[k]
                  - f_58 * di_99[k]
                  + f_54 * di_141[k]
                  - f_59 * di_146[k]
                  + f_60 * di_155[k];
    }

#pragma omp simd aligned(ab_x, dh_4, dh_11, dh_67, dh_74, dh_109, dh_116, di_4, di_11, di_88, \
                         di_95, di_144, di_151 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_61 * ab_x[k] * dh_4[k]
                  - f_61 * ab_x[k] * dh_11[k]
                  + f_61 * ab_x[k] * dh_67[k]
                  - f_61 * ab_x[k] * dh_74[k]
                  - f_62 * ab_x[k] * dh_109[k]
                  + f_62 * ab_x[k] * dh_116[k]
                  - f_61 * di_4[k]
                  + f_61 * di_11[k]
                  - f_61 * di_88[k]
                  + f_61 * di_95[k]
                  + f_62 * di_144[k]
                  - f_62 * di_151[k];
    }

#pragma omp simd aligned(ab_x, dh_1, dh_6, dh_8, dh_15, dh_17, dh_64, dh_69, dh_71, dh_78, \
                         dh_80, dh_106, dh_111, dh_113, dh_120, dh_122, di_1, di_6, di_8, \
                         di_15, di_17, di_85, di_90, di_92, di_99, di_101, di_141, di_146, \
                         di_148, di_155, di_157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_63 * ab_x[k] * dh_1[k]
                  - f_64 * ab_x[k] * dh_6[k]
                  + f_65 * ab_x[k] * dh_8[k]
                  + f_66 * ab_x[k] * dh_15[k]
                  - f_67 * ab_x[k] * dh_17[k]
                  - f_63 * ab_x[k] * dh_64[k]
                  - f_64 * ab_x[k] * dh_69[k]
                  + f_65 * ab_x[k] * dh_71[k]
                  + f_66 * ab_x[k] * dh_78[k]
                  - f_67 * ab_x[k] * dh_80[k]
                  + f_68 * ab_x[k] * dh_106[k]
                  + f_67 * ab_x[k] * dh_111[k]
                  - f_69 * ab_x[k] * dh_113[k]
                  - f_70 * ab_x[k] * dh_120[k]
                  + f_71 * ab_x[k] * dh_122[k]
                  + f_63 * di_1[k]
                  + f_64 * di_6[k]
                  - f_65 * di_8[k]
                  - f_66 * di_15[k]
                  + f_67 * di_17[k]
                  + f_63 * di_85[k]
                  + f_64 * di_90[k]
                  - f_65 * di_92[k]
                  - f_66 * di_99[k]
                  + f_67 * di_101[k]
                  - f_68 * di_141[k]
                  - f_67 * di_146[k]
                  + f_69 * di_148[k]
                  + f_70 * di_155[k]
                  - f_71 * di_157[k];
    }

#pragma omp simd aligned(ab_x, dh_4, dh_11, dh_13, dh_67, dh_74, dh_76, dh_109, dh_116, \
                         dh_118, di_4, di_11, di_13, di_88, di_95, di_97, di_144, di_151, \
                         di_153 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_72 * ab_x[k] * dh_4[k]
                  - f_72 * ab_x[k] * dh_11[k]
                  + f_73 * ab_x[k] * dh_13[k]
                  - f_72 * ab_x[k] * dh_67[k]
                  - f_72 * ab_x[k] * dh_74[k]
                  + f_73 * ab_x[k] * dh_76[k]
                  + f_74 * ab_x[k] * dh_109[k]
                  + f_74 * ab_x[k] * dh_116[k]
                  - f_75 * ab_x[k] * dh_118[k]
                  + f_72 * di_4[k]
                  + f_72 * di_11[k]
                  - f_73 * di_13[k]
                  + f_72 * di_88[k]
                  + f_72 * di_95[k]
                  - f_73 * di_97[k]
                  - f_74 * di_144[k]
                  - f_74 * di_151[k]
                  + f_75 * di_153[k];
    }

#pragma omp simd aligned(ab_x, dh_1, dh_6, dh_8, dh_15, dh_17, dh_19, dh_64, dh_69, dh_71, \
                         dh_78, dh_80, dh_82, dh_106, dh_111, dh_113, dh_120, dh_122, dh_124, \
                         di_1, di_6, di_8, di_15, di_17, di_19, di_85, di_90, di_92, di_99, \
                         di_101, di_103, di_141, di_146, di_148, di_155, di_157, \
                         di_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_76 * ab_x[k] * dh_1[k]
                  + f_77 * ab_x[k] * dh_6[k]
                  - f_78 * ab_x[k] * dh_8[k]
                  + f_76 * ab_x[k] * dh_15[k]
                  - f_78 * ab_x[k] * dh_17[k]
                  + f_31 * ab_x[k] * dh_19[k]
                  + f_76 * ab_x[k] * dh_64[k]
                  + f_77 * ab_x[k] * dh_69[k]
                  - f_78 * ab_x[k] * dh_71[k]
                  + f_76 * ab_x[k] * dh_78[k]
                  - f_78 * ab_x[k] * dh_80[k]
                  + f_31 * ab_x[k] * dh_82[k]
                  - f_79 * ab_x[k] * dh_106[k]
                  - f_31 * ab_x[k] * dh_111[k]
                  + f_80 * ab_x[k] * dh_113[k]
                  - f_79 * ab_x[k] * dh_120[k]
                  + f_80 * ab_x[k] * dh_122[k]
                  - f_81 * ab_x[k] * dh_124[k]
                  - f_76 * di_1[k]
                  - f_77 * di_6[k]
                  + f_78 * di_8[k]
                  - f_76 * di_15[k]
                  + f_78 * di_17[k]
                  - f_31 * di_19[k]
                  - f_76 * di_85[k]
                  - f_77 * di_90[k]
                  + f_78 * di_92[k]
                  - f_76 * di_99[k]
                  + f_78 * di_101[k]
                  - f_31 * di_103[k]
                  + f_79 * di_141[k]
                  + f_31 * di_146[k]
                  - f_80 * di_148[k]
                  + f_79 * di_155[k]
                  - f_80 * di_157[k]
                  + f_81 * di_159[k];
    }

#pragma omp simd aligned(ab_x, dh_2, dh_7, dh_9, dh_16, dh_18, dh_20, dh_65, dh_70, dh_72, \
                         dh_79, dh_81, dh_83, dh_107, dh_112, dh_114, dh_121, dh_123, dh_125, \
                         di_2, di_7, di_9, di_16, di_18, di_20, di_86, di_91, di_93, di_100, \
                         di_102, di_104, di_142, di_147, di_149, di_156, di_158, \
                         di_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_20 * ab_x[k] * dh_2[k]
                  + f_21 * ab_x[k] * dh_7[k]
                  - f_27 * ab_x[k] * dh_9[k]
                  + f_20 * ab_x[k] * dh_16[k]
                  - f_27 * ab_x[k] * dh_18[k]
                  + f_82 * ab_x[k] * dh_20[k]
                  + f_20 * ab_x[k] * dh_65[k]
                  + f_21 * ab_x[k] * dh_70[k]
                  - f_27 * ab_x[k] * dh_72[k]
                  + f_20 * ab_x[k] * dh_79[k]
                  - f_27 * ab_x[k] * dh_81[k]
                  + f_82 * ab_x[k] * dh_83[k]
                  - f_26 * ab_x[k] * dh_107[k]
                  - f_23 * ab_x[k] * dh_112[k]
                  + f_83 * ab_x[k] * dh_114[k]
                  - f_26 * ab_x[k] * dh_121[k]
                  + f_83 * ab_x[k] * dh_123[k]
                  - f_84 * ab_x[k] * dh_125[k]
                  - f_20 * di_2[k]
                  - f_21 * di_7[k]
                  + f_27 * di_9[k]
                  - f_20 * di_16[k]
                  + f_27 * di_18[k]
                  - f_82 * di_20[k]
                  - f_20 * di_86[k]
                  - f_21 * di_91[k]
                  + f_27 * di_93[k]
                  - f_20 * di_100[k]
                  + f_27 * di_102[k]
                  - f_82 * di_104[k]
                  + f_26 * di_142[k]
                  + f_23 * di_147[k]
                  - f_83 * di_149[k]
                  + f_26 * di_156[k]
                  - f_83 * di_158[k]
                  + f_84 * di_160[k];
    }

#pragma omp simd aligned(ab_x, dh_0, dh_3, dh_5, dh_10, dh_12, dh_14, dh_63, dh_66, dh_68, \
                         dh_73, dh_75, dh_77, dh_105, dh_108, dh_110, dh_115, dh_117, dh_119, \
                         di_0, di_3, di_5, di_10, di_12, di_14, di_84, di_87, di_89, di_94, \
                         di_96, di_98, di_140, di_143, di_145, di_150, di_152, \
                         di_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = f_76 * ab_x[k] * dh_0[k]
                  + f_77 * ab_x[k] * dh_3[k]
                  - f_78 * ab_x[k] * dh_5[k]
                  + f_76 * ab_x[k] * dh_10[k]
                  - f_78 * ab_x[k] * dh_12[k]
                  + f_31 * ab_x[k] * dh_14[k]
                  + f_76 * ab_x[k] * dh_63[k]
                  + f_77 * ab_x[k] * dh_66[k]
                  - f_78 * ab_x[k] * dh_68[k]
                  + f_76 * ab_x[k] * dh_73[k]
                  - f_78 * ab_x[k] * dh_75[k]
                  + f_31 * ab_x[k] * dh_77[k]
                  - f_79 * ab_x[k] * dh_105[k]
                  - f_31 * ab_x[k] * dh_108[k]
                  + f_80 * ab_x[k] * dh_110[k]
                  - f_79 * ab_x[k] * dh_115[k]
                  + f_80 * ab_x[k] * dh_117[k]
                  - f_81 * ab_x[k] * dh_119[k]
                  - f_76 * di_0[k]
                  - f_77 * di_3[k]
                  + f_78 * di_5[k]
                  - f_76 * di_10[k]
                  + f_78 * di_12[k]
                  - f_31 * di_14[k]
                  - f_76 * di_84[k]
                  - f_77 * di_87[k]
                  + f_78 * di_89[k]
                  - f_76 * di_94[k]
                  + f_78 * di_96[k]
                  - f_31 * di_98[k]
                  + f_79 * di_140[k]
                  + f_31 * di_143[k]
                  - f_80 * di_145[k]
                  + f_79 * di_150[k]
                  - f_80 * di_152[k]
                  + f_81 * di_154[k];
    }

#pragma omp simd aligned(ab_x, dh_2, dh_9, dh_16, dh_18, dh_65, dh_72, dh_79, dh_81, dh_107, \
                         dh_114, dh_121, dh_123, di_2, di_9, di_16, di_18, di_86, di_93, \
                         di_100, di_102, di_142, di_149, di_156, \
                         di_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_85 * ab_x[k] * dh_2[k]
                  + f_72 * ab_x[k] * dh_9[k]
                  + f_85 * ab_x[k] * dh_16[k]
                  - f_72 * ab_x[k] * dh_18[k]
                  - f_85 * ab_x[k] * dh_65[k]
                  + f_72 * ab_x[k] * dh_72[k]
                  + f_85 * ab_x[k] * dh_79[k]
                  - f_72 * ab_x[k] * dh_81[k]
                  + f_73 * ab_x[k] * dh_107[k]
                  - f_74 * ab_x[k] * dh_114[k]
                  - f_73 * ab_x[k] * dh_121[k]
                  + f_74 * ab_x[k] * dh_123[k]
                  + f_85 * di_2[k]
                  - f_72 * di_9[k]
                  - f_85 * di_16[k]
                  + f_72 * di_18[k]
                  + f_85 * di_86[k]
                  - f_72 * di_93[k]
                  - f_85 * di_100[k]
                  + f_72 * di_102[k]
                  - f_73 * di_142[k]
                  + f_74 * di_149[k]
                  + f_73 * di_156[k]
                  - f_74 * di_158[k];
    }

#pragma omp simd aligned(ab_x, dh_0, dh_3, dh_5, dh_10, dh_12, dh_63, dh_66, dh_68, dh_73, \
                         dh_75, dh_105, dh_108, dh_110, dh_115, dh_117, di_0, di_3, di_5, \
                         di_10, di_12, di_84, di_87, di_89, di_94, di_96, di_140, di_143, \
                         di_145, di_150, di_152 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_66 * ab_x[k] * dh_0[k]
                  + f_64 * ab_x[k] * dh_3[k]
                  + f_67 * ab_x[k] * dh_5[k]
                  + f_63 * ab_x[k] * dh_10[k]
                  - f_65 * ab_x[k] * dh_12[k]
                  - f_66 * ab_x[k] * dh_63[k]
                  + f_64 * ab_x[k] * dh_66[k]
                  + f_67 * ab_x[k] * dh_68[k]
                  + f_63 * ab_x[k] * dh_73[k]
                  - f_65 * ab_x[k] * dh_75[k]
                  + f_70 * ab_x[k] * dh_105[k]
                  - f_67 * ab_x[k] * dh_108[k]
                  - f_71 * ab_x[k] * dh_110[k]
                  - f_68 * ab_x[k] * dh_115[k]
                  + f_69 * ab_x[k] * dh_117[k]
                  + f_66 * di_0[k]
                  - f_64 * di_3[k]
                  - f_67 * di_5[k]
                  - f_63 * di_10[k]
                  + f_65 * di_12[k]
                  + f_66 * di_84[k]
                  - f_64 * di_87[k]
                  - f_67 * di_89[k]
                  - f_63 * di_94[k]
                  + f_65 * di_96[k]
                  - f_70 * di_140[k]
                  + f_67 * di_143[k]
                  + f_71 * di_145[k]
                  + f_68 * di_150[k]
                  - f_69 * di_152[k];
    }

#pragma omp simd aligned(ab_x, dh_2, dh_7, dh_16, dh_65, dh_70, dh_79, dh_107, dh_112, dh_121, \
                         di_2, di_7, di_16, di_86, di_91, di_100, di_142, di_147, \
                         di_156 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_86 * ab_x[k] * dh_2[k]
                  - f_87 * ab_x[k] * dh_7[k]
                  + f_86 * ab_x[k] * dh_16[k]
                  + f_86 * ab_x[k] * dh_65[k]
                  - f_87 * ab_x[k] * dh_70[k]
                  + f_86 * ab_x[k] * dh_79[k]
                  - f_61 * ab_x[k] * dh_107[k]
                  + f_88 * ab_x[k] * dh_112[k]
                  - f_61 * ab_x[k] * dh_121[k]
                  - f_86 * di_2[k]
                  + f_87 * di_7[k]
                  - f_86 * di_16[k]
                  - f_86 * di_86[k]
                  + f_87 * di_91[k]
                  - f_86 * di_100[k]
                  + f_61 * di_142[k]
                  - f_88 * di_147[k]
                  + f_61 * di_156[k];
    }

#pragma omp simd aligned(ab_x, dh_0, dh_3, dh_10, dh_63, dh_66, dh_73, dh_105, dh_108, dh_115, \
                         di_0, di_3, di_10, di_84, di_87, di_94, di_140, di_143, \
                         di_150 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_58 * ab_x[k] * dh_0[k]
                  - f_57 * ab_x[k] * dh_3[k]
                  + f_56 * ab_x[k] * dh_10[k]
                  + f_58 * ab_x[k] * dh_63[k]
                  - f_57 * ab_x[k] * dh_66[k]
                  + f_56 * ab_x[k] * dh_73[k]
                  - f_60 * ab_x[k] * dh_105[k]
                  + f_59 * ab_x[k] * dh_108[k]
                  - f_54 * ab_x[k] * dh_115[k]
                  - f_58 * di_0[k]
                  + f_57 * di_3[k]
                  - f_56 * di_10[k]
                  - f_58 * di_84[k]
                  + f_57 * di_87[k]
                  - f_56 * di_94[k]
                  + f_60 * di_140[k]
                  - f_59 * di_143[k]
                  + f_54 * di_150[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_43, dh_48, dh_57, dh_85, dh_90, dh_99, di_57, di_62, \
                         di_71, di_115, di_122, di_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_111 * ab_x[k] * dh_43[k]
                  + f_42 * ab_x[k] * dh_48[k]
                  - f_86 * ab_x[k] * dh_57[k]
                  + f_111 * ab_y[k] * dh_85[k]
                  - f_42 * ab_y[k] * dh_90[k]
                  + f_86 * ab_y[k] * dh_99[k]
                  + f_111 * di_57[k]
                  - f_42 * di_62[k]
                  + f_86 * di_71[k]
                  - f_111 * di_115[k]
                  + f_42 * di_122[k]
                  - f_86 * di_133[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_46, dh_53, dh_88, dh_95, di_60, di_67, di_119, \
                         di_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = -f_59 * ab_x[k] * dh_46[k]
                  + f_59 * ab_x[k] * dh_53[k]
                  + f_59 * ab_y[k] * dh_88[k]
                  - f_59 * ab_y[k] * dh_95[k]
                  + f_59 * di_60[k]
                  - f_59 * di_67[k]
                  - f_59 * di_119[k]
                  + f_59 * di_128[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_43, dh_48, dh_50, dh_57, dh_59, dh_85, dh_90, dh_92, \
                         dh_99, dh_101, di_57, di_62, di_64, di_71, di_73, di_115, di_122, \
                         di_124, di_133, di_135 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_112 * ab_x[k] * dh_43[k]
                  + f_37 * ab_x[k] * dh_48[k]
                  - f_17 * ab_x[k] * dh_50[k]
                  - f_113 * ab_x[k] * dh_57[k]
                  + f_19 * ab_x[k] * dh_59[k]
                  - f_112 * ab_y[k] * dh_85[k]
                  - f_37 * ab_y[k] * dh_90[k]
                  + f_17 * ab_y[k] * dh_92[k]
                  + f_113 * ab_y[k] * dh_99[k]
                  - f_19 * ab_y[k] * dh_101[k]
                  - f_112 * di_57[k]
                  - f_37 * di_62[k]
                  + f_17 * di_64[k]
                  + f_113 * di_71[k]
                  - f_19 * di_73[k]
                  + f_112 * di_115[k]
                  + f_37 * di_122[k]
                  - f_17 * di_124[k]
                  - f_113 * di_133[k]
                  + f_19 * di_135[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_46, dh_53, dh_55, dh_88, dh_95, dh_97, di_60, di_67, \
                         di_69, di_119, di_128, di_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_12 * ab_x[k] * dh_46[k]
                  + f_12 * ab_x[k] * dh_53[k]
                  - f_48 * ab_x[k] * dh_55[k]
                  - f_12 * ab_y[k] * dh_88[k]
                  - f_12 * ab_y[k] * dh_95[k]
                  + f_48 * ab_y[k] * dh_97[k]
                  - f_12 * di_60[k]
                  - f_12 * di_67[k]
                  + f_48 * di_69[k]
                  + f_12 * di_119[k]
                  + f_12 * di_128[k]
                  - f_48 * di_130[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_43, dh_48, dh_50, dh_57, dh_59, dh_61, dh_85, dh_90, \
                         dh_92, dh_99, dh_101, dh_103, di_57, di_62, di_64, di_71, di_73, \
                         di_75, di_115, di_122, di_124, di_133, di_135, \
                         di_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -0.9375 * ab_x[k] * dh_43[k]
                  - 1.875 * ab_x[k] * dh_48[k]
                  + 11.25 * ab_x[k] * dh_50[k]
                  - 0.9375 * ab_x[k] * dh_57[k]
                  + 11.25 * ab_x[k] * dh_59[k]
                  - 7.5 * ab_x[k] * dh_61[k]
                  + 0.9375 * ab_y[k] * dh_85[k]
                  + 1.875 * ab_y[k] * dh_90[k]
                  - 11.25 * ab_y[k] * dh_92[k]
                  + 0.9375 * ab_y[k] * dh_99[k]
                  - 11.25 * ab_y[k] * dh_101[k]
                  + 7.5 * ab_y[k] * dh_103[k]
                  + 0.9375 * di_57[k]
                  + 1.875 * di_62[k]
                  - 11.25 * di_64[k]
                  + 0.9375 * di_71[k]
                  - 11.25 * di_73[k]
                  + 7.5 * di_75[k]
                  - 0.9375 * di_115[k]
                  - 1.875 * di_122[k]
                  + 11.25 * di_124[k]
                  - 0.9375 * di_133[k]
                  + 11.25 * di_135[k]
                  - 7.5 * di_137[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_44, dh_49, dh_51, dh_58, dh_60, dh_62, dh_86, dh_91, \
                         dh_93, dh_100, dh_102, dh_104, di_58, di_63, di_65, di_72, di_74, \
                         di_76, di_116, di_123, di_125, di_134, di_136, \
                         di_138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_114 * ab_x[k] * dh_44[k]
                  - f_50 * ab_x[k] * dh_49[k]
                  + f_115 * ab_x[k] * dh_51[k]
                  - f_114 * ab_x[k] * dh_58[k]
                  + f_115 * ab_x[k] * dh_60[k]
                  - f_116 * ab_x[k] * dh_62[k]
                  + f_114 * ab_y[k] * dh_86[k]
                  + f_50 * ab_y[k] * dh_91[k]
                  - f_115 * ab_y[k] * dh_93[k]
                  + f_114 * ab_y[k] * dh_100[k]
                  - f_115 * ab_y[k] * dh_102[k]
                  + f_116 * ab_y[k] * dh_104[k]
                  + f_114 * di_58[k]
                  + f_50 * di_63[k]
                  - f_115 * di_65[k]
                  + f_114 * di_72[k]
                  - f_115 * di_74[k]
                  + f_116 * di_76[k]
                  - f_114 * di_116[k]
                  - f_50 * di_123[k]
                  + f_115 * di_125[k]
                  - f_114 * di_134[k]
                  + f_115 * di_136[k]
                  - f_116 * di_138[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_42, dh_45, dh_47, dh_52, dh_54, dh_56, dh_84, dh_87, \
                         dh_89, dh_94, dh_96, dh_98, di_56, di_59, di_61, di_66, di_68, di_70, \
                         di_113, di_118, di_120, di_127, di_129, \
                         di_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -0.9375 * ab_x[k] * dh_42[k]
                  - 1.875 * ab_x[k] * dh_45[k]
                  + 11.25 * ab_x[k] * dh_47[k]
                  - 0.9375 * ab_x[k] * dh_52[k]
                  + 11.25 * ab_x[k] * dh_54[k]
                  - 7.5 * ab_x[k] * dh_56[k]
                  + 0.9375 * ab_y[k] * dh_84[k]
                  + 1.875 * ab_y[k] * dh_87[k]
                  - 11.25 * ab_y[k] * dh_89[k]
                  + 0.9375 * ab_y[k] * dh_94[k]
                  - 11.25 * ab_y[k] * dh_96[k]
                  + 7.5 * ab_y[k] * dh_98[k]
                  + 0.9375 * di_56[k]
                  + 1.875 * di_59[k]
                  - 11.25 * di_61[k]
                  + 0.9375 * di_66[k]
                  - 11.25 * di_68[k]
                  + 7.5 * di_70[k]
                  - 0.9375 * di_113[k]
                  - 1.875 * di_118[k]
                  + 11.25 * di_120[k]
                  - 0.9375 * di_127[k]
                  + 11.25 * di_129[k]
                  - 7.5 * di_131[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_44, dh_51, dh_58, dh_60, dh_86, dh_93, dh_100, dh_102, \
                         di_58, di_65, di_72, di_74, di_116, di_125, di_134, \
                         di_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_117 * ab_x[k] * dh_44[k]
                  - f_12 * ab_x[k] * dh_51[k]
                  - f_117 * ab_x[k] * dh_58[k]
                  + f_12 * ab_x[k] * dh_60[k]
                  - f_117 * ab_y[k] * dh_86[k]
                  + f_12 * ab_y[k] * dh_93[k]
                  + f_117 * ab_y[k] * dh_100[k]
                  - f_12 * ab_y[k] * dh_102[k]
                  - f_117 * di_58[k]
                  + f_12 * di_65[k]
                  + f_117 * di_72[k]
                  - f_12 * di_74[k]
                  + f_117 * di_116[k]
                  - f_12 * di_125[k]
                  - f_117 * di_134[k]
                  + f_12 * di_136[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_42, dh_45, dh_47, dh_52, dh_54, dh_84, dh_87, dh_89, \
                         dh_94, dh_96, di_56, di_59, di_61, di_66, di_68, di_113, di_118, \
                         di_120, di_127, di_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_113 * ab_x[k] * dh_42[k]
                  - f_37 * ab_x[k] * dh_45[k]
                  - f_19 * ab_x[k] * dh_47[k]
                  - f_112 * ab_x[k] * dh_52[k]
                  + f_17 * ab_x[k] * dh_54[k]
                  - f_113 * ab_y[k] * dh_84[k]
                  + f_37 * ab_y[k] * dh_87[k]
                  + f_19 * ab_y[k] * dh_89[k]
                  + f_112 * ab_y[k] * dh_94[k]
                  - f_17 * ab_y[k] * dh_96[k]
                  - f_113 * di_56[k]
                  + f_37 * di_59[k]
                  + f_19 * di_61[k]
                  + f_112 * di_66[k]
                  - f_17 * di_68[k]
                  + f_113 * di_113[k]
                  - f_37 * di_118[k]
                  - f_19 * di_120[k]
                  - f_112 * di_127[k]
                  + f_17 * di_129[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_44, dh_49, dh_58, dh_86, dh_91, dh_100, di_58, di_63, \
                         di_72, di_116, di_123, di_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_57 * ab_x[k] * dh_44[k]
                  + f_118 * ab_x[k] * dh_49[k]
                  - f_57 * ab_x[k] * dh_58[k]
                  + f_57 * ab_y[k] * dh_86[k]
                  - f_118 * ab_y[k] * dh_91[k]
                  + f_57 * ab_y[k] * dh_100[k]
                  + f_57 * di_58[k]
                  - f_118 * di_63[k]
                  + f_57 * di_72[k]
                  - f_57 * di_116[k]
                  + f_118 * di_123[k]
                  - f_57 * di_134[k];
    }

#pragma omp simd aligned(ab_x, ab_y, dh_42, dh_45, dh_52, dh_84, dh_87, dh_94, di_56, di_59, \
                         di_66, di_113, di_118, di_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = -f_86 * ab_x[k] * dh_42[k]
                  + f_42 * ab_x[k] * dh_45[k]
                  - f_111 * ab_x[k] * dh_52[k]
                  + f_86 * ab_y[k] * dh_84[k]
                  - f_42 * ab_y[k] * dh_87[k]
                  + f_111 * ab_y[k] * dh_94[k]
                  + f_86 * di_56[k]
                  - f_42 * di_59[k]
                  + f_111 * di_66[k]
                  - f_86 * di_113[k]
                  + f_42 * di_118[k]
                  - f_111 * di_127[k];
    }

#pragma omp simd aligned(ab_x, dh_1, dh_6, dh_15, dh_64, dh_69, dh_78, di_1, di_6, di_15, \
                         di_85, di_90, di_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_3 * ab_x[k] * dh_1[k]
                  + f_4 * ab_x[k] * dh_6[k]
                  - f_5 * ab_x[k] * dh_15[k]
                  + f_0 * ab_x[k] * dh_64[k]
                  - f_1 * ab_x[k] * dh_69[k]
                  + f_2 * ab_x[k] * dh_78[k]
                  + f_3 * di_1[k]
                  - f_4 * di_6[k]
                  + f_5 * di_15[k]
                  - f_0 * di_85[k]
                  + f_1 * di_90[k]
                  - f_2 * di_99[k];
    }

#pragma omp simd aligned(ab_x, dh_4, dh_11, dh_67, dh_74, di_4, di_11, di_88, \
                         di_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_7 * ab_x[k] * dh_4[k]
                  + f_7 * ab_x[k] * dh_11[k]
                  + f_6 * ab_x[k] * dh_67[k]
                  - f_6 * ab_x[k] * dh_74[k]
                  + f_7 * di_4[k]
                  - f_7 * di_11[k]
                  - f_6 * di_88[k]
                  + f_6 * di_95[k];
    }

#pragma omp simd aligned(ab_x, dh_1, dh_6, dh_8, dh_15, dh_17, dh_64, dh_69, dh_71, dh_78, \
                         dh_80, di_1, di_6, di_8, di_15, di_17, di_85, di_90, di_92, di_99, \
                         di_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_11 * ab_x[k] * dh_1[k]
                  + f_13 * ab_x[k] * dh_6[k]
                  - f_12 * ab_x[k] * dh_8[k]
                  - f_14 * ab_x[k] * dh_15[k]
                  + f_15 * ab_x[k] * dh_17[k]
                  - f_8 * ab_x[k] * dh_64[k]
                  - f_9 * ab_x[k] * dh_69[k]
                  + f_10 * ab_x[k] * dh_71[k]
                  + f_11 * ab_x[k] * dh_78[k]
                  - f_12 * ab_x[k] * dh_80[k]
                  - f_11 * di_1[k]
                  - f_13 * di_6[k]
                  + f_12 * di_8[k]
                  + f_14 * di_15[k]
                  - f_15 * di_17[k]
                  + f_8 * di_85[k]
                  + f_9 * di_90[k]
                  - f_10 * di_92[k]
                  - f_11 * di_99[k]
                  + f_12 * di_101[k];
    }

#pragma omp simd aligned(ab_x, dh_4, dh_11, dh_13, dh_67, dh_74, dh_76, di_4, di_11, di_13, \
                         di_88, di_95, di_97 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_18 * ab_x[k] * dh_4[k]
                  + f_18 * ab_x[k] * dh_11[k]
                  - f_19 * ab_x[k] * dh_13[k]
                  - f_16 * ab_x[k] * dh_67[k]
                  - f_16 * ab_x[k] * dh_74[k]
                  + f_17 * ab_x[k] * dh_76[k]
                  - f_18 * di_4[k]
                  - f_18 * di_11[k]
                  + f_19 * di_13[k]
                  + f_16 * di_88[k]
                  + f_16 * di_95[k]
                  - f_17 * di_97[k];
    }

#pragma omp simd aligned(ab_x, dh_1, dh_6, dh_8, dh_15, dh_17, dh_19, dh_64, dh_69, dh_71, \
                         dh_78, dh_80, dh_82, di_1, di_6, di_8, di_15, di_17, di_19, di_85, \
                         di_90, di_92, di_99, di_101, di_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_24 * ab_x[k] * dh_1[k]
                  - f_25 * ab_x[k] * dh_6[k]
                  + f_26 * ab_x[k] * dh_8[k]
                  - f_24 * ab_x[k] * dh_15[k]
                  + f_26 * ab_x[k] * dh_17[k]
                  - f_27 * ab_x[k] * dh_19[k]
                  + f_20 * ab_x[k] * dh_64[k]
                  + f_21 * ab_x[k] * dh_69[k]
                  - f_22 * ab_x[k] * dh_71[k]
                  + f_20 * ab_x[k] * dh_78[k]
                  - f_22 * ab_x[k] * dh_80[k]
                  + f_23 * ab_x[k] * dh_82[k]
                  + f_24 * di_1[k]
                  + f_25 * di_6[k]
                  - f_26 * di_8[k]
                  + f_24 * di_15[k]
                  - f_26 * di_17[k]
                  + f_27 * di_19[k]
                  - f_20 * di_85[k]
                  - f_21 * di_90[k]
                  + f_22 * di_92[k]
                  - f_20 * di_99[k]
                  + f_22 * di_101[k]
                  - f_23 * di_103[k];
    }

#pragma omp simd aligned(ab_x, dh_2, dh_7, dh_9, dh_16, dh_18, dh_20, dh_65, dh_70, dh_72, \
                         dh_79, dh_81, dh_83, di_2, di_7, di_9, di_16, di_18, di_20, di_86, \
                         di_91, di_93, di_100, di_102, di_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_32 * ab_x[k] * dh_2[k]
                  - f_33 * ab_x[k] * dh_7[k]
                  + f_34 * ab_x[k] * dh_9[k]
                  - f_32 * ab_x[k] * dh_16[k]
                  + f_34 * ab_x[k] * dh_18[k]
                  - f_35 * ab_x[k] * dh_20[k]
                  + f_28 * ab_x[k] * dh_65[k]
                  + f_29 * ab_x[k] * dh_70[k]
                  - f_30 * ab_x[k] * dh_72[k]
                  + f_28 * ab_x[k] * dh_79[k]
                  - f_30 * ab_x[k] * dh_81[k]
                  + f_31 * ab_x[k] * dh_83[k]
                  + f_32 * di_2[k]
                  + f_33 * di_7[k]
                  - f_34 * di_9[k]
                  + f_32 * di_16[k]
                  - f_34 * di_18[k]
                  + f_35 * di_20[k]
                  - f_28 * di_86[k]
                  - f_29 * di_91[k]
                  + f_30 * di_93[k]
                  - f_28 * di_100[k]
                  + f_30 * di_102[k]
                  - f_31 * di_104[k];
    }

#pragma omp simd aligned(ab_x, dh_0, dh_3, dh_5, dh_10, dh_12, dh_14, dh_63, dh_66, dh_68, \
                         dh_73, dh_75, dh_77, di_0, di_3, di_5, di_10, di_12, di_14, di_84, \
                         di_87, di_89, di_94, di_96, di_98 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_24 * ab_x[k] * dh_0[k]
                  - f_25 * ab_x[k] * dh_3[k]
                  + f_26 * ab_x[k] * dh_5[k]
                  - f_24 * ab_x[k] * dh_10[k]
                  + f_26 * ab_x[k] * dh_12[k]
                  - f_27 * ab_x[k] * dh_14[k]
                  + f_20 * ab_x[k] * dh_63[k]
                  + f_21 * ab_x[k] * dh_66[k]
                  - f_22 * ab_x[k] * dh_68[k]
                  + f_20 * ab_x[k] * dh_73[k]
                  - f_22 * ab_x[k] * dh_75[k]
                  + f_23 * ab_x[k] * dh_77[k]
                  + f_24 * di_0[k]
                  + f_25 * di_3[k]
                  - f_26 * di_5[k]
                  + f_24 * di_10[k]
                  - f_26 * di_12[k]
                  + f_27 * di_14[k]
                  - f_20 * di_84[k]
                  - f_21 * di_87[k]
                  + f_22 * di_89[k]
                  - f_20 * di_94[k]
                  + f_22 * di_96[k]
                  - f_23 * di_98[k];
    }

#pragma omp simd aligned(ab_x, dh_2, dh_9, dh_16, dh_18, dh_65, dh_72, dh_79, dh_81, di_2, \
                         di_9, di_16, di_18, di_86, di_93, di_100, \
                         di_102 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_37 * ab_x[k] * dh_2[k]
                  - f_18 * ab_x[k] * dh_9[k]
                  - f_37 * ab_x[k] * dh_16[k]
                  + f_18 * ab_x[k] * dh_18[k]
                  - f_36 * ab_x[k] * dh_65[k]
                  + f_16 * ab_x[k] * dh_72[k]
                  + f_36 * ab_x[k] * dh_79[k]
                  - f_16 * ab_x[k] * dh_81[k]
                  - f_37 * di_2[k]
                  + f_18 * di_9[k]
                  + f_37 * di_16[k]
                  - f_18 * di_18[k]
                  + f_36 * di_86[k]
                  - f_16 * di_93[k]
                  - f_36 * di_100[k]
                  + f_16 * di_102[k];
    }

#pragma omp simd aligned(ab_x, dh_0, dh_3, dh_5, dh_10, dh_12, dh_63, dh_66, dh_68, dh_73, \
                         dh_75, di_0, di_3, di_5, di_10, di_12, di_84, di_87, di_89, di_94, \
                         di_96 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_14 * ab_x[k] * dh_0[k]
                  - f_13 * ab_x[k] * dh_3[k]
                  - f_15 * ab_x[k] * dh_5[k]
                  - f_11 * ab_x[k] * dh_10[k]
                  + f_12 * ab_x[k] * dh_12[k]
                  - f_11 * ab_x[k] * dh_63[k]
                  + f_9 * ab_x[k] * dh_66[k]
                  + f_12 * ab_x[k] * dh_68[k]
                  + f_8 * ab_x[k] * dh_73[k]
                  - f_10 * ab_x[k] * dh_75[k]
                  - f_14 * di_0[k]
                  + f_13 * di_3[k]
                  + f_15 * di_5[k]
                  + f_11 * di_10[k]
                  - f_12 * di_12[k]
                  + f_11 * di_84[k]
                  - f_9 * di_87[k]
                  - f_12 * di_89[k]
                  - f_8 * di_94[k]
                  + f_10 * di_96[k];
    }

#pragma omp simd aligned(ab_x, dh_2, dh_7, dh_16, dh_65, dh_70, dh_79, di_2, di_7, di_16, \
                         di_86, di_91, di_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_40 * ab_x[k] * dh_2[k]
                  + f_41 * ab_x[k] * dh_7[k]
                  - f_40 * ab_x[k] * dh_16[k]
                  + f_38 * ab_x[k] * dh_65[k]
                  - f_39 * ab_x[k] * dh_70[k]
                  + f_38 * ab_x[k] * dh_79[k]
                  + f_40 * di_2[k]
                  - f_41 * di_7[k]
                  + f_40 * di_16[k]
                  - f_38 * di_86[k]
                  + f_39 * di_91[k]
                  - f_38 * di_100[k];
    }

#pragma omp simd aligned(ab_x, dh_0, dh_3, dh_10, dh_63, dh_66, dh_73, di_0, di_3, di_10, \
                         di_84, di_87, di_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_5 * ab_x[k] * dh_0[k]
                  + f_4 * ab_x[k] * dh_3[k]
                  - f_3 * ab_x[k] * dh_10[k]
                  + f_2 * ab_x[k] * dh_63[k]
                  - f_1 * ab_x[k] * dh_66[k]
                  + f_0 * ab_x[k] * dh_73[k]
                  + f_5 * di_0[k]
                  - f_4 * di_3[k]
                  + f_3 * di_10[k]
                  - f_2 * di_84[k]
                  + f_1 * di_87[k]
                  - f_0 * di_94[k];
    }
}

auto
compute_hrr_fh(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t dh, const size_t di, const size_t nmax) -> void
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
    auto *t_150 = buffer.data(target + 150);
    auto *t_151 = buffer.data(target + 151);
    auto *t_152 = buffer.data(target + 152);
    auto *t_153 = buffer.data(target + 153);
    auto *t_154 = buffer.data(target + 154);
    auto *t_155 = buffer.data(target + 155);
    auto *t_156 = buffer.data(target + 156);
    auto *t_157 = buffer.data(target + 157);
    auto *t_158 = buffer.data(target + 158);
    auto *t_159 = buffer.data(target + 159);
    auto *t_160 = buffer.data(target + 160);
    auto *t_161 = buffer.data(target + 161);
    auto *t_162 = buffer.data(target + 162);
    auto *t_163 = buffer.data(target + 163);
    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);
    auto *t_168 = buffer.data(target + 168);
    auto *t_169 = buffer.data(target + 169);
    auto *t_170 = buffer.data(target + 170);
    auto *t_171 = buffer.data(target + 171);
    auto *t_172 = buffer.data(target + 172);
    auto *t_173 = buffer.data(target + 173);
    auto *t_174 = buffer.data(target + 174);
    auto *t_175 = buffer.data(target + 175);
    auto *t_176 = buffer.data(target + 176);
    auto *t_177 = buffer.data(target + 177);
    auto *t_178 = buffer.data(target + 178);
    auto *t_179 = buffer.data(target + 179);
    auto *t_180 = buffer.data(target + 180);
    auto *t_181 = buffer.data(target + 181);
    auto *t_182 = buffer.data(target + 182);
    auto *t_183 = buffer.data(target + 183);
    auto *t_184 = buffer.data(target + 184);
    auto *t_185 = buffer.data(target + 185);
    auto *t_186 = buffer.data(target + 186);
    auto *t_187 = buffer.data(target + 187);
    auto *t_188 = buffer.data(target + 188);
    auto *t_189 = buffer.data(target + 189);
    auto *t_190 = buffer.data(target + 190);
    auto *t_191 = buffer.data(target + 191);
    auto *t_192 = buffer.data(target + 192);
    auto *t_193 = buffer.data(target + 193);
    auto *t_194 = buffer.data(target + 194);
    auto *t_195 = buffer.data(target + 195);
    auto *t_196 = buffer.data(target + 196);
    auto *t_197 = buffer.data(target + 197);
    auto *t_198 = buffer.data(target + 198);
    auto *t_199 = buffer.data(target + 199);
    auto *t_200 = buffer.data(target + 200);
    auto *t_201 = buffer.data(target + 201);
    auto *t_202 = buffer.data(target + 202);
    auto *t_203 = buffer.data(target + 203);
    auto *t_204 = buffer.data(target + 204);
    auto *t_205 = buffer.data(target + 205);
    auto *t_206 = buffer.data(target + 206);
    auto *t_207 = buffer.data(target + 207);
    auto *t_208 = buffer.data(target + 208);
    auto *t_209 = buffer.data(target + 209);

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
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
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_41 = buffer.data(dh + 41);
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
    const auto *dh_57 = buffer.data(dh + 57);
    const auto *dh_58 = buffer.data(dh + 58);
    const auto *dh_59 = buffer.data(dh + 59);
    const auto *dh_60 = buffer.data(dh + 60);
    const auto *dh_61 = buffer.data(dh + 61);
    const auto *dh_62 = buffer.data(dh + 62);
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
    const auto *dh_83 = buffer.data(dh + 83);
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
    const auto *dh_104 = buffer.data(dh + 104);
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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_57 = buffer.data(di + 57);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_63 = buffer.data(di + 63);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_66 = buffer.data(di + 66);
    const auto *di_67 = buffer.data(di + 67);
    const auto *di_68 = buffer.data(di + 68);
    const auto *di_69 = buffer.data(di + 69);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_71 = buffer.data(di + 71);
    const auto *di_72 = buffer.data(di + 72);
    const auto *di_73 = buffer.data(di + 73);
    const auto *di_74 = buffer.data(di + 74);
    const auto *di_75 = buffer.data(di + 75);
    const auto *di_76 = buffer.data(di + 76);
    const auto *di_84 = buffer.data(di + 84);
    const auto *di_85 = buffer.data(di + 85);
    const auto *di_86 = buffer.data(di + 86);
    const auto *di_87 = buffer.data(di + 87);
    const auto *di_88 = buffer.data(di + 88);
    const auto *di_89 = buffer.data(di + 89);
    const auto *di_90 = buffer.data(di + 90);
    const auto *di_91 = buffer.data(di + 91);
    const auto *di_92 = buffer.data(di + 92);
    const auto *di_93 = buffer.data(di + 93);
    const auto *di_94 = buffer.data(di + 94);
    const auto *di_95 = buffer.data(di + 95);
    const auto *di_96 = buffer.data(di + 96);
    const auto *di_97 = buffer.data(di + 97);
    const auto *di_98 = buffer.data(di + 98);
    const auto *di_99 = buffer.data(di + 99);
    const auto *di_100 = buffer.data(di + 100);
    const auto *di_101 = buffer.data(di + 101);
    const auto *di_102 = buffer.data(di + 102);
    const auto *di_103 = buffer.data(di + 103);
    const auto *di_104 = buffer.data(di + 104);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_106 = buffer.data(di + 106);
    const auto *di_107 = buffer.data(di + 107);
    const auto *di_108 = buffer.data(di + 108);
    const auto *di_109 = buffer.data(di + 109);
    const auto *di_110 = buffer.data(di + 110);
    const auto *di_112 = buffer.data(di + 112);
    const auto *di_113 = buffer.data(di + 113);
    const auto *di_114 = buffer.data(di + 114);
    const auto *di_115 = buffer.data(di + 115);
    const auto *di_116 = buffer.data(di + 116);
    const auto *di_117 = buffer.data(di + 117);
    const auto *di_118 = buffer.data(di + 118);
    const auto *di_119 = buffer.data(di + 119);
    const auto *di_120 = buffer.data(di + 120);
    const auto *di_121 = buffer.data(di + 121);
    const auto *di_122 = buffer.data(di + 122);
    const auto *di_123 = buffer.data(di + 123);
    const auto *di_124 = buffer.data(di + 124);
    const auto *di_125 = buffer.data(di + 125);
    const auto *di_126 = buffer.data(di + 126);
    const auto *di_127 = buffer.data(di + 127);
    const auto *di_128 = buffer.data(di + 128);
    const auto *di_129 = buffer.data(di + 129);
    const auto *di_130 = buffer.data(di + 130);
    const auto *di_131 = buffer.data(di + 131);
    const auto *di_132 = buffer.data(di + 132);
    const auto *di_133 = buffer.data(di + 133);
    const auto *di_134 = buffer.data(di + 134);
    const auto *di_135 = buffer.data(di + 135);
    const auto *di_136 = buffer.data(di + 136);
    const auto *di_137 = buffer.data(di + 137);
    const auto *di_138 = buffer.data(di + 138);
    const auto *di_140 = buffer.data(di + 140);
    const auto *di_141 = buffer.data(di + 141);
    const auto *di_142 = buffer.data(di + 142);
    const auto *di_143 = buffer.data(di + 143);
    const auto *di_144 = buffer.data(di + 144);
    const auto *di_145 = buffer.data(di + 145);
    const auto *di_146 = buffer.data(di + 146);
    const auto *di_147 = buffer.data(di + 147);
    const auto *di_148 = buffer.data(di + 148);
    const auto *di_149 = buffer.data(di + 149);
    const auto *di_150 = buffer.data(di + 150);
    const auto *di_151 = buffer.data(di + 151);
    const auto *di_152 = buffer.data(di + 152);
    const auto *di_153 = buffer.data(di + 153);
    const auto *di_154 = buffer.data(di + 154);
    const auto *di_155 = buffer.data(di + 155);
    const auto *di_156 = buffer.data(di + 156);
    const auto *di_157 = buffer.data(di + 157);
    const auto *di_158 = buffer.data(di + 158);
    const auto *di_159 = buffer.data(di + 159);
    const auto *di_160 = buffer.data(di + 160);
    const auto *di_161 = buffer.data(di + 161);
    const auto *di_162 = buffer.data(di + 162);
    const auto *di_163 = buffer.data(di + 163);
    const auto *di_164 = buffer.data(di + 164);
    const auto *di_165 = buffer.data(di + 165);
    const auto *di_166 = buffer.data(di + 166);
    const auto *di_167 = buffer.data(di + 167);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, dh_0, dh_1, dh_2, dh_3, dh_4, di_0, \
                         di_1, di_2, di_3, di_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * dh_0[k]
                 + di_0[k];

        t_1[k] = -ab_x[k] * dh_1[k]
                 + di_1[k];

        t_2[k] = -ab_x[k] * dh_2[k]
                 + di_2[k];

        t_3[k] = -ab_x[k] * dh_3[k]
                 + di_3[k];

        t_4[k] = -ab_x[k] * dh_4[k]
                 + di_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, dh_5, dh_6, dh_7, dh_8, dh_9, di_5, \
                         di_6, di_7, di_8, di_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * dh_5[k]
                 + di_5[k];

        t_6[k] = -ab_x[k] * dh_6[k]
                 + di_6[k];

        t_7[k] = -ab_x[k] * dh_7[k]
                 + di_7[k];

        t_8[k] = -ab_x[k] * dh_8[k]
                 + di_8[k];

        t_9[k] = -ab_x[k] * dh_9[k]
                 + di_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, dh_10, dh_11, dh_12, dh_13, \
                         dh_14, di_10, di_11, di_12, di_13, di_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * dh_10[k]
                  + di_10[k];

        t_11[k] = -ab_x[k] * dh_11[k]
                  + di_11[k];

        t_12[k] = -ab_x[k] * dh_12[k]
                  + di_12[k];

        t_13[k] = -ab_x[k] * dh_13[k]
                  + di_13[k];

        t_14[k] = -ab_x[k] * dh_14[k]
                  + di_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, dh_15, dh_16, dh_17, dh_18, \
                         dh_19, di_15, di_16, di_17, di_18, di_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * dh_15[k]
                  + di_15[k];

        t_16[k] = -ab_x[k] * dh_16[k]
                  + di_16[k];

        t_17[k] = -ab_x[k] * dh_17[k]
                  + di_17[k];

        t_18[k] = -ab_x[k] * dh_18[k]
                  + di_18[k];

        t_19[k] = -ab_x[k] * dh_19[k]
                  + di_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, dh_20, dh_21, dh_22, dh_23, \
                         dh_24, di_20, di_28, di_29, di_30, di_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * dh_20[k]
                  + di_20[k];

        t_21[k] = -ab_x[k] * dh_21[k]
                  + di_28[k];

        t_22[k] = -ab_x[k] * dh_22[k]
                  + di_29[k];

        t_23[k] = -ab_x[k] * dh_23[k]
                  + di_30[k];

        t_24[k] = -ab_x[k] * dh_24[k]
                  + di_31[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, dh_25, dh_26, dh_27, dh_28, \
                         dh_29, di_32, di_33, di_34, di_35, di_36 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * dh_25[k]
                  + di_32[k];

        t_26[k] = -ab_x[k] * dh_26[k]
                  + di_33[k];

        t_27[k] = -ab_x[k] * dh_27[k]
                  + di_34[k];

        t_28[k] = -ab_x[k] * dh_28[k]
                  + di_35[k];

        t_29[k] = -ab_x[k] * dh_29[k]
                  + di_36[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, dh_30, dh_31, dh_32, dh_33, \
                         dh_34, di_37, di_38, di_39, di_40, di_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * dh_30[k]
                  + di_37[k];

        t_31[k] = -ab_x[k] * dh_31[k]
                  + di_38[k];

        t_32[k] = -ab_x[k] * dh_32[k]
                  + di_39[k];

        t_33[k] = -ab_x[k] * dh_33[k]
                  + di_40[k];

        t_34[k] = -ab_x[k] * dh_34[k]
                  + di_41[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, dh_35, dh_36, dh_37, dh_38, \
                         dh_39, di_42, di_43, di_44, di_45, di_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * dh_35[k]
                  + di_42[k];

        t_36[k] = -ab_x[k] * dh_36[k]
                  + di_43[k];

        t_37[k] = -ab_x[k] * dh_37[k]
                  + di_44[k];

        t_38[k] = -ab_x[k] * dh_38[k]
                  + di_45[k];

        t_39[k] = -ab_x[k] * dh_39[k]
                  + di_46[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, dh_40, dh_41, dh_42, dh_43, \
                         dh_44, di_47, di_48, di_56, di_57, di_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * dh_40[k]
                  + di_47[k];

        t_41[k] = -ab_x[k] * dh_41[k]
                  + di_48[k];

        t_42[k] = -ab_x[k] * dh_42[k]
                  + di_56[k];

        t_43[k] = -ab_x[k] * dh_43[k]
                  + di_57[k];

        t_44[k] = -ab_x[k] * dh_44[k]
                  + di_58[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, dh_45, dh_46, dh_47, dh_48, \
                         dh_49, di_59, di_60, di_61, di_62, di_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_x[k] * dh_45[k]
                  + di_59[k];

        t_46[k] = -ab_x[k] * dh_46[k]
                  + di_60[k];

        t_47[k] = -ab_x[k] * dh_47[k]
                  + di_61[k];

        t_48[k] = -ab_x[k] * dh_48[k]
                  + di_62[k];

        t_49[k] = -ab_x[k] * dh_49[k]
                  + di_63[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, dh_50, dh_51, dh_52, dh_53, \
                         dh_54, di_64, di_65, di_66, di_67, di_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_x[k] * dh_50[k]
                  + di_64[k];

        t_51[k] = -ab_x[k] * dh_51[k]
                  + di_65[k];

        t_52[k] = -ab_x[k] * dh_52[k]
                  + di_66[k];

        t_53[k] = -ab_x[k] * dh_53[k]
                  + di_67[k];

        t_54[k] = -ab_x[k] * dh_54[k]
                  + di_68[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, dh_55, dh_56, dh_57, dh_58, \
                         dh_59, di_69, di_70, di_71, di_72, di_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_x[k] * dh_55[k]
                  + di_69[k];

        t_56[k] = -ab_x[k] * dh_56[k]
                  + di_70[k];

        t_57[k] = -ab_x[k] * dh_57[k]
                  + di_71[k];

        t_58[k] = -ab_x[k] * dh_58[k]
                  + di_72[k];

        t_59[k] = -ab_x[k] * dh_59[k]
                  + di_73[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, dh_60, dh_61, dh_62, dh_63, \
                         dh_64, di_74, di_75, di_76, di_84, di_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_x[k] * dh_60[k]
                  + di_74[k];

        t_61[k] = -ab_x[k] * dh_61[k]
                  + di_75[k];

        t_62[k] = -ab_x[k] * dh_62[k]
                  + di_76[k];

        t_63[k] = -ab_x[k] * dh_63[k]
                  + di_84[k];

        t_64[k] = -ab_x[k] * dh_64[k]
                  + di_85[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, dh_65, dh_66, dh_67, dh_68, \
                         dh_69, di_86, di_87, di_88, di_89, di_90 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = -ab_x[k] * dh_65[k]
                  + di_86[k];

        t_66[k] = -ab_x[k] * dh_66[k]
                  + di_87[k];

        t_67[k] = -ab_x[k] * dh_67[k]
                  + di_88[k];

        t_68[k] = -ab_x[k] * dh_68[k]
                  + di_89[k];

        t_69[k] = -ab_x[k] * dh_69[k]
                  + di_90[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, dh_70, dh_71, dh_72, dh_73, \
                         dh_74, di_91, di_92, di_93, di_94, di_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = -ab_x[k] * dh_70[k]
                  + di_91[k];

        t_71[k] = -ab_x[k] * dh_71[k]
                  + di_92[k];

        t_72[k] = -ab_x[k] * dh_72[k]
                  + di_93[k];

        t_73[k] = -ab_x[k] * dh_73[k]
                  + di_94[k];

        t_74[k] = -ab_x[k] * dh_74[k]
                  + di_95[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, dh_75, dh_76, dh_77, dh_78, \
                         dh_79, di_96, di_97, di_98, di_99, di_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = -ab_x[k] * dh_75[k]
                  + di_96[k];

        t_76[k] = -ab_x[k] * dh_76[k]
                  + di_97[k];

        t_77[k] = -ab_x[k] * dh_77[k]
                  + di_98[k];

        t_78[k] = -ab_x[k] * dh_78[k]
                  + di_99[k];

        t_79[k] = -ab_x[k] * dh_79[k]
                  + di_100[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, dh_80, dh_81, dh_82, dh_83, \
                         dh_84, di_101, di_102, di_103, di_104, \
                         di_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = -ab_x[k] * dh_80[k]
                  + di_101[k];

        t_81[k] = -ab_x[k] * dh_81[k]
                  + di_102[k];

        t_82[k] = -ab_x[k] * dh_82[k]
                  + di_103[k];

        t_83[k] = -ab_x[k] * dh_83[k]
                  + di_104[k];

        t_84[k] = -ab_x[k] * dh_84[k]
                  + di_112[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, dh_85, dh_86, dh_87, dh_88, \
                         dh_89, di_113, di_114, di_115, di_116, \
                         di_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = -ab_x[k] * dh_85[k]
                  + di_113[k];

        t_86[k] = -ab_x[k] * dh_86[k]
                  + di_114[k];

        t_87[k] = -ab_x[k] * dh_87[k]
                  + di_115[k];

        t_88[k] = -ab_x[k] * dh_88[k]
                  + di_116[k];

        t_89[k] = -ab_x[k] * dh_89[k]
                  + di_117[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, dh_90, dh_91, dh_92, dh_93, \
                         dh_94, di_118, di_119, di_120, di_121, \
                         di_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = -ab_x[k] * dh_90[k]
                  + di_118[k];

        t_91[k] = -ab_x[k] * dh_91[k]
                  + di_119[k];

        t_92[k] = -ab_x[k] * dh_92[k]
                  + di_120[k];

        t_93[k] = -ab_x[k] * dh_93[k]
                  + di_121[k];

        t_94[k] = -ab_x[k] * dh_94[k]
                  + di_122[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, dh_95, dh_96, dh_97, dh_98, \
                         dh_99, di_123, di_124, di_125, di_126, \
                         di_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = -ab_x[k] * dh_95[k]
                  + di_123[k];

        t_96[k] = -ab_x[k] * dh_96[k]
                  + di_124[k];

        t_97[k] = -ab_x[k] * dh_97[k]
                  + di_125[k];

        t_98[k] = -ab_x[k] * dh_98[k]
                  + di_126[k];

        t_99[k] = -ab_x[k] * dh_99[k]
                  + di_127[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, dh_100, dh_101, dh_102, \
                         dh_103, dh_104, di_128, di_129, di_130, di_131, \
                         di_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = -ab_x[k] * dh_100[k]
                   + di_128[k];

        t_101[k] = -ab_x[k] * dh_101[k]
                   + di_129[k];

        t_102[k] = -ab_x[k] * dh_102[k]
                   + di_130[k];

        t_103[k] = -ab_x[k] * dh_103[k]
                   + di_131[k];

        t_104[k] = -ab_x[k] * dh_104[k]
                   + di_132[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, dh_105, dh_106, dh_107, \
                         dh_108, dh_109, di_140, di_141, di_142, di_143, \
                         di_144 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = -ab_x[k] * dh_105[k]
                   + di_140[k];

        t_106[k] = -ab_x[k] * dh_106[k]
                   + di_141[k];

        t_107[k] = -ab_x[k] * dh_107[k]
                   + di_142[k];

        t_108[k] = -ab_x[k] * dh_108[k]
                   + di_143[k];

        t_109[k] = -ab_x[k] * dh_109[k]
                   + di_144[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, dh_110, dh_111, dh_112, \
                         dh_113, dh_114, di_145, di_146, di_147, di_148, \
                         di_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = -ab_x[k] * dh_110[k]
                   + di_145[k];

        t_111[k] = -ab_x[k] * dh_111[k]
                   + di_146[k];

        t_112[k] = -ab_x[k] * dh_112[k]
                   + di_147[k];

        t_113[k] = -ab_x[k] * dh_113[k]
                   + di_148[k];

        t_114[k] = -ab_x[k] * dh_114[k]
                   + di_149[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, dh_115, dh_116, dh_117, \
                         dh_118, dh_119, di_150, di_151, di_152, di_153, \
                         di_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = -ab_x[k] * dh_115[k]
                   + di_150[k];

        t_116[k] = -ab_x[k] * dh_116[k]
                   + di_151[k];

        t_117[k] = -ab_x[k] * dh_117[k]
                   + di_152[k];

        t_118[k] = -ab_x[k] * dh_118[k]
                   + di_153[k];

        t_119[k] = -ab_x[k] * dh_119[k]
                   + di_154[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, dh_120, dh_121, dh_122, \
                         dh_123, dh_124, di_155, di_156, di_157, di_158, \
                         di_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = -ab_x[k] * dh_120[k]
                   + di_155[k];

        t_121[k] = -ab_x[k] * dh_121[k]
                   + di_156[k];

        t_122[k] = -ab_x[k] * dh_122[k]
                   + di_157[k];

        t_123[k] = -ab_x[k] * dh_123[k]
                   + di_158[k];

        t_124[k] = -ab_x[k] * dh_124[k]
                   + di_159[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, ab_x, ab_y, dh_63, dh_64, dh_65, dh_125, \
                         di_85, di_87, di_88, di_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = -ab_x[k] * dh_125[k]
                   + di_160[k];

        t_126[k] = -ab_y[k] * dh_63[k]
                   + di_85[k];

        t_127[k] = -ab_y[k] * dh_64[k]
                   + di_87[k];

        t_128[k] = -ab_y[k] * dh_65[k]
                   + di_88[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, ab_y, dh_66, dh_67, dh_68, dh_69, \
                         dh_70, di_90, di_91, di_92, di_94, di_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_129[k] = -ab_y[k] * dh_66[k]
                   + di_90[k];

        t_130[k] = -ab_y[k] * dh_67[k]
                   + di_91[k];

        t_131[k] = -ab_y[k] * dh_68[k]
                   + di_92[k];

        t_132[k] = -ab_y[k] * dh_69[k]
                   + di_94[k];

        t_133[k] = -ab_y[k] * dh_70[k]
                   + di_95[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, ab_y, dh_71, dh_72, dh_73, dh_74, \
                         dh_75, di_96, di_97, di_99, di_100, di_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_134[k] = -ab_y[k] * dh_71[k]
                   + di_96[k];

        t_135[k] = -ab_y[k] * dh_72[k]
                   + di_97[k];

        t_136[k] = -ab_y[k] * dh_73[k]
                   + di_99[k];

        t_137[k] = -ab_y[k] * dh_74[k]
                   + di_100[k];

        t_138[k] = -ab_y[k] * dh_75[k]
                   + di_101[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_y, dh_76, dh_77, dh_78, dh_79, \
                         dh_80, di_102, di_103, di_105, di_106, \
                         di_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_139[k] = -ab_y[k] * dh_76[k]
                   + di_102[k];

        t_140[k] = -ab_y[k] * dh_77[k]
                   + di_103[k];

        t_141[k] = -ab_y[k] * dh_78[k]
                   + di_105[k];

        t_142[k] = -ab_y[k] * dh_79[k]
                   + di_106[k];

        t_143[k] = -ab_y[k] * dh_80[k]
                   + di_107[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, ab_y, dh_81, dh_82, dh_83, dh_84, \
                         dh_85, di_108, di_109, di_110, di_113, \
                         di_115 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_144[k] = -ab_y[k] * dh_81[k]
                   + di_108[k];

        t_145[k] = -ab_y[k] * dh_82[k]
                   + di_109[k];

        t_146[k] = -ab_y[k] * dh_83[k]
                   + di_110[k];

        t_147[k] = -ab_y[k] * dh_84[k]
                   + di_113[k];

        t_148[k] = -ab_y[k] * dh_85[k]
                   + di_115[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, ab_y, dh_86, dh_87, dh_88, dh_89, \
                         dh_90, di_116, di_118, di_119, di_120, \
                         di_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_149[k] = -ab_y[k] * dh_86[k]
                   + di_116[k];

        t_150[k] = -ab_y[k] * dh_87[k]
                   + di_118[k];

        t_151[k] = -ab_y[k] * dh_88[k]
                   + di_119[k];

        t_152[k] = -ab_y[k] * dh_89[k]
                   + di_120[k];

        t_153[k] = -ab_y[k] * dh_90[k]
                   + di_122[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, ab_y, dh_91, dh_92, dh_93, dh_94, \
                         dh_95, di_123, di_124, di_125, di_127, \
                         di_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_154[k] = -ab_y[k] * dh_91[k]
                   + di_123[k];

        t_155[k] = -ab_y[k] * dh_92[k]
                   + di_124[k];

        t_156[k] = -ab_y[k] * dh_93[k]
                   + di_125[k];

        t_157[k] = -ab_y[k] * dh_94[k]
                   + di_127[k];

        t_158[k] = -ab_y[k] * dh_95[k]
                   + di_128[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, ab_y, dh_96, dh_97, dh_98, dh_99, \
                         dh_100, di_129, di_130, di_131, di_133, \
                         di_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_159[k] = -ab_y[k] * dh_96[k]
                   + di_129[k];

        t_160[k] = -ab_y[k] * dh_97[k]
                   + di_130[k];

        t_161[k] = -ab_y[k] * dh_98[k]
                   + di_131[k];

        t_162[k] = -ab_y[k] * dh_99[k]
                   + di_133[k];

        t_163[k] = -ab_y[k] * dh_100[k]
                   + di_134[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, ab_y, dh_101, dh_102, dh_103, \
                         dh_104, dh_105, di_135, di_136, di_137, di_138, \
                         di_141 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_164[k] = -ab_y[k] * dh_101[k]
                   + di_135[k];

        t_165[k] = -ab_y[k] * dh_102[k]
                   + di_136[k];

        t_166[k] = -ab_y[k] * dh_103[k]
                   + di_137[k];

        t_167[k] = -ab_y[k] * dh_104[k]
                   + di_138[k];

        t_168[k] = -ab_y[k] * dh_105[k]
                   + di_141[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, ab_y, dh_106, dh_107, dh_108, \
                         dh_109, dh_110, di_143, di_144, di_146, di_147, \
                         di_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_169[k] = -ab_y[k] * dh_106[k]
                   + di_143[k];

        t_170[k] = -ab_y[k] * dh_107[k]
                   + di_144[k];

        t_171[k] = -ab_y[k] * dh_108[k]
                   + di_146[k];

        t_172[k] = -ab_y[k] * dh_109[k]
                   + di_147[k];

        t_173[k] = -ab_y[k] * dh_110[k]
                   + di_148[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, ab_y, dh_111, dh_112, dh_113, \
                         dh_114, dh_115, di_150, di_151, di_152, di_153, \
                         di_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_174[k] = -ab_y[k] * dh_111[k]
                   + di_150[k];

        t_175[k] = -ab_y[k] * dh_112[k]
                   + di_151[k];

        t_176[k] = -ab_y[k] * dh_113[k]
                   + di_152[k];

        t_177[k] = -ab_y[k] * dh_114[k]
                   + di_153[k];

        t_178[k] = -ab_y[k] * dh_115[k]
                   + di_155[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, ab_y, dh_116, dh_117, dh_118, \
                         dh_119, dh_120, di_156, di_157, di_158, di_159, \
                         di_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_179[k] = -ab_y[k] * dh_116[k]
                   + di_156[k];

        t_180[k] = -ab_y[k] * dh_117[k]
                   + di_157[k];

        t_181[k] = -ab_y[k] * dh_118[k]
                   + di_158[k];

        t_182[k] = -ab_y[k] * dh_119[k]
                   + di_159[k];

        t_183[k] = -ab_y[k] * dh_120[k]
                   + di_161[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, ab_y, dh_121, dh_122, dh_123, \
                         dh_124, dh_125, di_162, di_163, di_164, di_165, \
                         di_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_184[k] = -ab_y[k] * dh_121[k]
                   + di_162[k];

        t_185[k] = -ab_y[k] * dh_122[k]
                   + di_163[k];

        t_186[k] = -ab_y[k] * dh_123[k]
                   + di_164[k];

        t_187[k] = -ab_y[k] * dh_124[k]
                   + di_165[k];

        t_188[k] = -ab_y[k] * dh_125[k]
                   + di_166[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, ab_z, dh_105, dh_106, dh_107, \
                         dh_108, dh_109, di_142, di_144, di_145, di_147, \
                         di_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_189[k] = -ab_z[k] * dh_105[k]
                   + di_142[k];

        t_190[k] = -ab_z[k] * dh_106[k]
                   + di_144[k];

        t_191[k] = -ab_z[k] * dh_107[k]
                   + di_145[k];

        t_192[k] = -ab_z[k] * dh_108[k]
                   + di_147[k];

        t_193[k] = -ab_z[k] * dh_109[k]
                   + di_148[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, ab_z, dh_110, dh_111, dh_112, \
                         dh_113, dh_114, di_149, di_151, di_152, di_153, \
                         di_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_194[k] = -ab_z[k] * dh_110[k]
                   + di_149[k];

        t_195[k] = -ab_z[k] * dh_111[k]
                   + di_151[k];

        t_196[k] = -ab_z[k] * dh_112[k]
                   + di_152[k];

        t_197[k] = -ab_z[k] * dh_113[k]
                   + di_153[k];

        t_198[k] = -ab_z[k] * dh_114[k]
                   + di_154[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, ab_z, dh_115, dh_116, dh_117, \
                         dh_118, dh_119, di_156, di_157, di_158, di_159, \
                         di_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_199[k] = -ab_z[k] * dh_115[k]
                   + di_156[k];

        t_200[k] = -ab_z[k] * dh_116[k]
                   + di_157[k];

        t_201[k] = -ab_z[k] * dh_117[k]
                   + di_158[k];

        t_202[k] = -ab_z[k] * dh_118[k]
                   + di_159[k];

        t_203[k] = -ab_z[k] * dh_119[k]
                   + di_160[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, ab_z, dh_120, dh_121, dh_122, \
                         dh_123, dh_124, di_162, di_163, di_164, di_165, \
                         di_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_204[k] = -ab_z[k] * dh_120[k]
                   + di_162[k];

        t_205[k] = -ab_z[k] * dh_121[k]
                   + di_163[k];

        t_206[k] = -ab_z[k] * dh_122[k]
                   + di_164[k];

        t_207[k] = -ab_z[k] * dh_123[k]
                   + di_165[k];

        t_208[k] = -ab_z[k] * dh_124[k]
                   + di_166[k];
    }

#pragma omp simd aligned(t_209, ab_z, dh_125, di_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_209[k] = -ab_z[k] * dh_125[k]
                   + di_167[k];
    }
}

}  // namespace simdovl
