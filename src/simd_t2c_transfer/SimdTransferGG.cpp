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


#include "SimdTransferGG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hrr_gg_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t fg, const size_t fh,
                   const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 13.125 * std::sqrt(2.0);
    const auto f_1 = 4.375 * std::sqrt(2.0);
    const auto f_2 = 1.25 * std::sqrt(7.0);
    const auto f_3 = 7.5 * std::sqrt(7.0);
    const auto f_4 = 1.875 * std::sqrt(14.0);
    const auto f_5 = 2.5 * std::sqrt(14.0);
    const auto f_6 = 0.1875 * std::sqrt(35.0);
    const auto f_7 = 0.375 * std::sqrt(35.0);
    const auto f_8 = 1.5 * std::sqrt(35.0);
    const auto f_9 = 0.5 * std::sqrt(35.0);
    const auto f_10 = 0.625 * std::sqrt(7.0);
    const auto f_11 = 3.75 * std::sqrt(7.0);
    const auto f_12 = 11.25 * std::sqrt(14.0);
    const auto f_13 = 0.625 * std::sqrt(14.0);
    const auto f_14 = 3.75 * std::sqrt(14.0);
    const auto f_15 = 5.625 * std::sqrt(7.0);
    const auto f_16 = 1.875 * std::sqrt(7.0);
    const auto f_17 = 2.5 * std::sqrt(7.0);
    const auto f_18 = 0.28125 * std::sqrt(70.0);
    const auto f_19 = 0.5625 * std::sqrt(70.0);
    const auto f_20 = 2.25 * std::sqrt(70.0);
    const auto f_21 = 0.75 * std::sqrt(70.0);
    const auto f_22 = 0.09375 * std::sqrt(70.0);
    const auto f_23 = 0.1875 * std::sqrt(70.0);
    const auto f_24 = 0.25 * std::sqrt(70.0);
    const auto f_25 = 0.9375 * std::sqrt(14.0);
    const auto f_26 = 5.625 * std::sqrt(14.0);
    const auto f_27 = 0.3125 * std::sqrt(14.0);
    const auto f_28 = 3.28125 * std::sqrt(2.0);
    const auto f_29 = 19.6875 * std::sqrt(2.0);
    const auto f_30 = 1.09375 * std::sqrt(2.0);
    const auto f_31 = 6.5625 * std::sqrt(2.0);
    const auto f_32 = 1.875 * std::sqrt(2.0);
    const auto f_33 = 2.5 * std::sqrt(2.0);
    const auto f_34 = 11.25 * std::sqrt(2.0);
    const auto f_35 = 15.0 * std::sqrt(2.0);
    const auto f_36 = 0.1875 * std::sqrt(5.0);
    const auto f_37 = 0.375 * std::sqrt(5.0);
    const auto f_38 = 1.5 * std::sqrt(5.0);
    const auto f_39 = 0.5 * std::sqrt(5.0);
    const auto f_40 = 1.125 * std::sqrt(5.0);
    const auto f_41 = 2.25 * std::sqrt(5.0);
    const auto f_42 = 9.0 * std::sqrt(5.0);
    const auto f_43 = 3.0 * std::sqrt(5.0);
    const auto f_44 = 0.3125 * std::sqrt(7.0);
    const auto f_45 = 11.25 * std::sqrt(7.0);
    const auto f_46 = 0.28125 * std::sqrt(10.0);
    const auto f_47 = 0.5625 * std::sqrt(10.0);
    const auto f_48 = 2.25 * std::sqrt(10.0);
    const auto f_49 = 0.75 * std::sqrt(10.0);
    const auto f_50 = 0.375 * std::sqrt(10.0);
    const auto f_51 = 3.0 * std::sqrt(10.0);
    const auto f_52 = std::sqrt(10.0);
    const auto f_53 = 0.9375 * std::sqrt(2.0);
    const auto f_54 = 5.625 * std::sqrt(2.0);
    const auto f_55 = 1.25 * std::sqrt(2.0);
    const auto f_56 = 7.5 * std::sqrt(2.0);
    const auto f_57 = 0.46875 * std::sqrt(14.0);
    const auto f_58 = 2.8125 * std::sqrt(14.0);
    const auto f_59 = 0.09375 * std::sqrt(5.0);
    const auto f_60 = 0.5625 * std::sqrt(5.0);
    const auto f_61 = 0.75 * std::sqrt(5.0);
    const auto f_62 = 4.5 * std::sqrt(5.0);
    const auto f_63 = 0.25 * std::sqrt(5.0);
    const auto f_64 = 0.046875 * std::sqrt(35.0);
    const auto f_65 = 0.28125 * std::sqrt(35.0);
    const auto f_66 = 0.09375 * std::sqrt(35.0);
    const auto f_67 = 0.5625 * std::sqrt(35.0);
    const auto f_68 = 2.25 * std::sqrt(35.0);
    const auto f_69 = 0.125 * std::sqrt(35.0);
    const auto f_70 = 0.75 * std::sqrt(35.0);
    const auto f_71 = 0.15625 * std::sqrt(7.0);
    const auto f_72 = 0.9375 * std::sqrt(7.0);

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
    auto *g_77 = values + 77 * nvalues;
    auto *g_78 = values + 78 * nvalues;
    auto *g_79 = values + 79 * nvalues;
    auto *g_80 = values + 80 * nvalues;

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_61 = buffer.data(fg + 61);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_64 = buffer.data(fg + 64);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_66 = buffer.data(fg + 66);
    const auto *fg_67 = buffer.data(fg + 67);
    const auto *fg_68 = buffer.data(fg + 68);
    const auto *fg_69 = buffer.data(fg + 69);
    const auto *fg_70 = buffer.data(fg + 70);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_74 = buffer.data(fg + 74);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_76 = buffer.data(fg + 76);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_78 = buffer.data(fg + 78);
    const auto *fg_79 = buffer.data(fg + 79);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_81 = buffer.data(fg + 81);
    const auto *fg_82 = buffer.data(fg + 82);
    const auto *fg_83 = buffer.data(fg + 83);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_85 = buffer.data(fg + 85);
    const auto *fg_86 = buffer.data(fg + 86);
    const auto *fg_87 = buffer.data(fg + 87);
    const auto *fg_88 = buffer.data(fg + 88);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_92 = buffer.data(fg + 92);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_94 = buffer.data(fg + 94);
    const auto *fg_95 = buffer.data(fg + 95);
    const auto *fg_96 = buffer.data(fg + 96);
    const auto *fg_97 = buffer.data(fg + 97);
    const auto *fg_98 = buffer.data(fg + 98);
    const auto *fg_99 = buffer.data(fg + 99);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_101 = buffer.data(fg + 101);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_103 = buffer.data(fg + 103);
    const auto *fg_104 = buffer.data(fg + 104);
    const auto *fg_105 = buffer.data(fg + 105);
    const auto *fg_106 = buffer.data(fg + 106);
    const auto *fg_107 = buffer.data(fg + 107);
    const auto *fg_108 = buffer.data(fg + 108);
    const auto *fg_109 = buffer.data(fg + 109);
    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_111 = buffer.data(fg + 111);
    const auto *fg_112 = buffer.data(fg + 112);
    const auto *fg_113 = buffer.data(fg + 113);
    const auto *fg_114 = buffer.data(fg + 114);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_120 = buffer.data(fg + 120);
    const auto *fg_121 = buffer.data(fg + 121);
    const auto *fg_122 = buffer.data(fg + 122);
    const auto *fg_123 = buffer.data(fg + 123);
    const auto *fg_124 = buffer.data(fg + 124);
    const auto *fg_125 = buffer.data(fg + 125);
    const auto *fg_126 = buffer.data(fg + 126);
    const auto *fg_127 = buffer.data(fg + 127);
    const auto *fg_128 = buffer.data(fg + 128);
    const auto *fg_129 = buffer.data(fg + 129);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_134 = buffer.data(fg + 134);
    const auto *fg_135 = buffer.data(fg + 135);
    const auto *fg_136 = buffer.data(fg + 136);
    const auto *fg_137 = buffer.data(fg + 137);
    const auto *fg_138 = buffer.data(fg + 138);
    const auto *fg_139 = buffer.data(fg + 139);
    const auto *fg_140 = buffer.data(fg + 140);
    const auto *fg_141 = buffer.data(fg + 141);
    const auto *fg_142 = buffer.data(fg + 142);
    const auto *fg_143 = buffer.data(fg + 143);
    const auto *fg_144 = buffer.data(fg + 144);
    const auto *fg_145 = buffer.data(fg + 145);
    const auto *fg_146 = buffer.data(fg + 146);
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_148 = buffer.data(fg + 148);
    const auto *fg_149 = buffer.data(fg + 149);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_48 = buffer.data(fh + 48);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_53 = buffer.data(fh + 53);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_55 = buffer.data(fh + 55);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_70 = buffer.data(fh + 70);
    const auto *fh_71 = buffer.data(fh + 71);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_74 = buffer.data(fh + 74);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_76 = buffer.data(fh + 76);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_84 = buffer.data(fh + 84);
    const auto *fh_85 = buffer.data(fh + 85);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_88 = buffer.data(fh + 88);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_91 = buffer.data(fh + 91);
    const auto *fh_92 = buffer.data(fh + 92);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_94 = buffer.data(fh + 94);
    const auto *fh_95 = buffer.data(fh + 95);
    const auto *fh_96 = buffer.data(fh + 96);
    const auto *fh_97 = buffer.data(fh + 97);
    const auto *fh_98 = buffer.data(fh + 98);
    const auto *fh_105 = buffer.data(fh + 105);
    const auto *fh_106 = buffer.data(fh + 106);
    const auto *fh_107 = buffer.data(fh + 107);
    const auto *fh_108 = buffer.data(fh + 108);
    const auto *fh_109 = buffer.data(fh + 109);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_111 = buffer.data(fh + 111);
    const auto *fh_112 = buffer.data(fh + 112);
    const auto *fh_113 = buffer.data(fh + 113);
    const auto *fh_114 = buffer.data(fh + 114);
    const auto *fh_115 = buffer.data(fh + 115);
    const auto *fh_116 = buffer.data(fh + 116);
    const auto *fh_117 = buffer.data(fh + 117);
    const auto *fh_118 = buffer.data(fh + 118);
    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_126 = buffer.data(fh + 126);
    const auto *fh_127 = buffer.data(fh + 127);
    const auto *fh_128 = buffer.data(fh + 128);
    const auto *fh_129 = buffer.data(fh + 129);
    const auto *fh_130 = buffer.data(fh + 130);
    const auto *fh_131 = buffer.data(fh + 131);
    const auto *fh_132 = buffer.data(fh + 132);
    const auto *fh_133 = buffer.data(fh + 133);
    const auto *fh_134 = buffer.data(fh + 134);
    const auto *fh_135 = buffer.data(fh + 135);
    const auto *fh_136 = buffer.data(fh + 136);
    const auto *fh_137 = buffer.data(fh + 137);
    const auto *fh_138 = buffer.data(fh + 138);
    const auto *fh_139 = buffer.data(fh + 139);
    const auto *fh_140 = buffer.data(fh + 140);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_142 = buffer.data(fh + 142);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_144 = buffer.data(fh + 144);
    const auto *fh_145 = buffer.data(fh + 145);
    const auto *fh_147 = buffer.data(fh + 147);
    const auto *fh_148 = buffer.data(fh + 148);
    const auto *fh_149 = buffer.data(fh + 149);
    const auto *fh_150 = buffer.data(fh + 150);
    const auto *fh_151 = buffer.data(fh + 151);
    const auto *fh_152 = buffer.data(fh + 152);
    const auto *fh_153 = buffer.data(fh + 153);
    const auto *fh_154 = buffer.data(fh + 154);
    const auto *fh_155 = buffer.data(fh + 155);
    const auto *fh_156 = buffer.data(fh + 156);
    const auto *fh_157 = buffer.data(fh + 157);
    const auto *fh_158 = buffer.data(fh + 158);
    const auto *fh_159 = buffer.data(fh + 159);
    const auto *fh_160 = buffer.data(fh + 160);
    const auto *fh_161 = buffer.data(fh + 161);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_163 = buffer.data(fh + 163);
    const auto *fh_164 = buffer.data(fh + 164);
    const auto *fh_165 = buffer.data(fh + 165);
    const auto *fh_166 = buffer.data(fh + 166);
    const auto *fh_168 = buffer.data(fh + 168);
    const auto *fh_169 = buffer.data(fh + 169);
    const auto *fh_170 = buffer.data(fh + 170);
    const auto *fh_171 = buffer.data(fh + 171);
    const auto *fh_172 = buffer.data(fh + 172);
    const auto *fh_173 = buffer.data(fh + 173);
    const auto *fh_174 = buffer.data(fh + 174);
    const auto *fh_175 = buffer.data(fh + 175);
    const auto *fh_176 = buffer.data(fh + 176);
    const auto *fh_177 = buffer.data(fh + 177);
    const auto *fh_178 = buffer.data(fh + 178);
    const auto *fh_179 = buffer.data(fh + 179);
    const auto *fh_180 = buffer.data(fh + 180);
    const auto *fh_181 = buffer.data(fh + 181);
    const auto *fh_182 = buffer.data(fh + 182);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_184 = buffer.data(fh + 184);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_189 = buffer.data(fh + 189);
    const auto *fh_190 = buffer.data(fh + 190);
    const auto *fh_191 = buffer.data(fh + 191);
    const auto *fh_192 = buffer.data(fh + 192);
    const auto *fh_193 = buffer.data(fh + 193);
    const auto *fh_194 = buffer.data(fh + 194);
    const auto *fh_195 = buffer.data(fh + 195);
    const auto *fh_196 = buffer.data(fh + 196);
    const auto *fh_197 = buffer.data(fh + 197);
    const auto *fh_198 = buffer.data(fh + 198);
    const auto *fh_199 = buffer.data(fh + 199);
    const auto *fh_200 = buffer.data(fh + 200);
    const auto *fh_201 = buffer.data(fh + 201);
    const auto *fh_202 = buffer.data(fh + 202);
    const auto *fh_203 = buffer.data(fh + 203);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_205 = buffer.data(fh + 205);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_208 = buffer.data(fh + 208);
    const auto *fh_209 = buffer.data(fh + 209);

#pragma omp simd aligned(ab_x, fg_16, fg_21, fg_91, fg_96, fh_22, fh_27, fh_127, \
                         fh_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -8.75 * ab_x[k] * fg_16[k]
                 + 8.75 * ab_x[k] * fg_21[k]
                 + 8.75 * ab_x[k] * fg_91[k]
                 - 8.75 * ab_x[k] * fg_96[k]
                 + 8.75 * fh_22[k]
                 - 8.75 * fh_27[k]
                 - 8.75 * fh_127[k]
                 + 8.75 * fh_132[k];
    }

#pragma omp simd aligned(ab_x, fg_19, fg_26, fg_94, fg_101, fh_25, fh_32, fh_130, \
                         fh_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = -f_0 * ab_x[k] * fg_19[k]
                 + f_1 * ab_x[k] * fg_26[k]
                 + f_0 * ab_x[k] * fg_94[k]
                 - f_1 * ab_x[k] * fg_101[k]
                 + f_0 * fh_25[k]
                 - f_1 * fh_32[k]
                 - f_0 * fh_130[k]
                 + f_1 * fh_137[k];
    }

#pragma omp simd aligned(ab_x, fg_16, fg_21, fg_23, fg_91, fg_96, fg_98, fh_22, fh_27, fh_29, \
                         fh_127, fh_132, fh_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = f_2 * ab_x[k] * fg_16[k]
                 + f_2 * ab_x[k] * fg_21[k]
                 - f_3 * ab_x[k] * fg_23[k]
                 - f_2 * ab_x[k] * fg_91[k]
                 - f_2 * ab_x[k] * fg_96[k]
                 + f_3 * ab_x[k] * fg_98[k]
                 - f_2 * fh_22[k]
                 - f_2 * fh_27[k]
                 + f_3 * fh_29[k]
                 + f_2 * fh_127[k]
                 + f_2 * fh_132[k]
                 - f_3 * fh_134[k];
    }

#pragma omp simd aligned(ab_x, fg_19, fg_26, fg_28, fg_94, fg_101, fg_103, fh_25, fh_32, \
                         fh_34, fh_130, fh_137, fh_139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_4 * ab_x[k] * fg_19[k]
                 + f_4 * ab_x[k] * fg_26[k]
                 - f_5 * ab_x[k] * fg_28[k]
                 - f_4 * ab_x[k] * fg_94[k]
                 - f_4 * ab_x[k] * fg_101[k]
                 + f_5 * ab_x[k] * fg_103[k]
                 - f_4 * fh_25[k]
                 - f_4 * fh_32[k]
                 + f_5 * fh_34[k]
                 + f_4 * fh_130[k]
                 + f_4 * fh_137[k]
                 - f_5 * fh_139[k];
    }

#pragma omp simd aligned(ab_x, fg_15, fg_18, fg_20, fg_25, fg_27, fg_29, fg_90, fg_93, fg_95, \
                         fg_100, fg_102, fg_104, fh_21, fh_24, fh_26, fh_31, fh_33, fh_35, \
                         fh_126, fh_129, fh_131, fh_136, fh_138, \
                         fh_140 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_6 * ab_x[k] * fg_15[k]
                 - f_7 * ab_x[k] * fg_18[k]
                 + f_8 * ab_x[k] * fg_20[k]
                 - f_6 * ab_x[k] * fg_25[k]
                 + f_8 * ab_x[k] * fg_27[k]
                 - f_9 * ab_x[k] * fg_29[k]
                 + f_6 * ab_x[k] * fg_90[k]
                 + f_7 * ab_x[k] * fg_93[k]
                 - f_8 * ab_x[k] * fg_95[k]
                 + f_6 * ab_x[k] * fg_100[k]
                 - f_8 * ab_x[k] * fg_102[k]
                 + f_9 * ab_x[k] * fg_104[k]
                 + f_6 * fh_21[k]
                 + f_7 * fh_24[k]
                 - f_8 * fh_26[k]
                 + f_6 * fh_31[k]
                 - f_8 * fh_33[k]
                 + f_9 * fh_35[k]
                 - f_6 * fh_126[k]
                 - f_7 * fh_129[k]
                 + f_8 * fh_131[k]
                 - f_6 * fh_136[k]
                 + f_8 * fh_138[k]
                 - f_9 * fh_140[k];
    }

#pragma omp simd aligned(ab_x, fg_17, fg_22, fg_24, fg_92, fg_97, fg_99, fh_23, fh_28, fh_30, \
                         fh_128, fh_133, fh_135 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_4 * ab_x[k] * fg_17[k]
                 + f_4 * ab_x[k] * fg_22[k]
                 - f_5 * ab_x[k] * fg_24[k]
                 - f_4 * ab_x[k] * fg_92[k]
                 - f_4 * ab_x[k] * fg_97[k]
                 + f_5 * ab_x[k] * fg_99[k]
                 - f_4 * fh_23[k]
                 - f_4 * fh_28[k]
                 + f_5 * fh_30[k]
                 + f_4 * fh_128[k]
                 + f_4 * fh_133[k]
                 - f_5 * fh_135[k];
    }

#pragma omp simd aligned(ab_x, fg_15, fg_20, fg_25, fg_27, fg_90, fg_95, fg_100, fg_102, \
                         fh_21, fh_26, fh_31, fh_33, fh_126, fh_131, fh_136, \
                         fh_138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_10 * ab_x[k] * fg_15[k]
                 - f_11 * ab_x[k] * fg_20[k]
                 - f_10 * ab_x[k] * fg_25[k]
                 + f_11 * ab_x[k] * fg_27[k]
                 - f_10 * ab_x[k] * fg_90[k]
                 + f_11 * ab_x[k] * fg_95[k]
                 + f_10 * ab_x[k] * fg_100[k]
                 - f_11 * ab_x[k] * fg_102[k]
                 - f_10 * fh_21[k]
                 + f_11 * fh_26[k]
                 + f_10 * fh_31[k]
                 - f_11 * fh_33[k]
                 + f_10 * fh_126[k]
                 - f_11 * fh_131[k]
                 - f_10 * fh_136[k]
                 + f_11 * fh_138[k];
    }

#pragma omp simd aligned(ab_x, fg_17, fg_22, fg_92, fg_97, fh_23, fh_28, fh_128, \
                         fh_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_1 * ab_x[k] * fg_17[k]
                 + f_0 * ab_x[k] * fg_22[k]
                 + f_1 * ab_x[k] * fg_92[k]
                 - f_0 * ab_x[k] * fg_97[k]
                 + f_1 * fh_23[k]
                 - f_0 * fh_28[k]
                 - f_1 * fh_128[k]
                 + f_0 * fh_133[k];
    }

#pragma omp simd aligned(ab_x, fg_15, fg_18, fg_25, fg_90, fg_93, fg_100, fh_21, fh_24, fh_31, \
                         fh_126, fh_129, fh_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -2.1875 * ab_x[k] * fg_15[k]
                 + 13.125 * ab_x[k] * fg_18[k]
                 - 2.1875 * ab_x[k] * fg_25[k]
                 + 2.1875 * ab_x[k] * fg_90[k]
                 - 13.125 * ab_x[k] * fg_93[k]
                 + 2.1875 * ab_x[k] * fg_100[k]
                 + 2.1875 * fh_21[k]
                 - 13.125 * fh_24[k]
                 + 2.1875 * fh_31[k]
                 - 2.1875 * fh_126[k]
                 + 13.125 * fh_129[k]
                 - 2.1875 * fh_136[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_61, fg_66, fg_106, fg_111, fh_85, fh_90, fh_150, \
                         fh_157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_0 * ab_x[k] * fg_61[k]
                 + f_0 * ab_x[k] * fg_66[k]
                 + f_1 * ab_y[k] * fg_106[k]
                 - f_1 * ab_y[k] * fg_111[k]
                 + f_0 * fh_85[k]
                 - f_0 * fh_90[k]
                 - f_1 * fh_150[k]
                 + f_1 * fh_157[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_64, fg_71, fg_109, fg_116, fh_88, fh_95, fh_154, \
                         fh_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -39.375 * ab_x[k] * fg_64[k]
                  + 13.125 * ab_x[k] * fg_71[k]
                  + 13.125 * ab_y[k] * fg_109[k]
                  - 4.375 * ab_y[k] * fg_116[k]
                  + 39.375 * fh_88[k]
                  - 13.125 * fh_95[k]
                  - 13.125 * fh_154[k]
                  + 4.375 * fh_163[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_61, fg_66, fg_68, fg_106, fg_111, fg_113, fh_85, \
                         fh_90, fh_92, fh_150, fh_157, fh_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_4 * ab_x[k] * fg_61[k]
                  + f_4 * ab_x[k] * fg_66[k]
                  - f_12 * ab_x[k] * fg_68[k]
                  - f_13 * ab_y[k] * fg_106[k]
                  - f_13 * ab_y[k] * fg_111[k]
                  + f_14 * ab_y[k] * fg_113[k]
                  - f_4 * fh_85[k]
                  - f_4 * fh_90[k]
                  + f_12 * fh_92[k]
                  + f_13 * fh_150[k]
                  + f_13 * fh_157[k]
                  - f_14 * fh_159[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_64, fg_71, fg_73, fg_109, fg_116, fg_118, fh_88, \
                         fh_95, fh_97, fh_154, fh_163, fh_165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_15 * ab_x[k] * fg_64[k]
                  + f_15 * ab_x[k] * fg_71[k]
                  - f_3 * ab_x[k] * fg_73[k]
                  - f_16 * ab_y[k] * fg_109[k]
                  - f_16 * ab_y[k] * fg_116[k]
                  + f_17 * ab_y[k] * fg_118[k]
                  - f_15 * fh_88[k]
                  - f_15 * fh_95[k]
                  + f_3 * fh_97[k]
                  + f_16 * fh_154[k]
                  + f_16 * fh_163[k]
                  - f_17 * fh_165[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_60, fg_63, fg_65, fg_70, fg_72, fg_74, fg_105, fg_108, \
                         fg_110, fg_115, fg_117, fg_119, fh_84, fh_87, fh_89, fh_94, fh_96, \
                         fh_98, fh_148, fh_153, fh_155, fh_162, fh_164, \
                         fh_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_18 * ab_x[k] * fg_60[k]
                  - f_19 * ab_x[k] * fg_63[k]
                  + f_20 * ab_x[k] * fg_65[k]
                  - f_18 * ab_x[k] * fg_70[k]
                  + f_20 * ab_x[k] * fg_72[k]
                  - f_21 * ab_x[k] * fg_74[k]
                  + f_22 * ab_y[k] * fg_105[k]
                  + f_23 * ab_y[k] * fg_108[k]
                  - f_21 * ab_y[k] * fg_110[k]
                  + f_22 * ab_y[k] * fg_115[k]
                  - f_21 * ab_y[k] * fg_117[k]
                  + f_24 * ab_y[k] * fg_119[k]
                  + f_18 * fh_84[k]
                  + f_19 * fh_87[k]
                  - f_20 * fh_89[k]
                  + f_18 * fh_94[k]
                  - f_20 * fh_96[k]
                  + f_21 * fh_98[k]
                  - f_22 * fh_148[k]
                  - f_23 * fh_153[k]
                  + f_21 * fh_155[k]
                  - f_22 * fh_162[k]
                  + f_21 * fh_164[k]
                  - f_24 * fh_166[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_62, fg_67, fg_69, fg_107, fg_112, fg_114, fh_86, \
                         fh_91, fh_93, fh_151, fh_158, fh_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_15 * ab_x[k] * fg_62[k]
                  + f_15 * ab_x[k] * fg_67[k]
                  - f_3 * ab_x[k] * fg_69[k]
                  - f_16 * ab_y[k] * fg_107[k]
                  - f_16 * ab_y[k] * fg_112[k]
                  + f_17 * ab_y[k] * fg_114[k]
                  - f_15 * fh_86[k]
                  - f_15 * fh_91[k]
                  + f_3 * fh_93[k]
                  + f_16 * fh_151[k]
                  + f_16 * fh_158[k]
                  - f_17 * fh_160[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_60, fg_65, fg_70, fg_72, fg_105, fg_110, fg_115, \
                         fg_117, fh_84, fh_89, fh_94, fh_96, fh_148, fh_155, fh_162, \
                         fh_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_25 * ab_x[k] * fg_60[k]
                  - f_26 * ab_x[k] * fg_65[k]
                  - f_25 * ab_x[k] * fg_70[k]
                  + f_26 * ab_x[k] * fg_72[k]
                  - f_27 * ab_y[k] * fg_105[k]
                  + f_4 * ab_y[k] * fg_110[k]
                  + f_27 * ab_y[k] * fg_115[k]
                  - f_4 * ab_y[k] * fg_117[k]
                  - f_25 * fh_84[k]
                  + f_26 * fh_89[k]
                  + f_25 * fh_94[k]
                  - f_26 * fh_96[k]
                  + f_27 * fh_148[k]
                  - f_4 * fh_155[k]
                  - f_27 * fh_162[k]
                  + f_4 * fh_164[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_62, fg_67, fg_107, fg_112, fh_86, fh_91, fh_151, \
                         fh_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -13.125 * ab_x[k] * fg_62[k]
                  + 39.375 * ab_x[k] * fg_67[k]
                  + 4.375 * ab_y[k] * fg_107[k]
                  - 13.125 * ab_y[k] * fg_112[k]
                  + 13.125 * fh_86[k]
                  - 39.375 * fh_91[k]
                  - 4.375 * fh_151[k]
                  + 13.125 * fh_158[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_60, fg_63, fg_70, fg_105, fg_108, fg_115, fh_84, \
                         fh_87, fh_94, fh_148, fh_153, fh_162 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_28 * ab_x[k] * fg_60[k]
                  + f_29 * ab_x[k] * fg_63[k]
                  - f_28 * ab_x[k] * fg_70[k]
                  + f_30 * ab_y[k] * fg_105[k]
                  - f_31 * ab_y[k] * fg_108[k]
                  + f_30 * ab_y[k] * fg_115[k]
                  + f_28 * fh_84[k]
                  - f_29 * fh_87[k]
                  + f_28 * fh_94[k]
                  - f_30 * fh_148[k]
                  + f_31 * fh_153[k]
                  - f_30 * fh_162[k];
    }

#pragma omp simd aligned(ab_x, fg_16, fg_21, fg_91, fg_96, fg_121, fg_126, fh_22, fh_27, \
                         fh_127, fh_132, fh_169, fh_174 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_2 * ab_x[k] * fg_16[k]
                  - f_2 * ab_x[k] * fg_21[k]
                  + f_2 * ab_x[k] * fg_91[k]
                  - f_2 * ab_x[k] * fg_96[k]
                  - f_3 * ab_x[k] * fg_121[k]
                  + f_3 * ab_x[k] * fg_126[k]
                  - f_2 * fh_22[k]
                  + f_2 * fh_27[k]
                  - f_2 * fh_127[k]
                  + f_2 * fh_132[k]
                  + f_3 * fh_169[k]
                  - f_3 * fh_174[k];
    }

#pragma omp simd aligned(ab_x, fg_19, fg_26, fg_94, fg_101, fg_124, fg_131, fh_25, fh_32, \
                         fh_130, fh_137, fh_172, fh_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_4 * ab_x[k] * fg_19[k]
                  - f_13 * ab_x[k] * fg_26[k]
                  + f_4 * ab_x[k] * fg_94[k]
                  - f_13 * ab_x[k] * fg_101[k]
                  - f_12 * ab_x[k] * fg_124[k]
                  + f_14 * ab_x[k] * fg_131[k]
                  - f_4 * fh_25[k]
                  + f_13 * fh_32[k]
                  - f_4 * fh_130[k]
                  + f_13 * fh_137[k]
                  + f_12 * fh_172[k]
                  - f_14 * fh_179[k];
    }

#pragma omp simd aligned(ab_x, fg_16, fg_21, fg_23, fg_91, fg_96, fg_98, fg_121, fg_126, \
                         fg_128, fh_22, fh_27, fh_29, fh_127, fh_132, fh_134, fh_169, fh_174, \
                         fh_176 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -1.25 * ab_x[k] * fg_16[k]
                  - 1.25 * ab_x[k] * fg_21[k]
                  + 7.5 * ab_x[k] * fg_23[k]
                  - 1.25 * ab_x[k] * fg_91[k]
                  - 1.25 * ab_x[k] * fg_96[k]
                  + 7.5 * ab_x[k] * fg_98[k]
                  + 7.5 * ab_x[k] * fg_121[k]
                  + 7.5 * ab_x[k] * fg_126[k]
                  - 45.0 * ab_x[k] * fg_128[k]
                  + 1.25 * fh_22[k]
                  + 1.25 * fh_27[k]
                  - 7.5 * fh_29[k]
                  + 1.25 * fh_127[k]
                  + 1.25 * fh_132[k]
                  - 7.5 * fh_134[k]
                  - 7.5 * fh_169[k]
                  - 7.5 * fh_174[k]
                  + 45.0 * fh_176[k];
    }

#pragma omp simd aligned(ab_x, fg_19, fg_26, fg_28, fg_94, fg_101, fg_103, fg_124, fg_131, \
                         fg_133, fh_25, fh_32, fh_34, fh_130, fh_137, fh_139, fh_172, fh_179, \
                         fh_181 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_32 * ab_x[k] * fg_19[k]
                  - f_32 * ab_x[k] * fg_26[k]
                  + f_33 * ab_x[k] * fg_28[k]
                  - f_32 * ab_x[k] * fg_94[k]
                  - f_32 * ab_x[k] * fg_101[k]
                  + f_33 * ab_x[k] * fg_103[k]
                  + f_34 * ab_x[k] * fg_124[k]
                  + f_34 * ab_x[k] * fg_131[k]
                  - f_35 * ab_x[k] * fg_133[k]
                  + f_32 * fh_25[k]
                  + f_32 * fh_32[k]
                  - f_33 * fh_34[k]
                  + f_32 * fh_130[k]
                  + f_32 * fh_137[k]
                  - f_33 * fh_139[k]
                  - f_34 * fh_172[k]
                  - f_34 * fh_179[k]
                  + f_35 * fh_181[k];
    }

#pragma omp simd aligned(ab_x, fg_15, fg_18, fg_20, fg_25, fg_27, fg_29, fg_90, fg_93, fg_95, \
                         fg_100, fg_102, fg_104, fg_120, fg_123, fg_125, fg_130, fg_132, \
                         fg_134, fh_21, fh_24, fh_26, fh_31, fh_33, fh_35, fh_126, fh_129, \
                         fh_131, fh_136, fh_138, fh_140, fh_168, fh_171, fh_173, fh_178, \
                         fh_180, fh_182 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_36 * ab_x[k] * fg_15[k]
                  + f_37 * ab_x[k] * fg_18[k]
                  - f_38 * ab_x[k] * fg_20[k]
                  + f_36 * ab_x[k] * fg_25[k]
                  - f_38 * ab_x[k] * fg_27[k]
                  + f_39 * ab_x[k] * fg_29[k]
                  + f_36 * ab_x[k] * fg_90[k]
                  + f_37 * ab_x[k] * fg_93[k]
                  - f_38 * ab_x[k] * fg_95[k]
                  + f_36 * ab_x[k] * fg_100[k]
                  - f_38 * ab_x[k] * fg_102[k]
                  + f_39 * ab_x[k] * fg_104[k]
                  - f_40 * ab_x[k] * fg_120[k]
                  - f_41 * ab_x[k] * fg_123[k]
                  + f_42 * ab_x[k] * fg_125[k]
                  - f_40 * ab_x[k] * fg_130[k]
                  + f_42 * ab_x[k] * fg_132[k]
                  - f_43 * ab_x[k] * fg_134[k]
                  - f_36 * fh_21[k]
                  - f_37 * fh_24[k]
                  + f_38 * fh_26[k]
                  - f_36 * fh_31[k]
                  + f_38 * fh_33[k]
                  - f_39 * fh_35[k]
                  - f_36 * fh_126[k]
                  - f_37 * fh_129[k]
                  + f_38 * fh_131[k]
                  - f_36 * fh_136[k]
                  + f_38 * fh_138[k]
                  - f_39 * fh_140[k]
                  + f_40 * fh_168[k]
                  + f_41 * fh_171[k]
                  - f_42 * fh_173[k]
                  + f_40 * fh_178[k]
                  - f_42 * fh_180[k]
                  + f_43 * fh_182[k];
    }

#pragma omp simd aligned(ab_x, fg_17, fg_22, fg_24, fg_92, fg_97, fg_99, fg_122, fg_127, \
                         fg_129, fh_23, fh_28, fh_30, fh_128, fh_133, fh_135, fh_170, fh_175, \
                         fh_177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_32 * ab_x[k] * fg_17[k]
                  - f_32 * ab_x[k] * fg_22[k]
                  + f_33 * ab_x[k] * fg_24[k]
                  - f_32 * ab_x[k] * fg_92[k]
                  - f_32 * ab_x[k] * fg_97[k]
                  + f_33 * ab_x[k] * fg_99[k]
                  + f_34 * ab_x[k] * fg_122[k]
                  + f_34 * ab_x[k] * fg_127[k]
                  - f_35 * ab_x[k] * fg_129[k]
                  + f_32 * fh_23[k]
                  + f_32 * fh_28[k]
                  - f_33 * fh_30[k]
                  + f_32 * fh_128[k]
                  + f_32 * fh_133[k]
                  - f_33 * fh_135[k]
                  - f_34 * fh_170[k]
                  - f_34 * fh_175[k]
                  + f_35 * fh_177[k];
    }

#pragma omp simd aligned(ab_x, fg_15, fg_20, fg_25, fg_27, fg_90, fg_95, fg_100, fg_102, \
                         fg_120, fg_125, fg_130, fg_132, fh_21, fh_26, fh_31, fh_33, fh_126, \
                         fh_131, fh_136, fh_138, fh_168, fh_173, fh_178, \
                         fh_180 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -0.625 * ab_x[k] * fg_15[k]
                  + 3.75 * ab_x[k] * fg_20[k]
                  + 0.625 * ab_x[k] * fg_25[k]
                  - 3.75 * ab_x[k] * fg_27[k]
                  - 0.625 * ab_x[k] * fg_90[k]
                  + 3.75 * ab_x[k] * fg_95[k]
                  + 0.625 * ab_x[k] * fg_100[k]
                  - 3.75 * ab_x[k] * fg_102[k]
                  + 3.75 * ab_x[k] * fg_120[k]
                  - 22.5 * ab_x[k] * fg_125[k]
                  - 3.75 * ab_x[k] * fg_130[k]
                  + 22.5 * ab_x[k] * fg_132[k]
                  + 0.625 * fh_21[k]
                  - 3.75 * fh_26[k]
                  - 0.625 * fh_31[k]
                  + 3.75 * fh_33[k]
                  + 0.625 * fh_126[k]
                  - 3.75 * fh_131[k]
                  - 0.625 * fh_136[k]
                  + 3.75 * fh_138[k]
                  - 3.75 * fh_168[k]
                  + 22.5 * fh_173[k]
                  + 3.75 * fh_178[k]
                  - 22.5 * fh_180[k];
    }

#pragma omp simd aligned(ab_x, fg_17, fg_22, fg_92, fg_97, fg_122, fg_127, fh_23, fh_28, \
                         fh_128, fh_133, fh_170, fh_175 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_13 * ab_x[k] * fg_17[k]
                  - f_4 * ab_x[k] * fg_22[k]
                  + f_13 * ab_x[k] * fg_92[k]
                  - f_4 * ab_x[k] * fg_97[k]
                  - f_14 * ab_x[k] * fg_122[k]
                  + f_12 * ab_x[k] * fg_127[k]
                  - f_13 * fh_23[k]
                  + f_4 * fh_28[k]
                  - f_13 * fh_128[k]
                  + f_4 * fh_133[k]
                  + f_14 * fh_170[k]
                  - f_12 * fh_175[k];
    }

#pragma omp simd aligned(ab_x, fg_15, fg_18, fg_25, fg_90, fg_93, fg_100, fg_120, fg_123, \
                         fg_130, fh_21, fh_24, fh_31, fh_126, fh_129, fh_136, fh_168, fh_171, \
                         fh_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_44 * ab_x[k] * fg_15[k]
                  - f_16 * ab_x[k] * fg_18[k]
                  + f_44 * ab_x[k] * fg_25[k]
                  + f_44 * ab_x[k] * fg_90[k]
                  - f_16 * ab_x[k] * fg_93[k]
                  + f_44 * ab_x[k] * fg_100[k]
                  - f_16 * ab_x[k] * fg_120[k]
                  + f_45 * ab_x[k] * fg_123[k]
                  - f_16 * ab_x[k] * fg_130[k]
                  - f_44 * fh_21[k]
                  + f_16 * fh_24[k]
                  - f_44 * fh_31[k]
                  - f_44 * fh_126[k]
                  + f_16 * fh_129[k]
                  - f_44 * fh_136[k]
                  + f_16 * fh_168[k]
                  - f_45 * fh_171[k]
                  + f_16 * fh_178[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_61, fg_66, fg_106, fg_111, fg_136, fg_141, fh_85, \
                         fh_90, fh_150, fh_157, fh_192, fh_199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = f_4 * ab_x[k] * fg_61[k]
                  - f_4 * ab_x[k] * fg_66[k]
                  + f_4 * ab_y[k] * fg_106[k]
                  - f_4 * ab_y[k] * fg_111[k]
                  - f_5 * ab_y[k] * fg_136[k]
                  + f_5 * ab_y[k] * fg_141[k]
                  - f_4 * fh_85[k]
                  + f_4 * fh_90[k]
                  - f_4 * fh_150[k]
                  + f_4 * fh_157[k]
                  + f_5 * fh_192[k]
                  - f_5 * fh_199[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_64, fg_71, fg_109, fg_116, fg_139, fg_146, fh_88, \
                         fh_95, fh_154, fh_163, fh_196, fh_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_15 * ab_x[k] * fg_64[k]
                  - f_16 * ab_x[k] * fg_71[k]
                  + f_15 * ab_y[k] * fg_109[k]
                  - f_16 * ab_y[k] * fg_116[k]
                  - f_3 * ab_y[k] * fg_139[k]
                  + f_17 * ab_y[k] * fg_146[k]
                  - f_15 * fh_88[k]
                  + f_16 * fh_95[k]
                  - f_15 * fh_154[k]
                  + f_16 * fh_163[k]
                  + f_3 * fh_196[k]
                  - f_17 * fh_205[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_61, fg_66, fg_68, fg_106, fg_111, fg_113, fg_136, \
                         fg_141, fg_143, fh_85, fh_90, fh_92, fh_150, fh_157, fh_159, fh_192, \
                         fh_199, fh_201 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = -f_32 * ab_x[k] * fg_61[k]
                  - f_32 * ab_x[k] * fg_66[k]
                  + f_34 * ab_x[k] * fg_68[k]
                  - f_32 * ab_y[k] * fg_106[k]
                  - f_32 * ab_y[k] * fg_111[k]
                  + f_34 * ab_y[k] * fg_113[k]
                  + f_33 * ab_y[k] * fg_136[k]
                  + f_33 * ab_y[k] * fg_141[k]
                  - f_35 * ab_y[k] * fg_143[k]
                  + f_32 * fh_85[k]
                  + f_32 * fh_90[k]
                  - f_34 * fh_92[k]
                  + f_32 * fh_150[k]
                  + f_32 * fh_157[k]
                  - f_34 * fh_159[k]
                  - f_33 * fh_192[k]
                  - f_33 * fh_199[k]
                  + f_35 * fh_201[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_64, fg_71, fg_73, fg_109, fg_116, fg_118, fg_139, \
                         fg_146, fg_148, fh_88, fh_95, fh_97, fh_154, fh_163, fh_165, fh_196, \
                         fh_205, fh_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -5.625 * ab_x[k] * fg_64[k]
                  - 5.625 * ab_x[k] * fg_71[k]
                  + 7.5 * ab_x[k] * fg_73[k]
                  - 5.625 * ab_y[k] * fg_109[k]
                  - 5.625 * ab_y[k] * fg_116[k]
                  + 7.5 * ab_y[k] * fg_118[k]
                  + 7.5 * ab_y[k] * fg_139[k]
                  + 7.5 * ab_y[k] * fg_146[k]
                  - 10.0 * ab_y[k] * fg_148[k]
                  + 5.625 * fh_88[k]
                  + 5.625 * fh_95[k]
                  - 7.5 * fh_97[k]
                  + 5.625 * fh_154[k]
                  + 5.625 * fh_163[k]
                  - 7.5 * fh_165[k]
                  - 7.5 * fh_196[k]
                  - 7.5 * fh_205[k]
                  + 10.0 * fh_207[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_60, fg_63, fg_65, fg_70, fg_72, fg_74, fg_105, fg_108, \
                         fg_110, fg_115, fg_117, fg_119, fg_135, fg_138, fg_140, fg_145, \
                         fg_147, fg_149, fh_84, fh_87, fh_89, fh_94, fh_96, fh_98, fh_148, \
                         fh_153, fh_155, fh_162, fh_164, fh_166, fh_190, fh_195, fh_197, \
                         fh_204, fh_206, fh_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = f_46 * ab_x[k] * fg_60[k]
                  + f_47 * ab_x[k] * fg_63[k]
                  - f_48 * ab_x[k] * fg_65[k]
                  + f_46 * ab_x[k] * fg_70[k]
                  - f_48 * ab_x[k] * fg_72[k]
                  + f_49 * ab_x[k] * fg_74[k]
                  + f_46 * ab_y[k] * fg_105[k]
                  + f_47 * ab_y[k] * fg_108[k]
                  - f_48 * ab_y[k] * fg_110[k]
                  + f_46 * ab_y[k] * fg_115[k]
                  - f_48 * ab_y[k] * fg_117[k]
                  + f_49 * ab_y[k] * fg_119[k]
                  - f_50 * ab_y[k] * fg_135[k]
                  - f_49 * ab_y[k] * fg_138[k]
                  + f_51 * ab_y[k] * fg_140[k]
                  - f_50 * ab_y[k] * fg_145[k]
                  + f_51 * ab_y[k] * fg_147[k]
                  - f_52 * ab_y[k] * fg_149[k]
                  - f_46 * fh_84[k]
                  - f_47 * fh_87[k]
                  + f_48 * fh_89[k]
                  - f_46 * fh_94[k]
                  + f_48 * fh_96[k]
                  - f_49 * fh_98[k]
                  - f_46 * fh_148[k]
                  - f_47 * fh_153[k]
                  + f_48 * fh_155[k]
                  - f_46 * fh_162[k]
                  + f_48 * fh_164[k]
                  - f_49 * fh_166[k]
                  + f_50 * fh_190[k]
                  + f_49 * fh_195[k]
                  - f_51 * fh_197[k]
                  + f_50 * fh_204[k]
                  - f_51 * fh_206[k]
                  + f_52 * fh_208[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_62, fg_67, fg_69, fg_107, fg_112, fg_114, fg_137, \
                         fg_142, fg_144, fh_86, fh_91, fh_93, fh_151, fh_158, fh_160, fh_193, \
                         fh_200, fh_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -5.625 * ab_x[k] * fg_62[k]
                  - 5.625 * ab_x[k] * fg_67[k]
                  + 7.5 * ab_x[k] * fg_69[k]
                  - 5.625 * ab_y[k] * fg_107[k]
                  - 5.625 * ab_y[k] * fg_112[k]
                  + 7.5 * ab_y[k] * fg_114[k]
                  + 7.5 * ab_y[k] * fg_137[k]
                  + 7.5 * ab_y[k] * fg_142[k]
                  - 10.0 * ab_y[k] * fg_144[k]
                  + 5.625 * fh_86[k]
                  + 5.625 * fh_91[k]
                  - 7.5 * fh_93[k]
                  + 5.625 * fh_151[k]
                  + 5.625 * fh_158[k]
                  - 7.5 * fh_160[k]
                  - 7.5 * fh_193[k]
                  - 7.5 * fh_200[k]
                  + 10.0 * fh_202[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_60, fg_65, fg_70, fg_72, fg_105, fg_110, fg_115, \
                         fg_117, fg_135, fg_140, fg_145, fg_147, fh_84, fh_89, fh_94, fh_96, \
                         fh_148, fh_155, fh_162, fh_164, fh_190, fh_197, fh_204, \
                         fh_206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_53 * ab_x[k] * fg_60[k]
                  + f_54 * ab_x[k] * fg_65[k]
                  + f_53 * ab_x[k] * fg_70[k]
                  - f_54 * ab_x[k] * fg_72[k]
                  - f_53 * ab_y[k] * fg_105[k]
                  + f_54 * ab_y[k] * fg_110[k]
                  + f_53 * ab_y[k] * fg_115[k]
                  - f_54 * ab_y[k] * fg_117[k]
                  + f_55 * ab_y[k] * fg_135[k]
                  - f_56 * ab_y[k] * fg_140[k]
                  - f_55 * ab_y[k] * fg_145[k]
                  + f_56 * ab_y[k] * fg_147[k]
                  + f_53 * fh_84[k]
                  - f_54 * fh_89[k]
                  - f_53 * fh_94[k]
                  + f_54 * fh_96[k]
                  + f_53 * fh_148[k]
                  - f_54 * fh_155[k]
                  - f_53 * fh_162[k]
                  + f_54 * fh_164[k]
                  - f_55 * fh_190[k]
                  + f_56 * fh_197[k]
                  + f_55 * fh_204[k]
                  - f_56 * fh_206[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_62, fg_67, fg_107, fg_112, fg_137, fg_142, fh_86, \
                         fh_91, fh_151, fh_158, fh_193, fh_200 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_16 * ab_x[k] * fg_62[k]
                  - f_15 * ab_x[k] * fg_67[k]
                  + f_16 * ab_y[k] * fg_107[k]
                  - f_15 * ab_y[k] * fg_112[k]
                  - f_17 * ab_y[k] * fg_137[k]
                  + f_3 * ab_y[k] * fg_142[k]
                  - f_16 * fh_86[k]
                  + f_15 * fh_91[k]
                  - f_16 * fh_151[k]
                  + f_15 * fh_158[k]
                  + f_17 * fh_193[k]
                  - f_3 * fh_200[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_60, fg_63, fg_70, fg_105, fg_108, fg_115, fg_135, \
                         fg_138, fg_145, fh_84, fh_87, fh_94, fh_148, fh_153, fh_162, fh_190, \
                         fh_195, fh_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_57 * ab_x[k] * fg_60[k]
                  - f_58 * ab_x[k] * fg_63[k]
                  + f_57 * ab_x[k] * fg_70[k]
                  + f_57 * ab_y[k] * fg_105[k]
                  - f_58 * ab_y[k] * fg_108[k]
                  + f_57 * ab_y[k] * fg_115[k]
                  - f_13 * ab_y[k] * fg_135[k]
                  + f_14 * ab_y[k] * fg_138[k]
                  - f_13 * ab_y[k] * fg_145[k]
                  - f_57 * fh_84[k]
                  + f_58 * fh_87[k]
                  - f_57 * fh_94[k]
                  - f_57 * fh_148[k]
                  + f_58 * fh_153[k]
                  - f_57 * fh_162[k]
                  + f_13 * fh_190[k]
                  - f_14 * fh_195[k]
                  + f_13 * fh_204[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_1, fg_6, fg_46, fg_51, fg_76, fg_81, fg_91, \
                         fg_96, fg_121, fg_126, fg_136, fg_141, fh_1, fh_6, fh_64, fh_69, \
                         fh_106, fh_111, fh_129, fh_136, fh_171, fh_178, fh_193, \
                         fh_200 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = -f_6 * ab_x[k] * fg_1[k]
                  + f_6 * ab_x[k] * fg_6[k]
                  - f_7 * ab_x[k] * fg_46[k]
                  + f_7 * ab_x[k] * fg_51[k]
                  + f_8 * ab_x[k] * fg_76[k]
                  - f_8 * ab_x[k] * fg_81[k]
                  - f_6 * ab_y[k] * fg_91[k]
                  + f_6 * ab_y[k] * fg_96[k]
                  + f_8 * ab_y[k] * fg_121[k]
                  - f_8 * ab_y[k] * fg_126[k]
                  - f_9 * ab_z[k] * fg_136[k]
                  + f_9 * ab_z[k] * fg_141[k]
                  + f_6 * fh_1[k]
                  - f_6 * fh_6[k]
                  + f_7 * fh_64[k]
                  - f_7 * fh_69[k]
                  - f_8 * fh_106[k]
                  + f_8 * fh_111[k]
                  + f_6 * fh_129[k]
                  - f_6 * fh_136[k]
                  - f_8 * fh_171[k]
                  + f_8 * fh_178[k]
                  + f_9 * fh_193[k]
                  - f_9 * fh_200[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_4, fg_11, fg_49, fg_56, fg_79, fg_86, fg_94, \
                         fg_101, fg_124, fg_131, fg_139, fg_146, fh_4, fh_11, fh_67, fh_74, \
                         fh_109, fh_116, fh_133, fh_142, fh_175, fh_184, fh_197, \
                         fh_206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_18 * ab_x[k] * fg_4[k]
                  + f_22 * ab_x[k] * fg_11[k]
                  - f_19 * ab_x[k] * fg_49[k]
                  + f_23 * ab_x[k] * fg_56[k]
                  + f_20 * ab_x[k] * fg_79[k]
                  - f_21 * ab_x[k] * fg_86[k]
                  - f_18 * ab_y[k] * fg_94[k]
                  + f_22 * ab_y[k] * fg_101[k]
                  + f_20 * ab_y[k] * fg_124[k]
                  - f_21 * ab_y[k] * fg_131[k]
                  - f_21 * ab_z[k] * fg_139[k]
                  + f_24 * ab_z[k] * fg_146[k]
                  + f_18 * fh_4[k]
                  - f_22 * fh_11[k]
                  + f_19 * fh_67[k]
                  - f_23 * fh_74[k]
                  - f_20 * fh_109[k]
                  + f_21 * fh_116[k]
                  + f_18 * fh_133[k]
                  - f_22 * fh_142[k]
                  - f_20 * fh_175[k]
                  + f_21 * fh_184[k]
                  + f_21 * fh_197[k]
                  - f_24 * fh_206[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_1, fg_6, fg_8, fg_46, fg_51, fg_53, fg_76, \
                         fg_81, fg_83, fg_91, fg_96, fg_98, fg_121, fg_126, fg_128, fg_136, \
                         fg_141, fg_143, fh_1, fh_6, fh_8, fh_64, fh_69, fh_71, fh_106, \
                         fh_111, fh_113, fh_129, fh_136, fh_138, fh_171, fh_178, fh_180, \
                         fh_193, fh_200, fh_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_36 * ab_x[k] * fg_1[k]
                  + f_36 * ab_x[k] * fg_6[k]
                  - f_40 * ab_x[k] * fg_8[k]
                  + f_37 * ab_x[k] * fg_46[k]
                  + f_37 * ab_x[k] * fg_51[k]
                  - f_41 * ab_x[k] * fg_53[k]
                  - f_38 * ab_x[k] * fg_76[k]
                  - f_38 * ab_x[k] * fg_81[k]
                  + f_42 * ab_x[k] * fg_83[k]
                  + f_36 * ab_y[k] * fg_91[k]
                  + f_36 * ab_y[k] * fg_96[k]
                  - f_40 * ab_y[k] * fg_98[k]
                  - f_38 * ab_y[k] * fg_121[k]
                  - f_38 * ab_y[k] * fg_126[k]
                  + f_42 * ab_y[k] * fg_128[k]
                  + f_39 * ab_z[k] * fg_136[k]
                  + f_39 * ab_z[k] * fg_141[k]
                  - f_43 * ab_z[k] * fg_143[k]
                  - f_36 * fh_1[k]
                  - f_36 * fh_6[k]
                  + f_40 * fh_8[k]
                  - f_37 * fh_64[k]
                  - f_37 * fh_69[k]
                  + f_41 * fh_71[k]
                  + f_38 * fh_106[k]
                  + f_38 * fh_111[k]
                  - f_42 * fh_113[k]
                  - f_36 * fh_129[k]
                  - f_36 * fh_136[k]
                  + f_40 * fh_138[k]
                  + f_38 * fh_171[k]
                  + f_38 * fh_178[k]
                  - f_42 * fh_180[k]
                  - f_39 * fh_193[k]
                  - f_39 * fh_200[k]
                  + f_43 * fh_202[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_4, fg_11, fg_13, fg_49, fg_56, fg_58, fg_79, \
                         fg_86, fg_88, fg_94, fg_101, fg_103, fg_124, fg_131, fg_133, fg_139, \
                         fg_146, fg_148, fh_4, fh_11, fh_13, fh_67, fh_74, fh_76, fh_109, \
                         fh_116, fh_118, fh_133, fh_142, fh_144, fh_175, fh_184, fh_186, \
                         fh_197, fh_206, fh_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_46 * ab_x[k] * fg_4[k]
                  + f_46 * ab_x[k] * fg_11[k]
                  - f_50 * ab_x[k] * fg_13[k]
                  + f_47 * ab_x[k] * fg_49[k]
                  + f_47 * ab_x[k] * fg_56[k]
                  - f_49 * ab_x[k] * fg_58[k]
                  - f_48 * ab_x[k] * fg_79[k]
                  - f_48 * ab_x[k] * fg_86[k]
                  + f_51 * ab_x[k] * fg_88[k]
                  + f_46 * ab_y[k] * fg_94[k]
                  + f_46 * ab_y[k] * fg_101[k]
                  - f_50 * ab_y[k] * fg_103[k]
                  - f_48 * ab_y[k] * fg_124[k]
                  - f_48 * ab_y[k] * fg_131[k]
                  + f_51 * ab_y[k] * fg_133[k]
                  + f_49 * ab_z[k] * fg_139[k]
                  + f_49 * ab_z[k] * fg_146[k]
                  - f_52 * ab_z[k] * fg_148[k]
                  - f_46 * fh_4[k]
                  - f_46 * fh_11[k]
                  + f_50 * fh_13[k]
                  - f_47 * fh_67[k]
                  - f_47 * fh_74[k]
                  + f_49 * fh_76[k]
                  + f_48 * fh_109[k]
                  + f_48 * fh_116[k]
                  - f_51 * fh_118[k]
                  - f_46 * fh_133[k]
                  - f_46 * fh_142[k]
                  + f_50 * fh_144[k]
                  + f_48 * fh_175[k]
                  + f_48 * fh_184[k]
                  - f_51 * fh_186[k]
                  - f_49 * fh_197[k]
                  - f_49 * fh_206[k]
                  + f_52 * fh_208[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_0, fg_3, fg_5, fg_10, fg_12, fg_14, fg_45, \
                         fg_48, fg_50, fg_55, fg_57, fg_59, fg_75, fg_78, fg_80, fg_85, fg_87, \
                         fg_89, fg_90, fg_93, fg_95, fg_100, fg_102, fg_104, fg_120, fg_123, \
                         fg_125, fg_130, fg_132, fg_134, fg_135, fg_138, fg_140, fg_145, \
                         fg_147, fg_149, fh_0, fh_3, fh_5, fh_10, fh_12, fh_14, fh_63, fh_66, \
                         fh_68, fh_73, fh_75, fh_77, fh_105, fh_108, fh_110, fh_115, fh_117, \
                         fh_119, fh_127, fh_132, fh_134, fh_141, fh_143, fh_145, fh_169, \
                         fh_174, fh_176, fh_183, fh_185, fh_187, fh_191, fh_196, fh_198, \
                         fh_205, fh_207, fh_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -0.140625 * ab_x[k] * fg_0[k]
                  - 0.28125 * ab_x[k] * fg_3[k]
                  + 1.125 * ab_x[k] * fg_5[k]
                  - 0.140625 * ab_x[k] * fg_10[k]
                  + 1.125 * ab_x[k] * fg_12[k]
                  - 0.375 * ab_x[k] * fg_14[k]
                  - 0.28125 * ab_x[k] * fg_45[k]
                  - 0.5625 * ab_x[k] * fg_48[k]
                  + 2.25 * ab_x[k] * fg_50[k]
                  - 0.28125 * ab_x[k] * fg_55[k]
                  + 2.25 * ab_x[k] * fg_57[k]
                  - 0.75 * ab_x[k] * fg_59[k]
                  + 1.125 * ab_x[k] * fg_75[k]
                  + 2.25 * ab_x[k] * fg_78[k]
                  - 9.0 * ab_x[k] * fg_80[k]
                  + 1.125 * ab_x[k] * fg_85[k]
                  - 9.0 * ab_x[k] * fg_87[k]
                  + 3.0 * ab_x[k] * fg_89[k]
                  - 0.140625 * ab_y[k] * fg_90[k]
                  - 0.28125 * ab_y[k] * fg_93[k]
                  + 1.125 * ab_y[k] * fg_95[k]
                  - 0.140625 * ab_y[k] * fg_100[k]
                  + 1.125 * ab_y[k] * fg_102[k]
                  - 0.375 * ab_y[k] * fg_104[k]
                  + 1.125 * ab_y[k] * fg_120[k]
                  + 2.25 * ab_y[k] * fg_123[k]
                  - 9.0 * ab_y[k] * fg_125[k]
                  + 1.125 * ab_y[k] * fg_130[k]
                  - 9.0 * ab_y[k] * fg_132[k]
                  + 3.0 * ab_y[k] * fg_134[k]
                  - 0.375 * ab_z[k] * fg_135[k]
                  - 0.75 * ab_z[k] * fg_138[k]
                  + 3.0 * ab_z[k] * fg_140[k]
                  - 0.375 * ab_z[k] * fg_145[k]
                  + 3.0 * ab_z[k] * fg_147[k]
                  - ab_z[k] * fg_149[k]
                  + 0.140625 * fh_0[k]
                  + 0.28125 * fh_3[k]
                  - 1.125 * fh_5[k]
                  + 0.140625 * fh_10[k]
                  - 1.125 * fh_12[k]
                  + 0.375 * fh_14[k]
                  + 0.28125 * fh_63[k]
                  + 0.5625 * fh_66[k]
                  - 2.25 * fh_68[k]
                  + 0.28125 * fh_73[k]
                  - 2.25 * fh_75[k]
                  + 0.75 * fh_77[k]
                  - 1.125 * fh_105[k]
                  - 2.25 * fh_108[k]
                  + 9.0 * fh_110[k]
                  - 1.125 * fh_115[k]
                  + 9.0 * fh_117[k]
                  - 3.0 * fh_119[k]
                  + 0.140625 * fh_127[k]
                  + 0.28125 * fh_132[k]
                  - 1.125 * fh_134[k]
                  + 0.140625 * fh_141[k]
                  - 1.125 * fh_143[k]
                  + 0.375 * fh_145[k]
                  - 1.125 * fh_169[k]
                  - 2.25 * fh_174[k]
                  + 9.0 * fh_176[k]
                  - 1.125 * fh_183[k]
                  + 9.0 * fh_185[k]
                  - 3.0 * fh_187[k]
                  + 0.375 * fh_191[k]
                  + 0.75 * fh_196[k]
                  - 3.0 * fh_198[k]
                  + 0.375 * fh_205[k]
                  - 3.0 * fh_207[k]
                  + fh_209[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_2, fg_7, fg_9, fg_47, fg_52, fg_54, fg_77, \
                         fg_82, fg_84, fg_92, fg_97, fg_99, fg_122, fg_127, fg_129, fg_137, \
                         fg_142, fg_144, fh_2, fh_7, fh_9, fh_65, fh_70, fh_72, fh_107, \
                         fh_112, fh_114, fh_130, fh_137, fh_139, fh_172, fh_179, fh_181, \
                         fh_194, fh_201, fh_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_46 * ab_x[k] * fg_2[k]
                  + f_46 * ab_x[k] * fg_7[k]
                  - f_50 * ab_x[k] * fg_9[k]
                  + f_47 * ab_x[k] * fg_47[k]
                  + f_47 * ab_x[k] * fg_52[k]
                  - f_49 * ab_x[k] * fg_54[k]
                  - f_48 * ab_x[k] * fg_77[k]
                  - f_48 * ab_x[k] * fg_82[k]
                  + f_51 * ab_x[k] * fg_84[k]
                  + f_46 * ab_y[k] * fg_92[k]
                  + f_46 * ab_y[k] * fg_97[k]
                  - f_50 * ab_y[k] * fg_99[k]
                  - f_48 * ab_y[k] * fg_122[k]
                  - f_48 * ab_y[k] * fg_127[k]
                  + f_51 * ab_y[k] * fg_129[k]
                  + f_49 * ab_z[k] * fg_137[k]
                  + f_49 * ab_z[k] * fg_142[k]
                  - f_52 * ab_z[k] * fg_144[k]
                  - f_46 * fh_2[k]
                  - f_46 * fh_7[k]
                  + f_50 * fh_9[k]
                  - f_47 * fh_65[k]
                  - f_47 * fh_70[k]
                  + f_49 * fh_72[k]
                  + f_48 * fh_107[k]
                  + f_48 * fh_112[k]
                  - f_51 * fh_114[k]
                  - f_46 * fh_130[k]
                  - f_46 * fh_137[k]
                  + f_50 * fh_139[k]
                  + f_48 * fh_172[k]
                  + f_48 * fh_179[k]
                  - f_51 * fh_181[k]
                  - f_49 * fh_194[k]
                  - f_49 * fh_201[k]
                  + f_52 * fh_203[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_0, fg_5, fg_10, fg_12, fg_45, fg_50, fg_55, \
                         fg_57, fg_75, fg_80, fg_85, fg_87, fg_90, fg_95, fg_100, fg_102, \
                         fg_120, fg_125, fg_130, fg_132, fg_135, fg_140, fg_145, fg_147, fh_0, \
                         fh_5, fh_10, fh_12, fh_63, fh_68, fh_73, fh_75, fh_105, fh_110, \
                         fh_115, fh_117, fh_127, fh_134, fh_141, fh_143, fh_169, fh_176, \
                         fh_183, fh_185, fh_191, fh_198, fh_205, \
                         fh_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_59 * ab_x[k] * fg_0[k]
                  - f_60 * ab_x[k] * fg_5[k]
                  - f_59 * ab_x[k] * fg_10[k]
                  + f_60 * ab_x[k] * fg_12[k]
                  + f_36 * ab_x[k] * fg_45[k]
                  - f_40 * ab_x[k] * fg_50[k]
                  - f_36 * ab_x[k] * fg_55[k]
                  + f_40 * ab_x[k] * fg_57[k]
                  - f_61 * ab_x[k] * fg_75[k]
                  + f_62 * ab_x[k] * fg_80[k]
                  + f_61 * ab_x[k] * fg_85[k]
                  - f_62 * ab_x[k] * fg_87[k]
                  + f_59 * ab_y[k] * fg_90[k]
                  - f_60 * ab_y[k] * fg_95[k]
                  - f_59 * ab_y[k] * fg_100[k]
                  + f_60 * ab_y[k] * fg_102[k]
                  - f_61 * ab_y[k] * fg_120[k]
                  + f_62 * ab_y[k] * fg_125[k]
                  + f_61 * ab_y[k] * fg_130[k]
                  - f_62 * ab_y[k] * fg_132[k]
                  + f_63 * ab_z[k] * fg_135[k]
                  - f_38 * ab_z[k] * fg_140[k]
                  - f_63 * ab_z[k] * fg_145[k]
                  + f_38 * ab_z[k] * fg_147[k]
                  - f_59 * fh_0[k]
                  + f_60 * fh_5[k]
                  + f_59 * fh_10[k]
                  - f_60 * fh_12[k]
                  - f_36 * fh_63[k]
                  + f_40 * fh_68[k]
                  + f_36 * fh_73[k]
                  - f_40 * fh_75[k]
                  + f_61 * fh_105[k]
                  - f_62 * fh_110[k]
                  - f_61 * fh_115[k]
                  + f_62 * fh_117[k]
                  - f_59 * fh_127[k]
                  + f_60 * fh_134[k]
                  + f_59 * fh_141[k]
                  - f_60 * fh_143[k]
                  + f_61 * fh_169[k]
                  - f_62 * fh_176[k]
                  - f_61 * fh_183[k]
                  + f_62 * fh_185[k]
                  - f_63 * fh_191[k]
                  + f_38 * fh_198[k]
                  + f_63 * fh_205[k]
                  - f_38 * fh_207[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_2, fg_7, fg_47, fg_52, fg_77, fg_82, fg_92, \
                         fg_97, fg_122, fg_127, fg_137, fg_142, fh_2, fh_7, fh_65, fh_70, \
                         fh_107, fh_112, fh_130, fh_137, fh_172, fh_179, fh_194, \
                         fh_201 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_22 * ab_x[k] * fg_2[k]
                  + f_18 * ab_x[k] * fg_7[k]
                  - f_23 * ab_x[k] * fg_47[k]
                  + f_19 * ab_x[k] * fg_52[k]
                  + f_21 * ab_x[k] * fg_77[k]
                  - f_20 * ab_x[k] * fg_82[k]
                  - f_22 * ab_y[k] * fg_92[k]
                  + f_18 * ab_y[k] * fg_97[k]
                  + f_21 * ab_y[k] * fg_122[k]
                  - f_20 * ab_y[k] * fg_127[k]
                  - f_24 * ab_z[k] * fg_137[k]
                  + f_21 * ab_z[k] * fg_142[k]
                  + f_22 * fh_2[k]
                  - f_18 * fh_7[k]
                  + f_23 * fh_65[k]
                  - f_19 * fh_70[k]
                  - f_21 * fh_107[k]
                  + f_20 * fh_112[k]
                  + f_22 * fh_130[k]
                  - f_18 * fh_137[k]
                  - f_21 * fh_172[k]
                  + f_20 * fh_179[k]
                  + f_24 * fh_194[k]
                  - f_21 * fh_201[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_0, fg_3, fg_10, fg_45, fg_48, fg_55, fg_75, \
                         fg_78, fg_85, fg_90, fg_93, fg_100, fg_120, fg_123, fg_130, fg_135, \
                         fg_138, fg_145, fh_0, fh_3, fh_10, fh_63, fh_66, fh_73, fh_105, \
                         fh_108, fh_115, fh_127, fh_132, fh_141, fh_169, fh_174, fh_183, \
                         fh_191, fh_196, fh_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_64 * ab_x[k] * fg_0[k]
                  + f_65 * ab_x[k] * fg_3[k]
                  - f_64 * ab_x[k] * fg_10[k]
                  - f_66 * ab_x[k] * fg_45[k]
                  + f_67 * ab_x[k] * fg_48[k]
                  - f_66 * ab_x[k] * fg_55[k]
                  + f_7 * ab_x[k] * fg_75[k]
                  - f_68 * ab_x[k] * fg_78[k]
                  + f_7 * ab_x[k] * fg_85[k]
                  - f_64 * ab_y[k] * fg_90[k]
                  + f_65 * ab_y[k] * fg_93[k]
                  - f_64 * ab_y[k] * fg_100[k]
                  + f_7 * ab_y[k] * fg_120[k]
                  - f_68 * ab_y[k] * fg_123[k]
                  + f_7 * ab_y[k] * fg_130[k]
                  - f_69 * ab_z[k] * fg_135[k]
                  + f_70 * ab_z[k] * fg_138[k]
                  - f_69 * ab_z[k] * fg_145[k]
                  + f_64 * fh_0[k]
                  - f_65 * fh_3[k]
                  + f_64 * fh_10[k]
                  + f_66 * fh_63[k]
                  - f_67 * fh_66[k]
                  + f_66 * fh_73[k]
                  - f_7 * fh_105[k]
                  + f_68 * fh_108[k]
                  - f_7 * fh_115[k]
                  + f_64 * fh_127[k]
                  - f_65 * fh_132[k]
                  + f_64 * fh_141[k]
                  - f_7 * fh_169[k]
                  + f_68 * fh_174[k]
                  - f_7 * fh_183[k]
                  + f_69 * fh_191[k]
                  - f_70 * fh_196[k]
                  + f_69 * fh_205[k];
    }

#pragma omp simd aligned(ab_x, fg_31, fg_36, fg_106, fg_111, fg_136, fg_141, fh_43, fh_48, \
                         fh_148, fh_153, fh_190, fh_195 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_4 * ab_x[k] * fg_31[k]
                  - f_4 * ab_x[k] * fg_36[k]
                  + f_4 * ab_x[k] * fg_106[k]
                  - f_4 * ab_x[k] * fg_111[k]
                  - f_5 * ab_x[k] * fg_136[k]
                  + f_5 * ab_x[k] * fg_141[k]
                  - f_4 * fh_43[k]
                  + f_4 * fh_48[k]
                  - f_4 * fh_148[k]
                  + f_4 * fh_153[k]
                  + f_5 * fh_190[k]
                  - f_5 * fh_195[k];
    }

#pragma omp simd aligned(ab_x, fg_34, fg_41, fg_109, fg_116, fg_139, fg_146, fh_46, fh_53, \
                         fh_151, fh_158, fh_193, fh_200 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_15 * ab_x[k] * fg_34[k]
                  - f_16 * ab_x[k] * fg_41[k]
                  + f_15 * ab_x[k] * fg_109[k]
                  - f_16 * ab_x[k] * fg_116[k]
                  - f_3 * ab_x[k] * fg_139[k]
                  + f_17 * ab_x[k] * fg_146[k]
                  - f_15 * fh_46[k]
                  + f_16 * fh_53[k]
                  - f_15 * fh_151[k]
                  + f_16 * fh_158[k]
                  + f_3 * fh_193[k]
                  - f_17 * fh_200[k];
    }

#pragma omp simd aligned(ab_x, fg_31, fg_36, fg_38, fg_106, fg_111, fg_113, fg_136, fg_141, \
                         fg_143, fh_43, fh_48, fh_50, fh_148, fh_153, fh_155, fh_190, fh_195, \
                         fh_197 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_32 * ab_x[k] * fg_31[k]
                  - f_32 * ab_x[k] * fg_36[k]
                  + f_34 * ab_x[k] * fg_38[k]
                  - f_32 * ab_x[k] * fg_106[k]
                  - f_32 * ab_x[k] * fg_111[k]
                  + f_34 * ab_x[k] * fg_113[k]
                  + f_33 * ab_x[k] * fg_136[k]
                  + f_33 * ab_x[k] * fg_141[k]
                  - f_35 * ab_x[k] * fg_143[k]
                  + f_32 * fh_43[k]
                  + f_32 * fh_48[k]
                  - f_34 * fh_50[k]
                  + f_32 * fh_148[k]
                  + f_32 * fh_153[k]
                  - f_34 * fh_155[k]
                  - f_33 * fh_190[k]
                  - f_33 * fh_195[k]
                  + f_35 * fh_197[k];
    }

#pragma omp simd aligned(ab_x, fg_34, fg_41, fg_43, fg_109, fg_116, fg_118, fg_139, fg_146, \
                         fg_148, fh_46, fh_53, fh_55, fh_151, fh_158, fh_160, fh_193, fh_200, \
                         fh_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -5.625 * ab_x[k] * fg_34[k]
                  - 5.625 * ab_x[k] * fg_41[k]
                  + 7.5 * ab_x[k] * fg_43[k]
                  - 5.625 * ab_x[k] * fg_109[k]
                  - 5.625 * ab_x[k] * fg_116[k]
                  + 7.5 * ab_x[k] * fg_118[k]
                  + 7.5 * ab_x[k] * fg_139[k]
                  + 7.5 * ab_x[k] * fg_146[k]
                  - 10.0 * ab_x[k] * fg_148[k]
                  + 5.625 * fh_46[k]
                  + 5.625 * fh_53[k]
                  - 7.5 * fh_55[k]
                  + 5.625 * fh_151[k]
                  + 5.625 * fh_158[k]
                  - 7.5 * fh_160[k]
                  - 7.5 * fh_193[k]
                  - 7.5 * fh_200[k]
                  + 10.0 * fh_202[k];
    }

#pragma omp simd aligned(ab_x, fg_30, fg_33, fg_35, fg_40, fg_42, fg_44, fg_105, fg_108, \
                         fg_110, fg_115, fg_117, fg_119, fg_135, fg_138, fg_140, fg_145, \
                         fg_147, fg_149, fh_42, fh_45, fh_47, fh_52, fh_54, fh_56, fh_147, \
                         fh_150, fh_152, fh_157, fh_159, fh_161, fh_189, fh_192, fh_194, \
                         fh_199, fh_201, fh_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_46 * ab_x[k] * fg_30[k]
                  + f_47 * ab_x[k] * fg_33[k]
                  - f_48 * ab_x[k] * fg_35[k]
                  + f_46 * ab_x[k] * fg_40[k]
                  - f_48 * ab_x[k] * fg_42[k]
                  + f_49 * ab_x[k] * fg_44[k]
                  + f_46 * ab_x[k] * fg_105[k]
                  + f_47 * ab_x[k] * fg_108[k]
                  - f_48 * ab_x[k] * fg_110[k]
                  + f_46 * ab_x[k] * fg_115[k]
                  - f_48 * ab_x[k] * fg_117[k]
                  + f_49 * ab_x[k] * fg_119[k]
                  - f_50 * ab_x[k] * fg_135[k]
                  - f_49 * ab_x[k] * fg_138[k]
                  + f_51 * ab_x[k] * fg_140[k]
                  - f_50 * ab_x[k] * fg_145[k]
                  + f_51 * ab_x[k] * fg_147[k]
                  - f_52 * ab_x[k] * fg_149[k]
                  - f_46 * fh_42[k]
                  - f_47 * fh_45[k]
                  + f_48 * fh_47[k]
                  - f_46 * fh_52[k]
                  + f_48 * fh_54[k]
                  - f_49 * fh_56[k]
                  - f_46 * fh_147[k]
                  - f_47 * fh_150[k]
                  + f_48 * fh_152[k]
                  - f_46 * fh_157[k]
                  + f_48 * fh_159[k]
                  - f_49 * fh_161[k]
                  + f_50 * fh_189[k]
                  + f_49 * fh_192[k]
                  - f_51 * fh_194[k]
                  + f_50 * fh_199[k]
                  - f_51 * fh_201[k]
                  + f_52 * fh_203[k];
    }

#pragma omp simd aligned(ab_x, fg_32, fg_37, fg_39, fg_107, fg_112, fg_114, fg_137, fg_142, \
                         fg_144, fh_44, fh_49, fh_51, fh_149, fh_154, fh_156, fh_191, fh_196, \
                         fh_198 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -5.625 * ab_x[k] * fg_32[k]
                  - 5.625 * ab_x[k] * fg_37[k]
                  + 7.5 * ab_x[k] * fg_39[k]
                  - 5.625 * ab_x[k] * fg_107[k]
                  - 5.625 * ab_x[k] * fg_112[k]
                  + 7.5 * ab_x[k] * fg_114[k]
                  + 7.5 * ab_x[k] * fg_137[k]
                  + 7.5 * ab_x[k] * fg_142[k]
                  - 10.0 * ab_x[k] * fg_144[k]
                  + 5.625 * fh_44[k]
                  + 5.625 * fh_49[k]
                  - 7.5 * fh_51[k]
                  + 5.625 * fh_149[k]
                  + 5.625 * fh_154[k]
                  - 7.5 * fh_156[k]
                  - 7.5 * fh_191[k]
                  - 7.5 * fh_196[k]
                  + 10.0 * fh_198[k];
    }

#pragma omp simd aligned(ab_x, fg_30, fg_35, fg_40, fg_42, fg_105, fg_110, fg_115, fg_117, \
                         fg_135, fg_140, fg_145, fg_147, fh_42, fh_47, fh_52, fh_54, fh_147, \
                         fh_152, fh_157, fh_159, fh_189, fh_194, fh_199, \
                         fh_201 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_53 * ab_x[k] * fg_30[k]
                  + f_54 * ab_x[k] * fg_35[k]
                  + f_53 * ab_x[k] * fg_40[k]
                  - f_54 * ab_x[k] * fg_42[k]
                  - f_53 * ab_x[k] * fg_105[k]
                  + f_54 * ab_x[k] * fg_110[k]
                  + f_53 * ab_x[k] * fg_115[k]
                  - f_54 * ab_x[k] * fg_117[k]
                  + f_55 * ab_x[k] * fg_135[k]
                  - f_56 * ab_x[k] * fg_140[k]
                  - f_55 * ab_x[k] * fg_145[k]
                  + f_56 * ab_x[k] * fg_147[k]
                  + f_53 * fh_42[k]
                  - f_54 * fh_47[k]
                  - f_53 * fh_52[k]
                  + f_54 * fh_54[k]
                  + f_53 * fh_147[k]
                  - f_54 * fh_152[k]
                  - f_53 * fh_157[k]
                  + f_54 * fh_159[k]
                  - f_55 * fh_189[k]
                  + f_56 * fh_194[k]
                  + f_55 * fh_199[k]
                  - f_56 * fh_201[k];
    }

#pragma omp simd aligned(ab_x, fg_32, fg_37, fg_107, fg_112, fg_137, fg_142, fh_44, fh_49, \
                         fh_149, fh_154, fh_191, fh_196 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_16 * ab_x[k] * fg_32[k]
                  - f_15 * ab_x[k] * fg_37[k]
                  + f_16 * ab_x[k] * fg_107[k]
                  - f_15 * ab_x[k] * fg_112[k]
                  - f_17 * ab_x[k] * fg_137[k]
                  + f_3 * ab_x[k] * fg_142[k]
                  - f_16 * fh_44[k]
                  + f_15 * fh_49[k]
                  - f_16 * fh_149[k]
                  + f_15 * fh_154[k]
                  + f_17 * fh_191[k]
                  - f_3 * fh_196[k];
    }

#pragma omp simd aligned(ab_x, fg_30, fg_33, fg_40, fg_105, fg_108, fg_115, fg_135, fg_138, \
                         fg_145, fh_42, fh_45, fh_52, fh_147, fh_150, fh_157, fh_189, fh_192, \
                         fh_199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_57 * ab_x[k] * fg_30[k]
                  - f_58 * ab_x[k] * fg_33[k]
                  + f_57 * ab_x[k] * fg_40[k]
                  + f_57 * ab_x[k] * fg_105[k]
                  - f_58 * ab_x[k] * fg_108[k]
                  + f_57 * ab_x[k] * fg_115[k]
                  - f_13 * ab_x[k] * fg_135[k]
                  + f_14 * ab_x[k] * fg_138[k]
                  - f_13 * ab_x[k] * fg_145[k]
                  - f_57 * fh_42[k]
                  + f_58 * fh_45[k]
                  - f_57 * fh_52[k]
                  - f_57 * fh_147[k]
                  + f_58 * fh_150[k]
                  - f_57 * fh_157[k]
                  + f_13 * fh_189[k]
                  - f_14 * fh_192[k]
                  + f_13 * fh_199[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_1, fg_6, fg_76, fg_81, fg_91, fg_96, fg_121, fg_126, \
                         fh_1, fh_6, fh_106, fh_111, fh_129, fh_136, fh_171, \
                         fh_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_10 * ab_x[k] * fg_1[k]
                  - f_10 * ab_x[k] * fg_6[k]
                  - f_11 * ab_x[k] * fg_76[k]
                  + f_11 * ab_x[k] * fg_81[k]
                  - f_10 * ab_y[k] * fg_91[k]
                  + f_10 * ab_y[k] * fg_96[k]
                  + f_11 * ab_y[k] * fg_121[k]
                  - f_11 * ab_y[k] * fg_126[k]
                  - f_10 * fh_1[k]
                  + f_10 * fh_6[k]
                  + f_11 * fh_106[k]
                  - f_11 * fh_111[k]
                  + f_10 * fh_129[k]
                  - f_10 * fh_136[k]
                  - f_11 * fh_171[k]
                  + f_11 * fh_178[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_4, fg_11, fg_79, fg_86, fg_94, fg_101, fg_124, fg_131, \
                         fh_4, fh_11, fh_109, fh_116, fh_133, fh_142, fh_175, \
                         fh_184 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_25 * ab_x[k] * fg_4[k]
                  - f_27 * ab_x[k] * fg_11[k]
                  - f_26 * ab_x[k] * fg_79[k]
                  + f_4 * ab_x[k] * fg_86[k]
                  - f_25 * ab_y[k] * fg_94[k]
                  + f_27 * ab_y[k] * fg_101[k]
                  + f_26 * ab_y[k] * fg_124[k]
                  - f_4 * ab_y[k] * fg_131[k]
                  - f_25 * fh_4[k]
                  + f_27 * fh_11[k]
                  + f_26 * fh_109[k]
                  - f_4 * fh_116[k]
                  + f_25 * fh_133[k]
                  - f_27 * fh_142[k]
                  - f_26 * fh_175[k]
                  + f_4 * fh_184[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_1, fg_6, fg_8, fg_76, fg_81, fg_83, fg_91, fg_96, \
                         fg_98, fg_121, fg_126, fg_128, fh_1, fh_6, fh_8, fh_106, fh_111, \
                         fh_113, fh_129, fh_136, fh_138, fh_171, fh_178, \
                         fh_180 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = -0.625 * ab_x[k] * fg_1[k]
                  - 0.625 * ab_x[k] * fg_6[k]
                  + 3.75 * ab_x[k] * fg_8[k]
                  + 3.75 * ab_x[k] * fg_76[k]
                  + 3.75 * ab_x[k] * fg_81[k]
                  - 22.5 * ab_x[k] * fg_83[k]
                  + 0.625 * ab_y[k] * fg_91[k]
                  + 0.625 * ab_y[k] * fg_96[k]
                  - 3.75 * ab_y[k] * fg_98[k]
                  - 3.75 * ab_y[k] * fg_121[k]
                  - 3.75 * ab_y[k] * fg_126[k]
                  + 22.5 * ab_y[k] * fg_128[k]
                  + 0.625 * fh_1[k]
                  + 0.625 * fh_6[k]
                  - 3.75 * fh_8[k]
                  - 3.75 * fh_106[k]
                  - 3.75 * fh_111[k]
                  + 22.5 * fh_113[k]
                  - 0.625 * fh_129[k]
                  - 0.625 * fh_136[k]
                  + 3.75 * fh_138[k]
                  + 3.75 * fh_171[k]
                  + 3.75 * fh_178[k]
                  - 22.5 * fh_180[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_4, fg_11, fg_13, fg_79, fg_86, fg_88, fg_94, fg_101, \
                         fg_103, fg_124, fg_131, fg_133, fh_4, fh_11, fh_13, fh_109, fh_116, \
                         fh_118, fh_133, fh_142, fh_144, fh_175, fh_184, \
                         fh_186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_53 * ab_x[k] * fg_4[k]
                  - f_53 * ab_x[k] * fg_11[k]
                  + f_55 * ab_x[k] * fg_13[k]
                  + f_54 * ab_x[k] * fg_79[k]
                  + f_54 * ab_x[k] * fg_86[k]
                  - f_56 * ab_x[k] * fg_88[k]
                  + f_53 * ab_y[k] * fg_94[k]
                  + f_53 * ab_y[k] * fg_101[k]
                  - f_55 * ab_y[k] * fg_103[k]
                  - f_54 * ab_y[k] * fg_124[k]
                  - f_54 * ab_y[k] * fg_131[k]
                  + f_56 * ab_y[k] * fg_133[k]
                  + f_53 * fh_4[k]
                  + f_53 * fh_11[k]
                  - f_55 * fh_13[k]
                  - f_54 * fh_109[k]
                  - f_54 * fh_116[k]
                  + f_56 * fh_118[k]
                  - f_53 * fh_133[k]
                  - f_53 * fh_142[k]
                  + f_55 * fh_144[k]
                  + f_54 * fh_175[k]
                  + f_54 * fh_184[k]
                  - f_56 * fh_186[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_0, fg_3, fg_5, fg_10, fg_12, fg_14, fg_75, fg_78, \
                         fg_80, fg_85, fg_87, fg_89, fg_90, fg_93, fg_95, fg_100, fg_102, \
                         fg_104, fg_120, fg_123, fg_125, fg_130, fg_132, fg_134, fh_0, fh_3, \
                         fh_5, fh_10, fh_12, fh_14, fh_105, fh_108, fh_110, fh_115, fh_117, \
                         fh_119, fh_127, fh_132, fh_134, fh_141, fh_143, fh_145, fh_169, \
                         fh_174, fh_176, fh_183, fh_185, fh_187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_59 * ab_x[k] * fg_0[k]
                  + f_36 * ab_x[k] * fg_3[k]
                  - f_61 * ab_x[k] * fg_5[k]
                  + f_59 * ab_x[k] * fg_10[k]
                  - f_61 * ab_x[k] * fg_12[k]
                  + f_63 * ab_x[k] * fg_14[k]
                  - f_60 * ab_x[k] * fg_75[k]
                  - f_40 * ab_x[k] * fg_78[k]
                  + f_62 * ab_x[k] * fg_80[k]
                  - f_60 * ab_x[k] * fg_85[k]
                  + f_62 * ab_x[k] * fg_87[k]
                  - f_38 * ab_x[k] * fg_89[k]
                  - f_59 * ab_y[k] * fg_90[k]
                  - f_36 * ab_y[k] * fg_93[k]
                  + f_61 * ab_y[k] * fg_95[k]
                  - f_59 * ab_y[k] * fg_100[k]
                  + f_61 * ab_y[k] * fg_102[k]
                  - f_63 * ab_y[k] * fg_104[k]
                  + f_60 * ab_y[k] * fg_120[k]
                  + f_40 * ab_y[k] * fg_123[k]
                  - f_62 * ab_y[k] * fg_125[k]
                  + f_60 * ab_y[k] * fg_130[k]
                  - f_62 * ab_y[k] * fg_132[k]
                  + f_38 * ab_y[k] * fg_134[k]
                  - f_59 * fh_0[k]
                  - f_36 * fh_3[k]
                  + f_61 * fh_5[k]
                  - f_59 * fh_10[k]
                  + f_61 * fh_12[k]
                  - f_63 * fh_14[k]
                  + f_60 * fh_105[k]
                  + f_40 * fh_108[k]
                  - f_62 * fh_110[k]
                  + f_60 * fh_115[k]
                  - f_62 * fh_117[k]
                  + f_38 * fh_119[k]
                  + f_59 * fh_127[k]
                  + f_36 * fh_132[k]
                  - f_61 * fh_134[k]
                  + f_59 * fh_141[k]
                  - f_61 * fh_143[k]
                  + f_63 * fh_145[k]
                  - f_60 * fh_169[k]
                  - f_40 * fh_174[k]
                  + f_62 * fh_176[k]
                  - f_60 * fh_183[k]
                  + f_62 * fh_185[k]
                  - f_38 * fh_187[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_2, fg_7, fg_9, fg_77, fg_82, fg_84, fg_92, fg_97, \
                         fg_99, fg_122, fg_127, fg_129, fh_2, fh_7, fh_9, fh_107, fh_112, \
                         fh_114, fh_130, fh_137, fh_139, fh_172, fh_179, \
                         fh_181 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_53 * ab_x[k] * fg_2[k]
                  - f_53 * ab_x[k] * fg_7[k]
                  + f_55 * ab_x[k] * fg_9[k]
                  + f_54 * ab_x[k] * fg_77[k]
                  + f_54 * ab_x[k] * fg_82[k]
                  - f_56 * ab_x[k] * fg_84[k]
                  + f_53 * ab_y[k] * fg_92[k]
                  + f_53 * ab_y[k] * fg_97[k]
                  - f_55 * ab_y[k] * fg_99[k]
                  - f_54 * ab_y[k] * fg_122[k]
                  - f_54 * ab_y[k] * fg_127[k]
                  + f_56 * ab_y[k] * fg_129[k]
                  + f_53 * fh_2[k]
                  + f_53 * fh_7[k]
                  - f_55 * fh_9[k]
                  - f_54 * fh_107[k]
                  - f_54 * fh_112[k]
                  + f_56 * fh_114[k]
                  - f_53 * fh_130[k]
                  - f_53 * fh_137[k]
                  + f_55 * fh_139[k]
                  + f_54 * fh_172[k]
                  + f_54 * fh_179[k]
                  - f_56 * fh_181[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_0, fg_5, fg_10, fg_12, fg_75, fg_80, fg_85, fg_87, \
                         fg_90, fg_95, fg_100, fg_102, fg_120, fg_125, fg_130, fg_132, fh_0, \
                         fh_5, fh_10, fh_12, fh_105, fh_110, fh_115, fh_117, fh_127, fh_134, \
                         fh_141, fh_143, fh_169, fh_176, fh_183, \
                         fh_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -0.3125 * ab_x[k] * fg_0[k]
                  + 1.875 * ab_x[k] * fg_5[k]
                  + 0.3125 * ab_x[k] * fg_10[k]
                  - 1.875 * ab_x[k] * fg_12[k]
                  + 1.875 * ab_x[k] * fg_75[k]
                  - 11.25 * ab_x[k] * fg_80[k]
                  - 1.875 * ab_x[k] * fg_85[k]
                  + 11.25 * ab_x[k] * fg_87[k]
                  + 0.3125 * ab_y[k] * fg_90[k]
                  - 1.875 * ab_y[k] * fg_95[k]
                  - 0.3125 * ab_y[k] * fg_100[k]
                  + 1.875 * ab_y[k] * fg_102[k]
                  - 1.875 * ab_y[k] * fg_120[k]
                  + 11.25 * ab_y[k] * fg_125[k]
                  + 1.875 * ab_y[k] * fg_130[k]
                  - 11.25 * ab_y[k] * fg_132[k]
                  + 0.3125 * fh_0[k]
                  - 1.875 * fh_5[k]
                  - 0.3125 * fh_10[k]
                  + 1.875 * fh_12[k]
                  - 1.875 * fh_105[k]
                  + 11.25 * fh_110[k]
                  + 1.875 * fh_115[k]
                  - 11.25 * fh_117[k]
                  - 0.3125 * fh_127[k]
                  + 1.875 * fh_134[k]
                  + 0.3125 * fh_141[k]
                  - 1.875 * fh_143[k]
                  + 1.875 * fh_169[k]
                  - 11.25 * fh_176[k]
                  - 1.875 * fh_183[k]
                  + 11.25 * fh_185[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_2, fg_7, fg_77, fg_82, fg_92, fg_97, fg_122, fg_127, \
                         fh_2, fh_7, fh_107, fh_112, fh_130, fh_137, fh_172, \
                         fh_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_27 * ab_x[k] * fg_2[k]
                  - f_25 * ab_x[k] * fg_7[k]
                  - f_4 * ab_x[k] * fg_77[k]
                  + f_26 * ab_x[k] * fg_82[k]
                  - f_27 * ab_y[k] * fg_92[k]
                  + f_25 * ab_y[k] * fg_97[k]
                  + f_4 * ab_y[k] * fg_122[k]
                  - f_26 * ab_y[k] * fg_127[k]
                  - f_27 * fh_2[k]
                  + f_25 * fh_7[k]
                  + f_4 * fh_107[k]
                  - f_26 * fh_112[k]
                  + f_27 * fh_130[k]
                  - f_25 * fh_137[k]
                  - f_4 * fh_172[k]
                  + f_26 * fh_179[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_0, fg_3, fg_10, fg_75, fg_78, fg_85, fg_90, fg_93, \
                         fg_100, fg_120, fg_123, fg_130, fh_0, fh_3, fh_10, fh_105, fh_108, \
                         fh_115, fh_127, fh_132, fh_141, fh_169, fh_174, \
                         fh_183 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_71 * ab_x[k] * fg_0[k]
                  - f_72 * ab_x[k] * fg_3[k]
                  + f_71 * ab_x[k] * fg_10[k]
                  - f_72 * ab_x[k] * fg_75[k]
                  + f_15 * ab_x[k] * fg_78[k]
                  - f_72 * ab_x[k] * fg_85[k]
                  - f_71 * ab_y[k] * fg_90[k]
                  + f_72 * ab_y[k] * fg_93[k]
                  - f_71 * ab_y[k] * fg_100[k]
                  + f_72 * ab_y[k] * fg_120[k]
                  - f_15 * ab_y[k] * fg_123[k]
                  + f_72 * ab_y[k] * fg_130[k]
                  - f_71 * fh_0[k]
                  + f_72 * fh_3[k]
                  - f_71 * fh_10[k]
                  + f_72 * fh_105[k]
                  - f_15 * fh_108[k]
                  + f_72 * fh_115[k]
                  + f_71 * fh_127[k]
                  - f_72 * fh_132[k]
                  + f_71 * fh_141[k]
                  - f_72 * fh_169[k]
                  + f_15 * fh_174[k]
                  - f_72 * fh_183[k];
    }

#pragma omp simd aligned(ab_x, fg_31, fg_36, fg_106, fg_111, fh_43, fh_48, fh_148, \
                         fh_153 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_1 * ab_x[k] * fg_31[k]
                  + f_1 * ab_x[k] * fg_36[k]
                  + f_0 * ab_x[k] * fg_106[k]
                  - f_0 * ab_x[k] * fg_111[k]
                  + f_1 * fh_43[k]
                  - f_1 * fh_48[k]
                  - f_0 * fh_148[k]
                  + f_0 * fh_153[k];
    }

#pragma omp simd aligned(ab_x, fg_34, fg_41, fg_109, fg_116, fh_46, fh_53, fh_151, \
                         fh_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -13.125 * ab_x[k] * fg_34[k]
                  + 4.375 * ab_x[k] * fg_41[k]
                  + 39.375 * ab_x[k] * fg_109[k]
                  - 13.125 * ab_x[k] * fg_116[k]
                  + 13.125 * fh_46[k]
                  - 4.375 * fh_53[k]
                  - 39.375 * fh_151[k]
                  + 13.125 * fh_158[k];
    }

#pragma omp simd aligned(ab_x, fg_31, fg_36, fg_38, fg_106, fg_111, fg_113, fh_43, fh_48, \
                         fh_50, fh_148, fh_153, fh_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_13 * ab_x[k] * fg_31[k]
                  + f_13 * ab_x[k] * fg_36[k]
                  - f_14 * ab_x[k] * fg_38[k]
                  - f_4 * ab_x[k] * fg_106[k]
                  - f_4 * ab_x[k] * fg_111[k]
                  + f_12 * ab_x[k] * fg_113[k]
                  - f_13 * fh_43[k]
                  - f_13 * fh_48[k]
                  + f_14 * fh_50[k]
                  + f_4 * fh_148[k]
                  + f_4 * fh_153[k]
                  - f_12 * fh_155[k];
    }

#pragma omp simd aligned(ab_x, fg_34, fg_41, fg_43, fg_109, fg_116, fg_118, fh_46, fh_53, \
                         fh_55, fh_151, fh_158, fh_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = f_16 * ab_x[k] * fg_34[k]
                  + f_16 * ab_x[k] * fg_41[k]
                  - f_17 * ab_x[k] * fg_43[k]
                  - f_15 * ab_x[k] * fg_109[k]
                  - f_15 * ab_x[k] * fg_116[k]
                  + f_3 * ab_x[k] * fg_118[k]
                  - f_16 * fh_46[k]
                  - f_16 * fh_53[k]
                  + f_17 * fh_55[k]
                  + f_15 * fh_151[k]
                  + f_15 * fh_158[k]
                  - f_3 * fh_160[k];
    }

#pragma omp simd aligned(ab_x, fg_30, fg_33, fg_35, fg_40, fg_42, fg_44, fg_105, fg_108, \
                         fg_110, fg_115, fg_117, fg_119, fh_42, fh_45, fh_47, fh_52, fh_54, \
                         fh_56, fh_147, fh_150, fh_152, fh_157, fh_159, \
                         fh_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_22 * ab_x[k] * fg_30[k]
                  - f_23 * ab_x[k] * fg_33[k]
                  + f_21 * ab_x[k] * fg_35[k]
                  - f_22 * ab_x[k] * fg_40[k]
                  + f_21 * ab_x[k] * fg_42[k]
                  - f_24 * ab_x[k] * fg_44[k]
                  + f_18 * ab_x[k] * fg_105[k]
                  + f_19 * ab_x[k] * fg_108[k]
                  - f_20 * ab_x[k] * fg_110[k]
                  + f_18 * ab_x[k] * fg_115[k]
                  - f_20 * ab_x[k] * fg_117[k]
                  + f_21 * ab_x[k] * fg_119[k]
                  + f_22 * fh_42[k]
                  + f_23 * fh_45[k]
                  - f_21 * fh_47[k]
                  + f_22 * fh_52[k]
                  - f_21 * fh_54[k]
                  + f_24 * fh_56[k]
                  - f_18 * fh_147[k]
                  - f_19 * fh_150[k]
                  + f_20 * fh_152[k]
                  - f_18 * fh_157[k]
                  + f_20 * fh_159[k]
                  - f_21 * fh_161[k];
    }

#pragma omp simd aligned(ab_x, fg_32, fg_37, fg_39, fg_107, fg_112, fg_114, fh_44, fh_49, \
                         fh_51, fh_149, fh_154, fh_156 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_16 * ab_x[k] * fg_32[k]
                  + f_16 * ab_x[k] * fg_37[k]
                  - f_17 * ab_x[k] * fg_39[k]
                  - f_15 * ab_x[k] * fg_107[k]
                  - f_15 * ab_x[k] * fg_112[k]
                  + f_3 * ab_x[k] * fg_114[k]
                  - f_16 * fh_44[k]
                  - f_16 * fh_49[k]
                  + f_17 * fh_51[k]
                  + f_15 * fh_149[k]
                  + f_15 * fh_154[k]
                  - f_3 * fh_156[k];
    }

#pragma omp simd aligned(ab_x, fg_30, fg_35, fg_40, fg_42, fg_105, fg_110, fg_115, fg_117, \
                         fh_42, fh_47, fh_52, fh_54, fh_147, fh_152, fh_157, \
                         fh_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_27 * ab_x[k] * fg_30[k]
                  - f_4 * ab_x[k] * fg_35[k]
                  - f_27 * ab_x[k] * fg_40[k]
                  + f_4 * ab_x[k] * fg_42[k]
                  - f_25 * ab_x[k] * fg_105[k]
                  + f_26 * ab_x[k] * fg_110[k]
                  + f_25 * ab_x[k] * fg_115[k]
                  - f_26 * ab_x[k] * fg_117[k]
                  - f_27 * fh_42[k]
                  + f_4 * fh_47[k]
                  + f_27 * fh_52[k]
                  - f_4 * fh_54[k]
                  + f_25 * fh_147[k]
                  - f_26 * fh_152[k]
                  - f_25 * fh_157[k]
                  + f_26 * fh_159[k];
    }

#pragma omp simd aligned(ab_x, fg_32, fg_37, fg_107, fg_112, fh_44, fh_49, fh_149, \
                         fh_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -4.375 * ab_x[k] * fg_32[k]
                  + 13.125 * ab_x[k] * fg_37[k]
                  + 13.125 * ab_x[k] * fg_107[k]
                  - 39.375 * ab_x[k] * fg_112[k]
                  + 4.375 * fh_44[k]
                  - 13.125 * fh_49[k]
                  - 13.125 * fh_149[k]
                  + 39.375 * fh_154[k];
    }

#pragma omp simd aligned(ab_x, fg_30, fg_33, fg_40, fg_105, fg_108, fg_115, fh_42, fh_45, \
                         fh_52, fh_147, fh_150, fh_157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_30 * ab_x[k] * fg_30[k]
                  + f_31 * ab_x[k] * fg_33[k]
                  - f_30 * ab_x[k] * fg_40[k]
                  + f_28 * ab_x[k] * fg_105[k]
                  - f_29 * ab_x[k] * fg_108[k]
                  + f_28 * ab_x[k] * fg_115[k]
                  + f_30 * fh_42[k]
                  - f_31 * fh_45[k]
                  + f_30 * fh_52[k]
                  - f_28 * fh_147[k]
                  + f_29 * fh_150[k]
                  - f_28 * fh_157[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_1, fg_6, fg_46, fg_51, fg_91, fg_96, fh_1, fh_6, \
                         fh_64, fh_69, fh_129, fh_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -2.1875 * ab_x[k] * fg_1[k]
                  + 2.1875 * ab_x[k] * fg_6[k]
                  + 13.125 * ab_x[k] * fg_46[k]
                  - 13.125 * ab_x[k] * fg_51[k]
                  - 2.1875 * ab_y[k] * fg_91[k]
                  + 2.1875 * ab_y[k] * fg_96[k]
                  + 2.1875 * fh_1[k]
                  - 2.1875 * fh_6[k]
                  - 13.125 * fh_64[k]
                  + 13.125 * fh_69[k]
                  + 2.1875 * fh_129[k]
                  - 2.1875 * fh_136[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_4, fg_11, fg_49, fg_56, fg_94, fg_101, fh_4, fh_11, \
                         fh_67, fh_74, fh_133, fh_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -f_28 * ab_x[k] * fg_4[k]
                  + f_30 * ab_x[k] * fg_11[k]
                  + f_29 * ab_x[k] * fg_49[k]
                  - f_31 * ab_x[k] * fg_56[k]
                  - f_28 * ab_y[k] * fg_94[k]
                  + f_30 * ab_y[k] * fg_101[k]
                  + f_28 * fh_4[k]
                  - f_30 * fh_11[k]
                  - f_29 * fh_67[k]
                  + f_31 * fh_74[k]
                  + f_28 * fh_133[k]
                  - f_30 * fh_142[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_1, fg_6, fg_8, fg_46, fg_51, fg_53, fg_91, fg_96, \
                         fg_98, fh_1, fh_6, fh_8, fh_64, fh_69, fh_71, fh_129, fh_136, \
                         fh_138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_44 * ab_x[k] * fg_1[k]
                  + f_44 * ab_x[k] * fg_6[k]
                  - f_16 * ab_x[k] * fg_8[k]
                  - f_16 * ab_x[k] * fg_46[k]
                  - f_16 * ab_x[k] * fg_51[k]
                  + f_45 * ab_x[k] * fg_53[k]
                  + f_44 * ab_y[k] * fg_91[k]
                  + f_44 * ab_y[k] * fg_96[k]
                  - f_16 * ab_y[k] * fg_98[k]
                  - f_44 * fh_1[k]
                  - f_44 * fh_6[k]
                  + f_16 * fh_8[k]
                  + f_16 * fh_64[k]
                  + f_16 * fh_69[k]
                  - f_45 * fh_71[k]
                  - f_44 * fh_129[k]
                  - f_44 * fh_136[k]
                  + f_16 * fh_138[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_4, fg_11, fg_13, fg_49, fg_56, fg_58, fg_94, fg_101, \
                         fg_103, fh_4, fh_11, fh_13, fh_67, fh_74, fh_76, fh_133, fh_142, \
                         fh_144 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_57 * ab_x[k] * fg_4[k]
                  + f_57 * ab_x[k] * fg_11[k]
                  - f_13 * ab_x[k] * fg_13[k]
                  - f_58 * ab_x[k] * fg_49[k]
                  - f_58 * ab_x[k] * fg_56[k]
                  + f_14 * ab_x[k] * fg_58[k]
                  + f_57 * ab_y[k] * fg_94[k]
                  + f_57 * ab_y[k] * fg_101[k]
                  - f_13 * ab_y[k] * fg_103[k]
                  - f_57 * fh_4[k]
                  - f_57 * fh_11[k]
                  + f_13 * fh_13[k]
                  + f_58 * fh_67[k]
                  + f_58 * fh_74[k]
                  - f_14 * fh_76[k]
                  - f_57 * fh_133[k]
                  - f_57 * fh_142[k]
                  + f_13 * fh_144[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_0, fg_3, fg_5, fg_10, fg_12, fg_14, fg_45, fg_48, \
                         fg_50, fg_55, fg_57, fg_59, fg_90, fg_93, fg_95, fg_100, fg_102, \
                         fg_104, fh_0, fh_3, fh_5, fh_10, fh_12, fh_14, fh_63, fh_66, fh_68, \
                         fh_73, fh_75, fh_77, fh_127, fh_132, fh_134, fh_141, fh_143, \
                         fh_145 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_64 * ab_x[k] * fg_0[k]
                  - f_66 * ab_x[k] * fg_3[k]
                  + f_7 * ab_x[k] * fg_5[k]
                  - f_64 * ab_x[k] * fg_10[k]
                  + f_7 * ab_x[k] * fg_12[k]
                  - f_69 * ab_x[k] * fg_14[k]
                  + f_65 * ab_x[k] * fg_45[k]
                  + f_67 * ab_x[k] * fg_48[k]
                  - f_68 * ab_x[k] * fg_50[k]
                  + f_65 * ab_x[k] * fg_55[k]
                  - f_68 * ab_x[k] * fg_57[k]
                  + f_70 * ab_x[k] * fg_59[k]
                  - f_64 * ab_y[k] * fg_90[k]
                  - f_66 * ab_y[k] * fg_93[k]
                  + f_7 * ab_y[k] * fg_95[k]
                  - f_64 * ab_y[k] * fg_100[k]
                  + f_7 * ab_y[k] * fg_102[k]
                  - f_69 * ab_y[k] * fg_104[k]
                  + f_64 * fh_0[k]
                  + f_66 * fh_3[k]
                  - f_7 * fh_5[k]
                  + f_64 * fh_10[k]
                  - f_7 * fh_12[k]
                  + f_69 * fh_14[k]
                  - f_65 * fh_63[k]
                  - f_67 * fh_66[k]
                  + f_68 * fh_68[k]
                  - f_65 * fh_73[k]
                  + f_68 * fh_75[k]
                  - f_70 * fh_77[k]
                  + f_64 * fh_127[k]
                  + f_66 * fh_132[k]
                  - f_7 * fh_134[k]
                  + f_64 * fh_141[k]
                  - f_7 * fh_143[k]
                  + f_69 * fh_145[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_2, fg_7, fg_9, fg_47, fg_52, fg_54, fg_92, fg_97, \
                         fg_99, fh_2, fh_7, fh_9, fh_65, fh_70, fh_72, fh_130, fh_137, \
                         fh_139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_57 * ab_x[k] * fg_2[k]
                  + f_57 * ab_x[k] * fg_7[k]
                  - f_13 * ab_x[k] * fg_9[k]
                  - f_58 * ab_x[k] * fg_47[k]
                  - f_58 * ab_x[k] * fg_52[k]
                  + f_14 * ab_x[k] * fg_54[k]
                  + f_57 * ab_y[k] * fg_92[k]
                  + f_57 * ab_y[k] * fg_97[k]
                  - f_13 * ab_y[k] * fg_99[k]
                  - f_57 * fh_2[k]
                  - f_57 * fh_7[k]
                  + f_13 * fh_9[k]
                  + f_58 * fh_65[k]
                  + f_58 * fh_70[k]
                  - f_14 * fh_72[k]
                  - f_57 * fh_130[k]
                  - f_57 * fh_137[k]
                  + f_13 * fh_139[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_0, fg_5, fg_10, fg_12, fg_45, fg_50, fg_55, fg_57, \
                         fg_90, fg_95, fg_100, fg_102, fh_0, fh_5, fh_10, fh_12, fh_63, fh_68, \
                         fh_73, fh_75, fh_127, fh_134, fh_141, fh_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = f_71 * ab_x[k] * fg_0[k]
                  - f_72 * ab_x[k] * fg_5[k]
                  - f_71 * ab_x[k] * fg_10[k]
                  + f_72 * ab_x[k] * fg_12[k]
                  - f_72 * ab_x[k] * fg_45[k]
                  + f_15 * ab_x[k] * fg_50[k]
                  + f_72 * ab_x[k] * fg_55[k]
                  - f_15 * ab_x[k] * fg_57[k]
                  + f_71 * ab_y[k] * fg_90[k]
                  - f_72 * ab_y[k] * fg_95[k]
                  - f_71 * ab_y[k] * fg_100[k]
                  + f_72 * ab_y[k] * fg_102[k]
                  - f_71 * fh_0[k]
                  + f_72 * fh_5[k]
                  + f_71 * fh_10[k]
                  - f_72 * fh_12[k]
                  + f_72 * fh_63[k]
                  - f_15 * fh_68[k]
                  - f_72 * fh_73[k]
                  + f_15 * fh_75[k]
                  - f_71 * fh_127[k]
                  + f_72 * fh_134[k]
                  + f_71 * fh_141[k]
                  - f_72 * fh_143[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_2, fg_7, fg_47, fg_52, fg_92, fg_97, fh_2, fh_7, \
                         fh_65, fh_70, fh_130, fh_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -f_30 * ab_x[k] * fg_2[k]
                  + f_28 * ab_x[k] * fg_7[k]
                  + f_31 * ab_x[k] * fg_47[k]
                  - f_29 * ab_x[k] * fg_52[k]
                  - f_30 * ab_y[k] * fg_92[k]
                  + f_28 * ab_y[k] * fg_97[k]
                  + f_30 * fh_2[k]
                  - f_28 * fh_7[k]
                  - f_31 * fh_65[k]
                  + f_29 * fh_70[k]
                  + f_30 * fh_130[k]
                  - f_28 * fh_137[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_0, fg_3, fg_10, fg_45, fg_48, fg_55, fg_90, fg_93, \
                         fg_100, fh_0, fh_3, fh_10, fh_63, fh_66, fh_73, fh_127, fh_132, \
                         fh_141 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -0.546875 * ab_x[k] * fg_0[k]
                  + 3.28125 * ab_x[k] * fg_3[k]
                  - 0.546875 * ab_x[k] * fg_10[k]
                  + 3.28125 * ab_x[k] * fg_45[k]
                  - 19.6875 * ab_x[k] * fg_48[k]
                  + 3.28125 * ab_x[k] * fg_55[k]
                  - 0.546875 * ab_y[k] * fg_90[k]
                  + 3.28125 * ab_y[k] * fg_93[k]
                  - 0.546875 * ab_y[k] * fg_100[k]
                  + 0.546875 * fh_0[k]
                  - 3.28125 * fh_3[k]
                  + 0.546875 * fh_10[k]
                  - 3.28125 * fh_63[k]
                  + 19.6875 * fh_66[k]
                  - 3.28125 * fh_73[k]
                  + 0.546875 * fh_127[k]
                  - 3.28125 * fh_132[k]
                  + 0.546875 * fh_141[k];
    }
}

auto
compute_hrr_gg_sph_tri(double *values, const size_t nvalues, CSimdMatrix &buffer,
                       const CSimdMatrix &coordinates, const size_t fg, const size_t fh,
                       const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 13.125 * std::sqrt(2.0);
    const auto f_1 = 4.375 * std::sqrt(2.0);
    const auto f_2 = 1.25 * std::sqrt(7.0);
    const auto f_3 = 7.5 * std::sqrt(7.0);
    const auto f_4 = 1.875 * std::sqrt(14.0);
    const auto f_5 = 2.5 * std::sqrt(14.0);
    const auto f_6 = 0.1875 * std::sqrt(35.0);
    const auto f_7 = 0.375 * std::sqrt(35.0);
    const auto f_8 = 1.5 * std::sqrt(35.0);
    const auto f_9 = 0.5 * std::sqrt(35.0);
    const auto f_10 = 0.625 * std::sqrt(7.0);
    const auto f_11 = 3.75 * std::sqrt(7.0);
    const auto f_12 = 11.25 * std::sqrt(14.0);
    const auto f_13 = 0.625 * std::sqrt(14.0);
    const auto f_14 = 3.75 * std::sqrt(14.0);
    const auto f_15 = 5.625 * std::sqrt(7.0);
    const auto f_16 = 1.875 * std::sqrt(7.0);
    const auto f_17 = 2.5 * std::sqrt(7.0);
    const auto f_18 = 0.28125 * std::sqrt(70.0);
    const auto f_19 = 0.5625 * std::sqrt(70.0);
    const auto f_20 = 2.25 * std::sqrt(70.0);
    const auto f_21 = 0.75 * std::sqrt(70.0);
    const auto f_22 = 0.09375 * std::sqrt(70.0);
    const auto f_23 = 0.1875 * std::sqrt(70.0);
    const auto f_24 = 0.25 * std::sqrt(70.0);
    const auto f_25 = 0.9375 * std::sqrt(14.0);
    const auto f_26 = 5.625 * std::sqrt(14.0);
    const auto f_27 = 0.3125 * std::sqrt(14.0);
    const auto f_28 = 3.28125 * std::sqrt(2.0);
    const auto f_29 = 19.6875 * std::sqrt(2.0);
    const auto f_30 = 1.09375 * std::sqrt(2.0);
    const auto f_31 = 6.5625 * std::sqrt(2.0);
    const auto f_32 = 1.875 * std::sqrt(2.0);
    const auto f_33 = 2.5 * std::sqrt(2.0);
    const auto f_34 = 11.25 * std::sqrt(2.0);
    const auto f_35 = 15.0 * std::sqrt(2.0);
    const auto f_36 = 0.1875 * std::sqrt(5.0);
    const auto f_37 = 0.375 * std::sqrt(5.0);
    const auto f_38 = 1.5 * std::sqrt(5.0);
    const auto f_39 = 0.5 * std::sqrt(5.0);
    const auto f_40 = 1.125 * std::sqrt(5.0);
    const auto f_41 = 2.25 * std::sqrt(5.0);
    const auto f_42 = 9.0 * std::sqrt(5.0);
    const auto f_43 = 3.0 * std::sqrt(5.0);
    const auto f_44 = 0.3125 * std::sqrt(7.0);
    const auto f_45 = 11.25 * std::sqrt(7.0);
    const auto f_46 = 0.28125 * std::sqrt(10.0);
    const auto f_47 = 0.5625 * std::sqrt(10.0);
    const auto f_48 = 2.25 * std::sqrt(10.0);
    const auto f_49 = 0.75 * std::sqrt(10.0);
    const auto f_50 = 0.375 * std::sqrt(10.0);
    const auto f_51 = 3.0 * std::sqrt(10.0);
    const auto f_52 = std::sqrt(10.0);
    const auto f_53 = 0.9375 * std::sqrt(2.0);
    const auto f_54 = 5.625 * std::sqrt(2.0);
    const auto f_55 = 1.25 * std::sqrt(2.0);
    const auto f_56 = 7.5 * std::sqrt(2.0);
    const auto f_57 = 0.46875 * std::sqrt(14.0);
    const auto f_58 = 2.8125 * std::sqrt(14.0);
    const auto f_59 = 0.09375 * std::sqrt(5.0);
    const auto f_60 = 0.5625 * std::sqrt(5.0);
    const auto f_61 = 0.75 * std::sqrt(5.0);
    const auto f_62 = 4.5 * std::sqrt(5.0);
    const auto f_63 = 0.25 * std::sqrt(5.0);
    const auto f_64 = 0.046875 * std::sqrt(35.0);
    const auto f_65 = 0.28125 * std::sqrt(35.0);
    const auto f_66 = 0.09375 * std::sqrt(35.0);
    const auto f_67 = 0.5625 * std::sqrt(35.0);
    const auto f_68 = 2.25 * std::sqrt(35.0);
    const auto f_69 = 0.125 * std::sqrt(35.0);
    const auto f_70 = 0.75 * std::sqrt(35.0);
    const auto f_71 = 0.15625 * std::sqrt(7.0);
    const auto f_72 = 0.9375 * std::sqrt(7.0);

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
    auto *g_77 = values + 77 * nvalues;
    auto *g_78 = values + 78 * nvalues;
    auto *g_79 = values + 79 * nvalues;
    auto *g_80 = values + 80 * nvalues;

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_61 = buffer.data(fg + 61);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_64 = buffer.data(fg + 64);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_66 = buffer.data(fg + 66);
    const auto *fg_67 = buffer.data(fg + 67);
    const auto *fg_68 = buffer.data(fg + 68);
    const auto *fg_69 = buffer.data(fg + 69);
    const auto *fg_70 = buffer.data(fg + 70);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_74 = buffer.data(fg + 74);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_78 = buffer.data(fg + 78);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_82 = buffer.data(fg + 82);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_85 = buffer.data(fg + 85);
    const auto *fg_87 = buffer.data(fg + 87);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_92 = buffer.data(fg + 92);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_94 = buffer.data(fg + 94);
    const auto *fg_95 = buffer.data(fg + 95);
    const auto *fg_96 = buffer.data(fg + 96);
    const auto *fg_97 = buffer.data(fg + 97);
    const auto *fg_98 = buffer.data(fg + 98);
    const auto *fg_99 = buffer.data(fg + 99);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_101 = buffer.data(fg + 101);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_103 = buffer.data(fg + 103);
    const auto *fg_104 = buffer.data(fg + 104);
    const auto *fg_105 = buffer.data(fg + 105);
    const auto *fg_106 = buffer.data(fg + 106);
    const auto *fg_107 = buffer.data(fg + 107);
    const auto *fg_108 = buffer.data(fg + 108);
    const auto *fg_109 = buffer.data(fg + 109);
    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_111 = buffer.data(fg + 111);
    const auto *fg_112 = buffer.data(fg + 112);
    const auto *fg_113 = buffer.data(fg + 113);
    const auto *fg_114 = buffer.data(fg + 114);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_120 = buffer.data(fg + 120);
    const auto *fg_121 = buffer.data(fg + 121);
    const auto *fg_122 = buffer.data(fg + 122);
    const auto *fg_123 = buffer.data(fg + 123);
    const auto *fg_124 = buffer.data(fg + 124);
    const auto *fg_125 = buffer.data(fg + 125);
    const auto *fg_126 = buffer.data(fg + 126);
    const auto *fg_127 = buffer.data(fg + 127);
    const auto *fg_128 = buffer.data(fg + 128);
    const auto *fg_129 = buffer.data(fg + 129);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_134 = buffer.data(fg + 134);
    const auto *fg_135 = buffer.data(fg + 135);
    const auto *fg_137 = buffer.data(fg + 137);
    const auto *fg_138 = buffer.data(fg + 138);
    const auto *fg_139 = buffer.data(fg + 139);
    const auto *fg_140 = buffer.data(fg + 140);
    const auto *fg_142 = buffer.data(fg + 142);
    const auto *fg_144 = buffer.data(fg + 144);
    const auto *fg_145 = buffer.data(fg + 145);
    const auto *fg_146 = buffer.data(fg + 146);
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_148 = buffer.data(fg + 148);
    const auto *fg_149 = buffer.data(fg + 149);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_70 = buffer.data(fh + 70);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_84 = buffer.data(fh + 84);
    const auto *fh_85 = buffer.data(fh + 85);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_88 = buffer.data(fh + 88);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_91 = buffer.data(fh + 91);
    const auto *fh_92 = buffer.data(fh + 92);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_94 = buffer.data(fh + 94);
    const auto *fh_95 = buffer.data(fh + 95);
    const auto *fh_96 = buffer.data(fh + 96);
    const auto *fh_97 = buffer.data(fh + 97);
    const auto *fh_98 = buffer.data(fh + 98);
    const auto *fh_105 = buffer.data(fh + 105);
    const auto *fh_107 = buffer.data(fh + 107);
    const auto *fh_108 = buffer.data(fh + 108);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_112 = buffer.data(fh + 112);
    const auto *fh_114 = buffer.data(fh + 114);
    const auto *fh_115 = buffer.data(fh + 115);
    const auto *fh_117 = buffer.data(fh + 117);
    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_126 = buffer.data(fh + 126);
    const auto *fh_127 = buffer.data(fh + 127);
    const auto *fh_128 = buffer.data(fh + 128);
    const auto *fh_129 = buffer.data(fh + 129);
    const auto *fh_130 = buffer.data(fh + 130);
    const auto *fh_131 = buffer.data(fh + 131);
    const auto *fh_132 = buffer.data(fh + 132);
    const auto *fh_133 = buffer.data(fh + 133);
    const auto *fh_134 = buffer.data(fh + 134);
    const auto *fh_135 = buffer.data(fh + 135);
    const auto *fh_136 = buffer.data(fh + 136);
    const auto *fh_137 = buffer.data(fh + 137);
    const auto *fh_138 = buffer.data(fh + 138);
    const auto *fh_139 = buffer.data(fh + 139);
    const auto *fh_140 = buffer.data(fh + 140);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_145 = buffer.data(fh + 145);
    const auto *fh_147 = buffer.data(fh + 147);
    const auto *fh_148 = buffer.data(fh + 148);
    const auto *fh_149 = buffer.data(fh + 149);
    const auto *fh_150 = buffer.data(fh + 150);
    const auto *fh_151 = buffer.data(fh + 151);
    const auto *fh_152 = buffer.data(fh + 152);
    const auto *fh_153 = buffer.data(fh + 153);
    const auto *fh_154 = buffer.data(fh + 154);
    const auto *fh_155 = buffer.data(fh + 155);
    const auto *fh_156 = buffer.data(fh + 156);
    const auto *fh_157 = buffer.data(fh + 157);
    const auto *fh_158 = buffer.data(fh + 158);
    const auto *fh_159 = buffer.data(fh + 159);
    const auto *fh_160 = buffer.data(fh + 160);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_163 = buffer.data(fh + 163);
    const auto *fh_164 = buffer.data(fh + 164);
    const auto *fh_165 = buffer.data(fh + 165);
    const auto *fh_166 = buffer.data(fh + 166);
    const auto *fh_168 = buffer.data(fh + 168);
    const auto *fh_169 = buffer.data(fh + 169);
    const auto *fh_170 = buffer.data(fh + 170);
    const auto *fh_171 = buffer.data(fh + 171);
    const auto *fh_172 = buffer.data(fh + 172);
    const auto *fh_173 = buffer.data(fh + 173);
    const auto *fh_174 = buffer.data(fh + 174);
    const auto *fh_175 = buffer.data(fh + 175);
    const auto *fh_176 = buffer.data(fh + 176);
    const auto *fh_177 = buffer.data(fh + 177);
    const auto *fh_178 = buffer.data(fh + 178);
    const auto *fh_179 = buffer.data(fh + 179);
    const auto *fh_180 = buffer.data(fh + 180);
    const auto *fh_181 = buffer.data(fh + 181);
    const auto *fh_182 = buffer.data(fh + 182);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_189 = buffer.data(fh + 189);
    const auto *fh_190 = buffer.data(fh + 190);
    const auto *fh_191 = buffer.data(fh + 191);
    const auto *fh_192 = buffer.data(fh + 192);
    const auto *fh_193 = buffer.data(fh + 193);
    const auto *fh_194 = buffer.data(fh + 194);
    const auto *fh_195 = buffer.data(fh + 195);
    const auto *fh_196 = buffer.data(fh + 196);
    const auto *fh_197 = buffer.data(fh + 197);
    const auto *fh_198 = buffer.data(fh + 198);
    const auto *fh_199 = buffer.data(fh + 199);
    const auto *fh_200 = buffer.data(fh + 200);
    const auto *fh_201 = buffer.data(fh + 201);
    const auto *fh_202 = buffer.data(fh + 202);
    const auto *fh_203 = buffer.data(fh + 203);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_205 = buffer.data(fh + 205);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_208 = buffer.data(fh + 208);
    const auto *fh_209 = buffer.data(fh + 209);

#pragma omp simd aligned(ab_x, fg_16, fg_21, fg_91, fg_96, fh_22, fh_27, fh_127, \
                         fh_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -8.75 * ab_x[k] * fg_16[k]
                 + 8.75 * ab_x[k] * fg_21[k]
                 + 8.75 * ab_x[k] * fg_91[k]
                 - 8.75 * ab_x[k] * fg_96[k]
                 + 8.75 * fh_22[k]
                 - 8.75 * fh_27[k]
                 - 8.75 * fh_127[k]
                 + 8.75 * fh_132[k];
    }

#pragma omp simd aligned(ab_x, fg_19, fg_26, fg_94, fg_101, fh_25, fh_32, fh_130, \
                         fh_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = -f_0 * ab_x[k] * fg_19[k]
                 + f_1 * ab_x[k] * fg_26[k]
                 + f_0 * ab_x[k] * fg_94[k]
                 - f_1 * ab_x[k] * fg_101[k]
                 + f_0 * fh_25[k]
                 - f_1 * fh_32[k]
                 - f_0 * fh_130[k]
                 + f_1 * fh_137[k];
        g_9[k] = g_1[k];
    }

#pragma omp simd aligned(ab_x, fg_16, fg_21, fg_23, fg_91, fg_96, fg_98, fh_22, fh_27, fh_29, \
                         fh_127, fh_132, fh_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = f_2 * ab_x[k] * fg_16[k]
                 + f_2 * ab_x[k] * fg_21[k]
                 - f_3 * ab_x[k] * fg_23[k]
                 - f_2 * ab_x[k] * fg_91[k]
                 - f_2 * ab_x[k] * fg_96[k]
                 + f_3 * ab_x[k] * fg_98[k]
                 - f_2 * fh_22[k]
                 - f_2 * fh_27[k]
                 + f_3 * fh_29[k]
                 + f_2 * fh_127[k]
                 + f_2 * fh_132[k]
                 - f_3 * fh_134[k];
        g_18[k] = g_2[k];
    }

#pragma omp simd aligned(ab_x, fg_19, fg_26, fg_28, fg_94, fg_101, fg_103, fh_25, fh_32, \
                         fh_34, fh_130, fh_137, fh_139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_4 * ab_x[k] * fg_19[k]
                 + f_4 * ab_x[k] * fg_26[k]
                 - f_5 * ab_x[k] * fg_28[k]
                 - f_4 * ab_x[k] * fg_94[k]
                 - f_4 * ab_x[k] * fg_101[k]
                 + f_5 * ab_x[k] * fg_103[k]
                 - f_4 * fh_25[k]
                 - f_4 * fh_32[k]
                 + f_5 * fh_34[k]
                 + f_4 * fh_130[k]
                 + f_4 * fh_137[k]
                 - f_5 * fh_139[k];
        g_27[k] = g_3[k];
    }

#pragma omp simd aligned(ab_x, fg_15, fg_18, fg_20, fg_25, fg_27, fg_29, fg_90, fg_93, fg_95, \
                         fg_100, fg_102, fg_104, fh_21, fh_24, fh_26, fh_31, fh_33, fh_35, \
                         fh_126, fh_129, fh_131, fh_136, fh_138, \
                         fh_140 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_6 * ab_x[k] * fg_15[k]
                 - f_7 * ab_x[k] * fg_18[k]
                 + f_8 * ab_x[k] * fg_20[k]
                 - f_6 * ab_x[k] * fg_25[k]
                 + f_8 * ab_x[k] * fg_27[k]
                 - f_9 * ab_x[k] * fg_29[k]
                 + f_6 * ab_x[k] * fg_90[k]
                 + f_7 * ab_x[k] * fg_93[k]
                 - f_8 * ab_x[k] * fg_95[k]
                 + f_6 * ab_x[k] * fg_100[k]
                 - f_8 * ab_x[k] * fg_102[k]
                 + f_9 * ab_x[k] * fg_104[k]
                 + f_6 * fh_21[k]
                 + f_7 * fh_24[k]
                 - f_8 * fh_26[k]
                 + f_6 * fh_31[k]
                 - f_8 * fh_33[k]
                 + f_9 * fh_35[k]
                 - f_6 * fh_126[k]
                 - f_7 * fh_129[k]
                 + f_8 * fh_131[k]
                 - f_6 * fh_136[k]
                 + f_8 * fh_138[k]
                 - f_9 * fh_140[k];
        g_36[k] = g_4[k];
    }

#pragma omp simd aligned(ab_x, fg_17, fg_22, fg_24, fg_92, fg_97, fg_99, fh_23, fh_28, fh_30, \
                         fh_128, fh_133, fh_135 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_4 * ab_x[k] * fg_17[k]
                 + f_4 * ab_x[k] * fg_22[k]
                 - f_5 * ab_x[k] * fg_24[k]
                 - f_4 * ab_x[k] * fg_92[k]
                 - f_4 * ab_x[k] * fg_97[k]
                 + f_5 * ab_x[k] * fg_99[k]
                 - f_4 * fh_23[k]
                 - f_4 * fh_28[k]
                 + f_5 * fh_30[k]
                 + f_4 * fh_128[k]
                 + f_4 * fh_133[k]
                 - f_5 * fh_135[k];
        g_45[k] = g_5[k];
    }

#pragma omp simd aligned(ab_x, fg_15, fg_20, fg_25, fg_27, fg_90, fg_95, fg_100, fg_102, \
                         fh_21, fh_26, fh_31, fh_33, fh_126, fh_131, fh_136, \
                         fh_138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_10 * ab_x[k] * fg_15[k]
                 - f_11 * ab_x[k] * fg_20[k]
                 - f_10 * ab_x[k] * fg_25[k]
                 + f_11 * ab_x[k] * fg_27[k]
                 - f_10 * ab_x[k] * fg_90[k]
                 + f_11 * ab_x[k] * fg_95[k]
                 + f_10 * ab_x[k] * fg_100[k]
                 - f_11 * ab_x[k] * fg_102[k]
                 - f_10 * fh_21[k]
                 + f_11 * fh_26[k]
                 + f_10 * fh_31[k]
                 - f_11 * fh_33[k]
                 + f_10 * fh_126[k]
                 - f_11 * fh_131[k]
                 - f_10 * fh_136[k]
                 + f_11 * fh_138[k];
        g_54[k] = g_6[k];
    }

#pragma omp simd aligned(ab_x, fg_17, fg_22, fg_92, fg_97, fh_23, fh_28, fh_128, \
                         fh_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_1 * ab_x[k] * fg_17[k]
                 + f_0 * ab_x[k] * fg_22[k]
                 + f_1 * ab_x[k] * fg_92[k]
                 - f_0 * ab_x[k] * fg_97[k]
                 + f_1 * fh_23[k]
                 - f_0 * fh_28[k]
                 - f_1 * fh_128[k]
                 + f_0 * fh_133[k];
        g_63[k] = g_7[k];
    }

#pragma omp simd aligned(ab_x, fg_15, fg_18, fg_25, fg_90, fg_93, fg_100, fh_21, fh_24, fh_31, \
                         fh_126, fh_129, fh_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -2.1875 * ab_x[k] * fg_15[k]
                 + 13.125 * ab_x[k] * fg_18[k]
                 - 2.1875 * ab_x[k] * fg_25[k]
                 + 2.1875 * ab_x[k] * fg_90[k]
                 - 13.125 * ab_x[k] * fg_93[k]
                 + 2.1875 * ab_x[k] * fg_100[k]
                 + 2.1875 * fh_21[k]
                 - 13.125 * fh_24[k]
                 + 2.1875 * fh_31[k]
                 - 2.1875 * fh_126[k]
                 + 13.125 * fh_129[k]
                 - 2.1875 * fh_136[k];
        g_72[k] = g_8[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_64, fg_71, fg_109, fg_116, fh_88, fh_95, fh_154, \
                         fh_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -39.375 * ab_x[k] * fg_64[k]
                  + 13.125 * ab_x[k] * fg_71[k]
                  + 13.125 * ab_y[k] * fg_109[k]
                  - 4.375 * ab_y[k] * fg_116[k]
                  + 39.375 * fh_88[k]
                  - 13.125 * fh_95[k]
                  - 13.125 * fh_154[k]
                  + 4.375 * fh_163[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_61, fg_66, fg_68, fg_106, fg_111, fg_113, fh_85, \
                         fh_90, fh_92, fh_150, fh_157, fh_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_4 * ab_x[k] * fg_61[k]
                  + f_4 * ab_x[k] * fg_66[k]
                  - f_12 * ab_x[k] * fg_68[k]
                  - f_13 * ab_y[k] * fg_106[k]
                  - f_13 * ab_y[k] * fg_111[k]
                  + f_14 * ab_y[k] * fg_113[k]
                  - f_4 * fh_85[k]
                  - f_4 * fh_90[k]
                  + f_12 * fh_92[k]
                  + f_13 * fh_150[k]
                  + f_13 * fh_157[k]
                  - f_14 * fh_159[k];
        g_19[k] = g_11[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_64, fg_71, fg_73, fg_109, fg_116, fg_118, fh_88, \
                         fh_95, fh_97, fh_154, fh_163, fh_165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_15 * ab_x[k] * fg_64[k]
                  + f_15 * ab_x[k] * fg_71[k]
                  - f_3 * ab_x[k] * fg_73[k]
                  - f_16 * ab_y[k] * fg_109[k]
                  - f_16 * ab_y[k] * fg_116[k]
                  + f_17 * ab_y[k] * fg_118[k]
                  - f_15 * fh_88[k]
                  - f_15 * fh_95[k]
                  + f_3 * fh_97[k]
                  + f_16 * fh_154[k]
                  + f_16 * fh_163[k]
                  - f_17 * fh_165[k];
        g_28[k] = g_12[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_60, fg_63, fg_65, fg_70, fg_72, fg_74, fg_105, fg_108, \
                         fg_110, fg_115, fg_117, fg_119, fh_84, fh_87, fh_89, fh_94, fh_96, \
                         fh_98, fh_148, fh_153, fh_155, fh_162, fh_164, \
                         fh_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_18 * ab_x[k] * fg_60[k]
                  - f_19 * ab_x[k] * fg_63[k]
                  + f_20 * ab_x[k] * fg_65[k]
                  - f_18 * ab_x[k] * fg_70[k]
                  + f_20 * ab_x[k] * fg_72[k]
                  - f_21 * ab_x[k] * fg_74[k]
                  + f_22 * ab_y[k] * fg_105[k]
                  + f_23 * ab_y[k] * fg_108[k]
                  - f_21 * ab_y[k] * fg_110[k]
                  + f_22 * ab_y[k] * fg_115[k]
                  - f_21 * ab_y[k] * fg_117[k]
                  + f_24 * ab_y[k] * fg_119[k]
                  + f_18 * fh_84[k]
                  + f_19 * fh_87[k]
                  - f_20 * fh_89[k]
                  + f_18 * fh_94[k]
                  - f_20 * fh_96[k]
                  + f_21 * fh_98[k]
                  - f_22 * fh_148[k]
                  - f_23 * fh_153[k]
                  + f_21 * fh_155[k]
                  - f_22 * fh_162[k]
                  + f_21 * fh_164[k]
                  - f_24 * fh_166[k];
        g_37[k] = g_13[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_62, fg_67, fg_69, fg_107, fg_112, fg_114, fh_86, \
                         fh_91, fh_93, fh_151, fh_158, fh_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_15 * ab_x[k] * fg_62[k]
                  + f_15 * ab_x[k] * fg_67[k]
                  - f_3 * ab_x[k] * fg_69[k]
                  - f_16 * ab_y[k] * fg_107[k]
                  - f_16 * ab_y[k] * fg_112[k]
                  + f_17 * ab_y[k] * fg_114[k]
                  - f_15 * fh_86[k]
                  - f_15 * fh_91[k]
                  + f_3 * fh_93[k]
                  + f_16 * fh_151[k]
                  + f_16 * fh_158[k]
                  - f_17 * fh_160[k];
        g_46[k] = g_14[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_60, fg_65, fg_70, fg_72, fg_105, fg_110, fg_115, \
                         fg_117, fh_84, fh_89, fh_94, fh_96, fh_148, fh_155, fh_162, \
                         fh_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_25 * ab_x[k] * fg_60[k]
                  - f_26 * ab_x[k] * fg_65[k]
                  - f_25 * ab_x[k] * fg_70[k]
                  + f_26 * ab_x[k] * fg_72[k]
                  - f_27 * ab_y[k] * fg_105[k]
                  + f_4 * ab_y[k] * fg_110[k]
                  + f_27 * ab_y[k] * fg_115[k]
                  - f_4 * ab_y[k] * fg_117[k]
                  - f_25 * fh_84[k]
                  + f_26 * fh_89[k]
                  + f_25 * fh_94[k]
                  - f_26 * fh_96[k]
                  + f_27 * fh_148[k]
                  - f_4 * fh_155[k]
                  - f_27 * fh_162[k]
                  + f_4 * fh_164[k];
        g_55[k] = g_15[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_62, fg_67, fg_107, fg_112, fh_86, fh_91, fh_151, \
                         fh_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -13.125 * ab_x[k] * fg_62[k]
                  + 39.375 * ab_x[k] * fg_67[k]
                  + 4.375 * ab_y[k] * fg_107[k]
                  - 13.125 * ab_y[k] * fg_112[k]
                  + 13.125 * fh_86[k]
                  - 39.375 * fh_91[k]
                  - 4.375 * fh_151[k]
                  + 13.125 * fh_158[k];
        g_64[k] = g_16[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_60, fg_63, fg_70, fg_105, fg_108, fg_115, fh_84, \
                         fh_87, fh_94, fh_148, fh_153, fh_162 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_28 * ab_x[k] * fg_60[k]
                  + f_29 * ab_x[k] * fg_63[k]
                  - f_28 * ab_x[k] * fg_70[k]
                  + f_30 * ab_y[k] * fg_105[k]
                  - f_31 * ab_y[k] * fg_108[k]
                  + f_30 * ab_y[k] * fg_115[k]
                  + f_28 * fh_84[k]
                  - f_29 * fh_87[k]
                  + f_28 * fh_94[k]
                  - f_30 * fh_148[k]
                  + f_31 * fh_153[k]
                  - f_30 * fh_162[k];
        g_73[k] = g_17[k];
    }

#pragma omp simd aligned(ab_x, fg_16, fg_21, fg_23, fg_91, fg_96, fg_98, fg_121, fg_126, \
                         fg_128, fh_22, fh_27, fh_29, fh_127, fh_132, fh_134, fh_169, fh_174, \
                         fh_176 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -1.25 * ab_x[k] * fg_16[k]
                  - 1.25 * ab_x[k] * fg_21[k]
                  + 7.5 * ab_x[k] * fg_23[k]
                  - 1.25 * ab_x[k] * fg_91[k]
                  - 1.25 * ab_x[k] * fg_96[k]
                  + 7.5 * ab_x[k] * fg_98[k]
                  + 7.5 * ab_x[k] * fg_121[k]
                  + 7.5 * ab_x[k] * fg_126[k]
                  - 45.0 * ab_x[k] * fg_128[k]
                  + 1.25 * fh_22[k]
                  + 1.25 * fh_27[k]
                  - 7.5 * fh_29[k]
                  + 1.25 * fh_127[k]
                  + 1.25 * fh_132[k]
                  - 7.5 * fh_134[k]
                  - 7.5 * fh_169[k]
                  - 7.5 * fh_174[k]
                  + 45.0 * fh_176[k];
    }

#pragma omp simd aligned(ab_x, fg_19, fg_26, fg_28, fg_94, fg_101, fg_103, fg_124, fg_131, \
                         fg_133, fh_25, fh_32, fh_34, fh_130, fh_137, fh_139, fh_172, fh_179, \
                         fh_181 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_32 * ab_x[k] * fg_19[k]
                  - f_32 * ab_x[k] * fg_26[k]
                  + f_33 * ab_x[k] * fg_28[k]
                  - f_32 * ab_x[k] * fg_94[k]
                  - f_32 * ab_x[k] * fg_101[k]
                  + f_33 * ab_x[k] * fg_103[k]
                  + f_34 * ab_x[k] * fg_124[k]
                  + f_34 * ab_x[k] * fg_131[k]
                  - f_35 * ab_x[k] * fg_133[k]
                  + f_32 * fh_25[k]
                  + f_32 * fh_32[k]
                  - f_33 * fh_34[k]
                  + f_32 * fh_130[k]
                  + f_32 * fh_137[k]
                  - f_33 * fh_139[k]
                  - f_34 * fh_172[k]
                  - f_34 * fh_179[k]
                  + f_35 * fh_181[k];
        g_29[k] = g_21[k];
    }

#pragma omp simd aligned(ab_x, fg_15, fg_18, fg_20, fg_25, fg_27, fg_29, fg_90, fg_93, fg_95, \
                         fg_100, fg_102, fg_104, fg_120, fg_123, fg_125, fg_130, fg_132, \
                         fg_134, fh_21, fh_24, fh_26, fh_31, fh_33, fh_35, fh_126, fh_129, \
                         fh_131, fh_136, fh_138, fh_140, fh_168, fh_171, fh_173, fh_178, \
                         fh_180, fh_182 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_36 * ab_x[k] * fg_15[k]
                  + f_37 * ab_x[k] * fg_18[k]
                  - f_38 * ab_x[k] * fg_20[k]
                  + f_36 * ab_x[k] * fg_25[k]
                  - f_38 * ab_x[k] * fg_27[k]
                  + f_39 * ab_x[k] * fg_29[k]
                  + f_36 * ab_x[k] * fg_90[k]
                  + f_37 * ab_x[k] * fg_93[k]
                  - f_38 * ab_x[k] * fg_95[k]
                  + f_36 * ab_x[k] * fg_100[k]
                  - f_38 * ab_x[k] * fg_102[k]
                  + f_39 * ab_x[k] * fg_104[k]
                  - f_40 * ab_x[k] * fg_120[k]
                  - f_41 * ab_x[k] * fg_123[k]
                  + f_42 * ab_x[k] * fg_125[k]
                  - f_40 * ab_x[k] * fg_130[k]
                  + f_42 * ab_x[k] * fg_132[k]
                  - f_43 * ab_x[k] * fg_134[k]
                  - f_36 * fh_21[k]
                  - f_37 * fh_24[k]
                  + f_38 * fh_26[k]
                  - f_36 * fh_31[k]
                  + f_38 * fh_33[k]
                  - f_39 * fh_35[k]
                  - f_36 * fh_126[k]
                  - f_37 * fh_129[k]
                  + f_38 * fh_131[k]
                  - f_36 * fh_136[k]
                  + f_38 * fh_138[k]
                  - f_39 * fh_140[k]
                  + f_40 * fh_168[k]
                  + f_41 * fh_171[k]
                  - f_42 * fh_173[k]
                  + f_40 * fh_178[k]
                  - f_42 * fh_180[k]
                  + f_43 * fh_182[k];
        g_38[k] = g_22[k];
    }

#pragma omp simd aligned(ab_x, fg_17, fg_22, fg_24, fg_92, fg_97, fg_99, fg_122, fg_127, \
                         fg_129, fh_23, fh_28, fh_30, fh_128, fh_133, fh_135, fh_170, fh_175, \
                         fh_177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_32 * ab_x[k] * fg_17[k]
                  - f_32 * ab_x[k] * fg_22[k]
                  + f_33 * ab_x[k] * fg_24[k]
                  - f_32 * ab_x[k] * fg_92[k]
                  - f_32 * ab_x[k] * fg_97[k]
                  + f_33 * ab_x[k] * fg_99[k]
                  + f_34 * ab_x[k] * fg_122[k]
                  + f_34 * ab_x[k] * fg_127[k]
                  - f_35 * ab_x[k] * fg_129[k]
                  + f_32 * fh_23[k]
                  + f_32 * fh_28[k]
                  - f_33 * fh_30[k]
                  + f_32 * fh_128[k]
                  + f_32 * fh_133[k]
                  - f_33 * fh_135[k]
                  - f_34 * fh_170[k]
                  - f_34 * fh_175[k]
                  + f_35 * fh_177[k];
        g_47[k] = g_23[k];
    }

#pragma omp simd aligned(ab_x, fg_15, fg_20, fg_25, fg_27, fg_90, fg_95, fg_100, fg_102, \
                         fg_120, fg_125, fg_130, fg_132, fh_21, fh_26, fh_31, fh_33, fh_126, \
                         fh_131, fh_136, fh_138, fh_168, fh_173, fh_178, \
                         fh_180 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -0.625 * ab_x[k] * fg_15[k]
                  + 3.75 * ab_x[k] * fg_20[k]
                  + 0.625 * ab_x[k] * fg_25[k]
                  - 3.75 * ab_x[k] * fg_27[k]
                  - 0.625 * ab_x[k] * fg_90[k]
                  + 3.75 * ab_x[k] * fg_95[k]
                  + 0.625 * ab_x[k] * fg_100[k]
                  - 3.75 * ab_x[k] * fg_102[k]
                  + 3.75 * ab_x[k] * fg_120[k]
                  - 22.5 * ab_x[k] * fg_125[k]
                  - 3.75 * ab_x[k] * fg_130[k]
                  + 22.5 * ab_x[k] * fg_132[k]
                  + 0.625 * fh_21[k]
                  - 3.75 * fh_26[k]
                  - 0.625 * fh_31[k]
                  + 3.75 * fh_33[k]
                  + 0.625 * fh_126[k]
                  - 3.75 * fh_131[k]
                  - 0.625 * fh_136[k]
                  + 3.75 * fh_138[k]
                  - 3.75 * fh_168[k]
                  + 22.5 * fh_173[k]
                  + 3.75 * fh_178[k]
                  - 22.5 * fh_180[k];
        g_56[k] = g_24[k];
    }

#pragma omp simd aligned(ab_x, fg_17, fg_22, fg_92, fg_97, fg_122, fg_127, fh_23, fh_28, \
                         fh_128, fh_133, fh_170, fh_175 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_13 * ab_x[k] * fg_17[k]
                  - f_4 * ab_x[k] * fg_22[k]
                  + f_13 * ab_x[k] * fg_92[k]
                  - f_4 * ab_x[k] * fg_97[k]
                  - f_14 * ab_x[k] * fg_122[k]
                  + f_12 * ab_x[k] * fg_127[k]
                  - f_13 * fh_23[k]
                  + f_4 * fh_28[k]
                  - f_13 * fh_128[k]
                  + f_4 * fh_133[k]
                  + f_14 * fh_170[k]
                  - f_12 * fh_175[k];
        g_65[k] = g_25[k];
    }

#pragma omp simd aligned(ab_x, fg_15, fg_18, fg_25, fg_90, fg_93, fg_100, fg_120, fg_123, \
                         fg_130, fh_21, fh_24, fh_31, fh_126, fh_129, fh_136, fh_168, fh_171, \
                         fh_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_44 * ab_x[k] * fg_15[k]
                  - f_16 * ab_x[k] * fg_18[k]
                  + f_44 * ab_x[k] * fg_25[k]
                  + f_44 * ab_x[k] * fg_90[k]
                  - f_16 * ab_x[k] * fg_93[k]
                  + f_44 * ab_x[k] * fg_100[k]
                  - f_16 * ab_x[k] * fg_120[k]
                  + f_45 * ab_x[k] * fg_123[k]
                  - f_16 * ab_x[k] * fg_130[k]
                  - f_44 * fh_21[k]
                  + f_16 * fh_24[k]
                  - f_44 * fh_31[k]
                  - f_44 * fh_126[k]
                  + f_16 * fh_129[k]
                  - f_44 * fh_136[k]
                  + f_16 * fh_168[k]
                  - f_45 * fh_171[k]
                  + f_16 * fh_178[k];
        g_74[k] = g_26[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_64, fg_71, fg_73, fg_109, fg_116, fg_118, fg_139, \
                         fg_146, fg_148, fh_88, fh_95, fh_97, fh_154, fh_163, fh_165, fh_196, \
                         fh_205, fh_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -5.625 * ab_x[k] * fg_64[k]
                  - 5.625 * ab_x[k] * fg_71[k]
                  + 7.5 * ab_x[k] * fg_73[k]
                  - 5.625 * ab_y[k] * fg_109[k]
                  - 5.625 * ab_y[k] * fg_116[k]
                  + 7.5 * ab_y[k] * fg_118[k]
                  + 7.5 * ab_y[k] * fg_139[k]
                  + 7.5 * ab_y[k] * fg_146[k]
                  - 10.0 * ab_y[k] * fg_148[k]
                  + 5.625 * fh_88[k]
                  + 5.625 * fh_95[k]
                  - 7.5 * fh_97[k]
                  + 5.625 * fh_154[k]
                  + 5.625 * fh_163[k]
                  - 7.5 * fh_165[k]
                  - 7.5 * fh_196[k]
                  - 7.5 * fh_205[k]
                  + 10.0 * fh_207[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_60, fg_63, fg_65, fg_70, fg_72, fg_74, fg_105, fg_108, \
                         fg_110, fg_115, fg_117, fg_119, fg_135, fg_138, fg_140, fg_145, \
                         fg_147, fg_149, fh_84, fh_87, fh_89, fh_94, fh_96, fh_98, fh_148, \
                         fh_153, fh_155, fh_162, fh_164, fh_166, fh_190, fh_195, fh_197, \
                         fh_204, fh_206, fh_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = f_46 * ab_x[k] * fg_60[k]
                  + f_47 * ab_x[k] * fg_63[k]
                  - f_48 * ab_x[k] * fg_65[k]
                  + f_46 * ab_x[k] * fg_70[k]
                  - f_48 * ab_x[k] * fg_72[k]
                  + f_49 * ab_x[k] * fg_74[k]
                  + f_46 * ab_y[k] * fg_105[k]
                  + f_47 * ab_y[k] * fg_108[k]
                  - f_48 * ab_y[k] * fg_110[k]
                  + f_46 * ab_y[k] * fg_115[k]
                  - f_48 * ab_y[k] * fg_117[k]
                  + f_49 * ab_y[k] * fg_119[k]
                  - f_50 * ab_y[k] * fg_135[k]
                  - f_49 * ab_y[k] * fg_138[k]
                  + f_51 * ab_y[k] * fg_140[k]
                  - f_50 * ab_y[k] * fg_145[k]
                  + f_51 * ab_y[k] * fg_147[k]
                  - f_52 * ab_y[k] * fg_149[k]
                  - f_46 * fh_84[k]
                  - f_47 * fh_87[k]
                  + f_48 * fh_89[k]
                  - f_46 * fh_94[k]
                  + f_48 * fh_96[k]
                  - f_49 * fh_98[k]
                  - f_46 * fh_148[k]
                  - f_47 * fh_153[k]
                  + f_48 * fh_155[k]
                  - f_46 * fh_162[k]
                  + f_48 * fh_164[k]
                  - f_49 * fh_166[k]
                  + f_50 * fh_190[k]
                  + f_49 * fh_195[k]
                  - f_51 * fh_197[k]
                  + f_50 * fh_204[k]
                  - f_51 * fh_206[k]
                  + f_52 * fh_208[k];
        g_39[k] = g_31[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_62, fg_67, fg_69, fg_107, fg_112, fg_114, fg_137, \
                         fg_142, fg_144, fh_86, fh_91, fh_93, fh_151, fh_158, fh_160, fh_193, \
                         fh_200, fh_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -5.625 * ab_x[k] * fg_62[k]
                  - 5.625 * ab_x[k] * fg_67[k]
                  + 7.5 * ab_x[k] * fg_69[k]
                  - 5.625 * ab_y[k] * fg_107[k]
                  - 5.625 * ab_y[k] * fg_112[k]
                  + 7.5 * ab_y[k] * fg_114[k]
                  + 7.5 * ab_y[k] * fg_137[k]
                  + 7.5 * ab_y[k] * fg_142[k]
                  - 10.0 * ab_y[k] * fg_144[k]
                  + 5.625 * fh_86[k]
                  + 5.625 * fh_91[k]
                  - 7.5 * fh_93[k]
                  + 5.625 * fh_151[k]
                  + 5.625 * fh_158[k]
                  - 7.5 * fh_160[k]
                  - 7.5 * fh_193[k]
                  - 7.5 * fh_200[k]
                  + 10.0 * fh_202[k];
        g_48[k] = g_32[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_60, fg_65, fg_70, fg_72, fg_105, fg_110, fg_115, \
                         fg_117, fg_135, fg_140, fg_145, fg_147, fh_84, fh_89, fh_94, fh_96, \
                         fh_148, fh_155, fh_162, fh_164, fh_190, fh_197, fh_204, \
                         fh_206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_53 * ab_x[k] * fg_60[k]
                  + f_54 * ab_x[k] * fg_65[k]
                  + f_53 * ab_x[k] * fg_70[k]
                  - f_54 * ab_x[k] * fg_72[k]
                  - f_53 * ab_y[k] * fg_105[k]
                  + f_54 * ab_y[k] * fg_110[k]
                  + f_53 * ab_y[k] * fg_115[k]
                  - f_54 * ab_y[k] * fg_117[k]
                  + f_55 * ab_y[k] * fg_135[k]
                  - f_56 * ab_y[k] * fg_140[k]
                  - f_55 * ab_y[k] * fg_145[k]
                  + f_56 * ab_y[k] * fg_147[k]
                  + f_53 * fh_84[k]
                  - f_54 * fh_89[k]
                  - f_53 * fh_94[k]
                  + f_54 * fh_96[k]
                  + f_53 * fh_148[k]
                  - f_54 * fh_155[k]
                  - f_53 * fh_162[k]
                  + f_54 * fh_164[k]
                  - f_55 * fh_190[k]
                  + f_56 * fh_197[k]
                  + f_55 * fh_204[k]
                  - f_56 * fh_206[k];
        g_57[k] = g_33[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_62, fg_67, fg_107, fg_112, fg_137, fg_142, fh_86, \
                         fh_91, fh_151, fh_158, fh_193, fh_200 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_16 * ab_x[k] * fg_62[k]
                  - f_15 * ab_x[k] * fg_67[k]
                  + f_16 * ab_y[k] * fg_107[k]
                  - f_15 * ab_y[k] * fg_112[k]
                  - f_17 * ab_y[k] * fg_137[k]
                  + f_3 * ab_y[k] * fg_142[k]
                  - f_16 * fh_86[k]
                  + f_15 * fh_91[k]
                  - f_16 * fh_151[k]
                  + f_15 * fh_158[k]
                  + f_17 * fh_193[k]
                  - f_3 * fh_200[k];
        g_66[k] = g_34[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_60, fg_63, fg_70, fg_105, fg_108, fg_115, fg_135, \
                         fg_138, fg_145, fh_84, fh_87, fh_94, fh_148, fh_153, fh_162, fh_190, \
                         fh_195, fh_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_57 * ab_x[k] * fg_60[k]
                  - f_58 * ab_x[k] * fg_63[k]
                  + f_57 * ab_x[k] * fg_70[k]
                  + f_57 * ab_y[k] * fg_105[k]
                  - f_58 * ab_y[k] * fg_108[k]
                  + f_57 * ab_y[k] * fg_115[k]
                  - f_13 * ab_y[k] * fg_135[k]
                  + f_14 * ab_y[k] * fg_138[k]
                  - f_13 * ab_y[k] * fg_145[k]
                  - f_57 * fh_84[k]
                  + f_58 * fh_87[k]
                  - f_57 * fh_94[k]
                  - f_57 * fh_148[k]
                  + f_58 * fh_153[k]
                  - f_57 * fh_162[k]
                  + f_13 * fh_190[k]
                  - f_14 * fh_195[k]
                  + f_13 * fh_204[k];
        g_75[k] = g_35[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_0, fg_3, fg_5, fg_10, fg_12, fg_14, fg_45, \
                         fg_48, fg_50, fg_55, fg_57, fg_59, fg_75, fg_78, fg_80, fg_85, fg_87, \
                         fg_89, fg_90, fg_93, fg_95, fg_100, fg_102, fg_104, fg_120, fg_123, \
                         fg_125, fg_130, fg_132, fg_134, fg_135, fg_138, fg_140, fg_145, \
                         fg_147, fg_149, fh_0, fh_3, fh_5, fh_10, fh_12, fh_14, fh_63, fh_66, \
                         fh_68, fh_73, fh_75, fh_77, fh_105, fh_108, fh_110, fh_115, fh_117, \
                         fh_119, fh_127, fh_132, fh_134, fh_141, fh_143, fh_145, fh_169, \
                         fh_174, fh_176, fh_183, fh_185, fh_187, fh_191, fh_196, fh_198, \
                         fh_205, fh_207, fh_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -0.140625 * ab_x[k] * fg_0[k]
                  - 0.28125 * ab_x[k] * fg_3[k]
                  + 1.125 * ab_x[k] * fg_5[k]
                  - 0.140625 * ab_x[k] * fg_10[k]
                  + 1.125 * ab_x[k] * fg_12[k]
                  - 0.375 * ab_x[k] * fg_14[k]
                  - 0.28125 * ab_x[k] * fg_45[k]
                  - 0.5625 * ab_x[k] * fg_48[k]
                  + 2.25 * ab_x[k] * fg_50[k]
                  - 0.28125 * ab_x[k] * fg_55[k]
                  + 2.25 * ab_x[k] * fg_57[k]
                  - 0.75 * ab_x[k] * fg_59[k]
                  + 1.125 * ab_x[k] * fg_75[k]
                  + 2.25 * ab_x[k] * fg_78[k]
                  - 9.0 * ab_x[k] * fg_80[k]
                  + 1.125 * ab_x[k] * fg_85[k]
                  - 9.0 * ab_x[k] * fg_87[k]
                  + 3.0 * ab_x[k] * fg_89[k]
                  - 0.140625 * ab_y[k] * fg_90[k]
                  - 0.28125 * ab_y[k] * fg_93[k]
                  + 1.125 * ab_y[k] * fg_95[k]
                  - 0.140625 * ab_y[k] * fg_100[k]
                  + 1.125 * ab_y[k] * fg_102[k]
                  - 0.375 * ab_y[k] * fg_104[k]
                  + 1.125 * ab_y[k] * fg_120[k]
                  + 2.25 * ab_y[k] * fg_123[k]
                  - 9.0 * ab_y[k] * fg_125[k]
                  + 1.125 * ab_y[k] * fg_130[k]
                  - 9.0 * ab_y[k] * fg_132[k]
                  + 3.0 * ab_y[k] * fg_134[k]
                  - 0.375 * ab_z[k] * fg_135[k]
                  - 0.75 * ab_z[k] * fg_138[k]
                  + 3.0 * ab_z[k] * fg_140[k]
                  - 0.375 * ab_z[k] * fg_145[k]
                  + 3.0 * ab_z[k] * fg_147[k]
                  - ab_z[k] * fg_149[k]
                  + 0.140625 * fh_0[k]
                  + 0.28125 * fh_3[k]
                  - 1.125 * fh_5[k]
                  + 0.140625 * fh_10[k]
                  - 1.125 * fh_12[k]
                  + 0.375 * fh_14[k]
                  + 0.28125 * fh_63[k]
                  + 0.5625 * fh_66[k]
                  - 2.25 * fh_68[k]
                  + 0.28125 * fh_73[k]
                  - 2.25 * fh_75[k]
                  + 0.75 * fh_77[k]
                  - 1.125 * fh_105[k]
                  - 2.25 * fh_108[k]
                  + 9.0 * fh_110[k]
                  - 1.125 * fh_115[k]
                  + 9.0 * fh_117[k]
                  - 3.0 * fh_119[k]
                  + 0.140625 * fh_127[k]
                  + 0.28125 * fh_132[k]
                  - 1.125 * fh_134[k]
                  + 0.140625 * fh_141[k]
                  - 1.125 * fh_143[k]
                  + 0.375 * fh_145[k]
                  - 1.125 * fh_169[k]
                  - 2.25 * fh_174[k]
                  + 9.0 * fh_176[k]
                  - 1.125 * fh_183[k]
                  + 9.0 * fh_185[k]
                  - 3.0 * fh_187[k]
                  + 0.375 * fh_191[k]
                  + 0.75 * fh_196[k]
                  - 3.0 * fh_198[k]
                  + 0.375 * fh_205[k]
                  - 3.0 * fh_207[k]
                  + fh_209[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_2, fg_7, fg_9, fg_47, fg_52, fg_54, fg_77, \
                         fg_82, fg_84, fg_92, fg_97, fg_99, fg_122, fg_127, fg_129, fg_137, \
                         fg_142, fg_144, fh_2, fh_7, fh_9, fh_65, fh_70, fh_72, fh_107, \
                         fh_112, fh_114, fh_130, fh_137, fh_139, fh_172, fh_179, fh_181, \
                         fh_194, fh_201, fh_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_46 * ab_x[k] * fg_2[k]
                  + f_46 * ab_x[k] * fg_7[k]
                  - f_50 * ab_x[k] * fg_9[k]
                  + f_47 * ab_x[k] * fg_47[k]
                  + f_47 * ab_x[k] * fg_52[k]
                  - f_49 * ab_x[k] * fg_54[k]
                  - f_48 * ab_x[k] * fg_77[k]
                  - f_48 * ab_x[k] * fg_82[k]
                  + f_51 * ab_x[k] * fg_84[k]
                  + f_46 * ab_y[k] * fg_92[k]
                  + f_46 * ab_y[k] * fg_97[k]
                  - f_50 * ab_y[k] * fg_99[k]
                  - f_48 * ab_y[k] * fg_122[k]
                  - f_48 * ab_y[k] * fg_127[k]
                  + f_51 * ab_y[k] * fg_129[k]
                  + f_49 * ab_z[k] * fg_137[k]
                  + f_49 * ab_z[k] * fg_142[k]
                  - f_52 * ab_z[k] * fg_144[k]
                  - f_46 * fh_2[k]
                  - f_46 * fh_7[k]
                  + f_50 * fh_9[k]
                  - f_47 * fh_65[k]
                  - f_47 * fh_70[k]
                  + f_49 * fh_72[k]
                  + f_48 * fh_107[k]
                  + f_48 * fh_112[k]
                  - f_51 * fh_114[k]
                  - f_46 * fh_130[k]
                  - f_46 * fh_137[k]
                  + f_50 * fh_139[k]
                  + f_48 * fh_172[k]
                  + f_48 * fh_179[k]
                  - f_51 * fh_181[k]
                  - f_49 * fh_194[k]
                  - f_49 * fh_201[k]
                  + f_52 * fh_203[k];
        g_49[k] = g_41[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_0, fg_5, fg_10, fg_12, fg_45, fg_50, fg_55, \
                         fg_57, fg_75, fg_80, fg_85, fg_87, fg_90, fg_95, fg_100, fg_102, \
                         fg_120, fg_125, fg_130, fg_132, fg_135, fg_140, fg_145, fg_147, fh_0, \
                         fh_5, fh_10, fh_12, fh_63, fh_68, fh_73, fh_75, fh_105, fh_110, \
                         fh_115, fh_117, fh_127, fh_134, fh_141, fh_143, fh_169, fh_176, \
                         fh_183, fh_185, fh_191, fh_198, fh_205, \
                         fh_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_59 * ab_x[k] * fg_0[k]
                  - f_60 * ab_x[k] * fg_5[k]
                  - f_59 * ab_x[k] * fg_10[k]
                  + f_60 * ab_x[k] * fg_12[k]
                  + f_36 * ab_x[k] * fg_45[k]
                  - f_40 * ab_x[k] * fg_50[k]
                  - f_36 * ab_x[k] * fg_55[k]
                  + f_40 * ab_x[k] * fg_57[k]
                  - f_61 * ab_x[k] * fg_75[k]
                  + f_62 * ab_x[k] * fg_80[k]
                  + f_61 * ab_x[k] * fg_85[k]
                  - f_62 * ab_x[k] * fg_87[k]
                  + f_59 * ab_y[k] * fg_90[k]
                  - f_60 * ab_y[k] * fg_95[k]
                  - f_59 * ab_y[k] * fg_100[k]
                  + f_60 * ab_y[k] * fg_102[k]
                  - f_61 * ab_y[k] * fg_120[k]
                  + f_62 * ab_y[k] * fg_125[k]
                  + f_61 * ab_y[k] * fg_130[k]
                  - f_62 * ab_y[k] * fg_132[k]
                  + f_63 * ab_z[k] * fg_135[k]
                  - f_38 * ab_z[k] * fg_140[k]
                  - f_63 * ab_z[k] * fg_145[k]
                  + f_38 * ab_z[k] * fg_147[k]
                  - f_59 * fh_0[k]
                  + f_60 * fh_5[k]
                  + f_59 * fh_10[k]
                  - f_60 * fh_12[k]
                  - f_36 * fh_63[k]
                  + f_40 * fh_68[k]
                  + f_36 * fh_73[k]
                  - f_40 * fh_75[k]
                  + f_61 * fh_105[k]
                  - f_62 * fh_110[k]
                  - f_61 * fh_115[k]
                  + f_62 * fh_117[k]
                  - f_59 * fh_127[k]
                  + f_60 * fh_134[k]
                  + f_59 * fh_141[k]
                  - f_60 * fh_143[k]
                  + f_61 * fh_169[k]
                  - f_62 * fh_176[k]
                  - f_61 * fh_183[k]
                  + f_62 * fh_185[k]
                  - f_63 * fh_191[k]
                  + f_38 * fh_198[k]
                  + f_63 * fh_205[k]
                  - f_38 * fh_207[k];
        g_58[k] = g_42[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_2, fg_7, fg_47, fg_52, fg_77, fg_82, fg_92, \
                         fg_97, fg_122, fg_127, fg_137, fg_142, fh_2, fh_7, fh_65, fh_70, \
                         fh_107, fh_112, fh_130, fh_137, fh_172, fh_179, fh_194, \
                         fh_201 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_22 * ab_x[k] * fg_2[k]
                  + f_18 * ab_x[k] * fg_7[k]
                  - f_23 * ab_x[k] * fg_47[k]
                  + f_19 * ab_x[k] * fg_52[k]
                  + f_21 * ab_x[k] * fg_77[k]
                  - f_20 * ab_x[k] * fg_82[k]
                  - f_22 * ab_y[k] * fg_92[k]
                  + f_18 * ab_y[k] * fg_97[k]
                  + f_21 * ab_y[k] * fg_122[k]
                  - f_20 * ab_y[k] * fg_127[k]
                  - f_24 * ab_z[k] * fg_137[k]
                  + f_21 * ab_z[k] * fg_142[k]
                  + f_22 * fh_2[k]
                  - f_18 * fh_7[k]
                  + f_23 * fh_65[k]
                  - f_19 * fh_70[k]
                  - f_21 * fh_107[k]
                  + f_20 * fh_112[k]
                  + f_22 * fh_130[k]
                  - f_18 * fh_137[k]
                  - f_21 * fh_172[k]
                  + f_20 * fh_179[k]
                  + f_24 * fh_194[k]
                  - f_21 * fh_201[k];
        g_67[k] = g_43[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, fg_0, fg_3, fg_10, fg_45, fg_48, fg_55, fg_75, \
                         fg_78, fg_85, fg_90, fg_93, fg_100, fg_120, fg_123, fg_130, fg_135, \
                         fg_138, fg_145, fh_0, fh_3, fh_10, fh_63, fh_66, fh_73, fh_105, \
                         fh_108, fh_115, fh_127, fh_132, fh_141, fh_169, fh_174, fh_183, \
                         fh_191, fh_196, fh_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_64 * ab_x[k] * fg_0[k]
                  + f_65 * ab_x[k] * fg_3[k]
                  - f_64 * ab_x[k] * fg_10[k]
                  - f_66 * ab_x[k] * fg_45[k]
                  + f_67 * ab_x[k] * fg_48[k]
                  - f_66 * ab_x[k] * fg_55[k]
                  + f_7 * ab_x[k] * fg_75[k]
                  - f_68 * ab_x[k] * fg_78[k]
                  + f_7 * ab_x[k] * fg_85[k]
                  - f_64 * ab_y[k] * fg_90[k]
                  + f_65 * ab_y[k] * fg_93[k]
                  - f_64 * ab_y[k] * fg_100[k]
                  + f_7 * ab_y[k] * fg_120[k]
                  - f_68 * ab_y[k] * fg_123[k]
                  + f_7 * ab_y[k] * fg_130[k]
                  - f_69 * ab_z[k] * fg_135[k]
                  + f_70 * ab_z[k] * fg_138[k]
                  - f_69 * ab_z[k] * fg_145[k]
                  + f_64 * fh_0[k]
                  - f_65 * fh_3[k]
                  + f_64 * fh_10[k]
                  + f_66 * fh_63[k]
                  - f_67 * fh_66[k]
                  + f_66 * fh_73[k]
                  - f_7 * fh_105[k]
                  + f_68 * fh_108[k]
                  - f_7 * fh_115[k]
                  + f_64 * fh_127[k]
                  - f_65 * fh_132[k]
                  + f_64 * fh_141[k]
                  - f_7 * fh_169[k]
                  + f_68 * fh_174[k]
                  - f_7 * fh_183[k]
                  + f_69 * fh_191[k]
                  - f_70 * fh_196[k]
                  + f_69 * fh_205[k];
        g_76[k] = g_44[k];
    }

#pragma omp simd aligned(ab_x, fg_32, fg_37, fg_39, fg_107, fg_112, fg_114, fg_137, fg_142, \
                         fg_144, fh_44, fh_49, fh_51, fh_149, fh_154, fh_156, fh_191, fh_196, \
                         fh_198 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -5.625 * ab_x[k] * fg_32[k]
                  - 5.625 * ab_x[k] * fg_37[k]
                  + 7.5 * ab_x[k] * fg_39[k]
                  - 5.625 * ab_x[k] * fg_107[k]
                  - 5.625 * ab_x[k] * fg_112[k]
                  + 7.5 * ab_x[k] * fg_114[k]
                  + 7.5 * ab_x[k] * fg_137[k]
                  + 7.5 * ab_x[k] * fg_142[k]
                  - 10.0 * ab_x[k] * fg_144[k]
                  + 5.625 * fh_44[k]
                  + 5.625 * fh_49[k]
                  - 7.5 * fh_51[k]
                  + 5.625 * fh_149[k]
                  + 5.625 * fh_154[k]
                  - 7.5 * fh_156[k]
                  - 7.5 * fh_191[k]
                  - 7.5 * fh_196[k]
                  + 10.0 * fh_198[k];
    }

#pragma omp simd aligned(ab_x, fg_30, fg_35, fg_40, fg_42, fg_105, fg_110, fg_115, fg_117, \
                         fg_135, fg_140, fg_145, fg_147, fh_42, fh_47, fh_52, fh_54, fh_147, \
                         fh_152, fh_157, fh_159, fh_189, fh_194, fh_199, \
                         fh_201 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_53 * ab_x[k] * fg_30[k]
                  + f_54 * ab_x[k] * fg_35[k]
                  + f_53 * ab_x[k] * fg_40[k]
                  - f_54 * ab_x[k] * fg_42[k]
                  - f_53 * ab_x[k] * fg_105[k]
                  + f_54 * ab_x[k] * fg_110[k]
                  + f_53 * ab_x[k] * fg_115[k]
                  - f_54 * ab_x[k] * fg_117[k]
                  + f_55 * ab_x[k] * fg_135[k]
                  - f_56 * ab_x[k] * fg_140[k]
                  - f_55 * ab_x[k] * fg_145[k]
                  + f_56 * ab_x[k] * fg_147[k]
                  + f_53 * fh_42[k]
                  - f_54 * fh_47[k]
                  - f_53 * fh_52[k]
                  + f_54 * fh_54[k]
                  + f_53 * fh_147[k]
                  - f_54 * fh_152[k]
                  - f_53 * fh_157[k]
                  + f_54 * fh_159[k]
                  - f_55 * fh_189[k]
                  + f_56 * fh_194[k]
                  + f_55 * fh_199[k]
                  - f_56 * fh_201[k];
        g_59[k] = g_51[k];
    }

#pragma omp simd aligned(ab_x, fg_32, fg_37, fg_107, fg_112, fg_137, fg_142, fh_44, fh_49, \
                         fh_149, fh_154, fh_191, fh_196 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_16 * ab_x[k] * fg_32[k]
                  - f_15 * ab_x[k] * fg_37[k]
                  + f_16 * ab_x[k] * fg_107[k]
                  - f_15 * ab_x[k] * fg_112[k]
                  - f_17 * ab_x[k] * fg_137[k]
                  + f_3 * ab_x[k] * fg_142[k]
                  - f_16 * fh_44[k]
                  + f_15 * fh_49[k]
                  - f_16 * fh_149[k]
                  + f_15 * fh_154[k]
                  + f_17 * fh_191[k]
                  - f_3 * fh_196[k];
        g_68[k] = g_52[k];
    }

#pragma omp simd aligned(ab_x, fg_30, fg_33, fg_40, fg_105, fg_108, fg_115, fg_135, fg_138, \
                         fg_145, fh_42, fh_45, fh_52, fh_147, fh_150, fh_157, fh_189, fh_192, \
                         fh_199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_57 * ab_x[k] * fg_30[k]
                  - f_58 * ab_x[k] * fg_33[k]
                  + f_57 * ab_x[k] * fg_40[k]
                  + f_57 * ab_x[k] * fg_105[k]
                  - f_58 * ab_x[k] * fg_108[k]
                  + f_57 * ab_x[k] * fg_115[k]
                  - f_13 * ab_x[k] * fg_135[k]
                  + f_14 * ab_x[k] * fg_138[k]
                  - f_13 * ab_x[k] * fg_145[k]
                  - f_57 * fh_42[k]
                  + f_58 * fh_45[k]
                  - f_57 * fh_52[k]
                  - f_57 * fh_147[k]
                  + f_58 * fh_150[k]
                  - f_57 * fh_157[k]
                  + f_13 * fh_189[k]
                  - f_14 * fh_192[k]
                  + f_13 * fh_199[k];
        g_77[k] = g_53[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_0, fg_5, fg_10, fg_12, fg_75, fg_80, fg_85, fg_87, \
                         fg_90, fg_95, fg_100, fg_102, fg_120, fg_125, fg_130, fg_132, fh_0, \
                         fh_5, fh_10, fh_12, fh_105, fh_110, fh_115, fh_117, fh_127, fh_134, \
                         fh_141, fh_143, fh_169, fh_176, fh_183, \
                         fh_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -0.3125 * ab_x[k] * fg_0[k]
                  + 1.875 * ab_x[k] * fg_5[k]
                  + 0.3125 * ab_x[k] * fg_10[k]
                  - 1.875 * ab_x[k] * fg_12[k]
                  + 1.875 * ab_x[k] * fg_75[k]
                  - 11.25 * ab_x[k] * fg_80[k]
                  - 1.875 * ab_x[k] * fg_85[k]
                  + 11.25 * ab_x[k] * fg_87[k]
                  + 0.3125 * ab_y[k] * fg_90[k]
                  - 1.875 * ab_y[k] * fg_95[k]
                  - 0.3125 * ab_y[k] * fg_100[k]
                  + 1.875 * ab_y[k] * fg_102[k]
                  - 1.875 * ab_y[k] * fg_120[k]
                  + 11.25 * ab_y[k] * fg_125[k]
                  + 1.875 * ab_y[k] * fg_130[k]
                  - 11.25 * ab_y[k] * fg_132[k]
                  + 0.3125 * fh_0[k]
                  - 1.875 * fh_5[k]
                  - 0.3125 * fh_10[k]
                  + 1.875 * fh_12[k]
                  - 1.875 * fh_105[k]
                  + 11.25 * fh_110[k]
                  + 1.875 * fh_115[k]
                  - 11.25 * fh_117[k]
                  - 0.3125 * fh_127[k]
                  + 1.875 * fh_134[k]
                  + 0.3125 * fh_141[k]
                  - 1.875 * fh_143[k]
                  + 1.875 * fh_169[k]
                  - 11.25 * fh_176[k]
                  - 1.875 * fh_183[k]
                  + 11.25 * fh_185[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_2, fg_7, fg_77, fg_82, fg_92, fg_97, fg_122, fg_127, \
                         fh_2, fh_7, fh_107, fh_112, fh_130, fh_137, fh_172, \
                         fh_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_27 * ab_x[k] * fg_2[k]
                  - f_25 * ab_x[k] * fg_7[k]
                  - f_4 * ab_x[k] * fg_77[k]
                  + f_26 * ab_x[k] * fg_82[k]
                  - f_27 * ab_y[k] * fg_92[k]
                  + f_25 * ab_y[k] * fg_97[k]
                  + f_4 * ab_y[k] * fg_122[k]
                  - f_26 * ab_y[k] * fg_127[k]
                  - f_27 * fh_2[k]
                  + f_25 * fh_7[k]
                  + f_4 * fh_107[k]
                  - f_26 * fh_112[k]
                  + f_27 * fh_130[k]
                  - f_25 * fh_137[k]
                  - f_4 * fh_172[k]
                  + f_26 * fh_179[k];
        g_69[k] = g_61[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_0, fg_3, fg_10, fg_75, fg_78, fg_85, fg_90, fg_93, \
                         fg_100, fg_120, fg_123, fg_130, fh_0, fh_3, fh_10, fh_105, fh_108, \
                         fh_115, fh_127, fh_132, fh_141, fh_169, fh_174, \
                         fh_183 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_71 * ab_x[k] * fg_0[k]
                  - f_72 * ab_x[k] * fg_3[k]
                  + f_71 * ab_x[k] * fg_10[k]
                  - f_72 * ab_x[k] * fg_75[k]
                  + f_15 * ab_x[k] * fg_78[k]
                  - f_72 * ab_x[k] * fg_85[k]
                  - f_71 * ab_y[k] * fg_90[k]
                  + f_72 * ab_y[k] * fg_93[k]
                  - f_71 * ab_y[k] * fg_100[k]
                  + f_72 * ab_y[k] * fg_120[k]
                  - f_15 * ab_y[k] * fg_123[k]
                  + f_72 * ab_y[k] * fg_130[k]
                  - f_71 * fh_0[k]
                  + f_72 * fh_3[k]
                  - f_71 * fh_10[k]
                  + f_72 * fh_105[k]
                  - f_15 * fh_108[k]
                  + f_72 * fh_115[k]
                  + f_71 * fh_127[k]
                  - f_72 * fh_132[k]
                  + f_71 * fh_141[k]
                  - f_72 * fh_169[k]
                  + f_15 * fh_174[k]
                  - f_72 * fh_183[k];
        g_78[k] = g_62[k];
    }

#pragma omp simd aligned(ab_x, fg_32, fg_37, fg_107, fg_112, fh_44, fh_49, fh_149, \
                         fh_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -4.375 * ab_x[k] * fg_32[k]
                  + 13.125 * ab_x[k] * fg_37[k]
                  + 13.125 * ab_x[k] * fg_107[k]
                  - 39.375 * ab_x[k] * fg_112[k]
                  + 4.375 * fh_44[k]
                  - 13.125 * fh_49[k]
                  - 13.125 * fh_149[k]
                  + 39.375 * fh_154[k];
    }

#pragma omp simd aligned(ab_x, fg_30, fg_33, fg_40, fg_105, fg_108, fg_115, fh_42, fh_45, \
                         fh_52, fh_147, fh_150, fh_157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_30 * ab_x[k] * fg_30[k]
                  + f_31 * ab_x[k] * fg_33[k]
                  - f_30 * ab_x[k] * fg_40[k]
                  + f_28 * ab_x[k] * fg_105[k]
                  - f_29 * ab_x[k] * fg_108[k]
                  + f_28 * ab_x[k] * fg_115[k]
                  + f_30 * fh_42[k]
                  - f_31 * fh_45[k]
                  + f_30 * fh_52[k]
                  - f_28 * fh_147[k]
                  + f_29 * fh_150[k]
                  - f_28 * fh_157[k];
        g_79[k] = g_71[k];
    }

#pragma omp simd aligned(ab_x, ab_y, fg_0, fg_3, fg_10, fg_45, fg_48, fg_55, fg_90, fg_93, \
                         fg_100, fh_0, fh_3, fh_10, fh_63, fh_66, fh_73, fh_127, fh_132, \
                         fh_141 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -0.546875 * ab_x[k] * fg_0[k]
                  + 3.28125 * ab_x[k] * fg_3[k]
                  - 0.546875 * ab_x[k] * fg_10[k]
                  + 3.28125 * ab_x[k] * fg_45[k]
                  - 19.6875 * ab_x[k] * fg_48[k]
                  + 3.28125 * ab_x[k] * fg_55[k]
                  - 0.546875 * ab_y[k] * fg_90[k]
                  + 3.28125 * ab_y[k] * fg_93[k]
                  - 0.546875 * ab_y[k] * fg_100[k]
                  + 0.546875 * fh_0[k]
                  - 3.28125 * fh_3[k]
                  + 0.546875 * fh_10[k]
                  - 3.28125 * fh_63[k]
                  + 19.6875 * fh_66[k]
                  - 3.28125 * fh_73[k]
                  + 0.546875 * fh_127[k]
                  - 3.28125 * fh_132[k]
                  + 0.546875 * fh_141[k];
    }
}

}  // namespace simdovl
