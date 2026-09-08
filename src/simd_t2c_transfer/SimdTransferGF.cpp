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

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_gf(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t gd, const size_t hd, const size_t nmax) -> void
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, gd_0, gd_1, gd_2, gd_3, gd_4, hd_0, \
                         hd_1, hd_2, hd_3, hd_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * gd_0[k]
                 + hd_0[k];

        t_1[k] = ab_x[k] * gd_1[k]
                 + hd_1[k];

        t_2[k] = ab_x[k] * gd_2[k]
                 + hd_2[k];

        t_3[k] = ab_x[k] * gd_3[k]
                 + hd_3[k];

        t_4[k] = ab_x[k] * gd_4[k]
                 + hd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ab_y, ab_z, gd_3, gd_4, gd_5, hd_5, \
                         hd_9, hd_10, hd_11, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = ab_x[k] * gd_5[k]
                 + hd_5[k];

        t_6[k] = ab_y[k] * gd_3[k]
                 + hd_9[k];

        t_7[k] = ab_y[k] * gd_4[k]
                 + hd_10[k];

        t_8[k] = ab_y[k] * gd_5[k]
                 + hd_11[k];

        t_9[k] = ab_z[k] * gd_5[k]
                 + hd_17[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, gd_6, gd_7, gd_8, gd_9, gd_10, \
                         hd_6, hd_7, hd_8, hd_9, hd_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = ab_x[k] * gd_6[k]
                  + hd_6[k];

        t_11[k] = ab_x[k] * gd_7[k]
                  + hd_7[k];

        t_12[k] = ab_x[k] * gd_8[k]
                  + hd_8[k];

        t_13[k] = ab_x[k] * gd_9[k]
                  + hd_9[k];

        t_14[k] = ab_x[k] * gd_10[k]
                  + hd_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, gd_9, gd_10, gd_11, \
                         hd_11, hd_21, hd_22, hd_23, hd_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = ab_x[k] * gd_11[k]
                  + hd_11[k];

        t_16[k] = ab_y[k] * gd_9[k]
                  + hd_21[k];

        t_17[k] = ab_y[k] * gd_10[k]
                  + hd_22[k];

        t_18[k] = ab_y[k] * gd_11[k]
                  + hd_23[k];

        t_19[k] = ab_z[k] * gd_11[k]
                  + hd_29[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, gd_12, gd_13, gd_14, gd_15, \
                         gd_16, hd_12, hd_13, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = ab_x[k] * gd_12[k]
                  + hd_12[k];

        t_21[k] = ab_x[k] * gd_13[k]
                  + hd_13[k];

        t_22[k] = ab_x[k] * gd_14[k]
                  + hd_14[k];

        t_23[k] = ab_x[k] * gd_15[k]
                  + hd_15[k];

        t_24[k] = ab_x[k] * gd_16[k]
                  + hd_16[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ab_y, ab_z, gd_15, gd_16, gd_17, \
                         hd_17, hd_27, hd_28, hd_29, hd_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = ab_x[k] * gd_17[k]
                  + hd_17[k];

        t_26[k] = ab_y[k] * gd_15[k]
                  + hd_27[k];

        t_27[k] = ab_y[k] * gd_16[k]
                  + hd_28[k];

        t_28[k] = ab_y[k] * gd_17[k]
                  + hd_29[k];

        t_29[k] = ab_z[k] * gd_17[k]
                  + hd_35[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, gd_18, gd_19, gd_20, gd_21, \
                         gd_22, hd_18, hd_19, hd_20, hd_21, hd_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = ab_x[k] * gd_18[k]
                  + hd_18[k];

        t_31[k] = ab_x[k] * gd_19[k]
                  + hd_19[k];

        t_32[k] = ab_x[k] * gd_20[k]
                  + hd_20[k];

        t_33[k] = ab_x[k] * gd_21[k]
                  + hd_21[k];

        t_34[k] = ab_x[k] * gd_22[k]
                  + hd_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, ab_y, ab_z, gd_21, gd_22, gd_23, \
                         hd_23, hd_39, hd_40, hd_41, hd_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = ab_x[k] * gd_23[k]
                  + hd_23[k];

        t_36[k] = ab_y[k] * gd_21[k]
                  + hd_39[k];

        t_37[k] = ab_y[k] * gd_22[k]
                  + hd_40[k];

        t_38[k] = ab_y[k] * gd_23[k]
                  + hd_41[k];

        t_39[k] = ab_z[k] * gd_23[k]
                  + hd_47[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, gd_24, gd_25, gd_26, gd_27, \
                         gd_28, hd_24, hd_25, hd_26, hd_27, hd_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = ab_x[k] * gd_24[k]
                  + hd_24[k];

        t_41[k] = ab_x[k] * gd_25[k]
                  + hd_25[k];

        t_42[k] = ab_x[k] * gd_26[k]
                  + hd_26[k];

        t_43[k] = ab_x[k] * gd_27[k]
                  + hd_27[k];

        t_44[k] = ab_x[k] * gd_28[k]
                  + hd_28[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ab_y, ab_z, gd_27, gd_28, gd_29, \
                         hd_29, hd_45, hd_46, hd_47, hd_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = ab_x[k] * gd_29[k]
                  + hd_29[k];

        t_46[k] = ab_y[k] * gd_27[k]
                  + hd_45[k];

        t_47[k] = ab_y[k] * gd_28[k]
                  + hd_46[k];

        t_48[k] = ab_y[k] * gd_29[k]
                  + hd_47[k];

        t_49[k] = ab_z[k] * gd_29[k]
                  + hd_53[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, gd_30, gd_31, gd_32, gd_33, \
                         gd_34, hd_30, hd_31, hd_32, hd_33, hd_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = ab_x[k] * gd_30[k]
                  + hd_30[k];

        t_51[k] = ab_x[k] * gd_31[k]
                  + hd_31[k];

        t_52[k] = ab_x[k] * gd_32[k]
                  + hd_32[k];

        t_53[k] = ab_x[k] * gd_33[k]
                  + hd_33[k];

        t_54[k] = ab_x[k] * gd_34[k]
                  + hd_34[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, gd_33, gd_34, gd_35, \
                         hd_35, hd_51, hd_52, hd_53, hd_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = ab_x[k] * gd_35[k]
                  + hd_35[k];

        t_56[k] = ab_y[k] * gd_33[k]
                  + hd_51[k];

        t_57[k] = ab_y[k] * gd_34[k]
                  + hd_52[k];

        t_58[k] = ab_y[k] * gd_35[k]
                  + hd_53[k];

        t_59[k] = ab_z[k] * gd_35[k]
                  + hd_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, gd_36, gd_37, gd_38, gd_39, \
                         gd_40, hd_36, hd_37, hd_38, hd_39, hd_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = ab_x[k] * gd_36[k]
                  + hd_36[k];

        t_61[k] = ab_x[k] * gd_37[k]
                  + hd_37[k];

        t_62[k] = ab_x[k] * gd_38[k]
                  + hd_38[k];

        t_63[k] = ab_x[k] * gd_39[k]
                  + hd_39[k];

        t_64[k] = ab_x[k] * gd_40[k]
                  + hd_40[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, ab_y, ab_z, gd_39, gd_40, gd_41, \
                         hd_41, hd_63, hd_64, hd_65, hd_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = ab_x[k] * gd_41[k]
                  + hd_41[k];

        t_66[k] = ab_y[k] * gd_39[k]
                  + hd_63[k];

        t_67[k] = ab_y[k] * gd_40[k]
                  + hd_64[k];

        t_68[k] = ab_y[k] * gd_41[k]
                  + hd_65[k];

        t_69[k] = ab_z[k] * gd_41[k]
                  + hd_71[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, gd_42, gd_43, gd_44, gd_45, \
                         gd_46, hd_42, hd_43, hd_44, hd_45, hd_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = ab_x[k] * gd_42[k]
                  + hd_42[k];

        t_71[k] = ab_x[k] * gd_43[k]
                  + hd_43[k];

        t_72[k] = ab_x[k] * gd_44[k]
                  + hd_44[k];

        t_73[k] = ab_x[k] * gd_45[k]
                  + hd_45[k];

        t_74[k] = ab_x[k] * gd_46[k]
                  + hd_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, ab_y, ab_z, gd_45, gd_46, gd_47, \
                         hd_47, hd_69, hd_70, hd_71, hd_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = ab_x[k] * gd_47[k]
                  + hd_47[k];

        t_76[k] = ab_y[k] * gd_45[k]
                  + hd_69[k];

        t_77[k] = ab_y[k] * gd_46[k]
                  + hd_70[k];

        t_78[k] = ab_y[k] * gd_47[k]
                  + hd_71[k];

        t_79[k] = ab_z[k] * gd_47[k]
                  + hd_77[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, gd_48, gd_49, gd_50, gd_51, \
                         gd_52, hd_48, hd_49, hd_50, hd_51, hd_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = ab_x[k] * gd_48[k]
                  + hd_48[k];

        t_81[k] = ab_x[k] * gd_49[k]
                  + hd_49[k];

        t_82[k] = ab_x[k] * gd_50[k]
                  + hd_50[k];

        t_83[k] = ab_x[k] * gd_51[k]
                  + hd_51[k];

        t_84[k] = ab_x[k] * gd_52[k]
                  + hd_52[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, gd_51, gd_52, gd_53, \
                         hd_53, hd_75, hd_76, hd_77, hd_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = ab_x[k] * gd_53[k]
                  + hd_53[k];

        t_86[k] = ab_y[k] * gd_51[k]
                  + hd_75[k];

        t_87[k] = ab_y[k] * gd_52[k]
                  + hd_76[k];

        t_88[k] = ab_y[k] * gd_53[k]
                  + hd_77[k];

        t_89[k] = ab_z[k] * gd_53[k]
                  + hd_83[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, gd_54, gd_55, gd_56, gd_57, \
                         gd_58, hd_54, hd_55, hd_56, hd_57, hd_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = ab_x[k] * gd_54[k]
                  + hd_54[k];

        t_91[k] = ab_x[k] * gd_55[k]
                  + hd_55[k];

        t_92[k] = ab_x[k] * gd_56[k]
                  + hd_56[k];

        t_93[k] = ab_x[k] * gd_57[k]
                  + hd_57[k];

        t_94[k] = ab_x[k] * gd_58[k]
                  + hd_58[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, ab_y, ab_z, gd_57, gd_58, gd_59, \
                         hd_59, hd_81, hd_82, hd_83, hd_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = ab_x[k] * gd_59[k]
                  + hd_59[k];

        t_96[k] = ab_y[k] * gd_57[k]
                  + hd_81[k];

        t_97[k] = ab_y[k] * gd_58[k]
                  + hd_82[k];

        t_98[k] = ab_y[k] * gd_59[k]
                  + hd_83[k];

        t_99[k] = ab_z[k] * gd_59[k]
                  + hd_89[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, gd_60, gd_61, gd_62, gd_63, \
                         gd_64, hd_60, hd_61, hd_62, hd_63, hd_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = ab_x[k] * gd_60[k]
                   + hd_60[k];

        t_101[k] = ab_x[k] * gd_61[k]
                   + hd_61[k];

        t_102[k] = ab_x[k] * gd_62[k]
                   + hd_62[k];

        t_103[k] = ab_x[k] * gd_63[k]
                   + hd_63[k];

        t_104[k] = ab_x[k] * gd_64[k]
                   + hd_64[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, gd_63, gd_64, \
                         gd_65, hd_65, hd_93, hd_94, hd_95, hd_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = ab_x[k] * gd_65[k]
                   + hd_65[k];

        t_106[k] = ab_y[k] * gd_63[k]
                   + hd_93[k];

        t_107[k] = ab_y[k] * gd_64[k]
                   + hd_94[k];

        t_108[k] = ab_y[k] * gd_65[k]
                   + hd_95[k];

        t_109[k] = ab_z[k] * gd_65[k]
                   + hd_101[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, gd_66, gd_67, gd_68, gd_69, \
                         gd_70, hd_66, hd_67, hd_68, hd_69, hd_70 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = ab_x[k] * gd_66[k]
                   + hd_66[k];

        t_111[k] = ab_x[k] * gd_67[k]
                   + hd_67[k];

        t_112[k] = ab_x[k] * gd_68[k]
                   + hd_68[k];

        t_113[k] = ab_x[k] * gd_69[k]
                   + hd_69[k];

        t_114[k] = ab_x[k] * gd_70[k]
                   + hd_70[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, ab_y, ab_z, gd_69, gd_70, \
                         gd_71, hd_71, hd_99, hd_100, hd_101, hd_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = ab_x[k] * gd_71[k]
                   + hd_71[k];

        t_116[k] = ab_y[k] * gd_69[k]
                   + hd_99[k];

        t_117[k] = ab_y[k] * gd_70[k]
                   + hd_100[k];

        t_118[k] = ab_y[k] * gd_71[k]
                   + hd_101[k];

        t_119[k] = ab_z[k] * gd_71[k]
                   + hd_107[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, gd_72, gd_73, gd_74, gd_75, \
                         gd_76, hd_72, hd_73, hd_74, hd_75, hd_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = ab_x[k] * gd_72[k]
                   + hd_72[k];

        t_121[k] = ab_x[k] * gd_73[k]
                   + hd_73[k];

        t_122[k] = ab_x[k] * gd_74[k]
                   + hd_74[k];

        t_123[k] = ab_x[k] * gd_75[k]
                   + hd_75[k];

        t_124[k] = ab_x[k] * gd_76[k]
                   + hd_76[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, ab_y, ab_z, gd_75, gd_76, \
                         gd_77, hd_77, hd_105, hd_106, hd_107, hd_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = ab_x[k] * gd_77[k]
                   + hd_77[k];

        t_126[k] = ab_y[k] * gd_75[k]
                   + hd_105[k];

        t_127[k] = ab_y[k] * gd_76[k]
                   + hd_106[k];

        t_128[k] = ab_y[k] * gd_77[k]
                   + hd_107[k];

        t_129[k] = ab_z[k] * gd_77[k]
                   + hd_113[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, gd_78, gd_79, gd_80, gd_81, \
                         gd_82, hd_78, hd_79, hd_80, hd_81, hd_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_130[k] = ab_x[k] * gd_78[k]
                   + hd_78[k];

        t_131[k] = ab_x[k] * gd_79[k]
                   + hd_79[k];

        t_132[k] = ab_x[k] * gd_80[k]
                   + hd_80[k];

        t_133[k] = ab_x[k] * gd_81[k]
                   + hd_81[k];

        t_134[k] = ab_x[k] * gd_82[k]
                   + hd_82[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, ab_y, ab_z, gd_81, gd_82, \
                         gd_83, hd_83, hd_111, hd_112, hd_113, hd_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = ab_x[k] * gd_83[k]
                   + hd_83[k];

        t_136[k] = ab_y[k] * gd_81[k]
                   + hd_111[k];

        t_137[k] = ab_y[k] * gd_82[k]
                   + hd_112[k];

        t_138[k] = ab_y[k] * gd_83[k]
                   + hd_113[k];

        t_139[k] = ab_z[k] * gd_83[k]
                   + hd_119[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, gd_84, gd_85, gd_86, gd_87, \
                         gd_88, hd_84, hd_85, hd_86, hd_87, hd_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_140[k] = ab_x[k] * gd_84[k]
                   + hd_84[k];

        t_141[k] = ab_x[k] * gd_85[k]
                   + hd_85[k];

        t_142[k] = ab_x[k] * gd_86[k]
                   + hd_86[k];

        t_143[k] = ab_x[k] * gd_87[k]
                   + hd_87[k];

        t_144[k] = ab_x[k] * gd_88[k]
                   + hd_88[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, gd_87, gd_88, \
                         gd_89, hd_89, hd_117, hd_118, hd_119, hd_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_145[k] = ab_x[k] * gd_89[k]
                   + hd_89[k];

        t_146[k] = ab_y[k] * gd_87[k]
                   + hd_117[k];

        t_147[k] = ab_y[k] * gd_88[k]
                   + hd_118[k];

        t_148[k] = ab_y[k] * gd_89[k]
                   + hd_119[k];

        t_149[k] = ab_z[k] * gd_89[k]
                   + hd_125[k];
    }
}

}  // namespace simdtrf
