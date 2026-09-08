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


#include "SimdTransferID.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_id(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t ip, const size_t kp, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_30 = buffer.data(ip + 30);
    const auto *ip_31 = buffer.data(ip + 31);
    const auto *ip_32 = buffer.data(ip + 32);
    const auto *ip_33 = buffer.data(ip + 33);
    const auto *ip_34 = buffer.data(ip + 34);
    const auto *ip_35 = buffer.data(ip + 35);
    const auto *ip_36 = buffer.data(ip + 36);
    const auto *ip_37 = buffer.data(ip + 37);
    const auto *ip_38 = buffer.data(ip + 38);
    const auto *ip_39 = buffer.data(ip + 39);
    const auto *ip_40 = buffer.data(ip + 40);
    const auto *ip_41 = buffer.data(ip + 41);
    const auto *ip_42 = buffer.data(ip + 42);
    const auto *ip_43 = buffer.data(ip + 43);
    const auto *ip_44 = buffer.data(ip + 44);
    const auto *ip_45 = buffer.data(ip + 45);
    const auto *ip_46 = buffer.data(ip + 46);
    const auto *ip_47 = buffer.data(ip + 47);
    const auto *ip_48 = buffer.data(ip + 48);
    const auto *ip_49 = buffer.data(ip + 49);
    const auto *ip_50 = buffer.data(ip + 50);
    const auto *ip_51 = buffer.data(ip + 51);
    const auto *ip_52 = buffer.data(ip + 52);
    const auto *ip_53 = buffer.data(ip + 53);
    const auto *ip_54 = buffer.data(ip + 54);
    const auto *ip_55 = buffer.data(ip + 55);
    const auto *ip_56 = buffer.data(ip + 56);
    const auto *ip_57 = buffer.data(ip + 57);
    const auto *ip_58 = buffer.data(ip + 58);
    const auto *ip_59 = buffer.data(ip + 59);
    const auto *ip_60 = buffer.data(ip + 60);
    const auto *ip_61 = buffer.data(ip + 61);
    const auto *ip_62 = buffer.data(ip + 62);
    const auto *ip_63 = buffer.data(ip + 63);
    const auto *ip_64 = buffer.data(ip + 64);
    const auto *ip_65 = buffer.data(ip + 65);
    const auto *ip_66 = buffer.data(ip + 66);
    const auto *ip_67 = buffer.data(ip + 67);
    const auto *ip_68 = buffer.data(ip + 68);
    const auto *ip_69 = buffer.data(ip + 69);
    const auto *ip_70 = buffer.data(ip + 70);
    const auto *ip_71 = buffer.data(ip + 71);
    const auto *ip_72 = buffer.data(ip + 72);
    const auto *ip_73 = buffer.data(ip + 73);
    const auto *ip_74 = buffer.data(ip + 74);
    const auto *ip_75 = buffer.data(ip + 75);
    const auto *ip_76 = buffer.data(ip + 76);
    const auto *ip_77 = buffer.data(ip + 77);
    const auto *ip_78 = buffer.data(ip + 78);
    const auto *ip_79 = buffer.data(ip + 79);
    const auto *ip_80 = buffer.data(ip + 80);
    const auto *ip_81 = buffer.data(ip + 81);
    const auto *ip_82 = buffer.data(ip + 82);
    const auto *ip_83 = buffer.data(ip + 83);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_34 = buffer.data(kp + 34);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_39 = buffer.data(kp + 39);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_48 = buffer.data(kp + 48);
    const auto *kp_49 = buffer.data(kp + 49);
    const auto *kp_50 = buffer.data(kp + 50);
    const auto *kp_51 = buffer.data(kp + 51);
    const auto *kp_52 = buffer.data(kp + 52);
    const auto *kp_53 = buffer.data(kp + 53);
    const auto *kp_54 = buffer.data(kp + 54);
    const auto *kp_55 = buffer.data(kp + 55);
    const auto *kp_56 = buffer.data(kp + 56);
    const auto *kp_57 = buffer.data(kp + 57);
    const auto *kp_58 = buffer.data(kp + 58);
    const auto *kp_59 = buffer.data(kp + 59);
    const auto *kp_60 = buffer.data(kp + 60);
    const auto *kp_61 = buffer.data(kp + 61);
    const auto *kp_62 = buffer.data(kp + 62);
    const auto *kp_63 = buffer.data(kp + 63);
    const auto *kp_64 = buffer.data(kp + 64);
    const auto *kp_65 = buffer.data(kp + 65);
    const auto *kp_66 = buffer.data(kp + 66);
    const auto *kp_67 = buffer.data(kp + 67);
    const auto *kp_68 = buffer.data(kp + 68);
    const auto *kp_69 = buffer.data(kp + 69);
    const auto *kp_70 = buffer.data(kp + 70);
    const auto *kp_71 = buffer.data(kp + 71);
    const auto *kp_72 = buffer.data(kp + 72);
    const auto *kp_73 = buffer.data(kp + 73);
    const auto *kp_74 = buffer.data(kp + 74);
    const auto *kp_75 = buffer.data(kp + 75);
    const auto *kp_76 = buffer.data(kp + 76);
    const auto *kp_77 = buffer.data(kp + 77);
    const auto *kp_78 = buffer.data(kp + 78);
    const auto *kp_79 = buffer.data(kp + 79);
    const auto *kp_80 = buffer.data(kp + 80);
    const auto *kp_81 = buffer.data(kp + 81);
    const auto *kp_82 = buffer.data(kp + 82);
    const auto *kp_83 = buffer.data(kp + 83);
    const auto *kp_85 = buffer.data(kp + 85);
    const auto *kp_86 = buffer.data(kp + 86);
    const auto *kp_88 = buffer.data(kp + 88);
    const auto *kp_89 = buffer.data(kp + 89);
    const auto *kp_91 = buffer.data(kp + 91);
    const auto *kp_92 = buffer.data(kp + 92);
    const auto *kp_94 = buffer.data(kp + 94);
    const auto *kp_95 = buffer.data(kp + 95);
    const auto *kp_97 = buffer.data(kp + 97);
    const auto *kp_98 = buffer.data(kp + 98);
    const auto *kp_100 = buffer.data(kp + 100);
    const auto *kp_101 = buffer.data(kp + 101);
    const auto *kp_103 = buffer.data(kp + 103);
    const auto *kp_104 = buffer.data(kp + 104);
    const auto *kp_107 = buffer.data(kp + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, ip_0, ip_1, ip_2, kp_0, kp_1, \
                         kp_2, kp_4, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * ip_0[k]
                 + kp_0[k];

        t_1[k] = ab_x[k] * ip_1[k]
                 + kp_1[k];

        t_2[k] = ab_x[k] * ip_2[k]
                 + kp_2[k];

        t_3[k] = ab_y[k] * ip_1[k]
                 + kp_4[k];

        t_4[k] = ab_y[k] * ip_2[k]
                 + kp_5[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, ip_2, ip_3, ip_4, ip_5, kp_3, kp_4, \
                         kp_5, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = ab_z[k] * ip_2[k]
                 + kp_8[k];

        t_6[k] = ab_x[k] * ip_3[k]
                 + kp_3[k];

        t_7[k] = ab_x[k] * ip_4[k]
                 + kp_4[k];

        t_8[k] = ab_x[k] * ip_5[k]
                 + kp_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, ip_4, ip_5, ip_6, kp_6, \
                         kp_10, kp_11, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_9[k] = ab_y[k] * ip_4[k]
                 + kp_10[k];

        t_10[k] = ab_y[k] * ip_5[k]
                  + kp_11[k];

        t_11[k] = ab_z[k] * ip_5[k]
                  + kp_14[k];

        t_12[k] = ab_x[k] * ip_6[k]
                  + kp_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, ip_7, ip_8, kp_7, \
                         kp_8, kp_13, kp_14, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_13[k] = ab_x[k] * ip_7[k]
                  + kp_7[k];

        t_14[k] = ab_x[k] * ip_8[k]
                  + kp_8[k];

        t_15[k] = ab_y[k] * ip_7[k]
                  + kp_13[k];

        t_16[k] = ab_y[k] * ip_8[k]
                  + kp_14[k];

        t_17[k] = ab_z[k] * ip_8[k]
                  + kp_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, ip_9, ip_10, ip_11, kp_9, \
                         kp_10, kp_11, kp_19, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_18[k] = ab_x[k] * ip_9[k]
                  + kp_9[k];

        t_19[k] = ab_x[k] * ip_10[k]
                  + kp_10[k];

        t_20[k] = ab_x[k] * ip_11[k]
                  + kp_11[k];

        t_21[k] = ab_y[k] * ip_10[k]
                  + kp_19[k];

        t_22[k] = ab_y[k] * ip_11[k]
                  + kp_20[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, ip_11, ip_12, ip_13, ip_14, \
                         kp_12, kp_13, kp_14, kp_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_23[k] = ab_z[k] * ip_11[k]
                  + kp_23[k];

        t_24[k] = ab_x[k] * ip_12[k]
                  + kp_12[k];

        t_25[k] = ab_x[k] * ip_13[k]
                  + kp_13[k];

        t_26[k] = ab_x[k] * ip_14[k]
                  + kp_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, ip_13, ip_14, ip_15, kp_15, \
                         kp_22, kp_23, kp_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_27[k] = ab_y[k] * ip_13[k]
                  + kp_22[k];

        t_28[k] = ab_y[k] * ip_14[k]
                  + kp_23[k];

        t_29[k] = ab_z[k] * ip_14[k]
                  + kp_26[k];

        t_30[k] = ab_x[k] * ip_15[k]
                  + kp_15[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, ip_16, ip_17, kp_16, \
                         kp_17, kp_25, kp_26, kp_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_31[k] = ab_x[k] * ip_16[k]
                  + kp_16[k];

        t_32[k] = ab_x[k] * ip_17[k]
                  + kp_17[k];

        t_33[k] = ab_y[k] * ip_16[k]
                  + kp_25[k];

        t_34[k] = ab_y[k] * ip_17[k]
                  + kp_26[k];

        t_35[k] = ab_z[k] * ip_17[k]
                  + kp_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, ip_18, ip_19, ip_20, kp_18, \
                         kp_19, kp_20, kp_31, kp_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_36[k] = ab_x[k] * ip_18[k]
                  + kp_18[k];

        t_37[k] = ab_x[k] * ip_19[k]
                  + kp_19[k];

        t_38[k] = ab_x[k] * ip_20[k]
                  + kp_20[k];

        t_39[k] = ab_y[k] * ip_19[k]
                  + kp_31[k];

        t_40[k] = ab_y[k] * ip_20[k]
                  + kp_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_z, ip_20, ip_21, ip_22, ip_23, \
                         kp_21, kp_22, kp_23, kp_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_41[k] = ab_z[k] * ip_20[k]
                  + kp_35[k];

        t_42[k] = ab_x[k] * ip_21[k]
                  + kp_21[k];

        t_43[k] = ab_x[k] * ip_22[k]
                  + kp_22[k];

        t_44[k] = ab_x[k] * ip_23[k]
                  + kp_23[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, ip_22, ip_23, ip_24, kp_24, \
                         kp_34, kp_35, kp_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = ab_y[k] * ip_22[k]
                  + kp_34[k];

        t_46[k] = ab_y[k] * ip_23[k]
                  + kp_35[k];

        t_47[k] = ab_z[k] * ip_23[k]
                  + kp_38[k];

        t_48[k] = ab_x[k] * ip_24[k]
                  + kp_24[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, ip_25, ip_26, kp_25, \
                         kp_26, kp_37, kp_38, kp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_49[k] = ab_x[k] * ip_25[k]
                  + kp_25[k];

        t_50[k] = ab_x[k] * ip_26[k]
                  + kp_26[k];

        t_51[k] = ab_y[k] * ip_25[k]
                  + kp_37[k];

        t_52[k] = ab_y[k] * ip_26[k]
                  + kp_38[k];

        t_53[k] = ab_z[k] * ip_26[k]
                  + kp_41[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, ip_27, ip_28, ip_29, kp_27, \
                         kp_28, kp_29, kp_40, kp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_54[k] = ab_x[k] * ip_27[k]
                  + kp_27[k];

        t_55[k] = ab_x[k] * ip_28[k]
                  + kp_28[k];

        t_56[k] = ab_x[k] * ip_29[k]
                  + kp_29[k];

        t_57[k] = ab_y[k] * ip_28[k]
                  + kp_40[k];

        t_58[k] = ab_y[k] * ip_29[k]
                  + kp_41[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, ab_x, ab_z, ip_29, ip_30, ip_31, ip_32, \
                         kp_30, kp_31, kp_32, kp_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_59[k] = ab_z[k] * ip_29[k]
                  + kp_44[k];

        t_60[k] = ab_x[k] * ip_30[k]
                  + kp_30[k];

        t_61[k] = ab_x[k] * ip_31[k]
                  + kp_31[k];

        t_62[k] = ab_x[k] * ip_32[k]
                  + kp_32[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, ab_x, ab_y, ab_z, ip_31, ip_32, ip_33, kp_33, \
                         kp_46, kp_47, kp_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_63[k] = ab_y[k] * ip_31[k]
                  + kp_46[k];

        t_64[k] = ab_y[k] * ip_32[k]
                  + kp_47[k];

        t_65[k] = ab_z[k] * ip_32[k]
                  + kp_50[k];

        t_66[k] = ab_x[k] * ip_33[k]
                  + kp_33[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, ab_x, ab_y, ab_z, ip_34, ip_35, kp_34, \
                         kp_35, kp_49, kp_50, kp_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_67[k] = ab_x[k] * ip_34[k]
                  + kp_34[k];

        t_68[k] = ab_x[k] * ip_35[k]
                  + kp_35[k];

        t_69[k] = ab_y[k] * ip_34[k]
                  + kp_49[k];

        t_70[k] = ab_y[k] * ip_35[k]
                  + kp_50[k];

        t_71[k] = ab_z[k] * ip_35[k]
                  + kp_53[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, ip_36, ip_37, ip_38, kp_36, \
                         kp_37, kp_38, kp_52, kp_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_72[k] = ab_x[k] * ip_36[k]
                  + kp_36[k];

        t_73[k] = ab_x[k] * ip_37[k]
                  + kp_37[k];

        t_74[k] = ab_x[k] * ip_38[k]
                  + kp_38[k];

        t_75[k] = ab_y[k] * ip_37[k]
                  + kp_52[k];

        t_76[k] = ab_y[k] * ip_38[k]
                  + kp_53[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, ab_x, ab_z, ip_38, ip_39, ip_40, ip_41, \
                         kp_39, kp_40, kp_41, kp_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_77[k] = ab_z[k] * ip_38[k]
                  + kp_56[k];

        t_78[k] = ab_x[k] * ip_39[k]
                  + kp_39[k];

        t_79[k] = ab_x[k] * ip_40[k]
                  + kp_40[k];

        t_80[k] = ab_x[k] * ip_41[k]
                  + kp_41[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, ab_x, ab_y, ab_z, ip_40, ip_41, ip_42, kp_42, \
                         kp_55, kp_56, kp_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_81[k] = ab_y[k] * ip_40[k]
                  + kp_55[k];

        t_82[k] = ab_y[k] * ip_41[k]
                  + kp_56[k];

        t_83[k] = ab_z[k] * ip_41[k]
                  + kp_59[k];

        t_84[k] = ab_x[k] * ip_42[k]
                  + kp_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, ip_43, ip_44, kp_43, \
                         kp_44, kp_58, kp_59, kp_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = ab_x[k] * ip_43[k]
                  + kp_43[k];

        t_86[k] = ab_x[k] * ip_44[k]
                  + kp_44[k];

        t_87[k] = ab_y[k] * ip_43[k]
                  + kp_58[k];

        t_88[k] = ab_y[k] * ip_44[k]
                  + kp_59[k];

        t_89[k] = ab_z[k] * ip_44[k]
                  + kp_62[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, ab_y, ip_45, ip_46, ip_47, kp_45, \
                         kp_46, kp_47, kp_64, kp_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = ab_x[k] * ip_45[k]
                  + kp_45[k];

        t_91[k] = ab_x[k] * ip_46[k]
                  + kp_46[k];

        t_92[k] = ab_x[k] * ip_47[k]
                  + kp_47[k];

        t_93[k] = ab_y[k] * ip_46[k]
                  + kp_64[k];

        t_94[k] = ab_y[k] * ip_47[k]
                  + kp_65[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, ab_x, ab_z, ip_47, ip_48, ip_49, ip_50, \
                         kp_48, kp_49, kp_50, kp_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = ab_z[k] * ip_47[k]
                  + kp_68[k];

        t_96[k] = ab_x[k] * ip_48[k]
                  + kp_48[k];

        t_97[k] = ab_x[k] * ip_49[k]
                  + kp_49[k];

        t_98[k] = ab_x[k] * ip_50[k]
                  + kp_50[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, ab_x, ab_y, ab_z, ip_49, ip_50, ip_51, \
                         kp_51, kp_67, kp_68, kp_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_99[k] = ab_y[k] * ip_49[k]
                  + kp_67[k];

        t_100[k] = ab_y[k] * ip_50[k]
                   + kp_68[k];

        t_101[k] = ab_z[k] * ip_50[k]
                   + kp_71[k];

        t_102[k] = ab_x[k] * ip_51[k]
                   + kp_51[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, ab_x, ab_y, ab_z, ip_52, ip_53, \
                         kp_52, kp_53, kp_70, kp_71, kp_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_103[k] = ab_x[k] * ip_52[k]
                   + kp_52[k];

        t_104[k] = ab_x[k] * ip_53[k]
                   + kp_53[k];

        t_105[k] = ab_y[k] * ip_52[k]
                   + kp_70[k];

        t_106[k] = ab_y[k] * ip_53[k]
                   + kp_71[k];

        t_107[k] = ab_z[k] * ip_53[k]
                   + kp_74[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, ab_x, ab_y, ip_54, ip_55, ip_56, \
                         kp_54, kp_55, kp_56, kp_73, kp_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_108[k] = ab_x[k] * ip_54[k]
                   + kp_54[k];

        t_109[k] = ab_x[k] * ip_55[k]
                   + kp_55[k];

        t_110[k] = ab_x[k] * ip_56[k]
                   + kp_56[k];

        t_111[k] = ab_y[k] * ip_55[k]
                   + kp_73[k];

        t_112[k] = ab_y[k] * ip_56[k]
                   + kp_74[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, ab_x, ab_z, ip_56, ip_57, ip_58, ip_59, \
                         kp_57, kp_58, kp_59, kp_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_113[k] = ab_z[k] * ip_56[k]
                   + kp_77[k];

        t_114[k] = ab_x[k] * ip_57[k]
                   + kp_57[k];

        t_115[k] = ab_x[k] * ip_58[k]
                   + kp_58[k];

        t_116[k] = ab_x[k] * ip_59[k]
                   + kp_59[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, ab_x, ab_y, ab_z, ip_58, ip_59, ip_60, \
                         kp_60, kp_76, kp_77, kp_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_117[k] = ab_y[k] * ip_58[k]
                   + kp_76[k];

        t_118[k] = ab_y[k] * ip_59[k]
                   + kp_77[k];

        t_119[k] = ab_z[k] * ip_59[k]
                   + kp_80[k];

        t_120[k] = ab_x[k] * ip_60[k]
                   + kp_60[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, ab_x, ab_y, ab_z, ip_61, ip_62, \
                         kp_61, kp_62, kp_79, kp_80, kp_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_121[k] = ab_x[k] * ip_61[k]
                   + kp_61[k];

        t_122[k] = ab_x[k] * ip_62[k]
                   + kp_62[k];

        t_123[k] = ab_y[k] * ip_61[k]
                   + kp_79[k];

        t_124[k] = ab_y[k] * ip_62[k]
                   + kp_80[k];

        t_125[k] = ab_z[k] * ip_62[k]
                   + kp_83[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, ab_x, ab_y, ip_63, ip_64, ip_65, \
                         kp_63, kp_64, kp_65, kp_85, kp_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_126[k] = ab_x[k] * ip_63[k]
                   + kp_63[k];

        t_127[k] = ab_x[k] * ip_64[k]
                   + kp_64[k];

        t_128[k] = ab_x[k] * ip_65[k]
                   + kp_65[k];

        t_129[k] = ab_y[k] * ip_64[k]
                   + kp_85[k];

        t_130[k] = ab_y[k] * ip_65[k]
                   + kp_86[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, ab_x, ab_z, ip_65, ip_66, ip_67, ip_68, \
                         kp_66, kp_67, kp_68, kp_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_131[k] = ab_z[k] * ip_65[k]
                   + kp_89[k];

        t_132[k] = ab_x[k] * ip_66[k]
                   + kp_66[k];

        t_133[k] = ab_x[k] * ip_67[k]
                   + kp_67[k];

        t_134[k] = ab_x[k] * ip_68[k]
                   + kp_68[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, ab_x, ab_y, ab_z, ip_67, ip_68, ip_69, \
                         kp_69, kp_88, kp_89, kp_92 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = ab_y[k] * ip_67[k]
                   + kp_88[k];

        t_136[k] = ab_y[k] * ip_68[k]
                   + kp_89[k];

        t_137[k] = ab_z[k] * ip_68[k]
                   + kp_92[k];

        t_138[k] = ab_x[k] * ip_69[k]
                   + kp_69[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_x, ab_y, ab_z, ip_70, ip_71, \
                         kp_70, kp_71, kp_91, kp_92, kp_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_139[k] = ab_x[k] * ip_70[k]
                   + kp_70[k];

        t_140[k] = ab_x[k] * ip_71[k]
                   + kp_71[k];

        t_141[k] = ab_y[k] * ip_70[k]
                   + kp_91[k];

        t_142[k] = ab_y[k] * ip_71[k]
                   + kp_92[k];

        t_143[k] = ab_z[k] * ip_71[k]
                   + kp_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, ab_x, ab_y, ip_72, ip_73, ip_74, \
                         kp_72, kp_73, kp_74, kp_94, kp_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_144[k] = ab_x[k] * ip_72[k]
                   + kp_72[k];

        t_145[k] = ab_x[k] * ip_73[k]
                   + kp_73[k];

        t_146[k] = ab_x[k] * ip_74[k]
                   + kp_74[k];

        t_147[k] = ab_y[k] * ip_73[k]
                   + kp_94[k];

        t_148[k] = ab_y[k] * ip_74[k]
                   + kp_95[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, ab_x, ab_z, ip_74, ip_75, ip_76, ip_77, \
                         kp_75, kp_76, kp_77, kp_98 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_149[k] = ab_z[k] * ip_74[k]
                   + kp_98[k];

        t_150[k] = ab_x[k] * ip_75[k]
                   + kp_75[k];

        t_151[k] = ab_x[k] * ip_76[k]
                   + kp_76[k];

        t_152[k] = ab_x[k] * ip_77[k]
                   + kp_77[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, ab_x, ab_y, ab_z, ip_76, ip_77, ip_78, \
                         kp_78, kp_97, kp_98, kp_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_153[k] = ab_y[k] * ip_76[k]
                   + kp_97[k];

        t_154[k] = ab_y[k] * ip_77[k]
                   + kp_98[k];

        t_155[k] = ab_z[k] * ip_77[k]
                   + kp_101[k];

        t_156[k] = ab_x[k] * ip_78[k]
                   + kp_78[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, ab_x, ab_y, ab_z, ip_79, ip_80, \
                         kp_79, kp_80, kp_100, kp_101, kp_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_157[k] = ab_x[k] * ip_79[k]
                   + kp_79[k];

        t_158[k] = ab_x[k] * ip_80[k]
                   + kp_80[k];

        t_159[k] = ab_y[k] * ip_79[k]
                   + kp_100[k];

        t_160[k] = ab_y[k] * ip_80[k]
                   + kp_101[k];

        t_161[k] = ab_z[k] * ip_80[k]
                   + kp_104[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, ab_x, ab_y, ip_81, ip_82, ip_83, \
                         kp_81, kp_82, kp_83, kp_103, kp_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_162[k] = ab_x[k] * ip_81[k]
                   + kp_81[k];

        t_163[k] = ab_x[k] * ip_82[k]
                   + kp_82[k];

        t_164[k] = ab_x[k] * ip_83[k]
                   + kp_83[k];

        t_165[k] = ab_y[k] * ip_82[k]
                   + kp_103[k];

        t_166[k] = ab_y[k] * ip_83[k]
                   + kp_104[k];
    }

#pragma omp simd aligned(t_167, ab_z, ip_83, kp_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_167[k] = ab_z[k] * ip_83[k]
                   + kp_107[k];
    }
}

}  // namespace simdtrf
