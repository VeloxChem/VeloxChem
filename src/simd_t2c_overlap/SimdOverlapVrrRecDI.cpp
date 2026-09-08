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


#include "SimdOverlapVrrRecDI.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_di_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t ph, const size_t pi, const size_t dg,
                          const size_t dh, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);
    const auto *ph_9 = buffer.data(ph + 9);
    const auto *ph_10 = buffer.data(ph + 10);
    const auto *ph_12 = buffer.data(ph + 12);
    const auto *ph_13 = buffer.data(ph + 13);
    const auto *ph_14 = buffer.data(ph + 14);
    const auto *ph_16 = buffer.data(ph + 16);
    const auto *ph_18 = buffer.data(ph + 18);
    const auto *ph_19 = buffer.data(ph + 19);
    const auto *ph_20 = buffer.data(ph + 20);
    const auto *ph_21 = buffer.data(ph + 21);
    const auto *ph_24 = buffer.data(ph + 24);
    const auto *ph_27 = buffer.data(ph + 27);
    const auto *ph_30 = buffer.data(ph + 30);
    const auto *ph_31 = buffer.data(ph + 31);
    const auto *ph_32 = buffer.data(ph + 32);
    const auto *ph_33 = buffer.data(ph + 33);
    const auto *ph_34 = buffer.data(ph + 34);
    const auto *ph_35 = buffer.data(ph + 35);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_16 = buffer.data(pi + 16);
    const auto *pi_17 = buffer.data(pi + 17);
    const auto *pi_18 = buffer.data(pi + 18);
    const auto *pi_19 = buffer.data(pi + 19);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_21 = buffer.data(pi + 21);
    const auto *pi_23 = buffer.data(pi + 23);
    const auto *pi_26 = buffer.data(pi + 26);
    const auto *pi_30 = buffer.data(pi + 30);
    const auto *pi_31 = buffer.data(pi + 31);
    const auto *pi_32 = buffer.data(pi + 32);
    const auto *pi_33 = buffer.data(pi + 33);
    const auto *pi_34 = buffer.data(pi + 34);
    const auto *pi_35 = buffer.data(pi + 35);
    const auto *pi_36 = buffer.data(pi + 36);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
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
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, ph_0, dg_0, dh_0, \
                         dh_1, dh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ph_0[k]
                 + f_1 * dg_0[k]
                 + pb_x[k] * dh_0[k];

        t_1[k] = pb_y[k] * dh_0[k];

        t_2[k] = pb_z[k] * dh_0[k];

        t_3[k] = f_2 * dg_0[k]
                 + pb_y[k] * dh_1[k];

        t_4[k] = pb_y[k] * dh_2[k];

        t_5[k] = f_2 * dg_0[k]
                 + pb_z[k] * dh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_y, pb_z, dg_1, dg_2, dg_3, dh_3, \
                         dh_4, dh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * dg_1[k]
                 + pb_y[k] * dh_3[k];

        t_7[k] = pb_z[k] * dh_3[k];

        t_8[k] = pb_y[k] * dh_4[k];

        t_9[k] = f_0 * dg_2[k]
                 + pb_z[k] * dh_4[k];

        t_10[k] = f_3 * dg_3[k]
                  + pb_y[k] * dh_5[k];

        t_11[k] = pb_z[k] * dh_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, pb_z, ph_7, dg_4, dh_6, \
                         dh_7, dh_8, dh_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * dg_4[k]
                  + pb_y[k] * dh_6[k];

        t_13[k] = pb_y[k] * dh_7[k];

        t_14[k] = f_3 * dg_4[k]
                  + pb_z[k] * dh_7[k];

        t_15[k] = f_0 * ph_7[k]
                  + pb_x[k] * dh_10[k];

        t_16[k] = pb_z[k] * dh_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pb_x, pb_y, ph_8, ph_9, ph_10, dg_5, \
                         dh_9, dh_10, dh_11, dh_12, dh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * ph_8[k]
                  + pb_x[k] * dh_11[k];

        t_18[k] = f_0 * ph_9[k]
                  + pb_x[k] * dh_12[k];

        t_19[k] = pb_y[k] * dh_9[k];

        t_20[k] = f_0 * ph_10[k]
                  + pb_x[k] * dh_14[k];

        t_21[k] = f_1 * dg_5[k]
                  + pb_y[k] * dh_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pb_y, pb_z, dg_6, dg_7, dg_8, \
                         dh_10, dh_11, dh_12, dh_13, dh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_z[k] * dh_10[k];

        t_23[k] = f_3 * dg_6[k]
                  + pb_y[k] * dh_11[k];

        t_24[k] = f_0 * dg_7[k]
                  + pb_y[k] * dh_12[k];

        t_25[k] = f_2 * dg_8[k]
                  + pb_y[k] * dh_13[k];

        t_26[k] = pb_y[k] * dh_14[k];

        t_27[k] = f_1 * dg_8[k]
                  + pb_z[k] * dh_14[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_x, pa_y, pb_y, pb_z, ph_0, ph_12, \
                         pi_0, pi_10, dh_15, dh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_y[k] * pi_0[k];

        t_29[k] = f_2 * ph_0[k]
                  + pb_y[k] * dh_15[k];

        t_30[k] = pb_z[k] * dh_15[k];

        t_31[k] = f_4 * ph_12[k]
                  + pa_x[k] * pi_10[k];

        t_32[k] = pb_z[k] * dh_16[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_x, pa_y, pb_y, pb_z, ph_4, ph_13, \
                         pi_2, pi_4, pi_11, dh_17, dh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_y[k] * pi_2[k];

        t_34[k] = f_3 * ph_13[k]
                  + pa_x[k] * pi_11[k];

        t_35[k] = pb_z[k] * dh_17[k];

        t_36[k] = f_2 * ph_4[k]
                  + pb_y[k] * dh_18[k];

        t_37[k] = pa_y[k] * pi_4[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pb_y, pb_z, ph_6, ph_14, pi_12, dg_10, \
                         dh_19, dh_20, dh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * ph_14[k]
                  + pa_x[k] * pi_12[k];

        t_39[k] = pb_z[k] * dh_19[k];

        t_40[k] = f_2 * dg_10[k]
                  + pb_z[k] * dh_20[k];

        t_41[k] = f_2 * ph_6[k]
                  + pb_y[k] * dh_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_y, pb_x, pb_z, ph_16, ph_18, ph_19, \
                         pi_6, dh_22, dh_23, dh_24, dh_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * pi_6[k];

        t_43[k] = f_2 * ph_16[k]
                  + pb_x[k] * dh_23[k];

        t_44[k] = pb_z[k] * dh_22[k];

        t_45[k] = f_2 * ph_18[k]
                  + pb_x[k] * dh_24[k];

        t_46[k] = f_2 * ph_19[k]
                  + pb_x[k] * dh_25[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pa_x, pa_y, pb_x, pb_z, ph_20, pi_8, \
                         pi_14, pi_15, dh_23, dh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_2 * ph_20[k]
                  + pb_x[k] * dh_26[k];

        t_48[k] = pa_y[k] * pi_8[k];

        t_49[k] = pa_x[k] * pi_14[k];

        t_50[k] = pb_z[k] * dh_23[k];

        t_51[k] = pa_x[k] * pi_15[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, pa_x, pa_z, pb_y, pi_0, pi_16, \
                         pi_17, pi_18, pi_19, dh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pa_x[k] * pi_16[k];

        t_53[k] = pa_x[k] * pi_17[k];

        t_54[k] = pa_x[k] * pi_18[k];

        t_55[k] = pa_x[k] * pi_19[k];

        t_56[k] = pa_z[k] * pi_0[k];

        t_57[k] = pb_y[k] * dh_27[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_x, pa_z, pb_y, pb_z, ph_0, ph_24, \
                         pi_1, pi_3, pi_23, dh_27, dh_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_2 * ph_0[k]
                  + pb_z[k] * dh_27[k];

        t_59[k] = pa_z[k] * pi_1[k];

        t_60[k] = pb_y[k] * dh_28[k];

        t_61[k] = f_4 * ph_24[k]
                  + pa_x[k] * pi_23[k];

        t_62[k] = pa_z[k] * pi_3[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_x, pa_z, pb_y, ph_27, pi_5, pi_26, \
                         dg_13, dg_14, dh_29, dh_30, dh_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_2 * dg_13[k]
                  + pb_y[k] * dh_29[k];

        t_64[k] = pb_y[k] * dh_30[k];

        t_65[k] = f_3 * ph_27[k]
                  + pa_x[k] * pi_26[k];

        t_66[k] = pa_z[k] * pi_5[k];

        t_67[k] = f_0 * dg_14[k]
                  + pb_y[k] * dh_31[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_x, pa_z, pb_y, ph_30, pi_7, pi_30, dg_15, \
                         dh_32, dh_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * dg_15[k]
                  + pb_y[k] * dh_32[k];

        t_69[k] = pb_y[k] * dh_33[k];

        t_70[k] = f_0 * ph_30[k]
                  + pa_x[k] * pi_30[k];

        t_71[k] = pa_z[k] * pi_7[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pb_x, pb_y, ph_31, ph_32, ph_33, ph_35, \
                         dh_34, dh_35, dh_36, dh_37, dh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_2 * ph_31[k]
                  + pb_x[k] * dh_35[k];

        t_73[k] = f_2 * ph_32[k]
                  + pb_x[k] * dh_36[k];

        t_74[k] = f_2 * ph_33[k]
                  + pb_x[k] * dh_37[k];

        t_75[k] = pb_y[k] * dh_34[k];

        t_76[k] = f_2 * ph_35[k]
                  + pb_x[k] * dh_38[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, t_83, pa_x, pb_y, pi_31, pi_32, \
                         pi_33, pi_34, pi_35, pi_36, dh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_x[k] * pi_31[k];

        t_78[k] = pa_x[k] * pi_32[k];

        t_79[k] = pa_x[k] * pi_33[k];

        t_80[k] = pa_x[k] * pi_34[k];

        t_81[k] = pa_x[k] * pi_35[k];

        t_82[k] = pb_y[k] * dh_38[k];

        t_83[k] = pa_x[k] * pi_36[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, pb_x, pb_z, dg_17, dg_18, dg_19, \
                         dg_20, dh_39, dh_40, dh_41, dh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * dg_17[k]
                  + pb_x[k] * dh_39[k];

        t_85[k] = f_4 * dg_18[k]
                  + pb_x[k] * dh_40[k];

        t_86[k] = pb_z[k] * dh_39[k];

        t_87[k] = f_3 * dg_19[k]
                  + pb_x[k] * dh_41[k];

        t_88[k] = pb_z[k] * dh_40[k];

        t_89[k] = f_3 * dg_20[k]
                  + pb_x[k] * dh_42[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pb_x, pb_z, dg_21, dg_22, dg_23, dg_24, \
                         dh_41, dh_43, dh_44, dh_45, dh_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * dg_21[k]
                  + pb_x[k] * dh_43[k];

        t_91[k] = pb_z[k] * dh_41[k];

        t_92[k] = f_0 * dg_22[k]
                  + pb_x[k] * dh_44[k];

        t_93[k] = f_0 * dg_23[k]
                  + pb_x[k] * dh_45[k];

        t_94[k] = f_2 * dg_24[k]
                  + pb_x[k] * dh_46[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pb_x, pb_z, dg_26, dg_27, dg_28, dh_43, \
                         dh_47, dh_48, dh_49, dh_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pb_z[k] * dh_43[k];

        t_96[k] = f_2 * dg_26[k]
                  + pb_x[k] * dh_47[k];

        t_97[k] = f_2 * dg_27[k]
                  + pb_x[k] * dh_48[k];

        t_98[k] = f_2 * dg_28[k]
                  + pb_x[k] * dh_49[k];

        t_99[k] = pb_x[k] * dh_50[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, t_105, pb_x, pb_y, ph_16, dg_24, \
                         dh_50, dh_51, dh_52, dh_53, dh_54, dh_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pb_x[k] * dh_51[k];

        t_101[k] = pb_x[k] * dh_52[k];

        t_102[k] = pb_x[k] * dh_53[k];

        t_103[k] = pb_x[k] * dh_54[k];

        t_104[k] = pb_x[k] * dh_55[k];

        t_105[k] = f_0 * ph_16[k]
                   + f_1 * dg_24[k]
                   + pb_y[k] * dh_50[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, pb_y, pb_z, ph_21, dg_24, dg_25, \
                         dg_26, dh_50, dh_51, dh_52, dh_53, dh_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = pb_z[k] * dh_50[k];

        t_107[k] = f_2 * dg_24[k]
                   + pb_z[k] * dh_51[k];

        t_108[k] = f_0 * dg_25[k]
                   + pb_z[k] * dh_52[k];

        t_109[k] = f_3 * dg_26[k]
                   + pb_z[k] * dh_53[k];

        t_110[k] = f_0 * ph_21[k]
                   + pb_y[k] * dh_55[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, pa_y, pa_z, pb_z, pi_9, pi_10, \
                         pi_20, pi_21, dg_28, dh_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_1 * dg_28[k]
                   + pb_z[k] * dh_55[k];

        t_112[k] = pa_y[k] * pi_20[k];

        t_113[k] = pa_z[k] * pi_9[k];

        t_114[k] = pa_y[k] * pi_21[k];

        t_115[k] = pa_z[k] * pi_10[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, pa_y, pa_z, pb_x, pi_11, pi_23, \
                         dg_29, dg_30, dg_31, dh_56, dh_57, dh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_3 * dg_29[k]
                   + pb_x[k] * dh_56[k];

        t_117[k] = pa_y[k] * pi_23[k];

        t_118[k] = pa_z[k] * pi_11[k];

        t_119[k] = f_0 * dg_30[k]
                   + pb_x[k] * dh_57[k];

        t_120[k] = f_0 * dg_31[k]
                   + pb_x[k] * dh_58[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pa_y, pa_z, pb_x, pi_12, pi_26, \
                         dg_33, dg_34, dg_35, dh_59, dh_60, dh_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * pi_26[k];

        t_122[k] = pa_z[k] * pi_12[k];

        t_123[k] = f_2 * dg_33[k]
                   + pb_x[k] * dh_59[k];

        t_124[k] = f_2 * dg_34[k]
                   + pb_x[k] * dh_60[k];

        t_125[k] = f_2 * dg_35[k]
                   + pb_x[k] * dh_61[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, t_131, t_132, pa_y, pb_x, pi_30, \
                         dh_62, dh_63, dh_64, dh_65, dh_66, dh_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_y[k] * pi_30[k];

        t_127[k] = pb_x[k] * dh_62[k];

        t_128[k] = pb_x[k] * dh_63[k];

        t_129[k] = pb_x[k] * dh_64[k];

        t_130[k] = pb_x[k] * dh_65[k];

        t_131[k] = pb_x[k] * dh_66[k];

        t_132[k] = pb_x[k] * dh_67[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_y, pa_z, pb_z, ph_16, ph_32, ph_33, \
                         pi_14, pi_33, pi_34, dh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * pi_14[k];

        t_134[k] = f_2 * ph_16[k]
                   + pb_z[k] * dh_62[k];

        t_135[k] = f_4 * ph_32[k]
                   + pa_y[k] * pi_33[k];

        t_136[k] = f_3 * ph_33[k]
                   + pa_y[k] * pi_34[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, pa_y, pb_x, pb_y, ph_34, ph_35, \
                         pi_35, pi_36, dg_37, dh_67, dh_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_0 * ph_34[k]
                   + pa_y[k] * pi_35[k];

        t_138[k] = f_2 * ph_35[k]
                   + pb_y[k] * dh_67[k];

        t_139[k] = pa_y[k] * pi_36[k];

        t_140[k] = f_1 * dg_37[k]
                   + pb_x[k] * dh_68[k];

        t_141[k] = pb_y[k] * dh_68[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, pb_x, pb_y, dg_38, dg_39, dg_40, \
                         dg_41, dh_69, dh_70, dh_71, dh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_4 * dg_38[k]
                   + pb_x[k] * dh_69[k];

        t_143[k] = f_3 * dg_39[k]
                   + pb_x[k] * dh_70[k];

        t_144[k] = pb_y[k] * dh_69[k];

        t_145[k] = f_3 * dg_40[k]
                   + pb_x[k] * dh_71[k];

        t_146[k] = f_0 * dg_41[k]
                   + pb_x[k] * dh_72[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, pb_x, pb_y, dg_42, dg_43, dg_44, \
                         dg_45, dh_71, dh_73, dh_74, dh_75, dh_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_0 * dg_42[k]
                   + pb_x[k] * dh_73[k];

        t_148[k] = pb_y[k] * dh_71[k];

        t_149[k] = f_0 * dg_43[k]
                   + pb_x[k] * dh_74[k];

        t_150[k] = f_2 * dg_44[k]
                   + pb_x[k] * dh_75[k];

        t_151[k] = f_2 * dg_45[k]
                   + pb_x[k] * dh_76[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, t_157, pb_x, pb_y, dg_46, dg_48, \
                         dh_74, dh_77, dh_78, dh_79, dh_80, dh_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_2 * dg_46[k]
                   + pb_x[k] * dh_77[k];

        t_153[k] = pb_y[k] * dh_74[k];

        t_154[k] = f_2 * dg_48[k]
                   + pb_x[k] * dh_78[k];

        t_155[k] = pb_x[k] * dh_79[k];

        t_156[k] = pb_x[k] * dh_80[k];

        t_157[k] = pb_x[k] * dh_81[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pb_x, pb_y, dg_44, dg_45, dh_79, \
                         dh_80, dh_82, dh_83, dh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = pb_x[k] * dh_82[k];

        t_159[k] = pb_x[k] * dh_83[k];

        t_160[k] = pb_x[k] * dh_84[k];

        t_161[k] = f_1 * dg_44[k]
                   + pb_y[k] * dh_79[k];

        t_162[k] = f_4 * dg_45[k]
                   + pb_y[k] * dh_80[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pb_y, pb_z, ph_35, dg_46, dg_47, \
                         dg_48, dh_81, dh_82, dh_83, dh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_3 * dg_46[k]
                   + pb_y[k] * dh_81[k];

        t_164[k] = f_0 * dg_47[k]
                   + pb_y[k] * dh_82[k];

        t_165[k] = f_2 * dg_48[k]
                   + pb_y[k] * dh_83[k];

        t_166[k] = pb_y[k] * dh_84[k];

        t_167[k] = f_0 * ph_35[k]
                   + f_1 * dg_48[k]
                   + pb_z[k] * dh_84[k];
    }
}

auto
compute_prim_di_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t ph, const size_t pi, const size_t dg,
                          const size_t dh, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);
    const auto *ph_12 = buffer.data(ph + 12);
    const auto *ph_13 = buffer.data(ph + 13);
    const auto *ph_15 = buffer.data(ph + 15);
    const auto *ph_17 = buffer.data(ph + 17);
    const auto *ph_19 = buffer.data(ph + 19);
    const auto *ph_20 = buffer.data(ph + 20);
    const auto *ph_21 = buffer.data(ph + 21);
    const auto *ph_22 = buffer.data(ph + 22);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_8 = buffer.data(dg + 8);
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
    const auto *dg_26 = buffer.data(dg + 26);
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

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
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
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_55 = buffer.data(dh + 55);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_58 = buffer.data(dh + 58);
    const auto *dh_59 = buffer.data(dh + 59);
    const auto *dh_61 = buffer.data(dh + 61);
    const auto *dh_62 = buffer.data(dh + 62);
    const auto *dh_63 = buffer.data(dh + 63);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ph_0, dg_0, dg_1, dh_0, \
                         dh_1, dh_2, dh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ph_0[k]
                 + f_1 * dg_0[k]
                 + pb_x[k] * dh_0[k];

        t_1[k] = f_2 * dg_0[k]
                 + pb_y[k] * dh_1[k];

        t_2[k] = f_2 * dg_0[k]
                 + pb_z[k] * dh_2[k];

        t_3[k] = f_0 * dg_1[k]
                 + pb_y[k] * dh_3[k];

        t_4[k] = pb_z[k] * dh_3[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, dg_2, dg_3, dg_4, dh_4, dh_5, \
                         dh_7, dh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_2[k]
                 + pb_z[k] * dh_4[k];

        t_6[k] = f_3 * dg_3[k]
                 + pb_y[k] * dh_5[k];

        t_7[k] = pb_z[k] * dh_5[k];

        t_8[k] = f_2 * dg_4[k]
                 + pb_y[k] * dh_7[k];

        t_9[k] = f_3 * dg_4[k]
                 + pb_z[k] * dh_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_x, pb_y, pb_z, ph_0, ph_3, ph_4, \
                         dg_5, dg_8, dh_9, dh_12, dh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * ph_3[k]
                  + pb_x[k] * dh_9[k];

        t_11[k] = f_0 * ph_4[k]
                  + pb_x[k] * dh_12[k];

        t_12[k] = f_1 * dg_5[k]
                  + pb_y[k] * dh_9[k];

        t_13[k] = f_1 * dg_8[k]
                  + pb_z[k] * dh_12[k];

        t_14[k] = f_2 * ph_0[k]
                  + pb_y[k] * dh_13[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_x, pb_x, ph_5, ph_6, ph_7, ph_8, \
                         pi_1, pi_2, pi_3, pi_4, dh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * ph_5[k]
                  + pa_x[k] * pi_1[k];

        t_16[k] = f_3 * ph_6[k]
                  + pa_x[k] * pi_2[k];

        t_17[k] = f_0 * ph_7[k]
                  + pa_x[k] * pi_3[k];

        t_18[k] = f_2 * ph_8[k]
                  + pb_x[k] * dh_18[k];

        t_19[k] = pa_x[k] * pi_4[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pa_z, pb_z, ph_0, ph_13, ph_15, pi_0, \
                         pi_5, pi_6, dh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * pi_0[k];

        t_21[k] = f_2 * ph_0[k]
                  + pb_z[k] * dh_22[k];

        t_22[k] = f_4 * ph_13[k]
                  + pa_x[k] * pi_5[k];

        t_23[k] = f_3 * ph_15[k]
                  + pa_x[k] * pi_6[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_x, pb_x, ph_17, ph_22, pi_8, pi_12, \
                         dg_13, dg_14, dh_29, dh_30, dh_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * ph_17[k]
                  + pa_x[k] * pi_8[k];

        t_25[k] = f_2 * ph_22[k]
                  + pb_x[k] * dh_29[k];

        t_26[k] = pa_x[k] * pi_12[k];

        t_27[k] = f_1 * dg_13[k]
                  + pb_x[k] * dh_30[k];

        t_28[k] = f_4 * dg_14[k]
                  + pb_x[k] * dh_31[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pb_x, pb_z, dg_15, dg_16, dg_17, dh_31, \
                         dh_32, dh_34, dh_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_3 * dg_15[k]
                  + pb_x[k] * dh_32[k];

        t_30[k] = pb_z[k] * dh_31[k];

        t_31[k] = f_3 * dg_16[k]
                  + pb_x[k] * dh_34[k];

        t_32[k] = f_0 * dg_17[k]
                  + pb_x[k] * dh_35[k];

        t_33[k] = pb_z[k] * dh_32[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pb_x, pb_z, dg_18, dg_19, dg_20, dg_22, \
                         dh_35, dh_37, dh_38, dh_39, dh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_0 * dg_18[k]
                  + pb_x[k] * dh_37[k];

        t_35[k] = f_0 * dg_19[k]
                  + pb_x[k] * dh_38[k];

        t_36[k] = f_2 * dg_20[k]
                  + pb_x[k] * dh_39[k];

        t_37[k] = pb_z[k] * dh_35[k];

        t_38[k] = f_2 * dg_22[k]
                  + pb_x[k] * dh_40[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pb_x, pb_y, pb_z, ph_8, dg_20, dg_23, dg_24, \
                         dh_41, dh_42, dh_43, dh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_2 * dg_23[k]
                  + pb_x[k] * dh_41[k];

        t_40[k] = f_2 * dg_24[k]
                  + pb_x[k] * dh_42[k];

        t_41[k] = f_0 * ph_8[k]
                  + f_1 * dg_20[k]
                  + pb_y[k] * dh_43[k];

        t_42[k] = f_2 * dg_20[k]
                  + pb_z[k] * dh_44[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pb_y, pb_z, ph_12, dg_21, dg_22, dg_24, \
                         dh_45, dh_46, dh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_0 * dg_21[k]
                  + pb_z[k] * dh_45[k];

        t_44[k] = f_3 * dg_22[k]
                  + pb_z[k] * dh_46[k];

        t_45[k] = f_0 * ph_12[k]
                  + pb_y[k] * dh_48[k];

        t_46[k] = f_1 * dg_24[k]
                  + pb_z[k] * dh_48[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pa_z, pb_x, pb_z, ph_8, ph_19, pi_4, \
                         pi_9, dg_26, dh_49, dh_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_2 * dg_26[k]
                  + pb_x[k] * dh_49[k];

        t_48[k] = pa_z[k] * pi_4[k];

        t_49[k] = f_2 * ph_8[k]
                  + pb_z[k] * dh_50[k];

        t_50[k] = f_4 * ph_19[k]
                  + pa_y[k] * pi_9[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pb_y, ph_20, ph_21, ph_22, pi_10, \
                         pi_11, pi_12, dh_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * ph_20[k]
                  + pa_y[k] * pi_10[k];

        t_52[k] = f_0 * ph_21[k]
                  + pa_y[k] * pi_11[k];

        t_53[k] = f_2 * ph_22[k]
                  + pb_y[k] * dh_55[k];

        t_54[k] = pa_y[k] * pi_12[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, pb_x, pb_y, dg_28, dg_29, dg_30, \
                         dg_31, dh_56, dh_58, dh_59, dh_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * dg_28[k]
                  + pb_x[k] * dh_56[k];

        t_56[k] = pb_y[k] * dh_56[k];

        t_57[k] = f_4 * dg_29[k]
                  + pb_x[k] * dh_58[k];

        t_58[k] = f_3 * dg_30[k]
                  + pb_x[k] * dh_59[k];

        t_59[k] = pb_y[k] * dh_58[k];

        t_60[k] = f_3 * dg_31[k]
                  + pb_x[k] * dh_61[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pb_x, pb_y, dg_32, dg_33, dg_34, dg_35, \
                         dh_61, dh_62, dh_63, dh_65, dh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * dg_32[k]
                  + pb_x[k] * dh_62[k];

        t_62[k] = f_0 * dg_33[k]
                  + pb_x[k] * dh_63[k];

        t_63[k] = pb_y[k] * dh_61[k];

        t_64[k] = f_0 * dg_34[k]
                  + pb_x[k] * dh_65[k];

        t_65[k] = f_2 * dg_35[k]
                  + pb_x[k] * dh_66[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pb_x, pb_y, dg_35, dg_36, dg_37, dg_39, \
                         dh_65, dh_67, dh_68, dh_69, dh_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_2 * dg_36[k]
                  + pb_x[k] * dh_67[k];

        t_67[k] = f_2 * dg_37[k]
                  + pb_x[k] * dh_68[k];

        t_68[k] = pb_y[k] * dh_65[k];

        t_69[k] = f_2 * dg_39[k]
                  + pb_x[k] * dh_69[k];

        t_70[k] = f_1 * dg_35[k]
                  + pb_y[k] * dh_70[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pb_y, dg_36, dg_37, dg_38, dg_39, dh_71, \
                         dh_72, dh_73, dh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_4 * dg_36[k]
                  + pb_y[k] * dh_71[k];

        t_72[k] = f_3 * dg_37[k]
                  + pb_y[k] * dh_72[k];

        t_73[k] = f_0 * dg_38[k]
                  + pb_y[k] * dh_73[k];

        t_74[k] = f_2 * dg_39[k]
                  + pb_y[k] * dh_74[k];
    }

#pragma omp simd aligned(t_75, pb_z, ph_22, dg_39, dh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * ph_22[k]
                  + f_1 * dg_39[k]
                  + pb_z[k] * dh_75[k];
    }
}

auto
compute_prim_di_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t ph, const size_t pi, const size_t dg,
                          const size_t dh, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);
    const auto *ph_9 = buffer.data(ph + 9);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_55 = buffer.data(dh + 55);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, ph_0, dg_0, dg_1, dh_0, dh_1, \
                         dh_2, dh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ph_0[k]
                 + f_1 * dg_0[k]
                 + pb_x[k] * dh_0[k];

        t_1[k] = f_2 * dg_0[k]
                 + pb_y[k] * dh_1[k];

        t_2[k] = f_2 * dg_0[k]
                 + pb_z[k] * dh_2[k];

        t_3[k] = f_0 * dg_1[k]
                 + pb_y[k] * dh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_y, pb_z, pi_0, dg_2, dg_3, dg_4, \
                         dh_4, dh_5, dh_7, dh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_0 * dg_2[k]
                 + pb_z[k] * dh_4[k];

        t_5[k] = f_3 * dg_3[k]
                 + pb_y[k] * dh_5[k];

        t_6[k] = f_2 * dg_4[k]
                 + pb_y[k] * dh_7[k];

        t_7[k] = f_3 * dg_4[k]
                 + pb_z[k] * dh_8[k];

        t_8[k] = pa_x[k] * pi_0[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_x, pb_x, pi_4, dg_10, dg_11, dg_13, \
                         dg_15, dh_20, dh_21, dh_23, dh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_x[k] * pi_4[k];

        t_10[k] = f_1 * dg_10[k]
                  + pb_x[k] * dh_20[k];

        t_11[k] = f_3 * dg_11[k]
                  + pb_x[k] * dh_21[k];

        t_12[k] = f_0 * dg_13[k]
                  + pb_x[k] * dh_23[k];

        t_13[k] = f_2 * dg_15[k]
                  + pb_x[k] * dh_26[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_x, pb_y, pb_z, ph_3, dg_15, dg_16, dg_17, \
                         dh_27, dh_28, dh_29, dh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_2 * dg_17[k]
                  + pb_x[k] * dh_27[k];

        t_15[k] = f_0 * ph_3[k]
                  + f_1 * dg_15[k]
                  + pb_y[k] * dh_28[k];

        t_16[k] = f_2 * dg_15[k]
                  + pb_z[k] * dh_29[k];

        t_17[k] = f_0 * dg_16[k]
                  + pb_z[k] * dh_30[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_y, pb_z, ph_6, ph_7, ph_8, pi_1, \
                         pi_2, pi_3, pi_4, dg_17, dh_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * dg_17[k]
                  + pb_z[k] * dh_31[k];

        t_19[k] = f_4 * ph_6[k]
                  + pa_y[k] * pi_1[k];

        t_20[k] = f_3 * ph_7[k]
                  + pa_y[k] * pi_2[k];

        t_21[k] = f_0 * ph_8[k]
                  + pa_y[k] * pi_3[k];

        t_22[k] = pa_y[k] * pi_4[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pb_x, dg_21, dg_23, dg_24, dg_25, \
                         dg_26, dh_40, dh_42, dh_43, dh_44, dh_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * dg_21[k]
                  + pb_x[k] * dh_40[k];

        t_24[k] = f_3 * dg_23[k]
                  + pb_x[k] * dh_42[k];

        t_25[k] = f_3 * dg_24[k]
                  + pb_x[k] * dh_43[k];

        t_26[k] = f_0 * dg_25[k]
                  + pb_x[k] * dh_44[k];

        t_27[k] = f_0 * dg_26[k]
                  + pb_x[k] * dh_46[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pb_x, pb_y, dg_27, dg_29, dg_31, dh_47, \
                         dh_48, dh_49, dh_50, dh_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_2 * dg_27[k]
                  + pb_x[k] * dh_47[k];

        t_29[k] = f_2 * dg_29[k]
                  + pb_x[k] * dh_48[k];

        t_30[k] = f_2 * dg_31[k]
                  + pb_x[k] * dh_49[k];

        t_31[k] = f_1 * dg_27[k]
                  + pb_y[k] * dh_50[k];

        t_32[k] = f_3 * dg_29[k]
                  + pb_y[k] * dh_52[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_y, pb_z, ph_9, dg_30, dg_31, dh_53, dh_54, \
                         dh_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * dg_30[k]
                  + pb_y[k] * dh_53[k];

        t_34[k] = f_2 * dg_31[k]
                  + pb_y[k] * dh_54[k];

        t_35[k] = f_0 * ph_9[k]
                  + f_1 * dg_31[k]
                  + pb_z[k] * dh_55[k];
    }
}

auto
compute_prim_di_overlap_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t ph, const size_t pi, const size_t dg,
                          const size_t dh, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_z, ph_0, dg_0, dg_2, dg_4, dh_0, dh_1, \
                         dh_3, dh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ph_0[k]
                 + f_1 * dg_0[k]
                 + pb_x[k] * dh_0[k];

        t_1[k] = f_2 * dg_0[k]
                 + pb_z[k] * dh_1[k];

        t_2[k] = f_0 * dg_2[k]
                 + pb_z[k] * dh_3[k];

        t_3[k] = f_3 * dg_4[k]
                 + pb_z[k] * dh_5[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_y, pb_z, ph_1, pi_0, pi_4, dg_12, \
                         dg_13, dh_11, dh_12, dh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_x[k] * pi_0[k];

        t_5[k] = pa_x[k] * pi_4[k];

        t_6[k] = f_0 * ph_1[k]
                 + f_1 * dg_12[k]
                 + pb_y[k] * dh_11[k];

        t_7[k] = f_2 * dg_12[k]
                 + pb_z[k] * dh_12[k];

        t_8[k] = f_0 * dg_13[k]
                 + pb_z[k] * dh_13[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_z, ph_2, ph_3, ph_4, pi_1, \
                         pi_2, pi_3, pi_4, dg_14, dh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * dg_14[k]
                 + pb_z[k] * dh_14[k];

        t_10[k] = f_4 * ph_2[k]
                  + pa_y[k] * pi_1[k];

        t_11[k] = f_3 * ph_3[k]
                  + pa_y[k] * pi_2[k];

        t_12[k] = f_0 * ph_4[k]
                  + pa_y[k] * pi_3[k];

        t_13[k] = pa_y[k] * pi_4[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_y, dg_24, dg_25, dg_26, dg_27, dh_23, \
                         dh_24, dh_25, dh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * dg_24[k]
                  + pb_y[k] * dh_23[k];

        t_15[k] = f_3 * dg_25[k]
                  + pb_y[k] * dh_24[k];

        t_16[k] = f_0 * dg_26[k]
                  + pb_y[k] * dh_25[k];

        t_17[k] = f_2 * dg_27[k]
                  + pb_y[k] * dh_26[k];
    }

#pragma omp simd aligned(t_18, pb_z, ph_5, dg_27, dh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * ph_5[k]
                  + f_1 * dg_27[k]
                  + pb_z[k] * dh_27[k];
    }
}

auto
compute_prim_di_overlap_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t ph, const size_t pi, const size_t dg,
                          const size_t dh, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_z, ph_0, dg_0, dg_1, dg_2, dh_0, dh_1, \
                         dh_2, dh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ph_0[k]
                 + f_1 * dg_0[k]
                 + pb_x[k] * dh_0[k];

        t_1[k] = f_2 * dg_0[k]
                 + pb_z[k] * dh_1[k];

        t_2[k] = f_0 * dg_1[k]
                 + pb_z[k] * dh_2[k];

        t_3[k] = f_3 * dg_2[k]
                 + pb_z[k] * dh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_y, pb_z, ph_1, pi_0, pi_4, dg_8, \
                         dg_9, dh_8, dh_9, dh_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_x[k] * pi_0[k];

        t_5[k] = pa_x[k] * pi_4[k];

        t_6[k] = f_0 * ph_1[k]
                 + f_1 * dg_8[k]
                 + pb_y[k] * dh_8[k];

        t_7[k] = f_2 * dg_8[k]
                 + pb_z[k] * dh_9[k];

        t_8[k] = f_0 * dg_9[k]
                 + pb_z[k] * dh_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_z, ph_2, ph_3, ph_4, pi_1, \
                         pi_2, pi_3, pi_4, dg_10, dh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * dg_10[k]
                 + pb_z[k] * dh_11[k];

        t_10[k] = f_4 * ph_2[k]
                  + pa_y[k] * pi_1[k];

        t_11[k] = f_3 * ph_3[k]
                  + pa_y[k] * pi_2[k];

        t_12[k] = f_0 * ph_4[k]
                  + pa_y[k] * pi_3[k];

        t_13[k] = pa_y[k] * pi_4[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_y, dg_15, dg_16, dg_17, dg_18, dh_17, \
                         dh_18, dh_19, dh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * dg_15[k]
                  + pb_y[k] * dh_17[k];

        t_15[k] = f_3 * dg_16[k]
                  + pb_y[k] * dh_18[k];

        t_16[k] = f_0 * dg_17[k]
                  + pb_y[k] * dh_19[k];

        t_17[k] = f_2 * dg_18[k]
                  + pb_y[k] * dh_20[k];
    }

#pragma omp simd aligned(t_18, pb_z, ph_5, dg_18, dh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * ph_5[k]
                  + f_1 * dg_18[k]
                  + pb_z[k] * dh_21[k];
    }
}

}  // namespace simdovl
