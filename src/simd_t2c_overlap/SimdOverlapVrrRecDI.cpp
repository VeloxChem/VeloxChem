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

static auto
compute_prim_di_overlap_0_piece0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t ph, const size_t pi,
                                 const size_t dg, const size_t dh, const size_t ncols,
                                 const double p) -> void
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_9 = buffer.data(ph + 9);
    const auto *ph_15 = buffer.data(ph + 15);
    const auto *ph_17 = buffer.data(ph + 17);
    const auto *ph_18 = buffer.data(ph + 18);
    const auto *ph_20 = buffer.data(ph + 20);
    const auto *ph_24 = buffer.data(ph + 24);
    const auto *ph_27 = buffer.data(ph + 27);
    const auto *ph_31 = buffer.data(ph + 31);
    const auto *ph_36 = buffer.data(ph + 36);
    const auto *ph_38 = buffer.data(ph + 38);
    const auto *ph_39 = buffer.data(ph + 39);
    const auto *ph_40 = buffer.data(ph + 40);
    const auto *ph_41 = buffer.data(ph + 41);
    const auto *ph_47 = buffer.data(ph + 47);
    const auto *ph_51 = buffer.data(ph + 51);
    const auto *ph_56 = buffer.data(ph + 56);
    const auto *ph_58 = buffer.data(ph + 58);
    const auto *ph_59 = buffer.data(ph + 59);
    const auto *ph_60 = buffer.data(ph + 60);
    const auto *ph_61 = buffer.data(ph + 61);
    const auto *ph_62 = buffer.data(ph + 62);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_29 = buffer.data(pi + 29);
    const auto *pi_31 = buffer.data(pi + 31);
    const auto *pi_34 = buffer.data(pi + 34);
    const auto *pi_38 = buffer.data(pi + 38);
    const auto *pi_49 = buffer.data(pi + 49);
    const auto *pi_51 = buffer.data(pi + 51);
    const auto *pi_52 = buffer.data(pi + 52);
    const auto *pi_53 = buffer.data(pi + 53);
    const auto *pi_54 = buffer.data(pi + 54);
    const auto *pi_55 = buffer.data(pi + 55);
    const auto *pi_56 = buffer.data(pi + 56);
    const auto *pi_58 = buffer.data(pi + 58);
    const auto *pi_61 = buffer.data(pi + 61);
    const auto *pi_65 = buffer.data(pi + 65);
    const auto *pi_70 = buffer.data(pi + 70);
    const auto *pi_77 = buffer.data(pi + 77);
    const auto *pi_78 = buffer.data(pi + 78);
    const auto *pi_79 = buffer.data(pi + 79);
    const auto *pi_80 = buffer.data(pi + 80);
    const auto *pi_81 = buffer.data(pi + 81);
    const auto *pi_83 = buffer.data(pi + 83);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_53 = buffer.data(dg + 53);
    const auto *dg_54 = buffer.data(dg + 54);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_56 = buffer.data(dg + 56);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_58 = buffer.data(dg + 58);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_64 = buffer.data(dg + 64);
    const auto *dg_67 = buffer.data(dg + 67);
    const auto *dg_68 = buffer.data(dg + 68);
    const auto *dg_71 = buffer.data(dg + 71);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_73 = buffer.data(dg + 73);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_77 = buffer.data(dg + 77);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_80 = buffer.data(dg + 80);
    const auto *dg_81 = buffer.data(dg + 81);
    const auto *dg_82 = buffer.data(dg + 82);
    const auto *dg_84 = buffer.data(dg + 84);
    const auto *dg_85 = buffer.data(dg + 85);
    const auto *dg_86 = buffer.data(dg + 86);
    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_89 = buffer.data(dg + 89);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_58 = buffer.data(dh + 58);
    const auto *dh_59 = buffer.data(dh + 59);
    const auto *dh_60 = buffer.data(dh + 60);
    const auto *dh_62 = buffer.data(dh + 62);
    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_64 = buffer.data(dh + 64);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_68 = buffer.data(dh + 68);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_71 = buffer.data(dh + 71);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_73 = buffer.data(dh + 73);
    const auto *dh_75 = buffer.data(dh + 75);
    const auto *dh_76 = buffer.data(dh + 76);
    const auto *dh_77 = buffer.data(dh + 77);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_79 = buffer.data(dh + 79);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_82 = buffer.data(dh + 82);
    const auto *dh_83 = buffer.data(dh + 83);
    const auto *dh_88 = buffer.data(dh + 88);
    const auto *dh_91 = buffer.data(dh + 91);
    const auto *dh_92 = buffer.data(dh + 92);
    const auto *dh_95 = buffer.data(dh + 95);
    const auto *dh_96 = buffer.data(dh + 96);
    const auto *dh_97 = buffer.data(dh + 97);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_104 = buffer.data(dh + 104);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_107 = buffer.data(dh + 107);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_111 = buffer.data(dh + 111);
    const auto *dh_112 = buffer.data(dh + 112);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_115 = buffer.data(dh + 115);
    const auto *dh_116 = buffer.data(dh + 116);
    const auto *dh_117 = buffer.data(dh + 117);
    const auto *dh_119 = buffer.data(dh + 119);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_121 = buffer.data(dh + 121);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_124 = buffer.data(dh + 124);
    const auto *dh_125 = buffer.data(dh + 125);

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
                         dh_5, dh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * dg_1[k]
                 + pb_y[k] * dh_3[k];

        t_7[k] = pb_z[k] * dh_3[k];

        t_8[k] = pb_y[k] * dh_5[k];

        t_9[k] = f_0 * dg_2[k]
                 + pb_z[k] * dh_5[k];

        t_10[k] = f_3 * dg_3[k]
                  + pb_y[k] * dh_6[k];

        t_11[k] = pb_z[k] * dh_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, pb_z, ph_15, dg_5, dh_8, \
                         dh_9, dh_10, dh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * dg_5[k]
                  + pb_y[k] * dh_8[k];

        t_13[k] = pb_y[k] * dh_9[k];

        t_14[k] = f_3 * dg_5[k]
                  + pb_z[k] * dh_9[k];

        t_15[k] = f_0 * ph_15[k]
                  + pb_x[k] * dh_15[k];

        t_16[k] = pb_z[k] * dh_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pb_x, pb_y, ph_17, ph_18, ph_20, dg_10, \
                         dh_14, dh_15, dh_17, dh_18, dh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * ph_17[k]
                  + pb_x[k] * dh_17[k];

        t_18[k] = f_0 * ph_18[k]
                  + pb_x[k] * dh_18[k];

        t_19[k] = pb_y[k] * dh_14[k];

        t_20[k] = f_0 * ph_20[k]
                  + pb_x[k] * dh_20[k];

        t_21[k] = f_1 * dg_10[k]
                  + pb_y[k] * dh_15[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pb_y, pb_z, dg_12, dg_13, dg_14, \
                         dh_15, dh_17, dh_18, dh_19, dh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_z[k] * dh_15[k];

        t_23[k] = f_3 * dg_12[k]
                  + pb_y[k] * dh_17[k];

        t_24[k] = f_0 * dg_13[k]
                  + pb_y[k] * dh_18[k];

        t_25[k] = f_2 * dg_14[k]
                  + pb_y[k] * dh_19[k];

        t_26[k] = pb_y[k] * dh_20[k];

        t_27[k] = f_1 * dg_14[k]
                  + pb_z[k] * dh_20[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_x, pa_y, pb_y, pb_z, ph_0, ph_24, \
                         pi_0, pi_31, dh_21, dh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_y[k] * pi_0[k];

        t_29[k] = f_2 * ph_0[k]
                  + pb_y[k] * dh_21[k];

        t_30[k] = pb_z[k] * dh_21[k];

        t_31[k] = f_4 * ph_24[k]
                  + pa_x[k] * pi_31[k];

        t_32[k] = pb_z[k] * dh_22[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_x, pa_y, pb_y, pb_z, ph_5, ph_27, \
                         pi_5, pi_9, pi_34, dh_24, dh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_y[k] * pi_5[k];

        t_34[k] = f_3 * ph_27[k]
                  + pa_x[k] * pi_34[k];

        t_35[k] = pb_z[k] * dh_24[k];

        t_36[k] = f_2 * ph_5[k]
                  + pb_y[k] * dh_26[k];

        t_37[k] = pa_y[k] * pi_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pb_y, pb_z, ph_9, ph_31, pi_38, dg_18, \
                         dh_27, dh_28, dh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * ph_31[k]
                  + pa_x[k] * pi_38[k];

        t_39[k] = pb_z[k] * dh_27[k];

        t_40[k] = f_2 * dg_18[k]
                  + pb_z[k] * dh_28[k];

        t_41[k] = f_2 * ph_9[k]
                  + pb_y[k] * dh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_y, pb_x, pb_z, ph_36, ph_38, ph_39, \
                         pi_14, dh_31, dh_36, dh_38, dh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * pi_14[k];

        t_43[k] = f_2 * ph_36[k]
                  + pb_x[k] * dh_36[k];

        t_44[k] = pb_z[k] * dh_31[k];

        t_45[k] = f_2 * ph_38[k]
                  + pb_x[k] * dh_38[k];

        t_46[k] = f_2 * ph_39[k]
                  + pb_x[k] * dh_39[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pa_x, pa_y, pb_x, pb_z, ph_40, pi_20, \
                         pi_49, pi_51, dh_36, dh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_2 * ph_40[k]
                  + pb_x[k] * dh_40[k];

        t_48[k] = pa_y[k] * pi_20[k];

        t_49[k] = pa_x[k] * pi_49[k];

        t_50[k] = pb_z[k] * dh_36[k];

        t_51[k] = pa_x[k] * pi_51[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, pa_x, pa_z, pb_y, pi_0, pi_52, \
                         pi_53, pi_54, pi_55, dh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pa_x[k] * pi_52[k];

        t_53[k] = pa_x[k] * pi_53[k];

        t_54[k] = pa_x[k] * pi_54[k];

        t_55[k] = pa_x[k] * pi_55[k];

        t_56[k] = pa_z[k] * pi_0[k];

        t_57[k] = pb_y[k] * dh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_x, pa_z, pb_y, pb_z, ph_0, ph_47, \
                         pi_3, pi_6, pi_61, dh_42, dh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_2 * ph_0[k]
                  + pb_z[k] * dh_42[k];

        t_59[k] = pa_z[k] * pi_3[k];

        t_60[k] = pb_y[k] * dh_44[k];

        t_61[k] = f_4 * ph_47[k]
                  + pa_x[k] * pi_61[k];

        t_62[k] = pa_z[k] * pi_6[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_x, pa_z, pb_y, ph_51, pi_10, pi_65, \
                         dg_32, dg_34, dh_46, dh_47, dh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_2 * dg_32[k]
                  + pb_y[k] * dh_46[k];

        t_64[k] = pb_y[k] * dh_47[k];

        t_65[k] = f_3 * ph_51[k]
                  + pa_x[k] * pi_65[k];

        t_66[k] = pa_z[k] * pi_10[k];

        t_67[k] = f_0 * dg_34[k]
                  + pb_y[k] * dh_49[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_x, pa_z, pb_y, ph_56, pi_15, pi_70, dg_35, \
                         dh_50, dh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * dg_35[k]
                  + pb_y[k] * dh_50[k];

        t_69[k] = pb_y[k] * dh_51[k];

        t_70[k] = f_0 * ph_56[k]
                  + pa_x[k] * pi_70[k];

        t_71[k] = pa_z[k] * pi_15[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pb_x, pb_y, ph_58, ph_59, ph_60, ph_62, \
                         dh_56, dh_58, dh_59, dh_60, dh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_2 * ph_58[k]
                  + pb_x[k] * dh_58[k];

        t_73[k] = f_2 * ph_59[k]
                  + pb_x[k] * dh_59[k];

        t_74[k] = f_2 * ph_60[k]
                  + pb_x[k] * dh_60[k];

        t_75[k] = pb_y[k] * dh_56[k];

        t_76[k] = f_2 * ph_62[k]
                  + pb_x[k] * dh_62[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, t_83, pa_x, pb_y, pi_77, pi_78, \
                         pi_79, pi_80, pi_81, pi_83, dh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_x[k] * pi_77[k];

        t_78[k] = pa_x[k] * pi_78[k];

        t_79[k] = pa_x[k] * pi_79[k];

        t_80[k] = pa_x[k] * pi_80[k];

        t_81[k] = pa_x[k] * pi_81[k];

        t_82[k] = pb_y[k] * dh_62[k];

        t_83[k] = pa_x[k] * pi_83[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, pb_x, pb_z, dg_45, dg_46, dg_48, \
                         dg_50, dh_63, dh_64, dh_66, dh_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * dg_45[k]
                  + pb_x[k] * dh_63[k];

        t_85[k] = f_4 * dg_46[k]
                  + pb_x[k] * dh_64[k];

        t_86[k] = pb_z[k] * dh_63[k];

        t_87[k] = f_3 * dg_48[k]
                  + pb_x[k] * dh_66[k];

        t_88[k] = pb_z[k] * dh_64[k];

        t_89[k] = f_3 * dg_50[k]
                  + pb_x[k] * dh_68[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pb_x, pb_z, dg_51, dg_53, dg_54, dg_55, \
                         dh_66, dh_69, dh_71, dh_72, dh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * dg_51[k]
                  + pb_x[k] * dh_69[k];

        t_91[k] = pb_z[k] * dh_66[k];

        t_92[k] = f_0 * dg_53[k]
                  + pb_x[k] * dh_71[k];

        t_93[k] = f_0 * dg_54[k]
                  + pb_x[k] * dh_72[k];

        t_94[k] = f_2 * dg_55[k]
                  + pb_x[k] * dh_73[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pb_x, pb_z, dg_57, dg_58, dg_59, dh_69, \
                         dh_75, dh_76, dh_77, dh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pb_z[k] * dh_69[k];

        t_96[k] = f_2 * dg_57[k]
                  + pb_x[k] * dh_75[k];

        t_97[k] = f_2 * dg_58[k]
                  + pb_x[k] * dh_76[k];

        t_98[k] = f_2 * dg_59[k]
                  + pb_x[k] * dh_77[k];

        t_99[k] = pb_x[k] * dh_78[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, t_105, pb_x, pb_y, ph_36, dg_55, \
                         dh_78, dh_79, dh_80, dh_81, dh_82, dh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pb_x[k] * dh_79[k];

        t_101[k] = pb_x[k] * dh_80[k];

        t_102[k] = pb_x[k] * dh_81[k];

        t_103[k] = pb_x[k] * dh_82[k];

        t_104[k] = pb_x[k] * dh_83[k];

        t_105[k] = f_0 * ph_36[k]
                   + f_1 * dg_55[k]
                   + pb_y[k] * dh_78[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, pb_y, pb_z, ph_41, dg_55, dg_56, \
                         dg_57, dh_78, dh_79, dh_80, dh_81, dh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = pb_z[k] * dh_78[k];

        t_107[k] = f_2 * dg_55[k]
                   + pb_z[k] * dh_79[k];

        t_108[k] = f_0 * dg_56[k]
                   + pb_z[k] * dh_80[k];

        t_109[k] = f_3 * dg_57[k]
                   + pb_z[k] * dh_81[k];

        t_110[k] = f_0 * ph_41[k]
                   + pb_y[k] * dh_83[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, pa_y, pa_z, pb_z, pi_29, pi_31, \
                         pi_56, pi_58, dg_59, dh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_1 * dg_59[k]
                   + pb_z[k] * dh_83[k];

        t_112[k] = pa_y[k] * pi_56[k];

        t_113[k] = pa_z[k] * pi_29[k];

        t_114[k] = pa_y[k] * pi_58[k];

        t_115[k] = pa_z[k] * pi_31[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, pa_y, pa_z, pb_x, pi_34, pi_61, \
                         dg_64, dg_67, dg_68, dh_88, dh_91, dh_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_3 * dg_64[k]
                   + pb_x[k] * dh_88[k];

        t_117[k] = pa_y[k] * pi_61[k];

        t_118[k] = pa_z[k] * pi_34[k];

        t_119[k] = f_0 * dg_67[k]
                   + pb_x[k] * dh_91[k];

        t_120[k] = f_0 * dg_68[k]
                   + pb_x[k] * dh_92[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pa_y, pa_z, pb_x, pi_38, pi_65, \
                         dg_71, dg_72, dg_73, dh_95, dh_96, dh_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * pi_65[k];

        t_122[k] = pa_z[k] * pi_38[k];

        t_123[k] = f_2 * dg_71[k]
                   + pb_x[k] * dh_95[k];

        t_124[k] = f_2 * dg_72[k]
                   + pb_x[k] * dh_96[k];

        t_125[k] = f_2 * dg_73[k]
                   + pb_x[k] * dh_97[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, t_131, t_132, pa_y, pb_x, pi_70, \
                         dh_99, dh_100, dh_101, dh_102, dh_103, \
                         dh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_y[k] * pi_70[k];

        t_127[k] = pb_x[k] * dh_99[k];

        t_128[k] = pb_x[k] * dh_100[k];

        t_129[k] = pb_x[k] * dh_101[k];

        t_130[k] = pb_x[k] * dh_102[k];

        t_131[k] = pb_x[k] * dh_103[k];

        t_132[k] = pb_x[k] * dh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_y, pa_z, pb_z, ph_36, ph_59, ph_60, \
                         pi_49, pi_79, pi_80, dh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * pi_49[k];

        t_134[k] = f_2 * ph_36[k]
                   + pb_z[k] * dh_99[k];

        t_135[k] = f_4 * ph_59[k]
                   + pa_y[k] * pi_79[k];

        t_136[k] = f_3 * ph_60[k]
                   + pa_y[k] * pi_80[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, pa_y, pb_x, pb_y, ph_61, ph_62, \
                         pi_81, pi_83, dg_75, dh_104, dh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_0 * ph_61[k]
                   + pa_y[k] * pi_81[k];

        t_138[k] = f_2 * ph_62[k]
                   + pb_y[k] * dh_104[k];

        t_139[k] = pa_y[k] * pi_83[k];

        t_140[k] = f_1 * dg_75[k]
                   + pb_x[k] * dh_105[k];

        t_141[k] = pb_y[k] * dh_105[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, pb_x, pb_y, dg_77, dg_78, dg_80, \
                         dg_81, dh_107, dh_108, dh_110, dh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_4 * dg_77[k]
                   + pb_x[k] * dh_107[k];

        t_143[k] = f_3 * dg_78[k]
                   + pb_x[k] * dh_108[k];

        t_144[k] = pb_y[k] * dh_107[k];

        t_145[k] = f_3 * dg_80[k]
                   + pb_x[k] * dh_110[k];

        t_146[k] = f_0 * dg_81[k]
                   + pb_x[k] * dh_111[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, pb_x, pb_y, dg_82, dg_84, dg_85, \
                         dg_86, dh_110, dh_112, dh_114, dh_115, \
                         dh_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_0 * dg_82[k]
                   + pb_x[k] * dh_112[k];

        t_148[k] = pb_y[k] * dh_110[k];

        t_149[k] = f_0 * dg_84[k]
                   + pb_x[k] * dh_114[k];

        t_150[k] = f_2 * dg_85[k]
                   + pb_x[k] * dh_115[k];

        t_151[k] = f_2 * dg_86[k]
                   + pb_x[k] * dh_116[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, t_157, pb_x, pb_y, dg_87, dg_89, \
                         dh_114, dh_117, dh_119, dh_120, dh_121, \
                         dh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_2 * dg_87[k]
                   + pb_x[k] * dh_117[k];

        t_153[k] = pb_y[k] * dh_114[k];

        t_154[k] = f_2 * dg_89[k]
                   + pb_x[k] * dh_119[k];

        t_155[k] = pb_x[k] * dh_120[k];

        t_156[k] = pb_x[k] * dh_121[k];

        t_157[k] = pb_x[k] * dh_122[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pb_x, pb_y, dg_85, dg_86, dh_120, \
                         dh_121, dh_123, dh_124, dh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = pb_x[k] * dh_123[k];

        t_159[k] = pb_x[k] * dh_124[k];

        t_160[k] = pb_x[k] * dh_125[k];

        t_161[k] = f_1 * dg_85[k]
                   + pb_y[k] * dh_120[k];

        t_162[k] = f_4 * dg_86[k]
                   + pb_y[k] * dh_121[k];
    }
}

static auto
compute_prim_di_overlap_0_piece1(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t ph, const size_t dg, const size_t dh,
                                 const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;

    auto *t_163 = buffer.data(target + 163);
    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_62 = buffer.data(ph + 62);

    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_88 = buffer.data(dg + 88);
    const auto *dg_89 = buffer.data(dg + 89);

    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_124 = buffer.data(dh + 124);
    const auto *dh_125 = buffer.data(dh + 125);

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pb_y, pb_z, ph_62, dg_87, dg_88, \
                         dg_89, dh_122, dh_123, dh_124, dh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_3 * dg_87[k]
                   + pb_y[k] * dh_122[k];

        t_164[k] = f_0 * dg_88[k]
                   + pb_y[k] * dh_123[k];

        t_165[k] = f_2 * dg_89[k]
                   + pb_y[k] * dh_124[k];

        t_166[k] = pb_y[k] * dh_125[k];

        t_167[k] = f_0 * ph_62[k]
                   + f_1 * dg_89[k]
                   + pb_z[k] * dh_125[k];
    }
}

auto
compute_prim_di_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t ph, const size_t pi, const size_t dg,
                          const size_t dh, const size_t ncols, const double p) -> void
{
    compute_prim_di_overlap_0_piece0(buffer, target, pa, pb, ph, pi, dg, dh, ncols, p);

    compute_prim_di_overlap_0_piece1(buffer, target, pb, ph, dg, dh, ncols, p);
}

}  // namespace simdovl
