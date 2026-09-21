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


#include "SimdElectronRepulsionVrrRecMP.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_mp_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ls, const size_t lp,
                                     const size_t ms, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 3.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 3.0 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 2.5 / p;
    const auto f_7 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ls_0 = buffer.data(ls + 0);
    const auto *ls_1 = buffer.data(ls + 1);
    const auto *ls_2 = buffer.data(ls + 2);
    const auto *ls_3 = buffer.data(ls + 3);
    const auto *ls_5 = buffer.data(ls + 5);
    const auto *ls_6 = buffer.data(ls + 6);
    const auto *ls_7 = buffer.data(ls + 7);
    const auto *ls_8 = buffer.data(ls + 8);
    const auto *ls_9 = buffer.data(ls + 9);
    const auto *ls_10 = buffer.data(ls + 10);
    const auto *ls_11 = buffer.data(ls + 11);
    const auto *ls_12 = buffer.data(ls + 12);
    const auto *ls_13 = buffer.data(ls + 13);
    const auto *ls_14 = buffer.data(ls + 14);
    const auto *ls_15 = buffer.data(ls + 15);
    const auto *ls_16 = buffer.data(ls + 16);
    const auto *ls_17 = buffer.data(ls + 17);
    const auto *ls_18 = buffer.data(ls + 18);
    const auto *ls_19 = buffer.data(ls + 19);
    const auto *ls_20 = buffer.data(ls + 20);
    const auto *ls_21 = buffer.data(ls + 21);
    const auto *ls_22 = buffer.data(ls + 22);
    const auto *ls_23 = buffer.data(ls + 23);
    const auto *ls_24 = buffer.data(ls + 24);
    const auto *ls_25 = buffer.data(ls + 25);
    const auto *ls_26 = buffer.data(ls + 26);
    const auto *ls_27 = buffer.data(ls + 27);
    const auto *ls_28 = buffer.data(ls + 28);
    const auto *ls_30 = buffer.data(ls + 30);
    const auto *ls_31 = buffer.data(ls + 31);
    const auto *ls_32 = buffer.data(ls + 32);
    const auto *ls_33 = buffer.data(ls + 33);
    const auto *ls_35 = buffer.data(ls + 35);
    const auto *ls_36 = buffer.data(ls + 36);
    const auto *ls_37 = buffer.data(ls + 37);
    const auto *ls_38 = buffer.data(ls + 38);
    const auto *ls_39 = buffer.data(ls + 39);
    const auto *ls_40 = buffer.data(ls + 40);
    const auto *ls_41 = buffer.data(ls + 41);
    const auto *ls_42 = buffer.data(ls + 42);
    const auto *ls_43 = buffer.data(ls + 43);
    const auto *ls_44 = buffer.data(ls + 44);

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_27 = buffer.data(lp + 27);
    const auto *lp_29 = buffer.data(lp + 29);
    const auto *lp_30 = buffer.data(lp + 30);
    const auto *lp_31 = buffer.data(lp + 31);
    const auto *lp_42 = buffer.data(lp + 42);
    const auto *lp_44 = buffer.data(lp + 44);
    const auto *lp_45 = buffer.data(lp + 45);
    const auto *lp_46 = buffer.data(lp + 46);
    const auto *lp_60 = buffer.data(lp + 60);
    const auto *lp_62 = buffer.data(lp + 62);
    const auto *lp_63 = buffer.data(lp + 63);
    const auto *lp_64 = buffer.data(lp + 64);
    const auto *lp_81 = buffer.data(lp + 81);
    const auto *lp_83 = buffer.data(lp + 83);
    const auto *lp_84 = buffer.data(lp + 84);
    const auto *lp_105 = buffer.data(lp + 105);
    const auto *lp_109 = buffer.data(lp + 109);
    const auto *lp_112 = buffer.data(lp + 112);
    const auto *lp_113 = buffer.data(lp + 113);
    const auto *lp_115 = buffer.data(lp + 115);
    const auto *lp_116 = buffer.data(lp + 116);
    const auto *lp_118 = buffer.data(lp + 118);
    const auto *lp_119 = buffer.data(lp + 119);
    const auto *lp_121 = buffer.data(lp + 121);
    const auto *lp_122 = buffer.data(lp + 122);
    const auto *lp_124 = buffer.data(lp + 124);
    const auto *lp_125 = buffer.data(lp + 125);
    const auto *lp_127 = buffer.data(lp + 127);
    const auto *lp_128 = buffer.data(lp + 128);
    const auto *lp_130 = buffer.data(lp + 130);
    const auto *lp_131 = buffer.data(lp + 131);
    const auto *lp_134 = buffer.data(lp + 134);

    const auto *ms_0 = buffer.data(ms + 0);
    const auto *ms_1 = buffer.data(ms + 1);
    const auto *ms_2 = buffer.data(ms + 2);
    const auto *ms_3 = buffer.data(ms + 3);
    const auto *ms_5 = buffer.data(ms + 5);
    const auto *ms_6 = buffer.data(ms + 6);
    const auto *ms_7 = buffer.data(ms + 7);
    const auto *ms_8 = buffer.data(ms + 8);
    const auto *ms_9 = buffer.data(ms + 9);
    const auto *ms_10 = buffer.data(ms + 10);
    const auto *ms_11 = buffer.data(ms + 11);
    const auto *ms_12 = buffer.data(ms + 12);
    const auto *ms_13 = buffer.data(ms + 13);
    const auto *ms_14 = buffer.data(ms + 14);
    const auto *ms_15 = buffer.data(ms + 15);
    const auto *ms_16 = buffer.data(ms + 16);
    const auto *ms_17 = buffer.data(ms + 17);
    const auto *ms_18 = buffer.data(ms + 18);
    const auto *ms_19 = buffer.data(ms + 19);
    const auto *ms_20 = buffer.data(ms + 20);
    const auto *ms_21 = buffer.data(ms + 21);
    const auto *ms_22 = buffer.data(ms + 22);
    const auto *ms_23 = buffer.data(ms + 23);
    const auto *ms_24 = buffer.data(ms + 24);
    const auto *ms_25 = buffer.data(ms + 25);
    const auto *ms_26 = buffer.data(ms + 26);
    const auto *ms_27 = buffer.data(ms + 27);
    const auto *ms_28 = buffer.data(ms + 28);
    const auto *ms_29 = buffer.data(ms + 29);
    const auto *ms_30 = buffer.data(ms + 30);
    const auto *ms_31 = buffer.data(ms + 31);
    const auto *ms_32 = buffer.data(ms + 32);
    const auto *ms_33 = buffer.data(ms + 33);
    const auto *ms_34 = buffer.data(ms + 34);
    const auto *ms_35 = buffer.data(ms + 35);
    const auto *ms_36 = buffer.data(ms + 36);
    const auto *ms_38 = buffer.data(ms + 38);
    const auto *ms_39 = buffer.data(ms + 39);
    const auto *ms_40 = buffer.data(ms + 40);
    const auto *ms_41 = buffer.data(ms + 41);
    const auto *ms_42 = buffer.data(ms + 42);
    const auto *ms_44 = buffer.data(ms + 44);
    const auto *ms_45 = buffer.data(ms + 45);
    const auto *ms_46 = buffer.data(ms + 46);
    const auto *ms_47 = buffer.data(ms + 47);
    const auto *ms_48 = buffer.data(ms + 48);
    const auto *ms_49 = buffer.data(ms + 49);
    const auto *ms_50 = buffer.data(ms + 50);
    const auto *ms_51 = buffer.data(ms + 51);
    const auto *ms_52 = buffer.data(ms + 52);
    const auto *ms_53 = buffer.data(ms + 53);
    const auto *ms_54 = buffer.data(ms + 54);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pa_y, pa_z, pb_x, pb_y, pb_z, \
                         ls_0, lp_0, ms_0, ms_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ls_0[k]
                 + pb_x[k] * ms_0[k];

        t_1[k] = pb_y[k] * ms_0[k];

        t_2[k] = pb_z[k] * ms_0[k];

        t_3[k] = pa_y[k] * lp_0[k];

        t_4[k] = f_1 * ls_0[k]
                 + pb_y[k] * ms_1[k];

        t_5[k] = pb_z[k] * ms_1[k];

        t_6[k] = pa_z[k] * lp_0[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pa_y, pb_x, pb_y, pb_z, ls_0, ls_1, \
                         ls_3, lp_6, ms_2, ms_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_y[k] * ms_2[k];

        t_8[k] = f_1 * ls_0[k]
                 + pb_z[k] * ms_2[k];

        t_9[k] = f_2 * ls_3[k]
                 + pb_x[k] * ms_3[k];

        t_10[k] = f_3 * ls_1[k]
                  + pb_y[k] * ms_3[k];

        t_11[k] = pb_z[k] * ms_3[k];

        t_12[k] = pa_y[k] * lp_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_y, pa_z, pb_x, pb_y, pb_z, ls_2, \
                         ls_5, lp_4, lp_8, ms_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_z[k] * lp_4[k];

        t_14[k] = pa_y[k] * lp_8[k];

        t_15[k] = f_2 * ls_5[k]
                  + pb_x[k] * ms_5[k];

        t_16[k] = pb_y[k] * ms_5[k];

        t_17[k] = f_3 * ls_2[k]
                  + pb_z[k] * ms_5[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, pa_z, pb_x, pb_y, pb_z, ls_3, \
                         ls_6, lp_9, lp_10, ms_6, ms_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_4 * ls_6[k]
                  + pb_x[k] * ms_6[k];

        t_19[k] = f_5 * ls_3[k]
                  + pb_y[k] * ms_6[k];

        t_20[k] = pb_z[k] * ms_6[k];

        t_21[k] = pa_z[k] * lp_9[k];

        t_22[k] = pa_z[k] * lp_10[k];

        t_23[k] = f_1 * ls_3[k]
                  + pb_z[k] * ms_7[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, pa_y, pb_x, pb_y, pb_z, ls_5, \
                         ls_9, lp_15, lp_17, ms_8, ms_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * lp_15[k];

        t_25[k] = f_1 * ls_5[k]
                  + pb_y[k] * ms_8[k];

        t_26[k] = pa_y[k] * lp_17[k];

        t_27[k] = f_4 * ls_9[k]
                  + pb_x[k] * ms_9[k];

        t_28[k] = pb_y[k] * ms_9[k];

        t_29[k] = f_5 * ls_5[k]
                  + pb_z[k] * ms_9[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pa_z, pb_x, pb_y, pb_z, ls_6, \
                         ls_10, lp_18, lp_19, ms_10, ms_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_6 * ls_10[k]
                  + pb_x[k] * ms_10[k];

        t_31[k] = f_7 * ls_6[k]
                  + pb_y[k] * ms_10[k];

        t_32[k] = pb_z[k] * ms_10[k];

        t_33[k] = pa_z[k] * lp_18[k];

        t_34[k] = pa_z[k] * lp_19[k];

        t_35[k] = f_1 * ls_6[k]
                  + pb_z[k] * ms_11[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_y, pb_x, pb_y, pb_z, ls_7, ls_8, \
                         ls_9, ls_12, lp_27, ms_12, ms_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_6 * ls_12[k]
                  + pb_x[k] * ms_12[k];

        t_37[k] = f_3 * ls_8[k]
                  + pb_y[k] * ms_12[k];

        t_38[k] = f_3 * ls_7[k]
                  + pb_z[k] * ms_12[k];

        t_39[k] = pa_y[k] * lp_27[k];

        t_40[k] = f_1 * ls_9[k]
                  + pb_y[k] * ms_13[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_y, pb_x, pb_y, pb_z, ls_9, ls_14, \
                         ls_15, lp_29, ms_14, ms_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pa_y[k] * lp_29[k];

        t_42[k] = f_6 * ls_14[k]
                  + pb_x[k] * ms_14[k];

        t_43[k] = pb_y[k] * ms_14[k];

        t_44[k] = f_7 * ls_9[k]
                  + pb_z[k] * ms_14[k];

        t_45[k] = f_7 * ls_15[k]
                  + pb_x[k] * ms_15[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_z, pb_y, pb_z, ls_10, lp_30, lp_31, \
                         ms_15, ms_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_6 * ls_10[k]
                  + pb_y[k] * ms_15[k];

        t_47[k] = pb_z[k] * ms_15[k];

        t_48[k] = pa_z[k] * lp_30[k];

        t_49[k] = pa_z[k] * lp_31[k];

        t_50[k] = f_1 * ls_10[k]
                  + pb_z[k] * ms_16[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, pb_x, pb_y, pb_z, ls_11, ls_12, \
                         ls_13, ls_17, ls_18, ms_17, ms_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_7 * ls_17[k]
                  + pb_x[k] * ms_17[k];

        t_52[k] = f_5 * ls_12[k]
                  + pb_y[k] * ms_17[k];

        t_53[k] = f_3 * ls_11[k]
                  + pb_z[k] * ms_17[k];

        t_54[k] = f_7 * ls_18[k]
                  + pb_x[k] * ms_18[k];

        t_55[k] = f_3 * ls_13[k]
                  + pb_y[k] * ms_18[k];

        t_56[k] = f_5 * ls_12[k]
                  + pb_z[k] * ms_18[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, pb_z, ls_14, \
                         ls_20, lp_42, lp_44, ms_19, ms_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_y[k] * lp_42[k];

        t_58[k] = f_1 * ls_14[k]
                  + pb_y[k] * ms_19[k];

        t_59[k] = pa_y[k] * lp_44[k];

        t_60[k] = f_7 * ls_20[k]
                  + pb_x[k] * ms_20[k];

        t_61[k] = pb_y[k] * ms_20[k];

        t_62[k] = f_6 * ls_14[k]
                  + pb_z[k] * ms_20[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, t_68, pa_z, pb_x, pb_y, pb_z, ls_15, \
                         ls_21, lp_45, lp_46, ms_21, ms_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_5 * ls_21[k]
                  + pb_x[k] * ms_21[k];

        t_64[k] = f_4 * ls_15[k]
                  + pb_y[k] * ms_21[k];

        t_65[k] = pb_z[k] * ms_21[k];

        t_66[k] = pa_z[k] * lp_45[k];

        t_67[k] = pa_z[k] * lp_46[k];

        t_68[k] = f_1 * ls_15[k]
                  + pb_z[k] * ms_22[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, t_74, pb_x, pb_y, pb_z, ls_16, ls_17, \
                         ls_18, ls_23, ls_24, ms_23, ms_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_5 * ls_23[k]
                  + pb_x[k] * ms_23[k];

        t_70[k] = f_7 * ls_17[k]
                  + pb_y[k] * ms_23[k];

        t_71[k] = f_3 * ls_16[k]
                  + pb_z[k] * ms_23[k];

        t_72[k] = f_5 * ls_24[k]
                  + pb_x[k] * ms_24[k];

        t_73[k] = f_5 * ls_18[k]
                  + pb_y[k] * ms_24[k];

        t_74[k] = f_5 * ls_17[k]
                  + pb_z[k] * ms_24[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, pb_z, ls_18, ls_19, \
                         ls_20, ls_25, lp_60, ms_25, ms_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_5 * ls_25[k]
                  + pb_x[k] * ms_25[k];

        t_76[k] = f_3 * ls_19[k]
                  + pb_y[k] * ms_25[k];

        t_77[k] = f_7 * ls_18[k]
                  + pb_z[k] * ms_25[k];

        t_78[k] = pa_y[k] * lp_60[k];

        t_79[k] = f_1 * ls_20[k]
                  + pb_y[k] * ms_26[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pa_y, pb_x, pb_y, pb_z, ls_20, ls_27, \
                         ls_28, lp_62, ms_27, ms_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_y[k] * lp_62[k];

        t_81[k] = f_5 * ls_27[k]
                  + pb_x[k] * ms_27[k];

        t_82[k] = pb_y[k] * ms_27[k];

        t_83[k] = f_4 * ls_20[k]
                  + pb_z[k] * ms_27[k];

        t_84[k] = f_3 * ls_28[k]
                  + pb_x[k] * ms_28[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_z, pb_y, pb_z, ls_21, lp_63, lp_64, \
                         ms_28, ms_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_2 * ls_21[k]
                  + pb_y[k] * ms_28[k];

        t_86[k] = pb_z[k] * ms_28[k];

        t_87[k] = pa_z[k] * lp_63[k];

        t_88[k] = pa_z[k] * lp_64[k];

        t_89[k] = f_1 * ls_21[k]
                  + pb_z[k] * ms_29[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, ls_22, ls_23, \
                         ls_24, ls_30, ls_31, ms_30, ms_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_3 * ls_30[k]
                  + pb_x[k] * ms_30[k];

        t_91[k] = f_6 * ls_23[k]
                  + pb_y[k] * ms_30[k];

        t_92[k] = f_3 * ls_22[k]
                  + pb_z[k] * ms_30[k];

        t_93[k] = f_3 * ls_31[k]
                  + pb_x[k] * ms_31[k];

        t_94[k] = f_7 * ls_24[k]
                  + pb_y[k] * ms_31[k];

        t_95[k] = f_5 * ls_23[k]
                  + pb_z[k] * ms_31[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, t_101, pb_x, pb_y, pb_z, ls_24, ls_25, \
                         ls_26, ls_32, ls_33, ms_32, ms_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_3 * ls_32[k]
                  + pb_x[k] * ms_32[k];

        t_97[k] = f_5 * ls_25[k]
                  + pb_y[k] * ms_32[k];

        t_98[k] = f_7 * ls_24[k]
                  + pb_z[k] * ms_32[k];

        t_99[k] = f_3 * ls_33[k]
                  + pb_x[k] * ms_33[k];

        t_100[k] = f_3 * ls_26[k]
                   + pb_y[k] * ms_33[k];

        t_101[k] = f_6 * ls_25[k]
                   + pb_z[k] * ms_33[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, t_107, pa_y, pb_x, pb_y, pb_z, \
                         ls_27, ls_35, lp_81, lp_83, ms_34, ms_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = pa_y[k] * lp_81[k];

        t_103[k] = f_1 * ls_27[k]
                   + pb_y[k] * ms_34[k];

        t_104[k] = pa_y[k] * lp_83[k];

        t_105[k] = f_3 * ls_35[k]
                   + pb_x[k] * ms_35[k];

        t_106[k] = pb_y[k] * ms_35[k];

        t_107[k] = f_2 * ls_27[k]
                   + pb_z[k] * ms_35[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, t_113, pa_x, pa_z, pb_x, pb_z, \
                         ls_36, lp_84, lp_109, lp_112, lp_113, ms_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_1 * ls_36[k]
                   + pb_x[k] * ms_36[k];

        t_109[k] = pa_x[k] * lp_109[k];

        t_110[k] = pb_z[k] * ms_36[k];

        t_111[k] = pa_z[k] * lp_84[k];

        t_112[k] = pa_x[k] * lp_112[k];

        t_113[k] = pa_x[k] * lp_113[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, pa_x, pb_x, ls_38, ls_39, \
                         lp_115, lp_116, lp_118, lp_119, ms_38, ms_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_1 * ls_38[k]
                   + pb_x[k] * ms_38[k];

        t_115[k] = pa_x[k] * lp_115[k];

        t_116[k] = pa_x[k] * lp_116[k];

        t_117[k] = f_1 * ls_39[k]
                   + pb_x[k] * ms_39[k];

        t_118[k] = pa_x[k] * lp_118[k];

        t_119[k] = pa_x[k] * lp_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, pa_x, pb_x, ls_40, ls_41, \
                         lp_121, lp_122, lp_124, lp_125, ms_40, ms_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_1 * ls_40[k]
                   + pb_x[k] * ms_40[k];

        t_121[k] = pa_x[k] * lp_121[k];

        t_122[k] = pa_x[k] * lp_122[k];

        t_123[k] = f_1 * ls_41[k]
                   + pb_x[k] * ms_41[k];

        t_124[k] = pa_x[k] * lp_124[k];

        t_125[k] = pa_x[k] * lp_125[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, t_131, pa_x, pa_y, pb_x, ls_42, \
                         lp_105, lp_127, lp_128, lp_130, lp_131, \
                         ms_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_1 * ls_42[k]
                   + pb_x[k] * ms_42[k];

        t_127[k] = pa_x[k] * lp_127[k];

        t_128[k] = pa_x[k] * lp_128[k];

        t_129[k] = pa_y[k] * lp_105[k];

        t_130[k] = pa_x[k] * lp_130[k];

        t_131[k] = pa_x[k] * lp_131[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, t_137, pa_x, pb_x, pb_y, pb_z, \
                         ls_36, ls_44, lp_134, ms_44, ms_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_1 * ls_44[k]
                   + pb_x[k] * ms_44[k];

        t_133[k] = pb_y[k] * ms_44[k];

        t_134[k] = pa_x[k] * lp_134[k];

        t_135[k] = pb_x[k] * ms_45[k];

        t_136[k] = f_0 * ls_36[k]
                   + pb_y[k] * ms_45[k];

        t_137[k] = pb_z[k] * ms_45[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, t_143, pa_z, pb_x, pb_y, pb_z, \
                         ls_36, ls_37, ls_38, lp_109, ms_46, ms_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = pb_x[k] * ms_46[k];

        t_139[k] = pa_z[k] * lp_109[k];

        t_140[k] = f_1 * ls_36[k]
                   + pb_z[k] * ms_46[k];

        t_141[k] = pb_x[k] * ms_47[k];

        t_142[k] = f_2 * ls_38[k]
                   + pb_y[k] * ms_47[k];

        t_143[k] = f_3 * ls_37[k]
                   + pb_z[k] * ms_47[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, t_149, t_150, pb_x, pb_y, pb_z, \
                         ls_38, ls_39, ls_40, ms_48, ms_49, ms_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pb_x[k] * ms_48[k];

        t_145[k] = f_4 * ls_39[k]
                   + pb_y[k] * ms_48[k];

        t_146[k] = f_5 * ls_38[k]
                   + pb_z[k] * ms_48[k];

        t_147[k] = pb_x[k] * ms_49[k];

        t_148[k] = f_6 * ls_40[k]
                   + pb_y[k] * ms_49[k];

        t_149[k] = f_7 * ls_39[k]
                   + pb_z[k] * ms_49[k];

        t_150[k] = pb_x[k] * ms_50[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, t_156, pb_x, pb_y, pb_z, ls_40, \
                         ls_41, ls_42, ms_50, ms_51, ms_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_7 * ls_41[k]
                   + pb_y[k] * ms_50[k];

        t_152[k] = f_6 * ls_40[k]
                   + pb_z[k] * ms_50[k];

        t_153[k] = pb_x[k] * ms_51[k];

        t_154[k] = f_5 * ls_42[k]
                   + pb_y[k] * ms_51[k];

        t_155[k] = f_4 * ls_41[k]
                   + pb_z[k] * ms_51[k];

        t_156[k] = pb_x[k] * ms_52[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pa_y, pb_x, pb_y, pb_z, ls_42, \
                         ls_43, ls_44, lp_134, ms_52, ms_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_3 * ls_43[k]
                   + pb_y[k] * ms_52[k];

        t_158[k] = f_2 * ls_42[k]
                   + pb_z[k] * ms_52[k];

        t_159[k] = pb_x[k] * ms_53[k];

        t_160[k] = f_1 * ls_44[k]
                   + pb_y[k] * ms_53[k];

        t_161[k] = pa_y[k] * lp_134[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pb_x, pb_y, pb_z, ls_44, \
                         ms_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pb_x[k] * ms_54[k];

        t_163[k] = pb_y[k] * ms_54[k];

        t_164[k] = f_0 * ls_44[k]
                   + pb_z[k] * ms_54[k];
    }
}

}  // namespace simdt2ceri
