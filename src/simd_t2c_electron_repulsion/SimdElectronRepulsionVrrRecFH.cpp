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


#include "SimdElectronRepulsionVrrRecFH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_fh_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t dg,
                                            const size_t dh, const size_t ff0, const size_t ff1,
                                            const size_t fg, const size_t ncols,
                                            const double alpha, const double beta,
                                            const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_47 = buffer.data(dg + 47);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_54 = buffer.data(dg + 54);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_58 = buffer.data(dg + 58);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_62 = buffer.data(dg + 62);
    const auto *dg_71 = buffer.data(dg + 71);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_73 = buffer.data(dg + 73);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_80 = buffer.data(dg + 80);
    const auto *dg_81 = buffer.data(dg + 81);
    const auto *dg_84 = buffer.data(dg + 84);
    const auto *dg_85 = buffer.data(dg + 85);
    const auto *dg_86 = buffer.data(dg + 86);
    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_89 = buffer.data(dg + 89);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_64 = buffer.data(dh + 64);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_68 = buffer.data(dh + 68);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_82 = buffer.data(dh + 82);
    const auto *dh_83 = buffer.data(dh + 83);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_104 = buffer.data(dh + 104);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_111 = buffer.data(dh + 111);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_121 = buffer.data(dh + 121);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_125 = buffer.data(dh + 125);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_9 = buffer.data(ff0 + 9);
    const auto *ff0_60 = buffer.data(ff0 + 60);
    const auto *ff0_63 = buffer.data(ff0 + 63);
    const auto *ff0_65 = buffer.data(ff0 + 65);
    const auto *ff0_66 = buffer.data(ff0 + 66);
    const auto *ff0_67 = buffer.data(ff0 + 67);
    const auto *ff0_69 = buffer.data(ff0 + 69);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_9 = buffer.data(ff1 + 9);
    const auto *ff1_60 = buffer.data(ff1 + 60);
    const auto *ff1_63 = buffer.data(ff1 + 63);
    const auto *ff1_65 = buffer.data(ff1 + 65);
    const auto *ff1_66 = buffer.data(ff1 + 66);
    const auto *ff1_67 = buffer.data(ff1 + 67);
    const auto *ff1_69 = buffer.data(ff1 + 69);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_78 = buffer.data(fg + 78);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_85 = buffer.data(fg + 85);
    const auto *fg_86 = buffer.data(fg + 86);
    const auto *fg_87 = buffer.data(fg + 87);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_95 = buffer.data(fg + 95);
    const auto *fg_96 = buffer.data(fg + 96);
    const auto *fg_99 = buffer.data(fg + 99);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_101 = buffer.data(fg + 101);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_103 = buffer.data(fg + 103);
    const auto *fg_104 = buffer.data(fg + 104);
    const auto *fg_105 = buffer.data(fg + 105);
    const auto *fg_107 = buffer.data(fg + 107);
    const auto *fg_108 = buffer.data(fg + 108);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, \
                         fg_0, fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = pb_y[k] * fg_2[k];

        t_5[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, dg_10, ff0_1, ff0_2, \
                         ff1_1, ff1_2, fg_3, fg_5, fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_7[k] = pb_z[k] * fg_3[k];

        t_8[k] = pb_y[k] * fg_5[k];

        t_9[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_5[k];

        t_10[k] = f_0 * dg_10[k]
                  + pb_x[k] * fg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pb_x, pb_y, pb_z, dg_12, dg_14, fg_6, fg_9, \
                         fg_12, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * fg_6[k];

        t_12[k] = f_0 * dg_12[k]
                  + pb_x[k] * fg_12[k];

        t_13[k] = pb_y[k] * fg_9[k];

        t_14[k] = f_0 * dg_14[k]
                  + pb_x[k] * fg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pb_z, ff0_6, ff0_8, ff0_9, ff1_6, \
                         ff1_8, ff1_9, fg_10, fg_12, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * ff0_6[k]
                  - f_2 * ff1_6[k]
                  + pb_y[k] * fg_10[k];

        t_16[k] = pb_z[k] * fg_10[k];

        t_17[k] = f_5 * ff0_8[k]
                  - f_6 * ff1_8[k]
                  + pb_y[k] * fg_12[k];

        t_18[k] = f_3 * ff0_9[k]
                  - f_4 * ff1_9[k]
                  + pb_y[k] * fg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pb_y, pb_z, dg_0, dh_0, ff0_9, \
                         ff1_9, fg_14, fg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * fg_14[k];

        t_20[k] = f_1 * ff0_9[k]
                  - f_2 * ff1_9[k]
                  + pb_z[k] * fg_14[k];

        t_21[k] = pa_y[k] * dh_0[k];

        t_22[k] = f_7 * dg_0[k]
                  + pb_y[k] * fg_15[k];

        t_23[k] = pb_z[k] * fg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pb_z, dg_1, dg_3, dh_3, dh_5, \
                         dh_6, fg_16, fg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * dg_1[k]
                  + pa_y[k] * dh_3[k];

        t_25[k] = pb_z[k] * fg_16[k];

        t_26[k] = pa_y[k] * dh_5[k];

        t_27[k] = f_0 * dg_3[k]
                  + pa_y[k] * dh_6[k];

        t_28[k] = pb_z[k] * fg_18[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, dg_5, dg_25, dh_9, \
                         fg_20, fg_21, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_7 * dg_5[k]
                  + pb_y[k] * fg_20[k];

        t_30[k] = pa_y[k] * dh_9[k];

        t_31[k] = f_8 * dg_25[k]
                  + pb_x[k] * fg_25[k];

        t_32[k] = pb_z[k] * fg_21[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pb_x, pb_z, dg_10, dg_27, dg_28, \
                         dh_14, dh_15, fg_25, fg_27, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_8 * dg_27[k]
                  + pb_x[k] * fg_27[k];

        t_34[k] = f_8 * dg_28[k]
                  + pb_x[k] * fg_28[k];

        t_35[k] = pa_y[k] * dh_14[k];

        t_36[k] = f_9 * dg_10[k]
                  + pa_y[k] * dh_15[k];

        t_37[k] = pb_z[k] * fg_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, dg_12, dg_13, dg_14, \
                         dh_0, dh_17, dh_18, dh_20, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * dg_12[k]
                  + pa_y[k] * dh_17[k];

        t_39[k] = f_8 * dg_13[k]
                  + pa_y[k] * dh_18[k];

        t_40[k] = f_7 * dg_14[k]
                  + pb_y[k] * fg_29[k];

        t_41[k] = pa_y[k] * dh_20[k];

        t_42[k] = pa_z[k] * dh_0[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, t_48, pa_z, pb_y, pb_z, dg_0, dg_2, \
                         dh_3, dh_5, dh_6, fg_30, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * fg_30[k];

        t_44[k] = f_7 * dg_0[k]
                  + pb_z[k] * fg_30[k];

        t_45[k] = pa_z[k] * dh_3[k];

        t_46[k] = pb_y[k] * fg_32[k];

        t_47[k] = f_8 * dg_2[k]
                  + pa_z[k] * dh_5[k];

        t_48[k] = pa_z[k] * dh_6[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_z, pb_y, pb_z, dg_3, dg_5, dh_9, dh_10, \
                         fg_33, fg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_7 * dg_3[k]
                  + pb_z[k] * fg_33[k];

        t_50[k] = pb_y[k] * fg_35[k];

        t_51[k] = f_0 * dg_5[k]
                  + pa_z[k] * dh_9[k];

        t_52[k] = pa_z[k] * dh_10[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, dg_41, dg_42, dg_44, \
                         dh_15, fg_39, fg_41, fg_42, fg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_8 * dg_41[k]
                  + pb_x[k] * fg_41[k];

        t_54[k] = f_8 * dg_42[k]
                  + pb_x[k] * fg_42[k];

        t_55[k] = pb_y[k] * fg_39[k];

        t_56[k] = f_8 * dg_44[k]
                  + pb_x[k] * fg_44[k];

        t_57[k] = pa_z[k] * dh_15[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pb_y, pb_z, dg_10, dg_11, dg_12, dh_17, \
                         dh_18, fg_40, fg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_7 * dg_10[k]
                  + pb_z[k] * fg_40[k];

        t_59[k] = f_8 * dg_11[k]
                  + pa_z[k] * dh_17[k];

        t_60[k] = f_0 * dg_12[k]
                  + pa_z[k] * dh_18[k];

        t_61[k] = pb_y[k] * fg_44[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pa_z, pb_y, pb_z, dg_14, dg_15, dg_45, \
                         dh_20, dh_63, fg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_9 * dg_14[k]
                  + pa_z[k] * dh_20[k];

        t_63[k] = f_9 * dg_45[k]
                  + pa_x[k] * dh_63[k];

        t_64[k] = f_8 * dg_15[k]
                  + pb_y[k] * fg_45[k];

        t_65[k] = pb_z[k] * fg_45[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_x, pb_z, dg_48, dg_50, dg_51, dh_66, \
                         dh_68, dh_69, fg_46, fg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_0 * dg_48[k]
                  + pa_x[k] * dh_66[k];

        t_67[k] = pb_z[k] * fg_46[k];

        t_68[k] = f_0 * dg_50[k]
                  + pa_x[k] * dh_68[k];

        t_69[k] = f_8 * dg_51[k]
                  + pa_x[k] * dh_69[k];

        t_70[k] = pb_z[k] * fg_48[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_x, pb_x, pb_y, pb_z, dg_20, dg_54, dg_55, \
                         dh_72, fg_50, fg_51, fg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_8 * dg_20[k]
                  + pb_y[k] * fg_50[k];

        t_72[k] = f_8 * dg_54[k]
                  + pa_x[k] * dh_72[k];

        t_73[k] = f_7 * dg_55[k]
                  + pb_x[k] * fg_55[k];

        t_74[k] = pb_z[k] * fg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pa_x, pb_x, pb_z, dg_57, dg_58, dg_59, \
                         dh_78, fg_55, fg_57, fg_58, fg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_7 * dg_57[k]
                  + pb_x[k] * fg_57[k];

        t_76[k] = f_7 * dg_58[k]
                  + pb_x[k] * fg_58[k];

        t_77[k] = f_7 * dg_59[k]
                  + pb_x[k] * fg_59[k];

        t_78[k] = pa_x[k] * dh_78[k];

        t_79[k] = pb_z[k] * fg_55[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, pa_x, pa_y, pa_z, dh_22, dh_42, \
                         dh_80, dh_81, dh_82, dh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_x[k] * dh_80[k];

        t_81[k] = pa_x[k] * dh_81[k];

        t_82[k] = pa_x[k] * dh_82[k];

        t_83[k] = pa_x[k] * dh_83[k];

        t_84[k] = pa_y[k] * dh_42[k];

        t_85[k] = pa_z[k] * dh_22[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, pa_y, pa_z, pb_y, dg_32, dh_24, dh_27, \
                         dh_44, dh_47, fg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_y[k] * dh_44[k];

        t_87[k] = pa_z[k] * dh_24[k];

        t_88[k] = f_7 * dg_32[k]
                  + pb_y[k] * fg_62[k];

        t_89[k] = pa_y[k] * dh_47[k];

        t_90[k] = pa_z[k] * dh_27[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pa_z, pb_y, pb_z, dg_18, dg_35, dh_31, \
                         dh_51, fg_63, fg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_7 * dg_18[k]
                  + pb_z[k] * fg_63[k];

        t_92[k] = f_7 * dg_35[k]
                  + pb_y[k] * fg_65[k];

        t_93[k] = pa_y[k] * dh_51[k];

        t_94[k] = pa_z[k] * dh_31[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pa_x, pa_y, pb_x, dg_71, dg_72, dg_73, \
                         dh_56, dh_99, fg_71, fg_72, fg_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_7 * dg_71[k]
                  + pb_x[k] * fg_71[k];

        t_96[k] = f_7 * dg_72[k]
                  + pb_x[k] * fg_72[k];

        t_97[k] = f_7 * dg_73[k]
                  + pb_x[k] * fg_73[k];

        t_98[k] = pa_y[k] * dh_56[k];

        t_99[k] = pa_x[k] * dh_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, t_105, pa_x, dg_75, dh_100, \
                         dh_101, dh_102, dh_103, dh_104, dh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_x[k] * dh_100[k];

        t_101[k] = pa_x[k] * dh_101[k];

        t_102[k] = pa_x[k] * dh_102[k];

        t_103[k] = pa_x[k] * dh_103[k];

        t_104[k] = pa_x[k] * dh_104[k];

        t_105[k] = f_9 * dg_75[k]
                   + pa_x[k] * dh_105[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, pa_x, pb_y, pb_z, dg_30, dg_78, \
                         dg_80, dh_108, dh_110, fg_75, fg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = pb_y[k] * fg_75[k];

        t_107[k] = f_8 * dg_30[k]
                   + pb_z[k] * fg_75[k];

        t_108[k] = f_0 * dg_78[k]
                   + pa_x[k] * dh_108[k];

        t_109[k] = pb_y[k] * fg_77[k];

        t_110[k] = f_0 * dg_80[k]
                   + pa_x[k] * dh_110[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_x, pb_y, pb_z, dg_33, dg_81, dg_84, \
                         dh_111, dh_114, fg_78, fg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_8 * dg_81[k]
                   + pa_x[k] * dh_111[k];

        t_112[k] = f_8 * dg_33[k]
                   + pb_z[k] * fg_78[k];

        t_113[k] = pb_y[k] * fg_80[k];

        t_114[k] = f_8 * dg_84[k]
                   + pa_x[k] * dh_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pb_x, pb_y, dg_85, dg_86, dg_87, \
                         dg_89, fg_84, fg_85, fg_86, fg_87, fg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_7 * dg_85[k]
                   + pb_x[k] * fg_85[k];

        t_116[k] = f_7 * dg_86[k]
                   + pb_x[k] * fg_86[k];

        t_117[k] = f_7 * dg_87[k]
                   + pb_x[k] * fg_87[k];

        t_118[k] = pb_y[k] * fg_84[k];

        t_119[k] = f_7 * dg_89[k]
                   + pb_x[k] * fg_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, pa_x, pb_y, dh_120, dh_121, \
                         dh_122, dh_123, dh_125, fg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_x[k] * dh_120[k];

        t_121[k] = pa_x[k] * dh_121[k];

        t_122[k] = pa_x[k] * dh_122[k];

        t_123[k] = pa_x[k] * dh_123[k];

        t_124[k] = pb_y[k] * fg_89[k];

        t_125[k] = pa_x[k] * dh_125[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, pb_x, pb_y, pb_z, dg_45, ff0_60, \
                         ff0_63, ff1_60, ff1_63, fg_90, fg_91, fg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_1 * ff0_60[k]
                   - f_2 * ff1_60[k]
                   + pb_x[k] * fg_90[k];

        t_127[k] = f_0 * dg_45[k]
                   + pb_y[k] * fg_90[k];

        t_128[k] = pb_z[k] * fg_90[k];

        t_129[k] = f_5 * ff0_63[k]
                   - f_6 * ff1_63[k]
                   + pb_x[k] * fg_93[k];

        t_130[k] = pb_z[k] * fg_91[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pb_x, pb_y, pb_z, dg_50, ff0_65, ff0_66, \
                         ff1_65, ff1_66, fg_93, fg_95, fg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_5 * ff0_65[k]
                   - f_6 * ff1_65[k]
                   + pb_x[k] * fg_95[k];

        t_132[k] = f_3 * ff0_66[k]
                   - f_4 * ff1_66[k]
                   + pb_x[k] * fg_96[k];

        t_133[k] = pb_z[k] * fg_93[k];

        t_134[k] = f_0 * dg_50[k]
                   + pb_y[k] * fg_95[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, t_140, pb_x, ff0_69, ff1_69, \
                         fg_99, fg_100, fg_101, fg_102, fg_103, \
                         fg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_3 * ff0_69[k]
                   - f_4 * ff1_69[k]
                   + pb_x[k] * fg_99[k];

        t_136[k] = pb_x[k] * fg_100[k];

        t_137[k] = pb_x[k] * fg_101[k];

        t_138[k] = pb_x[k] * fg_102[k];

        t_139[k] = pb_x[k] * fg_103[k];

        t_140[k] = pb_x[k] * fg_104[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pb_y, pb_z, dg_55, ff0_66, ff0_67, \
                         ff1_66, ff1_67, fg_100, fg_101, fg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_0 * dg_55[k]
                   + f_1 * ff0_66[k]
                   - f_2 * ff1_66[k]
                   + pb_y[k] * fg_100[k];

        t_142[k] = pb_z[k] * fg_100[k];

        t_143[k] = f_3 * ff0_66[k]
                   - f_4 * ff1_66[k]
                   + pb_z[k] * fg_101[k];

        t_144[k] = f_5 * ff0_67[k]
                   - f_6 * ff1_67[k]
                   + pb_z[k] * fg_102[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pa_z, pb_y, pb_z, dg_45, dg_59, \
                         dh_63, dh_64, ff0_69, ff1_69, fg_104, fg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_0 * dg_59[k]
                   + pb_y[k] * fg_104[k];

        t_146[k] = f_1 * ff0_69[k]
                   - f_2 * ff1_69[k]
                   + pb_z[k] * fg_104[k];

        t_147[k] = pa_z[k] * dh_63[k];

        t_148[k] = pa_z[k] * dh_64[k];

        t_149[k] = f_7 * dg_45[k]
                   + pb_z[k] * fg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pa_z, pb_y, pb_z, dg_47, dg_48, \
                         dg_62, dh_66, dh_68, dh_69, fg_107, fg_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * dh_66[k];

        t_151[k] = f_8 * dg_62[k]
                   + pb_y[k] * fg_107[k];

        t_152[k] = f_8 * dg_47[k]
                   + pa_z[k] * dh_68[k];

        t_153[k] = pa_z[k] * dh_69[k];

        t_154[k] = f_7 * dg_48[k]
                   + pb_z[k] * fg_108[k];
    }
}

static auto
compute_prim_fh_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t dg,
                                            const size_t dh, const size_t ff0, const size_t ff1,
                                            const size_t fg, const size_t ncols,
                                            const double alpha, const double beta,
                                            const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 2.5 / p;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_56 = buffer.data(dg + 56);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_63 = buffer.data(dg + 63);
    const auto *dg_65 = buffer.data(dg + 65);
    const auto *dg_70 = buffer.data(dg + 70);
    const auto *dg_74 = buffer.data(dg + 74);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_76 = buffer.data(dg + 76);
    const auto *dg_77 = buffer.data(dg + 77);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_80 = buffer.data(dg + 80);
    const auto *dg_85 = buffer.data(dg + 85);
    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_88 = buffer.data(dg + 88);
    const auto *dg_89 = buffer.data(dg + 89);

    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_83 = buffer.data(dh + 83);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_107 = buffer.data(dh + 107);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_111 = buffer.data(dh + 111);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_125 = buffer.data(dh + 125);

    const auto *ff0_90 = buffer.data(ff0 + 90);
    const auto *ff0_93 = buffer.data(ff0 + 93);
    const auto *ff0_95 = buffer.data(ff0 + 95);
    const auto *ff0_96 = buffer.data(ff0 + 96);
    const auto *ff0_98 = buffer.data(ff0 + 98);
    const auto *ff0_99 = buffer.data(ff0 + 99);

    const auto *ff1_90 = buffer.data(ff1 + 90);
    const auto *ff1_93 = buffer.data(ff1 + 93);
    const auto *ff1_95 = buffer.data(ff1 + 95);
    const auto *ff1_96 = buffer.data(ff1 + 96);
    const auto *ff1_98 = buffer.data(ff1 + 98);
    const auto *ff1_99 = buffer.data(ff1 + 99);

    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_120 = buffer.data(fg + 120);
    const auto *fg_122 = buffer.data(fg + 122);
    const auto *fg_123 = buffer.data(fg + 123);
    const auto *fg_125 = buffer.data(fg + 125);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_134 = buffer.data(fg + 134);
    const auto *fg_135 = buffer.data(fg + 135);
    const auto *fg_137 = buffer.data(fg + 137);
    const auto *fg_138 = buffer.data(fg + 138);
    const auto *fg_140 = buffer.data(fg + 140);
    const auto *fg_141 = buffer.data(fg + 141);
    const auto *fg_144 = buffer.data(fg + 144);
    const auto *fg_145 = buffer.data(fg + 145);
    const auto *fg_146 = buffer.data(fg + 146);
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_148 = buffer.data(fg + 148);
    const auto *fg_149 = buffer.data(fg + 149);

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pa_z, pb_x, pb_y, dg_50, dg_65, \
                         dh_72, fg_110, fg_115, fg_116, fg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_8 * dg_65[k]
                   + pb_y[k] * fg_110[k];

        t_156[k] = f_0 * dg_50[k]
                   + pa_z[k] * dh_72[k];

        t_157[k] = pb_x[k] * fg_115[k];

        t_158[k] = pb_x[k] * fg_116[k];

        t_159[k] = pb_x[k] * fg_117[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, pa_z, pb_x, pb_z, dg_55, dg_56, \
                         dh_78, dh_80, fg_115, fg_118, fg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pb_x[k] * fg_118[k];

        t_161[k] = pb_x[k] * fg_119[k];

        t_162[k] = pa_z[k] * dh_78[k];

        t_163[k] = f_7 * dg_55[k]
                   + pb_z[k] * fg_115[k];

        t_164[k] = f_8 * dg_56[k]
                   + pa_z[k] * dh_80[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pa_y, pa_z, pb_y, dg_57, dg_59, dg_74, \
                         dh_81, dh_83, dh_105, fg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_0 * dg_57[k]
                   + pa_z[k] * dh_81[k];

        t_166[k] = f_8 * dg_74[k]
                   + pb_y[k] * fg_119[k];

        t_167[k] = f_9 * dg_59[k]
                   + pa_z[k] * dh_83[k];

        t_168[k] = pa_y[k] * dh_105[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, pa_y, pb_y, dg_75, dg_76, dg_77, \
                         dh_107, dh_108, dh_110, fg_120, fg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_7 * dg_75[k]
                   + pb_y[k] * fg_120[k];

        t_170[k] = pa_y[k] * dh_107[k];

        t_171[k] = f_8 * dg_76[k]
                   + pa_y[k] * dh_108[k];

        t_172[k] = f_7 * dg_77[k]
                   + pb_y[k] * fg_122[k];

        t_173[k] = pa_y[k] * dh_110[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pa_y, pb_y, pb_z, dg_63, dg_78, dg_80, \
                         dh_111, dh_114, fg_123, fg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_0 * dg_78[k]
                   + pa_y[k] * dh_111[k];

        t_175[k] = f_8 * dg_63[k]
                   + pb_z[k] * fg_123[k];

        t_176[k] = f_7 * dg_80[k]
                   + pb_y[k] * fg_125[k];

        t_177[k] = pa_y[k] * dh_114[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, t_183, pa_y, pb_x, dg_85, dh_120, \
                         fg_130, fg_131, fg_132, fg_133, fg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = pb_x[k] * fg_130[k];

        t_179[k] = pb_x[k] * fg_131[k];

        t_180[k] = pb_x[k] * fg_132[k];

        t_181[k] = pb_x[k] * fg_133[k];

        t_182[k] = pb_x[k] * fg_134[k];

        t_183[k] = f_9 * dg_85[k]
                   + pa_y[k] * dh_120[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pb_y, pb_z, dg_70, dg_87, dg_88, \
                         dg_89, dh_122, dh_123, fg_130, fg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_8 * dg_70[k]
                   + pb_z[k] * fg_130[k];

        t_185[k] = f_0 * dg_87[k]
                   + pa_y[k] * dh_122[k];

        t_186[k] = f_8 * dg_88[k]
                   + pa_y[k] * dh_123[k];

        t_187[k] = f_7 * dg_89[k]
                   + pb_y[k] * fg_134[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pb_x, pb_y, pb_z, dg_75, dh_125, \
                         ff0_90, ff1_90, fg_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pa_y[k] * dh_125[k];

        t_189[k] = f_1 * ff0_90[k]
                   - f_2 * ff1_90[k]
                   + pb_x[k] * fg_135[k];

        t_190[k] = pb_y[k] * fg_135[k];

        t_191[k] = f_0 * dg_75[k]
                   + pb_z[k] * fg_135[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pb_x, pb_y, ff0_93, ff0_95, ff0_96, \
                         ff1_93, ff1_95, ff1_96, fg_137, fg_138, fg_140, \
                         fg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_5 * ff0_93[k]
                   - f_6 * ff1_93[k]
                   + pb_x[k] * fg_138[k];

        t_193[k] = pb_y[k] * fg_137[k];

        t_194[k] = f_5 * ff0_95[k]
                   - f_6 * ff1_95[k]
                   + pb_x[k] * fg_140[k];

        t_195[k] = f_3 * ff0_96[k]
                   - f_4 * ff1_96[k]
                   + pb_x[k] * fg_141[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pb_x, pb_y, pb_z, dg_78, ff0_99, \
                         ff1_99, fg_138, fg_140, fg_144, fg_145, \
                         fg_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_0 * dg_78[k]
                   + pb_z[k] * fg_138[k];

        t_197[k] = pb_y[k] * fg_140[k];

        t_198[k] = f_3 * ff0_99[k]
                   - f_4 * ff1_99[k]
                   + pb_x[k] * fg_144[k];

        t_199[k] = pb_x[k] * fg_145[k];

        t_200[k] = pb_x[k] * fg_146[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, pb_x, pb_y, pb_z, dg_85, ff0_96, \
                         ff1_96, fg_145, fg_147, fg_148, fg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pb_x[k] * fg_147[k];

        t_202[k] = pb_x[k] * fg_148[k];

        t_203[k] = pb_x[k] * fg_149[k];

        t_204[k] = f_1 * ff0_96[k]
                   - f_2 * ff1_96[k]
                   + pb_y[k] * fg_145[k];

        t_205[k] = f_0 * dg_85[k]
                   + pb_z[k] * fg_145[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pb_y, pb_z, dg_89, ff0_98, ff0_99, \
                         ff1_98, ff1_99, fg_147, fg_148, fg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_5 * ff0_98[k]
                   - f_6 * ff1_98[k]
                   + pb_y[k] * fg_147[k];

        t_207[k] = f_3 * ff0_99[k]
                   - f_4 * ff1_99[k]
                   + pb_y[k] * fg_148[k];

        t_208[k] = pb_y[k] * fg_149[k];

        t_209[k] = f_0 * dg_89[k]
                   + f_1 * ff0_99[k]
                   - f_2 * ff1_99[k]
                   + pb_z[k] * fg_149[k];
    }
}

auto
compute_prim_fh_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg, const size_t dh,
                                     const size_t ff0, const size_t ff1, const size_t fg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    compute_prim_fh_electron_repulsion_0_piece0(buffer, target, pa, pb, dg, dh, ff0, ff1, fg,
                                                ncols, alpha, beta, p);

    compute_prim_fh_electron_repulsion_0_piece1(buffer, target, pa, pb, dg, dh, ff0, ff1, fg,
                                                ncols, alpha, beta, p);
}

}  // namespace simdt2ceri
