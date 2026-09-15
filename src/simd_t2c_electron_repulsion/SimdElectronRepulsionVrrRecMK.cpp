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


#include "SimdElectronRepulsionVrrRecMK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_mk_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kk0,
                                            const size_t kk1, const size_t li, const size_t lk,
                                            const size_t mh0, const size_t mh1, const size_t mi,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_16 = 4.0 / p;
    const auto f_17 = 3.5 / p;
    const auto f_18 = 0.5 / alpha;
    const auto f_19 = 0.5 * beta / (alpha * p);
    const auto f_20 = 3.0 / alpha;
    const auto f_21 = 3.0 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kk0_0 = buffer.data(kk0 + 0);
    const auto *kk0_136 = buffer.data(kk0 + 136);

    const auto *kk1_0 = buffer.data(kk1 + 0);
    const auto *kk1_136 = buffer.data(kk1 + 136);

    const auto *li_0 = buffer.data(li + 0);
    const auto *li_1 = buffer.data(li + 1);
    const auto *li_2 = buffer.data(li + 2);
    const auto *li_3 = buffer.data(li + 3);
    const auto *li_5 = buffer.data(li + 5);
    const auto *li_6 = buffer.data(li + 6);
    const auto *li_7 = buffer.data(li + 7);
    const auto *li_8 = buffer.data(li + 8);
    const auto *li_9 = buffer.data(li + 9);
    const auto *li_10 = buffer.data(li + 10);
    const auto *li_11 = buffer.data(li + 11);
    const auto *li_12 = buffer.data(li + 12);
    const auto *li_13 = buffer.data(li + 13);
    const auto *li_14 = buffer.data(li + 14);
    const auto *li_21 = buffer.data(li + 21);
    const auto *li_22 = buffer.data(li + 22);
    const auto *li_23 = buffer.data(li + 23);
    const auto *li_24 = buffer.data(li + 24);
    const auto *li_25 = buffer.data(li + 25);
    const auto *li_26 = buffer.data(li + 26);
    const auto *li_27 = buffer.data(li + 27);
    const auto *li_28 = buffer.data(li + 28);
    const auto *li_33 = buffer.data(li + 33);
    const auto *li_37 = buffer.data(li + 37);
    const auto *li_42 = buffer.data(li + 42);
    const auto *li_49 = buffer.data(li + 49);
    const auto *li_51 = buffer.data(li + 51);
    const auto *li_52 = buffer.data(li + 52);
    const auto *li_53 = buffer.data(li + 53);
    const auto *li_54 = buffer.data(li + 54);
    const auto *li_78 = buffer.data(li + 78);
    const auto *li_79 = buffer.data(li + 79);
    const auto *li_80 = buffer.data(li + 80);
    const auto *li_81 = buffer.data(li + 81);
    const auto *li_83 = buffer.data(li + 83);
    const auto *li_87 = buffer.data(li + 87);
    const auto *li_90 = buffer.data(li + 90);
    const auto *li_94 = buffer.data(li + 94);
    const auto *li_99 = buffer.data(li + 99);
    const auto *li_105 = buffer.data(li + 105);
    const auto *li_107 = buffer.data(li + 107);
    const auto *li_108 = buffer.data(li + 108);
    const auto *li_109 = buffer.data(li + 109);
    const auto *li_110 = buffer.data(li + 110);
    const auto *li_111 = buffer.data(li + 111);

    const auto *lk_0 = buffer.data(lk + 0);
    const auto *lk_3 = buffer.data(lk + 3);
    const auto *lk_5 = buffer.data(lk + 5);
    const auto *lk_6 = buffer.data(lk + 6);
    const auto *lk_9 = buffer.data(lk + 9);
    const auto *lk_10 = buffer.data(lk + 10);
    const auto *lk_12 = buffer.data(lk + 12);
    const auto *lk_14 = buffer.data(lk + 14);
    const auto *lk_15 = buffer.data(lk + 15);
    const auto *lk_17 = buffer.data(lk + 17);
    const auto *lk_18 = buffer.data(lk + 18);
    const auto *lk_20 = buffer.data(lk + 20);
    const auto *lk_21 = buffer.data(lk + 21);
    const auto *lk_27 = buffer.data(lk + 27);
    const auto *lk_28 = buffer.data(lk + 28);
    const auto *lk_30 = buffer.data(lk + 30);
    const auto *lk_31 = buffer.data(lk + 31);
    const auto *lk_32 = buffer.data(lk + 32);
    const auto *lk_33 = buffer.data(lk + 33);
    const auto *lk_35 = buffer.data(lk + 35);
    const auto *lk_36 = buffer.data(lk + 36);
    const auto *lk_136 = buffer.data(lk + 136);

    const auto *mh0_0 = buffer.data(mh0 + 0);
    const auto *mh0_1 = buffer.data(mh0 + 1);
    const auto *mh0_2 = buffer.data(mh0 + 2);
    const auto *mh0_3 = buffer.data(mh0 + 3);
    const auto *mh0_5 = buffer.data(mh0 + 5);
    const auto *mh0_6 = buffer.data(mh0 + 6);
    const auto *mh0_8 = buffer.data(mh0 + 8);
    const auto *mh0_9 = buffer.data(mh0 + 9);
    const auto *mh0_15 = buffer.data(mh0 + 15);
    const auto *mh0_17 = buffer.data(mh0 + 17);
    const auto *mh0_18 = buffer.data(mh0 + 18);
    const auto *mh0_19 = buffer.data(mh0 + 19);
    const auto *mh0_20 = buffer.data(mh0 + 20);
    const auto *mh0_63 = buffer.data(mh0 + 63);
    const auto *mh0_65 = buffer.data(mh0 + 65);
    const auto *mh0_66 = buffer.data(mh0 + 66);
    const auto *mh0_68 = buffer.data(mh0 + 68);
    const auto *mh0_69 = buffer.data(mh0 + 69);
    const auto *mh0_70 = buffer.data(mh0 + 70);
    const auto *mh0_72 = buffer.data(mh0 + 72);
    const auto *mh0_73 = buffer.data(mh0 + 73);
    const auto *mh0_78 = buffer.data(mh0 + 78);
    const auto *mh0_79 = buffer.data(mh0 + 79);
    const auto *mh0_80 = buffer.data(mh0 + 80);

    const auto *mh1_0 = buffer.data(mh1 + 0);
    const auto *mh1_1 = buffer.data(mh1 + 1);
    const auto *mh1_2 = buffer.data(mh1 + 2);
    const auto *mh1_3 = buffer.data(mh1 + 3);
    const auto *mh1_5 = buffer.data(mh1 + 5);
    const auto *mh1_6 = buffer.data(mh1 + 6);
    const auto *mh1_8 = buffer.data(mh1 + 8);
    const auto *mh1_9 = buffer.data(mh1 + 9);
    const auto *mh1_15 = buffer.data(mh1 + 15);
    const auto *mh1_17 = buffer.data(mh1 + 17);
    const auto *mh1_18 = buffer.data(mh1 + 18);
    const auto *mh1_19 = buffer.data(mh1 + 19);
    const auto *mh1_20 = buffer.data(mh1 + 20);
    const auto *mh1_63 = buffer.data(mh1 + 63);
    const auto *mh1_65 = buffer.data(mh1 + 65);
    const auto *mh1_66 = buffer.data(mh1 + 66);
    const auto *mh1_68 = buffer.data(mh1 + 68);
    const auto *mh1_69 = buffer.data(mh1 + 69);
    const auto *mh1_70 = buffer.data(mh1 + 70);
    const auto *mh1_72 = buffer.data(mh1 + 72);
    const auto *mh1_73 = buffer.data(mh1 + 73);
    const auto *mh1_78 = buffer.data(mh1 + 78);
    const auto *mh1_79 = buffer.data(mh1 + 79);
    const auto *mh1_80 = buffer.data(mh1 + 80);

    const auto *mi_0 = buffer.data(mi + 0);
    const auto *mi_1 = buffer.data(mi + 1);
    const auto *mi_2 = buffer.data(mi + 2);
    const auto *mi_3 = buffer.data(mi + 3);
    const auto *mi_5 = buffer.data(mi + 5);
    const auto *mi_6 = buffer.data(mi + 6);
    const auto *mi_8 = buffer.data(mi + 8);
    const auto *mi_9 = buffer.data(mi + 9);
    const auto *mi_10 = buffer.data(mi + 10);
    const auto *mi_12 = buffer.data(mi + 12);
    const auto *mi_13 = buffer.data(mi + 13);
    const auto *mi_14 = buffer.data(mi + 14);
    const auto *mi_15 = buffer.data(mi + 15);
    const auto *mi_20 = buffer.data(mi + 20);
    const auto *mi_21 = buffer.data(mi + 21);
    const auto *mi_23 = buffer.data(mi + 23);
    const auto *mi_24 = buffer.data(mi + 24);
    const auto *mi_25 = buffer.data(mi + 25);
    const auto *mi_26 = buffer.data(mi + 26);
    const auto *mi_27 = buffer.data(mi + 27);
    const auto *mi_28 = buffer.data(mi + 28);
    const auto *mi_29 = buffer.data(mi + 29);
    const auto *mi_31 = buffer.data(mi + 31);
    const auto *mi_33 = buffer.data(mi + 33);
    const auto *mi_34 = buffer.data(mi + 34);
    const auto *mi_37 = buffer.data(mi + 37);
    const auto *mi_38 = buffer.data(mi + 38);
    const auto *mi_42 = buffer.data(mi + 42);
    const auto *mi_43 = buffer.data(mi + 43);
    const auto *mi_49 = buffer.data(mi + 49);
    const auto *mi_51 = buffer.data(mi + 51);
    const auto *mi_52 = buffer.data(mi + 52);
    const auto *mi_53 = buffer.data(mi + 53);
    const auto *mi_54 = buffer.data(mi + 54);
    const auto *mi_55 = buffer.data(mi + 55);
    const auto *mi_56 = buffer.data(mi + 56);
    const auto *mi_58 = buffer.data(mi + 58);
    const auto *mi_59 = buffer.data(mi + 59);
    const auto *mi_61 = buffer.data(mi + 61);
    const auto *mi_62 = buffer.data(mi + 62);
    const auto *mi_65 = buffer.data(mi + 65);
    const auto *mi_66 = buffer.data(mi + 66);
    const auto *mi_70 = buffer.data(mi + 70);
    const auto *mi_76 = buffer.data(mi + 76);
    const auto *mi_77 = buffer.data(mi + 77);
    const auto *mi_78 = buffer.data(mi + 78);
    const auto *mi_79 = buffer.data(mi + 79);
    const auto *mi_80 = buffer.data(mi + 80);
    const auto *mi_81 = buffer.data(mi + 81);
    const auto *mi_83 = buffer.data(mi + 83);
    const auto *mi_84 = buffer.data(mi + 84);
    const auto *mi_85 = buffer.data(mi + 85);
    const auto *mi_86 = buffer.data(mi + 86);
    const auto *mi_87 = buffer.data(mi + 87);
    const auto *mi_89 = buffer.data(mi + 89);
    const auto *mi_90 = buffer.data(mi + 90);
    const auto *mi_91 = buffer.data(mi + 91);
    const auto *mi_93 = buffer.data(mi + 93);
    const auto *mi_94 = buffer.data(mi + 94);
    const auto *mi_95 = buffer.data(mi + 95);
    const auto *mi_96 = buffer.data(mi + 96);
    const auto *mi_98 = buffer.data(mi + 98);
    const auto *mi_99 = buffer.data(mi + 99);
    const auto *mi_105 = buffer.data(mi + 105);
    const auto *mi_106 = buffer.data(mi + 106);
    const auto *mi_107 = buffer.data(mi + 107);
    const auto *mi_108 = buffer.data(mi + 108);
    const auto *mi_109 = buffer.data(mi + 109);
    const auto *mi_110 = buffer.data(mi + 110);
    const auto *mi_111 = buffer.data(mi + 111);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, li_0, mh0_0, mh1_0, \
                         mi_0, mi_1, mi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * li_0[k]
                 + f_1 * mh0_0[k]
                 - f_2 * mh1_0[k]
                 + pb_x[k] * mi_0[k];

        t_1[k] = pb_y[k] * mi_0[k];

        t_2[k] = pb_z[k] * mi_0[k];

        t_3[k] = f_3 * mh0_0[k]
                 - f_4 * mh1_0[k]
                 + pb_y[k] * mi_1[k];

        t_4[k] = pb_y[k] * mi_2[k];

        t_5[k] = f_3 * mh0_0[k]
                 - f_4 * mh1_0[k]
                 + pb_z[k] * mi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, mh0_1, mh0_2, mh0_3, mh1_1, \
                         mh1_2, mh1_3, mi_3, mi_5, mi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * mh0_1[k]
                 - f_6 * mh1_1[k]
                 + pb_y[k] * mi_3[k];

        t_7[k] = pb_z[k] * mi_3[k];

        t_8[k] = pb_y[k] * mi_5[k];

        t_9[k] = f_5 * mh0_2[k]
                 - f_6 * mh1_2[k]
                 + pb_z[k] * mi_5[k];

        t_10[k] = f_7 * mh0_3[k]
                  - f_8 * mh1_3[k]
                  + pb_y[k] * mi_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, mh0_5, mh0_6, mh1_5, \
                         mh1_6, mi_6, mi_8, mi_9, mi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * mi_6[k];

        t_12[k] = f_3 * mh0_5[k]
                  - f_4 * mh1_5[k]
                  + pb_y[k] * mi_8[k];

        t_13[k] = pb_y[k] * mi_9[k];

        t_14[k] = f_7 * mh0_5[k]
                  - f_8 * mh1_5[k]
                  + pb_z[k] * mi_9[k];

        t_15[k] = f_9 * mh0_6[k]
                  - f_10 * mh1_6[k]
                  + pb_y[k] * mi_10[k];

        t_16[k] = pb_z[k] * mi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, mh0_8, mh0_9, mh1_8, mh1_9, \
                         mi_12, mi_13, mi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * mh0_8[k]
                  - f_6 * mh1_8[k]
                  + pb_y[k] * mi_12[k];

        t_18[k] = f_3 * mh0_9[k]
                  - f_4 * mh1_9[k]
                  + pb_y[k] * mi_13[k];

        t_19[k] = pb_y[k] * mi_14[k];

        t_20[k] = f_9 * mh0_9[k]
                  - f_10 * mh1_9[k]
                  + pb_z[k] * mi_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_z, li_21, li_23, li_24, li_25, \
                         mi_15, mi_21, mi_23, mi_24, mi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * li_21[k]
                  + pb_x[k] * mi_21[k];

        t_22[k] = pb_z[k] * mi_15[k];

        t_23[k] = f_0 * li_23[k]
                  + pb_x[k] * mi_23[k];

        t_24[k] = f_0 * li_24[k]
                  + pb_x[k] * mi_24[k];

        t_25[k] = f_0 * li_25[k]
                  + pb_x[k] * mi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, li_27, mh0_15, mh1_15, \
                         mi_20, mi_21, mi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * mi_20[k];

        t_27[k] = f_0 * li_27[k]
                  + pb_x[k] * mi_27[k];

        t_28[k] = f_1 * mh0_15[k]
                  - f_2 * mh1_15[k]
                  + pb_y[k] * mi_21[k];

        t_29[k] = pb_z[k] * mi_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, mh0_17, mh0_18, mh0_19, mh1_17, mh1_18, \
                         mh1_19, mi_23, mi_24, mi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_9 * mh0_17[k]
                  - f_10 * mh1_17[k]
                  + pb_y[k] * mi_23[k];

        t_31[k] = f_7 * mh0_18[k]
                  - f_8 * mh1_18[k]
                  + pb_y[k] * mi_24[k];

        t_32[k] = f_5 * mh0_19[k]
                  - f_6 * mh1_19[k]
                  + pb_y[k] * mi_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, li_0, lk_0, \
                         mh0_20, mh1_20, mi_26, mi_27, mi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * mh0_20[k]
                  - f_4 * mh1_20[k]
                  + pb_y[k] * mi_26[k];

        t_34[k] = pb_y[k] * mi_27[k];

        t_35[k] = f_1 * mh0_20[k]
                  - f_2 * mh1_20[k]
                  + pb_z[k] * mi_27[k];

        t_36[k] = pa_y[k] * lk_0[k];

        t_37[k] = f_11 * li_0[k]
                  + pb_y[k] * mi_28[k];

        t_38[k] = pb_z[k] * mi_28[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pb_z, li_1, li_3, lk_3, lk_5, \
                         lk_6, mi_29, mi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * li_1[k]
                  + pa_y[k] * lk_3[k];

        t_40[k] = pb_z[k] * mi_29[k];

        t_41[k] = pa_y[k] * lk_5[k];

        t_42[k] = f_13 * li_3[k]
                  + pa_y[k] * lk_6[k];

        t_43[k] = pb_z[k] * mi_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_y, pb_y, pb_z, li_5, li_6, li_8, \
                         lk_9, lk_10, lk_12, mi_33, mi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_11 * li_5[k]
                  + pb_y[k] * mi_33[k];

        t_45[k] = pa_y[k] * lk_9[k];

        t_46[k] = f_14 * li_6[k]
                  + pa_y[k] * lk_10[k];

        t_47[k] = pb_z[k] * mi_34[k];

        t_48[k] = f_12 * li_8[k]
                  + pa_y[k] * lk_12[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_y, pb_y, pb_z, li_9, li_10, li_12, \
                         lk_14, lk_15, lk_17, mi_37, mi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_11 * li_9[k]
                  + pb_y[k] * mi_37[k];

        t_50[k] = pa_y[k] * lk_14[k];

        t_51[k] = f_15 * li_10[k]
                  + pa_y[k] * lk_15[k];

        t_52[k] = pb_z[k] * mi_38[k];

        t_53[k] = f_13 * li_12[k]
                  + pa_y[k] * lk_17[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pb_x, pb_y, li_13, li_14, li_49, lk_18, \
                         lk_20, mi_42, mi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * li_13[k]
                  + pa_y[k] * lk_18[k];

        t_55[k] = f_11 * li_14[k]
                  + pb_y[k] * mi_42[k];

        t_56[k] = pa_y[k] * lk_20[k];

        t_57[k] = f_16 * li_49[k]
                  + pb_x[k] * mi_49[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pb_z, li_51, li_52, li_53, li_54, \
                         mi_43, mi_51, mi_52, mi_53, mi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_z[k] * mi_43[k];

        t_59[k] = f_16 * li_51[k]
                  + pb_x[k] * mi_51[k];

        t_60[k] = f_16 * li_52[k]
                  + pb_x[k] * mi_52[k];

        t_61[k] = f_16 * li_53[k]
                  + pb_x[k] * mi_53[k];

        t_62[k] = f_16 * li_54[k]
                  + pb_x[k] * mi_54[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_z, li_21, li_23, li_24, lk_27, \
                         lk_28, lk_30, lk_31, mi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_y[k] * lk_27[k];

        t_64[k] = f_17 * li_21[k]
                  + pa_y[k] * lk_28[k];

        t_65[k] = pb_z[k] * mi_49[k];

        t_66[k] = f_15 * li_23[k]
                  + pa_y[k] * lk_30[k];

        t_67[k] = f_14 * li_24[k]
                  + pa_y[k] * lk_31[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pa_z, pb_y, li_25, li_26, li_27, \
                         lk_0, lk_32, lk_33, lk_35, mi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_13 * li_25[k]
                  + pa_y[k] * lk_32[k];

        t_69[k] = f_12 * li_26[k]
                  + pa_y[k] * lk_33[k];

        t_70[k] = f_11 * li_27[k]
                  + pb_y[k] * mi_55[k];

        t_71[k] = pa_y[k] * lk_35[k];

        t_72[k] = pa_z[k] * lk_0[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, li_0, li_2, \
                         lk_3, lk_5, lk_6, mi_56, mi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_y[k] * mi_56[k];

        t_74[k] = f_11 * li_0[k]
                  + pb_z[k] * mi_56[k];

        t_75[k] = pa_z[k] * lk_3[k];

        t_76[k] = pb_y[k] * mi_58[k];

        t_77[k] = f_12 * li_2[k]
                  + pa_z[k] * lk_5[k];

        t_78[k] = pa_z[k] * lk_6[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_z, pb_y, pb_z, li_3, li_5, li_6, \
                         lk_9, lk_10, mi_59, mi_61, mi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_11 * li_3[k]
                  + pb_z[k] * mi_59[k];

        t_80[k] = pb_y[k] * mi_61[k];

        t_81[k] = f_13 * li_5[k]
                  + pa_z[k] * lk_9[k];

        t_82[k] = pa_z[k] * lk_10[k];

        t_83[k] = f_11 * li_6[k]
                  + pb_z[k] * mi_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pa_z, pb_y, pb_z, li_7, li_9, li_10, \
                         lk_12, lk_14, lk_15, mi_65, mi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * li_7[k]
                  + pa_z[k] * lk_12[k];

        t_85[k] = pb_y[k] * mi_65[k];

        t_86[k] = f_14 * li_9[k]
                  + pa_z[k] * lk_14[k];

        t_87[k] = pa_z[k] * lk_15[k];

        t_88[k] = f_11 * li_10[k]
                  + pb_z[k] * mi_66[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, pa_z, pb_y, li_11, li_12, li_14, lk_17, \
                         lk_18, lk_20, lk_21, mi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_12 * li_11[k]
                  + pa_z[k] * lk_17[k];

        t_90[k] = f_13 * li_12[k]
                  + pa_z[k] * lk_18[k];

        t_91[k] = pb_y[k] * mi_70[k];

        t_92[k] = f_15 * li_14[k]
                  + pa_z[k] * lk_20[k];

        t_93[k] = pa_z[k] * lk_21[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pb_x, pb_y, li_78, li_79, li_80, li_81, \
                         mi_76, mi_78, mi_79, mi_80, mi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_16 * li_78[k]
                  + pb_x[k] * mi_78[k];

        t_95[k] = f_16 * li_79[k]
                  + pb_x[k] * mi_79[k];

        t_96[k] = f_16 * li_80[k]
                  + pb_x[k] * mi_80[k];

        t_97[k] = f_16 * li_81[k]
                  + pb_x[k] * mi_81[k];

        t_98[k] = pb_y[k] * mi_76[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_z, pb_x, pb_z, li_21, li_22, li_83, \
                         lk_28, lk_30, mi_77, mi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_16 * li_83[k]
                  + pb_x[k] * mi_83[k];

        t_100[k] = pa_z[k] * lk_28[k];

        t_101[k] = f_11 * li_21[k]
                   + pb_z[k] * mi_77[k];

        t_102[k] = f_12 * li_22[k]
                   + pa_z[k] * lk_30[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, li_23, li_24, li_25, \
                         li_27, lk_31, lk_32, lk_33, lk_35, mi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_13 * li_23[k]
                   + pa_z[k] * lk_31[k];

        t_104[k] = f_14 * li_24[k]
                   + pa_z[k] * lk_32[k];

        t_105[k] = f_15 * li_25[k]
                   + pa_z[k] * lk_33[k];

        t_106[k] = pb_y[k] * mi_83[k];

        t_107[k] = f_17 * li_27[k]
                   + pa_z[k] * lk_35[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_y, pb_y, pb_z, kk0_0, kk1_0, li_28, lk_36, \
                         mi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_18 * kk0_0[k]
                   - f_19 * kk1_0[k]
                   + pa_y[k] * lk_36[k];

        t_109[k] = f_12 * li_28[k]
                   + pb_y[k] * mi_84[k];

        t_110[k] = pb_z[k] * mi_84[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pb_x, pb_z, li_87, mh0_63, mh0_66, mh1_63, \
                         mh1_66, mi_85, mi_86, mi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_17 * li_87[k]
                   + f_9 * mh0_66[k]
                   - f_10 * mh1_66[k]
                   + pb_x[k] * mi_87[k];

        t_112[k] = pb_z[k] * mi_85[k];

        t_113[k] = f_3 * mh0_63[k]
                   - f_4 * mh1_63[k]
                   + pb_z[k] * mi_86[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, pb_y, pb_z, li_33, li_90, mh0_65, \
                         mh0_69, mh1_65, mh1_69, mi_87, mi_89, mi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_17 * li_90[k]
                   + f_7 * mh0_69[k]
                   - f_8 * mh1_69[k]
                   + pb_x[k] * mi_90[k];

        t_115[k] = pb_z[k] * mi_87[k];

        t_116[k] = f_12 * li_33[k]
                   + pb_y[k] * mi_89[k];

        t_117[k] = f_5 * mh0_65[k]
                   - f_6 * mh1_65[k]
                   + pb_z[k] * mi_89[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pb_x, pb_z, li_94, mh0_66, mh0_73, mh1_66, \
                         mh1_73, mi_90, mi_91, mi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_17 * li_94[k]
                   + f_5 * mh0_73[k]
                   - f_6 * mh1_73[k]
                   + pb_x[k] * mi_94[k];

        t_119[k] = pb_z[k] * mi_90[k];

        t_120[k] = f_3 * mh0_66[k]
                   - f_4 * mh1_66[k]
                   + pb_z[k] * mi_91[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_x, pb_y, pb_z, li_37, li_99, mh0_68, \
                         mh0_78, mh1_68, mh1_78, mi_93, mi_94, mi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_12 * li_37[k]
                   + pb_y[k] * mi_93[k];

        t_122[k] = f_7 * mh0_68[k]
                   - f_8 * mh1_68[k]
                   + pb_z[k] * mi_93[k];

        t_123[k] = f_17 * li_99[k]
                   + f_3 * mh0_78[k]
                   - f_4 * mh1_78[k]
                   + pb_x[k] * mi_99[k];

        t_124[k] = pb_z[k] * mi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pb_z, li_42, mh0_69, mh0_70, \
                         mh0_72, mh1_69, mh1_70, mh1_72, mi_95, mi_96, \
                         mi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * mh0_69[k]
                   - f_4 * mh1_69[k]
                   + pb_z[k] * mi_95[k];

        t_126[k] = f_5 * mh0_70[k]
                   - f_6 * mh1_70[k]
                   + pb_z[k] * mi_96[k];

        t_127[k] = f_12 * li_42[k]
                   + pb_y[k] * mi_98[k];

        t_128[k] = f_9 * mh0_72[k]
                   - f_10 * mh1_72[k]
                   + pb_z[k] * mi_98[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pb_x, pb_z, li_105, li_107, \
                         li_108, li_109, mi_99, mi_105, mi_107, mi_108, \
                         mi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_17 * li_105[k]
                   + pb_x[k] * mi_105[k];

        t_130[k] = pb_z[k] * mi_99[k];

        t_131[k] = f_17 * li_107[k]
                   + pb_x[k] * mi_107[k];

        t_132[k] = f_17 * li_108[k]
                   + pb_x[k] * mi_108[k];

        t_133[k] = f_17 * li_109[k]
                   + pb_x[k] * mi_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_x, pb_x, pb_z, kk0_136, kk1_136, \
                         li_110, li_111, lk_136, mi_105, mi_110, \
                         mi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_17 * li_110[k]
                   + pb_x[k] * mi_110[k];

        t_135[k] = f_17 * li_111[k]
                   + pb_x[k] * mi_111[k];

        t_136[k] = f_20 * kk0_136[k]
                   - f_21 * kk1_136[k]
                   + pa_x[k] * lk_136[k];

        t_137[k] = pb_z[k] * mi_105[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pb_z, mh0_78, mh0_79, mh0_80, mh1_78, mh1_79, \
                         mh1_80, mi_106, mi_107, mi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_3 * mh0_78[k]
                   - f_4 * mh1_78[k]
                   + pb_z[k] * mi_106[k];

        t_139[k] = f_5 * mh0_79[k]
                   - f_6 * mh1_79[k]
                   + pb_z[k] * mi_107[k];

        t_140[k] = f_7 * mh0_80[k]
                   - f_8 * mh1_80[k]
                   + pb_z[k] * mi_108[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kk0,
                                            const size_t kk1, const size_t li, const size_t lk,
                                            const size_t mh0, const size_t mh1, const size_t mi,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_17 = 3.5 / p;
    const auto f_18 = 0.5 / alpha;
    const auto f_19 = 0.5 * beta / (alpha * p);
    const auto f_20 = 3.0 / alpha;
    const auto f_21 = 3.0 * beta / (alpha * p);
    const auto f_22 = 1.0 / alpha;
    const auto f_23 = beta / (alpha * p);
    const auto f_24 = 3.0 / p;
    const auto f_25 = 2.5 / alpha;
    const auto f_26 = 2.5 * beta / (alpha * p);

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
    auto *t_210 = buffer.data(target + 210);
    auto *t_211 = buffer.data(target + 211);
    auto *t_212 = buffer.data(target + 212);
    auto *t_213 = buffer.data(target + 213);
    auto *t_214 = buffer.data(target + 214);
    auto *t_215 = buffer.data(target + 215);
    auto *t_216 = buffer.data(target + 216);
    auto *t_217 = buffer.data(target + 217);
    auto *t_218 = buffer.data(target + 218);
    auto *t_219 = buffer.data(target + 219);
    auto *t_220 = buffer.data(target + 220);
    auto *t_221 = buffer.data(target + 221);
    auto *t_222 = buffer.data(target + 222);
    auto *t_223 = buffer.data(target + 223);
    auto *t_224 = buffer.data(target + 224);
    auto *t_225 = buffer.data(target + 225);
    auto *t_226 = buffer.data(target + 226);
    auto *t_227 = buffer.data(target + 227);
    auto *t_228 = buffer.data(target + 228);
    auto *t_229 = buffer.data(target + 229);
    auto *t_230 = buffer.data(target + 230);
    auto *t_231 = buffer.data(target + 231);
    auto *t_232 = buffer.data(target + 232);
    auto *t_233 = buffer.data(target + 233);
    auto *t_234 = buffer.data(target + 234);
    auto *t_235 = buffer.data(target + 235);
    auto *t_236 = buffer.data(target + 236);
    auto *t_237 = buffer.data(target + 237);
    auto *t_238 = buffer.data(target + 238);
    auto *t_239 = buffer.data(target + 239);
    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);
    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);
    auto *t_248 = buffer.data(target + 248);
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);
    auto *t_256 = buffer.data(target + 256);
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);
    auto *t_259 = buffer.data(target + 259);
    auto *t_260 = buffer.data(target + 260);
    auto *t_261 = buffer.data(target + 261);
    auto *t_262 = buffer.data(target + 262);
    auto *t_263 = buffer.data(target + 263);
    auto *t_264 = buffer.data(target + 264);
    auto *t_265 = buffer.data(target + 265);
    auto *t_266 = buffer.data(target + 266);
    auto *t_267 = buffer.data(target + 267);
    auto *t_268 = buffer.data(target + 268);
    auto *t_269 = buffer.data(target + 269);
    auto *t_270 = buffer.data(target + 270);
    auto *t_271 = buffer.data(target + 271);
    auto *t_272 = buffer.data(target + 272);
    auto *t_273 = buffer.data(target + 273);
    auto *t_274 = buffer.data(target + 274);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kk0_0 = buffer.data(kk0 + 0);
    const auto *kk0_36 = buffer.data(kk0 + 36);
    const auto *kk0_215 = buffer.data(kk0 + 215);
    const auto *kk0_244 = buffer.data(kk0 + 244);

    const auto *kk1_0 = buffer.data(kk1 + 0);
    const auto *kk1_36 = buffer.data(kk1 + 36);
    const auto *kk1_215 = buffer.data(kk1 + 215);
    const auto *kk1_244 = buffer.data(kk1 + 244);

    const auto *li_31 = buffer.data(li + 31);
    const auto *li_34 = buffer.data(li + 34);
    const auto *li_38 = buffer.data(li + 38);
    const auto *li_49 = buffer.data(li + 49);
    const auto *li_55 = buffer.data(li + 55);
    const auto *li_56 = buffer.data(li + 56);
    const auto *li_58 = buffer.data(li + 58);
    const auto *li_59 = buffer.data(li + 59);
    const auto *li_61 = buffer.data(li + 61);
    const auto *li_62 = buffer.data(li + 62);
    const auto *li_64 = buffer.data(li + 64);
    const auto *li_65 = buffer.data(li + 65);
    const auto *li_66 = buffer.data(li + 66);
    const auto *li_68 = buffer.data(li + 68);
    const auto *li_69 = buffer.data(li + 69);
    const auto *li_70 = buffer.data(li + 70);
    const auto *li_77 = buffer.data(li + 77);
    const auto *li_79 = buffer.data(li + 79);
    const auto *li_80 = buffer.data(li + 80);
    const auto *li_81 = buffer.data(li + 81);
    const auto *li_82 = buffer.data(li + 82);
    const auto *li_83 = buffer.data(li + 83);
    const auto *li_84 = buffer.data(li + 84);
    const auto *li_86 = buffer.data(li + 86);
    const auto *li_87 = buffer.data(li + 87);
    const auto *li_89 = buffer.data(li + 89);
    const auto *li_90 = buffer.data(li + 90);
    const auto *li_91 = buffer.data(li + 91);
    const auto *li_93 = buffer.data(li + 93);
    const auto *li_94 = buffer.data(li + 94);
    const auto *li_95 = buffer.data(li + 95);
    const auto *li_96 = buffer.data(li + 96);
    const auto *li_98 = buffer.data(li + 98);
    const auto *li_111 = buffer.data(li + 111);
    const auto *li_114 = buffer.data(li + 114);
    const auto *li_117 = buffer.data(li + 117);
    const auto *li_121 = buffer.data(li + 121);
    const auto *li_126 = buffer.data(li + 126);
    const auto *li_134 = buffer.data(li + 134);
    const auto *li_135 = buffer.data(li + 135);
    const auto *li_136 = buffer.data(li + 136);
    const auto *li_137 = buffer.data(li + 137);
    const auto *li_138 = buffer.data(li + 138);
    const auto *li_145 = buffer.data(li + 145);
    const auto *li_149 = buffer.data(li + 149);
    const auto *li_154 = buffer.data(li + 154);
    const auto *li_160 = buffer.data(li + 160);
    const auto *li_161 = buffer.data(li + 161);
    const auto *li_162 = buffer.data(li + 162);
    const auto *li_163 = buffer.data(li + 163);
    const auto *li_164 = buffer.data(li + 164);
    const auto *li_165 = buffer.data(li + 165);
    const auto *li_167 = buffer.data(li + 167);
    const auto *li_171 = buffer.data(li + 171);
    const auto *li_174 = buffer.data(li + 174);
    const auto *li_178 = buffer.data(li + 178);
    const auto *li_183 = buffer.data(li + 183);
    const auto *li_189 = buffer.data(li + 189);
    const auto *li_191 = buffer.data(li + 191);
    const auto *li_192 = buffer.data(li + 192);
    const auto *li_193 = buffer.data(li + 193);
    const auto *li_194 = buffer.data(li + 194);
    const auto *li_195 = buffer.data(li + 195);
    const auto *li_218 = buffer.data(li + 218);

    const auto *lk_37 = buffer.data(lk + 37);
    const auto *lk_39 = buffer.data(lk + 39);
    const auto *lk_42 = buffer.data(lk + 42);
    const auto *lk_46 = buffer.data(lk + 46);
    const auto *lk_51 = buffer.data(lk + 51);
    const auto *lk_57 = buffer.data(lk + 57);
    const auto *lk_64 = buffer.data(lk + 64);
    const auto *lk_72 = buffer.data(lk + 72);
    const auto *lk_74 = buffer.data(lk + 74);
    const auto *lk_77 = buffer.data(lk + 77);
    const auto *lk_81 = buffer.data(lk + 81);
    const auto *lk_84 = buffer.data(lk + 84);
    const auto *lk_86 = buffer.data(lk + 86);
    const auto *lk_89 = buffer.data(lk + 89);
    const auto *lk_90 = buffer.data(lk + 90);
    const auto *lk_92 = buffer.data(lk + 92);
    const auto *lk_99 = buffer.data(lk + 99);
    const auto *lk_102 = buffer.data(lk + 102);
    const auto *lk_103 = buffer.data(lk + 103);
    const auto *lk_104 = buffer.data(lk + 104);
    const auto *lk_105 = buffer.data(lk + 105);
    const auto *lk_107 = buffer.data(lk + 107);
    const auto *lk_108 = buffer.data(lk + 108);
    const auto *lk_109 = buffer.data(lk + 109);
    const auto *lk_111 = buffer.data(lk + 111);
    const auto *lk_113 = buffer.data(lk + 113);
    const auto *lk_114 = buffer.data(lk + 114);
    const auto *lk_117 = buffer.data(lk + 117);
    const auto *lk_118 = buffer.data(lk + 118);
    const auto *lk_120 = buffer.data(lk + 120);
    const auto *lk_122 = buffer.data(lk + 122);
    const auto *lk_123 = buffer.data(lk + 123);
    const auto *lk_125 = buffer.data(lk + 125);
    const auto *lk_126 = buffer.data(lk + 126);
    const auto *lk_128 = buffer.data(lk + 128);
    const auto *lk_129 = buffer.data(lk + 129);
    const auto *lk_215 = buffer.data(lk + 215);
    const auto *lk_244 = buffer.data(lk + 244);

    const auto *mh0_81 = buffer.data(mh0 + 81);
    const auto *mh0_83 = buffer.data(mh0 + 83);
    const auto *mh0_105 = buffer.data(mh0 + 105);
    const auto *mh0_106 = buffer.data(mh0 + 106);
    const auto *mh0_108 = buffer.data(mh0 + 108);
    const auto *mh0_110 = buffer.data(mh0 + 110);
    const auto *mh0_111 = buffer.data(mh0 + 111);
    const auto *mh0_113 = buffer.data(mh0 + 113);
    const auto *mh0_114 = buffer.data(mh0 + 114);
    const auto *mh0_119 = buffer.data(mh0 + 119);
    const auto *mh0_120 = buffer.data(mh0 + 120);
    const auto *mh0_122 = buffer.data(mh0 + 122);
    const auto *mh0_123 = buffer.data(mh0 + 123);
    const auto *mh0_124 = buffer.data(mh0 + 124);
    const auto *mh0_125 = buffer.data(mh0 + 125);
    const auto *mh0_126 = buffer.data(mh0 + 126);
    const auto *mh0_128 = buffer.data(mh0 + 128);
    const auto *mh0_129 = buffer.data(mh0 + 129);
    const auto *mh0_131 = buffer.data(mh0 + 131);
    const auto *mh0_132 = buffer.data(mh0 + 132);
    const auto *mh0_133 = buffer.data(mh0 + 133);
    const auto *mh0_135 = buffer.data(mh0 + 135);
    const auto *mh0_136 = buffer.data(mh0 + 136);
    const auto *mh0_141 = buffer.data(mh0 + 141);
    const auto *mh0_142 = buffer.data(mh0 + 142);
    const auto *mh0_143 = buffer.data(mh0 + 143);
    const auto *mh0_144 = buffer.data(mh0 + 144);
    const auto *mh0_146 = buffer.data(mh0 + 146);

    const auto *mh1_81 = buffer.data(mh1 + 81);
    const auto *mh1_83 = buffer.data(mh1 + 83);
    const auto *mh1_105 = buffer.data(mh1 + 105);
    const auto *mh1_106 = buffer.data(mh1 + 106);
    const auto *mh1_108 = buffer.data(mh1 + 108);
    const auto *mh1_110 = buffer.data(mh1 + 110);
    const auto *mh1_111 = buffer.data(mh1 + 111);
    const auto *mh1_113 = buffer.data(mh1 + 113);
    const auto *mh1_114 = buffer.data(mh1 + 114);
    const auto *mh1_119 = buffer.data(mh1 + 119);
    const auto *mh1_120 = buffer.data(mh1 + 120);
    const auto *mh1_122 = buffer.data(mh1 + 122);
    const auto *mh1_123 = buffer.data(mh1 + 123);
    const auto *mh1_124 = buffer.data(mh1 + 124);
    const auto *mh1_125 = buffer.data(mh1 + 125);
    const auto *mh1_126 = buffer.data(mh1 + 126);
    const auto *mh1_128 = buffer.data(mh1 + 128);
    const auto *mh1_129 = buffer.data(mh1 + 129);
    const auto *mh1_131 = buffer.data(mh1 + 131);
    const auto *mh1_132 = buffer.data(mh1 + 132);
    const auto *mh1_133 = buffer.data(mh1 + 133);
    const auto *mh1_135 = buffer.data(mh1 + 135);
    const auto *mh1_136 = buffer.data(mh1 + 136);
    const auto *mh1_141 = buffer.data(mh1 + 141);
    const auto *mh1_142 = buffer.data(mh1 + 142);
    const auto *mh1_143 = buffer.data(mh1 + 143);
    const auto *mh1_144 = buffer.data(mh1 + 144);
    const auto *mh1_146 = buffer.data(mh1 + 146);

    const auto *mi_109 = buffer.data(mi + 109);
    const auto *mi_111 = buffer.data(mi + 111);
    const auto *mi_114 = buffer.data(mi + 114);
    const auto *mi_115 = buffer.data(mi + 115);
    const auto *mi_117 = buffer.data(mi + 117);
    const auto *mi_118 = buffer.data(mi + 118);
    const auto *mi_121 = buffer.data(mi + 121);
    const auto *mi_122 = buffer.data(mi + 122);
    const auto *mi_126 = buffer.data(mi + 126);
    const auto *mi_133 = buffer.data(mi + 133);
    const auto *mi_134 = buffer.data(mi + 134);
    const auto *mi_135 = buffer.data(mi + 135);
    const auto *mi_136 = buffer.data(mi + 136);
    const auto *mi_137 = buffer.data(mi + 137);
    const auto *mi_138 = buffer.data(mi + 138);
    const auto *mi_139 = buffer.data(mi + 139);
    const auto *mi_140 = buffer.data(mi + 140);
    const auto *mi_141 = buffer.data(mi + 141);
    const auto *mi_142 = buffer.data(mi + 142);
    const auto *mi_143 = buffer.data(mi + 143);
    const auto *mi_145 = buffer.data(mi + 145);
    const auto *mi_146 = buffer.data(mi + 146);
    const auto *mi_148 = buffer.data(mi + 148);
    const auto *mi_149 = buffer.data(mi + 149);
    const auto *mi_150 = buffer.data(mi + 150);
    const auto *mi_152 = buffer.data(mi + 152);
    const auto *mi_153 = buffer.data(mi + 153);
    const auto *mi_154 = buffer.data(mi + 154);
    const auto *mi_160 = buffer.data(mi + 160);
    const auto *mi_161 = buffer.data(mi + 161);
    const auto *mi_162 = buffer.data(mi + 162);
    const auto *mi_163 = buffer.data(mi + 163);
    const auto *mi_164 = buffer.data(mi + 164);
    const auto *mi_165 = buffer.data(mi + 165);
    const auto *mi_166 = buffer.data(mi + 166);
    const auto *mi_167 = buffer.data(mi + 167);
    const auto *mi_168 = buffer.data(mi + 168);
    const auto *mi_169 = buffer.data(mi + 169);
    const auto *mi_170 = buffer.data(mi + 170);
    const auto *mi_171 = buffer.data(mi + 171);
    const auto *mi_173 = buffer.data(mi + 173);
    const auto *mi_174 = buffer.data(mi + 174);
    const auto *mi_175 = buffer.data(mi + 175);
    const auto *mi_177 = buffer.data(mi + 177);
    const auto *mi_178 = buffer.data(mi + 178);
    const auto *mi_179 = buffer.data(mi + 179);
    const auto *mi_180 = buffer.data(mi + 180);
    const auto *mi_182 = buffer.data(mi + 182);
    const auto *mi_183 = buffer.data(mi + 183);
    const auto *mi_189 = buffer.data(mi + 189);
    const auto *mi_190 = buffer.data(mi + 190);
    const auto *mi_191 = buffer.data(mi + 191);
    const auto *mi_192 = buffer.data(mi + 192);
    const auto *mi_193 = buffer.data(mi + 193);
    const auto *mi_194 = buffer.data(mi + 194);
    const auto *mi_195 = buffer.data(mi + 195);
    const auto *mi_196 = buffer.data(mi + 196);
    const auto *mi_198 = buffer.data(mi + 198);
    const auto *mi_199 = buffer.data(mi + 199);
    const auto *mi_201 = buffer.data(mi + 201);
    const auto *mi_202 = buffer.data(mi + 202);
    const auto *mi_205 = buffer.data(mi + 205);
    const auto *mi_206 = buffer.data(mi + 206);
    const auto *mi_210 = buffer.data(mi + 210);
    const auto *mi_218 = buffer.data(mi + 218);

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pb_y, pb_z, li_55, lk_72, mh0_81, \
                         mh0_83, mh1_81, mh1_83, mi_109, mi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_9 * mh0_81[k]
                   - f_10 * mh1_81[k]
                   + pb_z[k] * mi_109[k];

        t_142[k] = f_12 * li_55[k]
                   + pb_y[k] * mi_111[k];

        t_143[k] = f_1 * mh0_83[k]
                   - f_2 * mh1_83[k]
                   + pb_z[k] * mi_111[k];

        t_144[k] = pa_y[k] * lk_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, t_150, pa_y, pa_z, pb_y, li_58, \
                         lk_37, lk_39, lk_42, lk_74, lk_77, mi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = pa_z[k] * lk_37[k];

        t_146[k] = pa_y[k] * lk_74[k];

        t_147[k] = pa_z[k] * lk_39[k];

        t_148[k] = f_11 * li_58[k]
                   + pb_y[k] * mi_114[k];

        t_149[k] = pa_y[k] * lk_77[k];

        t_150[k] = pa_z[k] * lk_42[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pa_y, pa_z, pb_y, pb_z, li_31, li_61, \
                         lk_46, lk_81, mi_115, mi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_11 * li_31[k]
                   + pb_z[k] * mi_115[k];

        t_152[k] = f_11 * li_61[k]
                   + pb_y[k] * mi_117[k];

        t_153[k] = pa_y[k] * lk_81[k];

        t_154[k] = pa_z[k] * lk_46[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_y, pb_y, pb_z, li_34, li_64, li_65, \
                         lk_84, lk_86, mi_118, mi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_11 * li_34[k]
                   + pb_z[k] * mi_118[k];

        t_156[k] = f_12 * li_64[k]
                   + pa_y[k] * lk_84[k];

        t_157[k] = f_11 * li_65[k]
                   + pb_y[k] * mi_121[k];

        t_158[k] = pa_y[k] * lk_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_y, pa_z, pb_z, li_38, li_68, li_69, \
                         lk_51, lk_89, lk_90, mi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * lk_51[k];

        t_160[k] = f_11 * li_38[k]
                   + pb_z[k] * mi_122[k];

        t_161[k] = f_13 * li_68[k]
                   + pa_y[k] * lk_89[k];

        t_162[k] = f_12 * li_69[k]
                   + pa_y[k] * lk_90[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_y, pa_z, pb_x, pb_y, li_70, li_134, \
                         lk_57, lk_92, mi_126, mi_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_11 * li_70[k]
                   + pb_y[k] * mi_126[k];

        t_164[k] = pa_y[k] * lk_92[k];

        t_165[k] = pa_z[k] * lk_57[k];

        t_166[k] = f_17 * li_134[k]
                   + pb_x[k] * mi_134[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pa_y, pb_x, li_135, li_136, \
                         li_137, li_138, lk_99, mi_135, mi_136, mi_137, \
                         mi_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_17 * li_135[k]
                   + pb_x[k] * mi_135[k];

        t_168[k] = f_17 * li_136[k]
                   + pb_x[k] * mi_136[k];

        t_169[k] = f_17 * li_137[k]
                   + pb_x[k] * mi_137[k];

        t_170[k] = f_17 * li_138[k]
                   + pb_x[k] * mi_138[k];

        t_171[k] = pa_y[k] * lk_99[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pa_z, pb_z, li_49, li_79, li_80, \
                         lk_64, lk_102, lk_103, mi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pa_z[k] * lk_64[k];

        t_173[k] = f_11 * li_49[k]
                   + pb_z[k] * mi_133[k];

        t_174[k] = f_15 * li_79[k]
                   + pa_y[k] * lk_102[k];

        t_175[k] = f_14 * li_80[k]
                   + pa_y[k] * lk_103[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_y, li_81, li_82, li_83, lk_104, \
                         lk_105, lk_107, mi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * li_81[k]
                   + pa_y[k] * lk_104[k];

        t_177[k] = f_12 * li_82[k]
                   + pa_y[k] * lk_105[k];

        t_178[k] = f_11 * li_83[k]
                   + pb_y[k] * mi_139[k];

        t_179[k] = pa_y[k] * lk_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_z, pb_y, pb_z, kk0_0, kk1_0, li_56, \
                         lk_72, mh0_105, mh1_105, mi_140, mi_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_18 * kk0_0[k]
                   - f_19 * kk1_0[k]
                   + pa_z[k] * lk_72[k];

        t_181[k] = pb_y[k] * mi_140[k];

        t_182[k] = f_12 * li_56[k]
                   + pb_z[k] * mi_140[k];

        t_183[k] = f_3 * mh0_105[k]
                   - f_4 * mh1_105[k]
                   + pb_y[k] * mi_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pb_x, pb_y, pb_z, li_59, li_145, mh0_106, \
                         mh0_110, mh1_106, mh1_110, mi_142, mi_143, \
                         mi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = pb_y[k] * mi_142[k];

        t_185[k] = f_17 * li_145[k]
                   + f_9 * mh0_110[k]
                   - f_10 * mh1_110[k]
                   + pb_x[k] * mi_145[k];

        t_186[k] = f_5 * mh0_106[k]
                   - f_6 * mh1_106[k]
                   + pb_y[k] * mi_143[k];

        t_187[k] = f_12 * li_59[k]
                   + pb_z[k] * mi_143[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pb_x, pb_y, pb_z, li_62, li_149, mh0_108, \
                         mh0_114, mh1_108, mh1_114, mi_145, mi_146, \
                         mi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pb_y[k] * mi_145[k];

        t_189[k] = f_17 * li_149[k]
                   + f_7 * mh0_114[k]
                   - f_8 * mh1_114[k]
                   + pb_x[k] * mi_149[k];

        t_190[k] = f_7 * mh0_108[k]
                   - f_8 * mh1_108[k]
                   + pb_y[k] * mi_146[k];

        t_191[k] = f_12 * li_62[k]
                   + pb_z[k] * mi_146[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_x, pb_y, li_154, mh0_110, mh0_119, mh1_110, \
                         mh1_119, mi_148, mi_149, mi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_3 * mh0_110[k]
                   - f_4 * mh1_110[k]
                   + pb_y[k] * mi_148[k];

        t_193[k] = pb_y[k] * mi_149[k];

        t_194[k] = f_17 * li_154[k]
                   + f_5 * mh0_119[k]
                   - f_6 * mh1_119[k]
                   + pb_x[k] * mi_154[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pb_y, pb_z, li_66, mh0_111, mh0_113, \
                         mh0_114, mh1_111, mh1_113, mh1_114, mi_150, mi_152, \
                         mi_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_9 * mh0_111[k]
                   - f_10 * mh1_111[k]
                   + pb_y[k] * mi_150[k];

        t_196[k] = f_12 * li_66[k]
                   + pb_z[k] * mi_150[k];

        t_197[k] = f_5 * mh0_113[k]
                   - f_6 * mh1_113[k]
                   + pb_y[k] * mi_152[k];

        t_198[k] = f_3 * mh0_114[k]
                   - f_4 * mh1_114[k]
                   + pb_y[k] * mi_153[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pb_x, pb_y, li_160, li_161, li_162, \
                         mh0_125, mh1_125, mi_154, mi_160, mi_161, \
                         mi_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = pb_y[k] * mi_154[k];

        t_200[k] = f_17 * li_160[k]
                   + f_3 * mh0_125[k]
                   - f_4 * mh1_125[k]
                   + pb_x[k] * mi_160[k];

        t_201[k] = f_17 * li_161[k]
                   + pb_x[k] * mi_161[k];

        t_202[k] = f_17 * li_162[k]
                   + pb_x[k] * mi_162[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pb_x, pb_y, li_163, li_164, \
                         li_165, li_167, mi_160, mi_163, mi_164, mi_165, \
                         mi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_17 * li_163[k]
                   + pb_x[k] * mi_163[k];

        t_204[k] = f_17 * li_164[k]
                   + pb_x[k] * mi_164[k];

        t_205[k] = f_17 * li_165[k]
                   + pb_x[k] * mi_165[k];

        t_206[k] = pb_y[k] * mi_160[k];

        t_207[k] = f_17 * li_167[k]
                   + pb_x[k] * mi_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_y, pb_z, li_77, mh0_120, mh0_122, \
                         mh0_123, mh1_120, mh1_122, mh1_123, mi_161, mi_163, \
                         mi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * mh0_120[k]
                   - f_2 * mh1_120[k]
                   + pb_y[k] * mi_161[k];

        t_209[k] = f_12 * li_77[k]
                   + pb_z[k] * mi_161[k];

        t_210[k] = f_9 * mh0_122[k]
                   - f_10 * mh1_122[k]
                   + pb_y[k] * mi_163[k];

        t_211[k] = f_7 * mh0_123[k]
                   - f_8 * mh1_123[k]
                   + pb_y[k] * mi_164[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pa_x, pb_y, kk0_215, kk1_215, lk_215, \
                         mh0_124, mh0_125, mh1_124, mh1_125, mi_165, mi_166, \
                         mi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_5 * mh0_124[k]
                   - f_6 * mh1_124[k]
                   + pb_y[k] * mi_165[k];

        t_213[k] = f_3 * mh0_125[k]
                   - f_4 * mh1_125[k]
                   + pb_y[k] * mi_166[k];

        t_214[k] = pb_y[k] * mi_167[k];

        t_215[k] = f_20 * kk0_215[k]
                   - f_21 * kk1_215[k]
                   + pa_x[k] * lk_215[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, pa_y, pb_y, pb_z, kk0_36, kk1_36, li_84, lk_108, \
                         mi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_22 * kk0_36[k]
                   - f_23 * kk1_36[k]
                   + pa_y[k] * lk_108[k];

        t_217[k] = f_13 * li_84[k]
                   + pb_y[k] * mi_168[k];

        t_218[k] = pb_z[k] * mi_168[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pb_x, pb_z, li_171, mh0_126, mh0_129, mh1_126, \
                         mh1_129, mi_169, mi_170, mi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_24 * li_171[k]
                   + f_9 * mh0_129[k]
                   - f_10 * mh1_129[k]
                   + pb_x[k] * mi_171[k];

        t_220[k] = pb_z[k] * mi_169[k];

        t_221[k] = f_3 * mh0_126[k]
                   - f_4 * mh1_126[k]
                   + pb_z[k] * mi_170[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pb_x, pb_y, pb_z, li_89, li_174, mh0_128, \
                         mh0_132, mh1_128, mh1_132, mi_171, mi_173, \
                         mi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_24 * li_174[k]
                   + f_7 * mh0_132[k]
                   - f_8 * mh1_132[k]
                   + pb_x[k] * mi_174[k];

        t_223[k] = pb_z[k] * mi_171[k];

        t_224[k] = f_13 * li_89[k]
                   + pb_y[k] * mi_173[k];

        t_225[k] = f_5 * mh0_128[k]
                   - f_6 * mh1_128[k]
                   + pb_z[k] * mi_173[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pb_x, pb_z, li_178, mh0_129, mh0_136, mh1_129, \
                         mh1_136, mi_174, mi_175, mi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_24 * li_178[k]
                   + f_5 * mh0_136[k]
                   - f_6 * mh1_136[k]
                   + pb_x[k] * mi_178[k];

        t_227[k] = pb_z[k] * mi_174[k];

        t_228[k] = f_3 * mh0_129[k]
                   - f_4 * mh1_129[k]
                   + pb_z[k] * mi_175[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_x, pb_y, pb_z, li_93, li_183, mh0_131, \
                         mh0_141, mh1_131, mh1_141, mi_177, mi_178, \
                         mi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_13 * li_93[k]
                   + pb_y[k] * mi_177[k];

        t_230[k] = f_7 * mh0_131[k]
                   - f_8 * mh1_131[k]
                   + pb_z[k] * mi_177[k];

        t_231[k] = f_24 * li_183[k]
                   + f_3 * mh0_141[k]
                   - f_4 * mh1_141[k]
                   + pb_x[k] * mi_183[k];

        t_232[k] = pb_z[k] * mi_178[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pb_y, pb_z, li_98, mh0_132, mh0_133, \
                         mh0_135, mh1_132, mh1_133, mh1_135, mi_179, mi_180, \
                         mi_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * mh0_132[k]
                   - f_4 * mh1_132[k]
                   + pb_z[k] * mi_179[k];

        t_234[k] = f_5 * mh0_133[k]
                   - f_6 * mh1_133[k]
                   + pb_z[k] * mi_180[k];

        t_235[k] = f_13 * li_98[k]
                   + pb_y[k] * mi_182[k];

        t_236[k] = f_9 * mh0_135[k]
                   - f_10 * mh1_135[k]
                   + pb_z[k] * mi_182[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, pb_x, pb_z, li_189, li_191, \
                         li_192, li_193, mi_183, mi_189, mi_191, mi_192, \
                         mi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_24 * li_189[k]
                   + pb_x[k] * mi_189[k];

        t_238[k] = pb_z[k] * mi_183[k];

        t_239[k] = f_24 * li_191[k]
                   + pb_x[k] * mi_191[k];

        t_240[k] = f_24 * li_192[k]
                   + pb_x[k] * mi_192[k];

        t_241[k] = f_24 * li_193[k]
                   + pb_x[k] * mi_193[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pa_x, pb_x, pb_z, kk0_244, kk1_244, \
                         li_194, li_195, lk_244, mi_189, mi_194, \
                         mi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_24 * li_194[k]
                   + pb_x[k] * mi_194[k];

        t_243[k] = f_24 * li_195[k]
                   + pb_x[k] * mi_195[k];

        t_244[k] = f_25 * kk0_244[k]
                   - f_26 * kk1_244[k]
                   + pa_x[k] * lk_244[k];

        t_245[k] = pb_z[k] * mi_189[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pb_z, mh0_141, mh0_142, mh0_143, mh1_141, \
                         mh1_142, mh1_143, mi_190, mi_191, mi_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_3 * mh0_141[k]
                   - f_4 * mh1_141[k]
                   + pb_z[k] * mi_190[k];

        t_247[k] = f_5 * mh0_142[k]
                   - f_6 * mh1_142[k]
                   + pb_z[k] * mi_191[k];

        t_248[k] = f_7 * mh0_143[k]
                   - f_8 * mh1_143[k]
                   + pb_z[k] * mi_192[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, pa_z, pb_y, pb_z, li_111, lk_108, \
                         mh0_144, mh0_146, mh1_144, mh1_146, mi_193, \
                         mi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_9 * mh0_144[k]
                   - f_10 * mh1_144[k]
                   + pb_z[k] * mi_193[k];

        t_250[k] = f_13 * li_111[k]
                   + pb_y[k] * mi_195[k];

        t_251[k] = f_1 * mh0_146[k]
                   - f_2 * mh1_146[k]
                   + pb_z[k] * mi_195[k];

        t_252[k] = pa_z[k] * lk_108[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, pa_z, pb_y, pb_z, li_84, li_86, \
                         li_114, lk_109, lk_111, lk_113, mi_196, \
                         mi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = pa_z[k] * lk_109[k];

        t_254[k] = f_11 * li_84[k]
                   + pb_z[k] * mi_196[k];

        t_255[k] = pa_z[k] * lk_111[k];

        t_256[k] = f_12 * li_114[k]
                   + pb_y[k] * mi_198[k];

        t_257[k] = f_12 * li_86[k]
                   + pa_z[k] * lk_113[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, pa_z, pb_y, pb_z, li_87, li_89, \
                         li_117, lk_114, lk_117, lk_118, mi_199, \
                         mi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = pa_z[k] * lk_114[k];

        t_259[k] = f_11 * li_87[k]
                   + pb_z[k] * mi_199[k];

        t_260[k] = f_12 * li_117[k]
                   + pb_y[k] * mi_201[k];

        t_261[k] = f_13 * li_89[k]
                   + pa_z[k] * lk_117[k];

        t_262[k] = pa_z[k] * lk_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, pa_z, pb_y, pb_z, li_90, li_91, li_93, \
                         li_121, lk_120, lk_122, mi_202, mi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_11 * li_90[k]
                   + pb_z[k] * mi_202[k];

        t_264[k] = f_12 * li_91[k]
                   + pa_z[k] * lk_120[k];

        t_265[k] = f_12 * li_121[k]
                   + pb_y[k] * mi_205[k];

        t_266[k] = f_14 * li_93[k]
                   + pa_z[k] * lk_122[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pa_z, pb_z, li_94, li_95, li_96, lk_123, \
                         lk_125, lk_126, mi_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = pa_z[k] * lk_123[k];

        t_268[k] = f_11 * li_94[k]
                   + pb_z[k] * mi_206[k];

        t_269[k] = f_12 * li_95[k]
                   + pa_z[k] * lk_125[k];

        t_270[k] = f_13 * li_96[k]
                   + pa_z[k] * lk_126[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pa_z, pb_x, pb_y, li_98, li_126, li_218, \
                         lk_128, lk_129, mi_210, mi_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_12 * li_126[k]
                   + pb_y[k] * mi_210[k];

        t_272[k] = f_15 * li_98[k]
                   + pa_z[k] * lk_128[k];

        t_273[k] = pa_z[k] * lk_129[k];

        t_274[k] = f_24 * li_218[k]
                   + pb_x[k] * mi_218[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kk0,
                                            const size_t kk1, const size_t li, const size_t lk,
                                            const size_t mh0, const size_t mh1, const size_t mi,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_17 = 3.5 / p;
    const auto f_22 = 1.0 / alpha;
    const auto f_23 = beta / (alpha * p);
    const auto f_24 = 3.0 / p;
    const auto f_25 = 2.5 / alpha;
    const auto f_26 = 2.5 * beta / (alpha * p);
    const auto f_27 = 1.5 / alpha;
    const auto f_28 = 1.5 * beta / (alpha * p);
    const auto f_29 = 2.0 / alpha;
    const auto f_30 = 2.0 * beta / (alpha * p);

    auto *t_275 = buffer.data(target + 275);
    auto *t_276 = buffer.data(target + 276);
    auto *t_277 = buffer.data(target + 277);
    auto *t_278 = buffer.data(target + 278);
    auto *t_279 = buffer.data(target + 279);
    auto *t_280 = buffer.data(target + 280);
    auto *t_281 = buffer.data(target + 281);
    auto *t_282 = buffer.data(target + 282);
    auto *t_283 = buffer.data(target + 283);
    auto *t_284 = buffer.data(target + 284);
    auto *t_285 = buffer.data(target + 285);
    auto *t_286 = buffer.data(target + 286);
    auto *t_287 = buffer.data(target + 287);
    auto *t_288 = buffer.data(target + 288);
    auto *t_289 = buffer.data(target + 289);
    auto *t_290 = buffer.data(target + 290);
    auto *t_291 = buffer.data(target + 291);
    auto *t_292 = buffer.data(target + 292);
    auto *t_293 = buffer.data(target + 293);
    auto *t_294 = buffer.data(target + 294);
    auto *t_295 = buffer.data(target + 295);
    auto *t_296 = buffer.data(target + 296);
    auto *t_297 = buffer.data(target + 297);
    auto *t_298 = buffer.data(target + 298);
    auto *t_299 = buffer.data(target + 299);
    auto *t_300 = buffer.data(target + 300);
    auto *t_301 = buffer.data(target + 301);
    auto *t_302 = buffer.data(target + 302);
    auto *t_303 = buffer.data(target + 303);
    auto *t_304 = buffer.data(target + 304);
    auto *t_305 = buffer.data(target + 305);
    auto *t_306 = buffer.data(target + 306);
    auto *t_307 = buffer.data(target + 307);
    auto *t_308 = buffer.data(target + 308);
    auto *t_309 = buffer.data(target + 309);
    auto *t_310 = buffer.data(target + 310);
    auto *t_311 = buffer.data(target + 311);
    auto *t_312 = buffer.data(target + 312);
    auto *t_313 = buffer.data(target + 313);
    auto *t_314 = buffer.data(target + 314);
    auto *t_315 = buffer.data(target + 315);
    auto *t_316 = buffer.data(target + 316);
    auto *t_317 = buffer.data(target + 317);
    auto *t_318 = buffer.data(target + 318);
    auto *t_319 = buffer.data(target + 319);
    auto *t_320 = buffer.data(target + 320);
    auto *t_321 = buffer.data(target + 321);
    auto *t_322 = buffer.data(target + 322);
    auto *t_323 = buffer.data(target + 323);
    auto *t_324 = buffer.data(target + 324);
    auto *t_325 = buffer.data(target + 325);
    auto *t_326 = buffer.data(target + 326);
    auto *t_327 = buffer.data(target + 327);
    auto *t_328 = buffer.data(target + 328);
    auto *t_329 = buffer.data(target + 329);
    auto *t_330 = buffer.data(target + 330);
    auto *t_331 = buffer.data(target + 331);
    auto *t_332 = buffer.data(target + 332);
    auto *t_333 = buffer.data(target + 333);
    auto *t_334 = buffer.data(target + 334);
    auto *t_335 = buffer.data(target + 335);
    auto *t_336 = buffer.data(target + 336);
    auto *t_337 = buffer.data(target + 337);
    auto *t_338 = buffer.data(target + 338);
    auto *t_339 = buffer.data(target + 339);
    auto *t_340 = buffer.data(target + 340);
    auto *t_341 = buffer.data(target + 341);
    auto *t_342 = buffer.data(target + 342);
    auto *t_343 = buffer.data(target + 343);
    auto *t_344 = buffer.data(target + 344);
    auto *t_345 = buffer.data(target + 345);
    auto *t_346 = buffer.data(target + 346);
    auto *t_347 = buffer.data(target + 347);
    auto *t_348 = buffer.data(target + 348);
    auto *t_349 = buffer.data(target + 349);
    auto *t_350 = buffer.data(target + 350);
    auto *t_351 = buffer.data(target + 351);
    auto *t_352 = buffer.data(target + 352);
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);
    auto *t_357 = buffer.data(target + 357);
    auto *t_358 = buffer.data(target + 358);
    auto *t_359 = buffer.data(target + 359);
    auto *t_360 = buffer.data(target + 360);
    auto *t_361 = buffer.data(target + 361);
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);
    auto *t_378 = buffer.data(target + 378);
    auto *t_379 = buffer.data(target + 379);
    auto *t_380 = buffer.data(target + 380);
    auto *t_381 = buffer.data(target + 381);
    auto *t_382 = buffer.data(target + 382);
    auto *t_383 = buffer.data(target + 383);
    auto *t_384 = buffer.data(target + 384);
    auto *t_385 = buffer.data(target + 385);
    auto *t_386 = buffer.data(target + 386);
    auto *t_387 = buffer.data(target + 387);
    auto *t_388 = buffer.data(target + 388);
    auto *t_389 = buffer.data(target + 389);
    auto *t_390 = buffer.data(target + 390);
    auto *t_391 = buffer.data(target + 391);
    auto *t_392 = buffer.data(target + 392);
    auto *t_393 = buffer.data(target + 393);
    auto *t_394 = buffer.data(target + 394);
    auto *t_395 = buffer.data(target + 395);
    auto *t_396 = buffer.data(target + 396);
    auto *t_397 = buffer.data(target + 397);
    auto *t_398 = buffer.data(target + 398);
    auto *t_399 = buffer.data(target + 399);
    auto *t_400 = buffer.data(target + 400);
    auto *t_401 = buffer.data(target + 401);
    auto *t_402 = buffer.data(target + 402);
    auto *t_403 = buffer.data(target + 403);
    auto *t_404 = buffer.data(target + 404);
    auto *t_405 = buffer.data(target + 405);
    auto *t_406 = buffer.data(target + 406);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kk0_72 = buffer.data(kk0 + 72);
    const auto *kk0_108 = buffer.data(kk0 + 108);
    const auto *kk0_359 = buffer.data(kk0 + 359);
    const auto *kk0_388 = buffer.data(kk0 + 388);

    const auto *kk1_72 = buffer.data(kk1 + 72);
    const auto *kk1_108 = buffer.data(kk1 + 108);
    const auto *kk1_359 = buffer.data(kk1 + 359);
    const auto *kk1_388 = buffer.data(kk1 + 388);

    const auto *li_105 = buffer.data(li + 105);
    const auto *li_106 = buffer.data(li + 106);
    const auto *li_107 = buffer.data(li + 107);
    const auto *li_108 = buffer.data(li + 108);
    const auto *li_109 = buffer.data(li + 109);
    const auto *li_111 = buffer.data(li + 111);
    const auto *li_115 = buffer.data(li + 115);
    const auto *li_118 = buffer.data(li + 118);
    const auto *li_122 = buffer.data(li + 122);
    const auto *li_133 = buffer.data(li + 133);
    const auto *li_139 = buffer.data(li + 139);
    const auto *li_140 = buffer.data(li + 140);
    const auto *li_141 = buffer.data(li + 141);
    const auto *li_142 = buffer.data(li + 142);
    const auto *li_143 = buffer.data(li + 143);
    const auto *li_145 = buffer.data(li + 145);
    const auto *li_146 = buffer.data(li + 146);
    const auto *li_148 = buffer.data(li + 148);
    const auto *li_149 = buffer.data(li + 149);
    const auto *li_150 = buffer.data(li + 150);
    const auto *li_152 = buffer.data(li + 152);
    const auto *li_153 = buffer.data(li + 153);
    const auto *li_154 = buffer.data(li + 154);
    const auto *li_161 = buffer.data(li + 161);
    const auto *li_163 = buffer.data(li + 163);
    const auto *li_164 = buffer.data(li + 164);
    const auto *li_165 = buffer.data(li + 165);
    const auto *li_166 = buffer.data(li + 166);
    const auto *li_167 = buffer.data(li + 167);
    const auto *li_168 = buffer.data(li + 168);
    const auto *li_170 = buffer.data(li + 170);
    const auto *li_171 = buffer.data(li + 171);
    const auto *li_173 = buffer.data(li + 173);
    const auto *li_177 = buffer.data(li + 177);
    const auto *li_182 = buffer.data(li + 182);
    const auto *li_195 = buffer.data(li + 195);
    const auto *li_198 = buffer.data(li + 198);
    const auto *li_201 = buffer.data(li + 201);
    const auto *li_219 = buffer.data(li + 219);
    const auto *li_220 = buffer.data(li + 220);
    const auto *li_221 = buffer.data(li + 221);
    const auto *li_222 = buffer.data(li + 222);
    const auto *li_223 = buffer.data(li + 223);
    const auto *li_245 = buffer.data(li + 245);
    const auto *li_246 = buffer.data(li + 246);
    const auto *li_247 = buffer.data(li + 247);
    const auto *li_248 = buffer.data(li + 248);
    const auto *li_249 = buffer.data(li + 249);
    const auto *li_250 = buffer.data(li + 250);
    const auto *li_257 = buffer.data(li + 257);
    const auto *li_261 = buffer.data(li + 261);
    const auto *li_266 = buffer.data(li + 266);
    const auto *li_272 = buffer.data(li + 272);
    const auto *li_273 = buffer.data(li + 273);
    const auto *li_274 = buffer.data(li + 274);
    const auto *li_275 = buffer.data(li + 275);
    const auto *li_276 = buffer.data(li + 276);
    const auto *li_277 = buffer.data(li + 277);
    const auto *li_279 = buffer.data(li + 279);
    const auto *li_283 = buffer.data(li + 283);
    const auto *li_286 = buffer.data(li + 286);
    const auto *li_290 = buffer.data(li + 290);
    const auto *li_295 = buffer.data(li + 295);
    const auto *li_301 = buffer.data(li + 301);
    const auto *li_303 = buffer.data(li + 303);
    const auto *li_304 = buffer.data(li + 304);
    const auto *li_305 = buffer.data(li + 305);
    const auto *li_306 = buffer.data(li + 306);
    const auto *li_307 = buffer.data(li + 307);

    const auto *lk_136 = buffer.data(lk + 136);
    const auto *lk_138 = buffer.data(lk + 138);
    const auto *lk_139 = buffer.data(lk + 139);
    const auto *lk_140 = buffer.data(lk + 140);
    const auto *lk_141 = buffer.data(lk + 141);
    const auto *lk_143 = buffer.data(lk + 143);
    const auto *lk_180 = buffer.data(lk + 180);
    const auto *lk_182 = buffer.data(lk + 182);
    const auto *lk_183 = buffer.data(lk + 183);
    const auto *lk_185 = buffer.data(lk + 185);
    const auto *lk_186 = buffer.data(lk + 186);
    const auto *lk_189 = buffer.data(lk + 189);
    const auto *lk_190 = buffer.data(lk + 190);
    const auto *lk_192 = buffer.data(lk + 192);
    const auto *lk_194 = buffer.data(lk + 194);
    const auto *lk_195 = buffer.data(lk + 195);
    const auto *lk_197 = buffer.data(lk + 197);
    const auto *lk_198 = buffer.data(lk + 198);
    const auto *lk_200 = buffer.data(lk + 200);
    const auto *lk_207 = buffer.data(lk + 207);
    const auto *lk_208 = buffer.data(lk + 208);
    const auto *lk_210 = buffer.data(lk + 210);
    const auto *lk_211 = buffer.data(lk + 211);
    const auto *lk_212 = buffer.data(lk + 212);
    const auto *lk_213 = buffer.data(lk + 213);
    const auto *lk_215 = buffer.data(lk + 215);
    const auto *lk_216 = buffer.data(lk + 216);
    const auto *lk_217 = buffer.data(lk + 217);
    const auto *lk_219 = buffer.data(lk + 219);
    const auto *lk_221 = buffer.data(lk + 221);
    const auto *lk_222 = buffer.data(lk + 222);
    const auto *lk_225 = buffer.data(lk + 225);
    const auto *lk_226 = buffer.data(lk + 226);
    const auto *lk_359 = buffer.data(lk + 359);
    const auto *lk_388 = buffer.data(lk + 388);

    const auto *mh0_189 = buffer.data(mh0 + 189);
    const auto *mh0_190 = buffer.data(mh0 + 190);
    const auto *mh0_192 = buffer.data(mh0 + 192);
    const auto *mh0_194 = buffer.data(mh0 + 194);
    const auto *mh0_195 = buffer.data(mh0 + 195);
    const auto *mh0_197 = buffer.data(mh0 + 197);
    const auto *mh0_198 = buffer.data(mh0 + 198);
    const auto *mh0_203 = buffer.data(mh0 + 203);
    const auto *mh0_204 = buffer.data(mh0 + 204);
    const auto *mh0_206 = buffer.data(mh0 + 206);
    const auto *mh0_207 = buffer.data(mh0 + 207);
    const auto *mh0_208 = buffer.data(mh0 + 208);
    const auto *mh0_209 = buffer.data(mh0 + 209);
    const auto *mh0_210 = buffer.data(mh0 + 210);
    const auto *mh0_212 = buffer.data(mh0 + 212);
    const auto *mh0_213 = buffer.data(mh0 + 213);
    const auto *mh0_215 = buffer.data(mh0 + 215);
    const auto *mh0_216 = buffer.data(mh0 + 216);
    const auto *mh0_217 = buffer.data(mh0 + 217);
    const auto *mh0_219 = buffer.data(mh0 + 219);
    const auto *mh0_220 = buffer.data(mh0 + 220);
    const auto *mh0_225 = buffer.data(mh0 + 225);
    const auto *mh0_226 = buffer.data(mh0 + 226);
    const auto *mh0_227 = buffer.data(mh0 + 227);
    const auto *mh0_228 = buffer.data(mh0 + 228);
    const auto *mh0_230 = buffer.data(mh0 + 230);

    const auto *mh1_189 = buffer.data(mh1 + 189);
    const auto *mh1_190 = buffer.data(mh1 + 190);
    const auto *mh1_192 = buffer.data(mh1 + 192);
    const auto *mh1_194 = buffer.data(mh1 + 194);
    const auto *mh1_195 = buffer.data(mh1 + 195);
    const auto *mh1_197 = buffer.data(mh1 + 197);
    const auto *mh1_198 = buffer.data(mh1 + 198);
    const auto *mh1_203 = buffer.data(mh1 + 203);
    const auto *mh1_204 = buffer.data(mh1 + 204);
    const auto *mh1_206 = buffer.data(mh1 + 206);
    const auto *mh1_207 = buffer.data(mh1 + 207);
    const auto *mh1_208 = buffer.data(mh1 + 208);
    const auto *mh1_209 = buffer.data(mh1 + 209);
    const auto *mh1_210 = buffer.data(mh1 + 210);
    const auto *mh1_212 = buffer.data(mh1 + 212);
    const auto *mh1_213 = buffer.data(mh1 + 213);
    const auto *mh1_215 = buffer.data(mh1 + 215);
    const auto *mh1_216 = buffer.data(mh1 + 216);
    const auto *mh1_217 = buffer.data(mh1 + 217);
    const auto *mh1_219 = buffer.data(mh1 + 219);
    const auto *mh1_220 = buffer.data(mh1 + 220);
    const auto *mh1_225 = buffer.data(mh1 + 225);
    const auto *mh1_226 = buffer.data(mh1 + 226);
    const auto *mh1_227 = buffer.data(mh1 + 227);
    const auto *mh1_228 = buffer.data(mh1 + 228);
    const auto *mh1_230 = buffer.data(mh1 + 230);

    const auto *mi_217 = buffer.data(mi + 217);
    const auto *mi_219 = buffer.data(mi + 219);
    const auto *mi_220 = buffer.data(mi + 220);
    const auto *mi_221 = buffer.data(mi + 221);
    const auto *mi_222 = buffer.data(mi + 222);
    const auto *mi_223 = buffer.data(mi + 223);
    const auto *mi_224 = buffer.data(mi + 224);
    const auto *mi_226 = buffer.data(mi + 226);
    const auto *mi_227 = buffer.data(mi + 227);
    const auto *mi_229 = buffer.data(mi + 229);
    const auto *mi_230 = buffer.data(mi + 230);
    const auto *mi_233 = buffer.data(mi + 233);
    const auto *mi_234 = buffer.data(mi + 234);
    const auto *mi_238 = buffer.data(mi + 238);
    const auto *mi_245 = buffer.data(mi + 245);
    const auto *mi_246 = buffer.data(mi + 246);
    const auto *mi_247 = buffer.data(mi + 247);
    const auto *mi_248 = buffer.data(mi + 248);
    const auto *mi_249 = buffer.data(mi + 249);
    const auto *mi_250 = buffer.data(mi + 250);
    const auto *mi_251 = buffer.data(mi + 251);
    const auto *mi_252 = buffer.data(mi + 252);
    const auto *mi_253 = buffer.data(mi + 253);
    const auto *mi_254 = buffer.data(mi + 254);
    const auto *mi_255 = buffer.data(mi + 255);
    const auto *mi_257 = buffer.data(mi + 257);
    const auto *mi_258 = buffer.data(mi + 258);
    const auto *mi_260 = buffer.data(mi + 260);
    const auto *mi_261 = buffer.data(mi + 261);
    const auto *mi_262 = buffer.data(mi + 262);
    const auto *mi_264 = buffer.data(mi + 264);
    const auto *mi_265 = buffer.data(mi + 265);
    const auto *mi_266 = buffer.data(mi + 266);
    const auto *mi_272 = buffer.data(mi + 272);
    const auto *mi_273 = buffer.data(mi + 273);
    const auto *mi_274 = buffer.data(mi + 274);
    const auto *mi_275 = buffer.data(mi + 275);
    const auto *mi_276 = buffer.data(mi + 276);
    const auto *mi_277 = buffer.data(mi + 277);
    const auto *mi_278 = buffer.data(mi + 278);
    const auto *mi_279 = buffer.data(mi + 279);
    const auto *mi_280 = buffer.data(mi + 280);
    const auto *mi_281 = buffer.data(mi + 281);
    const auto *mi_282 = buffer.data(mi + 282);
    const auto *mi_283 = buffer.data(mi + 283);
    const auto *mi_285 = buffer.data(mi + 285);
    const auto *mi_286 = buffer.data(mi + 286);
    const auto *mi_287 = buffer.data(mi + 287);
    const auto *mi_289 = buffer.data(mi + 289);
    const auto *mi_290 = buffer.data(mi + 290);
    const auto *mi_291 = buffer.data(mi + 291);
    const auto *mi_292 = buffer.data(mi + 292);
    const auto *mi_294 = buffer.data(mi + 294);
    const auto *mi_295 = buffer.data(mi + 295);
    const auto *mi_301 = buffer.data(mi + 301);
    const auto *mi_302 = buffer.data(mi + 302);
    const auto *mi_303 = buffer.data(mi + 303);
    const auto *mi_304 = buffer.data(mi + 304);
    const auto *mi_305 = buffer.data(mi + 305);
    const auto *mi_306 = buffer.data(mi + 306);
    const auto *mi_307 = buffer.data(mi + 307);
    const auto *mi_308 = buffer.data(mi + 308);
    const auto *mi_310 = buffer.data(mi + 310);
    const auto *mi_311 = buffer.data(mi + 311);
    const auto *mi_313 = buffer.data(mi + 313);

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, pb_x, li_219, li_220, li_221, \
                         li_222, li_223, mi_219, mi_220, mi_221, mi_222, \
                         mi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_24 * li_219[k]
                   + pb_x[k] * mi_219[k];

        t_276[k] = f_24 * li_220[k]
                   + pb_x[k] * mi_220[k];

        t_277[k] = f_24 * li_221[k]
                   + pb_x[k] * mi_221[k];

        t_278[k] = f_24 * li_222[k]
                   + pb_x[k] * mi_222[k];

        t_279[k] = f_24 * li_223[k]
                   + pb_x[k] * mi_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, pa_z, pb_z, li_105, li_106, \
                         li_107, li_108, lk_136, lk_138, lk_139, lk_140, \
                         mi_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pa_z[k] * lk_136[k];

        t_281[k] = f_11 * li_105[k]
                   + pb_z[k] * mi_217[k];

        t_282[k] = f_12 * li_106[k]
                   + pa_z[k] * lk_138[k];

        t_283[k] = f_13 * li_107[k]
                   + pa_z[k] * lk_139[k];

        t_284[k] = f_14 * li_108[k]
                   + pa_z[k] * lk_140[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, pa_y, pa_z, pb_y, li_109, li_111, li_139, \
                         lk_141, lk_143, lk_180, mi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_15 * li_109[k]
                   + pa_z[k] * lk_141[k];

        t_286[k] = f_12 * li_139[k]
                   + pb_y[k] * mi_223[k];

        t_287[k] = f_17 * li_111[k]
                   + pa_z[k] * lk_143[k];

        t_288[k] = pa_y[k] * lk_180[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, t_293, pa_y, pb_y, li_140, li_141, \
                         li_142, lk_182, lk_183, lk_185, mi_224, \
                         mi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_11 * li_140[k]
                   + pb_y[k] * mi_224[k];

        t_290[k] = pa_y[k] * lk_182[k];

        t_291[k] = f_12 * li_141[k]
                   + pa_y[k] * lk_183[k];

        t_292[k] = f_11 * li_142[k]
                   + pb_y[k] * mi_226[k];

        t_293[k] = pa_y[k] * lk_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_y, pb_y, pb_z, li_115, li_143, li_145, \
                         lk_186, lk_189, mi_227, mi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_13 * li_143[k]
                   + pa_y[k] * lk_186[k];

        t_295[k] = f_12 * li_115[k]
                   + pb_z[k] * mi_227[k];

        t_296[k] = f_11 * li_145[k]
                   + pb_y[k] * mi_229[k];

        t_297[k] = pa_y[k] * lk_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pa_y, pb_y, pb_z, li_118, li_146, li_148, \
                         li_149, lk_190, lk_192, mi_230, mi_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_14 * li_146[k]
                   + pa_y[k] * lk_190[k];

        t_299[k] = f_12 * li_118[k]
                   + pb_z[k] * mi_230[k];

        t_300[k] = f_12 * li_148[k]
                   + pa_y[k] * lk_192[k];

        t_301[k] = f_11 * li_149[k]
                   + pb_y[k] * mi_233[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, pa_y, pb_z, li_122, li_150, \
                         li_152, li_153, lk_194, lk_195, lk_197, lk_198, \
                         mi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = pa_y[k] * lk_194[k];

        t_303[k] = f_15 * li_150[k]
                   + pa_y[k] * lk_195[k];

        t_304[k] = f_12 * li_122[k]
                   + pb_z[k] * mi_234[k];

        t_305[k] = f_13 * li_152[k]
                   + pa_y[k] * lk_197[k];

        t_306[k] = f_12 * li_153[k]
                   + pa_y[k] * lk_198[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pa_y, pb_x, pb_y, li_154, li_245, li_246, \
                         lk_200, mi_238, mi_245, mi_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_11 * li_154[k]
                   + pb_y[k] * mi_238[k];

        t_308[k] = pa_y[k] * lk_200[k];

        t_309[k] = f_24 * li_245[k]
                   + pb_x[k] * mi_245[k];

        t_310[k] = f_24 * li_246[k]
                   + pb_x[k] * mi_246[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, t_315, pa_y, pb_x, li_247, li_248, \
                         li_249, li_250, lk_207, mi_247, mi_248, mi_249, \
                         mi_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_24 * li_247[k]
                   + pb_x[k] * mi_247[k];

        t_312[k] = f_24 * li_248[k]
                   + pb_x[k] * mi_248[k];

        t_313[k] = f_24 * li_249[k]
                   + pb_x[k] * mi_249[k];

        t_314[k] = f_24 * li_250[k]
                   + pb_x[k] * mi_250[k];

        t_315[k] = pa_y[k] * lk_207[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pa_y, pb_z, li_133, li_161, li_163, \
                         li_164, lk_208, lk_210, lk_211, mi_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_17 * li_161[k]
                   + pa_y[k] * lk_208[k];

        t_317[k] = f_12 * li_133[k]
                   + pb_z[k] * mi_245[k];

        t_318[k] = f_15 * li_163[k]
                   + pa_y[k] * lk_210[k];

        t_319[k] = f_14 * li_164[k]
                   + pa_y[k] * lk_211[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_y, pb_y, li_165, li_166, li_167, \
                         lk_212, lk_213, lk_215, mi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_13 * li_165[k]
                   + pa_y[k] * lk_212[k];

        t_321[k] = f_12 * li_166[k]
                   + pa_y[k] * lk_213[k];

        t_322[k] = f_11 * li_167[k]
                   + pb_y[k] * mi_251[k];

        t_323[k] = pa_y[k] * lk_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pa_z, pb_y, pb_z, kk0_72, kk1_72, li_140, \
                         lk_180, mh0_189, mh1_189, mi_252, mi_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_22 * kk0_72[k]
                   - f_23 * kk1_72[k]
                   + pa_z[k] * lk_180[k];

        t_325[k] = pb_y[k] * mi_252[k];

        t_326[k] = f_13 * li_140[k]
                   + pb_z[k] * mi_252[k];

        t_327[k] = f_3 * mh0_189[k]
                   - f_4 * mh1_189[k]
                   + pb_y[k] * mi_253[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pb_x, pb_y, pb_z, li_143, li_257, \
                         mh0_190, mh0_194, mh1_190, mh1_194, mi_254, mi_255, \
                         mi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = pb_y[k] * mi_254[k];

        t_329[k] = f_24 * li_257[k]
                   + f_9 * mh0_194[k]
                   - f_10 * mh1_194[k]
                   + pb_x[k] * mi_257[k];

        t_330[k] = f_5 * mh0_190[k]
                   - f_6 * mh1_190[k]
                   + pb_y[k] * mi_255[k];

        t_331[k] = f_13 * li_143[k]
                   + pb_z[k] * mi_255[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pb_x, pb_y, pb_z, li_146, li_261, \
                         mh0_192, mh0_198, mh1_192, mh1_198, mi_257, mi_258, \
                         mi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pb_y[k] * mi_257[k];

        t_333[k] = f_24 * li_261[k]
                   + f_7 * mh0_198[k]
                   - f_8 * mh1_198[k]
                   + pb_x[k] * mi_261[k];

        t_334[k] = f_7 * mh0_192[k]
                   - f_8 * mh1_192[k]
                   + pb_y[k] * mi_258[k];

        t_335[k] = f_13 * li_146[k]
                   + pb_z[k] * mi_258[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pb_x, pb_y, li_266, mh0_194, mh0_203, mh1_194, \
                         mh1_203, mi_260, mi_261, mi_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_3 * mh0_194[k]
                   - f_4 * mh1_194[k]
                   + pb_y[k] * mi_260[k];

        t_337[k] = pb_y[k] * mi_261[k];

        t_338[k] = f_24 * li_266[k]
                   + f_5 * mh0_203[k]
                   - f_6 * mh1_203[k]
                   + pb_x[k] * mi_266[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, pb_y, pb_z, li_150, mh0_195, mh0_197, \
                         mh0_198, mh1_195, mh1_197, mh1_198, mi_262, mi_264, \
                         mi_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_9 * mh0_195[k]
                   - f_10 * mh1_195[k]
                   + pb_y[k] * mi_262[k];

        t_340[k] = f_13 * li_150[k]
                   + pb_z[k] * mi_262[k];

        t_341[k] = f_5 * mh0_197[k]
                   - f_6 * mh1_197[k]
                   + pb_y[k] * mi_264[k];

        t_342[k] = f_3 * mh0_198[k]
                   - f_4 * mh1_198[k]
                   + pb_y[k] * mi_265[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, pb_x, pb_y, li_272, li_273, li_274, \
                         mh0_209, mh1_209, mi_266, mi_272, mi_273, \
                         mi_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = pb_y[k] * mi_266[k];

        t_344[k] = f_24 * li_272[k]
                   + f_3 * mh0_209[k]
                   - f_4 * mh1_209[k]
                   + pb_x[k] * mi_272[k];

        t_345[k] = f_24 * li_273[k]
                   + pb_x[k] * mi_273[k];

        t_346[k] = f_24 * li_274[k]
                   + pb_x[k] * mi_274[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, pb_x, pb_y, li_275, li_276, \
                         li_277, li_279, mi_272, mi_275, mi_276, mi_277, \
                         mi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_24 * li_275[k]
                   + pb_x[k] * mi_275[k];

        t_348[k] = f_24 * li_276[k]
                   + pb_x[k] * mi_276[k];

        t_349[k] = f_24 * li_277[k]
                   + pb_x[k] * mi_277[k];

        t_350[k] = pb_y[k] * mi_272[k];

        t_351[k] = f_24 * li_279[k]
                   + pb_x[k] * mi_279[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pb_y, pb_z, li_161, mh0_204, mh0_206, \
                         mh0_207, mh1_204, mh1_206, mh1_207, mi_273, mi_275, \
                         mi_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_1 * mh0_204[k]
                   - f_2 * mh1_204[k]
                   + pb_y[k] * mi_273[k];

        t_353[k] = f_13 * li_161[k]
                   + pb_z[k] * mi_273[k];

        t_354[k] = f_9 * mh0_206[k]
                   - f_10 * mh1_206[k]
                   + pb_y[k] * mi_275[k];

        t_355[k] = f_7 * mh0_207[k]
                   - f_8 * mh1_207[k]
                   + pb_y[k] * mi_276[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pa_x, pb_y, kk0_359, kk1_359, lk_359, \
                         mh0_208, mh0_209, mh1_208, mh1_209, mi_277, mi_278, \
                         mi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_5 * mh0_208[k]
                   - f_6 * mh1_208[k]
                   + pb_y[k] * mi_277[k];

        t_357[k] = f_3 * mh0_209[k]
                   - f_4 * mh1_209[k]
                   + pb_y[k] * mi_278[k];

        t_358[k] = pb_y[k] * mi_279[k];

        t_359[k] = f_25 * kk0_359[k]
                   - f_26 * kk1_359[k]
                   + pa_x[k] * lk_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pa_y, pb_y, pb_z, kk0_108, kk1_108, li_168, \
                         lk_216, mi_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_27 * kk0_108[k]
                   - f_28 * kk1_108[k]
                   + pa_y[k] * lk_216[k];

        t_361[k] = f_14 * li_168[k]
                   + pb_y[k] * mi_280[k];

        t_362[k] = pb_z[k] * mi_280[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pb_x, pb_z, li_283, mh0_210, mh0_213, mh1_210, \
                         mh1_213, mi_281, mi_282, mi_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_15 * li_283[k]
                   + f_9 * mh0_213[k]
                   - f_10 * mh1_213[k]
                   + pb_x[k] * mi_283[k];

        t_364[k] = pb_z[k] * mi_281[k];

        t_365[k] = f_3 * mh0_210[k]
                   - f_4 * mh1_210[k]
                   + pb_z[k] * mi_282[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_x, pb_y, pb_z, li_173, li_286, \
                         mh0_212, mh0_216, mh1_212, mh1_216, mi_283, mi_285, \
                         mi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_15 * li_286[k]
                   + f_7 * mh0_216[k]
                   - f_8 * mh1_216[k]
                   + pb_x[k] * mi_286[k];

        t_367[k] = pb_z[k] * mi_283[k];

        t_368[k] = f_14 * li_173[k]
                   + pb_y[k] * mi_285[k];

        t_369[k] = f_5 * mh0_212[k]
                   - f_6 * mh1_212[k]
                   + pb_z[k] * mi_285[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pb_x, pb_z, li_290, mh0_213, mh0_220, mh1_213, \
                         mh1_220, mi_286, mi_287, mi_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_15 * li_290[k]
                   + f_5 * mh0_220[k]
                   - f_6 * mh1_220[k]
                   + pb_x[k] * mi_290[k];

        t_371[k] = pb_z[k] * mi_286[k];

        t_372[k] = f_3 * mh0_213[k]
                   - f_4 * mh1_213[k]
                   + pb_z[k] * mi_287[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pb_x, pb_y, pb_z, li_177, li_295, \
                         mh0_215, mh0_225, mh1_215, mh1_225, mi_289, mi_290, \
                         mi_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_14 * li_177[k]
                   + pb_y[k] * mi_289[k];

        t_374[k] = f_7 * mh0_215[k]
                   - f_8 * mh1_215[k]
                   + pb_z[k] * mi_289[k];

        t_375[k] = f_15 * li_295[k]
                   + f_3 * mh0_225[k]
                   - f_4 * mh1_225[k]
                   + pb_x[k] * mi_295[k];

        t_376[k] = pb_z[k] * mi_290[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pb_y, pb_z, li_182, mh0_216, mh0_217, \
                         mh0_219, mh1_216, mh1_217, mh1_219, mi_291, mi_292, \
                         mi_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_3 * mh0_216[k]
                   - f_4 * mh1_216[k]
                   + pb_z[k] * mi_291[k];

        t_378[k] = f_5 * mh0_217[k]
                   - f_6 * mh1_217[k]
                   + pb_z[k] * mi_292[k];

        t_379[k] = f_14 * li_182[k]
                   + pb_y[k] * mi_294[k];

        t_380[k] = f_9 * mh0_219[k]
                   - f_10 * mh1_219[k]
                   + pb_z[k] * mi_294[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, t_385, pb_x, pb_z, li_301, li_303, \
                         li_304, li_305, mi_295, mi_301, mi_303, mi_304, \
                         mi_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_15 * li_301[k]
                   + pb_x[k] * mi_301[k];

        t_382[k] = pb_z[k] * mi_295[k];

        t_383[k] = f_15 * li_303[k]
                   + pb_x[k] * mi_303[k];

        t_384[k] = f_15 * li_304[k]
                   + pb_x[k] * mi_304[k];

        t_385[k] = f_15 * li_305[k]
                   + pb_x[k] * mi_305[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pa_x, pb_x, pb_z, kk0_388, kk1_388, \
                         li_306, li_307, lk_388, mi_301, mi_306, \
                         mi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_15 * li_306[k]
                   + pb_x[k] * mi_306[k];

        t_387[k] = f_15 * li_307[k]
                   + pb_x[k] * mi_307[k];

        t_388[k] = f_29 * kk0_388[k]
                   - f_30 * kk1_388[k]
                   + pa_x[k] * lk_388[k];

        t_389[k] = pb_z[k] * mi_301[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, pb_z, mh0_225, mh0_226, mh0_227, mh1_225, \
                         mh1_226, mh1_227, mi_302, mi_303, mi_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_3 * mh0_225[k]
                   - f_4 * mh1_225[k]
                   + pb_z[k] * mi_302[k];

        t_391[k] = f_5 * mh0_226[k]
                   - f_6 * mh1_226[k]
                   + pb_z[k] * mi_303[k];

        t_392[k] = f_7 * mh0_227[k]
                   - f_8 * mh1_227[k]
                   + pb_z[k] * mi_304[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, pa_z, pb_y, pb_z, li_195, lk_216, \
                         mh0_228, mh0_230, mh1_228, mh1_230, mi_305, \
                         mi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = f_9 * mh0_228[k]
                   - f_10 * mh1_228[k]
                   + pb_z[k] * mi_305[k];

        t_394[k] = f_14 * li_195[k]
                   + pb_y[k] * mi_307[k];

        t_395[k] = f_1 * mh0_230[k]
                   - f_2 * mh1_230[k]
                   + pb_z[k] * mi_307[k];

        t_396[k] = pa_z[k] * lk_216[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, pa_z, pb_y, pb_z, li_168, li_170, \
                         li_198, lk_217, lk_219, lk_221, mi_308, \
                         mi_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = pa_z[k] * lk_217[k];

        t_398[k] = f_11 * li_168[k]
                   + pb_z[k] * mi_308[k];

        t_399[k] = pa_z[k] * lk_219[k];

        t_400[k] = f_13 * li_198[k]
                   + pb_y[k] * mi_310[k];

        t_401[k] = f_12 * li_170[k]
                   + pa_z[k] * lk_221[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, pa_z, pb_y, pb_z, li_171, li_173, \
                         li_201, lk_222, lk_225, lk_226, mi_311, \
                         mi_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pa_z[k] * lk_222[k];

        t_403[k] = f_11 * li_171[k]
                   + pb_z[k] * mi_311[k];

        t_404[k] = f_13 * li_201[k]
                   + pb_y[k] * mi_313[k];

        t_405[k] = f_13 * li_173[k]
                   + pa_z[k] * lk_225[k];

        t_406[k] = pa_z[k] * lk_226[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kk0,
                                            const size_t kk1, const size_t li, const size_t lk,
                                            const size_t mh0, const size_t mh1, const size_t mi,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_17 = 3.5 / p;
    const auto f_18 = 0.5 / alpha;
    const auto f_19 = 0.5 * beta / (alpha * p);
    const auto f_27 = 1.5 / alpha;
    const auto f_28 = 1.5 * beta / (alpha * p);
    const auto f_29 = 2.0 / alpha;
    const auto f_30 = 2.0 * beta / (alpha * p);

    auto *t_407 = buffer.data(target + 407);
    auto *t_408 = buffer.data(target + 408);
    auto *t_409 = buffer.data(target + 409);
    auto *t_410 = buffer.data(target + 410);
    auto *t_411 = buffer.data(target + 411);
    auto *t_412 = buffer.data(target + 412);
    auto *t_413 = buffer.data(target + 413);
    auto *t_414 = buffer.data(target + 414);
    auto *t_415 = buffer.data(target + 415);
    auto *t_416 = buffer.data(target + 416);
    auto *t_417 = buffer.data(target + 417);
    auto *t_418 = buffer.data(target + 418);
    auto *t_419 = buffer.data(target + 419);
    auto *t_420 = buffer.data(target + 420);
    auto *t_421 = buffer.data(target + 421);
    auto *t_422 = buffer.data(target + 422);
    auto *t_423 = buffer.data(target + 423);
    auto *t_424 = buffer.data(target + 424);
    auto *t_425 = buffer.data(target + 425);
    auto *t_426 = buffer.data(target + 426);
    auto *t_427 = buffer.data(target + 427);
    auto *t_428 = buffer.data(target + 428);
    auto *t_429 = buffer.data(target + 429);
    auto *t_430 = buffer.data(target + 430);
    auto *t_431 = buffer.data(target + 431);
    auto *t_432 = buffer.data(target + 432);
    auto *t_433 = buffer.data(target + 433);
    auto *t_434 = buffer.data(target + 434);
    auto *t_435 = buffer.data(target + 435);
    auto *t_436 = buffer.data(target + 436);
    auto *t_437 = buffer.data(target + 437);
    auto *t_438 = buffer.data(target + 438);
    auto *t_439 = buffer.data(target + 439);
    auto *t_440 = buffer.data(target + 440);
    auto *t_441 = buffer.data(target + 441);
    auto *t_442 = buffer.data(target + 442);
    auto *t_443 = buffer.data(target + 443);
    auto *t_444 = buffer.data(target + 444);
    auto *t_445 = buffer.data(target + 445);
    auto *t_446 = buffer.data(target + 446);
    auto *t_447 = buffer.data(target + 447);
    auto *t_448 = buffer.data(target + 448);
    auto *t_449 = buffer.data(target + 449);
    auto *t_450 = buffer.data(target + 450);
    auto *t_451 = buffer.data(target + 451);
    auto *t_452 = buffer.data(target + 452);
    auto *t_453 = buffer.data(target + 453);
    auto *t_454 = buffer.data(target + 454);
    auto *t_455 = buffer.data(target + 455);
    auto *t_456 = buffer.data(target + 456);
    auto *t_457 = buffer.data(target + 457);
    auto *t_458 = buffer.data(target + 458);
    auto *t_459 = buffer.data(target + 459);
    auto *t_460 = buffer.data(target + 460);
    auto *t_461 = buffer.data(target + 461);
    auto *t_462 = buffer.data(target + 462);
    auto *t_463 = buffer.data(target + 463);
    auto *t_464 = buffer.data(target + 464);
    auto *t_465 = buffer.data(target + 465);
    auto *t_466 = buffer.data(target + 466);
    auto *t_467 = buffer.data(target + 467);
    auto *t_468 = buffer.data(target + 468);
    auto *t_469 = buffer.data(target + 469);
    auto *t_470 = buffer.data(target + 470);
    auto *t_471 = buffer.data(target + 471);
    auto *t_472 = buffer.data(target + 472);
    auto *t_473 = buffer.data(target + 473);
    auto *t_474 = buffer.data(target + 474);
    auto *t_475 = buffer.data(target + 475);
    auto *t_476 = buffer.data(target + 476);
    auto *t_477 = buffer.data(target + 477);
    auto *t_478 = buffer.data(target + 478);
    auto *t_479 = buffer.data(target + 479);
    auto *t_480 = buffer.data(target + 480);
    auto *t_481 = buffer.data(target + 481);
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);
    auto *t_495 = buffer.data(target + 495);
    auto *t_496 = buffer.data(target + 496);
    auto *t_497 = buffer.data(target + 497);
    auto *t_498 = buffer.data(target + 498);
    auto *t_499 = buffer.data(target + 499);
    auto *t_500 = buffer.data(target + 500);
    auto *t_501 = buffer.data(target + 501);
    auto *t_502 = buffer.data(target + 502);
    auto *t_503 = buffer.data(target + 503);
    auto *t_504 = buffer.data(target + 504);
    auto *t_505 = buffer.data(target + 505);
    auto *t_506 = buffer.data(target + 506);
    auto *t_507 = buffer.data(target + 507);
    auto *t_508 = buffer.data(target + 508);
    auto *t_509 = buffer.data(target + 509);
    auto *t_510 = buffer.data(target + 510);
    auto *t_511 = buffer.data(target + 511);
    auto *t_512 = buffer.data(target + 512);
    auto *t_513 = buffer.data(target + 513);
    auto *t_514 = buffer.data(target + 514);
    auto *t_515 = buffer.data(target + 515);
    auto *t_516 = buffer.data(target + 516);
    auto *t_517 = buffer.data(target + 517);
    auto *t_518 = buffer.data(target + 518);
    auto *t_519 = buffer.data(target + 519);
    auto *t_520 = buffer.data(target + 520);
    auto *t_521 = buffer.data(target + 521);
    auto *t_522 = buffer.data(target + 522);
    auto *t_523 = buffer.data(target + 523);
    auto *t_524 = buffer.data(target + 524);
    auto *t_525 = buffer.data(target + 525);
    auto *t_526 = buffer.data(target + 526);
    auto *t_527 = buffer.data(target + 527);
    auto *t_528 = buffer.data(target + 528);
    auto *t_529 = buffer.data(target + 529);
    auto *t_530 = buffer.data(target + 530);
    auto *t_531 = buffer.data(target + 531);
    auto *t_532 = buffer.data(target + 532);
    auto *t_533 = buffer.data(target + 533);
    auto *t_534 = buffer.data(target + 534);
    auto *t_535 = buffer.data(target + 535);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kk0_111 = buffer.data(kk0 + 111);
    const auto *kk0_114 = buffer.data(kk0 + 114);
    const auto *kk0_118 = buffer.data(kk0 + 118);
    const auto *kk0_123 = buffer.data(kk0 + 123);
    const auto *kk0_180 = buffer.data(kk0 + 180);
    const auto *kk0_185 = buffer.data(kk0 + 185);
    const auto *kk0_189 = buffer.data(kk0 + 189);
    const auto *kk0_194 = buffer.data(kk0 + 194);
    const auto *kk0_200 = buffer.data(kk0 + 200);
    const auto *kk0_460 = buffer.data(kk0 + 460);
    const auto *kk0_462 = buffer.data(kk0 + 462);
    const auto *kk0_463 = buffer.data(kk0 + 463);
    const auto *kk0_464 = buffer.data(kk0 + 464);
    const auto *kk0_465 = buffer.data(kk0 + 465);
    const auto *kk0_467 = buffer.data(kk0 + 467);

    const auto *kk1_111 = buffer.data(kk1 + 111);
    const auto *kk1_114 = buffer.data(kk1 + 114);
    const auto *kk1_118 = buffer.data(kk1 + 118);
    const auto *kk1_123 = buffer.data(kk1 + 123);
    const auto *kk1_180 = buffer.data(kk1 + 180);
    const auto *kk1_185 = buffer.data(kk1 + 185);
    const auto *kk1_189 = buffer.data(kk1 + 189);
    const auto *kk1_194 = buffer.data(kk1 + 194);
    const auto *kk1_200 = buffer.data(kk1 + 200);
    const auto *kk1_460 = buffer.data(kk1 + 460);
    const auto *kk1_462 = buffer.data(kk1 + 462);
    const auto *kk1_463 = buffer.data(kk1 + 463);
    const auto *kk1_464 = buffer.data(kk1 + 464);
    const auto *kk1_465 = buffer.data(kk1 + 465);
    const auto *kk1_467 = buffer.data(kk1 + 467);

    const auto *li_174 = buffer.data(li + 174);
    const auto *li_175 = buffer.data(li + 175);
    const auto *li_177 = buffer.data(li + 177);
    const auto *li_178 = buffer.data(li + 178);
    const auto *li_179 = buffer.data(li + 179);
    const auto *li_180 = buffer.data(li + 180);
    const auto *li_182 = buffer.data(li + 182);
    const auto *li_189 = buffer.data(li + 189);
    const auto *li_190 = buffer.data(li + 190);
    const auto *li_191 = buffer.data(li + 191);
    const auto *li_192 = buffer.data(li + 192);
    const auto *li_193 = buffer.data(li + 193);
    const auto *li_195 = buffer.data(li + 195);
    const auto *li_196 = buffer.data(li + 196);
    const auto *li_199 = buffer.data(li + 199);
    const auto *li_202 = buffer.data(li + 202);
    const auto *li_205 = buffer.data(li + 205);
    const auto *li_206 = buffer.data(li + 206);
    const auto *li_210 = buffer.data(li + 210);
    const auto *li_217 = buffer.data(li + 217);
    const auto *li_223 = buffer.data(li + 223);
    const auto *li_224 = buffer.data(li + 224);
    const auto *li_226 = buffer.data(li + 226);
    const auto *li_227 = buffer.data(li + 227);
    const auto *li_229 = buffer.data(li + 229);
    const auto *li_230 = buffer.data(li + 230);
    const auto *li_233 = buffer.data(li + 233);
    const auto *li_234 = buffer.data(li + 234);
    const auto *li_238 = buffer.data(li + 238);
    const auto *li_245 = buffer.data(li + 245);
    const auto *li_251 = buffer.data(li + 251);
    const auto *li_252 = buffer.data(li + 252);
    const auto *li_253 = buffer.data(li + 253);
    const auto *li_254 = buffer.data(li + 254);
    const auto *li_255 = buffer.data(li + 255);
    const auto *li_257 = buffer.data(li + 257);
    const auto *li_258 = buffer.data(li + 258);
    const auto *li_260 = buffer.data(li + 260);
    const auto *li_261 = buffer.data(li + 261);
    const auto *li_262 = buffer.data(li + 262);
    const auto *li_264 = buffer.data(li + 264);
    const auto *li_265 = buffer.data(li + 265);
    const auto *li_266 = buffer.data(li + 266);
    const auto *li_273 = buffer.data(li + 273);
    const auto *li_275 = buffer.data(li + 275);
    const auto *li_276 = buffer.data(li + 276);
    const auto *li_277 = buffer.data(li + 277);
    const auto *li_278 = buffer.data(li + 278);
    const auto *li_279 = buffer.data(li + 279);
    const auto *li_330 = buffer.data(li + 330);
    const auto *li_331 = buffer.data(li + 331);
    const auto *li_332 = buffer.data(li + 332);
    const auto *li_333 = buffer.data(li + 333);
    const auto *li_334 = buffer.data(li + 334);
    const auto *li_335 = buffer.data(li + 335);
    const auto *li_348 = buffer.data(li + 348);
    const auto *li_353 = buffer.data(li + 353);
    const auto *li_354 = buffer.data(li + 354);
    const auto *li_357 = buffer.data(li + 357);
    const auto *li_358 = buffer.data(li + 358);
    const auto *li_359 = buffer.data(li + 359);
    const auto *li_360 = buffer.data(li + 360);
    const auto *li_361 = buffer.data(li + 361);
    const auto *li_362 = buffer.data(li + 362);
    const auto *li_363 = buffer.data(li + 363);
    const auto *li_385 = buffer.data(li + 385);
    const auto *li_386 = buffer.data(li + 386);
    const auto *li_387 = buffer.data(li + 387);
    const auto *li_388 = buffer.data(li + 388);
    const auto *li_389 = buffer.data(li + 389);
    const auto *li_390 = buffer.data(li + 390);
    const auto *li_397 = buffer.data(li + 397);
    const auto *li_401 = buffer.data(li + 401);
    const auto *li_406 = buffer.data(li + 406);
    const auto *li_412 = buffer.data(li + 412);
    const auto *li_413 = buffer.data(li + 413);
    const auto *li_414 = buffer.data(li + 414);
    const auto *li_415 = buffer.data(li + 415);
    const auto *li_416 = buffer.data(li + 416);
    const auto *li_417 = buffer.data(li + 417);
    const auto *li_419 = buffer.data(li + 419);

    const auto *lk_228 = buffer.data(lk + 228);
    const auto *lk_230 = buffer.data(lk + 230);
    const auto *lk_231 = buffer.data(lk + 231);
    const auto *lk_233 = buffer.data(lk + 233);
    const auto *lk_234 = buffer.data(lk + 234);
    const auto *lk_236 = buffer.data(lk + 236);
    const auto *lk_237 = buffer.data(lk + 237);
    const auto *lk_244 = buffer.data(lk + 244);
    const auto *lk_246 = buffer.data(lk + 246);
    const auto *lk_247 = buffer.data(lk + 247);
    const auto *lk_248 = buffer.data(lk + 248);
    const auto *lk_249 = buffer.data(lk + 249);
    const auto *lk_251 = buffer.data(lk + 251);
    const auto *lk_255 = buffer.data(lk + 255);
    const auto *lk_258 = buffer.data(lk + 258);
    const auto *lk_262 = buffer.data(lk + 262);
    const auto *lk_267 = buffer.data(lk + 267);
    const auto *lk_288 = buffer.data(lk + 288);
    const auto *lk_293 = buffer.data(lk + 293);
    const auto *lk_297 = buffer.data(lk + 297);
    const auto *lk_302 = buffer.data(lk + 302);
    const auto *lk_308 = buffer.data(lk + 308);
    const auto *lk_324 = buffer.data(lk + 324);
    const auto *lk_326 = buffer.data(lk + 326);
    const auto *lk_327 = buffer.data(lk + 327);
    const auto *lk_329 = buffer.data(lk + 329);
    const auto *lk_330 = buffer.data(lk + 330);
    const auto *lk_333 = buffer.data(lk + 333);
    const auto *lk_334 = buffer.data(lk + 334);
    const auto *lk_336 = buffer.data(lk + 336);
    const auto *lk_338 = buffer.data(lk + 338);
    const auto *lk_339 = buffer.data(lk + 339);
    const auto *lk_341 = buffer.data(lk + 341);
    const auto *lk_342 = buffer.data(lk + 342);
    const auto *lk_344 = buffer.data(lk + 344);
    const auto *lk_351 = buffer.data(lk + 351);
    const auto *lk_352 = buffer.data(lk + 352);
    const auto *lk_354 = buffer.data(lk + 354);
    const auto *lk_355 = buffer.data(lk + 355);
    const auto *lk_356 = buffer.data(lk + 356);
    const auto *lk_357 = buffer.data(lk + 357);
    const auto *lk_359 = buffer.data(lk + 359);
    const auto *lk_460 = buffer.data(lk + 460);
    const auto *lk_462 = buffer.data(lk + 462);
    const auto *lk_463 = buffer.data(lk + 463);
    const auto *lk_464 = buffer.data(lk + 464);
    const auto *lk_465 = buffer.data(lk + 465);
    const auto *lk_467 = buffer.data(lk + 467);

    const auto *mh0_264 = buffer.data(mh0 + 264);
    const auto *mh0_269 = buffer.data(mh0 + 269);
    const auto *mh0_270 = buffer.data(mh0 + 270);
    const auto *mh0_294 = buffer.data(mh0 + 294);
    const auto *mh0_295 = buffer.data(mh0 + 295);
    const auto *mh0_297 = buffer.data(mh0 + 297);
    const auto *mh0_299 = buffer.data(mh0 + 299);
    const auto *mh0_300 = buffer.data(mh0 + 300);
    const auto *mh0_302 = buffer.data(mh0 + 302);
    const auto *mh0_303 = buffer.data(mh0 + 303);
    const auto *mh0_308 = buffer.data(mh0 + 308);
    const auto *mh0_309 = buffer.data(mh0 + 309);
    const auto *mh0_311 = buffer.data(mh0 + 311);
    const auto *mh0_312 = buffer.data(mh0 + 312);
    const auto *mh0_314 = buffer.data(mh0 + 314);

    const auto *mh1_264 = buffer.data(mh1 + 264);
    const auto *mh1_269 = buffer.data(mh1 + 269);
    const auto *mh1_270 = buffer.data(mh1 + 270);
    const auto *mh1_294 = buffer.data(mh1 + 294);
    const auto *mh1_295 = buffer.data(mh1 + 295);
    const auto *mh1_297 = buffer.data(mh1 + 297);
    const auto *mh1_299 = buffer.data(mh1 + 299);
    const auto *mh1_300 = buffer.data(mh1 + 300);
    const auto *mh1_302 = buffer.data(mh1 + 302);
    const auto *mh1_303 = buffer.data(mh1 + 303);
    const auto *mh1_308 = buffer.data(mh1 + 308);
    const auto *mh1_309 = buffer.data(mh1 + 309);
    const auto *mh1_311 = buffer.data(mh1 + 311);
    const auto *mh1_312 = buffer.data(mh1 + 312);
    const auto *mh1_314 = buffer.data(mh1 + 314);

    const auto *mi_314 = buffer.data(mi + 314);
    const auto *mi_317 = buffer.data(mi + 317);
    const auto *mi_318 = buffer.data(mi + 318);
    const auto *mi_322 = buffer.data(mi + 322);
    const auto *mi_329 = buffer.data(mi + 329);
    const auto *mi_330 = buffer.data(mi + 330);
    const auto *mi_331 = buffer.data(mi + 331);
    const auto *mi_332 = buffer.data(mi + 332);
    const auto *mi_333 = buffer.data(mi + 333);
    const auto *mi_334 = buffer.data(mi + 334);
    const auto *mi_335 = buffer.data(mi + 335);
    const auto *mi_336 = buffer.data(mi + 336);
    const auto *mi_338 = buffer.data(mi + 338);
    const auto *mi_339 = buffer.data(mi + 339);
    const auto *mi_341 = buffer.data(mi + 341);
    const auto *mi_342 = buffer.data(mi + 342);
    const auto *mi_345 = buffer.data(mi + 345);
    const auto *mi_346 = buffer.data(mi + 346);
    const auto *mi_348 = buffer.data(mi + 348);
    const auto *mi_350 = buffer.data(mi + 350);
    const auto *mi_353 = buffer.data(mi + 353);
    const auto *mi_354 = buffer.data(mi + 354);
    const auto *mi_357 = buffer.data(mi + 357);
    const auto *mi_358 = buffer.data(mi + 358);
    const auto *mi_359 = buffer.data(mi + 359);
    const auto *mi_360 = buffer.data(mi + 360);
    const auto *mi_361 = buffer.data(mi + 361);
    const auto *mi_362 = buffer.data(mi + 362);
    const auto *mi_363 = buffer.data(mi + 363);
    const auto *mi_364 = buffer.data(mi + 364);
    const auto *mi_366 = buffer.data(mi + 366);
    const auto *mi_367 = buffer.data(mi + 367);
    const auto *mi_369 = buffer.data(mi + 369);
    const auto *mi_370 = buffer.data(mi + 370);
    const auto *mi_373 = buffer.data(mi + 373);
    const auto *mi_374 = buffer.data(mi + 374);
    const auto *mi_378 = buffer.data(mi + 378);
    const auto *mi_385 = buffer.data(mi + 385);
    const auto *mi_386 = buffer.data(mi + 386);
    const auto *mi_387 = buffer.data(mi + 387);
    const auto *mi_388 = buffer.data(mi + 388);
    const auto *mi_389 = buffer.data(mi + 389);
    const auto *mi_390 = buffer.data(mi + 390);
    const auto *mi_391 = buffer.data(mi + 391);
    const auto *mi_392 = buffer.data(mi + 392);
    const auto *mi_393 = buffer.data(mi + 393);
    const auto *mi_394 = buffer.data(mi + 394);
    const auto *mi_395 = buffer.data(mi + 395);
    const auto *mi_397 = buffer.data(mi + 397);
    const auto *mi_398 = buffer.data(mi + 398);
    const auto *mi_400 = buffer.data(mi + 400);
    const auto *mi_401 = buffer.data(mi + 401);
    const auto *mi_402 = buffer.data(mi + 402);
    const auto *mi_404 = buffer.data(mi + 404);
    const auto *mi_405 = buffer.data(mi + 405);
    const auto *mi_406 = buffer.data(mi + 406);
    const auto *mi_412 = buffer.data(mi + 412);
    const auto *mi_413 = buffer.data(mi + 413);
    const auto *mi_414 = buffer.data(mi + 414);
    const auto *mi_415 = buffer.data(mi + 415);
    const auto *mi_416 = buffer.data(mi + 416);
    const auto *mi_417 = buffer.data(mi + 417);
    const auto *mi_419 = buffer.data(mi + 419);

#pragma omp simd aligned(t_407, t_408, t_409, t_410, pa_z, pb_y, pb_z, li_174, li_175, li_177, \
                         li_205, lk_228, lk_230, mi_314, mi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_11 * li_174[k]
                   + pb_z[k] * mi_314[k];

        t_408[k] = f_12 * li_175[k]
                   + pa_z[k] * lk_228[k];

        t_409[k] = f_13 * li_205[k]
                   + pb_y[k] * mi_317[k];

        t_410[k] = f_14 * li_177[k]
                   + pa_z[k] * lk_230[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, pa_z, pb_z, li_178, li_179, li_180, \
                         lk_231, lk_233, lk_234, mi_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = pa_z[k] * lk_231[k];

        t_412[k] = f_11 * li_178[k]
                   + pb_z[k] * mi_318[k];

        t_413[k] = f_12 * li_179[k]
                   + pa_z[k] * lk_233[k];

        t_414[k] = f_13 * li_180[k]
                   + pa_z[k] * lk_234[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, pa_z, pb_x, pb_y, li_182, li_210, li_330, \
                         lk_236, lk_237, mi_322, mi_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_13 * li_210[k]
                   + pb_y[k] * mi_322[k];

        t_416[k] = f_15 * li_182[k]
                   + pa_z[k] * lk_236[k];

        t_417[k] = pa_z[k] * lk_237[k];

        t_418[k] = f_15 * li_330[k]
                   + pb_x[k] * mi_330[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, t_423, pb_x, li_331, li_332, li_333, \
                         li_334, li_335, mi_331, mi_332, mi_333, mi_334, \
                         mi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_15 * li_331[k]
                   + pb_x[k] * mi_331[k];

        t_420[k] = f_15 * li_332[k]
                   + pb_x[k] * mi_332[k];

        t_421[k] = f_15 * li_333[k]
                   + pb_x[k] * mi_333[k];

        t_422[k] = f_15 * li_334[k]
                   + pb_x[k] * mi_334[k];

        t_423[k] = f_15 * li_335[k]
                   + pb_x[k] * mi_335[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, t_428, pa_z, pb_z, li_189, li_190, \
                         li_191, li_192, lk_244, lk_246, lk_247, lk_248, \
                         mi_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pa_z[k] * lk_244[k];

        t_425[k] = f_11 * li_189[k]
                   + pb_z[k] * mi_329[k];

        t_426[k] = f_12 * li_190[k]
                   + pa_z[k] * lk_246[k];

        t_427[k] = f_13 * li_191[k]
                   + pa_z[k] * lk_247[k];

        t_428[k] = f_14 * li_192[k]
                   + pa_z[k] * lk_248[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pa_y, pa_z, pb_y, kk0_180, kk1_180, \
                         li_193, li_195, li_223, lk_249, lk_251, lk_288, \
                         mi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_15 * li_193[k]
                   + pa_z[k] * lk_249[k];

        t_430[k] = f_13 * li_223[k]
                   + pb_y[k] * mi_335[k];

        t_431[k] = f_17 * li_195[k]
                   + pa_z[k] * lk_251[k];

        t_432[k] = f_18 * kk0_180[k]
                   - f_19 * kk1_180[k]
                   + pa_y[k] * lk_288[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pa_z, pb_y, pb_z, kk0_111, kk1_111, \
                         li_196, li_224, li_226, lk_255, mi_336, \
                         mi_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_12 * li_224[k]
                   + pb_y[k] * mi_336[k];

        t_434[k] = f_12 * li_196[k]
                   + pb_z[k] * mi_336[k];

        t_435[k] = f_18 * kk0_111[k]
                   - f_19 * kk1_111[k]
                   + pa_z[k] * lk_255[k];

        t_436[k] = f_12 * li_226[k]
                   + pb_y[k] * mi_338[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pa_y, pa_z, pb_z, kk0_114, kk0_185, kk1_114, \
                         kk1_185, li_199, lk_258, lk_293, mi_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_18 * kk0_185[k]
                   - f_19 * kk1_185[k]
                   + pa_y[k] * lk_293[k];

        t_438[k] = f_18 * kk0_114[k]
                   - f_19 * kk1_114[k]
                   + pa_z[k] * lk_258[k];

        t_439[k] = f_12 * li_199[k]
                   + pb_z[k] * mi_339[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pa_y, pa_z, pb_y, kk0_118, kk0_189, kk1_118, \
                         kk1_189, li_229, lk_262, lk_297, mi_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_12 * li_229[k]
                   + pb_y[k] * mi_341[k];

        t_441[k] = f_18 * kk0_189[k]
                   - f_19 * kk1_189[k]
                   + pa_y[k] * lk_297[k];

        t_442[k] = f_18 * kk0_118[k]
                   - f_19 * kk1_118[k]
                   + pa_z[k] * lk_262[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pb_x, pb_y, pb_z, li_202, li_233, li_348, \
                         mh0_264, mh1_264, mi_342, mi_345, mi_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_12 * li_202[k]
                   + pb_z[k] * mi_342[k];

        t_444[k] = f_15 * li_348[k]
                   + f_5 * mh0_264[k]
                   - f_6 * mh1_264[k]
                   + pb_x[k] * mi_348[k];

        t_445[k] = f_12 * li_233[k]
                   + pb_y[k] * mi_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pa_y, pa_z, pb_z, kk0_123, kk0_194, kk1_123, \
                         kk1_194, li_206, lk_267, lk_302, mi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_18 * kk0_194[k]
                   - f_19 * kk1_194[k]
                   + pa_y[k] * lk_302[k];

        t_447[k] = f_18 * kk0_123[k]
                   - f_19 * kk1_123[k]
                   + pa_z[k] * lk_267[k];

        t_448[k] = f_12 * li_206[k]
                   + pb_z[k] * mi_346[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pb_x, pb_y, li_238, li_353, li_354, mh0_269, \
                         mh0_270, mh1_269, mh1_270, mi_350, mi_353, \
                         mi_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_15 * li_353[k]
                   + f_3 * mh0_269[k]
                   - f_4 * mh1_269[k]
                   + pb_x[k] * mi_353[k];

        t_450[k] = f_15 * li_354[k]
                   + f_3 * mh0_270[k]
                   - f_4 * mh1_270[k]
                   + pb_x[k] * mi_354[k];

        t_451[k] = f_12 * li_238[k]
                   + pb_y[k] * mi_350[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pa_y, pb_x, kk0_200, kk1_200, li_357, \
                         li_358, li_359, lk_308, mi_357, mi_358, \
                         mi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_18 * kk0_200[k]
                   - f_19 * kk1_200[k]
                   + pa_y[k] * lk_308[k];

        t_453[k] = f_15 * li_357[k]
                   + pb_x[k] * mi_357[k];

        t_454[k] = f_15 * li_358[k]
                   + pb_x[k] * mi_358[k];

        t_455[k] = f_15 * li_359[k]
                   + pb_x[k] * mi_359[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pb_x, li_360, li_361, li_362, li_363, \
                         mi_360, mi_361, mi_362, mi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_15 * li_360[k]
                   + pb_x[k] * mi_360[k];

        t_457[k] = f_15 * li_361[k]
                   + pb_x[k] * mi_361[k];

        t_458[k] = f_15 * li_362[k]
                   + pb_x[k] * mi_362[k];

        t_459[k] = f_15 * li_363[k]
                   + pb_x[k] * mi_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pa_x, pb_z, kk0_460, kk0_462, kk1_460, kk1_462, \
                         li_217, lk_460, lk_462, mi_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_29 * kk0_460[k]
                   - f_30 * kk1_460[k]
                   + pa_x[k] * lk_460[k];

        t_461[k] = f_12 * li_217[k]
                   + pb_z[k] * mi_357[k];

        t_462[k] = f_29 * kk0_462[k]
                   - f_30 * kk1_462[k]
                   + pa_x[k] * lk_462[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pa_x, kk0_463, kk0_464, kk0_465, kk1_463, \
                         kk1_464, kk1_465, lk_463, lk_464, lk_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_29 * kk0_463[k]
                   - f_30 * kk1_463[k]
                   + pa_x[k] * lk_463[k];

        t_464[k] = f_29 * kk0_464[k]
                   - f_30 * kk1_464[k]
                   + pa_x[k] * lk_464[k];

        t_465[k] = f_29 * kk0_465[k]
                   - f_30 * kk1_465[k]
                   + pa_x[k] * lk_465[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_x, pa_y, pb_y, kk0_467, kk1_467, \
                         li_251, li_252, lk_324, lk_467, mi_363, \
                         mi_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_12 * li_251[k]
                   + pb_y[k] * mi_363[k];

        t_467[k] = f_29 * kk0_467[k]
                   - f_30 * kk1_467[k]
                   + pa_x[k] * lk_467[k];

        t_468[k] = pa_y[k] * lk_324[k];

        t_469[k] = f_11 * li_252[k]
                   + pb_y[k] * mi_364[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, pa_y, pb_y, li_253, li_254, \
                         li_255, lk_326, lk_327, lk_329, lk_330, \
                         mi_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = pa_y[k] * lk_326[k];

        t_471[k] = f_12 * li_253[k]
                   + pa_y[k] * lk_327[k];

        t_472[k] = f_11 * li_254[k]
                   + pb_y[k] * mi_366[k];

        t_473[k] = pa_y[k] * lk_329[k];

        t_474[k] = f_13 * li_255[k]
                   + pa_y[k] * lk_330[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_y, pb_y, pb_z, li_227, li_257, li_258, \
                         lk_333, lk_334, mi_367, mi_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_13 * li_227[k]
                   + pb_z[k] * mi_367[k];

        t_476[k] = f_11 * li_257[k]
                   + pb_y[k] * mi_369[k];

        t_477[k] = pa_y[k] * lk_333[k];

        t_478[k] = f_14 * li_258[k]
                   + pa_y[k] * lk_334[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, pa_y, pb_y, pb_z, li_230, li_260, li_261, \
                         lk_336, lk_338, mi_370, mi_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_13 * li_230[k]
                   + pb_z[k] * mi_370[k];

        t_480[k] = f_12 * li_260[k]
                   + pa_y[k] * lk_336[k];

        t_481[k] = f_11 * li_261[k]
                   + pb_y[k] * mi_373[k];

        t_482[k] = pa_y[k] * lk_338[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pa_y, pb_z, li_234, li_262, li_264, \
                         li_265, lk_339, lk_341, lk_342, mi_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_15 * li_262[k]
                   + pa_y[k] * lk_339[k];

        t_484[k] = f_13 * li_234[k]
                   + pb_z[k] * mi_374[k];

        t_485[k] = f_13 * li_264[k]
                   + pa_y[k] * lk_341[k];

        t_486[k] = f_12 * li_265[k]
                   + pa_y[k] * lk_342[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, pa_y, pb_x, pb_y, li_266, li_385, li_386, \
                         lk_344, mi_378, mi_385, mi_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_11 * li_266[k]
                   + pb_y[k] * mi_378[k];

        t_488[k] = pa_y[k] * lk_344[k];

        t_489[k] = f_15 * li_385[k]
                   + pb_x[k] * mi_385[k];

        t_490[k] = f_15 * li_386[k]
                   + pb_x[k] * mi_386[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pa_y, pb_x, li_387, li_388, \
                         li_389, li_390, lk_351, mi_387, mi_388, mi_389, \
                         mi_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_15 * li_387[k]
                   + pb_x[k] * mi_387[k];

        t_492[k] = f_15 * li_388[k]
                   + pb_x[k] * mi_388[k];

        t_493[k] = f_15 * li_389[k]
                   + pb_x[k] * mi_389[k];

        t_494[k] = f_15 * li_390[k]
                   + pb_x[k] * mi_390[k];

        t_495[k] = pa_y[k] * lk_351[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pa_y, pb_z, li_245, li_273, li_275, \
                         li_276, lk_352, lk_354, lk_355, mi_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_17 * li_273[k]
                   + pa_y[k] * lk_352[k];

        t_497[k] = f_13 * li_245[k]
                   + pb_z[k] * mi_385[k];

        t_498[k] = f_15 * li_275[k]
                   + pa_y[k] * lk_354[k];

        t_499[k] = f_14 * li_276[k]
                   + pa_y[k] * lk_355[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pa_y, pb_y, li_277, li_278, li_279, \
                         lk_356, lk_357, lk_359, mi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_13 * li_277[k]
                   + pa_y[k] * lk_356[k];

        t_501[k] = f_12 * li_278[k]
                   + pa_y[k] * lk_357[k];

        t_502[k] = f_11 * li_279[k]
                   + pb_y[k] * mi_391[k];

        t_503[k] = pa_y[k] * lk_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pa_z, pb_y, pb_z, kk0_180, kk1_180, \
                         li_252, lk_324, mh0_294, mh1_294, mi_392, \
                         mi_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_27 * kk0_180[k]
                   - f_28 * kk1_180[k]
                   + pa_z[k] * lk_324[k];

        t_505[k] = pb_y[k] * mi_392[k];

        t_506[k] = f_14 * li_252[k]
                   + pb_z[k] * mi_392[k];

        t_507[k] = f_3 * mh0_294[k]
                   - f_4 * mh1_294[k]
                   + pb_y[k] * mi_393[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pb_x, pb_y, pb_z, li_255, li_397, \
                         mh0_295, mh0_299, mh1_295, mh1_299, mi_394, mi_395, \
                         mi_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = pb_y[k] * mi_394[k];

        t_509[k] = f_15 * li_397[k]
                   + f_9 * mh0_299[k]
                   - f_10 * mh1_299[k]
                   + pb_x[k] * mi_397[k];

        t_510[k] = f_5 * mh0_295[k]
                   - f_6 * mh1_295[k]
                   + pb_y[k] * mi_395[k];

        t_511[k] = f_14 * li_255[k]
                   + pb_z[k] * mi_395[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pb_x, pb_y, pb_z, li_258, li_401, \
                         mh0_297, mh0_303, mh1_297, mh1_303, mi_397, mi_398, \
                         mi_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pb_y[k] * mi_397[k];

        t_513[k] = f_15 * li_401[k]
                   + f_7 * mh0_303[k]
                   - f_8 * mh1_303[k]
                   + pb_x[k] * mi_401[k];

        t_514[k] = f_7 * mh0_297[k]
                   - f_8 * mh1_297[k]
                   + pb_y[k] * mi_398[k];

        t_515[k] = f_14 * li_258[k]
                   + pb_z[k] * mi_398[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, pb_x, pb_y, li_406, mh0_299, mh0_308, mh1_299, \
                         mh1_308, mi_400, mi_401, mi_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_3 * mh0_299[k]
                   - f_4 * mh1_299[k]
                   + pb_y[k] * mi_400[k];

        t_517[k] = pb_y[k] * mi_401[k];

        t_518[k] = f_15 * li_406[k]
                   + f_5 * mh0_308[k]
                   - f_6 * mh1_308[k]
                   + pb_x[k] * mi_406[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, pb_y, pb_z, li_262, mh0_300, mh0_302, \
                         mh0_303, mh1_300, mh1_302, mh1_303, mi_402, mi_404, \
                         mi_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_9 * mh0_300[k]
                   - f_10 * mh1_300[k]
                   + pb_y[k] * mi_402[k];

        t_520[k] = f_14 * li_262[k]
                   + pb_z[k] * mi_402[k];

        t_521[k] = f_5 * mh0_302[k]
                   - f_6 * mh1_302[k]
                   + pb_y[k] * mi_404[k];

        t_522[k] = f_3 * mh0_303[k]
                   - f_4 * mh1_303[k]
                   + pb_y[k] * mi_405[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pb_x, pb_y, li_412, li_413, li_414, \
                         mh0_314, mh1_314, mi_406, mi_412, mi_413, \
                         mi_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = pb_y[k] * mi_406[k];

        t_524[k] = f_15 * li_412[k]
                   + f_3 * mh0_314[k]
                   - f_4 * mh1_314[k]
                   + pb_x[k] * mi_412[k];

        t_525[k] = f_15 * li_413[k]
                   + pb_x[k] * mi_413[k];

        t_526[k] = f_15 * li_414[k]
                   + pb_x[k] * mi_414[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, pb_x, pb_y, li_415, li_416, \
                         li_417, li_419, mi_412, mi_415, mi_416, mi_417, \
                         mi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_15 * li_415[k]
                   + pb_x[k] * mi_415[k];

        t_528[k] = f_15 * li_416[k]
                   + pb_x[k] * mi_416[k];

        t_529[k] = f_15 * li_417[k]
                   + pb_x[k] * mi_417[k];

        t_530[k] = pb_y[k] * mi_412[k];

        t_531[k] = f_15 * li_419[k]
                   + pb_x[k] * mi_419[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pb_y, pb_z, li_273, mh0_309, mh0_311, \
                         mh0_312, mh1_309, mh1_311, mh1_312, mi_413, mi_415, \
                         mi_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_1 * mh0_309[k]
                   - f_2 * mh1_309[k]
                   + pb_y[k] * mi_413[k];

        t_533[k] = f_14 * li_273[k]
                   + pb_z[k] * mi_413[k];

        t_534[k] = f_9 * mh0_311[k]
                   - f_10 * mh1_311[k]
                   + pb_y[k] * mi_415[k];

        t_535[k] = f_7 * mh0_312[k]
                   - f_8 * mh1_312[k]
                   + pb_y[k] * mi_416[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kk0,
                                            const size_t kk1, const size_t li, const size_t lk,
                                            const size_t mh0, const size_t mh1, const size_t mi,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_17 = 3.5 / p;
    const auto f_18 = 0.5 / alpha;
    const auto f_19 = 0.5 * beta / (alpha * p);
    const auto f_22 = 1.0 / alpha;
    const auto f_23 = beta / (alpha * p);
    const auto f_27 = 1.5 / alpha;
    const auto f_28 = 1.5 * beta / (alpha * p);
    const auto f_29 = 2.0 / alpha;
    const auto f_30 = 2.0 * beta / (alpha * p);

    auto *t_536 = buffer.data(target + 536);
    auto *t_537 = buffer.data(target + 537);
    auto *t_538 = buffer.data(target + 538);
    auto *t_539 = buffer.data(target + 539);
    auto *t_540 = buffer.data(target + 540);
    auto *t_541 = buffer.data(target + 541);
    auto *t_542 = buffer.data(target + 542);
    auto *t_543 = buffer.data(target + 543);
    auto *t_544 = buffer.data(target + 544);
    auto *t_545 = buffer.data(target + 545);
    auto *t_546 = buffer.data(target + 546);
    auto *t_547 = buffer.data(target + 547);
    auto *t_548 = buffer.data(target + 548);
    auto *t_549 = buffer.data(target + 549);
    auto *t_550 = buffer.data(target + 550);
    auto *t_551 = buffer.data(target + 551);
    auto *t_552 = buffer.data(target + 552);
    auto *t_553 = buffer.data(target + 553);
    auto *t_554 = buffer.data(target + 554);
    auto *t_555 = buffer.data(target + 555);
    auto *t_556 = buffer.data(target + 556);
    auto *t_557 = buffer.data(target + 557);
    auto *t_558 = buffer.data(target + 558);
    auto *t_559 = buffer.data(target + 559);
    auto *t_560 = buffer.data(target + 560);
    auto *t_561 = buffer.data(target + 561);
    auto *t_562 = buffer.data(target + 562);
    auto *t_563 = buffer.data(target + 563);
    auto *t_564 = buffer.data(target + 564);
    auto *t_565 = buffer.data(target + 565);
    auto *t_566 = buffer.data(target + 566);
    auto *t_567 = buffer.data(target + 567);
    auto *t_568 = buffer.data(target + 568);
    auto *t_569 = buffer.data(target + 569);
    auto *t_570 = buffer.data(target + 570);
    auto *t_571 = buffer.data(target + 571);
    auto *t_572 = buffer.data(target + 572);
    auto *t_573 = buffer.data(target + 573);
    auto *t_574 = buffer.data(target + 574);
    auto *t_575 = buffer.data(target + 575);
    auto *t_576 = buffer.data(target + 576);
    auto *t_577 = buffer.data(target + 577);
    auto *t_578 = buffer.data(target + 578);
    auto *t_579 = buffer.data(target + 579);
    auto *t_580 = buffer.data(target + 580);
    auto *t_581 = buffer.data(target + 581);
    auto *t_582 = buffer.data(target + 582);
    auto *t_583 = buffer.data(target + 583);
    auto *t_584 = buffer.data(target + 584);
    auto *t_585 = buffer.data(target + 585);
    auto *t_586 = buffer.data(target + 586);
    auto *t_587 = buffer.data(target + 587);
    auto *t_588 = buffer.data(target + 588);
    auto *t_589 = buffer.data(target + 589);
    auto *t_590 = buffer.data(target + 590);
    auto *t_591 = buffer.data(target + 591);
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);
    auto *t_604 = buffer.data(target + 604);
    auto *t_605 = buffer.data(target + 605);
    auto *t_606 = buffer.data(target + 606);
    auto *t_607 = buffer.data(target + 607);
    auto *t_608 = buffer.data(target + 608);
    auto *t_609 = buffer.data(target + 609);
    auto *t_610 = buffer.data(target + 610);
    auto *t_611 = buffer.data(target + 611);
    auto *t_612 = buffer.data(target + 612);
    auto *t_613 = buffer.data(target + 613);
    auto *t_614 = buffer.data(target + 614);
    auto *t_615 = buffer.data(target + 615);
    auto *t_616 = buffer.data(target + 616);
    auto *t_617 = buffer.data(target + 617);
    auto *t_618 = buffer.data(target + 618);
    auto *t_619 = buffer.data(target + 619);
    auto *t_620 = buffer.data(target + 620);
    auto *t_621 = buffer.data(target + 621);
    auto *t_622 = buffer.data(target + 622);
    auto *t_623 = buffer.data(target + 623);
    auto *t_624 = buffer.data(target + 624);
    auto *t_625 = buffer.data(target + 625);
    auto *t_626 = buffer.data(target + 626);
    auto *t_627 = buffer.data(target + 627);
    auto *t_628 = buffer.data(target + 628);
    auto *t_629 = buffer.data(target + 629);
    auto *t_630 = buffer.data(target + 630);
    auto *t_631 = buffer.data(target + 631);
    auto *t_632 = buffer.data(target + 632);
    auto *t_633 = buffer.data(target + 633);
    auto *t_634 = buffer.data(target + 634);
    auto *t_635 = buffer.data(target + 635);
    auto *t_636 = buffer.data(target + 636);
    auto *t_637 = buffer.data(target + 637);
    auto *t_638 = buffer.data(target + 638);
    auto *t_639 = buffer.data(target + 639);
    auto *t_640 = buffer.data(target + 640);
    auto *t_641 = buffer.data(target + 641);
    auto *t_642 = buffer.data(target + 642);
    auto *t_643 = buffer.data(target + 643);
    auto *t_644 = buffer.data(target + 644);
    auto *t_645 = buffer.data(target + 645);
    auto *t_646 = buffer.data(target + 646);
    auto *t_647 = buffer.data(target + 647);
    auto *t_648 = buffer.data(target + 648);
    auto *t_649 = buffer.data(target + 649);
    auto *t_650 = buffer.data(target + 650);
    auto *t_651 = buffer.data(target + 651);
    auto *t_652 = buffer.data(target + 652);
    auto *t_653 = buffer.data(target + 653);
    auto *t_654 = buffer.data(target + 654);
    auto *t_655 = buffer.data(target + 655);
    auto *t_656 = buffer.data(target + 656);
    auto *t_657 = buffer.data(target + 657);
    auto *t_658 = buffer.data(target + 658);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kk0_216 = buffer.data(kk0 + 216);
    const auto *kk0_219 = buffer.data(kk0 + 219);
    const auto *kk0_222 = buffer.data(kk0 + 222);
    const auto *kk0_226 = buffer.data(kk0 + 226);
    const auto *kk0_231 = buffer.data(kk0 + 231);
    const auto *kk0_255 = buffer.data(kk0 + 255);
    const auto *kk0_258 = buffer.data(kk0 + 258);
    const auto *kk0_262 = buffer.data(kk0 + 262);
    const auto *kk0_288 = buffer.data(kk0 + 288);
    const auto *kk0_293 = buffer.data(kk0 + 293);
    const auto *kk0_297 = buffer.data(kk0 + 297);
    const auto *kk0_302 = buffer.data(kk0 + 302);
    const auto *kk0_308 = buffer.data(kk0 + 308);
    const auto *kk0_324 = buffer.data(kk0 + 324);
    const auto *kk0_329 = buffer.data(kk0 + 329);
    const auto *kk0_333 = buffer.data(kk0 + 333);
    const auto *kk0_539 = buffer.data(kk0 + 539);
    const auto *kk0_568 = buffer.data(kk0 + 568);
    const auto *kk0_640 = buffer.data(kk0 + 640);
    const auto *kk0_642 = buffer.data(kk0 + 642);
    const auto *kk0_643 = buffer.data(kk0 + 643);
    const auto *kk0_644 = buffer.data(kk0 + 644);
    const auto *kk0_645 = buffer.data(kk0 + 645);
    const auto *kk0_647 = buffer.data(kk0 + 647);

    const auto *kk1_216 = buffer.data(kk1 + 216);
    const auto *kk1_219 = buffer.data(kk1 + 219);
    const auto *kk1_222 = buffer.data(kk1 + 222);
    const auto *kk1_226 = buffer.data(kk1 + 226);
    const auto *kk1_231 = buffer.data(kk1 + 231);
    const auto *kk1_255 = buffer.data(kk1 + 255);
    const auto *kk1_258 = buffer.data(kk1 + 258);
    const auto *kk1_262 = buffer.data(kk1 + 262);
    const auto *kk1_288 = buffer.data(kk1 + 288);
    const auto *kk1_293 = buffer.data(kk1 + 293);
    const auto *kk1_297 = buffer.data(kk1 + 297);
    const auto *kk1_302 = buffer.data(kk1 + 302);
    const auto *kk1_308 = buffer.data(kk1 + 308);
    const auto *kk1_324 = buffer.data(kk1 + 324);
    const auto *kk1_329 = buffer.data(kk1 + 329);
    const auto *kk1_333 = buffer.data(kk1 + 333);
    const auto *kk1_539 = buffer.data(kk1 + 539);
    const auto *kk1_568 = buffer.data(kk1 + 568);
    const auto *kk1_640 = buffer.data(kk1 + 640);
    const auto *kk1_642 = buffer.data(kk1 + 642);
    const auto *kk1_643 = buffer.data(kk1 + 643);
    const auto *kk1_644 = buffer.data(kk1 + 644);
    const auto *kk1_645 = buffer.data(kk1 + 645);
    const auto *kk1_647 = buffer.data(kk1 + 647);

    const auto *li_280 = buffer.data(li + 280);
    const auto *li_282 = buffer.data(li + 282);
    const auto *li_283 = buffer.data(li + 283);
    const auto *li_285 = buffer.data(li + 285);
    const auto *li_286 = buffer.data(li + 286);
    const auto *li_287 = buffer.data(li + 287);
    const auto *li_289 = buffer.data(li + 289);
    const auto *li_290 = buffer.data(li + 290);
    const auto *li_291 = buffer.data(li + 291);
    const auto *li_292 = buffer.data(li + 292);
    const auto *li_294 = buffer.data(li + 294);
    const auto *li_301 = buffer.data(li + 301);
    const auto *li_302 = buffer.data(li + 302);
    const auto *li_303 = buffer.data(li + 303);
    const auto *li_304 = buffer.data(li + 304);
    const auto *li_305 = buffer.data(li + 305);
    const auto *li_307 = buffer.data(li + 307);
    const auto *li_308 = buffer.data(li + 308);
    const auto *li_310 = buffer.data(li + 310);
    const auto *li_311 = buffer.data(li + 311);
    const auto *li_313 = buffer.data(li + 313);
    const auto *li_314 = buffer.data(li + 314);
    const auto *li_317 = buffer.data(li + 317);
    const auto *li_318 = buffer.data(li + 318);
    const auto *li_322 = buffer.data(li + 322);
    const auto *li_329 = buffer.data(li + 329);
    const auto *li_335 = buffer.data(li + 335);
    const auto *li_336 = buffer.data(li + 336);
    const auto *li_338 = buffer.data(li + 338);
    const auto *li_339 = buffer.data(li + 339);
    const auto *li_341 = buffer.data(li + 341);
    const auto *li_345 = buffer.data(li + 345);
    const auto *li_350 = buffer.data(li + 350);
    const auto *li_363 = buffer.data(li + 363);
    const auto *li_364 = buffer.data(li + 364);
    const auto *li_366 = buffer.data(li + 366);
    const auto *li_369 = buffer.data(li + 369);
    const auto *li_423 = buffer.data(li + 423);
    const auto *li_426 = buffer.data(li + 426);
    const auto *li_430 = buffer.data(li + 430);
    const auto *li_435 = buffer.data(li + 435);
    const auto *li_441 = buffer.data(li + 441);
    const auto *li_443 = buffer.data(li + 443);
    const auto *li_444 = buffer.data(li + 444);
    const auto *li_445 = buffer.data(li + 445);
    const auto *li_446 = buffer.data(li + 446);
    const auto *li_447 = buffer.data(li + 447);
    const auto *li_470 = buffer.data(li + 470);
    const auto *li_471 = buffer.data(li + 471);
    const auto *li_472 = buffer.data(li + 472);
    const auto *li_473 = buffer.data(li + 473);
    const auto *li_474 = buffer.data(li + 474);
    const auto *li_475 = buffer.data(li + 475);
    const auto *li_488 = buffer.data(li + 488);
    const auto *li_493 = buffer.data(li + 493);
    const auto *li_494 = buffer.data(li + 494);
    const auto *li_497 = buffer.data(li + 497);
    const auto *li_498 = buffer.data(li + 498);
    const auto *li_499 = buffer.data(li + 499);
    const auto *li_500 = buffer.data(li + 500);
    const auto *li_501 = buffer.data(li + 501);
    const auto *li_502 = buffer.data(li + 502);
    const auto *li_503 = buffer.data(li + 503);

    const auto *lk_360 = buffer.data(lk + 360);
    const auto *lk_361 = buffer.data(lk + 361);
    const auto *lk_363 = buffer.data(lk + 363);
    const auto *lk_365 = buffer.data(lk + 365);
    const auto *lk_366 = buffer.data(lk + 366);
    const auto *lk_369 = buffer.data(lk + 369);
    const auto *lk_370 = buffer.data(lk + 370);
    const auto *lk_372 = buffer.data(lk + 372);
    const auto *lk_374 = buffer.data(lk + 374);
    const auto *lk_375 = buffer.data(lk + 375);
    const auto *lk_377 = buffer.data(lk + 377);
    const auto *lk_378 = buffer.data(lk + 378);
    const auto *lk_380 = buffer.data(lk + 380);
    const auto *lk_381 = buffer.data(lk + 381);
    const auto *lk_388 = buffer.data(lk + 388);
    const auto *lk_390 = buffer.data(lk + 390);
    const auto *lk_391 = buffer.data(lk + 391);
    const auto *lk_392 = buffer.data(lk + 392);
    const auto *lk_393 = buffer.data(lk + 393);
    const auto *lk_395 = buffer.data(lk + 395);
    const auto *lk_399 = buffer.data(lk + 399);
    const auto *lk_402 = buffer.data(lk + 402);
    const auto *lk_406 = buffer.data(lk + 406);
    const auto *lk_411 = buffer.data(lk + 411);
    const auto *lk_432 = buffer.data(lk + 432);
    const auto *lk_435 = buffer.data(lk + 435);
    const auto *lk_437 = buffer.data(lk + 437);
    const auto *lk_438 = buffer.data(lk + 438);
    const auto *lk_441 = buffer.data(lk + 441);
    const auto *lk_442 = buffer.data(lk + 442);
    const auto *lk_446 = buffer.data(lk + 446);
    const auto *lk_452 = buffer.data(lk + 452);
    const auto *lk_468 = buffer.data(lk + 468);
    const auto *lk_473 = buffer.data(lk + 473);
    const auto *lk_477 = buffer.data(lk + 477);
    const auto *lk_539 = buffer.data(lk + 539);
    const auto *lk_568 = buffer.data(lk + 568);
    const auto *lk_640 = buffer.data(lk + 640);
    const auto *lk_642 = buffer.data(lk + 642);
    const auto *lk_643 = buffer.data(lk + 643);
    const auto *lk_644 = buffer.data(lk + 644);
    const auto *lk_645 = buffer.data(lk + 645);
    const auto *lk_647 = buffer.data(lk + 647);

    const auto *mh0_313 = buffer.data(mh0 + 313);
    const auto *mh0_314 = buffer.data(mh0 + 314);
    const auto *mh0_315 = buffer.data(mh0 + 315);
    const auto *mh0_317 = buffer.data(mh0 + 317);
    const auto *mh0_318 = buffer.data(mh0 + 318);
    const auto *mh0_320 = buffer.data(mh0 + 320);
    const auto *mh0_321 = buffer.data(mh0 + 321);
    const auto *mh0_322 = buffer.data(mh0 + 322);
    const auto *mh0_324 = buffer.data(mh0 + 324);
    const auto *mh0_325 = buffer.data(mh0 + 325);
    const auto *mh0_330 = buffer.data(mh0 + 330);
    const auto *mh0_331 = buffer.data(mh0 + 331);
    const auto *mh0_332 = buffer.data(mh0 + 332);
    const auto *mh0_333 = buffer.data(mh0 + 333);
    const auto *mh0_335 = buffer.data(mh0 + 335);
    const auto *mh0_369 = buffer.data(mh0 + 369);
    const auto *mh0_374 = buffer.data(mh0 + 374);
    const auto *mh0_375 = buffer.data(mh0 + 375);

    const auto *mh1_313 = buffer.data(mh1 + 313);
    const auto *mh1_314 = buffer.data(mh1 + 314);
    const auto *mh1_315 = buffer.data(mh1 + 315);
    const auto *mh1_317 = buffer.data(mh1 + 317);
    const auto *mh1_318 = buffer.data(mh1 + 318);
    const auto *mh1_320 = buffer.data(mh1 + 320);
    const auto *mh1_321 = buffer.data(mh1 + 321);
    const auto *mh1_322 = buffer.data(mh1 + 322);
    const auto *mh1_324 = buffer.data(mh1 + 324);
    const auto *mh1_325 = buffer.data(mh1 + 325);
    const auto *mh1_330 = buffer.data(mh1 + 330);
    const auto *mh1_331 = buffer.data(mh1 + 331);
    const auto *mh1_332 = buffer.data(mh1 + 332);
    const auto *mh1_333 = buffer.data(mh1 + 333);
    const auto *mh1_335 = buffer.data(mh1 + 335);
    const auto *mh1_369 = buffer.data(mh1 + 369);
    const auto *mh1_374 = buffer.data(mh1 + 374);
    const auto *mh1_375 = buffer.data(mh1 + 375);

    const auto *mi_417 = buffer.data(mi + 417);
    const auto *mi_418 = buffer.data(mi + 418);
    const auto *mi_419 = buffer.data(mi + 419);
    const auto *mi_420 = buffer.data(mi + 420);
    const auto *mi_421 = buffer.data(mi + 421);
    const auto *mi_422 = buffer.data(mi + 422);
    const auto *mi_423 = buffer.data(mi + 423);
    const auto *mi_425 = buffer.data(mi + 425);
    const auto *mi_426 = buffer.data(mi + 426);
    const auto *mi_427 = buffer.data(mi + 427);
    const auto *mi_429 = buffer.data(mi + 429);
    const auto *mi_430 = buffer.data(mi + 430);
    const auto *mi_431 = buffer.data(mi + 431);
    const auto *mi_432 = buffer.data(mi + 432);
    const auto *mi_434 = buffer.data(mi + 434);
    const auto *mi_435 = buffer.data(mi + 435);
    const auto *mi_441 = buffer.data(mi + 441);
    const auto *mi_442 = buffer.data(mi + 442);
    const auto *mi_443 = buffer.data(mi + 443);
    const auto *mi_444 = buffer.data(mi + 444);
    const auto *mi_445 = buffer.data(mi + 445);
    const auto *mi_446 = buffer.data(mi + 446);
    const auto *mi_447 = buffer.data(mi + 447);
    const auto *mi_448 = buffer.data(mi + 448);
    const auto *mi_450 = buffer.data(mi + 450);
    const auto *mi_451 = buffer.data(mi + 451);
    const auto *mi_453 = buffer.data(mi + 453);
    const auto *mi_454 = buffer.data(mi + 454);
    const auto *mi_457 = buffer.data(mi + 457);
    const auto *mi_458 = buffer.data(mi + 458);
    const auto *mi_462 = buffer.data(mi + 462);
    const auto *mi_469 = buffer.data(mi + 469);
    const auto *mi_470 = buffer.data(mi + 470);
    const auto *mi_471 = buffer.data(mi + 471);
    const auto *mi_472 = buffer.data(mi + 472);
    const auto *mi_473 = buffer.data(mi + 473);
    const auto *mi_474 = buffer.data(mi + 474);
    const auto *mi_475 = buffer.data(mi + 475);
    const auto *mi_476 = buffer.data(mi + 476);
    const auto *mi_478 = buffer.data(mi + 478);
    const auto *mi_479 = buffer.data(mi + 479);
    const auto *mi_481 = buffer.data(mi + 481);
    const auto *mi_482 = buffer.data(mi + 482);
    const auto *mi_485 = buffer.data(mi + 485);
    const auto *mi_486 = buffer.data(mi + 486);
    const auto *mi_488 = buffer.data(mi + 488);
    const auto *mi_490 = buffer.data(mi + 490);
    const auto *mi_493 = buffer.data(mi + 493);
    const auto *mi_494 = buffer.data(mi + 494);
    const auto *mi_497 = buffer.data(mi + 497);
    const auto *mi_498 = buffer.data(mi + 498);
    const auto *mi_499 = buffer.data(mi + 499);
    const auto *mi_500 = buffer.data(mi + 500);
    const auto *mi_501 = buffer.data(mi + 501);
    const auto *mi_502 = buffer.data(mi + 502);
    const auto *mi_503 = buffer.data(mi + 503);
    const auto *mi_504 = buffer.data(mi + 504);
    const auto *mi_506 = buffer.data(mi + 506);
    const auto *mi_507 = buffer.data(mi + 507);
    const auto *mi_509 = buffer.data(mi + 509);

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pa_x, pb_y, kk0_539, kk1_539, lk_539, \
                         mh0_313, mh0_314, mh1_313, mh1_314, mi_417, mi_418, \
                         mi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_5 * mh0_313[k]
                   - f_6 * mh1_313[k]
                   + pb_y[k] * mi_417[k];

        t_537[k] = f_3 * mh0_314[k]
                   - f_4 * mh1_314[k]
                   + pb_y[k] * mi_418[k];

        t_538[k] = pb_y[k] * mi_419[k];

        t_539[k] = f_29 * kk0_539[k]
                   - f_30 * kk1_539[k]
                   + pa_x[k] * lk_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, pa_y, pb_y, pb_z, kk0_216, kk1_216, li_280, \
                         lk_360, mi_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_29 * kk0_216[k]
                   - f_30 * kk1_216[k]
                   + pa_y[k] * lk_360[k];

        t_541[k] = f_15 * li_280[k]
                   + pb_y[k] * mi_420[k];

        t_542[k] = pb_z[k] * mi_420[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, pb_x, pb_z, li_423, mh0_315, mh0_318, mh1_315, \
                         mh1_318, mi_421, mi_422, mi_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_14 * li_423[k]
                   + f_9 * mh0_318[k]
                   - f_10 * mh1_318[k]
                   + pb_x[k] * mi_423[k];

        t_544[k] = pb_z[k] * mi_421[k];

        t_545[k] = f_3 * mh0_315[k]
                   - f_4 * mh1_315[k]
                   + pb_z[k] * mi_422[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, pb_x, pb_y, pb_z, li_285, li_426, \
                         mh0_317, mh0_321, mh1_317, mh1_321, mi_423, mi_425, \
                         mi_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = f_14 * li_426[k]
                   + f_7 * mh0_321[k]
                   - f_8 * mh1_321[k]
                   + pb_x[k] * mi_426[k];

        t_547[k] = pb_z[k] * mi_423[k];

        t_548[k] = f_15 * li_285[k]
                   + pb_y[k] * mi_425[k];

        t_549[k] = f_5 * mh0_317[k]
                   - f_6 * mh1_317[k]
                   + pb_z[k] * mi_425[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, pb_x, pb_z, li_430, mh0_318, mh0_325, mh1_318, \
                         mh1_325, mi_426, mi_427, mi_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_14 * li_430[k]
                   + f_5 * mh0_325[k]
                   - f_6 * mh1_325[k]
                   + pb_x[k] * mi_430[k];

        t_551[k] = pb_z[k] * mi_426[k];

        t_552[k] = f_3 * mh0_318[k]
                   - f_4 * mh1_318[k]
                   + pb_z[k] * mi_427[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, pb_x, pb_y, pb_z, li_289, li_435, \
                         mh0_320, mh0_330, mh1_320, mh1_330, mi_429, mi_430, \
                         mi_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_15 * li_289[k]
                   + pb_y[k] * mi_429[k];

        t_554[k] = f_7 * mh0_320[k]
                   - f_8 * mh1_320[k]
                   + pb_z[k] * mi_429[k];

        t_555[k] = f_14 * li_435[k]
                   + f_3 * mh0_330[k]
                   - f_4 * mh1_330[k]
                   + pb_x[k] * mi_435[k];

        t_556[k] = pb_z[k] * mi_430[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, pb_y, pb_z, li_294, mh0_321, mh0_322, \
                         mh0_324, mh1_321, mh1_322, mh1_324, mi_431, mi_432, \
                         mi_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_3 * mh0_321[k]
                   - f_4 * mh1_321[k]
                   + pb_z[k] * mi_431[k];

        t_558[k] = f_5 * mh0_322[k]
                   - f_6 * mh1_322[k]
                   + pb_z[k] * mi_432[k];

        t_559[k] = f_15 * li_294[k]
                   + pb_y[k] * mi_434[k];

        t_560[k] = f_9 * mh0_324[k]
                   - f_10 * mh1_324[k]
                   + pb_z[k] * mi_434[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, t_564, t_565, pb_x, pb_z, li_441, li_443, \
                         li_444, li_445, mi_435, mi_441, mi_443, mi_444, \
                         mi_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_14 * li_441[k]
                   + pb_x[k] * mi_441[k];

        t_562[k] = pb_z[k] * mi_435[k];

        t_563[k] = f_14 * li_443[k]
                   + pb_x[k] * mi_443[k];

        t_564[k] = f_14 * li_444[k]
                   + pb_x[k] * mi_444[k];

        t_565[k] = f_14 * li_445[k]
                   + pb_x[k] * mi_445[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pa_x, pb_x, pb_z, kk0_568, kk1_568, \
                         li_446, li_447, lk_568, mi_441, mi_446, \
                         mi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_14 * li_446[k]
                   + pb_x[k] * mi_446[k];

        t_567[k] = f_14 * li_447[k]
                   + pb_x[k] * mi_447[k];

        t_568[k] = f_27 * kk0_568[k]
                   - f_28 * kk1_568[k]
                   + pa_x[k] * lk_568[k];

        t_569[k] = pb_z[k] * mi_441[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, pb_z, mh0_330, mh0_331, mh0_332, mh1_330, \
                         mh1_331, mh1_332, mi_442, mi_443, mi_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_3 * mh0_330[k]
                   - f_4 * mh1_330[k]
                   + pb_z[k] * mi_442[k];

        t_571[k] = f_5 * mh0_331[k]
                   - f_6 * mh1_331[k]
                   + pb_z[k] * mi_443[k];

        t_572[k] = f_7 * mh0_332[k]
                   - f_8 * mh1_332[k]
                   + pb_z[k] * mi_444[k];
    }

#pragma omp simd aligned(t_573, t_574, t_575, t_576, pa_z, pb_y, pb_z, li_307, lk_360, \
                         mh0_333, mh0_335, mh1_333, mh1_335, mi_445, \
                         mi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_573[k] = f_9 * mh0_333[k]
                   - f_10 * mh1_333[k]
                   + pb_z[k] * mi_445[k];

        t_574[k] = f_15 * li_307[k]
                   + pb_y[k] * mi_447[k];

        t_575[k] = f_1 * mh0_335[k]
                   - f_2 * mh1_335[k]
                   + pb_z[k] * mi_447[k];

        t_576[k] = pa_z[k] * lk_360[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, pa_z, pb_y, pb_z, li_280, li_282, \
                         li_310, lk_361, lk_363, lk_365, mi_448, \
                         mi_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = pa_z[k] * lk_361[k];

        t_578[k] = f_11 * li_280[k]
                   + pb_z[k] * mi_448[k];

        t_579[k] = pa_z[k] * lk_363[k];

        t_580[k] = f_14 * li_310[k]
                   + pb_y[k] * mi_450[k];

        t_581[k] = f_12 * li_282[k]
                   + pa_z[k] * lk_365[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, pa_z, pb_y, pb_z, li_283, li_285, \
                         li_313, lk_366, lk_369, lk_370, mi_451, \
                         mi_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = pa_z[k] * lk_366[k];

        t_583[k] = f_11 * li_283[k]
                   + pb_z[k] * mi_451[k];

        t_584[k] = f_14 * li_313[k]
                   + pb_y[k] * mi_453[k];

        t_585[k] = f_13 * li_285[k]
                   + pa_z[k] * lk_369[k];

        t_586[k] = pa_z[k] * lk_370[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, t_590, pa_z, pb_y, pb_z, li_286, li_287, li_289, \
                         li_317, lk_372, lk_374, mi_454, mi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_11 * li_286[k]
                   + pb_z[k] * mi_454[k];

        t_588[k] = f_12 * li_287[k]
                   + pa_z[k] * lk_372[k];

        t_589[k] = f_14 * li_317[k]
                   + pb_y[k] * mi_457[k];

        t_590[k] = f_14 * li_289[k]
                   + pa_z[k] * lk_374[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, pa_z, pb_z, li_290, li_291, li_292, \
                         lk_375, lk_377, lk_378, mi_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = pa_z[k] * lk_375[k];

        t_592[k] = f_11 * li_290[k]
                   + pb_z[k] * mi_458[k];

        t_593[k] = f_12 * li_291[k]
                   + pa_z[k] * lk_377[k];

        t_594[k] = f_13 * li_292[k]
                   + pa_z[k] * lk_378[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, pa_z, pb_x, pb_y, li_294, li_322, li_470, \
                         lk_380, lk_381, mi_462, mi_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_14 * li_322[k]
                   + pb_y[k] * mi_462[k];

        t_596[k] = f_15 * li_294[k]
                   + pa_z[k] * lk_380[k];

        t_597[k] = pa_z[k] * lk_381[k];

        t_598[k] = f_14 * li_470[k]
                   + pb_x[k] * mi_470[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, t_603, pb_x, li_471, li_472, li_473, \
                         li_474, li_475, mi_471, mi_472, mi_473, mi_474, \
                         mi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = f_14 * li_471[k]
                   + pb_x[k] * mi_471[k];

        t_600[k] = f_14 * li_472[k]
                   + pb_x[k] * mi_472[k];

        t_601[k] = f_14 * li_473[k]
                   + pb_x[k] * mi_473[k];

        t_602[k] = f_14 * li_474[k]
                   + pb_x[k] * mi_474[k];

        t_603[k] = f_14 * li_475[k]
                   + pb_x[k] * mi_475[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, pa_z, pb_z, li_301, li_302, \
                         li_303, li_304, lk_388, lk_390, lk_391, lk_392, \
                         mi_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = pa_z[k] * lk_388[k];

        t_605[k] = f_11 * li_301[k]
                   + pb_z[k] * mi_469[k];

        t_606[k] = f_12 * li_302[k]
                   + pa_z[k] * lk_390[k];

        t_607[k] = f_13 * li_303[k]
                   + pa_z[k] * lk_391[k];

        t_608[k] = f_14 * li_304[k]
                   + pa_z[k] * lk_392[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, pa_y, pa_z, pb_y, kk0_288, kk1_288, \
                         li_305, li_307, li_335, lk_393, lk_395, lk_432, \
                         mi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_15 * li_305[k]
                   + pa_z[k] * lk_393[k];

        t_610[k] = f_14 * li_335[k]
                   + pb_y[k] * mi_475[k];

        t_611[k] = f_17 * li_307[k]
                   + pa_z[k] * lk_395[k];

        t_612[k] = f_22 * kk0_288[k]
                   - f_23 * kk1_288[k]
                   + pa_y[k] * lk_432[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pa_z, pb_y, pb_z, kk0_219, kk1_219, \
                         li_308, li_336, li_338, lk_399, mi_476, \
                         mi_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_13 * li_336[k]
                   + pb_y[k] * mi_476[k];

        t_614[k] = f_12 * li_308[k]
                   + pb_z[k] * mi_476[k];

        t_615[k] = f_18 * kk0_219[k]
                   - f_19 * kk1_219[k]
                   + pa_z[k] * lk_399[k];

        t_616[k] = f_13 * li_338[k]
                   + pb_y[k] * mi_478[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, pa_y, pa_z, pb_z, kk0_222, kk0_293, kk1_222, \
                         kk1_293, li_311, lk_402, lk_437, mi_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_22 * kk0_293[k]
                   - f_23 * kk1_293[k]
                   + pa_y[k] * lk_437[k];

        t_618[k] = f_18 * kk0_222[k]
                   - f_19 * kk1_222[k]
                   + pa_z[k] * lk_402[k];

        t_619[k] = f_12 * li_311[k]
                   + pb_z[k] * mi_479[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pa_y, pa_z, pb_y, kk0_226, kk0_297, kk1_226, \
                         kk1_297, li_341, lk_406, lk_441, mi_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_13 * li_341[k]
                   + pb_y[k] * mi_481[k];

        t_621[k] = f_22 * kk0_297[k]
                   - f_23 * kk1_297[k]
                   + pa_y[k] * lk_441[k];

        t_622[k] = f_18 * kk0_226[k]
                   - f_19 * kk1_226[k]
                   + pa_z[k] * lk_406[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pb_x, pb_y, pb_z, li_314, li_345, li_488, \
                         mh0_369, mh1_369, mi_482, mi_485, mi_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_12 * li_314[k]
                   + pb_z[k] * mi_482[k];

        t_624[k] = f_14 * li_488[k]
                   + f_5 * mh0_369[k]
                   - f_6 * mh1_369[k]
                   + pb_x[k] * mi_488[k];

        t_625[k] = f_13 * li_345[k]
                   + pb_y[k] * mi_485[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pa_y, pa_z, pb_z, kk0_231, kk0_302, kk1_231, \
                         kk1_302, li_318, lk_411, lk_446, mi_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_22 * kk0_302[k]
                   - f_23 * kk1_302[k]
                   + pa_y[k] * lk_446[k];

        t_627[k] = f_18 * kk0_231[k]
                   - f_19 * kk1_231[k]
                   + pa_z[k] * lk_411[k];

        t_628[k] = f_12 * li_318[k]
                   + pb_z[k] * mi_486[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pb_x, pb_y, li_350, li_493, li_494, mh0_374, \
                         mh0_375, mh1_374, mh1_375, mi_490, mi_493, \
                         mi_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_14 * li_493[k]
                   + f_3 * mh0_374[k]
                   - f_4 * mh1_374[k]
                   + pb_x[k] * mi_493[k];

        t_630[k] = f_14 * li_494[k]
                   + f_3 * mh0_375[k]
                   - f_4 * mh1_375[k]
                   + pb_x[k] * mi_494[k];

        t_631[k] = f_13 * li_350[k]
                   + pb_y[k] * mi_490[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, pa_y, pb_x, kk0_308, kk1_308, li_497, \
                         li_498, li_499, lk_452, mi_497, mi_498, \
                         mi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_22 * kk0_308[k]
                   - f_23 * kk1_308[k]
                   + pa_y[k] * lk_452[k];

        t_633[k] = f_14 * li_497[k]
                   + pb_x[k] * mi_497[k];

        t_634[k] = f_14 * li_498[k]
                   + pb_x[k] * mi_498[k];

        t_635[k] = f_14 * li_499[k]
                   + pb_x[k] * mi_499[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pb_x, li_500, li_501, li_502, li_503, \
                         mi_500, mi_501, mi_502, mi_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_14 * li_500[k]
                   + pb_x[k] * mi_500[k];

        t_637[k] = f_14 * li_501[k]
                   + pb_x[k] * mi_501[k];

        t_638[k] = f_14 * li_502[k]
                   + pb_x[k] * mi_502[k];

        t_639[k] = f_14 * li_503[k]
                   + pb_x[k] * mi_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pa_x, pb_z, kk0_640, kk0_642, kk1_640, kk1_642, \
                         li_329, lk_640, lk_642, mi_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_27 * kk0_640[k]
                   - f_28 * kk1_640[k]
                   + pa_x[k] * lk_640[k];

        t_641[k] = f_12 * li_329[k]
                   + pb_z[k] * mi_497[k];

        t_642[k] = f_27 * kk0_642[k]
                   - f_28 * kk1_642[k]
                   + pa_x[k] * lk_642[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pa_x, kk0_643, kk0_644, kk0_645, kk1_643, \
                         kk1_644, kk1_645, lk_643, lk_644, lk_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_27 * kk0_643[k]
                   - f_28 * kk1_643[k]
                   + pa_x[k] * lk_643[k];

        t_644[k] = f_27 * kk0_644[k]
                   - f_28 * kk1_644[k]
                   + pa_x[k] * lk_644[k];

        t_645[k] = f_27 * kk0_645[k]
                   - f_28 * kk1_645[k]
                   + pa_x[k] * lk_645[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pa_x, pa_y, pb_y, kk0_324, kk0_647, kk1_324, \
                         kk1_647, li_363, lk_468, lk_647, mi_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_13 * li_363[k]
                   + pb_y[k] * mi_503[k];

        t_647[k] = f_27 * kk0_647[k]
                   - f_28 * kk1_647[k]
                   + pa_x[k] * lk_647[k];

        t_648[k] = f_18 * kk0_324[k]
                   - f_19 * kk1_324[k]
                   + pa_y[k] * lk_468[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pa_z, pb_y, pb_z, kk0_255, kk1_255, \
                         li_336, li_364, li_366, lk_435, mi_504, \
                         mi_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_12 * li_364[k]
                   + pb_y[k] * mi_504[k];

        t_650[k] = f_13 * li_336[k]
                   + pb_z[k] * mi_504[k];

        t_651[k] = f_22 * kk0_255[k]
                   - f_23 * kk1_255[k]
                   + pa_z[k] * lk_435[k];

        t_652[k] = f_12 * li_366[k]
                   + pb_y[k] * mi_506[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pa_y, pa_z, pb_z, kk0_258, kk0_329, kk1_258, \
                         kk1_329, li_339, lk_438, lk_473, mi_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_18 * kk0_329[k]
                   - f_19 * kk1_329[k]
                   + pa_y[k] * lk_473[k];

        t_654[k] = f_22 * kk0_258[k]
                   - f_23 * kk1_258[k]
                   + pa_z[k] * lk_438[k];

        t_655[k] = f_13 * li_339[k]
                   + pb_z[k] * mi_507[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pa_y, pa_z, pb_y, kk0_262, kk0_333, kk1_262, \
                         kk1_333, li_369, lk_442, lk_477, mi_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_12 * li_369[k]
                   + pb_y[k] * mi_509[k];

        t_657[k] = f_18 * kk0_333[k]
                   - f_19 * kk1_333[k]
                   + pa_y[k] * lk_477[k];

        t_658[k] = f_22 * kk0_262[k]
                   - f_23 * kk1_262[k]
                   + pa_z[k] * lk_442[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kk0,
                                            const size_t kk1, const size_t li, const size_t lk,
                                            const size_t mh0, const size_t mh1, const size_t mi,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_17 = 3.5 / p;
    const auto f_18 = 0.5 / alpha;
    const auto f_19 = 0.5 * beta / (alpha * p);
    const auto f_22 = 1.0 / alpha;
    const auto f_23 = beta / (alpha * p);
    const auto f_24 = 3.0 / p;
    const auto f_25 = 2.5 / alpha;
    const auto f_26 = 2.5 * beta / (alpha * p);
    const auto f_27 = 1.5 / alpha;
    const auto f_28 = 1.5 * beta / (alpha * p);
    const auto f_29 = 2.0 / alpha;
    const auto f_30 = 2.0 * beta / (alpha * p);

    auto *t_659 = buffer.data(target + 659);
    auto *t_660 = buffer.data(target + 660);
    auto *t_661 = buffer.data(target + 661);
    auto *t_662 = buffer.data(target + 662);
    auto *t_663 = buffer.data(target + 663);
    auto *t_664 = buffer.data(target + 664);
    auto *t_665 = buffer.data(target + 665);
    auto *t_666 = buffer.data(target + 666);
    auto *t_667 = buffer.data(target + 667);
    auto *t_668 = buffer.data(target + 668);
    auto *t_669 = buffer.data(target + 669);
    auto *t_670 = buffer.data(target + 670);
    auto *t_671 = buffer.data(target + 671);
    auto *t_672 = buffer.data(target + 672);
    auto *t_673 = buffer.data(target + 673);
    auto *t_674 = buffer.data(target + 674);
    auto *t_675 = buffer.data(target + 675);
    auto *t_676 = buffer.data(target + 676);
    auto *t_677 = buffer.data(target + 677);
    auto *t_678 = buffer.data(target + 678);
    auto *t_679 = buffer.data(target + 679);
    auto *t_680 = buffer.data(target + 680);
    auto *t_681 = buffer.data(target + 681);
    auto *t_682 = buffer.data(target + 682);
    auto *t_683 = buffer.data(target + 683);
    auto *t_684 = buffer.data(target + 684);
    auto *t_685 = buffer.data(target + 685);
    auto *t_686 = buffer.data(target + 686);
    auto *t_687 = buffer.data(target + 687);
    auto *t_688 = buffer.data(target + 688);
    auto *t_689 = buffer.data(target + 689);
    auto *t_690 = buffer.data(target + 690);
    auto *t_691 = buffer.data(target + 691);
    auto *t_692 = buffer.data(target + 692);
    auto *t_693 = buffer.data(target + 693);
    auto *t_694 = buffer.data(target + 694);
    auto *t_695 = buffer.data(target + 695);
    auto *t_696 = buffer.data(target + 696);
    auto *t_697 = buffer.data(target + 697);
    auto *t_698 = buffer.data(target + 698);
    auto *t_699 = buffer.data(target + 699);
    auto *t_700 = buffer.data(target + 700);
    auto *t_701 = buffer.data(target + 701);
    auto *t_702 = buffer.data(target + 702);
    auto *t_703 = buffer.data(target + 703);
    auto *t_704 = buffer.data(target + 704);
    auto *t_705 = buffer.data(target + 705);
    auto *t_706 = buffer.data(target + 706);
    auto *t_707 = buffer.data(target + 707);
    auto *t_708 = buffer.data(target + 708);
    auto *t_709 = buffer.data(target + 709);
    auto *t_710 = buffer.data(target + 710);
    auto *t_711 = buffer.data(target + 711);
    auto *t_712 = buffer.data(target + 712);
    auto *t_713 = buffer.data(target + 713);
    auto *t_714 = buffer.data(target + 714);
    auto *t_715 = buffer.data(target + 715);
    auto *t_716 = buffer.data(target + 716);
    auto *t_717 = buffer.data(target + 717);
    auto *t_718 = buffer.data(target + 718);
    auto *t_719 = buffer.data(target + 719);
    auto *t_720 = buffer.data(target + 720);
    auto *t_721 = buffer.data(target + 721);
    auto *t_722 = buffer.data(target + 722);
    auto *t_723 = buffer.data(target + 723);
    auto *t_724 = buffer.data(target + 724);
    auto *t_725 = buffer.data(target + 725);
    auto *t_726 = buffer.data(target + 726);
    auto *t_727 = buffer.data(target + 727);
    auto *t_728 = buffer.data(target + 728);
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);
    auto *t_731 = buffer.data(target + 731);
    auto *t_732 = buffer.data(target + 732);
    auto *t_733 = buffer.data(target + 733);
    auto *t_734 = buffer.data(target + 734);
    auto *t_735 = buffer.data(target + 735);
    auto *t_736 = buffer.data(target + 736);
    auto *t_737 = buffer.data(target + 737);
    auto *t_738 = buffer.data(target + 738);
    auto *t_739 = buffer.data(target + 739);
    auto *t_740 = buffer.data(target + 740);
    auto *t_741 = buffer.data(target + 741);
    auto *t_742 = buffer.data(target + 742);
    auto *t_743 = buffer.data(target + 743);
    auto *t_744 = buffer.data(target + 744);
    auto *t_745 = buffer.data(target + 745);
    auto *t_746 = buffer.data(target + 746);
    auto *t_747 = buffer.data(target + 747);
    auto *t_748 = buffer.data(target + 748);
    auto *t_749 = buffer.data(target + 749);
    auto *t_750 = buffer.data(target + 750);
    auto *t_751 = buffer.data(target + 751);
    auto *t_752 = buffer.data(target + 752);
    auto *t_753 = buffer.data(target + 753);
    auto *t_754 = buffer.data(target + 754);
    auto *t_755 = buffer.data(target + 755);
    auto *t_756 = buffer.data(target + 756);
    auto *t_757 = buffer.data(target + 757);
    auto *t_758 = buffer.data(target + 758);
    auto *t_759 = buffer.data(target + 759);
    auto *t_760 = buffer.data(target + 760);
    auto *t_761 = buffer.data(target + 761);
    auto *t_762 = buffer.data(target + 762);
    auto *t_763 = buffer.data(target + 763);
    auto *t_764 = buffer.data(target + 764);
    auto *t_765 = buffer.data(target + 765);
    auto *t_766 = buffer.data(target + 766);
    auto *t_767 = buffer.data(target + 767);
    auto *t_768 = buffer.data(target + 768);
    auto *t_769 = buffer.data(target + 769);
    auto *t_770 = buffer.data(target + 770);
    auto *t_771 = buffer.data(target + 771);
    auto *t_772 = buffer.data(target + 772);
    auto *t_773 = buffer.data(target + 773);
    auto *t_774 = buffer.data(target + 774);
    auto *t_775 = buffer.data(target + 775);
    auto *t_776 = buffer.data(target + 776);
    auto *t_777 = buffer.data(target + 777);
    auto *t_778 = buffer.data(target + 778);
    auto *t_779 = buffer.data(target + 779);
    auto *t_780 = buffer.data(target + 780);
    auto *t_781 = buffer.data(target + 781);
    auto *t_782 = buffer.data(target + 782);
    auto *t_783 = buffer.data(target + 783);
    auto *t_784 = buffer.data(target + 784);
    auto *t_785 = buffer.data(target + 785);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kk0_267 = buffer.data(kk0 + 267);
    const auto *kk0_324 = buffer.data(kk0 + 324);
    const auto *kk0_338 = buffer.data(kk0 + 338);
    const auto *kk0_344 = buffer.data(kk0 + 344);
    const auto *kk0_360 = buffer.data(kk0 + 360);
    const auto *kk0_676 = buffer.data(kk0 + 676);
    const auto *kk0_678 = buffer.data(kk0 + 678);
    const auto *kk0_679 = buffer.data(kk0 + 679);
    const auto *kk0_680 = buffer.data(kk0 + 680);
    const auto *kk0_681 = buffer.data(kk0 + 681);
    const auto *kk0_683 = buffer.data(kk0 + 683);
    const auto *kk0_755 = buffer.data(kk0 + 755);
    const auto *kk0_784 = buffer.data(kk0 + 784);

    const auto *kk1_267 = buffer.data(kk1 + 267);
    const auto *kk1_324 = buffer.data(kk1 + 324);
    const auto *kk1_338 = buffer.data(kk1 + 338);
    const auto *kk1_344 = buffer.data(kk1 + 344);
    const auto *kk1_360 = buffer.data(kk1 + 360);
    const auto *kk1_676 = buffer.data(kk1 + 676);
    const auto *kk1_678 = buffer.data(kk1 + 678);
    const auto *kk1_679 = buffer.data(kk1 + 679);
    const auto *kk1_680 = buffer.data(kk1 + 680);
    const auto *kk1_681 = buffer.data(kk1 + 681);
    const auto *kk1_683 = buffer.data(kk1 + 683);
    const auto *kk1_755 = buffer.data(kk1 + 755);
    const auto *kk1_784 = buffer.data(kk1 + 784);

    const auto *li_342 = buffer.data(li + 342);
    const auto *li_346 = buffer.data(li + 346);
    const auto *li_357 = buffer.data(li + 357);
    const auto *li_367 = buffer.data(li + 367);
    const auto *li_370 = buffer.data(li + 370);
    const auto *li_373 = buffer.data(li + 373);
    const auto *li_374 = buffer.data(li + 374);
    const auto *li_378 = buffer.data(li + 378);
    const auto *li_385 = buffer.data(li + 385);
    const auto *li_391 = buffer.data(li + 391);
    const auto *li_392 = buffer.data(li + 392);
    const auto *li_393 = buffer.data(li + 393);
    const auto *li_394 = buffer.data(li + 394);
    const auto *li_395 = buffer.data(li + 395);
    const auto *li_397 = buffer.data(li + 397);
    const auto *li_398 = buffer.data(li + 398);
    const auto *li_400 = buffer.data(li + 400);
    const auto *li_401 = buffer.data(li + 401);
    const auto *li_402 = buffer.data(li + 402);
    const auto *li_404 = buffer.data(li + 404);
    const auto *li_405 = buffer.data(li + 405);
    const auto *li_406 = buffer.data(li + 406);
    const auto *li_413 = buffer.data(li + 413);
    const auto *li_415 = buffer.data(li + 415);
    const auto *li_416 = buffer.data(li + 416);
    const auto *li_417 = buffer.data(li + 417);
    const auto *li_418 = buffer.data(li + 418);
    const auto *li_419 = buffer.data(li + 419);
    const auto *li_420 = buffer.data(li + 420);
    const auto *li_425 = buffer.data(li + 425);
    const auto *li_429 = buffer.data(li + 429);
    const auto *li_434 = buffer.data(li + 434);
    const auto *li_516 = buffer.data(li + 516);
    const auto *li_521 = buffer.data(li + 521);
    const auto *li_522 = buffer.data(li + 522);
    const auto *li_525 = buffer.data(li + 525);
    const auto *li_526 = buffer.data(li + 526);
    const auto *li_527 = buffer.data(li + 527);
    const auto *li_528 = buffer.data(li + 528);
    const auto *li_529 = buffer.data(li + 529);
    const auto *li_530 = buffer.data(li + 530);
    const auto *li_531 = buffer.data(li + 531);
    const auto *li_553 = buffer.data(li + 553);
    const auto *li_554 = buffer.data(li + 554);
    const auto *li_555 = buffer.data(li + 555);
    const auto *li_556 = buffer.data(li + 556);
    const auto *li_557 = buffer.data(li + 557);
    const auto *li_558 = buffer.data(li + 558);
    const auto *li_565 = buffer.data(li + 565);
    const auto *li_569 = buffer.data(li + 569);
    const auto *li_574 = buffer.data(li + 574);
    const auto *li_580 = buffer.data(li + 580);
    const auto *li_581 = buffer.data(li + 581);
    const auto *li_582 = buffer.data(li + 582);
    const auto *li_583 = buffer.data(li + 583);
    const auto *li_584 = buffer.data(li + 584);
    const auto *li_585 = buffer.data(li + 585);
    const auto *li_587 = buffer.data(li + 587);
    const auto *li_591 = buffer.data(li + 591);
    const auto *li_594 = buffer.data(li + 594);
    const auto *li_598 = buffer.data(li + 598);
    const auto *li_603 = buffer.data(li + 603);
    const auto *li_609 = buffer.data(li + 609);
    const auto *li_611 = buffer.data(li + 611);
    const auto *li_612 = buffer.data(li + 612);
    const auto *li_613 = buffer.data(li + 613);
    const auto *li_614 = buffer.data(li + 614);
    const auto *li_615 = buffer.data(li + 615);

    const auto *lk_447 = buffer.data(lk + 447);
    const auto *lk_482 = buffer.data(lk + 482);
    const auto *lk_488 = buffer.data(lk + 488);
    const auto *lk_504 = buffer.data(lk + 504);
    const auto *lk_506 = buffer.data(lk + 506);
    const auto *lk_507 = buffer.data(lk + 507);
    const auto *lk_509 = buffer.data(lk + 509);
    const auto *lk_510 = buffer.data(lk + 510);
    const auto *lk_513 = buffer.data(lk + 513);
    const auto *lk_514 = buffer.data(lk + 514);
    const auto *lk_516 = buffer.data(lk + 516);
    const auto *lk_518 = buffer.data(lk + 518);
    const auto *lk_519 = buffer.data(lk + 519);
    const auto *lk_521 = buffer.data(lk + 521);
    const auto *lk_522 = buffer.data(lk + 522);
    const auto *lk_524 = buffer.data(lk + 524);
    const auto *lk_531 = buffer.data(lk + 531);
    const auto *lk_532 = buffer.data(lk + 532);
    const auto *lk_534 = buffer.data(lk + 534);
    const auto *lk_535 = buffer.data(lk + 535);
    const auto *lk_536 = buffer.data(lk + 536);
    const auto *lk_537 = buffer.data(lk + 537);
    const auto *lk_539 = buffer.data(lk + 539);
    const auto *lk_540 = buffer.data(lk + 540);
    const auto *lk_676 = buffer.data(lk + 676);
    const auto *lk_678 = buffer.data(lk + 678);
    const auto *lk_679 = buffer.data(lk + 679);
    const auto *lk_680 = buffer.data(lk + 680);
    const auto *lk_681 = buffer.data(lk + 681);
    const auto *lk_683 = buffer.data(lk + 683);
    const auto *lk_755 = buffer.data(lk + 755);
    const auto *lk_784 = buffer.data(lk + 784);

    const auto *mh0_390 = buffer.data(mh0 + 390);
    const auto *mh0_395 = buffer.data(mh0 + 395);
    const auto *mh0_396 = buffer.data(mh0 + 396);
    const auto *mh0_420 = buffer.data(mh0 + 420);
    const auto *mh0_421 = buffer.data(mh0 + 421);
    const auto *mh0_423 = buffer.data(mh0 + 423);
    const auto *mh0_425 = buffer.data(mh0 + 425);
    const auto *mh0_426 = buffer.data(mh0 + 426);
    const auto *mh0_428 = buffer.data(mh0 + 428);
    const auto *mh0_429 = buffer.data(mh0 + 429);
    const auto *mh0_434 = buffer.data(mh0 + 434);
    const auto *mh0_435 = buffer.data(mh0 + 435);
    const auto *mh0_437 = buffer.data(mh0 + 437);
    const auto *mh0_438 = buffer.data(mh0 + 438);
    const auto *mh0_439 = buffer.data(mh0 + 439);
    const auto *mh0_440 = buffer.data(mh0 + 440);
    const auto *mh0_441 = buffer.data(mh0 + 441);
    const auto *mh0_443 = buffer.data(mh0 + 443);
    const auto *mh0_444 = buffer.data(mh0 + 444);
    const auto *mh0_446 = buffer.data(mh0 + 446);
    const auto *mh0_447 = buffer.data(mh0 + 447);
    const auto *mh0_448 = buffer.data(mh0 + 448);
    const auto *mh0_450 = buffer.data(mh0 + 450);
    const auto *mh0_451 = buffer.data(mh0 + 451);
    const auto *mh0_456 = buffer.data(mh0 + 456);

    const auto *mh1_390 = buffer.data(mh1 + 390);
    const auto *mh1_395 = buffer.data(mh1 + 395);
    const auto *mh1_396 = buffer.data(mh1 + 396);
    const auto *mh1_420 = buffer.data(mh1 + 420);
    const auto *mh1_421 = buffer.data(mh1 + 421);
    const auto *mh1_423 = buffer.data(mh1 + 423);
    const auto *mh1_425 = buffer.data(mh1 + 425);
    const auto *mh1_426 = buffer.data(mh1 + 426);
    const auto *mh1_428 = buffer.data(mh1 + 428);
    const auto *mh1_429 = buffer.data(mh1 + 429);
    const auto *mh1_434 = buffer.data(mh1 + 434);
    const auto *mh1_435 = buffer.data(mh1 + 435);
    const auto *mh1_437 = buffer.data(mh1 + 437);
    const auto *mh1_438 = buffer.data(mh1 + 438);
    const auto *mh1_439 = buffer.data(mh1 + 439);
    const auto *mh1_440 = buffer.data(mh1 + 440);
    const auto *mh1_441 = buffer.data(mh1 + 441);
    const auto *mh1_443 = buffer.data(mh1 + 443);
    const auto *mh1_444 = buffer.data(mh1 + 444);
    const auto *mh1_446 = buffer.data(mh1 + 446);
    const auto *mh1_447 = buffer.data(mh1 + 447);
    const auto *mh1_448 = buffer.data(mh1 + 448);
    const auto *mh1_450 = buffer.data(mh1 + 450);
    const auto *mh1_451 = buffer.data(mh1 + 451);
    const auto *mh1_456 = buffer.data(mh1 + 456);

    const auto *mi_510 = buffer.data(mi + 510);
    const auto *mi_513 = buffer.data(mi + 513);
    const auto *mi_514 = buffer.data(mi + 514);
    const auto *mi_516 = buffer.data(mi + 516);
    const auto *mi_518 = buffer.data(mi + 518);
    const auto *mi_521 = buffer.data(mi + 521);
    const auto *mi_522 = buffer.data(mi + 522);
    const auto *mi_525 = buffer.data(mi + 525);
    const auto *mi_526 = buffer.data(mi + 526);
    const auto *mi_527 = buffer.data(mi + 527);
    const auto *mi_528 = buffer.data(mi + 528);
    const auto *mi_529 = buffer.data(mi + 529);
    const auto *mi_530 = buffer.data(mi + 530);
    const auto *mi_531 = buffer.data(mi + 531);
    const auto *mi_532 = buffer.data(mi + 532);
    const auto *mi_534 = buffer.data(mi + 534);
    const auto *mi_535 = buffer.data(mi + 535);
    const auto *mi_537 = buffer.data(mi + 537);
    const auto *mi_538 = buffer.data(mi + 538);
    const auto *mi_541 = buffer.data(mi + 541);
    const auto *mi_542 = buffer.data(mi + 542);
    const auto *mi_546 = buffer.data(mi + 546);
    const auto *mi_553 = buffer.data(mi + 553);
    const auto *mi_554 = buffer.data(mi + 554);
    const auto *mi_555 = buffer.data(mi + 555);
    const auto *mi_556 = buffer.data(mi + 556);
    const auto *mi_557 = buffer.data(mi + 557);
    const auto *mi_558 = buffer.data(mi + 558);
    const auto *mi_559 = buffer.data(mi + 559);
    const auto *mi_560 = buffer.data(mi + 560);
    const auto *mi_561 = buffer.data(mi + 561);
    const auto *mi_562 = buffer.data(mi + 562);
    const auto *mi_563 = buffer.data(mi + 563);
    const auto *mi_565 = buffer.data(mi + 565);
    const auto *mi_566 = buffer.data(mi + 566);
    const auto *mi_568 = buffer.data(mi + 568);
    const auto *mi_569 = buffer.data(mi + 569);
    const auto *mi_570 = buffer.data(mi + 570);
    const auto *mi_572 = buffer.data(mi + 572);
    const auto *mi_573 = buffer.data(mi + 573);
    const auto *mi_574 = buffer.data(mi + 574);
    const auto *mi_580 = buffer.data(mi + 580);
    const auto *mi_581 = buffer.data(mi + 581);
    const auto *mi_582 = buffer.data(mi + 582);
    const auto *mi_583 = buffer.data(mi + 583);
    const auto *mi_584 = buffer.data(mi + 584);
    const auto *mi_585 = buffer.data(mi + 585);
    const auto *mi_586 = buffer.data(mi + 586);
    const auto *mi_587 = buffer.data(mi + 587);
    const auto *mi_588 = buffer.data(mi + 588);
    const auto *mi_589 = buffer.data(mi + 589);
    const auto *mi_590 = buffer.data(mi + 590);
    const auto *mi_591 = buffer.data(mi + 591);
    const auto *mi_593 = buffer.data(mi + 593);
    const auto *mi_594 = buffer.data(mi + 594);
    const auto *mi_595 = buffer.data(mi + 595);
    const auto *mi_597 = buffer.data(mi + 597);
    const auto *mi_598 = buffer.data(mi + 598);
    const auto *mi_599 = buffer.data(mi + 599);
    const auto *mi_600 = buffer.data(mi + 600);
    const auto *mi_602 = buffer.data(mi + 602);
    const auto *mi_603 = buffer.data(mi + 603);
    const auto *mi_609 = buffer.data(mi + 609);
    const auto *mi_611 = buffer.data(mi + 611);
    const auto *mi_612 = buffer.data(mi + 612);
    const auto *mi_613 = buffer.data(mi + 613);
    const auto *mi_614 = buffer.data(mi + 614);
    const auto *mi_615 = buffer.data(mi + 615);

#pragma omp simd aligned(t_659, t_660, t_661, pb_x, pb_y, pb_z, li_342, li_373, li_516, \
                         mh0_390, mh1_390, mi_510, mi_513, mi_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_13 * li_342[k]
                   + pb_z[k] * mi_510[k];

        t_660[k] = f_14 * li_516[k]
                   + f_5 * mh0_390[k]
                   - f_6 * mh1_390[k]
                   + pb_x[k] * mi_516[k];

        t_661[k] = f_12 * li_373[k]
                   + pb_y[k] * mi_513[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pa_y, pa_z, pb_z, kk0_267, kk0_338, kk1_267, \
                         kk1_338, li_346, lk_447, lk_482, mi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_18 * kk0_338[k]
                   - f_19 * kk1_338[k]
                   + pa_y[k] * lk_482[k];

        t_663[k] = f_22 * kk0_267[k]
                   - f_23 * kk1_267[k]
                   + pa_z[k] * lk_447[k];

        t_664[k] = f_13 * li_346[k]
                   + pb_z[k] * mi_514[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, pb_x, pb_y, li_378, li_521, li_522, mh0_395, \
                         mh0_396, mh1_395, mh1_396, mi_518, mi_521, \
                         mi_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_14 * li_521[k]
                   + f_3 * mh0_395[k]
                   - f_4 * mh1_395[k]
                   + pb_x[k] * mi_521[k];

        t_666[k] = f_14 * li_522[k]
                   + f_3 * mh0_396[k]
                   - f_4 * mh1_396[k]
                   + pb_x[k] * mi_522[k];

        t_667[k] = f_12 * li_378[k]
                   + pb_y[k] * mi_518[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, pa_y, pb_x, kk0_344, kk1_344, li_525, \
                         li_526, li_527, lk_488, mi_525, mi_526, \
                         mi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_18 * kk0_344[k]
                   - f_19 * kk1_344[k]
                   + pa_y[k] * lk_488[k];

        t_669[k] = f_14 * li_525[k]
                   + pb_x[k] * mi_525[k];

        t_670[k] = f_14 * li_526[k]
                   + pb_x[k] * mi_526[k];

        t_671[k] = f_14 * li_527[k]
                   + pb_x[k] * mi_527[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pb_x, li_528, li_529, li_530, li_531, \
                         mi_528, mi_529, mi_530, mi_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_14 * li_528[k]
                   + pb_x[k] * mi_528[k];

        t_673[k] = f_14 * li_529[k]
                   + pb_x[k] * mi_529[k];

        t_674[k] = f_14 * li_530[k]
                   + pb_x[k] * mi_530[k];

        t_675[k] = f_14 * li_531[k]
                   + pb_x[k] * mi_531[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pa_x, pb_z, kk0_676, kk0_678, kk1_676, kk1_678, \
                         li_357, lk_676, lk_678, mi_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_27 * kk0_676[k]
                   - f_28 * kk1_676[k]
                   + pa_x[k] * lk_676[k];

        t_677[k] = f_13 * li_357[k]
                   + pb_z[k] * mi_525[k];

        t_678[k] = f_27 * kk0_678[k]
                   - f_28 * kk1_678[k]
                   + pa_x[k] * lk_678[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, pa_x, kk0_679, kk0_680, kk0_681, kk1_679, \
                         kk1_680, kk1_681, lk_679, lk_680, lk_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_27 * kk0_679[k]
                   - f_28 * kk1_679[k]
                   + pa_x[k] * lk_679[k];

        t_680[k] = f_27 * kk0_680[k]
                   - f_28 * kk1_680[k]
                   + pa_x[k] * lk_680[k];

        t_681[k] = f_27 * kk0_681[k]
                   - f_28 * kk1_681[k]
                   + pa_x[k] * lk_681[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, pa_x, pa_y, pb_y, kk0_683, kk1_683, \
                         li_391, li_392, lk_504, lk_683, mi_531, \
                         mi_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_12 * li_391[k]
                   + pb_y[k] * mi_531[k];

        t_683[k] = f_27 * kk0_683[k]
                   - f_28 * kk1_683[k]
                   + pa_x[k] * lk_683[k];

        t_684[k] = pa_y[k] * lk_504[k];

        t_685[k] = f_11 * li_392[k]
                   + pb_y[k] * mi_532[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, t_690, pa_y, pb_y, li_393, li_394, \
                         li_395, lk_506, lk_507, lk_509, lk_510, \
                         mi_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = pa_y[k] * lk_506[k];

        t_687[k] = f_12 * li_393[k]
                   + pa_y[k] * lk_507[k];

        t_688[k] = f_11 * li_394[k]
                   + pb_y[k] * mi_534[k];

        t_689[k] = pa_y[k] * lk_509[k];

        t_690[k] = f_13 * li_395[k]
                   + pa_y[k] * lk_510[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, t_694, pa_y, pb_y, pb_z, li_367, li_397, li_398, \
                         lk_513, lk_514, mi_535, mi_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_14 * li_367[k]
                   + pb_z[k] * mi_535[k];

        t_692[k] = f_11 * li_397[k]
                   + pb_y[k] * mi_537[k];

        t_693[k] = pa_y[k] * lk_513[k];

        t_694[k] = f_14 * li_398[k]
                   + pa_y[k] * lk_514[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, pa_y, pb_y, pb_z, li_370, li_400, li_401, \
                         lk_516, lk_518, mi_538, mi_541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_14 * li_370[k]
                   + pb_z[k] * mi_538[k];

        t_696[k] = f_12 * li_400[k]
                   + pa_y[k] * lk_516[k];

        t_697[k] = f_11 * li_401[k]
                   + pb_y[k] * mi_541[k];

        t_698[k] = pa_y[k] * lk_518[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, pa_y, pb_z, li_374, li_402, li_404, \
                         li_405, lk_519, lk_521, lk_522, mi_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_15 * li_402[k]
                   + pa_y[k] * lk_519[k];

        t_700[k] = f_14 * li_374[k]
                   + pb_z[k] * mi_542[k];

        t_701[k] = f_13 * li_404[k]
                   + pa_y[k] * lk_521[k];

        t_702[k] = f_12 * li_405[k]
                   + pa_y[k] * lk_522[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, t_706, pa_y, pb_x, pb_y, li_406, li_553, li_554, \
                         lk_524, mi_546, mi_553, mi_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_11 * li_406[k]
                   + pb_y[k] * mi_546[k];

        t_704[k] = pa_y[k] * lk_524[k];

        t_705[k] = f_14 * li_553[k]
                   + pb_x[k] * mi_553[k];

        t_706[k] = f_14 * li_554[k]
                   + pb_x[k] * mi_554[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, t_710, t_711, pa_y, pb_x, li_555, li_556, \
                         li_557, li_558, lk_531, mi_555, mi_556, mi_557, \
                         mi_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = f_14 * li_555[k]
                   + pb_x[k] * mi_555[k];

        t_708[k] = f_14 * li_556[k]
                   + pb_x[k] * mi_556[k];

        t_709[k] = f_14 * li_557[k]
                   + pb_x[k] * mi_557[k];

        t_710[k] = f_14 * li_558[k]
                   + pb_x[k] * mi_558[k];

        t_711[k] = pa_y[k] * lk_531[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pa_y, pb_z, li_385, li_413, li_415, \
                         li_416, lk_532, lk_534, lk_535, mi_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_17 * li_413[k]
                   + pa_y[k] * lk_532[k];

        t_713[k] = f_14 * li_385[k]
                   + pb_z[k] * mi_553[k];

        t_714[k] = f_15 * li_415[k]
                   + pa_y[k] * lk_534[k];

        t_715[k] = f_14 * li_416[k]
                   + pa_y[k] * lk_535[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pa_y, pb_y, li_417, li_418, li_419, \
                         lk_536, lk_537, lk_539, mi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_13 * li_417[k]
                   + pa_y[k] * lk_536[k];

        t_717[k] = f_12 * li_418[k]
                   + pa_y[k] * lk_537[k];

        t_718[k] = f_11 * li_419[k]
                   + pb_y[k] * mi_559[k];

        t_719[k] = pa_y[k] * lk_539[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pa_z, pb_y, pb_z, kk0_324, kk1_324, \
                         li_392, lk_504, mh0_420, mh1_420, mi_560, \
                         mi_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_29 * kk0_324[k]
                   - f_30 * kk1_324[k]
                   + pa_z[k] * lk_504[k];

        t_721[k] = pb_y[k] * mi_560[k];

        t_722[k] = f_15 * li_392[k]
                   + pb_z[k] * mi_560[k];

        t_723[k] = f_3 * mh0_420[k]
                   - f_4 * mh1_420[k]
                   + pb_y[k] * mi_561[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, pb_x, pb_y, pb_z, li_395, li_565, \
                         mh0_421, mh0_425, mh1_421, mh1_425, mi_562, mi_563, \
                         mi_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = pb_y[k] * mi_562[k];

        t_725[k] = f_14 * li_565[k]
                   + f_9 * mh0_425[k]
                   - f_10 * mh1_425[k]
                   + pb_x[k] * mi_565[k];

        t_726[k] = f_5 * mh0_421[k]
                   - f_6 * mh1_421[k]
                   + pb_y[k] * mi_563[k];

        t_727[k] = f_15 * li_395[k]
                   + pb_z[k] * mi_563[k];
    }

#pragma omp simd aligned(t_728, t_729, t_730, t_731, pb_x, pb_y, pb_z, li_398, li_569, \
                         mh0_423, mh0_429, mh1_423, mh1_429, mi_565, mi_566, \
                         mi_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_728[k] = pb_y[k] * mi_565[k];

        t_729[k] = f_14 * li_569[k]
                   + f_7 * mh0_429[k]
                   - f_8 * mh1_429[k]
                   + pb_x[k] * mi_569[k];

        t_730[k] = f_7 * mh0_423[k]
                   - f_8 * mh1_423[k]
                   + pb_y[k] * mi_566[k];

        t_731[k] = f_15 * li_398[k]
                   + pb_z[k] * mi_566[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pb_x, pb_y, li_574, mh0_425, mh0_434, mh1_425, \
                         mh1_434, mi_568, mi_569, mi_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_3 * mh0_425[k]
                   - f_4 * mh1_425[k]
                   + pb_y[k] * mi_568[k];

        t_733[k] = pb_y[k] * mi_569[k];

        t_734[k] = f_14 * li_574[k]
                   + f_5 * mh0_434[k]
                   - f_6 * mh1_434[k]
                   + pb_x[k] * mi_574[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, pb_y, pb_z, li_402, mh0_426, mh0_428, \
                         mh0_429, mh1_426, mh1_428, mh1_429, mi_570, mi_572, \
                         mi_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_9 * mh0_426[k]
                   - f_10 * mh1_426[k]
                   + pb_y[k] * mi_570[k];

        t_736[k] = f_15 * li_402[k]
                   + pb_z[k] * mi_570[k];

        t_737[k] = f_5 * mh0_428[k]
                   - f_6 * mh1_428[k]
                   + pb_y[k] * mi_572[k];

        t_738[k] = f_3 * mh0_429[k]
                   - f_4 * mh1_429[k]
                   + pb_y[k] * mi_573[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, pb_x, pb_y, li_580, li_581, li_582, \
                         mh0_440, mh1_440, mi_574, mi_580, mi_581, \
                         mi_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = pb_y[k] * mi_574[k];

        t_740[k] = f_14 * li_580[k]
                   + f_3 * mh0_440[k]
                   - f_4 * mh1_440[k]
                   + pb_x[k] * mi_580[k];

        t_741[k] = f_14 * li_581[k]
                   + pb_x[k] * mi_581[k];

        t_742[k] = f_14 * li_582[k]
                   + pb_x[k] * mi_582[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, t_747, pb_x, pb_y, li_583, li_584, \
                         li_585, li_587, mi_580, mi_583, mi_584, mi_585, \
                         mi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = f_14 * li_583[k]
                   + pb_x[k] * mi_583[k];

        t_744[k] = f_14 * li_584[k]
                   + pb_x[k] * mi_584[k];

        t_745[k] = f_14 * li_585[k]
                   + pb_x[k] * mi_585[k];

        t_746[k] = pb_y[k] * mi_580[k];

        t_747[k] = f_14 * li_587[k]
                   + pb_x[k] * mi_587[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pb_y, pb_z, li_413, mh0_435, mh0_437, \
                         mh0_438, mh1_435, mh1_437, mh1_438, mi_581, mi_583, \
                         mi_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_1 * mh0_435[k]
                   - f_2 * mh1_435[k]
                   + pb_y[k] * mi_581[k];

        t_749[k] = f_15 * li_413[k]
                   + pb_z[k] * mi_581[k];

        t_750[k] = f_9 * mh0_437[k]
                   - f_10 * mh1_437[k]
                   + pb_y[k] * mi_583[k];

        t_751[k] = f_7 * mh0_438[k]
                   - f_8 * mh1_438[k]
                   + pb_y[k] * mi_584[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pa_x, pb_y, kk0_755, kk1_755, lk_755, \
                         mh0_439, mh0_440, mh1_439, mh1_440, mi_585, mi_586, \
                         mi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_5 * mh0_439[k]
                   - f_6 * mh1_439[k]
                   + pb_y[k] * mi_585[k];

        t_753[k] = f_3 * mh0_440[k]
                   - f_4 * mh1_440[k]
                   + pb_y[k] * mi_586[k];

        t_754[k] = pb_y[k] * mi_587[k];

        t_755[k] = f_27 * kk0_755[k]
                   - f_28 * kk1_755[k]
                   + pa_x[k] * lk_755[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, pa_y, pb_y, pb_z, kk0_360, kk1_360, li_420, \
                         lk_540, mi_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_25 * kk0_360[k]
                   - f_26 * kk1_360[k]
                   + pa_y[k] * lk_540[k];

        t_757[k] = f_24 * li_420[k]
                   + pb_y[k] * mi_588[k];

        t_758[k] = pb_z[k] * mi_588[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pb_x, pb_z, li_591, mh0_441, mh0_444, mh1_441, \
                         mh1_444, mi_589, mi_590, mi_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_13 * li_591[k]
                   + f_9 * mh0_444[k]
                   - f_10 * mh1_444[k]
                   + pb_x[k] * mi_591[k];

        t_760[k] = pb_z[k] * mi_589[k];

        t_761[k] = f_3 * mh0_441[k]
                   - f_4 * mh1_441[k]
                   + pb_z[k] * mi_590[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, pb_x, pb_y, pb_z, li_425, li_594, \
                         mh0_443, mh0_447, mh1_443, mh1_447, mi_591, mi_593, \
                         mi_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_13 * li_594[k]
                   + f_7 * mh0_447[k]
                   - f_8 * mh1_447[k]
                   + pb_x[k] * mi_594[k];

        t_763[k] = pb_z[k] * mi_591[k];

        t_764[k] = f_24 * li_425[k]
                   + pb_y[k] * mi_593[k];

        t_765[k] = f_5 * mh0_443[k]
                   - f_6 * mh1_443[k]
                   + pb_z[k] * mi_593[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, pb_x, pb_z, li_598, mh0_444, mh0_451, mh1_444, \
                         mh1_451, mi_594, mi_595, mi_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = f_13 * li_598[k]
                   + f_5 * mh0_451[k]
                   - f_6 * mh1_451[k]
                   + pb_x[k] * mi_598[k];

        t_767[k] = pb_z[k] * mi_594[k];

        t_768[k] = f_3 * mh0_444[k]
                   - f_4 * mh1_444[k]
                   + pb_z[k] * mi_595[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, t_772, pb_x, pb_y, pb_z, li_429, li_603, \
                         mh0_446, mh0_456, mh1_446, mh1_456, mi_597, mi_598, \
                         mi_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = f_24 * li_429[k]
                   + pb_y[k] * mi_597[k];

        t_770[k] = f_7 * mh0_446[k]
                   - f_8 * mh1_446[k]
                   + pb_z[k] * mi_597[k];

        t_771[k] = f_13 * li_603[k]
                   + f_3 * mh0_456[k]
                   - f_4 * mh1_456[k]
                   + pb_x[k] * mi_603[k];

        t_772[k] = pb_z[k] * mi_598[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, t_776, pb_y, pb_z, li_434, mh0_447, mh0_448, \
                         mh0_450, mh1_447, mh1_448, mh1_450, mi_599, mi_600, \
                         mi_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_3 * mh0_447[k]
                   - f_4 * mh1_447[k]
                   + pb_z[k] * mi_599[k];

        t_774[k] = f_5 * mh0_448[k]
                   - f_6 * mh1_448[k]
                   + pb_z[k] * mi_600[k];

        t_775[k] = f_24 * li_434[k]
                   + pb_y[k] * mi_602[k];

        t_776[k] = f_9 * mh0_450[k]
                   - f_10 * mh1_450[k]
                   + pb_z[k] * mi_602[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, t_780, t_781, pb_x, pb_z, li_609, li_611, \
                         li_612, li_613, mi_603, mi_609, mi_611, mi_612, \
                         mi_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_13 * li_609[k]
                   + pb_x[k] * mi_609[k];

        t_778[k] = pb_z[k] * mi_603[k];

        t_779[k] = f_13 * li_611[k]
                   + pb_x[k] * mi_611[k];

        t_780[k] = f_13 * li_612[k]
                   + pb_x[k] * mi_612[k];

        t_781[k] = f_13 * li_613[k]
                   + pb_x[k] * mi_613[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, t_785, pa_x, pb_x, pb_z, kk0_784, kk1_784, \
                         li_614, li_615, lk_784, mi_609, mi_614, \
                         mi_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_13 * li_614[k]
                   + pb_x[k] * mi_614[k];

        t_783[k] = f_13 * li_615[k]
                   + pb_x[k] * mi_615[k];

        t_784[k] = f_22 * kk0_784[k]
                   - f_23 * kk1_784[k]
                   + pa_x[k] * lk_784[k];

        t_785[k] = pb_z[k] * mi_609[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece6(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kk0,
                                            const size_t kk1, const size_t li, const size_t lk,
                                            const size_t mh0, const size_t mh1, const size_t mi,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_17 = 3.5 / p;
    const auto f_18 = 0.5 / alpha;
    const auto f_19 = 0.5 * beta / (alpha * p);
    const auto f_22 = 1.0 / alpha;
    const auto f_23 = beta / (alpha * p);
    const auto f_24 = 3.0 / p;
    const auto f_27 = 1.5 / alpha;
    const auto f_28 = 1.5 * beta / (alpha * p);

    auto *t_786 = buffer.data(target + 786);
    auto *t_787 = buffer.data(target + 787);
    auto *t_788 = buffer.data(target + 788);
    auto *t_789 = buffer.data(target + 789);
    auto *t_790 = buffer.data(target + 790);
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);
    auto *t_794 = buffer.data(target + 794);
    auto *t_795 = buffer.data(target + 795);
    auto *t_796 = buffer.data(target + 796);
    auto *t_797 = buffer.data(target + 797);
    auto *t_798 = buffer.data(target + 798);
    auto *t_799 = buffer.data(target + 799);
    auto *t_800 = buffer.data(target + 800);
    auto *t_801 = buffer.data(target + 801);
    auto *t_802 = buffer.data(target + 802);
    auto *t_803 = buffer.data(target + 803);
    auto *t_804 = buffer.data(target + 804);
    auto *t_805 = buffer.data(target + 805);
    auto *t_806 = buffer.data(target + 806);
    auto *t_807 = buffer.data(target + 807);
    auto *t_808 = buffer.data(target + 808);
    auto *t_809 = buffer.data(target + 809);
    auto *t_810 = buffer.data(target + 810);
    auto *t_811 = buffer.data(target + 811);
    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);
    auto *t_815 = buffer.data(target + 815);
    auto *t_816 = buffer.data(target + 816);
    auto *t_817 = buffer.data(target + 817);
    auto *t_818 = buffer.data(target + 818);
    auto *t_819 = buffer.data(target + 819);
    auto *t_820 = buffer.data(target + 820);
    auto *t_821 = buffer.data(target + 821);
    auto *t_822 = buffer.data(target + 822);
    auto *t_823 = buffer.data(target + 823);
    auto *t_824 = buffer.data(target + 824);
    auto *t_825 = buffer.data(target + 825);
    auto *t_826 = buffer.data(target + 826);
    auto *t_827 = buffer.data(target + 827);
    auto *t_828 = buffer.data(target + 828);
    auto *t_829 = buffer.data(target + 829);
    auto *t_830 = buffer.data(target + 830);
    auto *t_831 = buffer.data(target + 831);
    auto *t_832 = buffer.data(target + 832);
    auto *t_833 = buffer.data(target + 833);
    auto *t_834 = buffer.data(target + 834);
    auto *t_835 = buffer.data(target + 835);
    auto *t_836 = buffer.data(target + 836);
    auto *t_837 = buffer.data(target + 837);
    auto *t_838 = buffer.data(target + 838);
    auto *t_839 = buffer.data(target + 839);
    auto *t_840 = buffer.data(target + 840);
    auto *t_841 = buffer.data(target + 841);
    auto *t_842 = buffer.data(target + 842);
    auto *t_843 = buffer.data(target + 843);
    auto *t_844 = buffer.data(target + 844);
    auto *t_845 = buffer.data(target + 845);
    auto *t_846 = buffer.data(target + 846);
    auto *t_847 = buffer.data(target + 847);
    auto *t_848 = buffer.data(target + 848);
    auto *t_849 = buffer.data(target + 849);
    auto *t_850 = buffer.data(target + 850);
    auto *t_851 = buffer.data(target + 851);
    auto *t_852 = buffer.data(target + 852);
    auto *t_853 = buffer.data(target + 853);
    auto *t_854 = buffer.data(target + 854);
    auto *t_855 = buffer.data(target + 855);
    auto *t_856 = buffer.data(target + 856);
    auto *t_857 = buffer.data(target + 857);
    auto *t_858 = buffer.data(target + 858);
    auto *t_859 = buffer.data(target + 859);
    auto *t_860 = buffer.data(target + 860);
    auto *t_861 = buffer.data(target + 861);
    auto *t_862 = buffer.data(target + 862);
    auto *t_863 = buffer.data(target + 863);
    auto *t_864 = buffer.data(target + 864);
    auto *t_865 = buffer.data(target + 865);
    auto *t_866 = buffer.data(target + 866);
    auto *t_867 = buffer.data(target + 867);
    auto *t_868 = buffer.data(target + 868);
    auto *t_869 = buffer.data(target + 869);
    auto *t_870 = buffer.data(target + 870);
    auto *t_871 = buffer.data(target + 871);
    auto *t_872 = buffer.data(target + 872);
    auto *t_873 = buffer.data(target + 873);
    auto *t_874 = buffer.data(target + 874);
    auto *t_875 = buffer.data(target + 875);
    auto *t_876 = buffer.data(target + 876);
    auto *t_877 = buffer.data(target + 877);
    auto *t_878 = buffer.data(target + 878);
    auto *t_879 = buffer.data(target + 879);
    auto *t_880 = buffer.data(target + 880);
    auto *t_881 = buffer.data(target + 881);
    auto *t_882 = buffer.data(target + 882);
    auto *t_883 = buffer.data(target + 883);
    auto *t_884 = buffer.data(target + 884);
    auto *t_885 = buffer.data(target + 885);
    auto *t_886 = buffer.data(target + 886);
    auto *t_887 = buffer.data(target + 887);
    auto *t_888 = buffer.data(target + 888);
    auto *t_889 = buffer.data(target + 889);
    auto *t_890 = buffer.data(target + 890);
    auto *t_891 = buffer.data(target + 891);
    auto *t_892 = buffer.data(target + 892);
    auto *t_893 = buffer.data(target + 893);
    auto *t_894 = buffer.data(target + 894);
    auto *t_895 = buffer.data(target + 895);
    auto *t_896 = buffer.data(target + 896);
    auto *t_897 = buffer.data(target + 897);
    auto *t_898 = buffer.data(target + 898);
    auto *t_899 = buffer.data(target + 899);
    auto *t_900 = buffer.data(target + 900);
    auto *t_901 = buffer.data(target + 901);
    auto *t_902 = buffer.data(target + 902);
    auto *t_903 = buffer.data(target + 903);
    auto *t_904 = buffer.data(target + 904);
    auto *t_905 = buffer.data(target + 905);
    auto *t_906 = buffer.data(target + 906);
    auto *t_907 = buffer.data(target + 907);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kk0_363 = buffer.data(kk0 + 363);
    const auto *kk0_366 = buffer.data(kk0 + 366);
    const auto *kk0_370 = buffer.data(kk0 + 370);
    const auto *kk0_375 = buffer.data(kk0 + 375);
    const auto *kk0_399 = buffer.data(kk0 + 399);
    const auto *kk0_402 = buffer.data(kk0 + 402);
    const auto *kk0_406 = buffer.data(kk0 + 406);
    const auto *kk0_411 = buffer.data(kk0 + 411);
    const auto *kk0_432 = buffer.data(kk0 + 432);
    const auto *kk0_435 = buffer.data(kk0 + 435);
    const auto *kk0_437 = buffer.data(kk0 + 437);
    const auto *kk0_438 = buffer.data(kk0 + 438);
    const auto *kk0_441 = buffer.data(kk0 + 441);
    const auto *kk0_446 = buffer.data(kk0 + 446);
    const auto *kk0_452 = buffer.data(kk0 + 452);
    const auto *kk0_468 = buffer.data(kk0 + 468);
    const auto *kk0_473 = buffer.data(kk0 + 473);
    const auto *kk0_477 = buffer.data(kk0 + 477);
    const auto *kk0_482 = buffer.data(kk0 + 482);
    const auto *kk0_488 = buffer.data(kk0 + 488);
    const auto *kk0_504 = buffer.data(kk0 + 504);
    const auto *kk0_509 = buffer.data(kk0 + 509);
    const auto *kk0_856 = buffer.data(kk0 + 856);
    const auto *kk0_858 = buffer.data(kk0 + 858);
    const auto *kk0_859 = buffer.data(kk0 + 859);
    const auto *kk0_860 = buffer.data(kk0 + 860);
    const auto *kk0_861 = buffer.data(kk0 + 861);
    const auto *kk0_863 = buffer.data(kk0 + 863);
    const auto *kk0_892 = buffer.data(kk0 + 892);
    const auto *kk0_894 = buffer.data(kk0 + 894);
    const auto *kk0_895 = buffer.data(kk0 + 895);
    const auto *kk0_896 = buffer.data(kk0 + 896);
    const auto *kk0_897 = buffer.data(kk0 + 897);
    const auto *kk0_899 = buffer.data(kk0 + 899);

    const auto *kk1_363 = buffer.data(kk1 + 363);
    const auto *kk1_366 = buffer.data(kk1 + 366);
    const auto *kk1_370 = buffer.data(kk1 + 370);
    const auto *kk1_375 = buffer.data(kk1 + 375);
    const auto *kk1_399 = buffer.data(kk1 + 399);
    const auto *kk1_402 = buffer.data(kk1 + 402);
    const auto *kk1_406 = buffer.data(kk1 + 406);
    const auto *kk1_411 = buffer.data(kk1 + 411);
    const auto *kk1_432 = buffer.data(kk1 + 432);
    const auto *kk1_435 = buffer.data(kk1 + 435);
    const auto *kk1_437 = buffer.data(kk1 + 437);
    const auto *kk1_438 = buffer.data(kk1 + 438);
    const auto *kk1_441 = buffer.data(kk1 + 441);
    const auto *kk1_446 = buffer.data(kk1 + 446);
    const auto *kk1_452 = buffer.data(kk1 + 452);
    const auto *kk1_468 = buffer.data(kk1 + 468);
    const auto *kk1_473 = buffer.data(kk1 + 473);
    const auto *kk1_477 = buffer.data(kk1 + 477);
    const auto *kk1_482 = buffer.data(kk1 + 482);
    const auto *kk1_488 = buffer.data(kk1 + 488);
    const auto *kk1_504 = buffer.data(kk1 + 504);
    const auto *kk1_509 = buffer.data(kk1 + 509);
    const auto *kk1_856 = buffer.data(kk1 + 856);
    const auto *kk1_858 = buffer.data(kk1 + 858);
    const auto *kk1_859 = buffer.data(kk1 + 859);
    const auto *kk1_860 = buffer.data(kk1 + 860);
    const auto *kk1_861 = buffer.data(kk1 + 861);
    const auto *kk1_863 = buffer.data(kk1 + 863);
    const auto *kk1_892 = buffer.data(kk1 + 892);
    const auto *kk1_894 = buffer.data(kk1 + 894);
    const auto *kk1_895 = buffer.data(kk1 + 895);
    const auto *kk1_896 = buffer.data(kk1 + 896);
    const auto *kk1_897 = buffer.data(kk1 + 897);
    const auto *kk1_899 = buffer.data(kk1 + 899);

    const auto *li_420 = buffer.data(li + 420);
    const auto *li_422 = buffer.data(li + 422);
    const auto *li_423 = buffer.data(li + 423);
    const auto *li_425 = buffer.data(li + 425);
    const auto *li_426 = buffer.data(li + 426);
    const auto *li_427 = buffer.data(li + 427);
    const auto *li_429 = buffer.data(li + 429);
    const auto *li_430 = buffer.data(li + 430);
    const auto *li_431 = buffer.data(li + 431);
    const auto *li_432 = buffer.data(li + 432);
    const auto *li_434 = buffer.data(li + 434);
    const auto *li_441 = buffer.data(li + 441);
    const auto *li_442 = buffer.data(li + 442);
    const auto *li_443 = buffer.data(li + 443);
    const auto *li_444 = buffer.data(li + 444);
    const auto *li_445 = buffer.data(li + 445);
    const auto *li_447 = buffer.data(li + 447);
    const auto *li_448 = buffer.data(li + 448);
    const auto *li_450 = buffer.data(li + 450);
    const auto *li_451 = buffer.data(li + 451);
    const auto *li_453 = buffer.data(li + 453);
    const auto *li_454 = buffer.data(li + 454);
    const auto *li_457 = buffer.data(li + 457);
    const auto *li_458 = buffer.data(li + 458);
    const auto *li_462 = buffer.data(li + 462);
    const auto *li_469 = buffer.data(li + 469);
    const auto *li_475 = buffer.data(li + 475);
    const auto *li_476 = buffer.data(li + 476);
    const auto *li_478 = buffer.data(li + 478);
    const auto *li_479 = buffer.data(li + 479);
    const auto *li_481 = buffer.data(li + 481);
    const auto *li_482 = buffer.data(li + 482);
    const auto *li_485 = buffer.data(li + 485);
    const auto *li_486 = buffer.data(li + 486);
    const auto *li_490 = buffer.data(li + 490);
    const auto *li_497 = buffer.data(li + 497);
    const auto *li_503 = buffer.data(li + 503);
    const auto *li_504 = buffer.data(li + 504);
    const auto *li_506 = buffer.data(li + 506);
    const auto *li_507 = buffer.data(li + 507);
    const auto *li_509 = buffer.data(li + 509);
    const auto *li_513 = buffer.data(li + 513);
    const auto *li_518 = buffer.data(li + 518);
    const auto *li_531 = buffer.data(li + 531);
    const auto *li_532 = buffer.data(li + 532);
    const auto *li_534 = buffer.data(li + 534);
    const auto *li_638 = buffer.data(li + 638);
    const auto *li_639 = buffer.data(li + 639);
    const auto *li_640 = buffer.data(li + 640);
    const auto *li_641 = buffer.data(li + 641);
    const auto *li_642 = buffer.data(li + 642);
    const auto *li_643 = buffer.data(li + 643);
    const auto *li_656 = buffer.data(li + 656);
    const auto *li_661 = buffer.data(li + 661);
    const auto *li_662 = buffer.data(li + 662);
    const auto *li_665 = buffer.data(li + 665);
    const auto *li_666 = buffer.data(li + 666);
    const auto *li_667 = buffer.data(li + 667);
    const auto *li_668 = buffer.data(li + 668);
    const auto *li_669 = buffer.data(li + 669);
    const auto *li_670 = buffer.data(li + 670);
    const auto *li_671 = buffer.data(li + 671);
    const auto *li_684 = buffer.data(li + 684);
    const auto *li_689 = buffer.data(li + 689);
    const auto *li_690 = buffer.data(li + 690);
    const auto *li_693 = buffer.data(li + 693);
    const auto *li_694 = buffer.data(li + 694);
    const auto *li_695 = buffer.data(li + 695);
    const auto *li_696 = buffer.data(li + 696);
    const auto *li_697 = buffer.data(li + 697);
    const auto *li_698 = buffer.data(li + 698);
    const auto *li_699 = buffer.data(li + 699);

    const auto *lk_540 = buffer.data(lk + 540);
    const auto *lk_541 = buffer.data(lk + 541);
    const auto *lk_543 = buffer.data(lk + 543);
    const auto *lk_545 = buffer.data(lk + 545);
    const auto *lk_546 = buffer.data(lk + 546);
    const auto *lk_549 = buffer.data(lk + 549);
    const auto *lk_550 = buffer.data(lk + 550);
    const auto *lk_552 = buffer.data(lk + 552);
    const auto *lk_554 = buffer.data(lk + 554);
    const auto *lk_555 = buffer.data(lk + 555);
    const auto *lk_557 = buffer.data(lk + 557);
    const auto *lk_558 = buffer.data(lk + 558);
    const auto *lk_560 = buffer.data(lk + 560);
    const auto *lk_561 = buffer.data(lk + 561);
    const auto *lk_568 = buffer.data(lk + 568);
    const auto *lk_570 = buffer.data(lk + 570);
    const auto *lk_571 = buffer.data(lk + 571);
    const auto *lk_572 = buffer.data(lk + 572);
    const auto *lk_573 = buffer.data(lk + 573);
    const auto *lk_575 = buffer.data(lk + 575);
    const auto *lk_579 = buffer.data(lk + 579);
    const auto *lk_582 = buffer.data(lk + 582);
    const auto *lk_586 = buffer.data(lk + 586);
    const auto *lk_591 = buffer.data(lk + 591);
    const auto *lk_612 = buffer.data(lk + 612);
    const auto *lk_615 = buffer.data(lk + 615);
    const auto *lk_617 = buffer.data(lk + 617);
    const auto *lk_618 = buffer.data(lk + 618);
    const auto *lk_621 = buffer.data(lk + 621);
    const auto *lk_622 = buffer.data(lk + 622);
    const auto *lk_626 = buffer.data(lk + 626);
    const auto *lk_627 = buffer.data(lk + 627);
    const auto *lk_632 = buffer.data(lk + 632);
    const auto *lk_648 = buffer.data(lk + 648);
    const auto *lk_651 = buffer.data(lk + 651);
    const auto *lk_653 = buffer.data(lk + 653);
    const auto *lk_654 = buffer.data(lk + 654);
    const auto *lk_657 = buffer.data(lk + 657);
    const auto *lk_662 = buffer.data(lk + 662);
    const auto *lk_668 = buffer.data(lk + 668);
    const auto *lk_684 = buffer.data(lk + 684);
    const auto *lk_689 = buffer.data(lk + 689);
    const auto *lk_856 = buffer.data(lk + 856);
    const auto *lk_858 = buffer.data(lk + 858);
    const auto *lk_859 = buffer.data(lk + 859);
    const auto *lk_860 = buffer.data(lk + 860);
    const auto *lk_861 = buffer.data(lk + 861);
    const auto *lk_863 = buffer.data(lk + 863);
    const auto *lk_892 = buffer.data(lk + 892);
    const auto *lk_894 = buffer.data(lk + 894);
    const auto *lk_895 = buffer.data(lk + 895);
    const auto *lk_896 = buffer.data(lk + 896);
    const auto *lk_897 = buffer.data(lk + 897);
    const auto *lk_899 = buffer.data(lk + 899);

    const auto *mh0_456 = buffer.data(mh0 + 456);
    const auto *mh0_457 = buffer.data(mh0 + 457);
    const auto *mh0_458 = buffer.data(mh0 + 458);
    const auto *mh0_459 = buffer.data(mh0 + 459);
    const auto *mh0_461 = buffer.data(mh0 + 461);
    const auto *mh0_495 = buffer.data(mh0 + 495);
    const auto *mh0_500 = buffer.data(mh0 + 500);
    const auto *mh0_501 = buffer.data(mh0 + 501);
    const auto *mh0_516 = buffer.data(mh0 + 516);
    const auto *mh0_521 = buffer.data(mh0 + 521);
    const auto *mh0_522 = buffer.data(mh0 + 522);

    const auto *mh1_456 = buffer.data(mh1 + 456);
    const auto *mh1_457 = buffer.data(mh1 + 457);
    const auto *mh1_458 = buffer.data(mh1 + 458);
    const auto *mh1_459 = buffer.data(mh1 + 459);
    const auto *mh1_461 = buffer.data(mh1 + 461);
    const auto *mh1_495 = buffer.data(mh1 + 495);
    const auto *mh1_500 = buffer.data(mh1 + 500);
    const auto *mh1_501 = buffer.data(mh1 + 501);
    const auto *mh1_516 = buffer.data(mh1 + 516);
    const auto *mh1_521 = buffer.data(mh1 + 521);
    const auto *mh1_522 = buffer.data(mh1 + 522);

    const auto *mi_610 = buffer.data(mi + 610);
    const auto *mi_611 = buffer.data(mi + 611);
    const auto *mi_612 = buffer.data(mi + 612);
    const auto *mi_613 = buffer.data(mi + 613);
    const auto *mi_615 = buffer.data(mi + 615);
    const auto *mi_616 = buffer.data(mi + 616);
    const auto *mi_618 = buffer.data(mi + 618);
    const auto *mi_619 = buffer.data(mi + 619);
    const auto *mi_621 = buffer.data(mi + 621);
    const auto *mi_622 = buffer.data(mi + 622);
    const auto *mi_625 = buffer.data(mi + 625);
    const auto *mi_626 = buffer.data(mi + 626);
    const auto *mi_630 = buffer.data(mi + 630);
    const auto *mi_637 = buffer.data(mi + 637);
    const auto *mi_638 = buffer.data(mi + 638);
    const auto *mi_639 = buffer.data(mi + 639);
    const auto *mi_640 = buffer.data(mi + 640);
    const auto *mi_641 = buffer.data(mi + 641);
    const auto *mi_642 = buffer.data(mi + 642);
    const auto *mi_643 = buffer.data(mi + 643);
    const auto *mi_644 = buffer.data(mi + 644);
    const auto *mi_646 = buffer.data(mi + 646);
    const auto *mi_647 = buffer.data(mi + 647);
    const auto *mi_649 = buffer.data(mi + 649);
    const auto *mi_650 = buffer.data(mi + 650);
    const auto *mi_653 = buffer.data(mi + 653);
    const auto *mi_654 = buffer.data(mi + 654);
    const auto *mi_656 = buffer.data(mi + 656);
    const auto *mi_658 = buffer.data(mi + 658);
    const auto *mi_661 = buffer.data(mi + 661);
    const auto *mi_662 = buffer.data(mi + 662);
    const auto *mi_665 = buffer.data(mi + 665);
    const auto *mi_666 = buffer.data(mi + 666);
    const auto *mi_667 = buffer.data(mi + 667);
    const auto *mi_668 = buffer.data(mi + 668);
    const auto *mi_669 = buffer.data(mi + 669);
    const auto *mi_670 = buffer.data(mi + 670);
    const auto *mi_671 = buffer.data(mi + 671);
    const auto *mi_672 = buffer.data(mi + 672);
    const auto *mi_674 = buffer.data(mi + 674);
    const auto *mi_675 = buffer.data(mi + 675);
    const auto *mi_677 = buffer.data(mi + 677);
    const auto *mi_678 = buffer.data(mi + 678);
    const auto *mi_681 = buffer.data(mi + 681);
    const auto *mi_682 = buffer.data(mi + 682);
    const auto *mi_684 = buffer.data(mi + 684);
    const auto *mi_686 = buffer.data(mi + 686);
    const auto *mi_689 = buffer.data(mi + 689);
    const auto *mi_690 = buffer.data(mi + 690);
    const auto *mi_693 = buffer.data(mi + 693);
    const auto *mi_694 = buffer.data(mi + 694);
    const auto *mi_695 = buffer.data(mi + 695);
    const auto *mi_696 = buffer.data(mi + 696);
    const auto *mi_697 = buffer.data(mi + 697);
    const auto *mi_698 = buffer.data(mi + 698);
    const auto *mi_699 = buffer.data(mi + 699);
    const auto *mi_700 = buffer.data(mi + 700);
    const auto *mi_702 = buffer.data(mi + 702);
    const auto *mi_703 = buffer.data(mi + 703);

#pragma omp simd aligned(t_786, t_787, t_788, pb_z, mh0_456, mh0_457, mh0_458, mh1_456, \
                         mh1_457, mh1_458, mi_610, mi_611, mi_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_3 * mh0_456[k]
                   - f_4 * mh1_456[k]
                   + pb_z[k] * mi_610[k];

        t_787[k] = f_5 * mh0_457[k]
                   - f_6 * mh1_457[k]
                   + pb_z[k] * mi_611[k];

        t_788[k] = f_7 * mh0_458[k]
                   - f_8 * mh1_458[k]
                   + pb_z[k] * mi_612[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, t_792, pa_z, pb_y, pb_z, li_447, lk_540, \
                         mh0_459, mh0_461, mh1_459, mh1_461, mi_613, \
                         mi_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = f_9 * mh0_459[k]
                   - f_10 * mh1_459[k]
                   + pb_z[k] * mi_613[k];

        t_790[k] = f_24 * li_447[k]
                   + pb_y[k] * mi_615[k];

        t_791[k] = f_1 * mh0_461[k]
                   - f_2 * mh1_461[k]
                   + pb_z[k] * mi_615[k];

        t_792[k] = pa_z[k] * lk_540[k];
    }

#pragma omp simd aligned(t_793, t_794, t_795, t_796, t_797, pa_z, pb_y, pb_z, li_420, li_422, \
                         li_450, lk_541, lk_543, lk_545, mi_616, \
                         mi_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_793[k] = pa_z[k] * lk_541[k];

        t_794[k] = f_11 * li_420[k]
                   + pb_z[k] * mi_616[k];

        t_795[k] = pa_z[k] * lk_543[k];

        t_796[k] = f_15 * li_450[k]
                   + pb_y[k] * mi_618[k];

        t_797[k] = f_12 * li_422[k]
                   + pa_z[k] * lk_545[k];
    }

#pragma omp simd aligned(t_798, t_799, t_800, t_801, t_802, pa_z, pb_y, pb_z, li_423, li_425, \
                         li_453, lk_546, lk_549, lk_550, mi_619, \
                         mi_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_798[k] = pa_z[k] * lk_546[k];

        t_799[k] = f_11 * li_423[k]
                   + pb_z[k] * mi_619[k];

        t_800[k] = f_15 * li_453[k]
                   + pb_y[k] * mi_621[k];

        t_801[k] = f_13 * li_425[k]
                   + pa_z[k] * lk_549[k];

        t_802[k] = pa_z[k] * lk_550[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, t_806, pa_z, pb_y, pb_z, li_426, li_427, li_429, \
                         li_457, lk_552, lk_554, mi_622, mi_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_11 * li_426[k]
                   + pb_z[k] * mi_622[k];

        t_804[k] = f_12 * li_427[k]
                   + pa_z[k] * lk_552[k];

        t_805[k] = f_15 * li_457[k]
                   + pb_y[k] * mi_625[k];

        t_806[k] = f_14 * li_429[k]
                   + pa_z[k] * lk_554[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, t_810, pa_z, pb_z, li_430, li_431, li_432, \
                         lk_555, lk_557, lk_558, mi_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = pa_z[k] * lk_555[k];

        t_808[k] = f_11 * li_430[k]
                   + pb_z[k] * mi_626[k];

        t_809[k] = f_12 * li_431[k]
                   + pa_z[k] * lk_557[k];

        t_810[k] = f_13 * li_432[k]
                   + pa_z[k] * lk_558[k];
    }

#pragma omp simd aligned(t_811, t_812, t_813, t_814, pa_z, pb_x, pb_y, li_434, li_462, li_638, \
                         lk_560, lk_561, mi_630, mi_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_811[k] = f_15 * li_462[k]
                   + pb_y[k] * mi_630[k];

        t_812[k] = f_15 * li_434[k]
                   + pa_z[k] * lk_560[k];

        t_813[k] = pa_z[k] * lk_561[k];

        t_814[k] = f_13 * li_638[k]
                   + pb_x[k] * mi_638[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, t_819, pb_x, li_639, li_640, li_641, \
                         li_642, li_643, mi_639, mi_640, mi_641, mi_642, \
                         mi_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_13 * li_639[k]
                   + pb_x[k] * mi_639[k];

        t_816[k] = f_13 * li_640[k]
                   + pb_x[k] * mi_640[k];

        t_817[k] = f_13 * li_641[k]
                   + pb_x[k] * mi_641[k];

        t_818[k] = f_13 * li_642[k]
                   + pb_x[k] * mi_642[k];

        t_819[k] = f_13 * li_643[k]
                   + pb_x[k] * mi_643[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, t_824, pa_z, pb_z, li_441, li_442, \
                         li_443, li_444, lk_568, lk_570, lk_571, lk_572, \
                         mi_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = pa_z[k] * lk_568[k];

        t_821[k] = f_11 * li_441[k]
                   + pb_z[k] * mi_637[k];

        t_822[k] = f_12 * li_442[k]
                   + pa_z[k] * lk_570[k];

        t_823[k] = f_13 * li_443[k]
                   + pa_z[k] * lk_571[k];

        t_824[k] = f_14 * li_444[k]
                   + pa_z[k] * lk_572[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, pa_y, pa_z, pb_y, kk0_432, kk1_432, \
                         li_445, li_447, li_475, lk_573, lk_575, lk_612, \
                         mi_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = f_15 * li_445[k]
                   + pa_z[k] * lk_573[k];

        t_826[k] = f_15 * li_475[k]
                   + pb_y[k] * mi_643[k];

        t_827[k] = f_17 * li_447[k]
                   + pa_z[k] * lk_575[k];

        t_828[k] = f_27 * kk0_432[k]
                   - f_28 * kk1_432[k]
                   + pa_y[k] * lk_612[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, pa_z, pb_y, pb_z, kk0_363, kk1_363, \
                         li_448, li_476, li_478, lk_579, mi_644, \
                         mi_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = f_14 * li_476[k]
                   + pb_y[k] * mi_644[k];

        t_830[k] = f_12 * li_448[k]
                   + pb_z[k] * mi_644[k];

        t_831[k] = f_18 * kk0_363[k]
                   - f_19 * kk1_363[k]
                   + pa_z[k] * lk_579[k];

        t_832[k] = f_14 * li_478[k]
                   + pb_y[k] * mi_646[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pa_y, pa_z, pb_z, kk0_366, kk0_437, kk1_366, \
                         kk1_437, li_451, lk_582, lk_617, mi_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_27 * kk0_437[k]
                   - f_28 * kk1_437[k]
                   + pa_y[k] * lk_617[k];

        t_834[k] = f_18 * kk0_366[k]
                   - f_19 * kk1_366[k]
                   + pa_z[k] * lk_582[k];

        t_835[k] = f_12 * li_451[k]
                   + pb_z[k] * mi_647[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, pa_y, pa_z, pb_y, kk0_370, kk0_441, kk1_370, \
                         kk1_441, li_481, lk_586, lk_621, mi_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_14 * li_481[k]
                   + pb_y[k] * mi_649[k];

        t_837[k] = f_27 * kk0_441[k]
                   - f_28 * kk1_441[k]
                   + pa_y[k] * lk_621[k];

        t_838[k] = f_18 * kk0_370[k]
                   - f_19 * kk1_370[k]
                   + pa_z[k] * lk_586[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, pb_x, pb_y, pb_z, li_454, li_485, li_656, \
                         mh0_495, mh1_495, mi_650, mi_653, mi_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = f_12 * li_454[k]
                   + pb_z[k] * mi_650[k];

        t_840[k] = f_13 * li_656[k]
                   + f_5 * mh0_495[k]
                   - f_6 * mh1_495[k]
                   + pb_x[k] * mi_656[k];

        t_841[k] = f_14 * li_485[k]
                   + pb_y[k] * mi_653[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, pa_y, pa_z, pb_z, kk0_375, kk0_446, kk1_375, \
                         kk1_446, li_458, lk_591, lk_626, mi_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_27 * kk0_446[k]
                   - f_28 * kk1_446[k]
                   + pa_y[k] * lk_626[k];

        t_843[k] = f_18 * kk0_375[k]
                   - f_19 * kk1_375[k]
                   + pa_z[k] * lk_591[k];

        t_844[k] = f_12 * li_458[k]
                   + pb_z[k] * mi_654[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pb_x, pb_y, li_490, li_661, li_662, mh0_500, \
                         mh0_501, mh1_500, mh1_501, mi_658, mi_661, \
                         mi_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_13 * li_661[k]
                   + f_3 * mh0_500[k]
                   - f_4 * mh1_500[k]
                   + pb_x[k] * mi_661[k];

        t_846[k] = f_13 * li_662[k]
                   + f_3 * mh0_501[k]
                   - f_4 * mh1_501[k]
                   + pb_x[k] * mi_662[k];

        t_847[k] = f_14 * li_490[k]
                   + pb_y[k] * mi_658[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, pa_y, pb_x, kk0_452, kk1_452, li_665, \
                         li_666, li_667, lk_632, mi_665, mi_666, \
                         mi_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_27 * kk0_452[k]
                   - f_28 * kk1_452[k]
                   + pa_y[k] * lk_632[k];

        t_849[k] = f_13 * li_665[k]
                   + pb_x[k] * mi_665[k];

        t_850[k] = f_13 * li_666[k]
                   + pb_x[k] * mi_666[k];

        t_851[k] = f_13 * li_667[k]
                   + pb_x[k] * mi_667[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pb_x, li_668, li_669, li_670, li_671, \
                         mi_668, mi_669, mi_670, mi_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_13 * li_668[k]
                   + pb_x[k] * mi_668[k];

        t_853[k] = f_13 * li_669[k]
                   + pb_x[k] * mi_669[k];

        t_854[k] = f_13 * li_670[k]
                   + pb_x[k] * mi_670[k];

        t_855[k] = f_13 * li_671[k]
                   + pb_x[k] * mi_671[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pa_x, pb_z, kk0_856, kk0_858, kk1_856, kk1_858, \
                         li_469, lk_856, lk_858, mi_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_22 * kk0_856[k]
                   - f_23 * kk1_856[k]
                   + pa_x[k] * lk_856[k];

        t_857[k] = f_12 * li_469[k]
                   + pb_z[k] * mi_665[k];

        t_858[k] = f_22 * kk0_858[k]
                   - f_23 * kk1_858[k]
                   + pa_x[k] * lk_858[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pa_x, kk0_859, kk0_860, kk0_861, kk1_859, \
                         kk1_860, kk1_861, lk_859, lk_860, lk_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_22 * kk0_859[k]
                   - f_23 * kk1_859[k]
                   + pa_x[k] * lk_859[k];

        t_860[k] = f_22 * kk0_860[k]
                   - f_23 * kk1_860[k]
                   + pa_x[k] * lk_860[k];

        t_861[k] = f_22 * kk0_861[k]
                   - f_23 * kk1_861[k]
                   + pa_x[k] * lk_861[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, pa_x, pa_y, pb_y, kk0_468, kk0_863, kk1_468, \
                         kk1_863, li_503, lk_648, lk_863, mi_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_14 * li_503[k]
                   + pb_y[k] * mi_671[k];

        t_863[k] = f_22 * kk0_863[k]
                   - f_23 * kk1_863[k]
                   + pa_x[k] * lk_863[k];

        t_864[k] = f_22 * kk0_468[k]
                   - f_23 * kk1_468[k]
                   + pa_y[k] * lk_648[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, pa_z, pb_y, pb_z, kk0_399, kk1_399, \
                         li_476, li_504, li_506, lk_615, mi_672, \
                         mi_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_13 * li_504[k]
                   + pb_y[k] * mi_672[k];

        t_866[k] = f_13 * li_476[k]
                   + pb_z[k] * mi_672[k];

        t_867[k] = f_22 * kk0_399[k]
                   - f_23 * kk1_399[k]
                   + pa_z[k] * lk_615[k];

        t_868[k] = f_13 * li_506[k]
                   + pb_y[k] * mi_674[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, pa_y, pa_z, pb_z, kk0_402, kk0_473, kk1_402, \
                         kk1_473, li_479, lk_618, lk_653, mi_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_22 * kk0_473[k]
                   - f_23 * kk1_473[k]
                   + pa_y[k] * lk_653[k];

        t_870[k] = f_22 * kk0_402[k]
                   - f_23 * kk1_402[k]
                   + pa_z[k] * lk_618[k];

        t_871[k] = f_13 * li_479[k]
                   + pb_z[k] * mi_675[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, pa_y, pa_z, pb_y, kk0_406, kk0_477, kk1_406, \
                         kk1_477, li_509, lk_622, lk_657, mi_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_13 * li_509[k]
                   + pb_y[k] * mi_677[k];

        t_873[k] = f_22 * kk0_477[k]
                   - f_23 * kk1_477[k]
                   + pa_y[k] * lk_657[k];

        t_874[k] = f_22 * kk0_406[k]
                   - f_23 * kk1_406[k]
                   + pa_z[k] * lk_622[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, pb_x, pb_y, pb_z, li_482, li_513, li_684, \
                         mh0_516, mh1_516, mi_678, mi_681, mi_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = f_13 * li_482[k]
                   + pb_z[k] * mi_678[k];

        t_876[k] = f_13 * li_684[k]
                   + f_5 * mh0_516[k]
                   - f_6 * mh1_516[k]
                   + pb_x[k] * mi_684[k];

        t_877[k] = f_13 * li_513[k]
                   + pb_y[k] * mi_681[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, pa_y, pa_z, pb_z, kk0_411, kk0_482, kk1_411, \
                         kk1_482, li_486, lk_627, lk_662, mi_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = f_22 * kk0_482[k]
                   - f_23 * kk1_482[k]
                   + pa_y[k] * lk_662[k];

        t_879[k] = f_22 * kk0_411[k]
                   - f_23 * kk1_411[k]
                   + pa_z[k] * lk_627[k];

        t_880[k] = f_13 * li_486[k]
                   + pb_z[k] * mi_682[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, pb_x, pb_y, li_518, li_689, li_690, mh0_521, \
                         mh0_522, mh1_521, mh1_522, mi_686, mi_689, \
                         mi_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = f_13 * li_689[k]
                   + f_3 * mh0_521[k]
                   - f_4 * mh1_521[k]
                   + pb_x[k] * mi_689[k];

        t_882[k] = f_13 * li_690[k]
                   + f_3 * mh0_522[k]
                   - f_4 * mh1_522[k]
                   + pb_x[k] * mi_690[k];

        t_883[k] = f_13 * li_518[k]
                   + pb_y[k] * mi_686[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, pa_y, pb_x, kk0_488, kk1_488, li_693, \
                         li_694, li_695, lk_668, mi_693, mi_694, \
                         mi_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_22 * kk0_488[k]
                   - f_23 * kk1_488[k]
                   + pa_y[k] * lk_668[k];

        t_885[k] = f_13 * li_693[k]
                   + pb_x[k] * mi_693[k];

        t_886[k] = f_13 * li_694[k]
                   + pb_x[k] * mi_694[k];

        t_887[k] = f_13 * li_695[k]
                   + pb_x[k] * mi_695[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, pb_x, li_696, li_697, li_698, li_699, \
                         mi_696, mi_697, mi_698, mi_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = f_13 * li_696[k]
                   + pb_x[k] * mi_696[k];

        t_889[k] = f_13 * li_697[k]
                   + pb_x[k] * mi_697[k];

        t_890[k] = f_13 * li_698[k]
                   + pb_x[k] * mi_698[k];

        t_891[k] = f_13 * li_699[k]
                   + pb_x[k] * mi_699[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, pa_x, pb_z, kk0_892, kk0_894, kk1_892, kk1_894, \
                         li_497, lk_892, lk_894, mi_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = f_22 * kk0_892[k]
                   - f_23 * kk1_892[k]
                   + pa_x[k] * lk_892[k];

        t_893[k] = f_13 * li_497[k]
                   + pb_z[k] * mi_693[k];

        t_894[k] = f_22 * kk0_894[k]
                   - f_23 * kk1_894[k]
                   + pa_x[k] * lk_894[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, pa_x, kk0_895, kk0_896, kk0_897, kk1_895, \
                         kk1_896, kk1_897, lk_895, lk_896, lk_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_22 * kk0_895[k]
                   - f_23 * kk1_895[k]
                   + pa_x[k] * lk_895[k];

        t_896[k] = f_22 * kk0_896[k]
                   - f_23 * kk1_896[k]
                   + pa_x[k] * lk_896[k];

        t_897[k] = f_22 * kk0_897[k]
                   - f_23 * kk1_897[k]
                   + pa_x[k] * lk_897[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, pa_x, pa_y, pb_y, kk0_504, kk0_899, kk1_504, \
                         kk1_899, li_531, lk_684, lk_899, mi_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_13 * li_531[k]
                   + pb_y[k] * mi_699[k];

        t_899[k] = f_22 * kk0_899[k]
                   - f_23 * kk1_899[k]
                   + pa_x[k] * lk_899[k];

        t_900[k] = f_18 * kk0_504[k]
                   - f_19 * kk1_504[k]
                   + pa_y[k] * lk_684[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, t_904, pa_z, pb_y, pb_z, kk0_435, kk1_435, \
                         li_504, li_532, li_534, lk_651, mi_700, \
                         mi_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = f_12 * li_532[k]
                   + pb_y[k] * mi_700[k];

        t_902[k] = f_14 * li_504[k]
                   + pb_z[k] * mi_700[k];

        t_903[k] = f_27 * kk0_435[k]
                   - f_28 * kk1_435[k]
                   + pa_z[k] * lk_651[k];

        t_904[k] = f_12 * li_534[k]
                   + pb_y[k] * mi_702[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pa_y, pa_z, pb_z, kk0_438, kk0_509, kk1_438, \
                         kk1_509, li_507, lk_654, lk_689, mi_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_18 * kk0_509[k]
                   - f_19 * kk1_509[k]
                   + pa_y[k] * lk_689[k];

        t_906[k] = f_27 * kk0_438[k]
                   - f_28 * kk1_438[k]
                   + pa_z[k] * lk_654[k];

        t_907[k] = f_14 * li_507[k]
                   + pb_z[k] * mi_703[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece7(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kk0,
                                            const size_t kk1, const size_t li, const size_t lk,
                                            const size_t mh0, const size_t mh1, const size_t mi,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_17 = 3.5 / p;
    const auto f_18 = 0.5 / alpha;
    const auto f_19 = 0.5 * beta / (alpha * p);
    const auto f_20 = 3.0 / alpha;
    const auto f_21 = 3.0 * beta / (alpha * p);
    const auto f_22 = 1.0 / alpha;
    const auto f_23 = beta / (alpha * p);
    const auto f_24 = 3.0 / p;
    const auto f_25 = 2.5 / alpha;
    const auto f_26 = 2.5 * beta / (alpha * p);
    const auto f_27 = 1.5 / alpha;
    const auto f_28 = 1.5 * beta / (alpha * p);

    auto *t_908 = buffer.data(target + 908);
    auto *t_909 = buffer.data(target + 909);
    auto *t_910 = buffer.data(target + 910);
    auto *t_911 = buffer.data(target + 911);
    auto *t_912 = buffer.data(target + 912);
    auto *t_913 = buffer.data(target + 913);
    auto *t_914 = buffer.data(target + 914);
    auto *t_915 = buffer.data(target + 915);
    auto *t_916 = buffer.data(target + 916);
    auto *t_917 = buffer.data(target + 917);
    auto *t_918 = buffer.data(target + 918);
    auto *t_919 = buffer.data(target + 919);
    auto *t_920 = buffer.data(target + 920);
    auto *t_921 = buffer.data(target + 921);
    auto *t_922 = buffer.data(target + 922);
    auto *t_923 = buffer.data(target + 923);
    auto *t_924 = buffer.data(target + 924);
    auto *t_925 = buffer.data(target + 925);
    auto *t_926 = buffer.data(target + 926);
    auto *t_927 = buffer.data(target + 927);
    auto *t_928 = buffer.data(target + 928);
    auto *t_929 = buffer.data(target + 929);
    auto *t_930 = buffer.data(target + 930);
    auto *t_931 = buffer.data(target + 931);
    auto *t_932 = buffer.data(target + 932);
    auto *t_933 = buffer.data(target + 933);
    auto *t_934 = buffer.data(target + 934);
    auto *t_935 = buffer.data(target + 935);
    auto *t_936 = buffer.data(target + 936);
    auto *t_937 = buffer.data(target + 937);
    auto *t_938 = buffer.data(target + 938);
    auto *t_939 = buffer.data(target + 939);
    auto *t_940 = buffer.data(target + 940);
    auto *t_941 = buffer.data(target + 941);
    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);
    auto *t_945 = buffer.data(target + 945);
    auto *t_946 = buffer.data(target + 946);
    auto *t_947 = buffer.data(target + 947);
    auto *t_948 = buffer.data(target + 948);
    auto *t_949 = buffer.data(target + 949);
    auto *t_950 = buffer.data(target + 950);
    auto *t_951 = buffer.data(target + 951);
    auto *t_952 = buffer.data(target + 952);
    auto *t_953 = buffer.data(target + 953);
    auto *t_954 = buffer.data(target + 954);
    auto *t_955 = buffer.data(target + 955);
    auto *t_956 = buffer.data(target + 956);
    auto *t_957 = buffer.data(target + 957);
    auto *t_958 = buffer.data(target + 958);
    auto *t_959 = buffer.data(target + 959);
    auto *t_960 = buffer.data(target + 960);
    auto *t_961 = buffer.data(target + 961);
    auto *t_962 = buffer.data(target + 962);
    auto *t_963 = buffer.data(target + 963);
    auto *t_964 = buffer.data(target + 964);
    auto *t_965 = buffer.data(target + 965);
    auto *t_966 = buffer.data(target + 966);
    auto *t_967 = buffer.data(target + 967);
    auto *t_968 = buffer.data(target + 968);
    auto *t_969 = buffer.data(target + 969);
    auto *t_970 = buffer.data(target + 970);
    auto *t_971 = buffer.data(target + 971);
    auto *t_972 = buffer.data(target + 972);
    auto *t_973 = buffer.data(target + 973);
    auto *t_974 = buffer.data(target + 974);
    auto *t_975 = buffer.data(target + 975);
    auto *t_976 = buffer.data(target + 976);
    auto *t_977 = buffer.data(target + 977);
    auto *t_978 = buffer.data(target + 978);
    auto *t_979 = buffer.data(target + 979);
    auto *t_980 = buffer.data(target + 980);
    auto *t_981 = buffer.data(target + 981);
    auto *t_982 = buffer.data(target + 982);
    auto *t_983 = buffer.data(target + 983);
    auto *t_984 = buffer.data(target + 984);
    auto *t_985 = buffer.data(target + 985);
    auto *t_986 = buffer.data(target + 986);
    auto *t_987 = buffer.data(target + 987);
    auto *t_988 = buffer.data(target + 988);
    auto *t_989 = buffer.data(target + 989);
    auto *t_990 = buffer.data(target + 990);
    auto *t_991 = buffer.data(target + 991);
    auto *t_992 = buffer.data(target + 992);
    auto *t_993 = buffer.data(target + 993);
    auto *t_994 = buffer.data(target + 994);
    auto *t_995 = buffer.data(target + 995);
    auto *t_996 = buffer.data(target + 996);
    auto *t_997 = buffer.data(target + 997);
    auto *t_998 = buffer.data(target + 998);
    auto *t_999 = buffer.data(target + 999);
    auto *t_1000 = buffer.data(target + 1000);
    auto *t_1001 = buffer.data(target + 1001);
    auto *t_1002 = buffer.data(target + 1002);
    auto *t_1003 = buffer.data(target + 1003);
    auto *t_1004 = buffer.data(target + 1004);
    auto *t_1005 = buffer.data(target + 1005);
    auto *t_1006 = buffer.data(target + 1006);
    auto *t_1007 = buffer.data(target + 1007);
    auto *t_1008 = buffer.data(target + 1008);
    auto *t_1009 = buffer.data(target + 1009);
    auto *t_1010 = buffer.data(target + 1010);
    auto *t_1011 = buffer.data(target + 1011);
    auto *t_1012 = buffer.data(target + 1012);
    auto *t_1013 = buffer.data(target + 1013);
    auto *t_1014 = buffer.data(target + 1014);
    auto *t_1015 = buffer.data(target + 1015);
    auto *t_1016 = buffer.data(target + 1016);
    auto *t_1017 = buffer.data(target + 1017);
    auto *t_1018 = buffer.data(target + 1018);
    auto *t_1019 = buffer.data(target + 1019);
    auto *t_1020 = buffer.data(target + 1020);
    auto *t_1021 = buffer.data(target + 1021);
    auto *t_1022 = buffer.data(target + 1022);
    auto *t_1023 = buffer.data(target + 1023);
    auto *t_1024 = buffer.data(target + 1024);
    auto *t_1025 = buffer.data(target + 1025);
    auto *t_1026 = buffer.data(target + 1026);
    auto *t_1027 = buffer.data(target + 1027);
    auto *t_1028 = buffer.data(target + 1028);
    auto *t_1029 = buffer.data(target + 1029);
    auto *t_1030 = buffer.data(target + 1030);
    auto *t_1031 = buffer.data(target + 1031);
    auto *t_1032 = buffer.data(target + 1032);
    auto *t_1033 = buffer.data(target + 1033);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kk0_442 = buffer.data(kk0 + 442);
    const auto *kk0_447 = buffer.data(kk0 + 447);
    const auto *kk0_504 = buffer.data(kk0 + 504);
    const auto *kk0_513 = buffer.data(kk0 + 513);
    const auto *kk0_518 = buffer.data(kk0 + 518);
    const auto *kk0_524 = buffer.data(kk0 + 524);
    const auto *kk0_540 = buffer.data(kk0 + 540);
    const auto *kk0_928 = buffer.data(kk0 + 928);
    const auto *kk0_930 = buffer.data(kk0 + 930);
    const auto *kk0_931 = buffer.data(kk0 + 931);
    const auto *kk0_932 = buffer.data(kk0 + 932);
    const auto *kk0_933 = buffer.data(kk0 + 933);
    const auto *kk0_935 = buffer.data(kk0 + 935);
    const auto *kk0_1007 = buffer.data(kk0 + 1007);

    const auto *kk1_442 = buffer.data(kk1 + 442);
    const auto *kk1_447 = buffer.data(kk1 + 447);
    const auto *kk1_504 = buffer.data(kk1 + 504);
    const auto *kk1_513 = buffer.data(kk1 + 513);
    const auto *kk1_518 = buffer.data(kk1 + 518);
    const auto *kk1_524 = buffer.data(kk1 + 524);
    const auto *kk1_540 = buffer.data(kk1 + 540);
    const auto *kk1_928 = buffer.data(kk1 + 928);
    const auto *kk1_930 = buffer.data(kk1 + 930);
    const auto *kk1_931 = buffer.data(kk1 + 931);
    const auto *kk1_932 = buffer.data(kk1 + 932);
    const auto *kk1_933 = buffer.data(kk1 + 933);
    const auto *kk1_935 = buffer.data(kk1 + 935);
    const auto *kk1_1007 = buffer.data(kk1 + 1007);

    const auto *li_510 = buffer.data(li + 510);
    const auto *li_514 = buffer.data(li + 514);
    const auto *li_525 = buffer.data(li + 525);
    const auto *li_535 = buffer.data(li + 535);
    const auto *li_537 = buffer.data(li + 537);
    const auto *li_538 = buffer.data(li + 538);
    const auto *li_541 = buffer.data(li + 541);
    const auto *li_542 = buffer.data(li + 542);
    const auto *li_546 = buffer.data(li + 546);
    const auto *li_553 = buffer.data(li + 553);
    const auto *li_559 = buffer.data(li + 559);
    const auto *li_560 = buffer.data(li + 560);
    const auto *li_561 = buffer.data(li + 561);
    const auto *li_562 = buffer.data(li + 562);
    const auto *li_563 = buffer.data(li + 563);
    const auto *li_565 = buffer.data(li + 565);
    const auto *li_566 = buffer.data(li + 566);
    const auto *li_568 = buffer.data(li + 568);
    const auto *li_569 = buffer.data(li + 569);
    const auto *li_570 = buffer.data(li + 570);
    const auto *li_572 = buffer.data(li + 572);
    const auto *li_573 = buffer.data(li + 573);
    const auto *li_574 = buffer.data(li + 574);
    const auto *li_581 = buffer.data(li + 581);
    const auto *li_583 = buffer.data(li + 583);
    const auto *li_584 = buffer.data(li + 584);
    const auto *li_585 = buffer.data(li + 585);
    const auto *li_586 = buffer.data(li + 586);
    const auto *li_587 = buffer.data(li + 587);
    const auto *li_588 = buffer.data(li + 588);
    const auto *li_593 = buffer.data(li + 593);
    const auto *li_597 = buffer.data(li + 597);
    const auto *li_602 = buffer.data(li + 602);
    const auto *li_712 = buffer.data(li + 712);
    const auto *li_717 = buffer.data(li + 717);
    const auto *li_718 = buffer.data(li + 718);
    const auto *li_721 = buffer.data(li + 721);
    const auto *li_722 = buffer.data(li + 722);
    const auto *li_723 = buffer.data(li + 723);
    const auto *li_724 = buffer.data(li + 724);
    const auto *li_725 = buffer.data(li + 725);
    const auto *li_726 = buffer.data(li + 726);
    const auto *li_727 = buffer.data(li + 727);
    const auto *li_749 = buffer.data(li + 749);
    const auto *li_750 = buffer.data(li + 750);
    const auto *li_751 = buffer.data(li + 751);
    const auto *li_752 = buffer.data(li + 752);
    const auto *li_753 = buffer.data(li + 753);
    const auto *li_754 = buffer.data(li + 754);
    const auto *li_761 = buffer.data(li + 761);
    const auto *li_765 = buffer.data(li + 765);
    const auto *li_770 = buffer.data(li + 770);
    const auto *li_776 = buffer.data(li + 776);
    const auto *li_777 = buffer.data(li + 777);
    const auto *li_778 = buffer.data(li + 778);
    const auto *li_779 = buffer.data(li + 779);
    const auto *li_780 = buffer.data(li + 780);
    const auto *li_781 = buffer.data(li + 781);
    const auto *li_783 = buffer.data(li + 783);
    const auto *li_787 = buffer.data(li + 787);
    const auto *li_790 = buffer.data(li + 790);
    const auto *li_794 = buffer.data(li + 794);
    const auto *li_799 = buffer.data(li + 799);
    const auto *li_805 = buffer.data(li + 805);
    const auto *li_807 = buffer.data(li + 807);
    const auto *li_808 = buffer.data(li + 808);
    const auto *li_809 = buffer.data(li + 809);

    const auto *lk_658 = buffer.data(lk + 658);
    const auto *lk_663 = buffer.data(lk + 663);
    const auto *lk_693 = buffer.data(lk + 693);
    const auto *lk_698 = buffer.data(lk + 698);
    const auto *lk_704 = buffer.data(lk + 704);
    const auto *lk_720 = buffer.data(lk + 720);
    const auto *lk_722 = buffer.data(lk + 722);
    const auto *lk_723 = buffer.data(lk + 723);
    const auto *lk_725 = buffer.data(lk + 725);
    const auto *lk_726 = buffer.data(lk + 726);
    const auto *lk_729 = buffer.data(lk + 729);
    const auto *lk_730 = buffer.data(lk + 730);
    const auto *lk_732 = buffer.data(lk + 732);
    const auto *lk_734 = buffer.data(lk + 734);
    const auto *lk_735 = buffer.data(lk + 735);
    const auto *lk_737 = buffer.data(lk + 737);
    const auto *lk_738 = buffer.data(lk + 738);
    const auto *lk_740 = buffer.data(lk + 740);
    const auto *lk_747 = buffer.data(lk + 747);
    const auto *lk_748 = buffer.data(lk + 748);
    const auto *lk_750 = buffer.data(lk + 750);
    const auto *lk_751 = buffer.data(lk + 751);
    const auto *lk_752 = buffer.data(lk + 752);
    const auto *lk_753 = buffer.data(lk + 753);
    const auto *lk_755 = buffer.data(lk + 755);
    const auto *lk_756 = buffer.data(lk + 756);
    const auto *lk_928 = buffer.data(lk + 928);
    const auto *lk_930 = buffer.data(lk + 930);
    const auto *lk_931 = buffer.data(lk + 931);
    const auto *lk_932 = buffer.data(lk + 932);
    const auto *lk_933 = buffer.data(lk + 933);
    const auto *lk_935 = buffer.data(lk + 935);
    const auto *lk_1007 = buffer.data(lk + 1007);

    const auto *mh0_537 = buffer.data(mh0 + 537);
    const auto *mh0_542 = buffer.data(mh0 + 542);
    const auto *mh0_543 = buffer.data(mh0 + 543);
    const auto *mh0_567 = buffer.data(mh0 + 567);
    const auto *mh0_568 = buffer.data(mh0 + 568);
    const auto *mh0_570 = buffer.data(mh0 + 570);
    const auto *mh0_572 = buffer.data(mh0 + 572);
    const auto *mh0_573 = buffer.data(mh0 + 573);
    const auto *mh0_575 = buffer.data(mh0 + 575);
    const auto *mh0_576 = buffer.data(mh0 + 576);
    const auto *mh0_581 = buffer.data(mh0 + 581);
    const auto *mh0_582 = buffer.data(mh0 + 582);
    const auto *mh0_584 = buffer.data(mh0 + 584);
    const auto *mh0_585 = buffer.data(mh0 + 585);
    const auto *mh0_586 = buffer.data(mh0 + 586);
    const auto *mh0_587 = buffer.data(mh0 + 587);
    const auto *mh0_588 = buffer.data(mh0 + 588);
    const auto *mh0_590 = buffer.data(mh0 + 590);
    const auto *mh0_591 = buffer.data(mh0 + 591);
    const auto *mh0_593 = buffer.data(mh0 + 593);
    const auto *mh0_594 = buffer.data(mh0 + 594);
    const auto *mh0_595 = buffer.data(mh0 + 595);
    const auto *mh0_597 = buffer.data(mh0 + 597);
    const auto *mh0_598 = buffer.data(mh0 + 598);
    const auto *mh0_603 = buffer.data(mh0 + 603);

    const auto *mh1_537 = buffer.data(mh1 + 537);
    const auto *mh1_542 = buffer.data(mh1 + 542);
    const auto *mh1_543 = buffer.data(mh1 + 543);
    const auto *mh1_567 = buffer.data(mh1 + 567);
    const auto *mh1_568 = buffer.data(mh1 + 568);
    const auto *mh1_570 = buffer.data(mh1 + 570);
    const auto *mh1_572 = buffer.data(mh1 + 572);
    const auto *mh1_573 = buffer.data(mh1 + 573);
    const auto *mh1_575 = buffer.data(mh1 + 575);
    const auto *mh1_576 = buffer.data(mh1 + 576);
    const auto *mh1_581 = buffer.data(mh1 + 581);
    const auto *mh1_582 = buffer.data(mh1 + 582);
    const auto *mh1_584 = buffer.data(mh1 + 584);
    const auto *mh1_585 = buffer.data(mh1 + 585);
    const auto *mh1_586 = buffer.data(mh1 + 586);
    const auto *mh1_587 = buffer.data(mh1 + 587);
    const auto *mh1_588 = buffer.data(mh1 + 588);
    const auto *mh1_590 = buffer.data(mh1 + 590);
    const auto *mh1_591 = buffer.data(mh1 + 591);
    const auto *mh1_593 = buffer.data(mh1 + 593);
    const auto *mh1_594 = buffer.data(mh1 + 594);
    const auto *mh1_595 = buffer.data(mh1 + 595);
    const auto *mh1_597 = buffer.data(mh1 + 597);
    const auto *mh1_598 = buffer.data(mh1 + 598);
    const auto *mh1_603 = buffer.data(mh1 + 603);

    const auto *mi_705 = buffer.data(mi + 705);
    const auto *mi_706 = buffer.data(mi + 706);
    const auto *mi_709 = buffer.data(mi + 709);
    const auto *mi_710 = buffer.data(mi + 710);
    const auto *mi_712 = buffer.data(mi + 712);
    const auto *mi_714 = buffer.data(mi + 714);
    const auto *mi_717 = buffer.data(mi + 717);
    const auto *mi_718 = buffer.data(mi + 718);
    const auto *mi_721 = buffer.data(mi + 721);
    const auto *mi_722 = buffer.data(mi + 722);
    const auto *mi_723 = buffer.data(mi + 723);
    const auto *mi_724 = buffer.data(mi + 724);
    const auto *mi_725 = buffer.data(mi + 725);
    const auto *mi_726 = buffer.data(mi + 726);
    const auto *mi_727 = buffer.data(mi + 727);
    const auto *mi_728 = buffer.data(mi + 728);
    const auto *mi_730 = buffer.data(mi + 730);
    const auto *mi_731 = buffer.data(mi + 731);
    const auto *mi_733 = buffer.data(mi + 733);
    const auto *mi_734 = buffer.data(mi + 734);
    const auto *mi_737 = buffer.data(mi + 737);
    const auto *mi_738 = buffer.data(mi + 738);
    const auto *mi_742 = buffer.data(mi + 742);
    const auto *mi_749 = buffer.data(mi + 749);
    const auto *mi_750 = buffer.data(mi + 750);
    const auto *mi_751 = buffer.data(mi + 751);
    const auto *mi_752 = buffer.data(mi + 752);
    const auto *mi_753 = buffer.data(mi + 753);
    const auto *mi_754 = buffer.data(mi + 754);
    const auto *mi_755 = buffer.data(mi + 755);
    const auto *mi_756 = buffer.data(mi + 756);
    const auto *mi_757 = buffer.data(mi + 757);
    const auto *mi_758 = buffer.data(mi + 758);
    const auto *mi_759 = buffer.data(mi + 759);
    const auto *mi_761 = buffer.data(mi + 761);
    const auto *mi_762 = buffer.data(mi + 762);
    const auto *mi_764 = buffer.data(mi + 764);
    const auto *mi_765 = buffer.data(mi + 765);
    const auto *mi_766 = buffer.data(mi + 766);
    const auto *mi_768 = buffer.data(mi + 768);
    const auto *mi_769 = buffer.data(mi + 769);
    const auto *mi_770 = buffer.data(mi + 770);
    const auto *mi_776 = buffer.data(mi + 776);
    const auto *mi_777 = buffer.data(mi + 777);
    const auto *mi_778 = buffer.data(mi + 778);
    const auto *mi_779 = buffer.data(mi + 779);
    const auto *mi_780 = buffer.data(mi + 780);
    const auto *mi_781 = buffer.data(mi + 781);
    const auto *mi_782 = buffer.data(mi + 782);
    const auto *mi_783 = buffer.data(mi + 783);
    const auto *mi_784 = buffer.data(mi + 784);
    const auto *mi_785 = buffer.data(mi + 785);
    const auto *mi_786 = buffer.data(mi + 786);
    const auto *mi_787 = buffer.data(mi + 787);
    const auto *mi_789 = buffer.data(mi + 789);
    const auto *mi_790 = buffer.data(mi + 790);
    const auto *mi_791 = buffer.data(mi + 791);
    const auto *mi_793 = buffer.data(mi + 793);
    const auto *mi_794 = buffer.data(mi + 794);
    const auto *mi_795 = buffer.data(mi + 795);
    const auto *mi_796 = buffer.data(mi + 796);
    const auto *mi_798 = buffer.data(mi + 798);
    const auto *mi_799 = buffer.data(mi + 799);
    const auto *mi_805 = buffer.data(mi + 805);
    const auto *mi_807 = buffer.data(mi + 807);
    const auto *mi_808 = buffer.data(mi + 808);
    const auto *mi_809 = buffer.data(mi + 809);

#pragma omp simd aligned(t_908, t_909, t_910, pa_y, pa_z, pb_y, kk0_442, kk0_513, kk1_442, \
                         kk1_513, li_537, lk_658, lk_693, mi_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_12 * li_537[k]
                   + pb_y[k] * mi_705[k];

        t_909[k] = f_18 * kk0_513[k]
                   - f_19 * kk1_513[k]
                   + pa_y[k] * lk_693[k];

        t_910[k] = f_27 * kk0_442[k]
                   - f_28 * kk1_442[k]
                   + pa_z[k] * lk_658[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, pb_x, pb_y, pb_z, li_510, li_541, li_712, \
                         mh0_537, mh1_537, mi_706, mi_709, mi_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_14 * li_510[k]
                   + pb_z[k] * mi_706[k];

        t_912[k] = f_13 * li_712[k]
                   + f_5 * mh0_537[k]
                   - f_6 * mh1_537[k]
                   + pb_x[k] * mi_712[k];

        t_913[k] = f_12 * li_541[k]
                   + pb_y[k] * mi_709[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, pa_y, pa_z, pb_z, kk0_447, kk0_518, kk1_447, \
                         kk1_518, li_514, lk_663, lk_698, mi_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_18 * kk0_518[k]
                   - f_19 * kk1_518[k]
                   + pa_y[k] * lk_698[k];

        t_915[k] = f_27 * kk0_447[k]
                   - f_28 * kk1_447[k]
                   + pa_z[k] * lk_663[k];

        t_916[k] = f_14 * li_514[k]
                   + pb_z[k] * mi_710[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, pb_x, pb_y, li_546, li_717, li_718, mh0_542, \
                         mh0_543, mh1_542, mh1_543, mi_714, mi_717, \
                         mi_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = f_13 * li_717[k]
                   + f_3 * mh0_542[k]
                   - f_4 * mh1_542[k]
                   + pb_x[k] * mi_717[k];

        t_918[k] = f_13 * li_718[k]
                   + f_3 * mh0_543[k]
                   - f_4 * mh1_543[k]
                   + pb_x[k] * mi_718[k];

        t_919[k] = f_12 * li_546[k]
                   + pb_y[k] * mi_714[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pa_y, pb_x, kk0_524, kk1_524, li_721, \
                         li_722, li_723, lk_704, mi_721, mi_722, \
                         mi_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = f_18 * kk0_524[k]
                   - f_19 * kk1_524[k]
                   + pa_y[k] * lk_704[k];

        t_921[k] = f_13 * li_721[k]
                   + pb_x[k] * mi_721[k];

        t_922[k] = f_13 * li_722[k]
                   + pb_x[k] * mi_722[k];

        t_923[k] = f_13 * li_723[k]
                   + pb_x[k] * mi_723[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pb_x, li_724, li_725, li_726, li_727, \
                         mi_724, mi_725, mi_726, mi_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_13 * li_724[k]
                   + pb_x[k] * mi_724[k];

        t_925[k] = f_13 * li_725[k]
                   + pb_x[k] * mi_725[k];

        t_926[k] = f_13 * li_726[k]
                   + pb_x[k] * mi_726[k];

        t_927[k] = f_13 * li_727[k]
                   + pb_x[k] * mi_727[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, pa_x, pb_z, kk0_928, kk0_930, kk1_928, kk1_930, \
                         li_525, lk_928, lk_930, mi_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_22 * kk0_928[k]
                   - f_23 * kk1_928[k]
                   + pa_x[k] * lk_928[k];

        t_929[k] = f_14 * li_525[k]
                   + pb_z[k] * mi_721[k];

        t_930[k] = f_22 * kk0_930[k]
                   - f_23 * kk1_930[k]
                   + pa_x[k] * lk_930[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, pa_x, kk0_931, kk0_932, kk0_933, kk1_931, \
                         kk1_932, kk1_933, lk_931, lk_932, lk_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_22 * kk0_931[k]
                   - f_23 * kk1_931[k]
                   + pa_x[k] * lk_931[k];

        t_932[k] = f_22 * kk0_932[k]
                   - f_23 * kk1_932[k]
                   + pa_x[k] * lk_932[k];

        t_933[k] = f_22 * kk0_933[k]
                   - f_23 * kk1_933[k]
                   + pa_x[k] * lk_933[k];
    }

#pragma omp simd aligned(t_934, t_935, t_936, t_937, pa_x, pa_y, pb_y, kk0_935, kk1_935, \
                         li_559, li_560, lk_720, lk_935, mi_727, \
                         mi_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_12 * li_559[k]
                   + pb_y[k] * mi_727[k];

        t_935[k] = f_22 * kk0_935[k]
                   - f_23 * kk1_935[k]
                   + pa_x[k] * lk_935[k];

        t_936[k] = pa_y[k] * lk_720[k];

        t_937[k] = f_11 * li_560[k]
                   + pb_y[k] * mi_728[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, t_942, pa_y, pb_y, li_561, li_562, \
                         li_563, lk_722, lk_723, lk_725, lk_726, \
                         mi_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = pa_y[k] * lk_722[k];

        t_939[k] = f_12 * li_561[k]
                   + pa_y[k] * lk_723[k];

        t_940[k] = f_11 * li_562[k]
                   + pb_y[k] * mi_730[k];

        t_941[k] = pa_y[k] * lk_725[k];

        t_942[k] = f_13 * li_563[k]
                   + pa_y[k] * lk_726[k];
    }

#pragma omp simd aligned(t_943, t_944, t_945, t_946, pa_y, pb_y, pb_z, li_535, li_565, li_566, \
                         lk_729, lk_730, mi_731, mi_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = f_15 * li_535[k]
                   + pb_z[k] * mi_731[k];

        t_944[k] = f_11 * li_565[k]
                   + pb_y[k] * mi_733[k];

        t_945[k] = pa_y[k] * lk_729[k];

        t_946[k] = f_14 * li_566[k]
                   + pa_y[k] * lk_730[k];
    }

#pragma omp simd aligned(t_947, t_948, t_949, t_950, pa_y, pb_y, pb_z, li_538, li_568, li_569, \
                         lk_732, lk_734, mi_734, mi_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_947[k] = f_15 * li_538[k]
                   + pb_z[k] * mi_734[k];

        t_948[k] = f_12 * li_568[k]
                   + pa_y[k] * lk_732[k];

        t_949[k] = f_11 * li_569[k]
                   + pb_y[k] * mi_737[k];

        t_950[k] = pa_y[k] * lk_734[k];
    }

#pragma omp simd aligned(t_951, t_952, t_953, t_954, pa_y, pb_z, li_542, li_570, li_572, \
                         li_573, lk_735, lk_737, lk_738, mi_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_951[k] = f_15 * li_570[k]
                   + pa_y[k] * lk_735[k];

        t_952[k] = f_15 * li_542[k]
                   + pb_z[k] * mi_738[k];

        t_953[k] = f_13 * li_572[k]
                   + pa_y[k] * lk_737[k];

        t_954[k] = f_12 * li_573[k]
                   + pa_y[k] * lk_738[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, t_958, pa_y, pb_x, pb_y, li_574, li_749, li_750, \
                         lk_740, mi_742, mi_749, mi_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = f_11 * li_574[k]
                   + pb_y[k] * mi_742[k];

        t_956[k] = pa_y[k] * lk_740[k];

        t_957[k] = f_13 * li_749[k]
                   + pb_x[k] * mi_749[k];

        t_958[k] = f_13 * li_750[k]
                   + pb_x[k] * mi_750[k];
    }

#pragma omp simd aligned(t_959, t_960, t_961, t_962, t_963, pa_y, pb_x, li_751, li_752, \
                         li_753, li_754, lk_747, mi_751, mi_752, mi_753, \
                         mi_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_959[k] = f_13 * li_751[k]
                   + pb_x[k] * mi_751[k];

        t_960[k] = f_13 * li_752[k]
                   + pb_x[k] * mi_752[k];

        t_961[k] = f_13 * li_753[k]
                   + pb_x[k] * mi_753[k];

        t_962[k] = f_13 * li_754[k]
                   + pb_x[k] * mi_754[k];

        t_963[k] = pa_y[k] * lk_747[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pa_y, pb_z, li_553, li_581, li_583, \
                         li_584, lk_748, lk_750, lk_751, mi_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = f_17 * li_581[k]
                   + pa_y[k] * lk_748[k];

        t_965[k] = f_15 * li_553[k]
                   + pb_z[k] * mi_749[k];

        t_966[k] = f_15 * li_583[k]
                   + pa_y[k] * lk_750[k];

        t_967[k] = f_14 * li_584[k]
                   + pa_y[k] * lk_751[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, pa_y, pb_y, li_585, li_586, li_587, \
                         lk_752, lk_753, lk_755, mi_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_13 * li_585[k]
                   + pa_y[k] * lk_752[k];

        t_969[k] = f_12 * li_586[k]
                   + pa_y[k] * lk_753[k];

        t_970[k] = f_11 * li_587[k]
                   + pb_y[k] * mi_755[k];

        t_971[k] = pa_y[k] * lk_755[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, pa_z, pb_y, pb_z, kk0_504, kk1_504, \
                         li_560, lk_720, mh0_567, mh1_567, mi_756, \
                         mi_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = f_25 * kk0_504[k]
                   - f_26 * kk1_504[k]
                   + pa_z[k] * lk_720[k];

        t_973[k] = pb_y[k] * mi_756[k];

        t_974[k] = f_24 * li_560[k]
                   + pb_z[k] * mi_756[k];

        t_975[k] = f_3 * mh0_567[k]
                   - f_4 * mh1_567[k]
                   + pb_y[k] * mi_757[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, t_979, pb_x, pb_y, pb_z, li_563, li_761, \
                         mh0_568, mh0_572, mh1_568, mh1_572, mi_758, mi_759, \
                         mi_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = pb_y[k] * mi_758[k];

        t_977[k] = f_13 * li_761[k]
                   + f_9 * mh0_572[k]
                   - f_10 * mh1_572[k]
                   + pb_x[k] * mi_761[k];

        t_978[k] = f_5 * mh0_568[k]
                   - f_6 * mh1_568[k]
                   + pb_y[k] * mi_759[k];

        t_979[k] = f_24 * li_563[k]
                   + pb_z[k] * mi_759[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, pb_x, pb_y, pb_z, li_566, li_765, \
                         mh0_570, mh0_576, mh1_570, mh1_576, mi_761, mi_762, \
                         mi_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = pb_y[k] * mi_761[k];

        t_981[k] = f_13 * li_765[k]
                   + f_7 * mh0_576[k]
                   - f_8 * mh1_576[k]
                   + pb_x[k] * mi_765[k];

        t_982[k] = f_7 * mh0_570[k]
                   - f_8 * mh1_570[k]
                   + pb_y[k] * mi_762[k];

        t_983[k] = f_24 * li_566[k]
                   + pb_z[k] * mi_762[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, pb_x, pb_y, li_770, mh0_572, mh0_581, mh1_572, \
                         mh1_581, mi_764, mi_765, mi_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_3 * mh0_572[k]
                   - f_4 * mh1_572[k]
                   + pb_y[k] * mi_764[k];

        t_985[k] = pb_y[k] * mi_765[k];

        t_986[k] = f_13 * li_770[k]
                   + f_5 * mh0_581[k]
                   - f_6 * mh1_581[k]
                   + pb_x[k] * mi_770[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, t_990, pb_y, pb_z, li_570, mh0_573, mh0_575, \
                         mh0_576, mh1_573, mh1_575, mh1_576, mi_766, mi_768, \
                         mi_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = f_9 * mh0_573[k]
                   - f_10 * mh1_573[k]
                   + pb_y[k] * mi_766[k];

        t_988[k] = f_24 * li_570[k]
                   + pb_z[k] * mi_766[k];

        t_989[k] = f_5 * mh0_575[k]
                   - f_6 * mh1_575[k]
                   + pb_y[k] * mi_768[k];

        t_990[k] = f_3 * mh0_576[k]
                   - f_4 * mh1_576[k]
                   + pb_y[k] * mi_769[k];
    }

#pragma omp simd aligned(t_991, t_992, t_993, t_994, pb_x, pb_y, li_776, li_777, li_778, \
                         mh0_587, mh1_587, mi_770, mi_776, mi_777, \
                         mi_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_991[k] = pb_y[k] * mi_770[k];

        t_992[k] = f_13 * li_776[k]
                   + f_3 * mh0_587[k]
                   - f_4 * mh1_587[k]
                   + pb_x[k] * mi_776[k];

        t_993[k] = f_13 * li_777[k]
                   + pb_x[k] * mi_777[k];

        t_994[k] = f_13 * li_778[k]
                   + pb_x[k] * mi_778[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, t_999, pb_x, pb_y, li_779, li_780, \
                         li_781, li_783, mi_776, mi_779, mi_780, mi_781, \
                         mi_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_13 * li_779[k]
                   + pb_x[k] * mi_779[k];

        t_996[k] = f_13 * li_780[k]
                   + pb_x[k] * mi_780[k];

        t_997[k] = f_13 * li_781[k]
                   + pb_x[k] * mi_781[k];

        t_998[k] = pb_y[k] * mi_776[k];

        t_999[k] = f_13 * li_783[k]
                   + pb_x[k] * mi_783[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pb_y, pb_z, li_581, mh0_582, mh0_584, \
                         mh0_585, mh1_582, mh1_584, mh1_585, mi_777, mi_779, \
                         mi_780 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = f_1 * mh0_582[k]
                    - f_2 * mh1_582[k]
                    + pb_y[k] * mi_777[k];

        t_1001[k] = f_24 * li_581[k]
                    + pb_z[k] * mi_777[k];

        t_1002[k] = f_9 * mh0_584[k]
                    - f_10 * mh1_584[k]
                    + pb_y[k] * mi_779[k];

        t_1003[k] = f_7 * mh0_585[k]
                    - f_8 * mh1_585[k]
                    + pb_y[k] * mi_780[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pa_x, pb_y, kk0_1007, kk1_1007, \
                         lk_1007, mh0_586, mh0_587, mh1_586, mh1_587, mi_781, mi_782, \
                         mi_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_5 * mh0_586[k]
                    - f_6 * mh1_586[k]
                    + pb_y[k] * mi_781[k];

        t_1005[k] = f_3 * mh0_587[k]
                    - f_4 * mh1_587[k]
                    + pb_y[k] * mi_782[k];

        t_1006[k] = pb_y[k] * mi_783[k];

        t_1007[k] = f_22 * kk0_1007[k]
                    - f_23 * kk1_1007[k]
                    + pa_x[k] * lk_1007[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, pa_y, pb_y, pb_z, kk0_540, kk1_540, li_588, \
                         lk_756, mi_784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_20 * kk0_540[k]
                    - f_21 * kk1_540[k]
                    + pa_y[k] * lk_756[k];

        t_1009[k] = f_17 * li_588[k]
                    + pb_y[k] * mi_784[k];

        t_1010[k] = pb_z[k] * mi_784[k];
    }

#pragma omp simd aligned(t_1011, t_1012, t_1013, pb_x, pb_z, li_787, mh0_588, mh0_591, \
                         mh1_588, mh1_591, mi_785, mi_786, mi_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1011[k] = f_12 * li_787[k]
                    + f_9 * mh0_591[k]
                    - f_10 * mh1_591[k]
                    + pb_x[k] * mi_787[k];

        t_1012[k] = pb_z[k] * mi_785[k];

        t_1013[k] = f_3 * mh0_588[k]
                    - f_4 * mh1_588[k]
                    + pb_z[k] * mi_786[k];
    }

#pragma omp simd aligned(t_1014, t_1015, t_1016, t_1017, pb_x, pb_y, pb_z, li_593, li_790, \
                         mh0_590, mh0_594, mh1_590, mh1_594, mi_787, mi_789, \
                         mi_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1014[k] = f_12 * li_790[k]
                    + f_7 * mh0_594[k]
                    - f_8 * mh1_594[k]
                    + pb_x[k] * mi_790[k];

        t_1015[k] = pb_z[k] * mi_787[k];

        t_1016[k] = f_17 * li_593[k]
                    + pb_y[k] * mi_789[k];

        t_1017[k] = f_5 * mh0_590[k]
                    - f_6 * mh1_590[k]
                    + pb_z[k] * mi_789[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, pb_x, pb_z, li_794, mh0_591, mh0_598, \
                         mh1_591, mh1_598, mi_790, mi_791, mi_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = f_12 * li_794[k]
                    + f_5 * mh0_598[k]
                    - f_6 * mh1_598[k]
                    + pb_x[k] * mi_794[k];

        t_1019[k] = pb_z[k] * mi_790[k];

        t_1020[k] = f_3 * mh0_591[k]
                    - f_4 * mh1_591[k]
                    + pb_z[k] * mi_791[k];
    }

#pragma omp simd aligned(t_1021, t_1022, t_1023, t_1024, pb_x, pb_y, pb_z, li_597, li_799, \
                         mh0_593, mh0_603, mh1_593, mh1_603, mi_793, mi_794, \
                         mi_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_17 * li_597[k]
                    + pb_y[k] * mi_793[k];

        t_1022[k] = f_7 * mh0_593[k]
                    - f_8 * mh1_593[k]
                    + pb_z[k] * mi_793[k];

        t_1023[k] = f_12 * li_799[k]
                    + f_3 * mh0_603[k]
                    - f_4 * mh1_603[k]
                    + pb_x[k] * mi_799[k];

        t_1024[k] = pb_z[k] * mi_794[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, t_1028, pb_y, pb_z, li_602, mh0_594, mh0_595, \
                         mh0_597, mh1_594, mh1_595, mh1_597, mi_795, mi_796, \
                         mi_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_3 * mh0_594[k]
                    - f_4 * mh1_594[k]
                    + pb_z[k] * mi_795[k];

        t_1026[k] = f_5 * mh0_595[k]
                    - f_6 * mh1_595[k]
                    + pb_z[k] * mi_796[k];

        t_1027[k] = f_17 * li_602[k]
                    + pb_y[k] * mi_798[k];

        t_1028[k] = f_9 * mh0_597[k]
                    - f_10 * mh1_597[k]
                    + pb_z[k] * mi_798[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, t_1032, t_1033, pb_x, pb_z, li_805, li_807, \
                         li_808, li_809, mi_799, mi_805, mi_807, mi_808, \
                         mi_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = f_12 * li_805[k]
                    + pb_x[k] * mi_805[k];

        t_1030[k] = pb_z[k] * mi_799[k];

        t_1031[k] = f_12 * li_807[k]
                    + pb_x[k] * mi_807[k];

        t_1032[k] = f_12 * li_808[k]
                    + pb_x[k] * mi_808[k];

        t_1033[k] = f_12 * li_809[k]
                    + pb_x[k] * mi_809[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece8(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kk0,
                                            const size_t kk1, const size_t li, const size_t lk,
                                            const size_t mh0, const size_t mh1, const size_t mi,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_17 = 3.5 / p;
    const auto f_18 = 0.5 / alpha;
    const auto f_19 = 0.5 * beta / (alpha * p);
    const auto f_22 = 1.0 / alpha;
    const auto f_23 = beta / (alpha * p);
    const auto f_24 = 3.0 / p;
    const auto f_27 = 1.5 / alpha;
    const auto f_28 = 1.5 * beta / (alpha * p);
    const auto f_29 = 2.0 / alpha;
    const auto f_30 = 2.0 * beta / (alpha * p);

    auto *t_1034 = buffer.data(target + 1034);
    auto *t_1035 = buffer.data(target + 1035);
    auto *t_1036 = buffer.data(target + 1036);
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);
    auto *t_1041 = buffer.data(target + 1041);
    auto *t_1042 = buffer.data(target + 1042);
    auto *t_1043 = buffer.data(target + 1043);
    auto *t_1044 = buffer.data(target + 1044);
    auto *t_1045 = buffer.data(target + 1045);
    auto *t_1046 = buffer.data(target + 1046);
    auto *t_1047 = buffer.data(target + 1047);
    auto *t_1048 = buffer.data(target + 1048);
    auto *t_1049 = buffer.data(target + 1049);
    auto *t_1050 = buffer.data(target + 1050);
    auto *t_1051 = buffer.data(target + 1051);
    auto *t_1052 = buffer.data(target + 1052);
    auto *t_1053 = buffer.data(target + 1053);
    auto *t_1054 = buffer.data(target + 1054);
    auto *t_1055 = buffer.data(target + 1055);
    auto *t_1056 = buffer.data(target + 1056);
    auto *t_1057 = buffer.data(target + 1057);
    auto *t_1058 = buffer.data(target + 1058);
    auto *t_1059 = buffer.data(target + 1059);
    auto *t_1060 = buffer.data(target + 1060);
    auto *t_1061 = buffer.data(target + 1061);
    auto *t_1062 = buffer.data(target + 1062);
    auto *t_1063 = buffer.data(target + 1063);
    auto *t_1064 = buffer.data(target + 1064);
    auto *t_1065 = buffer.data(target + 1065);
    auto *t_1066 = buffer.data(target + 1066);
    auto *t_1067 = buffer.data(target + 1067);
    auto *t_1068 = buffer.data(target + 1068);
    auto *t_1069 = buffer.data(target + 1069);
    auto *t_1070 = buffer.data(target + 1070);
    auto *t_1071 = buffer.data(target + 1071);
    auto *t_1072 = buffer.data(target + 1072);
    auto *t_1073 = buffer.data(target + 1073);
    auto *t_1074 = buffer.data(target + 1074);
    auto *t_1075 = buffer.data(target + 1075);
    auto *t_1076 = buffer.data(target + 1076);
    auto *t_1077 = buffer.data(target + 1077);
    auto *t_1078 = buffer.data(target + 1078);
    auto *t_1079 = buffer.data(target + 1079);
    auto *t_1080 = buffer.data(target + 1080);
    auto *t_1081 = buffer.data(target + 1081);
    auto *t_1082 = buffer.data(target + 1082);
    auto *t_1083 = buffer.data(target + 1083);
    auto *t_1084 = buffer.data(target + 1084);
    auto *t_1085 = buffer.data(target + 1085);
    auto *t_1086 = buffer.data(target + 1086);
    auto *t_1087 = buffer.data(target + 1087);
    auto *t_1088 = buffer.data(target + 1088);
    auto *t_1089 = buffer.data(target + 1089);
    auto *t_1090 = buffer.data(target + 1090);
    auto *t_1091 = buffer.data(target + 1091);
    auto *t_1092 = buffer.data(target + 1092);
    auto *t_1093 = buffer.data(target + 1093);
    auto *t_1094 = buffer.data(target + 1094);
    auto *t_1095 = buffer.data(target + 1095);
    auto *t_1096 = buffer.data(target + 1096);
    auto *t_1097 = buffer.data(target + 1097);
    auto *t_1098 = buffer.data(target + 1098);
    auto *t_1099 = buffer.data(target + 1099);
    auto *t_1100 = buffer.data(target + 1100);
    auto *t_1101 = buffer.data(target + 1101);
    auto *t_1102 = buffer.data(target + 1102);
    auto *t_1103 = buffer.data(target + 1103);
    auto *t_1104 = buffer.data(target + 1104);
    auto *t_1105 = buffer.data(target + 1105);
    auto *t_1106 = buffer.data(target + 1106);
    auto *t_1107 = buffer.data(target + 1107);
    auto *t_1108 = buffer.data(target + 1108);
    auto *t_1109 = buffer.data(target + 1109);
    auto *t_1110 = buffer.data(target + 1110);
    auto *t_1111 = buffer.data(target + 1111);
    auto *t_1112 = buffer.data(target + 1112);
    auto *t_1113 = buffer.data(target + 1113);
    auto *t_1114 = buffer.data(target + 1114);
    auto *t_1115 = buffer.data(target + 1115);
    auto *t_1116 = buffer.data(target + 1116);
    auto *t_1117 = buffer.data(target + 1117);
    auto *t_1118 = buffer.data(target + 1118);
    auto *t_1119 = buffer.data(target + 1119);
    auto *t_1120 = buffer.data(target + 1120);
    auto *t_1121 = buffer.data(target + 1121);
    auto *t_1122 = buffer.data(target + 1122);
    auto *t_1123 = buffer.data(target + 1123);
    auto *t_1124 = buffer.data(target + 1124);
    auto *t_1125 = buffer.data(target + 1125);
    auto *t_1126 = buffer.data(target + 1126);
    auto *t_1127 = buffer.data(target + 1127);
    auto *t_1128 = buffer.data(target + 1128);
    auto *t_1129 = buffer.data(target + 1129);
    auto *t_1130 = buffer.data(target + 1130);
    auto *t_1131 = buffer.data(target + 1131);
    auto *t_1132 = buffer.data(target + 1132);
    auto *t_1133 = buffer.data(target + 1133);
    auto *t_1134 = buffer.data(target + 1134);
    auto *t_1135 = buffer.data(target + 1135);
    auto *t_1136 = buffer.data(target + 1136);
    auto *t_1137 = buffer.data(target + 1137);
    auto *t_1138 = buffer.data(target + 1138);
    auto *t_1139 = buffer.data(target + 1139);
    auto *t_1140 = buffer.data(target + 1140);
    auto *t_1141 = buffer.data(target + 1141);
    auto *t_1142 = buffer.data(target + 1142);
    auto *t_1143 = buffer.data(target + 1143);
    auto *t_1144 = buffer.data(target + 1144);
    auto *t_1145 = buffer.data(target + 1145);
    auto *t_1146 = buffer.data(target + 1146);
    auto *t_1147 = buffer.data(target + 1147);
    auto *t_1148 = buffer.data(target + 1148);
    auto *t_1149 = buffer.data(target + 1149);
    auto *t_1150 = buffer.data(target + 1150);
    auto *t_1151 = buffer.data(target + 1151);
    auto *t_1152 = buffer.data(target + 1152);
    auto *t_1153 = buffer.data(target + 1153);
    auto *t_1154 = buffer.data(target + 1154);
    auto *t_1155 = buffer.data(target + 1155);
    auto *t_1156 = buffer.data(target + 1156);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kk0_543 = buffer.data(kk0 + 543);
    const auto *kk0_546 = buffer.data(kk0 + 546);
    const auto *kk0_550 = buffer.data(kk0 + 550);
    const auto *kk0_555 = buffer.data(kk0 + 555);
    const auto *kk0_579 = buffer.data(kk0 + 579);
    const auto *kk0_582 = buffer.data(kk0 + 582);
    const auto *kk0_586 = buffer.data(kk0 + 586);
    const auto *kk0_591 = buffer.data(kk0 + 591);
    const auto *kk0_612 = buffer.data(kk0 + 612);
    const auto *kk0_615 = buffer.data(kk0 + 615);
    const auto *kk0_617 = buffer.data(kk0 + 617);
    const auto *kk0_621 = buffer.data(kk0 + 621);
    const auto *kk0_626 = buffer.data(kk0 + 626);
    const auto *kk0_632 = buffer.data(kk0 + 632);
    const auto *kk0_648 = buffer.data(kk0 + 648);
    const auto *kk0_653 = buffer.data(kk0 + 653);
    const auto *kk0_657 = buffer.data(kk0 + 657);
    const auto *kk0_662 = buffer.data(kk0 + 662);
    const auto *kk0_668 = buffer.data(kk0 + 668);
    const auto *kk0_684 = buffer.data(kk0 + 684);
    const auto *kk0_1036 = buffer.data(kk0 + 1036);
    const auto *kk0_1108 = buffer.data(kk0 + 1108);
    const auto *kk0_1110 = buffer.data(kk0 + 1110);
    const auto *kk0_1111 = buffer.data(kk0 + 1111);
    const auto *kk0_1112 = buffer.data(kk0 + 1112);
    const auto *kk0_1113 = buffer.data(kk0 + 1113);
    const auto *kk0_1115 = buffer.data(kk0 + 1115);
    const auto *kk0_1144 = buffer.data(kk0 + 1144);
    const auto *kk0_1146 = buffer.data(kk0 + 1146);
    const auto *kk0_1147 = buffer.data(kk0 + 1147);
    const auto *kk0_1148 = buffer.data(kk0 + 1148);
    const auto *kk0_1149 = buffer.data(kk0 + 1149);
    const auto *kk0_1151 = buffer.data(kk0 + 1151);

    const auto *kk1_543 = buffer.data(kk1 + 543);
    const auto *kk1_546 = buffer.data(kk1 + 546);
    const auto *kk1_550 = buffer.data(kk1 + 550);
    const auto *kk1_555 = buffer.data(kk1 + 555);
    const auto *kk1_579 = buffer.data(kk1 + 579);
    const auto *kk1_582 = buffer.data(kk1 + 582);
    const auto *kk1_586 = buffer.data(kk1 + 586);
    const auto *kk1_591 = buffer.data(kk1 + 591);
    const auto *kk1_612 = buffer.data(kk1 + 612);
    const auto *kk1_615 = buffer.data(kk1 + 615);
    const auto *kk1_617 = buffer.data(kk1 + 617);
    const auto *kk1_621 = buffer.data(kk1 + 621);
    const auto *kk1_626 = buffer.data(kk1 + 626);
    const auto *kk1_632 = buffer.data(kk1 + 632);
    const auto *kk1_648 = buffer.data(kk1 + 648);
    const auto *kk1_653 = buffer.data(kk1 + 653);
    const auto *kk1_657 = buffer.data(kk1 + 657);
    const auto *kk1_662 = buffer.data(kk1 + 662);
    const auto *kk1_668 = buffer.data(kk1 + 668);
    const auto *kk1_684 = buffer.data(kk1 + 684);
    const auto *kk1_1036 = buffer.data(kk1 + 1036);
    const auto *kk1_1108 = buffer.data(kk1 + 1108);
    const auto *kk1_1110 = buffer.data(kk1 + 1110);
    const auto *kk1_1111 = buffer.data(kk1 + 1111);
    const auto *kk1_1112 = buffer.data(kk1 + 1112);
    const auto *kk1_1113 = buffer.data(kk1 + 1113);
    const auto *kk1_1115 = buffer.data(kk1 + 1115);
    const auto *kk1_1144 = buffer.data(kk1 + 1144);
    const auto *kk1_1146 = buffer.data(kk1 + 1146);
    const auto *kk1_1147 = buffer.data(kk1 + 1147);
    const auto *kk1_1148 = buffer.data(kk1 + 1148);
    const auto *kk1_1149 = buffer.data(kk1 + 1149);
    const auto *kk1_1151 = buffer.data(kk1 + 1151);

    const auto *li_588 = buffer.data(li + 588);
    const auto *li_590 = buffer.data(li + 590);
    const auto *li_591 = buffer.data(li + 591);
    const auto *li_593 = buffer.data(li + 593);
    const auto *li_594 = buffer.data(li + 594);
    const auto *li_595 = buffer.data(li + 595);
    const auto *li_597 = buffer.data(li + 597);
    const auto *li_598 = buffer.data(li + 598);
    const auto *li_599 = buffer.data(li + 599);
    const auto *li_600 = buffer.data(li + 600);
    const auto *li_602 = buffer.data(li + 602);
    const auto *li_609 = buffer.data(li + 609);
    const auto *li_610 = buffer.data(li + 610);
    const auto *li_611 = buffer.data(li + 611);
    const auto *li_612 = buffer.data(li + 612);
    const auto *li_613 = buffer.data(li + 613);
    const auto *li_615 = buffer.data(li + 615);
    const auto *li_616 = buffer.data(li + 616);
    const auto *li_618 = buffer.data(li + 618);
    const auto *li_619 = buffer.data(li + 619);
    const auto *li_621 = buffer.data(li + 621);
    const auto *li_622 = buffer.data(li + 622);
    const auto *li_625 = buffer.data(li + 625);
    const auto *li_626 = buffer.data(li + 626);
    const auto *li_630 = buffer.data(li + 630);
    const auto *li_637 = buffer.data(li + 637);
    const auto *li_643 = buffer.data(li + 643);
    const auto *li_644 = buffer.data(li + 644);
    const auto *li_646 = buffer.data(li + 646);
    const auto *li_647 = buffer.data(li + 647);
    const auto *li_649 = buffer.data(li + 649);
    const auto *li_650 = buffer.data(li + 650);
    const auto *li_653 = buffer.data(li + 653);
    const auto *li_654 = buffer.data(li + 654);
    const auto *li_658 = buffer.data(li + 658);
    const auto *li_665 = buffer.data(li + 665);
    const auto *li_671 = buffer.data(li + 671);
    const auto *li_672 = buffer.data(li + 672);
    const auto *li_674 = buffer.data(li + 674);
    const auto *li_677 = buffer.data(li + 677);
    const auto *li_681 = buffer.data(li + 681);
    const auto *li_686 = buffer.data(li + 686);
    const auto *li_699 = buffer.data(li + 699);
    const auto *li_700 = buffer.data(li + 700);
    const auto *li_702 = buffer.data(li + 702);
    const auto *li_810 = buffer.data(li + 810);
    const auto *li_811 = buffer.data(li + 811);
    const auto *li_834 = buffer.data(li + 834);
    const auto *li_835 = buffer.data(li + 835);
    const auto *li_836 = buffer.data(li + 836);
    const auto *li_837 = buffer.data(li + 837);
    const auto *li_838 = buffer.data(li + 838);
    const auto *li_839 = buffer.data(li + 839);
    const auto *li_852 = buffer.data(li + 852);
    const auto *li_857 = buffer.data(li + 857);
    const auto *li_858 = buffer.data(li + 858);
    const auto *li_861 = buffer.data(li + 861);
    const auto *li_862 = buffer.data(li + 862);
    const auto *li_863 = buffer.data(li + 863);
    const auto *li_864 = buffer.data(li + 864);
    const auto *li_865 = buffer.data(li + 865);
    const auto *li_866 = buffer.data(li + 866);
    const auto *li_867 = buffer.data(li + 867);
    const auto *li_880 = buffer.data(li + 880);
    const auto *li_885 = buffer.data(li + 885);
    const auto *li_886 = buffer.data(li + 886);
    const auto *li_889 = buffer.data(li + 889);
    const auto *li_890 = buffer.data(li + 890);
    const auto *li_891 = buffer.data(li + 891);
    const auto *li_892 = buffer.data(li + 892);
    const auto *li_893 = buffer.data(li + 893);
    const auto *li_894 = buffer.data(li + 894);
    const auto *li_895 = buffer.data(li + 895);

    const auto *lk_756 = buffer.data(lk + 756);
    const auto *lk_757 = buffer.data(lk + 757);
    const auto *lk_759 = buffer.data(lk + 759);
    const auto *lk_761 = buffer.data(lk + 761);
    const auto *lk_762 = buffer.data(lk + 762);
    const auto *lk_765 = buffer.data(lk + 765);
    const auto *lk_766 = buffer.data(lk + 766);
    const auto *lk_768 = buffer.data(lk + 768);
    const auto *lk_770 = buffer.data(lk + 770);
    const auto *lk_771 = buffer.data(lk + 771);
    const auto *lk_773 = buffer.data(lk + 773);
    const auto *lk_774 = buffer.data(lk + 774);
    const auto *lk_776 = buffer.data(lk + 776);
    const auto *lk_777 = buffer.data(lk + 777);
    const auto *lk_784 = buffer.data(lk + 784);
    const auto *lk_786 = buffer.data(lk + 786);
    const auto *lk_787 = buffer.data(lk + 787);
    const auto *lk_788 = buffer.data(lk + 788);
    const auto *lk_789 = buffer.data(lk + 789);
    const auto *lk_791 = buffer.data(lk + 791);
    const auto *lk_795 = buffer.data(lk + 795);
    const auto *lk_798 = buffer.data(lk + 798);
    const auto *lk_802 = buffer.data(lk + 802);
    const auto *lk_807 = buffer.data(lk + 807);
    const auto *lk_828 = buffer.data(lk + 828);
    const auto *lk_831 = buffer.data(lk + 831);
    const auto *lk_833 = buffer.data(lk + 833);
    const auto *lk_834 = buffer.data(lk + 834);
    const auto *lk_837 = buffer.data(lk + 837);
    const auto *lk_838 = buffer.data(lk + 838);
    const auto *lk_842 = buffer.data(lk + 842);
    const auto *lk_843 = buffer.data(lk + 843);
    const auto *lk_848 = buffer.data(lk + 848);
    const auto *lk_864 = buffer.data(lk + 864);
    const auto *lk_867 = buffer.data(lk + 867);
    const auto *lk_869 = buffer.data(lk + 869);
    const auto *lk_873 = buffer.data(lk + 873);
    const auto *lk_878 = buffer.data(lk + 878);
    const auto *lk_884 = buffer.data(lk + 884);
    const auto *lk_900 = buffer.data(lk + 900);
    const auto *lk_1036 = buffer.data(lk + 1036);
    const auto *lk_1108 = buffer.data(lk + 1108);
    const auto *lk_1110 = buffer.data(lk + 1110);
    const auto *lk_1111 = buffer.data(lk + 1111);
    const auto *lk_1112 = buffer.data(lk + 1112);
    const auto *lk_1113 = buffer.data(lk + 1113);
    const auto *lk_1115 = buffer.data(lk + 1115);
    const auto *lk_1144 = buffer.data(lk + 1144);
    const auto *lk_1146 = buffer.data(lk + 1146);
    const auto *lk_1147 = buffer.data(lk + 1147);
    const auto *lk_1148 = buffer.data(lk + 1148);
    const auto *lk_1149 = buffer.data(lk + 1149);
    const auto *lk_1151 = buffer.data(lk + 1151);

    const auto *mh0_603 = buffer.data(mh0 + 603);
    const auto *mh0_604 = buffer.data(mh0 + 604);
    const auto *mh0_605 = buffer.data(mh0 + 605);
    const auto *mh0_606 = buffer.data(mh0 + 606);
    const auto *mh0_608 = buffer.data(mh0 + 608);
    const auto *mh0_642 = buffer.data(mh0 + 642);
    const auto *mh0_647 = buffer.data(mh0 + 647);
    const auto *mh0_648 = buffer.data(mh0 + 648);
    const auto *mh0_663 = buffer.data(mh0 + 663);
    const auto *mh0_668 = buffer.data(mh0 + 668);
    const auto *mh0_669 = buffer.data(mh0 + 669);

    const auto *mh1_603 = buffer.data(mh1 + 603);
    const auto *mh1_604 = buffer.data(mh1 + 604);
    const auto *mh1_605 = buffer.data(mh1 + 605);
    const auto *mh1_606 = buffer.data(mh1 + 606);
    const auto *mh1_608 = buffer.data(mh1 + 608);
    const auto *mh1_642 = buffer.data(mh1 + 642);
    const auto *mh1_647 = buffer.data(mh1 + 647);
    const auto *mh1_648 = buffer.data(mh1 + 648);
    const auto *mh1_663 = buffer.data(mh1 + 663);
    const auto *mh1_668 = buffer.data(mh1 + 668);
    const auto *mh1_669 = buffer.data(mh1 + 669);

    const auto *mi_805 = buffer.data(mi + 805);
    const auto *mi_806 = buffer.data(mi + 806);
    const auto *mi_807 = buffer.data(mi + 807);
    const auto *mi_808 = buffer.data(mi + 808);
    const auto *mi_809 = buffer.data(mi + 809);
    const auto *mi_810 = buffer.data(mi + 810);
    const auto *mi_811 = buffer.data(mi + 811);
    const auto *mi_812 = buffer.data(mi + 812);
    const auto *mi_814 = buffer.data(mi + 814);
    const auto *mi_815 = buffer.data(mi + 815);
    const auto *mi_817 = buffer.data(mi + 817);
    const auto *mi_818 = buffer.data(mi + 818);
    const auto *mi_821 = buffer.data(mi + 821);
    const auto *mi_822 = buffer.data(mi + 822);
    const auto *mi_826 = buffer.data(mi + 826);
    const auto *mi_833 = buffer.data(mi + 833);
    const auto *mi_834 = buffer.data(mi + 834);
    const auto *mi_835 = buffer.data(mi + 835);
    const auto *mi_836 = buffer.data(mi + 836);
    const auto *mi_837 = buffer.data(mi + 837);
    const auto *mi_838 = buffer.data(mi + 838);
    const auto *mi_839 = buffer.data(mi + 839);
    const auto *mi_840 = buffer.data(mi + 840);
    const auto *mi_842 = buffer.data(mi + 842);
    const auto *mi_843 = buffer.data(mi + 843);
    const auto *mi_845 = buffer.data(mi + 845);
    const auto *mi_846 = buffer.data(mi + 846);
    const auto *mi_849 = buffer.data(mi + 849);
    const auto *mi_850 = buffer.data(mi + 850);
    const auto *mi_852 = buffer.data(mi + 852);
    const auto *mi_854 = buffer.data(mi + 854);
    const auto *mi_857 = buffer.data(mi + 857);
    const auto *mi_858 = buffer.data(mi + 858);
    const auto *mi_861 = buffer.data(mi + 861);
    const auto *mi_862 = buffer.data(mi + 862);
    const auto *mi_863 = buffer.data(mi + 863);
    const auto *mi_864 = buffer.data(mi + 864);
    const auto *mi_865 = buffer.data(mi + 865);
    const auto *mi_866 = buffer.data(mi + 866);
    const auto *mi_867 = buffer.data(mi + 867);
    const auto *mi_868 = buffer.data(mi + 868);
    const auto *mi_870 = buffer.data(mi + 870);
    const auto *mi_871 = buffer.data(mi + 871);
    const auto *mi_873 = buffer.data(mi + 873);
    const auto *mi_874 = buffer.data(mi + 874);
    const auto *mi_877 = buffer.data(mi + 877);
    const auto *mi_878 = buffer.data(mi + 878);
    const auto *mi_880 = buffer.data(mi + 880);
    const auto *mi_882 = buffer.data(mi + 882);
    const auto *mi_885 = buffer.data(mi + 885);
    const auto *mi_886 = buffer.data(mi + 886);
    const auto *mi_889 = buffer.data(mi + 889);
    const auto *mi_890 = buffer.data(mi + 890);
    const auto *mi_891 = buffer.data(mi + 891);
    const auto *mi_892 = buffer.data(mi + 892);
    const auto *mi_893 = buffer.data(mi + 893);
    const auto *mi_894 = buffer.data(mi + 894);
    const auto *mi_895 = buffer.data(mi + 895);
    const auto *mi_896 = buffer.data(mi + 896);
    const auto *mi_898 = buffer.data(mi + 898);

#pragma omp simd aligned(t_1034, t_1035, t_1036, t_1037, pa_x, pb_x, pb_z, kk0_1036, kk1_1036, \
                         li_810, li_811, lk_1036, mi_805, mi_810, \
                         mi_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = f_12 * li_810[k]
                    + pb_x[k] * mi_810[k];

        t_1035[k] = f_12 * li_811[k]
                    + pb_x[k] * mi_811[k];

        t_1036[k] = f_18 * kk0_1036[k]
                    - f_19 * kk1_1036[k]
                    + pa_x[k] * lk_1036[k];

        t_1037[k] = pb_z[k] * mi_805[k];
    }

#pragma omp simd aligned(t_1038, t_1039, t_1040, pb_z, mh0_603, mh0_604, mh0_605, mh1_603, \
                         mh1_604, mh1_605, mi_806, mi_807, mi_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1038[k] = f_3 * mh0_603[k]
                    - f_4 * mh1_603[k]
                    + pb_z[k] * mi_806[k];

        t_1039[k] = f_5 * mh0_604[k]
                    - f_6 * mh1_604[k]
                    + pb_z[k] * mi_807[k];

        t_1040[k] = f_7 * mh0_605[k]
                    - f_8 * mh1_605[k]
                    + pb_z[k] * mi_808[k];
    }

#pragma omp simd aligned(t_1041, t_1042, t_1043, t_1044, pa_z, pb_y, pb_z, li_615, lk_756, \
                         mh0_606, mh0_608, mh1_606, mh1_608, mi_809, \
                         mi_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1041[k] = f_9 * mh0_606[k]
                    - f_10 * mh1_606[k]
                    + pb_z[k] * mi_809[k];

        t_1042[k] = f_17 * li_615[k]
                    + pb_y[k] * mi_811[k];

        t_1043[k] = f_1 * mh0_608[k]
                    - f_2 * mh1_608[k]
                    + pb_z[k] * mi_811[k];

        t_1044[k] = pa_z[k] * lk_756[k];
    }

#pragma omp simd aligned(t_1045, t_1046, t_1047, t_1048, t_1049, pa_z, pb_y, pb_z, li_588, \
                         li_590, li_618, lk_757, lk_759, lk_761, mi_812, \
                         mi_814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1045[k] = pa_z[k] * lk_757[k];

        t_1046[k] = f_11 * li_588[k]
                    + pb_z[k] * mi_812[k];

        t_1047[k] = pa_z[k] * lk_759[k];

        t_1048[k] = f_24 * li_618[k]
                    + pb_y[k] * mi_814[k];

        t_1049[k] = f_12 * li_590[k]
                    + pa_z[k] * lk_761[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, t_1053, t_1054, pa_z, pb_y, pb_z, li_591, \
                         li_593, li_621, lk_762, lk_765, lk_766, mi_815, \
                         mi_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = pa_z[k] * lk_762[k];

        t_1051[k] = f_11 * li_591[k]
                    + pb_z[k] * mi_815[k];

        t_1052[k] = f_24 * li_621[k]
                    + pb_y[k] * mi_817[k];

        t_1053[k] = f_13 * li_593[k]
                    + pa_z[k] * lk_765[k];

        t_1054[k] = pa_z[k] * lk_766[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, t_1058, pa_z, pb_y, pb_z, li_594, li_595, \
                         li_597, li_625, lk_768, lk_770, mi_818, \
                         mi_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_11 * li_594[k]
                    + pb_z[k] * mi_818[k];

        t_1056[k] = f_12 * li_595[k]
                    + pa_z[k] * lk_768[k];

        t_1057[k] = f_24 * li_625[k]
                    + pb_y[k] * mi_821[k];

        t_1058[k] = f_14 * li_597[k]
                    + pa_z[k] * lk_770[k];
    }

#pragma omp simd aligned(t_1059, t_1060, t_1061, t_1062, pa_z, pb_z, li_598, li_599, li_600, \
                         lk_771, lk_773, lk_774, mi_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = pa_z[k] * lk_771[k];

        t_1060[k] = f_11 * li_598[k]
                    + pb_z[k] * mi_822[k];

        t_1061[k] = f_12 * li_599[k]
                    + pa_z[k] * lk_773[k];

        t_1062[k] = f_13 * li_600[k]
                    + pa_z[k] * lk_774[k];
    }

#pragma omp simd aligned(t_1063, t_1064, t_1065, t_1066, pa_z, pb_x, pb_y, li_602, li_630, \
                         li_834, lk_776, lk_777, mi_826, mi_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1063[k] = f_24 * li_630[k]
                    + pb_y[k] * mi_826[k];

        t_1064[k] = f_15 * li_602[k]
                    + pa_z[k] * lk_776[k];

        t_1065[k] = pa_z[k] * lk_777[k];

        t_1066[k] = f_12 * li_834[k]
                    + pb_x[k] * mi_834[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, t_1070, t_1071, pb_x, li_835, li_836, li_837, \
                         li_838, li_839, mi_835, mi_836, mi_837, mi_838, \
                         mi_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = f_12 * li_835[k]
                    + pb_x[k] * mi_835[k];

        t_1068[k] = f_12 * li_836[k]
                    + pb_x[k] * mi_836[k];

        t_1069[k] = f_12 * li_837[k]
                    + pb_x[k] * mi_837[k];

        t_1070[k] = f_12 * li_838[k]
                    + pb_x[k] * mi_838[k];

        t_1071[k] = f_12 * li_839[k]
                    + pb_x[k] * mi_839[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, t_1075, t_1076, pa_z, pb_z, li_609, li_610, \
                         li_611, li_612, lk_784, lk_786, lk_787, lk_788, \
                         mi_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = pa_z[k] * lk_784[k];

        t_1073[k] = f_11 * li_609[k]
                    + pb_z[k] * mi_833[k];

        t_1074[k] = f_12 * li_610[k]
                    + pa_z[k] * lk_786[k];

        t_1075[k] = f_13 * li_611[k]
                    + pa_z[k] * lk_787[k];

        t_1076[k] = f_14 * li_612[k]
                    + pa_z[k] * lk_788[k];
    }

#pragma omp simd aligned(t_1077, t_1078, t_1079, t_1080, pa_y, pa_z, pb_y, kk0_612, kk1_612, \
                         li_613, li_615, li_643, lk_789, lk_791, lk_828, \
                         mi_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1077[k] = f_15 * li_613[k]
                    + pa_z[k] * lk_789[k];

        t_1078[k] = f_24 * li_643[k]
                    + pb_y[k] * mi_839[k];

        t_1079[k] = f_17 * li_615[k]
                    + pa_z[k] * lk_791[k];

        t_1080[k] = f_29 * kk0_612[k]
                    - f_30 * kk1_612[k]
                    + pa_y[k] * lk_828[k];
    }

#pragma omp simd aligned(t_1081, t_1082, t_1083, t_1084, pa_z, pb_y, pb_z, kk0_543, kk1_543, \
                         li_616, li_644, li_646, lk_795, mi_840, \
                         mi_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1081[k] = f_15 * li_644[k]
                    + pb_y[k] * mi_840[k];

        t_1082[k] = f_12 * li_616[k]
                    + pb_z[k] * mi_840[k];

        t_1083[k] = f_18 * kk0_543[k]
                    - f_19 * kk1_543[k]
                    + pa_z[k] * lk_795[k];

        t_1084[k] = f_15 * li_646[k]
                    + pb_y[k] * mi_842[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pa_y, pa_z, pb_z, kk0_546, kk0_617, kk1_546, \
                         kk1_617, li_619, lk_798, lk_833, mi_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = f_29 * kk0_617[k]
                    - f_30 * kk1_617[k]
                    + pa_y[k] * lk_833[k];

        t_1086[k] = f_18 * kk0_546[k]
                    - f_19 * kk1_546[k]
                    + pa_z[k] * lk_798[k];

        t_1087[k] = f_12 * li_619[k]
                    + pb_z[k] * mi_843[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pa_y, pa_z, pb_y, kk0_550, kk0_621, kk1_550, \
                         kk1_621, li_649, lk_802, lk_837, mi_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_15 * li_649[k]
                    + pb_y[k] * mi_845[k];

        t_1089[k] = f_29 * kk0_621[k]
                    - f_30 * kk1_621[k]
                    + pa_y[k] * lk_837[k];

        t_1090[k] = f_18 * kk0_550[k]
                    - f_19 * kk1_550[k]
                    + pa_z[k] * lk_802[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pb_x, pb_y, pb_z, li_622, li_653, li_852, \
                         mh0_642, mh1_642, mi_846, mi_849, mi_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_12 * li_622[k]
                    + pb_z[k] * mi_846[k];

        t_1092[k] = f_12 * li_852[k]
                    + f_5 * mh0_642[k]
                    - f_6 * mh1_642[k]
                    + pb_x[k] * mi_852[k];

        t_1093[k] = f_15 * li_653[k]
                    + pb_y[k] * mi_849[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pa_y, pa_z, pb_z, kk0_555, kk0_626, kk1_555, \
                         kk1_626, li_626, lk_807, lk_842, mi_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = f_29 * kk0_626[k]
                    - f_30 * kk1_626[k]
                    + pa_y[k] * lk_842[k];

        t_1095[k] = f_18 * kk0_555[k]
                    - f_19 * kk1_555[k]
                    + pa_z[k] * lk_807[k];

        t_1096[k] = f_12 * li_626[k]
                    + pb_z[k] * mi_850[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pb_x, pb_y, li_658, li_857, li_858, mh0_647, \
                         mh0_648, mh1_647, mh1_648, mi_854, mi_857, \
                         mi_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_12 * li_857[k]
                    + f_3 * mh0_647[k]
                    - f_4 * mh1_647[k]
                    + pb_x[k] * mi_857[k];

        t_1098[k] = f_12 * li_858[k]
                    + f_3 * mh0_648[k]
                    - f_4 * mh1_648[k]
                    + pb_x[k] * mi_858[k];

        t_1099[k] = f_15 * li_658[k]
                    + pb_y[k] * mi_854[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, t_1103, pa_y, pb_x, kk0_632, kk1_632, li_861, \
                         li_862, li_863, lk_848, mi_861, mi_862, \
                         mi_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_29 * kk0_632[k]
                    - f_30 * kk1_632[k]
                    + pa_y[k] * lk_848[k];

        t_1101[k] = f_12 * li_861[k]
                    + pb_x[k] * mi_861[k];

        t_1102[k] = f_12 * li_862[k]
                    + pb_x[k] * mi_862[k];

        t_1103[k] = f_12 * li_863[k]
                    + pb_x[k] * mi_863[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, t_1107, pb_x, li_864, li_865, li_866, li_867, \
                         mi_864, mi_865, mi_866, mi_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_12 * li_864[k]
                    + pb_x[k] * mi_864[k];

        t_1105[k] = f_12 * li_865[k]
                    + pb_x[k] * mi_865[k];

        t_1106[k] = f_12 * li_866[k]
                    + pb_x[k] * mi_866[k];

        t_1107[k] = f_12 * li_867[k]
                    + pb_x[k] * mi_867[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, pa_x, pb_z, kk0_1108, kk0_1110, kk1_1108, \
                         kk1_1110, li_637, lk_1108, lk_1110, mi_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_18 * kk0_1108[k]
                    - f_19 * kk1_1108[k]
                    + pa_x[k] * lk_1108[k];

        t_1109[k] = f_12 * li_637[k]
                    + pb_z[k] * mi_861[k];

        t_1110[k] = f_18 * kk0_1110[k]
                    - f_19 * kk1_1110[k]
                    + pa_x[k] * lk_1110[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, pa_x, kk0_1111, kk0_1112, kk0_1113, kk1_1111, \
                         kk1_1112, kk1_1113, lk_1111, lk_1112, \
                         lk_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_18 * kk0_1111[k]
                    - f_19 * kk1_1111[k]
                    + pa_x[k] * lk_1111[k];

        t_1112[k] = f_18 * kk0_1112[k]
                    - f_19 * kk1_1112[k]
                    + pa_x[k] * lk_1112[k];

        t_1113[k] = f_18 * kk0_1113[k]
                    - f_19 * kk1_1113[k]
                    + pa_x[k] * lk_1113[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, pa_x, pa_y, pb_y, kk0_648, kk0_1115, kk1_648, \
                         kk1_1115, li_671, lk_864, lk_1115, mi_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_15 * li_671[k]
                    + pb_y[k] * mi_867[k];

        t_1115[k] = f_18 * kk0_1115[k]
                    - f_19 * kk1_1115[k]
                    + pa_x[k] * lk_1115[k];

        t_1116[k] = f_27 * kk0_648[k]
                    - f_28 * kk1_648[k]
                    + pa_y[k] * lk_864[k];
    }

#pragma omp simd aligned(t_1117, t_1118, t_1119, t_1120, pa_z, pb_y, pb_z, kk0_579, kk1_579, \
                         li_644, li_672, li_674, lk_831, mi_868, \
                         mi_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = f_14 * li_672[k]
                    + pb_y[k] * mi_868[k];

        t_1118[k] = f_13 * li_644[k]
                    + pb_z[k] * mi_868[k];

        t_1119[k] = f_22 * kk0_579[k]
                    - f_23 * kk1_579[k]
                    + pa_z[k] * lk_831[k];

        t_1120[k] = f_14 * li_674[k]
                    + pb_y[k] * mi_870[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pa_y, pa_z, pb_z, kk0_582, kk0_653, kk1_582, \
                         kk1_653, li_647, lk_834, lk_869, mi_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_27 * kk0_653[k]
                    - f_28 * kk1_653[k]
                    + pa_y[k] * lk_869[k];

        t_1122[k] = f_22 * kk0_582[k]
                    - f_23 * kk1_582[k]
                    + pa_z[k] * lk_834[k];

        t_1123[k] = f_13 * li_647[k]
                    + pb_z[k] * mi_871[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, pa_y, pa_z, pb_y, kk0_586, kk0_657, kk1_586, \
                         kk1_657, li_677, lk_838, lk_873, mi_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_14 * li_677[k]
                    + pb_y[k] * mi_873[k];

        t_1125[k] = f_27 * kk0_657[k]
                    - f_28 * kk1_657[k]
                    + pa_y[k] * lk_873[k];

        t_1126[k] = f_22 * kk0_586[k]
                    - f_23 * kk1_586[k]
                    + pa_z[k] * lk_838[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pb_x, pb_y, pb_z, li_650, li_681, li_880, \
                         mh0_663, mh1_663, mi_874, mi_877, mi_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_13 * li_650[k]
                    + pb_z[k] * mi_874[k];

        t_1128[k] = f_12 * li_880[k]
                    + f_5 * mh0_663[k]
                    - f_6 * mh1_663[k]
                    + pb_x[k] * mi_880[k];

        t_1129[k] = f_14 * li_681[k]
                    + pb_y[k] * mi_877[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pa_y, pa_z, pb_z, kk0_591, kk0_662, kk1_591, \
                         kk1_662, li_654, lk_843, lk_878, mi_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = f_27 * kk0_662[k]
                    - f_28 * kk1_662[k]
                    + pa_y[k] * lk_878[k];

        t_1131[k] = f_22 * kk0_591[k]
                    - f_23 * kk1_591[k]
                    + pa_z[k] * lk_843[k];

        t_1132[k] = f_13 * li_654[k]
                    + pb_z[k] * mi_878[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, pb_x, pb_y, li_686, li_885, li_886, mh0_668, \
                         mh0_669, mh1_668, mh1_669, mi_882, mi_885, \
                         mi_886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_12 * li_885[k]
                    + f_3 * mh0_668[k]
                    - f_4 * mh1_668[k]
                    + pb_x[k] * mi_885[k];

        t_1134[k] = f_12 * li_886[k]
                    + f_3 * mh0_669[k]
                    - f_4 * mh1_669[k]
                    + pb_x[k] * mi_886[k];

        t_1135[k] = f_14 * li_686[k]
                    + pb_y[k] * mi_882[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, t_1139, pa_y, pb_x, kk0_668, kk1_668, li_889, \
                         li_890, li_891, lk_884, mi_889, mi_890, \
                         mi_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_27 * kk0_668[k]
                    - f_28 * kk1_668[k]
                    + pa_y[k] * lk_884[k];

        t_1137[k] = f_12 * li_889[k]
                    + pb_x[k] * mi_889[k];

        t_1138[k] = f_12 * li_890[k]
                    + pb_x[k] * mi_890[k];

        t_1139[k] = f_12 * li_891[k]
                    + pb_x[k] * mi_891[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, t_1143, pb_x, li_892, li_893, li_894, li_895, \
                         mi_892, mi_893, mi_894, mi_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_12 * li_892[k]
                    + pb_x[k] * mi_892[k];

        t_1141[k] = f_12 * li_893[k]
                    + pb_x[k] * mi_893[k];

        t_1142[k] = f_12 * li_894[k]
                    + pb_x[k] * mi_894[k];

        t_1143[k] = f_12 * li_895[k]
                    + pb_x[k] * mi_895[k];
    }

#pragma omp simd aligned(t_1144, t_1145, t_1146, pa_x, pb_z, kk0_1144, kk0_1146, kk1_1144, \
                         kk1_1146, li_665, lk_1144, lk_1146, mi_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1144[k] = f_18 * kk0_1144[k]
                    - f_19 * kk1_1144[k]
                    + pa_x[k] * lk_1144[k];

        t_1145[k] = f_13 * li_665[k]
                    + pb_z[k] * mi_889[k];

        t_1146[k] = f_18 * kk0_1146[k]
                    - f_19 * kk1_1146[k]
                    + pa_x[k] * lk_1146[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, pa_x, kk0_1147, kk0_1148, kk0_1149, kk1_1147, \
                         kk1_1148, kk1_1149, lk_1147, lk_1148, \
                         lk_1149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = f_18 * kk0_1147[k]
                    - f_19 * kk1_1147[k]
                    + pa_x[k] * lk_1147[k];

        t_1148[k] = f_18 * kk0_1148[k]
                    - f_19 * kk1_1148[k]
                    + pa_x[k] * lk_1148[k];

        t_1149[k] = f_18 * kk0_1149[k]
                    - f_19 * kk1_1149[k]
                    + pa_x[k] * lk_1149[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, pa_x, pa_y, pb_y, kk0_684, kk0_1151, kk1_684, \
                         kk1_1151, li_699, lk_900, lk_1151, mi_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = f_14 * li_699[k]
                    + pb_y[k] * mi_895[k];

        t_1151[k] = f_18 * kk0_1151[k]
                    - f_19 * kk1_1151[k]
                    + pa_x[k] * lk_1151[k];

        t_1152[k] = f_22 * kk0_684[k]
                    - f_23 * kk1_684[k]
                    + pa_y[k] * lk_900[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, t_1156, pa_z, pb_y, pb_z, kk0_615, kk1_615, \
                         li_672, li_700, li_702, lk_867, mi_896, \
                         mi_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_13 * li_700[k]
                    + pb_y[k] * mi_896[k];

        t_1154[k] = f_14 * li_672[k]
                    + pb_z[k] * mi_896[k];

        t_1155[k] = f_27 * kk0_615[k]
                    - f_28 * kk1_615[k]
                    + pa_z[k] * lk_867[k];

        t_1156[k] = f_13 * li_702[k]
                    + pb_y[k] * mi_898[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece9(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kk0,
                                            const size_t kk1, const size_t li, const size_t lk,
                                            const size_t mh0, const size_t mh1, const size_t mi,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_17 = 3.5 / p;
    const auto f_18 = 0.5 / alpha;
    const auto f_19 = 0.5 * beta / (alpha * p);
    const auto f_20 = 3.0 / alpha;
    const auto f_21 = 3.0 * beta / (alpha * p);
    const auto f_22 = 1.0 / alpha;
    const auto f_23 = beta / (alpha * p);
    const auto f_24 = 3.0 / p;
    const auto f_27 = 1.5 / alpha;
    const auto f_28 = 1.5 * beta / (alpha * p);
    const auto f_29 = 2.0 / alpha;
    const auto f_30 = 2.0 * beta / (alpha * p);

    auto *t_1157 = buffer.data(target + 1157);
    auto *t_1158 = buffer.data(target + 1158);
    auto *t_1159 = buffer.data(target + 1159);
    auto *t_1160 = buffer.data(target + 1160);
    auto *t_1161 = buffer.data(target + 1161);
    auto *t_1162 = buffer.data(target + 1162);
    auto *t_1163 = buffer.data(target + 1163);
    auto *t_1164 = buffer.data(target + 1164);
    auto *t_1165 = buffer.data(target + 1165);
    auto *t_1166 = buffer.data(target + 1166);
    auto *t_1167 = buffer.data(target + 1167);
    auto *t_1168 = buffer.data(target + 1168);
    auto *t_1169 = buffer.data(target + 1169);
    auto *t_1170 = buffer.data(target + 1170);
    auto *t_1171 = buffer.data(target + 1171);
    auto *t_1172 = buffer.data(target + 1172);
    auto *t_1173 = buffer.data(target + 1173);
    auto *t_1174 = buffer.data(target + 1174);
    auto *t_1175 = buffer.data(target + 1175);
    auto *t_1176 = buffer.data(target + 1176);
    auto *t_1177 = buffer.data(target + 1177);
    auto *t_1178 = buffer.data(target + 1178);
    auto *t_1179 = buffer.data(target + 1179);
    auto *t_1180 = buffer.data(target + 1180);
    auto *t_1181 = buffer.data(target + 1181);
    auto *t_1182 = buffer.data(target + 1182);
    auto *t_1183 = buffer.data(target + 1183);
    auto *t_1184 = buffer.data(target + 1184);
    auto *t_1185 = buffer.data(target + 1185);
    auto *t_1186 = buffer.data(target + 1186);
    auto *t_1187 = buffer.data(target + 1187);
    auto *t_1188 = buffer.data(target + 1188);
    auto *t_1189 = buffer.data(target + 1189);
    auto *t_1190 = buffer.data(target + 1190);
    auto *t_1191 = buffer.data(target + 1191);
    auto *t_1192 = buffer.data(target + 1192);
    auto *t_1193 = buffer.data(target + 1193);
    auto *t_1194 = buffer.data(target + 1194);
    auto *t_1195 = buffer.data(target + 1195);
    auto *t_1196 = buffer.data(target + 1196);
    auto *t_1197 = buffer.data(target + 1197);
    auto *t_1198 = buffer.data(target + 1198);
    auto *t_1199 = buffer.data(target + 1199);
    auto *t_1200 = buffer.data(target + 1200);
    auto *t_1201 = buffer.data(target + 1201);
    auto *t_1202 = buffer.data(target + 1202);
    auto *t_1203 = buffer.data(target + 1203);
    auto *t_1204 = buffer.data(target + 1204);
    auto *t_1205 = buffer.data(target + 1205);
    auto *t_1206 = buffer.data(target + 1206);
    auto *t_1207 = buffer.data(target + 1207);
    auto *t_1208 = buffer.data(target + 1208);
    auto *t_1209 = buffer.data(target + 1209);
    auto *t_1210 = buffer.data(target + 1210);
    auto *t_1211 = buffer.data(target + 1211);
    auto *t_1212 = buffer.data(target + 1212);
    auto *t_1213 = buffer.data(target + 1213);
    auto *t_1214 = buffer.data(target + 1214);
    auto *t_1215 = buffer.data(target + 1215);
    auto *t_1216 = buffer.data(target + 1216);
    auto *t_1217 = buffer.data(target + 1217);
    auto *t_1218 = buffer.data(target + 1218);
    auto *t_1219 = buffer.data(target + 1219);
    auto *t_1220 = buffer.data(target + 1220);
    auto *t_1221 = buffer.data(target + 1221);
    auto *t_1222 = buffer.data(target + 1222);
    auto *t_1223 = buffer.data(target + 1223);
    auto *t_1224 = buffer.data(target + 1224);
    auto *t_1225 = buffer.data(target + 1225);
    auto *t_1226 = buffer.data(target + 1226);
    auto *t_1227 = buffer.data(target + 1227);
    auto *t_1228 = buffer.data(target + 1228);
    auto *t_1229 = buffer.data(target + 1229);
    auto *t_1230 = buffer.data(target + 1230);
    auto *t_1231 = buffer.data(target + 1231);
    auto *t_1232 = buffer.data(target + 1232);
    auto *t_1233 = buffer.data(target + 1233);
    auto *t_1234 = buffer.data(target + 1234);
    auto *t_1235 = buffer.data(target + 1235);
    auto *t_1236 = buffer.data(target + 1236);
    auto *t_1237 = buffer.data(target + 1237);
    auto *t_1238 = buffer.data(target + 1238);
    auto *t_1239 = buffer.data(target + 1239);
    auto *t_1240 = buffer.data(target + 1240);
    auto *t_1241 = buffer.data(target + 1241);
    auto *t_1242 = buffer.data(target + 1242);
    auto *t_1243 = buffer.data(target + 1243);
    auto *t_1244 = buffer.data(target + 1244);
    auto *t_1245 = buffer.data(target + 1245);
    auto *t_1246 = buffer.data(target + 1246);
    auto *t_1247 = buffer.data(target + 1247);
    auto *t_1248 = buffer.data(target + 1248);
    auto *t_1249 = buffer.data(target + 1249);
    auto *t_1250 = buffer.data(target + 1250);
    auto *t_1251 = buffer.data(target + 1251);
    auto *t_1252 = buffer.data(target + 1252);
    auto *t_1253 = buffer.data(target + 1253);
    auto *t_1254 = buffer.data(target + 1254);
    auto *t_1255 = buffer.data(target + 1255);
    auto *t_1256 = buffer.data(target + 1256);
    auto *t_1257 = buffer.data(target + 1257);
    auto *t_1258 = buffer.data(target + 1258);
    auto *t_1259 = buffer.data(target + 1259);
    auto *t_1260 = buffer.data(target + 1260);
    auto *t_1261 = buffer.data(target + 1261);
    auto *t_1262 = buffer.data(target + 1262);
    auto *t_1263 = buffer.data(target + 1263);
    auto *t_1264 = buffer.data(target + 1264);
    auto *t_1265 = buffer.data(target + 1265);
    auto *t_1266 = buffer.data(target + 1266);
    auto *t_1267 = buffer.data(target + 1267);
    auto *t_1268 = buffer.data(target + 1268);
    auto *t_1269 = buffer.data(target + 1269);
    auto *t_1270 = buffer.data(target + 1270);
    auto *t_1271 = buffer.data(target + 1271);
    auto *t_1272 = buffer.data(target + 1272);
    auto *t_1273 = buffer.data(target + 1273);
    auto *t_1274 = buffer.data(target + 1274);
    auto *t_1275 = buffer.data(target + 1275);
    auto *t_1276 = buffer.data(target + 1276);
    auto *t_1277 = buffer.data(target + 1277);
    auto *t_1278 = buffer.data(target + 1278);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kk0_618 = buffer.data(kk0 + 618);
    const auto *kk0_622 = buffer.data(kk0 + 622);
    const auto *kk0_627 = buffer.data(kk0 + 627);
    const auto *kk0_651 = buffer.data(kk0 + 651);
    const auto *kk0_654 = buffer.data(kk0 + 654);
    const auto *kk0_658 = buffer.data(kk0 + 658);
    const auto *kk0_663 = buffer.data(kk0 + 663);
    const auto *kk0_689 = buffer.data(kk0 + 689);
    const auto *kk0_693 = buffer.data(kk0 + 693);
    const auto *kk0_698 = buffer.data(kk0 + 698);
    const auto *kk0_704 = buffer.data(kk0 + 704);
    const auto *kk0_720 = buffer.data(kk0 + 720);
    const auto *kk0_725 = buffer.data(kk0 + 725);
    const auto *kk0_729 = buffer.data(kk0 + 729);
    const auto *kk0_734 = buffer.data(kk0 + 734);
    const auto *kk0_740 = buffer.data(kk0 + 740);
    const auto *kk0_1180 = buffer.data(kk0 + 1180);
    const auto *kk0_1182 = buffer.data(kk0 + 1182);
    const auto *kk0_1183 = buffer.data(kk0 + 1183);
    const auto *kk0_1184 = buffer.data(kk0 + 1184);
    const auto *kk0_1185 = buffer.data(kk0 + 1185);
    const auto *kk0_1187 = buffer.data(kk0 + 1187);
    const auto *kk0_1216 = buffer.data(kk0 + 1216);
    const auto *kk0_1218 = buffer.data(kk0 + 1218);
    const auto *kk0_1219 = buffer.data(kk0 + 1219);
    const auto *kk0_1220 = buffer.data(kk0 + 1220);
    const auto *kk0_1221 = buffer.data(kk0 + 1221);
    const auto *kk0_1223 = buffer.data(kk0 + 1223);

    const auto *kk1_618 = buffer.data(kk1 + 618);
    const auto *kk1_622 = buffer.data(kk1 + 622);
    const auto *kk1_627 = buffer.data(kk1 + 627);
    const auto *kk1_651 = buffer.data(kk1 + 651);
    const auto *kk1_654 = buffer.data(kk1 + 654);
    const auto *kk1_658 = buffer.data(kk1 + 658);
    const auto *kk1_663 = buffer.data(kk1 + 663);
    const auto *kk1_689 = buffer.data(kk1 + 689);
    const auto *kk1_693 = buffer.data(kk1 + 693);
    const auto *kk1_698 = buffer.data(kk1 + 698);
    const auto *kk1_704 = buffer.data(kk1 + 704);
    const auto *kk1_720 = buffer.data(kk1 + 720);
    const auto *kk1_725 = buffer.data(kk1 + 725);
    const auto *kk1_729 = buffer.data(kk1 + 729);
    const auto *kk1_734 = buffer.data(kk1 + 734);
    const auto *kk1_740 = buffer.data(kk1 + 740);
    const auto *kk1_1180 = buffer.data(kk1 + 1180);
    const auto *kk1_1182 = buffer.data(kk1 + 1182);
    const auto *kk1_1183 = buffer.data(kk1 + 1183);
    const auto *kk1_1184 = buffer.data(kk1 + 1184);
    const auto *kk1_1185 = buffer.data(kk1 + 1185);
    const auto *kk1_1187 = buffer.data(kk1 + 1187);
    const auto *kk1_1216 = buffer.data(kk1 + 1216);
    const auto *kk1_1218 = buffer.data(kk1 + 1218);
    const auto *kk1_1219 = buffer.data(kk1 + 1219);
    const auto *kk1_1220 = buffer.data(kk1 + 1220);
    const auto *kk1_1221 = buffer.data(kk1 + 1221);
    const auto *kk1_1223 = buffer.data(kk1 + 1223);

    const auto *li_675 = buffer.data(li + 675);
    const auto *li_678 = buffer.data(li + 678);
    const auto *li_682 = buffer.data(li + 682);
    const auto *li_693 = buffer.data(li + 693);
    const auto *li_700 = buffer.data(li + 700);
    const auto *li_703 = buffer.data(li + 703);
    const auto *li_705 = buffer.data(li + 705);
    const auto *li_706 = buffer.data(li + 706);
    const auto *li_709 = buffer.data(li + 709);
    const auto *li_710 = buffer.data(li + 710);
    const auto *li_714 = buffer.data(li + 714);
    const auto *li_721 = buffer.data(li + 721);
    const auto *li_727 = buffer.data(li + 727);
    const auto *li_728 = buffer.data(li + 728);
    const auto *li_730 = buffer.data(li + 730);
    const auto *li_731 = buffer.data(li + 731);
    const auto *li_733 = buffer.data(li + 733);
    const auto *li_734 = buffer.data(li + 734);
    const auto *li_737 = buffer.data(li + 737);
    const auto *li_738 = buffer.data(li + 738);
    const auto *li_742 = buffer.data(li + 742);
    const auto *li_749 = buffer.data(li + 749);
    const auto *li_755 = buffer.data(li + 755);
    const auto *li_756 = buffer.data(li + 756);
    const auto *li_757 = buffer.data(li + 757);
    const auto *li_758 = buffer.data(li + 758);
    const auto *li_759 = buffer.data(li + 759);
    const auto *li_761 = buffer.data(li + 761);
    const auto *li_762 = buffer.data(li + 762);
    const auto *li_764 = buffer.data(li + 764);
    const auto *li_765 = buffer.data(li + 765);
    const auto *li_766 = buffer.data(li + 766);
    const auto *li_768 = buffer.data(li + 768);
    const auto *li_769 = buffer.data(li + 769);
    const auto *li_770 = buffer.data(li + 770);
    const auto *li_777 = buffer.data(li + 777);
    const auto *li_779 = buffer.data(li + 779);
    const auto *li_780 = buffer.data(li + 780);
    const auto *li_781 = buffer.data(li + 781);
    const auto *li_782 = buffer.data(li + 782);
    const auto *li_783 = buffer.data(li + 783);
    const auto *li_908 = buffer.data(li + 908);
    const auto *li_913 = buffer.data(li + 913);
    const auto *li_914 = buffer.data(li + 914);
    const auto *li_917 = buffer.data(li + 917);
    const auto *li_918 = buffer.data(li + 918);
    const auto *li_919 = buffer.data(li + 919);
    const auto *li_920 = buffer.data(li + 920);
    const auto *li_921 = buffer.data(li + 921);
    const auto *li_922 = buffer.data(li + 922);
    const auto *li_923 = buffer.data(li + 923);
    const auto *li_936 = buffer.data(li + 936);
    const auto *li_941 = buffer.data(li + 941);
    const auto *li_942 = buffer.data(li + 942);
    const auto *li_945 = buffer.data(li + 945);
    const auto *li_946 = buffer.data(li + 946);
    const auto *li_947 = buffer.data(li + 947);
    const auto *li_948 = buffer.data(li + 948);
    const auto *li_949 = buffer.data(li + 949);
    const auto *li_950 = buffer.data(li + 950);
    const auto *li_951 = buffer.data(li + 951);
    const auto *li_973 = buffer.data(li + 973);
    const auto *li_974 = buffer.data(li + 974);
    const auto *li_975 = buffer.data(li + 975);
    const auto *li_976 = buffer.data(li + 976);
    const auto *li_977 = buffer.data(li + 977);
    const auto *li_978 = buffer.data(li + 978);
    const auto *li_985 = buffer.data(li + 985);
    const auto *li_989 = buffer.data(li + 989);
    const auto *li_994 = buffer.data(li + 994);

    const auto *lk_870 = buffer.data(lk + 870);
    const auto *lk_874 = buffer.data(lk + 874);
    const auto *lk_879 = buffer.data(lk + 879);
    const auto *lk_903 = buffer.data(lk + 903);
    const auto *lk_905 = buffer.data(lk + 905);
    const auto *lk_906 = buffer.data(lk + 906);
    const auto *lk_909 = buffer.data(lk + 909);
    const auto *lk_910 = buffer.data(lk + 910);
    const auto *lk_914 = buffer.data(lk + 914);
    const auto *lk_915 = buffer.data(lk + 915);
    const auto *lk_920 = buffer.data(lk + 920);
    const auto *lk_936 = buffer.data(lk + 936);
    const auto *lk_941 = buffer.data(lk + 941);
    const auto *lk_945 = buffer.data(lk + 945);
    const auto *lk_950 = buffer.data(lk + 950);
    const auto *lk_956 = buffer.data(lk + 956);
    const auto *lk_972 = buffer.data(lk + 972);
    const auto *lk_974 = buffer.data(lk + 974);
    const auto *lk_975 = buffer.data(lk + 975);
    const auto *lk_977 = buffer.data(lk + 977);
    const auto *lk_978 = buffer.data(lk + 978);
    const auto *lk_981 = buffer.data(lk + 981);
    const auto *lk_982 = buffer.data(lk + 982);
    const auto *lk_984 = buffer.data(lk + 984);
    const auto *lk_986 = buffer.data(lk + 986);
    const auto *lk_987 = buffer.data(lk + 987);
    const auto *lk_989 = buffer.data(lk + 989);
    const auto *lk_990 = buffer.data(lk + 990);
    const auto *lk_992 = buffer.data(lk + 992);
    const auto *lk_999 = buffer.data(lk + 999);
    const auto *lk_1000 = buffer.data(lk + 1000);
    const auto *lk_1002 = buffer.data(lk + 1002);
    const auto *lk_1003 = buffer.data(lk + 1003);
    const auto *lk_1004 = buffer.data(lk + 1004);
    const auto *lk_1005 = buffer.data(lk + 1005);
    const auto *lk_1007 = buffer.data(lk + 1007);
    const auto *lk_1180 = buffer.data(lk + 1180);
    const auto *lk_1182 = buffer.data(lk + 1182);
    const auto *lk_1183 = buffer.data(lk + 1183);
    const auto *lk_1184 = buffer.data(lk + 1184);
    const auto *lk_1185 = buffer.data(lk + 1185);
    const auto *lk_1187 = buffer.data(lk + 1187);
    const auto *lk_1216 = buffer.data(lk + 1216);
    const auto *lk_1218 = buffer.data(lk + 1218);
    const auto *lk_1219 = buffer.data(lk + 1219);
    const auto *lk_1220 = buffer.data(lk + 1220);
    const auto *lk_1221 = buffer.data(lk + 1221);
    const auto *lk_1223 = buffer.data(lk + 1223);

    const auto *mh0_684 = buffer.data(mh0 + 684);
    const auto *mh0_689 = buffer.data(mh0 + 689);
    const auto *mh0_690 = buffer.data(mh0 + 690);
    const auto *mh0_705 = buffer.data(mh0 + 705);
    const auto *mh0_710 = buffer.data(mh0 + 710);
    const auto *mh0_711 = buffer.data(mh0 + 711);
    const auto *mh0_735 = buffer.data(mh0 + 735);
    const auto *mh0_736 = buffer.data(mh0 + 736);
    const auto *mh0_738 = buffer.data(mh0 + 738);
    const auto *mh0_740 = buffer.data(mh0 + 740);
    const auto *mh0_741 = buffer.data(mh0 + 741);
    const auto *mh0_743 = buffer.data(mh0 + 743);
    const auto *mh0_744 = buffer.data(mh0 + 744);
    const auto *mh0_749 = buffer.data(mh0 + 749);

    const auto *mh1_684 = buffer.data(mh1 + 684);
    const auto *mh1_689 = buffer.data(mh1 + 689);
    const auto *mh1_690 = buffer.data(mh1 + 690);
    const auto *mh1_705 = buffer.data(mh1 + 705);
    const auto *mh1_710 = buffer.data(mh1 + 710);
    const auto *mh1_711 = buffer.data(mh1 + 711);
    const auto *mh1_735 = buffer.data(mh1 + 735);
    const auto *mh1_736 = buffer.data(mh1 + 736);
    const auto *mh1_738 = buffer.data(mh1 + 738);
    const auto *mh1_740 = buffer.data(mh1 + 740);
    const auto *mh1_741 = buffer.data(mh1 + 741);
    const auto *mh1_743 = buffer.data(mh1 + 743);
    const auto *mh1_744 = buffer.data(mh1 + 744);
    const auto *mh1_749 = buffer.data(mh1 + 749);

    const auto *mi_899 = buffer.data(mi + 899);
    const auto *mi_901 = buffer.data(mi + 901);
    const auto *mi_902 = buffer.data(mi + 902);
    const auto *mi_905 = buffer.data(mi + 905);
    const auto *mi_906 = buffer.data(mi + 906);
    const auto *mi_908 = buffer.data(mi + 908);
    const auto *mi_910 = buffer.data(mi + 910);
    const auto *mi_913 = buffer.data(mi + 913);
    const auto *mi_914 = buffer.data(mi + 914);
    const auto *mi_917 = buffer.data(mi + 917);
    const auto *mi_918 = buffer.data(mi + 918);
    const auto *mi_919 = buffer.data(mi + 919);
    const auto *mi_920 = buffer.data(mi + 920);
    const auto *mi_921 = buffer.data(mi + 921);
    const auto *mi_922 = buffer.data(mi + 922);
    const auto *mi_923 = buffer.data(mi + 923);
    const auto *mi_924 = buffer.data(mi + 924);
    const auto *mi_926 = buffer.data(mi + 926);
    const auto *mi_927 = buffer.data(mi + 927);
    const auto *mi_929 = buffer.data(mi + 929);
    const auto *mi_930 = buffer.data(mi + 930);
    const auto *mi_933 = buffer.data(mi + 933);
    const auto *mi_934 = buffer.data(mi + 934);
    const auto *mi_936 = buffer.data(mi + 936);
    const auto *mi_938 = buffer.data(mi + 938);
    const auto *mi_941 = buffer.data(mi + 941);
    const auto *mi_942 = buffer.data(mi + 942);
    const auto *mi_945 = buffer.data(mi + 945);
    const auto *mi_946 = buffer.data(mi + 946);
    const auto *mi_947 = buffer.data(mi + 947);
    const auto *mi_948 = buffer.data(mi + 948);
    const auto *mi_949 = buffer.data(mi + 949);
    const auto *mi_950 = buffer.data(mi + 950);
    const auto *mi_951 = buffer.data(mi + 951);
    const auto *mi_952 = buffer.data(mi + 952);
    const auto *mi_954 = buffer.data(mi + 954);
    const auto *mi_955 = buffer.data(mi + 955);
    const auto *mi_957 = buffer.data(mi + 957);
    const auto *mi_958 = buffer.data(mi + 958);
    const auto *mi_961 = buffer.data(mi + 961);
    const auto *mi_962 = buffer.data(mi + 962);
    const auto *mi_966 = buffer.data(mi + 966);
    const auto *mi_973 = buffer.data(mi + 973);
    const auto *mi_974 = buffer.data(mi + 974);
    const auto *mi_975 = buffer.data(mi + 975);
    const auto *mi_976 = buffer.data(mi + 976);
    const auto *mi_977 = buffer.data(mi + 977);
    const auto *mi_978 = buffer.data(mi + 978);
    const auto *mi_979 = buffer.data(mi + 979);
    const auto *mi_980 = buffer.data(mi + 980);
    const auto *mi_981 = buffer.data(mi + 981);
    const auto *mi_982 = buffer.data(mi + 982);
    const auto *mi_983 = buffer.data(mi + 983);
    const auto *mi_985 = buffer.data(mi + 985);
    const auto *mi_986 = buffer.data(mi + 986);
    const auto *mi_988 = buffer.data(mi + 988);
    const auto *mi_989 = buffer.data(mi + 989);
    const auto *mi_990 = buffer.data(mi + 990);
    const auto *mi_992 = buffer.data(mi + 992);
    const auto *mi_993 = buffer.data(mi + 993);
    const auto *mi_994 = buffer.data(mi + 994);

#pragma omp simd aligned(t_1157, t_1158, t_1159, pa_y, pa_z, pb_z, kk0_618, kk0_689, kk1_618, \
                         kk1_689, li_675, lk_870, lk_905, mi_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_22 * kk0_689[k]
                    - f_23 * kk1_689[k]
                    + pa_y[k] * lk_905[k];

        t_1158[k] = f_27 * kk0_618[k]
                    - f_28 * kk1_618[k]
                    + pa_z[k] * lk_870[k];

        t_1159[k] = f_14 * li_675[k]
                    + pb_z[k] * mi_899[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, pa_y, pa_z, pb_y, kk0_622, kk0_693, kk1_622, \
                         kk1_693, li_705, lk_874, lk_909, mi_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_13 * li_705[k]
                    + pb_y[k] * mi_901[k];

        t_1161[k] = f_22 * kk0_693[k]
                    - f_23 * kk1_693[k]
                    + pa_y[k] * lk_909[k];

        t_1162[k] = f_27 * kk0_622[k]
                    - f_28 * kk1_622[k]
                    + pa_z[k] * lk_874[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, pb_x, pb_y, pb_z, li_678, li_709, li_908, \
                         mh0_684, mh1_684, mi_902, mi_905, mi_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_14 * li_678[k]
                    + pb_z[k] * mi_902[k];

        t_1164[k] = f_12 * li_908[k]
                    + f_5 * mh0_684[k]
                    - f_6 * mh1_684[k]
                    + pb_x[k] * mi_908[k];

        t_1165[k] = f_13 * li_709[k]
                    + pb_y[k] * mi_905[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, pa_y, pa_z, pb_z, kk0_627, kk0_698, kk1_627, \
                         kk1_698, li_682, lk_879, lk_914, mi_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_22 * kk0_698[k]
                    - f_23 * kk1_698[k]
                    + pa_y[k] * lk_914[k];

        t_1167[k] = f_27 * kk0_627[k]
                    - f_28 * kk1_627[k]
                    + pa_z[k] * lk_879[k];

        t_1168[k] = f_14 * li_682[k]
                    + pb_z[k] * mi_906[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, pb_x, pb_y, li_714, li_913, li_914, mh0_689, \
                         mh0_690, mh1_689, mh1_690, mi_910, mi_913, \
                         mi_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_12 * li_913[k]
                    + f_3 * mh0_689[k]
                    - f_4 * mh1_689[k]
                    + pb_x[k] * mi_913[k];

        t_1170[k] = f_12 * li_914[k]
                    + f_3 * mh0_690[k]
                    - f_4 * mh1_690[k]
                    + pb_x[k] * mi_914[k];

        t_1171[k] = f_13 * li_714[k]
                    + pb_y[k] * mi_910[k];
    }

#pragma omp simd aligned(t_1172, t_1173, t_1174, t_1175, pa_y, pb_x, kk0_704, kk1_704, li_917, \
                         li_918, li_919, lk_920, mi_917, mi_918, \
                         mi_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1172[k] = f_22 * kk0_704[k]
                    - f_23 * kk1_704[k]
                    + pa_y[k] * lk_920[k];

        t_1173[k] = f_12 * li_917[k]
                    + pb_x[k] * mi_917[k];

        t_1174[k] = f_12 * li_918[k]
                    + pb_x[k] * mi_918[k];

        t_1175[k] = f_12 * li_919[k]
                    + pb_x[k] * mi_919[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, t_1179, pb_x, li_920, li_921, li_922, li_923, \
                         mi_920, mi_921, mi_922, mi_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = f_12 * li_920[k]
                    + pb_x[k] * mi_920[k];

        t_1177[k] = f_12 * li_921[k]
                    + pb_x[k] * mi_921[k];

        t_1178[k] = f_12 * li_922[k]
                    + pb_x[k] * mi_922[k];

        t_1179[k] = f_12 * li_923[k]
                    + pb_x[k] * mi_923[k];
    }

#pragma omp simd aligned(t_1180, t_1181, t_1182, pa_x, pb_z, kk0_1180, kk0_1182, kk1_1180, \
                         kk1_1182, li_693, lk_1180, lk_1182, mi_917 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1180[k] = f_18 * kk0_1180[k]
                    - f_19 * kk1_1180[k]
                    + pa_x[k] * lk_1180[k];

        t_1181[k] = f_14 * li_693[k]
                    + pb_z[k] * mi_917[k];

        t_1182[k] = f_18 * kk0_1182[k]
                    - f_19 * kk1_1182[k]
                    + pa_x[k] * lk_1182[k];
    }

#pragma omp simd aligned(t_1183, t_1184, t_1185, pa_x, kk0_1183, kk0_1184, kk0_1185, kk1_1183, \
                         kk1_1184, kk1_1185, lk_1183, lk_1184, \
                         lk_1185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1183[k] = f_18 * kk0_1183[k]
                    - f_19 * kk1_1183[k]
                    + pa_x[k] * lk_1183[k];

        t_1184[k] = f_18 * kk0_1184[k]
                    - f_19 * kk1_1184[k]
                    + pa_x[k] * lk_1184[k];

        t_1185[k] = f_18 * kk0_1185[k]
                    - f_19 * kk1_1185[k]
                    + pa_x[k] * lk_1185[k];
    }

#pragma omp simd aligned(t_1186, t_1187, t_1188, pa_x, pa_y, pb_y, kk0_720, kk0_1187, kk1_720, \
                         kk1_1187, li_727, lk_936, lk_1187, mi_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1186[k] = f_13 * li_727[k]
                    + pb_y[k] * mi_923[k];

        t_1187[k] = f_18 * kk0_1187[k]
                    - f_19 * kk1_1187[k]
                    + pa_x[k] * lk_1187[k];

        t_1188[k] = f_18 * kk0_720[k]
                    - f_19 * kk1_720[k]
                    + pa_y[k] * lk_936[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, t_1192, pa_z, pb_y, pb_z, kk0_651, kk1_651, \
                         li_700, li_728, li_730, lk_903, mi_924, \
                         mi_926 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_12 * li_728[k]
                    + pb_y[k] * mi_924[k];

        t_1190[k] = f_15 * li_700[k]
                    + pb_z[k] * mi_924[k];

        t_1191[k] = f_29 * kk0_651[k]
                    - f_30 * kk1_651[k]
                    + pa_z[k] * lk_903[k];

        t_1192[k] = f_12 * li_730[k]
                    + pb_y[k] * mi_926[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, pa_y, pa_z, pb_z, kk0_654, kk0_725, kk1_654, \
                         kk1_725, li_703, lk_906, lk_941, mi_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = f_18 * kk0_725[k]
                    - f_19 * kk1_725[k]
                    + pa_y[k] * lk_941[k];

        t_1194[k] = f_29 * kk0_654[k]
                    - f_30 * kk1_654[k]
                    + pa_z[k] * lk_906[k];

        t_1195[k] = f_15 * li_703[k]
                    + pb_z[k] * mi_927[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, pa_y, pa_z, pb_y, kk0_658, kk0_729, kk1_658, \
                         kk1_729, li_733, lk_910, lk_945, mi_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_12 * li_733[k]
                    + pb_y[k] * mi_929[k];

        t_1197[k] = f_18 * kk0_729[k]
                    - f_19 * kk1_729[k]
                    + pa_y[k] * lk_945[k];

        t_1198[k] = f_29 * kk0_658[k]
                    - f_30 * kk1_658[k]
                    + pa_z[k] * lk_910[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pb_x, pb_y, pb_z, li_706, li_737, li_936, \
                         mh0_705, mh1_705, mi_930, mi_933, mi_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_15 * li_706[k]
                    + pb_z[k] * mi_930[k];

        t_1200[k] = f_12 * li_936[k]
                    + f_5 * mh0_705[k]
                    - f_6 * mh1_705[k]
                    + pb_x[k] * mi_936[k];

        t_1201[k] = f_12 * li_737[k]
                    + pb_y[k] * mi_933[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, pa_y, pa_z, pb_z, kk0_663, kk0_734, kk1_663, \
                         kk1_734, li_710, lk_915, lk_950, mi_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_18 * kk0_734[k]
                    - f_19 * kk1_734[k]
                    + pa_y[k] * lk_950[k];

        t_1203[k] = f_29 * kk0_663[k]
                    - f_30 * kk1_663[k]
                    + pa_z[k] * lk_915[k];

        t_1204[k] = f_15 * li_710[k]
                    + pb_z[k] * mi_934[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, pb_x, pb_y, li_742, li_941, li_942, mh0_710, \
                         mh0_711, mh1_710, mh1_711, mi_938, mi_941, \
                         mi_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = f_12 * li_941[k]
                    + f_3 * mh0_710[k]
                    - f_4 * mh1_710[k]
                    + pb_x[k] * mi_941[k];

        t_1206[k] = f_12 * li_942[k]
                    + f_3 * mh0_711[k]
                    - f_4 * mh1_711[k]
                    + pb_x[k] * mi_942[k];

        t_1207[k] = f_12 * li_742[k]
                    + pb_y[k] * mi_938[k];
    }

#pragma omp simd aligned(t_1208, t_1209, t_1210, t_1211, pa_y, pb_x, kk0_740, kk1_740, li_945, \
                         li_946, li_947, lk_956, mi_945, mi_946, \
                         mi_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1208[k] = f_18 * kk0_740[k]
                    - f_19 * kk1_740[k]
                    + pa_y[k] * lk_956[k];

        t_1209[k] = f_12 * li_945[k]
                    + pb_x[k] * mi_945[k];

        t_1210[k] = f_12 * li_946[k]
                    + pb_x[k] * mi_946[k];

        t_1211[k] = f_12 * li_947[k]
                    + pb_x[k] * mi_947[k];
    }

#pragma omp simd aligned(t_1212, t_1213, t_1214, t_1215, pb_x, li_948, li_949, li_950, li_951, \
                         mi_948, mi_949, mi_950, mi_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1212[k] = f_12 * li_948[k]
                    + pb_x[k] * mi_948[k];

        t_1213[k] = f_12 * li_949[k]
                    + pb_x[k] * mi_949[k];

        t_1214[k] = f_12 * li_950[k]
                    + pb_x[k] * mi_950[k];

        t_1215[k] = f_12 * li_951[k]
                    + pb_x[k] * mi_951[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, pa_x, pb_z, kk0_1216, kk0_1218, kk1_1216, \
                         kk1_1218, li_721, lk_1216, lk_1218, mi_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = f_18 * kk0_1216[k]
                    - f_19 * kk1_1216[k]
                    + pa_x[k] * lk_1216[k];

        t_1217[k] = f_15 * li_721[k]
                    + pb_z[k] * mi_945[k];

        t_1218[k] = f_18 * kk0_1218[k]
                    - f_19 * kk1_1218[k]
                    + pa_x[k] * lk_1218[k];
    }

#pragma omp simd aligned(t_1219, t_1220, t_1221, pa_x, kk0_1219, kk0_1220, kk0_1221, kk1_1219, \
                         kk1_1220, kk1_1221, lk_1219, lk_1220, \
                         lk_1221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1219[k] = f_18 * kk0_1219[k]
                    - f_19 * kk1_1219[k]
                    + pa_x[k] * lk_1219[k];

        t_1220[k] = f_18 * kk0_1220[k]
                    - f_19 * kk1_1220[k]
                    + pa_x[k] * lk_1220[k];

        t_1221[k] = f_18 * kk0_1221[k]
                    - f_19 * kk1_1221[k]
                    + pa_x[k] * lk_1221[k];
    }

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, pa_x, pa_y, pb_y, kk0_1223, kk1_1223, \
                         li_755, li_756, lk_972, lk_1223, mi_951, \
                         mi_952 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_12 * li_755[k]
                    + pb_y[k] * mi_951[k];

        t_1223[k] = f_18 * kk0_1223[k]
                    - f_19 * kk1_1223[k]
                    + pa_x[k] * lk_1223[k];

        t_1224[k] = pa_y[k] * lk_972[k];

        t_1225[k] = f_11 * li_756[k]
                    + pb_y[k] * mi_952[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, t_1229, t_1230, pa_y, pb_y, li_757, li_758, \
                         li_759, lk_974, lk_975, lk_977, lk_978, \
                         mi_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = pa_y[k] * lk_974[k];

        t_1227[k] = f_12 * li_757[k]
                    + pa_y[k] * lk_975[k];

        t_1228[k] = f_11 * li_758[k]
                    + pb_y[k] * mi_954[k];

        t_1229[k] = pa_y[k] * lk_977[k];

        t_1230[k] = f_13 * li_759[k]
                    + pa_y[k] * lk_978[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, t_1234, pa_y, pb_y, pb_z, li_731, li_761, \
                         li_762, lk_981, lk_982, mi_955, mi_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = f_24 * li_731[k]
                    + pb_z[k] * mi_955[k];

        t_1232[k] = f_11 * li_761[k]
                    + pb_y[k] * mi_957[k];

        t_1233[k] = pa_y[k] * lk_981[k];

        t_1234[k] = f_14 * li_762[k]
                    + pa_y[k] * lk_982[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, t_1238, pa_y, pb_y, pb_z, li_734, li_764, \
                         li_765, lk_984, lk_986, mi_958, mi_961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_24 * li_734[k]
                    + pb_z[k] * mi_958[k];

        t_1236[k] = f_12 * li_764[k]
                    + pa_y[k] * lk_984[k];

        t_1237[k] = f_11 * li_765[k]
                    + pb_y[k] * mi_961[k];

        t_1238[k] = pa_y[k] * lk_986[k];
    }

#pragma omp simd aligned(t_1239, t_1240, t_1241, t_1242, pa_y, pb_z, li_738, li_766, li_768, \
                         li_769, lk_987, lk_989, lk_990, mi_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1239[k] = f_15 * li_766[k]
                    + pa_y[k] * lk_987[k];

        t_1240[k] = f_24 * li_738[k]
                    + pb_z[k] * mi_962[k];

        t_1241[k] = f_13 * li_768[k]
                    + pa_y[k] * lk_989[k];

        t_1242[k] = f_12 * li_769[k]
                    + pa_y[k] * lk_990[k];
    }

#pragma omp simd aligned(t_1243, t_1244, t_1245, t_1246, pa_y, pb_x, pb_y, li_770, li_973, \
                         li_974, lk_992, mi_966, mi_973, mi_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1243[k] = f_11 * li_770[k]
                    + pb_y[k] * mi_966[k];

        t_1244[k] = pa_y[k] * lk_992[k];

        t_1245[k] = f_12 * li_973[k]
                    + pb_x[k] * mi_973[k];

        t_1246[k] = f_12 * li_974[k]
                    + pb_x[k] * mi_974[k];
    }

#pragma omp simd aligned(t_1247, t_1248, t_1249, t_1250, t_1251, pa_y, pb_x, li_975, li_976, \
                         li_977, li_978, lk_999, mi_975, mi_976, mi_977, \
                         mi_978 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1247[k] = f_12 * li_975[k]
                    + pb_x[k] * mi_975[k];

        t_1248[k] = f_12 * li_976[k]
                    + pb_x[k] * mi_976[k];

        t_1249[k] = f_12 * li_977[k]
                    + pb_x[k] * mi_977[k];

        t_1250[k] = f_12 * li_978[k]
                    + pb_x[k] * mi_978[k];

        t_1251[k] = pa_y[k] * lk_999[k];
    }

#pragma omp simd aligned(t_1252, t_1253, t_1254, t_1255, pa_y, pb_z, li_749, li_777, li_779, \
                         li_780, lk_1000, lk_1002, lk_1003, mi_973 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1252[k] = f_17 * li_777[k]
                    + pa_y[k] * lk_1000[k];

        t_1253[k] = f_24 * li_749[k]
                    + pb_z[k] * mi_973[k];

        t_1254[k] = f_15 * li_779[k]
                    + pa_y[k] * lk_1002[k];

        t_1255[k] = f_14 * li_780[k]
                    + pa_y[k] * lk_1003[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, t_1259, pa_y, pb_y, li_781, li_782, li_783, \
                         lk_1004, lk_1005, lk_1007, mi_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = f_13 * li_781[k]
                    + pa_y[k] * lk_1004[k];

        t_1257[k] = f_12 * li_782[k]
                    + pa_y[k] * lk_1005[k];

        t_1258[k] = f_11 * li_783[k]
                    + pb_y[k] * mi_979[k];

        t_1259[k] = pa_y[k] * lk_1007[k];
    }

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, pa_z, pb_y, pb_z, kk0_720, kk1_720, \
                         li_756, lk_972, mh0_735, mh1_735, mi_980, \
                         mi_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = f_20 * kk0_720[k]
                    - f_21 * kk1_720[k]
                    + pa_z[k] * lk_972[k];

        t_1261[k] = pb_y[k] * mi_980[k];

        t_1262[k] = f_17 * li_756[k]
                    + pb_z[k] * mi_980[k];

        t_1263[k] = f_3 * mh0_735[k]
                    - f_4 * mh1_735[k]
                    + pb_y[k] * mi_981[k];
    }

#pragma omp simd aligned(t_1264, t_1265, t_1266, t_1267, pb_x, pb_y, pb_z, li_759, li_985, \
                         mh0_736, mh0_740, mh1_736, mh1_740, mi_982, mi_983, \
                         mi_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1264[k] = pb_y[k] * mi_982[k];

        t_1265[k] = f_12 * li_985[k]
                    + f_9 * mh0_740[k]
                    - f_10 * mh1_740[k]
                    + pb_x[k] * mi_985[k];

        t_1266[k] = f_5 * mh0_736[k]
                    - f_6 * mh1_736[k]
                    + pb_y[k] * mi_983[k];

        t_1267[k] = f_17 * li_759[k]
                    + pb_z[k] * mi_983[k];
    }

#pragma omp simd aligned(t_1268, t_1269, t_1270, t_1271, pb_x, pb_y, pb_z, li_762, li_989, \
                         mh0_738, mh0_744, mh1_738, mh1_744, mi_985, mi_986, \
                         mi_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1268[k] = pb_y[k] * mi_985[k];

        t_1269[k] = f_12 * li_989[k]
                    + f_7 * mh0_744[k]
                    - f_8 * mh1_744[k]
                    + pb_x[k] * mi_989[k];

        t_1270[k] = f_7 * mh0_738[k]
                    - f_8 * mh1_738[k]
                    + pb_y[k] * mi_986[k];

        t_1271[k] = f_17 * li_762[k]
                    + pb_z[k] * mi_986[k];
    }

#pragma omp simd aligned(t_1272, t_1273, t_1274, pb_x, pb_y, li_994, mh0_740, mh0_749, \
                         mh1_740, mh1_749, mi_988, mi_989, mi_994 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1272[k] = f_3 * mh0_740[k]
                    - f_4 * mh1_740[k]
                    + pb_y[k] * mi_988[k];

        t_1273[k] = pb_y[k] * mi_989[k];

        t_1274[k] = f_12 * li_994[k]
                    + f_5 * mh0_749[k]
                    - f_6 * mh1_749[k]
                    + pb_x[k] * mi_994[k];
    }

#pragma omp simd aligned(t_1275, t_1276, t_1277, t_1278, pb_y, pb_z, li_766, mh0_741, mh0_743, \
                         mh0_744, mh1_741, mh1_743, mh1_744, mi_990, mi_992, \
                         mi_993 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1275[k] = f_9 * mh0_741[k]
                    - f_10 * mh1_741[k]
                    + pb_y[k] * mi_990[k];

        t_1276[k] = f_17 * li_766[k]
                    + pb_z[k] * mi_990[k];

        t_1277[k] = f_5 * mh0_743[k]
                    - f_6 * mh1_743[k]
                    + pb_y[k] * mi_992[k];

        t_1278[k] = f_3 * mh0_744[k]
                    - f_4 * mh1_744[k]
                    + pb_y[k] * mi_993[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece10(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t kk0,
                                             const size_t kk1, const size_t li, const size_t lk,
                                             const size_t mh0, const size_t mh1, const size_t mi,
                                             const size_t ncols, const double alpha,
                                             const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_16 = 4.0 / p;
    const auto f_17 = 3.5 / p;
    const auto f_18 = 0.5 / alpha;
    const auto f_19 = 0.5 * beta / (alpha * p);
    const auto f_24 = 3.0 / p;

    auto *t_1279 = buffer.data(target + 1279);
    auto *t_1280 = buffer.data(target + 1280);
    auto *t_1281 = buffer.data(target + 1281);
    auto *t_1282 = buffer.data(target + 1282);
    auto *t_1283 = buffer.data(target + 1283);
    auto *t_1284 = buffer.data(target + 1284);
    auto *t_1285 = buffer.data(target + 1285);
    auto *t_1286 = buffer.data(target + 1286);
    auto *t_1287 = buffer.data(target + 1287);
    auto *t_1288 = buffer.data(target + 1288);
    auto *t_1289 = buffer.data(target + 1289);
    auto *t_1290 = buffer.data(target + 1290);
    auto *t_1291 = buffer.data(target + 1291);
    auto *t_1292 = buffer.data(target + 1292);
    auto *t_1293 = buffer.data(target + 1293);
    auto *t_1294 = buffer.data(target + 1294);
    auto *t_1295 = buffer.data(target + 1295);
    auto *t_1296 = buffer.data(target + 1296);
    auto *t_1297 = buffer.data(target + 1297);
    auto *t_1298 = buffer.data(target + 1298);
    auto *t_1299 = buffer.data(target + 1299);
    auto *t_1300 = buffer.data(target + 1300);
    auto *t_1301 = buffer.data(target + 1301);
    auto *t_1302 = buffer.data(target + 1302);
    auto *t_1303 = buffer.data(target + 1303);
    auto *t_1304 = buffer.data(target + 1304);
    auto *t_1305 = buffer.data(target + 1305);
    auto *t_1306 = buffer.data(target + 1306);
    auto *t_1307 = buffer.data(target + 1307);
    auto *t_1308 = buffer.data(target + 1308);
    auto *t_1309 = buffer.data(target + 1309);
    auto *t_1310 = buffer.data(target + 1310);
    auto *t_1311 = buffer.data(target + 1311);
    auto *t_1312 = buffer.data(target + 1312);
    auto *t_1313 = buffer.data(target + 1313);
    auto *t_1314 = buffer.data(target + 1314);
    auto *t_1315 = buffer.data(target + 1315);
    auto *t_1316 = buffer.data(target + 1316);
    auto *t_1317 = buffer.data(target + 1317);
    auto *t_1318 = buffer.data(target + 1318);
    auto *t_1319 = buffer.data(target + 1319);
    auto *t_1320 = buffer.data(target + 1320);
    auto *t_1321 = buffer.data(target + 1321);
    auto *t_1322 = buffer.data(target + 1322);
    auto *t_1323 = buffer.data(target + 1323);
    auto *t_1324 = buffer.data(target + 1324);
    auto *t_1325 = buffer.data(target + 1325);
    auto *t_1326 = buffer.data(target + 1326);
    auto *t_1327 = buffer.data(target + 1327);
    auto *t_1328 = buffer.data(target + 1328);
    auto *t_1329 = buffer.data(target + 1329);
    auto *t_1330 = buffer.data(target + 1330);
    auto *t_1331 = buffer.data(target + 1331);
    auto *t_1332 = buffer.data(target + 1332);
    auto *t_1333 = buffer.data(target + 1333);
    auto *t_1334 = buffer.data(target + 1334);
    auto *t_1335 = buffer.data(target + 1335);
    auto *t_1336 = buffer.data(target + 1336);
    auto *t_1337 = buffer.data(target + 1337);
    auto *t_1338 = buffer.data(target + 1338);
    auto *t_1339 = buffer.data(target + 1339);
    auto *t_1340 = buffer.data(target + 1340);
    auto *t_1341 = buffer.data(target + 1341);
    auto *t_1342 = buffer.data(target + 1342);
    auto *t_1343 = buffer.data(target + 1343);
    auto *t_1344 = buffer.data(target + 1344);
    auto *t_1345 = buffer.data(target + 1345);
    auto *t_1346 = buffer.data(target + 1346);
    auto *t_1347 = buffer.data(target + 1347);
    auto *t_1348 = buffer.data(target + 1348);
    auto *t_1349 = buffer.data(target + 1349);
    auto *t_1350 = buffer.data(target + 1350);
    auto *t_1351 = buffer.data(target + 1351);
    auto *t_1352 = buffer.data(target + 1352);
    auto *t_1353 = buffer.data(target + 1353);
    auto *t_1354 = buffer.data(target + 1354);
    auto *t_1355 = buffer.data(target + 1355);
    auto *t_1356 = buffer.data(target + 1356);
    auto *t_1357 = buffer.data(target + 1357);
    auto *t_1358 = buffer.data(target + 1358);
    auto *t_1359 = buffer.data(target + 1359);
    auto *t_1360 = buffer.data(target + 1360);
    auto *t_1361 = buffer.data(target + 1361);
    auto *t_1362 = buffer.data(target + 1362);
    auto *t_1363 = buffer.data(target + 1363);
    auto *t_1364 = buffer.data(target + 1364);
    auto *t_1365 = buffer.data(target + 1365);
    auto *t_1366 = buffer.data(target + 1366);
    auto *t_1367 = buffer.data(target + 1367);
    auto *t_1368 = buffer.data(target + 1368);
    auto *t_1369 = buffer.data(target + 1369);
    auto *t_1370 = buffer.data(target + 1370);
    auto *t_1371 = buffer.data(target + 1371);
    auto *t_1372 = buffer.data(target + 1372);
    auto *t_1373 = buffer.data(target + 1373);
    auto *t_1374 = buffer.data(target + 1374);
    auto *t_1375 = buffer.data(target + 1375);
    auto *t_1376 = buffer.data(target + 1376);
    auto *t_1377 = buffer.data(target + 1377);
    auto *t_1378 = buffer.data(target + 1378);
    auto *t_1379 = buffer.data(target + 1379);
    auto *t_1380 = buffer.data(target + 1380);
    auto *t_1381 = buffer.data(target + 1381);
    auto *t_1382 = buffer.data(target + 1382);
    auto *t_1383 = buffer.data(target + 1383);
    auto *t_1384 = buffer.data(target + 1384);
    auto *t_1385 = buffer.data(target + 1385);
    auto *t_1386 = buffer.data(target + 1386);
    auto *t_1387 = buffer.data(target + 1387);
    auto *t_1388 = buffer.data(target + 1388);
    auto *t_1389 = buffer.data(target + 1389);
    auto *t_1390 = buffer.data(target + 1390);
    auto *t_1391 = buffer.data(target + 1391);
    auto *t_1392 = buffer.data(target + 1392);
    auto *t_1393 = buffer.data(target + 1393);
    auto *t_1394 = buffer.data(target + 1394);
    auto *t_1395 = buffer.data(target + 1395);
    auto *t_1396 = buffer.data(target + 1396);
    auto *t_1397 = buffer.data(target + 1397);
    auto *t_1398 = buffer.data(target + 1398);
    auto *t_1399 = buffer.data(target + 1399);
    auto *t_1400 = buffer.data(target + 1400);
    auto *t_1401 = buffer.data(target + 1401);
    auto *t_1402 = buffer.data(target + 1402);
    auto *t_1403 = buffer.data(target + 1403);
    auto *t_1404 = buffer.data(target + 1404);
    auto *t_1405 = buffer.data(target + 1405);
    auto *t_1406 = buffer.data(target + 1406);
    auto *t_1407 = buffer.data(target + 1407);
    auto *t_1408 = buffer.data(target + 1408);
    auto *t_1409 = buffer.data(target + 1409);
    auto *t_1410 = buffer.data(target + 1410);
    auto *t_1411 = buffer.data(target + 1411);
    auto *t_1412 = buffer.data(target + 1412);
    auto *t_1413 = buffer.data(target + 1413);
    auto *t_1414 = buffer.data(target + 1414);
    auto *t_1415 = buffer.data(target + 1415);
    auto *t_1416 = buffer.data(target + 1416);
    auto *t_1417 = buffer.data(target + 1417);
    auto *t_1418 = buffer.data(target + 1418);
    auto *t_1419 = buffer.data(target + 1419);
    auto *t_1420 = buffer.data(target + 1420);
    auto *t_1421 = buffer.data(target + 1421);
    auto *t_1422 = buffer.data(target + 1422);
    auto *t_1423 = buffer.data(target + 1423);
    auto *t_1424 = buffer.data(target + 1424);
    auto *t_1425 = buffer.data(target + 1425);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kk0_1295 = buffer.data(kk0 + 1295);

    const auto *kk1_1295 = buffer.data(kk1 + 1295);

    const auto *li_777 = buffer.data(li + 777);
    const auto *li_784 = buffer.data(li + 784);
    const auto *li_787 = buffer.data(li + 787);
    const auto *li_789 = buffer.data(li + 789);
    const auto *li_790 = buffer.data(li + 790);
    const auto *li_793 = buffer.data(li + 793);
    const auto *li_794 = buffer.data(li + 794);
    const auto *li_798 = buffer.data(li + 798);
    const auto *li_812 = buffer.data(li + 812);
    const auto *li_814 = buffer.data(li + 814);
    const auto *li_815 = buffer.data(li + 815);
    const auto *li_817 = buffer.data(li + 817);
    const auto *li_818 = buffer.data(li + 818);
    const auto *li_821 = buffer.data(li + 821);
    const auto *li_822 = buffer.data(li + 822);
    const auto *li_826 = buffer.data(li + 826);
    const auto *li_840 = buffer.data(li + 840);
    const auto *li_842 = buffer.data(li + 842);
    const auto *li_843 = buffer.data(li + 843);
    const auto *li_845 = buffer.data(li + 845);
    const auto *li_846 = buffer.data(li + 846);
    const auto *li_849 = buffer.data(li + 849);
    const auto *li_850 = buffer.data(li + 850);
    const auto *li_854 = buffer.data(li + 854);
    const auto *li_868 = buffer.data(li + 868);
    const auto *li_870 = buffer.data(li + 870);
    const auto *li_873 = buffer.data(li + 873);
    const auto *li_877 = buffer.data(li + 877);
    const auto *li_882 = buffer.data(li + 882);
    const auto *li_1000 = buffer.data(li + 1000);
    const auto *li_1001 = buffer.data(li + 1001);
    const auto *li_1002 = buffer.data(li + 1002);
    const auto *li_1003 = buffer.data(li + 1003);
    const auto *li_1004 = buffer.data(li + 1004);
    const auto *li_1005 = buffer.data(li + 1005);
    const auto *li_1007 = buffer.data(li + 1007);
    const auto *li_1008 = buffer.data(li + 1008);
    const auto *li_1011 = buffer.data(li + 1011);
    const auto *li_1013 = buffer.data(li + 1013);
    const auto *li_1014 = buffer.data(li + 1014);
    const auto *li_1017 = buffer.data(li + 1017);
    const auto *li_1018 = buffer.data(li + 1018);
    const auto *li_1020 = buffer.data(li + 1020);
    const auto *li_1022 = buffer.data(li + 1022);
    const auto *li_1023 = buffer.data(li + 1023);
    const auto *li_1025 = buffer.data(li + 1025);
    const auto *li_1026 = buffer.data(li + 1026);
    const auto *li_1028 = buffer.data(li + 1028);
    const auto *li_1029 = buffer.data(li + 1029);
    const auto *li_1031 = buffer.data(li + 1031);
    const auto *li_1032 = buffer.data(li + 1032);
    const auto *li_1033 = buffer.data(li + 1033);
    const auto *li_1034 = buffer.data(li + 1034);
    const auto *li_1035 = buffer.data(li + 1035);
    const auto *li_1041 = buffer.data(li + 1041);
    const auto *li_1045 = buffer.data(li + 1045);
    const auto *li_1048 = buffer.data(li + 1048);
    const auto *li_1050 = buffer.data(li + 1050);
    const auto *li_1053 = buffer.data(li + 1053);
    const auto *li_1054 = buffer.data(li + 1054);
    const auto *li_1056 = buffer.data(li + 1056);
    const auto *li_1058 = buffer.data(li + 1058);
    const auto *li_1059 = buffer.data(li + 1059);
    const auto *li_1060 = buffer.data(li + 1060);
    const auto *li_1061 = buffer.data(li + 1061);
    const auto *li_1062 = buffer.data(li + 1062);
    const auto *li_1063 = buffer.data(li + 1063);
    const auto *li_1064 = buffer.data(li + 1064);
    const auto *li_1067 = buffer.data(li + 1067);
    const auto *li_1069 = buffer.data(li + 1069);
    const auto *li_1070 = buffer.data(li + 1070);
    const auto *li_1073 = buffer.data(li + 1073);
    const auto *li_1074 = buffer.data(li + 1074);
    const auto *li_1076 = buffer.data(li + 1076);
    const auto *li_1078 = buffer.data(li + 1078);
    const auto *li_1079 = buffer.data(li + 1079);
    const auto *li_1081 = buffer.data(li + 1081);
    const auto *li_1082 = buffer.data(li + 1082);
    const auto *li_1084 = buffer.data(li + 1084);
    const auto *li_1085 = buffer.data(li + 1085);
    const auto *li_1086 = buffer.data(li + 1086);
    const auto *li_1087 = buffer.data(li + 1087);
    const auto *li_1088 = buffer.data(li + 1088);
    const auto *li_1089 = buffer.data(li + 1089);
    const auto *li_1090 = buffer.data(li + 1090);
    const auto *li_1091 = buffer.data(li + 1091);
    const auto *li_1092 = buffer.data(li + 1092);
    const auto *li_1095 = buffer.data(li + 1095);
    const auto *li_1097 = buffer.data(li + 1097);
    const auto *li_1098 = buffer.data(li + 1098);
    const auto *li_1101 = buffer.data(li + 1101);
    const auto *li_1102 = buffer.data(li + 1102);
    const auto *li_1104 = buffer.data(li + 1104);
    const auto *li_1106 = buffer.data(li + 1106);
    const auto *li_1107 = buffer.data(li + 1107);
    const auto *li_1109 = buffer.data(li + 1109);
    const auto *li_1110 = buffer.data(li + 1110);
    const auto *li_1112 = buffer.data(li + 1112);
    const auto *li_1113 = buffer.data(li + 1113);

    const auto *lk_1008 = buffer.data(lk + 1008);
    const auto *lk_1009 = buffer.data(lk + 1009);
    const auto *lk_1011 = buffer.data(lk + 1011);
    const auto *lk_1014 = buffer.data(lk + 1014);
    const auto *lk_1018 = buffer.data(lk + 1018);
    const auto *lk_1023 = buffer.data(lk + 1023);
    const auto *lk_1029 = buffer.data(lk + 1029);
    const auto *lk_1295 = buffer.data(lk + 1295);
    const auto *lk_1296 = buffer.data(lk + 1296);
    const auto *lk_1299 = buffer.data(lk + 1299);
    const auto *lk_1301 = buffer.data(lk + 1301);
    const auto *lk_1302 = buffer.data(lk + 1302);
    const auto *lk_1305 = buffer.data(lk + 1305);
    const auto *lk_1306 = buffer.data(lk + 1306);
    const auto *lk_1308 = buffer.data(lk + 1308);
    const auto *lk_1310 = buffer.data(lk + 1310);
    const auto *lk_1311 = buffer.data(lk + 1311);
    const auto *lk_1313 = buffer.data(lk + 1313);
    const auto *lk_1314 = buffer.data(lk + 1314);
    const auto *lk_1316 = buffer.data(lk + 1316);
    const auto *lk_1324 = buffer.data(lk + 1324);
    const auto *lk_1326 = buffer.data(lk + 1326);
    const auto *lk_1327 = buffer.data(lk + 1327);
    const auto *lk_1328 = buffer.data(lk + 1328);
    const auto *lk_1329 = buffer.data(lk + 1329);
    const auto *lk_1330 = buffer.data(lk + 1330);
    const auto *lk_1331 = buffer.data(lk + 1331);
    const auto *lk_1337 = buffer.data(lk + 1337);
    const auto *lk_1341 = buffer.data(lk + 1341);
    const auto *lk_1344 = buffer.data(lk + 1344);
    const auto *lk_1346 = buffer.data(lk + 1346);
    const auto *lk_1349 = buffer.data(lk + 1349);
    const auto *lk_1350 = buffer.data(lk + 1350);
    const auto *lk_1352 = buffer.data(lk + 1352);
    const auto *lk_1360 = buffer.data(lk + 1360);
    const auto *lk_1361 = buffer.data(lk + 1361);
    const auto *lk_1362 = buffer.data(lk + 1362);
    const auto *lk_1363 = buffer.data(lk + 1363);
    const auto *lk_1364 = buffer.data(lk + 1364);
    const auto *lk_1365 = buffer.data(lk + 1365);
    const auto *lk_1366 = buffer.data(lk + 1366);
    const auto *lk_1367 = buffer.data(lk + 1367);
    const auto *lk_1368 = buffer.data(lk + 1368);
    const auto *lk_1371 = buffer.data(lk + 1371);
    const auto *lk_1373 = buffer.data(lk + 1373);
    const auto *lk_1374 = buffer.data(lk + 1374);
    const auto *lk_1377 = buffer.data(lk + 1377);
    const auto *lk_1378 = buffer.data(lk + 1378);
    const auto *lk_1380 = buffer.data(lk + 1380);
    const auto *lk_1382 = buffer.data(lk + 1382);
    const auto *lk_1383 = buffer.data(lk + 1383);
    const auto *lk_1385 = buffer.data(lk + 1385);
    const auto *lk_1386 = buffer.data(lk + 1386);
    const auto *lk_1388 = buffer.data(lk + 1388);
    const auto *lk_1396 = buffer.data(lk + 1396);
    const auto *lk_1397 = buffer.data(lk + 1397);
    const auto *lk_1398 = buffer.data(lk + 1398);
    const auto *lk_1399 = buffer.data(lk + 1399);
    const auto *lk_1400 = buffer.data(lk + 1400);
    const auto *lk_1401 = buffer.data(lk + 1401);
    const auto *lk_1402 = buffer.data(lk + 1402);
    const auto *lk_1403 = buffer.data(lk + 1403);
    const auto *lk_1404 = buffer.data(lk + 1404);
    const auto *lk_1407 = buffer.data(lk + 1407);
    const auto *lk_1409 = buffer.data(lk + 1409);
    const auto *lk_1410 = buffer.data(lk + 1410);
    const auto *lk_1413 = buffer.data(lk + 1413);
    const auto *lk_1414 = buffer.data(lk + 1414);
    const auto *lk_1416 = buffer.data(lk + 1416);
    const auto *lk_1418 = buffer.data(lk + 1418);
    const auto *lk_1419 = buffer.data(lk + 1419);
    const auto *lk_1421 = buffer.data(lk + 1421);
    const auto *lk_1422 = buffer.data(lk + 1422);
    const auto *lk_1424 = buffer.data(lk + 1424);

    const auto *mh0_750 = buffer.data(mh0 + 750);
    const auto *mh0_752 = buffer.data(mh0 + 752);
    const auto *mh0_753 = buffer.data(mh0 + 753);
    const auto *mh0_754 = buffer.data(mh0 + 754);
    const auto *mh0_755 = buffer.data(mh0 + 755);

    const auto *mh1_750 = buffer.data(mh1 + 750);
    const auto *mh1_752 = buffer.data(mh1 + 752);
    const auto *mh1_753 = buffer.data(mh1 + 753);
    const auto *mh1_754 = buffer.data(mh1 + 754);
    const auto *mh1_755 = buffer.data(mh1 + 755);

    const auto *mi_994 = buffer.data(mi + 994);
    const auto *mi_1000 = buffer.data(mi + 1000);
    const auto *mi_1001 = buffer.data(mi + 1001);
    const auto *mi_1002 = buffer.data(mi + 1002);
    const auto *mi_1003 = buffer.data(mi + 1003);
    const auto *mi_1004 = buffer.data(mi + 1004);
    const auto *mi_1005 = buffer.data(mi + 1005);
    const auto *mi_1006 = buffer.data(mi + 1006);
    const auto *mi_1007 = buffer.data(mi + 1007);
    const auto *mi_1008 = buffer.data(mi + 1008);
    const auto *mi_1009 = buffer.data(mi + 1009);
    const auto *mi_1011 = buffer.data(mi + 1011);
    const auto *mi_1013 = buffer.data(mi + 1013);
    const auto *mi_1014 = buffer.data(mi + 1014);
    const auto *mi_1017 = buffer.data(mi + 1017);
    const auto *mi_1018 = buffer.data(mi + 1018);
    const auto *mi_1022 = buffer.data(mi + 1022);
    const auto *mi_1023 = buffer.data(mi + 1023);
    const auto *mi_1029 = buffer.data(mi + 1029);
    const auto *mi_1031 = buffer.data(mi + 1031);
    const auto *mi_1032 = buffer.data(mi + 1032);
    const auto *mi_1033 = buffer.data(mi + 1033);
    const auto *mi_1034 = buffer.data(mi + 1034);
    const auto *mi_1035 = buffer.data(mi + 1035);
    const auto *mi_1036 = buffer.data(mi + 1036);
    const auto *mi_1038 = buffer.data(mi + 1038);
    const auto *mi_1039 = buffer.data(mi + 1039);
    const auto *mi_1041 = buffer.data(mi + 1041);
    const auto *mi_1042 = buffer.data(mi + 1042);
    const auto *mi_1045 = buffer.data(mi + 1045);
    const auto *mi_1046 = buffer.data(mi + 1046);
    const auto *mi_1050 = buffer.data(mi + 1050);
    const auto *mi_1058 = buffer.data(mi + 1058);
    const auto *mi_1059 = buffer.data(mi + 1059);
    const auto *mi_1060 = buffer.data(mi + 1060);
    const auto *mi_1061 = buffer.data(mi + 1061);
    const auto *mi_1062 = buffer.data(mi + 1062);
    const auto *mi_1063 = buffer.data(mi + 1063);
    const auto *mi_1064 = buffer.data(mi + 1064);
    const auto *mi_1066 = buffer.data(mi + 1066);
    const auto *mi_1067 = buffer.data(mi + 1067);
    const auto *mi_1069 = buffer.data(mi + 1069);
    const auto *mi_1070 = buffer.data(mi + 1070);
    const auto *mi_1073 = buffer.data(mi + 1073);
    const auto *mi_1074 = buffer.data(mi + 1074);
    const auto *mi_1078 = buffer.data(mi + 1078);
    const auto *mi_1085 = buffer.data(mi + 1085);
    const auto *mi_1086 = buffer.data(mi + 1086);
    const auto *mi_1087 = buffer.data(mi + 1087);
    const auto *mi_1088 = buffer.data(mi + 1088);
    const auto *mi_1089 = buffer.data(mi + 1089);
    const auto *mi_1090 = buffer.data(mi + 1090);
    const auto *mi_1091 = buffer.data(mi + 1091);
    const auto *mi_1092 = buffer.data(mi + 1092);
    const auto *mi_1094 = buffer.data(mi + 1094);
    const auto *mi_1095 = buffer.data(mi + 1095);
    const auto *mi_1097 = buffer.data(mi + 1097);
    const auto *mi_1098 = buffer.data(mi + 1098);
    const auto *mi_1101 = buffer.data(mi + 1101);
    const auto *mi_1102 = buffer.data(mi + 1102);
    const auto *mi_1106 = buffer.data(mi + 1106);
    const auto *mi_1113 = buffer.data(mi + 1113);

#pragma omp simd aligned(t_1279, t_1280, t_1281, t_1282, pb_x, pb_y, li_1000, li_1001, \
                         li_1002, mh0_755, mh1_755, mi_994, mi_1000, mi_1001, \
                         mi_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1279[k] = pb_y[k] * mi_994[k];

        t_1280[k] = f_12 * li_1000[k]
                    + f_3 * mh0_755[k]
                    - f_4 * mh1_755[k]
                    + pb_x[k] * mi_1000[k];

        t_1281[k] = f_12 * li_1001[k]
                    + pb_x[k] * mi_1001[k];

        t_1282[k] = f_12 * li_1002[k]
                    + pb_x[k] * mi_1002[k];
    }

#pragma omp simd aligned(t_1283, t_1284, t_1285, t_1286, t_1287, pb_x, pb_y, li_1003, li_1004, \
                         li_1005, li_1007, mi_1000, mi_1003, mi_1004, mi_1005, \
                         mi_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1283[k] = f_12 * li_1003[k]
                    + pb_x[k] * mi_1003[k];

        t_1284[k] = f_12 * li_1004[k]
                    + pb_x[k] * mi_1004[k];

        t_1285[k] = f_12 * li_1005[k]
                    + pb_x[k] * mi_1005[k];

        t_1286[k] = pb_y[k] * mi_1000[k];

        t_1287[k] = f_12 * li_1007[k]
                    + pb_x[k] * mi_1007[k];
    }

#pragma omp simd aligned(t_1288, t_1289, t_1290, t_1291, pb_y, pb_z, li_777, mh0_750, mh0_752, \
                         mh0_753, mh1_750, mh1_752, mh1_753, mi_1001, mi_1003, \
                         mi_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1288[k] = f_1 * mh0_750[k]
                    - f_2 * mh1_750[k]
                    + pb_y[k] * mi_1001[k];

        t_1289[k] = f_17 * li_777[k]
                    + pb_z[k] * mi_1001[k];

        t_1290[k] = f_9 * mh0_752[k]
                    - f_10 * mh1_752[k]
                    + pb_y[k] * mi_1003[k];

        t_1291[k] = f_7 * mh0_753[k]
                    - f_8 * mh1_753[k]
                    + pb_y[k] * mi_1004[k];
    }

#pragma omp simd aligned(t_1292, t_1293, t_1294, t_1295, pa_x, pb_y, kk0_1295, kk1_1295, \
                         lk_1295, mh0_754, mh0_755, mh1_754, mh1_755, mi_1005, mi_1006, \
                         mi_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1292[k] = f_5 * mh0_754[k]
                    - f_6 * mh1_754[k]
                    + pb_y[k] * mi_1005[k];

        t_1293[k] = f_3 * mh0_755[k]
                    - f_4 * mh1_755[k]
                    + pb_y[k] * mi_1006[k];

        t_1294[k] = pb_y[k] * mi_1007[k];

        t_1295[k] = f_18 * kk0_1295[k]
                    - f_19 * kk1_1295[k]
                    + pa_x[k] * lk_1295[k];
    }

#pragma omp simd aligned(t_1296, t_1297, t_1298, t_1299, t_1300, pa_x, pb_y, pb_z, li_784, \
                         li_1008, li_1011, lk_1296, lk_1299, mi_1008, \
                         mi_1009 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1296[k] = f_17 * li_1008[k]
                    + pa_x[k] * lk_1296[k];

        t_1297[k] = f_16 * li_784[k]
                    + pb_y[k] * mi_1008[k];

        t_1298[k] = pb_z[k] * mi_1008[k];

        t_1299[k] = f_15 * li_1011[k]
                    + pa_x[k] * lk_1299[k];

        t_1300[k] = pb_z[k] * mi_1009[k];
    }

#pragma omp simd aligned(t_1301, t_1302, t_1303, t_1304, pa_x, pb_y, pb_z, li_789, li_1013, \
                         li_1014, lk_1301, lk_1302, mi_1011, mi_1013 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1301[k] = f_15 * li_1013[k]
                    + pa_x[k] * lk_1301[k];

        t_1302[k] = f_14 * li_1014[k]
                    + pa_x[k] * lk_1302[k];

        t_1303[k] = pb_z[k] * mi_1011[k];

        t_1304[k] = f_16 * li_789[k]
                    + pb_y[k] * mi_1013[k];
    }

#pragma omp simd aligned(t_1305, t_1306, t_1307, t_1308, pa_x, pb_z, li_1017, li_1018, \
                         li_1020, lk_1305, lk_1306, lk_1308, mi_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1305[k] = f_14 * li_1017[k]
                    + pa_x[k] * lk_1305[k];

        t_1306[k] = f_13 * li_1018[k]
                    + pa_x[k] * lk_1306[k];

        t_1307[k] = pb_z[k] * mi_1014[k];

        t_1308[k] = f_13 * li_1020[k]
                    + pa_x[k] * lk_1308[k];
    }

#pragma omp simd aligned(t_1309, t_1310, t_1311, t_1312, pa_x, pb_y, pb_z, li_793, li_1022, \
                         li_1023, lk_1310, lk_1311, mi_1017, mi_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1309[k] = f_16 * li_793[k]
                    + pb_y[k] * mi_1017[k];

        t_1310[k] = f_13 * li_1022[k]
                    + pa_x[k] * lk_1310[k];

        t_1311[k] = f_12 * li_1023[k]
                    + pa_x[k] * lk_1311[k];

        t_1312[k] = pb_z[k] * mi_1018[k];
    }

#pragma omp simd aligned(t_1313, t_1314, t_1315, t_1316, pa_x, pb_y, li_798, li_1025, li_1026, \
                         li_1028, lk_1313, lk_1314, lk_1316, mi_1022 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1313[k] = f_12 * li_1025[k]
                    + pa_x[k] * lk_1313[k];

        t_1314[k] = f_12 * li_1026[k]
                    + pa_x[k] * lk_1314[k];

        t_1315[k] = f_16 * li_798[k]
                    + pb_y[k] * mi_1022[k];

        t_1316[k] = f_12 * li_1028[k]
                    + pa_x[k] * lk_1316[k];
    }

#pragma omp simd aligned(t_1317, t_1318, t_1319, t_1320, t_1321, pb_x, pb_z, li_1029, li_1031, \
                         li_1032, li_1033, mi_1023, mi_1029, mi_1031, mi_1032, \
                         mi_1033 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1317[k] = f_11 * li_1029[k]
                    + pb_x[k] * mi_1029[k];

        t_1318[k] = pb_z[k] * mi_1023[k];

        t_1319[k] = f_11 * li_1031[k]
                    + pb_x[k] * mi_1031[k];

        t_1320[k] = f_11 * li_1032[k]
                    + pb_x[k] * mi_1032[k];

        t_1321[k] = f_11 * li_1033[k]
                    + pb_x[k] * mi_1033[k];
    }

#pragma omp simd aligned(t_1322, t_1323, t_1324, t_1325, t_1326, pa_x, pb_x, pb_z, li_1034, \
                         li_1035, lk_1324, lk_1326, mi_1029, mi_1034, \
                         mi_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1322[k] = f_11 * li_1034[k]
                    + pb_x[k] * mi_1034[k];

        t_1323[k] = f_11 * li_1035[k]
                    + pb_x[k] * mi_1035[k];

        t_1324[k] = pa_x[k] * lk_1324[k];

        t_1325[k] = pb_z[k] * mi_1029[k];

        t_1326[k] = pa_x[k] * lk_1326[k];
    }

#pragma omp simd aligned(t_1327, t_1328, t_1329, t_1330, t_1331, t_1332, t_1333, pa_x, pa_z, \
                         lk_1008, lk_1009, lk_1327, lk_1328, lk_1329, lk_1330, \
                         lk_1331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1327[k] = pa_x[k] * lk_1327[k];

        t_1328[k] = pa_x[k] * lk_1328[k];

        t_1329[k] = pa_x[k] * lk_1329[k];

        t_1330[k] = pa_x[k] * lk_1330[k];

        t_1331[k] = pa_x[k] * lk_1331[k];

        t_1332[k] = pa_z[k] * lk_1008[k];

        t_1333[k] = pa_z[k] * lk_1009[k];
    }

#pragma omp simd aligned(t_1334, t_1335, t_1336, t_1337, pa_x, pa_z, pb_y, pb_z, li_784, \
                         li_814, li_1041, lk_1011, lk_1337, mi_1036, \
                         mi_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1334[k] = f_11 * li_784[k]
                    + pb_z[k] * mi_1036[k];

        t_1335[k] = pa_z[k] * lk_1011[k];

        t_1336[k] = f_17 * li_814[k]
                    + pb_y[k] * mi_1038[k];

        t_1337[k] = f_15 * li_1041[k]
                    + pa_x[k] * lk_1337[k];
    }

#pragma omp simd aligned(t_1338, t_1339, t_1340, t_1341, pa_x, pa_z, pb_y, pb_z, li_787, \
                         li_817, li_1045, lk_1014, lk_1341, mi_1039, \
                         mi_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1338[k] = pa_z[k] * lk_1014[k];

        t_1339[k] = f_11 * li_787[k]
                    + pb_z[k] * mi_1039[k];

        t_1340[k] = f_17 * li_817[k]
                    + pb_y[k] * mi_1041[k];

        t_1341[k] = f_14 * li_1045[k]
                    + pa_x[k] * lk_1341[k];
    }

#pragma omp simd aligned(t_1342, t_1343, t_1344, t_1345, pa_x, pa_z, pb_y, pb_z, li_790, \
                         li_821, li_1048, lk_1018, lk_1344, mi_1042, \
                         mi_1045 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1342[k] = pa_z[k] * lk_1018[k];

        t_1343[k] = f_11 * li_790[k]
                    + pb_z[k] * mi_1042[k];

        t_1344[k] = f_13 * li_1048[k]
                    + pa_x[k] * lk_1344[k];

        t_1345[k] = f_17 * li_821[k]
                    + pb_y[k] * mi_1045[k];
    }

#pragma omp simd aligned(t_1346, t_1347, t_1348, t_1349, pa_x, pa_z, pb_z, li_794, li_1050, \
                         li_1053, lk_1023, lk_1346, lk_1349, mi_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1346[k] = f_13 * li_1050[k]
                    + pa_x[k] * lk_1346[k];

        t_1347[k] = pa_z[k] * lk_1023[k];

        t_1348[k] = f_11 * li_794[k]
                    + pb_z[k] * mi_1046[k];

        t_1349[k] = f_12 * li_1053[k]
                    + pa_x[k] * lk_1349[k];
    }

#pragma omp simd aligned(t_1350, t_1351, t_1352, t_1353, pa_x, pa_z, pb_y, li_826, li_1054, \
                         li_1056, lk_1029, lk_1350, lk_1352, mi_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1350[k] = f_12 * li_1054[k]
                    + pa_x[k] * lk_1350[k];

        t_1351[k] = f_17 * li_826[k]
                    + pb_y[k] * mi_1050[k];

        t_1352[k] = f_12 * li_1056[k]
                    + pa_x[k] * lk_1352[k];

        t_1353[k] = pa_z[k] * lk_1029[k];
    }

#pragma omp simd aligned(t_1354, t_1355, t_1356, t_1357, t_1358, pb_x, li_1058, li_1059, \
                         li_1060, li_1061, li_1062, mi_1058, mi_1059, mi_1060, mi_1061, \
                         mi_1062 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1354[k] = f_11 * li_1058[k]
                    + pb_x[k] * mi_1058[k];

        t_1355[k] = f_11 * li_1059[k]
                    + pb_x[k] * mi_1059[k];

        t_1356[k] = f_11 * li_1060[k]
                    + pb_x[k] * mi_1060[k];

        t_1357[k] = f_11 * li_1061[k]
                    + pb_x[k] * mi_1061[k];

        t_1358[k] = f_11 * li_1062[k]
                    + pb_x[k] * mi_1062[k];
    }

#pragma omp simd aligned(t_1359, t_1360, t_1361, t_1362, t_1363, t_1364, pa_x, pb_x, li_1063, \
                         lk_1360, lk_1361, lk_1362, lk_1363, lk_1364, \
                         mi_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1359[k] = f_11 * li_1063[k]
                    + pb_x[k] * mi_1063[k];

        t_1360[k] = pa_x[k] * lk_1360[k];

        t_1361[k] = pa_x[k] * lk_1361[k];

        t_1362[k] = pa_x[k] * lk_1362[k];

        t_1363[k] = pa_x[k] * lk_1363[k];

        t_1364[k] = pa_x[k] * lk_1364[k];
    }

#pragma omp simd aligned(t_1365, t_1366, t_1367, t_1368, t_1369, pa_x, pb_y, li_840, li_1064, \
                         lk_1365, lk_1366, lk_1367, lk_1368, mi_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1365[k] = pa_x[k] * lk_1365[k];

        t_1366[k] = pa_x[k] * lk_1366[k];

        t_1367[k] = pa_x[k] * lk_1367[k];

        t_1368[k] = f_17 * li_1064[k]
                    + pa_x[k] * lk_1368[k];

        t_1369[k] = f_24 * li_840[k]
                    + pb_y[k] * mi_1064[k];
    }

#pragma omp simd aligned(t_1370, t_1371, t_1372, t_1373, pa_x, pb_y, pb_z, li_812, li_842, \
                         li_1067, li_1069, lk_1371, lk_1373, mi_1064, \
                         mi_1066 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1370[k] = f_12 * li_812[k]
                    + pb_z[k] * mi_1064[k];

        t_1371[k] = f_15 * li_1067[k]
                    + pa_x[k] * lk_1371[k];

        t_1372[k] = f_24 * li_842[k]
                    + pb_y[k] * mi_1066[k];

        t_1373[k] = f_15 * li_1069[k]
                    + pa_x[k] * lk_1373[k];
    }

#pragma omp simd aligned(t_1374, t_1375, t_1376, t_1377, pa_x, pb_y, pb_z, li_815, li_845, \
                         li_1070, li_1073, lk_1374, lk_1377, mi_1067, \
                         mi_1069 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1374[k] = f_14 * li_1070[k]
                    + pa_x[k] * lk_1374[k];

        t_1375[k] = f_12 * li_815[k]
                    + pb_z[k] * mi_1067[k];

        t_1376[k] = f_24 * li_845[k]
                    + pb_y[k] * mi_1069[k];

        t_1377[k] = f_14 * li_1073[k]
                    + pa_x[k] * lk_1377[k];
    }

#pragma omp simd aligned(t_1378, t_1379, t_1380, t_1381, pa_x, pb_y, pb_z, li_818, li_849, \
                         li_1074, li_1076, lk_1378, lk_1380, mi_1070, \
                         mi_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1378[k] = f_13 * li_1074[k]
                    + pa_x[k] * lk_1378[k];

        t_1379[k] = f_12 * li_818[k]
                    + pb_z[k] * mi_1070[k];

        t_1380[k] = f_13 * li_1076[k]
                    + pa_x[k] * lk_1380[k];

        t_1381[k] = f_24 * li_849[k]
                    + pb_y[k] * mi_1073[k];
    }

#pragma omp simd aligned(t_1382, t_1383, t_1384, t_1385, pa_x, pb_z, li_822, li_1078, li_1079, \
                         li_1081, lk_1382, lk_1383, lk_1385, mi_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1382[k] = f_13 * li_1078[k]
                    + pa_x[k] * lk_1382[k];

        t_1383[k] = f_12 * li_1079[k]
                    + pa_x[k] * lk_1383[k];

        t_1384[k] = f_12 * li_822[k]
                    + pb_z[k] * mi_1074[k];

        t_1385[k] = f_12 * li_1081[k]
                    + pa_x[k] * lk_1385[k];
    }

#pragma omp simd aligned(t_1386, t_1387, t_1388, t_1389, pa_x, pb_x, pb_y, li_854, li_1082, \
                         li_1084, li_1085, lk_1386, lk_1388, mi_1078, \
                         mi_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1386[k] = f_12 * li_1082[k]
                    + pa_x[k] * lk_1386[k];

        t_1387[k] = f_24 * li_854[k]
                    + pb_y[k] * mi_1078[k];

        t_1388[k] = f_12 * li_1084[k]
                    + pa_x[k] * lk_1388[k];

        t_1389[k] = f_11 * li_1085[k]
                    + pb_x[k] * mi_1085[k];
    }

#pragma omp simd aligned(t_1390, t_1391, t_1392, t_1393, t_1394, pb_x, li_1086, li_1087, \
                         li_1088, li_1089, li_1090, mi_1086, mi_1087, mi_1088, mi_1089, \
                         mi_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1390[k] = f_11 * li_1086[k]
                    + pb_x[k] * mi_1086[k];

        t_1391[k] = f_11 * li_1087[k]
                    + pb_x[k] * mi_1087[k];

        t_1392[k] = f_11 * li_1088[k]
                    + pb_x[k] * mi_1088[k];

        t_1393[k] = f_11 * li_1089[k]
                    + pb_x[k] * mi_1089[k];

        t_1394[k] = f_11 * li_1090[k]
                    + pb_x[k] * mi_1090[k];
    }

#pragma omp simd aligned(t_1395, t_1396, t_1397, t_1398, t_1399, t_1400, pa_x, pb_x, li_1091, \
                         lk_1396, lk_1397, lk_1398, lk_1399, lk_1400, \
                         mi_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1395[k] = f_11 * li_1091[k]
                    + pb_x[k] * mi_1091[k];

        t_1396[k] = pa_x[k] * lk_1396[k];

        t_1397[k] = pa_x[k] * lk_1397[k];

        t_1398[k] = pa_x[k] * lk_1398[k];

        t_1399[k] = pa_x[k] * lk_1399[k];

        t_1400[k] = pa_x[k] * lk_1400[k];
    }

#pragma omp simd aligned(t_1401, t_1402, t_1403, t_1404, t_1405, pa_x, pb_y, li_868, li_1092, \
                         lk_1401, lk_1402, lk_1403, lk_1404, mi_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1401[k] = pa_x[k] * lk_1401[k];

        t_1402[k] = pa_x[k] * lk_1402[k];

        t_1403[k] = pa_x[k] * lk_1403[k];

        t_1404[k] = f_17 * li_1092[k]
                    + pa_x[k] * lk_1404[k];

        t_1405[k] = f_15 * li_868[k]
                    + pb_y[k] * mi_1092[k];
    }

#pragma omp simd aligned(t_1406, t_1407, t_1408, t_1409, pa_x, pb_y, pb_z, li_840, li_870, \
                         li_1095, li_1097, lk_1407, lk_1409, mi_1092, \
                         mi_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1406[k] = f_13 * li_840[k]
                    + pb_z[k] * mi_1092[k];

        t_1407[k] = f_15 * li_1095[k]
                    + pa_x[k] * lk_1407[k];

        t_1408[k] = f_15 * li_870[k]
                    + pb_y[k] * mi_1094[k];

        t_1409[k] = f_15 * li_1097[k]
                    + pa_x[k] * lk_1409[k];
    }

#pragma omp simd aligned(t_1410, t_1411, t_1412, t_1413, pa_x, pb_y, pb_z, li_843, li_873, \
                         li_1098, li_1101, lk_1410, lk_1413, mi_1095, \
                         mi_1097 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1410[k] = f_14 * li_1098[k]
                    + pa_x[k] * lk_1410[k];

        t_1411[k] = f_13 * li_843[k]
                    + pb_z[k] * mi_1095[k];

        t_1412[k] = f_15 * li_873[k]
                    + pb_y[k] * mi_1097[k];

        t_1413[k] = f_14 * li_1101[k]
                    + pa_x[k] * lk_1413[k];
    }

#pragma omp simd aligned(t_1414, t_1415, t_1416, t_1417, pa_x, pb_y, pb_z, li_846, li_877, \
                         li_1102, li_1104, lk_1414, lk_1416, mi_1098, \
                         mi_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1414[k] = f_13 * li_1102[k]
                    + pa_x[k] * lk_1414[k];

        t_1415[k] = f_13 * li_846[k]
                    + pb_z[k] * mi_1098[k];

        t_1416[k] = f_13 * li_1104[k]
                    + pa_x[k] * lk_1416[k];

        t_1417[k] = f_15 * li_877[k]
                    + pb_y[k] * mi_1101[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, t_1421, pa_x, pb_z, li_850, li_1106, li_1107, \
                         li_1109, lk_1418, lk_1419, lk_1421, mi_1102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = f_13 * li_1106[k]
                    + pa_x[k] * lk_1418[k];

        t_1419[k] = f_12 * li_1107[k]
                    + pa_x[k] * lk_1419[k];

        t_1420[k] = f_13 * li_850[k]
                    + pb_z[k] * mi_1102[k];

        t_1421[k] = f_12 * li_1109[k]
                    + pa_x[k] * lk_1421[k];
    }

#pragma omp simd aligned(t_1422, t_1423, t_1424, t_1425, pa_x, pb_x, pb_y, li_882, li_1110, \
                         li_1112, li_1113, lk_1422, lk_1424, mi_1106, \
                         mi_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1422[k] = f_12 * li_1110[k]
                    + pa_x[k] * lk_1422[k];

        t_1423[k] = f_15 * li_882[k]
                    + pb_y[k] * mi_1106[k];

        t_1424[k] = f_12 * li_1112[k]
                    + pa_x[k] * lk_1424[k];

        t_1425[k] = f_11 * li_1113[k]
                    + pb_x[k] * mi_1113[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece11(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t li,
                                             const size_t lk, const size_t mi,
                                             const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_17 = 3.5 / p;
    const auto f_24 = 3.0 / p;

    auto *t_1426 = buffer.data(target + 1426);
    auto *t_1427 = buffer.data(target + 1427);
    auto *t_1428 = buffer.data(target + 1428);
    auto *t_1429 = buffer.data(target + 1429);
    auto *t_1430 = buffer.data(target + 1430);
    auto *t_1431 = buffer.data(target + 1431);
    auto *t_1432 = buffer.data(target + 1432);
    auto *t_1433 = buffer.data(target + 1433);
    auto *t_1434 = buffer.data(target + 1434);
    auto *t_1435 = buffer.data(target + 1435);
    auto *t_1436 = buffer.data(target + 1436);
    auto *t_1437 = buffer.data(target + 1437);
    auto *t_1438 = buffer.data(target + 1438);
    auto *t_1439 = buffer.data(target + 1439);
    auto *t_1440 = buffer.data(target + 1440);
    auto *t_1441 = buffer.data(target + 1441);
    auto *t_1442 = buffer.data(target + 1442);
    auto *t_1443 = buffer.data(target + 1443);
    auto *t_1444 = buffer.data(target + 1444);
    auto *t_1445 = buffer.data(target + 1445);
    auto *t_1446 = buffer.data(target + 1446);
    auto *t_1447 = buffer.data(target + 1447);
    auto *t_1448 = buffer.data(target + 1448);
    auto *t_1449 = buffer.data(target + 1449);
    auto *t_1450 = buffer.data(target + 1450);
    auto *t_1451 = buffer.data(target + 1451);
    auto *t_1452 = buffer.data(target + 1452);
    auto *t_1453 = buffer.data(target + 1453);
    auto *t_1454 = buffer.data(target + 1454);
    auto *t_1455 = buffer.data(target + 1455);
    auto *t_1456 = buffer.data(target + 1456);
    auto *t_1457 = buffer.data(target + 1457);
    auto *t_1458 = buffer.data(target + 1458);
    auto *t_1459 = buffer.data(target + 1459);
    auto *t_1460 = buffer.data(target + 1460);
    auto *t_1461 = buffer.data(target + 1461);
    auto *t_1462 = buffer.data(target + 1462);
    auto *t_1463 = buffer.data(target + 1463);
    auto *t_1464 = buffer.data(target + 1464);
    auto *t_1465 = buffer.data(target + 1465);
    auto *t_1466 = buffer.data(target + 1466);
    auto *t_1467 = buffer.data(target + 1467);
    auto *t_1468 = buffer.data(target + 1468);
    auto *t_1469 = buffer.data(target + 1469);
    auto *t_1470 = buffer.data(target + 1470);
    auto *t_1471 = buffer.data(target + 1471);
    auto *t_1472 = buffer.data(target + 1472);
    auto *t_1473 = buffer.data(target + 1473);
    auto *t_1474 = buffer.data(target + 1474);
    auto *t_1475 = buffer.data(target + 1475);
    auto *t_1476 = buffer.data(target + 1476);
    auto *t_1477 = buffer.data(target + 1477);
    auto *t_1478 = buffer.data(target + 1478);
    auto *t_1479 = buffer.data(target + 1479);
    auto *t_1480 = buffer.data(target + 1480);
    auto *t_1481 = buffer.data(target + 1481);
    auto *t_1482 = buffer.data(target + 1482);
    auto *t_1483 = buffer.data(target + 1483);
    auto *t_1484 = buffer.data(target + 1484);
    auto *t_1485 = buffer.data(target + 1485);
    auto *t_1486 = buffer.data(target + 1486);
    auto *t_1487 = buffer.data(target + 1487);
    auto *t_1488 = buffer.data(target + 1488);
    auto *t_1489 = buffer.data(target + 1489);
    auto *t_1490 = buffer.data(target + 1490);
    auto *t_1491 = buffer.data(target + 1491);
    auto *t_1492 = buffer.data(target + 1492);
    auto *t_1493 = buffer.data(target + 1493);
    auto *t_1494 = buffer.data(target + 1494);
    auto *t_1495 = buffer.data(target + 1495);
    auto *t_1496 = buffer.data(target + 1496);
    auto *t_1497 = buffer.data(target + 1497);
    auto *t_1498 = buffer.data(target + 1498);
    auto *t_1499 = buffer.data(target + 1499);
    auto *t_1500 = buffer.data(target + 1500);
    auto *t_1501 = buffer.data(target + 1501);
    auto *t_1502 = buffer.data(target + 1502);
    auto *t_1503 = buffer.data(target + 1503);
    auto *t_1504 = buffer.data(target + 1504);
    auto *t_1505 = buffer.data(target + 1505);
    auto *t_1506 = buffer.data(target + 1506);
    auto *t_1507 = buffer.data(target + 1507);
    auto *t_1508 = buffer.data(target + 1508);
    auto *t_1509 = buffer.data(target + 1509);
    auto *t_1510 = buffer.data(target + 1510);
    auto *t_1511 = buffer.data(target + 1511);
    auto *t_1512 = buffer.data(target + 1512);
    auto *t_1513 = buffer.data(target + 1513);
    auto *t_1514 = buffer.data(target + 1514);
    auto *t_1515 = buffer.data(target + 1515);
    auto *t_1516 = buffer.data(target + 1516);
    auto *t_1517 = buffer.data(target + 1517);
    auto *t_1518 = buffer.data(target + 1518);
    auto *t_1519 = buffer.data(target + 1519);
    auto *t_1520 = buffer.data(target + 1520);
    auto *t_1521 = buffer.data(target + 1521);
    auto *t_1522 = buffer.data(target + 1522);
    auto *t_1523 = buffer.data(target + 1523);
    auto *t_1524 = buffer.data(target + 1524);
    auto *t_1525 = buffer.data(target + 1525);
    auto *t_1526 = buffer.data(target + 1526);
    auto *t_1527 = buffer.data(target + 1527);
    auto *t_1528 = buffer.data(target + 1528);
    auto *t_1529 = buffer.data(target + 1529);
    auto *t_1530 = buffer.data(target + 1530);
    auto *t_1531 = buffer.data(target + 1531);
    auto *t_1532 = buffer.data(target + 1532);
    auto *t_1533 = buffer.data(target + 1533);
    auto *t_1534 = buffer.data(target + 1534);
    auto *t_1535 = buffer.data(target + 1535);
    auto *t_1536 = buffer.data(target + 1536);
    auto *t_1537 = buffer.data(target + 1537);
    auto *t_1538 = buffer.data(target + 1538);
    auto *t_1539 = buffer.data(target + 1539);
    auto *t_1540 = buffer.data(target + 1540);
    auto *t_1541 = buffer.data(target + 1541);
    auto *t_1542 = buffer.data(target + 1542);
    auto *t_1543 = buffer.data(target + 1543);
    auto *t_1544 = buffer.data(target + 1544);
    auto *t_1545 = buffer.data(target + 1545);
    auto *t_1546 = buffer.data(target + 1546);
    auto *t_1547 = buffer.data(target + 1547);
    auto *t_1548 = buffer.data(target + 1548);
    auto *t_1549 = buffer.data(target + 1549);
    auto *t_1550 = buffer.data(target + 1550);
    auto *t_1551 = buffer.data(target + 1551);
    auto *t_1552 = buffer.data(target + 1552);
    auto *t_1553 = buffer.data(target + 1553);
    auto *t_1554 = buffer.data(target + 1554);
    auto *t_1555 = buffer.data(target + 1555);
    auto *t_1556 = buffer.data(target + 1556);
    auto *t_1557 = buffer.data(target + 1557);
    auto *t_1558 = buffer.data(target + 1558);
    auto *t_1559 = buffer.data(target + 1559);
    auto *t_1560 = buffer.data(target + 1560);
    auto *t_1561 = buffer.data(target + 1561);
    auto *t_1562 = buffer.data(target + 1562);
    auto *t_1563 = buffer.data(target + 1563);
    auto *t_1564 = buffer.data(target + 1564);
    auto *t_1565 = buffer.data(target + 1565);
    auto *t_1566 = buffer.data(target + 1566);
    auto *t_1567 = buffer.data(target + 1567);
    auto *t_1568 = buffer.data(target + 1568);
    auto *t_1569 = buffer.data(target + 1569);
    auto *t_1570 = buffer.data(target + 1570);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *li_868 = buffer.data(li + 868);
    const auto *li_871 = buffer.data(li + 871);
    const auto *li_874 = buffer.data(li + 874);
    const auto *li_878 = buffer.data(li + 878);
    const auto *li_896 = buffer.data(li + 896);
    const auto *li_898 = buffer.data(li + 898);
    const auto *li_899 = buffer.data(li + 899);
    const auto *li_901 = buffer.data(li + 901);
    const auto *li_902 = buffer.data(li + 902);
    const auto *li_905 = buffer.data(li + 905);
    const auto *li_906 = buffer.data(li + 906);
    const auto *li_910 = buffer.data(li + 910);
    const auto *li_924 = buffer.data(li + 924);
    const auto *li_926 = buffer.data(li + 926);
    const auto *li_927 = buffer.data(li + 927);
    const auto *li_929 = buffer.data(li + 929);
    const auto *li_930 = buffer.data(li + 930);
    const auto *li_933 = buffer.data(li + 933);
    const auto *li_934 = buffer.data(li + 934);
    const auto *li_938 = buffer.data(li + 938);
    const auto *li_952 = buffer.data(li + 952);
    const auto *li_954 = buffer.data(li + 954);
    const auto *li_955 = buffer.data(li + 955);
    const auto *li_957 = buffer.data(li + 957);
    const auto *li_958 = buffer.data(li + 958);
    const auto *li_961 = buffer.data(li + 961);
    const auto *li_962 = buffer.data(li + 962);
    const auto *li_966 = buffer.data(li + 966);
    const auto *li_980 = buffer.data(li + 980);
    const auto *li_982 = buffer.data(li + 982);
    const auto *li_985 = buffer.data(li + 985);
    const auto *li_989 = buffer.data(li + 989);
    const auto *li_994 = buffer.data(li + 994);
    const auto *li_1114 = buffer.data(li + 1114);
    const auto *li_1115 = buffer.data(li + 1115);
    const auto *li_1116 = buffer.data(li + 1116);
    const auto *li_1117 = buffer.data(li + 1117);
    const auto *li_1118 = buffer.data(li + 1118);
    const auto *li_1119 = buffer.data(li + 1119);
    const auto *li_1120 = buffer.data(li + 1120);
    const auto *li_1123 = buffer.data(li + 1123);
    const auto *li_1125 = buffer.data(li + 1125);
    const auto *li_1126 = buffer.data(li + 1126);
    const auto *li_1129 = buffer.data(li + 1129);
    const auto *li_1130 = buffer.data(li + 1130);
    const auto *li_1132 = buffer.data(li + 1132);
    const auto *li_1134 = buffer.data(li + 1134);
    const auto *li_1135 = buffer.data(li + 1135);
    const auto *li_1137 = buffer.data(li + 1137);
    const auto *li_1138 = buffer.data(li + 1138);
    const auto *li_1140 = buffer.data(li + 1140);
    const auto *li_1141 = buffer.data(li + 1141);
    const auto *li_1142 = buffer.data(li + 1142);
    const auto *li_1143 = buffer.data(li + 1143);
    const auto *li_1144 = buffer.data(li + 1144);
    const auto *li_1145 = buffer.data(li + 1145);
    const auto *li_1146 = buffer.data(li + 1146);
    const auto *li_1147 = buffer.data(li + 1147);
    const auto *li_1148 = buffer.data(li + 1148);
    const auto *li_1151 = buffer.data(li + 1151);
    const auto *li_1153 = buffer.data(li + 1153);
    const auto *li_1154 = buffer.data(li + 1154);
    const auto *li_1157 = buffer.data(li + 1157);
    const auto *li_1158 = buffer.data(li + 1158);
    const auto *li_1160 = buffer.data(li + 1160);
    const auto *li_1162 = buffer.data(li + 1162);
    const auto *li_1163 = buffer.data(li + 1163);
    const auto *li_1165 = buffer.data(li + 1165);
    const auto *li_1166 = buffer.data(li + 1166);
    const auto *li_1168 = buffer.data(li + 1168);
    const auto *li_1169 = buffer.data(li + 1169);
    const auto *li_1170 = buffer.data(li + 1170);
    const auto *li_1171 = buffer.data(li + 1171);
    const auto *li_1172 = buffer.data(li + 1172);
    const auto *li_1173 = buffer.data(li + 1173);
    const auto *li_1174 = buffer.data(li + 1174);
    const auto *li_1175 = buffer.data(li + 1175);
    const auto *li_1176 = buffer.data(li + 1176);
    const auto *li_1179 = buffer.data(li + 1179);
    const auto *li_1181 = buffer.data(li + 1181);
    const auto *li_1182 = buffer.data(li + 1182);
    const auto *li_1185 = buffer.data(li + 1185);
    const auto *li_1186 = buffer.data(li + 1186);
    const auto *li_1188 = buffer.data(li + 1188);
    const auto *li_1190 = buffer.data(li + 1190);
    const auto *li_1191 = buffer.data(li + 1191);
    const auto *li_1193 = buffer.data(li + 1193);
    const auto *li_1194 = buffer.data(li + 1194);
    const auto *li_1196 = buffer.data(li + 1196);
    const auto *li_1197 = buffer.data(li + 1197);
    const auto *li_1198 = buffer.data(li + 1198);
    const auto *li_1199 = buffer.data(li + 1199);
    const auto *li_1200 = buffer.data(li + 1200);
    const auto *li_1201 = buffer.data(li + 1201);
    const auto *li_1202 = buffer.data(li + 1202);
    const auto *li_1203 = buffer.data(li + 1203);
    const auto *li_1207 = buffer.data(li + 1207);
    const auto *li_1210 = buffer.data(li + 1210);
    const auto *li_1214 = buffer.data(li + 1214);
    const auto *li_1216 = buffer.data(li + 1216);
    const auto *li_1219 = buffer.data(li + 1219);
    const auto *li_1221 = buffer.data(li + 1221);
    const auto *li_1222 = buffer.data(li + 1222);
    const auto *li_1225 = buffer.data(li + 1225);
    const auto *li_1226 = buffer.data(li + 1226);

    const auto *lk_1260 = buffer.data(lk + 1260);
    const auto *lk_1262 = buffer.data(lk + 1262);
    const auto *lk_1265 = buffer.data(lk + 1265);
    const auto *lk_1269 = buffer.data(lk + 1269);
    const auto *lk_1274 = buffer.data(lk + 1274);
    const auto *lk_1280 = buffer.data(lk + 1280);
    const auto *lk_1432 = buffer.data(lk + 1432);
    const auto *lk_1433 = buffer.data(lk + 1433);
    const auto *lk_1434 = buffer.data(lk + 1434);
    const auto *lk_1435 = buffer.data(lk + 1435);
    const auto *lk_1436 = buffer.data(lk + 1436);
    const auto *lk_1437 = buffer.data(lk + 1437);
    const auto *lk_1438 = buffer.data(lk + 1438);
    const auto *lk_1439 = buffer.data(lk + 1439);
    const auto *lk_1440 = buffer.data(lk + 1440);
    const auto *lk_1443 = buffer.data(lk + 1443);
    const auto *lk_1445 = buffer.data(lk + 1445);
    const auto *lk_1446 = buffer.data(lk + 1446);
    const auto *lk_1449 = buffer.data(lk + 1449);
    const auto *lk_1450 = buffer.data(lk + 1450);
    const auto *lk_1452 = buffer.data(lk + 1452);
    const auto *lk_1454 = buffer.data(lk + 1454);
    const auto *lk_1455 = buffer.data(lk + 1455);
    const auto *lk_1457 = buffer.data(lk + 1457);
    const auto *lk_1458 = buffer.data(lk + 1458);
    const auto *lk_1460 = buffer.data(lk + 1460);
    const auto *lk_1468 = buffer.data(lk + 1468);
    const auto *lk_1469 = buffer.data(lk + 1469);
    const auto *lk_1470 = buffer.data(lk + 1470);
    const auto *lk_1471 = buffer.data(lk + 1471);
    const auto *lk_1472 = buffer.data(lk + 1472);
    const auto *lk_1473 = buffer.data(lk + 1473);
    const auto *lk_1474 = buffer.data(lk + 1474);
    const auto *lk_1475 = buffer.data(lk + 1475);
    const auto *lk_1476 = buffer.data(lk + 1476);
    const auto *lk_1479 = buffer.data(lk + 1479);
    const auto *lk_1481 = buffer.data(lk + 1481);
    const auto *lk_1482 = buffer.data(lk + 1482);
    const auto *lk_1485 = buffer.data(lk + 1485);
    const auto *lk_1486 = buffer.data(lk + 1486);
    const auto *lk_1488 = buffer.data(lk + 1488);
    const auto *lk_1490 = buffer.data(lk + 1490);
    const auto *lk_1491 = buffer.data(lk + 1491);
    const auto *lk_1493 = buffer.data(lk + 1493);
    const auto *lk_1494 = buffer.data(lk + 1494);
    const auto *lk_1496 = buffer.data(lk + 1496);
    const auto *lk_1504 = buffer.data(lk + 1504);
    const auto *lk_1505 = buffer.data(lk + 1505);
    const auto *lk_1506 = buffer.data(lk + 1506);
    const auto *lk_1507 = buffer.data(lk + 1507);
    const auto *lk_1508 = buffer.data(lk + 1508);
    const auto *lk_1509 = buffer.data(lk + 1509);
    const auto *lk_1510 = buffer.data(lk + 1510);
    const auto *lk_1511 = buffer.data(lk + 1511);
    const auto *lk_1512 = buffer.data(lk + 1512);
    const auto *lk_1515 = buffer.data(lk + 1515);
    const auto *lk_1517 = buffer.data(lk + 1517);
    const auto *lk_1518 = buffer.data(lk + 1518);
    const auto *lk_1521 = buffer.data(lk + 1521);
    const auto *lk_1522 = buffer.data(lk + 1522);
    const auto *lk_1524 = buffer.data(lk + 1524);
    const auto *lk_1526 = buffer.data(lk + 1526);
    const auto *lk_1527 = buffer.data(lk + 1527);
    const auto *lk_1529 = buffer.data(lk + 1529);
    const auto *lk_1530 = buffer.data(lk + 1530);
    const auto *lk_1532 = buffer.data(lk + 1532);
    const auto *lk_1540 = buffer.data(lk + 1540);
    const auto *lk_1541 = buffer.data(lk + 1541);
    const auto *lk_1542 = buffer.data(lk + 1542);
    const auto *lk_1543 = buffer.data(lk + 1543);
    const auto *lk_1544 = buffer.data(lk + 1544);
    const auto *lk_1545 = buffer.data(lk + 1545);
    const auto *lk_1546 = buffer.data(lk + 1546);
    const auto *lk_1547 = buffer.data(lk + 1547);
    const auto *lk_1551 = buffer.data(lk + 1551);
    const auto *lk_1554 = buffer.data(lk + 1554);
    const auto *lk_1558 = buffer.data(lk + 1558);
    const auto *lk_1560 = buffer.data(lk + 1560);
    const auto *lk_1563 = buffer.data(lk + 1563);
    const auto *lk_1565 = buffer.data(lk + 1565);
    const auto *lk_1566 = buffer.data(lk + 1566);

    const auto *mi_1114 = buffer.data(mi + 1114);
    const auto *mi_1115 = buffer.data(mi + 1115);
    const auto *mi_1116 = buffer.data(mi + 1116);
    const auto *mi_1117 = buffer.data(mi + 1117);
    const auto *mi_1118 = buffer.data(mi + 1118);
    const auto *mi_1119 = buffer.data(mi + 1119);
    const auto *mi_1120 = buffer.data(mi + 1120);
    const auto *mi_1122 = buffer.data(mi + 1122);
    const auto *mi_1123 = buffer.data(mi + 1123);
    const auto *mi_1125 = buffer.data(mi + 1125);
    const auto *mi_1126 = buffer.data(mi + 1126);
    const auto *mi_1129 = buffer.data(mi + 1129);
    const auto *mi_1130 = buffer.data(mi + 1130);
    const auto *mi_1134 = buffer.data(mi + 1134);
    const auto *mi_1141 = buffer.data(mi + 1141);
    const auto *mi_1142 = buffer.data(mi + 1142);
    const auto *mi_1143 = buffer.data(mi + 1143);
    const auto *mi_1144 = buffer.data(mi + 1144);
    const auto *mi_1145 = buffer.data(mi + 1145);
    const auto *mi_1146 = buffer.data(mi + 1146);
    const auto *mi_1147 = buffer.data(mi + 1147);
    const auto *mi_1148 = buffer.data(mi + 1148);
    const auto *mi_1150 = buffer.data(mi + 1150);
    const auto *mi_1151 = buffer.data(mi + 1151);
    const auto *mi_1153 = buffer.data(mi + 1153);
    const auto *mi_1154 = buffer.data(mi + 1154);
    const auto *mi_1157 = buffer.data(mi + 1157);
    const auto *mi_1158 = buffer.data(mi + 1158);
    const auto *mi_1162 = buffer.data(mi + 1162);
    const auto *mi_1169 = buffer.data(mi + 1169);
    const auto *mi_1170 = buffer.data(mi + 1170);
    const auto *mi_1171 = buffer.data(mi + 1171);
    const auto *mi_1172 = buffer.data(mi + 1172);
    const auto *mi_1173 = buffer.data(mi + 1173);
    const auto *mi_1174 = buffer.data(mi + 1174);
    const auto *mi_1175 = buffer.data(mi + 1175);
    const auto *mi_1176 = buffer.data(mi + 1176);
    const auto *mi_1178 = buffer.data(mi + 1178);
    const auto *mi_1179 = buffer.data(mi + 1179);
    const auto *mi_1181 = buffer.data(mi + 1181);
    const auto *mi_1182 = buffer.data(mi + 1182);
    const auto *mi_1185 = buffer.data(mi + 1185);
    const auto *mi_1186 = buffer.data(mi + 1186);
    const auto *mi_1190 = buffer.data(mi + 1190);
    const auto *mi_1197 = buffer.data(mi + 1197);
    const auto *mi_1198 = buffer.data(mi + 1198);
    const auto *mi_1199 = buffer.data(mi + 1199);
    const auto *mi_1200 = buffer.data(mi + 1200);
    const auto *mi_1201 = buffer.data(mi + 1201);
    const auto *mi_1202 = buffer.data(mi + 1202);
    const auto *mi_1203 = buffer.data(mi + 1203);
    const auto *mi_1204 = buffer.data(mi + 1204);
    const auto *mi_1206 = buffer.data(mi + 1206);
    const auto *mi_1207 = buffer.data(mi + 1207);
    const auto *mi_1209 = buffer.data(mi + 1209);
    const auto *mi_1210 = buffer.data(mi + 1210);
    const auto *mi_1213 = buffer.data(mi + 1213);
    const auto *mi_1214 = buffer.data(mi + 1214);
    const auto *mi_1218 = buffer.data(mi + 1218);
    const auto *mi_1225 = buffer.data(mi + 1225);
    const auto *mi_1226 = buffer.data(mi + 1226);

#pragma omp simd aligned(t_1426, t_1427, t_1428, t_1429, t_1430, pb_x, li_1114, li_1115, \
                         li_1116, li_1117, li_1118, mi_1114, mi_1115, mi_1116, mi_1117, \
                         mi_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1426[k] = f_11 * li_1114[k]
                    + pb_x[k] * mi_1114[k];

        t_1427[k] = f_11 * li_1115[k]
                    + pb_x[k] * mi_1115[k];

        t_1428[k] = f_11 * li_1116[k]
                    + pb_x[k] * mi_1116[k];

        t_1429[k] = f_11 * li_1117[k]
                    + pb_x[k] * mi_1117[k];

        t_1430[k] = f_11 * li_1118[k]
                    + pb_x[k] * mi_1118[k];
    }

#pragma omp simd aligned(t_1431, t_1432, t_1433, t_1434, t_1435, t_1436, pa_x, pb_x, li_1119, \
                         lk_1432, lk_1433, lk_1434, lk_1435, lk_1436, \
                         mi_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1431[k] = f_11 * li_1119[k]
                    + pb_x[k] * mi_1119[k];

        t_1432[k] = pa_x[k] * lk_1432[k];

        t_1433[k] = pa_x[k] * lk_1433[k];

        t_1434[k] = pa_x[k] * lk_1434[k];

        t_1435[k] = pa_x[k] * lk_1435[k];

        t_1436[k] = pa_x[k] * lk_1436[k];
    }

#pragma omp simd aligned(t_1437, t_1438, t_1439, t_1440, t_1441, pa_x, pb_y, li_896, li_1120, \
                         lk_1437, lk_1438, lk_1439, lk_1440, mi_1120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1437[k] = pa_x[k] * lk_1437[k];

        t_1438[k] = pa_x[k] * lk_1438[k];

        t_1439[k] = pa_x[k] * lk_1439[k];

        t_1440[k] = f_17 * li_1120[k]
                    + pa_x[k] * lk_1440[k];

        t_1441[k] = f_14 * li_896[k]
                    + pb_y[k] * mi_1120[k];
    }

#pragma omp simd aligned(t_1442, t_1443, t_1444, t_1445, pa_x, pb_y, pb_z, li_868, li_898, \
                         li_1123, li_1125, lk_1443, lk_1445, mi_1120, \
                         mi_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1442[k] = f_14 * li_868[k]
                    + pb_z[k] * mi_1120[k];

        t_1443[k] = f_15 * li_1123[k]
                    + pa_x[k] * lk_1443[k];

        t_1444[k] = f_14 * li_898[k]
                    + pb_y[k] * mi_1122[k];

        t_1445[k] = f_15 * li_1125[k]
                    + pa_x[k] * lk_1445[k];
    }

#pragma omp simd aligned(t_1446, t_1447, t_1448, t_1449, pa_x, pb_y, pb_z, li_871, li_901, \
                         li_1126, li_1129, lk_1446, lk_1449, mi_1123, \
                         mi_1125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1446[k] = f_14 * li_1126[k]
                    + pa_x[k] * lk_1446[k];

        t_1447[k] = f_14 * li_871[k]
                    + pb_z[k] * mi_1123[k];

        t_1448[k] = f_14 * li_901[k]
                    + pb_y[k] * mi_1125[k];

        t_1449[k] = f_14 * li_1129[k]
                    + pa_x[k] * lk_1449[k];
    }

#pragma omp simd aligned(t_1450, t_1451, t_1452, t_1453, pa_x, pb_y, pb_z, li_874, li_905, \
                         li_1130, li_1132, lk_1450, lk_1452, mi_1126, \
                         mi_1129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1450[k] = f_13 * li_1130[k]
                    + pa_x[k] * lk_1450[k];

        t_1451[k] = f_14 * li_874[k]
                    + pb_z[k] * mi_1126[k];

        t_1452[k] = f_13 * li_1132[k]
                    + pa_x[k] * lk_1452[k];

        t_1453[k] = f_14 * li_905[k]
                    + pb_y[k] * mi_1129[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, t_1457, pa_x, pb_z, li_878, li_1134, li_1135, \
                         li_1137, lk_1454, lk_1455, lk_1457, mi_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_13 * li_1134[k]
                    + pa_x[k] * lk_1454[k];

        t_1455[k] = f_12 * li_1135[k]
                    + pa_x[k] * lk_1455[k];

        t_1456[k] = f_14 * li_878[k]
                    + pb_z[k] * mi_1130[k];

        t_1457[k] = f_12 * li_1137[k]
                    + pa_x[k] * lk_1457[k];
    }

#pragma omp simd aligned(t_1458, t_1459, t_1460, t_1461, pa_x, pb_x, pb_y, li_910, li_1138, \
                         li_1140, li_1141, lk_1458, lk_1460, mi_1134, \
                         mi_1141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1458[k] = f_12 * li_1138[k]
                    + pa_x[k] * lk_1458[k];

        t_1459[k] = f_14 * li_910[k]
                    + pb_y[k] * mi_1134[k];

        t_1460[k] = f_12 * li_1140[k]
                    + pa_x[k] * lk_1460[k];

        t_1461[k] = f_11 * li_1141[k]
                    + pb_x[k] * mi_1141[k];
    }

#pragma omp simd aligned(t_1462, t_1463, t_1464, t_1465, t_1466, pb_x, li_1142, li_1143, \
                         li_1144, li_1145, li_1146, mi_1142, mi_1143, mi_1144, mi_1145, \
                         mi_1146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1462[k] = f_11 * li_1142[k]
                    + pb_x[k] * mi_1142[k];

        t_1463[k] = f_11 * li_1143[k]
                    + pb_x[k] * mi_1143[k];

        t_1464[k] = f_11 * li_1144[k]
                    + pb_x[k] * mi_1144[k];

        t_1465[k] = f_11 * li_1145[k]
                    + pb_x[k] * mi_1145[k];

        t_1466[k] = f_11 * li_1146[k]
                    + pb_x[k] * mi_1146[k];
    }

#pragma omp simd aligned(t_1467, t_1468, t_1469, t_1470, t_1471, t_1472, pa_x, pb_x, li_1147, \
                         lk_1468, lk_1469, lk_1470, lk_1471, lk_1472, \
                         mi_1147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1467[k] = f_11 * li_1147[k]
                    + pb_x[k] * mi_1147[k];

        t_1468[k] = pa_x[k] * lk_1468[k];

        t_1469[k] = pa_x[k] * lk_1469[k];

        t_1470[k] = pa_x[k] * lk_1470[k];

        t_1471[k] = pa_x[k] * lk_1471[k];

        t_1472[k] = pa_x[k] * lk_1472[k];
    }

#pragma omp simd aligned(t_1473, t_1474, t_1475, t_1476, t_1477, pa_x, pb_y, li_924, li_1148, \
                         lk_1473, lk_1474, lk_1475, lk_1476, mi_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1473[k] = pa_x[k] * lk_1473[k];

        t_1474[k] = pa_x[k] * lk_1474[k];

        t_1475[k] = pa_x[k] * lk_1475[k];

        t_1476[k] = f_17 * li_1148[k]
                    + pa_x[k] * lk_1476[k];

        t_1477[k] = f_13 * li_924[k]
                    + pb_y[k] * mi_1148[k];
    }

#pragma omp simd aligned(t_1478, t_1479, t_1480, t_1481, pa_x, pb_y, pb_z, li_896, li_926, \
                         li_1151, li_1153, lk_1479, lk_1481, mi_1148, \
                         mi_1150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1478[k] = f_15 * li_896[k]
                    + pb_z[k] * mi_1148[k];

        t_1479[k] = f_15 * li_1151[k]
                    + pa_x[k] * lk_1479[k];

        t_1480[k] = f_13 * li_926[k]
                    + pb_y[k] * mi_1150[k];

        t_1481[k] = f_15 * li_1153[k]
                    + pa_x[k] * lk_1481[k];
    }

#pragma omp simd aligned(t_1482, t_1483, t_1484, t_1485, pa_x, pb_y, pb_z, li_899, li_929, \
                         li_1154, li_1157, lk_1482, lk_1485, mi_1151, \
                         mi_1153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1482[k] = f_14 * li_1154[k]
                    + pa_x[k] * lk_1482[k];

        t_1483[k] = f_15 * li_899[k]
                    + pb_z[k] * mi_1151[k];

        t_1484[k] = f_13 * li_929[k]
                    + pb_y[k] * mi_1153[k];

        t_1485[k] = f_14 * li_1157[k]
                    + pa_x[k] * lk_1485[k];
    }

#pragma omp simd aligned(t_1486, t_1487, t_1488, t_1489, pa_x, pb_y, pb_z, li_902, li_933, \
                         li_1158, li_1160, lk_1486, lk_1488, mi_1154, \
                         mi_1157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1486[k] = f_13 * li_1158[k]
                    + pa_x[k] * lk_1486[k];

        t_1487[k] = f_15 * li_902[k]
                    + pb_z[k] * mi_1154[k];

        t_1488[k] = f_13 * li_1160[k]
                    + pa_x[k] * lk_1488[k];

        t_1489[k] = f_13 * li_933[k]
                    + pb_y[k] * mi_1157[k];
    }

#pragma omp simd aligned(t_1490, t_1491, t_1492, t_1493, pa_x, pb_z, li_906, li_1162, li_1163, \
                         li_1165, lk_1490, lk_1491, lk_1493, mi_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1490[k] = f_13 * li_1162[k]
                    + pa_x[k] * lk_1490[k];

        t_1491[k] = f_12 * li_1163[k]
                    + pa_x[k] * lk_1491[k];

        t_1492[k] = f_15 * li_906[k]
                    + pb_z[k] * mi_1158[k];

        t_1493[k] = f_12 * li_1165[k]
                    + pa_x[k] * lk_1493[k];
    }

#pragma omp simd aligned(t_1494, t_1495, t_1496, t_1497, pa_x, pb_x, pb_y, li_938, li_1166, \
                         li_1168, li_1169, lk_1494, lk_1496, mi_1162, \
                         mi_1169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1494[k] = f_12 * li_1166[k]
                    + pa_x[k] * lk_1494[k];

        t_1495[k] = f_13 * li_938[k]
                    + pb_y[k] * mi_1162[k];

        t_1496[k] = f_12 * li_1168[k]
                    + pa_x[k] * lk_1496[k];

        t_1497[k] = f_11 * li_1169[k]
                    + pb_x[k] * mi_1169[k];
    }

#pragma omp simd aligned(t_1498, t_1499, t_1500, t_1501, t_1502, pb_x, li_1170, li_1171, \
                         li_1172, li_1173, li_1174, mi_1170, mi_1171, mi_1172, mi_1173, \
                         mi_1174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1498[k] = f_11 * li_1170[k]
                    + pb_x[k] * mi_1170[k];

        t_1499[k] = f_11 * li_1171[k]
                    + pb_x[k] * mi_1171[k];

        t_1500[k] = f_11 * li_1172[k]
                    + pb_x[k] * mi_1172[k];

        t_1501[k] = f_11 * li_1173[k]
                    + pb_x[k] * mi_1173[k];

        t_1502[k] = f_11 * li_1174[k]
                    + pb_x[k] * mi_1174[k];
    }

#pragma omp simd aligned(t_1503, t_1504, t_1505, t_1506, t_1507, t_1508, pa_x, pb_x, li_1175, \
                         lk_1504, lk_1505, lk_1506, lk_1507, lk_1508, \
                         mi_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1503[k] = f_11 * li_1175[k]
                    + pb_x[k] * mi_1175[k];

        t_1504[k] = pa_x[k] * lk_1504[k];

        t_1505[k] = pa_x[k] * lk_1505[k];

        t_1506[k] = pa_x[k] * lk_1506[k];

        t_1507[k] = pa_x[k] * lk_1507[k];

        t_1508[k] = pa_x[k] * lk_1508[k];
    }

#pragma omp simd aligned(t_1509, t_1510, t_1511, t_1512, t_1513, pa_x, pb_y, li_952, li_1176, \
                         lk_1509, lk_1510, lk_1511, lk_1512, mi_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1509[k] = pa_x[k] * lk_1509[k];

        t_1510[k] = pa_x[k] * lk_1510[k];

        t_1511[k] = pa_x[k] * lk_1511[k];

        t_1512[k] = f_17 * li_1176[k]
                    + pa_x[k] * lk_1512[k];

        t_1513[k] = f_12 * li_952[k]
                    + pb_y[k] * mi_1176[k];
    }

#pragma omp simd aligned(t_1514, t_1515, t_1516, t_1517, pa_x, pb_y, pb_z, li_924, li_954, \
                         li_1179, li_1181, lk_1515, lk_1517, mi_1176, \
                         mi_1178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1514[k] = f_24 * li_924[k]
                    + pb_z[k] * mi_1176[k];

        t_1515[k] = f_15 * li_1179[k]
                    + pa_x[k] * lk_1515[k];

        t_1516[k] = f_12 * li_954[k]
                    + pb_y[k] * mi_1178[k];

        t_1517[k] = f_15 * li_1181[k]
                    + pa_x[k] * lk_1517[k];
    }

#pragma omp simd aligned(t_1518, t_1519, t_1520, t_1521, pa_x, pb_y, pb_z, li_927, li_957, \
                         li_1182, li_1185, lk_1518, lk_1521, mi_1179, \
                         mi_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1518[k] = f_14 * li_1182[k]
                    + pa_x[k] * lk_1518[k];

        t_1519[k] = f_24 * li_927[k]
                    + pb_z[k] * mi_1179[k];

        t_1520[k] = f_12 * li_957[k]
                    + pb_y[k] * mi_1181[k];

        t_1521[k] = f_14 * li_1185[k]
                    + pa_x[k] * lk_1521[k];
    }

#pragma omp simd aligned(t_1522, t_1523, t_1524, t_1525, pa_x, pb_y, pb_z, li_930, li_961, \
                         li_1186, li_1188, lk_1522, lk_1524, mi_1182, \
                         mi_1185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1522[k] = f_13 * li_1186[k]
                    + pa_x[k] * lk_1522[k];

        t_1523[k] = f_24 * li_930[k]
                    + pb_z[k] * mi_1182[k];

        t_1524[k] = f_13 * li_1188[k]
                    + pa_x[k] * lk_1524[k];

        t_1525[k] = f_12 * li_961[k]
                    + pb_y[k] * mi_1185[k];
    }

#pragma omp simd aligned(t_1526, t_1527, t_1528, t_1529, pa_x, pb_z, li_934, li_1190, li_1191, \
                         li_1193, lk_1526, lk_1527, lk_1529, mi_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1526[k] = f_13 * li_1190[k]
                    + pa_x[k] * lk_1526[k];

        t_1527[k] = f_12 * li_1191[k]
                    + pa_x[k] * lk_1527[k];

        t_1528[k] = f_24 * li_934[k]
                    + pb_z[k] * mi_1186[k];

        t_1529[k] = f_12 * li_1193[k]
                    + pa_x[k] * lk_1529[k];
    }

#pragma omp simd aligned(t_1530, t_1531, t_1532, t_1533, pa_x, pb_x, pb_y, li_966, li_1194, \
                         li_1196, li_1197, lk_1530, lk_1532, mi_1190, \
                         mi_1197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1530[k] = f_12 * li_1194[k]
                    + pa_x[k] * lk_1530[k];

        t_1531[k] = f_12 * li_966[k]
                    + pb_y[k] * mi_1190[k];

        t_1532[k] = f_12 * li_1196[k]
                    + pa_x[k] * lk_1532[k];

        t_1533[k] = f_11 * li_1197[k]
                    + pb_x[k] * mi_1197[k];
    }

#pragma omp simd aligned(t_1534, t_1535, t_1536, t_1537, t_1538, pb_x, li_1198, li_1199, \
                         li_1200, li_1201, li_1202, mi_1198, mi_1199, mi_1200, mi_1201, \
                         mi_1202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1534[k] = f_11 * li_1198[k]
                    + pb_x[k] * mi_1198[k];

        t_1535[k] = f_11 * li_1199[k]
                    + pb_x[k] * mi_1199[k];

        t_1536[k] = f_11 * li_1200[k]
                    + pb_x[k] * mi_1200[k];

        t_1537[k] = f_11 * li_1201[k]
                    + pb_x[k] * mi_1201[k];

        t_1538[k] = f_11 * li_1202[k]
                    + pb_x[k] * mi_1202[k];
    }

#pragma omp simd aligned(t_1539, t_1540, t_1541, t_1542, t_1543, t_1544, pa_x, pb_x, li_1203, \
                         lk_1540, lk_1541, lk_1542, lk_1543, lk_1544, \
                         mi_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1539[k] = f_11 * li_1203[k]
                    + pb_x[k] * mi_1203[k];

        t_1540[k] = pa_x[k] * lk_1540[k];

        t_1541[k] = pa_x[k] * lk_1541[k];

        t_1542[k] = pa_x[k] * lk_1542[k];

        t_1543[k] = pa_x[k] * lk_1543[k];

        t_1544[k] = pa_x[k] * lk_1544[k];
    }

#pragma omp simd aligned(t_1545, t_1546, t_1547, t_1548, t_1549, t_1550, pa_x, pa_y, pb_y, \
                         li_980, lk_1260, lk_1262, lk_1545, lk_1546, lk_1547, \
                         mi_1204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1545[k] = pa_x[k] * lk_1545[k];

        t_1546[k] = pa_x[k] * lk_1546[k];

        t_1547[k] = pa_x[k] * lk_1547[k];

        t_1548[k] = pa_y[k] * lk_1260[k];

        t_1549[k] = f_11 * li_980[k]
                    + pb_y[k] * mi_1204[k];

        t_1550[k] = pa_y[k] * lk_1262[k];
    }

#pragma omp simd aligned(t_1551, t_1552, t_1553, t_1554, pa_x, pa_y, pb_y, li_982, li_1207, \
                         li_1210, lk_1265, lk_1551, lk_1554, mi_1206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1551[k] = f_15 * li_1207[k]
                    + pa_x[k] * lk_1551[k];

        t_1552[k] = f_11 * li_982[k]
                    + pb_y[k] * mi_1206[k];

        t_1553[k] = pa_y[k] * lk_1265[k];

        t_1554[k] = f_14 * li_1210[k]
                    + pa_x[k] * lk_1554[k];
    }

#pragma omp simd aligned(t_1555, t_1556, t_1557, t_1558, pa_x, pa_y, pb_y, pb_z, li_955, \
                         li_985, li_1214, lk_1269, lk_1558, mi_1207, \
                         mi_1209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1555[k] = f_17 * li_955[k]
                    + pb_z[k] * mi_1207[k];

        t_1556[k] = f_11 * li_985[k]
                    + pb_y[k] * mi_1209[k];

        t_1557[k] = pa_y[k] * lk_1269[k];

        t_1558[k] = f_13 * li_1214[k]
                    + pa_x[k] * lk_1558[k];
    }

#pragma omp simd aligned(t_1559, t_1560, t_1561, t_1562, pa_x, pa_y, pb_y, pb_z, li_958, \
                         li_989, li_1216, lk_1274, lk_1560, mi_1210, \
                         mi_1213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1559[k] = f_17 * li_958[k]
                    + pb_z[k] * mi_1210[k];

        t_1560[k] = f_13 * li_1216[k]
                    + pa_x[k] * lk_1560[k];

        t_1561[k] = f_11 * li_989[k]
                    + pb_y[k] * mi_1213[k];

        t_1562[k] = pa_y[k] * lk_1274[k];
    }

#pragma omp simd aligned(t_1563, t_1564, t_1565, t_1566, pa_x, pb_z, li_962, li_1219, li_1221, \
                         li_1222, lk_1563, lk_1565, lk_1566, mi_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1563[k] = f_12 * li_1219[k]
                    + pa_x[k] * lk_1563[k];

        t_1564[k] = f_17 * li_962[k]
                    + pb_z[k] * mi_1214[k];

        t_1565[k] = f_12 * li_1221[k]
                    + pa_x[k] * lk_1565[k];

        t_1566[k] = f_12 * li_1222[k]
                    + pa_x[k] * lk_1566[k];
    }

#pragma omp simd aligned(t_1567, t_1568, t_1569, t_1570, pa_y, pb_x, pb_y, li_994, li_1225, \
                         li_1226, lk_1280, mi_1218, mi_1225, mi_1226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1567[k] = f_11 * li_994[k]
                    + pb_y[k] * mi_1218[k];

        t_1568[k] = pa_y[k] * lk_1280[k];

        t_1569[k] = f_11 * li_1225[k]
                    + pb_x[k] * mi_1225[k];

        t_1570[k] = f_11 * li_1226[k]
                    + pb_x[k] * mi_1226[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece12(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t li,
                                             const size_t lk, const size_t mh0, const size_t mh1,
                                             const size_t mi, const size_t ncols,
                                             const double alpha, const double beta,
                                             const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_16 = 4.0 / p;
    const auto f_17 = 3.5 / p;

    auto *t_1571 = buffer.data(target + 1571);
    auto *t_1572 = buffer.data(target + 1572);
    auto *t_1573 = buffer.data(target + 1573);
    auto *t_1574 = buffer.data(target + 1574);
    auto *t_1575 = buffer.data(target + 1575);
    auto *t_1576 = buffer.data(target + 1576);
    auto *t_1577 = buffer.data(target + 1577);
    auto *t_1578 = buffer.data(target + 1578);
    auto *t_1579 = buffer.data(target + 1579);
    auto *t_1580 = buffer.data(target + 1580);
    auto *t_1581 = buffer.data(target + 1581);
    auto *t_1582 = buffer.data(target + 1582);
    auto *t_1583 = buffer.data(target + 1583);
    auto *t_1584 = buffer.data(target + 1584);
    auto *t_1585 = buffer.data(target + 1585);
    auto *t_1586 = buffer.data(target + 1586);
    auto *t_1587 = buffer.data(target + 1587);
    auto *t_1588 = buffer.data(target + 1588);
    auto *t_1589 = buffer.data(target + 1589);
    auto *t_1590 = buffer.data(target + 1590);
    auto *t_1591 = buffer.data(target + 1591);
    auto *t_1592 = buffer.data(target + 1592);
    auto *t_1593 = buffer.data(target + 1593);
    auto *t_1594 = buffer.data(target + 1594);
    auto *t_1595 = buffer.data(target + 1595);
    auto *t_1596 = buffer.data(target + 1596);
    auto *t_1597 = buffer.data(target + 1597);
    auto *t_1598 = buffer.data(target + 1598);
    auto *t_1599 = buffer.data(target + 1599);
    auto *t_1600 = buffer.data(target + 1600);
    auto *t_1601 = buffer.data(target + 1601);
    auto *t_1602 = buffer.data(target + 1602);
    auto *t_1603 = buffer.data(target + 1603);
    auto *t_1604 = buffer.data(target + 1604);
    auto *t_1605 = buffer.data(target + 1605);
    auto *t_1606 = buffer.data(target + 1606);
    auto *t_1607 = buffer.data(target + 1607);
    auto *t_1608 = buffer.data(target + 1608);
    auto *t_1609 = buffer.data(target + 1609);
    auto *t_1610 = buffer.data(target + 1610);
    auto *t_1611 = buffer.data(target + 1611);
    auto *t_1612 = buffer.data(target + 1612);
    auto *t_1613 = buffer.data(target + 1613);
    auto *t_1614 = buffer.data(target + 1614);
    auto *t_1615 = buffer.data(target + 1615);
    auto *t_1616 = buffer.data(target + 1616);
    auto *t_1617 = buffer.data(target + 1617);
    auto *t_1618 = buffer.data(target + 1618);
    auto *t_1619 = buffer.data(target + 1619);
    auto *t_1620 = buffer.data(target + 1620);
    auto *t_1621 = buffer.data(target + 1621);
    auto *t_1622 = buffer.data(target + 1622);
    auto *t_1623 = buffer.data(target + 1623);
    auto *t_1624 = buffer.data(target + 1624);
    auto *t_1625 = buffer.data(target + 1625);
    auto *t_1626 = buffer.data(target + 1626);
    auto *t_1627 = buffer.data(target + 1627);
    auto *t_1628 = buffer.data(target + 1628);
    auto *t_1629 = buffer.data(target + 1629);
    auto *t_1630 = buffer.data(target + 1630);
    auto *t_1631 = buffer.data(target + 1631);
    auto *t_1632 = buffer.data(target + 1632);
    auto *t_1633 = buffer.data(target + 1633);
    auto *t_1634 = buffer.data(target + 1634);
    auto *t_1635 = buffer.data(target + 1635);
    auto *t_1636 = buffer.data(target + 1636);
    auto *t_1637 = buffer.data(target + 1637);
    auto *t_1638 = buffer.data(target + 1638);
    auto *t_1639 = buffer.data(target + 1639);
    auto *t_1640 = buffer.data(target + 1640);
    auto *t_1641 = buffer.data(target + 1641);
    auto *t_1642 = buffer.data(target + 1642);
    auto *t_1643 = buffer.data(target + 1643);
    auto *t_1644 = buffer.data(target + 1644);
    auto *t_1645 = buffer.data(target + 1645);
    auto *t_1646 = buffer.data(target + 1646);
    auto *t_1647 = buffer.data(target + 1647);
    auto *t_1648 = buffer.data(target + 1648);
    auto *t_1649 = buffer.data(target + 1649);
    auto *t_1650 = buffer.data(target + 1650);
    auto *t_1651 = buffer.data(target + 1651);
    auto *t_1652 = buffer.data(target + 1652);
    auto *t_1653 = buffer.data(target + 1653);
    auto *t_1654 = buffer.data(target + 1654);
    auto *t_1655 = buffer.data(target + 1655);
    auto *t_1656 = buffer.data(target + 1656);
    auto *t_1657 = buffer.data(target + 1657);
    auto *t_1658 = buffer.data(target + 1658);
    auto *t_1659 = buffer.data(target + 1659);
    auto *t_1660 = buffer.data(target + 1660);
    auto *t_1661 = buffer.data(target + 1661);
    auto *t_1662 = buffer.data(target + 1662);
    auto *t_1663 = buffer.data(target + 1663);
    auto *t_1664 = buffer.data(target + 1664);
    auto *t_1665 = buffer.data(target + 1665);
    auto *t_1666 = buffer.data(target + 1666);
    auto *t_1667 = buffer.data(target + 1667);
    auto *t_1668 = buffer.data(target + 1668);
    auto *t_1669 = buffer.data(target + 1669);
    auto *t_1670 = buffer.data(target + 1670);
    auto *t_1671 = buffer.data(target + 1671);
    auto *t_1672 = buffer.data(target + 1672);
    auto *t_1673 = buffer.data(target + 1673);
    auto *t_1674 = buffer.data(target + 1674);
    auto *t_1675 = buffer.data(target + 1675);
    auto *t_1676 = buffer.data(target + 1676);
    auto *t_1677 = buffer.data(target + 1677);
    auto *t_1678 = buffer.data(target + 1678);
    auto *t_1679 = buffer.data(target + 1679);
    auto *t_1680 = buffer.data(target + 1680);
    auto *t_1681 = buffer.data(target + 1681);
    auto *t_1682 = buffer.data(target + 1682);
    auto *t_1683 = buffer.data(target + 1683);
    auto *t_1684 = buffer.data(target + 1684);
    auto *t_1685 = buffer.data(target + 1685);
    auto *t_1686 = buffer.data(target + 1686);
    auto *t_1687 = buffer.data(target + 1687);
    auto *t_1688 = buffer.data(target + 1688);
    auto *t_1689 = buffer.data(target + 1689);
    auto *t_1690 = buffer.data(target + 1690);
    auto *t_1691 = buffer.data(target + 1691);
    auto *t_1692 = buffer.data(target + 1692);
    auto *t_1693 = buffer.data(target + 1693);
    auto *t_1694 = buffer.data(target + 1694);
    auto *t_1695 = buffer.data(target + 1695);
    auto *t_1696 = buffer.data(target + 1696);
    auto *t_1697 = buffer.data(target + 1697);
    auto *t_1698 = buffer.data(target + 1698);
    auto *t_1699 = buffer.data(target + 1699);
    auto *t_1700 = buffer.data(target + 1700);
    auto *t_1701 = buffer.data(target + 1701);
    auto *t_1702 = buffer.data(target + 1702);
    auto *t_1703 = buffer.data(target + 1703);
    auto *t_1704 = buffer.data(target + 1704);
    auto *t_1705 = buffer.data(target + 1705);
    auto *t_1706 = buffer.data(target + 1706);
    auto *t_1707 = buffer.data(target + 1707);
    auto *t_1708 = buffer.data(target + 1708);
    auto *t_1709 = buffer.data(target + 1709);
    auto *t_1710 = buffer.data(target + 1710);
    auto *t_1711 = buffer.data(target + 1711);
    auto *t_1712 = buffer.data(target + 1712);
    auto *t_1713 = buffer.data(target + 1713);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *li_980 = buffer.data(li + 980);
    const auto *li_983 = buffer.data(li + 983);
    const auto *li_986 = buffer.data(li + 986);
    const auto *li_990 = buffer.data(li + 990);
    const auto *li_1008 = buffer.data(li + 1008);
    const auto *li_1010 = buffer.data(li + 1010);
    const auto *li_1011 = buffer.data(li + 1011);
    const auto *li_1013 = buffer.data(li + 1013);
    const auto *li_1014 = buffer.data(li + 1014);
    const auto *li_1015 = buffer.data(li + 1015);
    const auto *li_1017 = buffer.data(li + 1017);
    const auto *li_1018 = buffer.data(li + 1018);
    const auto *li_1019 = buffer.data(li + 1019);
    const auto *li_1020 = buffer.data(li + 1020);
    const auto *li_1022 = buffer.data(li + 1022);
    const auto *li_1029 = buffer.data(li + 1029);
    const auto *li_1030 = buffer.data(li + 1030);
    const auto *li_1031 = buffer.data(li + 1031);
    const auto *li_1032 = buffer.data(li + 1032);
    const auto *li_1033 = buffer.data(li + 1033);
    const auto *li_1035 = buffer.data(li + 1035);
    const auto *li_1036 = buffer.data(li + 1036);
    const auto *li_1038 = buffer.data(li + 1038);
    const auto *li_1039 = buffer.data(li + 1039);
    const auto *li_1041 = buffer.data(li + 1041);
    const auto *li_1042 = buffer.data(li + 1042);
    const auto *li_1045 = buffer.data(li + 1045);
    const auto *li_1046 = buffer.data(li + 1046);
    const auto *li_1050 = buffer.data(li + 1050);
    const auto *li_1063 = buffer.data(li + 1063);
    const auto *li_1064 = buffer.data(li + 1064);
    const auto *li_1066 = buffer.data(li + 1066);
    const auto *li_1069 = buffer.data(li + 1069);
    const auto *li_1073 = buffer.data(li + 1073);
    const auto *li_1078 = buffer.data(li + 1078);
    const auto *li_1227 = buffer.data(li + 1227);
    const auto *li_1228 = buffer.data(li + 1228);
    const auto *li_1229 = buffer.data(li + 1229);
    const auto *li_1230 = buffer.data(li + 1230);
    const auto *li_1232 = buffer.data(li + 1232);
    const auto *li_1235 = buffer.data(li + 1235);
    const auto *li_1237 = buffer.data(li + 1237);
    const auto *li_1238 = buffer.data(li + 1238);
    const auto *li_1241 = buffer.data(li + 1241);
    const auto *li_1242 = buffer.data(li + 1242);
    const auto *li_1244 = buffer.data(li + 1244);
    const auto *li_1246 = buffer.data(li + 1246);
    const auto *li_1247 = buffer.data(li + 1247);
    const auto *li_1249 = buffer.data(li + 1249);
    const auto *li_1250 = buffer.data(li + 1250);
    const auto *li_1252 = buffer.data(li + 1252);
    const auto *li_1253 = buffer.data(li + 1253);
    const auto *li_1254 = buffer.data(li + 1254);
    const auto *li_1255 = buffer.data(li + 1255);
    const auto *li_1256 = buffer.data(li + 1256);
    const auto *li_1257 = buffer.data(li + 1257);
    const auto *li_1259 = buffer.data(li + 1259);

    const auto *lk_1287 = buffer.data(lk + 1287);
    const auto *lk_1296 = buffer.data(lk + 1296);
    const auto *lk_1297 = buffer.data(lk + 1297);
    const auto *lk_1299 = buffer.data(lk + 1299);
    const auto *lk_1301 = buffer.data(lk + 1301);
    const auto *lk_1302 = buffer.data(lk + 1302);
    const auto *lk_1305 = buffer.data(lk + 1305);
    const auto *lk_1306 = buffer.data(lk + 1306);
    const auto *lk_1308 = buffer.data(lk + 1308);
    const auto *lk_1310 = buffer.data(lk + 1310);
    const auto *lk_1311 = buffer.data(lk + 1311);
    const auto *lk_1313 = buffer.data(lk + 1313);
    const auto *lk_1314 = buffer.data(lk + 1314);
    const auto *lk_1316 = buffer.data(lk + 1316);
    const auto *lk_1324 = buffer.data(lk + 1324);
    const auto *lk_1326 = buffer.data(lk + 1326);
    const auto *lk_1327 = buffer.data(lk + 1327);
    const auto *lk_1328 = buffer.data(lk + 1328);
    const auto *lk_1329 = buffer.data(lk + 1329);
    const auto *lk_1331 = buffer.data(lk + 1331);
    const auto *lk_1576 = buffer.data(lk + 1576);
    const auto *lk_1577 = buffer.data(lk + 1577);
    const auto *lk_1578 = buffer.data(lk + 1578);
    const auto *lk_1579 = buffer.data(lk + 1579);
    const auto *lk_1580 = buffer.data(lk + 1580);
    const auto *lk_1581 = buffer.data(lk + 1581);
    const auto *lk_1582 = buffer.data(lk + 1582);
    const auto *lk_1583 = buffer.data(lk + 1583);
    const auto *lk_1584 = buffer.data(lk + 1584);
    const auto *lk_1587 = buffer.data(lk + 1587);
    const auto *lk_1589 = buffer.data(lk + 1589);
    const auto *lk_1590 = buffer.data(lk + 1590);
    const auto *lk_1593 = buffer.data(lk + 1593);
    const auto *lk_1594 = buffer.data(lk + 1594);
    const auto *lk_1596 = buffer.data(lk + 1596);
    const auto *lk_1598 = buffer.data(lk + 1598);
    const auto *lk_1599 = buffer.data(lk + 1599);
    const auto *lk_1601 = buffer.data(lk + 1601);
    const auto *lk_1602 = buffer.data(lk + 1602);
    const auto *lk_1604 = buffer.data(lk + 1604);
    const auto *lk_1612 = buffer.data(lk + 1612);
    const auto *lk_1613 = buffer.data(lk + 1613);
    const auto *lk_1614 = buffer.data(lk + 1614);
    const auto *lk_1615 = buffer.data(lk + 1615);
    const auto *lk_1616 = buffer.data(lk + 1616);
    const auto *lk_1617 = buffer.data(lk + 1617);
    const auto *lk_1619 = buffer.data(lk + 1619);

    const auto *mh0_945 = buffer.data(mh0 + 945);
    const auto *mh0_948 = buffer.data(mh0 + 948);
    const auto *mh0_950 = buffer.data(mh0 + 950);
    const auto *mh0_951 = buffer.data(mh0 + 951);
    const auto *mh0_954 = buffer.data(mh0 + 954);
    const auto *mh0_955 = buffer.data(mh0 + 955);
    const auto *mh0_957 = buffer.data(mh0 + 957);
    const auto *mh0_959 = buffer.data(mh0 + 959);
    const auto *mh0_960 = buffer.data(mh0 + 960);
    const auto *mh0_961 = buffer.data(mh0 + 961);
    const auto *mh0_962 = buffer.data(mh0 + 962);
    const auto *mh0_963 = buffer.data(mh0 + 963);
    const auto *mh0_965 = buffer.data(mh0 + 965);
    const auto *mh0_987 = buffer.data(mh0 + 987);
    const auto *mh0_990 = buffer.data(mh0 + 990);
    const auto *mh0_992 = buffer.data(mh0 + 992);
    const auto *mh0_993 = buffer.data(mh0 + 993);
    const auto *mh0_996 = buffer.data(mh0 + 996);
    const auto *mh0_997 = buffer.data(mh0 + 997);
    const auto *mh0_999 = buffer.data(mh0 + 999);
    const auto *mh0_1001 = buffer.data(mh0 + 1001);
    const auto *mh0_1002 = buffer.data(mh0 + 1002);
    const auto *mh0_1004 = buffer.data(mh0 + 1004);
    const auto *mh0_1005 = buffer.data(mh0 + 1005);
    const auto *mh0_1007 = buffer.data(mh0 + 1007);

    const auto *mh1_945 = buffer.data(mh1 + 945);
    const auto *mh1_948 = buffer.data(mh1 + 948);
    const auto *mh1_950 = buffer.data(mh1 + 950);
    const auto *mh1_951 = buffer.data(mh1 + 951);
    const auto *mh1_954 = buffer.data(mh1 + 954);
    const auto *mh1_955 = buffer.data(mh1 + 955);
    const auto *mh1_957 = buffer.data(mh1 + 957);
    const auto *mh1_959 = buffer.data(mh1 + 959);
    const auto *mh1_960 = buffer.data(mh1 + 960);
    const auto *mh1_961 = buffer.data(mh1 + 961);
    const auto *mh1_962 = buffer.data(mh1 + 962);
    const auto *mh1_963 = buffer.data(mh1 + 963);
    const auto *mh1_965 = buffer.data(mh1 + 965);
    const auto *mh1_987 = buffer.data(mh1 + 987);
    const auto *mh1_990 = buffer.data(mh1 + 990);
    const auto *mh1_992 = buffer.data(mh1 + 992);
    const auto *mh1_993 = buffer.data(mh1 + 993);
    const auto *mh1_996 = buffer.data(mh1 + 996);
    const auto *mh1_997 = buffer.data(mh1 + 997);
    const auto *mh1_999 = buffer.data(mh1 + 999);
    const auto *mh1_1001 = buffer.data(mh1 + 1001);
    const auto *mh1_1002 = buffer.data(mh1 + 1002);
    const auto *mh1_1004 = buffer.data(mh1 + 1004);
    const auto *mh1_1005 = buffer.data(mh1 + 1005);
    const auto *mh1_1007 = buffer.data(mh1 + 1007);

    const auto *mi_1227 = buffer.data(mi + 1227);
    const auto *mi_1228 = buffer.data(mi + 1228);
    const auto *mi_1229 = buffer.data(mi + 1229);
    const auto *mi_1230 = buffer.data(mi + 1230);
    const auto *mi_1232 = buffer.data(mi + 1232);
    const auto *mi_1234 = buffer.data(mi + 1234);
    const auto *mi_1235 = buffer.data(mi + 1235);
    const auto *mi_1237 = buffer.data(mi + 1237);
    const auto *mi_1238 = buffer.data(mi + 1238);
    const auto *mi_1241 = buffer.data(mi + 1241);
    const auto *mi_1242 = buffer.data(mi + 1242);
    const auto *mi_1246 = buffer.data(mi + 1246);
    const auto *mi_1252 = buffer.data(mi + 1252);
    const auto *mi_1253 = buffer.data(mi + 1253);
    const auto *mi_1254 = buffer.data(mi + 1254);
    const auto *mi_1255 = buffer.data(mi + 1255);
    const auto *mi_1256 = buffer.data(mi + 1256);
    const auto *mi_1257 = buffer.data(mi + 1257);
    const auto *mi_1259 = buffer.data(mi + 1259);
    const auto *mi_1260 = buffer.data(mi + 1260);
    const auto *mi_1261 = buffer.data(mi + 1261);
    const auto *mi_1263 = buffer.data(mi + 1263);
    const auto *mi_1265 = buffer.data(mi + 1265);
    const auto *mi_1266 = buffer.data(mi + 1266);
    const auto *mi_1269 = buffer.data(mi + 1269);
    const auto *mi_1270 = buffer.data(mi + 1270);
    const auto *mi_1272 = buffer.data(mi + 1272);
    const auto *mi_1274 = buffer.data(mi + 1274);
    const auto *mi_1275 = buffer.data(mi + 1275);
    const auto *mi_1277 = buffer.data(mi + 1277);
    const auto *mi_1278 = buffer.data(mi + 1278);
    const auto *mi_1280 = buffer.data(mi + 1280);
    const auto *mi_1281 = buffer.data(mi + 1281);
    const auto *mi_1282 = buffer.data(mi + 1282);
    const auto *mi_1283 = buffer.data(mi + 1283);
    const auto *mi_1284 = buffer.data(mi + 1284);
    const auto *mi_1285 = buffer.data(mi + 1285);
    const auto *mi_1286 = buffer.data(mi + 1286);
    const auto *mi_1287 = buffer.data(mi + 1287);
    const auto *mi_1288 = buffer.data(mi + 1288);
    const auto *mi_1290 = buffer.data(mi + 1290);
    const auto *mi_1291 = buffer.data(mi + 1291);
    const auto *mi_1293 = buffer.data(mi + 1293);
    const auto *mi_1294 = buffer.data(mi + 1294);
    const auto *mi_1297 = buffer.data(mi + 1297);
    const auto *mi_1298 = buffer.data(mi + 1298);
    const auto *mi_1302 = buffer.data(mi + 1302);
    const auto *mi_1309 = buffer.data(mi + 1309);
    const auto *mi_1310 = buffer.data(mi + 1310);
    const auto *mi_1311 = buffer.data(mi + 1311);
    const auto *mi_1312 = buffer.data(mi + 1312);
    const auto *mi_1313 = buffer.data(mi + 1313);
    const auto *mi_1314 = buffer.data(mi + 1314);
    const auto *mi_1315 = buffer.data(mi + 1315);
    const auto *mi_1316 = buffer.data(mi + 1316);
    const auto *mi_1318 = buffer.data(mi + 1318);
    const auto *mi_1319 = buffer.data(mi + 1319);
    const auto *mi_1321 = buffer.data(mi + 1321);
    const auto *mi_1322 = buffer.data(mi + 1322);
    const auto *mi_1325 = buffer.data(mi + 1325);
    const auto *mi_1326 = buffer.data(mi + 1326);
    const auto *mi_1328 = buffer.data(mi + 1328);
    const auto *mi_1330 = buffer.data(mi + 1330);
    const auto *mi_1331 = buffer.data(mi + 1331);
    const auto *mi_1333 = buffer.data(mi + 1333);
    const auto *mi_1334 = buffer.data(mi + 1334);
    const auto *mi_1336 = buffer.data(mi + 1336);
    const auto *mi_1337 = buffer.data(mi + 1337);

#pragma omp simd aligned(t_1571, t_1572, t_1573, t_1574, t_1575, pa_y, pb_x, li_1227, li_1228, \
                         li_1229, li_1230, lk_1287, mi_1227, mi_1228, mi_1229, \
                         mi_1230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1571[k] = f_11 * li_1227[k]
                    + pb_x[k] * mi_1227[k];

        t_1572[k] = f_11 * li_1228[k]
                    + pb_x[k] * mi_1228[k];

        t_1573[k] = f_11 * li_1229[k]
                    + pb_x[k] * mi_1229[k];

        t_1574[k] = f_11 * li_1230[k]
                    + pb_x[k] * mi_1230[k];

        t_1575[k] = pa_y[k] * lk_1287[k];
    }

#pragma omp simd aligned(t_1576, t_1577, t_1578, t_1579, t_1580, t_1581, t_1582, pa_x, \
                         lk_1576, lk_1577, lk_1578, lk_1579, lk_1580, lk_1581, \
                         lk_1582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1576[k] = pa_x[k] * lk_1576[k];

        t_1577[k] = pa_x[k] * lk_1577[k];

        t_1578[k] = pa_x[k] * lk_1578[k];

        t_1579[k] = pa_x[k] * lk_1579[k];

        t_1580[k] = pa_x[k] * lk_1580[k];

        t_1581[k] = pa_x[k] * lk_1581[k];

        t_1582[k] = pa_x[k] * lk_1582[k];
    }

#pragma omp simd aligned(t_1583, t_1584, t_1585, t_1586, t_1587, pa_x, pb_y, pb_z, li_980, \
                         li_1232, li_1235, lk_1583, lk_1584, lk_1587, \
                         mi_1232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1583[k] = pa_x[k] * lk_1583[k];

        t_1584[k] = f_17 * li_1232[k]
                    + pa_x[k] * lk_1584[k];

        t_1585[k] = pb_y[k] * mi_1232[k];

        t_1586[k] = f_16 * li_980[k]
                    + pb_z[k] * mi_1232[k];

        t_1587[k] = f_15 * li_1235[k]
                    + pa_x[k] * lk_1587[k];
    }

#pragma omp simd aligned(t_1588, t_1589, t_1590, t_1591, t_1592, pa_x, pb_y, pb_z, li_983, \
                         li_1237, li_1238, lk_1589, lk_1590, mi_1234, mi_1235, \
                         mi_1237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1588[k] = pb_y[k] * mi_1234[k];

        t_1589[k] = f_15 * li_1237[k]
                    + pa_x[k] * lk_1589[k];

        t_1590[k] = f_14 * li_1238[k]
                    + pa_x[k] * lk_1590[k];

        t_1591[k] = f_16 * li_983[k]
                    + pb_z[k] * mi_1235[k];

        t_1592[k] = pb_y[k] * mi_1237[k];
    }

#pragma omp simd aligned(t_1593, t_1594, t_1595, t_1596, pa_x, pb_z, li_986, li_1241, li_1242, \
                         li_1244, lk_1593, lk_1594, lk_1596, mi_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1593[k] = f_14 * li_1241[k]
                    + pa_x[k] * lk_1593[k];

        t_1594[k] = f_13 * li_1242[k]
                    + pa_x[k] * lk_1594[k];

        t_1595[k] = f_16 * li_986[k]
                    + pb_z[k] * mi_1238[k];

        t_1596[k] = f_13 * li_1244[k]
                    + pa_x[k] * lk_1596[k];
    }

#pragma omp simd aligned(t_1597, t_1598, t_1599, t_1600, pa_x, pb_y, pb_z, li_990, li_1246, \
                         li_1247, lk_1598, lk_1599, mi_1241, mi_1242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1597[k] = pb_y[k] * mi_1241[k];

        t_1598[k] = f_13 * li_1246[k]
                    + pa_x[k] * lk_1598[k];

        t_1599[k] = f_12 * li_1247[k]
                    + pa_x[k] * lk_1599[k];

        t_1600[k] = f_16 * li_990[k]
                    + pb_z[k] * mi_1242[k];
    }

#pragma omp simd aligned(t_1601, t_1602, t_1603, t_1604, pa_x, pb_y, li_1249, li_1250, \
                         li_1252, lk_1601, lk_1602, lk_1604, mi_1246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1601[k] = f_12 * li_1249[k]
                    + pa_x[k] * lk_1601[k];

        t_1602[k] = f_12 * li_1250[k]
                    + pa_x[k] * lk_1602[k];

        t_1603[k] = pb_y[k] * mi_1246[k];

        t_1604[k] = f_12 * li_1252[k]
                    + pa_x[k] * lk_1604[k];
    }

#pragma omp simd aligned(t_1605, t_1606, t_1607, t_1608, t_1609, pb_x, li_1253, li_1254, \
                         li_1255, li_1256, li_1257, mi_1253, mi_1254, mi_1255, mi_1256, \
                         mi_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1605[k] = f_11 * li_1253[k]
                    + pb_x[k] * mi_1253[k];

        t_1606[k] = f_11 * li_1254[k]
                    + pb_x[k] * mi_1254[k];

        t_1607[k] = f_11 * li_1255[k]
                    + pb_x[k] * mi_1255[k];

        t_1608[k] = f_11 * li_1256[k]
                    + pb_x[k] * mi_1256[k];

        t_1609[k] = f_11 * li_1257[k]
                    + pb_x[k] * mi_1257[k];
    }

#pragma omp simd aligned(t_1610, t_1611, t_1612, t_1613, t_1614, t_1615, pa_x, pb_x, pb_y, \
                         li_1259, lk_1612, lk_1613, lk_1614, lk_1615, mi_1252, \
                         mi_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1610[k] = pb_y[k] * mi_1252[k];

        t_1611[k] = f_11 * li_1259[k]
                    + pb_x[k] * mi_1259[k];

        t_1612[k] = pa_x[k] * lk_1612[k];

        t_1613[k] = pa_x[k] * lk_1613[k];

        t_1614[k] = pa_x[k] * lk_1614[k];

        t_1615[k] = pa_x[k] * lk_1615[k];
    }

#pragma omp simd aligned(t_1616, t_1617, t_1618, t_1619, t_1620, pa_x, pb_x, pb_y, lk_1616, \
                         lk_1617, lk_1619, mh0_945, mh1_945, mi_1259, \
                         mi_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1616[k] = pa_x[k] * lk_1616[k];

        t_1617[k] = pa_x[k] * lk_1617[k];

        t_1618[k] = pb_y[k] * mi_1259[k];

        t_1619[k] = pa_x[k] * lk_1619[k];

        t_1620[k] = f_1 * mh0_945[k]
                    - f_2 * mh1_945[k]
                    + pb_x[k] * mi_1260[k];
    }

#pragma omp simd aligned(t_1621, t_1622, t_1623, t_1624, pb_x, pb_y, pb_z, li_1008, mh0_948, \
                         mh1_948, mi_1260, mi_1261, mi_1263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1621[k] = f_0 * li_1008[k]
                    + pb_y[k] * mi_1260[k];

        t_1622[k] = pb_z[k] * mi_1260[k];

        t_1623[k] = f_9 * mh0_948[k]
                    - f_10 * mh1_948[k]
                    + pb_x[k] * mi_1263[k];

        t_1624[k] = pb_z[k] * mi_1261[k];
    }

#pragma omp simd aligned(t_1625, t_1626, t_1627, t_1628, pb_x, pb_y, pb_z, li_1013, mh0_950, \
                         mh0_951, mh1_950, mh1_951, mi_1263, mi_1265, \
                         mi_1266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1625[k] = f_9 * mh0_950[k]
                    - f_10 * mh1_950[k]
                    + pb_x[k] * mi_1265[k];

        t_1626[k] = f_7 * mh0_951[k]
                    - f_8 * mh1_951[k]
                    + pb_x[k] * mi_1266[k];

        t_1627[k] = pb_z[k] * mi_1263[k];

        t_1628[k] = f_0 * li_1013[k]
                    + pb_y[k] * mi_1265[k];
    }

#pragma omp simd aligned(t_1629, t_1630, t_1631, t_1632, pb_x, pb_z, mh0_954, mh0_955, \
                         mh0_957, mh1_954, mh1_955, mh1_957, mi_1266, mi_1269, mi_1270, \
                         mi_1272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1629[k] = f_7 * mh0_954[k]
                    - f_8 * mh1_954[k]
                    + pb_x[k] * mi_1269[k];

        t_1630[k] = f_5 * mh0_955[k]
                    - f_6 * mh1_955[k]
                    + pb_x[k] * mi_1270[k];

        t_1631[k] = pb_z[k] * mi_1266[k];

        t_1632[k] = f_5 * mh0_957[k]
                    - f_6 * mh1_957[k]
                    + pb_x[k] * mi_1272[k];
    }

#pragma omp simd aligned(t_1633, t_1634, t_1635, t_1636, pb_x, pb_y, pb_z, li_1017, mh0_959, \
                         mh0_960, mh1_959, mh1_960, mi_1269, mi_1270, mi_1274, \
                         mi_1275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1633[k] = f_0 * li_1017[k]
                    + pb_y[k] * mi_1269[k];

        t_1634[k] = f_5 * mh0_959[k]
                    - f_6 * mh1_959[k]
                    + pb_x[k] * mi_1274[k];

        t_1635[k] = f_3 * mh0_960[k]
                    - f_4 * mh1_960[k]
                    + pb_x[k] * mi_1275[k];

        t_1636[k] = pb_z[k] * mi_1270[k];
    }

#pragma omp simd aligned(t_1637, t_1638, t_1639, pb_x, pb_y, li_1022, mh0_962, mh0_963, \
                         mh1_962, mh1_963, mi_1274, mi_1277, mi_1278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1637[k] = f_3 * mh0_962[k]
                    - f_4 * mh1_962[k]
                    + pb_x[k] * mi_1277[k];

        t_1638[k] = f_3 * mh0_963[k]
                    - f_4 * mh1_963[k]
                    + pb_x[k] * mi_1278[k];

        t_1639[k] = f_0 * li_1022[k]
                    + pb_y[k] * mi_1274[k];
    }

#pragma omp simd aligned(t_1640, t_1641, t_1642, t_1643, t_1644, t_1645, pb_x, mh0_965, \
                         mh1_965, mi_1280, mi_1281, mi_1282, mi_1283, mi_1284, \
                         mi_1285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1640[k] = f_3 * mh0_965[k]
                    - f_4 * mh1_965[k]
                    + pb_x[k] * mi_1280[k];

        t_1641[k] = pb_x[k] * mi_1281[k];

        t_1642[k] = pb_x[k] * mi_1282[k];

        t_1643[k] = pb_x[k] * mi_1283[k];

        t_1644[k] = pb_x[k] * mi_1284[k];

        t_1645[k] = pb_x[k] * mi_1285[k];
    }

#pragma omp simd aligned(t_1646, t_1647, t_1648, t_1649, t_1650, pb_x, pb_y, pb_z, li_1029, \
                         mh0_960, mh1_960, mi_1281, mi_1282, mi_1286, \
                         mi_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1646[k] = pb_x[k] * mi_1286[k];

        t_1647[k] = pb_x[k] * mi_1287[k];

        t_1648[k] = f_0 * li_1029[k]
                    + f_1 * mh0_960[k]
                    - f_2 * mh1_960[k]
                    + pb_y[k] * mi_1281[k];

        t_1649[k] = pb_z[k] * mi_1281[k];

        t_1650[k] = f_3 * mh0_960[k]
                    - f_4 * mh1_960[k]
                    + pb_z[k] * mi_1282[k];
    }

#pragma omp simd aligned(t_1651, t_1652, t_1653, pb_z, mh0_961, mh0_962, mh0_963, mh1_961, \
                         mh1_962, mh1_963, mi_1283, mi_1284, mi_1285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1651[k] = f_5 * mh0_961[k]
                    - f_6 * mh1_961[k]
                    + pb_z[k] * mi_1283[k];

        t_1652[k] = f_7 * mh0_962[k]
                    - f_8 * mh1_962[k]
                    + pb_z[k] * mi_1284[k];

        t_1653[k] = f_9 * mh0_963[k]
                    - f_10 * mh1_963[k]
                    + pb_z[k] * mi_1285[k];
    }

#pragma omp simd aligned(t_1654, t_1655, t_1656, t_1657, t_1658, pa_z, pb_y, pb_z, li_1008, \
                         li_1035, lk_1296, lk_1297, mh0_965, mh1_965, mi_1287, \
                         mi_1288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1654[k] = f_0 * li_1035[k]
                    + pb_y[k] * mi_1287[k];

        t_1655[k] = f_1 * mh0_965[k]
                    - f_2 * mh1_965[k]
                    + pb_z[k] * mi_1287[k];

        t_1656[k] = pa_z[k] * lk_1296[k];

        t_1657[k] = pa_z[k] * lk_1297[k];

        t_1658[k] = f_11 * li_1008[k]
                    + pb_z[k] * mi_1288[k];
    }

#pragma omp simd aligned(t_1659, t_1660, t_1661, t_1662, t_1663, pa_z, pb_y, pb_z, li_1010, \
                         li_1011, li_1038, lk_1299, lk_1301, lk_1302, mi_1290, \
                         mi_1291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1659[k] = pa_z[k] * lk_1299[k];

        t_1660[k] = f_16 * li_1038[k]
                    + pb_y[k] * mi_1290[k];

        t_1661[k] = f_12 * li_1010[k]
                    + pa_z[k] * lk_1301[k];

        t_1662[k] = pa_z[k] * lk_1302[k];

        t_1663[k] = f_11 * li_1011[k]
                    + pb_z[k] * mi_1291[k];
    }

#pragma omp simd aligned(t_1664, t_1665, t_1666, t_1667, pa_z, pb_y, pb_z, li_1013, li_1014, \
                         li_1041, lk_1305, lk_1306, mi_1293, mi_1294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1664[k] = f_16 * li_1041[k]
                    + pb_y[k] * mi_1293[k];

        t_1665[k] = f_13 * li_1013[k]
                    + pa_z[k] * lk_1305[k];

        t_1666[k] = pa_z[k] * lk_1306[k];

        t_1667[k] = f_11 * li_1014[k]
                    + pb_z[k] * mi_1294[k];
    }

#pragma omp simd aligned(t_1668, t_1669, t_1670, t_1671, pa_z, pb_y, li_1015, li_1017, \
                         li_1045, lk_1308, lk_1310, lk_1311, mi_1297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1668[k] = f_12 * li_1015[k]
                    + pa_z[k] * lk_1308[k];

        t_1669[k] = f_16 * li_1045[k]
                    + pb_y[k] * mi_1297[k];

        t_1670[k] = f_14 * li_1017[k]
                    + pa_z[k] * lk_1310[k];

        t_1671[k] = pa_z[k] * lk_1311[k];
    }

#pragma omp simd aligned(t_1672, t_1673, t_1674, t_1675, pa_z, pb_y, pb_z, li_1018, li_1019, \
                         li_1020, li_1050, lk_1313, lk_1314, mi_1298, \
                         mi_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1672[k] = f_11 * li_1018[k]
                    + pb_z[k] * mi_1298[k];

        t_1673[k] = f_12 * li_1019[k]
                    + pa_z[k] * lk_1313[k];

        t_1674[k] = f_13 * li_1020[k]
                    + pa_z[k] * lk_1314[k];

        t_1675[k] = f_16 * li_1050[k]
                    + pb_y[k] * mi_1302[k];
    }

#pragma omp simd aligned(t_1676, t_1677, t_1678, t_1679, t_1680, t_1681, pa_z, pb_x, li_1022, \
                         lk_1316, mi_1309, mi_1310, mi_1311, mi_1312, \
                         mi_1313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1676[k] = f_15 * li_1022[k]
                    + pa_z[k] * lk_1316[k];

        t_1677[k] = pb_x[k] * mi_1309[k];

        t_1678[k] = pb_x[k] * mi_1310[k];

        t_1679[k] = pb_x[k] * mi_1311[k];

        t_1680[k] = pb_x[k] * mi_1312[k];

        t_1681[k] = pb_x[k] * mi_1313[k];
    }

#pragma omp simd aligned(t_1682, t_1683, t_1684, t_1685, t_1686, pa_z, pb_x, pb_z, li_1029, \
                         li_1030, lk_1324, lk_1326, mi_1309, mi_1314, \
                         mi_1315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1682[k] = pb_x[k] * mi_1314[k];

        t_1683[k] = pb_x[k] * mi_1315[k];

        t_1684[k] = pa_z[k] * lk_1324[k];

        t_1685[k] = f_11 * li_1029[k]
                    + pb_z[k] * mi_1309[k];

        t_1686[k] = f_12 * li_1030[k]
                    + pa_z[k] * lk_1326[k];
    }

#pragma omp simd aligned(t_1687, t_1688, t_1689, t_1690, pa_z, pb_y, li_1031, li_1032, \
                         li_1033, li_1063, lk_1327, lk_1328, lk_1329, \
                         mi_1315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1687[k] = f_13 * li_1031[k]
                    + pa_z[k] * lk_1327[k];

        t_1688[k] = f_14 * li_1032[k]
                    + pa_z[k] * lk_1328[k];

        t_1689[k] = f_15 * li_1033[k]
                    + pa_z[k] * lk_1329[k];

        t_1690[k] = f_16 * li_1063[k]
                    + pb_y[k] * mi_1315[k];
    }

#pragma omp simd aligned(t_1691, t_1692, t_1693, t_1694, pa_z, pb_x, pb_y, pb_z, li_1035, \
                         li_1036, li_1064, lk_1331, mh0_987, mh1_987, \
                         mi_1316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1691[k] = f_17 * li_1035[k]
                    + pa_z[k] * lk_1331[k];

        t_1692[k] = f_1 * mh0_987[k]
                    - f_2 * mh1_987[k]
                    + pb_x[k] * mi_1316[k];

        t_1693[k] = f_17 * li_1064[k]
                    + pb_y[k] * mi_1316[k];

        t_1694[k] = f_12 * li_1036[k]
                    + pb_z[k] * mi_1316[k];
    }

#pragma omp simd aligned(t_1695, t_1696, t_1697, pb_x, pb_y, li_1066, mh0_990, mh0_992, \
                         mh1_990, mh1_992, mi_1318, mi_1319, mi_1321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1695[k] = f_9 * mh0_990[k]
                    - f_10 * mh1_990[k]
                    + pb_x[k] * mi_1319[k];

        t_1696[k] = f_17 * li_1066[k]
                    + pb_y[k] * mi_1318[k];

        t_1697[k] = f_9 * mh0_992[k]
                    - f_10 * mh1_992[k]
                    + pb_x[k] * mi_1321[k];
    }

#pragma omp simd aligned(t_1698, t_1699, t_1700, pb_x, pb_y, pb_z, li_1039, li_1069, mh0_993, \
                         mh1_993, mi_1319, mi_1321, mi_1322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1698[k] = f_7 * mh0_993[k]
                    - f_8 * mh1_993[k]
                    + pb_x[k] * mi_1322[k];

        t_1699[k] = f_12 * li_1039[k]
                    + pb_z[k] * mi_1319[k];

        t_1700[k] = f_17 * li_1069[k]
                    + pb_y[k] * mi_1321[k];
    }

#pragma omp simd aligned(t_1701, t_1702, t_1703, pb_x, pb_z, li_1042, mh0_996, mh0_997, \
                         mh1_996, mh1_997, mi_1322, mi_1325, mi_1326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1701[k] = f_7 * mh0_996[k]
                    - f_8 * mh1_996[k]
                    + pb_x[k] * mi_1325[k];

        t_1702[k] = f_5 * mh0_997[k]
                    - f_6 * mh1_997[k]
                    + pb_x[k] * mi_1326[k];

        t_1703[k] = f_12 * li_1042[k]
                    + pb_z[k] * mi_1322[k];
    }

#pragma omp simd aligned(t_1704, t_1705, t_1706, pb_x, pb_y, li_1073, mh0_999, mh0_1001, \
                         mh1_999, mh1_1001, mi_1325, mi_1328, mi_1330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1704[k] = f_5 * mh0_999[k]
                    - f_6 * mh1_999[k]
                    + pb_x[k] * mi_1328[k];

        t_1705[k] = f_17 * li_1073[k]
                    + pb_y[k] * mi_1325[k];

        t_1706[k] = f_5 * mh0_1001[k]
                    - f_6 * mh1_1001[k]
                    + pb_x[k] * mi_1330[k];
    }

#pragma omp simd aligned(t_1707, t_1708, t_1709, pb_x, pb_z, li_1046, mh0_1002, mh0_1004, \
                         mh1_1002, mh1_1004, mi_1326, mi_1331, \
                         mi_1333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1707[k] = f_3 * mh0_1002[k]
                    - f_4 * mh1_1002[k]
                    + pb_x[k] * mi_1331[k];

        t_1708[k] = f_12 * li_1046[k]
                    + pb_z[k] * mi_1326[k];

        t_1709[k] = f_3 * mh0_1004[k]
                    - f_4 * mh1_1004[k]
                    + pb_x[k] * mi_1333[k];
    }

#pragma omp simd aligned(t_1710, t_1711, t_1712, t_1713, pb_x, pb_y, li_1078, mh0_1005, \
                         mh0_1007, mh1_1005, mh1_1007, mi_1330, mi_1334, mi_1336, \
                         mi_1337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1710[k] = f_3 * mh0_1005[k]
                    - f_4 * mh1_1005[k]
                    + pb_x[k] * mi_1334[k];

        t_1711[k] = f_17 * li_1078[k]
                    + pb_y[k] * mi_1330[k];

        t_1712[k] = f_3 * mh0_1007[k]
                    - f_4 * mh1_1007[k]
                    + pb_x[k] * mi_1336[k];

        t_1713[k] = pb_x[k] * mi_1337[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece13(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t kk0,
                                             const size_t kk1, const size_t li, const size_t lk,
                                             const size_t mh0, const size_t mh1, const size_t mi,
                                             const size_t ncols, const double alpha,
                                             const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_17 = 3.5 / p;
    const auto f_18 = 0.5 / alpha;
    const auto f_19 = 0.5 * beta / (alpha * p);
    const auto f_20 = 3.0 / alpha;
    const auto f_21 = 3.0 * beta / (alpha * p);
    const auto f_22 = 1.0 / alpha;
    const auto f_23 = beta / (alpha * p);
    const auto f_24 = 3.0 / p;
    const auto f_25 = 2.5 / alpha;
    const auto f_26 = 2.5 * beta / (alpha * p);
    const auto f_27 = 1.5 / alpha;
    const auto f_28 = 1.5 * beta / (alpha * p);
    const auto f_29 = 2.0 / alpha;
    const auto f_30 = 2.0 * beta / (alpha * p);

    auto *t_1714 = buffer.data(target + 1714);
    auto *t_1715 = buffer.data(target + 1715);
    auto *t_1716 = buffer.data(target + 1716);
    auto *t_1717 = buffer.data(target + 1717);
    auto *t_1718 = buffer.data(target + 1718);
    auto *t_1719 = buffer.data(target + 1719);
    auto *t_1720 = buffer.data(target + 1720);
    auto *t_1721 = buffer.data(target + 1721);
    auto *t_1722 = buffer.data(target + 1722);
    auto *t_1723 = buffer.data(target + 1723);
    auto *t_1724 = buffer.data(target + 1724);
    auto *t_1725 = buffer.data(target + 1725);
    auto *t_1726 = buffer.data(target + 1726);
    auto *t_1727 = buffer.data(target + 1727);
    auto *t_1728 = buffer.data(target + 1728);
    auto *t_1729 = buffer.data(target + 1729);
    auto *t_1730 = buffer.data(target + 1730);
    auto *t_1731 = buffer.data(target + 1731);
    auto *t_1732 = buffer.data(target + 1732);
    auto *t_1733 = buffer.data(target + 1733);
    auto *t_1734 = buffer.data(target + 1734);
    auto *t_1735 = buffer.data(target + 1735);
    auto *t_1736 = buffer.data(target + 1736);
    auto *t_1737 = buffer.data(target + 1737);
    auto *t_1738 = buffer.data(target + 1738);
    auto *t_1739 = buffer.data(target + 1739);
    auto *t_1740 = buffer.data(target + 1740);
    auto *t_1741 = buffer.data(target + 1741);
    auto *t_1742 = buffer.data(target + 1742);
    auto *t_1743 = buffer.data(target + 1743);
    auto *t_1744 = buffer.data(target + 1744);
    auto *t_1745 = buffer.data(target + 1745);
    auto *t_1746 = buffer.data(target + 1746);
    auto *t_1747 = buffer.data(target + 1747);
    auto *t_1748 = buffer.data(target + 1748);
    auto *t_1749 = buffer.data(target + 1749);
    auto *t_1750 = buffer.data(target + 1750);
    auto *t_1751 = buffer.data(target + 1751);
    auto *t_1752 = buffer.data(target + 1752);
    auto *t_1753 = buffer.data(target + 1753);
    auto *t_1754 = buffer.data(target + 1754);
    auto *t_1755 = buffer.data(target + 1755);
    auto *t_1756 = buffer.data(target + 1756);
    auto *t_1757 = buffer.data(target + 1757);
    auto *t_1758 = buffer.data(target + 1758);
    auto *t_1759 = buffer.data(target + 1759);
    auto *t_1760 = buffer.data(target + 1760);
    auto *t_1761 = buffer.data(target + 1761);
    auto *t_1762 = buffer.data(target + 1762);
    auto *t_1763 = buffer.data(target + 1763);
    auto *t_1764 = buffer.data(target + 1764);
    auto *t_1765 = buffer.data(target + 1765);
    auto *t_1766 = buffer.data(target + 1766);
    auto *t_1767 = buffer.data(target + 1767);
    auto *t_1768 = buffer.data(target + 1768);
    auto *t_1769 = buffer.data(target + 1769);
    auto *t_1770 = buffer.data(target + 1770);
    auto *t_1771 = buffer.data(target + 1771);
    auto *t_1772 = buffer.data(target + 1772);
    auto *t_1773 = buffer.data(target + 1773);
    auto *t_1774 = buffer.data(target + 1774);
    auto *t_1775 = buffer.data(target + 1775);
    auto *t_1776 = buffer.data(target + 1776);
    auto *t_1777 = buffer.data(target + 1777);
    auto *t_1778 = buffer.data(target + 1778);
    auto *t_1779 = buffer.data(target + 1779);
    auto *t_1780 = buffer.data(target + 1780);
    auto *t_1781 = buffer.data(target + 1781);
    auto *t_1782 = buffer.data(target + 1782);
    auto *t_1783 = buffer.data(target + 1783);
    auto *t_1784 = buffer.data(target + 1784);
    auto *t_1785 = buffer.data(target + 1785);
    auto *t_1786 = buffer.data(target + 1786);
    auto *t_1787 = buffer.data(target + 1787);
    auto *t_1788 = buffer.data(target + 1788);
    auto *t_1789 = buffer.data(target + 1789);
    auto *t_1790 = buffer.data(target + 1790);
    auto *t_1791 = buffer.data(target + 1791);
    auto *t_1792 = buffer.data(target + 1792);
    auto *t_1793 = buffer.data(target + 1793);
    auto *t_1794 = buffer.data(target + 1794);
    auto *t_1795 = buffer.data(target + 1795);
    auto *t_1796 = buffer.data(target + 1796);
    auto *t_1797 = buffer.data(target + 1797);
    auto *t_1798 = buffer.data(target + 1798);
    auto *t_1799 = buffer.data(target + 1799);
    auto *t_1800 = buffer.data(target + 1800);
    auto *t_1801 = buffer.data(target + 1801);
    auto *t_1802 = buffer.data(target + 1802);
    auto *t_1803 = buffer.data(target + 1803);
    auto *t_1804 = buffer.data(target + 1804);
    auto *t_1805 = buffer.data(target + 1805);
    auto *t_1806 = buffer.data(target + 1806);
    auto *t_1807 = buffer.data(target + 1807);
    auto *t_1808 = buffer.data(target + 1808);
    auto *t_1809 = buffer.data(target + 1809);
    auto *t_1810 = buffer.data(target + 1810);
    auto *t_1811 = buffer.data(target + 1811);
    auto *t_1812 = buffer.data(target + 1812);
    auto *t_1813 = buffer.data(target + 1813);
    auto *t_1814 = buffer.data(target + 1814);
    auto *t_1815 = buffer.data(target + 1815);
    auto *t_1816 = buffer.data(target + 1816);
    auto *t_1817 = buffer.data(target + 1817);
    auto *t_1818 = buffer.data(target + 1818);
    auto *t_1819 = buffer.data(target + 1819);
    auto *t_1820 = buffer.data(target + 1820);
    auto *t_1821 = buffer.data(target + 1821);
    auto *t_1822 = buffer.data(target + 1822);
    auto *t_1823 = buffer.data(target + 1823);
    auto *t_1824 = buffer.data(target + 1824);
    auto *t_1825 = buffer.data(target + 1825);
    auto *t_1826 = buffer.data(target + 1826);
    auto *t_1827 = buffer.data(target + 1827);
    auto *t_1828 = buffer.data(target + 1828);
    auto *t_1829 = buffer.data(target + 1829);
    auto *t_1830 = buffer.data(target + 1830);
    auto *t_1831 = buffer.data(target + 1831);
    auto *t_1832 = buffer.data(target + 1832);
    auto *t_1833 = buffer.data(target + 1833);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kk0_1036 = buffer.data(kk0 + 1036);
    const auto *kk0_1072 = buffer.data(kk0 + 1072);
    const auto *kk0_1108 = buffer.data(kk0 + 1108);
    const auto *kk0_1115 = buffer.data(kk0 + 1115);
    const auto *kk0_1144 = buffer.data(kk0 + 1144);
    const auto *kk0_1151 = buffer.data(kk0 + 1151);
    const auto *kk0_1187 = buffer.data(kk0 + 1187);

    const auto *kk1_1036 = buffer.data(kk1 + 1036);
    const auto *kk1_1072 = buffer.data(kk1 + 1072);
    const auto *kk1_1108 = buffer.data(kk1 + 1108);
    const auto *kk1_1115 = buffer.data(kk1 + 1115);
    const auto *kk1_1144 = buffer.data(kk1 + 1144);
    const auto *kk1_1151 = buffer.data(kk1 + 1151);
    const auto *kk1_1187 = buffer.data(kk1 + 1187);

    const auto *li_1057 = buffer.data(li + 1057);
    const auto *li_1064 = buffer.data(li + 1064);
    const auto *li_1067 = buffer.data(li + 1067);
    const auto *li_1070 = buffer.data(li + 1070);
    const auto *li_1074 = buffer.data(li + 1074);
    const auto *li_1085 = buffer.data(li + 1085);
    const auto *li_1087 = buffer.data(li + 1087);
    const auto *li_1088 = buffer.data(li + 1088);
    const auto *li_1089 = buffer.data(li + 1089);
    const auto *li_1090 = buffer.data(li + 1090);
    const auto *li_1091 = buffer.data(li + 1091);
    const auto *li_1092 = buffer.data(li + 1092);
    const auto *li_1094 = buffer.data(li + 1094);
    const auto *li_1095 = buffer.data(li + 1095);
    const auto *li_1097 = buffer.data(li + 1097);
    const auto *li_1098 = buffer.data(li + 1098);
    const auto *li_1101 = buffer.data(li + 1101);
    const auto *li_1102 = buffer.data(li + 1102);
    const auto *li_1106 = buffer.data(li + 1106);
    const auto *li_1113 = buffer.data(li + 1113);
    const auto *li_1115 = buffer.data(li + 1115);
    const auto *li_1116 = buffer.data(li + 1116);
    const auto *li_1117 = buffer.data(li + 1117);
    const auto *li_1118 = buffer.data(li + 1118);
    const auto *li_1119 = buffer.data(li + 1119);
    const auto *li_1120 = buffer.data(li + 1120);
    const auto *li_1122 = buffer.data(li + 1122);
    const auto *li_1123 = buffer.data(li + 1123);
    const auto *li_1125 = buffer.data(li + 1125);
    const auto *li_1126 = buffer.data(li + 1126);
    const auto *li_1129 = buffer.data(li + 1129);
    const auto *li_1130 = buffer.data(li + 1130);
    const auto *li_1134 = buffer.data(li + 1134);
    const auto *li_1141 = buffer.data(li + 1141);
    const auto *li_1143 = buffer.data(li + 1143);
    const auto *li_1144 = buffer.data(li + 1144);
    const auto *li_1145 = buffer.data(li + 1145);
    const auto *li_1146 = buffer.data(li + 1146);
    const auto *li_1147 = buffer.data(li + 1147);
    const auto *li_1148 = buffer.data(li + 1148);
    const auto *li_1150 = buffer.data(li + 1150);
    const auto *li_1153 = buffer.data(li + 1153);
    const auto *li_1157 = buffer.data(li + 1157);
    const auto *li_1162 = buffer.data(li + 1162);
    const auto *li_1171 = buffer.data(li + 1171);
    const auto *li_1172 = buffer.data(li + 1172);
    const auto *li_1173 = buffer.data(li + 1173);
    const auto *li_1174 = buffer.data(li + 1174);

    const auto *lk_1360 = buffer.data(lk + 1360);
    const auto *lk_1396 = buffer.data(lk + 1396);
    const auto *lk_1403 = buffer.data(lk + 1403);
    const auto *lk_1432 = buffer.data(lk + 1432);
    const auto *lk_1439 = buffer.data(lk + 1439);
    const auto *lk_1468 = buffer.data(lk + 1468);
    const auto *lk_1475 = buffer.data(lk + 1475);

    const auto *mh0_1004 = buffer.data(mh0 + 1004);
    const auto *mh0_1005 = buffer.data(mh0 + 1005);
    const auto *mh0_1006 = buffer.data(mh0 + 1006);
    const auto *mh0_1007 = buffer.data(mh0 + 1007);
    const auto *mh0_1008 = buffer.data(mh0 + 1008);
    const auto *mh0_1011 = buffer.data(mh0 + 1011);
    const auto *mh0_1013 = buffer.data(mh0 + 1013);
    const auto *mh0_1014 = buffer.data(mh0 + 1014);
    const auto *mh0_1017 = buffer.data(mh0 + 1017);
    const auto *mh0_1018 = buffer.data(mh0 + 1018);
    const auto *mh0_1020 = buffer.data(mh0 + 1020);
    const auto *mh0_1022 = buffer.data(mh0 + 1022);
    const auto *mh0_1023 = buffer.data(mh0 + 1023);
    const auto *mh0_1025 = buffer.data(mh0 + 1025);
    const auto *mh0_1026 = buffer.data(mh0 + 1026);
    const auto *mh0_1027 = buffer.data(mh0 + 1027);
    const auto *mh0_1028 = buffer.data(mh0 + 1028);
    const auto *mh0_1029 = buffer.data(mh0 + 1029);
    const auto *mh0_1032 = buffer.data(mh0 + 1032);
    const auto *mh0_1034 = buffer.data(mh0 + 1034);
    const auto *mh0_1035 = buffer.data(mh0 + 1035);
    const auto *mh0_1038 = buffer.data(mh0 + 1038);
    const auto *mh0_1039 = buffer.data(mh0 + 1039);
    const auto *mh0_1041 = buffer.data(mh0 + 1041);
    const auto *mh0_1043 = buffer.data(mh0 + 1043);
    const auto *mh0_1044 = buffer.data(mh0 + 1044);
    const auto *mh0_1046 = buffer.data(mh0 + 1046);
    const auto *mh0_1047 = buffer.data(mh0 + 1047);
    const auto *mh0_1048 = buffer.data(mh0 + 1048);
    const auto *mh0_1049 = buffer.data(mh0 + 1049);
    const auto *mh0_1050 = buffer.data(mh0 + 1050);
    const auto *mh0_1053 = buffer.data(mh0 + 1053);
    const auto *mh0_1055 = buffer.data(mh0 + 1055);
    const auto *mh0_1056 = buffer.data(mh0 + 1056);
    const auto *mh0_1059 = buffer.data(mh0 + 1059);
    const auto *mh0_1060 = buffer.data(mh0 + 1060);
    const auto *mh0_1062 = buffer.data(mh0 + 1062);
    const auto *mh0_1064 = buffer.data(mh0 + 1064);
    const auto *mh0_1065 = buffer.data(mh0 + 1065);
    const auto *mh0_1067 = buffer.data(mh0 + 1067);
    const auto *mh0_1068 = buffer.data(mh0 + 1068);
    const auto *mh0_1069 = buffer.data(mh0 + 1069);
    const auto *mh0_1070 = buffer.data(mh0 + 1070);

    const auto *mh1_1004 = buffer.data(mh1 + 1004);
    const auto *mh1_1005 = buffer.data(mh1 + 1005);
    const auto *mh1_1006 = buffer.data(mh1 + 1006);
    const auto *mh1_1007 = buffer.data(mh1 + 1007);
    const auto *mh1_1008 = buffer.data(mh1 + 1008);
    const auto *mh1_1011 = buffer.data(mh1 + 1011);
    const auto *mh1_1013 = buffer.data(mh1 + 1013);
    const auto *mh1_1014 = buffer.data(mh1 + 1014);
    const auto *mh1_1017 = buffer.data(mh1 + 1017);
    const auto *mh1_1018 = buffer.data(mh1 + 1018);
    const auto *mh1_1020 = buffer.data(mh1 + 1020);
    const auto *mh1_1022 = buffer.data(mh1 + 1022);
    const auto *mh1_1023 = buffer.data(mh1 + 1023);
    const auto *mh1_1025 = buffer.data(mh1 + 1025);
    const auto *mh1_1026 = buffer.data(mh1 + 1026);
    const auto *mh1_1027 = buffer.data(mh1 + 1027);
    const auto *mh1_1028 = buffer.data(mh1 + 1028);
    const auto *mh1_1029 = buffer.data(mh1 + 1029);
    const auto *mh1_1032 = buffer.data(mh1 + 1032);
    const auto *mh1_1034 = buffer.data(mh1 + 1034);
    const auto *mh1_1035 = buffer.data(mh1 + 1035);
    const auto *mh1_1038 = buffer.data(mh1 + 1038);
    const auto *mh1_1039 = buffer.data(mh1 + 1039);
    const auto *mh1_1041 = buffer.data(mh1 + 1041);
    const auto *mh1_1043 = buffer.data(mh1 + 1043);
    const auto *mh1_1044 = buffer.data(mh1 + 1044);
    const auto *mh1_1046 = buffer.data(mh1 + 1046);
    const auto *mh1_1047 = buffer.data(mh1 + 1047);
    const auto *mh1_1048 = buffer.data(mh1 + 1048);
    const auto *mh1_1049 = buffer.data(mh1 + 1049);
    const auto *mh1_1050 = buffer.data(mh1 + 1050);
    const auto *mh1_1053 = buffer.data(mh1 + 1053);
    const auto *mh1_1055 = buffer.data(mh1 + 1055);
    const auto *mh1_1056 = buffer.data(mh1 + 1056);
    const auto *mh1_1059 = buffer.data(mh1 + 1059);
    const auto *mh1_1060 = buffer.data(mh1 + 1060);
    const auto *mh1_1062 = buffer.data(mh1 + 1062);
    const auto *mh1_1064 = buffer.data(mh1 + 1064);
    const auto *mh1_1065 = buffer.data(mh1 + 1065);
    const auto *mh1_1067 = buffer.data(mh1 + 1067);
    const auto *mh1_1068 = buffer.data(mh1 + 1068);
    const auto *mh1_1069 = buffer.data(mh1 + 1069);
    const auto *mh1_1070 = buffer.data(mh1 + 1070);

    const auto *mi_1337 = buffer.data(mi + 1337);
    const auto *mi_1338 = buffer.data(mi + 1338);
    const auto *mi_1339 = buffer.data(mi + 1339);
    const auto *mi_1340 = buffer.data(mi + 1340);
    const auto *mi_1341 = buffer.data(mi + 1341);
    const auto *mi_1342 = buffer.data(mi + 1342);
    const auto *mi_1343 = buffer.data(mi + 1343);
    const auto *mi_1344 = buffer.data(mi + 1344);
    const auto *mi_1346 = buffer.data(mi + 1346);
    const auto *mi_1347 = buffer.data(mi + 1347);
    const auto *mi_1349 = buffer.data(mi + 1349);
    const auto *mi_1350 = buffer.data(mi + 1350);
    const auto *mi_1353 = buffer.data(mi + 1353);
    const auto *mi_1354 = buffer.data(mi + 1354);
    const auto *mi_1356 = buffer.data(mi + 1356);
    const auto *mi_1358 = buffer.data(mi + 1358);
    const auto *mi_1359 = buffer.data(mi + 1359);
    const auto *mi_1361 = buffer.data(mi + 1361);
    const auto *mi_1362 = buffer.data(mi + 1362);
    const auto *mi_1364 = buffer.data(mi + 1364);
    const auto *mi_1365 = buffer.data(mi + 1365);
    const auto *mi_1366 = buffer.data(mi + 1366);
    const auto *mi_1367 = buffer.data(mi + 1367);
    const auto *mi_1368 = buffer.data(mi + 1368);
    const auto *mi_1369 = buffer.data(mi + 1369);
    const auto *mi_1370 = buffer.data(mi + 1370);
    const auto *mi_1371 = buffer.data(mi + 1371);
    const auto *mi_1372 = buffer.data(mi + 1372);
    const auto *mi_1374 = buffer.data(mi + 1374);
    const auto *mi_1375 = buffer.data(mi + 1375);
    const auto *mi_1377 = buffer.data(mi + 1377);
    const auto *mi_1378 = buffer.data(mi + 1378);
    const auto *mi_1381 = buffer.data(mi + 1381);
    const auto *mi_1382 = buffer.data(mi + 1382);
    const auto *mi_1384 = buffer.data(mi + 1384);
    const auto *mi_1386 = buffer.data(mi + 1386);
    const auto *mi_1387 = buffer.data(mi + 1387);
    const auto *mi_1389 = buffer.data(mi + 1389);
    const auto *mi_1390 = buffer.data(mi + 1390);
    const auto *mi_1392 = buffer.data(mi + 1392);
    const auto *mi_1393 = buffer.data(mi + 1393);
    const auto *mi_1394 = buffer.data(mi + 1394);
    const auto *mi_1395 = buffer.data(mi + 1395);
    const auto *mi_1396 = buffer.data(mi + 1396);
    const auto *mi_1397 = buffer.data(mi + 1397);
    const auto *mi_1398 = buffer.data(mi + 1398);
    const auto *mi_1399 = buffer.data(mi + 1399);
    const auto *mi_1400 = buffer.data(mi + 1400);
    const auto *mi_1402 = buffer.data(mi + 1402);
    const auto *mi_1403 = buffer.data(mi + 1403);
    const auto *mi_1405 = buffer.data(mi + 1405);
    const auto *mi_1406 = buffer.data(mi + 1406);
    const auto *mi_1409 = buffer.data(mi + 1409);
    const auto *mi_1410 = buffer.data(mi + 1410);
    const auto *mi_1412 = buffer.data(mi + 1412);
    const auto *mi_1414 = buffer.data(mi + 1414);
    const auto *mi_1415 = buffer.data(mi + 1415);
    const auto *mi_1417 = buffer.data(mi + 1417);
    const auto *mi_1418 = buffer.data(mi + 1418);
    const auto *mi_1420 = buffer.data(mi + 1420);
    const auto *mi_1421 = buffer.data(mi + 1421);
    const auto *mi_1422 = buffer.data(mi + 1422);
    const auto *mi_1423 = buffer.data(mi + 1423);
    const auto *mi_1424 = buffer.data(mi + 1424);
    const auto *mi_1425 = buffer.data(mi + 1425);
    const auto *mi_1426 = buffer.data(mi + 1426);
    const auto *mi_1427 = buffer.data(mi + 1427);

#pragma omp simd aligned(t_1714, t_1715, t_1716, t_1717, t_1718, t_1719, pb_x, mi_1338, \
                         mi_1339, mi_1340, mi_1341, mi_1342, mi_1343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1714[k] = pb_x[k] * mi_1338[k];

        t_1715[k] = pb_x[k] * mi_1339[k];

        t_1716[k] = pb_x[k] * mi_1340[k];

        t_1717[k] = pb_x[k] * mi_1341[k];

        t_1718[k] = pb_x[k] * mi_1342[k];

        t_1719[k] = pb_x[k] * mi_1343[k];
    }

#pragma omp simd aligned(t_1720, t_1721, t_1722, pa_z, pb_y, pb_z, kk0_1036, kk1_1036, \
                         li_1057, li_1087, lk_1360, mh0_1004, mh1_1004, mi_1337, \
                         mi_1339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1720[k] = f_18 * kk0_1036[k]
                    - f_19 * kk1_1036[k]
                    + pa_z[k] * lk_1360[k];

        t_1721[k] = f_12 * li_1057[k]
                    + pb_z[k] * mi_1337[k];

        t_1722[k] = f_17 * li_1087[k]
                    + f_9 * mh0_1004[k]
                    - f_10 * mh1_1004[k]
                    + pb_y[k] * mi_1339[k];
    }

#pragma omp simd aligned(t_1723, t_1724, t_1725, pb_y, li_1088, li_1089, li_1090, mh0_1005, \
                         mh0_1006, mh0_1007, mh1_1005, mh1_1006, mh1_1007, mi_1340, mi_1341, \
                         mi_1342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1723[k] = f_17 * li_1088[k]
                    + f_7 * mh0_1005[k]
                    - f_8 * mh1_1005[k]
                    + pb_y[k] * mi_1340[k];

        t_1724[k] = f_17 * li_1089[k]
                    + f_5 * mh0_1006[k]
                    - f_6 * mh1_1006[k]
                    + pb_y[k] * mi_1341[k];

        t_1725[k] = f_17 * li_1090[k]
                    + f_3 * mh0_1007[k]
                    - f_4 * mh1_1007[k]
                    + pb_y[k] * mi_1342[k];
    }

#pragma omp simd aligned(t_1726, t_1727, t_1728, t_1729, pa_y, pb_x, pb_y, kk0_1115, kk1_1115, \
                         li_1091, li_1092, lk_1403, mh0_1008, mh1_1008, mi_1343, \
                         mi_1344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1726[k] = f_17 * li_1091[k]
                    + pb_y[k] * mi_1343[k];

        t_1727[k] = f_20 * kk0_1115[k]
                    - f_21 * kk1_1115[k]
                    + pa_y[k] * lk_1403[k];

        t_1728[k] = f_1 * mh0_1008[k]
                    - f_2 * mh1_1008[k]
                    + pb_x[k] * mi_1344[k];

        t_1729[k] = f_24 * li_1092[k]
                    + pb_y[k] * mi_1344[k];
    }

#pragma omp simd aligned(t_1730, t_1731, t_1732, pb_x, pb_y, pb_z, li_1064, li_1094, mh0_1011, \
                         mh1_1011, mi_1344, mi_1346, mi_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1730[k] = f_13 * li_1064[k]
                    + pb_z[k] * mi_1344[k];

        t_1731[k] = f_9 * mh0_1011[k]
                    - f_10 * mh1_1011[k]
                    + pb_x[k] * mi_1347[k];

        t_1732[k] = f_24 * li_1094[k]
                    + pb_y[k] * mi_1346[k];
    }

#pragma omp simd aligned(t_1733, t_1734, t_1735, t_1736, pb_x, pb_y, pb_z, li_1067, li_1097, \
                         mh0_1013, mh0_1014, mh1_1013, mh1_1014, mi_1347, mi_1349, \
                         mi_1350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1733[k] = f_9 * mh0_1013[k]
                    - f_10 * mh1_1013[k]
                    + pb_x[k] * mi_1349[k];

        t_1734[k] = f_7 * mh0_1014[k]
                    - f_8 * mh1_1014[k]
                    + pb_x[k] * mi_1350[k];

        t_1735[k] = f_13 * li_1067[k]
                    + pb_z[k] * mi_1347[k];

        t_1736[k] = f_24 * li_1097[k]
                    + pb_y[k] * mi_1349[k];
    }

#pragma omp simd aligned(t_1737, t_1738, t_1739, pb_x, pb_z, li_1070, mh0_1017, mh0_1018, \
                         mh1_1017, mh1_1018, mi_1350, mi_1353, \
                         mi_1354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1737[k] = f_7 * mh0_1017[k]
                    - f_8 * mh1_1017[k]
                    + pb_x[k] * mi_1353[k];

        t_1738[k] = f_5 * mh0_1018[k]
                    - f_6 * mh1_1018[k]
                    + pb_x[k] * mi_1354[k];

        t_1739[k] = f_13 * li_1070[k]
                    + pb_z[k] * mi_1350[k];
    }

#pragma omp simd aligned(t_1740, t_1741, t_1742, pb_x, pb_y, li_1101, mh0_1020, mh0_1022, \
                         mh1_1020, mh1_1022, mi_1353, mi_1356, \
                         mi_1358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1740[k] = f_5 * mh0_1020[k]
                    - f_6 * mh1_1020[k]
                    + pb_x[k] * mi_1356[k];

        t_1741[k] = f_24 * li_1101[k]
                    + pb_y[k] * mi_1353[k];

        t_1742[k] = f_5 * mh0_1022[k]
                    - f_6 * mh1_1022[k]
                    + pb_x[k] * mi_1358[k];
    }

#pragma omp simd aligned(t_1743, t_1744, t_1745, pb_x, pb_z, li_1074, mh0_1023, mh0_1025, \
                         mh1_1023, mh1_1025, mi_1354, mi_1359, \
                         mi_1361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1743[k] = f_3 * mh0_1023[k]
                    - f_4 * mh1_1023[k]
                    + pb_x[k] * mi_1359[k];

        t_1744[k] = f_13 * li_1074[k]
                    + pb_z[k] * mi_1354[k];

        t_1745[k] = f_3 * mh0_1025[k]
                    - f_4 * mh1_1025[k]
                    + pb_x[k] * mi_1361[k];
    }

#pragma omp simd aligned(t_1746, t_1747, t_1748, t_1749, pb_x, pb_y, li_1106, mh0_1026, \
                         mh0_1028, mh1_1026, mh1_1028, mi_1358, mi_1362, mi_1364, \
                         mi_1365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1746[k] = f_3 * mh0_1026[k]
                    - f_4 * mh1_1026[k]
                    + pb_x[k] * mi_1362[k];

        t_1747[k] = f_24 * li_1106[k]
                    + pb_y[k] * mi_1358[k];

        t_1748[k] = f_3 * mh0_1028[k]
                    - f_4 * mh1_1028[k]
                    + pb_x[k] * mi_1364[k];

        t_1749[k] = pb_x[k] * mi_1365[k];
    }

#pragma omp simd aligned(t_1750, t_1751, t_1752, t_1753, t_1754, t_1755, pb_x, mi_1366, \
                         mi_1367, mi_1368, mi_1369, mi_1370, mi_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1750[k] = pb_x[k] * mi_1366[k];

        t_1751[k] = pb_x[k] * mi_1367[k];

        t_1752[k] = pb_x[k] * mi_1368[k];

        t_1753[k] = pb_x[k] * mi_1369[k];

        t_1754[k] = pb_x[k] * mi_1370[k];

        t_1755[k] = pb_x[k] * mi_1371[k];
    }

#pragma omp simd aligned(t_1756, t_1757, t_1758, pa_z, pb_y, pb_z, kk0_1072, kk1_1072, \
                         li_1085, li_1115, lk_1396, mh0_1025, mh1_1025, mi_1365, \
                         mi_1367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1756[k] = f_22 * kk0_1072[k]
                    - f_23 * kk1_1072[k]
                    + pa_z[k] * lk_1396[k];

        t_1757[k] = f_13 * li_1085[k]
                    + pb_z[k] * mi_1365[k];

        t_1758[k] = f_24 * li_1115[k]
                    + f_9 * mh0_1025[k]
                    - f_10 * mh1_1025[k]
                    + pb_y[k] * mi_1367[k];
    }

#pragma omp simd aligned(t_1759, t_1760, t_1761, pb_y, li_1116, li_1117, li_1118, mh0_1026, \
                         mh0_1027, mh0_1028, mh1_1026, mh1_1027, mh1_1028, mi_1368, mi_1369, \
                         mi_1370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1759[k] = f_24 * li_1116[k]
                    + f_7 * mh0_1026[k]
                    - f_8 * mh1_1026[k]
                    + pb_y[k] * mi_1368[k];

        t_1760[k] = f_24 * li_1117[k]
                    + f_5 * mh0_1027[k]
                    - f_6 * mh1_1027[k]
                    + pb_y[k] * mi_1369[k];

        t_1761[k] = f_24 * li_1118[k]
                    + f_3 * mh0_1028[k]
                    - f_4 * mh1_1028[k]
                    + pb_y[k] * mi_1370[k];
    }

#pragma omp simd aligned(t_1762, t_1763, t_1764, t_1765, pa_y, pb_x, pb_y, kk0_1151, kk1_1151, \
                         li_1119, li_1120, lk_1439, mh0_1029, mh1_1029, mi_1371, \
                         mi_1372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1762[k] = f_24 * li_1119[k]
                    + pb_y[k] * mi_1371[k];

        t_1763[k] = f_25 * kk0_1151[k]
                    - f_26 * kk1_1151[k]
                    + pa_y[k] * lk_1439[k];

        t_1764[k] = f_1 * mh0_1029[k]
                    - f_2 * mh1_1029[k]
                    + pb_x[k] * mi_1372[k];

        t_1765[k] = f_15 * li_1120[k]
                    + pb_y[k] * mi_1372[k];
    }

#pragma omp simd aligned(t_1766, t_1767, t_1768, pb_x, pb_y, pb_z, li_1092, li_1122, mh0_1032, \
                         mh1_1032, mi_1372, mi_1374, mi_1375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1766[k] = f_14 * li_1092[k]
                    + pb_z[k] * mi_1372[k];

        t_1767[k] = f_9 * mh0_1032[k]
                    - f_10 * mh1_1032[k]
                    + pb_x[k] * mi_1375[k];

        t_1768[k] = f_15 * li_1122[k]
                    + pb_y[k] * mi_1374[k];
    }

#pragma omp simd aligned(t_1769, t_1770, t_1771, t_1772, pb_x, pb_y, pb_z, li_1095, li_1125, \
                         mh0_1034, mh0_1035, mh1_1034, mh1_1035, mi_1375, mi_1377, \
                         mi_1378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1769[k] = f_9 * mh0_1034[k]
                    - f_10 * mh1_1034[k]
                    + pb_x[k] * mi_1377[k];

        t_1770[k] = f_7 * mh0_1035[k]
                    - f_8 * mh1_1035[k]
                    + pb_x[k] * mi_1378[k];

        t_1771[k] = f_14 * li_1095[k]
                    + pb_z[k] * mi_1375[k];

        t_1772[k] = f_15 * li_1125[k]
                    + pb_y[k] * mi_1377[k];
    }

#pragma omp simd aligned(t_1773, t_1774, t_1775, pb_x, pb_z, li_1098, mh0_1038, mh0_1039, \
                         mh1_1038, mh1_1039, mi_1378, mi_1381, \
                         mi_1382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1773[k] = f_7 * mh0_1038[k]
                    - f_8 * mh1_1038[k]
                    + pb_x[k] * mi_1381[k];

        t_1774[k] = f_5 * mh0_1039[k]
                    - f_6 * mh1_1039[k]
                    + pb_x[k] * mi_1382[k];

        t_1775[k] = f_14 * li_1098[k]
                    + pb_z[k] * mi_1378[k];
    }

#pragma omp simd aligned(t_1776, t_1777, t_1778, pb_x, pb_y, li_1129, mh0_1041, mh0_1043, \
                         mh1_1041, mh1_1043, mi_1381, mi_1384, \
                         mi_1386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1776[k] = f_5 * mh0_1041[k]
                    - f_6 * mh1_1041[k]
                    + pb_x[k] * mi_1384[k];

        t_1777[k] = f_15 * li_1129[k]
                    + pb_y[k] * mi_1381[k];

        t_1778[k] = f_5 * mh0_1043[k]
                    - f_6 * mh1_1043[k]
                    + pb_x[k] * mi_1386[k];
    }

#pragma omp simd aligned(t_1779, t_1780, t_1781, pb_x, pb_z, li_1102, mh0_1044, mh0_1046, \
                         mh1_1044, mh1_1046, mi_1382, mi_1387, \
                         mi_1389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1779[k] = f_3 * mh0_1044[k]
                    - f_4 * mh1_1044[k]
                    + pb_x[k] * mi_1387[k];

        t_1780[k] = f_14 * li_1102[k]
                    + pb_z[k] * mi_1382[k];

        t_1781[k] = f_3 * mh0_1046[k]
                    - f_4 * mh1_1046[k]
                    + pb_x[k] * mi_1389[k];
    }

#pragma omp simd aligned(t_1782, t_1783, t_1784, t_1785, pb_x, pb_y, li_1134, mh0_1047, \
                         mh0_1049, mh1_1047, mh1_1049, mi_1386, mi_1390, mi_1392, \
                         mi_1393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1782[k] = f_3 * mh0_1047[k]
                    - f_4 * mh1_1047[k]
                    + pb_x[k] * mi_1390[k];

        t_1783[k] = f_15 * li_1134[k]
                    + pb_y[k] * mi_1386[k];

        t_1784[k] = f_3 * mh0_1049[k]
                    - f_4 * mh1_1049[k]
                    + pb_x[k] * mi_1392[k];

        t_1785[k] = pb_x[k] * mi_1393[k];
    }

#pragma omp simd aligned(t_1786, t_1787, t_1788, t_1789, t_1790, t_1791, pb_x, mi_1394, \
                         mi_1395, mi_1396, mi_1397, mi_1398, mi_1399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1786[k] = pb_x[k] * mi_1394[k];

        t_1787[k] = pb_x[k] * mi_1395[k];

        t_1788[k] = pb_x[k] * mi_1396[k];

        t_1789[k] = pb_x[k] * mi_1397[k];

        t_1790[k] = pb_x[k] * mi_1398[k];

        t_1791[k] = pb_x[k] * mi_1399[k];
    }

#pragma omp simd aligned(t_1792, t_1793, t_1794, pa_z, pb_y, pb_z, kk0_1108, kk1_1108, \
                         li_1113, li_1143, lk_1432, mh0_1046, mh1_1046, mi_1393, \
                         mi_1395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1792[k] = f_27 * kk0_1108[k]
                    - f_28 * kk1_1108[k]
                    + pa_z[k] * lk_1432[k];

        t_1793[k] = f_14 * li_1113[k]
                    + pb_z[k] * mi_1393[k];

        t_1794[k] = f_15 * li_1143[k]
                    + f_9 * mh0_1046[k]
                    - f_10 * mh1_1046[k]
                    + pb_y[k] * mi_1395[k];
    }

#pragma omp simd aligned(t_1795, t_1796, t_1797, pb_y, li_1144, li_1145, li_1146, mh0_1047, \
                         mh0_1048, mh0_1049, mh1_1047, mh1_1048, mh1_1049, mi_1396, mi_1397, \
                         mi_1398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1795[k] = f_15 * li_1144[k]
                    + f_7 * mh0_1047[k]
                    - f_8 * mh1_1047[k]
                    + pb_y[k] * mi_1396[k];

        t_1796[k] = f_15 * li_1145[k]
                    + f_5 * mh0_1048[k]
                    - f_6 * mh1_1048[k]
                    + pb_y[k] * mi_1397[k];

        t_1797[k] = f_15 * li_1146[k]
                    + f_3 * mh0_1049[k]
                    - f_4 * mh1_1049[k]
                    + pb_y[k] * mi_1398[k];
    }

#pragma omp simd aligned(t_1798, t_1799, t_1800, t_1801, pa_y, pb_x, pb_y, kk0_1187, kk1_1187, \
                         li_1147, li_1148, lk_1475, mh0_1050, mh1_1050, mi_1399, \
                         mi_1400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1798[k] = f_15 * li_1147[k]
                    + pb_y[k] * mi_1399[k];

        t_1799[k] = f_29 * kk0_1187[k]
                    - f_30 * kk1_1187[k]
                    + pa_y[k] * lk_1475[k];

        t_1800[k] = f_1 * mh0_1050[k]
                    - f_2 * mh1_1050[k]
                    + pb_x[k] * mi_1400[k];

        t_1801[k] = f_14 * li_1148[k]
                    + pb_y[k] * mi_1400[k];
    }

#pragma omp simd aligned(t_1802, t_1803, t_1804, pb_x, pb_y, pb_z, li_1120, li_1150, mh0_1053, \
                         mh1_1053, mi_1400, mi_1402, mi_1403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1802[k] = f_15 * li_1120[k]
                    + pb_z[k] * mi_1400[k];

        t_1803[k] = f_9 * mh0_1053[k]
                    - f_10 * mh1_1053[k]
                    + pb_x[k] * mi_1403[k];

        t_1804[k] = f_14 * li_1150[k]
                    + pb_y[k] * mi_1402[k];
    }

#pragma omp simd aligned(t_1805, t_1806, t_1807, t_1808, pb_x, pb_y, pb_z, li_1123, li_1153, \
                         mh0_1055, mh0_1056, mh1_1055, mh1_1056, mi_1403, mi_1405, \
                         mi_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1805[k] = f_9 * mh0_1055[k]
                    - f_10 * mh1_1055[k]
                    + pb_x[k] * mi_1405[k];

        t_1806[k] = f_7 * mh0_1056[k]
                    - f_8 * mh1_1056[k]
                    + pb_x[k] * mi_1406[k];

        t_1807[k] = f_15 * li_1123[k]
                    + pb_z[k] * mi_1403[k];

        t_1808[k] = f_14 * li_1153[k]
                    + pb_y[k] * mi_1405[k];
    }

#pragma omp simd aligned(t_1809, t_1810, t_1811, pb_x, pb_z, li_1126, mh0_1059, mh0_1060, \
                         mh1_1059, mh1_1060, mi_1406, mi_1409, \
                         mi_1410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1809[k] = f_7 * mh0_1059[k]
                    - f_8 * mh1_1059[k]
                    + pb_x[k] * mi_1409[k];

        t_1810[k] = f_5 * mh0_1060[k]
                    - f_6 * mh1_1060[k]
                    + pb_x[k] * mi_1410[k];

        t_1811[k] = f_15 * li_1126[k]
                    + pb_z[k] * mi_1406[k];
    }

#pragma omp simd aligned(t_1812, t_1813, t_1814, pb_x, pb_y, li_1157, mh0_1062, mh0_1064, \
                         mh1_1062, mh1_1064, mi_1409, mi_1412, \
                         mi_1414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1812[k] = f_5 * mh0_1062[k]
                    - f_6 * mh1_1062[k]
                    + pb_x[k] * mi_1412[k];

        t_1813[k] = f_14 * li_1157[k]
                    + pb_y[k] * mi_1409[k];

        t_1814[k] = f_5 * mh0_1064[k]
                    - f_6 * mh1_1064[k]
                    + pb_x[k] * mi_1414[k];
    }

#pragma omp simd aligned(t_1815, t_1816, t_1817, pb_x, pb_z, li_1130, mh0_1065, mh0_1067, \
                         mh1_1065, mh1_1067, mi_1410, mi_1415, \
                         mi_1417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1815[k] = f_3 * mh0_1065[k]
                    - f_4 * mh1_1065[k]
                    + pb_x[k] * mi_1415[k];

        t_1816[k] = f_15 * li_1130[k]
                    + pb_z[k] * mi_1410[k];

        t_1817[k] = f_3 * mh0_1067[k]
                    - f_4 * mh1_1067[k]
                    + pb_x[k] * mi_1417[k];
    }

#pragma omp simd aligned(t_1818, t_1819, t_1820, t_1821, pb_x, pb_y, li_1162, mh0_1068, \
                         mh0_1070, mh1_1068, mh1_1070, mi_1414, mi_1418, mi_1420, \
                         mi_1421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1818[k] = f_3 * mh0_1068[k]
                    - f_4 * mh1_1068[k]
                    + pb_x[k] * mi_1418[k];

        t_1819[k] = f_14 * li_1162[k]
                    + pb_y[k] * mi_1414[k];

        t_1820[k] = f_3 * mh0_1070[k]
                    - f_4 * mh1_1070[k]
                    + pb_x[k] * mi_1420[k];

        t_1821[k] = pb_x[k] * mi_1421[k];
    }

#pragma omp simd aligned(t_1822, t_1823, t_1824, t_1825, t_1826, t_1827, pb_x, mi_1422, \
                         mi_1423, mi_1424, mi_1425, mi_1426, mi_1427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1822[k] = pb_x[k] * mi_1422[k];

        t_1823[k] = pb_x[k] * mi_1423[k];

        t_1824[k] = pb_x[k] * mi_1424[k];

        t_1825[k] = pb_x[k] * mi_1425[k];

        t_1826[k] = pb_x[k] * mi_1426[k];

        t_1827[k] = pb_x[k] * mi_1427[k];
    }

#pragma omp simd aligned(t_1828, t_1829, t_1830, pa_z, pb_y, pb_z, kk0_1144, kk1_1144, \
                         li_1141, li_1171, lk_1468, mh0_1067, mh1_1067, mi_1421, \
                         mi_1423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1828[k] = f_29 * kk0_1144[k]
                    - f_30 * kk1_1144[k]
                    + pa_z[k] * lk_1468[k];

        t_1829[k] = f_15 * li_1141[k]
                    + pb_z[k] * mi_1421[k];

        t_1830[k] = f_14 * li_1171[k]
                    + f_9 * mh0_1067[k]
                    - f_10 * mh1_1067[k]
                    + pb_y[k] * mi_1423[k];
    }

#pragma omp simd aligned(t_1831, t_1832, t_1833, pb_y, li_1172, li_1173, li_1174, mh0_1068, \
                         mh0_1069, mh0_1070, mh1_1068, mh1_1069, mh1_1070, mi_1424, mi_1425, \
                         mi_1426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1831[k] = f_14 * li_1172[k]
                    + f_7 * mh0_1068[k]
                    - f_8 * mh1_1068[k]
                    + pb_y[k] * mi_1424[k];

        t_1832[k] = f_14 * li_1173[k]
                    + f_5 * mh0_1069[k]
                    - f_6 * mh1_1069[k]
                    + pb_y[k] * mi_1425[k];

        t_1833[k] = f_14 * li_1174[k]
                    + f_3 * mh0_1070[k]
                    - f_4 * mh1_1070[k]
                    + pb_y[k] * mi_1426[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece14(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t kk0,
                                             const size_t kk1, const size_t li, const size_t lk,
                                             const size_t mh0, const size_t mh1, const size_t mi,
                                             const size_t ncols, const double alpha,
                                             const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;
    const auto f_16 = 4.0 / p;
    const auto f_17 = 3.5 / p;
    const auto f_18 = 0.5 / alpha;
    const auto f_19 = 0.5 * beta / (alpha * p);
    const auto f_20 = 3.0 / alpha;
    const auto f_21 = 3.0 * beta / (alpha * p);
    const auto f_22 = 1.0 / alpha;
    const auto f_23 = beta / (alpha * p);
    const auto f_24 = 3.0 / p;
    const auto f_25 = 2.5 / alpha;
    const auto f_26 = 2.5 * beta / (alpha * p);
    const auto f_27 = 1.5 / alpha;
    const auto f_28 = 1.5 * beta / (alpha * p);

    auto *t_1834 = buffer.data(target + 1834);
    auto *t_1835 = buffer.data(target + 1835);
    auto *t_1836 = buffer.data(target + 1836);
    auto *t_1837 = buffer.data(target + 1837);
    auto *t_1838 = buffer.data(target + 1838);
    auto *t_1839 = buffer.data(target + 1839);
    auto *t_1840 = buffer.data(target + 1840);
    auto *t_1841 = buffer.data(target + 1841);
    auto *t_1842 = buffer.data(target + 1842);
    auto *t_1843 = buffer.data(target + 1843);
    auto *t_1844 = buffer.data(target + 1844);
    auto *t_1845 = buffer.data(target + 1845);
    auto *t_1846 = buffer.data(target + 1846);
    auto *t_1847 = buffer.data(target + 1847);
    auto *t_1848 = buffer.data(target + 1848);
    auto *t_1849 = buffer.data(target + 1849);
    auto *t_1850 = buffer.data(target + 1850);
    auto *t_1851 = buffer.data(target + 1851);
    auto *t_1852 = buffer.data(target + 1852);
    auto *t_1853 = buffer.data(target + 1853);
    auto *t_1854 = buffer.data(target + 1854);
    auto *t_1855 = buffer.data(target + 1855);
    auto *t_1856 = buffer.data(target + 1856);
    auto *t_1857 = buffer.data(target + 1857);
    auto *t_1858 = buffer.data(target + 1858);
    auto *t_1859 = buffer.data(target + 1859);
    auto *t_1860 = buffer.data(target + 1860);
    auto *t_1861 = buffer.data(target + 1861);
    auto *t_1862 = buffer.data(target + 1862);
    auto *t_1863 = buffer.data(target + 1863);
    auto *t_1864 = buffer.data(target + 1864);
    auto *t_1865 = buffer.data(target + 1865);
    auto *t_1866 = buffer.data(target + 1866);
    auto *t_1867 = buffer.data(target + 1867);
    auto *t_1868 = buffer.data(target + 1868);
    auto *t_1869 = buffer.data(target + 1869);
    auto *t_1870 = buffer.data(target + 1870);
    auto *t_1871 = buffer.data(target + 1871);
    auto *t_1872 = buffer.data(target + 1872);
    auto *t_1873 = buffer.data(target + 1873);
    auto *t_1874 = buffer.data(target + 1874);
    auto *t_1875 = buffer.data(target + 1875);
    auto *t_1876 = buffer.data(target + 1876);
    auto *t_1877 = buffer.data(target + 1877);
    auto *t_1878 = buffer.data(target + 1878);
    auto *t_1879 = buffer.data(target + 1879);
    auto *t_1880 = buffer.data(target + 1880);
    auto *t_1881 = buffer.data(target + 1881);
    auto *t_1882 = buffer.data(target + 1882);
    auto *t_1883 = buffer.data(target + 1883);
    auto *t_1884 = buffer.data(target + 1884);
    auto *t_1885 = buffer.data(target + 1885);
    auto *t_1886 = buffer.data(target + 1886);
    auto *t_1887 = buffer.data(target + 1887);
    auto *t_1888 = buffer.data(target + 1888);
    auto *t_1889 = buffer.data(target + 1889);
    auto *t_1890 = buffer.data(target + 1890);
    auto *t_1891 = buffer.data(target + 1891);
    auto *t_1892 = buffer.data(target + 1892);
    auto *t_1893 = buffer.data(target + 1893);
    auto *t_1894 = buffer.data(target + 1894);
    auto *t_1895 = buffer.data(target + 1895);
    auto *t_1896 = buffer.data(target + 1896);
    auto *t_1897 = buffer.data(target + 1897);
    auto *t_1898 = buffer.data(target + 1898);
    auto *t_1899 = buffer.data(target + 1899);
    auto *t_1900 = buffer.data(target + 1900);
    auto *t_1901 = buffer.data(target + 1901);
    auto *t_1902 = buffer.data(target + 1902);
    auto *t_1903 = buffer.data(target + 1903);
    auto *t_1904 = buffer.data(target + 1904);
    auto *t_1905 = buffer.data(target + 1905);
    auto *t_1906 = buffer.data(target + 1906);
    auto *t_1907 = buffer.data(target + 1907);
    auto *t_1908 = buffer.data(target + 1908);
    auto *t_1909 = buffer.data(target + 1909);
    auto *t_1910 = buffer.data(target + 1910);
    auto *t_1911 = buffer.data(target + 1911);
    auto *t_1912 = buffer.data(target + 1912);
    auto *t_1913 = buffer.data(target + 1913);
    auto *t_1914 = buffer.data(target + 1914);
    auto *t_1915 = buffer.data(target + 1915);
    auto *t_1916 = buffer.data(target + 1916);
    auto *t_1917 = buffer.data(target + 1917);
    auto *t_1918 = buffer.data(target + 1918);
    auto *t_1919 = buffer.data(target + 1919);
    auto *t_1920 = buffer.data(target + 1920);
    auto *t_1921 = buffer.data(target + 1921);
    auto *t_1922 = buffer.data(target + 1922);
    auto *t_1923 = buffer.data(target + 1923);
    auto *t_1924 = buffer.data(target + 1924);
    auto *t_1925 = buffer.data(target + 1925);
    auto *t_1926 = buffer.data(target + 1926);
    auto *t_1927 = buffer.data(target + 1927);
    auto *t_1928 = buffer.data(target + 1928);
    auto *t_1929 = buffer.data(target + 1929);
    auto *t_1930 = buffer.data(target + 1930);
    auto *t_1931 = buffer.data(target + 1931);
    auto *t_1932 = buffer.data(target + 1932);
    auto *t_1933 = buffer.data(target + 1933);
    auto *t_1934 = buffer.data(target + 1934);
    auto *t_1935 = buffer.data(target + 1935);
    auto *t_1936 = buffer.data(target + 1936);
    auto *t_1937 = buffer.data(target + 1937);
    auto *t_1938 = buffer.data(target + 1938);
    auto *t_1939 = buffer.data(target + 1939);
    auto *t_1940 = buffer.data(target + 1940);
    auto *t_1941 = buffer.data(target + 1941);
    auto *t_1942 = buffer.data(target + 1942);
    auto *t_1943 = buffer.data(target + 1943);
    auto *t_1944 = buffer.data(target + 1944);
    auto *t_1945 = buffer.data(target + 1945);
    auto *t_1946 = buffer.data(target + 1946);
    auto *t_1947 = buffer.data(target + 1947);
    auto *t_1948 = buffer.data(target + 1948);
    auto *t_1949 = buffer.data(target + 1949);
    auto *t_1950 = buffer.data(target + 1950);
    auto *t_1951 = buffer.data(target + 1951);
    auto *t_1952 = buffer.data(target + 1952);
    auto *t_1953 = buffer.data(target + 1953);
    auto *t_1954 = buffer.data(target + 1954);
    auto *t_1955 = buffer.data(target + 1955);
    auto *t_1956 = buffer.data(target + 1956);
    auto *t_1957 = buffer.data(target + 1957);
    auto *t_1958 = buffer.data(target + 1958);
    auto *t_1959 = buffer.data(target + 1959);
    auto *t_1960 = buffer.data(target + 1960);
    auto *t_1961 = buffer.data(target + 1961);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kk0_1180 = buffer.data(kk0 + 1180);
    const auto *kk0_1216 = buffer.data(kk0 + 1216);
    const auto *kk0_1223 = buffer.data(kk0 + 1223);
    const auto *kk0_1259 = buffer.data(kk0 + 1259);
    const auto *kk0_1295 = buffer.data(kk0 + 1295);

    const auto *kk1_1180 = buffer.data(kk1 + 1180);
    const auto *kk1_1216 = buffer.data(kk1 + 1216);
    const auto *kk1_1223 = buffer.data(kk1 + 1223);
    const auto *kk1_1259 = buffer.data(kk1 + 1259);
    const auto *kk1_1295 = buffer.data(kk1 + 1295);

    const auto *li_1148 = buffer.data(li + 1148);
    const auto *li_1151 = buffer.data(li + 1151);
    const auto *li_1154 = buffer.data(li + 1154);
    const auto *li_1158 = buffer.data(li + 1158);
    const auto *li_1169 = buffer.data(li + 1169);
    const auto *li_1175 = buffer.data(li + 1175);
    const auto *li_1176 = buffer.data(li + 1176);
    const auto *li_1178 = buffer.data(li + 1178);
    const auto *li_1179 = buffer.data(li + 1179);
    const auto *li_1181 = buffer.data(li + 1181);
    const auto *li_1182 = buffer.data(li + 1182);
    const auto *li_1185 = buffer.data(li + 1185);
    const auto *li_1186 = buffer.data(li + 1186);
    const auto *li_1190 = buffer.data(li + 1190);
    const auto *li_1197 = buffer.data(li + 1197);
    const auto *li_1199 = buffer.data(li + 1199);
    const auto *li_1200 = buffer.data(li + 1200);
    const auto *li_1201 = buffer.data(li + 1201);
    const auto *li_1202 = buffer.data(li + 1202);
    const auto *li_1203 = buffer.data(li + 1203);
    const auto *li_1204 = buffer.data(li + 1204);
    const auto *li_1206 = buffer.data(li + 1206);
    const auto *li_1207 = buffer.data(li + 1207);
    const auto *li_1209 = buffer.data(li + 1209);
    const auto *li_1210 = buffer.data(li + 1210);
    const auto *li_1213 = buffer.data(li + 1213);
    const auto *li_1214 = buffer.data(li + 1214);
    const auto *li_1218 = buffer.data(li + 1218);
    const auto *li_1225 = buffer.data(li + 1225);
    const auto *li_1227 = buffer.data(li + 1227);
    const auto *li_1228 = buffer.data(li + 1228);
    const auto *li_1229 = buffer.data(li + 1229);
    const auto *li_1230 = buffer.data(li + 1230);
    const auto *li_1231 = buffer.data(li + 1231);
    const auto *li_1232 = buffer.data(li + 1232);
    const auto *li_1233 = buffer.data(li + 1233);
    const auto *li_1234 = buffer.data(li + 1234);
    const auto *li_1235 = buffer.data(li + 1235);
    const auto *li_1237 = buffer.data(li + 1237);
    const auto *li_1238 = buffer.data(li + 1238);
    const auto *li_1240 = buffer.data(li + 1240);
    const auto *li_1241 = buffer.data(li + 1241);
    const auto *li_1242 = buffer.data(li + 1242);
    const auto *li_1244 = buffer.data(li + 1244);
    const auto *li_1245 = buffer.data(li + 1245);
    const auto *li_1246 = buffer.data(li + 1246);
    const auto *li_1253 = buffer.data(li + 1253);
    const auto *li_1255 = buffer.data(li + 1255);
    const auto *li_1256 = buffer.data(li + 1256);
    const auto *li_1257 = buffer.data(li + 1257);
    const auto *li_1258 = buffer.data(li + 1258);
    const auto *li_1259 = buffer.data(li + 1259);

    const auto *lk_1504 = buffer.data(lk + 1504);
    const auto *lk_1511 = buffer.data(lk + 1511);
    const auto *lk_1540 = buffer.data(lk + 1540);
    const auto *lk_1547 = buffer.data(lk + 1547);
    const auto *lk_1583 = buffer.data(lk + 1583);
    const auto *lk_1584 = buffer.data(lk + 1584);
    const auto *lk_1586 = buffer.data(lk + 1586);
    const auto *lk_1587 = buffer.data(lk + 1587);
    const auto *lk_1589 = buffer.data(lk + 1589);
    const auto *lk_1590 = buffer.data(lk + 1590);
    const auto *lk_1593 = buffer.data(lk + 1593);
    const auto *lk_1594 = buffer.data(lk + 1594);
    const auto *lk_1596 = buffer.data(lk + 1596);
    const auto *lk_1598 = buffer.data(lk + 1598);
    const auto *lk_1599 = buffer.data(lk + 1599);
    const auto *lk_1601 = buffer.data(lk + 1601);
    const auto *lk_1602 = buffer.data(lk + 1602);
    const auto *lk_1604 = buffer.data(lk + 1604);
    const auto *lk_1612 = buffer.data(lk + 1612);
    const auto *lk_1614 = buffer.data(lk + 1614);
    const auto *lk_1615 = buffer.data(lk + 1615);
    const auto *lk_1616 = buffer.data(lk + 1616);
    const auto *lk_1617 = buffer.data(lk + 1617);
    const auto *lk_1619 = buffer.data(lk + 1619);

    const auto *mh0_1071 = buffer.data(mh0 + 1071);
    const auto *mh0_1074 = buffer.data(mh0 + 1074);
    const auto *mh0_1076 = buffer.data(mh0 + 1076);
    const auto *mh0_1077 = buffer.data(mh0 + 1077);
    const auto *mh0_1080 = buffer.data(mh0 + 1080);
    const auto *mh0_1081 = buffer.data(mh0 + 1081);
    const auto *mh0_1083 = buffer.data(mh0 + 1083);
    const auto *mh0_1085 = buffer.data(mh0 + 1085);
    const auto *mh0_1086 = buffer.data(mh0 + 1086);
    const auto *mh0_1088 = buffer.data(mh0 + 1088);
    const auto *mh0_1089 = buffer.data(mh0 + 1089);
    const auto *mh0_1090 = buffer.data(mh0 + 1090);
    const auto *mh0_1091 = buffer.data(mh0 + 1091);
    const auto *mh0_1092 = buffer.data(mh0 + 1092);
    const auto *mh0_1095 = buffer.data(mh0 + 1095);
    const auto *mh0_1097 = buffer.data(mh0 + 1097);
    const auto *mh0_1098 = buffer.data(mh0 + 1098);
    const auto *mh0_1101 = buffer.data(mh0 + 1101);
    const auto *mh0_1102 = buffer.data(mh0 + 1102);
    const auto *mh0_1104 = buffer.data(mh0 + 1104);
    const auto *mh0_1106 = buffer.data(mh0 + 1106);
    const auto *mh0_1107 = buffer.data(mh0 + 1107);
    const auto *mh0_1109 = buffer.data(mh0 + 1109);
    const auto *mh0_1110 = buffer.data(mh0 + 1110);
    const auto *mh0_1111 = buffer.data(mh0 + 1111);
    const auto *mh0_1112 = buffer.data(mh0 + 1112);
    const auto *mh0_1134 = buffer.data(mh0 + 1134);
    const auto *mh0_1137 = buffer.data(mh0 + 1137);
    const auto *mh0_1139 = buffer.data(mh0 + 1139);
    const auto *mh0_1140 = buffer.data(mh0 + 1140);
    const auto *mh0_1143 = buffer.data(mh0 + 1143);
    const auto *mh0_1144 = buffer.data(mh0 + 1144);
    const auto *mh0_1146 = buffer.data(mh0 + 1146);
    const auto *mh0_1148 = buffer.data(mh0 + 1148);
    const auto *mh0_1149 = buffer.data(mh0 + 1149);
    const auto *mh0_1151 = buffer.data(mh0 + 1151);

    const auto *mh1_1071 = buffer.data(mh1 + 1071);
    const auto *mh1_1074 = buffer.data(mh1 + 1074);
    const auto *mh1_1076 = buffer.data(mh1 + 1076);
    const auto *mh1_1077 = buffer.data(mh1 + 1077);
    const auto *mh1_1080 = buffer.data(mh1 + 1080);
    const auto *mh1_1081 = buffer.data(mh1 + 1081);
    const auto *mh1_1083 = buffer.data(mh1 + 1083);
    const auto *mh1_1085 = buffer.data(mh1 + 1085);
    const auto *mh1_1086 = buffer.data(mh1 + 1086);
    const auto *mh1_1088 = buffer.data(mh1 + 1088);
    const auto *mh1_1089 = buffer.data(mh1 + 1089);
    const auto *mh1_1090 = buffer.data(mh1 + 1090);
    const auto *mh1_1091 = buffer.data(mh1 + 1091);
    const auto *mh1_1092 = buffer.data(mh1 + 1092);
    const auto *mh1_1095 = buffer.data(mh1 + 1095);
    const auto *mh1_1097 = buffer.data(mh1 + 1097);
    const auto *mh1_1098 = buffer.data(mh1 + 1098);
    const auto *mh1_1101 = buffer.data(mh1 + 1101);
    const auto *mh1_1102 = buffer.data(mh1 + 1102);
    const auto *mh1_1104 = buffer.data(mh1 + 1104);
    const auto *mh1_1106 = buffer.data(mh1 + 1106);
    const auto *mh1_1107 = buffer.data(mh1 + 1107);
    const auto *mh1_1109 = buffer.data(mh1 + 1109);
    const auto *mh1_1110 = buffer.data(mh1 + 1110);
    const auto *mh1_1111 = buffer.data(mh1 + 1111);
    const auto *mh1_1112 = buffer.data(mh1 + 1112);
    const auto *mh1_1134 = buffer.data(mh1 + 1134);
    const auto *mh1_1137 = buffer.data(mh1 + 1137);
    const auto *mh1_1139 = buffer.data(mh1 + 1139);
    const auto *mh1_1140 = buffer.data(mh1 + 1140);
    const auto *mh1_1143 = buffer.data(mh1 + 1143);
    const auto *mh1_1144 = buffer.data(mh1 + 1144);
    const auto *mh1_1146 = buffer.data(mh1 + 1146);
    const auto *mh1_1148 = buffer.data(mh1 + 1148);
    const auto *mh1_1149 = buffer.data(mh1 + 1149);
    const auto *mh1_1151 = buffer.data(mh1 + 1151);

    const auto *mi_1427 = buffer.data(mi + 1427);
    const auto *mi_1428 = buffer.data(mi + 1428);
    const auto *mi_1430 = buffer.data(mi + 1430);
    const auto *mi_1431 = buffer.data(mi + 1431);
    const auto *mi_1433 = buffer.data(mi + 1433);
    const auto *mi_1434 = buffer.data(mi + 1434);
    const auto *mi_1437 = buffer.data(mi + 1437);
    const auto *mi_1438 = buffer.data(mi + 1438);
    const auto *mi_1440 = buffer.data(mi + 1440);
    const auto *mi_1442 = buffer.data(mi + 1442);
    const auto *mi_1443 = buffer.data(mi + 1443);
    const auto *mi_1445 = buffer.data(mi + 1445);
    const auto *mi_1446 = buffer.data(mi + 1446);
    const auto *mi_1448 = buffer.data(mi + 1448);
    const auto *mi_1449 = buffer.data(mi + 1449);
    const auto *mi_1450 = buffer.data(mi + 1450);
    const auto *mi_1451 = buffer.data(mi + 1451);
    const auto *mi_1452 = buffer.data(mi + 1452);
    const auto *mi_1453 = buffer.data(mi + 1453);
    const auto *mi_1454 = buffer.data(mi + 1454);
    const auto *mi_1455 = buffer.data(mi + 1455);
    const auto *mi_1456 = buffer.data(mi + 1456);
    const auto *mi_1458 = buffer.data(mi + 1458);
    const auto *mi_1459 = buffer.data(mi + 1459);
    const auto *mi_1461 = buffer.data(mi + 1461);
    const auto *mi_1462 = buffer.data(mi + 1462);
    const auto *mi_1465 = buffer.data(mi + 1465);
    const auto *mi_1466 = buffer.data(mi + 1466);
    const auto *mi_1468 = buffer.data(mi + 1468);
    const auto *mi_1470 = buffer.data(mi + 1470);
    const auto *mi_1471 = buffer.data(mi + 1471);
    const auto *mi_1473 = buffer.data(mi + 1473);
    const auto *mi_1474 = buffer.data(mi + 1474);
    const auto *mi_1476 = buffer.data(mi + 1476);
    const auto *mi_1477 = buffer.data(mi + 1477);
    const auto *mi_1478 = buffer.data(mi + 1478);
    const auto *mi_1479 = buffer.data(mi + 1479);
    const auto *mi_1480 = buffer.data(mi + 1480);
    const auto *mi_1481 = buffer.data(mi + 1481);
    const auto *mi_1482 = buffer.data(mi + 1482);
    const auto *mi_1483 = buffer.data(mi + 1483);
    const auto *mi_1484 = buffer.data(mi + 1484);
    const auto *mi_1486 = buffer.data(mi + 1486);
    const auto *mi_1487 = buffer.data(mi + 1487);
    const auto *mi_1489 = buffer.data(mi + 1489);
    const auto *mi_1490 = buffer.data(mi + 1490);
    const auto *mi_1493 = buffer.data(mi + 1493);
    const auto *mi_1494 = buffer.data(mi + 1494);
    const auto *mi_1498 = buffer.data(mi + 1498);
    const auto *mi_1505 = buffer.data(mi + 1505);
    const auto *mi_1506 = buffer.data(mi + 1506);
    const auto *mi_1507 = buffer.data(mi + 1507);
    const auto *mi_1508 = buffer.data(mi + 1508);
    const auto *mi_1509 = buffer.data(mi + 1509);
    const auto *mi_1510 = buffer.data(mi + 1510);
    const auto *mi_1511 = buffer.data(mi + 1511);
    const auto *mi_1512 = buffer.data(mi + 1512);
    const auto *mi_1514 = buffer.data(mi + 1514);
    const auto *mi_1515 = buffer.data(mi + 1515);
    const auto *mi_1517 = buffer.data(mi + 1517);
    const auto *mi_1518 = buffer.data(mi + 1518);
    const auto *mi_1521 = buffer.data(mi + 1521);
    const auto *mi_1522 = buffer.data(mi + 1522);
    const auto *mi_1524 = buffer.data(mi + 1524);
    const auto *mi_1526 = buffer.data(mi + 1526);
    const auto *mi_1527 = buffer.data(mi + 1527);
    const auto *mi_1529 = buffer.data(mi + 1529);

#pragma omp simd aligned(t_1834, t_1835, t_1836, t_1837, pa_y, pb_x, pb_y, kk0_1223, kk1_1223, \
                         li_1175, li_1176, lk_1511, mh0_1071, mh1_1071, mi_1427, \
                         mi_1428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1834[k] = f_14 * li_1175[k]
                    + pb_y[k] * mi_1427[k];

        t_1835[k] = f_27 * kk0_1223[k]
                    - f_28 * kk1_1223[k]
                    + pa_y[k] * lk_1511[k];

        t_1836[k] = f_1 * mh0_1071[k]
                    - f_2 * mh1_1071[k]
                    + pb_x[k] * mi_1428[k];

        t_1837[k] = f_13 * li_1176[k]
                    + pb_y[k] * mi_1428[k];
    }

#pragma omp simd aligned(t_1838, t_1839, t_1840, pb_x, pb_y, pb_z, li_1148, li_1178, mh0_1074, \
                         mh1_1074, mi_1428, mi_1430, mi_1431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1838[k] = f_24 * li_1148[k]
                    + pb_z[k] * mi_1428[k];

        t_1839[k] = f_9 * mh0_1074[k]
                    - f_10 * mh1_1074[k]
                    + pb_x[k] * mi_1431[k];

        t_1840[k] = f_13 * li_1178[k]
                    + pb_y[k] * mi_1430[k];
    }

#pragma omp simd aligned(t_1841, t_1842, t_1843, t_1844, pb_x, pb_y, pb_z, li_1151, li_1181, \
                         mh0_1076, mh0_1077, mh1_1076, mh1_1077, mi_1431, mi_1433, \
                         mi_1434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1841[k] = f_9 * mh0_1076[k]
                    - f_10 * mh1_1076[k]
                    + pb_x[k] * mi_1433[k];

        t_1842[k] = f_7 * mh0_1077[k]
                    - f_8 * mh1_1077[k]
                    + pb_x[k] * mi_1434[k];

        t_1843[k] = f_24 * li_1151[k]
                    + pb_z[k] * mi_1431[k];

        t_1844[k] = f_13 * li_1181[k]
                    + pb_y[k] * mi_1433[k];
    }

#pragma omp simd aligned(t_1845, t_1846, t_1847, pb_x, pb_z, li_1154, mh0_1080, mh0_1081, \
                         mh1_1080, mh1_1081, mi_1434, mi_1437, \
                         mi_1438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1845[k] = f_7 * mh0_1080[k]
                    - f_8 * mh1_1080[k]
                    + pb_x[k] * mi_1437[k];

        t_1846[k] = f_5 * mh0_1081[k]
                    - f_6 * mh1_1081[k]
                    + pb_x[k] * mi_1438[k];

        t_1847[k] = f_24 * li_1154[k]
                    + pb_z[k] * mi_1434[k];
    }

#pragma omp simd aligned(t_1848, t_1849, t_1850, pb_x, pb_y, li_1185, mh0_1083, mh0_1085, \
                         mh1_1083, mh1_1085, mi_1437, mi_1440, \
                         mi_1442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1848[k] = f_5 * mh0_1083[k]
                    - f_6 * mh1_1083[k]
                    + pb_x[k] * mi_1440[k];

        t_1849[k] = f_13 * li_1185[k]
                    + pb_y[k] * mi_1437[k];

        t_1850[k] = f_5 * mh0_1085[k]
                    - f_6 * mh1_1085[k]
                    + pb_x[k] * mi_1442[k];
    }

#pragma omp simd aligned(t_1851, t_1852, t_1853, pb_x, pb_z, li_1158, mh0_1086, mh0_1088, \
                         mh1_1086, mh1_1088, mi_1438, mi_1443, \
                         mi_1445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1851[k] = f_3 * mh0_1086[k]
                    - f_4 * mh1_1086[k]
                    + pb_x[k] * mi_1443[k];

        t_1852[k] = f_24 * li_1158[k]
                    + pb_z[k] * mi_1438[k];

        t_1853[k] = f_3 * mh0_1088[k]
                    - f_4 * mh1_1088[k]
                    + pb_x[k] * mi_1445[k];
    }

#pragma omp simd aligned(t_1854, t_1855, t_1856, t_1857, pb_x, pb_y, li_1190, mh0_1089, \
                         mh0_1091, mh1_1089, mh1_1091, mi_1442, mi_1446, mi_1448, \
                         mi_1449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1854[k] = f_3 * mh0_1089[k]
                    - f_4 * mh1_1089[k]
                    + pb_x[k] * mi_1446[k];

        t_1855[k] = f_13 * li_1190[k]
                    + pb_y[k] * mi_1442[k];

        t_1856[k] = f_3 * mh0_1091[k]
                    - f_4 * mh1_1091[k]
                    + pb_x[k] * mi_1448[k];

        t_1857[k] = pb_x[k] * mi_1449[k];
    }

#pragma omp simd aligned(t_1858, t_1859, t_1860, t_1861, t_1862, t_1863, pb_x, mi_1450, \
                         mi_1451, mi_1452, mi_1453, mi_1454, mi_1455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1858[k] = pb_x[k] * mi_1450[k];

        t_1859[k] = pb_x[k] * mi_1451[k];

        t_1860[k] = pb_x[k] * mi_1452[k];

        t_1861[k] = pb_x[k] * mi_1453[k];

        t_1862[k] = pb_x[k] * mi_1454[k];

        t_1863[k] = pb_x[k] * mi_1455[k];
    }

#pragma omp simd aligned(t_1864, t_1865, t_1866, pa_z, pb_y, pb_z, kk0_1180, kk1_1180, \
                         li_1169, li_1199, lk_1504, mh0_1088, mh1_1088, mi_1449, \
                         mi_1451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1864[k] = f_25 * kk0_1180[k]
                    - f_26 * kk1_1180[k]
                    + pa_z[k] * lk_1504[k];

        t_1865[k] = f_24 * li_1169[k]
                    + pb_z[k] * mi_1449[k];

        t_1866[k] = f_13 * li_1199[k]
                    + f_9 * mh0_1088[k]
                    - f_10 * mh1_1088[k]
                    + pb_y[k] * mi_1451[k];
    }

#pragma omp simd aligned(t_1867, t_1868, t_1869, pb_y, li_1200, li_1201, li_1202, mh0_1089, \
                         mh0_1090, mh0_1091, mh1_1089, mh1_1090, mh1_1091, mi_1452, mi_1453, \
                         mi_1454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1867[k] = f_13 * li_1200[k]
                    + f_7 * mh0_1089[k]
                    - f_8 * mh1_1089[k]
                    + pb_y[k] * mi_1452[k];

        t_1868[k] = f_13 * li_1201[k]
                    + f_5 * mh0_1090[k]
                    - f_6 * mh1_1090[k]
                    + pb_y[k] * mi_1453[k];

        t_1869[k] = f_13 * li_1202[k]
                    + f_3 * mh0_1091[k]
                    - f_4 * mh1_1091[k]
                    + pb_y[k] * mi_1454[k];
    }

#pragma omp simd aligned(t_1870, t_1871, t_1872, t_1873, pa_y, pb_x, pb_y, kk0_1259, kk1_1259, \
                         li_1203, li_1204, lk_1547, mh0_1092, mh1_1092, mi_1455, \
                         mi_1456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1870[k] = f_13 * li_1203[k]
                    + pb_y[k] * mi_1455[k];

        t_1871[k] = f_22 * kk0_1259[k]
                    - f_23 * kk1_1259[k]
                    + pa_y[k] * lk_1547[k];

        t_1872[k] = f_1 * mh0_1092[k]
                    - f_2 * mh1_1092[k]
                    + pb_x[k] * mi_1456[k];

        t_1873[k] = f_12 * li_1204[k]
                    + pb_y[k] * mi_1456[k];
    }

#pragma omp simd aligned(t_1874, t_1875, t_1876, pb_x, pb_y, pb_z, li_1176, li_1206, mh0_1095, \
                         mh1_1095, mi_1456, mi_1458, mi_1459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1874[k] = f_17 * li_1176[k]
                    + pb_z[k] * mi_1456[k];

        t_1875[k] = f_9 * mh0_1095[k]
                    - f_10 * mh1_1095[k]
                    + pb_x[k] * mi_1459[k];

        t_1876[k] = f_12 * li_1206[k]
                    + pb_y[k] * mi_1458[k];
    }

#pragma omp simd aligned(t_1877, t_1878, t_1879, t_1880, pb_x, pb_y, pb_z, li_1179, li_1209, \
                         mh0_1097, mh0_1098, mh1_1097, mh1_1098, mi_1459, mi_1461, \
                         mi_1462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1877[k] = f_9 * mh0_1097[k]
                    - f_10 * mh1_1097[k]
                    + pb_x[k] * mi_1461[k];

        t_1878[k] = f_7 * mh0_1098[k]
                    - f_8 * mh1_1098[k]
                    + pb_x[k] * mi_1462[k];

        t_1879[k] = f_17 * li_1179[k]
                    + pb_z[k] * mi_1459[k];

        t_1880[k] = f_12 * li_1209[k]
                    + pb_y[k] * mi_1461[k];
    }

#pragma omp simd aligned(t_1881, t_1882, t_1883, pb_x, pb_z, li_1182, mh0_1101, mh0_1102, \
                         mh1_1101, mh1_1102, mi_1462, mi_1465, \
                         mi_1466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1881[k] = f_7 * mh0_1101[k]
                    - f_8 * mh1_1101[k]
                    + pb_x[k] * mi_1465[k];

        t_1882[k] = f_5 * mh0_1102[k]
                    - f_6 * mh1_1102[k]
                    + pb_x[k] * mi_1466[k];

        t_1883[k] = f_17 * li_1182[k]
                    + pb_z[k] * mi_1462[k];
    }

#pragma omp simd aligned(t_1884, t_1885, t_1886, pb_x, pb_y, li_1213, mh0_1104, mh0_1106, \
                         mh1_1104, mh1_1106, mi_1465, mi_1468, \
                         mi_1470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1884[k] = f_5 * mh0_1104[k]
                    - f_6 * mh1_1104[k]
                    + pb_x[k] * mi_1468[k];

        t_1885[k] = f_12 * li_1213[k]
                    + pb_y[k] * mi_1465[k];

        t_1886[k] = f_5 * mh0_1106[k]
                    - f_6 * mh1_1106[k]
                    + pb_x[k] * mi_1470[k];
    }

#pragma omp simd aligned(t_1887, t_1888, t_1889, pb_x, pb_z, li_1186, mh0_1107, mh0_1109, \
                         mh1_1107, mh1_1109, mi_1466, mi_1471, \
                         mi_1473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1887[k] = f_3 * mh0_1107[k]
                    - f_4 * mh1_1107[k]
                    + pb_x[k] * mi_1471[k];

        t_1888[k] = f_17 * li_1186[k]
                    + pb_z[k] * mi_1466[k];

        t_1889[k] = f_3 * mh0_1109[k]
                    - f_4 * mh1_1109[k]
                    + pb_x[k] * mi_1473[k];
    }

#pragma omp simd aligned(t_1890, t_1891, t_1892, t_1893, pb_x, pb_y, li_1218, mh0_1110, \
                         mh0_1112, mh1_1110, mh1_1112, mi_1470, mi_1474, mi_1476, \
                         mi_1477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1890[k] = f_3 * mh0_1110[k]
                    - f_4 * mh1_1110[k]
                    + pb_x[k] * mi_1474[k];

        t_1891[k] = f_12 * li_1218[k]
                    + pb_y[k] * mi_1470[k];

        t_1892[k] = f_3 * mh0_1112[k]
                    - f_4 * mh1_1112[k]
                    + pb_x[k] * mi_1476[k];

        t_1893[k] = pb_x[k] * mi_1477[k];
    }

#pragma omp simd aligned(t_1894, t_1895, t_1896, t_1897, t_1898, t_1899, pb_x, mi_1478, \
                         mi_1479, mi_1480, mi_1481, mi_1482, mi_1483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1894[k] = pb_x[k] * mi_1478[k];

        t_1895[k] = pb_x[k] * mi_1479[k];

        t_1896[k] = pb_x[k] * mi_1480[k];

        t_1897[k] = pb_x[k] * mi_1481[k];

        t_1898[k] = pb_x[k] * mi_1482[k];

        t_1899[k] = pb_x[k] * mi_1483[k];
    }

#pragma omp simd aligned(t_1900, t_1901, t_1902, pa_z, pb_y, pb_z, kk0_1216, kk1_1216, \
                         li_1197, li_1227, lk_1540, mh0_1109, mh1_1109, mi_1477, \
                         mi_1479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1900[k] = f_20 * kk0_1216[k]
                    - f_21 * kk1_1216[k]
                    + pa_z[k] * lk_1540[k];

        t_1901[k] = f_17 * li_1197[k]
                    + pb_z[k] * mi_1477[k];

        t_1902[k] = f_12 * li_1227[k]
                    + f_9 * mh0_1109[k]
                    - f_10 * mh1_1109[k]
                    + pb_y[k] * mi_1479[k];
    }

#pragma omp simd aligned(t_1903, t_1904, t_1905, pb_y, li_1228, li_1229, li_1230, mh0_1110, \
                         mh0_1111, mh0_1112, mh1_1110, mh1_1111, mh1_1112, mi_1480, mi_1481, \
                         mi_1482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1903[k] = f_12 * li_1228[k]
                    + f_7 * mh0_1110[k]
                    - f_8 * mh1_1110[k]
                    + pb_y[k] * mi_1480[k];

        t_1904[k] = f_12 * li_1229[k]
                    + f_5 * mh0_1111[k]
                    - f_6 * mh1_1111[k]
                    + pb_y[k] * mi_1481[k];

        t_1905[k] = f_12 * li_1230[k]
                    + f_3 * mh0_1112[k]
                    - f_4 * mh1_1112[k]
                    + pb_y[k] * mi_1482[k];
    }

#pragma omp simd aligned(t_1906, t_1907, t_1908, t_1909, t_1910, pa_y, pb_y, kk0_1295, \
                         kk1_1295, li_1231, li_1232, lk_1583, lk_1584, lk_1586, mi_1483, \
                         mi_1484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1906[k] = f_12 * li_1231[k]
                    + pb_y[k] * mi_1483[k];

        t_1907[k] = f_18 * kk0_1295[k]
                    - f_19 * kk1_1295[k]
                    + pa_y[k] * lk_1583[k];

        t_1908[k] = pa_y[k] * lk_1584[k];

        t_1909[k] = f_11 * li_1232[k]
                    + pb_y[k] * mi_1484[k];

        t_1910[k] = pa_y[k] * lk_1586[k];
    }

#pragma omp simd aligned(t_1911, t_1912, t_1913, t_1914, pa_y, pb_y, li_1233, li_1234, \
                         li_1235, lk_1587, lk_1589, lk_1590, mi_1486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1911[k] = f_12 * li_1233[k]
                    + pa_y[k] * lk_1587[k];

        t_1912[k] = f_11 * li_1234[k]
                    + pb_y[k] * mi_1486[k];

        t_1913[k] = pa_y[k] * lk_1589[k];

        t_1914[k] = f_13 * li_1235[k]
                    + pa_y[k] * lk_1590[k];
    }

#pragma omp simd aligned(t_1915, t_1916, t_1917, t_1918, pa_y, pb_y, pb_z, li_1207, li_1237, \
                         li_1238, lk_1593, lk_1594, mi_1487, mi_1489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1915[k] = f_16 * li_1207[k]
                    + pb_z[k] * mi_1487[k];

        t_1916[k] = f_11 * li_1237[k]
                    + pb_y[k] * mi_1489[k];

        t_1917[k] = pa_y[k] * lk_1593[k];

        t_1918[k] = f_14 * li_1238[k]
                    + pa_y[k] * lk_1594[k];
    }

#pragma omp simd aligned(t_1919, t_1920, t_1921, t_1922, pa_y, pb_y, pb_z, li_1210, li_1240, \
                         li_1241, lk_1596, lk_1598, mi_1490, mi_1493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1919[k] = f_16 * li_1210[k]
                    + pb_z[k] * mi_1490[k];

        t_1920[k] = f_12 * li_1240[k]
                    + pa_y[k] * lk_1596[k];

        t_1921[k] = f_11 * li_1241[k]
                    + pb_y[k] * mi_1493[k];

        t_1922[k] = pa_y[k] * lk_1598[k];
    }

#pragma omp simd aligned(t_1923, t_1924, t_1925, t_1926, pa_y, pb_z, li_1214, li_1242, \
                         li_1244, li_1245, lk_1599, lk_1601, lk_1602, \
                         mi_1494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1923[k] = f_15 * li_1242[k]
                    + pa_y[k] * lk_1599[k];

        t_1924[k] = f_16 * li_1214[k]
                    + pb_z[k] * mi_1494[k];

        t_1925[k] = f_13 * li_1244[k]
                    + pa_y[k] * lk_1601[k];

        t_1926[k] = f_12 * li_1245[k]
                    + pa_y[k] * lk_1602[k];
    }

#pragma omp simd aligned(t_1927, t_1928, t_1929, t_1930, t_1931, t_1932, pa_y, pb_x, pb_y, \
                         li_1246, lk_1604, mi_1498, mi_1505, mi_1506, mi_1507, \
                         mi_1508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1927[k] = f_11 * li_1246[k]
                    + pb_y[k] * mi_1498[k];

        t_1928[k] = pa_y[k] * lk_1604[k];

        t_1929[k] = pb_x[k] * mi_1505[k];

        t_1930[k] = pb_x[k] * mi_1506[k];

        t_1931[k] = pb_x[k] * mi_1507[k];

        t_1932[k] = pb_x[k] * mi_1508[k];
    }

#pragma omp simd aligned(t_1933, t_1934, t_1935, t_1936, t_1937, pa_y, pb_x, pb_z, li_1225, \
                         li_1253, lk_1612, mi_1505, mi_1509, mi_1510, \
                         mi_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1933[k] = pb_x[k] * mi_1509[k];

        t_1934[k] = pb_x[k] * mi_1510[k];

        t_1935[k] = pb_x[k] * mi_1511[k];

        t_1936[k] = f_17 * li_1253[k]
                    + pa_y[k] * lk_1612[k];

        t_1937[k] = f_16 * li_1225[k]
                    + pb_z[k] * mi_1505[k];
    }

#pragma omp simd aligned(t_1938, t_1939, t_1940, t_1941, pa_y, li_1255, li_1256, li_1257, \
                         li_1258, lk_1614, lk_1615, lk_1616, lk_1617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1938[k] = f_15 * li_1255[k]
                    + pa_y[k] * lk_1614[k];

        t_1939[k] = f_14 * li_1256[k]
                    + pa_y[k] * lk_1615[k];

        t_1940[k] = f_13 * li_1257[k]
                    + pa_y[k] * lk_1616[k];

        t_1941[k] = f_12 * li_1258[k]
                    + pa_y[k] * lk_1617[k];
    }

#pragma omp simd aligned(t_1942, t_1943, t_1944, t_1945, t_1946, pa_y, pb_x, pb_y, pb_z, \
                         li_1232, li_1259, lk_1619, mh0_1134, mh1_1134, mi_1511, \
                         mi_1512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1942[k] = f_11 * li_1259[k]
                    + pb_y[k] * mi_1511[k];

        t_1943[k] = pa_y[k] * lk_1619[k];

        t_1944[k] = f_1 * mh0_1134[k]
                    - f_2 * mh1_1134[k]
                    + pb_x[k] * mi_1512[k];

        t_1945[k] = pb_y[k] * mi_1512[k];

        t_1946[k] = f_0 * li_1232[k]
                    + pb_z[k] * mi_1512[k];
    }

#pragma omp simd aligned(t_1947, t_1948, t_1949, t_1950, pb_x, pb_y, mh0_1137, mh0_1139, \
                         mh0_1140, mh1_1137, mh1_1139, mh1_1140, mi_1514, mi_1515, mi_1517, \
                         mi_1518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1947[k] = f_9 * mh0_1137[k]
                    - f_10 * mh1_1137[k]
                    + pb_x[k] * mi_1515[k];

        t_1948[k] = pb_y[k] * mi_1514[k];

        t_1949[k] = f_9 * mh0_1139[k]
                    - f_10 * mh1_1139[k]
                    + pb_x[k] * mi_1517[k];

        t_1950[k] = f_7 * mh0_1140[k]
                    - f_8 * mh1_1140[k]
                    + pb_x[k] * mi_1518[k];
    }

#pragma omp simd aligned(t_1951, t_1952, t_1953, t_1954, pb_x, pb_y, pb_z, li_1235, mh0_1143, \
                         mh0_1144, mh1_1143, mh1_1144, mi_1515, mi_1517, mi_1521, \
                         mi_1522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1951[k] = f_0 * li_1235[k]
                    + pb_z[k] * mi_1515[k];

        t_1952[k] = pb_y[k] * mi_1517[k];

        t_1953[k] = f_7 * mh0_1143[k]
                    - f_8 * mh1_1143[k]
                    + pb_x[k] * mi_1521[k];

        t_1954[k] = f_5 * mh0_1144[k]
                    - f_6 * mh1_1144[k]
                    + pb_x[k] * mi_1522[k];
    }

#pragma omp simd aligned(t_1955, t_1956, t_1957, t_1958, pb_x, pb_y, pb_z, li_1238, mh0_1146, \
                         mh0_1148, mh1_1146, mh1_1148, mi_1518, mi_1521, mi_1524, \
                         mi_1526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1955[k] = f_0 * li_1238[k]
                    + pb_z[k] * mi_1518[k];

        t_1956[k] = f_5 * mh0_1146[k]
                    - f_6 * mh1_1146[k]
                    + pb_x[k] * mi_1524[k];

        t_1957[k] = pb_y[k] * mi_1521[k];

        t_1958[k] = f_5 * mh0_1148[k]
                    - f_6 * mh1_1148[k]
                    + pb_x[k] * mi_1526[k];
    }

#pragma omp simd aligned(t_1959, t_1960, t_1961, pb_x, pb_z, li_1242, mh0_1149, mh0_1151, \
                         mh1_1149, mh1_1151, mi_1522, mi_1527, \
                         mi_1529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1959[k] = f_3 * mh0_1149[k]
                    - f_4 * mh1_1149[k]
                    + pb_x[k] * mi_1527[k];

        t_1960[k] = f_0 * li_1242[k]
                    + pb_z[k] * mi_1522[k];

        t_1961[k] = f_3 * mh0_1151[k]
                    - f_4 * mh1_1151[k]
                    + pb_x[k] * mi_1529[k];
    }
}

static auto
compute_prim_mk_electron_repulsion_0_piece15(CSimdMatrix &buffer, const size_t target,
                                             const size_t pb, const size_t li, const size_t mh0,
                                             const size_t mh1, const size_t mi,
                                             const size_t ncols, const double alpha,
                                             const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);

    auto *t_1962 = buffer.data(target + 1962);
    auto *t_1963 = buffer.data(target + 1963);
    auto *t_1964 = buffer.data(target + 1964);
    auto *t_1965 = buffer.data(target + 1965);
    auto *t_1966 = buffer.data(target + 1966);
    auto *t_1967 = buffer.data(target + 1967);
    auto *t_1968 = buffer.data(target + 1968);
    auto *t_1969 = buffer.data(target + 1969);
    auto *t_1970 = buffer.data(target + 1970);
    auto *t_1971 = buffer.data(target + 1971);
    auto *t_1972 = buffer.data(target + 1972);
    auto *t_1973 = buffer.data(target + 1973);
    auto *t_1974 = buffer.data(target + 1974);
    auto *t_1975 = buffer.data(target + 1975);
    auto *t_1976 = buffer.data(target + 1976);
    auto *t_1977 = buffer.data(target + 1977);
    auto *t_1978 = buffer.data(target + 1978);
    auto *t_1979 = buffer.data(target + 1979);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *li_1253 = buffer.data(li + 1253);
    const auto *li_1259 = buffer.data(li + 1259);

    const auto *mh0_1149 = buffer.data(mh0 + 1149);
    const auto *mh0_1151 = buffer.data(mh0 + 1151);
    const auto *mh0_1152 = buffer.data(mh0 + 1152);
    const auto *mh0_1153 = buffer.data(mh0 + 1153);
    const auto *mh0_1154 = buffer.data(mh0 + 1154);

    const auto *mh1_1149 = buffer.data(mh1 + 1149);
    const auto *mh1_1151 = buffer.data(mh1 + 1151);
    const auto *mh1_1152 = buffer.data(mh1 + 1152);
    const auto *mh1_1153 = buffer.data(mh1 + 1153);
    const auto *mh1_1154 = buffer.data(mh1 + 1154);

    const auto *mi_1526 = buffer.data(mi + 1526);
    const auto *mi_1530 = buffer.data(mi + 1530);
    const auto *mi_1532 = buffer.data(mi + 1532);
    const auto *mi_1533 = buffer.data(mi + 1533);
    const auto *mi_1534 = buffer.data(mi + 1534);
    const auto *mi_1535 = buffer.data(mi + 1535);
    const auto *mi_1536 = buffer.data(mi + 1536);
    const auto *mi_1537 = buffer.data(mi + 1537);
    const auto *mi_1538 = buffer.data(mi + 1538);
    const auto *mi_1539 = buffer.data(mi + 1539);

#pragma omp simd aligned(t_1962, t_1963, t_1964, t_1965, t_1966, pb_x, pb_y, mh0_1152, \
                         mh0_1154, mh1_1152, mh1_1154, mi_1526, mi_1530, mi_1532, mi_1533, \
                         mi_1534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1962[k] = f_3 * mh0_1152[k]
                    - f_4 * mh1_1152[k]
                    + pb_x[k] * mi_1530[k];

        t_1963[k] = pb_y[k] * mi_1526[k];

        t_1964[k] = f_3 * mh0_1154[k]
                    - f_4 * mh1_1154[k]
                    + pb_x[k] * mi_1532[k];

        t_1965[k] = pb_x[k] * mi_1533[k];

        t_1966[k] = pb_x[k] * mi_1534[k];
    }

#pragma omp simd aligned(t_1967, t_1968, t_1969, t_1970, t_1971, t_1972, pb_x, pb_y, mh0_1149, \
                         mh1_1149, mi_1533, mi_1535, mi_1536, mi_1537, mi_1538, \
                         mi_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1967[k] = pb_x[k] * mi_1535[k];

        t_1968[k] = pb_x[k] * mi_1536[k];

        t_1969[k] = pb_x[k] * mi_1537[k];

        t_1970[k] = pb_x[k] * mi_1538[k];

        t_1971[k] = pb_x[k] * mi_1539[k];

        t_1972[k] = f_1 * mh0_1149[k]
                    - f_2 * mh1_1149[k]
                    + pb_y[k] * mi_1533[k];
    }

#pragma omp simd aligned(t_1973, t_1974, t_1975, pb_y, pb_z, li_1253, mh0_1151, mh0_1152, \
                         mh1_1151, mh1_1152, mi_1533, mi_1535, \
                         mi_1536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1973[k] = f_0 * li_1253[k]
                    + pb_z[k] * mi_1533[k];

        t_1974[k] = f_9 * mh0_1151[k]
                    - f_10 * mh1_1151[k]
                    + pb_y[k] * mi_1535[k];

        t_1975[k] = f_7 * mh0_1152[k]
                    - f_8 * mh1_1152[k]
                    + pb_y[k] * mi_1536[k];
    }

#pragma omp simd aligned(t_1976, t_1977, t_1978, t_1979, pb_y, pb_z, li_1259, mh0_1153, \
                         mh0_1154, mh1_1153, mh1_1154, mi_1537, mi_1538, \
                         mi_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1976[k] = f_5 * mh0_1153[k]
                    - f_6 * mh1_1153[k]
                    + pb_y[k] * mi_1537[k];

        t_1977[k] = f_3 * mh0_1154[k]
                    - f_4 * mh1_1154[k]
                    + pb_y[k] * mi_1538[k];

        t_1978[k] = pb_y[k] * mi_1539[k];

        t_1979[k] = f_0 * li_1259[k]
                    + f_1 * mh0_1154[k]
                    - f_2 * mh1_1154[k]
                    + pb_z[k] * mi_1539[k];
    }
}

auto
compute_prim_mk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t kk0, const size_t kk1,
                                     const size_t li, const size_t lk, const size_t mh0,
                                     const size_t mh1, const size_t mi, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    compute_prim_mk_electron_repulsion_0_piece0(buffer, target, pa, pb, kk0, kk1, li, lk, mh0,
                                                mh1, mi, ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece1(buffer, target, pa, pb, kk0, kk1, li, lk, mh0,
                                                mh1, mi, ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece2(buffer, target, pa, pb, kk0, kk1, li, lk, mh0,
                                                mh1, mi, ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece3(buffer, target, pa, pb, kk0, kk1, li, lk, mh0,
                                                mh1, mi, ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece4(buffer, target, pa, pb, kk0, kk1, li, lk, mh0,
                                                mh1, mi, ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece5(buffer, target, pa, pb, kk0, kk1, li, lk, mh0,
                                                mh1, mi, ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece6(buffer, target, pa, pb, kk0, kk1, li, lk, mh0,
                                                mh1, mi, ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece7(buffer, target, pa, pb, kk0, kk1, li, lk, mh0,
                                                mh1, mi, ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece8(buffer, target, pa, pb, kk0, kk1, li, lk, mh0,
                                                mh1, mi, ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece9(buffer, target, pa, pb, kk0, kk1, li, lk, mh0,
                                                mh1, mi, ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece10(buffer, target, pa, pb, kk0, kk1, li, lk, mh0,
                                                 mh1, mi, ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece11(buffer, target, pa, pb, li, lk, mi, ncols, p);

    compute_prim_mk_electron_repulsion_0_piece12(buffer, target, pa, pb, li, lk, mh0, mh1, mi,
                                                 ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece13(buffer, target, pa, pb, kk0, kk1, li, lk, mh0,
                                                 mh1, mi, ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece14(buffer, target, pa, pb, kk0, kk1, li, lk, mh0,
                                                 mh1, mi, ncols, alpha, beta, p);

    compute_prim_mk_electron_repulsion_0_piece15(buffer, target, pb, li, mh0, mh1, mi, ncols,
                                                 alpha, beta, p);
}

}  // namespace simdt2ceri
