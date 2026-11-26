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


#include "SimdElectronRepulsionVrrRecMH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_mh_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kh0,
                                            const size_t kh1, const size_t lg, const size_t lh,
                                            const size_t mf0, const size_t mf1, const size_t mg,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 4.0 / p;
    const auto f_11 = 2.5 / p;
    const auto f_12 = 0.5 / alpha;
    const auto f_13 = 0.5 * beta / (alpha * p);
    const auto f_14 = 3.5 / p;
    const auto f_15 = 3.0 / alpha;
    const auto f_16 = 3.0 * beta / (alpha * p);
    const auto f_17 = 1.0 / alpha;
    const auto f_18 = beta / (alpha * p);
    const auto f_19 = 3.0 / p;

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

    const auto *kh0_0 = buffer.data(kh0 + 0);
    const auto *kh0_21 = buffer.data(kh0 + 21);
    const auto *kh0_78 = buffer.data(kh0 + 78);
    const auto *kh0_125 = buffer.data(kh0 + 125);

    const auto *kh1_0 = buffer.data(kh1 + 0);
    const auto *kh1_21 = buffer.data(kh1 + 21);
    const auto *kh1_78 = buffer.data(kh1 + 78);
    const auto *kh1_125 = buffer.data(kh1 + 125);

    const auto *lg_0 = buffer.data(lg + 0);
    const auto *lg_1 = buffer.data(lg + 1);
    const auto *lg_2 = buffer.data(lg + 2);
    const auto *lg_3 = buffer.data(lg + 3);
    const auto *lg_5 = buffer.data(lg + 5);
    const auto *lg_10 = buffer.data(lg + 10);
    const auto *lg_11 = buffer.data(lg + 11);
    const auto *lg_12 = buffer.data(lg + 12);
    const auto *lg_13 = buffer.data(lg + 13);
    const auto *lg_14 = buffer.data(lg + 14);
    const auto *lg_15 = buffer.data(lg + 15);
    const auto *lg_18 = buffer.data(lg + 18);
    const auto *lg_20 = buffer.data(lg + 20);
    const auto *lg_25 = buffer.data(lg + 25);
    const auto *lg_27 = buffer.data(lg + 27);
    const auto *lg_28 = buffer.data(lg + 28);
    const auto *lg_29 = buffer.data(lg + 29);
    const auto *lg_30 = buffer.data(lg + 30);
    const auto *lg_32 = buffer.data(lg + 32);
    const auto *lg_33 = buffer.data(lg + 33);
    const auto *lg_35 = buffer.data(lg + 35);
    const auto *lg_40 = buffer.data(lg + 40);
    const auto *lg_41 = buffer.data(lg + 41);
    const auto *lg_42 = buffer.data(lg + 42);
    const auto *lg_43 = buffer.data(lg + 43);
    const auto *lg_44 = buffer.data(lg + 44);
    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_51 = buffer.data(lg + 51);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_58 = buffer.data(lg + 58);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_71 = buffer.data(lg + 71);
    const auto *lg_72 = buffer.data(lg + 72);
    const auto *lg_73 = buffer.data(lg + 73);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_84 = buffer.data(lg + 84);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_86 = buffer.data(lg + 86);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_96 = buffer.data(lg + 96);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_103 = buffer.data(lg + 103);
    const auto *lg_104 = buffer.data(lg + 104);

    const auto *lh_0 = buffer.data(lh + 0);
    const auto *lh_3 = buffer.data(lh + 3);
    const auto *lh_5 = buffer.data(lh + 5);
    const auto *lh_6 = buffer.data(lh + 6);
    const auto *lh_9 = buffer.data(lh + 9);
    const auto *lh_10 = buffer.data(lh + 10);
    const auto *lh_14 = buffer.data(lh + 14);
    const auto *lh_15 = buffer.data(lh + 15);
    const auto *lh_17 = buffer.data(lh + 17);
    const auto *lh_18 = buffer.data(lh + 18);
    const auto *lh_20 = buffer.data(lh + 20);
    const auto *lh_21 = buffer.data(lh + 21);
    const auto *lh_22 = buffer.data(lh + 22);
    const auto *lh_24 = buffer.data(lh + 24);
    const auto *lh_27 = buffer.data(lh + 27);
    const auto *lh_31 = buffer.data(lh + 31);
    const auto *lh_36 = buffer.data(lh + 36);
    const auto *lh_42 = buffer.data(lh + 42);
    const auto *lh_44 = buffer.data(lh + 44);
    const auto *lh_47 = buffer.data(lh + 47);
    const auto *lh_51 = buffer.data(lh + 51);
    const auto *lh_56 = buffer.data(lh + 56);
    const auto *lh_59 = buffer.data(lh + 59);
    const auto *lh_60 = buffer.data(lh + 60);
    const auto *lh_62 = buffer.data(lh + 62);
    const auto *lh_63 = buffer.data(lh + 63);
    const auto *lh_78 = buffer.data(lh + 78);
    const auto *lh_125 = buffer.data(lh + 125);

    const auto *mf0_0 = buffer.data(mf0 + 0);
    const auto *mf0_1 = buffer.data(mf0 + 1);
    const auto *mf0_2 = buffer.data(mf0 + 2);
    const auto *mf0_6 = buffer.data(mf0 + 6);
    const auto *mf0_8 = buffer.data(mf0 + 8);
    const auto *mf0_9 = buffer.data(mf0 + 9);
    const auto *mf0_30 = buffer.data(mf0 + 30);
    const auto *mf0_32 = buffer.data(mf0 + 32);
    const auto *mf0_33 = buffer.data(mf0 + 33);
    const auto *mf0_36 = buffer.data(mf0 + 36);
    const auto *mf0_37 = buffer.data(mf0 + 37);
    const auto *mf0_39 = buffer.data(mf0 + 39);
    const auto *mf0_50 = buffer.data(mf0 + 50);
    const auto *mf0_51 = buffer.data(mf0 + 51);
    const auto *mf0_55 = buffer.data(mf0 + 55);
    const auto *mf0_56 = buffer.data(mf0 + 56);
    const auto *mf0_58 = buffer.data(mf0 + 58);
    const auto *mf0_59 = buffer.data(mf0 + 59);
    const auto *mf0_60 = buffer.data(mf0 + 60);
    const auto *mf0_62 = buffer.data(mf0 + 62);
    const auto *mf0_63 = buffer.data(mf0 + 63);
    const auto *mf0_66 = buffer.data(mf0 + 66);

    const auto *mf1_0 = buffer.data(mf1 + 0);
    const auto *mf1_1 = buffer.data(mf1 + 1);
    const auto *mf1_2 = buffer.data(mf1 + 2);
    const auto *mf1_6 = buffer.data(mf1 + 6);
    const auto *mf1_8 = buffer.data(mf1 + 8);
    const auto *mf1_9 = buffer.data(mf1 + 9);
    const auto *mf1_30 = buffer.data(mf1 + 30);
    const auto *mf1_32 = buffer.data(mf1 + 32);
    const auto *mf1_33 = buffer.data(mf1 + 33);
    const auto *mf1_36 = buffer.data(mf1 + 36);
    const auto *mf1_37 = buffer.data(mf1 + 37);
    const auto *mf1_39 = buffer.data(mf1 + 39);
    const auto *mf1_50 = buffer.data(mf1 + 50);
    const auto *mf1_51 = buffer.data(mf1 + 51);
    const auto *mf1_55 = buffer.data(mf1 + 55);
    const auto *mf1_56 = buffer.data(mf1 + 56);
    const auto *mf1_58 = buffer.data(mf1 + 58);
    const auto *mf1_59 = buffer.data(mf1 + 59);
    const auto *mf1_60 = buffer.data(mf1 + 60);
    const auto *mf1_62 = buffer.data(mf1 + 62);
    const auto *mf1_63 = buffer.data(mf1 + 63);
    const auto *mf1_66 = buffer.data(mf1 + 66);

    const auto *mg_0 = buffer.data(mg + 0);
    const auto *mg_1 = buffer.data(mg + 1);
    const auto *mg_2 = buffer.data(mg + 2);
    const auto *mg_3 = buffer.data(mg + 3);
    const auto *mg_5 = buffer.data(mg + 5);
    const auto *mg_6 = buffer.data(mg + 6);
    const auto *mg_9 = buffer.data(mg + 9);
    const auto *mg_10 = buffer.data(mg + 10);
    const auto *mg_12 = buffer.data(mg + 12);
    const auto *mg_13 = buffer.data(mg + 13);
    const auto *mg_14 = buffer.data(mg + 14);
    const auto *mg_15 = buffer.data(mg + 15);
    const auto *mg_16 = buffer.data(mg + 16);
    const auto *mg_18 = buffer.data(mg + 18);
    const auto *mg_20 = buffer.data(mg + 20);
    const auto *mg_21 = buffer.data(mg + 21);
    const auto *mg_25 = buffer.data(mg + 25);
    const auto *mg_27 = buffer.data(mg + 27);
    const auto *mg_28 = buffer.data(mg + 28);
    const auto *mg_29 = buffer.data(mg + 29);
    const auto *mg_30 = buffer.data(mg + 30);
    const auto *mg_32 = buffer.data(mg + 32);
    const auto *mg_33 = buffer.data(mg + 33);
    const auto *mg_35 = buffer.data(mg + 35);
    const auto *mg_39 = buffer.data(mg + 39);
    const auto *mg_40 = buffer.data(mg + 40);
    const auto *mg_41 = buffer.data(mg + 41);
    const auto *mg_42 = buffer.data(mg + 42);
    const auto *mg_44 = buffer.data(mg + 44);
    const auto *mg_45 = buffer.data(mg + 45);
    const auto *mg_46 = buffer.data(mg + 46);
    const auto *mg_47 = buffer.data(mg + 47);
    const auto *mg_48 = buffer.data(mg + 48);
    const auto *mg_50 = buffer.data(mg + 50);
    const auto *mg_51 = buffer.data(mg + 51);
    const auto *mg_55 = buffer.data(mg + 55);
    const auto *mg_56 = buffer.data(mg + 56);
    const auto *mg_57 = buffer.data(mg + 57);
    const auto *mg_58 = buffer.data(mg + 58);
    const auto *mg_59 = buffer.data(mg + 59);
    const auto *mg_62 = buffer.data(mg + 62);
    const auto *mg_63 = buffer.data(mg + 63);
    const auto *mg_65 = buffer.data(mg + 65);
    const auto *mg_70 = buffer.data(mg + 70);
    const auto *mg_71 = buffer.data(mg + 71);
    const auto *mg_72 = buffer.data(mg + 72);
    const auto *mg_73 = buffer.data(mg + 73);
    const auto *mg_74 = buffer.data(mg + 74);
    const auto *mg_75 = buffer.data(mg + 75);
    const auto *mg_76 = buffer.data(mg + 76);
    const auto *mg_77 = buffer.data(mg + 77);
    const auto *mg_78 = buffer.data(mg + 78);
    const auto *mg_80 = buffer.data(mg + 80);
    const auto *mg_84 = buffer.data(mg + 84);
    const auto *mg_85 = buffer.data(mg + 85);
    const auto *mg_86 = buffer.data(mg + 86);
    const auto *mg_87 = buffer.data(mg + 87);
    const auto *mg_88 = buffer.data(mg + 88);
    const auto *mg_89 = buffer.data(mg + 89);
    const auto *mg_90 = buffer.data(mg + 90);
    const auto *mg_91 = buffer.data(mg + 91);
    const auto *mg_92 = buffer.data(mg + 92);
    const auto *mg_93 = buffer.data(mg + 93);
    const auto *mg_95 = buffer.data(mg + 95);
    const auto *mg_96 = buffer.data(mg + 96);
    const auto *mg_100 = buffer.data(mg + 100);
    const auto *mg_102 = buffer.data(mg + 102);
    const auto *mg_103 = buffer.data(mg + 103);
    const auto *mg_104 = buffer.data(mg + 104);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, lg_0, mf0_0, mf1_0, \
                         mg_0, mg_1, mg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lg_0[k]
                 + f_1 * mf0_0[k]
                 - f_2 * mf1_0[k]
                 + pb_x[k] * mg_0[k];

        t_1[k] = pb_y[k] * mg_0[k];

        t_2[k] = pb_z[k] * mg_0[k];

        t_3[k] = f_3 * mf0_0[k]
                 - f_4 * mf1_0[k]
                 + pb_y[k] * mg_1[k];

        t_4[k] = pb_y[k] * mg_2[k];

        t_5[k] = f_3 * mf0_0[k]
                 - f_4 * mf1_0[k]
                 + pb_z[k] * mg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, lg_10, mf0_1, mf0_2, \
                         mf1_1, mf1_2, mg_3, mg_5, mg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * mf0_1[k]
                 - f_6 * mf1_1[k]
                 + pb_y[k] * mg_3[k];

        t_7[k] = pb_z[k] * mg_3[k];

        t_8[k] = pb_y[k] * mg_5[k];

        t_9[k] = f_5 * mf0_2[k]
                 - f_6 * mf1_2[k]
                 + pb_z[k] * mg_5[k];

        t_10[k] = f_0 * lg_10[k]
                  + pb_x[k] * mg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pb_x, pb_y, pb_z, lg_12, lg_14, mg_6, mg_9, \
                         mg_12, mg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * mg_6[k];

        t_12[k] = f_0 * lg_12[k]
                  + pb_x[k] * mg_12[k];

        t_13[k] = pb_y[k] * mg_9[k];

        t_14[k] = f_0 * lg_14[k]
                  + pb_x[k] * mg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pb_z, mf0_6, mf0_8, mf0_9, mf1_6, \
                         mf1_8, mf1_9, mg_10, mg_12, mg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * mf0_6[k]
                  - f_2 * mf1_6[k]
                  + pb_y[k] * mg_10[k];

        t_16[k] = pb_z[k] * mg_10[k];

        t_17[k] = f_5 * mf0_8[k]
                  - f_6 * mf1_8[k]
                  + pb_y[k] * mg_12[k];

        t_18[k] = f_3 * mf0_9[k]
                  - f_4 * mf1_9[k]
                  + pb_y[k] * mg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pb_y, pb_z, lg_0, lh_0, mf0_9, \
                         mf1_9, mg_14, mg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * mg_14[k];

        t_20[k] = f_1 * mf0_9[k]
                  - f_2 * mf1_9[k]
                  + pb_z[k] * mg_14[k];

        t_21[k] = pa_y[k] * lh_0[k];

        t_22[k] = f_7 * lg_0[k]
                  + pb_y[k] * mg_15[k];

        t_23[k] = pb_z[k] * mg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pb_z, lg_1, lg_3, lh_3, lh_5, \
                         lh_6, mg_16, mg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * lg_1[k]
                  + pa_y[k] * lh_3[k];

        t_25[k] = pb_z[k] * mg_16[k];

        t_26[k] = pa_y[k] * lh_5[k];

        t_27[k] = f_9 * lg_3[k]
                  + pa_y[k] * lh_6[k];

        t_28[k] = pb_z[k] * mg_18[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, lg_5, lg_25, lh_9, \
                         mg_20, mg_21, mg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_7 * lg_5[k]
                  + pb_y[k] * mg_20[k];

        t_30[k] = pa_y[k] * lh_9[k];

        t_31[k] = f_10 * lg_25[k]
                  + pb_x[k] * mg_25[k];

        t_32[k] = pb_z[k] * mg_21[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pb_x, pb_z, lg_10, lg_27, lg_28, \
                         lh_14, lh_15, mg_25, mg_27, mg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_10 * lg_27[k]
                  + pb_x[k] * mg_27[k];

        t_34[k] = f_10 * lg_28[k]
                  + pb_x[k] * mg_28[k];

        t_35[k] = pa_y[k] * lh_14[k];

        t_36[k] = f_11 * lg_10[k]
                  + pa_y[k] * lh_15[k];

        t_37[k] = pb_z[k] * mg_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, lg_12, lg_13, lg_14, \
                         lh_0, lh_17, lh_18, lh_20, mg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_9 * lg_12[k]
                  + pa_y[k] * lh_17[k];

        t_39[k] = f_8 * lg_13[k]
                  + pa_y[k] * lh_18[k];

        t_40[k] = f_7 * lg_14[k]
                  + pb_y[k] * mg_29[k];

        t_41[k] = pa_y[k] * lh_20[k];

        t_42[k] = pa_z[k] * lh_0[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, t_48, pa_z, pb_y, pb_z, lg_0, lg_2, \
                         lh_3, lh_5, lh_6, mg_30, mg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * mg_30[k];

        t_44[k] = f_7 * lg_0[k]
                  + pb_z[k] * mg_30[k];

        t_45[k] = pa_z[k] * lh_3[k];

        t_46[k] = pb_y[k] * mg_32[k];

        t_47[k] = f_8 * lg_2[k]
                  + pa_z[k] * lh_5[k];

        t_48[k] = pa_z[k] * lh_6[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_z, pb_y, pb_z, lg_3, lg_5, lh_9, lh_10, \
                         mg_33, mg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_7 * lg_3[k]
                  + pb_z[k] * mg_33[k];

        t_50[k] = pb_y[k] * mg_35[k];

        t_51[k] = f_9 * lg_5[k]
                  + pa_z[k] * lh_9[k];

        t_52[k] = pa_z[k] * lh_10[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, lg_41, lg_42, lg_44, \
                         lh_15, mg_39, mg_41, mg_42, mg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * lg_41[k]
                  + pb_x[k] * mg_41[k];

        t_54[k] = f_10 * lg_42[k]
                  + pb_x[k] * mg_42[k];

        t_55[k] = pb_y[k] * mg_39[k];

        t_56[k] = f_10 * lg_44[k]
                  + pb_x[k] * mg_44[k];

        t_57[k] = pa_z[k] * lh_15[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pb_y, pb_z, lg_10, lg_11, lg_12, lh_17, \
                         lh_18, mg_40, mg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_7 * lg_10[k]
                  + pb_z[k] * mg_40[k];

        t_59[k] = f_8 * lg_11[k]
                  + pa_z[k] * lh_17[k];

        t_60[k] = f_9 * lg_12[k]
                  + pa_z[k] * lh_18[k];

        t_61[k] = pb_y[k] * mg_44[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, pb_z, kh0_0, kh1_0, lg_14, \
                         lg_15, lh_20, lh_21, mg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_11 * lg_14[k]
                  + pa_z[k] * lh_20[k];

        t_63[k] = f_12 * kh0_0[k]
                  - f_13 * kh1_0[k]
                  + pa_y[k] * lh_21[k];

        t_64[k] = f_8 * lg_15[k]
                  + pb_y[k] * mg_45[k];

        t_65[k] = pb_z[k] * mg_45[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, lg_48, mf0_30, mf0_33, mf1_30, mf1_33, \
                         mg_46, mg_47, mg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_14 * lg_48[k]
                  + f_5 * mf0_33[k]
                  - f_6 * mf1_33[k]
                  + pb_x[k] * mg_48[k];

        t_67[k] = pb_z[k] * mg_46[k];

        t_68[k] = f_3 * mf0_30[k]
                  - f_4 * mf1_30[k]
                  + pb_z[k] * mg_47[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_x, pb_y, pb_z, lg_20, lg_51, mf0_32, \
                         mf0_36, mf1_32, mf1_36, mg_48, mg_50, mg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_14 * lg_51[k]
                  + f_3 * mf0_36[k]
                  - f_4 * mf1_36[k]
                  + pb_x[k] * mg_51[k];

        t_70[k] = pb_z[k] * mg_48[k];

        t_71[k] = f_8 * lg_20[k]
                  + pb_y[k] * mg_50[k];

        t_72[k] = f_5 * mf0_32[k]
                  - f_6 * mf1_32[k]
                  + pb_z[k] * mg_50[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_x, pb_z, lg_55, lg_57, lg_58, lg_59, \
                         mg_51, mg_55, mg_57, mg_58, mg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_14 * lg_55[k]
                  + pb_x[k] * mg_55[k];

        t_74[k] = pb_z[k] * mg_51[k];

        t_75[k] = f_14 * lg_57[k]
                  + pb_x[k] * mg_57[k];

        t_76[k] = f_14 * lg_58[k]
                  + pb_x[k] * mg_58[k];

        t_77[k] = f_14 * lg_59[k]
                  + pb_x[k] * mg_59[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_z, kh0_78, kh1_78, lh_78, mf0_36, \
                         mf0_37, mf1_36, mf1_37, mg_55, mg_56, mg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_15 * kh0_78[k]
                  - f_16 * kh1_78[k]
                  + pa_x[k] * lh_78[k];

        t_79[k] = pb_z[k] * mg_55[k];

        t_80[k] = f_3 * mf0_36[k]
                  - f_4 * mf1_36[k]
                  + pb_z[k] * mg_56[k];

        t_81[k] = f_5 * mf0_37[k]
                  - f_6 * mf1_37[k]
                  + pb_z[k] * mg_57[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pa_y, pa_z, pb_y, pb_z, lg_29, lh_22, \
                         lh_42, lh_44, mf0_39, mf1_39, mg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_8 * lg_29[k]
                  + pb_y[k] * mg_59[k];

        t_83[k] = f_1 * mf0_39[k]
                  - f_2 * mf1_39[k]
                  + pb_z[k] * mg_59[k];

        t_84[k] = pa_y[k] * lh_42[k];

        t_85[k] = pa_z[k] * lh_22[k];

        t_86[k] = pa_y[k] * lh_44[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, pb_z, lg_18, lg_32, \
                         lh_24, lh_27, lh_47, mg_62, mg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * lh_24[k];

        t_88[k] = f_7 * lg_32[k]
                  + pb_y[k] * mg_62[k];

        t_89[k] = pa_y[k] * lh_47[k];

        t_90[k] = pa_z[k] * lh_27[k];

        t_91[k] = f_7 * lg_18[k]
                  + pb_z[k] * mg_63[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, pb_y, lg_35, lg_71, lh_31, \
                         lh_51, mg_65, mg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_7 * lg_35[k]
                  + pb_y[k] * mg_65[k];

        t_93[k] = pa_y[k] * lh_51[k];

        t_94[k] = pa_z[k] * lh_31[k];

        t_95[k] = f_14 * lg_71[k]
                  + pb_x[k] * mg_71[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_y, pa_z, pb_x, lg_72, lg_73, lh_36, lh_56, \
                         mg_72, mg_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_14 * lg_72[k]
                  + pb_x[k] * mg_72[k];

        t_97[k] = f_14 * lg_73[k]
                  + pb_x[k] * mg_73[k];

        t_98[k] = pa_y[k] * lh_56[k];

        t_99[k] = pa_z[k] * lh_36[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_y, pb_z, lg_25, lg_42, lg_43, \
                         lg_44, lh_59, lh_60, mg_70, mg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_7 * lg_25[k]
                   + pb_z[k] * mg_70[k];

        t_101[k] = f_9 * lg_42[k]
                   + pa_y[k] * lh_59[k];

        t_102[k] = f_8 * lg_43[k]
                   + pa_y[k] * lh_60[k];

        t_103[k] = f_7 * lg_44[k]
                   + pb_y[k] * mg_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_y, pa_z, pb_y, pb_z, kh0_0, kh1_0, \
                         lg_30, lh_42, lh_62, mg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * lh_62[k];

        t_105[k] = f_12 * kh0_0[k]
                   - f_13 * kh1_0[k]
                   + pa_z[k] * lh_42[k];

        t_106[k] = pb_y[k] * mg_75[k];

        t_107[k] = f_8 * lg_30[k]
                   + pb_z[k] * mg_75[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pb_x, pb_y, lg_80, mf0_50, mf0_55, mf1_50, \
                         mf1_55, mg_76, mg_77, mg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_3 * mf0_50[k]
                   - f_4 * mf1_50[k]
                   + pb_y[k] * mg_76[k];

        t_109[k] = pb_y[k] * mg_77[k];

        t_110[k] = f_14 * lg_80[k]
                   + f_5 * mf0_55[k]
                   - f_6 * mf1_55[k]
                   + pb_x[k] * mg_80[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pb_x, pb_y, pb_z, lg_33, lg_84, mf0_51, \
                         mf0_59, mf1_51, mf1_59, mg_78, mg_80, mg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_5 * mf0_51[k]
                   - f_6 * mf1_51[k]
                   + pb_y[k] * mg_78[k];

        t_112[k] = f_8 * lg_33[k]
                   + pb_z[k] * mg_78[k];

        t_113[k] = pb_y[k] * mg_80[k];

        t_114[k] = f_14 * lg_84[k]
                   + f_3 * mf0_59[k]
                   - f_4 * mf1_59[k]
                   + pb_x[k] * mg_84[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pb_x, pb_y, lg_85, lg_86, lg_87, \
                         lg_89, mg_84, mg_85, mg_86, mg_87, mg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_14 * lg_85[k]
                   + pb_x[k] * mg_85[k];

        t_116[k] = f_14 * lg_86[k]
                   + pb_x[k] * mg_86[k];

        t_117[k] = f_14 * lg_87[k]
                   + pb_x[k] * mg_87[k];

        t_118[k] = pb_y[k] * mg_84[k];

        t_119[k] = f_14 * lg_89[k]
                   + pb_x[k] * mg_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pb_z, lg_40, mf0_56, mf0_58, \
                         mf0_59, mf1_56, mf1_58, mf1_59, mg_85, mg_87, \
                         mg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_1 * mf0_56[k]
                   - f_2 * mf1_56[k]
                   + pb_y[k] * mg_85[k];

        t_121[k] = f_8 * lg_40[k]
                   + pb_z[k] * mg_85[k];

        t_122[k] = f_5 * mf0_58[k]
                   - f_6 * mf1_58[k]
                   + pb_y[k] * mg_87[k];

        t_123[k] = f_3 * mf0_59[k]
                   - f_4 * mf1_59[k]
                   + pb_y[k] * mg_88[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_x, pa_y, pb_y, kh0_21, kh0_125, \
                         kh1_21, kh1_125, lg_45, lh_63, lh_125, mg_89, \
                         mg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * mg_89[k];

        t_125[k] = f_15 * kh0_125[k]
                   - f_16 * kh1_125[k]
                   + pa_x[k] * lh_125[k];

        t_126[k] = f_17 * kh0_21[k]
                   - f_18 * kh1_21[k]
                   + pa_y[k] * lh_63[k];

        t_127[k] = f_9 * lg_45[k]
                   + pb_y[k] * mg_90[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pb_x, pb_z, lg_93, mf0_60, mf0_63, \
                         mf1_60, mf1_63, mg_90, mg_91, mg_92, mg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = pb_z[k] * mg_90[k];

        t_129[k] = f_19 * lg_93[k]
                   + f_5 * mf0_63[k]
                   - f_6 * mf1_63[k]
                   + pb_x[k] * mg_93[k];

        t_130[k] = pb_z[k] * mg_91[k];

        t_131[k] = f_3 * mf0_60[k]
                   - f_4 * mf1_60[k]
                   + pb_z[k] * mg_92[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, pb_y, pb_z, lg_50, lg_96, mf0_62, \
                         mf0_66, mf1_62, mf1_66, mg_93, mg_95, mg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_19 * lg_96[k]
                   + f_3 * mf0_66[k]
                   - f_4 * mf1_66[k]
                   + pb_x[k] * mg_96[k];

        t_133[k] = pb_z[k] * mg_93[k];

        t_134[k] = f_9 * lg_50[k]
                   + pb_y[k] * mg_95[k];

        t_135[k] = f_5 * mf0_62[k]
                   - f_6 * mf1_62[k]
                   + pb_z[k] * mg_95[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, pb_x, pb_z, lg_100, lg_102, \
                         lg_103, lg_104, mg_96, mg_100, mg_102, mg_103, \
                         mg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_19 * lg_100[k]
                   + pb_x[k] * mg_100[k];

        t_137[k] = pb_z[k] * mg_96[k];

        t_138[k] = f_19 * lg_102[k]
                   + pb_x[k] * mg_102[k];

        t_139[k] = f_19 * lg_103[k]
                   + pb_x[k] * mg_103[k];

        t_140[k] = f_19 * lg_104[k]
                   + pb_x[k] * mg_104[k];
    }
}

static auto
compute_prim_mh_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kh0,
                                            const size_t kh1, const size_t lg, const size_t lh,
                                            const size_t mf0, const size_t mf1, const size_t mg,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_11 = 2.5 / p;
    const auto f_12 = 0.5 / alpha;
    const auto f_13 = 0.5 * beta / (alpha * p);
    const auto f_17 = 1.0 / alpha;
    const auto f_18 = beta / (alpha * p);
    const auto f_19 = 3.0 / p;
    const auto f_20 = 2.5 / alpha;
    const auto f_21 = 2.5 * beta / (alpha * p);
    const auto f_22 = 1.5 / alpha;
    const auto f_23 = 1.5 * beta / (alpha * p);
    const auto f_24 = 2.0 / p;
    const auto f_25 = 2.0 / alpha;
    const auto f_26 = 2.0 * beta / (alpha * p);

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
    auto *t_275 = buffer.data(target + 275);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kh0_42 = buffer.data(kh0 + 42);
    const auto *kh0_63 = buffer.data(kh0 + 63);
    const auto *kh0_66 = buffer.data(kh0 + 66);
    const auto *kh0_69 = buffer.data(kh0 + 69);
    const auto *kh0_105 = buffer.data(kh0 + 105);
    const auto *kh0_110 = buffer.data(kh0 + 110);
    const auto *kh0_114 = buffer.data(kh0 + 114);
    const auto *kh0_141 = buffer.data(kh0 + 141);
    const auto *kh0_209 = buffer.data(kh0 + 209);
    const auto *kh0_225 = buffer.data(kh0 + 225);
    const auto *kh0_267 = buffer.data(kh0 + 267);
    const auto *kh0_269 = buffer.data(kh0 + 269);
    const auto *kh0_270 = buffer.data(kh0 + 270);
    const auto *kh0_272 = buffer.data(kh0 + 272);

    const auto *kh1_42 = buffer.data(kh1 + 42);
    const auto *kh1_63 = buffer.data(kh1 + 63);
    const auto *kh1_66 = buffer.data(kh1 + 66);
    const auto *kh1_69 = buffer.data(kh1 + 69);
    const auto *kh1_105 = buffer.data(kh1 + 105);
    const auto *kh1_110 = buffer.data(kh1 + 110);
    const auto *kh1_114 = buffer.data(kh1 + 114);
    const auto *kh1_141 = buffer.data(kh1 + 141);
    const auto *kh1_209 = buffer.data(kh1 + 209);
    const auto *kh1_225 = buffer.data(kh1 + 225);
    const auto *kh1_267 = buffer.data(kh1 + 267);
    const auto *kh1_269 = buffer.data(kh1 + 269);
    const auto *kh1_270 = buffer.data(kh1 + 270);
    const auto *kh1_272 = buffer.data(kh1 + 272);

    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_47 = buffer.data(lg + 47);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_56 = buffer.data(lg + 56);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_62 = buffer.data(lg + 62);
    const auto *lg_63 = buffer.data(lg + 63);
    const auto *lg_65 = buffer.data(lg + 65);
    const auto *lg_70 = buffer.data(lg + 70);
    const auto *lg_74 = buffer.data(lg + 74);
    const auto *lg_75 = buffer.data(lg + 75);
    const auto *lg_76 = buffer.data(lg + 76);
    const auto *lg_77 = buffer.data(lg + 77);
    const auto *lg_78 = buffer.data(lg + 78);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_88 = buffer.data(lg + 88);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_90 = buffer.data(lg + 90);
    const auto *lg_92 = buffer.data(lg + 92);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_95 = buffer.data(lg + 95);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_101 = buffer.data(lg + 101);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_104 = buffer.data(lg + 104);
    const auto *lg_105 = buffer.data(lg + 105);
    const auto *lg_107 = buffer.data(lg + 107);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_110 = buffer.data(lg + 110);
    const auto *lg_115 = buffer.data(lg + 115);
    const auto *lg_116 = buffer.data(lg + 116);
    const auto *lg_117 = buffer.data(lg + 117);
    const auto *lg_118 = buffer.data(lg + 118);
    const auto *lg_119 = buffer.data(lg + 119);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_122 = buffer.data(lg + 122);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_131 = buffer.data(lg + 131);
    const auto *lg_132 = buffer.data(lg + 132);
    const auto *lg_133 = buffer.data(lg + 133);
    const auto *lg_134 = buffer.data(lg + 134);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_144 = buffer.data(lg + 144);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_146 = buffer.data(lg + 146);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_149 = buffer.data(lg + 149);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_156 = buffer.data(lg + 156);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_163 = buffer.data(lg + 163);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_176 = buffer.data(lg + 176);
    const auto *lg_177 = buffer.data(lg + 177);
    const auto *lg_178 = buffer.data(lg + 178);
    const auto *lg_179 = buffer.data(lg + 179);
    const auto *lg_190 = buffer.data(lg + 190);
    const auto *lg_191 = buffer.data(lg + 191);
    const auto *lg_192 = buffer.data(lg + 192);
    const auto *lg_193 = buffer.data(lg + 193);
    const auto *lg_194 = buffer.data(lg + 194);

    const auto *lh_63 = buffer.data(lh + 63);
    const auto *lh_64 = buffer.data(lh + 64);
    const auto *lh_66 = buffer.data(lh + 66);
    const auto *lh_68 = buffer.data(lh + 68);
    const auto *lh_69 = buffer.data(lh + 69);
    const auto *lh_72 = buffer.data(lh + 72);
    const auto *lh_73 = buffer.data(lh + 73);
    const auto *lh_78 = buffer.data(lh + 78);
    const auto *lh_80 = buffer.data(lh + 80);
    const auto *lh_81 = buffer.data(lh + 81);
    const auto *lh_83 = buffer.data(lh + 83);
    const auto *lh_105 = buffer.data(lh + 105);
    const auto *lh_107 = buffer.data(lh + 107);
    const auto *lh_108 = buffer.data(lh + 108);
    const auto *lh_110 = buffer.data(lh + 110);
    const auto *lh_111 = buffer.data(lh + 111);
    const auto *lh_114 = buffer.data(lh + 114);
    const auto *lh_119 = buffer.data(lh + 119);
    const auto *lh_120 = buffer.data(lh + 120);
    const auto *lh_122 = buffer.data(lh + 122);
    const auto *lh_123 = buffer.data(lh + 123);
    const auto *lh_125 = buffer.data(lh + 125);
    const auto *lh_126 = buffer.data(lh + 126);
    const auto *lh_127 = buffer.data(lh + 127);
    const auto *lh_129 = buffer.data(lh + 129);
    const auto *lh_131 = buffer.data(lh + 131);
    const auto *lh_132 = buffer.data(lh + 132);
    const auto *lh_135 = buffer.data(lh + 135);
    const auto *lh_136 = buffer.data(lh + 136);
    const auto *lh_141 = buffer.data(lh + 141);
    const auto *lh_143 = buffer.data(lh + 143);
    const auto *lh_144 = buffer.data(lh + 144);
    const auto *lh_146 = buffer.data(lh + 146);
    const auto *lh_150 = buffer.data(lh + 150);
    const auto *lh_153 = buffer.data(lh + 153);
    const auto *lh_168 = buffer.data(lh + 168);
    const auto *lh_173 = buffer.data(lh + 173);
    const auto *lh_177 = buffer.data(lh + 177);
    const auto *lh_189 = buffer.data(lh + 189);
    const auto *lh_191 = buffer.data(lh + 191);
    const auto *lh_209 = buffer.data(lh + 209);
    const auto *lh_225 = buffer.data(lh + 225);
    const auto *lh_267 = buffer.data(lh + 267);
    const auto *lh_269 = buffer.data(lh + 269);
    const auto *lh_270 = buffer.data(lh + 270);
    const auto *lh_272 = buffer.data(lh + 272);

    const auto *mf0_66 = buffer.data(mf0 + 66);
    const auto *mf0_67 = buffer.data(mf0 + 67);
    const auto *mf0_69 = buffer.data(mf0 + 69);
    const auto *mf0_90 = buffer.data(mf0 + 90);
    const auto *mf0_91 = buffer.data(mf0 + 91);
    const auto *mf0_95 = buffer.data(mf0 + 95);
    const auto *mf0_96 = buffer.data(mf0 + 96);
    const auto *mf0_98 = buffer.data(mf0 + 98);
    const auto *mf0_99 = buffer.data(mf0 + 99);
    const auto *mf0_100 = buffer.data(mf0 + 100);
    const auto *mf0_102 = buffer.data(mf0 + 102);
    const auto *mf0_103 = buffer.data(mf0 + 103);
    const auto *mf0_106 = buffer.data(mf0 + 106);
    const auto *mf0_107 = buffer.data(mf0 + 107);
    const auto *mf0_109 = buffer.data(mf0 + 109);

    const auto *mf1_66 = buffer.data(mf1 + 66);
    const auto *mf1_67 = buffer.data(mf1 + 67);
    const auto *mf1_69 = buffer.data(mf1 + 69);
    const auto *mf1_90 = buffer.data(mf1 + 90);
    const auto *mf1_91 = buffer.data(mf1 + 91);
    const auto *mf1_95 = buffer.data(mf1 + 95);
    const auto *mf1_96 = buffer.data(mf1 + 96);
    const auto *mf1_98 = buffer.data(mf1 + 98);
    const auto *mf1_99 = buffer.data(mf1 + 99);
    const auto *mf1_100 = buffer.data(mf1 + 100);
    const auto *mf1_102 = buffer.data(mf1 + 102);
    const auto *mf1_103 = buffer.data(mf1 + 103);
    const auto *mf1_106 = buffer.data(mf1 + 106);
    const auto *mf1_107 = buffer.data(mf1 + 107);
    const auto *mf1_109 = buffer.data(mf1 + 109);

    const auto *mg_100 = buffer.data(mg + 100);
    const auto *mg_101 = buffer.data(mg + 101);
    const auto *mg_102 = buffer.data(mg + 102);
    const auto *mg_104 = buffer.data(mg + 104);
    const auto *mg_105 = buffer.data(mg + 105);
    const auto *mg_107 = buffer.data(mg + 107);
    const auto *mg_108 = buffer.data(mg + 108);
    const auto *mg_110 = buffer.data(mg + 110);
    const auto *mg_115 = buffer.data(mg + 115);
    const auto *mg_116 = buffer.data(mg + 116);
    const auto *mg_117 = buffer.data(mg + 117);
    const auto *mg_118 = buffer.data(mg + 118);
    const auto *mg_119 = buffer.data(mg + 119);
    const auto *mg_120 = buffer.data(mg + 120);
    const auto *mg_122 = buffer.data(mg + 122);
    const auto *mg_123 = buffer.data(mg + 123);
    const auto *mg_125 = buffer.data(mg + 125);
    const auto *mg_130 = buffer.data(mg + 130);
    const auto *mg_131 = buffer.data(mg + 131);
    const auto *mg_132 = buffer.data(mg + 132);
    const auto *mg_133 = buffer.data(mg + 133);
    const auto *mg_134 = buffer.data(mg + 134);
    const auto *mg_135 = buffer.data(mg + 135);
    const auto *mg_136 = buffer.data(mg + 136);
    const auto *mg_137 = buffer.data(mg + 137);
    const auto *mg_138 = buffer.data(mg + 138);
    const auto *mg_140 = buffer.data(mg + 140);
    const auto *mg_144 = buffer.data(mg + 144);
    const auto *mg_145 = buffer.data(mg + 145);
    const auto *mg_146 = buffer.data(mg + 146);
    const auto *mg_147 = buffer.data(mg + 147);
    const auto *mg_148 = buffer.data(mg + 148);
    const auto *mg_149 = buffer.data(mg + 149);
    const auto *mg_150 = buffer.data(mg + 150);
    const auto *mg_151 = buffer.data(mg + 151);
    const auto *mg_152 = buffer.data(mg + 152);
    const auto *mg_153 = buffer.data(mg + 153);
    const auto *mg_155 = buffer.data(mg + 155);
    const auto *mg_156 = buffer.data(mg + 156);
    const auto *mg_160 = buffer.data(mg + 160);
    const auto *mg_161 = buffer.data(mg + 161);
    const auto *mg_162 = buffer.data(mg + 162);
    const auto *mg_163 = buffer.data(mg + 163);
    const auto *mg_164 = buffer.data(mg + 164);
    const auto *mg_165 = buffer.data(mg + 165);
    const auto *mg_167 = buffer.data(mg + 167);
    const auto *mg_168 = buffer.data(mg + 168);
    const auto *mg_170 = buffer.data(mg + 170);
    const auto *mg_175 = buffer.data(mg + 175);
    const auto *mg_176 = buffer.data(mg + 176);
    const auto *mg_177 = buffer.data(mg + 177);
    const auto *mg_178 = buffer.data(mg + 178);
    const auto *mg_179 = buffer.data(mg + 179);
    const auto *mg_180 = buffer.data(mg + 180);
    const auto *mg_182 = buffer.data(mg + 182);
    const auto *mg_183 = buffer.data(mg + 183);
    const auto *mg_185 = buffer.data(mg + 185);
    const auto *mg_190 = buffer.data(mg + 190);
    const auto *mg_191 = buffer.data(mg + 191);
    const auto *mg_192 = buffer.data(mg + 192);
    const auto *mg_193 = buffer.data(mg + 193);
    const auto *mg_194 = buffer.data(mg + 194);
    const auto *mg_195 = buffer.data(mg + 195);

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_x, pb_z, kh0_141, kh1_141, lh_141, \
                         mf0_66, mf0_67, mf1_66, mf1_67, mg_100, mg_101, \
                         mg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_20 * kh0_141[k]
                   - f_21 * kh1_141[k]
                   + pa_x[k] * lh_141[k];

        t_142[k] = pb_z[k] * mg_100[k];

        t_143[k] = f_3 * mf0_66[k]
                   - f_4 * mf1_66[k]
                   + pb_z[k] * mg_101[k];

        t_144[k] = f_5 * mf0_67[k]
                   - f_6 * mf1_67[k]
                   + pb_z[k] * mg_102[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pa_z, pb_y, pb_z, lg_45, lg_59, \
                         lh_63, lh_64, mf0_69, mf1_69, mg_104, mg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_9 * lg_59[k]
                   + pb_y[k] * mg_104[k];

        t_146[k] = f_1 * mf0_69[k]
                   - f_2 * mf1_69[k]
                   + pb_z[k] * mg_104[k];

        t_147[k] = pa_z[k] * lh_63[k];

        t_148[k] = pa_z[k] * lh_64[k];

        t_149[k] = f_7 * lg_45[k]
                   + pb_z[k] * mg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pa_z, pb_y, pb_z, lg_47, lg_48, \
                         lg_62, lh_66, lh_68, lh_69, mg_107, mg_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * lh_66[k];

        t_151[k] = f_8 * lg_62[k]
                   + pb_y[k] * mg_107[k];

        t_152[k] = f_8 * lg_47[k]
                   + pa_z[k] * lh_68[k];

        t_153[k] = pa_z[k] * lh_69[k];

        t_154[k] = f_7 * lg_48[k]
                   + pb_z[k] * mg_108[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_z, pb_x, pb_y, lg_50, lg_65, lg_116, \
                         lh_72, lh_73, mg_110, mg_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_8 * lg_65[k]
                   + pb_y[k] * mg_110[k];

        t_156[k] = f_9 * lg_50[k]
                   + pa_z[k] * lh_72[k];

        t_157[k] = pa_z[k] * lh_73[k];

        t_158[k] = f_19 * lg_116[k]
                   + pb_x[k] * mg_116[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_z, pb_x, lg_117, lg_118, lg_119, \
                         lh_78, mg_117, mg_118, mg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_19 * lg_117[k]
                   + pb_x[k] * mg_117[k];

        t_160[k] = f_19 * lg_118[k]
                   + pb_x[k] * mg_118[k];

        t_161[k] = f_19 * lg_119[k]
                   + pb_x[k] * mg_119[k];

        t_162[k] = pa_z[k] * lh_78[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_z, pb_y, pb_z, lg_55, lg_56, lg_57, \
                         lg_74, lh_80, lh_81, mg_115, mg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_7 * lg_55[k]
                   + pb_z[k] * mg_115[k];

        t_164[k] = f_8 * lg_56[k]
                   + pa_z[k] * lh_80[k];

        t_165[k] = f_9 * lg_57[k]
                   + pa_z[k] * lh_81[k];

        t_166[k] = f_8 * lg_74[k]
                   + pb_y[k] * mg_119[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pa_y, pa_z, pb_y, lg_59, lg_75, \
                         lg_76, lh_83, lh_105, lh_107, lh_108, mg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_11 * lg_59[k]
                   + pa_z[k] * lh_83[k];

        t_168[k] = pa_y[k] * lh_105[k];

        t_169[k] = f_7 * lg_75[k]
                   + pb_y[k] * mg_120[k];

        t_170[k] = pa_y[k] * lh_107[k];

        t_171[k] = f_8 * lg_76[k]
                   + pa_y[k] * lh_108[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pb_y, pb_z, lg_63, lg_77, lg_78, \
                         lh_110, lh_111, mg_122, mg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_7 * lg_77[k]
                   + pb_y[k] * mg_122[k];

        t_173[k] = pa_y[k] * lh_110[k];

        t_174[k] = f_9 * lg_78[k]
                   + pa_y[k] * lh_111[k];

        t_175[k] = f_8 * lg_63[k]
                   + pb_z[k] * mg_123[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_x, pb_y, lg_80, lg_130, lg_131, \
                         lh_114, mg_125, mg_130, mg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_7 * lg_80[k]
                   + pb_y[k] * mg_125[k];

        t_177[k] = pa_y[k] * lh_114[k];

        t_178[k] = f_19 * lg_130[k]
                   + pb_x[k] * mg_130[k];

        t_179[k] = f_19 * lg_131[k]
                   + pb_x[k] * mg_131[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pb_x, lg_85, lg_132, lg_133, \
                         lh_119, lh_120, mg_132, mg_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_19 * lg_132[k]
                   + pb_x[k] * mg_132[k];

        t_181[k] = f_19 * lg_133[k]
                   + pb_x[k] * mg_133[k];

        t_182[k] = pa_y[k] * lh_119[k];

        t_183[k] = f_11 * lg_85[k]
                   + pa_y[k] * lh_120[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pb_y, pb_z, lg_70, lg_87, lg_88, \
                         lg_89, lh_122, lh_123, mg_130, mg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_8 * lg_70[k]
                   + pb_z[k] * mg_130[k];

        t_185[k] = f_9 * lg_87[k]
                   + pa_y[k] * lh_122[k];

        t_186[k] = f_8 * lg_88[k]
                   + pa_y[k] * lh_123[k];

        t_187[k] = f_7 * lg_89[k]
                   + pb_y[k] * mg_134[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pa_z, pb_y, pb_z, kh0_42, kh1_42, \
                         lg_75, lh_105, lh_125, mg_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pa_y[k] * lh_125[k];

        t_189[k] = f_17 * kh0_42[k]
                   - f_18 * kh1_42[k]
                   + pa_z[k] * lh_105[k];

        t_190[k] = pb_y[k] * mg_135[k];

        t_191[k] = f_9 * lg_75[k]
                   + pb_z[k] * mg_135[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_x, pb_y, lg_140, mf0_90, mf0_95, mf1_90, \
                         mf1_95, mg_136, mg_137, mg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_3 * mf0_90[k]
                   - f_4 * mf1_90[k]
                   + pb_y[k] * mg_136[k];

        t_193[k] = pb_y[k] * mg_137[k];

        t_194[k] = f_19 * lg_140[k]
                   + f_5 * mf0_95[k]
                   - f_6 * mf1_95[k]
                   + pb_x[k] * mg_140[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pb_x, pb_y, pb_z, lg_78, lg_144, mf0_91, \
                         mf0_99, mf1_91, mf1_99, mg_138, mg_140, \
                         mg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_5 * mf0_91[k]
                   - f_6 * mf1_91[k]
                   + pb_y[k] * mg_138[k];

        t_196[k] = f_9 * lg_78[k]
                   + pb_z[k] * mg_138[k];

        t_197[k] = pb_y[k] * mg_140[k];

        t_198[k] = f_19 * lg_144[k]
                   + f_3 * mf0_99[k]
                   - f_4 * mf1_99[k]
                   + pb_x[k] * mg_144[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, pb_x, pb_y, lg_145, lg_146, \
                         lg_147, lg_149, mg_144, mg_145, mg_146, mg_147, \
                         mg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_19 * lg_145[k]
                   + pb_x[k] * mg_145[k];

        t_200[k] = f_19 * lg_146[k]
                   + pb_x[k] * mg_146[k];

        t_201[k] = f_19 * lg_147[k]
                   + pb_x[k] * mg_147[k];

        t_202[k] = pb_y[k] * mg_144[k];

        t_203[k] = f_19 * lg_149[k]
                   + pb_x[k] * mg_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_y, pb_z, lg_85, mf0_96, mf0_98, \
                         mf0_99, mf1_96, mf1_98, mf1_99, mg_145, mg_147, \
                         mg_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * mf0_96[k]
                   - f_2 * mf1_96[k]
                   + pb_y[k] * mg_145[k];

        t_205[k] = f_9 * lg_85[k]
                   + pb_z[k] * mg_145[k];

        t_206[k] = f_5 * mf0_98[k]
                   - f_6 * mf1_98[k]
                   + pb_y[k] * mg_147[k];

        t_207[k] = f_3 * mf0_99[k]
                   - f_4 * mf1_99[k]
                   + pb_y[k] * mg_148[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_x, pa_y, pb_y, kh0_63, kh0_209, \
                         kh1_63, kh1_209, lg_90, lh_126, lh_209, mg_149, \
                         mg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pb_y[k] * mg_149[k];

        t_209[k] = f_20 * kh0_209[k]
                   - f_21 * kh1_209[k]
                   + pa_x[k] * lh_209[k];

        t_210[k] = f_22 * kh0_63[k]
                   - f_23 * kh1_63[k]
                   + pa_y[k] * lh_126[k];

        t_211[k] = f_24 * lg_90[k]
                   + pb_y[k] * mg_150[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pb_x, pb_z, lg_153, mf0_100, mf0_103, \
                         mf1_100, mf1_103, mg_150, mg_151, mg_152, \
                         mg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pb_z[k] * mg_150[k];

        t_213[k] = f_11 * lg_153[k]
                   + f_5 * mf0_103[k]
                   - f_6 * mf1_103[k]
                   + pb_x[k] * mg_153[k];

        t_214[k] = pb_z[k] * mg_151[k];

        t_215[k] = f_3 * mf0_100[k]
                   - f_4 * mf1_100[k]
                   + pb_z[k] * mg_152[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pb_x, pb_y, pb_z, lg_95, lg_156, mf0_102, \
                         mf0_106, mf1_102, mf1_106, mg_153, mg_155, \
                         mg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_11 * lg_156[k]
                   + f_3 * mf0_106[k]
                   - f_4 * mf1_106[k]
                   + pb_x[k] * mg_156[k];

        t_217[k] = pb_z[k] * mg_153[k];

        t_218[k] = f_24 * lg_95[k]
                   + pb_y[k] * mg_155[k];

        t_219[k] = f_5 * mf0_102[k]
                   - f_6 * mf1_102[k]
                   + pb_z[k] * mg_155[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, pb_x, pb_z, lg_160, lg_162, \
                         lg_163, lg_164, mg_156, mg_160, mg_162, mg_163, \
                         mg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_11 * lg_160[k]
                   + pb_x[k] * mg_160[k];

        t_221[k] = pb_z[k] * mg_156[k];

        t_222[k] = f_11 * lg_162[k]
                   + pb_x[k] * mg_162[k];

        t_223[k] = f_11 * lg_163[k]
                   + pb_x[k] * mg_163[k];

        t_224[k] = f_11 * lg_164[k]
                   + pb_x[k] * mg_164[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_x, pb_z, kh0_225, kh1_225, lh_225, \
                         mf0_106, mf0_107, mf1_106, mf1_107, mg_160, mg_161, \
                         mg_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_25 * kh0_225[k]
                   - f_26 * kh1_225[k]
                   + pa_x[k] * lh_225[k];

        t_226[k] = pb_z[k] * mg_160[k];

        t_227[k] = f_3 * mf0_106[k]
                   - f_4 * mf1_106[k]
                   + pb_z[k] * mg_161[k];

        t_228[k] = f_5 * mf0_107[k]
                   - f_6 * mf1_107[k]
                   + pb_z[k] * mg_162[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, pa_z, pb_y, pb_z, lg_90, lg_104, \
                         lh_126, lh_127, mf0_109, mf1_109, mg_164, \
                         mg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_24 * lg_104[k]
                   + pb_y[k] * mg_164[k];

        t_230[k] = f_1 * mf0_109[k]
                   - f_2 * mf1_109[k]
                   + pb_z[k] * mg_164[k];

        t_231[k] = pa_z[k] * lh_126[k];

        t_232[k] = pa_z[k] * lh_127[k];

        t_233[k] = f_7 * lg_90[k]
                   + pb_z[k] * mg_165[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pa_z, pb_y, pb_z, lg_92, lg_93, \
                         lg_107, lh_129, lh_131, lh_132, mg_167, \
                         mg_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_z[k] * lh_129[k];

        t_235[k] = f_9 * lg_107[k]
                   + pb_y[k] * mg_167[k];

        t_236[k] = f_8 * lg_92[k]
                   + pa_z[k] * lh_131[k];

        t_237[k] = pa_z[k] * lh_132[k];

        t_238[k] = f_7 * lg_93[k]
                   + pb_z[k] * mg_168[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pa_z, pb_x, pb_y, lg_95, lg_110, lg_176, \
                         lh_135, lh_136, mg_170, mg_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_9 * lg_110[k]
                   + pb_y[k] * mg_170[k];

        t_240[k] = f_9 * lg_95[k]
                   + pa_z[k] * lh_135[k];

        t_241[k] = pa_z[k] * lh_136[k];

        t_242[k] = f_11 * lg_176[k]
                   + pb_x[k] * mg_176[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pa_z, pb_x, lg_177, lg_178, lg_179, \
                         lh_141, mg_177, mg_178, mg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_11 * lg_177[k]
                   + pb_x[k] * mg_177[k];

        t_244[k] = f_11 * lg_178[k]
                   + pb_x[k] * mg_178[k];

        t_245[k] = f_11 * lg_179[k]
                   + pb_x[k] * mg_179[k];

        t_246[k] = pa_z[k] * lh_141[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pa_z, pb_y, pb_z, lg_100, lg_101, lg_102, \
                         lg_119, lh_143, lh_144, mg_175, mg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_7 * lg_100[k]
                   + pb_z[k] * mg_175[k];

        t_248[k] = f_8 * lg_101[k]
                   + pa_z[k] * lh_143[k];

        t_249[k] = f_9 * lg_102[k]
                   + pa_z[k] * lh_144[k];

        t_250[k] = f_9 * lg_119[k]
                   + pb_y[k] * mg_179[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_y, pa_z, pb_y, pb_z, kh0_105, kh1_105, \
                         lg_104, lg_105, lg_120, lh_146, lh_168, \
                         mg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_11 * lg_104[k]
                   + pa_z[k] * lh_146[k];

        t_252[k] = f_12 * kh0_105[k]
                   - f_13 * kh1_105[k]
                   + pa_y[k] * lh_168[k];

        t_253[k] = f_8 * lg_120[k]
                   + pb_y[k] * mg_180[k];

        t_254[k] = f_8 * lg_105[k]
                   + pb_z[k] * mg_180[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pa_y, pa_z, pb_y, kh0_66, kh0_110, kh1_66, \
                         kh1_110, lg_122, lh_150, lh_173, mg_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_12 * kh0_66[k]
                   - f_13 * kh1_66[k]
                   + pa_z[k] * lh_150[k];

        t_256[k] = f_8 * lg_122[k]
                   + pb_y[k] * mg_182[k];

        t_257[k] = f_12 * kh0_110[k]
                   - f_13 * kh1_110[k]
                   + pa_y[k] * lh_173[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pa_z, pb_y, pb_z, kh0_69, kh1_69, lg_108, \
                         lg_125, lh_153, mg_183, mg_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_12 * kh0_69[k]
                   - f_13 * kh1_69[k]
                   + pa_z[k] * lh_153[k];

        t_259[k] = f_8 * lg_108[k]
                   + pb_z[k] * mg_183[k];

        t_260[k] = f_8 * lg_125[k]
                   + pb_y[k] * mg_185[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pa_y, pb_x, kh0_114, kh1_114, lg_190, \
                         lg_191, lg_192, lh_177, mg_190, mg_191, \
                         mg_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_12 * kh0_114[k]
                   - f_13 * kh1_114[k]
                   + pa_y[k] * lh_177[k];

        t_262[k] = f_11 * lg_190[k]
                   + pb_x[k] * mg_190[k];

        t_263[k] = f_11 * lg_191[k]
                   + pb_x[k] * mg_191[k];

        t_264[k] = f_11 * lg_192[k]
                   + pb_x[k] * mg_192[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pa_x, pb_x, pb_z, kh0_267, kh1_267, \
                         lg_115, lg_193, lg_194, lh_267, mg_190, mg_193, \
                         mg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_11 * lg_193[k]
                   + pb_x[k] * mg_193[k];

        t_266[k] = f_11 * lg_194[k]
                   + pb_x[k] * mg_194[k];

        t_267[k] = f_25 * kh0_267[k]
                   - f_26 * kh1_267[k]
                   + pa_x[k] * lh_267[k];

        t_268[k] = f_8 * lg_115[k]
                   + pb_z[k] * mg_190[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pa_x, pb_y, kh0_269, kh0_270, kh1_269, kh1_270, \
                         lg_134, lh_269, lh_270, mg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_25 * kh0_269[k]
                   - f_26 * kh1_269[k]
                   + pa_x[k] * lh_269[k];

        t_270[k] = f_25 * kh0_270[k]
                   - f_26 * kh1_270[k]
                   + pa_x[k] * lh_270[k];

        t_271[k] = f_8 * lg_134[k]
                   + pb_y[k] * mg_194[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pa_x, pa_y, pb_y, kh0_272, kh1_272, \
                         lg_135, lh_189, lh_191, lh_272, mg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_25 * kh0_272[k]
                   - f_26 * kh1_272[k]
                   + pa_x[k] * lh_272[k];

        t_273[k] = pa_y[k] * lh_189[k];

        t_274[k] = f_7 * lg_135[k]
                   + pb_y[k] * mg_195[k];

        t_275[k] = pa_y[k] * lh_191[k];
    }
}

static auto
compute_prim_mh_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kh0,
                                            const size_t kh1, const size_t lg, const size_t lh,
                                            const size_t mf0, const size_t mf1, const size_t mg,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_11 = 2.5 / p;
    const auto f_12 = 0.5 / alpha;
    const auto f_13 = 0.5 * beta / (alpha * p);
    const auto f_17 = 1.0 / alpha;
    const auto f_18 = beta / (alpha * p);
    const auto f_22 = 1.5 / alpha;
    const auto f_23 = 1.5 * beta / (alpha * p);
    const auto f_24 = 2.0 / p;
    const auto f_25 = 2.0 / alpha;
    const auto f_26 = 2.0 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kh0_105 = buffer.data(kh0 + 105);
    const auto *kh0_126 = buffer.data(kh0 + 126);
    const auto *kh0_129 = buffer.data(kh0 + 129);
    const auto *kh0_132 = buffer.data(kh0 + 132);
    const auto *kh0_150 = buffer.data(kh0 + 150);
    const auto *kh0_153 = buffer.data(kh0 + 153);
    const auto *kh0_168 = buffer.data(kh0 + 168);
    const auto *kh0_173 = buffer.data(kh0 + 173);
    const auto *kh0_177 = buffer.data(kh0 + 177);
    const auto *kh0_189 = buffer.data(kh0 + 189);
    const auto *kh0_194 = buffer.data(kh0 + 194);
    const auto *kh0_198 = buffer.data(kh0 + 198);
    const auto *kh0_314 = buffer.data(kh0 + 314);
    const auto *kh0_330 = buffer.data(kh0 + 330);
    const auto *kh0_372 = buffer.data(kh0 + 372);
    const auto *kh0_374 = buffer.data(kh0 + 374);
    const auto *kh0_375 = buffer.data(kh0 + 375);
    const auto *kh0_377 = buffer.data(kh0 + 377);
    const auto *kh0_393 = buffer.data(kh0 + 393);
    const auto *kh0_395 = buffer.data(kh0 + 395);
    const auto *kh0_396 = buffer.data(kh0 + 396);
    const auto *kh0_398 = buffer.data(kh0 + 398);

    const auto *kh1_105 = buffer.data(kh1 + 105);
    const auto *kh1_126 = buffer.data(kh1 + 126);
    const auto *kh1_129 = buffer.data(kh1 + 129);
    const auto *kh1_132 = buffer.data(kh1 + 132);
    const auto *kh1_150 = buffer.data(kh1 + 150);
    const auto *kh1_153 = buffer.data(kh1 + 153);
    const auto *kh1_168 = buffer.data(kh1 + 168);
    const auto *kh1_173 = buffer.data(kh1 + 173);
    const auto *kh1_177 = buffer.data(kh1 + 177);
    const auto *kh1_189 = buffer.data(kh1 + 189);
    const auto *kh1_194 = buffer.data(kh1 + 194);
    const auto *kh1_198 = buffer.data(kh1 + 198);
    const auto *kh1_314 = buffer.data(kh1 + 314);
    const auto *kh1_330 = buffer.data(kh1 + 330);
    const auto *kh1_372 = buffer.data(kh1 + 372);
    const auto *kh1_374 = buffer.data(kh1 + 374);
    const auto *kh1_375 = buffer.data(kh1 + 375);
    const auto *kh1_377 = buffer.data(kh1 + 377);
    const auto *kh1_393 = buffer.data(kh1 + 393);
    const auto *kh1_395 = buffer.data(kh1 + 395);
    const auto *kh1_396 = buffer.data(kh1 + 396);
    const auto *kh1_398 = buffer.data(kh1 + 398);

    const auto *lg_123 = buffer.data(lg + 123);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_136 = buffer.data(lg + 136);
    const auto *lg_137 = buffer.data(lg + 137);
    const auto *lg_138 = buffer.data(lg + 138);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_148 = buffer.data(lg + 148);
    const auto *lg_149 = buffer.data(lg + 149);
    const auto *lg_150 = buffer.data(lg + 150);
    const auto *lg_152 = buffer.data(lg + 152);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_155 = buffer.data(lg + 155);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_161 = buffer.data(lg + 161);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_165 = buffer.data(lg + 165);
    const auto *lg_167 = buffer.data(lg + 167);
    const auto *lg_168 = buffer.data(lg + 168);
    const auto *lg_170 = buffer.data(lg + 170);
    const auto *lg_175 = buffer.data(lg + 175);
    const auto *lg_179 = buffer.data(lg + 179);
    const auto *lg_180 = buffer.data(lg + 180);
    const auto *lg_182 = buffer.data(lg + 182);
    const auto *lg_183 = buffer.data(lg + 183);
    const auto *lg_185 = buffer.data(lg + 185);
    const auto *lg_190 = buffer.data(lg + 190);
    const auto *lg_194 = buffer.data(lg + 194);
    const auto *lg_195 = buffer.data(lg + 195);
    const auto *lg_197 = buffer.data(lg + 197);
    const auto *lg_200 = buffer.data(lg + 200);
    const auto *lg_205 = buffer.data(lg + 205);
    const auto *lg_206 = buffer.data(lg + 206);
    const auto *lg_207 = buffer.data(lg + 207);
    const auto *lg_208 = buffer.data(lg + 208);
    const auto *lg_209 = buffer.data(lg + 209);
    const auto *lg_210 = buffer.data(lg + 210);
    const auto *lg_211 = buffer.data(lg + 211);
    const auto *lg_212 = buffer.data(lg + 212);
    const auto *lg_213 = buffer.data(lg + 213);
    const auto *lg_215 = buffer.data(lg + 215);
    const auto *lg_219 = buffer.data(lg + 219);
    const auto *lg_220 = buffer.data(lg + 220);
    const auto *lg_221 = buffer.data(lg + 221);
    const auto *lg_222 = buffer.data(lg + 222);
    const auto *lg_224 = buffer.data(lg + 224);
    const auto *lg_228 = buffer.data(lg + 228);
    const auto *lg_231 = buffer.data(lg + 231);
    const auto *lg_235 = buffer.data(lg + 235);
    const auto *lg_237 = buffer.data(lg + 237);
    const auto *lg_238 = buffer.data(lg + 238);
    const auto *lg_239 = buffer.data(lg + 239);
    const auto *lg_251 = buffer.data(lg + 251);
    const auto *lg_252 = buffer.data(lg + 252);
    const auto *lg_253 = buffer.data(lg + 253);
    const auto *lg_254 = buffer.data(lg + 254);
    const auto *lg_265 = buffer.data(lg + 265);
    const auto *lg_266 = buffer.data(lg + 266);
    const auto *lg_267 = buffer.data(lg + 267);
    const auto *lg_268 = buffer.data(lg + 268);
    const auto *lg_269 = buffer.data(lg + 269);
    const auto *lg_280 = buffer.data(lg + 280);
    const auto *lg_281 = buffer.data(lg + 281);
    const auto *lg_282 = buffer.data(lg + 282);
    const auto *lg_283 = buffer.data(lg + 283);
    const auto *lg_284 = buffer.data(lg + 284);

    const auto *lh_189 = buffer.data(lh + 189);
    const auto *lh_192 = buffer.data(lh + 192);
    const auto *lh_194 = buffer.data(lh + 194);
    const auto *lh_195 = buffer.data(lh + 195);
    const auto *lh_198 = buffer.data(lh + 198);
    const auto *lh_203 = buffer.data(lh + 203);
    const auto *lh_204 = buffer.data(lh + 204);
    const auto *lh_206 = buffer.data(lh + 206);
    const auto *lh_207 = buffer.data(lh + 207);
    const auto *lh_209 = buffer.data(lh + 209);
    const auto *lh_210 = buffer.data(lh + 210);
    const auto *lh_211 = buffer.data(lh + 211);
    const auto *lh_213 = buffer.data(lh + 213);
    const auto *lh_215 = buffer.data(lh + 215);
    const auto *lh_216 = buffer.data(lh + 216);
    const auto *lh_219 = buffer.data(lh + 219);
    const auto *lh_220 = buffer.data(lh + 220);
    const auto *lh_225 = buffer.data(lh + 225);
    const auto *lh_227 = buffer.data(lh + 227);
    const auto *lh_228 = buffer.data(lh + 228);
    const auto *lh_230 = buffer.data(lh + 230);
    const auto *lh_234 = buffer.data(lh + 234);
    const auto *lh_237 = buffer.data(lh + 237);
    const auto *lh_252 = buffer.data(lh + 252);
    const auto *lh_255 = buffer.data(lh + 255);
    const auto *lh_257 = buffer.data(lh + 257);
    const auto *lh_258 = buffer.data(lh + 258);
    const auto *lh_261 = buffer.data(lh + 261);
    const auto *lh_273 = buffer.data(lh + 273);
    const auto *lh_278 = buffer.data(lh + 278);
    const auto *lh_282 = buffer.data(lh + 282);
    const auto *lh_294 = buffer.data(lh + 294);
    const auto *lh_296 = buffer.data(lh + 296);
    const auto *lh_297 = buffer.data(lh + 297);
    const auto *lh_299 = buffer.data(lh + 299);
    const auto *lh_300 = buffer.data(lh + 300);
    const auto *lh_314 = buffer.data(lh + 314);
    const auto *lh_330 = buffer.data(lh + 330);
    const auto *lh_372 = buffer.data(lh + 372);
    const auto *lh_374 = buffer.data(lh + 374);
    const auto *lh_375 = buffer.data(lh + 375);
    const auto *lh_377 = buffer.data(lh + 377);
    const auto *lh_393 = buffer.data(lh + 393);
    const auto *lh_395 = buffer.data(lh + 395);
    const auto *lh_396 = buffer.data(lh + 396);
    const auto *lh_398 = buffer.data(lh + 398);

    const auto *mf0_140 = buffer.data(mf0 + 140);
    const auto *mf0_141 = buffer.data(mf0 + 141);
    const auto *mf0_145 = buffer.data(mf0 + 145);
    const auto *mf0_146 = buffer.data(mf0 + 146);
    const auto *mf0_148 = buffer.data(mf0 + 148);
    const auto *mf0_149 = buffer.data(mf0 + 149);
    const auto *mf0_150 = buffer.data(mf0 + 150);
    const auto *mf0_152 = buffer.data(mf0 + 152);
    const auto *mf0_153 = buffer.data(mf0 + 153);
    const auto *mf0_156 = buffer.data(mf0 + 156);
    const auto *mf0_157 = buffer.data(mf0 + 157);
    const auto *mf0_159 = buffer.data(mf0 + 159);

    const auto *mf1_140 = buffer.data(mf1 + 140);
    const auto *mf1_141 = buffer.data(mf1 + 141);
    const auto *mf1_145 = buffer.data(mf1 + 145);
    const auto *mf1_146 = buffer.data(mf1 + 146);
    const auto *mf1_148 = buffer.data(mf1 + 148);
    const auto *mf1_149 = buffer.data(mf1 + 149);
    const auto *mf1_150 = buffer.data(mf1 + 150);
    const auto *mf1_152 = buffer.data(mf1 + 152);
    const auto *mf1_153 = buffer.data(mf1 + 153);
    const auto *mf1_156 = buffer.data(mf1 + 156);
    const auto *mf1_157 = buffer.data(mf1 + 157);
    const auto *mf1_159 = buffer.data(mf1 + 159);

    const auto *mg_197 = buffer.data(mg + 197);
    const auto *mg_198 = buffer.data(mg + 198);
    const auto *mg_200 = buffer.data(mg + 200);
    const auto *mg_205 = buffer.data(mg + 205);
    const auto *mg_206 = buffer.data(mg + 206);
    const auto *mg_207 = buffer.data(mg + 207);
    const auto *mg_208 = buffer.data(mg + 208);
    const auto *mg_209 = buffer.data(mg + 209);
    const auto *mg_210 = buffer.data(mg + 210);
    const auto *mg_211 = buffer.data(mg + 211);
    const auto *mg_212 = buffer.data(mg + 212);
    const auto *mg_213 = buffer.data(mg + 213);
    const auto *mg_215 = buffer.data(mg + 215);
    const auto *mg_219 = buffer.data(mg + 219);
    const auto *mg_220 = buffer.data(mg + 220);
    const auto *mg_221 = buffer.data(mg + 221);
    const auto *mg_222 = buffer.data(mg + 222);
    const auto *mg_223 = buffer.data(mg + 223);
    const auto *mg_224 = buffer.data(mg + 224);
    const auto *mg_225 = buffer.data(mg + 225);
    const auto *mg_226 = buffer.data(mg + 226);
    const auto *mg_227 = buffer.data(mg + 227);
    const auto *mg_228 = buffer.data(mg + 228);
    const auto *mg_230 = buffer.data(mg + 230);
    const auto *mg_231 = buffer.data(mg + 231);
    const auto *mg_235 = buffer.data(mg + 235);
    const auto *mg_236 = buffer.data(mg + 236);
    const auto *mg_237 = buffer.data(mg + 237);
    const auto *mg_238 = buffer.data(mg + 238);
    const auto *mg_239 = buffer.data(mg + 239);
    const auto *mg_240 = buffer.data(mg + 240);
    const auto *mg_242 = buffer.data(mg + 242);
    const auto *mg_243 = buffer.data(mg + 243);
    const auto *mg_245 = buffer.data(mg + 245);
    const auto *mg_250 = buffer.data(mg + 250);
    const auto *mg_251 = buffer.data(mg + 251);
    const auto *mg_252 = buffer.data(mg + 252);
    const auto *mg_253 = buffer.data(mg + 253);
    const auto *mg_254 = buffer.data(mg + 254);
    const auto *mg_255 = buffer.data(mg + 255);
    const auto *mg_257 = buffer.data(mg + 257);
    const auto *mg_258 = buffer.data(mg + 258);
    const auto *mg_260 = buffer.data(mg + 260);
    const auto *mg_265 = buffer.data(mg + 265);
    const auto *mg_266 = buffer.data(mg + 266);
    const auto *mg_267 = buffer.data(mg + 267);
    const auto *mg_268 = buffer.data(mg + 268);
    const auto *mg_269 = buffer.data(mg + 269);
    const auto *mg_270 = buffer.data(mg + 270);
    const auto *mg_272 = buffer.data(mg + 272);
    const auto *mg_273 = buffer.data(mg + 273);
    const auto *mg_275 = buffer.data(mg + 275);
    const auto *mg_280 = buffer.data(mg + 280);
    const auto *mg_281 = buffer.data(mg + 281);
    const auto *mg_282 = buffer.data(mg + 282);
    const auto *mg_283 = buffer.data(mg + 283);
    const auto *mg_284 = buffer.data(mg + 284);
    const auto *mg_285 = buffer.data(mg + 285);
    const auto *mg_287 = buffer.data(mg + 287);

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pa_y, pb_y, lg_136, lg_137, lg_138, \
                         lh_192, lh_194, lh_195, mg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_8 * lg_136[k]
                   + pa_y[k] * lh_192[k];

        t_277[k] = f_7 * lg_137[k]
                   + pb_y[k] * mg_197[k];

        t_278[k] = pa_y[k] * lh_194[k];

        t_279[k] = f_9 * lg_138[k]
                   + pa_y[k] * lh_195[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pa_y, pb_x, pb_y, pb_z, lg_123, lg_140, \
                         lg_205, lh_198, mg_198, mg_200, mg_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_9 * lg_123[k]
                   + pb_z[k] * mg_198[k];

        t_281[k] = f_7 * lg_140[k]
                   + pb_y[k] * mg_200[k];

        t_282[k] = pa_y[k] * lh_198[k];

        t_283[k] = f_11 * lg_205[k]
                   + pb_x[k] * mg_205[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, pa_y, pb_x, lg_145, lg_206, \
                         lg_207, lg_208, lh_203, lh_204, mg_206, mg_207, \
                         mg_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_11 * lg_206[k]
                   + pb_x[k] * mg_206[k];

        t_285[k] = f_11 * lg_207[k]
                   + pb_x[k] * mg_207[k];

        t_286[k] = f_11 * lg_208[k]
                   + pb_x[k] * mg_208[k];

        t_287[k] = pa_y[k] * lh_203[k];

        t_288[k] = f_11 * lg_145[k]
                   + pa_y[k] * lh_204[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pa_y, pb_y, pb_z, lg_130, lg_147, lg_148, \
                         lg_149, lh_206, lh_207, mg_205, mg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_9 * lg_130[k]
                   + pb_z[k] * mg_205[k];

        t_290[k] = f_9 * lg_147[k]
                   + pa_y[k] * lh_206[k];

        t_291[k] = f_8 * lg_148[k]
                   + pa_y[k] * lh_207[k];

        t_292[k] = f_7 * lg_149[k]
                   + pb_y[k] * mg_209[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_y, pa_z, pb_y, pb_z, kh0_105, kh1_105, \
                         lg_135, lh_189, lh_209, mg_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = pa_y[k] * lh_209[k];

        t_294[k] = f_22 * kh0_105[k]
                   - f_23 * kh1_105[k]
                   + pa_z[k] * lh_189[k];

        t_295[k] = pb_y[k] * mg_210[k];

        t_296[k] = f_24 * lg_135[k]
                   + pb_z[k] * mg_210[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_x, pb_y, lg_215, mf0_140, mf0_145, mf1_140, \
                         mf1_145, mg_211, mg_212, mg_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_3 * mf0_140[k]
                   - f_4 * mf1_140[k]
                   + pb_y[k] * mg_211[k];

        t_298[k] = pb_y[k] * mg_212[k];

        t_299[k] = f_11 * lg_215[k]
                   + f_5 * mf0_145[k]
                   - f_6 * mf1_145[k]
                   + pb_x[k] * mg_215[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pb_x, pb_y, pb_z, lg_138, lg_219, \
                         mf0_141, mf0_149, mf1_141, mf1_149, mg_213, mg_215, \
                         mg_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_5 * mf0_141[k]
                   - f_6 * mf1_141[k]
                   + pb_y[k] * mg_213[k];

        t_301[k] = f_24 * lg_138[k]
                   + pb_z[k] * mg_213[k];

        t_302[k] = pb_y[k] * mg_215[k];

        t_303[k] = f_11 * lg_219[k]
                   + f_3 * mf0_149[k]
                   - f_4 * mf1_149[k]
                   + pb_x[k] * mg_219[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pb_x, pb_y, lg_220, lg_221, \
                         lg_222, lg_224, mg_219, mg_220, mg_221, mg_222, \
                         mg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_11 * lg_220[k]
                   + pb_x[k] * mg_220[k];

        t_305[k] = f_11 * lg_221[k]
                   + pb_x[k] * mg_221[k];

        t_306[k] = f_11 * lg_222[k]
                   + pb_x[k] * mg_222[k];

        t_307[k] = pb_y[k] * mg_219[k];

        t_308[k] = f_11 * lg_224[k]
                   + pb_x[k] * mg_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pb_y, pb_z, lg_145, mf0_146, mf0_148, \
                         mf0_149, mf1_146, mf1_148, mf1_149, mg_220, mg_222, \
                         mg_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_1 * mf0_146[k]
                   - f_2 * mf1_146[k]
                   + pb_y[k] * mg_220[k];

        t_310[k] = f_24 * lg_145[k]
                   + pb_z[k] * mg_220[k];

        t_311[k] = f_5 * mf0_148[k]
                   - f_6 * mf1_148[k]
                   + pb_y[k] * mg_222[k];

        t_312[k] = f_3 * mf0_149[k]
                   - f_4 * mf1_149[k]
                   + pb_y[k] * mg_223[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pa_x, pa_y, pb_y, kh0_126, kh0_314, \
                         kh1_126, kh1_314, lg_150, lh_210, lh_314, mg_224, \
                         mg_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = pb_y[k] * mg_224[k];

        t_314[k] = f_25 * kh0_314[k]
                   - f_26 * kh1_314[k]
                   + pa_x[k] * lh_314[k];

        t_315[k] = f_25 * kh0_126[k]
                   - f_26 * kh1_126[k]
                   + pa_y[k] * lh_210[k];

        t_316[k] = f_11 * lg_150[k]
                   + pb_y[k] * mg_225[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, pb_x, pb_z, lg_228, mf0_150, mf0_153, \
                         mf1_150, mf1_153, mg_225, mg_226, mg_227, \
                         mg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = pb_z[k] * mg_225[k];

        t_318[k] = f_24 * lg_228[k]
                   + f_5 * mf0_153[k]
                   - f_6 * mf1_153[k]
                   + pb_x[k] * mg_228[k];

        t_319[k] = pb_z[k] * mg_226[k];

        t_320[k] = f_3 * mf0_150[k]
                   - f_4 * mf1_150[k]
                   + pb_z[k] * mg_227[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pb_x, pb_y, pb_z, lg_155, lg_231, \
                         mf0_152, mf0_156, mf1_152, mf1_156, mg_228, mg_230, \
                         mg_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_24 * lg_231[k]
                   + f_3 * mf0_156[k]
                   - f_4 * mf1_156[k]
                   + pb_x[k] * mg_231[k];

        t_322[k] = pb_z[k] * mg_228[k];

        t_323[k] = f_11 * lg_155[k]
                   + pb_y[k] * mg_230[k];

        t_324[k] = f_5 * mf0_152[k]
                   - f_6 * mf1_152[k]
                   + pb_z[k] * mg_230[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, pb_x, pb_z, lg_235, lg_237, \
                         lg_238, lg_239, mg_231, mg_235, mg_237, mg_238, \
                         mg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_24 * lg_235[k]
                   + pb_x[k] * mg_235[k];

        t_326[k] = pb_z[k] * mg_231[k];

        t_327[k] = f_24 * lg_237[k]
                   + pb_x[k] * mg_237[k];

        t_328[k] = f_24 * lg_238[k]
                   + pb_x[k] * mg_238[k];

        t_329[k] = f_24 * lg_239[k]
                   + pb_x[k] * mg_239[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pa_x, pb_z, kh0_330, kh1_330, lh_330, \
                         mf0_156, mf0_157, mf1_156, mf1_157, mg_235, mg_236, \
                         mg_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_22 * kh0_330[k]
                   - f_23 * kh1_330[k]
                   + pa_x[k] * lh_330[k];

        t_331[k] = pb_z[k] * mg_235[k];

        t_332[k] = f_3 * mf0_156[k]
                   - f_4 * mf1_156[k]
                   + pb_z[k] * mg_236[k];

        t_333[k] = f_5 * mf0_157[k]
                   - f_6 * mf1_157[k]
                   + pb_z[k] * mg_237[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, pa_z, pb_y, pb_z, lg_150, lg_164, \
                         lh_210, lh_211, mf0_159, mf1_159, mg_239, \
                         mg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_11 * lg_164[k]
                   + pb_y[k] * mg_239[k];

        t_335[k] = f_1 * mf0_159[k]
                   - f_2 * mf1_159[k]
                   + pb_z[k] * mg_239[k];

        t_336[k] = pa_z[k] * lh_210[k];

        t_337[k] = pa_z[k] * lh_211[k];

        t_338[k] = f_7 * lg_150[k]
                   + pb_z[k] * mg_240[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, t_343, pa_z, pb_y, pb_z, lg_152, lg_153, \
                         lg_167, lh_213, lh_215, lh_216, mg_242, \
                         mg_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pa_z[k] * lh_213[k];

        t_340[k] = f_24 * lg_167[k]
                   + pb_y[k] * mg_242[k];

        t_341[k] = f_8 * lg_152[k]
                   + pa_z[k] * lh_215[k];

        t_342[k] = pa_z[k] * lh_216[k];

        t_343[k] = f_7 * lg_153[k]
                   + pb_z[k] * mg_243[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pa_z, pb_x, pb_y, lg_155, lg_170, lg_251, \
                         lh_219, lh_220, mg_245, mg_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_24 * lg_170[k]
                   + pb_y[k] * mg_245[k];

        t_345[k] = f_9 * lg_155[k]
                   + pa_z[k] * lh_219[k];

        t_346[k] = pa_z[k] * lh_220[k];

        t_347[k] = f_24 * lg_251[k]
                   + pb_x[k] * mg_251[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pa_z, pb_x, lg_252, lg_253, lg_254, \
                         lh_225, mg_252, mg_253, mg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_24 * lg_252[k]
                   + pb_x[k] * mg_252[k];

        t_349[k] = f_24 * lg_253[k]
                   + pb_x[k] * mg_253[k];

        t_350[k] = f_24 * lg_254[k]
                   + pb_x[k] * mg_254[k];

        t_351[k] = pa_z[k] * lh_225[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pa_z, pb_y, pb_z, lg_160, lg_161, lg_162, \
                         lg_179, lh_227, lh_228, mg_250, mg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_7 * lg_160[k]
                   + pb_z[k] * mg_250[k];

        t_353[k] = f_8 * lg_161[k]
                   + pa_z[k] * lh_227[k];

        t_354[k] = f_9 * lg_162[k]
                   + pa_z[k] * lh_228[k];

        t_355[k] = f_24 * lg_179[k]
                   + pb_y[k] * mg_254[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pa_y, pa_z, pb_y, pb_z, kh0_168, kh1_168, \
                         lg_164, lg_165, lg_180, lh_230, lh_252, \
                         mg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_11 * lg_164[k]
                   + pa_z[k] * lh_230[k];

        t_357[k] = f_17 * kh0_168[k]
                   - f_18 * kh1_168[k]
                   + pa_y[k] * lh_252[k];

        t_358[k] = f_9 * lg_180[k]
                   + pb_y[k] * mg_255[k];

        t_359[k] = f_8 * lg_165[k]
                   + pb_z[k] * mg_255[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pa_y, pa_z, pb_y, kh0_129, kh0_173, kh1_129, \
                         kh1_173, lg_182, lh_234, lh_257, mg_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_12 * kh0_129[k]
                   - f_13 * kh1_129[k]
                   + pa_z[k] * lh_234[k];

        t_361[k] = f_9 * lg_182[k]
                   + pb_y[k] * mg_257[k];

        t_362[k] = f_17 * kh0_173[k]
                   - f_18 * kh1_173[k]
                   + pa_y[k] * lh_257[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pa_z, pb_y, pb_z, kh0_132, kh1_132, lg_168, \
                         lg_185, lh_237, mg_258, mg_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_12 * kh0_132[k]
                   - f_13 * kh1_132[k]
                   + pa_z[k] * lh_237[k];

        t_364[k] = f_8 * lg_168[k]
                   + pb_z[k] * mg_258[k];

        t_365[k] = f_9 * lg_185[k]
                   + pb_y[k] * mg_260[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_y, pb_x, kh0_177, kh1_177, lg_265, \
                         lg_266, lg_267, lh_261, mg_265, mg_266, \
                         mg_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_17 * kh0_177[k]
                   - f_18 * kh1_177[k]
                   + pa_y[k] * lh_261[k];

        t_367[k] = f_24 * lg_265[k]
                   + pb_x[k] * mg_265[k];

        t_368[k] = f_24 * lg_266[k]
                   + pb_x[k] * mg_266[k];

        t_369[k] = f_24 * lg_267[k]
                   + pb_x[k] * mg_267[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_x, pb_x, pb_z, kh0_372, kh1_372, \
                         lg_175, lg_268, lg_269, lh_372, mg_265, mg_268, \
                         mg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_24 * lg_268[k]
                   + pb_x[k] * mg_268[k];

        t_371[k] = f_24 * lg_269[k]
                   + pb_x[k] * mg_269[k];

        t_372[k] = f_22 * kh0_372[k]
                   - f_23 * kh1_372[k]
                   + pa_x[k] * lh_372[k];

        t_373[k] = f_8 * lg_175[k]
                   + pb_z[k] * mg_265[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pa_x, pb_y, kh0_374, kh0_375, kh1_374, kh1_375, \
                         lg_194, lh_374, lh_375, mg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = f_22 * kh0_374[k]
                   - f_23 * kh1_374[k]
                   + pa_x[k] * lh_374[k];

        t_375[k] = f_22 * kh0_375[k]
                   - f_23 * kh1_375[k]
                   + pa_x[k] * lh_375[k];

        t_376[k] = f_9 * lg_194[k]
                   + pb_y[k] * mg_269[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, pa_x, pa_y, pb_y, kh0_189, kh0_377, kh1_189, \
                         kh1_377, lg_195, lh_273, lh_377, mg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_22 * kh0_377[k]
                   - f_23 * kh1_377[k]
                   + pa_x[k] * lh_377[k];

        t_378[k] = f_12 * kh0_189[k]
                   - f_13 * kh1_189[k]
                   + pa_y[k] * lh_273[k];

        t_379[k] = f_8 * lg_195[k]
                   + pb_y[k] * mg_270[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, pa_z, pb_y, pb_z, kh0_150, kh1_150, lg_180, \
                         lg_197, lh_255, mg_270, mg_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_9 * lg_180[k]
                   + pb_z[k] * mg_270[k];

        t_381[k] = f_17 * kh0_150[k]
                   - f_18 * kh1_150[k]
                   + pa_z[k] * lh_255[k];

        t_382[k] = f_8 * lg_197[k]
                   + pb_y[k] * mg_272[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pa_y, pa_z, pb_z, kh0_153, kh0_194, kh1_153, \
                         kh1_194, lg_183, lh_258, lh_278, mg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_12 * kh0_194[k]
                   - f_13 * kh1_194[k]
                   + pa_y[k] * lh_278[k];

        t_384[k] = f_17 * kh0_153[k]
                   - f_18 * kh1_153[k]
                   + pa_z[k] * lh_258[k];

        t_385[k] = f_9 * lg_183[k]
                   + pb_z[k] * mg_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pa_y, pb_x, pb_y, kh0_198, kh1_198, \
                         lg_200, lg_280, lg_281, lh_282, mg_275, mg_280, \
                         mg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_8 * lg_200[k]
                   + pb_y[k] * mg_275[k];

        t_387[k] = f_12 * kh0_198[k]
                   - f_13 * kh1_198[k]
                   + pa_y[k] * lh_282[k];

        t_388[k] = f_24 * lg_280[k]
                   + pb_x[k] * mg_280[k];

        t_389[k] = f_24 * lg_281[k]
                   + pb_x[k] * mg_281[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pa_x, pb_x, kh0_393, kh1_393, lg_282, \
                         lg_283, lg_284, lh_393, mg_282, mg_283, \
                         mg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_24 * lg_282[k]
                   + pb_x[k] * mg_282[k];

        t_391[k] = f_24 * lg_283[k]
                   + pb_x[k] * mg_283[k];

        t_392[k] = f_24 * lg_284[k]
                   + pb_x[k] * mg_284[k];

        t_393[k] = f_22 * kh0_393[k]
                   - f_23 * kh1_393[k]
                   + pa_x[k] * lh_393[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pa_x, pb_z, kh0_395, kh0_396, kh1_395, kh1_396, \
                         lg_190, lh_395, lh_396, mg_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_9 * lg_190[k]
                   + pb_z[k] * mg_280[k];

        t_395[k] = f_22 * kh0_395[k]
                   - f_23 * kh1_395[k]
                   + pa_x[k] * lh_395[k];

        t_396[k] = f_22 * kh0_396[k]
                   - f_23 * kh1_396[k]
                   + pa_x[k] * lh_396[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pa_x, pa_y, pb_y, kh0_398, kh1_398, \
                         lg_209, lg_210, lh_294, lh_398, mg_284, \
                         mg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_8 * lg_209[k]
                   + pb_y[k] * mg_284[k];

        t_398[k] = f_22 * kh0_398[k]
                   - f_23 * kh1_398[k]
                   + pa_x[k] * lh_398[k];

        t_399[k] = pa_y[k] * lh_294[k];

        t_400[k] = f_7 * lg_210[k]
                   + pb_y[k] * mg_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, t_405, pa_y, pb_y, lg_211, lg_212, \
                         lg_213, lh_296, lh_297, lh_299, lh_300, \
                         mg_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = pa_y[k] * lh_296[k];

        t_402[k] = f_8 * lg_211[k]
                   + pa_y[k] * lh_297[k];

        t_403[k] = f_7 * lg_212[k]
                   + pb_y[k] * mg_287[k];

        t_404[k] = pa_y[k] * lh_299[k];

        t_405[k] = f_9 * lg_213[k]
                   + pa_y[k] * lh_300[k];
    }
}

static auto
compute_prim_mh_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kh0,
                                            const size_t kh1, const size_t lg, const size_t lh,
                                            const size_t mf0, const size_t mf1, const size_t mg,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_11 = 2.5 / p;
    const auto f_12 = 0.5 / alpha;
    const auto f_13 = 0.5 * beta / (alpha * p);
    const auto f_17 = 1.0 / alpha;
    const auto f_18 = beta / (alpha * p);
    const auto f_19 = 3.0 / p;
    const auto f_20 = 2.5 / alpha;
    const auto f_21 = 2.5 * beta / (alpha * p);
    const auto f_22 = 1.5 / alpha;
    const auto f_23 = 1.5 * beta / (alpha * p);
    const auto f_24 = 2.0 / p;
    const auto f_25 = 2.0 / alpha;
    const auto f_26 = 2.0 * beta / (alpha * p);

    auto *t_406 = buffer.data(target + 406);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kh0_189 = buffer.data(kh0 + 189);
    const auto *kh0_210 = buffer.data(kh0 + 210);
    const auto *kh0_213 = buffer.data(kh0 + 213);
    const auto *kh0_216 = buffer.data(kh0 + 216);
    const auto *kh0_234 = buffer.data(kh0 + 234);
    const auto *kh0_237 = buffer.data(kh0 + 237);
    const auto *kh0_252 = buffer.data(kh0 + 252);
    const auto *kh0_255 = buffer.data(kh0 + 255);
    const auto *kh0_257 = buffer.data(kh0 + 257);
    const auto *kh0_258 = buffer.data(kh0 + 258);
    const auto *kh0_261 = buffer.data(kh0 + 261);
    const auto *kh0_273 = buffer.data(kh0 + 273);
    const auto *kh0_278 = buffer.data(kh0 + 278);
    const auto *kh0_282 = buffer.data(kh0 + 282);
    const auto *kh0_294 = buffer.data(kh0 + 294);
    const auto *kh0_299 = buffer.data(kh0 + 299);
    const auto *kh0_440 = buffer.data(kh0 + 440);
    const auto *kh0_456 = buffer.data(kh0 + 456);
    const auto *kh0_498 = buffer.data(kh0 + 498);
    const auto *kh0_500 = buffer.data(kh0 + 500);
    const auto *kh0_501 = buffer.data(kh0 + 501);
    const auto *kh0_503 = buffer.data(kh0 + 503);
    const auto *kh0_519 = buffer.data(kh0 + 519);
    const auto *kh0_521 = buffer.data(kh0 + 521);
    const auto *kh0_522 = buffer.data(kh0 + 522);
    const auto *kh0_524 = buffer.data(kh0 + 524);

    const auto *kh1_189 = buffer.data(kh1 + 189);
    const auto *kh1_210 = buffer.data(kh1 + 210);
    const auto *kh1_213 = buffer.data(kh1 + 213);
    const auto *kh1_216 = buffer.data(kh1 + 216);
    const auto *kh1_234 = buffer.data(kh1 + 234);
    const auto *kh1_237 = buffer.data(kh1 + 237);
    const auto *kh1_252 = buffer.data(kh1 + 252);
    const auto *kh1_255 = buffer.data(kh1 + 255);
    const auto *kh1_257 = buffer.data(kh1 + 257);
    const auto *kh1_258 = buffer.data(kh1 + 258);
    const auto *kh1_261 = buffer.data(kh1 + 261);
    const auto *kh1_273 = buffer.data(kh1 + 273);
    const auto *kh1_278 = buffer.data(kh1 + 278);
    const auto *kh1_282 = buffer.data(kh1 + 282);
    const auto *kh1_294 = buffer.data(kh1 + 294);
    const auto *kh1_299 = buffer.data(kh1 + 299);
    const auto *kh1_440 = buffer.data(kh1 + 440);
    const auto *kh1_456 = buffer.data(kh1 + 456);
    const auto *kh1_498 = buffer.data(kh1 + 498);
    const auto *kh1_500 = buffer.data(kh1 + 500);
    const auto *kh1_501 = buffer.data(kh1 + 501);
    const auto *kh1_503 = buffer.data(kh1 + 503);
    const auto *kh1_519 = buffer.data(kh1 + 519);
    const auto *kh1_521 = buffer.data(kh1 + 521);
    const auto *kh1_522 = buffer.data(kh1 + 522);
    const auto *kh1_524 = buffer.data(kh1 + 524);

    const auto *lg_198 = buffer.data(lg + 198);
    const auto *lg_205 = buffer.data(lg + 205);
    const auto *lg_210 = buffer.data(lg + 210);
    const auto *lg_213 = buffer.data(lg + 213);
    const auto *lg_215 = buffer.data(lg + 215);
    const auto *lg_220 = buffer.data(lg + 220);
    const auto *lg_222 = buffer.data(lg + 222);
    const auto *lg_223 = buffer.data(lg + 223);
    const auto *lg_224 = buffer.data(lg + 224);
    const auto *lg_225 = buffer.data(lg + 225);
    const auto *lg_227 = buffer.data(lg + 227);
    const auto *lg_228 = buffer.data(lg + 228);
    const auto *lg_230 = buffer.data(lg + 230);
    const auto *lg_235 = buffer.data(lg + 235);
    const auto *lg_236 = buffer.data(lg + 236);
    const auto *lg_237 = buffer.data(lg + 237);
    const auto *lg_239 = buffer.data(lg + 239);
    const auto *lg_240 = buffer.data(lg + 240);
    const auto *lg_242 = buffer.data(lg + 242);
    const auto *lg_243 = buffer.data(lg + 243);
    const auto *lg_245 = buffer.data(lg + 245);
    const auto *lg_250 = buffer.data(lg + 250);
    const auto *lg_254 = buffer.data(lg + 254);
    const auto *lg_255 = buffer.data(lg + 255);
    const auto *lg_257 = buffer.data(lg + 257);
    const auto *lg_258 = buffer.data(lg + 258);
    const auto *lg_260 = buffer.data(lg + 260);
    const auto *lg_265 = buffer.data(lg + 265);
    const auto *lg_269 = buffer.data(lg + 269);
    const auto *lg_270 = buffer.data(lg + 270);
    const auto *lg_272 = buffer.data(lg + 272);
    const auto *lg_273 = buffer.data(lg + 273);
    const auto *lg_275 = buffer.data(lg + 275);
    const auto *lg_284 = buffer.data(lg + 284);
    const auto *lg_285 = buffer.data(lg + 285);
    const auto *lg_287 = buffer.data(lg + 287);
    const auto *lg_295 = buffer.data(lg + 295);
    const auto *lg_296 = buffer.data(lg + 296);
    const auto *lg_297 = buffer.data(lg + 297);
    const auto *lg_298 = buffer.data(lg + 298);
    const auto *lg_305 = buffer.data(lg + 305);
    const auto *lg_309 = buffer.data(lg + 309);
    const auto *lg_310 = buffer.data(lg + 310);
    const auto *lg_311 = buffer.data(lg + 311);
    const auto *lg_312 = buffer.data(lg + 312);
    const auto *lg_314 = buffer.data(lg + 314);
    const auto *lg_318 = buffer.data(lg + 318);
    const auto *lg_321 = buffer.data(lg + 321);
    const auto *lg_325 = buffer.data(lg + 325);
    const auto *lg_327 = buffer.data(lg + 327);
    const auto *lg_328 = buffer.data(lg + 328);
    const auto *lg_329 = buffer.data(lg + 329);
    const auto *lg_341 = buffer.data(lg + 341);
    const auto *lg_342 = buffer.data(lg + 342);
    const auto *lg_343 = buffer.data(lg + 343);
    const auto *lg_344 = buffer.data(lg + 344);
    const auto *lg_355 = buffer.data(lg + 355);
    const auto *lg_356 = buffer.data(lg + 356);
    const auto *lg_357 = buffer.data(lg + 357);
    const auto *lg_358 = buffer.data(lg + 358);
    const auto *lg_359 = buffer.data(lg + 359);
    const auto *lg_370 = buffer.data(lg + 370);
    const auto *lg_371 = buffer.data(lg + 371);
    const auto *lg_372 = buffer.data(lg + 372);
    const auto *lg_373 = buffer.data(lg + 373);
    const auto *lg_374 = buffer.data(lg + 374);

    const auto *lh_294 = buffer.data(lh + 294);
    const auto *lh_303 = buffer.data(lh + 303);
    const auto *lh_308 = buffer.data(lh + 308);
    const auto *lh_309 = buffer.data(lh + 309);
    const auto *lh_311 = buffer.data(lh + 311);
    const auto *lh_312 = buffer.data(lh + 312);
    const auto *lh_314 = buffer.data(lh + 314);
    const auto *lh_315 = buffer.data(lh + 315);
    const auto *lh_316 = buffer.data(lh + 316);
    const auto *lh_318 = buffer.data(lh + 318);
    const auto *lh_320 = buffer.data(lh + 320);
    const auto *lh_321 = buffer.data(lh + 321);
    const auto *lh_324 = buffer.data(lh + 324);
    const auto *lh_325 = buffer.data(lh + 325);
    const auto *lh_330 = buffer.data(lh + 330);
    const auto *lh_332 = buffer.data(lh + 332);
    const auto *lh_333 = buffer.data(lh + 333);
    const auto *lh_335 = buffer.data(lh + 335);
    const auto *lh_339 = buffer.data(lh + 339);
    const auto *lh_342 = buffer.data(lh + 342);
    const auto *lh_357 = buffer.data(lh + 357);
    const auto *lh_360 = buffer.data(lh + 360);
    const auto *lh_362 = buffer.data(lh + 362);
    const auto *lh_363 = buffer.data(lh + 363);
    const auto *lh_366 = buffer.data(lh + 366);
    const auto *lh_378 = buffer.data(lh + 378);
    const auto *lh_381 = buffer.data(lh + 381);
    const auto *lh_383 = buffer.data(lh + 383);
    const auto *lh_384 = buffer.data(lh + 384);
    const auto *lh_387 = buffer.data(lh + 387);
    const auto *lh_399 = buffer.data(lh + 399);
    const auto *lh_404 = buffer.data(lh + 404);
    const auto *lh_440 = buffer.data(lh + 440);
    const auto *lh_456 = buffer.data(lh + 456);
    const auto *lh_498 = buffer.data(lh + 498);
    const auto *lh_500 = buffer.data(lh + 500);
    const auto *lh_501 = buffer.data(lh + 501);
    const auto *lh_503 = buffer.data(lh + 503);
    const auto *lh_519 = buffer.data(lh + 519);
    const auto *lh_521 = buffer.data(lh + 521);
    const auto *lh_522 = buffer.data(lh + 522);
    const auto *lh_524 = buffer.data(lh + 524);

    const auto *mf0_200 = buffer.data(mf0 + 200);
    const auto *mf0_201 = buffer.data(mf0 + 201);
    const auto *mf0_205 = buffer.data(mf0 + 205);
    const auto *mf0_206 = buffer.data(mf0 + 206);
    const auto *mf0_208 = buffer.data(mf0 + 208);
    const auto *mf0_209 = buffer.data(mf0 + 209);
    const auto *mf0_210 = buffer.data(mf0 + 210);
    const auto *mf0_212 = buffer.data(mf0 + 212);
    const auto *mf0_213 = buffer.data(mf0 + 213);
    const auto *mf0_216 = buffer.data(mf0 + 216);
    const auto *mf0_217 = buffer.data(mf0 + 217);
    const auto *mf0_219 = buffer.data(mf0 + 219);

    const auto *mf1_200 = buffer.data(mf1 + 200);
    const auto *mf1_201 = buffer.data(mf1 + 201);
    const auto *mf1_205 = buffer.data(mf1 + 205);
    const auto *mf1_206 = buffer.data(mf1 + 206);
    const auto *mf1_208 = buffer.data(mf1 + 208);
    const auto *mf1_209 = buffer.data(mf1 + 209);
    const auto *mf1_210 = buffer.data(mf1 + 210);
    const auto *mf1_212 = buffer.data(mf1 + 212);
    const auto *mf1_213 = buffer.data(mf1 + 213);
    const auto *mf1_216 = buffer.data(mf1 + 216);
    const auto *mf1_217 = buffer.data(mf1 + 217);
    const auto *mf1_219 = buffer.data(mf1 + 219);

    const auto *mg_288 = buffer.data(mg + 288);
    const auto *mg_290 = buffer.data(mg + 290);
    const auto *mg_295 = buffer.data(mg + 295);
    const auto *mg_296 = buffer.data(mg + 296);
    const auto *mg_297 = buffer.data(mg + 297);
    const auto *mg_298 = buffer.data(mg + 298);
    const auto *mg_299 = buffer.data(mg + 299);
    const auto *mg_300 = buffer.data(mg + 300);
    const auto *mg_301 = buffer.data(mg + 301);
    const auto *mg_302 = buffer.data(mg + 302);
    const auto *mg_303 = buffer.data(mg + 303);
    const auto *mg_305 = buffer.data(mg + 305);
    const auto *mg_309 = buffer.data(mg + 309);
    const auto *mg_310 = buffer.data(mg + 310);
    const auto *mg_311 = buffer.data(mg + 311);
    const auto *mg_312 = buffer.data(mg + 312);
    const auto *mg_313 = buffer.data(mg + 313);
    const auto *mg_314 = buffer.data(mg + 314);
    const auto *mg_315 = buffer.data(mg + 315);
    const auto *mg_316 = buffer.data(mg + 316);
    const auto *mg_317 = buffer.data(mg + 317);
    const auto *mg_318 = buffer.data(mg + 318);
    const auto *mg_320 = buffer.data(mg + 320);
    const auto *mg_321 = buffer.data(mg + 321);
    const auto *mg_325 = buffer.data(mg + 325);
    const auto *mg_326 = buffer.data(mg + 326);
    const auto *mg_327 = buffer.data(mg + 327);
    const auto *mg_328 = buffer.data(mg + 328);
    const auto *mg_329 = buffer.data(mg + 329);
    const auto *mg_330 = buffer.data(mg + 330);
    const auto *mg_332 = buffer.data(mg + 332);
    const auto *mg_333 = buffer.data(mg + 333);
    const auto *mg_335 = buffer.data(mg + 335);
    const auto *mg_340 = buffer.data(mg + 340);
    const auto *mg_341 = buffer.data(mg + 341);
    const auto *mg_342 = buffer.data(mg + 342);
    const auto *mg_343 = buffer.data(mg + 343);
    const auto *mg_344 = buffer.data(mg + 344);
    const auto *mg_345 = buffer.data(mg + 345);
    const auto *mg_347 = buffer.data(mg + 347);
    const auto *mg_348 = buffer.data(mg + 348);
    const auto *mg_350 = buffer.data(mg + 350);
    const auto *mg_355 = buffer.data(mg + 355);
    const auto *mg_356 = buffer.data(mg + 356);
    const auto *mg_357 = buffer.data(mg + 357);
    const auto *mg_358 = buffer.data(mg + 358);
    const auto *mg_359 = buffer.data(mg + 359);
    const auto *mg_360 = buffer.data(mg + 360);
    const auto *mg_362 = buffer.data(mg + 362);
    const auto *mg_363 = buffer.data(mg + 363);
    const auto *mg_365 = buffer.data(mg + 365);
    const auto *mg_370 = buffer.data(mg + 370);
    const auto *mg_371 = buffer.data(mg + 371);
    const auto *mg_372 = buffer.data(mg + 372);
    const auto *mg_373 = buffer.data(mg + 373);
    const auto *mg_374 = buffer.data(mg + 374);
    const auto *mg_375 = buffer.data(mg + 375);
    const auto *mg_377 = buffer.data(mg + 377);
    const auto *mg_378 = buffer.data(mg + 378);

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pa_y, pb_x, pb_y, pb_z, lg_198, lg_215, \
                         lg_295, lh_303, mg_288, mg_290, mg_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = f_24 * lg_198[k]
                   + pb_z[k] * mg_288[k];

        t_407[k] = f_7 * lg_215[k]
                   + pb_y[k] * mg_290[k];

        t_408[k] = pa_y[k] * lh_303[k];

        t_409[k] = f_24 * lg_295[k]
                   + pb_x[k] * mg_295[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, pa_y, pb_x, lg_220, lg_296, \
                         lg_297, lg_298, lh_308, lh_309, mg_296, mg_297, \
                         mg_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_24 * lg_296[k]
                   + pb_x[k] * mg_296[k];

        t_411[k] = f_24 * lg_297[k]
                   + pb_x[k] * mg_297[k];

        t_412[k] = f_24 * lg_298[k]
                   + pb_x[k] * mg_298[k];

        t_413[k] = pa_y[k] * lh_308[k];

        t_414[k] = f_11 * lg_220[k]
                   + pa_y[k] * lh_309[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, pa_y, pb_y, pb_z, lg_205, lg_222, lg_223, \
                         lg_224, lh_311, lh_312, mg_295, mg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_24 * lg_205[k]
                   + pb_z[k] * mg_295[k];

        t_416[k] = f_9 * lg_222[k]
                   + pa_y[k] * lh_311[k];

        t_417[k] = f_8 * lg_223[k]
                   + pa_y[k] * lh_312[k];

        t_418[k] = f_7 * lg_224[k]
                   + pb_y[k] * mg_299[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, pa_y, pa_z, pb_y, pb_z, kh0_189, kh1_189, \
                         lg_210, lh_294, lh_314, mg_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = pa_y[k] * lh_314[k];

        t_420[k] = f_25 * kh0_189[k]
                   - f_26 * kh1_189[k]
                   + pa_z[k] * lh_294[k];

        t_421[k] = pb_y[k] * mg_300[k];

        t_422[k] = f_11 * lg_210[k]
                   + pb_z[k] * mg_300[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pb_x, pb_y, lg_305, mf0_200, mf0_205, mf1_200, \
                         mf1_205, mg_301, mg_302, mg_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_3 * mf0_200[k]
                   - f_4 * mf1_200[k]
                   + pb_y[k] * mg_301[k];

        t_424[k] = pb_y[k] * mg_302[k];

        t_425[k] = f_24 * lg_305[k]
                   + f_5 * mf0_205[k]
                   - f_6 * mf1_205[k]
                   + pb_x[k] * mg_305[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, pb_x, pb_y, pb_z, lg_213, lg_309, \
                         mf0_201, mf0_209, mf1_201, mf1_209, mg_303, mg_305, \
                         mg_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_5 * mf0_201[k]
                   - f_6 * mf1_201[k]
                   + pb_y[k] * mg_303[k];

        t_427[k] = f_11 * lg_213[k]
                   + pb_z[k] * mg_303[k];

        t_428[k] = pb_y[k] * mg_305[k];

        t_429[k] = f_24 * lg_309[k]
                   + f_3 * mf0_209[k]
                   - f_4 * mf1_209[k]
                   + pb_x[k] * mg_309[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, pb_x, pb_y, lg_310, lg_311, \
                         lg_312, lg_314, mg_309, mg_310, mg_311, mg_312, \
                         mg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_24 * lg_310[k]
                   + pb_x[k] * mg_310[k];

        t_431[k] = f_24 * lg_311[k]
                   + pb_x[k] * mg_311[k];

        t_432[k] = f_24 * lg_312[k]
                   + pb_x[k] * mg_312[k];

        t_433[k] = pb_y[k] * mg_309[k];

        t_434[k] = f_24 * lg_314[k]
                   + pb_x[k] * mg_314[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, pb_y, pb_z, lg_220, mf0_206, mf0_208, \
                         mf0_209, mf1_206, mf1_208, mf1_209, mg_310, mg_312, \
                         mg_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_1 * mf0_206[k]
                   - f_2 * mf1_206[k]
                   + pb_y[k] * mg_310[k];

        t_436[k] = f_11 * lg_220[k]
                   + pb_z[k] * mg_310[k];

        t_437[k] = f_5 * mf0_208[k]
                   - f_6 * mf1_208[k]
                   + pb_y[k] * mg_312[k];

        t_438[k] = f_3 * mf0_209[k]
                   - f_4 * mf1_209[k]
                   + pb_y[k] * mg_313[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, pa_x, pa_y, pb_y, kh0_210, kh0_440, \
                         kh1_210, kh1_440, lg_225, lh_315, lh_440, mg_314, \
                         mg_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = pb_y[k] * mg_314[k];

        t_440[k] = f_22 * kh0_440[k]
                   - f_23 * kh1_440[k]
                   + pa_x[k] * lh_440[k];

        t_441[k] = f_20 * kh0_210[k]
                   - f_21 * kh1_210[k]
                   + pa_y[k] * lh_315[k];

        t_442[k] = f_19 * lg_225[k]
                   + pb_y[k] * mg_315[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, pb_x, pb_z, lg_318, mf0_210, mf0_213, \
                         mf1_210, mf1_213, mg_315, mg_316, mg_317, \
                         mg_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = pb_z[k] * mg_315[k];

        t_444[k] = f_9 * lg_318[k]
                   + f_5 * mf0_213[k]
                   - f_6 * mf1_213[k]
                   + pb_x[k] * mg_318[k];

        t_445[k] = pb_z[k] * mg_316[k];

        t_446[k] = f_3 * mf0_210[k]
                   - f_4 * mf1_210[k]
                   + pb_z[k] * mg_317[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, pb_x, pb_y, pb_z, lg_230, lg_321, \
                         mf0_212, mf0_216, mf1_212, mf1_216, mg_318, mg_320, \
                         mg_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_9 * lg_321[k]
                   + f_3 * mf0_216[k]
                   - f_4 * mf1_216[k]
                   + pb_x[k] * mg_321[k];

        t_448[k] = pb_z[k] * mg_318[k];

        t_449[k] = f_19 * lg_230[k]
                   + pb_y[k] * mg_320[k];

        t_450[k] = f_5 * mf0_212[k]
                   - f_6 * mf1_212[k]
                   + pb_z[k] * mg_320[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, t_455, pb_x, pb_z, lg_325, lg_327, \
                         lg_328, lg_329, mg_321, mg_325, mg_327, mg_328, \
                         mg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_9 * lg_325[k]
                   + pb_x[k] * mg_325[k];

        t_452[k] = pb_z[k] * mg_321[k];

        t_453[k] = f_9 * lg_327[k]
                   + pb_x[k] * mg_327[k];

        t_454[k] = f_9 * lg_328[k]
                   + pb_x[k] * mg_328[k];

        t_455[k] = f_9 * lg_329[k]
                   + pb_x[k] * mg_329[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pa_x, pb_z, kh0_456, kh1_456, lh_456, \
                         mf0_216, mf0_217, mf1_216, mf1_217, mg_325, mg_326, \
                         mg_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_17 * kh0_456[k]
                   - f_18 * kh1_456[k]
                   + pa_x[k] * lh_456[k];

        t_457[k] = pb_z[k] * mg_325[k];

        t_458[k] = f_3 * mf0_216[k]
                   - f_4 * mf1_216[k]
                   + pb_z[k] * mg_326[k];

        t_459[k] = f_5 * mf0_217[k]
                   - f_6 * mf1_217[k]
                   + pb_z[k] * mg_327[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, pa_z, pb_y, pb_z, lg_225, lg_239, \
                         lh_315, lh_316, mf0_219, mf1_219, mg_329, \
                         mg_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_19 * lg_239[k]
                   + pb_y[k] * mg_329[k];

        t_461[k] = f_1 * mf0_219[k]
                   - f_2 * mf1_219[k]
                   + pb_z[k] * mg_329[k];

        t_462[k] = pa_z[k] * lh_315[k];

        t_463[k] = pa_z[k] * lh_316[k];

        t_464[k] = f_7 * lg_225[k]
                   + pb_z[k] * mg_330[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, pa_z, pb_y, pb_z, lg_227, lg_228, \
                         lg_242, lh_318, lh_320, lh_321, mg_332, \
                         mg_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pa_z[k] * lh_318[k];

        t_466[k] = f_11 * lg_242[k]
                   + pb_y[k] * mg_332[k];

        t_467[k] = f_8 * lg_227[k]
                   + pa_z[k] * lh_320[k];

        t_468[k] = pa_z[k] * lh_321[k];

        t_469[k] = f_7 * lg_228[k]
                   + pb_z[k] * mg_333[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pa_z, pb_x, pb_y, lg_230, lg_245, lg_341, \
                         lh_324, lh_325, mg_335, mg_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_11 * lg_245[k]
                   + pb_y[k] * mg_335[k];

        t_471[k] = f_9 * lg_230[k]
                   + pa_z[k] * lh_324[k];

        t_472[k] = pa_z[k] * lh_325[k];

        t_473[k] = f_9 * lg_341[k]
                   + pb_x[k] * mg_341[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pa_z, pb_x, lg_342, lg_343, lg_344, \
                         lh_330, mg_342, mg_343, mg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_9 * lg_342[k]
                   + pb_x[k] * mg_342[k];

        t_475[k] = f_9 * lg_343[k]
                   + pb_x[k] * mg_343[k];

        t_476[k] = f_9 * lg_344[k]
                   + pb_x[k] * mg_344[k];

        t_477[k] = pa_z[k] * lh_330[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, pa_z, pb_y, pb_z, lg_235, lg_236, lg_237, \
                         lg_254, lh_332, lh_333, mg_340, mg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_7 * lg_235[k]
                   + pb_z[k] * mg_340[k];

        t_479[k] = f_8 * lg_236[k]
                   + pa_z[k] * lh_332[k];

        t_480[k] = f_9 * lg_237[k]
                   + pa_z[k] * lh_333[k];

        t_481[k] = f_11 * lg_254[k]
                   + pb_y[k] * mg_344[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, pa_y, pa_z, pb_y, pb_z, kh0_252, kh1_252, \
                         lg_239, lg_240, lg_255, lh_335, lh_357, \
                         mg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_11 * lg_239[k]
                   + pa_z[k] * lh_335[k];

        t_483[k] = f_22 * kh0_252[k]
                   - f_23 * kh1_252[k]
                   + pa_y[k] * lh_357[k];

        t_484[k] = f_24 * lg_255[k]
                   + pb_y[k] * mg_345[k];

        t_485[k] = f_8 * lg_240[k]
                   + pb_z[k] * mg_345[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, pa_y, pa_z, pb_y, kh0_213, kh0_257, kh1_213, \
                         kh1_257, lg_257, lh_339, lh_362, mg_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_12 * kh0_213[k]
                   - f_13 * kh1_213[k]
                   + pa_z[k] * lh_339[k];

        t_487[k] = f_24 * lg_257[k]
                   + pb_y[k] * mg_347[k];

        t_488[k] = f_22 * kh0_257[k]
                   - f_23 * kh1_257[k]
                   + pa_y[k] * lh_362[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pa_z, pb_y, pb_z, kh0_216, kh1_216, lg_243, \
                         lg_260, lh_342, mg_348, mg_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_12 * kh0_216[k]
                   - f_13 * kh1_216[k]
                   + pa_z[k] * lh_342[k];

        t_490[k] = f_8 * lg_243[k]
                   + pb_z[k] * mg_348[k];

        t_491[k] = f_24 * lg_260[k]
                   + pb_y[k] * mg_350[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pa_y, pb_x, kh0_261, kh1_261, lg_355, \
                         lg_356, lg_357, lh_366, mg_355, mg_356, \
                         mg_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_22 * kh0_261[k]
                   - f_23 * kh1_261[k]
                   + pa_y[k] * lh_366[k];

        t_493[k] = f_9 * lg_355[k]
                   + pb_x[k] * mg_355[k];

        t_494[k] = f_9 * lg_356[k]
                   + pb_x[k] * mg_356[k];

        t_495[k] = f_9 * lg_357[k]
                   + pb_x[k] * mg_357[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pa_x, pb_x, pb_z, kh0_498, kh1_498, \
                         lg_250, lg_358, lg_359, lh_498, mg_355, mg_358, \
                         mg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_9 * lg_358[k]
                   + pb_x[k] * mg_358[k];

        t_497[k] = f_9 * lg_359[k]
                   + pb_x[k] * mg_359[k];

        t_498[k] = f_17 * kh0_498[k]
                   - f_18 * kh1_498[k]
                   + pa_x[k] * lh_498[k];

        t_499[k] = f_8 * lg_250[k]
                   + pb_z[k] * mg_355[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, pa_x, pb_y, kh0_500, kh0_501, kh1_500, kh1_501, \
                         lg_269, lh_500, lh_501, mg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_17 * kh0_500[k]
                   - f_18 * kh1_500[k]
                   + pa_x[k] * lh_500[k];

        t_501[k] = f_17 * kh0_501[k]
                   - f_18 * kh1_501[k]
                   + pa_x[k] * lh_501[k];

        t_502[k] = f_24 * lg_269[k]
                   + pb_y[k] * mg_359[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, pa_x, pa_y, pb_y, kh0_273, kh0_503, kh1_273, \
                         kh1_503, lg_270, lh_378, lh_503, mg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_17 * kh0_503[k]
                   - f_18 * kh1_503[k]
                   + pa_x[k] * lh_503[k];

        t_504[k] = f_17 * kh0_273[k]
                   - f_18 * kh1_273[k]
                   + pa_y[k] * lh_378[k];

        t_505[k] = f_9 * lg_270[k]
                   + pb_y[k] * mg_360[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pa_z, pb_y, pb_z, kh0_234, kh1_234, lg_255, \
                         lg_272, lh_360, mg_360, mg_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_9 * lg_255[k]
                   + pb_z[k] * mg_360[k];

        t_507[k] = f_17 * kh0_234[k]
                   - f_18 * kh1_234[k]
                   + pa_z[k] * lh_360[k];

        t_508[k] = f_9 * lg_272[k]
                   + pb_y[k] * mg_362[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pa_y, pa_z, pb_z, kh0_237, kh0_278, kh1_237, \
                         kh1_278, lg_258, lh_363, lh_383, mg_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_17 * kh0_278[k]
                   - f_18 * kh1_278[k]
                   + pa_y[k] * lh_383[k];

        t_510[k] = f_17 * kh0_237[k]
                   - f_18 * kh1_237[k]
                   + pa_z[k] * lh_363[k];

        t_511[k] = f_9 * lg_258[k]
                   + pb_z[k] * mg_363[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pa_y, pb_x, pb_y, kh0_282, kh1_282, \
                         lg_275, lg_370, lg_371, lh_387, mg_365, mg_370, \
                         mg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_9 * lg_275[k]
                   + pb_y[k] * mg_365[k];

        t_513[k] = f_17 * kh0_282[k]
                   - f_18 * kh1_282[k]
                   + pa_y[k] * lh_387[k];

        t_514[k] = f_9 * lg_370[k]
                   + pb_x[k] * mg_370[k];

        t_515[k] = f_9 * lg_371[k]
                   + pb_x[k] * mg_371[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pa_x, pb_x, kh0_519, kh1_519, lg_372, \
                         lg_373, lg_374, lh_519, mg_372, mg_373, \
                         mg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_9 * lg_372[k]
                   + pb_x[k] * mg_372[k];

        t_517[k] = f_9 * lg_373[k]
                   + pb_x[k] * mg_373[k];

        t_518[k] = f_9 * lg_374[k]
                   + pb_x[k] * mg_374[k];

        t_519[k] = f_17 * kh0_519[k]
                   - f_18 * kh1_519[k]
                   + pa_x[k] * lh_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pa_x, pb_z, kh0_521, kh0_522, kh1_521, kh1_522, \
                         lg_265, lh_521, lh_522, mg_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_9 * lg_265[k]
                   + pb_z[k] * mg_370[k];

        t_521[k] = f_17 * kh0_521[k]
                   - f_18 * kh1_521[k]
                   + pa_x[k] * lh_521[k];

        t_522[k] = f_17 * kh0_522[k]
                   - f_18 * kh1_522[k]
                   + pa_x[k] * lh_522[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pa_x, pa_y, pb_y, kh0_294, kh0_524, kh1_294, \
                         kh1_524, lg_284, lh_399, lh_524, mg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_9 * lg_284[k]
                   + pb_y[k] * mg_374[k];

        t_524[k] = f_17 * kh0_524[k]
                   - f_18 * kh1_524[k]
                   + pa_x[k] * lh_524[k];

        t_525[k] = f_12 * kh0_294[k]
                   - f_13 * kh1_294[k]
                   + pa_y[k] * lh_399[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pa_z, pb_y, pb_z, kh0_255, kh1_255, \
                         lg_270, lg_285, lg_287, lh_381, mg_375, \
                         mg_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_8 * lg_285[k]
                   + pb_y[k] * mg_375[k];

        t_527[k] = f_24 * lg_270[k]
                   + pb_z[k] * mg_375[k];

        t_528[k] = f_22 * kh0_255[k]
                   - f_23 * kh1_255[k]
                   + pa_z[k] * lh_381[k];

        t_529[k] = f_8 * lg_287[k]
                   + pb_y[k] * mg_377[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pa_y, pa_z, pb_z, kh0_258, kh0_299, kh1_258, \
                         kh1_299, lg_273, lh_384, lh_404, mg_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_12 * kh0_299[k]
                   - f_13 * kh1_299[k]
                   + pa_y[k] * lh_404[k];

        t_531[k] = f_22 * kh0_258[k]
                   - f_23 * kh1_258[k]
                   + pa_z[k] * lh_384[k];

        t_532[k] = f_24 * lg_273[k]
                   + pb_z[k] * mg_378[k];
    }
}

static auto
compute_prim_mh_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kh0,
                                            const size_t kh1, const size_t lg, const size_t lh,
                                            const size_t mf0, const size_t mf1, const size_t mg,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_11 = 2.5 / p;
    const auto f_12 = 0.5 / alpha;
    const auto f_13 = 0.5 * beta / (alpha * p);
    const auto f_14 = 3.5 / p;
    const auto f_15 = 3.0 / alpha;
    const auto f_16 = 3.0 * beta / (alpha * p);
    const auto f_17 = 1.0 / alpha;
    const auto f_18 = beta / (alpha * p);
    const auto f_19 = 3.0 / p;
    const auto f_20 = 2.5 / alpha;
    const auto f_21 = 2.5 * beta / (alpha * p);
    const auto f_22 = 1.5 / alpha;
    const auto f_23 = 1.5 * beta / (alpha * p);
    const auto f_24 = 2.0 / p;
    const auto f_25 = 2.0 / alpha;
    const auto f_26 = 2.0 * beta / (alpha * p);

    auto *t_533 = buffer.data(target + 533);
    auto *t_534 = buffer.data(target + 534);
    auto *t_535 = buffer.data(target + 535);
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
    auto *t_659 = buffer.data(target + 659);
    auto *t_660 = buffer.data(target + 660);
    auto *t_661 = buffer.data(target + 661);
    auto *t_662 = buffer.data(target + 662);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kh0_294 = buffer.data(kh0 + 294);
    const auto *kh0_303 = buffer.data(kh0 + 303);
    const auto *kh0_315 = buffer.data(kh0 + 315);
    const auto *kh0_318 = buffer.data(kh0 + 318);
    const auto *kh0_321 = buffer.data(kh0 + 321);
    const auto *kh0_339 = buffer.data(kh0 + 339);
    const auto *kh0_342 = buffer.data(kh0 + 342);
    const auto *kh0_357 = buffer.data(kh0 + 357);
    const auto *kh0_362 = buffer.data(kh0 + 362);
    const auto *kh0_366 = buffer.data(kh0 + 366);
    const auto *kh0_378 = buffer.data(kh0 + 378);
    const auto *kh0_383 = buffer.data(kh0 + 383);
    const auto *kh0_387 = buffer.data(kh0 + 387);
    const auto *kh0_540 = buffer.data(kh0 + 540);
    const auto *kh0_542 = buffer.data(kh0 + 542);
    const auto *kh0_543 = buffer.data(kh0 + 543);
    const auto *kh0_545 = buffer.data(kh0 + 545);
    const auto *kh0_587 = buffer.data(kh0 + 587);
    const auto *kh0_603 = buffer.data(kh0 + 603);
    const auto *kh0_645 = buffer.data(kh0 + 645);
    const auto *kh0_647 = buffer.data(kh0 + 647);
    const auto *kh0_648 = buffer.data(kh0 + 648);
    const auto *kh0_650 = buffer.data(kh0 + 650);

    const auto *kh1_294 = buffer.data(kh1 + 294);
    const auto *kh1_303 = buffer.data(kh1 + 303);
    const auto *kh1_315 = buffer.data(kh1 + 315);
    const auto *kh1_318 = buffer.data(kh1 + 318);
    const auto *kh1_321 = buffer.data(kh1 + 321);
    const auto *kh1_339 = buffer.data(kh1 + 339);
    const auto *kh1_342 = buffer.data(kh1 + 342);
    const auto *kh1_357 = buffer.data(kh1 + 357);
    const auto *kh1_362 = buffer.data(kh1 + 362);
    const auto *kh1_366 = buffer.data(kh1 + 366);
    const auto *kh1_378 = buffer.data(kh1 + 378);
    const auto *kh1_383 = buffer.data(kh1 + 383);
    const auto *kh1_387 = buffer.data(kh1 + 387);
    const auto *kh1_540 = buffer.data(kh1 + 540);
    const auto *kh1_542 = buffer.data(kh1 + 542);
    const auto *kh1_543 = buffer.data(kh1 + 543);
    const auto *kh1_545 = buffer.data(kh1 + 545);
    const auto *kh1_587 = buffer.data(kh1 + 587);
    const auto *kh1_603 = buffer.data(kh1 + 603);
    const auto *kh1_645 = buffer.data(kh1 + 645);
    const auto *kh1_647 = buffer.data(kh1 + 647);
    const auto *kh1_648 = buffer.data(kh1 + 648);
    const auto *kh1_650 = buffer.data(kh1 + 650);

    const auto *lg_280 = buffer.data(lg + 280);
    const auto *lg_288 = buffer.data(lg + 288);
    const auto *lg_290 = buffer.data(lg + 290);
    const auto *lg_295 = buffer.data(lg + 295);
    const auto *lg_299 = buffer.data(lg + 299);
    const auto *lg_300 = buffer.data(lg + 300);
    const auto *lg_301 = buffer.data(lg + 301);
    const auto *lg_302 = buffer.data(lg + 302);
    const auto *lg_303 = buffer.data(lg + 303);
    const auto *lg_305 = buffer.data(lg + 305);
    const auto *lg_310 = buffer.data(lg + 310);
    const auto *lg_312 = buffer.data(lg + 312);
    const auto *lg_313 = buffer.data(lg + 313);
    const auto *lg_314 = buffer.data(lg + 314);
    const auto *lg_315 = buffer.data(lg + 315);
    const auto *lg_317 = buffer.data(lg + 317);
    const auto *lg_318 = buffer.data(lg + 318);
    const auto *lg_320 = buffer.data(lg + 320);
    const auto *lg_325 = buffer.data(lg + 325);
    const auto *lg_326 = buffer.data(lg + 326);
    const auto *lg_327 = buffer.data(lg + 327);
    const auto *lg_329 = buffer.data(lg + 329);
    const auto *lg_330 = buffer.data(lg + 330);
    const auto *lg_332 = buffer.data(lg + 332);
    const auto *lg_333 = buffer.data(lg + 333);
    const auto *lg_335 = buffer.data(lg + 335);
    const auto *lg_340 = buffer.data(lg + 340);
    const auto *lg_344 = buffer.data(lg + 344);
    const auto *lg_345 = buffer.data(lg + 345);
    const auto *lg_347 = buffer.data(lg + 347);
    const auto *lg_348 = buffer.data(lg + 348);
    const auto *lg_350 = buffer.data(lg + 350);
    const auto *lg_359 = buffer.data(lg + 359);
    const auto *lg_360 = buffer.data(lg + 360);
    const auto *lg_362 = buffer.data(lg + 362);
    const auto *lg_365 = buffer.data(lg + 365);
    const auto *lg_385 = buffer.data(lg + 385);
    const auto *lg_386 = buffer.data(lg + 386);
    const auto *lg_387 = buffer.data(lg + 387);
    const auto *lg_388 = buffer.data(lg + 388);
    const auto *lg_389 = buffer.data(lg + 389);
    const auto *lg_400 = buffer.data(lg + 400);
    const auto *lg_401 = buffer.data(lg + 401);
    const auto *lg_402 = buffer.data(lg + 402);
    const auto *lg_403 = buffer.data(lg + 403);
    const auto *lg_410 = buffer.data(lg + 410);
    const auto *lg_414 = buffer.data(lg + 414);
    const auto *lg_415 = buffer.data(lg + 415);
    const auto *lg_416 = buffer.data(lg + 416);
    const auto *lg_417 = buffer.data(lg + 417);
    const auto *lg_419 = buffer.data(lg + 419);
    const auto *lg_423 = buffer.data(lg + 423);
    const auto *lg_426 = buffer.data(lg + 426);
    const auto *lg_430 = buffer.data(lg + 430);
    const auto *lg_432 = buffer.data(lg + 432);
    const auto *lg_433 = buffer.data(lg + 433);
    const auto *lg_434 = buffer.data(lg + 434);
    const auto *lg_446 = buffer.data(lg + 446);
    const auto *lg_447 = buffer.data(lg + 447);
    const auto *lg_448 = buffer.data(lg + 448);
    const auto *lg_449 = buffer.data(lg + 449);
    const auto *lg_460 = buffer.data(lg + 460);
    const auto *lg_461 = buffer.data(lg + 461);
    const auto *lg_462 = buffer.data(lg + 462);
    const auto *lg_463 = buffer.data(lg + 463);
    const auto *lg_464 = buffer.data(lg + 464);
    const auto *lg_475 = buffer.data(lg + 475);
    const auto *lg_476 = buffer.data(lg + 476);

    const auto *lh_408 = buffer.data(lh + 408);
    const auto *lh_420 = buffer.data(lh + 420);
    const auto *lh_422 = buffer.data(lh + 422);
    const auto *lh_423 = buffer.data(lh + 423);
    const auto *lh_425 = buffer.data(lh + 425);
    const auto *lh_426 = buffer.data(lh + 426);
    const auto *lh_429 = buffer.data(lh + 429);
    const auto *lh_434 = buffer.data(lh + 434);
    const auto *lh_435 = buffer.data(lh + 435);
    const auto *lh_437 = buffer.data(lh + 437);
    const auto *lh_438 = buffer.data(lh + 438);
    const auto *lh_440 = buffer.data(lh + 440);
    const auto *lh_441 = buffer.data(lh + 441);
    const auto *lh_442 = buffer.data(lh + 442);
    const auto *lh_444 = buffer.data(lh + 444);
    const auto *lh_446 = buffer.data(lh + 446);
    const auto *lh_447 = buffer.data(lh + 447);
    const auto *lh_450 = buffer.data(lh + 450);
    const auto *lh_451 = buffer.data(lh + 451);
    const auto *lh_456 = buffer.data(lh + 456);
    const auto *lh_458 = buffer.data(lh + 458);
    const auto *lh_459 = buffer.data(lh + 459);
    const auto *lh_461 = buffer.data(lh + 461);
    const auto *lh_465 = buffer.data(lh + 465);
    const auto *lh_468 = buffer.data(lh + 468);
    const auto *lh_483 = buffer.data(lh + 483);
    const auto *lh_486 = buffer.data(lh + 486);
    const auto *lh_488 = buffer.data(lh + 488);
    const auto *lh_489 = buffer.data(lh + 489);
    const auto *lh_492 = buffer.data(lh + 492);
    const auto *lh_504 = buffer.data(lh + 504);
    const auto *lh_509 = buffer.data(lh + 509);
    const auto *lh_513 = buffer.data(lh + 513);
    const auto *lh_540 = buffer.data(lh + 540);
    const auto *lh_542 = buffer.data(lh + 542);
    const auto *lh_543 = buffer.data(lh + 543);
    const auto *lh_545 = buffer.data(lh + 545);
    const auto *lh_587 = buffer.data(lh + 587);
    const auto *lh_603 = buffer.data(lh + 603);
    const auto *lh_645 = buffer.data(lh + 645);
    const auto *lh_647 = buffer.data(lh + 647);
    const auto *lh_648 = buffer.data(lh + 648);
    const auto *lh_650 = buffer.data(lh + 650);

    const auto *mf0_270 = buffer.data(mf0 + 270);
    const auto *mf0_271 = buffer.data(mf0 + 271);
    const auto *mf0_275 = buffer.data(mf0 + 275);
    const auto *mf0_276 = buffer.data(mf0 + 276);
    const auto *mf0_278 = buffer.data(mf0 + 278);
    const auto *mf0_279 = buffer.data(mf0 + 279);
    const auto *mf0_280 = buffer.data(mf0 + 280);
    const auto *mf0_282 = buffer.data(mf0 + 282);
    const auto *mf0_283 = buffer.data(mf0 + 283);
    const auto *mf0_286 = buffer.data(mf0 + 286);
    const auto *mf0_287 = buffer.data(mf0 + 287);
    const auto *mf0_289 = buffer.data(mf0 + 289);

    const auto *mf1_270 = buffer.data(mf1 + 270);
    const auto *mf1_271 = buffer.data(mf1 + 271);
    const auto *mf1_275 = buffer.data(mf1 + 275);
    const auto *mf1_276 = buffer.data(mf1 + 276);
    const auto *mf1_278 = buffer.data(mf1 + 278);
    const auto *mf1_279 = buffer.data(mf1 + 279);
    const auto *mf1_280 = buffer.data(mf1 + 280);
    const auto *mf1_282 = buffer.data(mf1 + 282);
    const auto *mf1_283 = buffer.data(mf1 + 283);
    const auto *mf1_286 = buffer.data(mf1 + 286);
    const auto *mf1_287 = buffer.data(mf1 + 287);
    const auto *mf1_289 = buffer.data(mf1 + 289);

    const auto *mg_380 = buffer.data(mg + 380);
    const auto *mg_385 = buffer.data(mg + 385);
    const auto *mg_386 = buffer.data(mg + 386);
    const auto *mg_387 = buffer.data(mg + 387);
    const auto *mg_388 = buffer.data(mg + 388);
    const auto *mg_389 = buffer.data(mg + 389);
    const auto *mg_390 = buffer.data(mg + 390);
    const auto *mg_392 = buffer.data(mg + 392);
    const auto *mg_393 = buffer.data(mg + 393);
    const auto *mg_395 = buffer.data(mg + 395);
    const auto *mg_400 = buffer.data(mg + 400);
    const auto *mg_401 = buffer.data(mg + 401);
    const auto *mg_402 = buffer.data(mg + 402);
    const auto *mg_403 = buffer.data(mg + 403);
    const auto *mg_404 = buffer.data(mg + 404);
    const auto *mg_405 = buffer.data(mg + 405);
    const auto *mg_406 = buffer.data(mg + 406);
    const auto *mg_407 = buffer.data(mg + 407);
    const auto *mg_408 = buffer.data(mg + 408);
    const auto *mg_410 = buffer.data(mg + 410);
    const auto *mg_414 = buffer.data(mg + 414);
    const auto *mg_415 = buffer.data(mg + 415);
    const auto *mg_416 = buffer.data(mg + 416);
    const auto *mg_417 = buffer.data(mg + 417);
    const auto *mg_418 = buffer.data(mg + 418);
    const auto *mg_419 = buffer.data(mg + 419);
    const auto *mg_420 = buffer.data(mg + 420);
    const auto *mg_421 = buffer.data(mg + 421);
    const auto *mg_422 = buffer.data(mg + 422);
    const auto *mg_423 = buffer.data(mg + 423);
    const auto *mg_425 = buffer.data(mg + 425);
    const auto *mg_426 = buffer.data(mg + 426);
    const auto *mg_430 = buffer.data(mg + 430);
    const auto *mg_431 = buffer.data(mg + 431);
    const auto *mg_432 = buffer.data(mg + 432);
    const auto *mg_433 = buffer.data(mg + 433);
    const auto *mg_434 = buffer.data(mg + 434);
    const auto *mg_435 = buffer.data(mg + 435);
    const auto *mg_437 = buffer.data(mg + 437);
    const auto *mg_438 = buffer.data(mg + 438);
    const auto *mg_440 = buffer.data(mg + 440);
    const auto *mg_445 = buffer.data(mg + 445);
    const auto *mg_446 = buffer.data(mg + 446);
    const auto *mg_447 = buffer.data(mg + 447);
    const auto *mg_448 = buffer.data(mg + 448);
    const auto *mg_449 = buffer.data(mg + 449);
    const auto *mg_450 = buffer.data(mg + 450);
    const auto *mg_452 = buffer.data(mg + 452);
    const auto *mg_453 = buffer.data(mg + 453);
    const auto *mg_455 = buffer.data(mg + 455);
    const auto *mg_460 = buffer.data(mg + 460);
    const auto *mg_461 = buffer.data(mg + 461);
    const auto *mg_462 = buffer.data(mg + 462);
    const auto *mg_463 = buffer.data(mg + 463);
    const auto *mg_464 = buffer.data(mg + 464);
    const auto *mg_465 = buffer.data(mg + 465);
    const auto *mg_467 = buffer.data(mg + 467);
    const auto *mg_468 = buffer.data(mg + 468);
    const auto *mg_470 = buffer.data(mg + 470);
    const auto *mg_475 = buffer.data(mg + 475);
    const auto *mg_476 = buffer.data(mg + 476);

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pa_y, pb_x, pb_y, kh0_303, kh1_303, \
                         lg_290, lg_385, lg_386, lh_408, mg_380, mg_385, \
                         mg_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_8 * lg_290[k]
                   + pb_y[k] * mg_380[k];

        t_534[k] = f_12 * kh0_303[k]
                   - f_13 * kh1_303[k]
                   + pa_y[k] * lh_408[k];

        t_535[k] = f_9 * lg_385[k]
                   + pb_x[k] * mg_385[k];

        t_536[k] = f_9 * lg_386[k]
                   + pb_x[k] * mg_386[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pa_x, pb_x, kh0_540, kh1_540, lg_387, \
                         lg_388, lg_389, lh_540, mg_387, mg_388, \
                         mg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_9 * lg_387[k]
                   + pb_x[k] * mg_387[k];

        t_538[k] = f_9 * lg_388[k]
                   + pb_x[k] * mg_388[k];

        t_539[k] = f_9 * lg_389[k]
                   + pb_x[k] * mg_389[k];

        t_540[k] = f_17 * kh0_540[k]
                   - f_18 * kh1_540[k]
                   + pa_x[k] * lh_540[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pa_x, pb_z, kh0_542, kh0_543, kh1_542, kh1_543, \
                         lg_280, lh_542, lh_543, mg_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_24 * lg_280[k]
                   + pb_z[k] * mg_385[k];

        t_542[k] = f_17 * kh0_542[k]
                   - f_18 * kh1_542[k]
                   + pa_x[k] * lh_542[k];

        t_543[k] = f_17 * kh0_543[k]
                   - f_18 * kh1_543[k]
                   + pa_x[k] * lh_543[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pa_x, pa_y, pb_y, kh0_545, kh1_545, \
                         lg_299, lg_300, lh_420, lh_545, mg_389, \
                         mg_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_8 * lg_299[k]
                   + pb_y[k] * mg_389[k];

        t_545[k] = f_17 * kh0_545[k]
                   - f_18 * kh1_545[k]
                   + pa_x[k] * lh_545[k];

        t_546[k] = pa_y[k] * lh_420[k];

        t_547[k] = f_7 * lg_300[k]
                   + pb_y[k] * mg_390[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pa_y, pb_y, lg_301, lg_302, \
                         lg_303, lh_422, lh_423, lh_425, lh_426, \
                         mg_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = pa_y[k] * lh_422[k];

        t_549[k] = f_8 * lg_301[k]
                   + pa_y[k] * lh_423[k];

        t_550[k] = f_7 * lg_302[k]
                   + pb_y[k] * mg_392[k];

        t_551[k] = pa_y[k] * lh_425[k];

        t_552[k] = f_9 * lg_303[k]
                   + pa_y[k] * lh_426[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, pa_y, pb_x, pb_y, pb_z, lg_288, lg_305, \
                         lg_400, lh_429, mg_393, mg_395, mg_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_11 * lg_288[k]
                   + pb_z[k] * mg_393[k];

        t_554[k] = f_7 * lg_305[k]
                   + pb_y[k] * mg_395[k];

        t_555[k] = pa_y[k] * lh_429[k];

        t_556[k] = f_9 * lg_400[k]
                   + pb_x[k] * mg_400[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, pa_y, pb_x, lg_310, lg_401, \
                         lg_402, lg_403, lh_434, lh_435, mg_401, mg_402, \
                         mg_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_9 * lg_401[k]
                   + pb_x[k] * mg_401[k];

        t_558[k] = f_9 * lg_402[k]
                   + pb_x[k] * mg_402[k];

        t_559[k] = f_9 * lg_403[k]
                   + pb_x[k] * mg_403[k];

        t_560[k] = pa_y[k] * lh_434[k];

        t_561[k] = f_11 * lg_310[k]
                   + pa_y[k] * lh_435[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, pa_y, pb_y, pb_z, lg_295, lg_312, lg_313, \
                         lg_314, lh_437, lh_438, mg_400, mg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_11 * lg_295[k]
                   + pb_z[k] * mg_400[k];

        t_563[k] = f_9 * lg_312[k]
                   + pa_y[k] * lh_437[k];

        t_564[k] = f_8 * lg_313[k]
                   + pa_y[k] * lh_438[k];

        t_565[k] = f_7 * lg_314[k]
                   + pb_y[k] * mg_404[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pa_y, pa_z, pb_y, pb_z, kh0_294, kh1_294, \
                         lg_300, lh_420, lh_440, mg_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = pa_y[k] * lh_440[k];

        t_567[k] = f_20 * kh0_294[k]
                   - f_21 * kh1_294[k]
                   + pa_z[k] * lh_420[k];

        t_568[k] = pb_y[k] * mg_405[k];

        t_569[k] = f_19 * lg_300[k]
                   + pb_z[k] * mg_405[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, pb_x, pb_y, lg_410, mf0_270, mf0_275, mf1_270, \
                         mf1_275, mg_406, mg_407, mg_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_3 * mf0_270[k]
                   - f_4 * mf1_270[k]
                   + pb_y[k] * mg_406[k];

        t_571[k] = pb_y[k] * mg_407[k];

        t_572[k] = f_9 * lg_410[k]
                   + f_5 * mf0_275[k]
                   - f_6 * mf1_275[k]
                   + pb_x[k] * mg_410[k];
    }

#pragma omp simd aligned(t_573, t_574, t_575, t_576, pb_x, pb_y, pb_z, lg_303, lg_414, \
                         mf0_271, mf0_279, mf1_271, mf1_279, mg_408, mg_410, \
                         mg_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_573[k] = f_5 * mf0_271[k]
                   - f_6 * mf1_271[k]
                   + pb_y[k] * mg_408[k];

        t_574[k] = f_19 * lg_303[k]
                   + pb_z[k] * mg_408[k];

        t_575[k] = pb_y[k] * mg_410[k];

        t_576[k] = f_9 * lg_414[k]
                   + f_3 * mf0_279[k]
                   - f_4 * mf1_279[k]
                   + pb_x[k] * mg_414[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, pb_x, pb_y, lg_415, lg_416, \
                         lg_417, lg_419, mg_414, mg_415, mg_416, mg_417, \
                         mg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = f_9 * lg_415[k]
                   + pb_x[k] * mg_415[k];

        t_578[k] = f_9 * lg_416[k]
                   + pb_x[k] * mg_416[k];

        t_579[k] = f_9 * lg_417[k]
                   + pb_x[k] * mg_417[k];

        t_580[k] = pb_y[k] * mg_414[k];

        t_581[k] = f_9 * lg_419[k]
                   + pb_x[k] * mg_419[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pb_y, pb_z, lg_310, mf0_276, mf0_278, \
                         mf0_279, mf1_276, mf1_278, mf1_279, mg_415, mg_417, \
                         mg_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_1 * mf0_276[k]
                   - f_2 * mf1_276[k]
                   + pb_y[k] * mg_415[k];

        t_583[k] = f_19 * lg_310[k]
                   + pb_z[k] * mg_415[k];

        t_584[k] = f_5 * mf0_278[k]
                   - f_6 * mf1_278[k]
                   + pb_y[k] * mg_417[k];

        t_585[k] = f_3 * mf0_279[k]
                   - f_4 * mf1_279[k]
                   + pb_y[k] * mg_418[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, pa_x, pa_y, pb_y, kh0_315, kh0_587, \
                         kh1_315, kh1_587, lg_315, lh_441, lh_587, mg_419, \
                         mg_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = pb_y[k] * mg_419[k];

        t_587[k] = f_17 * kh0_587[k]
                   - f_18 * kh1_587[k]
                   + pa_x[k] * lh_587[k];

        t_588[k] = f_15 * kh0_315[k]
                   - f_16 * kh1_315[k]
                   + pa_y[k] * lh_441[k];

        t_589[k] = f_14 * lg_315[k]
                   + pb_y[k] * mg_420[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, pb_x, pb_z, lg_423, mf0_280, mf0_283, \
                         mf1_280, mf1_283, mg_420, mg_421, mg_422, \
                         mg_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = pb_z[k] * mg_420[k];

        t_591[k] = f_8 * lg_423[k]
                   + f_5 * mf0_283[k]
                   - f_6 * mf1_283[k]
                   + pb_x[k] * mg_423[k];

        t_592[k] = pb_z[k] * mg_421[k];

        t_593[k] = f_3 * mf0_280[k]
                   - f_4 * mf1_280[k]
                   + pb_z[k] * mg_422[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, pb_x, pb_y, pb_z, lg_320, lg_426, \
                         mf0_282, mf0_286, mf1_282, mf1_286, mg_423, mg_425, \
                         mg_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_8 * lg_426[k]
                   + f_3 * mf0_286[k]
                   - f_4 * mf1_286[k]
                   + pb_x[k] * mg_426[k];

        t_595[k] = pb_z[k] * mg_423[k];

        t_596[k] = f_14 * lg_320[k]
                   + pb_y[k] * mg_425[k];

        t_597[k] = f_5 * mf0_282[k]
                   - f_6 * mf1_282[k]
                   + pb_z[k] * mg_425[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, t_601, t_602, pb_x, pb_z, lg_430, lg_432, \
                         lg_433, lg_434, mg_426, mg_430, mg_432, mg_433, \
                         mg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_8 * lg_430[k]
                   + pb_x[k] * mg_430[k];

        t_599[k] = pb_z[k] * mg_426[k];

        t_600[k] = f_8 * lg_432[k]
                   + pb_x[k] * mg_432[k];

        t_601[k] = f_8 * lg_433[k]
                   + pb_x[k] * mg_433[k];

        t_602[k] = f_8 * lg_434[k]
                   + pb_x[k] * mg_434[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, t_606, pa_x, pb_z, kh0_603, kh1_603, lh_603, \
                         mf0_286, mf0_287, mf1_286, mf1_287, mg_430, mg_431, \
                         mg_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = f_12 * kh0_603[k]
                   - f_13 * kh1_603[k]
                   + pa_x[k] * lh_603[k];

        t_604[k] = pb_z[k] * mg_430[k];

        t_605[k] = f_3 * mf0_286[k]
                   - f_4 * mf1_286[k]
                   + pb_z[k] * mg_431[k];

        t_606[k] = f_5 * mf0_287[k]
                   - f_6 * mf1_287[k]
                   + pb_z[k] * mg_432[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, t_611, pa_z, pb_y, pb_z, lg_315, lg_329, \
                         lh_441, lh_442, mf0_289, mf1_289, mg_434, \
                         mg_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_14 * lg_329[k]
                   + pb_y[k] * mg_434[k];

        t_608[k] = f_1 * mf0_289[k]
                   - f_2 * mf1_289[k]
                   + pb_z[k] * mg_434[k];

        t_609[k] = pa_z[k] * lh_441[k];

        t_610[k] = pa_z[k] * lh_442[k];

        t_611[k] = f_7 * lg_315[k]
                   + pb_z[k] * mg_435[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, t_616, pa_z, pb_y, pb_z, lg_317, lg_318, \
                         lg_332, lh_444, lh_446, lh_447, mg_437, \
                         mg_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = pa_z[k] * lh_444[k];

        t_613[k] = f_19 * lg_332[k]
                   + pb_y[k] * mg_437[k];

        t_614[k] = f_8 * lg_317[k]
                   + pa_z[k] * lh_446[k];

        t_615[k] = pa_z[k] * lh_447[k];

        t_616[k] = f_7 * lg_318[k]
                   + pb_z[k] * mg_438[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, pa_z, pb_x, pb_y, lg_320, lg_335, lg_446, \
                         lh_450, lh_451, mg_440, mg_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_19 * lg_335[k]
                   + pb_y[k] * mg_440[k];

        t_618[k] = f_9 * lg_320[k]
                   + pa_z[k] * lh_450[k];

        t_619[k] = pa_z[k] * lh_451[k];

        t_620[k] = f_8 * lg_446[k]
                   + pb_x[k] * mg_446[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, t_624, pa_z, pb_x, lg_447, lg_448, lg_449, \
                         lh_456, mg_447, mg_448, mg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_8 * lg_447[k]
                   + pb_x[k] * mg_447[k];

        t_622[k] = f_8 * lg_448[k]
                   + pb_x[k] * mg_448[k];

        t_623[k] = f_8 * lg_449[k]
                   + pb_x[k] * mg_449[k];

        t_624[k] = pa_z[k] * lh_456[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, pa_z, pb_y, pb_z, lg_325, lg_326, lg_327, \
                         lg_344, lh_458, lh_459, mg_445, mg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = f_7 * lg_325[k]
                   + pb_z[k] * mg_445[k];

        t_626[k] = f_8 * lg_326[k]
                   + pa_z[k] * lh_458[k];

        t_627[k] = f_9 * lg_327[k]
                   + pa_z[k] * lh_459[k];

        t_628[k] = f_19 * lg_344[k]
                   + pb_y[k] * mg_449[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, pa_y, pa_z, pb_y, pb_z, kh0_357, kh1_357, \
                         lg_329, lg_330, lg_345, lh_461, lh_483, \
                         mg_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_11 * lg_329[k]
                   + pa_z[k] * lh_461[k];

        t_630[k] = f_25 * kh0_357[k]
                   - f_26 * kh1_357[k]
                   + pa_y[k] * lh_483[k];

        t_631[k] = f_11 * lg_345[k]
                   + pb_y[k] * mg_450[k];

        t_632[k] = f_8 * lg_330[k]
                   + pb_z[k] * mg_450[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pa_y, pa_z, pb_y, kh0_318, kh0_362, kh1_318, \
                         kh1_362, lg_347, lh_465, lh_488, mg_452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_12 * kh0_318[k]
                   - f_13 * kh1_318[k]
                   + pa_z[k] * lh_465[k];

        t_634[k] = f_11 * lg_347[k]
                   + pb_y[k] * mg_452[k];

        t_635[k] = f_25 * kh0_362[k]
                   - f_26 * kh1_362[k]
                   + pa_y[k] * lh_488[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pa_z, pb_y, pb_z, kh0_321, kh1_321, lg_333, \
                         lg_350, lh_468, mg_453, mg_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_12 * kh0_321[k]
                   - f_13 * kh1_321[k]
                   + pa_z[k] * lh_468[k];

        t_637[k] = f_8 * lg_333[k]
                   + pb_z[k] * mg_453[k];

        t_638[k] = f_11 * lg_350[k]
                   + pb_y[k] * mg_455[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, pa_y, pb_x, kh0_366, kh1_366, lg_460, \
                         lg_461, lg_462, lh_492, mg_460, mg_461, \
                         mg_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_25 * kh0_366[k]
                   - f_26 * kh1_366[k]
                   + pa_y[k] * lh_492[k];

        t_640[k] = f_8 * lg_460[k]
                   + pb_x[k] * mg_460[k];

        t_641[k] = f_8 * lg_461[k]
                   + pb_x[k] * mg_461[k];

        t_642[k] = f_8 * lg_462[k]
                   + pb_x[k] * mg_462[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, t_646, pa_x, pb_x, pb_z, kh0_645, kh1_645, \
                         lg_340, lg_463, lg_464, lh_645, mg_460, mg_463, \
                         mg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_8 * lg_463[k]
                   + pb_x[k] * mg_463[k];

        t_644[k] = f_8 * lg_464[k]
                   + pb_x[k] * mg_464[k];

        t_645[k] = f_12 * kh0_645[k]
                   - f_13 * kh1_645[k]
                   + pa_x[k] * lh_645[k];

        t_646[k] = f_8 * lg_340[k]
                   + pb_z[k] * mg_460[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pa_x, pb_y, kh0_647, kh0_648, kh1_647, kh1_648, \
                         lg_359, lh_647, lh_648, mg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_12 * kh0_647[k]
                   - f_13 * kh1_647[k]
                   + pa_x[k] * lh_647[k];

        t_648[k] = f_12 * kh0_648[k]
                   - f_13 * kh1_648[k]
                   + pa_x[k] * lh_648[k];

        t_649[k] = f_11 * lg_359[k]
                   + pb_y[k] * mg_464[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pa_x, pa_y, pb_y, kh0_378, kh0_650, kh1_378, \
                         kh1_650, lg_360, lh_504, lh_650, mg_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_12 * kh0_650[k]
                   - f_13 * kh1_650[k]
                   + pa_x[k] * lh_650[k];

        t_651[k] = f_22 * kh0_378[k]
                   - f_23 * kh1_378[k]
                   + pa_y[k] * lh_504[k];

        t_652[k] = f_24 * lg_360[k]
                   + pb_y[k] * mg_465[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pa_z, pb_y, pb_z, kh0_339, kh1_339, lg_345, \
                         lg_362, lh_486, mg_465, mg_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_9 * lg_345[k]
                   + pb_z[k] * mg_465[k];

        t_654[k] = f_17 * kh0_339[k]
                   - f_18 * kh1_339[k]
                   + pa_z[k] * lh_486[k];

        t_655[k] = f_24 * lg_362[k]
                   + pb_y[k] * mg_467[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pa_y, pa_z, pb_z, kh0_342, kh0_383, kh1_342, \
                         kh1_383, lg_348, lh_489, lh_509, mg_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_22 * kh0_383[k]
                   - f_23 * kh1_383[k]
                   + pa_y[k] * lh_509[k];

        t_657[k] = f_17 * kh0_342[k]
                   - f_18 * kh1_342[k]
                   + pa_z[k] * lh_489[k];

        t_658[k] = f_9 * lg_348[k]
                   + pb_z[k] * mg_468[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, pa_y, pb_x, pb_y, kh0_387, kh1_387, \
                         lg_365, lg_475, lg_476, lh_513, mg_470, mg_475, \
                         mg_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_24 * lg_365[k]
                   + pb_y[k] * mg_470[k];

        t_660[k] = f_22 * kh0_387[k]
                   - f_23 * kh1_387[k]
                   + pa_y[k] * lh_513[k];

        t_661[k] = f_8 * lg_475[k]
                   + pb_x[k] * mg_475[k];

        t_662[k] = f_8 * lg_476[k]
                   + pb_x[k] * mg_476[k];
    }
}

static auto
compute_prim_mh_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kh0,
                                            const size_t kh1, const size_t lg, const size_t lh,
                                            const size_t mf0, const size_t mf1, const size_t mg,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 4.0 / p;
    const auto f_11 = 2.5 / p;
    const auto f_12 = 0.5 / alpha;
    const auto f_13 = 0.5 * beta / (alpha * p);
    const auto f_14 = 3.5 / p;
    const auto f_15 = 3.0 / alpha;
    const auto f_16 = 3.0 * beta / (alpha * p);
    const auto f_17 = 1.0 / alpha;
    const auto f_18 = beta / (alpha * p);
    const auto f_19 = 3.0 / p;
    const auto f_22 = 1.5 / alpha;
    const auto f_23 = 1.5 * beta / (alpha * p);
    const auto f_24 = 2.0 / p;
    const auto f_25 = 2.0 / alpha;
    const auto f_26 = 2.0 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kh0_360 = buffer.data(kh0 + 360);
    const auto *kh0_363 = buffer.data(kh0 + 363);
    const auto *kh0_381 = buffer.data(kh0 + 381);
    const auto *kh0_384 = buffer.data(kh0 + 384);
    const auto *kh0_399 = buffer.data(kh0 + 399);
    const auto *kh0_404 = buffer.data(kh0 + 404);
    const auto *kh0_408 = buffer.data(kh0 + 408);
    const auto *kh0_420 = buffer.data(kh0 + 420);
    const auto *kh0_425 = buffer.data(kh0 + 425);
    const auto *kh0_429 = buffer.data(kh0 + 429);
    const auto *kh0_666 = buffer.data(kh0 + 666);
    const auto *kh0_668 = buffer.data(kh0 + 668);
    const auto *kh0_669 = buffer.data(kh0 + 669);
    const auto *kh0_671 = buffer.data(kh0 + 671);
    const auto *kh0_687 = buffer.data(kh0 + 687);
    const auto *kh0_689 = buffer.data(kh0 + 689);
    const auto *kh0_690 = buffer.data(kh0 + 690);
    const auto *kh0_692 = buffer.data(kh0 + 692);
    const auto *kh0_708 = buffer.data(kh0 + 708);
    const auto *kh0_710 = buffer.data(kh0 + 710);
    const auto *kh0_711 = buffer.data(kh0 + 711);
    const auto *kh0_713 = buffer.data(kh0 + 713);
    const auto *kh0_755 = buffer.data(kh0 + 755);

    const auto *kh1_360 = buffer.data(kh1 + 360);
    const auto *kh1_363 = buffer.data(kh1 + 363);
    const auto *kh1_381 = buffer.data(kh1 + 381);
    const auto *kh1_384 = buffer.data(kh1 + 384);
    const auto *kh1_399 = buffer.data(kh1 + 399);
    const auto *kh1_404 = buffer.data(kh1 + 404);
    const auto *kh1_408 = buffer.data(kh1 + 408);
    const auto *kh1_420 = buffer.data(kh1 + 420);
    const auto *kh1_425 = buffer.data(kh1 + 425);
    const auto *kh1_429 = buffer.data(kh1 + 429);
    const auto *kh1_666 = buffer.data(kh1 + 666);
    const auto *kh1_668 = buffer.data(kh1 + 668);
    const auto *kh1_669 = buffer.data(kh1 + 669);
    const auto *kh1_671 = buffer.data(kh1 + 671);
    const auto *kh1_687 = buffer.data(kh1 + 687);
    const auto *kh1_689 = buffer.data(kh1 + 689);
    const auto *kh1_690 = buffer.data(kh1 + 690);
    const auto *kh1_692 = buffer.data(kh1 + 692);
    const auto *kh1_708 = buffer.data(kh1 + 708);
    const auto *kh1_710 = buffer.data(kh1 + 710);
    const auto *kh1_711 = buffer.data(kh1 + 711);
    const auto *kh1_713 = buffer.data(kh1 + 713);
    const auto *kh1_755 = buffer.data(kh1 + 755);

    const auto *lg_355 = buffer.data(lg + 355);
    const auto *lg_360 = buffer.data(lg + 360);
    const auto *lg_363 = buffer.data(lg + 363);
    const auto *lg_370 = buffer.data(lg + 370);
    const auto *lg_374 = buffer.data(lg + 374);
    const auto *lg_375 = buffer.data(lg + 375);
    const auto *lg_377 = buffer.data(lg + 377);
    const auto *lg_378 = buffer.data(lg + 378);
    const auto *lg_380 = buffer.data(lg + 380);
    const auto *lg_385 = buffer.data(lg + 385);
    const auto *lg_389 = buffer.data(lg + 389);
    const auto *lg_390 = buffer.data(lg + 390);
    const auto *lg_392 = buffer.data(lg + 392);
    const auto *lg_393 = buffer.data(lg + 393);
    const auto *lg_395 = buffer.data(lg + 395);
    const auto *lg_400 = buffer.data(lg + 400);
    const auto *lg_404 = buffer.data(lg + 404);
    const auto *lg_405 = buffer.data(lg + 405);
    const auto *lg_406 = buffer.data(lg + 406);
    const auto *lg_407 = buffer.data(lg + 407);
    const auto *lg_408 = buffer.data(lg + 408);
    const auto *lg_410 = buffer.data(lg + 410);
    const auto *lg_415 = buffer.data(lg + 415);
    const auto *lg_417 = buffer.data(lg + 417);
    const auto *lg_418 = buffer.data(lg + 418);
    const auto *lg_419 = buffer.data(lg + 419);
    const auto *lg_420 = buffer.data(lg + 420);
    const auto *lg_423 = buffer.data(lg + 423);
    const auto *lg_425 = buffer.data(lg + 425);
    const auto *lg_437 = buffer.data(lg + 437);
    const auto *lg_440 = buffer.data(lg + 440);
    const auto *lg_477 = buffer.data(lg + 477);
    const auto *lg_478 = buffer.data(lg + 478);
    const auto *lg_479 = buffer.data(lg + 479);
    const auto *lg_490 = buffer.data(lg + 490);
    const auto *lg_491 = buffer.data(lg + 491);
    const auto *lg_492 = buffer.data(lg + 492);
    const auto *lg_493 = buffer.data(lg + 493);
    const auto *lg_494 = buffer.data(lg + 494);
    const auto *lg_505 = buffer.data(lg + 505);
    const auto *lg_506 = buffer.data(lg + 506);
    const auto *lg_507 = buffer.data(lg + 507);
    const auto *lg_508 = buffer.data(lg + 508);
    const auto *lg_509 = buffer.data(lg + 509);
    const auto *lg_520 = buffer.data(lg + 520);
    const auto *lg_521 = buffer.data(lg + 521);
    const auto *lg_522 = buffer.data(lg + 522);
    const auto *lg_523 = buffer.data(lg + 523);
    const auto *lg_530 = buffer.data(lg + 530);
    const auto *lg_534 = buffer.data(lg + 534);
    const auto *lg_535 = buffer.data(lg + 535);
    const auto *lg_536 = buffer.data(lg + 536);
    const auto *lg_537 = buffer.data(lg + 537);
    const auto *lg_539 = buffer.data(lg + 539);
    const auto *lg_540 = buffer.data(lg + 540);
    const auto *lg_543 = buffer.data(lg + 543);
    const auto *lg_545 = buffer.data(lg + 545);
    const auto *lg_546 = buffer.data(lg + 546);
    const auto *lg_549 = buffer.data(lg + 549);
    const auto *lg_550 = buffer.data(lg + 550);
    const auto *lg_552 = buffer.data(lg + 552);
    const auto *lg_553 = buffer.data(lg + 553);
    const auto *lg_554 = buffer.data(lg + 554);
    const auto *lg_560 = buffer.data(lg + 560);
    const auto *lg_564 = buffer.data(lg + 564);
    const auto *lg_566 = buffer.data(lg + 566);
    const auto *lg_567 = buffer.data(lg + 567);
    const auto *lg_568 = buffer.data(lg + 568);
    const auto *lg_569 = buffer.data(lg + 569);
    const auto *lg_570 = buffer.data(lg + 570);

    const auto *lh_507 = buffer.data(lh + 507);
    const auto *lh_510 = buffer.data(lh + 510);
    const auto *lh_525 = buffer.data(lh + 525);
    const auto *lh_528 = buffer.data(lh + 528);
    const auto *lh_530 = buffer.data(lh + 530);
    const auto *lh_531 = buffer.data(lh + 531);
    const auto *lh_534 = buffer.data(lh + 534);
    const auto *lh_546 = buffer.data(lh + 546);
    const auto *lh_551 = buffer.data(lh + 551);
    const auto *lh_555 = buffer.data(lh + 555);
    const auto *lh_567 = buffer.data(lh + 567);
    const auto *lh_569 = buffer.data(lh + 569);
    const auto *lh_570 = buffer.data(lh + 570);
    const auto *lh_572 = buffer.data(lh + 572);
    const auto *lh_573 = buffer.data(lh + 573);
    const auto *lh_576 = buffer.data(lh + 576);
    const auto *lh_581 = buffer.data(lh + 581);
    const auto *lh_582 = buffer.data(lh + 582);
    const auto *lh_584 = buffer.data(lh + 584);
    const auto *lh_585 = buffer.data(lh + 585);
    const auto *lh_587 = buffer.data(lh + 587);
    const auto *lh_588 = buffer.data(lh + 588);
    const auto *lh_589 = buffer.data(lh + 589);
    const auto *lh_591 = buffer.data(lh + 591);
    const auto *lh_594 = buffer.data(lh + 594);
    const auto *lh_598 = buffer.data(lh + 598);
    const auto *lh_666 = buffer.data(lh + 666);
    const auto *lh_668 = buffer.data(lh + 668);
    const auto *lh_669 = buffer.data(lh + 669);
    const auto *lh_671 = buffer.data(lh + 671);
    const auto *lh_687 = buffer.data(lh + 687);
    const auto *lh_689 = buffer.data(lh + 689);
    const auto *lh_690 = buffer.data(lh + 690);
    const auto *lh_692 = buffer.data(lh + 692);
    const auto *lh_708 = buffer.data(lh + 708);
    const auto *lh_710 = buffer.data(lh + 710);
    const auto *lh_711 = buffer.data(lh + 711);
    const auto *lh_713 = buffer.data(lh + 713);
    const auto *lh_755 = buffer.data(lh + 755);
    const auto *lh_756 = buffer.data(lh + 756);
    const auto *lh_759 = buffer.data(lh + 759);
    const auto *lh_761 = buffer.data(lh + 761);
    const auto *lh_762 = buffer.data(lh + 762);
    const auto *lh_765 = buffer.data(lh + 765);
    const auto *lh_771 = buffer.data(lh + 771);
    const auto *lh_773 = buffer.data(lh + 773);
    const auto *lh_774 = buffer.data(lh + 774);
    const auto *lh_775 = buffer.data(lh + 775);
    const auto *lh_776 = buffer.data(lh + 776);
    const auto *lh_782 = buffer.data(lh + 782);
    const auto *lh_786 = buffer.data(lh + 786);
    const auto *lh_792 = buffer.data(lh + 792);
    const auto *lh_793 = buffer.data(lh + 793);
    const auto *lh_794 = buffer.data(lh + 794);
    const auto *lh_795 = buffer.data(lh + 795);
    const auto *lh_796 = buffer.data(lh + 796);
    const auto *lh_797 = buffer.data(lh + 797);
    const auto *lh_798 = buffer.data(lh + 798);

    const auto *mf0_350 = buffer.data(mf0 + 350);
    const auto *mf0_351 = buffer.data(mf0 + 351);
    const auto *mf0_355 = buffer.data(mf0 + 355);
    const auto *mf0_356 = buffer.data(mf0 + 356);
    const auto *mf0_358 = buffer.data(mf0 + 358);
    const auto *mf0_359 = buffer.data(mf0 + 359);

    const auto *mf1_350 = buffer.data(mf1 + 350);
    const auto *mf1_351 = buffer.data(mf1 + 351);
    const auto *mf1_355 = buffer.data(mf1 + 355);
    const auto *mf1_356 = buffer.data(mf1 + 356);
    const auto *mf1_358 = buffer.data(mf1 + 358);
    const auto *mf1_359 = buffer.data(mf1 + 359);

    const auto *mg_475 = buffer.data(mg + 475);
    const auto *mg_477 = buffer.data(mg + 477);
    const auto *mg_478 = buffer.data(mg + 478);
    const auto *mg_479 = buffer.data(mg + 479);
    const auto *mg_480 = buffer.data(mg + 480);
    const auto *mg_482 = buffer.data(mg + 482);
    const auto *mg_483 = buffer.data(mg + 483);
    const auto *mg_485 = buffer.data(mg + 485);
    const auto *mg_490 = buffer.data(mg + 490);
    const auto *mg_491 = buffer.data(mg + 491);
    const auto *mg_492 = buffer.data(mg + 492);
    const auto *mg_493 = buffer.data(mg + 493);
    const auto *mg_494 = buffer.data(mg + 494);
    const auto *mg_495 = buffer.data(mg + 495);
    const auto *mg_497 = buffer.data(mg + 497);
    const auto *mg_498 = buffer.data(mg + 498);
    const auto *mg_500 = buffer.data(mg + 500);
    const auto *mg_505 = buffer.data(mg + 505);
    const auto *mg_506 = buffer.data(mg + 506);
    const auto *mg_507 = buffer.data(mg + 507);
    const auto *mg_508 = buffer.data(mg + 508);
    const auto *mg_509 = buffer.data(mg + 509);
    const auto *mg_510 = buffer.data(mg + 510);
    const auto *mg_512 = buffer.data(mg + 512);
    const auto *mg_513 = buffer.data(mg + 513);
    const auto *mg_515 = buffer.data(mg + 515);
    const auto *mg_520 = buffer.data(mg + 520);
    const auto *mg_521 = buffer.data(mg + 521);
    const auto *mg_522 = buffer.data(mg + 522);
    const auto *mg_523 = buffer.data(mg + 523);
    const auto *mg_524 = buffer.data(mg + 524);
    const auto *mg_525 = buffer.data(mg + 525);
    const auto *mg_526 = buffer.data(mg + 526);
    const auto *mg_527 = buffer.data(mg + 527);
    const auto *mg_528 = buffer.data(mg + 528);
    const auto *mg_530 = buffer.data(mg + 530);
    const auto *mg_534 = buffer.data(mg + 534);
    const auto *mg_535 = buffer.data(mg + 535);
    const auto *mg_536 = buffer.data(mg + 536);
    const auto *mg_537 = buffer.data(mg + 537);
    const auto *mg_538 = buffer.data(mg + 538);
    const auto *mg_539 = buffer.data(mg + 539);
    const auto *mg_540 = buffer.data(mg + 540);
    const auto *mg_541 = buffer.data(mg + 541);
    const auto *mg_543 = buffer.data(mg + 543);
    const auto *mg_545 = buffer.data(mg + 545);
    const auto *mg_546 = buffer.data(mg + 546);
    const auto *mg_550 = buffer.data(mg + 550);
    const auto *mg_552 = buffer.data(mg + 552);
    const auto *mg_553 = buffer.data(mg + 553);
    const auto *mg_554 = buffer.data(mg + 554);
    const auto *mg_555 = buffer.data(mg + 555);
    const auto *mg_557 = buffer.data(mg + 557);
    const auto *mg_558 = buffer.data(mg + 558);
    const auto *mg_560 = buffer.data(mg + 560);
    const auto *mg_566 = buffer.data(mg + 566);
    const auto *mg_567 = buffer.data(mg + 567);
    const auto *mg_568 = buffer.data(mg + 568);
    const auto *mg_569 = buffer.data(mg + 569);

#pragma omp simd aligned(t_663, t_664, t_665, t_666, pa_x, pb_x, kh0_666, kh1_666, lg_477, \
                         lg_478, lg_479, lh_666, mg_477, mg_478, \
                         mg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_8 * lg_477[k]
                   + pb_x[k] * mg_477[k];

        t_664[k] = f_8 * lg_478[k]
                   + pb_x[k] * mg_478[k];

        t_665[k] = f_8 * lg_479[k]
                   + pb_x[k] * mg_479[k];

        t_666[k] = f_12 * kh0_666[k]
                   - f_13 * kh1_666[k]
                   + pa_x[k] * lh_666[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pa_x, pb_z, kh0_668, kh0_669, kh1_668, kh1_669, \
                         lg_355, lh_668, lh_669, mg_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_9 * lg_355[k]
                   + pb_z[k] * mg_475[k];

        t_668[k] = f_12 * kh0_668[k]
                   - f_13 * kh1_668[k]
                   + pa_x[k] * lh_668[k];

        t_669[k] = f_12 * kh0_669[k]
                   - f_13 * kh1_669[k]
                   + pa_x[k] * lh_669[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pa_x, pa_y, pb_y, kh0_399, kh0_671, kh1_399, \
                         kh1_671, lg_374, lh_525, lh_671, mg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_24 * lg_374[k]
                   + pb_y[k] * mg_479[k];

        t_671[k] = f_12 * kh0_671[k]
                   - f_13 * kh1_671[k]
                   + pa_x[k] * lh_671[k];

        t_672[k] = f_17 * kh0_399[k]
                   - f_18 * kh1_399[k]
                   + pa_y[k] * lh_525[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pa_z, pb_y, pb_z, kh0_360, kh1_360, \
                         lg_360, lg_375, lg_377, lh_507, mg_480, \
                         mg_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_9 * lg_375[k]
                   + pb_y[k] * mg_480[k];

        t_674[k] = f_24 * lg_360[k]
                   + pb_z[k] * mg_480[k];

        t_675[k] = f_22 * kh0_360[k]
                   - f_23 * kh1_360[k]
                   + pa_z[k] * lh_507[k];

        t_676[k] = f_9 * lg_377[k]
                   + pb_y[k] * mg_482[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pa_y, pa_z, pb_z, kh0_363, kh0_404, kh1_363, \
                         kh1_404, lg_363, lh_510, lh_530, mg_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_17 * kh0_404[k]
                   - f_18 * kh1_404[k]
                   + pa_y[k] * lh_530[k];

        t_678[k] = f_22 * kh0_363[k]
                   - f_23 * kh1_363[k]
                   + pa_z[k] * lh_510[k];

        t_679[k] = f_24 * lg_363[k]
                   + pb_z[k] * mg_483[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pa_y, pb_x, pb_y, kh0_408, kh1_408, \
                         lg_380, lg_490, lg_491, lh_534, mg_485, mg_490, \
                         mg_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_9 * lg_380[k]
                   + pb_y[k] * mg_485[k];

        t_681[k] = f_17 * kh0_408[k]
                   - f_18 * kh1_408[k]
                   + pa_y[k] * lh_534[k];

        t_682[k] = f_8 * lg_490[k]
                   + pb_x[k] * mg_490[k];

        t_683[k] = f_8 * lg_491[k]
                   + pb_x[k] * mg_491[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, pa_x, pb_x, kh0_687, kh1_687, lg_492, \
                         lg_493, lg_494, lh_687, mg_492, mg_493, \
                         mg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_8 * lg_492[k]
                   + pb_x[k] * mg_492[k];

        t_685[k] = f_8 * lg_493[k]
                   + pb_x[k] * mg_493[k];

        t_686[k] = f_8 * lg_494[k]
                   + pb_x[k] * mg_494[k];

        t_687[k] = f_12 * kh0_687[k]
                   - f_13 * kh1_687[k]
                   + pa_x[k] * lh_687[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, pa_x, pb_z, kh0_689, kh0_690, kh1_689, kh1_690, \
                         lg_370, lh_689, lh_690, mg_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_24 * lg_370[k]
                   + pb_z[k] * mg_490[k];

        t_689[k] = f_12 * kh0_689[k]
                   - f_13 * kh1_689[k]
                   + pa_x[k] * lh_689[k];

        t_690[k] = f_12 * kh0_690[k]
                   - f_13 * kh1_690[k]
                   + pa_x[k] * lh_690[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, pa_x, pa_y, pb_y, kh0_420, kh0_692, kh1_420, \
                         kh1_692, lg_389, lh_546, lh_692, mg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_9 * lg_389[k]
                   + pb_y[k] * mg_494[k];

        t_692[k] = f_12 * kh0_692[k]
                   - f_13 * kh1_692[k]
                   + pa_x[k] * lh_692[k];

        t_693[k] = f_12 * kh0_420[k]
                   - f_13 * kh1_420[k]
                   + pa_y[k] * lh_546[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pa_z, pb_y, pb_z, kh0_381, kh1_381, \
                         lg_375, lg_390, lg_392, lh_528, mg_495, \
                         mg_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_8 * lg_390[k]
                   + pb_y[k] * mg_495[k];

        t_695[k] = f_11 * lg_375[k]
                   + pb_z[k] * mg_495[k];

        t_696[k] = f_25 * kh0_381[k]
                   - f_26 * kh1_381[k]
                   + pa_z[k] * lh_528[k];

        t_697[k] = f_8 * lg_392[k]
                   + pb_y[k] * mg_497[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pa_y, pa_z, pb_z, kh0_384, kh0_425, kh1_384, \
                         kh1_425, lg_378, lh_531, lh_551, mg_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_12 * kh0_425[k]
                   - f_13 * kh1_425[k]
                   + pa_y[k] * lh_551[k];

        t_699[k] = f_25 * kh0_384[k]
                   - f_26 * kh1_384[k]
                   + pa_z[k] * lh_531[k];

        t_700[k] = f_11 * lg_378[k]
                   + pb_z[k] * mg_498[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pa_y, pb_x, pb_y, kh0_429, kh1_429, \
                         lg_395, lg_505, lg_506, lh_555, mg_500, mg_505, \
                         mg_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_8 * lg_395[k]
                   + pb_y[k] * mg_500[k];

        t_702[k] = f_12 * kh0_429[k]
                   - f_13 * kh1_429[k]
                   + pa_y[k] * lh_555[k];

        t_703[k] = f_8 * lg_505[k]
                   + pb_x[k] * mg_505[k];

        t_704[k] = f_8 * lg_506[k]
                   + pb_x[k] * mg_506[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, pa_x, pb_x, kh0_708, kh1_708, lg_507, \
                         lg_508, lg_509, lh_708, mg_507, mg_508, \
                         mg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_8 * lg_507[k]
                   + pb_x[k] * mg_507[k];

        t_706[k] = f_8 * lg_508[k]
                   + pb_x[k] * mg_508[k];

        t_707[k] = f_8 * lg_509[k]
                   + pb_x[k] * mg_509[k];

        t_708[k] = f_12 * kh0_708[k]
                   - f_13 * kh1_708[k]
                   + pa_x[k] * lh_708[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pa_x, pb_z, kh0_710, kh0_711, kh1_710, kh1_711, \
                         lg_385, lh_710, lh_711, mg_505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_11 * lg_385[k]
                   + pb_z[k] * mg_505[k];

        t_710[k] = f_12 * kh0_710[k]
                   - f_13 * kh1_710[k]
                   + pa_x[k] * lh_710[k];

        t_711[k] = f_12 * kh0_711[k]
                   - f_13 * kh1_711[k]
                   + pa_x[k] * lh_711[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pa_x, pa_y, pb_y, kh0_713, kh1_713, \
                         lg_404, lg_405, lh_567, lh_713, mg_509, \
                         mg_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_8 * lg_404[k]
                   + pb_y[k] * mg_509[k];

        t_713[k] = f_12 * kh0_713[k]
                   - f_13 * kh1_713[k]
                   + pa_x[k] * lh_713[k];

        t_714[k] = pa_y[k] * lh_567[k];

        t_715[k] = f_7 * lg_405[k]
                   + pb_y[k] * mg_510[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, t_720, pa_y, pb_y, lg_406, lg_407, \
                         lg_408, lh_569, lh_570, lh_572, lh_573, \
                         mg_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = pa_y[k] * lh_569[k];

        t_717[k] = f_8 * lg_406[k]
                   + pa_y[k] * lh_570[k];

        t_718[k] = f_7 * lg_407[k]
                   + pb_y[k] * mg_512[k];

        t_719[k] = pa_y[k] * lh_572[k];

        t_720[k] = f_9 * lg_408[k]
                   + pa_y[k] * lh_573[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, t_724, pa_y, pb_x, pb_y, pb_z, lg_393, lg_410, \
                         lg_520, lh_576, mg_513, mg_515, mg_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = f_19 * lg_393[k]
                   + pb_z[k] * mg_513[k];

        t_722[k] = f_7 * lg_410[k]
                   + pb_y[k] * mg_515[k];

        t_723[k] = pa_y[k] * lh_576[k];

        t_724[k] = f_8 * lg_520[k]
                   + pb_x[k] * mg_520[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, pa_y, pb_x, lg_415, lg_521, \
                         lg_522, lg_523, lh_581, lh_582, mg_521, mg_522, \
                         mg_523 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_8 * lg_521[k]
                   + pb_x[k] * mg_521[k];

        t_726[k] = f_8 * lg_522[k]
                   + pb_x[k] * mg_522[k];

        t_727[k] = f_8 * lg_523[k]
                   + pb_x[k] * mg_523[k];

        t_728[k] = pa_y[k] * lh_581[k];

        t_729[k] = f_11 * lg_415[k]
                   + pa_y[k] * lh_582[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pa_y, pb_y, pb_z, lg_400, lg_417, lg_418, \
                         lg_419, lh_584, lh_585, mg_520, mg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_19 * lg_400[k]
                   + pb_z[k] * mg_520[k];

        t_731[k] = f_9 * lg_417[k]
                   + pa_y[k] * lh_584[k];

        t_732[k] = f_8 * lg_418[k]
                   + pa_y[k] * lh_585[k];

        t_733[k] = f_7 * lg_419[k]
                   + pb_y[k] * mg_524[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, pa_y, pa_z, pb_y, pb_z, kh0_420, kh1_420, \
                         lg_405, lh_567, lh_587, mg_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pa_y[k] * lh_587[k];

        t_735[k] = f_15 * kh0_420[k]
                   - f_16 * kh1_420[k]
                   + pa_z[k] * lh_567[k];

        t_736[k] = pb_y[k] * mg_525[k];

        t_737[k] = f_14 * lg_405[k]
                   + pb_z[k] * mg_525[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pb_x, pb_y, lg_530, mf0_350, mf0_355, mf1_350, \
                         mf1_355, mg_526, mg_527, mg_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = f_3 * mf0_350[k]
                   - f_4 * mf1_350[k]
                   + pb_y[k] * mg_526[k];

        t_739[k] = pb_y[k] * mg_527[k];

        t_740[k] = f_8 * lg_530[k]
                   + f_5 * mf0_355[k]
                   - f_6 * mf1_355[k]
                   + pb_x[k] * mg_530[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pb_x, pb_y, pb_z, lg_408, lg_534, \
                         mf0_351, mf0_359, mf1_351, mf1_359, mg_528, mg_530, \
                         mg_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_5 * mf0_351[k]
                   - f_6 * mf1_351[k]
                   + pb_y[k] * mg_528[k];

        t_742[k] = f_14 * lg_408[k]
                   + pb_z[k] * mg_528[k];

        t_743[k] = pb_y[k] * mg_530[k];

        t_744[k] = f_8 * lg_534[k]
                   + f_3 * mf0_359[k]
                   - f_4 * mf1_359[k]
                   + pb_x[k] * mg_534[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, pb_x, pb_y, lg_535, lg_536, \
                         lg_537, lg_539, mg_534, mg_535, mg_536, mg_537, \
                         mg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_8 * lg_535[k]
                   + pb_x[k] * mg_535[k];

        t_746[k] = f_8 * lg_536[k]
                   + pb_x[k] * mg_536[k];

        t_747[k] = f_8 * lg_537[k]
                   + pb_x[k] * mg_537[k];

        t_748[k] = pb_y[k] * mg_534[k];

        t_749[k] = f_8 * lg_539[k]
                   + pb_x[k] * mg_539[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, pb_y, pb_z, lg_415, mf0_356, mf0_358, \
                         mf0_359, mf1_356, mf1_358, mf1_359, mg_535, mg_537, \
                         mg_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_1 * mf0_356[k]
                   - f_2 * mf1_356[k]
                   + pb_y[k] * mg_535[k];

        t_751[k] = f_14 * lg_415[k]
                   + pb_z[k] * mg_535[k];

        t_752[k] = f_5 * mf0_358[k]
                   - f_6 * mf1_358[k]
                   + pb_y[k] * mg_537[k];

        t_753[k] = f_3 * mf0_359[k]
                   - f_4 * mf1_359[k]
                   + pb_y[k] * mg_538[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, t_758, pa_x, pb_y, pb_z, kh0_755, \
                         kh1_755, lg_420, lg_540, lh_755, lh_756, mg_539, \
                         mg_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = pb_y[k] * mg_539[k];

        t_755[k] = f_12 * kh0_755[k]
                   - f_13 * kh1_755[k]
                   + pa_x[k] * lh_755[k];

        t_756[k] = f_11 * lg_540[k]
                   + pa_x[k] * lh_756[k];

        t_757[k] = f_10 * lg_420[k]
                   + pb_y[k] * mg_540[k];

        t_758[k] = pb_z[k] * mg_540[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, t_762, t_763, pa_x, pb_z, lg_543, lg_545, \
                         lg_546, lh_759, lh_761, lh_762, mg_541, \
                         mg_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_9 * lg_543[k]
                   + pa_x[k] * lh_759[k];

        t_760[k] = pb_z[k] * mg_541[k];

        t_761[k] = f_9 * lg_545[k]
                   + pa_x[k] * lh_761[k];

        t_762[k] = f_8 * lg_546[k]
                   + pa_x[k] * lh_762[k];

        t_763[k] = pb_z[k] * mg_543[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pa_x, pb_x, pb_y, pb_z, lg_425, lg_549, \
                         lg_550, lh_765, mg_545, mg_546, mg_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_10 * lg_425[k]
                   + pb_y[k] * mg_545[k];

        t_765[k] = f_8 * lg_549[k]
                   + pa_x[k] * lh_765[k];

        t_766[k] = f_7 * lg_550[k]
                   + pb_x[k] * mg_550[k];

        t_767[k] = pb_z[k] * mg_546[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, t_772, pa_x, pb_x, pb_z, lg_552, lg_553, \
                         lg_554, lh_771, mg_550, mg_552, mg_553, \
                         mg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_7 * lg_552[k]
                   + pb_x[k] * mg_552[k];

        t_769[k] = f_7 * lg_553[k]
                   + pb_x[k] * mg_553[k];

        t_770[k] = f_7 * lg_554[k]
                   + pb_x[k] * mg_554[k];

        t_771[k] = pa_x[k] * lh_771[k];

        t_772[k] = pb_z[k] * mg_550[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, t_776, t_777, t_778, pa_x, pa_z, lh_588, lh_589, \
                         lh_773, lh_774, lh_775, lh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = pa_x[k] * lh_773[k];

        t_774[k] = pa_x[k] * lh_774[k];

        t_775[k] = pa_x[k] * lh_775[k];

        t_776[k] = pa_x[k] * lh_776[k];

        t_777[k] = pa_z[k] * lh_588[k];

        t_778[k] = pa_z[k] * lh_589[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, t_782, pa_x, pa_z, pb_y, pb_z, lg_420, lg_437, \
                         lg_560, lh_591, lh_782, mg_555, mg_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_7 * lg_420[k]
                   + pb_z[k] * mg_555[k];

        t_780[k] = pa_z[k] * lh_591[k];

        t_781[k] = f_14 * lg_437[k]
                   + pb_y[k] * mg_557[k];

        t_782[k] = f_9 * lg_560[k]
                   + pa_x[k] * lh_782[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, t_786, pa_x, pa_z, pb_y, pb_z, lg_423, lg_440, \
                         lg_564, lh_594, lh_786, mg_558, mg_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = pa_z[k] * lh_594[k];

        t_784[k] = f_7 * lg_423[k]
                   + pb_z[k] * mg_558[k];

        t_785[k] = f_14 * lg_440[k]
                   + pb_y[k] * mg_560[k];

        t_786[k] = f_8 * lg_564[k]
                   + pa_x[k] * lh_786[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, t_790, t_791, pa_z, pb_x, lg_566, lg_567, \
                         lg_568, lg_569, lh_598, mg_566, mg_567, mg_568, \
                         mg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = pa_z[k] * lh_598[k];

        t_788[k] = f_7 * lg_566[k]
                   + pb_x[k] * mg_566[k];

        t_789[k] = f_7 * lg_567[k]
                   + pb_x[k] * mg_567[k];

        t_790[k] = f_7 * lg_568[k]
                   + pb_x[k] * mg_568[k];

        t_791[k] = f_7 * lg_569[k]
                   + pb_x[k] * mg_569[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, t_796, t_797, t_798, pa_x, lg_570, \
                         lh_792, lh_793, lh_794, lh_795, lh_796, lh_797, \
                         lh_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = pa_x[k] * lh_792[k];

        t_793[k] = pa_x[k] * lh_793[k];

        t_794[k] = pa_x[k] * lh_794[k];

        t_795[k] = pa_x[k] * lh_795[k];

        t_796[k] = pa_x[k] * lh_796[k];

        t_797[k] = pa_x[k] * lh_797[k];

        t_798[k] = f_11 * lg_570[k]
                   + pa_x[k] * lh_798[k];
    }
}

static auto
compute_prim_mh_electron_repulsion_0_piece6(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t lg,
                                            const size_t lh, const size_t mf0, const size_t mf1,
                                            const size_t mg, const size_t ncols,
                                            const double alpha, const double beta,
                                            const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 4.0 / p;
    const auto f_11 = 2.5 / p;
    const auto f_14 = 3.5 / p;
    const auto f_19 = 3.0 / p;
    const auto f_24 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *lg_435 = buffer.data(lg + 435);
    const auto *lg_438 = buffer.data(lg + 438);
    const auto *lg_450 = buffer.data(lg + 450);
    const auto *lg_452 = buffer.data(lg + 452);
    const auto *lg_453 = buffer.data(lg + 453);
    const auto *lg_455 = buffer.data(lg + 455);
    const auto *lg_465 = buffer.data(lg + 465);
    const auto *lg_467 = buffer.data(lg + 467);
    const auto *lg_468 = buffer.data(lg + 468);
    const auto *lg_470 = buffer.data(lg + 470);
    const auto *lg_480 = buffer.data(lg + 480);
    const auto *lg_482 = buffer.data(lg + 482);
    const auto *lg_483 = buffer.data(lg + 483);
    const auto *lg_485 = buffer.data(lg + 485);
    const auto *lg_495 = buffer.data(lg + 495);
    const auto *lg_497 = buffer.data(lg + 497);
    const auto *lg_498 = buffer.data(lg + 498);
    const auto *lg_500 = buffer.data(lg + 500);
    const auto *lg_510 = buffer.data(lg + 510);
    const auto *lg_512 = buffer.data(lg + 512);
    const auto *lg_513 = buffer.data(lg + 513);
    const auto *lg_515 = buffer.data(lg + 515);
    const auto *lg_525 = buffer.data(lg + 525);
    const auto *lg_527 = buffer.data(lg + 527);
    const auto *lg_528 = buffer.data(lg + 528);
    const auto *lg_530 = buffer.data(lg + 530);
    const auto *lg_540 = buffer.data(lg + 540);
    const auto *lg_573 = buffer.data(lg + 573);
    const auto *lg_575 = buffer.data(lg + 575);
    const auto *lg_576 = buffer.data(lg + 576);
    const auto *lg_579 = buffer.data(lg + 579);
    const auto *lg_580 = buffer.data(lg + 580);
    const auto *lg_581 = buffer.data(lg + 581);
    const auto *lg_582 = buffer.data(lg + 582);
    const auto *lg_583 = buffer.data(lg + 583);
    const auto *lg_584 = buffer.data(lg + 584);
    const auto *lg_585 = buffer.data(lg + 585);
    const auto *lg_588 = buffer.data(lg + 588);
    const auto *lg_590 = buffer.data(lg + 590);
    const auto *lg_591 = buffer.data(lg + 591);
    const auto *lg_594 = buffer.data(lg + 594);
    const auto *lg_595 = buffer.data(lg + 595);
    const auto *lg_596 = buffer.data(lg + 596);
    const auto *lg_597 = buffer.data(lg + 597);
    const auto *lg_598 = buffer.data(lg + 598);
    const auto *lg_599 = buffer.data(lg + 599);
    const auto *lg_600 = buffer.data(lg + 600);
    const auto *lg_603 = buffer.data(lg + 603);
    const auto *lg_605 = buffer.data(lg + 605);
    const auto *lg_606 = buffer.data(lg + 606);
    const auto *lg_609 = buffer.data(lg + 609);
    const auto *lg_610 = buffer.data(lg + 610);
    const auto *lg_611 = buffer.data(lg + 611);
    const auto *lg_612 = buffer.data(lg + 612);
    const auto *lg_613 = buffer.data(lg + 613);
    const auto *lg_614 = buffer.data(lg + 614);
    const auto *lg_615 = buffer.data(lg + 615);
    const auto *lg_618 = buffer.data(lg + 618);
    const auto *lg_620 = buffer.data(lg + 620);
    const auto *lg_621 = buffer.data(lg + 621);
    const auto *lg_624 = buffer.data(lg + 624);
    const auto *lg_625 = buffer.data(lg + 625);
    const auto *lg_626 = buffer.data(lg + 626);
    const auto *lg_627 = buffer.data(lg + 627);
    const auto *lg_628 = buffer.data(lg + 628);
    const auto *lg_629 = buffer.data(lg + 629);
    const auto *lg_630 = buffer.data(lg + 630);
    const auto *lg_633 = buffer.data(lg + 633);
    const auto *lg_635 = buffer.data(lg + 635);
    const auto *lg_636 = buffer.data(lg + 636);
    const auto *lg_639 = buffer.data(lg + 639);
    const auto *lg_640 = buffer.data(lg + 640);
    const auto *lg_641 = buffer.data(lg + 641);
    const auto *lg_642 = buffer.data(lg + 642);
    const auto *lg_643 = buffer.data(lg + 643);
    const auto *lg_644 = buffer.data(lg + 644);
    const auto *lg_648 = buffer.data(lg + 648);
    const auto *lg_651 = buffer.data(lg + 651);
    const auto *lg_655 = buffer.data(lg + 655);
    const auto *lg_656 = buffer.data(lg + 656);
    const auto *lg_657 = buffer.data(lg + 657);
    const auto *lg_658 = buffer.data(lg + 658);
    const auto *lg_660 = buffer.data(lg + 660);
    const auto *lg_663 = buffer.data(lg + 663);
    const auto *lg_665 = buffer.data(lg + 665);
    const auto *lg_666 = buffer.data(lg + 666);
    const auto *lg_669 = buffer.data(lg + 669);
    const auto *lg_670 = buffer.data(lg + 670);
    const auto *lg_671 = buffer.data(lg + 671);
    const auto *lg_672 = buffer.data(lg + 672);
    const auto *lg_674 = buffer.data(lg + 674);

    const auto *lh_735 = buffer.data(lh + 735);
    const auto *lh_737 = buffer.data(lh + 737);
    const auto *lh_740 = buffer.data(lh + 740);
    const auto *lh_744 = buffer.data(lh + 744);
    const auto *lh_749 = buffer.data(lh + 749);
    const auto *lh_801 = buffer.data(lh + 801);
    const auto *lh_803 = buffer.data(lh + 803);
    const auto *lh_804 = buffer.data(lh + 804);
    const auto *lh_807 = buffer.data(lh + 807);
    const auto *lh_813 = buffer.data(lh + 813);
    const auto *lh_814 = buffer.data(lh + 814);
    const auto *lh_815 = buffer.data(lh + 815);
    const auto *lh_816 = buffer.data(lh + 816);
    const auto *lh_817 = buffer.data(lh + 817);
    const auto *lh_818 = buffer.data(lh + 818);
    const auto *lh_819 = buffer.data(lh + 819);
    const auto *lh_822 = buffer.data(lh + 822);
    const auto *lh_824 = buffer.data(lh + 824);
    const auto *lh_825 = buffer.data(lh + 825);
    const auto *lh_828 = buffer.data(lh + 828);
    const auto *lh_834 = buffer.data(lh + 834);
    const auto *lh_835 = buffer.data(lh + 835);
    const auto *lh_836 = buffer.data(lh + 836);
    const auto *lh_837 = buffer.data(lh + 837);
    const auto *lh_838 = buffer.data(lh + 838);
    const auto *lh_839 = buffer.data(lh + 839);
    const auto *lh_840 = buffer.data(lh + 840);
    const auto *lh_843 = buffer.data(lh + 843);
    const auto *lh_845 = buffer.data(lh + 845);
    const auto *lh_846 = buffer.data(lh + 846);
    const auto *lh_849 = buffer.data(lh + 849);
    const auto *lh_855 = buffer.data(lh + 855);
    const auto *lh_856 = buffer.data(lh + 856);
    const auto *lh_857 = buffer.data(lh + 857);
    const auto *lh_858 = buffer.data(lh + 858);
    const auto *lh_859 = buffer.data(lh + 859);
    const auto *lh_860 = buffer.data(lh + 860);
    const auto *lh_861 = buffer.data(lh + 861);
    const auto *lh_864 = buffer.data(lh + 864);
    const auto *lh_866 = buffer.data(lh + 866);
    const auto *lh_867 = buffer.data(lh + 867);
    const auto *lh_870 = buffer.data(lh + 870);
    const auto *lh_876 = buffer.data(lh + 876);
    const auto *lh_877 = buffer.data(lh + 877);
    const auto *lh_878 = buffer.data(lh + 878);
    const auto *lh_879 = buffer.data(lh + 879);
    const auto *lh_880 = buffer.data(lh + 880);
    const auto *lh_881 = buffer.data(lh + 881);
    const auto *lh_882 = buffer.data(lh + 882);
    const auto *lh_885 = buffer.data(lh + 885);
    const auto *lh_887 = buffer.data(lh + 887);
    const auto *lh_888 = buffer.data(lh + 888);
    const auto *lh_891 = buffer.data(lh + 891);
    const auto *lh_897 = buffer.data(lh + 897);
    const auto *lh_898 = buffer.data(lh + 898);
    const auto *lh_899 = buffer.data(lh + 899);
    const auto *lh_900 = buffer.data(lh + 900);
    const auto *lh_901 = buffer.data(lh + 901);
    const auto *lh_902 = buffer.data(lh + 902);
    const auto *lh_906 = buffer.data(lh + 906);
    const auto *lh_909 = buffer.data(lh + 909);
    const auto *lh_918 = buffer.data(lh + 918);
    const auto *lh_919 = buffer.data(lh + 919);
    const auto *lh_920 = buffer.data(lh + 920);
    const auto *lh_921 = buffer.data(lh + 921);
    const auto *lh_922 = buffer.data(lh + 922);
    const auto *lh_923 = buffer.data(lh + 923);
    const auto *lh_924 = buffer.data(lh + 924);
    const auto *lh_927 = buffer.data(lh + 927);
    const auto *lh_929 = buffer.data(lh + 929);
    const auto *lh_930 = buffer.data(lh + 930);
    const auto *lh_933 = buffer.data(lh + 933);
    const auto *lh_939 = buffer.data(lh + 939);
    const auto *lh_940 = buffer.data(lh + 940);
    const auto *lh_941 = buffer.data(lh + 941);
    const auto *lh_942 = buffer.data(lh + 942);
    const auto *lh_944 = buffer.data(lh + 944);

    const auto *mf0_450 = buffer.data(mf0 + 450);
    const auto *mf0_453 = buffer.data(mf0 + 453);

    const auto *mf1_450 = buffer.data(mf1 + 450);
    const auto *mf1_453 = buffer.data(mf1 + 453);

    const auto *mg_570 = buffer.data(mg + 570);
    const auto *mg_572 = buffer.data(mg + 572);
    const auto *mg_573 = buffer.data(mg + 573);
    const auto *mg_575 = buffer.data(mg + 575);
    const auto *mg_580 = buffer.data(mg + 580);
    const auto *mg_581 = buffer.data(mg + 581);
    const auto *mg_582 = buffer.data(mg + 582);
    const auto *mg_583 = buffer.data(mg + 583);
    const auto *mg_584 = buffer.data(mg + 584);
    const auto *mg_585 = buffer.data(mg + 585);
    const auto *mg_587 = buffer.data(mg + 587);
    const auto *mg_588 = buffer.data(mg + 588);
    const auto *mg_590 = buffer.data(mg + 590);
    const auto *mg_595 = buffer.data(mg + 595);
    const auto *mg_596 = buffer.data(mg + 596);
    const auto *mg_597 = buffer.data(mg + 597);
    const auto *mg_598 = buffer.data(mg + 598);
    const auto *mg_599 = buffer.data(mg + 599);
    const auto *mg_600 = buffer.data(mg + 600);
    const auto *mg_602 = buffer.data(mg + 602);
    const auto *mg_603 = buffer.data(mg + 603);
    const auto *mg_605 = buffer.data(mg + 605);
    const auto *mg_610 = buffer.data(mg + 610);
    const auto *mg_611 = buffer.data(mg + 611);
    const auto *mg_612 = buffer.data(mg + 612);
    const auto *mg_613 = buffer.data(mg + 613);
    const auto *mg_614 = buffer.data(mg + 614);
    const auto *mg_615 = buffer.data(mg + 615);
    const auto *mg_617 = buffer.data(mg + 617);
    const auto *mg_618 = buffer.data(mg + 618);
    const auto *mg_620 = buffer.data(mg + 620);
    const auto *mg_625 = buffer.data(mg + 625);
    const auto *mg_626 = buffer.data(mg + 626);
    const auto *mg_627 = buffer.data(mg + 627);
    const auto *mg_628 = buffer.data(mg + 628);
    const auto *mg_629 = buffer.data(mg + 629);
    const auto *mg_630 = buffer.data(mg + 630);
    const auto *mg_632 = buffer.data(mg + 632);
    const auto *mg_633 = buffer.data(mg + 633);
    const auto *mg_635 = buffer.data(mg + 635);
    const auto *mg_640 = buffer.data(mg + 640);
    const auto *mg_641 = buffer.data(mg + 641);
    const auto *mg_642 = buffer.data(mg + 642);
    const auto *mg_643 = buffer.data(mg + 643);
    const auto *mg_644 = buffer.data(mg + 644);
    const auto *mg_645 = buffer.data(mg + 645);
    const auto *mg_647 = buffer.data(mg + 647);
    const auto *mg_648 = buffer.data(mg + 648);
    const auto *mg_650 = buffer.data(mg + 650);
    const auto *mg_655 = buffer.data(mg + 655);
    const auto *mg_656 = buffer.data(mg + 656);
    const auto *mg_657 = buffer.data(mg + 657);
    const auto *mg_658 = buffer.data(mg + 658);
    const auto *mg_660 = buffer.data(mg + 660);
    const auto *mg_662 = buffer.data(mg + 662);
    const auto *mg_663 = buffer.data(mg + 663);
    const auto *mg_665 = buffer.data(mg + 665);
    const auto *mg_669 = buffer.data(mg + 669);
    const auto *mg_670 = buffer.data(mg + 670);
    const auto *mg_671 = buffer.data(mg + 671);
    const auto *mg_672 = buffer.data(mg + 672);
    const auto *mg_674 = buffer.data(mg + 674);
    const auto *mg_675 = buffer.data(mg + 675);
    const auto *mg_676 = buffer.data(mg + 676);
    const auto *mg_678 = buffer.data(mg + 678);

#pragma omp simd aligned(t_799, t_800, t_801, t_802, pa_x, pb_y, pb_z, lg_435, lg_450, lg_452, \
                         lg_573, lh_801, mg_570, mg_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_799[k] = f_19 * lg_450[k]
                   + pb_y[k] * mg_570[k];

        t_800[k] = f_8 * lg_435[k]
                   + pb_z[k] * mg_570[k];

        t_801[k] = f_9 * lg_573[k]
                   + pa_x[k] * lh_801[k];

        t_802[k] = f_19 * lg_452[k]
                   + pb_y[k] * mg_572[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, t_806, pa_x, pb_y, pb_z, lg_438, lg_455, lg_575, \
                         lg_576, lh_803, lh_804, mg_573, mg_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_9 * lg_575[k]
                   + pa_x[k] * lh_803[k];

        t_804[k] = f_8 * lg_576[k]
                   + pa_x[k] * lh_804[k];

        t_805[k] = f_8 * lg_438[k]
                   + pb_z[k] * mg_573[k];

        t_806[k] = f_19 * lg_455[k]
                   + pb_y[k] * mg_575[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, t_810, pa_x, pb_x, lg_579, lg_580, lg_581, \
                         lg_582, lh_807, mg_580, mg_581, mg_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = f_8 * lg_579[k]
                   + pa_x[k] * lh_807[k];

        t_808[k] = f_7 * lg_580[k]
                   + pb_x[k] * mg_580[k];

        t_809[k] = f_7 * lg_581[k]
                   + pb_x[k] * mg_581[k];

        t_810[k] = f_7 * lg_582[k]
                   + pb_x[k] * mg_582[k];
    }

#pragma omp simd aligned(t_811, t_812, t_813, t_814, t_815, t_816, pa_x, pb_x, lg_583, lg_584, \
                         lh_813, lh_814, lh_815, lh_816, mg_583, \
                         mg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_811[k] = f_7 * lg_583[k]
                   + pb_x[k] * mg_583[k];

        t_812[k] = f_7 * lg_584[k]
                   + pb_x[k] * mg_584[k];

        t_813[k] = pa_x[k] * lh_813[k];

        t_814[k] = pa_x[k] * lh_814[k];

        t_815[k] = pa_x[k] * lh_815[k];

        t_816[k] = pa_x[k] * lh_816[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, t_820, t_821, pa_x, pb_y, pb_z, lg_450, lg_465, \
                         lg_585, lh_817, lh_818, lh_819, mg_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = pa_x[k] * lh_817[k];

        t_818[k] = pa_x[k] * lh_818[k];

        t_819[k] = f_11 * lg_585[k]
                   + pa_x[k] * lh_819[k];

        t_820[k] = f_11 * lg_465[k]
                   + pb_y[k] * mg_585[k];

        t_821[k] = f_9 * lg_450[k]
                   + pb_z[k] * mg_585[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, t_825, pa_x, pb_y, lg_467, lg_588, lg_590, \
                         lg_591, lh_822, lh_824, lh_825, mg_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = f_9 * lg_588[k]
                   + pa_x[k] * lh_822[k];

        t_823[k] = f_11 * lg_467[k]
                   + pb_y[k] * mg_587[k];

        t_824[k] = f_9 * lg_590[k]
                   + pa_x[k] * lh_824[k];

        t_825[k] = f_8 * lg_591[k]
                   + pa_x[k] * lh_825[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pa_x, pb_x, pb_y, pb_z, lg_453, lg_470, \
                         lg_594, lg_595, lh_828, mg_588, mg_590, \
                         mg_595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_9 * lg_453[k]
                   + pb_z[k] * mg_588[k];

        t_827[k] = f_11 * lg_470[k]
                   + pb_y[k] * mg_590[k];

        t_828[k] = f_8 * lg_594[k]
                   + pa_x[k] * lh_828[k];

        t_829[k] = f_7 * lg_595[k]
                   + pb_x[k] * mg_595[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, pa_x, pb_x, lg_596, lg_597, \
                         lg_598, lg_599, lh_834, mg_596, mg_597, mg_598, \
                         mg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_7 * lg_596[k]
                   + pb_x[k] * mg_596[k];

        t_831[k] = f_7 * lg_597[k]
                   + pb_x[k] * mg_597[k];

        t_832[k] = f_7 * lg_598[k]
                   + pb_x[k] * mg_598[k];

        t_833[k] = f_7 * lg_599[k]
                   + pb_x[k] * mg_599[k];

        t_834[k] = pa_x[k] * lh_834[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, t_838, t_839, t_840, pa_x, lg_600, lh_835, \
                         lh_836, lh_837, lh_838, lh_839, lh_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = pa_x[k] * lh_835[k];

        t_836[k] = pa_x[k] * lh_836[k];

        t_837[k] = pa_x[k] * lh_837[k];

        t_838[k] = pa_x[k] * lh_838[k];

        t_839[k] = pa_x[k] * lh_839[k];

        t_840[k] = f_11 * lg_600[k]
                   + pa_x[k] * lh_840[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, t_844, pa_x, pb_y, pb_z, lg_465, lg_480, lg_482, \
                         lg_603, lh_843, mg_600, mg_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_24 * lg_480[k]
                   + pb_y[k] * mg_600[k];

        t_842[k] = f_24 * lg_465[k]
                   + pb_z[k] * mg_600[k];

        t_843[k] = f_9 * lg_603[k]
                   + pa_x[k] * lh_843[k];

        t_844[k] = f_24 * lg_482[k]
                   + pb_y[k] * mg_602[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, pa_x, pb_y, pb_z, lg_468, lg_485, lg_605, \
                         lg_606, lh_845, lh_846, mg_603, mg_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_9 * lg_605[k]
                   + pa_x[k] * lh_845[k];

        t_846[k] = f_8 * lg_606[k]
                   + pa_x[k] * lh_846[k];

        t_847[k] = f_24 * lg_468[k]
                   + pb_z[k] * mg_603[k];

        t_848[k] = f_24 * lg_485[k]
                   + pb_y[k] * mg_605[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, t_852, pa_x, pb_x, lg_609, lg_610, lg_611, \
                         lg_612, lh_849, mg_610, mg_611, mg_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_8 * lg_609[k]
                   + pa_x[k] * lh_849[k];

        t_850[k] = f_7 * lg_610[k]
                   + pb_x[k] * mg_610[k];

        t_851[k] = f_7 * lg_611[k]
                   + pb_x[k] * mg_611[k];

        t_852[k] = f_7 * lg_612[k]
                   + pb_x[k] * mg_612[k];
    }

#pragma omp simd aligned(t_853, t_854, t_855, t_856, t_857, t_858, pa_x, pb_x, lg_613, lg_614, \
                         lh_855, lh_856, lh_857, lh_858, mg_613, \
                         mg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_853[k] = f_7 * lg_613[k]
                   + pb_x[k] * mg_613[k];

        t_854[k] = f_7 * lg_614[k]
                   + pb_x[k] * mg_614[k];

        t_855[k] = pa_x[k] * lh_855[k];

        t_856[k] = pa_x[k] * lh_856[k];

        t_857[k] = pa_x[k] * lh_857[k];

        t_858[k] = pa_x[k] * lh_858[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, t_862, t_863, pa_x, pb_y, pb_z, lg_480, lg_495, \
                         lg_615, lh_859, lh_860, lh_861, mg_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = pa_x[k] * lh_859[k];

        t_860[k] = pa_x[k] * lh_860[k];

        t_861[k] = f_11 * lg_615[k]
                   + pa_x[k] * lh_861[k];

        t_862[k] = f_9 * lg_495[k]
                   + pb_y[k] * mg_615[k];

        t_863[k] = f_11 * lg_480[k]
                   + pb_z[k] * mg_615[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, t_867, pa_x, pb_y, lg_497, lg_618, lg_620, \
                         lg_621, lh_864, lh_866, lh_867, mg_617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = f_9 * lg_618[k]
                   + pa_x[k] * lh_864[k];

        t_865[k] = f_9 * lg_497[k]
                   + pb_y[k] * mg_617[k];

        t_866[k] = f_9 * lg_620[k]
                   + pa_x[k] * lh_866[k];

        t_867[k] = f_8 * lg_621[k]
                   + pa_x[k] * lh_867[k];
    }

#pragma omp simd aligned(t_868, t_869, t_870, t_871, pa_x, pb_x, pb_y, pb_z, lg_483, lg_500, \
                         lg_624, lg_625, lh_870, mg_618, mg_620, \
                         mg_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_868[k] = f_11 * lg_483[k]
                   + pb_z[k] * mg_618[k];

        t_869[k] = f_9 * lg_500[k]
                   + pb_y[k] * mg_620[k];

        t_870[k] = f_8 * lg_624[k]
                   + pa_x[k] * lh_870[k];

        t_871[k] = f_7 * lg_625[k]
                   + pb_x[k] * mg_625[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, t_875, t_876, pa_x, pb_x, lg_626, lg_627, \
                         lg_628, lg_629, lh_876, mg_626, mg_627, mg_628, \
                         mg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_7 * lg_626[k]
                   + pb_x[k] * mg_626[k];

        t_873[k] = f_7 * lg_627[k]
                   + pb_x[k] * mg_627[k];

        t_874[k] = f_7 * lg_628[k]
                   + pb_x[k] * mg_628[k];

        t_875[k] = f_7 * lg_629[k]
                   + pb_x[k] * mg_629[k];

        t_876[k] = pa_x[k] * lh_876[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, t_880, t_881, t_882, pa_x, lg_630, lh_877, \
                         lh_878, lh_879, lh_880, lh_881, lh_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = pa_x[k] * lh_877[k];

        t_878[k] = pa_x[k] * lh_878[k];

        t_879[k] = pa_x[k] * lh_879[k];

        t_880[k] = pa_x[k] * lh_880[k];

        t_881[k] = pa_x[k] * lh_881[k];

        t_882[k] = f_11 * lg_630[k]
                   + pa_x[k] * lh_882[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, pa_x, pb_y, pb_z, lg_495, lg_510, lg_512, \
                         lg_633, lh_885, mg_630, mg_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = f_8 * lg_510[k]
                   + pb_y[k] * mg_630[k];

        t_884[k] = f_19 * lg_495[k]
                   + pb_z[k] * mg_630[k];

        t_885[k] = f_9 * lg_633[k]
                   + pa_x[k] * lh_885[k];

        t_886[k] = f_8 * lg_512[k]
                   + pb_y[k] * mg_632[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, t_890, pa_x, pb_y, pb_z, lg_498, lg_515, lg_635, \
                         lg_636, lh_887, lh_888, mg_633, mg_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = f_9 * lg_635[k]
                   + pa_x[k] * lh_887[k];

        t_888[k] = f_8 * lg_636[k]
                   + pa_x[k] * lh_888[k];

        t_889[k] = f_19 * lg_498[k]
                   + pb_z[k] * mg_633[k];

        t_890[k] = f_8 * lg_515[k]
                   + pb_y[k] * mg_635[k];
    }

#pragma omp simd aligned(t_891, t_892, t_893, t_894, pa_x, pb_x, lg_639, lg_640, lg_641, \
                         lg_642, lh_891, mg_640, mg_641, mg_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_891[k] = f_8 * lg_639[k]
                   + pa_x[k] * lh_891[k];

        t_892[k] = f_7 * lg_640[k]
                   + pb_x[k] * mg_640[k];

        t_893[k] = f_7 * lg_641[k]
                   + pb_x[k] * mg_641[k];

        t_894[k] = f_7 * lg_642[k]
                   + pb_x[k] * mg_642[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, t_899, t_900, pa_x, pb_x, lg_643, lg_644, \
                         lh_897, lh_898, lh_899, lh_900, mg_643, \
                         mg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_7 * lg_643[k]
                   + pb_x[k] * mg_643[k];

        t_896[k] = f_7 * lg_644[k]
                   + pb_x[k] * mg_644[k];

        t_897[k] = pa_x[k] * lh_897[k];

        t_898[k] = pa_x[k] * lh_898[k];

        t_899[k] = pa_x[k] * lh_899[k];

        t_900[k] = pa_x[k] * lh_900[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, t_904, t_905, pa_x, pa_y, pb_y, lg_525, lh_735, \
                         lh_737, lh_901, lh_902, mg_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = pa_x[k] * lh_901[k];

        t_902[k] = pa_x[k] * lh_902[k];

        t_903[k] = pa_y[k] * lh_735[k];

        t_904[k] = f_7 * lg_525[k]
                   + pb_y[k] * mg_645[k];

        t_905[k] = pa_y[k] * lh_737[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, t_909, pa_x, pa_y, pb_y, lg_527, lg_648, lg_651, \
                         lh_740, lh_906, lh_909, mg_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_9 * lg_648[k]
                   + pa_x[k] * lh_906[k];

        t_907[k] = f_7 * lg_527[k]
                   + pb_y[k] * mg_647[k];

        t_908[k] = pa_y[k] * lh_740[k];

        t_909[k] = f_8 * lg_651[k]
                   + pa_x[k] * lh_909[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, pa_y, pb_x, pb_y, pb_z, lg_513, lg_530, \
                         lg_655, lh_744, mg_648, mg_650, mg_655 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = f_14 * lg_513[k]
                   + pb_z[k] * mg_648[k];

        t_911[k] = f_7 * lg_530[k]
                   + pb_y[k] * mg_650[k];

        t_912[k] = pa_y[k] * lh_744[k];

        t_913[k] = f_7 * lg_655[k]
                   + pb_x[k] * mg_655[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, t_917, t_918, pa_x, pa_y, pb_x, lg_656, lg_657, \
                         lg_658, lh_749, lh_918, mg_656, mg_657, \
                         mg_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_7 * lg_656[k]
                   + pb_x[k] * mg_656[k];

        t_915[k] = f_7 * lg_657[k]
                   + pb_x[k] * mg_657[k];

        t_916[k] = f_7 * lg_658[k]
                   + pb_x[k] * mg_658[k];

        t_917[k] = pa_y[k] * lh_749[k];

        t_918[k] = pa_x[k] * lh_918[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, t_922, t_923, t_924, pa_x, lg_660, lh_919, \
                         lh_920, lh_921, lh_922, lh_923, lh_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = pa_x[k] * lh_919[k];

        t_920[k] = pa_x[k] * lh_920[k];

        t_921[k] = pa_x[k] * lh_921[k];

        t_922[k] = pa_x[k] * lh_922[k];

        t_923[k] = pa_x[k] * lh_923[k];

        t_924[k] = f_11 * lg_660[k]
                   + pa_x[k] * lh_924[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, t_929, pa_x, pb_y, pb_z, lg_525, lg_663, \
                         lg_665, lh_927, lh_929, mg_660, mg_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = pb_y[k] * mg_660[k];

        t_926[k] = f_10 * lg_525[k]
                   + pb_z[k] * mg_660[k];

        t_927[k] = f_9 * lg_663[k]
                   + pa_x[k] * lh_927[k];

        t_928[k] = pb_y[k] * mg_662[k];

        t_929[k] = f_9 * lg_665[k]
                   + pa_x[k] * lh_929[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, pa_x, pb_y, pb_z, lg_528, lg_666, lg_669, \
                         lh_930, lh_933, mg_663, mg_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_8 * lg_666[k]
                   + pa_x[k] * lh_930[k];

        t_931[k] = f_10 * lg_528[k]
                   + pb_z[k] * mg_663[k];

        t_932[k] = pb_y[k] * mg_665[k];

        t_933[k] = f_8 * lg_669[k]
                   + pa_x[k] * lh_933[k];
    }

#pragma omp simd aligned(t_934, t_935, t_936, t_937, t_938, pb_x, pb_y, lg_670, lg_671, \
                         lg_672, lg_674, mg_669, mg_670, mg_671, mg_672, \
                         mg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_7 * lg_670[k]
                   + pb_x[k] * mg_670[k];

        t_935[k] = f_7 * lg_671[k]
                   + pb_x[k] * mg_671[k];

        t_936[k] = f_7 * lg_672[k]
                   + pb_x[k] * mg_672[k];

        t_937[k] = pb_y[k] * mg_669[k];

        t_938[k] = f_7 * lg_674[k]
                   + pb_x[k] * mg_674[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, t_942, t_943, t_944, pa_x, pb_y, lh_939, lh_940, \
                         lh_941, lh_942, lh_944, mg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = pa_x[k] * lh_939[k];

        t_940[k] = pa_x[k] * lh_940[k];

        t_941[k] = pa_x[k] * lh_941[k];

        t_942[k] = pa_x[k] * lh_942[k];

        t_943[k] = pb_y[k] * mg_674[k];

        t_944[k] = pa_x[k] * lh_944[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, t_948, t_949, pb_x, pb_y, pb_z, lg_540, mf0_450, \
                         mf0_453, mf1_450, mf1_453, mg_675, mg_676, \
                         mg_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = f_1 * mf0_450[k]
                   - f_2 * mf1_450[k]
                   + pb_x[k] * mg_675[k];

        t_946[k] = f_0 * lg_540[k]
                   + pb_y[k] * mg_675[k];

        t_947[k] = pb_z[k] * mg_675[k];

        t_948[k] = f_5 * mf0_453[k]
                   - f_6 * mf1_453[k]
                   + pb_x[k] * mg_678[k];

        t_949[k] = pb_z[k] * mg_676[k];
    }
}

static auto
compute_prim_mh_electron_repulsion_0_piece7(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kh0,
                                            const size_t kh1, const size_t lg, const size_t lh,
                                            const size_t mf0, const size_t mf1, const size_t mg,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 4.0 / p;
    const auto f_11 = 2.5 / p;
    const auto f_12 = 0.5 / alpha;
    const auto f_13 = 0.5 * beta / (alpha * p);
    const auto f_14 = 3.5 / p;
    const auto f_15 = 3.0 / alpha;
    const auto f_16 = 3.0 * beta / (alpha * p);
    const auto f_17 = 1.0 / alpha;
    const auto f_18 = beta / (alpha * p);
    const auto f_19 = 3.0 / p;
    const auto f_20 = 2.5 / alpha;
    const auto f_21 = 2.5 * beta / (alpha * p);
    const auto f_22 = 1.5 / alpha;
    const auto f_23 = 1.5 * beta / (alpha * p);
    const auto f_24 = 2.0 / p;
    const auto f_25 = 2.0 / alpha;
    const auto f_26 = 2.0 * beta / (alpha * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kh0_603 = buffer.data(kh0 + 603);
    const auto *kh0_624 = buffer.data(kh0 + 624);
    const auto *kh0_645 = buffer.data(kh0 + 645);
    const auto *kh0_650 = buffer.data(kh0 + 650);
    const auto *kh0_666 = buffer.data(kh0 + 666);
    const auto *kh0_671 = buffer.data(kh0 + 671);
    const auto *kh0_692 = buffer.data(kh0 + 692);
    const auto *kh0_713 = buffer.data(kh0 + 713);

    const auto *kh1_603 = buffer.data(kh1 + 603);
    const auto *kh1_624 = buffer.data(kh1 + 624);
    const auto *kh1_645 = buffer.data(kh1 + 645);
    const auto *kh1_650 = buffer.data(kh1 + 650);
    const auto *kh1_666 = buffer.data(kh1 + 666);
    const auto *kh1_671 = buffer.data(kh1 + 671);
    const auto *kh1_692 = buffer.data(kh1 + 692);
    const auto *kh1_713 = buffer.data(kh1 + 713);

    const auto *lg_540 = buffer.data(lg + 540);
    const auto *lg_542 = buffer.data(lg + 542);
    const auto *lg_543 = buffer.data(lg + 543);
    const auto *lg_545 = buffer.data(lg + 545);
    const auto *lg_550 = buffer.data(lg + 550);
    const auto *lg_551 = buffer.data(lg + 551);
    const auto *lg_552 = buffer.data(lg + 552);
    const auto *lg_554 = buffer.data(lg + 554);
    const auto *lg_555 = buffer.data(lg + 555);
    const auto *lg_557 = buffer.data(lg + 557);
    const auto *lg_558 = buffer.data(lg + 558);
    const auto *lg_560 = buffer.data(lg + 560);
    const auto *lg_565 = buffer.data(lg + 565);
    const auto *lg_569 = buffer.data(lg + 569);
    const auto *lg_570 = buffer.data(lg + 570);
    const auto *lg_572 = buffer.data(lg + 572);
    const auto *lg_573 = buffer.data(lg + 573);
    const auto *lg_575 = buffer.data(lg + 575);
    const auto *lg_580 = buffer.data(lg + 580);
    const auto *lg_582 = buffer.data(lg + 582);
    const auto *lg_583 = buffer.data(lg + 583);
    const auto *lg_584 = buffer.data(lg + 584);
    const auto *lg_585 = buffer.data(lg + 585);
    const auto *lg_587 = buffer.data(lg + 587);
    const auto *lg_588 = buffer.data(lg + 588);
    const auto *lg_590 = buffer.data(lg + 590);
    const auto *lg_595 = buffer.data(lg + 595);
    const auto *lg_597 = buffer.data(lg + 597);
    const auto *lg_598 = buffer.data(lg + 598);
    const auto *lg_599 = buffer.data(lg + 599);
    const auto *lg_600 = buffer.data(lg + 600);
    const auto *lg_602 = buffer.data(lg + 602);
    const auto *lg_603 = buffer.data(lg + 603);
    const auto *lg_605 = buffer.data(lg + 605);
    const auto *lg_610 = buffer.data(lg + 610);
    const auto *lg_612 = buffer.data(lg + 612);
    const auto *lg_613 = buffer.data(lg + 613);
    const auto *lg_614 = buffer.data(lg + 614);
    const auto *lg_615 = buffer.data(lg + 615);
    const auto *lg_617 = buffer.data(lg + 617);
    const auto *lg_618 = buffer.data(lg + 618);
    const auto *lg_620 = buffer.data(lg + 620);
    const auto *lg_627 = buffer.data(lg + 627);
    const auto *lg_628 = buffer.data(lg + 628);
    const auto *lg_629 = buffer.data(lg + 629);
    const auto *lg_630 = buffer.data(lg + 630);
    const auto *lg_632 = buffer.data(lg + 632);
    const auto *lg_635 = buffer.data(lg + 635);

    const auto *lh_756 = buffer.data(lh + 756);
    const auto *lh_757 = buffer.data(lh + 757);
    const auto *lh_759 = buffer.data(lh + 759);
    const auto *lh_761 = buffer.data(lh + 761);
    const auto *lh_762 = buffer.data(lh + 762);
    const auto *lh_765 = buffer.data(lh + 765);
    const auto *lh_771 = buffer.data(lh + 771);
    const auto *lh_773 = buffer.data(lh + 773);
    const auto *lh_774 = buffer.data(lh + 774);
    const auto *lh_776 = buffer.data(lh + 776);
    const auto *lh_792 = buffer.data(lh + 792);
    const auto *lh_813 = buffer.data(lh + 813);
    const auto *lh_818 = buffer.data(lh + 818);
    const auto *lh_834 = buffer.data(lh + 834);
    const auto *lh_839 = buffer.data(lh + 839);
    const auto *lh_855 = buffer.data(lh + 855);
    const auto *lh_860 = buffer.data(lh + 860);
    const auto *lh_881 = buffer.data(lh + 881);

    const auto *mf0_455 = buffer.data(mf0 + 455);
    const auto *mf0_456 = buffer.data(mf0 + 456);
    const auto *mf0_457 = buffer.data(mf0 + 457);
    const auto *mf0_459 = buffer.data(mf0 + 459);
    const auto *mf0_470 = buffer.data(mf0 + 470);
    const auto *mf0_473 = buffer.data(mf0 + 473);
    const auto *mf0_475 = buffer.data(mf0 + 475);
    const auto *mf0_476 = buffer.data(mf0 + 476);
    const auto *mf0_478 = buffer.data(mf0 + 478);
    const auto *mf0_479 = buffer.data(mf0 + 479);
    const auto *mf0_480 = buffer.data(mf0 + 480);
    const auto *mf0_483 = buffer.data(mf0 + 483);
    const auto *mf0_485 = buffer.data(mf0 + 485);
    const auto *mf0_486 = buffer.data(mf0 + 486);
    const auto *mf0_488 = buffer.data(mf0 + 488);
    const auto *mf0_489 = buffer.data(mf0 + 489);
    const auto *mf0_490 = buffer.data(mf0 + 490);
    const auto *mf0_493 = buffer.data(mf0 + 493);
    const auto *mf0_495 = buffer.data(mf0 + 495);
    const auto *mf0_496 = buffer.data(mf0 + 496);
    const auto *mf0_498 = buffer.data(mf0 + 498);
    const auto *mf0_499 = buffer.data(mf0 + 499);
    const auto *mf0_500 = buffer.data(mf0 + 500);
    const auto *mf0_503 = buffer.data(mf0 + 503);
    const auto *mf0_505 = buffer.data(mf0 + 505);
    const auto *mf0_506 = buffer.data(mf0 + 506);
    const auto *mf0_508 = buffer.data(mf0 + 508);
    const auto *mf0_509 = buffer.data(mf0 + 509);
    const auto *mf0_510 = buffer.data(mf0 + 510);
    const auto *mf0_513 = buffer.data(mf0 + 513);
    const auto *mf0_515 = buffer.data(mf0 + 515);
    const auto *mf0_516 = buffer.data(mf0 + 516);

    const auto *mf1_455 = buffer.data(mf1 + 455);
    const auto *mf1_456 = buffer.data(mf1 + 456);
    const auto *mf1_457 = buffer.data(mf1 + 457);
    const auto *mf1_459 = buffer.data(mf1 + 459);
    const auto *mf1_470 = buffer.data(mf1 + 470);
    const auto *mf1_473 = buffer.data(mf1 + 473);
    const auto *mf1_475 = buffer.data(mf1 + 475);
    const auto *mf1_476 = buffer.data(mf1 + 476);
    const auto *mf1_478 = buffer.data(mf1 + 478);
    const auto *mf1_479 = buffer.data(mf1 + 479);
    const auto *mf1_480 = buffer.data(mf1 + 480);
    const auto *mf1_483 = buffer.data(mf1 + 483);
    const auto *mf1_485 = buffer.data(mf1 + 485);
    const auto *mf1_486 = buffer.data(mf1 + 486);
    const auto *mf1_488 = buffer.data(mf1 + 488);
    const auto *mf1_489 = buffer.data(mf1 + 489);
    const auto *mf1_490 = buffer.data(mf1 + 490);
    const auto *mf1_493 = buffer.data(mf1 + 493);
    const auto *mf1_495 = buffer.data(mf1 + 495);
    const auto *mf1_496 = buffer.data(mf1 + 496);
    const auto *mf1_498 = buffer.data(mf1 + 498);
    const auto *mf1_499 = buffer.data(mf1 + 499);
    const auto *mf1_500 = buffer.data(mf1 + 500);
    const auto *mf1_503 = buffer.data(mf1 + 503);
    const auto *mf1_505 = buffer.data(mf1 + 505);
    const auto *mf1_506 = buffer.data(mf1 + 506);
    const auto *mf1_508 = buffer.data(mf1 + 508);
    const auto *mf1_509 = buffer.data(mf1 + 509);
    const auto *mf1_510 = buffer.data(mf1 + 510);
    const auto *mf1_513 = buffer.data(mf1 + 513);
    const auto *mf1_515 = buffer.data(mf1 + 515);
    const auto *mf1_516 = buffer.data(mf1 + 516);

    const auto *mg_678 = buffer.data(mg + 678);
    const auto *mg_680 = buffer.data(mg + 680);
    const auto *mg_681 = buffer.data(mg + 681);
    const auto *mg_684 = buffer.data(mg + 684);
    const auto *mg_685 = buffer.data(mg + 685);
    const auto *mg_686 = buffer.data(mg + 686);
    const auto *mg_687 = buffer.data(mg + 687);
    const auto *mg_688 = buffer.data(mg + 688);
    const auto *mg_689 = buffer.data(mg + 689);
    const auto *mg_690 = buffer.data(mg + 690);
    const auto *mg_692 = buffer.data(mg + 692);
    const auto *mg_693 = buffer.data(mg + 693);
    const auto *mg_695 = buffer.data(mg + 695);
    const auto *mg_700 = buffer.data(mg + 700);
    const auto *mg_701 = buffer.data(mg + 701);
    const auto *mg_702 = buffer.data(mg + 702);
    const auto *mg_703 = buffer.data(mg + 703);
    const auto *mg_704 = buffer.data(mg + 704);
    const auto *mg_705 = buffer.data(mg + 705);
    const auto *mg_707 = buffer.data(mg + 707);
    const auto *mg_708 = buffer.data(mg + 708);
    const auto *mg_710 = buffer.data(mg + 710);
    const auto *mg_711 = buffer.data(mg + 711);
    const auto *mg_714 = buffer.data(mg + 714);
    const auto *mg_715 = buffer.data(mg + 715);
    const auto *mg_716 = buffer.data(mg + 716);
    const auto *mg_717 = buffer.data(mg + 717);
    const auto *mg_718 = buffer.data(mg + 718);
    const auto *mg_719 = buffer.data(mg + 719);
    const auto *mg_720 = buffer.data(mg + 720);
    const auto *mg_722 = buffer.data(mg + 722);
    const auto *mg_723 = buffer.data(mg + 723);
    const auto *mg_725 = buffer.data(mg + 725);
    const auto *mg_726 = buffer.data(mg + 726);
    const auto *mg_729 = buffer.data(mg + 729);
    const auto *mg_730 = buffer.data(mg + 730);
    const auto *mg_731 = buffer.data(mg + 731);
    const auto *mg_732 = buffer.data(mg + 732);
    const auto *mg_733 = buffer.data(mg + 733);
    const auto *mg_734 = buffer.data(mg + 734);
    const auto *mg_735 = buffer.data(mg + 735);
    const auto *mg_737 = buffer.data(mg + 737);
    const auto *mg_738 = buffer.data(mg + 738);
    const auto *mg_740 = buffer.data(mg + 740);
    const auto *mg_741 = buffer.data(mg + 741);
    const auto *mg_744 = buffer.data(mg + 744);
    const auto *mg_745 = buffer.data(mg + 745);
    const auto *mg_746 = buffer.data(mg + 746);
    const auto *mg_747 = buffer.data(mg + 747);
    const auto *mg_748 = buffer.data(mg + 748);
    const auto *mg_749 = buffer.data(mg + 749);
    const auto *mg_750 = buffer.data(mg + 750);
    const auto *mg_752 = buffer.data(mg + 752);
    const auto *mg_753 = buffer.data(mg + 753);
    const auto *mg_755 = buffer.data(mg + 755);
    const auto *mg_756 = buffer.data(mg + 756);
    const auto *mg_759 = buffer.data(mg + 759);
    const auto *mg_760 = buffer.data(mg + 760);
    const auto *mg_761 = buffer.data(mg + 761);
    const auto *mg_762 = buffer.data(mg + 762);
    const auto *mg_763 = buffer.data(mg + 763);
    const auto *mg_764 = buffer.data(mg + 764);
    const auto *mg_765 = buffer.data(mg + 765);
    const auto *mg_767 = buffer.data(mg + 767);
    const auto *mg_768 = buffer.data(mg + 768);
    const auto *mg_770 = buffer.data(mg + 770);
    const auto *mg_771 = buffer.data(mg + 771);

#pragma omp simd aligned(t_950, t_951, t_952, t_953, pb_x, pb_y, pb_z, lg_545, mf0_455, \
                         mf0_456, mf1_455, mf1_456, mg_678, mg_680, \
                         mg_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = f_5 * mf0_455[k]
                   - f_6 * mf1_455[k]
                   + pb_x[k] * mg_680[k];

        t_951[k] = f_3 * mf0_456[k]
                   - f_4 * mf1_456[k]
                   + pb_x[k] * mg_681[k];

        t_952[k] = pb_z[k] * mg_678[k];

        t_953[k] = f_0 * lg_545[k]
                   + pb_y[k] * mg_680[k];
    }

#pragma omp simd aligned(t_954, t_955, t_956, t_957, t_958, t_959, pb_x, mf0_459, mf1_459, \
                         mg_684, mg_685, mg_686, mg_687, mg_688, \
                         mg_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_954[k] = f_3 * mf0_459[k]
                   - f_4 * mf1_459[k]
                   + pb_x[k] * mg_684[k];

        t_955[k] = pb_x[k] * mg_685[k];

        t_956[k] = pb_x[k] * mg_686[k];

        t_957[k] = pb_x[k] * mg_687[k];

        t_958[k] = pb_x[k] * mg_688[k];

        t_959[k] = pb_x[k] * mg_689[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, t_963, pb_y, pb_z, lg_550, mf0_456, mf0_457, \
                         mf1_456, mf1_457, mg_685, mg_686, mg_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = f_0 * lg_550[k]
                   + f_1 * mf0_456[k]
                   - f_2 * mf1_456[k]
                   + pb_y[k] * mg_685[k];

        t_961[k] = pb_z[k] * mg_685[k];

        t_962[k] = f_3 * mf0_456[k]
                   - f_4 * mf1_456[k]
                   + pb_z[k] * mg_686[k];

        t_963[k] = f_5 * mf0_457[k]
                   - f_6 * mf1_457[k]
                   + pb_z[k] * mg_687[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, t_968, pa_z, pb_y, pb_z, lg_540, lg_554, \
                         lh_756, lh_757, mf0_459, mf1_459, mg_689, \
                         mg_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = f_0 * lg_554[k]
                   + pb_y[k] * mg_689[k];

        t_965[k] = f_1 * mf0_459[k]
                   - f_2 * mf1_459[k]
                   + pb_z[k] * mg_689[k];

        t_966[k] = pa_z[k] * lh_756[k];

        t_967[k] = pa_z[k] * lh_757[k];

        t_968[k] = f_7 * lg_540[k]
                   + pb_z[k] * mg_690[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, t_972, t_973, pa_z, pb_y, pb_z, lg_542, lg_543, \
                         lg_557, lh_759, lh_761, lh_762, mg_692, \
                         mg_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = pa_z[k] * lh_759[k];

        t_970[k] = f_10 * lg_557[k]
                   + pb_y[k] * mg_692[k];

        t_971[k] = f_8 * lg_542[k]
                   + pa_z[k] * lh_761[k];

        t_972[k] = pa_z[k] * lh_762[k];

        t_973[k] = f_7 * lg_543[k]
                   + pb_z[k] * mg_693[k];
    }

#pragma omp simd aligned(t_974, t_975, t_976, t_977, t_978, pa_z, pb_x, pb_y, lg_545, lg_560, \
                         lh_765, mg_695, mg_700, mg_701, mg_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_974[k] = f_10 * lg_560[k]
                   + pb_y[k] * mg_695[k];

        t_975[k] = f_9 * lg_545[k]
                   + pa_z[k] * lh_765[k];

        t_976[k] = pb_x[k] * mg_700[k];

        t_977[k] = pb_x[k] * mg_701[k];

        t_978[k] = pb_x[k] * mg_702[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, t_982, t_983, pa_z, pb_x, pb_z, lg_550, lg_551, \
                         lh_771, lh_773, mg_700, mg_703, mg_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = pb_x[k] * mg_703[k];

        t_980[k] = pb_x[k] * mg_704[k];

        t_981[k] = pa_z[k] * lh_771[k];

        t_982[k] = f_7 * lg_550[k]
                   + pb_z[k] * mg_700[k];

        t_983[k] = f_8 * lg_551[k]
                   + pa_z[k] * lh_773[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, t_987, pa_z, pb_x, pb_y, lg_552, lg_554, lg_569, \
                         lh_774, lh_776, mf0_470, mf1_470, mg_704, \
                         mg_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_9 * lg_552[k]
                   + pa_z[k] * lh_774[k];

        t_985[k] = f_10 * lg_569[k]
                   + pb_y[k] * mg_704[k];

        t_986[k] = f_11 * lg_554[k]
                   + pa_z[k] * lh_776[k];

        t_987[k] = f_1 * mf0_470[k]
                   - f_2 * mf1_470[k]
                   + pb_x[k] * mg_705[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, t_991, pb_x, pb_y, pb_z, lg_555, lg_570, lg_572, \
                         mf0_473, mf1_473, mg_705, mg_707, mg_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_14 * lg_570[k]
                   + pb_y[k] * mg_705[k];

        t_989[k] = f_8 * lg_555[k]
                   + pb_z[k] * mg_705[k];

        t_990[k] = f_5 * mf0_473[k]
                   - f_6 * mf1_473[k]
                   + pb_x[k] * mg_708[k];

        t_991[k] = f_14 * lg_572[k]
                   + pb_y[k] * mg_707[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, pb_x, pb_y, pb_z, lg_558, lg_575, \
                         mf0_475, mf0_476, mf1_475, mf1_476, mg_708, mg_710, \
                         mg_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_5 * mf0_475[k]
                   - f_6 * mf1_475[k]
                   + pb_x[k] * mg_710[k];

        t_993[k] = f_3 * mf0_476[k]
                   - f_4 * mf1_476[k]
                   + pb_x[k] * mg_711[k];

        t_994[k] = f_8 * lg_558[k]
                   + pb_z[k] * mg_708[k];

        t_995[k] = f_14 * lg_575[k]
                   + pb_y[k] * mg_710[k];
    }

#pragma omp simd aligned(t_996, t_997, t_998, t_999, t_1000, t_1001, pb_x, mf0_479, mf1_479, \
                         mg_714, mg_715, mg_716, mg_717, mg_718, \
                         mg_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_996[k] = f_3 * mf0_479[k]
                   - f_4 * mf1_479[k]
                   + pb_x[k] * mg_714[k];

        t_997[k] = pb_x[k] * mg_715[k];

        t_998[k] = pb_x[k] * mg_716[k];

        t_999[k] = pb_x[k] * mg_717[k];

        t_1000[k] = pb_x[k] * mg_718[k];

        t_1001[k] = pb_x[k] * mg_719[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, pa_z, pb_y, pb_z, kh0_603, kh1_603, lg_565, \
                         lg_582, lh_792, mf0_478, mf1_478, mg_715, \
                         mg_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_12 * kh0_603[k]
                    - f_13 * kh1_603[k]
                    + pa_z[k] * lh_792[k];

        t_1003[k] = f_8 * lg_565[k]
                    + pb_z[k] * mg_715[k];

        t_1004[k] = f_14 * lg_582[k]
                    + f_5 * mf0_478[k]
                    - f_6 * mf1_478[k]
                    + pb_y[k] * mg_717[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, pa_y, pb_y, kh0_650, kh1_650, lg_583, lg_584, \
                         lh_818, mf0_479, mf1_479, mg_718, mg_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = f_14 * lg_583[k]
                    + f_3 * mf0_479[k]
                    - f_4 * mf1_479[k]
                    + pb_y[k] * mg_718[k];

        t_1006[k] = f_14 * lg_584[k]
                    + pb_y[k] * mg_719[k];

        t_1007[k] = f_15 * kh0_650[k]
                    - f_16 * kh1_650[k]
                    + pa_y[k] * lh_818[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pb_x, pb_y, pb_z, lg_570, lg_585, \
                         mf0_480, mf0_483, mf1_480, mf1_483, mg_720, \
                         mg_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_1 * mf0_480[k]
                    - f_2 * mf1_480[k]
                    + pb_x[k] * mg_720[k];

        t_1009[k] = f_19 * lg_585[k]
                    + pb_y[k] * mg_720[k];

        t_1010[k] = f_9 * lg_570[k]
                    + pb_z[k] * mg_720[k];

        t_1011[k] = f_5 * mf0_483[k]
                    - f_6 * mf1_483[k]
                    + pb_x[k] * mg_723[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, pb_x, pb_y, lg_587, mf0_485, mf0_486, \
                         mf1_485, mf1_486, mg_722, mg_725, mg_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_19 * lg_587[k]
                    + pb_y[k] * mg_722[k];

        t_1013[k] = f_5 * mf0_485[k]
                    - f_6 * mf1_485[k]
                    + pb_x[k] * mg_725[k];

        t_1014[k] = f_3 * mf0_486[k]
                    - f_4 * mf1_486[k]
                    + pb_x[k] * mg_726[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, t_1018, pb_x, pb_y, pb_z, lg_573, lg_590, \
                         mf0_489, mf1_489, mg_723, mg_725, mg_729, \
                         mg_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = f_9 * lg_573[k]
                    + pb_z[k] * mg_723[k];

        t_1016[k] = f_19 * lg_590[k]
                    + pb_y[k] * mg_725[k];

        t_1017[k] = f_3 * mf0_489[k]
                    - f_4 * mf1_489[k]
                    + pb_x[k] * mg_729[k];

        t_1018[k] = pb_x[k] * mg_730[k];
    }

#pragma omp simd aligned(t_1019, t_1020, t_1021, t_1022, t_1023, pa_z, pb_x, kh0_624, kh1_624, \
                         lh_813, mg_731, mg_732, mg_733, mg_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1019[k] = pb_x[k] * mg_731[k];

        t_1020[k] = pb_x[k] * mg_732[k];

        t_1021[k] = pb_x[k] * mg_733[k];

        t_1022[k] = pb_x[k] * mg_734[k];

        t_1023[k] = f_17 * kh0_624[k]
                    - f_18 * kh1_624[k]
                    + pa_z[k] * lh_813[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, pb_y, pb_z, lg_580, lg_597, lg_598, mf0_488, \
                         mf0_489, mf1_488, mf1_489, mg_730, mg_732, \
                         mg_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_9 * lg_580[k]
                    + pb_z[k] * mg_730[k];

        t_1025[k] = f_19 * lg_597[k]
                    + f_5 * mf0_488[k]
                    - f_6 * mf1_488[k]
                    + pb_y[k] * mg_732[k];

        t_1026[k] = f_19 * lg_598[k]
                    + f_3 * mf0_489[k]
                    - f_4 * mf1_489[k]
                    + pb_y[k] * mg_733[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, t_1030, pa_y, pb_x, pb_y, kh0_671, kh1_671, \
                         lg_599, lg_600, lh_839, mf0_490, mf1_490, mg_734, \
                         mg_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_19 * lg_599[k]
                    + pb_y[k] * mg_734[k];

        t_1028[k] = f_20 * kh0_671[k]
                    - f_21 * kh1_671[k]
                    + pa_y[k] * lh_839[k];

        t_1029[k] = f_1 * mf0_490[k]
                    - f_2 * mf1_490[k]
                    + pb_x[k] * mg_735[k];

        t_1030[k] = f_11 * lg_600[k]
                    + pb_y[k] * mg_735[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, pb_x, pb_y, pb_z, lg_585, lg_602, mf0_493, \
                         mf1_493, mg_735, mg_737, mg_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = f_24 * lg_585[k]
                    + pb_z[k] * mg_735[k];

        t_1032[k] = f_5 * mf0_493[k]
                    - f_6 * mf1_493[k]
                    + pb_x[k] * mg_738[k];

        t_1033[k] = f_11 * lg_602[k]
                    + pb_y[k] * mg_737[k];
    }

#pragma omp simd aligned(t_1034, t_1035, t_1036, t_1037, pb_x, pb_y, pb_z, lg_588, lg_605, \
                         mf0_495, mf0_496, mf1_495, mf1_496, mg_738, mg_740, \
                         mg_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = f_5 * mf0_495[k]
                    - f_6 * mf1_495[k]
                    + pb_x[k] * mg_740[k];

        t_1035[k] = f_3 * mf0_496[k]
                    - f_4 * mf1_496[k]
                    + pb_x[k] * mg_741[k];

        t_1036[k] = f_24 * lg_588[k]
                    + pb_z[k] * mg_738[k];

        t_1037[k] = f_11 * lg_605[k]
                    + pb_y[k] * mg_740[k];
    }

#pragma omp simd aligned(t_1038, t_1039, t_1040, t_1041, t_1042, t_1043, pb_x, mf0_499, \
                         mf1_499, mg_744, mg_745, mg_746, mg_747, mg_748, \
                         mg_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1038[k] = f_3 * mf0_499[k]
                    - f_4 * mf1_499[k]
                    + pb_x[k] * mg_744[k];

        t_1039[k] = pb_x[k] * mg_745[k];

        t_1040[k] = pb_x[k] * mg_746[k];

        t_1041[k] = pb_x[k] * mg_747[k];

        t_1042[k] = pb_x[k] * mg_748[k];

        t_1043[k] = pb_x[k] * mg_749[k];
    }

#pragma omp simd aligned(t_1044, t_1045, t_1046, pa_z, pb_y, pb_z, kh0_645, kh1_645, lg_595, \
                         lg_612, lh_834, mf0_498, mf1_498, mg_745, \
                         mg_747 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1044[k] = f_22 * kh0_645[k]
                    - f_23 * kh1_645[k]
                    + pa_z[k] * lh_834[k];

        t_1045[k] = f_24 * lg_595[k]
                    + pb_z[k] * mg_745[k];

        t_1046[k] = f_11 * lg_612[k]
                    + f_5 * mf0_498[k]
                    - f_6 * mf1_498[k]
                    + pb_y[k] * mg_747[k];
    }

#pragma omp simd aligned(t_1047, t_1048, t_1049, pa_y, pb_y, kh0_692, kh1_692, lg_613, lg_614, \
                         lh_860, mf0_499, mf1_499, mg_748, mg_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1047[k] = f_11 * lg_613[k]
                    + f_3 * mf0_499[k]
                    - f_4 * mf1_499[k]
                    + pb_y[k] * mg_748[k];

        t_1048[k] = f_11 * lg_614[k]
                    + pb_y[k] * mg_749[k];

        t_1049[k] = f_25 * kh0_692[k]
                    - f_26 * kh1_692[k]
                    + pa_y[k] * lh_860[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, t_1053, pb_x, pb_y, pb_z, lg_600, lg_615, \
                         mf0_500, mf0_503, mf1_500, mf1_503, mg_750, \
                         mg_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = f_1 * mf0_500[k]
                    - f_2 * mf1_500[k]
                    + pb_x[k] * mg_750[k];

        t_1051[k] = f_24 * lg_615[k]
                    + pb_y[k] * mg_750[k];

        t_1052[k] = f_11 * lg_600[k]
                    + pb_z[k] * mg_750[k];

        t_1053[k] = f_5 * mf0_503[k]
                    - f_6 * mf1_503[k]
                    + pb_x[k] * mg_753[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, pb_x, pb_y, lg_617, mf0_505, mf0_506, \
                         mf1_505, mf1_506, mg_752, mg_755, mg_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = f_24 * lg_617[k]
                    + pb_y[k] * mg_752[k];

        t_1055[k] = f_5 * mf0_505[k]
                    - f_6 * mf1_505[k]
                    + pb_x[k] * mg_755[k];

        t_1056[k] = f_3 * mf0_506[k]
                    - f_4 * mf1_506[k]
                    + pb_x[k] * mg_756[k];
    }

#pragma omp simd aligned(t_1057, t_1058, t_1059, t_1060, pb_x, pb_y, pb_z, lg_603, lg_620, \
                         mf0_509, mf1_509, mg_753, mg_755, mg_759, \
                         mg_760 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1057[k] = f_11 * lg_603[k]
                    + pb_z[k] * mg_753[k];

        t_1058[k] = f_24 * lg_620[k]
                    + pb_y[k] * mg_755[k];

        t_1059[k] = f_3 * mf0_509[k]
                    - f_4 * mf1_509[k]
                    + pb_x[k] * mg_759[k];

        t_1060[k] = pb_x[k] * mg_760[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, t_1064, t_1065, pa_z, pb_x, kh0_666, kh1_666, \
                         lh_855, mg_761, mg_762, mg_763, mg_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = pb_x[k] * mg_761[k];

        t_1062[k] = pb_x[k] * mg_762[k];

        t_1063[k] = pb_x[k] * mg_763[k];

        t_1064[k] = pb_x[k] * mg_764[k];

        t_1065[k] = f_25 * kh0_666[k]
                    - f_26 * kh1_666[k]
                    + pa_z[k] * lh_855[k];
    }

#pragma omp simd aligned(t_1066, t_1067, t_1068, pb_y, pb_z, lg_610, lg_627, lg_628, mf0_508, \
                         mf0_509, mf1_508, mf1_509, mg_760, mg_762, \
                         mg_763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1066[k] = f_11 * lg_610[k]
                    + pb_z[k] * mg_760[k];

        t_1067[k] = f_24 * lg_627[k]
                    + f_5 * mf0_508[k]
                    - f_6 * mf1_508[k]
                    + pb_y[k] * mg_762[k];

        t_1068[k] = f_24 * lg_628[k]
                    + f_3 * mf0_509[k]
                    - f_4 * mf1_509[k]
                    + pb_y[k] * mg_763[k];
    }

#pragma omp simd aligned(t_1069, t_1070, t_1071, t_1072, pa_y, pb_x, pb_y, kh0_713, kh1_713, \
                         lg_629, lg_630, lh_881, mf0_510, mf1_510, mg_764, \
                         mg_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1069[k] = f_24 * lg_629[k]
                    + pb_y[k] * mg_764[k];

        t_1070[k] = f_22 * kh0_713[k]
                    - f_23 * kh1_713[k]
                    + pa_y[k] * lh_881[k];

        t_1071[k] = f_1 * mf0_510[k]
                    - f_2 * mf1_510[k]
                    + pb_x[k] * mg_765[k];

        t_1072[k] = f_9 * lg_630[k]
                    + pb_y[k] * mg_765[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pb_x, pb_y, pb_z, lg_615, lg_632, mf0_513, \
                         mf1_513, mg_765, mg_767, mg_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_19 * lg_615[k]
                    + pb_z[k] * mg_765[k];

        t_1074[k] = f_5 * mf0_513[k]
                    - f_6 * mf1_513[k]
                    + pb_x[k] * mg_768[k];

        t_1075[k] = f_9 * lg_632[k]
                    + pb_y[k] * mg_767[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, t_1079, pb_x, pb_y, pb_z, lg_618, lg_635, \
                         mf0_515, mf0_516, mf1_515, mf1_516, mg_768, mg_770, \
                         mg_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_5 * mf0_515[k]
                    - f_6 * mf1_515[k]
                    + pb_x[k] * mg_770[k];

        t_1077[k] = f_3 * mf0_516[k]
                    - f_4 * mf1_516[k]
                    + pb_x[k] * mg_771[k];

        t_1078[k] = f_19 * lg_618[k]
                    + pb_z[k] * mg_768[k];

        t_1079[k] = f_9 * lg_635[k]
                    + pb_y[k] * mg_770[k];
    }
}

static auto
compute_prim_mh_electron_repulsion_0_piece8(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kh0,
                                            const size_t kh1, const size_t lg, const size_t lh,
                                            const size_t mf0, const size_t mf1, const size_t mg,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 4.0 / p;
    const auto f_11 = 2.5 / p;
    const auto f_12 = 0.5 / alpha;
    const auto f_13 = 0.5 * beta / (alpha * p);
    const auto f_14 = 3.5 / p;
    const auto f_15 = 3.0 / alpha;
    const auto f_16 = 3.0 * beta / (alpha * p);
    const auto f_17 = 1.0 / alpha;
    const auto f_18 = beta / (alpha * p);
    const auto f_19 = 3.0 / p;
    const auto f_20 = 2.5 / alpha;
    const auto f_21 = 2.5 * beta / (alpha * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kh0_687 = buffer.data(kh0 + 687);
    const auto *kh0_708 = buffer.data(kh0 + 708);
    const auto *kh0_734 = buffer.data(kh0 + 734);
    const auto *kh0_755 = buffer.data(kh0 + 755);

    const auto *kh1_687 = buffer.data(kh1 + 687);
    const auto *kh1_708 = buffer.data(kh1 + 708);
    const auto *kh1_734 = buffer.data(kh1 + 734);
    const auto *kh1_755 = buffer.data(kh1 + 755);

    const auto *lg_625 = buffer.data(lg + 625);
    const auto *lg_630 = buffer.data(lg + 630);
    const auto *lg_633 = buffer.data(lg + 633);
    const auto *lg_640 = buffer.data(lg + 640);
    const auto *lg_642 = buffer.data(lg + 642);
    const auto *lg_643 = buffer.data(lg + 643);
    const auto *lg_644 = buffer.data(lg + 644);
    const auto *lg_645 = buffer.data(lg + 645);
    const auto *lg_647 = buffer.data(lg + 647);
    const auto *lg_648 = buffer.data(lg + 648);
    const auto *lg_650 = buffer.data(lg + 650);
    const auto *lg_655 = buffer.data(lg + 655);
    const auto *lg_657 = buffer.data(lg + 657);
    const auto *lg_658 = buffer.data(lg + 658);
    const auto *lg_659 = buffer.data(lg + 659);
    const auto *lg_660 = buffer.data(lg + 660);
    const auto *lg_661 = buffer.data(lg + 661);
    const auto *lg_662 = buffer.data(lg + 662);
    const auto *lg_663 = buffer.data(lg + 663);
    const auto *lg_665 = buffer.data(lg + 665);
    const auto *lg_670 = buffer.data(lg + 670);
    const auto *lg_672 = buffer.data(lg + 672);
    const auto *lg_673 = buffer.data(lg + 673);
    const auto *lg_674 = buffer.data(lg + 674);

    const auto *lh_876 = buffer.data(lh + 876);
    const auto *lh_897 = buffer.data(lh + 897);
    const auto *lh_902 = buffer.data(lh + 902);
    const auto *lh_923 = buffer.data(lh + 923);
    const auto *lh_924 = buffer.data(lh + 924);
    const auto *lh_926 = buffer.data(lh + 926);
    const auto *lh_927 = buffer.data(lh + 927);
    const auto *lh_929 = buffer.data(lh + 929);
    const auto *lh_930 = buffer.data(lh + 930);
    const auto *lh_933 = buffer.data(lh + 933);
    const auto *lh_939 = buffer.data(lh + 939);
    const auto *lh_941 = buffer.data(lh + 941);
    const auto *lh_942 = buffer.data(lh + 942);
    const auto *lh_944 = buffer.data(lh + 944);

    const auto *mf0_518 = buffer.data(mf0 + 518);
    const auto *mf0_519 = buffer.data(mf0 + 519);
    const auto *mf0_520 = buffer.data(mf0 + 520);
    const auto *mf0_523 = buffer.data(mf0 + 523);
    const auto *mf0_525 = buffer.data(mf0 + 525);
    const auto *mf0_526 = buffer.data(mf0 + 526);
    const auto *mf0_528 = buffer.data(mf0 + 528);
    const auto *mf0_529 = buffer.data(mf0 + 529);
    const auto *mf0_540 = buffer.data(mf0 + 540);
    const auto *mf0_543 = buffer.data(mf0 + 543);
    const auto *mf0_545 = buffer.data(mf0 + 545);
    const auto *mf0_546 = buffer.data(mf0 + 546);
    const auto *mf0_548 = buffer.data(mf0 + 548);
    const auto *mf0_549 = buffer.data(mf0 + 549);

    const auto *mf1_518 = buffer.data(mf1 + 518);
    const auto *mf1_519 = buffer.data(mf1 + 519);
    const auto *mf1_520 = buffer.data(mf1 + 520);
    const auto *mf1_523 = buffer.data(mf1 + 523);
    const auto *mf1_525 = buffer.data(mf1 + 525);
    const auto *mf1_526 = buffer.data(mf1 + 526);
    const auto *mf1_528 = buffer.data(mf1 + 528);
    const auto *mf1_529 = buffer.data(mf1 + 529);
    const auto *mf1_540 = buffer.data(mf1 + 540);
    const auto *mf1_543 = buffer.data(mf1 + 543);
    const auto *mf1_545 = buffer.data(mf1 + 545);
    const auto *mf1_546 = buffer.data(mf1 + 546);
    const auto *mf1_548 = buffer.data(mf1 + 548);
    const auto *mf1_549 = buffer.data(mf1 + 549);

    const auto *mg_774 = buffer.data(mg + 774);
    const auto *mg_775 = buffer.data(mg + 775);
    const auto *mg_776 = buffer.data(mg + 776);
    const auto *mg_777 = buffer.data(mg + 777);
    const auto *mg_778 = buffer.data(mg + 778);
    const auto *mg_779 = buffer.data(mg + 779);
    const auto *mg_780 = buffer.data(mg + 780);
    const auto *mg_782 = buffer.data(mg + 782);
    const auto *mg_783 = buffer.data(mg + 783);
    const auto *mg_785 = buffer.data(mg + 785);
    const auto *mg_786 = buffer.data(mg + 786);
    const auto *mg_789 = buffer.data(mg + 789);
    const auto *mg_790 = buffer.data(mg + 790);
    const auto *mg_791 = buffer.data(mg + 791);
    const auto *mg_792 = buffer.data(mg + 792);
    const auto *mg_793 = buffer.data(mg + 793);
    const auto *mg_794 = buffer.data(mg + 794);
    const auto *mg_795 = buffer.data(mg + 795);
    const auto *mg_797 = buffer.data(mg + 797);
    const auto *mg_798 = buffer.data(mg + 798);
    const auto *mg_800 = buffer.data(mg + 800);
    const auto *mg_805 = buffer.data(mg + 805);
    const auto *mg_806 = buffer.data(mg + 806);
    const auto *mg_807 = buffer.data(mg + 807);
    const auto *mg_808 = buffer.data(mg + 808);
    const auto *mg_809 = buffer.data(mg + 809);
    const auto *mg_810 = buffer.data(mg + 810);
    const auto *mg_812 = buffer.data(mg + 812);
    const auto *mg_813 = buffer.data(mg + 813);
    const auto *mg_815 = buffer.data(mg + 815);
    const auto *mg_816 = buffer.data(mg + 816);
    const auto *mg_819 = buffer.data(mg + 819);
    const auto *mg_820 = buffer.data(mg + 820);
    const auto *mg_821 = buffer.data(mg + 821);
    const auto *mg_822 = buffer.data(mg + 822);
    const auto *mg_823 = buffer.data(mg + 823);
    const auto *mg_824 = buffer.data(mg + 824);

#pragma omp simd aligned(t_1080, t_1081, t_1082, t_1083, t_1084, t_1085, pb_x, mf0_519, \
                         mf1_519, mg_774, mg_775, mg_776, mg_777, mg_778, \
                         mg_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = f_3 * mf0_519[k]
                    - f_4 * mf1_519[k]
                    + pb_x[k] * mg_774[k];

        t_1081[k] = pb_x[k] * mg_775[k];

        t_1082[k] = pb_x[k] * mg_776[k];

        t_1083[k] = pb_x[k] * mg_777[k];

        t_1084[k] = pb_x[k] * mg_778[k];

        t_1085[k] = pb_x[k] * mg_779[k];
    }

#pragma omp simd aligned(t_1086, t_1087, t_1088, pa_z, pb_y, pb_z, kh0_687, kh1_687, lg_625, \
                         lg_642, lh_876, mf0_518, mf1_518, mg_775, \
                         mg_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1086[k] = f_20 * kh0_687[k]
                    - f_21 * kh1_687[k]
                    + pa_z[k] * lh_876[k];

        t_1087[k] = f_19 * lg_625[k]
                    + pb_z[k] * mg_775[k];

        t_1088[k] = f_9 * lg_642[k]
                    + f_5 * mf0_518[k]
                    - f_6 * mf1_518[k]
                    + pb_y[k] * mg_777[k];
    }

#pragma omp simd aligned(t_1089, t_1090, t_1091, pa_y, pb_y, kh0_734, kh1_734, lg_643, lg_644, \
                         lh_902, mf0_519, mf1_519, mg_778, mg_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1089[k] = f_9 * lg_643[k]
                    + f_3 * mf0_519[k]
                    - f_4 * mf1_519[k]
                    + pb_y[k] * mg_778[k];

        t_1090[k] = f_9 * lg_644[k]
                    + pb_y[k] * mg_779[k];

        t_1091[k] = f_17 * kh0_734[k]
                    - f_18 * kh1_734[k]
                    + pa_y[k] * lh_902[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, t_1095, pb_x, pb_y, pb_z, lg_630, lg_645, \
                         mf0_520, mf0_523, mf1_520, mf1_523, mg_780, \
                         mg_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = f_1 * mf0_520[k]
                    - f_2 * mf1_520[k]
                    + pb_x[k] * mg_780[k];

        t_1093[k] = f_8 * lg_645[k]
                    + pb_y[k] * mg_780[k];

        t_1094[k] = f_14 * lg_630[k]
                    + pb_z[k] * mg_780[k];

        t_1095[k] = f_5 * mf0_523[k]
                    - f_6 * mf1_523[k]
                    + pb_x[k] * mg_783[k];
    }

#pragma omp simd aligned(t_1096, t_1097, t_1098, pb_x, pb_y, lg_647, mf0_525, mf0_526, \
                         mf1_525, mf1_526, mg_782, mg_785, mg_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1096[k] = f_8 * lg_647[k]
                    + pb_y[k] * mg_782[k];

        t_1097[k] = f_5 * mf0_525[k]
                    - f_6 * mf1_525[k]
                    + pb_x[k] * mg_785[k];

        t_1098[k] = f_3 * mf0_526[k]
                    - f_4 * mf1_526[k]
                    + pb_x[k] * mg_786[k];
    }

#pragma omp simd aligned(t_1099, t_1100, t_1101, t_1102, pb_x, pb_y, pb_z, lg_633, lg_650, \
                         mf0_529, mf1_529, mg_783, mg_785, mg_789, \
                         mg_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1099[k] = f_14 * lg_633[k]
                    + pb_z[k] * mg_783[k];

        t_1100[k] = f_8 * lg_650[k]
                    + pb_y[k] * mg_785[k];

        t_1101[k] = f_3 * mf0_529[k]
                    - f_4 * mf1_529[k]
                    + pb_x[k] * mg_789[k];

        t_1102[k] = pb_x[k] * mg_790[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, t_1106, t_1107, pa_z, pb_x, kh0_708, kh1_708, \
                         lh_897, mg_791, mg_792, mg_793, mg_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = pb_x[k] * mg_791[k];

        t_1104[k] = pb_x[k] * mg_792[k];

        t_1105[k] = pb_x[k] * mg_793[k];

        t_1106[k] = pb_x[k] * mg_794[k];

        t_1107[k] = f_15 * kh0_708[k]
                    - f_16 * kh1_708[k]
                    + pa_z[k] * lh_897[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, pb_y, pb_z, lg_640, lg_657, lg_658, mf0_528, \
                         mf0_529, mf1_528, mf1_529, mg_790, mg_792, \
                         mg_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_14 * lg_640[k]
                    + pb_z[k] * mg_790[k];

        t_1109[k] = f_8 * lg_657[k]
                    + f_5 * mf0_528[k]
                    - f_6 * mf1_528[k]
                    + pb_y[k] * mg_792[k];

        t_1110[k] = f_8 * lg_658[k]
                    + f_3 * mf0_529[k]
                    - f_4 * mf1_529[k]
                    + pb_y[k] * mg_793[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, t_1114, t_1115, pa_y, pb_y, kh0_755, kh1_755, \
                         lg_659, lg_660, lh_923, lh_924, lh_926, mg_794, \
                         mg_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_8 * lg_659[k]
                    + pb_y[k] * mg_794[k];

        t_1112[k] = f_12 * kh0_755[k]
                    - f_13 * kh1_755[k]
                    + pa_y[k] * lh_923[k];

        t_1113[k] = pa_y[k] * lh_924[k];

        t_1114[k] = f_7 * lg_660[k]
                    + pb_y[k] * mg_795[k];

        t_1115[k] = pa_y[k] * lh_926[k];
    }

#pragma omp simd aligned(t_1116, t_1117, t_1118, t_1119, pa_y, pb_y, lg_661, lg_662, lg_663, \
                         lh_927, lh_929, lh_930, mg_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1116[k] = f_8 * lg_661[k]
                    + pa_y[k] * lh_927[k];

        t_1117[k] = f_7 * lg_662[k]
                    + pb_y[k] * mg_797[k];

        t_1118[k] = pa_y[k] * lh_929[k];

        t_1119[k] = f_9 * lg_663[k]
                    + pa_y[k] * lh_930[k];
    }

#pragma omp simd aligned(t_1120, t_1121, t_1122, t_1123, t_1124, pa_y, pb_x, pb_y, pb_z, \
                         lg_648, lg_665, lh_933, mg_798, mg_800, mg_805, \
                         mg_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1120[k] = f_10 * lg_648[k]
                    + pb_z[k] * mg_798[k];

        t_1121[k] = f_7 * lg_665[k]
                    + pb_y[k] * mg_800[k];

        t_1122[k] = pa_y[k] * lh_933[k];

        t_1123[k] = pb_x[k] * mg_805[k];

        t_1124[k] = pb_x[k] * mg_806[k];
    }

#pragma omp simd aligned(t_1125, t_1126, t_1127, t_1128, t_1129, pa_y, pb_x, pb_z, lg_655, \
                         lg_670, lh_939, mg_805, mg_807, mg_808, \
                         mg_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1125[k] = pb_x[k] * mg_807[k];

        t_1126[k] = pb_x[k] * mg_808[k];

        t_1127[k] = pb_x[k] * mg_809[k];

        t_1128[k] = f_11 * lg_670[k]
                    + pa_y[k] * lh_939[k];

        t_1129[k] = f_10 * lg_655[k]
                    + pb_z[k] * mg_805[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, t_1133, pa_y, pb_y, lg_672, lg_673, lg_674, \
                         lh_941, lh_942, lh_944, mg_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = f_9 * lg_672[k]
                    + pa_y[k] * lh_941[k];

        t_1131[k] = f_8 * lg_673[k]
                    + pa_y[k] * lh_942[k];

        t_1132[k] = f_7 * lg_674[k]
                    + pb_y[k] * mg_809[k];

        t_1133[k] = pa_y[k] * lh_944[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, t_1138, pb_x, pb_y, pb_z, lg_660, \
                         mf0_540, mf0_543, mf1_540, mf1_543, mg_810, mg_812, \
                         mg_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = f_1 * mf0_540[k]
                    - f_2 * mf1_540[k]
                    + pb_x[k] * mg_810[k];

        t_1135[k] = pb_y[k] * mg_810[k];

        t_1136[k] = f_0 * lg_660[k]
                    + pb_z[k] * mg_810[k];

        t_1137[k] = f_5 * mf0_543[k]
                    - f_6 * mf1_543[k]
                    + pb_x[k] * mg_813[k];

        t_1138[k] = pb_y[k] * mg_812[k];
    }

#pragma omp simd aligned(t_1139, t_1140, t_1141, t_1142, pb_x, pb_y, pb_z, lg_663, mf0_545, \
                         mf0_546, mf1_545, mf1_546, mg_813, mg_815, \
                         mg_816 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = f_5 * mf0_545[k]
                    - f_6 * mf1_545[k]
                    + pb_x[k] * mg_815[k];

        t_1140[k] = f_3 * mf0_546[k]
                    - f_4 * mf1_546[k]
                    + pb_x[k] * mg_816[k];

        t_1141[k] = f_0 * lg_663[k]
                    + pb_z[k] * mg_813[k];

        t_1142[k] = pb_y[k] * mg_815[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, t_1146, t_1147, t_1148, pb_x, mf0_549, \
                         mf1_549, mg_819, mg_820, mg_821, mg_822, mg_823, \
                         mg_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = f_3 * mf0_549[k]
                    - f_4 * mf1_549[k]
                    + pb_x[k] * mg_819[k];

        t_1144[k] = pb_x[k] * mg_820[k];

        t_1145[k] = pb_x[k] * mg_821[k];

        t_1146[k] = pb_x[k] * mg_822[k];

        t_1147[k] = pb_x[k] * mg_823[k];

        t_1148[k] = pb_x[k] * mg_824[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, t_1152, pb_y, pb_z, lg_670, mf0_546, mf0_548, \
                         mf0_549, mf1_546, mf1_548, mf1_549, mg_820, mg_822, \
                         mg_823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_1 * mf0_546[k]
                    - f_2 * mf1_546[k]
                    + pb_y[k] * mg_820[k];

        t_1150[k] = f_0 * lg_670[k]
                    + pb_z[k] * mg_820[k];

        t_1151[k] = f_5 * mf0_548[k]
                    - f_6 * mf1_548[k]
                    + pb_y[k] * mg_822[k];

        t_1152[k] = f_3 * mf0_549[k]
                    - f_4 * mf1_549[k]
                    + pb_y[k] * mg_823[k];
    }

#pragma omp simd aligned(t_1153, t_1154, pb_y, pb_z, lg_674, mf0_549, mf1_549, \
                         mg_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = pb_y[k] * mg_824[k];

        t_1154[k] = f_0 * lg_674[k]
                    + f_1 * mf0_549[k]
                    - f_2 * mf1_549[k]
                    + pb_z[k] * mg_824[k];
    }
}

auto
compute_prim_mh_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t kh0, const size_t kh1,
                                     const size_t lg, const size_t lh, const size_t mf0,
                                     const size_t mf1, const size_t mg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    compute_prim_mh_electron_repulsion_0_piece0(buffer, target, pa, pb, kh0, kh1, lg, lh, mf0,
                                                mf1, mg, ncols, alpha, beta, p);

    compute_prim_mh_electron_repulsion_0_piece1(buffer, target, pa, pb, kh0, kh1, lg, lh, mf0,
                                                mf1, mg, ncols, alpha, beta, p);

    compute_prim_mh_electron_repulsion_0_piece2(buffer, target, pa, pb, kh0, kh1, lg, lh, mf0,
                                                mf1, mg, ncols, alpha, beta, p);

    compute_prim_mh_electron_repulsion_0_piece3(buffer, target, pa, pb, kh0, kh1, lg, lh, mf0,
                                                mf1, mg, ncols, alpha, beta, p);

    compute_prim_mh_electron_repulsion_0_piece4(buffer, target, pa, pb, kh0, kh1, lg, lh, mf0,
                                                mf1, mg, ncols, alpha, beta, p);

    compute_prim_mh_electron_repulsion_0_piece5(buffer, target, pa, pb, kh0, kh1, lg, lh, mf0,
                                                mf1, mg, ncols, alpha, beta, p);

    compute_prim_mh_electron_repulsion_0_piece6(buffer, target, pa, pb, lg, lh, mf0, mf1, mg,
                                                ncols, alpha, beta, p);

    compute_prim_mh_electron_repulsion_0_piece7(buffer, target, pa, pb, kh0, kh1, lg, lh, mf0,
                                                mf1, mg, ncols, alpha, beta, p);

    compute_prim_mh_electron_repulsion_0_piece8(buffer, target, pa, pb, kh0, kh1, lg, lh, mf0,
                                                mf1, mg, ncols, alpha, beta, p);
}

}  // namespace simdt2ceri
