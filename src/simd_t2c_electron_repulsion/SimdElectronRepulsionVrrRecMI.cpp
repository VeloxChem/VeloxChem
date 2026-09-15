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


#include "SimdElectronRepulsionVrrRecMI.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_mi_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t ki0,
                                            const size_t ki1, const size_t lh, const size_t li,
                                            const size_t mg0, const size_t mg1, const size_t mh,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_13 = 4.0 / p;
    const auto f_14 = 3.0 / p;
    const auto f_15 = 0.5 / alpha;
    const auto f_16 = 0.5 * beta / (alpha * p);
    const auto f_17 = 3.5 / p;
    const auto f_18 = 3.0 / alpha;
    const auto f_19 = 3.0 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ki0_0 = buffer.data(ki0 + 0);
    const auto *ki0_105 = buffer.data(ki0 + 105);

    const auto *ki1_0 = buffer.data(ki1 + 0);
    const auto *ki1_105 = buffer.data(ki1 + 105);

    const auto *lh_0 = buffer.data(lh + 0);
    const auto *lh_1 = buffer.data(lh + 1);
    const auto *lh_2 = buffer.data(lh + 2);
    const auto *lh_3 = buffer.data(lh + 3);
    const auto *lh_5 = buffer.data(lh + 5);
    const auto *lh_6 = buffer.data(lh + 6);
    const auto *lh_7 = buffer.data(lh + 7);
    const auto *lh_8 = buffer.data(lh + 8);
    const auto *lh_9 = buffer.data(lh + 9);
    const auto *lh_15 = buffer.data(lh + 15);
    const auto *lh_16 = buffer.data(lh + 16);
    const auto *lh_17 = buffer.data(lh + 17);
    const auto *lh_18 = buffer.data(lh + 18);
    const auto *lh_19 = buffer.data(lh + 19);
    const auto *lh_20 = buffer.data(lh + 20);
    const auto *lh_21 = buffer.data(lh + 21);
    const auto *lh_24 = buffer.data(lh + 24);
    const auto *lh_26 = buffer.data(lh + 26);
    const auto *lh_27 = buffer.data(lh + 27);
    const auto *lh_30 = buffer.data(lh + 30);
    const auto *lh_36 = buffer.data(lh + 36);
    const auto *lh_38 = buffer.data(lh + 38);
    const auto *lh_39 = buffer.data(lh + 39);
    const auto *lh_40 = buffer.data(lh + 40);
    const auto *lh_41 = buffer.data(lh + 41);
    const auto *lh_42 = buffer.data(lh + 42);
    const auto *lh_44 = buffer.data(lh + 44);
    const auto *lh_47 = buffer.data(lh + 47);
    const auto *lh_50 = buffer.data(lh + 50);
    const auto *lh_51 = buffer.data(lh + 51);
    const auto *lh_58 = buffer.data(lh + 58);
    const auto *lh_59 = buffer.data(lh + 59);
    const auto *lh_60 = buffer.data(lh + 60);
    const auto *lh_61 = buffer.data(lh + 61);
    const auto *lh_62 = buffer.data(lh + 62);
    const auto *lh_66 = buffer.data(lh + 66);
    const auto *lh_69 = buffer.data(lh + 69);
    const auto *lh_73 = buffer.data(lh + 73);
    const auto *lh_78 = buffer.data(lh + 78);
    const auto *lh_80 = buffer.data(lh + 80);
    const auto *lh_81 = buffer.data(lh + 81);
    const auto *lh_82 = buffer.data(lh + 82);
    const auto *lh_83 = buffer.data(lh + 83);
    const auto *lh_100 = buffer.data(lh + 100);
    const auto *lh_101 = buffer.data(lh + 101);
    const auto *lh_102 = buffer.data(lh + 102);
    const auto *lh_103 = buffer.data(lh + 103);

    const auto *li_0 = buffer.data(li + 0);
    const auto *li_3 = buffer.data(li + 3);
    const auto *li_5 = buffer.data(li + 5);
    const auto *li_6 = buffer.data(li + 6);
    const auto *li_9 = buffer.data(li + 9);
    const auto *li_10 = buffer.data(li + 10);
    const auto *li_12 = buffer.data(li + 12);
    const auto *li_14 = buffer.data(li + 14);
    const auto *li_15 = buffer.data(li + 15);
    const auto *li_20 = buffer.data(li + 20);
    const auto *li_21 = buffer.data(li + 21);
    const auto *li_23 = buffer.data(li + 23);
    const auto *li_24 = buffer.data(li + 24);
    const auto *li_25 = buffer.data(li + 25);
    const auto *li_27 = buffer.data(li + 27);
    const auto *li_28 = buffer.data(li + 28);
    const auto *li_29 = buffer.data(li + 29);
    const auto *li_31 = buffer.data(li + 31);
    const auto *li_34 = buffer.data(li + 34);
    const auto *li_38 = buffer.data(li + 38);
    const auto *li_43 = buffer.data(li + 43);
    const auto *li_49 = buffer.data(li + 49);
    const auto *li_56 = buffer.data(li + 56);
    const auto *li_58 = buffer.data(li + 58);
    const auto *li_61 = buffer.data(li + 61);
    const auto *li_65 = buffer.data(li + 65);
    const auto *li_68 = buffer.data(li + 68);
    const auto *li_70 = buffer.data(li + 70);
    const auto *li_76 = buffer.data(li + 76);
    const auto *li_79 = buffer.data(li + 79);
    const auto *li_80 = buffer.data(li + 80);
    const auto *li_81 = buffer.data(li + 81);
    const auto *li_83 = buffer.data(li + 83);
    const auto *li_105 = buffer.data(li + 105);

    const auto *mg0_0 = buffer.data(mg0 + 0);
    const auto *mg0_1 = buffer.data(mg0 + 1);
    const auto *mg0_2 = buffer.data(mg0 + 2);
    const auto *mg0_3 = buffer.data(mg0 + 3);
    const auto *mg0_5 = buffer.data(mg0 + 5);
    const auto *mg0_10 = buffer.data(mg0 + 10);
    const auto *mg0_12 = buffer.data(mg0 + 12);
    const auto *mg0_13 = buffer.data(mg0 + 13);
    const auto *mg0_14 = buffer.data(mg0 + 14);
    const auto *mg0_45 = buffer.data(mg0 + 45);
    const auto *mg0_47 = buffer.data(mg0 + 47);
    const auto *mg0_48 = buffer.data(mg0 + 48);
    const auto *mg0_50 = buffer.data(mg0 + 50);
    const auto *mg0_51 = buffer.data(mg0 + 51);
    const auto *mg0_55 = buffer.data(mg0 + 55);
    const auto *mg0_56 = buffer.data(mg0 + 56);
    const auto *mg0_57 = buffer.data(mg0 + 57);
    const auto *mg0_59 = buffer.data(mg0 + 59);
    const auto *mg0_75 = buffer.data(mg0 + 75);

    const auto *mg1_0 = buffer.data(mg1 + 0);
    const auto *mg1_1 = buffer.data(mg1 + 1);
    const auto *mg1_2 = buffer.data(mg1 + 2);
    const auto *mg1_3 = buffer.data(mg1 + 3);
    const auto *mg1_5 = buffer.data(mg1 + 5);
    const auto *mg1_10 = buffer.data(mg1 + 10);
    const auto *mg1_12 = buffer.data(mg1 + 12);
    const auto *mg1_13 = buffer.data(mg1 + 13);
    const auto *mg1_14 = buffer.data(mg1 + 14);
    const auto *mg1_45 = buffer.data(mg1 + 45);
    const auto *mg1_47 = buffer.data(mg1 + 47);
    const auto *mg1_48 = buffer.data(mg1 + 48);
    const auto *mg1_50 = buffer.data(mg1 + 50);
    const auto *mg1_51 = buffer.data(mg1 + 51);
    const auto *mg1_55 = buffer.data(mg1 + 55);
    const auto *mg1_56 = buffer.data(mg1 + 56);
    const auto *mg1_57 = buffer.data(mg1 + 57);
    const auto *mg1_59 = buffer.data(mg1 + 59);
    const auto *mg1_75 = buffer.data(mg1 + 75);

    const auto *mh_0 = buffer.data(mh + 0);
    const auto *mh_1 = buffer.data(mh + 1);
    const auto *mh_2 = buffer.data(mh + 2);
    const auto *mh_3 = buffer.data(mh + 3);
    const auto *mh_5 = buffer.data(mh + 5);
    const auto *mh_6 = buffer.data(mh + 6);
    const auto *mh_8 = buffer.data(mh + 8);
    const auto *mh_9 = buffer.data(mh + 9);
    const auto *mh_10 = buffer.data(mh + 10);
    const auto *mh_14 = buffer.data(mh + 14);
    const auto *mh_15 = buffer.data(mh + 15);
    const auto *mh_17 = buffer.data(mh + 17);
    const auto *mh_18 = buffer.data(mh + 18);
    const auto *mh_19 = buffer.data(mh + 19);
    const auto *mh_20 = buffer.data(mh + 20);
    const auto *mh_21 = buffer.data(mh + 21);
    const auto *mh_22 = buffer.data(mh + 22);
    const auto *mh_24 = buffer.data(mh + 24);
    const auto *mh_26 = buffer.data(mh + 26);
    const auto *mh_27 = buffer.data(mh + 27);
    const auto *mh_30 = buffer.data(mh + 30);
    const auto *mh_31 = buffer.data(mh + 31);
    const auto *mh_36 = buffer.data(mh + 36);
    const auto *mh_38 = buffer.data(mh + 38);
    const auto *mh_39 = buffer.data(mh + 39);
    const auto *mh_40 = buffer.data(mh + 40);
    const auto *mh_41 = buffer.data(mh + 41);
    const auto *mh_42 = buffer.data(mh + 42);
    const auto *mh_44 = buffer.data(mh + 44);
    const auto *mh_45 = buffer.data(mh + 45);
    const auto *mh_47 = buffer.data(mh + 47);
    const auto *mh_48 = buffer.data(mh + 48);
    const auto *mh_51 = buffer.data(mh + 51);
    const auto *mh_56 = buffer.data(mh + 56);
    const auto *mh_57 = buffer.data(mh + 57);
    const auto *mh_58 = buffer.data(mh + 58);
    const auto *mh_59 = buffer.data(mh + 59);
    const auto *mh_60 = buffer.data(mh + 60);
    const auto *mh_62 = buffer.data(mh + 62);
    const auto *mh_63 = buffer.data(mh + 63);
    const auto *mh_64 = buffer.data(mh + 64);
    const auto *mh_65 = buffer.data(mh + 65);
    const auto *mh_66 = buffer.data(mh + 66);
    const auto *mh_68 = buffer.data(mh + 68);
    const auto *mh_69 = buffer.data(mh + 69);
    const auto *mh_70 = buffer.data(mh + 70);
    const auto *mh_72 = buffer.data(mh + 72);
    const auto *mh_73 = buffer.data(mh + 73);
    const auto *mh_78 = buffer.data(mh + 78);
    const auto *mh_79 = buffer.data(mh + 79);
    const auto *mh_80 = buffer.data(mh + 80);
    const auto *mh_81 = buffer.data(mh + 81);
    const auto *mh_82 = buffer.data(mh + 82);
    const auto *mh_83 = buffer.data(mh + 83);
    const auto *mh_86 = buffer.data(mh + 86);
    const auto *mh_87 = buffer.data(mh + 87);
    const auto *mh_89 = buffer.data(mh + 89);
    const auto *mh_90 = buffer.data(mh + 90);
    const auto *mh_93 = buffer.data(mh + 93);
    const auto *mh_99 = buffer.data(mh + 99);
    const auto *mh_100 = buffer.data(mh + 100);
    const auto *mh_101 = buffer.data(mh + 101);
    const auto *mh_102 = buffer.data(mh + 102);
    const auto *mh_103 = buffer.data(mh + 103);
    const auto *mh_104 = buffer.data(mh + 104);
    const auto *mh_105 = buffer.data(mh + 105);
    const auto *mh_106 = buffer.data(mh + 106);
    const auto *mh_107 = buffer.data(mh + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, lh_0, mg0_0, mg1_0, \
                         mh_0, mh_1, mh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lh_0[k]
                 + f_1 * mg0_0[k]
                 - f_2 * mg1_0[k]
                 + pb_x[k] * mh_0[k];

        t_1[k] = pb_y[k] * mh_0[k];

        t_2[k] = pb_z[k] * mh_0[k];

        t_3[k] = f_3 * mg0_0[k]
                 - f_4 * mg1_0[k]
                 + pb_y[k] * mh_1[k];

        t_4[k] = pb_y[k] * mh_2[k];

        t_5[k] = f_3 * mg0_0[k]
                 - f_4 * mg1_0[k]
                 + pb_z[k] * mh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, mg0_1, mg0_2, mg0_3, mg1_1, \
                         mg1_2, mg1_3, mh_3, mh_5, mh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * mg0_1[k]
                 - f_6 * mg1_1[k]
                 + pb_y[k] * mh_3[k];

        t_7[k] = pb_z[k] * mh_3[k];

        t_8[k] = pb_y[k] * mh_5[k];

        t_9[k] = f_5 * mg0_2[k]
                 - f_6 * mg1_2[k]
                 + pb_z[k] * mh_5[k];

        t_10[k] = f_7 * mg0_3[k]
                  - f_8 * mg1_3[k]
                  + pb_y[k] * mh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pb_x, pb_y, pb_z, lh_15, mg0_5, mg1_5, \
                         mh_6, mh_8, mh_9, mh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * mh_6[k];

        t_12[k] = f_3 * mg0_5[k]
                  - f_4 * mg1_5[k]
                  + pb_y[k] * mh_8[k];

        t_13[k] = pb_y[k] * mh_9[k];

        t_14[k] = f_7 * mg0_5[k]
                  - f_8 * mg1_5[k]
                  + pb_z[k] * mh_9[k];

        t_15[k] = f_0 * lh_15[k]
                  + pb_x[k] * mh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pb_x, pb_y, pb_z, lh_17, lh_18, lh_20, \
                         mh_10, mh_14, mh_17, mh_18, mh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_z[k] * mh_10[k];

        t_17[k] = f_0 * lh_17[k]
                  + pb_x[k] * mh_17[k];

        t_18[k] = f_0 * lh_18[k]
                  + pb_x[k] * mh_18[k];

        t_19[k] = pb_y[k] * mh_14[k];

        t_20[k] = f_0 * lh_20[k]
                  + pb_x[k] * mh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pb_z, mg0_10, mg0_12, mg0_13, mg1_10, \
                         mg1_12, mg1_13, mh_15, mh_17, mh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * mg0_10[k]
                  - f_2 * mg1_10[k]
                  + pb_y[k] * mh_15[k];

        t_22[k] = pb_z[k] * mh_15[k];

        t_23[k] = f_7 * mg0_12[k]
                  - f_8 * mg1_12[k]
                  + pb_y[k] * mh_17[k];

        t_24[k] = f_5 * mg0_13[k]
                  - f_6 * mg1_13[k]
                  + pb_y[k] * mh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, pa_y, pb_y, pb_z, lh_0, li_0, \
                         mg0_14, mg1_14, mh_19, mh_20, mh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * mg0_14[k]
                  - f_4 * mg1_14[k]
                  + pb_y[k] * mh_19[k];

        t_26[k] = pb_y[k] * mh_20[k];

        t_27[k] = f_1 * mg0_14[k]
                  - f_2 * mg1_14[k]
                  + pb_z[k] * mh_20[k];

        t_28[k] = pa_y[k] * li_0[k];

        t_29[k] = f_9 * lh_0[k]
                  + pb_y[k] * mh_21[k];

        t_30[k] = pb_z[k] * mh_21[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_y, pb_z, lh_1, lh_3, li_3, li_5, \
                         li_6, mh_22, mh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * lh_1[k]
                  + pa_y[k] * li_3[k];

        t_32[k] = pb_z[k] * mh_22[k];

        t_33[k] = pa_y[k] * li_5[k];

        t_34[k] = f_11 * lh_3[k]
                  + pa_y[k] * li_6[k];

        t_35[k] = pb_z[k] * mh_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_y, pb_y, pb_z, lh_5, lh_6, lh_8, \
                         li_9, li_10, li_12, mh_26, mh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * lh_5[k]
                  + pb_y[k] * mh_26[k];

        t_37[k] = pa_y[k] * li_9[k];

        t_38[k] = f_12 * lh_6[k]
                  + pa_y[k] * li_10[k];

        t_39[k] = pb_z[k] * mh_27[k];

        t_40[k] = f_10 * lh_8[k]
                  + pa_y[k] * li_12[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_y, pb_x, pb_y, pb_z, lh_9, lh_36, li_14, \
                         mh_30, mh_31, mh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_9 * lh_9[k]
                  + pb_y[k] * mh_30[k];

        t_42[k] = pa_y[k] * li_14[k];

        t_43[k] = f_13 * lh_36[k]
                  + pb_x[k] * mh_36[k];

        t_44[k] = pb_z[k] * mh_31[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_y, pb_x, lh_15, lh_38, lh_39, lh_40, \
                         li_20, li_21, mh_38, mh_39, mh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_13 * lh_38[k]
                  + pb_x[k] * mh_38[k];

        t_46[k] = f_13 * lh_39[k]
                  + pb_x[k] * mh_39[k];

        t_47[k] = f_13 * lh_40[k]
                  + pb_x[k] * mh_40[k];

        t_48[k] = pa_y[k] * li_20[k];

        t_49[k] = f_14 * lh_15[k]
                  + pa_y[k] * li_21[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pb_z, lh_17, lh_18, lh_19, li_23, \
                         li_24, li_25, mh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pb_z[k] * mh_36[k];

        t_51[k] = f_12 * lh_17[k]
                  + pa_y[k] * li_23[k];

        t_52[k] = f_11 * lh_18[k]
                  + pa_y[k] * li_24[k];

        t_53[k] = f_10 * lh_19[k]
                  + pa_y[k] * li_25[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pa_y, pa_z, pb_y, pb_z, lh_0, lh_20, \
                         li_0, li_27, mh_41, mh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * lh_20[k]
                  + pb_y[k] * mh_41[k];

        t_55[k] = pa_y[k] * li_27[k];

        t_56[k] = pa_z[k] * li_0[k];

        t_57[k] = pb_y[k] * mh_42[k];

        t_58[k] = f_9 * lh_0[k]
                  + pb_z[k] * mh_42[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, t_63, pa_z, pb_y, pb_z, lh_2, lh_3, li_3, \
                         li_5, li_6, mh_44, mh_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_z[k] * li_3[k];

        t_60[k] = pb_y[k] * mh_44[k];

        t_61[k] = f_10 * lh_2[k]
                  + pa_z[k] * li_5[k];

        t_62[k] = pa_z[k] * li_6[k];

        t_63[k] = f_9 * lh_3[k]
                  + pb_z[k] * mh_45[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, pa_z, pb_y, pb_z, lh_5, lh_6, lh_7, \
                         li_9, li_10, li_12, mh_47, mh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_y[k] * mh_47[k];

        t_65[k] = f_11 * lh_5[k]
                  + pa_z[k] * li_9[k];

        t_66[k] = pa_z[k] * li_10[k];

        t_67[k] = f_9 * lh_6[k]
                  + pb_z[k] * mh_48[k];

        t_68[k] = f_10 * lh_7[k]
                  + pa_z[k] * li_12[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_z, pb_x, pb_y, lh_9, lh_58, lh_59, \
                         li_14, li_15, mh_51, mh_58, mh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_y[k] * mh_51[k];

        t_70[k] = f_12 * lh_9[k]
                  + pa_z[k] * li_14[k];

        t_71[k] = pa_z[k] * li_15[k];

        t_72[k] = f_13 * lh_58[k]
                  + pb_x[k] * mh_58[k];

        t_73[k] = f_13 * lh_59[k]
                  + pb_x[k] * mh_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pb_x, pb_y, lh_60, lh_62, li_21, mh_56, \
                         mh_60, mh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * lh_60[k]
                  + pb_x[k] * mh_60[k];

        t_75[k] = pb_y[k] * mh_56[k];

        t_76[k] = f_13 * lh_62[k]
                  + pb_x[k] * mh_62[k];

        t_77[k] = pa_z[k] * li_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_z, pb_z, lh_15, lh_16, lh_17, lh_18, \
                         li_23, li_24, li_25, mh_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_9 * lh_15[k]
                  + pb_z[k] * mh_57[k];

        t_79[k] = f_10 * lh_16[k]
                  + pa_z[k] * li_23[k];

        t_80[k] = f_11 * lh_17[k]
                  + pa_z[k] * li_24[k];

        t_81[k] = f_12 * lh_18[k]
                  + pa_z[k] * li_25[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_y, pa_z, pb_y, ki0_0, ki1_0, lh_20, lh_21, \
                         li_27, li_28, mh_62, mh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pb_y[k] * mh_62[k];

        t_83[k] = f_14 * lh_20[k]
                  + pa_z[k] * li_27[k];

        t_84[k] = f_15 * ki0_0[k]
                  - f_16 * ki1_0[k]
                  + pa_y[k] * li_28[k];

        t_85[k] = f_10 * lh_21[k]
                  + pb_y[k] * mh_63[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_x, pb_z, lh_66, mg0_45, mg0_48, mg1_45, \
                         mg1_48, mh_63, mh_64, mh_65, mh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_z[k] * mh_63[k];

        t_87[k] = f_17 * lh_66[k]
                  + f_7 * mg0_48[k]
                  - f_8 * mg1_48[k]
                  + pb_x[k] * mh_66[k];

        t_88[k] = pb_z[k] * mh_64[k];

        t_89[k] = f_3 * mg0_45[k]
                  - f_4 * mg1_45[k]
                  + pb_z[k] * mh_65[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pb_x, pb_y, pb_z, lh_26, lh_69, mg0_47, \
                         mg0_51, mg1_47, mg1_51, mh_66, mh_68, mh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_17 * lh_69[k]
                  + f_5 * mg0_51[k]
                  - f_6 * mg1_51[k]
                  + pb_x[k] * mh_69[k];

        t_91[k] = pb_z[k] * mh_66[k];

        t_92[k] = f_10 * lh_26[k]
                  + pb_y[k] * mh_68[k];

        t_93[k] = f_5 * mg0_47[k]
                  - f_6 * mg1_47[k]
                  + pb_z[k] * mh_68[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pb_x, pb_z, lh_73, mg0_48, mg0_55, mg1_48, mg1_55, \
                         mh_69, mh_70, mh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_17 * lh_73[k]
                  + f_3 * mg0_55[k]
                  - f_4 * mg1_55[k]
                  + pb_x[k] * mh_73[k];

        t_95[k] = pb_z[k] * mh_69[k];

        t_96[k] = f_3 * mg0_48[k]
                  - f_4 * mg1_48[k]
                  + pb_z[k] * mh_70[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_x, pb_y, pb_z, lh_30, lh_78, mg0_50, \
                         mg1_50, mh_72, mh_73, mh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_10 * lh_30[k]
                  + pb_y[k] * mh_72[k];

        t_98[k] = f_7 * mg0_50[k]
                  - f_8 * mg1_50[k]
                  + pb_z[k] * mh_72[k];

        t_99[k] = f_17 * lh_78[k]
                  + pb_x[k] * mh_78[k];

        t_100[k] = pb_z[k] * mh_73[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pb_x, lh_80, lh_81, lh_82, lh_83, mh_80, \
                         mh_81, mh_82, mh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_17 * lh_80[k]
                   + pb_x[k] * mh_80[k];

        t_102[k] = f_17 * lh_81[k]
                   + pb_x[k] * mh_81[k];

        t_103[k] = f_17 * lh_82[k]
                   + pb_x[k] * mh_82[k];

        t_104[k] = f_17 * lh_83[k]
                   + pb_x[k] * mh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_x, pb_z, ki0_105, ki1_105, li_105, \
                         mg0_55, mg0_56, mg1_55, mg1_56, mh_78, mh_79, \
                         mh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_18 * ki0_105[k]
                   - f_19 * ki1_105[k]
                   + pa_x[k] * li_105[k];

        t_106[k] = pb_z[k] * mh_78[k];

        t_107[k] = f_3 * mg0_55[k]
                   - f_4 * mg1_55[k]
                   + pb_z[k] * mh_79[k];

        t_108[k] = f_5 * mg0_56[k]
                   - f_6 * mg1_56[k]
                   + pb_z[k] * mh_80[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pb_y, pb_z, lh_41, li_56, mg0_57, \
                         mg0_59, mg1_57, mg1_59, mh_81, mh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_7 * mg0_57[k]
                   - f_8 * mg1_57[k]
                   + pb_z[k] * mh_81[k];

        t_110[k] = f_10 * lh_41[k]
                   + pb_y[k] * mh_83[k];

        t_111[k] = f_1 * mg0_59[k]
                   - f_2 * mg1_59[k]
                   + pb_z[k] * mh_83[k];

        t_112[k] = pa_y[k] * li_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, t_118, pa_y, pa_z, pb_y, lh_44, \
                         li_29, li_31, li_34, li_58, li_61, mh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_z[k] * li_29[k];

        t_114[k] = pa_y[k] * li_58[k];

        t_115[k] = pa_z[k] * li_31[k];

        t_116[k] = f_9 * lh_44[k]
                   + pb_y[k] * mh_86[k];

        t_117[k] = pa_y[k] * li_61[k];

        t_118[k] = pa_z[k] * li_34[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pa_y, pa_z, pb_y, pb_z, lh_24, lh_47, \
                         li_38, li_65, mh_87, mh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_9 * lh_24[k]
                   + pb_z[k] * mh_87[k];

        t_120[k] = f_9 * lh_47[k]
                   + pb_y[k] * mh_89[k];

        t_121[k] = pa_y[k] * li_65[k];

        t_122[k] = pa_z[k] * li_38[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_y, pb_y, pb_z, lh_27, lh_50, lh_51, \
                         li_68, li_70, mh_90, mh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_9 * lh_27[k]
                   + pb_z[k] * mh_90[k];

        t_124[k] = f_10 * lh_50[k]
                   + pa_y[k] * li_68[k];

        t_125[k] = f_9 * lh_51[k]
                   + pb_y[k] * mh_93[k];

        t_126[k] = pa_y[k] * li_70[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pa_z, pb_x, lh_100, lh_101, \
                         lh_102, lh_103, li_43, mh_100, mh_101, mh_102, \
                         mh_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = pa_z[k] * li_43[k];

        t_128[k] = f_17 * lh_100[k]
                   + pb_x[k] * mh_100[k];

        t_129[k] = f_17 * lh_101[k]
                   + pb_x[k] * mh_101[k];

        t_130[k] = f_17 * lh_102[k]
                   + pb_x[k] * mh_102[k];

        t_131[k] = f_17 * lh_103[k]
                   + pb_x[k] * mh_103[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, pa_y, pa_z, pb_z, lh_36, lh_59, \
                         lh_60, li_49, li_76, li_79, li_80, mh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_y[k] * li_76[k];

        t_133[k] = pa_z[k] * li_49[k];

        t_134[k] = f_9 * lh_36[k]
                   + pb_z[k] * mh_99[k];

        t_135[k] = f_12 * lh_59[k]
                   + pa_y[k] * li_79[k];

        t_136[k] = f_11 * lh_60[k]
                   + pa_y[k] * li_80[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pa_y, pa_z, pb_y, ki0_0, ki1_0, lh_61, \
                         lh_62, li_56, li_81, li_83, mh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_10 * lh_61[k]
                   + pa_y[k] * li_81[k];

        t_138[k] = f_9 * lh_62[k]
                   + pb_y[k] * mh_104[k];

        t_139[k] = pa_y[k] * li_83[k];

        t_140[k] = f_15 * ki0_0[k]
                   - f_16 * ki1_0[k]
                   + pa_z[k] * li_56[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pb_y, pb_z, lh_42, mg0_75, mg1_75, \
                         mh_105, mh_106, mh_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = pb_y[k] * mh_105[k];

        t_142[k] = f_10 * lh_42[k]
                   + pb_z[k] * mh_105[k];

        t_143[k] = f_3 * mg0_75[k]
                   - f_4 * mg1_75[k]
                   + pb_y[k] * mh_106[k];

        t_144[k] = pb_y[k] * mh_107[k];
    }
}

static auto
compute_prim_mi_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t ki0,
                                            const size_t ki1, const size_t lh, const size_t li,
                                            const size_t mg0, const size_t mg1, const size_t mh,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_14 = 3.0 / p;
    const auto f_17 = 3.5 / p;
    const auto f_18 = 3.0 / alpha;
    const auto f_19 = 3.0 * beta / (alpha * p);
    const auto f_20 = 1.0 / alpha;
    const auto f_21 = beta / (alpha * p);
    const auto f_22 = 2.5 / alpha;
    const auto f_23 = 2.5 * beta / (alpha * p);

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
    auto *t_276 = buffer.data(target + 276);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ki0_28 = buffer.data(ki0 + 28);
    const auto *ki0_56 = buffer.data(ki0 + 56);
    const auto *ki0_167 = buffer.data(ki0 + 167);
    const auto *ki0_189 = buffer.data(ki0 + 189);

    const auto *ki1_28 = buffer.data(ki1 + 28);
    const auto *ki1_56 = buffer.data(ki1 + 56);
    const auto *ki1_167 = buffer.data(ki1 + 167);
    const auto *ki1_189 = buffer.data(ki1 + 189);

    const auto *lh_45 = buffer.data(lh + 45);
    const auto *lh_48 = buffer.data(lh + 48);
    const auto *lh_57 = buffer.data(lh + 57);
    const auto *lh_63 = buffer.data(lh + 63);
    const auto *lh_65 = buffer.data(lh + 65);
    const auto *lh_66 = buffer.data(lh + 66);
    const auto *lh_68 = buffer.data(lh + 68);
    const auto *lh_69 = buffer.data(lh + 69);
    const auto *lh_70 = buffer.data(lh + 70);
    const auto *lh_72 = buffer.data(lh + 72);
    const auto *lh_78 = buffer.data(lh + 78);
    const auto *lh_79 = buffer.data(lh + 79);
    const auto *lh_80 = buffer.data(lh + 80);
    const auto *lh_81 = buffer.data(lh + 81);
    const auto *lh_83 = buffer.data(lh + 83);
    const auto *lh_86 = buffer.data(lh + 86);
    const auto *lh_87 = buffer.data(lh + 87);
    const auto *lh_89 = buffer.data(lh + 89);
    const auto *lh_90 = buffer.data(lh + 90);
    const auto *lh_93 = buffer.data(lh + 93);
    const auto *lh_99 = buffer.data(lh + 99);
    const auto *lh_104 = buffer.data(lh + 104);
    const auto *lh_105 = buffer.data(lh + 105);
    const auto *lh_106 = buffer.data(lh + 106);
    const auto *lh_107 = buffer.data(lh + 107);
    const auto *lh_108 = buffer.data(lh + 108);
    const auto *lh_110 = buffer.data(lh + 110);
    const auto *lh_111 = buffer.data(lh + 111);
    const auto *lh_113 = buffer.data(lh + 113);
    const auto *lh_114 = buffer.data(lh + 114);
    const auto *lh_119 = buffer.data(lh + 119);
    const auto *lh_120 = buffer.data(lh + 120);
    const auto *lh_121 = buffer.data(lh + 121);
    const auto *lh_122 = buffer.data(lh + 122);
    const auto *lh_123 = buffer.data(lh + 123);
    const auto *lh_124 = buffer.data(lh + 124);
    const auto *lh_125 = buffer.data(lh + 125);
    const auto *lh_129 = buffer.data(lh + 129);
    const auto *lh_132 = buffer.data(lh + 132);
    const auto *lh_136 = buffer.data(lh + 136);
    const auto *lh_141 = buffer.data(lh + 141);
    const auto *lh_143 = buffer.data(lh + 143);
    const auto *lh_144 = buffer.data(lh + 144);
    const auto *lh_145 = buffer.data(lh + 145);
    const auto *lh_146 = buffer.data(lh + 146);
    const auto *lh_163 = buffer.data(lh + 163);
    const auto *lh_164 = buffer.data(lh + 164);
    const auto *lh_165 = buffer.data(lh + 165);
    const auto *lh_166 = buffer.data(lh + 166);
    const auto *lh_167 = buffer.data(lh + 167);
    const auto *lh_183 = buffer.data(lh + 183);
    const auto *lh_184 = buffer.data(lh + 184);
    const auto *lh_185 = buffer.data(lh + 185);
    const auto *lh_186 = buffer.data(lh + 186);
    const auto *lh_187 = buffer.data(lh + 187);
    const auto *lh_194 = buffer.data(lh + 194);
    const auto *lh_198 = buffer.data(lh + 198);
    const auto *lh_203 = buffer.data(lh + 203);
    const auto *lh_204 = buffer.data(lh + 204);
    const auto *lh_205 = buffer.data(lh + 205);
    const auto *lh_206 = buffer.data(lh + 206);
    const auto *lh_207 = buffer.data(lh + 207);
    const auto *lh_209 = buffer.data(lh + 209);

    const auto *li_84 = buffer.data(li + 84);
    const auto *li_85 = buffer.data(li + 85);
    const auto *li_87 = buffer.data(li + 87);
    const auto *li_89 = buffer.data(li + 89);
    const auto *li_90 = buffer.data(li + 90);
    const auto *li_93 = buffer.data(li + 93);
    const auto *li_94 = buffer.data(li + 94);
    const auto *li_96 = buffer.data(li + 96);
    const auto *li_98 = buffer.data(li + 98);
    const auto *li_99 = buffer.data(li + 99);
    const auto *li_105 = buffer.data(li + 105);
    const auto *li_107 = buffer.data(li + 107);
    const auto *li_108 = buffer.data(li + 108);
    const auto *li_109 = buffer.data(li + 109);
    const auto *li_111 = buffer.data(li + 111);
    const auto *li_140 = buffer.data(li + 140);
    const auto *li_142 = buffer.data(li + 142);
    const auto *li_143 = buffer.data(li + 143);
    const auto *li_145 = buffer.data(li + 145);
    const auto *li_146 = buffer.data(li + 146);
    const auto *li_149 = buffer.data(li + 149);
    const auto *li_150 = buffer.data(li + 150);
    const auto *li_152 = buffer.data(li + 152);
    const auto *li_154 = buffer.data(li + 154);
    const auto *li_160 = buffer.data(li + 160);
    const auto *li_161 = buffer.data(li + 161);
    const auto *li_163 = buffer.data(li + 163);
    const auto *li_164 = buffer.data(li + 164);
    const auto *li_165 = buffer.data(li + 165);
    const auto *li_167 = buffer.data(li + 167);
    const auto *li_189 = buffer.data(li + 189);

    const auto *mg0_76 = buffer.data(mg0 + 76);
    const auto *mg0_78 = buffer.data(mg0 + 78);
    const auto *mg0_80 = buffer.data(mg0 + 80);
    const auto *mg0_84 = buffer.data(mg0 + 84);
    const auto *mg0_85 = buffer.data(mg0 + 85);
    const auto *mg0_87 = buffer.data(mg0 + 87);
    const auto *mg0_88 = buffer.data(mg0 + 88);
    const auto *mg0_89 = buffer.data(mg0 + 89);
    const auto *mg0_90 = buffer.data(mg0 + 90);
    const auto *mg0_92 = buffer.data(mg0 + 92);
    const auto *mg0_93 = buffer.data(mg0 + 93);
    const auto *mg0_95 = buffer.data(mg0 + 95);
    const auto *mg0_96 = buffer.data(mg0 + 96);
    const auto *mg0_100 = buffer.data(mg0 + 100);
    const auto *mg0_101 = buffer.data(mg0 + 101);
    const auto *mg0_102 = buffer.data(mg0 + 102);
    const auto *mg0_104 = buffer.data(mg0 + 104);
    const auto *mg0_135 = buffer.data(mg0 + 135);
    const auto *mg0_136 = buffer.data(mg0 + 136);
    const auto *mg0_138 = buffer.data(mg0 + 138);
    const auto *mg0_140 = buffer.data(mg0 + 140);
    const auto *mg0_144 = buffer.data(mg0 + 144);
    const auto *mg0_145 = buffer.data(mg0 + 145);
    const auto *mg0_147 = buffer.data(mg0 + 147);
    const auto *mg0_148 = buffer.data(mg0 + 148);
    const auto *mg0_149 = buffer.data(mg0 + 149);

    const auto *mg1_76 = buffer.data(mg1 + 76);
    const auto *mg1_78 = buffer.data(mg1 + 78);
    const auto *mg1_80 = buffer.data(mg1 + 80);
    const auto *mg1_84 = buffer.data(mg1 + 84);
    const auto *mg1_85 = buffer.data(mg1 + 85);
    const auto *mg1_87 = buffer.data(mg1 + 87);
    const auto *mg1_88 = buffer.data(mg1 + 88);
    const auto *mg1_89 = buffer.data(mg1 + 89);
    const auto *mg1_90 = buffer.data(mg1 + 90);
    const auto *mg1_92 = buffer.data(mg1 + 92);
    const auto *mg1_93 = buffer.data(mg1 + 93);
    const auto *mg1_95 = buffer.data(mg1 + 95);
    const auto *mg1_96 = buffer.data(mg1 + 96);
    const auto *mg1_100 = buffer.data(mg1 + 100);
    const auto *mg1_101 = buffer.data(mg1 + 101);
    const auto *mg1_102 = buffer.data(mg1 + 102);
    const auto *mg1_104 = buffer.data(mg1 + 104);
    const auto *mg1_135 = buffer.data(mg1 + 135);
    const auto *mg1_136 = buffer.data(mg1 + 136);
    const auto *mg1_138 = buffer.data(mg1 + 138);
    const auto *mg1_140 = buffer.data(mg1 + 140);
    const auto *mg1_144 = buffer.data(mg1 + 144);
    const auto *mg1_145 = buffer.data(mg1 + 145);
    const auto *mg1_147 = buffer.data(mg1 + 147);
    const auto *mg1_148 = buffer.data(mg1 + 148);
    const auto *mg1_149 = buffer.data(mg1 + 149);

    const auto *mh_108 = buffer.data(mh + 108);
    const auto *mh_110 = buffer.data(mh + 110);
    const auto *mh_111 = buffer.data(mh + 111);
    const auto *mh_113 = buffer.data(mh + 113);
    const auto *mh_114 = buffer.data(mh + 114);
    const auto *mh_119 = buffer.data(mh + 119);
    const auto *mh_120 = buffer.data(mh + 120);
    const auto *mh_121 = buffer.data(mh + 121);
    const auto *mh_122 = buffer.data(mh + 122);
    const auto *mh_123 = buffer.data(mh + 123);
    const auto *mh_124 = buffer.data(mh + 124);
    const auto *mh_125 = buffer.data(mh + 125);
    const auto *mh_126 = buffer.data(mh + 126);
    const auto *mh_127 = buffer.data(mh + 127);
    const auto *mh_128 = buffer.data(mh + 128);
    const auto *mh_129 = buffer.data(mh + 129);
    const auto *mh_131 = buffer.data(mh + 131);
    const auto *mh_132 = buffer.data(mh + 132);
    const auto *mh_133 = buffer.data(mh + 133);
    const auto *mh_135 = buffer.data(mh + 135);
    const auto *mh_136 = buffer.data(mh + 136);
    const auto *mh_141 = buffer.data(mh + 141);
    const auto *mh_142 = buffer.data(mh + 142);
    const auto *mh_143 = buffer.data(mh + 143);
    const auto *mh_144 = buffer.data(mh + 144);
    const auto *mh_145 = buffer.data(mh + 145);
    const auto *mh_146 = buffer.data(mh + 146);
    const auto *mh_147 = buffer.data(mh + 147);
    const auto *mh_149 = buffer.data(mh + 149);
    const auto *mh_150 = buffer.data(mh + 150);
    const auto *mh_152 = buffer.data(mh + 152);
    const auto *mh_153 = buffer.data(mh + 153);
    const auto *mh_156 = buffer.data(mh + 156);
    const auto *mh_162 = buffer.data(mh + 162);
    const auto *mh_163 = buffer.data(mh + 163);
    const auto *mh_164 = buffer.data(mh + 164);
    const auto *mh_165 = buffer.data(mh + 165);
    const auto *mh_166 = buffer.data(mh + 166);
    const auto *mh_167 = buffer.data(mh + 167);
    const auto *mh_168 = buffer.data(mh + 168);
    const auto *mh_170 = buffer.data(mh + 170);
    const auto *mh_171 = buffer.data(mh + 171);
    const auto *mh_173 = buffer.data(mh + 173);
    const auto *mh_174 = buffer.data(mh + 174);
    const auto *mh_177 = buffer.data(mh + 177);
    const auto *mh_183 = buffer.data(mh + 183);
    const auto *mh_184 = buffer.data(mh + 184);
    const auto *mh_185 = buffer.data(mh + 185);
    const auto *mh_186 = buffer.data(mh + 186);
    const auto *mh_187 = buffer.data(mh + 187);
    const auto *mh_188 = buffer.data(mh + 188);
    const auto *mh_189 = buffer.data(mh + 189);
    const auto *mh_190 = buffer.data(mh + 190);
    const auto *mh_191 = buffer.data(mh + 191);
    const auto *mh_192 = buffer.data(mh + 192);
    const auto *mh_194 = buffer.data(mh + 194);
    const auto *mh_195 = buffer.data(mh + 195);
    const auto *mh_197 = buffer.data(mh + 197);
    const auto *mh_198 = buffer.data(mh + 198);
    const auto *mh_203 = buffer.data(mh + 203);
    const auto *mh_204 = buffer.data(mh + 204);
    const auto *mh_205 = buffer.data(mh + 205);
    const auto *mh_206 = buffer.data(mh + 206);
    const auto *mh_207 = buffer.data(mh + 207);
    const auto *mh_209 = buffer.data(mh + 209);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pb_x, pb_y, pb_z, lh_45, lh_110, mg0_76, \
                         mg0_80, mg1_76, mg1_80, mh_108, mh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_17 * lh_110[k]
                   + f_7 * mg0_80[k]
                   - f_8 * mg1_80[k]
                   + pb_x[k] * mh_110[k];

        t_146[k] = f_5 * mg0_76[k]
                   - f_6 * mg1_76[k]
                   + pb_y[k] * mh_108[k];

        t_147[k] = f_10 * lh_45[k]
                   + pb_z[k] * mh_108[k];

        t_148[k] = pb_y[k] * mh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pb_x, pb_y, pb_z, lh_48, lh_114, mg0_78, mg0_84, \
                         mg1_78, mg1_84, mh_111, mh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_17 * lh_114[k]
                   + f_5 * mg0_84[k]
                   - f_6 * mg1_84[k]
                   + pb_x[k] * mh_114[k];

        t_150[k] = f_7 * mg0_78[k]
                   - f_8 * mg1_78[k]
                   + pb_y[k] * mh_111[k];

        t_151[k] = f_10 * lh_48[k]
                   + pb_z[k] * mh_111[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pb_y, lh_119, lh_120, mg0_80, \
                         mg0_89, mg1_80, mg1_89, mh_113, mh_114, mh_119, \
                         mh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_3 * mg0_80[k]
                   - f_4 * mg1_80[k]
                   + pb_y[k] * mh_113[k];

        t_153[k] = pb_y[k] * mh_114[k];

        t_154[k] = f_17 * lh_119[k]
                   + f_3 * mg0_89[k]
                   - f_4 * mg1_89[k]
                   + pb_x[k] * mh_119[k];

        t_155[k] = f_17 * lh_120[k]
                   + pb_x[k] * mh_120[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pb_x, pb_y, lh_121, lh_122, \
                         lh_123, lh_125, mh_119, mh_121, mh_122, mh_123, \
                         mh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_17 * lh_121[k]
                   + pb_x[k] * mh_121[k];

        t_157[k] = f_17 * lh_122[k]
                   + pb_x[k] * mh_122[k];

        t_158[k] = f_17 * lh_123[k]
                   + pb_x[k] * mh_123[k];

        t_159[k] = pb_y[k] * mh_119[k];

        t_160[k] = f_17 * lh_125[k]
                   + pb_x[k] * mh_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_y, pb_z, lh_57, mg0_85, mg0_87, \
                         mg0_88, mg1_85, mg1_87, mg1_88, mh_120, mh_122, \
                         mh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_1 * mg0_85[k]
                   - f_2 * mg1_85[k]
                   + pb_y[k] * mh_120[k];

        t_162[k] = f_10 * lh_57[k]
                   + pb_z[k] * mh_120[k];

        t_163[k] = f_7 * mg0_87[k]
                   - f_8 * mg1_87[k]
                   + pb_y[k] * mh_122[k];

        t_164[k] = f_5 * mg0_88[k]
                   - f_6 * mg1_88[k]
                   + pb_y[k] * mh_123[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pa_x, pb_y, ki0_167, ki1_167, li_167, mg0_89, \
                         mg1_89, mh_124, mh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_3 * mg0_89[k]
                   - f_4 * mg1_89[k]
                   + pb_y[k] * mh_124[k];

        t_166[k] = pb_y[k] * mh_125[k];

        t_167[k] = f_18 * ki0_167[k]
                   - f_19 * ki1_167[k]
                   + pa_x[k] * li_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_y, pb_y, pb_z, ki0_28, ki1_28, lh_63, li_84, \
                         mh_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_20 * ki0_28[k]
                   - f_21 * ki1_28[k]
                   + pa_y[k] * li_84[k];

        t_169[k] = f_11 * lh_63[k]
                   + pb_y[k] * mh_126[k];

        t_170[k] = pb_z[k] * mh_126[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pb_x, pb_z, lh_129, mg0_90, mg0_93, mg1_90, \
                         mg1_93, mh_127, mh_128, mh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_14 * lh_129[k]
                   + f_7 * mg0_93[k]
                   - f_8 * mg1_93[k]
                   + pb_x[k] * mh_129[k];

        t_172[k] = pb_z[k] * mh_127[k];

        t_173[k] = f_3 * mg0_90[k]
                   - f_4 * mg1_90[k]
                   + pb_z[k] * mh_128[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pb_x, pb_y, pb_z, lh_68, lh_132, mg0_92, \
                         mg0_96, mg1_92, mg1_96, mh_129, mh_131, \
                         mh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_14 * lh_132[k]
                   + f_5 * mg0_96[k]
                   - f_6 * mg1_96[k]
                   + pb_x[k] * mh_132[k];

        t_175[k] = pb_z[k] * mh_129[k];

        t_176[k] = f_11 * lh_68[k]
                   + pb_y[k] * mh_131[k];

        t_177[k] = f_5 * mg0_92[k]
                   - f_6 * mg1_92[k]
                   + pb_z[k] * mh_131[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pb_x, pb_z, lh_136, mg0_93, mg0_100, mg1_93, \
                         mg1_100, mh_132, mh_133, mh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * lh_136[k]
                   + f_3 * mg0_100[k]
                   - f_4 * mg1_100[k]
                   + pb_x[k] * mh_136[k];

        t_179[k] = pb_z[k] * mh_132[k];

        t_180[k] = f_3 * mg0_93[k]
                   - f_4 * mg1_93[k]
                   + pb_z[k] * mh_133[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_x, pb_y, pb_z, lh_72, lh_141, mg0_95, \
                         mg1_95, mh_135, mh_136, mh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_11 * lh_72[k]
                   + pb_y[k] * mh_135[k];

        t_182[k] = f_7 * mg0_95[k]
                   - f_8 * mg1_95[k]
                   + pb_z[k] * mh_135[k];

        t_183[k] = f_14 * lh_141[k]
                   + pb_x[k] * mh_141[k];

        t_184[k] = pb_z[k] * mh_136[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_x, lh_143, lh_144, lh_145, lh_146, \
                         mh_143, mh_144, mh_145, mh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_14 * lh_143[k]
                   + pb_x[k] * mh_143[k];

        t_186[k] = f_14 * lh_144[k]
                   + pb_x[k] * mh_144[k];

        t_187[k] = f_14 * lh_145[k]
                   + pb_x[k] * mh_145[k];

        t_188[k] = f_14 * lh_146[k]
                   + pb_x[k] * mh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_x, pb_z, ki0_189, ki1_189, li_189, \
                         mg0_100, mg0_101, mg1_100, mg1_101, mh_141, mh_142, \
                         mh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_22 * ki0_189[k]
                   - f_23 * ki1_189[k]
                   + pa_x[k] * li_189[k];

        t_190[k] = pb_z[k] * mh_141[k];

        t_191[k] = f_3 * mg0_100[k]
                   - f_4 * mg1_100[k]
                   + pb_z[k] * mh_142[k];

        t_192[k] = f_5 * mg0_101[k]
                   - f_6 * mg1_101[k]
                   + pb_z[k] * mh_143[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pb_y, pb_z, lh_83, li_84, mg0_102, \
                         mg0_104, mg1_102, mg1_104, mh_144, mh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_7 * mg0_102[k]
                   - f_8 * mg1_102[k]
                   + pb_z[k] * mh_144[k];

        t_194[k] = f_11 * lh_83[k]
                   + pb_y[k] * mh_146[k];

        t_195[k] = f_1 * mg0_104[k]
                   - f_2 * mg1_104[k]
                   + pb_z[k] * mh_146[k];

        t_196[k] = pa_z[k] * li_84[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, pa_z, pb_y, pb_z, lh_63, lh_65, \
                         lh_86, li_85, li_87, li_89, mh_147, mh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = pa_z[k] * li_85[k];

        t_198[k] = f_9 * lh_63[k]
                   + pb_z[k] * mh_147[k];

        t_199[k] = pa_z[k] * li_87[k];

        t_200[k] = f_10 * lh_86[k]
                   + pb_y[k] * mh_149[k];

        t_201[k] = f_10 * lh_65[k]
                   + pa_z[k] * li_89[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, pa_z, pb_y, pb_z, lh_66, lh_68, \
                         lh_89, li_90, li_93, li_94, mh_150, mh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = pa_z[k] * li_90[k];

        t_203[k] = f_9 * lh_66[k]
                   + pb_z[k] * mh_150[k];

        t_204[k] = f_10 * lh_89[k]
                   + pb_y[k] * mh_152[k];

        t_205[k] = f_11 * lh_68[k]
                   + pa_z[k] * li_93[k];

        t_206[k] = pa_z[k] * li_94[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, pa_z, pb_y, pb_z, lh_69, lh_70, lh_72, \
                         lh_93, li_96, li_98, mh_153, mh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_9 * lh_69[k]
                   + pb_z[k] * mh_153[k];

        t_208[k] = f_10 * lh_70[k]
                   + pa_z[k] * li_96[k];

        t_209[k] = f_10 * lh_93[k]
                   + pb_y[k] * mh_156[k];

        t_210[k] = f_12 * lh_72[k]
                   + pa_z[k] * li_98[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, pa_z, pb_x, lh_163, lh_164, \
                         lh_165, lh_166, li_99, mh_163, mh_164, mh_165, \
                         mh_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = pa_z[k] * li_99[k];

        t_212[k] = f_14 * lh_163[k]
                   + pb_x[k] * mh_163[k];

        t_213[k] = f_14 * lh_164[k]
                   + pb_x[k] * mh_164[k];

        t_214[k] = f_14 * lh_165[k]
                   + pb_x[k] * mh_165[k];

        t_215[k] = f_14 * lh_166[k]
                   + pb_x[k] * mh_166[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pa_z, pb_x, pb_z, lh_78, lh_79, lh_167, \
                         li_105, li_107, mh_162, mh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_14 * lh_167[k]
                   + pb_x[k] * mh_167[k];

        t_217[k] = pa_z[k] * li_105[k];

        t_218[k] = f_9 * lh_78[k]
                   + pb_z[k] * mh_162[k];

        t_219[k] = f_10 * lh_79[k]
                   + pa_z[k] * li_107[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pa_z, pb_y, lh_80, lh_81, lh_83, lh_104, \
                         li_108, li_109, li_111, mh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_11 * lh_80[k]
                   + pa_z[k] * li_108[k];

        t_221[k] = f_12 * lh_81[k]
                   + pa_z[k] * li_109[k];

        t_222[k] = f_10 * lh_104[k]
                   + pb_y[k] * mh_167[k];

        t_223[k] = f_14 * lh_83[k]
                   + pa_z[k] * li_111[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, pa_y, pb_y, lh_105, lh_106, \
                         lh_107, li_140, li_142, li_143, mh_168, \
                         mh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pa_y[k] * li_140[k];

        t_225[k] = f_9 * lh_105[k]
                   + pb_y[k] * mh_168[k];

        t_226[k] = pa_y[k] * li_142[k];

        t_227[k] = f_10 * lh_106[k]
                   + pa_y[k] * li_143[k];

        t_228[k] = f_9 * lh_107[k]
                   + pb_y[k] * mh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, pa_y, pb_y, pb_z, lh_87, lh_108, \
                         lh_110, li_145, li_146, li_149, mh_171, \
                         mh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * li_145[k];

        t_230[k] = f_11 * lh_108[k]
                   + pa_y[k] * li_146[k];

        t_231[k] = f_10 * lh_87[k]
                   + pb_z[k] * mh_171[k];

        t_232[k] = f_9 * lh_110[k]
                   + pb_y[k] * mh_173[k];

        t_233[k] = pa_y[k] * li_149[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pa_y, pb_y, pb_z, lh_90, lh_111, lh_113, \
                         lh_114, li_150, li_152, mh_174, mh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_12 * lh_111[k]
                   + pa_y[k] * li_150[k];

        t_235[k] = f_10 * lh_90[k]
                   + pb_z[k] * mh_174[k];

        t_236[k] = f_10 * lh_113[k]
                   + pa_y[k] * li_152[k];

        t_237[k] = f_9 * lh_114[k]
                   + pb_y[k] * mh_177[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, pa_y, pb_x, lh_183, lh_184, \
                         lh_185, lh_186, li_154, mh_183, mh_184, mh_185, \
                         mh_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = pa_y[k] * li_154[k];

        t_239[k] = f_14 * lh_183[k]
                   + pb_x[k] * mh_183[k];

        t_240[k] = f_14 * lh_184[k]
                   + pb_x[k] * mh_184[k];

        t_241[k] = f_14 * lh_185[k]
                   + pb_x[k] * mh_185[k];

        t_242[k] = f_14 * lh_186[k]
                   + pb_x[k] * mh_186[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pa_y, pb_x, pb_z, lh_99, lh_120, lh_187, \
                         li_160, li_161, mh_183, mh_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_14 * lh_187[k]
                   + pb_x[k] * mh_187[k];

        t_244[k] = pa_y[k] * li_160[k];

        t_245[k] = f_14 * lh_120[k]
                   + pa_y[k] * li_161[k];

        t_246[k] = f_10 * lh_99[k]
                   + pb_z[k] * mh_183[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, pa_y, pb_y, lh_122, lh_123, \
                         lh_124, lh_125, li_163, li_164, li_165, li_167, \
                         mh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_12 * lh_122[k]
                   + pa_y[k] * li_163[k];

        t_248[k] = f_11 * lh_123[k]
                   + pa_y[k] * li_164[k];

        t_249[k] = f_10 * lh_124[k]
                   + pa_y[k] * li_165[k];

        t_250[k] = f_9 * lh_125[k]
                   + pb_y[k] * mh_188[k];

        t_251[k] = pa_y[k] * li_167[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_z, pb_y, pb_z, ki0_56, ki1_56, lh_105, \
                         li_140, mg0_135, mg1_135, mh_189, mh_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_20 * ki0_56[k]
                   - f_21 * ki1_56[k]
                   + pa_z[k] * li_140[k];

        t_253[k] = pb_y[k] * mh_189[k];

        t_254[k] = f_11 * lh_105[k]
                   + pb_z[k] * mh_189[k];

        t_255[k] = f_3 * mg0_135[k]
                   - f_4 * mg1_135[k]
                   + pb_y[k] * mh_190[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pb_x, pb_y, pb_z, lh_108, lh_194, \
                         mg0_136, mg0_140, mg1_136, mg1_140, mh_191, mh_192, \
                         mh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = pb_y[k] * mh_191[k];

        t_257[k] = f_14 * lh_194[k]
                   + f_7 * mg0_140[k]
                   - f_8 * mg1_140[k]
                   + pb_x[k] * mh_194[k];

        t_258[k] = f_5 * mg0_136[k]
                   - f_6 * mg1_136[k]
                   + pb_y[k] * mh_192[k];

        t_259[k] = f_11 * lh_108[k]
                   + pb_z[k] * mh_192[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pb_x, pb_y, pb_z, lh_111, lh_198, \
                         mg0_138, mg0_144, mg1_138, mg1_144, mh_194, mh_195, \
                         mh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pb_y[k] * mh_194[k];

        t_261[k] = f_14 * lh_198[k]
                   + f_5 * mg0_144[k]
                   - f_6 * mg1_144[k]
                   + pb_x[k] * mh_198[k];

        t_262[k] = f_7 * mg0_138[k]
                   - f_8 * mg1_138[k]
                   + pb_y[k] * mh_195[k];

        t_263[k] = f_11 * lh_111[k]
                   + pb_z[k] * mh_195[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pb_x, pb_y, lh_203, lh_204, mg0_140, \
                         mg0_149, mg1_140, mg1_149, mh_197, mh_198, mh_203, \
                         mh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_3 * mg0_140[k]
                   - f_4 * mg1_140[k]
                   + pb_y[k] * mh_197[k];

        t_265[k] = pb_y[k] * mh_198[k];

        t_266[k] = f_14 * lh_203[k]
                   + f_3 * mg0_149[k]
                   - f_4 * mg1_149[k]
                   + pb_x[k] * mh_203[k];

        t_267[k] = f_14 * lh_204[k]
                   + pb_x[k] * mh_204[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, pb_x, pb_y, lh_205, lh_206, \
                         lh_207, lh_209, mh_203, mh_205, mh_206, mh_207, \
                         mh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_14 * lh_205[k]
                   + pb_x[k] * mh_205[k];

        t_269[k] = f_14 * lh_206[k]
                   + pb_x[k] * mh_206[k];

        t_270[k] = f_14 * lh_207[k]
                   + pb_x[k] * mh_207[k];

        t_271[k] = pb_y[k] * mh_203[k];

        t_272[k] = f_14 * lh_209[k]
                   + pb_x[k] * mh_209[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, pb_y, pb_z, lh_120, mg0_145, mg0_147, \
                         mg0_148, mg1_145, mg1_147, mg1_148, mh_204, mh_206, \
                         mh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_1 * mg0_145[k]
                   - f_2 * mg1_145[k]
                   + pb_y[k] * mh_204[k];

        t_274[k] = f_11 * lh_120[k]
                   + pb_z[k] * mh_204[k];

        t_275[k] = f_7 * mg0_147[k]
                   - f_8 * mg1_147[k]
                   + pb_y[k] * mh_206[k];

        t_276[k] = f_5 * mg0_148[k]
                   - f_6 * mg1_148[k]
                   + pb_y[k] * mh_207[k];
    }
}

static auto
compute_prim_mi_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t ki0,
                                            const size_t ki1, const size_t lh, const size_t li,
                                            const size_t mg0, const size_t mg1, const size_t mh,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_14 = 3.0 / p;
    const auto f_15 = 0.5 / alpha;
    const auto f_16 = 0.5 * beta / (alpha * p);
    const auto f_22 = 2.5 / alpha;
    const auto f_23 = 2.5 * beta / (alpha * p);
    const auto f_24 = 1.5 / alpha;
    const auto f_25 = 1.5 * beta / (alpha * p);
    const auto f_26 = 2.5 / p;
    const auto f_27 = 2.0 / alpha;
    const auto f_28 = 2.0 * beta / (alpha * p);

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
    auto *t_407 = buffer.data(target + 407);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ki0_84 = buffer.data(ki0 + 84);
    const auto *ki0_87 = buffer.data(ki0 + 87);
    const auto *ki0_90 = buffer.data(ki0 + 90);
    const auto *ki0_94 = buffer.data(ki0 + 94);
    const auto *ki0_140 = buffer.data(ki0 + 140);
    const auto *ki0_145 = buffer.data(ki0 + 145);
    const auto *ki0_149 = buffer.data(ki0 + 149);
    const auto *ki0_154 = buffer.data(ki0 + 154);
    const auto *ki0_279 = buffer.data(ki0 + 279);
    const auto *ki0_301 = buffer.data(ki0 + 301);
    const auto *ki0_357 = buffer.data(ki0 + 357);
    const auto *ki0_359 = buffer.data(ki0 + 359);
    const auto *ki0_360 = buffer.data(ki0 + 360);
    const auto *ki0_361 = buffer.data(ki0 + 361);
    const auto *ki0_363 = buffer.data(ki0 + 363);

    const auto *ki1_84 = buffer.data(ki1 + 84);
    const auto *ki1_87 = buffer.data(ki1 + 87);
    const auto *ki1_90 = buffer.data(ki1 + 90);
    const auto *ki1_94 = buffer.data(ki1 + 94);
    const auto *ki1_140 = buffer.data(ki1 + 140);
    const auto *ki1_145 = buffer.data(ki1 + 145);
    const auto *ki1_149 = buffer.data(ki1 + 149);
    const auto *ki1_154 = buffer.data(ki1 + 154);
    const auto *ki1_279 = buffer.data(ki1 + 279);
    const auto *ki1_301 = buffer.data(ki1 + 301);
    const auto *ki1_357 = buffer.data(ki1 + 357);
    const auto *ki1_359 = buffer.data(ki1 + 359);
    const auto *ki1_360 = buffer.data(ki1 + 360);
    const auto *ki1_361 = buffer.data(ki1 + 361);
    const auto *ki1_363 = buffer.data(ki1 + 363);

    const auto *lh_126 = buffer.data(lh + 126);
    const auto *lh_128 = buffer.data(lh + 128);
    const auto *lh_129 = buffer.data(lh + 129);
    const auto *lh_131 = buffer.data(lh + 131);
    const auto *lh_132 = buffer.data(lh + 132);
    const auto *lh_133 = buffer.data(lh + 133);
    const auto *lh_135 = buffer.data(lh + 135);
    const auto *lh_141 = buffer.data(lh + 141);
    const auto *lh_142 = buffer.data(lh + 142);
    const auto *lh_143 = buffer.data(lh + 143);
    const auto *lh_144 = buffer.data(lh + 144);
    const auto *lh_146 = buffer.data(lh + 146);
    const auto *lh_147 = buffer.data(lh + 147);
    const auto *lh_149 = buffer.data(lh + 149);
    const auto *lh_150 = buffer.data(lh + 150);
    const auto *lh_152 = buffer.data(lh + 152);
    const auto *lh_153 = buffer.data(lh + 153);
    const auto *lh_156 = buffer.data(lh + 156);
    const auto *lh_162 = buffer.data(lh + 162);
    const auto *lh_167 = buffer.data(lh + 167);
    const auto *lh_168 = buffer.data(lh + 168);
    const auto *lh_170 = buffer.data(lh + 170);
    const auto *lh_171 = buffer.data(lh + 171);
    const auto *lh_173 = buffer.data(lh + 173);
    const auto *lh_174 = buffer.data(lh + 174);
    const auto *lh_177 = buffer.data(lh + 177);
    const auto *lh_183 = buffer.data(lh + 183);
    const auto *lh_188 = buffer.data(lh + 188);
    const auto *lh_189 = buffer.data(lh + 189);
    const auto *lh_190 = buffer.data(lh + 190);
    const auto *lh_191 = buffer.data(lh + 191);
    const auto *lh_192 = buffer.data(lh + 192);
    const auto *lh_194 = buffer.data(lh + 194);
    const auto *lh_195 = buffer.data(lh + 195);
    const auto *lh_197 = buffer.data(lh + 197);
    const auto *lh_198 = buffer.data(lh + 198);
    const auto *lh_204 = buffer.data(lh + 204);
    const auto *lh_206 = buffer.data(lh + 206);
    const auto *lh_207 = buffer.data(lh + 207);
    const auto *lh_208 = buffer.data(lh + 208);
    const auto *lh_209 = buffer.data(lh + 209);
    const auto *lh_213 = buffer.data(lh + 213);
    const auto *lh_216 = buffer.data(lh + 216);
    const auto *lh_220 = buffer.data(lh + 220);
    const auto *lh_225 = buffer.data(lh + 225);
    const auto *lh_227 = buffer.data(lh + 227);
    const auto *lh_228 = buffer.data(lh + 228);
    const auto *lh_229 = buffer.data(lh + 229);
    const auto *lh_230 = buffer.data(lh + 230);
    const auto *lh_247 = buffer.data(lh + 247);
    const auto *lh_248 = buffer.data(lh + 248);
    const auto *lh_249 = buffer.data(lh + 249);
    const auto *lh_250 = buffer.data(lh + 250);
    const auto *lh_251 = buffer.data(lh + 251);
    const auto *lh_264 = buffer.data(lh + 264);
    const auto *lh_267 = buffer.data(lh + 267);
    const auto *lh_268 = buffer.data(lh + 268);
    const auto *lh_269 = buffer.data(lh + 269);
    const auto *lh_270 = buffer.data(lh + 270);
    const auto *lh_271 = buffer.data(lh + 271);
    const auto *lh_272 = buffer.data(lh + 272);
    const auto *lh_288 = buffer.data(lh + 288);
    const auto *lh_289 = buffer.data(lh + 289);
    const auto *lh_290 = buffer.data(lh + 290);
    const auto *lh_291 = buffer.data(lh + 291);
    const auto *lh_292 = buffer.data(lh + 292);
    const auto *lh_299 = buffer.data(lh + 299);
    const auto *lh_303 = buffer.data(lh + 303);
    const auto *lh_308 = buffer.data(lh + 308);
    const auto *lh_309 = buffer.data(lh + 309);

    const auto *li_168 = buffer.data(li + 168);
    const auto *li_169 = buffer.data(li + 169);
    const auto *li_171 = buffer.data(li + 171);
    const auto *li_173 = buffer.data(li + 173);
    const auto *li_174 = buffer.data(li + 174);
    const auto *li_177 = buffer.data(li + 177);
    const auto *li_178 = buffer.data(li + 178);
    const auto *li_180 = buffer.data(li + 180);
    const auto *li_182 = buffer.data(li + 182);
    const auto *li_183 = buffer.data(li + 183);
    const auto *li_189 = buffer.data(li + 189);
    const auto *li_191 = buffer.data(li + 191);
    const auto *li_192 = buffer.data(li + 192);
    const auto *li_193 = buffer.data(li + 193);
    const auto *li_195 = buffer.data(li + 195);
    const auto *li_199 = buffer.data(li + 199);
    const auto *li_202 = buffer.data(li + 202);
    const auto *li_206 = buffer.data(li + 206);
    const auto *li_224 = buffer.data(li + 224);
    const auto *li_229 = buffer.data(li + 229);
    const auto *li_233 = buffer.data(li + 233);
    const auto *li_238 = buffer.data(li + 238);
    const auto *li_252 = buffer.data(li + 252);
    const auto *li_254 = buffer.data(li + 254);
    const auto *li_255 = buffer.data(li + 255);
    const auto *li_257 = buffer.data(li + 257);
    const auto *li_258 = buffer.data(li + 258);
    const auto *li_261 = buffer.data(li + 261);
    const auto *li_262 = buffer.data(li + 262);
    const auto *li_264 = buffer.data(li + 264);
    const auto *li_266 = buffer.data(li + 266);
    const auto *li_272 = buffer.data(li + 272);
    const auto *li_273 = buffer.data(li + 273);
    const auto *li_275 = buffer.data(li + 275);
    const auto *li_276 = buffer.data(li + 276);
    const auto *li_277 = buffer.data(li + 277);
    const auto *li_279 = buffer.data(li + 279);
    const auto *li_301 = buffer.data(li + 301);
    const auto *li_357 = buffer.data(li + 357);
    const auto *li_359 = buffer.data(li + 359);
    const auto *li_360 = buffer.data(li + 360);
    const auto *li_361 = buffer.data(li + 361);
    const auto *li_363 = buffer.data(li + 363);

    const auto *mg0_149 = buffer.data(mg0 + 149);
    const auto *mg0_150 = buffer.data(mg0 + 150);
    const auto *mg0_152 = buffer.data(mg0 + 152);
    const auto *mg0_153 = buffer.data(mg0 + 153);
    const auto *mg0_155 = buffer.data(mg0 + 155);
    const auto *mg0_156 = buffer.data(mg0 + 156);
    const auto *mg0_160 = buffer.data(mg0 + 160);
    const auto *mg0_161 = buffer.data(mg0 + 161);
    const auto *mg0_162 = buffer.data(mg0 + 162);
    const auto *mg0_164 = buffer.data(mg0 + 164);
    const auto *mg0_192 = buffer.data(mg0 + 192);
    const auto *mg0_210 = buffer.data(mg0 + 210);
    const auto *mg0_211 = buffer.data(mg0 + 211);
    const auto *mg0_213 = buffer.data(mg0 + 213);
    const auto *mg0_215 = buffer.data(mg0 + 215);
    const auto *mg0_219 = buffer.data(mg0 + 219);
    const auto *mg0_224 = buffer.data(mg0 + 224);

    const auto *mg1_149 = buffer.data(mg1 + 149);
    const auto *mg1_150 = buffer.data(mg1 + 150);
    const auto *mg1_152 = buffer.data(mg1 + 152);
    const auto *mg1_153 = buffer.data(mg1 + 153);
    const auto *mg1_155 = buffer.data(mg1 + 155);
    const auto *mg1_156 = buffer.data(mg1 + 156);
    const auto *mg1_160 = buffer.data(mg1 + 160);
    const auto *mg1_161 = buffer.data(mg1 + 161);
    const auto *mg1_162 = buffer.data(mg1 + 162);
    const auto *mg1_164 = buffer.data(mg1 + 164);
    const auto *mg1_192 = buffer.data(mg1 + 192);
    const auto *mg1_210 = buffer.data(mg1 + 210);
    const auto *mg1_211 = buffer.data(mg1 + 211);
    const auto *mg1_213 = buffer.data(mg1 + 213);
    const auto *mg1_215 = buffer.data(mg1 + 215);
    const auto *mg1_219 = buffer.data(mg1 + 219);
    const auto *mg1_224 = buffer.data(mg1 + 224);

    const auto *mh_208 = buffer.data(mh + 208);
    const auto *mh_209 = buffer.data(mh + 209);
    const auto *mh_210 = buffer.data(mh + 210);
    const auto *mh_211 = buffer.data(mh + 211);
    const auto *mh_212 = buffer.data(mh + 212);
    const auto *mh_213 = buffer.data(mh + 213);
    const auto *mh_215 = buffer.data(mh + 215);
    const auto *mh_216 = buffer.data(mh + 216);
    const auto *mh_217 = buffer.data(mh + 217);
    const auto *mh_219 = buffer.data(mh + 219);
    const auto *mh_220 = buffer.data(mh + 220);
    const auto *mh_225 = buffer.data(mh + 225);
    const auto *mh_226 = buffer.data(mh + 226);
    const auto *mh_227 = buffer.data(mh + 227);
    const auto *mh_228 = buffer.data(mh + 228);
    const auto *mh_229 = buffer.data(mh + 229);
    const auto *mh_230 = buffer.data(mh + 230);
    const auto *mh_231 = buffer.data(mh + 231);
    const auto *mh_233 = buffer.data(mh + 233);
    const auto *mh_234 = buffer.data(mh + 234);
    const auto *mh_236 = buffer.data(mh + 236);
    const auto *mh_237 = buffer.data(mh + 237);
    const auto *mh_240 = buffer.data(mh + 240);
    const auto *mh_246 = buffer.data(mh + 246);
    const auto *mh_247 = buffer.data(mh + 247);
    const auto *mh_248 = buffer.data(mh + 248);
    const auto *mh_249 = buffer.data(mh + 249);
    const auto *mh_250 = buffer.data(mh + 250);
    const auto *mh_251 = buffer.data(mh + 251);
    const auto *mh_252 = buffer.data(mh + 252);
    const auto *mh_254 = buffer.data(mh + 254);
    const auto *mh_255 = buffer.data(mh + 255);
    const auto *mh_257 = buffer.data(mh + 257);
    const auto *mh_258 = buffer.data(mh + 258);
    const auto *mh_261 = buffer.data(mh + 261);
    const auto *mh_264 = buffer.data(mh + 264);
    const auto *mh_267 = buffer.data(mh + 267);
    const auto *mh_268 = buffer.data(mh + 268);
    const auto *mh_269 = buffer.data(mh + 269);
    const auto *mh_270 = buffer.data(mh + 270);
    const auto *mh_271 = buffer.data(mh + 271);
    const auto *mh_272 = buffer.data(mh + 272);
    const auto *mh_273 = buffer.data(mh + 273);
    const auto *mh_275 = buffer.data(mh + 275);
    const auto *mh_276 = buffer.data(mh + 276);
    const auto *mh_278 = buffer.data(mh + 278);
    const auto *mh_279 = buffer.data(mh + 279);
    const auto *mh_282 = buffer.data(mh + 282);
    const auto *mh_288 = buffer.data(mh + 288);
    const auto *mh_289 = buffer.data(mh + 289);
    const auto *mh_290 = buffer.data(mh + 290);
    const auto *mh_291 = buffer.data(mh + 291);
    const auto *mh_292 = buffer.data(mh + 292);
    const auto *mh_293 = buffer.data(mh + 293);
    const auto *mh_294 = buffer.data(mh + 294);
    const auto *mh_295 = buffer.data(mh + 295);
    const auto *mh_296 = buffer.data(mh + 296);
    const auto *mh_297 = buffer.data(mh + 297);
    const auto *mh_299 = buffer.data(mh + 299);
    const auto *mh_300 = buffer.data(mh + 300);
    const auto *mh_302 = buffer.data(mh + 302);
    const auto *mh_303 = buffer.data(mh + 303);
    const auto *mh_308 = buffer.data(mh + 308);
    const auto *mh_309 = buffer.data(mh + 309);

#pragma omp simd aligned(t_277, t_278, t_279, pa_x, pb_y, ki0_279, ki1_279, li_279, mg0_149, \
                         mg1_149, mh_208, mh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_3 * mg0_149[k]
                   - f_4 * mg1_149[k]
                   + pb_y[k] * mh_208[k];

        t_278[k] = pb_y[k] * mh_209[k];

        t_279[k] = f_22 * ki0_279[k]
                   - f_23 * ki1_279[k]
                   + pa_x[k] * li_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pa_y, pb_y, pb_z, ki0_84, ki1_84, lh_126, \
                         li_168, mh_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_24 * ki0_84[k]
                   - f_25 * ki1_84[k]
                   + pa_y[k] * li_168[k];

        t_281[k] = f_12 * lh_126[k]
                   + pb_y[k] * mh_210[k];

        t_282[k] = pb_z[k] * mh_210[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pb_x, pb_z, lh_213, mg0_150, mg0_153, mg1_150, \
                         mg1_153, mh_211, mh_212, mh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_26 * lh_213[k]
                   + f_7 * mg0_153[k]
                   - f_8 * mg1_153[k]
                   + pb_x[k] * mh_213[k];

        t_284[k] = pb_z[k] * mh_211[k];

        t_285[k] = f_3 * mg0_150[k]
                   - f_4 * mg1_150[k]
                   + pb_z[k] * mh_212[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_x, pb_y, pb_z, lh_131, lh_216, \
                         mg0_152, mg0_156, mg1_152, mg1_156, mh_213, mh_215, \
                         mh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_26 * lh_216[k]
                   + f_5 * mg0_156[k]
                   - f_6 * mg1_156[k]
                   + pb_x[k] * mh_216[k];

        t_287[k] = pb_z[k] * mh_213[k];

        t_288[k] = f_12 * lh_131[k]
                   + pb_y[k] * mh_215[k];

        t_289[k] = f_5 * mg0_152[k]
                   - f_6 * mg1_152[k]
                   + pb_z[k] * mh_215[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pb_x, pb_z, lh_220, mg0_153, mg0_160, mg1_153, \
                         mg1_160, mh_216, mh_217, mh_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_26 * lh_220[k]
                   + f_3 * mg0_160[k]
                   - f_4 * mg1_160[k]
                   + pb_x[k] * mh_220[k];

        t_291[k] = pb_z[k] * mh_216[k];

        t_292[k] = f_3 * mg0_153[k]
                   - f_4 * mg1_153[k]
                   + pb_z[k] * mh_217[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pb_x, pb_y, pb_z, lh_135, lh_225, \
                         mg0_155, mg1_155, mh_219, mh_220, mh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_12 * lh_135[k]
                   + pb_y[k] * mh_219[k];

        t_294[k] = f_7 * mg0_155[k]
                   - f_8 * mg1_155[k]
                   + pb_z[k] * mh_219[k];

        t_295[k] = f_26 * lh_225[k]
                   + pb_x[k] * mh_225[k];

        t_296[k] = pb_z[k] * mh_220[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, pb_x, lh_227, lh_228, lh_229, lh_230, \
                         mh_227, mh_228, mh_229, mh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_26 * lh_227[k]
                   + pb_x[k] * mh_227[k];

        t_298[k] = f_26 * lh_228[k]
                   + pb_x[k] * mh_228[k];

        t_299[k] = f_26 * lh_229[k]
                   + pb_x[k] * mh_229[k];

        t_300[k] = f_26 * lh_230[k]
                   + pb_x[k] * mh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pa_x, pb_z, ki0_301, ki1_301, li_301, \
                         mg0_160, mg0_161, mg1_160, mg1_161, mh_225, mh_226, \
                         mh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_27 * ki0_301[k]
                   - f_28 * ki1_301[k]
                   + pa_x[k] * li_301[k];

        t_302[k] = pb_z[k] * mh_225[k];

        t_303[k] = f_3 * mg0_160[k]
                   - f_4 * mg1_160[k]
                   + pb_z[k] * mh_226[k];

        t_304[k] = f_5 * mg0_161[k]
                   - f_6 * mg1_161[k]
                   + pb_z[k] * mh_227[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_z, pb_y, pb_z, lh_146, li_168, \
                         mg0_162, mg0_164, mg1_162, mg1_164, mh_228, \
                         mh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_7 * mg0_162[k]
                   - f_8 * mg1_162[k]
                   + pb_z[k] * mh_228[k];

        t_306[k] = f_12 * lh_146[k]
                   + pb_y[k] * mh_230[k];

        t_307[k] = f_1 * mg0_164[k]
                   - f_2 * mg1_164[k]
                   + pb_z[k] * mh_230[k];

        t_308[k] = pa_z[k] * li_168[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pa_z, pb_y, pb_z, lh_126, lh_128, \
                         lh_149, li_169, li_171, li_173, mh_231, \
                         mh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = pa_z[k] * li_169[k];

        t_310[k] = f_9 * lh_126[k]
                   + pb_z[k] * mh_231[k];

        t_311[k] = pa_z[k] * li_171[k];

        t_312[k] = f_11 * lh_149[k]
                   + pb_y[k] * mh_233[k];

        t_313[k] = f_10 * lh_128[k]
                   + pa_z[k] * li_173[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, t_318, pa_z, pb_y, pb_z, lh_129, lh_131, \
                         lh_152, li_174, li_177, li_178, mh_234, \
                         mh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = pa_z[k] * li_174[k];

        t_315[k] = f_9 * lh_129[k]
                   + pb_z[k] * mh_234[k];

        t_316[k] = f_11 * lh_152[k]
                   + pb_y[k] * mh_236[k];

        t_317[k] = f_11 * lh_131[k]
                   + pa_z[k] * li_177[k];

        t_318[k] = pa_z[k] * li_178[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, pa_z, pb_y, pb_z, lh_132, lh_133, lh_135, \
                         lh_156, li_180, li_182, mh_237, mh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_9 * lh_132[k]
                   + pb_z[k] * mh_237[k];

        t_320[k] = f_10 * lh_133[k]
                   + pa_z[k] * li_180[k];

        t_321[k] = f_11 * lh_156[k]
                   + pb_y[k] * mh_240[k];

        t_322[k] = f_12 * lh_135[k]
                   + pa_z[k] * li_182[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, pa_z, pb_x, lh_247, lh_248, \
                         lh_249, lh_250, li_183, mh_247, mh_248, mh_249, \
                         mh_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = pa_z[k] * li_183[k];

        t_324[k] = f_26 * lh_247[k]
                   + pb_x[k] * mh_247[k];

        t_325[k] = f_26 * lh_248[k]
                   + pb_x[k] * mh_248[k];

        t_326[k] = f_26 * lh_249[k]
                   + pb_x[k] * mh_249[k];

        t_327[k] = f_26 * lh_250[k]
                   + pb_x[k] * mh_250[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pa_z, pb_x, pb_z, lh_141, lh_142, lh_251, \
                         li_189, li_191, mh_246, mh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_26 * lh_251[k]
                   + pb_x[k] * mh_251[k];

        t_329[k] = pa_z[k] * li_189[k];

        t_330[k] = f_9 * lh_141[k]
                   + pb_z[k] * mh_246[k];

        t_331[k] = f_10 * lh_142[k]
                   + pa_z[k] * li_191[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pa_z, pb_y, lh_143, lh_144, lh_146, \
                         lh_167, li_192, li_193, li_195, mh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_11 * lh_143[k]
                   + pa_z[k] * li_192[k];

        t_333[k] = f_12 * lh_144[k]
                   + pa_z[k] * li_193[k];

        t_334[k] = f_11 * lh_167[k]
                   + pb_y[k] * mh_251[k];

        t_335[k] = f_14 * lh_146[k]
                   + pa_z[k] * li_195[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pa_y, pb_y, pb_z, ki0_140, ki1_140, lh_147, \
                         lh_168, li_224, mh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_15 * ki0_140[k]
                   - f_16 * ki1_140[k]
                   + pa_y[k] * li_224[k];

        t_337[k] = f_10 * lh_168[k]
                   + pb_y[k] * mh_252[k];

        t_338[k] = f_10 * lh_147[k]
                   + pb_z[k] * mh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_y, pa_z, pb_y, ki0_87, ki0_145, ki1_87, \
                         ki1_145, lh_170, li_199, li_229, mh_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_15 * ki0_87[k]
                   - f_16 * ki1_87[k]
                   + pa_z[k] * li_199[k];

        t_340[k] = f_10 * lh_170[k]
                   + pb_y[k] * mh_254[k];

        t_341[k] = f_15 * ki0_145[k]
                   - f_16 * ki1_145[k]
                   + pa_y[k] * li_229[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pa_z, pb_y, pb_z, ki0_90, ki1_90, lh_150, \
                         lh_173, li_202, mh_255, mh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_15 * ki0_90[k]
                   - f_16 * ki1_90[k]
                   + pa_z[k] * li_202[k];

        t_343[k] = f_10 * lh_150[k]
                   + pb_z[k] * mh_255[k];

        t_344[k] = f_10 * lh_173[k]
                   + pb_y[k] * mh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pa_y, pa_z, pb_z, ki0_94, ki0_149, ki1_94, \
                         ki1_149, lh_153, li_206, li_233, mh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_15 * ki0_149[k]
                   - f_16 * ki1_149[k]
                   + pa_y[k] * li_233[k];

        t_346[k] = f_15 * ki0_94[k]
                   - f_16 * ki1_94[k]
                   + pa_z[k] * li_206[k];

        t_347[k] = f_10 * lh_153[k]
                   + pb_z[k] * mh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pa_y, pb_x, pb_y, ki0_154, ki1_154, lh_177, \
                         lh_264, li_238, mg0_192, mg1_192, mh_261, \
                         mh_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_26 * lh_264[k]
                   + f_3 * mg0_192[k]
                   - f_4 * mg1_192[k]
                   + pb_x[k] * mh_264[k];

        t_349[k] = f_10 * lh_177[k]
                   + pb_y[k] * mh_261[k];

        t_350[k] = f_15 * ki0_154[k]
                   - f_16 * ki1_154[k]
                   + pa_y[k] * li_238[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pb_x, lh_267, lh_268, lh_269, \
                         lh_270, lh_271, mh_267, mh_268, mh_269, mh_270, \
                         mh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_26 * lh_267[k]
                   + pb_x[k] * mh_267[k];

        t_352[k] = f_26 * lh_268[k]
                   + pb_x[k] * mh_268[k];

        t_353[k] = f_26 * lh_269[k]
                   + pb_x[k] * mh_269[k];

        t_354[k] = f_26 * lh_270[k]
                   + pb_x[k] * mh_270[k];

        t_355[k] = f_26 * lh_271[k]
                   + pb_x[k] * mh_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pa_x, pb_x, pb_z, ki0_357, ki1_357, lh_162, \
                         lh_272, li_357, mh_267, mh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_26 * lh_272[k]
                   + pb_x[k] * mh_272[k];

        t_357[k] = f_27 * ki0_357[k]
                   - f_28 * ki1_357[k]
                   + pa_x[k] * li_357[k];

        t_358[k] = f_10 * lh_162[k]
                   + pb_z[k] * mh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pa_x, ki0_359, ki0_360, ki0_361, ki1_359, \
                         ki1_360, ki1_361, li_359, li_360, li_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_27 * ki0_359[k]
                   - f_28 * ki1_359[k]
                   + pa_x[k] * li_359[k];

        t_360[k] = f_27 * ki0_360[k]
                   - f_28 * ki1_360[k]
                   + pa_x[k] * li_360[k];

        t_361[k] = f_27 * ki0_361[k]
                   - f_28 * ki1_361[k]
                   + pa_x[k] * li_361[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pa_x, pa_y, pb_y, ki0_363, ki1_363, \
                         lh_188, lh_189, li_252, li_363, mh_272, \
                         mh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_10 * lh_188[k]
                   + pb_y[k] * mh_272[k];

        t_363[k] = f_27 * ki0_363[k]
                   - f_28 * ki1_363[k]
                   + pa_x[k] * li_363[k];

        t_364[k] = pa_y[k] * li_252[k];

        t_365[k] = f_9 * lh_189[k]
                   + pb_y[k] * mh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, t_370, pa_y, pb_y, lh_190, lh_191, \
                         lh_192, li_254, li_255, li_257, li_258, \
                         mh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = pa_y[k] * li_254[k];

        t_367[k] = f_10 * lh_190[k]
                   + pa_y[k] * li_255[k];

        t_368[k] = f_9 * lh_191[k]
                   + pb_y[k] * mh_275[k];

        t_369[k] = pa_y[k] * li_257[k];

        t_370[k] = f_11 * lh_192[k]
                   + pa_y[k] * li_258[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pa_y, pb_y, pb_z, lh_171, lh_194, lh_195, \
                         li_261, li_262, mh_276, mh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_11 * lh_171[k]
                   + pb_z[k] * mh_276[k];

        t_372[k] = f_9 * lh_194[k]
                   + pb_y[k] * mh_278[k];

        t_373[k] = pa_y[k] * li_261[k];

        t_374[k] = f_12 * lh_195[k]
                   + pa_y[k] * li_262[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, pa_y, pb_y, pb_z, lh_174, lh_197, lh_198, \
                         li_264, li_266, mh_279, mh_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_11 * lh_174[k]
                   + pb_z[k] * mh_279[k];

        t_376[k] = f_10 * lh_197[k]
                   + pa_y[k] * li_264[k];

        t_377[k] = f_9 * lh_198[k]
                   + pb_y[k] * mh_282[k];

        t_378[k] = pa_y[k] * li_266[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, t_383, pb_x, lh_288, lh_289, lh_290, \
                         lh_291, lh_292, mh_288, mh_289, mh_290, mh_291, \
                         mh_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_26 * lh_288[k]
                   + pb_x[k] * mh_288[k];

        t_380[k] = f_26 * lh_289[k]
                   + pb_x[k] * mh_289[k];

        t_381[k] = f_26 * lh_290[k]
                   + pb_x[k] * mh_290[k];

        t_382[k] = f_26 * lh_291[k]
                   + pb_x[k] * mh_291[k];

        t_383[k] = f_26 * lh_292[k]
                   + pb_x[k] * mh_292[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, t_388, pa_y, pb_z, lh_183, lh_204, \
                         lh_206, lh_207, li_272, li_273, li_275, li_276, \
                         mh_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = pa_y[k] * li_272[k];

        t_385[k] = f_14 * lh_204[k]
                   + pa_y[k] * li_273[k];

        t_386[k] = f_11 * lh_183[k]
                   + pb_z[k] * mh_288[k];

        t_387[k] = f_12 * lh_206[k]
                   + pa_y[k] * li_275[k];

        t_388[k] = f_11 * lh_207[k]
                   + pa_y[k] * li_276[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, pa_y, pa_z, pb_y, ki0_140, ki1_140, \
                         lh_208, lh_209, li_252, li_277, li_279, \
                         mh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_10 * lh_208[k]
                   + pa_y[k] * li_277[k];

        t_390[k] = f_9 * lh_209[k]
                   + pb_y[k] * mh_293[k];

        t_391[k] = pa_y[k] * li_279[k];

        t_392[k] = f_24 * ki0_140[k]
                   - f_25 * ki1_140[k]
                   + pa_z[k] * li_252[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, pb_y, pb_z, lh_189, mg0_210, mg1_210, \
                         mh_294, mh_295, mh_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = pb_y[k] * mh_294[k];

        t_394[k] = f_12 * lh_189[k]
                   + pb_z[k] * mh_294[k];

        t_395[k] = f_3 * mg0_210[k]
                   - f_4 * mg1_210[k]
                   + pb_y[k] * mh_295[k];

        t_396[k] = pb_y[k] * mh_296[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pb_x, pb_y, pb_z, lh_192, lh_299, \
                         mg0_211, mg0_215, mg1_211, mg1_215, mh_297, \
                         mh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_26 * lh_299[k]
                   + f_7 * mg0_215[k]
                   - f_8 * mg1_215[k]
                   + pb_x[k] * mh_299[k];

        t_398[k] = f_5 * mg0_211[k]
                   - f_6 * mg1_211[k]
                   + pb_y[k] * mh_297[k];

        t_399[k] = f_12 * lh_192[k]
                   + pb_z[k] * mh_297[k];

        t_400[k] = pb_y[k] * mh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pb_x, pb_y, pb_z, lh_195, lh_303, mg0_213, \
                         mg0_219, mg1_213, mg1_219, mh_300, mh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_26 * lh_303[k]
                   + f_5 * mg0_219[k]
                   - f_6 * mg1_219[k]
                   + pb_x[k] * mh_303[k];

        t_402[k] = f_7 * mg0_213[k]
                   - f_8 * mg1_213[k]
                   + pb_y[k] * mh_300[k];

        t_403[k] = f_12 * lh_195[k]
                   + pb_z[k] * mh_300[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pb_x, pb_y, lh_308, lh_309, mg0_215, \
                         mg0_224, mg1_215, mg1_224, mh_302, mh_303, mh_308, \
                         mh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_3 * mg0_215[k]
                   - f_4 * mg1_215[k]
                   + pb_y[k] * mh_302[k];

        t_405[k] = pb_y[k] * mh_303[k];

        t_406[k] = f_26 * lh_308[k]
                   + f_3 * mg0_224[k]
                   - f_4 * mg1_224[k]
                   + pb_x[k] * mh_308[k];

        t_407[k] = f_26 * lh_309[k]
                   + pb_x[k] * mh_309[k];
    }
}

static auto
compute_prim_mi_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t ki0,
                                            const size_t ki1, const size_t lh, const size_t li,
                                            const size_t mg0, const size_t mg1, const size_t mh,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_14 = 3.0 / p;
    const auto f_15 = 0.5 / alpha;
    const auto f_16 = 0.5 * beta / (alpha * p);
    const auto f_20 = 1.0 / alpha;
    const auto f_21 = beta / (alpha * p);
    const auto f_24 = 1.5 / alpha;
    const auto f_25 = 1.5 * beta / (alpha * p);
    const auto f_26 = 2.5 / p;
    const auto f_27 = 2.0 / alpha;
    const auto f_28 = 2.0 * beta / (alpha * p);

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

    const auto *ki0_168 = buffer.data(ki0 + 168);
    const auto *ki0_171 = buffer.data(ki0 + 171);
    const auto *ki0_174 = buffer.data(ki0 + 174);
    const auto *ki0_178 = buffer.data(ki0 + 178);
    const auto *ki0_199 = buffer.data(ki0 + 199);
    const auto *ki0_202 = buffer.data(ki0 + 202);
    const auto *ki0_206 = buffer.data(ki0 + 206);
    const auto *ki0_224 = buffer.data(ki0 + 224);
    const auto *ki0_229 = buffer.data(ki0 + 229);
    const auto *ki0_233 = buffer.data(ki0 + 233);
    const auto *ki0_238 = buffer.data(ki0 + 238);
    const auto *ki0_252 = buffer.data(ki0 + 252);
    const auto *ki0_257 = buffer.data(ki0 + 257);
    const auto *ki0_261 = buffer.data(ki0 + 261);
    const auto *ki0_266 = buffer.data(ki0 + 266);
    const auto *ki0_419 = buffer.data(ki0 + 419);
    const auto *ki0_441 = buffer.data(ki0 + 441);
    const auto *ki0_497 = buffer.data(ki0 + 497);
    const auto *ki0_499 = buffer.data(ki0 + 499);
    const auto *ki0_500 = buffer.data(ki0 + 500);
    const auto *ki0_501 = buffer.data(ki0 + 501);
    const auto *ki0_503 = buffer.data(ki0 + 503);
    const auto *ki0_525 = buffer.data(ki0 + 525);
    const auto *ki0_527 = buffer.data(ki0 + 527);
    const auto *ki0_528 = buffer.data(ki0 + 528);
    const auto *ki0_529 = buffer.data(ki0 + 529);
    const auto *ki0_531 = buffer.data(ki0 + 531);

    const auto *ki1_168 = buffer.data(ki1 + 168);
    const auto *ki1_171 = buffer.data(ki1 + 171);
    const auto *ki1_174 = buffer.data(ki1 + 174);
    const auto *ki1_178 = buffer.data(ki1 + 178);
    const auto *ki1_199 = buffer.data(ki1 + 199);
    const auto *ki1_202 = buffer.data(ki1 + 202);
    const auto *ki1_206 = buffer.data(ki1 + 206);
    const auto *ki1_224 = buffer.data(ki1 + 224);
    const auto *ki1_229 = buffer.data(ki1 + 229);
    const auto *ki1_233 = buffer.data(ki1 + 233);
    const auto *ki1_238 = buffer.data(ki1 + 238);
    const auto *ki1_252 = buffer.data(ki1 + 252);
    const auto *ki1_257 = buffer.data(ki1 + 257);
    const auto *ki1_261 = buffer.data(ki1 + 261);
    const auto *ki1_266 = buffer.data(ki1 + 266);
    const auto *ki1_419 = buffer.data(ki1 + 419);
    const auto *ki1_441 = buffer.data(ki1 + 441);
    const auto *ki1_497 = buffer.data(ki1 + 497);
    const auto *ki1_499 = buffer.data(ki1 + 499);
    const auto *ki1_500 = buffer.data(ki1 + 500);
    const auto *ki1_501 = buffer.data(ki1 + 501);
    const auto *ki1_503 = buffer.data(ki1 + 503);
    const auto *ki1_525 = buffer.data(ki1 + 525);
    const auto *ki1_527 = buffer.data(ki1 + 527);
    const auto *ki1_528 = buffer.data(ki1 + 528);
    const auto *ki1_529 = buffer.data(ki1 + 529);
    const auto *ki1_531 = buffer.data(ki1 + 531);

    const auto *lh_204 = buffer.data(lh + 204);
    const auto *lh_210 = buffer.data(lh + 210);
    const auto *lh_212 = buffer.data(lh + 212);
    const auto *lh_213 = buffer.data(lh + 213);
    const auto *lh_215 = buffer.data(lh + 215);
    const auto *lh_216 = buffer.data(lh + 216);
    const auto *lh_217 = buffer.data(lh + 217);
    const auto *lh_219 = buffer.data(lh + 219);
    const auto *lh_225 = buffer.data(lh + 225);
    const auto *lh_226 = buffer.data(lh + 226);
    const auto *lh_227 = buffer.data(lh + 227);
    const auto *lh_228 = buffer.data(lh + 228);
    const auto *lh_230 = buffer.data(lh + 230);
    const auto *lh_231 = buffer.data(lh + 231);
    const auto *lh_233 = buffer.data(lh + 233);
    const auto *lh_234 = buffer.data(lh + 234);
    const auto *lh_236 = buffer.data(lh + 236);
    const auto *lh_237 = buffer.data(lh + 237);
    const auto *lh_240 = buffer.data(lh + 240);
    const auto *lh_246 = buffer.data(lh + 246);
    const auto *lh_251 = buffer.data(lh + 251);
    const auto *lh_252 = buffer.data(lh + 252);
    const auto *lh_254 = buffer.data(lh + 254);
    const auto *lh_255 = buffer.data(lh + 255);
    const auto *lh_257 = buffer.data(lh + 257);
    const auto *lh_258 = buffer.data(lh + 258);
    const auto *lh_261 = buffer.data(lh + 261);
    const auto *lh_267 = buffer.data(lh + 267);
    const auto *lh_272 = buffer.data(lh + 272);
    const auto *lh_273 = buffer.data(lh + 273);
    const auto *lh_275 = buffer.data(lh + 275);
    const auto *lh_278 = buffer.data(lh + 278);
    const auto *lh_282 = buffer.data(lh + 282);
    const auto *lh_293 = buffer.data(lh + 293);
    const auto *lh_310 = buffer.data(lh + 310);
    const auto *lh_311 = buffer.data(lh + 311);
    const auto *lh_312 = buffer.data(lh + 312);
    const auto *lh_314 = buffer.data(lh + 314);
    const auto *lh_318 = buffer.data(lh + 318);
    const auto *lh_321 = buffer.data(lh + 321);
    const auto *lh_325 = buffer.data(lh + 325);
    const auto *lh_330 = buffer.data(lh + 330);
    const auto *lh_332 = buffer.data(lh + 332);
    const auto *lh_333 = buffer.data(lh + 333);
    const auto *lh_334 = buffer.data(lh + 334);
    const auto *lh_335 = buffer.data(lh + 335);
    const auto *lh_352 = buffer.data(lh + 352);
    const auto *lh_353 = buffer.data(lh + 353);
    const auto *lh_354 = buffer.data(lh + 354);
    const auto *lh_355 = buffer.data(lh + 355);
    const auto *lh_356 = buffer.data(lh + 356);
    const auto *lh_369 = buffer.data(lh + 369);
    const auto *lh_372 = buffer.data(lh + 372);
    const auto *lh_373 = buffer.data(lh + 373);
    const auto *lh_374 = buffer.data(lh + 374);
    const auto *lh_375 = buffer.data(lh + 375);
    const auto *lh_376 = buffer.data(lh + 376);
    const auto *lh_377 = buffer.data(lh + 377);
    const auto *lh_390 = buffer.data(lh + 390);
    const auto *lh_393 = buffer.data(lh + 393);
    const auto *lh_394 = buffer.data(lh + 394);
    const auto *lh_395 = buffer.data(lh + 395);
    const auto *lh_396 = buffer.data(lh + 396);
    const auto *lh_397 = buffer.data(lh + 397);
    const auto *lh_398 = buffer.data(lh + 398);

    const auto *li_280 = buffer.data(li + 280);
    const auto *li_281 = buffer.data(li + 281);
    const auto *li_283 = buffer.data(li + 283);
    const auto *li_285 = buffer.data(li + 285);
    const auto *li_286 = buffer.data(li + 286);
    const auto *li_289 = buffer.data(li + 289);
    const auto *li_290 = buffer.data(li + 290);
    const auto *li_292 = buffer.data(li + 292);
    const auto *li_294 = buffer.data(li + 294);
    const auto *li_295 = buffer.data(li + 295);
    const auto *li_301 = buffer.data(li + 301);
    const auto *li_303 = buffer.data(li + 303);
    const auto *li_304 = buffer.data(li + 304);
    const auto *li_305 = buffer.data(li + 305);
    const auto *li_307 = buffer.data(li + 307);
    const auto *li_311 = buffer.data(li + 311);
    const auto *li_314 = buffer.data(li + 314);
    const auto *li_318 = buffer.data(li + 318);
    const auto *li_336 = buffer.data(li + 336);
    const auto *li_339 = buffer.data(li + 339);
    const auto *li_341 = buffer.data(li + 341);
    const auto *li_342 = buffer.data(li + 342);
    const auto *li_345 = buffer.data(li + 345);
    const auto *li_346 = buffer.data(li + 346);
    const auto *li_350 = buffer.data(li + 350);
    const auto *li_364 = buffer.data(li + 364);
    const auto *li_369 = buffer.data(li + 369);
    const auto *li_373 = buffer.data(li + 373);
    const auto *li_378 = buffer.data(li + 378);
    const auto *li_392 = buffer.data(li + 392);
    const auto *li_419 = buffer.data(li + 419);
    const auto *li_441 = buffer.data(li + 441);
    const auto *li_497 = buffer.data(li + 497);
    const auto *li_499 = buffer.data(li + 499);
    const auto *li_500 = buffer.data(li + 500);
    const auto *li_501 = buffer.data(li + 501);
    const auto *li_503 = buffer.data(li + 503);
    const auto *li_525 = buffer.data(li + 525);
    const auto *li_527 = buffer.data(li + 527);
    const auto *li_528 = buffer.data(li + 528);
    const auto *li_529 = buffer.data(li + 529);
    const auto *li_531 = buffer.data(li + 531);

    const auto *mg0_220 = buffer.data(mg0 + 220);
    const auto *mg0_222 = buffer.data(mg0 + 222);
    const auto *mg0_223 = buffer.data(mg0 + 223);
    const auto *mg0_224 = buffer.data(mg0 + 224);
    const auto *mg0_225 = buffer.data(mg0 + 225);
    const auto *mg0_227 = buffer.data(mg0 + 227);
    const auto *mg0_228 = buffer.data(mg0 + 228);
    const auto *mg0_230 = buffer.data(mg0 + 230);
    const auto *mg0_231 = buffer.data(mg0 + 231);
    const auto *mg0_235 = buffer.data(mg0 + 235);
    const auto *mg0_236 = buffer.data(mg0 + 236);
    const auto *mg0_237 = buffer.data(mg0 + 237);
    const auto *mg0_239 = buffer.data(mg0 + 239);
    const auto *mg0_267 = buffer.data(mg0 + 267);
    const auto *mg0_282 = buffer.data(mg0 + 282);

    const auto *mg1_220 = buffer.data(mg1 + 220);
    const auto *mg1_222 = buffer.data(mg1 + 222);
    const auto *mg1_223 = buffer.data(mg1 + 223);
    const auto *mg1_224 = buffer.data(mg1 + 224);
    const auto *mg1_225 = buffer.data(mg1 + 225);
    const auto *mg1_227 = buffer.data(mg1 + 227);
    const auto *mg1_228 = buffer.data(mg1 + 228);
    const auto *mg1_230 = buffer.data(mg1 + 230);
    const auto *mg1_231 = buffer.data(mg1 + 231);
    const auto *mg1_235 = buffer.data(mg1 + 235);
    const auto *mg1_236 = buffer.data(mg1 + 236);
    const auto *mg1_237 = buffer.data(mg1 + 237);
    const auto *mg1_239 = buffer.data(mg1 + 239);
    const auto *mg1_267 = buffer.data(mg1 + 267);
    const auto *mg1_282 = buffer.data(mg1 + 282);

    const auto *mh_308 = buffer.data(mh + 308);
    const auto *mh_309 = buffer.data(mh + 309);
    const auto *mh_310 = buffer.data(mh + 310);
    const auto *mh_311 = buffer.data(mh + 311);
    const auto *mh_312 = buffer.data(mh + 312);
    const auto *mh_313 = buffer.data(mh + 313);
    const auto *mh_314 = buffer.data(mh + 314);
    const auto *mh_315 = buffer.data(mh + 315);
    const auto *mh_316 = buffer.data(mh + 316);
    const auto *mh_317 = buffer.data(mh + 317);
    const auto *mh_318 = buffer.data(mh + 318);
    const auto *mh_320 = buffer.data(mh + 320);
    const auto *mh_321 = buffer.data(mh + 321);
    const auto *mh_322 = buffer.data(mh + 322);
    const auto *mh_324 = buffer.data(mh + 324);
    const auto *mh_325 = buffer.data(mh + 325);
    const auto *mh_330 = buffer.data(mh + 330);
    const auto *mh_331 = buffer.data(mh + 331);
    const auto *mh_332 = buffer.data(mh + 332);
    const auto *mh_333 = buffer.data(mh + 333);
    const auto *mh_334 = buffer.data(mh + 334);
    const auto *mh_335 = buffer.data(mh + 335);
    const auto *mh_336 = buffer.data(mh + 336);
    const auto *mh_338 = buffer.data(mh + 338);
    const auto *mh_339 = buffer.data(mh + 339);
    const auto *mh_341 = buffer.data(mh + 341);
    const auto *mh_342 = buffer.data(mh + 342);
    const auto *mh_345 = buffer.data(mh + 345);
    const auto *mh_351 = buffer.data(mh + 351);
    const auto *mh_352 = buffer.data(mh + 352);
    const auto *mh_353 = buffer.data(mh + 353);
    const auto *mh_354 = buffer.data(mh + 354);
    const auto *mh_355 = buffer.data(mh + 355);
    const auto *mh_356 = buffer.data(mh + 356);
    const auto *mh_357 = buffer.data(mh + 357);
    const auto *mh_359 = buffer.data(mh + 359);
    const auto *mh_360 = buffer.data(mh + 360);
    const auto *mh_362 = buffer.data(mh + 362);
    const auto *mh_363 = buffer.data(mh + 363);
    const auto *mh_366 = buffer.data(mh + 366);
    const auto *mh_369 = buffer.data(mh + 369);
    const auto *mh_372 = buffer.data(mh + 372);
    const auto *mh_373 = buffer.data(mh + 373);
    const auto *mh_374 = buffer.data(mh + 374);
    const auto *mh_375 = buffer.data(mh + 375);
    const auto *mh_376 = buffer.data(mh + 376);
    const auto *mh_377 = buffer.data(mh + 377);
    const auto *mh_378 = buffer.data(mh + 378);
    const auto *mh_380 = buffer.data(mh + 380);
    const auto *mh_381 = buffer.data(mh + 381);
    const auto *mh_383 = buffer.data(mh + 383);
    const auto *mh_384 = buffer.data(mh + 384);
    const auto *mh_387 = buffer.data(mh + 387);
    const auto *mh_390 = buffer.data(mh + 390);
    const auto *mh_393 = buffer.data(mh + 393);
    const auto *mh_394 = buffer.data(mh + 394);
    const auto *mh_395 = buffer.data(mh + 395);
    const auto *mh_396 = buffer.data(mh + 396);
    const auto *mh_397 = buffer.data(mh + 397);
    const auto *mh_398 = buffer.data(mh + 398);

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, pb_x, pb_y, lh_310, lh_311, \
                         lh_312, lh_314, mh_308, mh_310, mh_311, mh_312, \
                         mh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_26 * lh_310[k]
                   + pb_x[k] * mh_310[k];

        t_409[k] = f_26 * lh_311[k]
                   + pb_x[k] * mh_311[k];

        t_410[k] = f_26 * lh_312[k]
                   + pb_x[k] * mh_312[k];

        t_411[k] = pb_y[k] * mh_308[k];

        t_412[k] = f_26 * lh_314[k]
                   + pb_x[k] * mh_314[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, pb_y, pb_z, lh_204, mg0_220, mg0_222, \
                         mg0_223, mg1_220, mg1_222, mg1_223, mh_309, mh_311, \
                         mh_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_1 * mg0_220[k]
                   - f_2 * mg1_220[k]
                   + pb_y[k] * mh_309[k];

        t_414[k] = f_12 * lh_204[k]
                   + pb_z[k] * mh_309[k];

        t_415[k] = f_7 * mg0_222[k]
                   - f_8 * mg1_222[k]
                   + pb_y[k] * mh_311[k];

        t_416[k] = f_5 * mg0_223[k]
                   - f_6 * mg1_223[k]
                   + pb_y[k] * mh_312[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pa_x, pb_y, ki0_419, ki1_419, li_419, mg0_224, \
                         mg1_224, mh_313, mh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_3 * mg0_224[k]
                   - f_4 * mg1_224[k]
                   + pb_y[k] * mh_313[k];

        t_418[k] = pb_y[k] * mh_314[k];

        t_419[k] = f_27 * ki0_419[k]
                   - f_28 * ki1_419[k]
                   + pa_x[k] * li_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pa_y, pb_y, pb_z, ki0_168, ki1_168, lh_210, \
                         li_280, mh_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_27 * ki0_168[k]
                   - f_28 * ki1_168[k]
                   + pa_y[k] * li_280[k];

        t_421[k] = f_26 * lh_210[k]
                   + pb_y[k] * mh_315[k];

        t_422[k] = pb_z[k] * mh_315[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pb_x, pb_z, lh_318, mg0_225, mg0_228, mg1_225, \
                         mg1_228, mh_316, mh_317, mh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_12 * lh_318[k]
                   + f_7 * mg0_228[k]
                   - f_8 * mg1_228[k]
                   + pb_x[k] * mh_318[k];

        t_424[k] = pb_z[k] * mh_316[k];

        t_425[k] = f_3 * mg0_225[k]
                   - f_4 * mg1_225[k]
                   + pb_z[k] * mh_317[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, pb_x, pb_y, pb_z, lh_215, lh_321, \
                         mg0_227, mg0_231, mg1_227, mg1_231, mh_318, mh_320, \
                         mh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_12 * lh_321[k]
                   + f_5 * mg0_231[k]
                   - f_6 * mg1_231[k]
                   + pb_x[k] * mh_321[k];

        t_427[k] = pb_z[k] * mh_318[k];

        t_428[k] = f_26 * lh_215[k]
                   + pb_y[k] * mh_320[k];

        t_429[k] = f_5 * mg0_227[k]
                   - f_6 * mg1_227[k]
                   + pb_z[k] * mh_320[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pb_x, pb_z, lh_325, mg0_228, mg0_235, mg1_228, \
                         mg1_235, mh_321, mh_322, mh_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_12 * lh_325[k]
                   + f_3 * mg0_235[k]
                   - f_4 * mg1_235[k]
                   + pb_x[k] * mh_325[k];

        t_431[k] = pb_z[k] * mh_321[k];

        t_432[k] = f_3 * mg0_228[k]
                   - f_4 * mg1_228[k]
                   + pb_z[k] * mh_322[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pb_x, pb_y, pb_z, lh_219, lh_330, \
                         mg0_230, mg1_230, mh_324, mh_325, mh_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_26 * lh_219[k]
                   + pb_y[k] * mh_324[k];

        t_434[k] = f_7 * mg0_230[k]
                   - f_8 * mg1_230[k]
                   + pb_z[k] * mh_324[k];

        t_435[k] = f_12 * lh_330[k]
                   + pb_x[k] * mh_330[k];

        t_436[k] = pb_z[k] * mh_325[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pb_x, lh_332, lh_333, lh_334, lh_335, \
                         mh_332, mh_333, mh_334, mh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_12 * lh_332[k]
                   + pb_x[k] * mh_332[k];

        t_438[k] = f_12 * lh_333[k]
                   + pb_x[k] * mh_333[k];

        t_439[k] = f_12 * lh_334[k]
                   + pb_x[k] * mh_334[k];

        t_440[k] = f_12 * lh_335[k]
                   + pb_x[k] * mh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pa_x, pb_z, ki0_441, ki1_441, li_441, \
                         mg0_235, mg0_236, mg1_235, mg1_236, mh_330, mh_331, \
                         mh_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_24 * ki0_441[k]
                   - f_25 * ki1_441[k]
                   + pa_x[k] * li_441[k];

        t_442[k] = pb_z[k] * mh_330[k];

        t_443[k] = f_3 * mg0_235[k]
                   - f_4 * mg1_235[k]
                   + pb_z[k] * mh_331[k];

        t_444[k] = f_5 * mg0_236[k]
                   - f_6 * mg1_236[k]
                   + pb_z[k] * mh_332[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pa_z, pb_y, pb_z, lh_230, li_280, \
                         mg0_237, mg0_239, mg1_237, mg1_239, mh_333, \
                         mh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_7 * mg0_237[k]
                   - f_8 * mg1_237[k]
                   + pb_z[k] * mh_333[k];

        t_446[k] = f_26 * lh_230[k]
                   + pb_y[k] * mh_335[k];

        t_447[k] = f_1 * mg0_239[k]
                   - f_2 * mg1_239[k]
                   + pb_z[k] * mh_335[k];

        t_448[k] = pa_z[k] * li_280[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, t_453, pa_z, pb_y, pb_z, lh_210, lh_212, \
                         lh_233, li_281, li_283, li_285, mh_336, \
                         mh_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = pa_z[k] * li_281[k];

        t_450[k] = f_9 * lh_210[k]
                   + pb_z[k] * mh_336[k];

        t_451[k] = pa_z[k] * li_283[k];

        t_452[k] = f_12 * lh_233[k]
                   + pb_y[k] * mh_338[k];

        t_453[k] = f_10 * lh_212[k]
                   + pa_z[k] * li_285[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, pa_z, pb_y, pb_z, lh_213, lh_215, \
                         lh_236, li_286, li_289, li_290, mh_339, \
                         mh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = pa_z[k] * li_286[k];

        t_455[k] = f_9 * lh_213[k]
                   + pb_z[k] * mh_339[k];

        t_456[k] = f_12 * lh_236[k]
                   + pb_y[k] * mh_341[k];

        t_457[k] = f_11 * lh_215[k]
                   + pa_z[k] * li_289[k];

        t_458[k] = pa_z[k] * li_290[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, t_462, pa_z, pb_y, pb_z, lh_216, lh_217, lh_219, \
                         lh_240, li_292, li_294, mh_342, mh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_9 * lh_216[k]
                   + pb_z[k] * mh_342[k];

        t_460[k] = f_10 * lh_217[k]
                   + pa_z[k] * li_292[k];

        t_461[k] = f_12 * lh_240[k]
                   + pb_y[k] * mh_345[k];

        t_462[k] = f_12 * lh_219[k]
                   + pa_z[k] * li_294[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, t_467, pa_z, pb_x, lh_352, lh_353, \
                         lh_354, lh_355, li_295, mh_352, mh_353, mh_354, \
                         mh_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = pa_z[k] * li_295[k];

        t_464[k] = f_12 * lh_352[k]
                   + pb_x[k] * mh_352[k];

        t_465[k] = f_12 * lh_353[k]
                   + pb_x[k] * mh_353[k];

        t_466[k] = f_12 * lh_354[k]
                   + pb_x[k] * mh_354[k];

        t_467[k] = f_12 * lh_355[k]
                   + pb_x[k] * mh_355[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, t_471, pa_z, pb_x, pb_z, lh_225, lh_226, lh_356, \
                         li_301, li_303, mh_351, mh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_12 * lh_356[k]
                   + pb_x[k] * mh_356[k];

        t_469[k] = pa_z[k] * li_301[k];

        t_470[k] = f_9 * lh_225[k]
                   + pb_z[k] * mh_351[k];

        t_471[k] = f_10 * lh_226[k]
                   + pa_z[k] * li_303[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, pa_z, pb_y, lh_227, lh_228, lh_230, \
                         lh_251, li_304, li_305, li_307, mh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = f_11 * lh_227[k]
                   + pa_z[k] * li_304[k];

        t_473[k] = f_12 * lh_228[k]
                   + pa_z[k] * li_305[k];

        t_474[k] = f_12 * lh_251[k]
                   + pb_y[k] * mh_356[k];

        t_475[k] = f_14 * lh_230[k]
                   + pa_z[k] * li_307[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pa_y, pb_y, pb_z, ki0_224, ki1_224, lh_231, \
                         lh_252, li_336, mh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_20 * ki0_224[k]
                   - f_21 * ki1_224[k]
                   + pa_y[k] * li_336[k];

        t_477[k] = f_11 * lh_252[k]
                   + pb_y[k] * mh_357[k];

        t_478[k] = f_10 * lh_231[k]
                   + pb_z[k] * mh_357[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pa_y, pa_z, pb_y, ki0_171, ki0_229, ki1_171, \
                         ki1_229, lh_254, li_311, li_341, mh_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_15 * ki0_171[k]
                   - f_16 * ki1_171[k]
                   + pa_z[k] * li_311[k];

        t_480[k] = f_11 * lh_254[k]
                   + pb_y[k] * mh_359[k];

        t_481[k] = f_20 * ki0_229[k]
                   - f_21 * ki1_229[k]
                   + pa_y[k] * li_341[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pa_z, pb_y, pb_z, ki0_174, ki1_174, lh_234, \
                         lh_257, li_314, mh_360, mh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_15 * ki0_174[k]
                   - f_16 * ki1_174[k]
                   + pa_z[k] * li_314[k];

        t_483[k] = f_10 * lh_234[k]
                   + pb_z[k] * mh_360[k];

        t_484[k] = f_11 * lh_257[k]
                   + pb_y[k] * mh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pa_y, pa_z, pb_z, ki0_178, ki0_233, ki1_178, \
                         ki1_233, lh_237, li_318, li_345, mh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_20 * ki0_233[k]
                   - f_21 * ki1_233[k]
                   + pa_y[k] * li_345[k];

        t_486[k] = f_15 * ki0_178[k]
                   - f_16 * ki1_178[k]
                   + pa_z[k] * li_318[k];

        t_487[k] = f_10 * lh_237[k]
                   + pb_z[k] * mh_363[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pa_y, pb_x, pb_y, ki0_238, ki1_238, lh_261, \
                         lh_369, li_350, mg0_267, mg1_267, mh_366, \
                         mh_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_12 * lh_369[k]
                   + f_3 * mg0_267[k]
                   - f_4 * mg1_267[k]
                   + pb_x[k] * mh_369[k];

        t_489[k] = f_11 * lh_261[k]
                   + pb_y[k] * mh_366[k];

        t_490[k] = f_20 * ki0_238[k]
                   - f_21 * ki1_238[k]
                   + pa_y[k] * li_350[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pb_x, lh_372, lh_373, lh_374, \
                         lh_375, lh_376, mh_372, mh_373, mh_374, mh_375, \
                         mh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_12 * lh_372[k]
                   + pb_x[k] * mh_372[k];

        t_492[k] = f_12 * lh_373[k]
                   + pb_x[k] * mh_373[k];

        t_493[k] = f_12 * lh_374[k]
                   + pb_x[k] * mh_374[k];

        t_494[k] = f_12 * lh_375[k]
                   + pb_x[k] * mh_375[k];

        t_495[k] = f_12 * lh_376[k]
                   + pb_x[k] * mh_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pa_x, pb_x, pb_z, ki0_497, ki1_497, lh_246, \
                         lh_377, li_497, mh_372, mh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_12 * lh_377[k]
                   + pb_x[k] * mh_377[k];

        t_497[k] = f_24 * ki0_497[k]
                   - f_25 * ki1_497[k]
                   + pa_x[k] * li_497[k];

        t_498[k] = f_10 * lh_246[k]
                   + pb_z[k] * mh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pa_x, ki0_499, ki0_500, ki0_501, ki1_499, \
                         ki1_500, ki1_501, li_499, li_500, li_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_24 * ki0_499[k]
                   - f_25 * ki1_499[k]
                   + pa_x[k] * li_499[k];

        t_500[k] = f_24 * ki0_500[k]
                   - f_25 * ki1_500[k]
                   + pa_x[k] * li_500[k];

        t_501[k] = f_24 * ki0_501[k]
                   - f_25 * ki1_501[k]
                   + pa_x[k] * li_501[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pa_x, pa_y, pb_y, ki0_252, ki0_503, ki1_252, \
                         ki1_503, lh_272, li_364, li_503, mh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_11 * lh_272[k]
                   + pb_y[k] * mh_377[k];

        t_503[k] = f_24 * ki0_503[k]
                   - f_25 * ki1_503[k]
                   + pa_x[k] * li_503[k];

        t_504[k] = f_15 * ki0_252[k]
                   - f_16 * ki1_252[k]
                   + pa_y[k] * li_364[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pa_z, pb_y, pb_z, ki0_199, ki1_199, \
                         lh_252, lh_273, lh_275, li_339, mh_378, \
                         mh_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_10 * lh_273[k]
                   + pb_y[k] * mh_378[k];

        t_506[k] = f_11 * lh_252[k]
                   + pb_z[k] * mh_378[k];

        t_507[k] = f_20 * ki0_199[k]
                   - f_21 * ki1_199[k]
                   + pa_z[k] * li_339[k];

        t_508[k] = f_10 * lh_275[k]
                   + pb_y[k] * mh_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pa_y, pa_z, pb_z, ki0_202, ki0_257, ki1_202, \
                         ki1_257, lh_255, li_342, li_369, mh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_15 * ki0_257[k]
                   - f_16 * ki1_257[k]
                   + pa_y[k] * li_369[k];

        t_510[k] = f_20 * ki0_202[k]
                   - f_21 * ki1_202[k]
                   + pa_z[k] * li_342[k];

        t_511[k] = f_11 * lh_255[k]
                   + pb_z[k] * mh_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pa_y, pa_z, pb_y, ki0_206, ki0_261, ki1_206, \
                         ki1_261, lh_278, li_346, li_373, mh_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_10 * lh_278[k]
                   + pb_y[k] * mh_383[k];

        t_513[k] = f_15 * ki0_261[k]
                   - f_16 * ki1_261[k]
                   + pa_y[k] * li_373[k];

        t_514[k] = f_20 * ki0_206[k]
                   - f_21 * ki1_206[k]
                   + pa_z[k] * li_346[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pb_x, pb_y, pb_z, lh_258, lh_282, lh_390, \
                         mg0_282, mg1_282, mh_384, mh_387, mh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_11 * lh_258[k]
                   + pb_z[k] * mh_384[k];

        t_516[k] = f_12 * lh_390[k]
                   + f_3 * mg0_282[k]
                   - f_4 * mg1_282[k]
                   + pb_x[k] * mh_390[k];

        t_517[k] = f_10 * lh_282[k]
                   + pb_y[k] * mh_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pa_y, pb_x, ki0_266, ki1_266, lh_393, \
                         lh_394, lh_395, li_378, mh_393, mh_394, \
                         mh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_15 * ki0_266[k]
                   - f_16 * ki1_266[k]
                   + pa_y[k] * li_378[k];

        t_519[k] = f_12 * lh_393[k]
                   + pb_x[k] * mh_393[k];

        t_520[k] = f_12 * lh_394[k]
                   + pb_x[k] * mh_394[k];

        t_521[k] = f_12 * lh_395[k]
                   + pb_x[k] * mh_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pa_x, pb_x, ki0_525, ki1_525, lh_396, \
                         lh_397, lh_398, li_525, mh_396, mh_397, \
                         mh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_12 * lh_396[k]
                   + pb_x[k] * mh_396[k];

        t_523[k] = f_12 * lh_397[k]
                   + pb_x[k] * mh_397[k];

        t_524[k] = f_12 * lh_398[k]
                   + pb_x[k] * mh_398[k];

        t_525[k] = f_24 * ki0_525[k]
                   - f_25 * ki1_525[k]
                   + pa_x[k] * li_525[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pa_x, pb_z, ki0_527, ki0_528, ki1_527, ki1_528, \
                         lh_267, li_527, li_528, mh_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_11 * lh_267[k]
                   + pb_z[k] * mh_393[k];

        t_527[k] = f_24 * ki0_527[k]
                   - f_25 * ki1_527[k]
                   + pa_x[k] * li_527[k];

        t_528[k] = f_24 * ki0_528[k]
                   - f_25 * ki1_528[k]
                   + pa_x[k] * li_528[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pa_x, pa_y, pb_y, ki0_529, ki0_531, \
                         ki1_529, ki1_531, lh_293, li_392, li_529, li_531, \
                         mh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_24 * ki0_529[k]
                   - f_25 * ki1_529[k]
                   + pa_x[k] * li_529[k];

        t_530[k] = f_10 * lh_293[k]
                   + pb_y[k] * mh_398[k];

        t_531[k] = f_24 * ki0_531[k]
                   - f_25 * ki1_531[k]
                   + pa_x[k] * li_531[k];

        t_532[k] = pa_y[k] * li_392[k];
    }
}

static auto
compute_prim_mi_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t ki0,
                                            const size_t ki1, const size_t lh, const size_t li,
                                            const size_t mg0, const size_t mg1, const size_t mh,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_14 = 3.0 / p;
    const auto f_15 = 0.5 / alpha;
    const auto f_16 = 0.5 * beta / (alpha * p);
    const auto f_20 = 1.0 / alpha;
    const auto f_21 = beta / (alpha * p);
    const auto f_22 = 2.5 / alpha;
    const auto f_23 = 2.5 * beta / (alpha * p);
    const auto f_24 = 1.5 / alpha;
    const auto f_25 = 1.5 * beta / (alpha * p);
    const auto f_26 = 2.5 / p;
    const auto f_27 = 2.0 / alpha;
    const auto f_28 = 2.0 * beta / (alpha * p);

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
    auto *t_663 = buffer.data(target + 663);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ki0_252 = buffer.data(ki0 + 252);
    const auto *ki0_280 = buffer.data(ki0 + 280);
    const auto *ki0_283 = buffer.data(ki0 + 283);
    const auto *ki0_286 = buffer.data(ki0 + 286);
    const auto *ki0_290 = buffer.data(ki0 + 290);
    const auto *ki0_336 = buffer.data(ki0 + 336);
    const auto *ki0_341 = buffer.data(ki0 + 341);
    const auto *ki0_345 = buffer.data(ki0 + 345);
    const auto *ki0_350 = buffer.data(ki0 + 350);
    const auto *ki0_587 = buffer.data(ki0 + 587);
    const auto *ki0_609 = buffer.data(ki0 + 609);

    const auto *ki1_252 = buffer.data(ki1 + 252);
    const auto *ki1_280 = buffer.data(ki1 + 280);
    const auto *ki1_283 = buffer.data(ki1 + 283);
    const auto *ki1_286 = buffer.data(ki1 + 286);
    const auto *ki1_290 = buffer.data(ki1 + 290);
    const auto *ki1_336 = buffer.data(ki1 + 336);
    const auto *ki1_341 = buffer.data(ki1 + 341);
    const auto *ki1_345 = buffer.data(ki1 + 345);
    const auto *ki1_350 = buffer.data(ki1 + 350);
    const auto *ki1_587 = buffer.data(ki1 + 587);
    const auto *ki1_609 = buffer.data(ki1 + 609);

    const auto *lh_276 = buffer.data(lh + 276);
    const auto *lh_279 = buffer.data(lh + 279);
    const auto *lh_288 = buffer.data(lh + 288);
    const auto *lh_294 = buffer.data(lh + 294);
    const auto *lh_295 = buffer.data(lh + 295);
    const auto *lh_296 = buffer.data(lh + 296);
    const auto *lh_297 = buffer.data(lh + 297);
    const auto *lh_299 = buffer.data(lh + 299);
    const auto *lh_300 = buffer.data(lh + 300);
    const auto *lh_302 = buffer.data(lh + 302);
    const auto *lh_303 = buffer.data(lh + 303);
    const auto *lh_309 = buffer.data(lh + 309);
    const auto *lh_311 = buffer.data(lh + 311);
    const auto *lh_312 = buffer.data(lh + 312);
    const auto *lh_313 = buffer.data(lh + 313);
    const auto *lh_314 = buffer.data(lh + 314);
    const auto *lh_315 = buffer.data(lh + 315);
    const auto *lh_317 = buffer.data(lh + 317);
    const auto *lh_318 = buffer.data(lh + 318);
    const auto *lh_320 = buffer.data(lh + 320);
    const auto *lh_321 = buffer.data(lh + 321);
    const auto *lh_322 = buffer.data(lh + 322);
    const auto *lh_324 = buffer.data(lh + 324);
    const auto *lh_330 = buffer.data(lh + 330);
    const auto *lh_331 = buffer.data(lh + 331);
    const auto *lh_332 = buffer.data(lh + 332);
    const auto *lh_333 = buffer.data(lh + 333);
    const auto *lh_335 = buffer.data(lh + 335);
    const auto *lh_336 = buffer.data(lh + 336);
    const auto *lh_338 = buffer.data(lh + 338);
    const auto *lh_339 = buffer.data(lh + 339);
    const auto *lh_341 = buffer.data(lh + 341);
    const auto *lh_342 = buffer.data(lh + 342);
    const auto *lh_345 = buffer.data(lh + 345);
    const auto *lh_356 = buffer.data(lh + 356);
    const auto *lh_357 = buffer.data(lh + 357);
    const auto *lh_359 = buffer.data(lh + 359);
    const auto *lh_362 = buffer.data(lh + 362);
    const auto *lh_366 = buffer.data(lh + 366);
    const auto *lh_414 = buffer.data(lh + 414);
    const auto *lh_415 = buffer.data(lh + 415);
    const auto *lh_416 = buffer.data(lh + 416);
    const auto *lh_417 = buffer.data(lh + 417);
    const auto *lh_418 = buffer.data(lh + 418);
    const auto *lh_425 = buffer.data(lh + 425);
    const auto *lh_429 = buffer.data(lh + 429);
    const auto *lh_434 = buffer.data(lh + 434);
    const auto *lh_435 = buffer.data(lh + 435);
    const auto *lh_436 = buffer.data(lh + 436);
    const auto *lh_437 = buffer.data(lh + 437);
    const auto *lh_438 = buffer.data(lh + 438);
    const auto *lh_440 = buffer.data(lh + 440);
    const auto *lh_444 = buffer.data(lh + 444);
    const auto *lh_447 = buffer.data(lh + 447);
    const auto *lh_451 = buffer.data(lh + 451);
    const auto *lh_456 = buffer.data(lh + 456);
    const auto *lh_458 = buffer.data(lh + 458);
    const auto *lh_459 = buffer.data(lh + 459);
    const auto *lh_460 = buffer.data(lh + 460);
    const auto *lh_461 = buffer.data(lh + 461);
    const auto *lh_478 = buffer.data(lh + 478);
    const auto *lh_479 = buffer.data(lh + 479);
    const auto *lh_480 = buffer.data(lh + 480);
    const auto *lh_481 = buffer.data(lh + 481);
    const auto *lh_482 = buffer.data(lh + 482);
    const auto *lh_495 = buffer.data(lh + 495);
    const auto *lh_498 = buffer.data(lh + 498);
    const auto *lh_499 = buffer.data(lh + 499);
    const auto *lh_500 = buffer.data(lh + 500);
    const auto *lh_501 = buffer.data(lh + 501);
    const auto *lh_502 = buffer.data(lh + 502);

    const auto *li_392 = buffer.data(li + 392);
    const auto *li_394 = buffer.data(li + 394);
    const auto *li_395 = buffer.data(li + 395);
    const auto *li_397 = buffer.data(li + 397);
    const auto *li_398 = buffer.data(li + 398);
    const auto *li_401 = buffer.data(li + 401);
    const auto *li_402 = buffer.data(li + 402);
    const auto *li_404 = buffer.data(li + 404);
    const auto *li_406 = buffer.data(li + 406);
    const auto *li_412 = buffer.data(li + 412);
    const auto *li_413 = buffer.data(li + 413);
    const auto *li_415 = buffer.data(li + 415);
    const auto *li_416 = buffer.data(li + 416);
    const auto *li_417 = buffer.data(li + 417);
    const auto *li_419 = buffer.data(li + 419);
    const auto *li_420 = buffer.data(li + 420);
    const auto *li_421 = buffer.data(li + 421);
    const auto *li_423 = buffer.data(li + 423);
    const auto *li_425 = buffer.data(li + 425);
    const auto *li_426 = buffer.data(li + 426);
    const auto *li_429 = buffer.data(li + 429);
    const auto *li_430 = buffer.data(li + 430);
    const auto *li_432 = buffer.data(li + 432);
    const auto *li_434 = buffer.data(li + 434);
    const auto *li_435 = buffer.data(li + 435);
    const auto *li_441 = buffer.data(li + 441);
    const auto *li_443 = buffer.data(li + 443);
    const auto *li_444 = buffer.data(li + 444);
    const auto *li_445 = buffer.data(li + 445);
    const auto *li_447 = buffer.data(li + 447);
    const auto *li_451 = buffer.data(li + 451);
    const auto *li_454 = buffer.data(li + 454);
    const auto *li_458 = buffer.data(li + 458);
    const auto *li_476 = buffer.data(li + 476);
    const auto *li_481 = buffer.data(li + 481);
    const auto *li_485 = buffer.data(li + 485);
    const auto *li_490 = buffer.data(li + 490);
    const auto *li_587 = buffer.data(li + 587);
    const auto *li_609 = buffer.data(li + 609);

    const auto *mg0_300 = buffer.data(mg0 + 300);
    const auto *mg0_301 = buffer.data(mg0 + 301);
    const auto *mg0_303 = buffer.data(mg0 + 303);
    const auto *mg0_305 = buffer.data(mg0 + 305);
    const auto *mg0_309 = buffer.data(mg0 + 309);
    const auto *mg0_310 = buffer.data(mg0 + 310);
    const auto *mg0_312 = buffer.data(mg0 + 312);
    const auto *mg0_313 = buffer.data(mg0 + 313);
    const auto *mg0_314 = buffer.data(mg0 + 314);
    const auto *mg0_315 = buffer.data(mg0 + 315);
    const auto *mg0_317 = buffer.data(mg0 + 317);
    const auto *mg0_318 = buffer.data(mg0 + 318);
    const auto *mg0_320 = buffer.data(mg0 + 320);
    const auto *mg0_321 = buffer.data(mg0 + 321);
    const auto *mg0_325 = buffer.data(mg0 + 325);
    const auto *mg0_326 = buffer.data(mg0 + 326);
    const auto *mg0_327 = buffer.data(mg0 + 327);
    const auto *mg0_329 = buffer.data(mg0 + 329);
    const auto *mg0_357 = buffer.data(mg0 + 357);

    const auto *mg1_300 = buffer.data(mg1 + 300);
    const auto *mg1_301 = buffer.data(mg1 + 301);
    const auto *mg1_303 = buffer.data(mg1 + 303);
    const auto *mg1_305 = buffer.data(mg1 + 305);
    const auto *mg1_309 = buffer.data(mg1 + 309);
    const auto *mg1_310 = buffer.data(mg1 + 310);
    const auto *mg1_312 = buffer.data(mg1 + 312);
    const auto *mg1_313 = buffer.data(mg1 + 313);
    const auto *mg1_314 = buffer.data(mg1 + 314);
    const auto *mg1_315 = buffer.data(mg1 + 315);
    const auto *mg1_317 = buffer.data(mg1 + 317);
    const auto *mg1_318 = buffer.data(mg1 + 318);
    const auto *mg1_320 = buffer.data(mg1 + 320);
    const auto *mg1_321 = buffer.data(mg1 + 321);
    const auto *mg1_325 = buffer.data(mg1 + 325);
    const auto *mg1_326 = buffer.data(mg1 + 326);
    const auto *mg1_327 = buffer.data(mg1 + 327);
    const auto *mg1_329 = buffer.data(mg1 + 329);
    const auto *mg1_357 = buffer.data(mg1 + 357);

    const auto *mh_399 = buffer.data(mh + 399);
    const auto *mh_401 = buffer.data(mh + 401);
    const auto *mh_402 = buffer.data(mh + 402);
    const auto *mh_404 = buffer.data(mh + 404);
    const auto *mh_405 = buffer.data(mh + 405);
    const auto *mh_408 = buffer.data(mh + 408);
    const auto *mh_414 = buffer.data(mh + 414);
    const auto *mh_415 = buffer.data(mh + 415);
    const auto *mh_416 = buffer.data(mh + 416);
    const auto *mh_417 = buffer.data(mh + 417);
    const auto *mh_418 = buffer.data(mh + 418);
    const auto *mh_419 = buffer.data(mh + 419);
    const auto *mh_420 = buffer.data(mh + 420);
    const auto *mh_421 = buffer.data(mh + 421);
    const auto *mh_422 = buffer.data(mh + 422);
    const auto *mh_423 = buffer.data(mh + 423);
    const auto *mh_425 = buffer.data(mh + 425);
    const auto *mh_426 = buffer.data(mh + 426);
    const auto *mh_428 = buffer.data(mh + 428);
    const auto *mh_429 = buffer.data(mh + 429);
    const auto *mh_434 = buffer.data(mh + 434);
    const auto *mh_435 = buffer.data(mh + 435);
    const auto *mh_436 = buffer.data(mh + 436);
    const auto *mh_437 = buffer.data(mh + 437);
    const auto *mh_438 = buffer.data(mh + 438);
    const auto *mh_439 = buffer.data(mh + 439);
    const auto *mh_440 = buffer.data(mh + 440);
    const auto *mh_441 = buffer.data(mh + 441);
    const auto *mh_442 = buffer.data(mh + 442);
    const auto *mh_443 = buffer.data(mh + 443);
    const auto *mh_444 = buffer.data(mh + 444);
    const auto *mh_446 = buffer.data(mh + 446);
    const auto *mh_447 = buffer.data(mh + 447);
    const auto *mh_448 = buffer.data(mh + 448);
    const auto *mh_450 = buffer.data(mh + 450);
    const auto *mh_451 = buffer.data(mh + 451);
    const auto *mh_456 = buffer.data(mh + 456);
    const auto *mh_457 = buffer.data(mh + 457);
    const auto *mh_458 = buffer.data(mh + 458);
    const auto *mh_459 = buffer.data(mh + 459);
    const auto *mh_460 = buffer.data(mh + 460);
    const auto *mh_461 = buffer.data(mh + 461);
    const auto *mh_462 = buffer.data(mh + 462);
    const auto *mh_464 = buffer.data(mh + 464);
    const auto *mh_465 = buffer.data(mh + 465);
    const auto *mh_467 = buffer.data(mh + 467);
    const auto *mh_468 = buffer.data(mh + 468);
    const auto *mh_471 = buffer.data(mh + 471);
    const auto *mh_477 = buffer.data(mh + 477);
    const auto *mh_478 = buffer.data(mh + 478);
    const auto *mh_479 = buffer.data(mh + 479);
    const auto *mh_480 = buffer.data(mh + 480);
    const auto *mh_481 = buffer.data(mh + 481);
    const auto *mh_482 = buffer.data(mh + 482);
    const auto *mh_483 = buffer.data(mh + 483);
    const auto *mh_485 = buffer.data(mh + 485);
    const auto *mh_486 = buffer.data(mh + 486);
    const auto *mh_488 = buffer.data(mh + 488);
    const auto *mh_489 = buffer.data(mh + 489);
    const auto *mh_492 = buffer.data(mh + 492);
    const auto *mh_495 = buffer.data(mh + 495);
    const auto *mh_498 = buffer.data(mh + 498);
    const auto *mh_499 = buffer.data(mh + 499);
    const auto *mh_500 = buffer.data(mh + 500);
    const auto *mh_501 = buffer.data(mh + 501);
    const auto *mh_502 = buffer.data(mh + 502);

#pragma omp simd aligned(t_533, t_534, t_535, t_536, t_537, pa_y, pb_y, lh_294, lh_295, \
                         lh_296, li_394, li_395, li_397, mh_399, \
                         mh_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_9 * lh_294[k]
                   + pb_y[k] * mh_399[k];

        t_534[k] = pa_y[k] * li_394[k];

        t_535[k] = f_10 * lh_295[k]
                   + pa_y[k] * li_395[k];

        t_536[k] = f_9 * lh_296[k]
                   + pb_y[k] * mh_401[k];

        t_537[k] = pa_y[k] * li_397[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, pa_y, pb_y, pb_z, lh_276, lh_297, lh_299, \
                         li_398, li_401, mh_402, mh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_11 * lh_297[k]
                   + pa_y[k] * li_398[k];

        t_539[k] = f_12 * lh_276[k]
                   + pb_z[k] * mh_402[k];

        t_540[k] = f_9 * lh_299[k]
                   + pb_y[k] * mh_404[k];

        t_541[k] = pa_y[k] * li_401[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, pa_y, pb_y, pb_z, lh_279, lh_300, lh_302, \
                         lh_303, li_402, li_404, mh_405, mh_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_12 * lh_300[k]
                   + pa_y[k] * li_402[k];

        t_543[k] = f_12 * lh_279[k]
                   + pb_z[k] * mh_405[k];

        t_544[k] = f_10 * lh_302[k]
                   + pa_y[k] * li_404[k];

        t_545[k] = f_9 * lh_303[k]
                   + pb_y[k] * mh_408[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, t_550, pa_y, pb_x, lh_414, lh_415, \
                         lh_416, lh_417, li_406, mh_414, mh_415, mh_416, \
                         mh_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = pa_y[k] * li_406[k];

        t_547[k] = f_12 * lh_414[k]
                   + pb_x[k] * mh_414[k];

        t_548[k] = f_12 * lh_415[k]
                   + pb_x[k] * mh_415[k];

        t_549[k] = f_12 * lh_416[k]
                   + pb_x[k] * mh_416[k];

        t_550[k] = f_12 * lh_417[k]
                   + pb_x[k] * mh_417[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, t_554, pa_y, pb_x, pb_z, lh_288, lh_309, lh_418, \
                         li_412, li_413, mh_414, mh_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_12 * lh_418[k]
                   + pb_x[k] * mh_418[k];

        t_552[k] = pa_y[k] * li_412[k];

        t_553[k] = f_14 * lh_309[k]
                   + pa_y[k] * li_413[k];

        t_554[k] = f_12 * lh_288[k]
                   + pb_z[k] * mh_414[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, pa_y, pb_y, lh_311, lh_312, \
                         lh_313, lh_314, li_415, li_416, li_417, li_419, \
                         mh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_12 * lh_311[k]
                   + pa_y[k] * li_415[k];

        t_556[k] = f_11 * lh_312[k]
                   + pa_y[k] * li_416[k];

        t_557[k] = f_10 * lh_313[k]
                   + pa_y[k] * li_417[k];

        t_558[k] = f_9 * lh_314[k]
                   + pb_y[k] * mh_419[k];

        t_559[k] = pa_y[k] * li_419[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, pa_z, pb_y, pb_z, ki0_252, ki1_252, \
                         lh_294, li_392, mg0_300, mg1_300, mh_420, \
                         mh_421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_27 * ki0_252[k]
                   - f_28 * ki1_252[k]
                   + pa_z[k] * li_392[k];

        t_561[k] = pb_y[k] * mh_420[k];

        t_562[k] = f_26 * lh_294[k]
                   + pb_z[k] * mh_420[k];

        t_563[k] = f_3 * mg0_300[k]
                   - f_4 * mg1_300[k]
                   + pb_y[k] * mh_421[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, t_567, pb_x, pb_y, pb_z, lh_297, lh_425, \
                         mg0_301, mg0_305, mg1_301, mg1_305, mh_422, mh_423, \
                         mh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = pb_y[k] * mh_422[k];

        t_565[k] = f_12 * lh_425[k]
                   + f_7 * mg0_305[k]
                   - f_8 * mg1_305[k]
                   + pb_x[k] * mh_425[k];

        t_566[k] = f_5 * mg0_301[k]
                   - f_6 * mg1_301[k]
                   + pb_y[k] * mh_423[k];

        t_567[k] = f_26 * lh_297[k]
                   + pb_z[k] * mh_423[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, t_571, pb_x, pb_y, pb_z, lh_300, lh_429, \
                         mg0_303, mg0_309, mg1_303, mg1_309, mh_425, mh_426, \
                         mh_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = pb_y[k] * mh_425[k];

        t_569[k] = f_12 * lh_429[k]
                   + f_5 * mg0_309[k]
                   - f_6 * mg1_309[k]
                   + pb_x[k] * mh_429[k];

        t_570[k] = f_7 * mg0_303[k]
                   - f_8 * mg1_303[k]
                   + pb_y[k] * mh_426[k];

        t_571[k] = f_26 * lh_300[k]
                   + pb_z[k] * mh_426[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pb_x, pb_y, lh_434, lh_435, mg0_305, \
                         mg0_314, mg1_305, mg1_314, mh_428, mh_429, mh_434, \
                         mh_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_3 * mg0_305[k]
                   - f_4 * mg1_305[k]
                   + pb_y[k] * mh_428[k];

        t_573[k] = pb_y[k] * mh_429[k];

        t_574[k] = f_12 * lh_434[k]
                   + f_3 * mg0_314[k]
                   - f_4 * mg1_314[k]
                   + pb_x[k] * mh_434[k];

        t_575[k] = f_12 * lh_435[k]
                   + pb_x[k] * mh_435[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, pb_x, pb_y, lh_436, lh_437, \
                         lh_438, lh_440, mh_434, mh_436, mh_437, mh_438, \
                         mh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_12 * lh_436[k]
                   + pb_x[k] * mh_436[k];

        t_577[k] = f_12 * lh_437[k]
                   + pb_x[k] * mh_437[k];

        t_578[k] = f_12 * lh_438[k]
                   + pb_x[k] * mh_438[k];

        t_579[k] = pb_y[k] * mh_434[k];

        t_580[k] = f_12 * lh_440[k]
                   + pb_x[k] * mh_440[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, t_584, pb_y, pb_z, lh_309, mg0_310, mg0_312, \
                         mg0_313, mg1_310, mg1_312, mg1_313, mh_435, mh_437, \
                         mh_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_1 * mg0_310[k]
                   - f_2 * mg1_310[k]
                   + pb_y[k] * mh_435[k];

        t_582[k] = f_26 * lh_309[k]
                   + pb_z[k] * mh_435[k];

        t_583[k] = f_7 * mg0_312[k]
                   - f_8 * mg1_312[k]
                   + pb_y[k] * mh_437[k];

        t_584[k] = f_5 * mg0_313[k]
                   - f_6 * mg1_313[k]
                   + pb_y[k] * mh_438[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, pa_x, pb_y, ki0_587, ki1_587, li_587, mg0_314, \
                         mg1_314, mh_439, mh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_3 * mg0_314[k]
                   - f_4 * mg1_314[k]
                   + pb_y[k] * mh_439[k];

        t_586[k] = pb_y[k] * mh_440[k];

        t_587[k] = f_24 * ki0_587[k]
                   - f_25 * ki1_587[k]
                   + pa_x[k] * li_587[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, pa_y, pb_y, pb_z, ki0_280, ki1_280, lh_315, \
                         li_420, mh_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_22 * ki0_280[k]
                   - f_23 * ki1_280[k]
                   + pa_y[k] * li_420[k];

        t_589[k] = f_14 * lh_315[k]
                   + pb_y[k] * mh_441[k];

        t_590[k] = pb_z[k] * mh_441[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, pb_x, pb_z, lh_444, mg0_315, mg0_318, mg1_315, \
                         mg1_318, mh_442, mh_443, mh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_11 * lh_444[k]
                   + f_7 * mg0_318[k]
                   - f_8 * mg1_318[k]
                   + pb_x[k] * mh_444[k];

        t_592[k] = pb_z[k] * mh_442[k];

        t_593[k] = f_3 * mg0_315[k]
                   - f_4 * mg1_315[k]
                   + pb_z[k] * mh_443[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, pb_x, pb_y, pb_z, lh_320, lh_447, \
                         mg0_317, mg0_321, mg1_317, mg1_321, mh_444, mh_446, \
                         mh_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_11 * lh_447[k]
                   + f_5 * mg0_321[k]
                   - f_6 * mg1_321[k]
                   + pb_x[k] * mh_447[k];

        t_595[k] = pb_z[k] * mh_444[k];

        t_596[k] = f_14 * lh_320[k]
                   + pb_y[k] * mh_446[k];

        t_597[k] = f_5 * mg0_317[k]
                   - f_6 * mg1_317[k]
                   + pb_z[k] * mh_446[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, pb_x, pb_z, lh_451, mg0_318, mg0_325, mg1_318, \
                         mg1_325, mh_447, mh_448, mh_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_11 * lh_451[k]
                   + f_3 * mg0_325[k]
                   - f_4 * mg1_325[k]
                   + pb_x[k] * mh_451[k];

        t_599[k] = pb_z[k] * mh_447[k];

        t_600[k] = f_3 * mg0_318[k]
                   - f_4 * mg1_318[k]
                   + pb_z[k] * mh_448[k];
    }

#pragma omp simd aligned(t_601, t_602, t_603, t_604, pb_x, pb_y, pb_z, lh_324, lh_456, \
                         mg0_320, mg1_320, mh_450, mh_451, mh_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = f_14 * lh_324[k]
                   + pb_y[k] * mh_450[k];

        t_602[k] = f_7 * mg0_320[k]
                   - f_8 * mg1_320[k]
                   + pb_z[k] * mh_450[k];

        t_603[k] = f_11 * lh_456[k]
                   + pb_x[k] * mh_456[k];

        t_604[k] = pb_z[k] * mh_451[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, pb_x, lh_458, lh_459, lh_460, lh_461, \
                         mh_458, mh_459, mh_460, mh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = f_11 * lh_458[k]
                   + pb_x[k] * mh_458[k];

        t_606[k] = f_11 * lh_459[k]
                   + pb_x[k] * mh_459[k];

        t_607[k] = f_11 * lh_460[k]
                   + pb_x[k] * mh_460[k];

        t_608[k] = f_11 * lh_461[k]
                   + pb_x[k] * mh_461[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, pa_x, pb_z, ki0_609, ki1_609, li_609, \
                         mg0_325, mg0_326, mg1_325, mg1_326, mh_456, mh_457, \
                         mh_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_20 * ki0_609[k]
                   - f_21 * ki1_609[k]
                   + pa_x[k] * li_609[k];

        t_610[k] = pb_z[k] * mh_456[k];

        t_611[k] = f_3 * mg0_325[k]
                   - f_4 * mg1_325[k]
                   + pb_z[k] * mh_457[k];

        t_612[k] = f_5 * mg0_326[k]
                   - f_6 * mg1_326[k]
                   + pb_z[k] * mh_458[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pa_z, pb_y, pb_z, lh_335, li_420, \
                         mg0_327, mg0_329, mg1_327, mg1_329, mh_459, \
                         mh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_7 * mg0_327[k]
                   - f_8 * mg1_327[k]
                   + pb_z[k] * mh_459[k];

        t_614[k] = f_14 * lh_335[k]
                   + pb_y[k] * mh_461[k];

        t_615[k] = f_1 * mg0_329[k]
                   - f_2 * mg1_329[k]
                   + pb_z[k] * mh_461[k];

        t_616[k] = pa_z[k] * li_420[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, t_621, pa_z, pb_y, pb_z, lh_315, lh_317, \
                         lh_338, li_421, li_423, li_425, mh_462, \
                         mh_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = pa_z[k] * li_421[k];

        t_618[k] = f_9 * lh_315[k]
                   + pb_z[k] * mh_462[k];

        t_619[k] = pa_z[k] * li_423[k];

        t_620[k] = f_26 * lh_338[k]
                   + pb_y[k] * mh_464[k];

        t_621[k] = f_10 * lh_317[k]
                   + pa_z[k] * li_425[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, t_626, pa_z, pb_y, pb_z, lh_318, lh_320, \
                         lh_341, li_426, li_429, li_430, mh_465, \
                         mh_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = pa_z[k] * li_426[k];

        t_623[k] = f_9 * lh_318[k]
                   + pb_z[k] * mh_465[k];

        t_624[k] = f_26 * lh_341[k]
                   + pb_y[k] * mh_467[k];

        t_625[k] = f_11 * lh_320[k]
                   + pa_z[k] * li_429[k];

        t_626[k] = pa_z[k] * li_430[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, t_630, pa_z, pb_y, pb_z, lh_321, lh_322, lh_324, \
                         lh_345, li_432, li_434, mh_468, mh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_9 * lh_321[k]
                   + pb_z[k] * mh_468[k];

        t_628[k] = f_10 * lh_322[k]
                   + pa_z[k] * li_432[k];

        t_629[k] = f_26 * lh_345[k]
                   + pb_y[k] * mh_471[k];

        t_630[k] = f_12 * lh_324[k]
                   + pa_z[k] * li_434[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, t_635, pa_z, pb_x, lh_478, lh_479, \
                         lh_480, lh_481, li_435, mh_478, mh_479, mh_480, \
                         mh_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = pa_z[k] * li_435[k];

        t_632[k] = f_11 * lh_478[k]
                   + pb_x[k] * mh_478[k];

        t_633[k] = f_11 * lh_479[k]
                   + pb_x[k] * mh_479[k];

        t_634[k] = f_11 * lh_480[k]
                   + pb_x[k] * mh_480[k];

        t_635[k] = f_11 * lh_481[k]
                   + pb_x[k] * mh_481[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pa_z, pb_x, pb_z, lh_330, lh_331, lh_482, \
                         li_441, li_443, mh_477, mh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_11 * lh_482[k]
                   + pb_x[k] * mh_482[k];

        t_637[k] = pa_z[k] * li_441[k];

        t_638[k] = f_9 * lh_330[k]
                   + pb_z[k] * mh_477[k];

        t_639[k] = f_10 * lh_331[k]
                   + pa_z[k] * li_443[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, pa_z, pb_y, lh_332, lh_333, lh_335, \
                         lh_356, li_444, li_445, li_447, mh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_11 * lh_332[k]
                   + pa_z[k] * li_444[k];

        t_641[k] = f_12 * lh_333[k]
                   + pa_z[k] * li_445[k];

        t_642[k] = f_26 * lh_356[k]
                   + pb_y[k] * mh_482[k];

        t_643[k] = f_14 * lh_335[k]
                   + pa_z[k] * li_447[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, pa_y, pb_y, pb_z, ki0_336, ki1_336, lh_336, \
                         lh_357, li_476, mh_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_24 * ki0_336[k]
                   - f_25 * ki1_336[k]
                   + pa_y[k] * li_476[k];

        t_645[k] = f_12 * lh_357[k]
                   + pb_y[k] * mh_483[k];

        t_646[k] = f_10 * lh_336[k]
                   + pb_z[k] * mh_483[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pa_y, pa_z, pb_y, ki0_283, ki0_341, ki1_283, \
                         ki1_341, lh_359, li_451, li_481, mh_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_15 * ki0_283[k]
                   - f_16 * ki1_283[k]
                   + pa_z[k] * li_451[k];

        t_648[k] = f_12 * lh_359[k]
                   + pb_y[k] * mh_485[k];

        t_649[k] = f_24 * ki0_341[k]
                   - f_25 * ki1_341[k]
                   + pa_y[k] * li_481[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pa_z, pb_y, pb_z, ki0_286, ki1_286, lh_339, \
                         lh_362, li_454, mh_486, mh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_15 * ki0_286[k]
                   - f_16 * ki1_286[k]
                   + pa_z[k] * li_454[k];

        t_651[k] = f_10 * lh_339[k]
                   + pb_z[k] * mh_486[k];

        t_652[k] = f_12 * lh_362[k]
                   + pb_y[k] * mh_488[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pa_y, pa_z, pb_z, ki0_290, ki0_345, ki1_290, \
                         ki1_345, lh_342, li_458, li_485, mh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_24 * ki0_345[k]
                   - f_25 * ki1_345[k]
                   + pa_y[k] * li_485[k];

        t_654[k] = f_15 * ki0_290[k]
                   - f_16 * ki1_290[k]
                   + pa_z[k] * li_458[k];

        t_655[k] = f_10 * lh_342[k]
                   + pb_z[k] * mh_489[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pa_y, pb_x, pb_y, ki0_350, ki1_350, lh_366, \
                         lh_495, li_490, mg0_357, mg1_357, mh_492, \
                         mh_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_11 * lh_495[k]
                   + f_3 * mg0_357[k]
                   - f_4 * mg1_357[k]
                   + pb_x[k] * mh_495[k];

        t_657[k] = f_12 * lh_366[k]
                   + pb_y[k] * mh_492[k];

        t_658[k] = f_24 * ki0_350[k]
                   - f_25 * ki1_350[k]
                   + pa_y[k] * li_490[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, pb_x, lh_498, lh_499, lh_500, \
                         lh_501, lh_502, mh_498, mh_499, mh_500, mh_501, \
                         mh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_11 * lh_498[k]
                   + pb_x[k] * mh_498[k];

        t_660[k] = f_11 * lh_499[k]
                   + pb_x[k] * mh_499[k];

        t_661[k] = f_11 * lh_500[k]
                   + pb_x[k] * mh_500[k];

        t_662[k] = f_11 * lh_501[k]
                   + pb_x[k] * mh_501[k];

        t_663[k] = f_11 * lh_502[k]
                   + pb_x[k] * mh_502[k];
    }
}

static auto
compute_prim_mi_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t ki0,
                                            const size_t ki1, const size_t lh, const size_t li,
                                            const size_t mg0, const size_t mg1, const size_t mh,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_14 = 3.0 / p;
    const auto f_15 = 0.5 / alpha;
    const auto f_16 = 0.5 * beta / (alpha * p);
    const auto f_17 = 3.5 / p;
    const auto f_18 = 3.0 / alpha;
    const auto f_19 = 3.0 * beta / (alpha * p);
    const auto f_20 = 1.0 / alpha;
    const auto f_21 = beta / (alpha * p);
    const auto f_22 = 2.5 / alpha;
    const auto f_23 = 2.5 * beta / (alpha * p);
    const auto f_24 = 1.5 / alpha;
    const auto f_25 = 1.5 * beta / (alpha * p);
    const auto f_26 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ki0_311 = buffer.data(ki0 + 311);
    const auto *ki0_314 = buffer.data(ki0 + 314);
    const auto *ki0_318 = buffer.data(ki0 + 318);
    const auto *ki0_339 = buffer.data(ki0 + 339);
    const auto *ki0_342 = buffer.data(ki0 + 342);
    const auto *ki0_346 = buffer.data(ki0 + 346);
    const auto *ki0_364 = buffer.data(ki0 + 364);
    const auto *ki0_369 = buffer.data(ki0 + 369);
    const auto *ki0_373 = buffer.data(ki0 + 373);
    const auto *ki0_378 = buffer.data(ki0 + 378);
    const auto *ki0_392 = buffer.data(ki0 + 392);
    const auto *ki0_397 = buffer.data(ki0 + 397);
    const auto *ki0_401 = buffer.data(ki0 + 401);
    const auto *ki0_406 = buffer.data(ki0 + 406);
    const auto *ki0_420 = buffer.data(ki0 + 420);
    const auto *ki0_665 = buffer.data(ki0 + 665);
    const auto *ki0_667 = buffer.data(ki0 + 667);
    const auto *ki0_668 = buffer.data(ki0 + 668);
    const auto *ki0_669 = buffer.data(ki0 + 669);
    const auto *ki0_671 = buffer.data(ki0 + 671);
    const auto *ki0_693 = buffer.data(ki0 + 693);
    const auto *ki0_695 = buffer.data(ki0 + 695);
    const auto *ki0_696 = buffer.data(ki0 + 696);
    const auto *ki0_697 = buffer.data(ki0 + 697);
    const auto *ki0_699 = buffer.data(ki0 + 699);
    const auto *ki0_721 = buffer.data(ki0 + 721);
    const auto *ki0_723 = buffer.data(ki0 + 723);
    const auto *ki0_724 = buffer.data(ki0 + 724);
    const auto *ki0_725 = buffer.data(ki0 + 725);
    const auto *ki0_727 = buffer.data(ki0 + 727);
    const auto *ki0_783 = buffer.data(ki0 + 783);

    const auto *ki1_311 = buffer.data(ki1 + 311);
    const auto *ki1_314 = buffer.data(ki1 + 314);
    const auto *ki1_318 = buffer.data(ki1 + 318);
    const auto *ki1_339 = buffer.data(ki1 + 339);
    const auto *ki1_342 = buffer.data(ki1 + 342);
    const auto *ki1_346 = buffer.data(ki1 + 346);
    const auto *ki1_364 = buffer.data(ki1 + 364);
    const auto *ki1_369 = buffer.data(ki1 + 369);
    const auto *ki1_373 = buffer.data(ki1 + 373);
    const auto *ki1_378 = buffer.data(ki1 + 378);
    const auto *ki1_392 = buffer.data(ki1 + 392);
    const auto *ki1_397 = buffer.data(ki1 + 397);
    const auto *ki1_401 = buffer.data(ki1 + 401);
    const auto *ki1_406 = buffer.data(ki1 + 406);
    const auto *ki1_420 = buffer.data(ki1 + 420);
    const auto *ki1_665 = buffer.data(ki1 + 665);
    const auto *ki1_667 = buffer.data(ki1 + 667);
    const auto *ki1_668 = buffer.data(ki1 + 668);
    const auto *ki1_669 = buffer.data(ki1 + 669);
    const auto *ki1_671 = buffer.data(ki1 + 671);
    const auto *ki1_693 = buffer.data(ki1 + 693);
    const auto *ki1_695 = buffer.data(ki1 + 695);
    const auto *ki1_696 = buffer.data(ki1 + 696);
    const auto *ki1_697 = buffer.data(ki1 + 697);
    const auto *ki1_699 = buffer.data(ki1 + 699);
    const auto *ki1_721 = buffer.data(ki1 + 721);
    const auto *ki1_723 = buffer.data(ki1 + 723);
    const auto *ki1_724 = buffer.data(ki1 + 724);
    const auto *ki1_725 = buffer.data(ki1 + 725);
    const auto *ki1_727 = buffer.data(ki1 + 727);
    const auto *ki1_783 = buffer.data(ki1 + 783);

    const auto *lh_351 = buffer.data(lh + 351);
    const auto *lh_357 = buffer.data(lh + 357);
    const auto *lh_360 = buffer.data(lh + 360);
    const auto *lh_363 = buffer.data(lh + 363);
    const auto *lh_372 = buffer.data(lh + 372);
    const auto *lh_377 = buffer.data(lh + 377);
    const auto *lh_378 = buffer.data(lh + 378);
    const auto *lh_380 = buffer.data(lh + 380);
    const auto *lh_381 = buffer.data(lh + 381);
    const auto *lh_383 = buffer.data(lh + 383);
    const auto *lh_384 = buffer.data(lh + 384);
    const auto *lh_387 = buffer.data(lh + 387);
    const auto *lh_393 = buffer.data(lh + 393);
    const auto *lh_398 = buffer.data(lh + 398);
    const auto *lh_399 = buffer.data(lh + 399);
    const auto *lh_401 = buffer.data(lh + 401);
    const auto *lh_402 = buffer.data(lh + 402);
    const auto *lh_404 = buffer.data(lh + 404);
    const auto *lh_405 = buffer.data(lh + 405);
    const auto *lh_408 = buffer.data(lh + 408);
    const auto *lh_414 = buffer.data(lh + 414);
    const auto *lh_419 = buffer.data(lh + 419);
    const auto *lh_420 = buffer.data(lh + 420);
    const auto *lh_421 = buffer.data(lh + 421);
    const auto *lh_422 = buffer.data(lh + 422);
    const auto *lh_423 = buffer.data(lh + 423);
    const auto *lh_425 = buffer.data(lh + 425);
    const auto *lh_426 = buffer.data(lh + 426);
    const auto *lh_428 = buffer.data(lh + 428);
    const auto *lh_429 = buffer.data(lh + 429);
    const auto *lh_435 = buffer.data(lh + 435);
    const auto *lh_437 = buffer.data(lh + 437);
    const auto *lh_438 = buffer.data(lh + 438);
    const auto *lh_439 = buffer.data(lh + 439);
    const auto *lh_440 = buffer.data(lh + 440);
    const auto *lh_441 = buffer.data(lh + 441);
    const auto *lh_503 = buffer.data(lh + 503);
    const auto *lh_516 = buffer.data(lh + 516);
    const auto *lh_519 = buffer.data(lh + 519);
    const auto *lh_520 = buffer.data(lh + 520);
    const auto *lh_521 = buffer.data(lh + 521);
    const auto *lh_522 = buffer.data(lh + 522);
    const auto *lh_523 = buffer.data(lh + 523);
    const auto *lh_524 = buffer.data(lh + 524);
    const auto *lh_537 = buffer.data(lh + 537);
    const auto *lh_540 = buffer.data(lh + 540);
    const auto *lh_541 = buffer.data(lh + 541);
    const auto *lh_542 = buffer.data(lh + 542);
    const auto *lh_543 = buffer.data(lh + 543);
    const auto *lh_544 = buffer.data(lh + 544);
    const auto *lh_545 = buffer.data(lh + 545);
    const auto *lh_561 = buffer.data(lh + 561);
    const auto *lh_562 = buffer.data(lh + 562);
    const auto *lh_563 = buffer.data(lh + 563);
    const auto *lh_564 = buffer.data(lh + 564);
    const auto *lh_565 = buffer.data(lh + 565);
    const auto *lh_572 = buffer.data(lh + 572);
    const auto *lh_576 = buffer.data(lh + 576);
    const auto *lh_581 = buffer.data(lh + 581);
    const auto *lh_582 = buffer.data(lh + 582);
    const auto *lh_583 = buffer.data(lh + 583);
    const auto *lh_584 = buffer.data(lh + 584);
    const auto *lh_585 = buffer.data(lh + 585);
    const auto *lh_587 = buffer.data(lh + 587);

    const auto *li_479 = buffer.data(li + 479);
    const auto *li_482 = buffer.data(li + 482);
    const auto *li_486 = buffer.data(li + 486);
    const auto *li_504 = buffer.data(li + 504);
    const auto *li_507 = buffer.data(li + 507);
    const auto *li_509 = buffer.data(li + 509);
    const auto *li_510 = buffer.data(li + 510);
    const auto *li_513 = buffer.data(li + 513);
    const auto *li_514 = buffer.data(li + 514);
    const auto *li_518 = buffer.data(li + 518);
    const auto *li_532 = buffer.data(li + 532);
    const auto *li_537 = buffer.data(li + 537);
    const auto *li_541 = buffer.data(li + 541);
    const auto *li_546 = buffer.data(li + 546);
    const auto *li_560 = buffer.data(li + 560);
    const auto *li_562 = buffer.data(li + 562);
    const auto *li_563 = buffer.data(li + 563);
    const auto *li_565 = buffer.data(li + 565);
    const auto *li_566 = buffer.data(li + 566);
    const auto *li_569 = buffer.data(li + 569);
    const auto *li_570 = buffer.data(li + 570);
    const auto *li_572 = buffer.data(li + 572);
    const auto *li_574 = buffer.data(li + 574);
    const auto *li_580 = buffer.data(li + 580);
    const auto *li_581 = buffer.data(li + 581);
    const auto *li_583 = buffer.data(li + 583);
    const auto *li_584 = buffer.data(li + 584);
    const auto *li_585 = buffer.data(li + 585);
    const auto *li_587 = buffer.data(li + 587);
    const auto *li_588 = buffer.data(li + 588);
    const auto *li_665 = buffer.data(li + 665);
    const auto *li_667 = buffer.data(li + 667);
    const auto *li_668 = buffer.data(li + 668);
    const auto *li_669 = buffer.data(li + 669);
    const auto *li_671 = buffer.data(li + 671);
    const auto *li_693 = buffer.data(li + 693);
    const auto *li_695 = buffer.data(li + 695);
    const auto *li_696 = buffer.data(li + 696);
    const auto *li_697 = buffer.data(li + 697);
    const auto *li_699 = buffer.data(li + 699);
    const auto *li_721 = buffer.data(li + 721);
    const auto *li_723 = buffer.data(li + 723);
    const auto *li_724 = buffer.data(li + 724);
    const auto *li_725 = buffer.data(li + 725);
    const auto *li_727 = buffer.data(li + 727);
    const auto *li_783 = buffer.data(li + 783);

    const auto *mg0_372 = buffer.data(mg0 + 372);
    const auto *mg0_387 = buffer.data(mg0 + 387);
    const auto *mg0_405 = buffer.data(mg0 + 405);
    const auto *mg0_406 = buffer.data(mg0 + 406);
    const auto *mg0_408 = buffer.data(mg0 + 408);
    const auto *mg0_410 = buffer.data(mg0 + 410);
    const auto *mg0_414 = buffer.data(mg0 + 414);
    const auto *mg0_415 = buffer.data(mg0 + 415);
    const auto *mg0_417 = buffer.data(mg0 + 417);
    const auto *mg0_418 = buffer.data(mg0 + 418);
    const auto *mg0_419 = buffer.data(mg0 + 419);

    const auto *mg1_372 = buffer.data(mg1 + 372);
    const auto *mg1_387 = buffer.data(mg1 + 387);
    const auto *mg1_405 = buffer.data(mg1 + 405);
    const auto *mg1_406 = buffer.data(mg1 + 406);
    const auto *mg1_408 = buffer.data(mg1 + 408);
    const auto *mg1_410 = buffer.data(mg1 + 410);
    const auto *mg1_414 = buffer.data(mg1 + 414);
    const auto *mg1_415 = buffer.data(mg1 + 415);
    const auto *mg1_417 = buffer.data(mg1 + 417);
    const auto *mg1_418 = buffer.data(mg1 + 418);
    const auto *mg1_419 = buffer.data(mg1 + 419);

    const auto *mh_498 = buffer.data(mh + 498);
    const auto *mh_503 = buffer.data(mh + 503);
    const auto *mh_504 = buffer.data(mh + 504);
    const auto *mh_506 = buffer.data(mh + 506);
    const auto *mh_507 = buffer.data(mh + 507);
    const auto *mh_509 = buffer.data(mh + 509);
    const auto *mh_510 = buffer.data(mh + 510);
    const auto *mh_513 = buffer.data(mh + 513);
    const auto *mh_516 = buffer.data(mh + 516);
    const auto *mh_519 = buffer.data(mh + 519);
    const auto *mh_520 = buffer.data(mh + 520);
    const auto *mh_521 = buffer.data(mh + 521);
    const auto *mh_522 = buffer.data(mh + 522);
    const auto *mh_523 = buffer.data(mh + 523);
    const auto *mh_524 = buffer.data(mh + 524);
    const auto *mh_525 = buffer.data(mh + 525);
    const auto *mh_527 = buffer.data(mh + 527);
    const auto *mh_528 = buffer.data(mh + 528);
    const auto *mh_530 = buffer.data(mh + 530);
    const auto *mh_531 = buffer.data(mh + 531);
    const auto *mh_534 = buffer.data(mh + 534);
    const auto *mh_537 = buffer.data(mh + 537);
    const auto *mh_540 = buffer.data(mh + 540);
    const auto *mh_541 = buffer.data(mh + 541);
    const auto *mh_542 = buffer.data(mh + 542);
    const auto *mh_543 = buffer.data(mh + 543);
    const auto *mh_544 = buffer.data(mh + 544);
    const auto *mh_545 = buffer.data(mh + 545);
    const auto *mh_546 = buffer.data(mh + 546);
    const auto *mh_548 = buffer.data(mh + 548);
    const auto *mh_549 = buffer.data(mh + 549);
    const auto *mh_551 = buffer.data(mh + 551);
    const auto *mh_552 = buffer.data(mh + 552);
    const auto *mh_555 = buffer.data(mh + 555);
    const auto *mh_561 = buffer.data(mh + 561);
    const auto *mh_562 = buffer.data(mh + 562);
    const auto *mh_563 = buffer.data(mh + 563);
    const auto *mh_564 = buffer.data(mh + 564);
    const auto *mh_565 = buffer.data(mh + 565);
    const auto *mh_566 = buffer.data(mh + 566);
    const auto *mh_567 = buffer.data(mh + 567);
    const auto *mh_568 = buffer.data(mh + 568);
    const auto *mh_569 = buffer.data(mh + 569);
    const auto *mh_570 = buffer.data(mh + 570);
    const auto *mh_572 = buffer.data(mh + 572);
    const auto *mh_573 = buffer.data(mh + 573);
    const auto *mh_575 = buffer.data(mh + 575);
    const auto *mh_576 = buffer.data(mh + 576);
    const auto *mh_581 = buffer.data(mh + 581);
    const auto *mh_582 = buffer.data(mh + 582);
    const auto *mh_583 = buffer.data(mh + 583);
    const auto *mh_584 = buffer.data(mh + 584);
    const auto *mh_585 = buffer.data(mh + 585);
    const auto *mh_586 = buffer.data(mh + 586);
    const auto *mh_587 = buffer.data(mh + 587);
    const auto *mh_588 = buffer.data(mh + 588);

#pragma omp simd aligned(t_664, t_665, t_666, pa_x, pb_x, pb_z, ki0_665, ki1_665, lh_351, \
                         lh_503, li_665, mh_498, mh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_11 * lh_503[k]
                   + pb_x[k] * mh_503[k];

        t_665[k] = f_20 * ki0_665[k]
                   - f_21 * ki1_665[k]
                   + pa_x[k] * li_665[k];

        t_666[k] = f_10 * lh_351[k]
                   + pb_z[k] * mh_498[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pa_x, ki0_667, ki0_668, ki0_669, ki1_667, \
                         ki1_668, ki1_669, li_667, li_668, li_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_20 * ki0_667[k]
                   - f_21 * ki1_667[k]
                   + pa_x[k] * li_667[k];

        t_668[k] = f_20 * ki0_668[k]
                   - f_21 * ki1_668[k]
                   + pa_x[k] * li_668[k];

        t_669[k] = f_20 * ki0_669[k]
                   - f_21 * ki1_669[k]
                   + pa_x[k] * li_669[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pa_x, pa_y, pb_y, ki0_364, ki0_671, ki1_364, \
                         ki1_671, lh_377, li_504, li_671, mh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_12 * lh_377[k]
                   + pb_y[k] * mh_503[k];

        t_671[k] = f_20 * ki0_671[k]
                   - f_21 * ki1_671[k]
                   + pa_x[k] * li_671[k];

        t_672[k] = f_20 * ki0_364[k]
                   - f_21 * ki1_364[k]
                   + pa_y[k] * li_504[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pa_z, pb_y, pb_z, ki0_311, ki1_311, \
                         lh_357, lh_378, lh_380, li_479, mh_504, \
                         mh_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_11 * lh_378[k]
                   + pb_y[k] * mh_504[k];

        t_674[k] = f_11 * lh_357[k]
                   + pb_z[k] * mh_504[k];

        t_675[k] = f_20 * ki0_311[k]
                   - f_21 * ki1_311[k]
                   + pa_z[k] * li_479[k];

        t_676[k] = f_11 * lh_380[k]
                   + pb_y[k] * mh_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pa_y, pa_z, pb_z, ki0_314, ki0_369, ki1_314, \
                         ki1_369, lh_360, li_482, li_509, mh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_20 * ki0_369[k]
                   - f_21 * ki1_369[k]
                   + pa_y[k] * li_509[k];

        t_678[k] = f_20 * ki0_314[k]
                   - f_21 * ki1_314[k]
                   + pa_z[k] * li_482[k];

        t_679[k] = f_11 * lh_360[k]
                   + pb_z[k] * mh_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pa_y, pa_z, pb_y, ki0_318, ki0_373, ki1_318, \
                         ki1_373, lh_383, li_486, li_513, mh_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_11 * lh_383[k]
                   + pb_y[k] * mh_509[k];

        t_681[k] = f_20 * ki0_373[k]
                   - f_21 * ki1_373[k]
                   + pa_y[k] * li_513[k];

        t_682[k] = f_20 * ki0_318[k]
                   - f_21 * ki1_318[k]
                   + pa_z[k] * li_486[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pb_x, pb_y, pb_z, lh_363, lh_387, lh_516, \
                         mg0_372, mg1_372, mh_510, mh_513, mh_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_11 * lh_363[k]
                   + pb_z[k] * mh_510[k];

        t_684[k] = f_11 * lh_516[k]
                   + f_3 * mg0_372[k]
                   - f_4 * mg1_372[k]
                   + pb_x[k] * mh_516[k];

        t_685[k] = f_11 * lh_387[k]
                   + pb_y[k] * mh_513[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pa_y, pb_x, ki0_378, ki1_378, lh_519, \
                         lh_520, lh_521, li_518, mh_519, mh_520, \
                         mh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_20 * ki0_378[k]
                   - f_21 * ki1_378[k]
                   + pa_y[k] * li_518[k];

        t_687[k] = f_11 * lh_519[k]
                   + pb_x[k] * mh_519[k];

        t_688[k] = f_11 * lh_520[k]
                   + pb_x[k] * mh_520[k];

        t_689[k] = f_11 * lh_521[k]
                   + pb_x[k] * mh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pa_x, pb_x, ki0_693, ki1_693, lh_522, \
                         lh_523, lh_524, li_693, mh_522, mh_523, \
                         mh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_11 * lh_522[k]
                   + pb_x[k] * mh_522[k];

        t_691[k] = f_11 * lh_523[k]
                   + pb_x[k] * mh_523[k];

        t_692[k] = f_11 * lh_524[k]
                   + pb_x[k] * mh_524[k];

        t_693[k] = f_20 * ki0_693[k]
                   - f_21 * ki1_693[k]
                   + pa_x[k] * li_693[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pa_x, pb_z, ki0_695, ki0_696, ki1_695, ki1_696, \
                         lh_372, li_695, li_696, mh_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_11 * lh_372[k]
                   + pb_z[k] * mh_519[k];

        t_695[k] = f_20 * ki0_695[k]
                   - f_21 * ki1_695[k]
                   + pa_x[k] * li_695[k];

        t_696[k] = f_20 * ki0_696[k]
                   - f_21 * ki1_696[k]
                   + pa_x[k] * li_696[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pa_x, pb_y, ki0_697, ki0_699, ki1_697, ki1_699, \
                         lh_398, li_697, li_699, mh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_20 * ki0_697[k]
                   - f_21 * ki1_697[k]
                   + pa_x[k] * li_697[k];

        t_698[k] = f_11 * lh_398[k]
                   + pb_y[k] * mh_524[k];

        t_699[k] = f_20 * ki0_699[k]
                   - f_21 * ki1_699[k]
                   + pa_x[k] * li_699[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pa_y, pb_y, pb_z, ki0_392, ki1_392, lh_378, \
                         lh_399, li_532, mh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_15 * ki0_392[k]
                   - f_16 * ki1_392[k]
                   + pa_y[k] * li_532[k];

        t_701[k] = f_10 * lh_399[k]
                   + pb_y[k] * mh_525[k];

        t_702[k] = f_12 * lh_378[k]
                   + pb_z[k] * mh_525[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, pa_y, pa_z, pb_y, ki0_339, ki0_397, ki1_339, \
                         ki1_397, lh_401, li_507, li_537, mh_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_24 * ki0_339[k]
                   - f_25 * ki1_339[k]
                   + pa_z[k] * li_507[k];

        t_704[k] = f_10 * lh_401[k]
                   + pb_y[k] * mh_527[k];

        t_705[k] = f_15 * ki0_397[k]
                   - f_16 * ki1_397[k]
                   + pa_y[k] * li_537[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, pa_z, pb_y, pb_z, ki0_342, ki1_342, lh_381, \
                         lh_404, li_510, mh_528, mh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_24 * ki0_342[k]
                   - f_25 * ki1_342[k]
                   + pa_z[k] * li_510[k];

        t_707[k] = f_12 * lh_381[k]
                   + pb_z[k] * mh_528[k];

        t_708[k] = f_10 * lh_404[k]
                   + pb_y[k] * mh_530[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pa_y, pa_z, pb_z, ki0_346, ki0_401, ki1_346, \
                         ki1_401, lh_384, li_514, li_541, mh_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_15 * ki0_401[k]
                   - f_16 * ki1_401[k]
                   + pa_y[k] * li_541[k];

        t_710[k] = f_24 * ki0_346[k]
                   - f_25 * ki1_346[k]
                   + pa_z[k] * li_514[k];

        t_711[k] = f_12 * lh_384[k]
                   + pb_z[k] * mh_531[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, pa_y, pb_x, pb_y, ki0_406, ki1_406, lh_408, \
                         lh_537, li_546, mg0_387, mg1_387, mh_534, \
                         mh_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_11 * lh_537[k]
                   + f_3 * mg0_387[k]
                   - f_4 * mg1_387[k]
                   + pb_x[k] * mh_537[k];

        t_713[k] = f_10 * lh_408[k]
                   + pb_y[k] * mh_534[k];

        t_714[k] = f_15 * ki0_406[k]
                   - f_16 * ki1_406[k]
                   + pa_y[k] * li_546[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, pb_x, lh_540, lh_541, lh_542, \
                         lh_543, lh_544, mh_540, mh_541, mh_542, mh_543, \
                         mh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_11 * lh_540[k]
                   + pb_x[k] * mh_540[k];

        t_716[k] = f_11 * lh_541[k]
                   + pb_x[k] * mh_541[k];

        t_717[k] = f_11 * lh_542[k]
                   + pb_x[k] * mh_542[k];

        t_718[k] = f_11 * lh_543[k]
                   + pb_x[k] * mh_543[k];

        t_719[k] = f_11 * lh_544[k]
                   + pb_x[k] * mh_544[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, pa_x, pb_x, pb_z, ki0_721, ki1_721, lh_393, \
                         lh_545, li_721, mh_540, mh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_11 * lh_545[k]
                   + pb_x[k] * mh_545[k];

        t_721[k] = f_20 * ki0_721[k]
                   - f_21 * ki1_721[k]
                   + pa_x[k] * li_721[k];

        t_722[k] = f_12 * lh_393[k]
                   + pb_z[k] * mh_540[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pa_x, ki0_723, ki0_724, ki0_725, ki1_723, \
                         ki1_724, ki1_725, li_723, li_724, li_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_20 * ki0_723[k]
                   - f_21 * ki1_723[k]
                   + pa_x[k] * li_723[k];

        t_724[k] = f_20 * ki0_724[k]
                   - f_21 * ki1_724[k]
                   + pa_x[k] * li_724[k];

        t_725[k] = f_20 * ki0_725[k]
                   - f_21 * ki1_725[k]
                   + pa_x[k] * li_725[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pa_x, pa_y, pb_y, ki0_727, ki1_727, \
                         lh_419, lh_420, li_560, li_727, mh_545, \
                         mh_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_10 * lh_419[k]
                   + pb_y[k] * mh_545[k];

        t_727[k] = f_20 * ki0_727[k]
                   - f_21 * ki1_727[k]
                   + pa_x[k] * li_727[k];

        t_728[k] = pa_y[k] * li_560[k];

        t_729[k] = f_9 * lh_420[k]
                   + pb_y[k] * mh_546[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, pa_y, pb_y, lh_421, lh_422, \
                         lh_423, li_562, li_563, li_565, li_566, \
                         mh_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = pa_y[k] * li_562[k];

        t_731[k] = f_10 * lh_421[k]
                   + pa_y[k] * li_563[k];

        t_732[k] = f_9 * lh_422[k]
                   + pb_y[k] * mh_548[k];

        t_733[k] = pa_y[k] * li_565[k];

        t_734[k] = f_11 * lh_423[k]
                   + pa_y[k] * li_566[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, pa_y, pb_y, pb_z, lh_402, lh_425, lh_426, \
                         li_569, li_570, mh_549, mh_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_26 * lh_402[k]
                   + pb_z[k] * mh_549[k];

        t_736[k] = f_9 * lh_425[k]
                   + pb_y[k] * mh_551[k];

        t_737[k] = pa_y[k] * li_569[k];

        t_738[k] = f_12 * lh_426[k]
                   + pa_y[k] * li_570[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, pa_y, pb_y, pb_z, lh_405, lh_428, lh_429, \
                         li_572, li_574, mh_552, mh_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_26 * lh_405[k]
                   + pb_z[k] * mh_552[k];

        t_740[k] = f_10 * lh_428[k]
                   + pa_y[k] * li_572[k];

        t_741[k] = f_9 * lh_429[k]
                   + pb_y[k] * mh_555[k];

        t_742[k] = pa_y[k] * li_574[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, t_747, pb_x, lh_561, lh_562, lh_563, \
                         lh_564, lh_565, mh_561, mh_562, mh_563, mh_564, \
                         mh_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = f_11 * lh_561[k]
                   + pb_x[k] * mh_561[k];

        t_744[k] = f_11 * lh_562[k]
                   + pb_x[k] * mh_562[k];

        t_745[k] = f_11 * lh_563[k]
                   + pb_x[k] * mh_563[k];

        t_746[k] = f_11 * lh_564[k]
                   + pb_x[k] * mh_564[k];

        t_747[k] = f_11 * lh_565[k]
                   + pb_x[k] * mh_565[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, t_752, pa_y, pb_z, lh_414, lh_435, \
                         lh_437, lh_438, li_580, li_581, li_583, li_584, \
                         mh_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = pa_y[k] * li_580[k];

        t_749[k] = f_14 * lh_435[k]
                   + pa_y[k] * li_581[k];

        t_750[k] = f_26 * lh_414[k]
                   + pb_z[k] * mh_561[k];

        t_751[k] = f_12 * lh_437[k]
                   + pa_y[k] * li_583[k];

        t_752[k] = f_11 * lh_438[k]
                   + pa_y[k] * li_584[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, t_756, pa_y, pa_z, pb_y, ki0_392, ki1_392, \
                         lh_439, lh_440, li_560, li_585, li_587, \
                         mh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_10 * lh_439[k]
                   + pa_y[k] * li_585[k];

        t_754[k] = f_9 * lh_440[k]
                   + pb_y[k] * mh_566[k];

        t_755[k] = pa_y[k] * li_587[k];

        t_756[k] = f_22 * ki0_392[k]
                   - f_23 * ki1_392[k]
                   + pa_z[k] * li_560[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, t_760, pb_y, pb_z, lh_420, mg0_405, mg1_405, \
                         mh_567, mh_568, mh_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = pb_y[k] * mh_567[k];

        t_758[k] = f_14 * lh_420[k]
                   + pb_z[k] * mh_567[k];

        t_759[k] = f_3 * mg0_405[k]
                   - f_4 * mg1_405[k]
                   + pb_y[k] * mh_568[k];

        t_760[k] = pb_y[k] * mh_569[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, t_764, pb_x, pb_y, pb_z, lh_423, lh_572, \
                         mg0_406, mg0_410, mg1_406, mg1_410, mh_570, \
                         mh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_11 * lh_572[k]
                   + f_7 * mg0_410[k]
                   - f_8 * mg1_410[k]
                   + pb_x[k] * mh_572[k];

        t_762[k] = f_5 * mg0_406[k]
                   - f_6 * mg1_406[k]
                   + pb_y[k] * mh_570[k];

        t_763[k] = f_14 * lh_423[k]
                   + pb_z[k] * mh_570[k];

        t_764[k] = pb_y[k] * mh_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pb_x, pb_y, pb_z, lh_426, lh_576, mg0_408, \
                         mg0_414, mg1_408, mg1_414, mh_573, mh_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_11 * lh_576[k]
                   + f_5 * mg0_414[k]
                   - f_6 * mg1_414[k]
                   + pb_x[k] * mh_576[k];

        t_766[k] = f_7 * mg0_408[k]
                   - f_8 * mg1_408[k]
                   + pb_y[k] * mh_573[k];

        t_767[k] = f_14 * lh_426[k]
                   + pb_z[k] * mh_573[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pb_x, pb_y, lh_581, lh_582, mg0_410, \
                         mg0_419, mg1_410, mg1_419, mh_575, mh_576, mh_581, \
                         mh_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_3 * mg0_410[k]
                   - f_4 * mg1_410[k]
                   + pb_y[k] * mh_575[k];

        t_769[k] = pb_y[k] * mh_576[k];

        t_770[k] = f_11 * lh_581[k]
                   + f_3 * mg0_419[k]
                   - f_4 * mg1_419[k]
                   + pb_x[k] * mh_581[k];

        t_771[k] = f_11 * lh_582[k]
                   + pb_x[k] * mh_582[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, t_776, pb_x, pb_y, lh_583, lh_584, \
                         lh_585, lh_587, mh_581, mh_583, mh_584, mh_585, \
                         mh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_11 * lh_583[k]
                   + pb_x[k] * mh_583[k];

        t_773[k] = f_11 * lh_584[k]
                   + pb_x[k] * mh_584[k];

        t_774[k] = f_11 * lh_585[k]
                   + pb_x[k] * mh_585[k];

        t_775[k] = pb_y[k] * mh_581[k];

        t_776[k] = f_11 * lh_587[k]
                   + pb_x[k] * mh_587[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, t_780, pb_y, pb_z, lh_435, mg0_415, mg0_417, \
                         mg0_418, mg1_415, mg1_417, mg1_418, mh_582, mh_584, \
                         mh_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_1 * mg0_415[k]
                   - f_2 * mg1_415[k]
                   + pb_y[k] * mh_582[k];

        t_778[k] = f_14 * lh_435[k]
                   + pb_z[k] * mh_582[k];

        t_779[k] = f_7 * mg0_417[k]
                   - f_8 * mg1_417[k]
                   + pb_y[k] * mh_584[k];

        t_780[k] = f_5 * mg0_418[k]
                   - f_6 * mg1_418[k]
                   + pb_y[k] * mh_585[k];
    }

#pragma omp simd aligned(t_781, t_782, t_783, pa_x, pb_y, ki0_783, ki1_783, li_783, mg0_419, \
                         mg1_419, mh_586, mh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_781[k] = f_3 * mg0_419[k]
                   - f_4 * mg1_419[k]
                   + pb_y[k] * mh_586[k];

        t_782[k] = pb_y[k] * mh_587[k];

        t_783[k] = f_20 * ki0_783[k]
                   - f_21 * ki1_783[k]
                   + pa_x[k] * li_783[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, pa_y, pb_y, pb_z, ki0_420, ki1_420, lh_441, \
                         li_588, mh_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_18 * ki0_420[k]
                   - f_19 * ki1_420[k]
                   + pa_y[k] * li_588[k];

        t_785[k] = f_17 * lh_441[k]
                   + pb_y[k] * mh_588[k];

        t_786[k] = pb_z[k] * mh_588[k];
    }
}

static auto
compute_prim_mi_electron_repulsion_0_piece6(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t ki0,
                                            const size_t ki1, const size_t lh, const size_t li,
                                            const size_t mg0, const size_t mg1, const size_t mh,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_14 = 3.0 / p;
    const auto f_15 = 0.5 / alpha;
    const auto f_16 = 0.5 * beta / (alpha * p);
    const auto f_17 = 3.5 / p;
    const auto f_20 = 1.0 / alpha;
    const auto f_21 = beta / (alpha * p);
    const auto f_24 = 1.5 / alpha;
    const auto f_25 = 1.5 * beta / (alpha * p);
    const auto f_26 = 2.5 / p;
    const auto f_27 = 2.0 / alpha;
    const auto f_28 = 2.0 * beta / (alpha * p);

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

    const auto *ki0_423 = buffer.data(ki0 + 423);
    const auto *ki0_426 = buffer.data(ki0 + 426);
    const auto *ki0_430 = buffer.data(ki0 + 430);
    const auto *ki0_451 = buffer.data(ki0 + 451);
    const auto *ki0_454 = buffer.data(ki0 + 454);
    const auto *ki0_458 = buffer.data(ki0 + 458);
    const auto *ki0_476 = buffer.data(ki0 + 476);
    const auto *ki0_479 = buffer.data(ki0 + 479);
    const auto *ki0_481 = buffer.data(ki0 + 481);
    const auto *ki0_482 = buffer.data(ki0 + 482);
    const auto *ki0_485 = buffer.data(ki0 + 485);
    const auto *ki0_486 = buffer.data(ki0 + 486);
    const auto *ki0_490 = buffer.data(ki0 + 490);
    const auto *ki0_504 = buffer.data(ki0 + 504);
    const auto *ki0_509 = buffer.data(ki0 + 509);
    const auto *ki0_513 = buffer.data(ki0 + 513);
    const auto *ki0_518 = buffer.data(ki0 + 518);
    const auto *ki0_532 = buffer.data(ki0 + 532);
    const auto *ki0_537 = buffer.data(ki0 + 537);
    const auto *ki0_541 = buffer.data(ki0 + 541);
    const auto *ki0_805 = buffer.data(ki0 + 805);
    const auto *ki0_861 = buffer.data(ki0 + 861);
    const auto *ki0_863 = buffer.data(ki0 + 863);
    const auto *ki0_864 = buffer.data(ki0 + 864);
    const auto *ki0_865 = buffer.data(ki0 + 865);
    const auto *ki0_867 = buffer.data(ki0 + 867);
    const auto *ki0_889 = buffer.data(ki0 + 889);
    const auto *ki0_891 = buffer.data(ki0 + 891);
    const auto *ki0_892 = buffer.data(ki0 + 892);
    const auto *ki0_893 = buffer.data(ki0 + 893);
    const auto *ki0_895 = buffer.data(ki0 + 895);

    const auto *ki1_423 = buffer.data(ki1 + 423);
    const auto *ki1_426 = buffer.data(ki1 + 426);
    const auto *ki1_430 = buffer.data(ki1 + 430);
    const auto *ki1_451 = buffer.data(ki1 + 451);
    const auto *ki1_454 = buffer.data(ki1 + 454);
    const auto *ki1_458 = buffer.data(ki1 + 458);
    const auto *ki1_476 = buffer.data(ki1 + 476);
    const auto *ki1_479 = buffer.data(ki1 + 479);
    const auto *ki1_481 = buffer.data(ki1 + 481);
    const auto *ki1_482 = buffer.data(ki1 + 482);
    const auto *ki1_485 = buffer.data(ki1 + 485);
    const auto *ki1_486 = buffer.data(ki1 + 486);
    const auto *ki1_490 = buffer.data(ki1 + 490);
    const auto *ki1_504 = buffer.data(ki1 + 504);
    const auto *ki1_509 = buffer.data(ki1 + 509);
    const auto *ki1_513 = buffer.data(ki1 + 513);
    const auto *ki1_518 = buffer.data(ki1 + 518);
    const auto *ki1_532 = buffer.data(ki1 + 532);
    const auto *ki1_537 = buffer.data(ki1 + 537);
    const auto *ki1_541 = buffer.data(ki1 + 541);
    const auto *ki1_805 = buffer.data(ki1 + 805);
    const auto *ki1_861 = buffer.data(ki1 + 861);
    const auto *ki1_863 = buffer.data(ki1 + 863);
    const auto *ki1_864 = buffer.data(ki1 + 864);
    const auto *ki1_865 = buffer.data(ki1 + 865);
    const auto *ki1_867 = buffer.data(ki1 + 867);
    const auto *ki1_889 = buffer.data(ki1 + 889);
    const auto *ki1_891 = buffer.data(ki1 + 891);
    const auto *ki1_892 = buffer.data(ki1 + 892);
    const auto *ki1_893 = buffer.data(ki1 + 893);
    const auto *ki1_895 = buffer.data(ki1 + 895);

    const auto *lh_441 = buffer.data(lh + 441);
    const auto *lh_443 = buffer.data(lh + 443);
    const auto *lh_444 = buffer.data(lh + 444);
    const auto *lh_446 = buffer.data(lh + 446);
    const auto *lh_447 = buffer.data(lh + 447);
    const auto *lh_448 = buffer.data(lh + 448);
    const auto *lh_450 = buffer.data(lh + 450);
    const auto *lh_456 = buffer.data(lh + 456);
    const auto *lh_457 = buffer.data(lh + 457);
    const auto *lh_458 = buffer.data(lh + 458);
    const auto *lh_459 = buffer.data(lh + 459);
    const auto *lh_461 = buffer.data(lh + 461);
    const auto *lh_462 = buffer.data(lh + 462);
    const auto *lh_464 = buffer.data(lh + 464);
    const auto *lh_465 = buffer.data(lh + 465);
    const auto *lh_467 = buffer.data(lh + 467);
    const auto *lh_468 = buffer.data(lh + 468);
    const auto *lh_471 = buffer.data(lh + 471);
    const auto *lh_477 = buffer.data(lh + 477);
    const auto *lh_482 = buffer.data(lh + 482);
    const auto *lh_483 = buffer.data(lh + 483);
    const auto *lh_485 = buffer.data(lh + 485);
    const auto *lh_486 = buffer.data(lh + 486);
    const auto *lh_488 = buffer.data(lh + 488);
    const auto *lh_489 = buffer.data(lh + 489);
    const auto *lh_492 = buffer.data(lh + 492);
    const auto *lh_498 = buffer.data(lh + 498);
    const auto *lh_503 = buffer.data(lh + 503);
    const auto *lh_504 = buffer.data(lh + 504);
    const auto *lh_506 = buffer.data(lh + 506);
    const auto *lh_507 = buffer.data(lh + 507);
    const auto *lh_509 = buffer.data(lh + 509);
    const auto *lh_510 = buffer.data(lh + 510);
    const auto *lh_513 = buffer.data(lh + 513);
    const auto *lh_524 = buffer.data(lh + 524);
    const auto *lh_525 = buffer.data(lh + 525);
    const auto *lh_527 = buffer.data(lh + 527);
    const auto *lh_530 = buffer.data(lh + 530);
    const auto *lh_591 = buffer.data(lh + 591);
    const auto *lh_594 = buffer.data(lh + 594);
    const auto *lh_598 = buffer.data(lh + 598);
    const auto *lh_603 = buffer.data(lh + 603);
    const auto *lh_605 = buffer.data(lh + 605);
    const auto *lh_606 = buffer.data(lh + 606);
    const auto *lh_607 = buffer.data(lh + 607);
    const auto *lh_608 = buffer.data(lh + 608);
    const auto *lh_625 = buffer.data(lh + 625);
    const auto *lh_626 = buffer.data(lh + 626);
    const auto *lh_627 = buffer.data(lh + 627);
    const auto *lh_628 = buffer.data(lh + 628);
    const auto *lh_629 = buffer.data(lh + 629);
    const auto *lh_642 = buffer.data(lh + 642);
    const auto *lh_645 = buffer.data(lh + 645);
    const auto *lh_646 = buffer.data(lh + 646);
    const auto *lh_647 = buffer.data(lh + 647);
    const auto *lh_648 = buffer.data(lh + 648);
    const auto *lh_649 = buffer.data(lh + 649);
    const auto *lh_650 = buffer.data(lh + 650);
    const auto *lh_663 = buffer.data(lh + 663);
    const auto *lh_666 = buffer.data(lh + 666);
    const auto *lh_667 = buffer.data(lh + 667);
    const auto *lh_668 = buffer.data(lh + 668);
    const auto *lh_669 = buffer.data(lh + 669);
    const auto *lh_670 = buffer.data(lh + 670);
    const auto *lh_671 = buffer.data(lh + 671);

    const auto *li_588 = buffer.data(li + 588);
    const auto *li_589 = buffer.data(li + 589);
    const auto *li_591 = buffer.data(li + 591);
    const auto *li_593 = buffer.data(li + 593);
    const auto *li_594 = buffer.data(li + 594);
    const auto *li_597 = buffer.data(li + 597);
    const auto *li_598 = buffer.data(li + 598);
    const auto *li_600 = buffer.data(li + 600);
    const auto *li_602 = buffer.data(li + 602);
    const auto *li_603 = buffer.data(li + 603);
    const auto *li_609 = buffer.data(li + 609);
    const auto *li_611 = buffer.data(li + 611);
    const auto *li_612 = buffer.data(li + 612);
    const auto *li_613 = buffer.data(li + 613);
    const auto *li_615 = buffer.data(li + 615);
    const auto *li_619 = buffer.data(li + 619);
    const auto *li_622 = buffer.data(li + 622);
    const auto *li_626 = buffer.data(li + 626);
    const auto *li_644 = buffer.data(li + 644);
    const auto *li_647 = buffer.data(li + 647);
    const auto *li_649 = buffer.data(li + 649);
    const auto *li_650 = buffer.data(li + 650);
    const auto *li_653 = buffer.data(li + 653);
    const auto *li_654 = buffer.data(li + 654);
    const auto *li_658 = buffer.data(li + 658);
    const auto *li_672 = buffer.data(li + 672);
    const auto *li_675 = buffer.data(li + 675);
    const auto *li_677 = buffer.data(li + 677);
    const auto *li_678 = buffer.data(li + 678);
    const auto *li_681 = buffer.data(li + 681);
    const auto *li_682 = buffer.data(li + 682);
    const auto *li_686 = buffer.data(li + 686);
    const auto *li_700 = buffer.data(li + 700);
    const auto *li_705 = buffer.data(li + 705);
    const auto *li_709 = buffer.data(li + 709);
    const auto *li_805 = buffer.data(li + 805);
    const auto *li_861 = buffer.data(li + 861);
    const auto *li_863 = buffer.data(li + 863);
    const auto *li_864 = buffer.data(li + 864);
    const auto *li_865 = buffer.data(li + 865);
    const auto *li_867 = buffer.data(li + 867);
    const auto *li_889 = buffer.data(li + 889);
    const auto *li_891 = buffer.data(li + 891);
    const auto *li_892 = buffer.data(li + 892);
    const auto *li_893 = buffer.data(li + 893);
    const auto *li_895 = buffer.data(li + 895);

    const auto *mg0_420 = buffer.data(mg0 + 420);
    const auto *mg0_422 = buffer.data(mg0 + 422);
    const auto *mg0_423 = buffer.data(mg0 + 423);
    const auto *mg0_425 = buffer.data(mg0 + 425);
    const auto *mg0_426 = buffer.data(mg0 + 426);
    const auto *mg0_430 = buffer.data(mg0 + 430);
    const auto *mg0_431 = buffer.data(mg0 + 431);
    const auto *mg0_432 = buffer.data(mg0 + 432);
    const auto *mg0_434 = buffer.data(mg0 + 434);
    const auto *mg0_462 = buffer.data(mg0 + 462);
    const auto *mg0_477 = buffer.data(mg0 + 477);

    const auto *mg1_420 = buffer.data(mg1 + 420);
    const auto *mg1_422 = buffer.data(mg1 + 422);
    const auto *mg1_423 = buffer.data(mg1 + 423);
    const auto *mg1_425 = buffer.data(mg1 + 425);
    const auto *mg1_426 = buffer.data(mg1 + 426);
    const auto *mg1_430 = buffer.data(mg1 + 430);
    const auto *mg1_431 = buffer.data(mg1 + 431);
    const auto *mg1_432 = buffer.data(mg1 + 432);
    const auto *mg1_434 = buffer.data(mg1 + 434);
    const auto *mg1_462 = buffer.data(mg1 + 462);
    const auto *mg1_477 = buffer.data(mg1 + 477);

    const auto *mh_589 = buffer.data(mh + 589);
    const auto *mh_590 = buffer.data(mh + 590);
    const auto *mh_591 = buffer.data(mh + 591);
    const auto *mh_593 = buffer.data(mh + 593);
    const auto *mh_594 = buffer.data(mh + 594);
    const auto *mh_595 = buffer.data(mh + 595);
    const auto *mh_597 = buffer.data(mh + 597);
    const auto *mh_598 = buffer.data(mh + 598);
    const auto *mh_603 = buffer.data(mh + 603);
    const auto *mh_604 = buffer.data(mh + 604);
    const auto *mh_605 = buffer.data(mh + 605);
    const auto *mh_606 = buffer.data(mh + 606);
    const auto *mh_607 = buffer.data(mh + 607);
    const auto *mh_608 = buffer.data(mh + 608);
    const auto *mh_609 = buffer.data(mh + 609);
    const auto *mh_611 = buffer.data(mh + 611);
    const auto *mh_612 = buffer.data(mh + 612);
    const auto *mh_614 = buffer.data(mh + 614);
    const auto *mh_615 = buffer.data(mh + 615);
    const auto *mh_618 = buffer.data(mh + 618);
    const auto *mh_624 = buffer.data(mh + 624);
    const auto *mh_625 = buffer.data(mh + 625);
    const auto *mh_626 = buffer.data(mh + 626);
    const auto *mh_627 = buffer.data(mh + 627);
    const auto *mh_628 = buffer.data(mh + 628);
    const auto *mh_629 = buffer.data(mh + 629);
    const auto *mh_630 = buffer.data(mh + 630);
    const auto *mh_632 = buffer.data(mh + 632);
    const auto *mh_633 = buffer.data(mh + 633);
    const auto *mh_635 = buffer.data(mh + 635);
    const auto *mh_636 = buffer.data(mh + 636);
    const auto *mh_639 = buffer.data(mh + 639);
    const auto *mh_642 = buffer.data(mh + 642);
    const auto *mh_645 = buffer.data(mh + 645);
    const auto *mh_646 = buffer.data(mh + 646);
    const auto *mh_647 = buffer.data(mh + 647);
    const auto *mh_648 = buffer.data(mh + 648);
    const auto *mh_649 = buffer.data(mh + 649);
    const auto *mh_650 = buffer.data(mh + 650);
    const auto *mh_651 = buffer.data(mh + 651);
    const auto *mh_653 = buffer.data(mh + 653);
    const auto *mh_654 = buffer.data(mh + 654);
    const auto *mh_656 = buffer.data(mh + 656);
    const auto *mh_657 = buffer.data(mh + 657);
    const auto *mh_660 = buffer.data(mh + 660);
    const auto *mh_663 = buffer.data(mh + 663);
    const auto *mh_666 = buffer.data(mh + 666);
    const auto *mh_667 = buffer.data(mh + 667);
    const auto *mh_668 = buffer.data(mh + 668);
    const auto *mh_669 = buffer.data(mh + 669);
    const auto *mh_670 = buffer.data(mh + 670);
    const auto *mh_671 = buffer.data(mh + 671);
    const auto *mh_672 = buffer.data(mh + 672);
    const auto *mh_674 = buffer.data(mh + 674);
    const auto *mh_675 = buffer.data(mh + 675);
    const auto *mh_677 = buffer.data(mh + 677);
    const auto *mh_678 = buffer.data(mh + 678);

#pragma omp simd aligned(t_787, t_788, t_789, pb_x, pb_z, lh_591, mg0_420, mg0_423, mg1_420, \
                         mg1_423, mh_589, mh_590, mh_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_10 * lh_591[k]
                   + f_7 * mg0_423[k]
                   - f_8 * mg1_423[k]
                   + pb_x[k] * mh_591[k];

        t_788[k] = pb_z[k] * mh_589[k];

        t_789[k] = f_3 * mg0_420[k]
                   - f_4 * mg1_420[k]
                   + pb_z[k] * mh_590[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pb_x, pb_y, pb_z, lh_446, lh_594, \
                         mg0_422, mg0_426, mg1_422, mg1_426, mh_591, mh_593, \
                         mh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_10 * lh_594[k]
                   + f_5 * mg0_426[k]
                   - f_6 * mg1_426[k]
                   + pb_x[k] * mh_594[k];

        t_791[k] = pb_z[k] * mh_591[k];

        t_792[k] = f_17 * lh_446[k]
                   + pb_y[k] * mh_593[k];

        t_793[k] = f_5 * mg0_422[k]
                   - f_6 * mg1_422[k]
                   + pb_z[k] * mh_593[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pb_x, pb_z, lh_598, mg0_423, mg0_430, mg1_423, \
                         mg1_430, mh_594, mh_595, mh_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_10 * lh_598[k]
                   + f_3 * mg0_430[k]
                   - f_4 * mg1_430[k]
                   + pb_x[k] * mh_598[k];

        t_795[k] = pb_z[k] * mh_594[k];

        t_796[k] = f_3 * mg0_423[k]
                   - f_4 * mg1_423[k]
                   + pb_z[k] * mh_595[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, t_800, pb_x, pb_y, pb_z, lh_450, lh_603, \
                         mg0_425, mg1_425, mh_597, mh_598, mh_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_17 * lh_450[k]
                   + pb_y[k] * mh_597[k];

        t_798[k] = f_7 * mg0_425[k]
                   - f_8 * mg1_425[k]
                   + pb_z[k] * mh_597[k];

        t_799[k] = f_10 * lh_603[k]
                   + pb_x[k] * mh_603[k];

        t_800[k] = pb_z[k] * mh_598[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, t_804, pb_x, lh_605, lh_606, lh_607, lh_608, \
                         mh_605, mh_606, mh_607, mh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = f_10 * lh_605[k]
                   + pb_x[k] * mh_605[k];

        t_802[k] = f_10 * lh_606[k]
                   + pb_x[k] * mh_606[k];

        t_803[k] = f_10 * lh_607[k]
                   + pb_x[k] * mh_607[k];

        t_804[k] = f_10 * lh_608[k]
                   + pb_x[k] * mh_608[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, pa_x, pb_z, ki0_805, ki1_805, li_805, \
                         mg0_430, mg0_431, mg1_430, mg1_431, mh_603, mh_604, \
                         mh_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_15 * ki0_805[k]
                   - f_16 * ki1_805[k]
                   + pa_x[k] * li_805[k];

        t_806[k] = pb_z[k] * mh_603[k];

        t_807[k] = f_3 * mg0_430[k]
                   - f_4 * mg1_430[k]
                   + pb_z[k] * mh_604[k];

        t_808[k] = f_5 * mg0_431[k]
                   - f_6 * mg1_431[k]
                   + pb_z[k] * mh_605[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, t_812, pa_z, pb_y, pb_z, lh_461, li_588, \
                         mg0_432, mg0_434, mg1_432, mg1_434, mh_606, \
                         mh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_7 * mg0_432[k]
                   - f_8 * mg1_432[k]
                   + pb_z[k] * mh_606[k];

        t_810[k] = f_17 * lh_461[k]
                   + pb_y[k] * mh_608[k];

        t_811[k] = f_1 * mg0_434[k]
                   - f_2 * mg1_434[k]
                   + pb_z[k] * mh_608[k];

        t_812[k] = pa_z[k] * li_588[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, t_816, t_817, pa_z, pb_y, pb_z, lh_441, lh_443, \
                         lh_464, li_589, li_591, li_593, mh_609, \
                         mh_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = pa_z[k] * li_589[k];

        t_814[k] = f_9 * lh_441[k]
                   + pb_z[k] * mh_609[k];

        t_815[k] = pa_z[k] * li_591[k];

        t_816[k] = f_14 * lh_464[k]
                   + pb_y[k] * mh_611[k];

        t_817[k] = f_10 * lh_443[k]
                   + pa_z[k] * li_593[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, t_821, t_822, pa_z, pb_y, pb_z, lh_444, lh_446, \
                         lh_467, li_594, li_597, li_598, mh_612, \
                         mh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = pa_z[k] * li_594[k];

        t_819[k] = f_9 * lh_444[k]
                   + pb_z[k] * mh_612[k];

        t_820[k] = f_14 * lh_467[k]
                   + pb_y[k] * mh_614[k];

        t_821[k] = f_11 * lh_446[k]
                   + pa_z[k] * li_597[k];

        t_822[k] = pa_z[k] * li_598[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, t_826, pa_z, pb_y, pb_z, lh_447, lh_448, lh_450, \
                         lh_471, li_600, li_602, mh_615, mh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_9 * lh_447[k]
                   + pb_z[k] * mh_615[k];

        t_824[k] = f_10 * lh_448[k]
                   + pa_z[k] * li_600[k];

        t_825[k] = f_14 * lh_471[k]
                   + pb_y[k] * mh_618[k];

        t_826[k] = f_12 * lh_450[k]
                   + pa_z[k] * li_602[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, t_830, t_831, pa_z, pb_x, lh_625, lh_626, \
                         lh_627, lh_628, li_603, mh_625, mh_626, mh_627, \
                         mh_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = pa_z[k] * li_603[k];

        t_828[k] = f_10 * lh_625[k]
                   + pb_x[k] * mh_625[k];

        t_829[k] = f_10 * lh_626[k]
                   + pb_x[k] * mh_626[k];

        t_830[k] = f_10 * lh_627[k]
                   + pb_x[k] * mh_627[k];

        t_831[k] = f_10 * lh_628[k]
                   + pb_x[k] * mh_628[k];
    }

#pragma omp simd aligned(t_832, t_833, t_834, t_835, pa_z, pb_x, pb_z, lh_456, lh_457, lh_629, \
                         li_609, li_611, mh_624, mh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_832[k] = f_10 * lh_629[k]
                   + pb_x[k] * mh_629[k];

        t_833[k] = pa_z[k] * li_609[k];

        t_834[k] = f_9 * lh_456[k]
                   + pb_z[k] * mh_624[k];

        t_835[k] = f_10 * lh_457[k]
                   + pa_z[k] * li_611[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, t_839, pa_z, pb_y, lh_458, lh_459, lh_461, \
                         lh_482, li_612, li_613, li_615, mh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_11 * lh_458[k]
                   + pa_z[k] * li_612[k];

        t_837[k] = f_12 * lh_459[k]
                   + pa_z[k] * li_613[k];

        t_838[k] = f_14 * lh_482[k]
                   + pb_y[k] * mh_629[k];

        t_839[k] = f_14 * lh_461[k]
                   + pa_z[k] * li_615[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pa_y, pb_y, pb_z, ki0_476, ki1_476, lh_462, \
                         lh_483, li_644, mh_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_27 * ki0_476[k]
                   - f_28 * ki1_476[k]
                   + pa_y[k] * li_644[k];

        t_841[k] = f_26 * lh_483[k]
                   + pb_y[k] * mh_630[k];

        t_842[k] = f_10 * lh_462[k]
                   + pb_z[k] * mh_630[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pa_y, pa_z, pb_y, ki0_423, ki0_481, ki1_423, \
                         ki1_481, lh_485, li_619, li_649, mh_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = f_15 * ki0_423[k]
                   - f_16 * ki1_423[k]
                   + pa_z[k] * li_619[k];

        t_844[k] = f_26 * lh_485[k]
                   + pb_y[k] * mh_632[k];

        t_845[k] = f_27 * ki0_481[k]
                   - f_28 * ki1_481[k]
                   + pa_y[k] * li_649[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, pa_z, pb_y, pb_z, ki0_426, ki1_426, lh_465, \
                         lh_488, li_622, mh_633, mh_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = f_15 * ki0_426[k]
                   - f_16 * ki1_426[k]
                   + pa_z[k] * li_622[k];

        t_847[k] = f_10 * lh_465[k]
                   + pb_z[k] * mh_633[k];

        t_848[k] = f_26 * lh_488[k]
                   + pb_y[k] * mh_635[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, pa_y, pa_z, pb_z, ki0_430, ki0_485, ki1_430, \
                         ki1_485, lh_468, li_626, li_653, mh_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_27 * ki0_485[k]
                   - f_28 * ki1_485[k]
                   + pa_y[k] * li_653[k];

        t_850[k] = f_15 * ki0_430[k]
                   - f_16 * ki1_430[k]
                   + pa_z[k] * li_626[k];

        t_851[k] = f_10 * lh_468[k]
                   + pb_z[k] * mh_636[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, pa_y, pb_x, pb_y, ki0_490, ki1_490, lh_492, \
                         lh_642, li_658, mg0_462, mg1_462, mh_639, \
                         mh_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_10 * lh_642[k]
                   + f_3 * mg0_462[k]
                   - f_4 * mg1_462[k]
                   + pb_x[k] * mh_642[k];

        t_853[k] = f_26 * lh_492[k]
                   + pb_y[k] * mh_639[k];

        t_854[k] = f_27 * ki0_490[k]
                   - f_28 * ki1_490[k]
                   + pa_y[k] * li_658[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, pb_x, lh_645, lh_646, lh_647, \
                         lh_648, lh_649, mh_645, mh_646, mh_647, mh_648, \
                         mh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = f_10 * lh_645[k]
                   + pb_x[k] * mh_645[k];

        t_856[k] = f_10 * lh_646[k]
                   + pb_x[k] * mh_646[k];

        t_857[k] = f_10 * lh_647[k]
                   + pb_x[k] * mh_647[k];

        t_858[k] = f_10 * lh_648[k]
                   + pb_x[k] * mh_648[k];

        t_859[k] = f_10 * lh_649[k]
                   + pb_x[k] * mh_649[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, pa_x, pb_x, pb_z, ki0_861, ki1_861, lh_477, \
                         lh_650, li_861, mh_645, mh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = f_10 * lh_650[k]
                   + pb_x[k] * mh_650[k];

        t_861[k] = f_15 * ki0_861[k]
                   - f_16 * ki1_861[k]
                   + pa_x[k] * li_861[k];

        t_862[k] = f_10 * lh_477[k]
                   + pb_z[k] * mh_645[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, pa_x, ki0_863, ki0_864, ki0_865, ki1_863, \
                         ki1_864, ki1_865, li_863, li_864, li_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = f_15 * ki0_863[k]
                   - f_16 * ki1_863[k]
                   + pa_x[k] * li_863[k];

        t_864[k] = f_15 * ki0_864[k]
                   - f_16 * ki1_864[k]
                   + pa_x[k] * li_864[k];

        t_865[k] = f_15 * ki0_865[k]
                   - f_16 * ki1_865[k]
                   + pa_x[k] * li_865[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pa_x, pa_y, pb_y, ki0_504, ki0_867, ki1_504, \
                         ki1_867, lh_503, li_672, li_867, mh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_26 * lh_503[k]
                   + pb_y[k] * mh_650[k];

        t_867[k] = f_15 * ki0_867[k]
                   - f_16 * ki1_867[k]
                   + pa_x[k] * li_867[k];

        t_868[k] = f_24 * ki0_504[k]
                   - f_25 * ki1_504[k]
                   + pa_y[k] * li_672[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pa_z, pb_y, pb_z, ki0_451, ki1_451, \
                         lh_483, lh_504, lh_506, li_647, mh_651, \
                         mh_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_12 * lh_504[k]
                   + pb_y[k] * mh_651[k];

        t_870[k] = f_11 * lh_483[k]
                   + pb_z[k] * mh_651[k];

        t_871[k] = f_20 * ki0_451[k]
                   - f_21 * ki1_451[k]
                   + pa_z[k] * li_647[k];

        t_872[k] = f_12 * lh_506[k]
                   + pb_y[k] * mh_653[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pa_y, pa_z, pb_z, ki0_454, ki0_509, ki1_454, \
                         ki1_509, lh_486, li_650, li_677, mh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_24 * ki0_509[k]
                   - f_25 * ki1_509[k]
                   + pa_y[k] * li_677[k];

        t_874[k] = f_20 * ki0_454[k]
                   - f_21 * ki1_454[k]
                   + pa_z[k] * li_650[k];

        t_875[k] = f_11 * lh_486[k]
                   + pb_z[k] * mh_654[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pa_y, pa_z, pb_y, ki0_458, ki0_513, ki1_458, \
                         ki1_513, lh_509, li_654, li_681, mh_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_12 * lh_509[k]
                   + pb_y[k] * mh_656[k];

        t_877[k] = f_24 * ki0_513[k]
                   - f_25 * ki1_513[k]
                   + pa_y[k] * li_681[k];

        t_878[k] = f_20 * ki0_458[k]
                   - f_21 * ki1_458[k]
                   + pa_z[k] * li_654[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pb_x, pb_y, pb_z, lh_489, lh_513, lh_663, \
                         mg0_477, mg1_477, mh_657, mh_660, mh_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_11 * lh_489[k]
                   + pb_z[k] * mh_657[k];

        t_880[k] = f_10 * lh_663[k]
                   + f_3 * mg0_477[k]
                   - f_4 * mg1_477[k]
                   + pb_x[k] * mh_663[k];

        t_881[k] = f_12 * lh_513[k]
                   + pb_y[k] * mh_660[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pa_y, pb_x, ki0_518, ki1_518, lh_666, \
                         lh_667, lh_668, li_686, mh_666, mh_667, \
                         mh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_24 * ki0_518[k]
                   - f_25 * ki1_518[k]
                   + pa_y[k] * li_686[k];

        t_883[k] = f_10 * lh_666[k]
                   + pb_x[k] * mh_666[k];

        t_884[k] = f_10 * lh_667[k]
                   + pb_x[k] * mh_667[k];

        t_885[k] = f_10 * lh_668[k]
                   + pb_x[k] * mh_668[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, pa_x, pb_x, ki0_889, ki1_889, lh_669, \
                         lh_670, lh_671, li_889, mh_669, mh_670, \
                         mh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_10 * lh_669[k]
                   + pb_x[k] * mh_669[k];

        t_887[k] = f_10 * lh_670[k]
                   + pb_x[k] * mh_670[k];

        t_888[k] = f_10 * lh_671[k]
                   + pb_x[k] * mh_671[k];

        t_889[k] = f_15 * ki0_889[k]
                   - f_16 * ki1_889[k]
                   + pa_x[k] * li_889[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, pa_x, pb_z, ki0_891, ki0_892, ki1_891, ki1_892, \
                         lh_498, li_891, li_892, mh_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_11 * lh_498[k]
                   + pb_z[k] * mh_666[k];

        t_891[k] = f_15 * ki0_891[k]
                   - f_16 * ki1_891[k]
                   + pa_x[k] * li_891[k];

        t_892[k] = f_15 * ki0_892[k]
                   - f_16 * ki1_892[k]
                   + pa_x[k] * li_892[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pa_x, pb_y, ki0_893, ki0_895, ki1_893, ki1_895, \
                         lh_524, li_893, li_895, mh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_15 * ki0_893[k]
                   - f_16 * ki1_893[k]
                   + pa_x[k] * li_893[k];

        t_894[k] = f_12 * lh_524[k]
                   + pb_y[k] * mh_671[k];

        t_895[k] = f_15 * ki0_895[k]
                   - f_16 * ki1_895[k]
                   + pa_x[k] * li_895[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pa_y, pb_y, pb_z, ki0_532, ki1_532, lh_504, \
                         lh_525, li_700, mh_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_20 * ki0_532[k]
                   - f_21 * ki1_532[k]
                   + pa_y[k] * li_700[k];

        t_897[k] = f_11 * lh_525[k]
                   + pb_y[k] * mh_672[k];

        t_898[k] = f_12 * lh_504[k]
                   + pb_z[k] * mh_672[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, pa_y, pa_z, pb_y, ki0_479, ki0_537, ki1_479, \
                         ki1_537, lh_527, li_675, li_705, mh_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = f_24 * ki0_479[k]
                   - f_25 * ki1_479[k]
                   + pa_z[k] * li_675[k];

        t_900[k] = f_11 * lh_527[k]
                   + pb_y[k] * mh_674[k];

        t_901[k] = f_20 * ki0_537[k]
                   - f_21 * ki1_537[k]
                   + pa_y[k] * li_705[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pa_z, pb_y, pb_z, ki0_482, ki1_482, lh_507, \
                         lh_530, li_678, mh_675, mh_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_24 * ki0_482[k]
                   - f_25 * ki1_482[k]
                   + pa_z[k] * li_678[k];

        t_903[k] = f_12 * lh_507[k]
                   + pb_z[k] * mh_675[k];

        t_904[k] = f_11 * lh_530[k]
                   + pb_y[k] * mh_677[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pa_y, pa_z, pb_z, ki0_486, ki0_541, ki1_486, \
                         ki1_541, lh_510, li_682, li_709, mh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_20 * ki0_541[k]
                   - f_21 * ki1_541[k]
                   + pa_y[k] * li_709[k];

        t_906[k] = f_24 * ki0_486[k]
                   - f_25 * ki1_486[k]
                   + pa_z[k] * li_682[k];

        t_907[k] = f_12 * lh_510[k]
                   + pb_z[k] * mh_678[k];
    }
}

static auto
compute_prim_mi_electron_repulsion_0_piece7(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t ki0,
                                            const size_t ki1, const size_t lh, const size_t li,
                                            const size_t mg0, const size_t mg1, const size_t mh,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_13 = 4.0 / p;
    const auto f_14 = 3.0 / p;
    const auto f_15 = 0.5 / alpha;
    const auto f_16 = 0.5 * beta / (alpha * p);
    const auto f_17 = 3.5 / p;
    const auto f_18 = 3.0 / alpha;
    const auto f_19 = 3.0 * beta / (alpha * p);
    const auto f_20 = 1.0 / alpha;
    const auto f_21 = beta / (alpha * p);
    const auto f_26 = 2.5 / p;
    const auto f_27 = 2.0 / alpha;
    const auto f_28 = 2.0 * beta / (alpha * p);

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
    auto *t_1034 = buffer.data(target + 1034);
    auto *t_1035 = buffer.data(target + 1035);
    auto *t_1036 = buffer.data(target + 1036);
    auto *t_1037 = buffer.data(target + 1037);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ki0_507 = buffer.data(ki0 + 507);
    const auto *ki0_510 = buffer.data(ki0 + 510);
    const auto *ki0_514 = buffer.data(ki0 + 514);
    const auto *ki0_546 = buffer.data(ki0 + 546);
    const auto *ki0_560 = buffer.data(ki0 + 560);
    const auto *ki0_565 = buffer.data(ki0 + 565);
    const auto *ki0_569 = buffer.data(ki0 + 569);
    const auto *ki0_574 = buffer.data(ki0 + 574);
    const auto *ki0_917 = buffer.data(ki0 + 917);
    const auto *ki0_919 = buffer.data(ki0 + 919);
    const auto *ki0_920 = buffer.data(ki0 + 920);
    const auto *ki0_921 = buffer.data(ki0 + 921);
    const auto *ki0_923 = buffer.data(ki0 + 923);
    const auto *ki0_945 = buffer.data(ki0 + 945);
    const auto *ki0_947 = buffer.data(ki0 + 947);
    const auto *ki0_948 = buffer.data(ki0 + 948);
    const auto *ki0_949 = buffer.data(ki0 + 949);
    const auto *ki0_951 = buffer.data(ki0 + 951);
    const auto *ki0_1007 = buffer.data(ki0 + 1007);

    const auto *ki1_507 = buffer.data(ki1 + 507);
    const auto *ki1_510 = buffer.data(ki1 + 510);
    const auto *ki1_514 = buffer.data(ki1 + 514);
    const auto *ki1_546 = buffer.data(ki1 + 546);
    const auto *ki1_560 = buffer.data(ki1 + 560);
    const auto *ki1_565 = buffer.data(ki1 + 565);
    const auto *ki1_569 = buffer.data(ki1 + 569);
    const auto *ki1_574 = buffer.data(ki1 + 574);
    const auto *ki1_917 = buffer.data(ki1 + 917);
    const auto *ki1_919 = buffer.data(ki1 + 919);
    const auto *ki1_920 = buffer.data(ki1 + 920);
    const auto *ki1_921 = buffer.data(ki1 + 921);
    const auto *ki1_923 = buffer.data(ki1 + 923);
    const auto *ki1_945 = buffer.data(ki1 + 945);
    const auto *ki1_947 = buffer.data(ki1 + 947);
    const auto *ki1_948 = buffer.data(ki1 + 948);
    const auto *ki1_949 = buffer.data(ki1 + 949);
    const auto *ki1_951 = buffer.data(ki1 + 951);
    const auto *ki1_1007 = buffer.data(ki1 + 1007);

    const auto *lh_519 = buffer.data(lh + 519);
    const auto *lh_525 = buffer.data(lh + 525);
    const auto *lh_528 = buffer.data(lh + 528);
    const auto *lh_531 = buffer.data(lh + 531);
    const auto *lh_534 = buffer.data(lh + 534);
    const auto *lh_540 = buffer.data(lh + 540);
    const auto *lh_545 = buffer.data(lh + 545);
    const auto *lh_546 = buffer.data(lh + 546);
    const auto *lh_548 = buffer.data(lh + 548);
    const auto *lh_549 = buffer.data(lh + 549);
    const auto *lh_551 = buffer.data(lh + 551);
    const auto *lh_552 = buffer.data(lh + 552);
    const auto *lh_555 = buffer.data(lh + 555);
    const auto *lh_561 = buffer.data(lh + 561);
    const auto *lh_566 = buffer.data(lh + 566);
    const auto *lh_567 = buffer.data(lh + 567);
    const auto *lh_568 = buffer.data(lh + 568);
    const auto *lh_569 = buffer.data(lh + 569);
    const auto *lh_570 = buffer.data(lh + 570);
    const auto *lh_572 = buffer.data(lh + 572);
    const auto *lh_573 = buffer.data(lh + 573);
    const auto *lh_575 = buffer.data(lh + 575);
    const auto *lh_576 = buffer.data(lh + 576);
    const auto *lh_582 = buffer.data(lh + 582);
    const auto *lh_584 = buffer.data(lh + 584);
    const auto *lh_585 = buffer.data(lh + 585);
    const auto *lh_586 = buffer.data(lh + 586);
    const auto *lh_587 = buffer.data(lh + 587);
    const auto *lh_588 = buffer.data(lh + 588);
    const auto *lh_593 = buffer.data(lh + 593);
    const auto *lh_597 = buffer.data(lh + 597);
    const auto *lh_684 = buffer.data(lh + 684);
    const auto *lh_687 = buffer.data(lh + 687);
    const auto *lh_688 = buffer.data(lh + 688);
    const auto *lh_689 = buffer.data(lh + 689);
    const auto *lh_690 = buffer.data(lh + 690);
    const auto *lh_691 = buffer.data(lh + 691);
    const auto *lh_692 = buffer.data(lh + 692);
    const auto *lh_705 = buffer.data(lh + 705);
    const auto *lh_708 = buffer.data(lh + 708);
    const auto *lh_709 = buffer.data(lh + 709);
    const auto *lh_710 = buffer.data(lh + 710);
    const auto *lh_711 = buffer.data(lh + 711);
    const auto *lh_712 = buffer.data(lh + 712);
    const auto *lh_713 = buffer.data(lh + 713);
    const auto *lh_729 = buffer.data(lh + 729);
    const auto *lh_730 = buffer.data(lh + 730);
    const auto *lh_731 = buffer.data(lh + 731);
    const auto *lh_732 = buffer.data(lh + 732);
    const auto *lh_733 = buffer.data(lh + 733);
    const auto *lh_740 = buffer.data(lh + 740);
    const auto *lh_744 = buffer.data(lh + 744);
    const auto *lh_749 = buffer.data(lh + 749);
    const auto *lh_750 = buffer.data(lh + 750);
    const auto *lh_751 = buffer.data(lh + 751);
    const auto *lh_752 = buffer.data(lh + 752);
    const auto *lh_753 = buffer.data(lh + 753);
    const auto *lh_755 = buffer.data(lh + 755);
    const auto *lh_756 = buffer.data(lh + 756);
    const auto *lh_759 = buffer.data(lh + 759);
    const auto *lh_761 = buffer.data(lh + 761);
    const auto *lh_762 = buffer.data(lh + 762);
    const auto *lh_765 = buffer.data(lh + 765);
    const auto *lh_766 = buffer.data(lh + 766);
    const auto *lh_768 = buffer.data(lh + 768);
    const auto *lh_770 = buffer.data(lh + 770);
    const auto *lh_771 = buffer.data(lh + 771);
    const auto *lh_773 = buffer.data(lh + 773);
    const auto *lh_774 = buffer.data(lh + 774);
    const auto *lh_775 = buffer.data(lh + 775);
    const auto *lh_776 = buffer.data(lh + 776);

    const auto *li_703 = buffer.data(li + 703);
    const auto *li_706 = buffer.data(li + 706);
    const auto *li_710 = buffer.data(li + 710);
    const auto *li_714 = buffer.data(li + 714);
    const auto *li_728 = buffer.data(li + 728);
    const auto *li_733 = buffer.data(li + 733);
    const auto *li_737 = buffer.data(li + 737);
    const auto *li_742 = buffer.data(li + 742);
    const auto *li_756 = buffer.data(li + 756);
    const auto *li_758 = buffer.data(li + 758);
    const auto *li_759 = buffer.data(li + 759);
    const auto *li_761 = buffer.data(li + 761);
    const auto *li_762 = buffer.data(li + 762);
    const auto *li_765 = buffer.data(li + 765);
    const auto *li_766 = buffer.data(li + 766);
    const auto *li_768 = buffer.data(li + 768);
    const auto *li_770 = buffer.data(li + 770);
    const auto *li_776 = buffer.data(li + 776);
    const auto *li_777 = buffer.data(li + 777);
    const auto *li_779 = buffer.data(li + 779);
    const auto *li_780 = buffer.data(li + 780);
    const auto *li_781 = buffer.data(li + 781);
    const auto *li_783 = buffer.data(li + 783);
    const auto *li_784 = buffer.data(li + 784);
    const auto *li_785 = buffer.data(li + 785);
    const auto *li_917 = buffer.data(li + 917);
    const auto *li_919 = buffer.data(li + 919);
    const auto *li_920 = buffer.data(li + 920);
    const auto *li_921 = buffer.data(li + 921);
    const auto *li_923 = buffer.data(li + 923);
    const auto *li_945 = buffer.data(li + 945);
    const auto *li_947 = buffer.data(li + 947);
    const auto *li_948 = buffer.data(li + 948);
    const auto *li_949 = buffer.data(li + 949);
    const auto *li_951 = buffer.data(li + 951);
    const auto *li_1007 = buffer.data(li + 1007);
    const auto *li_1008 = buffer.data(li + 1008);
    const auto *li_1011 = buffer.data(li + 1011);
    const auto *li_1013 = buffer.data(li + 1013);
    const auto *li_1014 = buffer.data(li + 1014);
    const auto *li_1017 = buffer.data(li + 1017);
    const auto *li_1018 = buffer.data(li + 1018);
    const auto *li_1020 = buffer.data(li + 1020);
    const auto *li_1022 = buffer.data(li + 1022);
    const auto *li_1029 = buffer.data(li + 1029);
    const auto *li_1031 = buffer.data(li + 1031);
    const auto *li_1032 = buffer.data(li + 1032);
    const auto *li_1033 = buffer.data(li + 1033);
    const auto *li_1034 = buffer.data(li + 1034);
    const auto *li_1035 = buffer.data(li + 1035);

    const auto *mg0_492 = buffer.data(mg0 + 492);
    const auto *mg0_507 = buffer.data(mg0 + 507);
    const auto *mg0_525 = buffer.data(mg0 + 525);
    const auto *mg0_526 = buffer.data(mg0 + 526);
    const auto *mg0_528 = buffer.data(mg0 + 528);
    const auto *mg0_530 = buffer.data(mg0 + 530);
    const auto *mg0_534 = buffer.data(mg0 + 534);
    const auto *mg0_535 = buffer.data(mg0 + 535);
    const auto *mg0_537 = buffer.data(mg0 + 537);
    const auto *mg0_538 = buffer.data(mg0 + 538);
    const auto *mg0_539 = buffer.data(mg0 + 539);

    const auto *mg1_492 = buffer.data(mg1 + 492);
    const auto *mg1_507 = buffer.data(mg1 + 507);
    const auto *mg1_525 = buffer.data(mg1 + 525);
    const auto *mg1_526 = buffer.data(mg1 + 526);
    const auto *mg1_528 = buffer.data(mg1 + 528);
    const auto *mg1_530 = buffer.data(mg1 + 530);
    const auto *mg1_534 = buffer.data(mg1 + 534);
    const auto *mg1_535 = buffer.data(mg1 + 535);
    const auto *mg1_537 = buffer.data(mg1 + 537);
    const auto *mg1_538 = buffer.data(mg1 + 538);
    const auto *mg1_539 = buffer.data(mg1 + 539);

    const auto *mh_681 = buffer.data(mh + 681);
    const auto *mh_684 = buffer.data(mh + 684);
    const auto *mh_687 = buffer.data(mh + 687);
    const auto *mh_688 = buffer.data(mh + 688);
    const auto *mh_689 = buffer.data(mh + 689);
    const auto *mh_690 = buffer.data(mh + 690);
    const auto *mh_691 = buffer.data(mh + 691);
    const auto *mh_692 = buffer.data(mh + 692);
    const auto *mh_693 = buffer.data(mh + 693);
    const auto *mh_695 = buffer.data(mh + 695);
    const auto *mh_696 = buffer.data(mh + 696);
    const auto *mh_698 = buffer.data(mh + 698);
    const auto *mh_699 = buffer.data(mh + 699);
    const auto *mh_702 = buffer.data(mh + 702);
    const auto *mh_705 = buffer.data(mh + 705);
    const auto *mh_708 = buffer.data(mh + 708);
    const auto *mh_709 = buffer.data(mh + 709);
    const auto *mh_710 = buffer.data(mh + 710);
    const auto *mh_711 = buffer.data(mh + 711);
    const auto *mh_712 = buffer.data(mh + 712);
    const auto *mh_713 = buffer.data(mh + 713);
    const auto *mh_714 = buffer.data(mh + 714);
    const auto *mh_716 = buffer.data(mh + 716);
    const auto *mh_717 = buffer.data(mh + 717);
    const auto *mh_719 = buffer.data(mh + 719);
    const auto *mh_720 = buffer.data(mh + 720);
    const auto *mh_723 = buffer.data(mh + 723);
    const auto *mh_729 = buffer.data(mh + 729);
    const auto *mh_730 = buffer.data(mh + 730);
    const auto *mh_731 = buffer.data(mh + 731);
    const auto *mh_732 = buffer.data(mh + 732);
    const auto *mh_733 = buffer.data(mh + 733);
    const auto *mh_734 = buffer.data(mh + 734);
    const auto *mh_735 = buffer.data(mh + 735);
    const auto *mh_736 = buffer.data(mh + 736);
    const auto *mh_737 = buffer.data(mh + 737);
    const auto *mh_738 = buffer.data(mh + 738);
    const auto *mh_740 = buffer.data(mh + 740);
    const auto *mh_741 = buffer.data(mh + 741);
    const auto *mh_743 = buffer.data(mh + 743);
    const auto *mh_744 = buffer.data(mh + 744);
    const auto *mh_749 = buffer.data(mh + 749);
    const auto *mh_750 = buffer.data(mh + 750);
    const auto *mh_751 = buffer.data(mh + 751);
    const auto *mh_752 = buffer.data(mh + 752);
    const auto *mh_753 = buffer.data(mh + 753);
    const auto *mh_754 = buffer.data(mh + 754);
    const auto *mh_755 = buffer.data(mh + 755);
    const auto *mh_756 = buffer.data(mh + 756);
    const auto *mh_757 = buffer.data(mh + 757);
    const auto *mh_759 = buffer.data(mh + 759);
    const auto *mh_761 = buffer.data(mh + 761);
    const auto *mh_762 = buffer.data(mh + 762);
    const auto *mh_765 = buffer.data(mh + 765);
    const auto *mh_766 = buffer.data(mh + 766);
    const auto *mh_771 = buffer.data(mh + 771);
    const auto *mh_773 = buffer.data(mh + 773);
    const auto *mh_774 = buffer.data(mh + 774);
    const auto *mh_775 = buffer.data(mh + 775);
    const auto *mh_776 = buffer.data(mh + 776);

#pragma omp simd aligned(t_908, t_909, t_910, pa_y, pb_x, pb_y, ki0_546, ki1_546, lh_534, \
                         lh_684, li_714, mg0_492, mg1_492, mh_681, \
                         mh_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_10 * lh_684[k]
                   + f_3 * mg0_492[k]
                   - f_4 * mg1_492[k]
                   + pb_x[k] * mh_684[k];

        t_909[k] = f_11 * lh_534[k]
                   + pb_y[k] * mh_681[k];

        t_910[k] = f_20 * ki0_546[k]
                   - f_21 * ki1_546[k]
                   + pa_y[k] * li_714[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, t_914, t_915, pb_x, lh_687, lh_688, lh_689, \
                         lh_690, lh_691, mh_687, mh_688, mh_689, mh_690, \
                         mh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_10 * lh_687[k]
                   + pb_x[k] * mh_687[k];

        t_912[k] = f_10 * lh_688[k]
                   + pb_x[k] * mh_688[k];

        t_913[k] = f_10 * lh_689[k]
                   + pb_x[k] * mh_689[k];

        t_914[k] = f_10 * lh_690[k]
                   + pb_x[k] * mh_690[k];

        t_915[k] = f_10 * lh_691[k]
                   + pb_x[k] * mh_691[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, pa_x, pb_x, pb_z, ki0_917, ki1_917, lh_519, \
                         lh_692, li_917, mh_687, mh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_10 * lh_692[k]
                   + pb_x[k] * mh_692[k];

        t_917[k] = f_15 * ki0_917[k]
                   - f_16 * ki1_917[k]
                   + pa_x[k] * li_917[k];

        t_918[k] = f_12 * lh_519[k]
                   + pb_z[k] * mh_687[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, pa_x, ki0_919, ki0_920, ki0_921, ki1_919, \
                         ki1_920, ki1_921, li_919, li_920, li_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = f_15 * ki0_919[k]
                   - f_16 * ki1_919[k]
                   + pa_x[k] * li_919[k];

        t_920[k] = f_15 * ki0_920[k]
                   - f_16 * ki1_920[k]
                   + pa_x[k] * li_920[k];

        t_921[k] = f_15 * ki0_921[k]
                   - f_16 * ki1_921[k]
                   + pa_x[k] * li_921[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, pa_x, pa_y, pb_y, ki0_560, ki0_923, ki1_560, \
                         ki1_923, lh_545, li_728, li_923, mh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_11 * lh_545[k]
                   + pb_y[k] * mh_692[k];

        t_923[k] = f_15 * ki0_923[k]
                   - f_16 * ki1_923[k]
                   + pa_x[k] * li_923[k];

        t_924[k] = f_15 * ki0_560[k]
                   - f_16 * ki1_560[k]
                   + pa_y[k] * li_728[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, pa_z, pb_y, pb_z, ki0_507, ki1_507, \
                         lh_525, lh_546, lh_548, li_703, mh_693, \
                         mh_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = f_10 * lh_546[k]
                   + pb_y[k] * mh_693[k];

        t_926[k] = f_26 * lh_525[k]
                   + pb_z[k] * mh_693[k];

        t_927[k] = f_27 * ki0_507[k]
                   - f_28 * ki1_507[k]
                   + pa_z[k] * li_703[k];

        t_928[k] = f_10 * lh_548[k]
                   + pb_y[k] * mh_695[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, pa_y, pa_z, pb_z, ki0_510, ki0_565, ki1_510, \
                         ki1_565, lh_528, li_706, li_733, mh_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = f_15 * ki0_565[k]
                   - f_16 * ki1_565[k]
                   + pa_y[k] * li_733[k];

        t_930[k] = f_27 * ki0_510[k]
                   - f_28 * ki1_510[k]
                   + pa_z[k] * li_706[k];

        t_931[k] = f_26 * lh_528[k]
                   + pb_z[k] * mh_696[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, pa_y, pa_z, pb_y, ki0_514, ki0_569, ki1_514, \
                         ki1_569, lh_551, li_710, li_737, mh_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_10 * lh_551[k]
                   + pb_y[k] * mh_698[k];

        t_933[k] = f_15 * ki0_569[k]
                   - f_16 * ki1_569[k]
                   + pa_y[k] * li_737[k];

        t_934[k] = f_27 * ki0_514[k]
                   - f_28 * ki1_514[k]
                   + pa_z[k] * li_710[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, pb_x, pb_y, pb_z, lh_531, lh_555, lh_705, \
                         mg0_507, mg1_507, mh_699, mh_702, mh_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_26 * lh_531[k]
                   + pb_z[k] * mh_699[k];

        t_936[k] = f_10 * lh_705[k]
                   + f_3 * mg0_507[k]
                   - f_4 * mg1_507[k]
                   + pb_x[k] * mh_705[k];

        t_937[k] = f_10 * lh_555[k]
                   + pb_y[k] * mh_702[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pa_y, pb_x, ki0_574, ki1_574, lh_708, \
                         lh_709, lh_710, li_742, mh_708, mh_709, \
                         mh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_15 * ki0_574[k]
                   - f_16 * ki1_574[k]
                   + pa_y[k] * li_742[k];

        t_939[k] = f_10 * lh_708[k]
                   + pb_x[k] * mh_708[k];

        t_940[k] = f_10 * lh_709[k]
                   + pb_x[k] * mh_709[k];

        t_941[k] = f_10 * lh_710[k]
                   + pb_x[k] * mh_710[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pa_x, pb_x, ki0_945, ki1_945, lh_711, \
                         lh_712, lh_713, li_945, mh_711, mh_712, \
                         mh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_10 * lh_711[k]
                   + pb_x[k] * mh_711[k];

        t_943[k] = f_10 * lh_712[k]
                   + pb_x[k] * mh_712[k];

        t_944[k] = f_10 * lh_713[k]
                   + pb_x[k] * mh_713[k];

        t_945[k] = f_15 * ki0_945[k]
                   - f_16 * ki1_945[k]
                   + pa_x[k] * li_945[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pa_x, pb_z, ki0_947, ki0_948, ki1_947, ki1_948, \
                         lh_540, li_947, li_948, mh_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_26 * lh_540[k]
                   + pb_z[k] * mh_708[k];

        t_947[k] = f_15 * ki0_947[k]
                   - f_16 * ki1_947[k]
                   + pa_x[k] * li_947[k];

        t_948[k] = f_15 * ki0_948[k]
                   - f_16 * ki1_948[k]
                   + pa_x[k] * li_948[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pa_x, pa_y, pb_y, ki0_949, ki0_951, \
                         ki1_949, ki1_951, lh_566, li_756, li_949, li_951, \
                         mh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_15 * ki0_949[k]
                   - f_16 * ki1_949[k]
                   + pa_x[k] * li_949[k];

        t_950[k] = f_10 * lh_566[k]
                   + pb_y[k] * mh_713[k];

        t_951[k] = f_15 * ki0_951[k]
                   - f_16 * ki1_951[k]
                   + pa_x[k] * li_951[k];

        t_952[k] = pa_y[k] * li_756[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, t_957, pa_y, pb_y, lh_567, lh_568, \
                         lh_569, li_758, li_759, li_761, mh_714, \
                         mh_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_9 * lh_567[k]
                   + pb_y[k] * mh_714[k];

        t_954[k] = pa_y[k] * li_758[k];

        t_955[k] = f_10 * lh_568[k]
                   + pa_y[k] * li_759[k];

        t_956[k] = f_9 * lh_569[k]
                   + pb_y[k] * mh_716[k];

        t_957[k] = pa_y[k] * li_761[k];
    }

#pragma omp simd aligned(t_958, t_959, t_960, t_961, pa_y, pb_y, pb_z, lh_549, lh_570, lh_572, \
                         li_762, li_765, mh_717, mh_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_958[k] = f_11 * lh_570[k]
                   + pa_y[k] * li_762[k];

        t_959[k] = f_14 * lh_549[k]
                   + pb_z[k] * mh_717[k];

        t_960[k] = f_9 * lh_572[k]
                   + pb_y[k] * mh_719[k];

        t_961[k] = pa_y[k] * li_765[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, pa_y, pb_y, pb_z, lh_552, lh_573, lh_575, \
                         lh_576, li_766, li_768, mh_720, mh_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_12 * lh_573[k]
                   + pa_y[k] * li_766[k];

        t_963[k] = f_14 * lh_552[k]
                   + pb_z[k] * mh_720[k];

        t_964[k] = f_10 * lh_575[k]
                   + pa_y[k] * li_768[k];

        t_965[k] = f_9 * lh_576[k]
                   + pb_y[k] * mh_723[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, t_969, t_970, pa_y, pb_x, lh_729, lh_730, \
                         lh_731, lh_732, li_770, mh_729, mh_730, mh_731, \
                         mh_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = pa_y[k] * li_770[k];

        t_967[k] = f_10 * lh_729[k]
                   + pb_x[k] * mh_729[k];

        t_968[k] = f_10 * lh_730[k]
                   + pb_x[k] * mh_730[k];

        t_969[k] = f_10 * lh_731[k]
                   + pb_x[k] * mh_731[k];

        t_970[k] = f_10 * lh_732[k]
                   + pb_x[k] * mh_732[k];
    }

#pragma omp simd aligned(t_971, t_972, t_973, t_974, pa_y, pb_x, pb_z, lh_561, lh_582, lh_733, \
                         li_776, li_777, mh_729, mh_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_971[k] = f_10 * lh_733[k]
                   + pb_x[k] * mh_733[k];

        t_972[k] = pa_y[k] * li_776[k];

        t_973[k] = f_14 * lh_582[k]
                   + pa_y[k] * li_777[k];

        t_974[k] = f_14 * lh_561[k]
                   + pb_z[k] * mh_729[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, t_979, pa_y, pb_y, lh_584, lh_585, \
                         lh_586, lh_587, li_779, li_780, li_781, li_783, \
                         mh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = f_12 * lh_584[k]
                   + pa_y[k] * li_779[k];

        t_976[k] = f_11 * lh_585[k]
                   + pa_y[k] * li_780[k];

        t_977[k] = f_10 * lh_586[k]
                   + pa_y[k] * li_781[k];

        t_978[k] = f_9 * lh_587[k]
                   + pb_y[k] * mh_734[k];

        t_979[k] = pa_y[k] * li_783[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, pa_z, pb_y, pb_z, ki0_560, ki1_560, \
                         lh_567, li_756, mg0_525, mg1_525, mh_735, \
                         mh_736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = f_18 * ki0_560[k]
                   - f_19 * ki1_560[k]
                   + pa_z[k] * li_756[k];

        t_981[k] = pb_y[k] * mh_735[k];

        t_982[k] = f_17 * lh_567[k]
                   + pb_z[k] * mh_735[k];

        t_983[k] = f_3 * mg0_525[k]
                   - f_4 * mg1_525[k]
                   + pb_y[k] * mh_736[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, t_987, pb_x, pb_y, pb_z, lh_570, lh_740, \
                         mg0_526, mg0_530, mg1_526, mg1_530, mh_737, mh_738, \
                         mh_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = pb_y[k] * mh_737[k];

        t_985[k] = f_10 * lh_740[k]
                   + f_7 * mg0_530[k]
                   - f_8 * mg1_530[k]
                   + pb_x[k] * mh_740[k];

        t_986[k] = f_5 * mg0_526[k]
                   - f_6 * mg1_526[k]
                   + pb_y[k] * mh_738[k];

        t_987[k] = f_17 * lh_570[k]
                   + pb_z[k] * mh_738[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, t_991, pb_x, pb_y, pb_z, lh_573, lh_744, \
                         mg0_528, mg0_534, mg1_528, mg1_534, mh_740, mh_741, \
                         mh_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = pb_y[k] * mh_740[k];

        t_989[k] = f_10 * lh_744[k]
                   + f_5 * mg0_534[k]
                   - f_6 * mg1_534[k]
                   + pb_x[k] * mh_744[k];

        t_990[k] = f_7 * mg0_528[k]
                   - f_8 * mg1_528[k]
                   + pb_y[k] * mh_741[k];

        t_991[k] = f_17 * lh_573[k]
                   + pb_z[k] * mh_741[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, pb_x, pb_y, lh_749, lh_750, mg0_530, \
                         mg0_539, mg1_530, mg1_539, mh_743, mh_744, mh_749, \
                         mh_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_3 * mg0_530[k]
                   - f_4 * mg1_530[k]
                   + pb_y[k] * mh_743[k];

        t_993[k] = pb_y[k] * mh_744[k];

        t_994[k] = f_10 * lh_749[k]
                   + f_3 * mg0_539[k]
                   - f_4 * mg1_539[k]
                   + pb_x[k] * mh_749[k];

        t_995[k] = f_10 * lh_750[k]
                   + pb_x[k] * mh_750[k];
    }

#pragma omp simd aligned(t_996, t_997, t_998, t_999, t_1000, pb_x, pb_y, lh_751, lh_752, \
                         lh_753, lh_755, mh_749, mh_751, mh_752, mh_753, \
                         mh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_996[k] = f_10 * lh_751[k]
                   + pb_x[k] * mh_751[k];

        t_997[k] = f_10 * lh_752[k]
                   + pb_x[k] * mh_752[k];

        t_998[k] = f_10 * lh_753[k]
                   + pb_x[k] * mh_753[k];

        t_999[k] = pb_y[k] * mh_749[k];

        t_1000[k] = f_10 * lh_755[k]
                    + pb_x[k] * mh_755[k];
    }

#pragma omp simd aligned(t_1001, t_1002, t_1003, t_1004, pb_y, pb_z, lh_582, mg0_535, mg0_537, \
                         mg0_538, mg1_535, mg1_537, mg1_538, mh_750, mh_752, \
                         mh_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1001[k] = f_1 * mg0_535[k]
                    - f_2 * mg1_535[k]
                    + pb_y[k] * mh_750[k];

        t_1002[k] = f_17 * lh_582[k]
                    + pb_z[k] * mh_750[k];

        t_1003[k] = f_7 * mg0_537[k]
                    - f_8 * mg1_537[k]
                    + pb_y[k] * mh_752[k];

        t_1004[k] = f_5 * mg0_538[k]
                    - f_6 * mg1_538[k]
                    + pb_y[k] * mh_753[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, t_1008, pa_x, pb_y, ki0_1007, ki1_1007, \
                         lh_756, li_1007, li_1008, mg0_539, mg1_539, mh_754, \
                         mh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = f_3 * mg0_539[k]
                    - f_4 * mg1_539[k]
                    + pb_y[k] * mh_754[k];

        t_1006[k] = pb_y[k] * mh_755[k];

        t_1007[k] = f_15 * ki0_1007[k]
                    - f_16 * ki1_1007[k]
                    + pa_x[k] * li_1007[k];

        t_1008[k] = f_14 * lh_756[k]
                    + pa_x[k] * li_1008[k];
    }

#pragma omp simd aligned(t_1009, t_1010, t_1011, t_1012, t_1013, pa_x, pb_y, pb_z, lh_588, \
                         lh_759, lh_761, li_1011, li_1013, mh_756, \
                         mh_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1009[k] = f_13 * lh_588[k]
                    + pb_y[k] * mh_756[k];

        t_1010[k] = pb_z[k] * mh_756[k];

        t_1011[k] = f_12 * lh_759[k]
                    + pa_x[k] * li_1011[k];

        t_1012[k] = pb_z[k] * mh_757[k];

        t_1013[k] = f_12 * lh_761[k]
                    + pa_x[k] * li_1013[k];
    }

#pragma omp simd aligned(t_1014, t_1015, t_1016, t_1017, pa_x, pb_y, pb_z, lh_593, lh_762, \
                         lh_765, li_1014, li_1017, mh_759, mh_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1014[k] = f_11 * lh_762[k]
                    + pa_x[k] * li_1014[k];

        t_1015[k] = pb_z[k] * mh_759[k];

        t_1016[k] = f_13 * lh_593[k]
                    + pb_y[k] * mh_761[k];

        t_1017[k] = f_11 * lh_765[k]
                    + pa_x[k] * li_1017[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, t_1021, pa_x, pb_y, pb_z, lh_597, lh_766, \
                         lh_768, li_1018, li_1020, mh_762, mh_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = f_10 * lh_766[k]
                    + pa_x[k] * li_1018[k];

        t_1019[k] = pb_z[k] * mh_762[k];

        t_1020[k] = f_10 * lh_768[k]
                    + pa_x[k] * li_1020[k];

        t_1021[k] = f_13 * lh_597[k]
                    + pb_y[k] * mh_765[k];
    }

#pragma omp simd aligned(t_1022, t_1023, t_1024, t_1025, pa_x, pb_x, pb_z, lh_770, lh_771, \
                         lh_773, li_1022, mh_766, mh_771, mh_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1022[k] = f_10 * lh_770[k]
                    + pa_x[k] * li_1022[k];

        t_1023[k] = f_9 * lh_771[k]
                    + pb_x[k] * mh_771[k];

        t_1024[k] = pb_z[k] * mh_766[k];

        t_1025[k] = f_9 * lh_773[k]
                    + pb_x[k] * mh_773[k];
    }

#pragma omp simd aligned(t_1026, t_1027, t_1028, t_1029, t_1030, pa_x, pb_x, pb_z, lh_774, \
                         lh_775, lh_776, li_1029, mh_771, mh_774, mh_775, \
                         mh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1026[k] = f_9 * lh_774[k]
                    + pb_x[k] * mh_774[k];

        t_1027[k] = f_9 * lh_775[k]
                    + pb_x[k] * mh_775[k];

        t_1028[k] = f_9 * lh_776[k]
                    + pb_x[k] * mh_776[k];

        t_1029[k] = pa_x[k] * li_1029[k];

        t_1030[k] = pb_z[k] * mh_771[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, t_1034, t_1035, t_1036, t_1037, pa_x, pa_z, \
                         li_784, li_785, li_1031, li_1032, li_1033, li_1034, \
                         li_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = pa_x[k] * li_1031[k];

        t_1032[k] = pa_x[k] * li_1032[k];

        t_1033[k] = pa_x[k] * li_1033[k];

        t_1034[k] = pa_x[k] * li_1034[k];

        t_1035[k] = pa_x[k] * li_1035[k];

        t_1036[k] = pa_z[k] * li_784[k];

        t_1037[k] = pa_z[k] * li_785[k];
    }
}

static auto
compute_prim_mi_electron_repulsion_0_piece8(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t lh,
                                            const size_t li, const size_t mh, const size_t ncols,
                                            const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_9 = 0.5 / p;
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_14 = 3.0 / p;
    const auto f_17 = 3.5 / p;
    const auto f_26 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *lh_588 = buffer.data(lh + 588);
    const auto *lh_591 = buffer.data(lh + 591);
    const auto *lh_594 = buffer.data(lh + 594);
    const auto *lh_609 = buffer.data(lh + 609);
    const auto *lh_611 = buffer.data(lh + 611);
    const auto *lh_612 = buffer.data(lh + 612);
    const auto *lh_614 = buffer.data(lh + 614);
    const auto *lh_615 = buffer.data(lh + 615);
    const auto *lh_618 = buffer.data(lh + 618);
    const auto *lh_630 = buffer.data(lh + 630);
    const auto *lh_632 = buffer.data(lh + 632);
    const auto *lh_633 = buffer.data(lh + 633);
    const auto *lh_635 = buffer.data(lh + 635);
    const auto *lh_636 = buffer.data(lh + 636);
    const auto *lh_639 = buffer.data(lh + 639);
    const auto *lh_651 = buffer.data(lh + 651);
    const auto *lh_653 = buffer.data(lh + 653);
    const auto *lh_654 = buffer.data(lh + 654);
    const auto *lh_656 = buffer.data(lh + 656);
    const auto *lh_657 = buffer.data(lh + 657);
    const auto *lh_660 = buffer.data(lh + 660);
    const auto *lh_672 = buffer.data(lh + 672);
    const auto *lh_674 = buffer.data(lh + 674);
    const auto *lh_675 = buffer.data(lh + 675);
    const auto *lh_677 = buffer.data(lh + 677);
    const auto *lh_678 = buffer.data(lh + 678);
    const auto *lh_681 = buffer.data(lh + 681);
    const auto *lh_693 = buffer.data(lh + 693);
    const auto *lh_695 = buffer.data(lh + 695);
    const auto *lh_696 = buffer.data(lh + 696);
    const auto *lh_698 = buffer.data(lh + 698);
    const auto *lh_702 = buffer.data(lh + 702);
    const auto *lh_714 = buffer.data(lh + 714);
    const auto *lh_716 = buffer.data(lh + 716);
    const auto *lh_719 = buffer.data(lh + 719);
    const auto *lh_782 = buffer.data(lh + 782);
    const auto *lh_786 = buffer.data(lh + 786);
    const auto *lh_789 = buffer.data(lh + 789);
    const auto *lh_791 = buffer.data(lh + 791);
    const auto *lh_793 = buffer.data(lh + 793);
    const auto *lh_794 = buffer.data(lh + 794);
    const auto *lh_795 = buffer.data(lh + 795);
    const auto *lh_796 = buffer.data(lh + 796);
    const auto *lh_797 = buffer.data(lh + 797);
    const auto *lh_798 = buffer.data(lh + 798);
    const auto *lh_801 = buffer.data(lh + 801);
    const auto *lh_803 = buffer.data(lh + 803);
    const auto *lh_804 = buffer.data(lh + 804);
    const auto *lh_807 = buffer.data(lh + 807);
    const auto *lh_808 = buffer.data(lh + 808);
    const auto *lh_810 = buffer.data(lh + 810);
    const auto *lh_812 = buffer.data(lh + 812);
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
    const auto *lh_829 = buffer.data(lh + 829);
    const auto *lh_831 = buffer.data(lh + 831);
    const auto *lh_833 = buffer.data(lh + 833);
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
    const auto *lh_850 = buffer.data(lh + 850);
    const auto *lh_852 = buffer.data(lh + 852);
    const auto *lh_854 = buffer.data(lh + 854);
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
    const auto *lh_871 = buffer.data(lh + 871);
    const auto *lh_873 = buffer.data(lh + 873);
    const auto *lh_875 = buffer.data(lh + 875);
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

    const auto *li_787 = buffer.data(li + 787);
    const auto *li_790 = buffer.data(li + 790);
    const auto *li_794 = buffer.data(li + 794);
    const auto *li_799 = buffer.data(li + 799);
    const auto *li_1041 = buffer.data(li + 1041);
    const auto *li_1045 = buffer.data(li + 1045);
    const auto *li_1048 = buffer.data(li + 1048);
    const auto *li_1050 = buffer.data(li + 1050);
    const auto *li_1057 = buffer.data(li + 1057);
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
    const auto *li_1113 = buffer.data(li + 1113);
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

    const auto *mh_777 = buffer.data(mh + 777);
    const auto *mh_779 = buffer.data(mh + 779);
    const auto *mh_780 = buffer.data(mh + 780);
    const auto *mh_782 = buffer.data(mh + 782);
    const auto *mh_783 = buffer.data(mh + 783);
    const auto *mh_786 = buffer.data(mh + 786);
    const auto *mh_793 = buffer.data(mh + 793);
    const auto *mh_794 = buffer.data(mh + 794);
    const auto *mh_795 = buffer.data(mh + 795);
    const auto *mh_796 = buffer.data(mh + 796);
    const auto *mh_797 = buffer.data(mh + 797);
    const auto *mh_798 = buffer.data(mh + 798);
    const auto *mh_800 = buffer.data(mh + 800);
    const auto *mh_801 = buffer.data(mh + 801);
    const auto *mh_803 = buffer.data(mh + 803);
    const auto *mh_804 = buffer.data(mh + 804);
    const auto *mh_807 = buffer.data(mh + 807);
    const auto *mh_813 = buffer.data(mh + 813);
    const auto *mh_814 = buffer.data(mh + 814);
    const auto *mh_815 = buffer.data(mh + 815);
    const auto *mh_816 = buffer.data(mh + 816);
    const auto *mh_817 = buffer.data(mh + 817);
    const auto *mh_818 = buffer.data(mh + 818);
    const auto *mh_819 = buffer.data(mh + 819);
    const auto *mh_821 = buffer.data(mh + 821);
    const auto *mh_822 = buffer.data(mh + 822);
    const auto *mh_824 = buffer.data(mh + 824);
    const auto *mh_825 = buffer.data(mh + 825);
    const auto *mh_828 = buffer.data(mh + 828);
    const auto *mh_834 = buffer.data(mh + 834);
    const auto *mh_835 = buffer.data(mh + 835);
    const auto *mh_836 = buffer.data(mh + 836);
    const auto *mh_837 = buffer.data(mh + 837);
    const auto *mh_838 = buffer.data(mh + 838);
    const auto *mh_839 = buffer.data(mh + 839);
    const auto *mh_840 = buffer.data(mh + 840);
    const auto *mh_842 = buffer.data(mh + 842);
    const auto *mh_843 = buffer.data(mh + 843);
    const auto *mh_845 = buffer.data(mh + 845);
    const auto *mh_846 = buffer.data(mh + 846);
    const auto *mh_849 = buffer.data(mh + 849);
    const auto *mh_855 = buffer.data(mh + 855);
    const auto *mh_856 = buffer.data(mh + 856);
    const auto *mh_857 = buffer.data(mh + 857);
    const auto *mh_858 = buffer.data(mh + 858);
    const auto *mh_859 = buffer.data(mh + 859);
    const auto *mh_860 = buffer.data(mh + 860);
    const auto *mh_861 = buffer.data(mh + 861);
    const auto *mh_863 = buffer.data(mh + 863);
    const auto *mh_864 = buffer.data(mh + 864);
    const auto *mh_866 = buffer.data(mh + 866);
    const auto *mh_867 = buffer.data(mh + 867);
    const auto *mh_870 = buffer.data(mh + 870);
    const auto *mh_876 = buffer.data(mh + 876);
    const auto *mh_877 = buffer.data(mh + 877);
    const auto *mh_878 = buffer.data(mh + 878);
    const auto *mh_879 = buffer.data(mh + 879);
    const auto *mh_880 = buffer.data(mh + 880);
    const auto *mh_881 = buffer.data(mh + 881);
    const auto *mh_882 = buffer.data(mh + 882);
    const auto *mh_884 = buffer.data(mh + 884);
    const auto *mh_885 = buffer.data(mh + 885);
    const auto *mh_887 = buffer.data(mh + 887);

#pragma omp simd aligned(t_1038, t_1039, t_1040, t_1041, pa_x, pa_z, pb_y, pb_z, lh_588, \
                         lh_611, lh_782, li_787, li_1041, mh_777, \
                         mh_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1038[k] = f_9 * lh_588[k]
                    + pb_z[k] * mh_777[k];

        t_1039[k] = pa_z[k] * li_787[k];

        t_1040[k] = f_17 * lh_611[k]
                    + pb_y[k] * mh_779[k];

        t_1041[k] = f_12 * lh_782[k]
                    + pa_x[k] * li_1041[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, t_1045, pa_x, pa_z, pb_y, pb_z, lh_591, \
                         lh_614, lh_786, li_790, li_1045, mh_780, \
                         mh_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = pa_z[k] * li_790[k];

        t_1043[k] = f_9 * lh_591[k]
                    + pb_z[k] * mh_780[k];

        t_1044[k] = f_17 * lh_614[k]
                    + pb_y[k] * mh_782[k];

        t_1045[k] = f_11 * lh_786[k]
                    + pa_x[k] * li_1045[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, t_1049, pa_x, pa_z, pb_y, pb_z, lh_594, \
                         lh_618, lh_789, li_794, li_1048, mh_783, \
                         mh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = pa_z[k] * li_794[k];

        t_1047[k] = f_9 * lh_594[k]
                    + pb_z[k] * mh_783[k];

        t_1048[k] = f_10 * lh_789[k]
                    + pa_x[k] * li_1048[k];

        t_1049[k] = f_17 * lh_618[k]
                    + pb_y[k] * mh_786[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, t_1053, pa_x, pa_z, pb_x, lh_791, lh_793, \
                         lh_794, li_799, li_1050, mh_793, mh_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = f_10 * lh_791[k]
                    + pa_x[k] * li_1050[k];

        t_1051[k] = pa_z[k] * li_799[k];

        t_1052[k] = f_9 * lh_793[k]
                    + pb_x[k] * mh_793[k];

        t_1053[k] = f_9 * lh_794[k]
                    + pb_x[k] * mh_794[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, t_1057, t_1058, pa_x, pb_x, lh_795, lh_796, \
                         lh_797, li_1057, li_1058, mh_795, mh_796, \
                         mh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = f_9 * lh_795[k]
                    + pb_x[k] * mh_795[k];

        t_1055[k] = f_9 * lh_796[k]
                    + pb_x[k] * mh_796[k];

        t_1056[k] = f_9 * lh_797[k]
                    + pb_x[k] * mh_797[k];

        t_1057[k] = pa_x[k] * li_1057[k];

        t_1058[k] = pa_x[k] * li_1058[k];
    }

#pragma omp simd aligned(t_1059, t_1060, t_1061, t_1062, t_1063, t_1064, pa_x, lh_798, \
                         li_1059, li_1060, li_1061, li_1062, li_1063, \
                         li_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = pa_x[k] * li_1059[k];

        t_1060[k] = pa_x[k] * li_1060[k];

        t_1061[k] = pa_x[k] * li_1061[k];

        t_1062[k] = pa_x[k] * li_1062[k];

        t_1063[k] = pa_x[k] * li_1063[k];

        t_1064[k] = f_14 * lh_798[k]
                    + pa_x[k] * li_1064[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, pa_x, pb_y, pb_z, lh_609, lh_630, \
                         lh_632, lh_801, li_1067, mh_798, mh_800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_14 * lh_630[k]
                    + pb_y[k] * mh_798[k];

        t_1066[k] = f_10 * lh_609[k]
                    + pb_z[k] * mh_798[k];

        t_1067[k] = f_12 * lh_801[k]
                    + pa_x[k] * li_1067[k];

        t_1068[k] = f_14 * lh_632[k]
                    + pb_y[k] * mh_800[k];
    }

#pragma omp simd aligned(t_1069, t_1070, t_1071, t_1072, pa_x, pb_y, pb_z, lh_612, lh_635, \
                         lh_803, lh_804, li_1069, li_1070, mh_801, \
                         mh_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1069[k] = f_12 * lh_803[k]
                    + pa_x[k] * li_1069[k];

        t_1070[k] = f_11 * lh_804[k]
                    + pa_x[k] * li_1070[k];

        t_1071[k] = f_10 * lh_612[k]
                    + pb_z[k] * mh_801[k];

        t_1072[k] = f_14 * lh_635[k]
                    + pb_y[k] * mh_803[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, t_1076, pa_x, pb_z, lh_615, lh_807, lh_808, \
                         lh_810, li_1073, li_1074, li_1076, mh_804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_11 * lh_807[k]
                    + pa_x[k] * li_1073[k];

        t_1074[k] = f_10 * lh_808[k]
                    + pa_x[k] * li_1074[k];

        t_1075[k] = f_10 * lh_615[k]
                    + pb_z[k] * mh_804[k];

        t_1076[k] = f_10 * lh_810[k]
                    + pa_x[k] * li_1076[k];
    }

#pragma omp simd aligned(t_1077, t_1078, t_1079, t_1080, pa_x, pb_x, pb_y, lh_639, lh_812, \
                         lh_813, lh_814, li_1078, mh_807, mh_813, \
                         mh_814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1077[k] = f_14 * lh_639[k]
                    + pb_y[k] * mh_807[k];

        t_1078[k] = f_10 * lh_812[k]
                    + pa_x[k] * li_1078[k];

        t_1079[k] = f_9 * lh_813[k]
                    + pb_x[k] * mh_813[k];

        t_1080[k] = f_9 * lh_814[k]
                    + pb_x[k] * mh_814[k];
    }

#pragma omp simd aligned(t_1081, t_1082, t_1083, t_1084, t_1085, pa_x, pb_x, lh_815, lh_816, \
                         lh_817, lh_818, li_1085, mh_815, mh_816, mh_817, \
                         mh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1081[k] = f_9 * lh_815[k]
                    + pb_x[k] * mh_815[k];

        t_1082[k] = f_9 * lh_816[k]
                    + pb_x[k] * mh_816[k];

        t_1083[k] = f_9 * lh_817[k]
                    + pb_x[k] * mh_817[k];

        t_1084[k] = f_9 * lh_818[k]
                    + pb_x[k] * mh_818[k];

        t_1085[k] = pa_x[k] * li_1085[k];
    }

#pragma omp simd aligned(t_1086, t_1087, t_1088, t_1089, t_1090, t_1091, t_1092, pa_x, lh_819, \
                         li_1086, li_1087, li_1088, li_1089, li_1090, li_1091, \
                         li_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1086[k] = pa_x[k] * li_1086[k];

        t_1087[k] = pa_x[k] * li_1087[k];

        t_1088[k] = pa_x[k] * li_1088[k];

        t_1089[k] = pa_x[k] * li_1089[k];

        t_1090[k] = pa_x[k] * li_1090[k];

        t_1091[k] = pa_x[k] * li_1091[k];

        t_1092[k] = f_14 * lh_819[k]
                    + pa_x[k] * li_1092[k];
    }

#pragma omp simd aligned(t_1093, t_1094, t_1095, t_1096, pa_x, pb_y, pb_z, lh_630, lh_651, \
                         lh_653, lh_822, li_1095, mh_819, mh_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1093[k] = f_26 * lh_651[k]
                    + pb_y[k] * mh_819[k];

        t_1094[k] = f_11 * lh_630[k]
                    + pb_z[k] * mh_819[k];

        t_1095[k] = f_12 * lh_822[k]
                    + pa_x[k] * li_1095[k];

        t_1096[k] = f_26 * lh_653[k]
                    + pb_y[k] * mh_821[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, t_1100, pa_x, pb_y, pb_z, lh_633, lh_656, \
                         lh_824, lh_825, li_1097, li_1098, mh_822, \
                         mh_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_12 * lh_824[k]
                    + pa_x[k] * li_1097[k];

        t_1098[k] = f_11 * lh_825[k]
                    + pa_x[k] * li_1098[k];

        t_1099[k] = f_11 * lh_633[k]
                    + pb_z[k] * mh_822[k];

        t_1100[k] = f_26 * lh_656[k]
                    + pb_y[k] * mh_824[k];
    }

#pragma omp simd aligned(t_1101, t_1102, t_1103, t_1104, pa_x, pb_z, lh_636, lh_828, lh_829, \
                         lh_831, li_1101, li_1102, li_1104, mh_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1101[k] = f_11 * lh_828[k]
                    + pa_x[k] * li_1101[k];

        t_1102[k] = f_10 * lh_829[k]
                    + pa_x[k] * li_1102[k];

        t_1103[k] = f_11 * lh_636[k]
                    + pb_z[k] * mh_825[k];

        t_1104[k] = f_10 * lh_831[k]
                    + pa_x[k] * li_1104[k];
    }

#pragma omp simd aligned(t_1105, t_1106, t_1107, t_1108, pa_x, pb_x, pb_y, lh_660, lh_833, \
                         lh_834, lh_835, li_1106, mh_828, mh_834, \
                         mh_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1105[k] = f_26 * lh_660[k]
                    + pb_y[k] * mh_828[k];

        t_1106[k] = f_10 * lh_833[k]
                    + pa_x[k] * li_1106[k];

        t_1107[k] = f_9 * lh_834[k]
                    + pb_x[k] * mh_834[k];

        t_1108[k] = f_9 * lh_835[k]
                    + pb_x[k] * mh_835[k];
    }

#pragma omp simd aligned(t_1109, t_1110, t_1111, t_1112, t_1113, pa_x, pb_x, lh_836, lh_837, \
                         lh_838, lh_839, li_1113, mh_836, mh_837, mh_838, \
                         mh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1109[k] = f_9 * lh_836[k]
                    + pb_x[k] * mh_836[k];

        t_1110[k] = f_9 * lh_837[k]
                    + pb_x[k] * mh_837[k];

        t_1111[k] = f_9 * lh_838[k]
                    + pb_x[k] * mh_838[k];

        t_1112[k] = f_9 * lh_839[k]
                    + pb_x[k] * mh_839[k];

        t_1113[k] = pa_x[k] * li_1113[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, t_1117, t_1118, t_1119, t_1120, pa_x, lh_840, \
                         li_1114, li_1115, li_1116, li_1117, li_1118, li_1119, \
                         li_1120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = pa_x[k] * li_1114[k];

        t_1115[k] = pa_x[k] * li_1115[k];

        t_1116[k] = pa_x[k] * li_1116[k];

        t_1117[k] = pa_x[k] * li_1117[k];

        t_1118[k] = pa_x[k] * li_1118[k];

        t_1119[k] = pa_x[k] * li_1119[k];

        t_1120[k] = f_14 * lh_840[k]
                    + pa_x[k] * li_1120[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, t_1124, pa_x, pb_y, pb_z, lh_651, lh_672, \
                         lh_674, lh_843, li_1123, mh_840, mh_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_12 * lh_672[k]
                    + pb_y[k] * mh_840[k];

        t_1122[k] = f_12 * lh_651[k]
                    + pb_z[k] * mh_840[k];

        t_1123[k] = f_12 * lh_843[k]
                    + pa_x[k] * li_1123[k];

        t_1124[k] = f_12 * lh_674[k]
                    + pb_y[k] * mh_842[k];
    }

#pragma omp simd aligned(t_1125, t_1126, t_1127, t_1128, pa_x, pb_y, pb_z, lh_654, lh_677, \
                         lh_845, lh_846, li_1125, li_1126, mh_843, \
                         mh_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1125[k] = f_12 * lh_845[k]
                    + pa_x[k] * li_1125[k];

        t_1126[k] = f_11 * lh_846[k]
                    + pa_x[k] * li_1126[k];

        t_1127[k] = f_12 * lh_654[k]
                    + pb_z[k] * mh_843[k];

        t_1128[k] = f_12 * lh_677[k]
                    + pb_y[k] * mh_845[k];
    }

#pragma omp simd aligned(t_1129, t_1130, t_1131, t_1132, pa_x, pb_z, lh_657, lh_849, lh_850, \
                         lh_852, li_1129, li_1130, li_1132, mh_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1129[k] = f_11 * lh_849[k]
                    + pa_x[k] * li_1129[k];

        t_1130[k] = f_10 * lh_850[k]
                    + pa_x[k] * li_1130[k];

        t_1131[k] = f_12 * lh_657[k]
                    + pb_z[k] * mh_846[k];

        t_1132[k] = f_10 * lh_852[k]
                    + pa_x[k] * li_1132[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, t_1136, pa_x, pb_x, pb_y, lh_681, lh_854, \
                         lh_855, lh_856, li_1134, mh_849, mh_855, \
                         mh_856 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_12 * lh_681[k]
                    + pb_y[k] * mh_849[k];

        t_1134[k] = f_10 * lh_854[k]
                    + pa_x[k] * li_1134[k];

        t_1135[k] = f_9 * lh_855[k]
                    + pb_x[k] * mh_855[k];

        t_1136[k] = f_9 * lh_856[k]
                    + pb_x[k] * mh_856[k];
    }

#pragma omp simd aligned(t_1137, t_1138, t_1139, t_1140, t_1141, pa_x, pb_x, lh_857, lh_858, \
                         lh_859, lh_860, li_1141, mh_857, mh_858, mh_859, \
                         mh_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1137[k] = f_9 * lh_857[k]
                    + pb_x[k] * mh_857[k];

        t_1138[k] = f_9 * lh_858[k]
                    + pb_x[k] * mh_858[k];

        t_1139[k] = f_9 * lh_859[k]
                    + pb_x[k] * mh_859[k];

        t_1140[k] = f_9 * lh_860[k]
                    + pb_x[k] * mh_860[k];

        t_1141[k] = pa_x[k] * li_1141[k];
    }

#pragma omp simd aligned(t_1142, t_1143, t_1144, t_1145, t_1146, t_1147, t_1148, pa_x, lh_861, \
                         li_1142, li_1143, li_1144, li_1145, li_1146, li_1147, \
                         li_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = pa_x[k] * li_1142[k];

        t_1143[k] = pa_x[k] * li_1143[k];

        t_1144[k] = pa_x[k] * li_1144[k];

        t_1145[k] = pa_x[k] * li_1145[k];

        t_1146[k] = pa_x[k] * li_1146[k];

        t_1147[k] = pa_x[k] * li_1147[k];

        t_1148[k] = f_14 * lh_861[k]
                    + pa_x[k] * li_1148[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, t_1152, pa_x, pb_y, pb_z, lh_672, lh_693, \
                         lh_695, lh_864, li_1151, mh_861, mh_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_11 * lh_693[k]
                    + pb_y[k] * mh_861[k];

        t_1150[k] = f_26 * lh_672[k]
                    + pb_z[k] * mh_861[k];

        t_1151[k] = f_12 * lh_864[k]
                    + pa_x[k] * li_1151[k];

        t_1152[k] = f_11 * lh_695[k]
                    + pb_y[k] * mh_863[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, t_1156, pa_x, pb_y, pb_z, lh_675, lh_698, \
                         lh_866, lh_867, li_1153, li_1154, mh_864, \
                         mh_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_12 * lh_866[k]
                    + pa_x[k] * li_1153[k];

        t_1154[k] = f_11 * lh_867[k]
                    + pa_x[k] * li_1154[k];

        t_1155[k] = f_26 * lh_675[k]
                    + pb_z[k] * mh_864[k];

        t_1156[k] = f_11 * lh_698[k]
                    + pb_y[k] * mh_866[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, t_1160, pa_x, pb_z, lh_678, lh_870, lh_871, \
                         lh_873, li_1157, li_1158, li_1160, mh_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_11 * lh_870[k]
                    + pa_x[k] * li_1157[k];

        t_1158[k] = f_10 * lh_871[k]
                    + pa_x[k] * li_1158[k];

        t_1159[k] = f_26 * lh_678[k]
                    + pb_z[k] * mh_867[k];

        t_1160[k] = f_10 * lh_873[k]
                    + pa_x[k] * li_1160[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, t_1164, pa_x, pb_x, pb_y, lh_702, lh_875, \
                         lh_876, lh_877, li_1162, mh_870, mh_876, \
                         mh_877 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = f_11 * lh_702[k]
                    + pb_y[k] * mh_870[k];

        t_1162[k] = f_10 * lh_875[k]
                    + pa_x[k] * li_1162[k];

        t_1163[k] = f_9 * lh_876[k]
                    + pb_x[k] * mh_876[k];

        t_1164[k] = f_9 * lh_877[k]
                    + pb_x[k] * mh_877[k];
    }

#pragma omp simd aligned(t_1165, t_1166, t_1167, t_1168, t_1169, pa_x, pb_x, lh_878, lh_879, \
                         lh_880, lh_881, li_1169, mh_878, mh_879, mh_880, \
                         mh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1165[k] = f_9 * lh_878[k]
                    + pb_x[k] * mh_878[k];

        t_1166[k] = f_9 * lh_879[k]
                    + pb_x[k] * mh_879[k];

        t_1167[k] = f_9 * lh_880[k]
                    + pb_x[k] * mh_880[k];

        t_1168[k] = f_9 * lh_881[k]
                    + pb_x[k] * mh_881[k];

        t_1169[k] = pa_x[k] * li_1169[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, t_1173, t_1174, t_1175, t_1176, pa_x, lh_882, \
                         li_1170, li_1171, li_1172, li_1173, li_1174, li_1175, \
                         li_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = pa_x[k] * li_1170[k];

        t_1171[k] = pa_x[k] * li_1171[k];

        t_1172[k] = pa_x[k] * li_1172[k];

        t_1173[k] = pa_x[k] * li_1173[k];

        t_1174[k] = pa_x[k] * li_1174[k];

        t_1175[k] = pa_x[k] * li_1175[k];

        t_1176[k] = f_14 * lh_882[k]
                    + pa_x[k] * li_1176[k];
    }

#pragma omp simd aligned(t_1177, t_1178, t_1179, t_1180, pa_x, pb_y, pb_z, lh_693, lh_714, \
                         lh_716, lh_885, li_1179, mh_882, mh_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1177[k] = f_10 * lh_714[k]
                    + pb_y[k] * mh_882[k];

        t_1178[k] = f_14 * lh_693[k]
                    + pb_z[k] * mh_882[k];

        t_1179[k] = f_12 * lh_885[k]
                    + pa_x[k] * li_1179[k];

        t_1180[k] = f_10 * lh_716[k]
                    + pb_y[k] * mh_884[k];
    }

#pragma omp simd aligned(t_1181, t_1182, t_1183, t_1184, pa_x, pb_y, pb_z, lh_696, lh_719, \
                         lh_887, lh_888, li_1181, li_1182, mh_885, \
                         mh_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1181[k] = f_12 * lh_887[k]
                    + pa_x[k] * li_1181[k];

        t_1182[k] = f_11 * lh_888[k]
                    + pa_x[k] * li_1182[k];

        t_1183[k] = f_14 * lh_696[k]
                    + pb_z[k] * mh_885[k];

        t_1184[k] = f_10 * lh_719[k]
                    + pb_y[k] * mh_887[k];
    }
}

static auto
compute_prim_mi_electron_repulsion_0_piece9(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t lh,
                                            const size_t li, const size_t mg0, const size_t mg1,
                                            const size_t mh, const size_t ncols,
                                            const double alpha, const double beta,
                                            const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_13 = 4.0 / p;
    const auto f_14 = 3.0 / p;
    const auto f_17 = 3.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *lh_699 = buffer.data(lh + 699);
    const auto *lh_717 = buffer.data(lh + 717);
    const auto *lh_720 = buffer.data(lh + 720);
    const auto *lh_723 = buffer.data(lh + 723);
    const auto *lh_735 = buffer.data(lh + 735);
    const auto *lh_737 = buffer.data(lh + 737);
    const auto *lh_738 = buffer.data(lh + 738);
    const auto *lh_740 = buffer.data(lh + 740);
    const auto *lh_741 = buffer.data(lh + 741);
    const auto *lh_744 = buffer.data(lh + 744);
    const auto *lh_756 = buffer.data(lh + 756);
    const auto *lh_758 = buffer.data(lh + 758);
    const auto *lh_759 = buffer.data(lh + 759);
    const auto *lh_761 = buffer.data(lh + 761);
    const auto *lh_762 = buffer.data(lh + 762);
    const auto *lh_763 = buffer.data(lh + 763);
    const auto *lh_765 = buffer.data(lh + 765);
    const auto *lh_771 = buffer.data(lh + 771);
    const auto *lh_772 = buffer.data(lh + 772);
    const auto *lh_773 = buffer.data(lh + 773);
    const auto *lh_774 = buffer.data(lh + 774);
    const auto *lh_776 = buffer.data(lh + 776);
    const auto *lh_777 = buffer.data(lh + 777);
    const auto *lh_779 = buffer.data(lh + 779);
    const auto *lh_780 = buffer.data(lh + 780);
    const auto *lh_782 = buffer.data(lh + 782);
    const auto *lh_783 = buffer.data(lh + 783);
    const auto *lh_786 = buffer.data(lh + 786);
    const auto *lh_797 = buffer.data(lh + 797);
    const auto *lh_798 = buffer.data(lh + 798);
    const auto *lh_800 = buffer.data(lh + 800);
    const auto *lh_803 = buffer.data(lh + 803);
    const auto *lh_807 = buffer.data(lh + 807);
    const auto *lh_891 = buffer.data(lh + 891);
    const auto *lh_892 = buffer.data(lh + 892);
    const auto *lh_894 = buffer.data(lh + 894);
    const auto *lh_896 = buffer.data(lh + 896);
    const auto *lh_897 = buffer.data(lh + 897);
    const auto *lh_898 = buffer.data(lh + 898);
    const auto *lh_899 = buffer.data(lh + 899);
    const auto *lh_900 = buffer.data(lh + 900);
    const auto *lh_901 = buffer.data(lh + 901);
    const auto *lh_902 = buffer.data(lh + 902);
    const auto *lh_906 = buffer.data(lh + 906);
    const auto *lh_909 = buffer.data(lh + 909);
    const auto *lh_913 = buffer.data(lh + 913);
    const auto *lh_915 = buffer.data(lh + 915);
    const auto *lh_918 = buffer.data(lh + 918);
    const auto *lh_919 = buffer.data(lh + 919);
    const auto *lh_920 = buffer.data(lh + 920);
    const auto *lh_921 = buffer.data(lh + 921);
    const auto *lh_922 = buffer.data(lh + 922);
    const auto *lh_924 = buffer.data(lh + 924);
    const auto *lh_927 = buffer.data(lh + 927);
    const auto *lh_929 = buffer.data(lh + 929);
    const auto *lh_930 = buffer.data(lh + 930);
    const auto *lh_933 = buffer.data(lh + 933);
    const auto *lh_934 = buffer.data(lh + 934);
    const auto *lh_936 = buffer.data(lh + 936);
    const auto *lh_938 = buffer.data(lh + 938);
    const auto *lh_939 = buffer.data(lh + 939);
    const auto *lh_940 = buffer.data(lh + 940);
    const auto *lh_941 = buffer.data(lh + 941);
    const auto *lh_942 = buffer.data(lh + 942);
    const auto *lh_944 = buffer.data(lh + 944);

    const auto *li_980 = buffer.data(li + 980);
    const auto *li_982 = buffer.data(li + 982);
    const auto *li_985 = buffer.data(li + 985);
    const auto *li_989 = buffer.data(li + 989);
    const auto *li_994 = buffer.data(li + 994);
    const auto *li_1000 = buffer.data(li + 1000);
    const auto *li_1008 = buffer.data(li + 1008);
    const auto *li_1009 = buffer.data(li + 1009);
    const auto *li_1011 = buffer.data(li + 1011);
    const auto *li_1013 = buffer.data(li + 1013);
    const auto *li_1014 = buffer.data(li + 1014);
    const auto *li_1017 = buffer.data(li + 1017);
    const auto *li_1018 = buffer.data(li + 1018);
    const auto *li_1020 = buffer.data(li + 1020);
    const auto *li_1022 = buffer.data(li + 1022);
    const auto *li_1029 = buffer.data(li + 1029);
    const auto *li_1031 = buffer.data(li + 1031);
    const auto *li_1032 = buffer.data(li + 1032);
    const auto *li_1033 = buffer.data(li + 1033);
    const auto *li_1035 = buffer.data(li + 1035);
    const auto *li_1185 = buffer.data(li + 1185);
    const auto *li_1186 = buffer.data(li + 1186);
    const auto *li_1188 = buffer.data(li + 1188);
    const auto *li_1190 = buffer.data(li + 1190);
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
    const auto *li_1225 = buffer.data(li + 1225);
    const auto *li_1226 = buffer.data(li + 1226);
    const auto *li_1227 = buffer.data(li + 1227);
    const auto *li_1228 = buffer.data(li + 1228);
    const auto *li_1229 = buffer.data(li + 1229);
    const auto *li_1230 = buffer.data(li + 1230);
    const auto *li_1231 = buffer.data(li + 1231);
    const auto *li_1232 = buffer.data(li + 1232);
    const auto *li_1235 = buffer.data(li + 1235);
    const auto *li_1237 = buffer.data(li + 1237);
    const auto *li_1238 = buffer.data(li + 1238);
    const auto *li_1241 = buffer.data(li + 1241);
    const auto *li_1242 = buffer.data(li + 1242);
    const auto *li_1244 = buffer.data(li + 1244);
    const auto *li_1246 = buffer.data(li + 1246);
    const auto *li_1253 = buffer.data(li + 1253);
    const auto *li_1254 = buffer.data(li + 1254);
    const auto *li_1255 = buffer.data(li + 1255);
    const auto *li_1256 = buffer.data(li + 1256);
    const auto *li_1257 = buffer.data(li + 1257);
    const auto *li_1259 = buffer.data(li + 1259);

    const auto *mg0_675 = buffer.data(mg0 + 675);
    const auto *mg0_678 = buffer.data(mg0 + 678);
    const auto *mg0_680 = buffer.data(mg0 + 680);
    const auto *mg0_681 = buffer.data(mg0 + 681);
    const auto *mg0_684 = buffer.data(mg0 + 684);
    const auto *mg0_685 = buffer.data(mg0 + 685);
    const auto *mg0_686 = buffer.data(mg0 + 686);
    const auto *mg0_687 = buffer.data(mg0 + 687);
    const auto *mg0_689 = buffer.data(mg0 + 689);
    const auto *mg0_705 = buffer.data(mg0 + 705);
    const auto *mg0_708 = buffer.data(mg0 + 708);
    const auto *mg0_710 = buffer.data(mg0 + 710);
    const auto *mg0_711 = buffer.data(mg0 + 711);
    const auto *mg0_714 = buffer.data(mg0 + 714);
    const auto *mg0_715 = buffer.data(mg0 + 715);
    const auto *mg0_717 = buffer.data(mg0 + 717);
    const auto *mg0_719 = buffer.data(mg0 + 719);

    const auto *mg1_675 = buffer.data(mg1 + 675);
    const auto *mg1_678 = buffer.data(mg1 + 678);
    const auto *mg1_680 = buffer.data(mg1 + 680);
    const auto *mg1_681 = buffer.data(mg1 + 681);
    const auto *mg1_684 = buffer.data(mg1 + 684);
    const auto *mg1_685 = buffer.data(mg1 + 685);
    const auto *mg1_686 = buffer.data(mg1 + 686);
    const auto *mg1_687 = buffer.data(mg1 + 687);
    const auto *mg1_689 = buffer.data(mg1 + 689);
    const auto *mg1_705 = buffer.data(mg1 + 705);
    const auto *mg1_708 = buffer.data(mg1 + 708);
    const auto *mg1_710 = buffer.data(mg1 + 710);
    const auto *mg1_711 = buffer.data(mg1 + 711);
    const auto *mg1_714 = buffer.data(mg1 + 714);
    const auto *mg1_715 = buffer.data(mg1 + 715);
    const auto *mg1_717 = buffer.data(mg1 + 717);
    const auto *mg1_719 = buffer.data(mg1 + 719);

    const auto *mh_888 = buffer.data(mh + 888);
    const auto *mh_891 = buffer.data(mh + 891);
    const auto *mh_897 = buffer.data(mh + 897);
    const auto *mh_898 = buffer.data(mh + 898);
    const auto *mh_899 = buffer.data(mh + 899);
    const auto *mh_900 = buffer.data(mh + 900);
    const auto *mh_901 = buffer.data(mh + 901);
    const auto *mh_902 = buffer.data(mh + 902);
    const auto *mh_903 = buffer.data(mh + 903);
    const auto *mh_905 = buffer.data(mh + 905);
    const auto *mh_906 = buffer.data(mh + 906);
    const auto *mh_908 = buffer.data(mh + 908);
    const auto *mh_909 = buffer.data(mh + 909);
    const auto *mh_912 = buffer.data(mh + 912);
    const auto *mh_918 = buffer.data(mh + 918);
    const auto *mh_919 = buffer.data(mh + 919);
    const auto *mh_920 = buffer.data(mh + 920);
    const auto *mh_921 = buffer.data(mh + 921);
    const auto *mh_922 = buffer.data(mh + 922);
    const auto *mh_924 = buffer.data(mh + 924);
    const auto *mh_926 = buffer.data(mh + 926);
    const auto *mh_927 = buffer.data(mh + 927);
    const auto *mh_929 = buffer.data(mh + 929);
    const auto *mh_930 = buffer.data(mh + 930);
    const auto *mh_933 = buffer.data(mh + 933);
    const auto *mh_938 = buffer.data(mh + 938);
    const auto *mh_939 = buffer.data(mh + 939);
    const auto *mh_940 = buffer.data(mh + 940);
    const auto *mh_941 = buffer.data(mh + 941);
    const auto *mh_942 = buffer.data(mh + 942);
    const auto *mh_944 = buffer.data(mh + 944);
    const auto *mh_945 = buffer.data(mh + 945);
    const auto *mh_946 = buffer.data(mh + 946);
    const auto *mh_948 = buffer.data(mh + 948);
    const auto *mh_950 = buffer.data(mh + 950);
    const auto *mh_951 = buffer.data(mh + 951);
    const auto *mh_954 = buffer.data(mh + 954);
    const auto *mh_955 = buffer.data(mh + 955);
    const auto *mh_957 = buffer.data(mh + 957);
    const auto *mh_959 = buffer.data(mh + 959);
    const auto *mh_960 = buffer.data(mh + 960);
    const auto *mh_961 = buffer.data(mh + 961);
    const auto *mh_962 = buffer.data(mh + 962);
    const auto *mh_963 = buffer.data(mh + 963);
    const auto *mh_964 = buffer.data(mh + 964);
    const auto *mh_965 = buffer.data(mh + 965);
    const auto *mh_966 = buffer.data(mh + 966);
    const auto *mh_968 = buffer.data(mh + 968);
    const auto *mh_969 = buffer.data(mh + 969);
    const auto *mh_971 = buffer.data(mh + 971);
    const auto *mh_972 = buffer.data(mh + 972);
    const auto *mh_975 = buffer.data(mh + 975);
    const auto *mh_981 = buffer.data(mh + 981);
    const auto *mh_982 = buffer.data(mh + 982);
    const auto *mh_983 = buffer.data(mh + 983);
    const auto *mh_984 = buffer.data(mh + 984);
    const auto *mh_985 = buffer.data(mh + 985);
    const auto *mh_986 = buffer.data(mh + 986);
    const auto *mh_987 = buffer.data(mh + 987);
    const auto *mh_989 = buffer.data(mh + 989);
    const auto *mh_990 = buffer.data(mh + 990);
    const auto *mh_992 = buffer.data(mh + 992);
    const auto *mh_993 = buffer.data(mh + 993);
    const auto *mh_996 = buffer.data(mh + 996);
    const auto *mh_997 = buffer.data(mh + 997);
    const auto *mh_999 = buffer.data(mh + 999);
    const auto *mh_1001 = buffer.data(mh + 1001);
    const auto *mh_1002 = buffer.data(mh + 1002);

#pragma omp simd aligned(t_1185, t_1186, t_1187, t_1188, pa_x, pb_z, lh_699, lh_891, lh_892, \
                         lh_894, li_1185, li_1186, li_1188, mh_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_11 * lh_891[k]
                    + pa_x[k] * li_1185[k];

        t_1186[k] = f_10 * lh_892[k]
                    + pa_x[k] * li_1186[k];

        t_1187[k] = f_14 * lh_699[k]
                    + pb_z[k] * mh_888[k];

        t_1188[k] = f_10 * lh_894[k]
                    + pa_x[k] * li_1188[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, t_1192, pa_x, pb_x, pb_y, lh_723, lh_896, \
                         lh_897, lh_898, li_1190, mh_891, mh_897, \
                         mh_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_10 * lh_723[k]
                    + pb_y[k] * mh_891[k];

        t_1190[k] = f_10 * lh_896[k]
                    + pa_x[k] * li_1190[k];

        t_1191[k] = f_9 * lh_897[k]
                    + pb_x[k] * mh_897[k];

        t_1192[k] = f_9 * lh_898[k]
                    + pb_x[k] * mh_898[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, t_1196, t_1197, pa_x, pb_x, lh_899, lh_900, \
                         lh_901, lh_902, li_1197, mh_899, mh_900, mh_901, \
                         mh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = f_9 * lh_899[k]
                    + pb_x[k] * mh_899[k];

        t_1194[k] = f_9 * lh_900[k]
                    + pb_x[k] * mh_900[k];

        t_1195[k] = f_9 * lh_901[k]
                    + pb_x[k] * mh_901[k];

        t_1196[k] = f_9 * lh_902[k]
                    + pb_x[k] * mh_902[k];

        t_1197[k] = pa_x[k] * li_1197[k];
    }

#pragma omp simd aligned(t_1198, t_1199, t_1200, t_1201, t_1202, t_1203, t_1204, pa_x, pa_y, \
                         li_980, li_1198, li_1199, li_1200, li_1201, li_1202, \
                         li_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1198[k] = pa_x[k] * li_1198[k];

        t_1199[k] = pa_x[k] * li_1199[k];

        t_1200[k] = pa_x[k] * li_1200[k];

        t_1201[k] = pa_x[k] * li_1201[k];

        t_1202[k] = pa_x[k] * li_1202[k];

        t_1203[k] = pa_x[k] * li_1203[k];

        t_1204[k] = pa_y[k] * li_980[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, t_1208, t_1209, pa_x, pa_y, pb_y, lh_735, \
                         lh_737, lh_906, li_982, li_985, li_1207, mh_903, \
                         mh_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = f_9 * lh_735[k]
                    + pb_y[k] * mh_903[k];

        t_1206[k] = pa_y[k] * li_982[k];

        t_1207[k] = f_12 * lh_906[k]
                    + pa_x[k] * li_1207[k];

        t_1208[k] = f_9 * lh_737[k]
                    + pb_y[k] * mh_905[k];

        t_1209[k] = pa_y[k] * li_985[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, pa_x, pa_y, pb_y, pb_z, lh_717, \
                         lh_740, lh_909, li_989, li_1210, mh_906, \
                         mh_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = f_11 * lh_909[k]
                    + pa_x[k] * li_1210[k];

        t_1211[k] = f_17 * lh_717[k]
                    + pb_z[k] * mh_906[k];

        t_1212[k] = f_9 * lh_740[k]
                    + pb_y[k] * mh_908[k];

        t_1213[k] = pa_y[k] * li_989[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, t_1217, pa_x, pb_y, pb_z, lh_720, lh_744, \
                         lh_913, lh_915, li_1214, li_1216, mh_909, \
                         mh_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = f_10 * lh_913[k]
                    + pa_x[k] * li_1214[k];

        t_1215[k] = f_17 * lh_720[k]
                    + pb_z[k] * mh_909[k];

        t_1216[k] = f_10 * lh_915[k]
                    + pa_x[k] * li_1216[k];

        t_1217[k] = f_9 * lh_744[k]
                    + pb_y[k] * mh_912[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, t_1221, t_1222, pa_y, pb_x, lh_918, lh_919, \
                         lh_920, lh_921, li_994, mh_918, mh_919, mh_920, \
                         mh_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = pa_y[k] * li_994[k];

        t_1219[k] = f_9 * lh_918[k]
                    + pb_x[k] * mh_918[k];

        t_1220[k] = f_9 * lh_919[k]
                    + pb_x[k] * mh_919[k];

        t_1221[k] = f_9 * lh_920[k]
                    + pb_x[k] * mh_920[k];

        t_1222[k] = f_9 * lh_921[k]
                    + pb_x[k] * mh_921[k];
    }

#pragma omp simd aligned(t_1223, t_1224, t_1225, t_1226, t_1227, t_1228, pa_x, pa_y, pb_x, \
                         lh_922, li_1000, li_1225, li_1226, li_1227, li_1228, \
                         mh_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1223[k] = f_9 * lh_922[k]
                    + pb_x[k] * mh_922[k];

        t_1224[k] = pa_y[k] * li_1000[k];

        t_1225[k] = pa_x[k] * li_1225[k];

        t_1226[k] = pa_x[k] * li_1226[k];

        t_1227[k] = pa_x[k] * li_1227[k];

        t_1228[k] = pa_x[k] * li_1228[k];
    }

#pragma omp simd aligned(t_1229, t_1230, t_1231, t_1232, t_1233, t_1234, pa_x, pb_y, pb_z, \
                         lh_735, lh_924, li_1229, li_1230, li_1231, li_1232, \
                         mh_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1229[k] = pa_x[k] * li_1229[k];

        t_1230[k] = pa_x[k] * li_1230[k];

        t_1231[k] = pa_x[k] * li_1231[k];

        t_1232[k] = f_14 * lh_924[k]
                    + pa_x[k] * li_1232[k];

        t_1233[k] = pb_y[k] * mh_924[k];

        t_1234[k] = f_13 * lh_735[k]
                    + pb_z[k] * mh_924[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, t_1238, pa_x, pb_y, lh_927, lh_929, lh_930, \
                         li_1235, li_1237, li_1238, mh_926 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_12 * lh_927[k]
                    + pa_x[k] * li_1235[k];

        t_1236[k] = pb_y[k] * mh_926[k];

        t_1237[k] = f_12 * lh_929[k]
                    + pa_x[k] * li_1237[k];

        t_1238[k] = f_11 * lh_930[k]
                    + pa_x[k] * li_1238[k];
    }

#pragma omp simd aligned(t_1239, t_1240, t_1241, t_1242, pa_x, pb_y, pb_z, lh_738, lh_933, \
                         lh_934, li_1241, li_1242, mh_927, mh_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1239[k] = f_13 * lh_738[k]
                    + pb_z[k] * mh_927[k];

        t_1240[k] = pb_y[k] * mh_929[k];

        t_1241[k] = f_11 * lh_933[k]
                    + pa_x[k] * li_1241[k];

        t_1242[k] = f_10 * lh_934[k]
                    + pa_x[k] * li_1242[k];
    }

#pragma omp simd aligned(t_1243, t_1244, t_1245, t_1246, pa_x, pb_y, pb_z, lh_741, lh_936, \
                         lh_938, li_1244, li_1246, mh_930, mh_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1243[k] = f_13 * lh_741[k]
                    + pb_z[k] * mh_930[k];

        t_1244[k] = f_10 * lh_936[k]
                    + pa_x[k] * li_1244[k];

        t_1245[k] = pb_y[k] * mh_933[k];

        t_1246[k] = f_10 * lh_938[k]
                    + pa_x[k] * li_1246[k];
    }

#pragma omp simd aligned(t_1247, t_1248, t_1249, t_1250, t_1251, pb_x, pb_y, lh_939, lh_940, \
                         lh_941, lh_942, mh_938, mh_939, mh_940, mh_941, \
                         mh_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1247[k] = f_9 * lh_939[k]
                    + pb_x[k] * mh_939[k];

        t_1248[k] = f_9 * lh_940[k]
                    + pb_x[k] * mh_940[k];

        t_1249[k] = f_9 * lh_941[k]
                    + pb_x[k] * mh_941[k];

        t_1250[k] = f_9 * lh_942[k]
                    + pb_x[k] * mh_942[k];

        t_1251[k] = pb_y[k] * mh_938[k];
    }

#pragma omp simd aligned(t_1252, t_1253, t_1254, t_1255, t_1256, t_1257, pa_x, pb_x, lh_944, \
                         li_1253, li_1254, li_1255, li_1256, li_1257, \
                         mh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1252[k] = f_9 * lh_944[k]
                    + pb_x[k] * mh_944[k];

        t_1253[k] = pa_x[k] * li_1253[k];

        t_1254[k] = pa_x[k] * li_1254[k];

        t_1255[k] = pa_x[k] * li_1255[k];

        t_1256[k] = pa_x[k] * li_1256[k];

        t_1257[k] = pa_x[k] * li_1257[k];
    }

#pragma omp simd aligned(t_1258, t_1259, t_1260, t_1261, t_1262, pa_x, pb_x, pb_y, pb_z, \
                         lh_756, li_1259, mg0_675, mg1_675, mh_944, \
                         mh_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1258[k] = pb_y[k] * mh_944[k];

        t_1259[k] = pa_x[k] * li_1259[k];

        t_1260[k] = f_1 * mg0_675[k]
                    - f_2 * mg1_675[k]
                    + pb_x[k] * mh_945[k];

        t_1261[k] = f_0 * lh_756[k]
                    + pb_y[k] * mh_945[k];

        t_1262[k] = pb_z[k] * mh_945[k];
    }

#pragma omp simd aligned(t_1263, t_1264, t_1265, t_1266, pb_x, pb_z, mg0_678, mg0_680, \
                         mg0_681, mg1_678, mg1_680, mg1_681, mh_946, mh_948, mh_950, \
                         mh_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1263[k] = f_7 * mg0_678[k]
                    - f_8 * mg1_678[k]
                    + pb_x[k] * mh_948[k];

        t_1264[k] = pb_z[k] * mh_946[k];

        t_1265[k] = f_7 * mg0_680[k]
                    - f_8 * mg1_680[k]
                    + pb_x[k] * mh_950[k];

        t_1266[k] = f_5 * mg0_681[k]
                    - f_6 * mg1_681[k]
                    + pb_x[k] * mh_951[k];
    }

#pragma omp simd aligned(t_1267, t_1268, t_1269, t_1270, pb_x, pb_y, pb_z, lh_761, mg0_684, \
                         mg0_685, mg1_684, mg1_685, mh_948, mh_950, mh_954, \
                         mh_955 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1267[k] = pb_z[k] * mh_948[k];

        t_1268[k] = f_0 * lh_761[k]
                    + pb_y[k] * mh_950[k];

        t_1269[k] = f_5 * mg0_684[k]
                    - f_6 * mg1_684[k]
                    + pb_x[k] * mh_954[k];

        t_1270[k] = f_3 * mg0_685[k]
                    - f_4 * mg1_685[k]
                    + pb_x[k] * mh_955[k];
    }

#pragma omp simd aligned(t_1271, t_1272, t_1273, t_1274, pb_x, pb_y, pb_z, lh_765, mg0_687, \
                         mg0_689, mg1_687, mg1_689, mh_951, mh_954, mh_957, \
                         mh_959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1271[k] = pb_z[k] * mh_951[k];

        t_1272[k] = f_3 * mg0_687[k]
                    - f_4 * mg1_687[k]
                    + pb_x[k] * mh_957[k];

        t_1273[k] = f_0 * lh_765[k]
                    + pb_y[k] * mh_954[k];

        t_1274[k] = f_3 * mg0_689[k]
                    - f_4 * mg1_689[k]
                    + pb_x[k] * mh_959[k];
    }

#pragma omp simd aligned(t_1275, t_1276, t_1277, t_1278, t_1279, t_1280, pb_x, mh_960, mh_961, \
                         mh_962, mh_963, mh_964, mh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1275[k] = pb_x[k] * mh_960[k];

        t_1276[k] = pb_x[k] * mh_961[k];

        t_1277[k] = pb_x[k] * mh_962[k];

        t_1278[k] = pb_x[k] * mh_963[k];

        t_1279[k] = pb_x[k] * mh_964[k];

        t_1280[k] = pb_x[k] * mh_965[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, t_1284, pb_y, pb_z, lh_771, mg0_685, mg0_686, \
                         mg1_685, mg1_686, mh_960, mh_961, mh_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = f_0 * lh_771[k]
                    + f_1 * mg0_685[k]
                    - f_2 * mg1_685[k]
                    + pb_y[k] * mh_960[k];

        t_1282[k] = pb_z[k] * mh_960[k];

        t_1283[k] = f_3 * mg0_685[k]
                    - f_4 * mg1_685[k]
                    + pb_z[k] * mh_961[k];

        t_1284[k] = f_5 * mg0_686[k]
                    - f_6 * mg1_686[k]
                    + pb_z[k] * mh_962[k];
    }

#pragma omp simd aligned(t_1285, t_1286, t_1287, t_1288, pa_z, pb_y, pb_z, lh_776, li_1008, \
                         mg0_687, mg0_689, mg1_687, mg1_689, mh_963, \
                         mh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1285[k] = f_7 * mg0_687[k]
                    - f_8 * mg1_687[k]
                    + pb_z[k] * mh_963[k];

        t_1286[k] = f_0 * lh_776[k]
                    + pb_y[k] * mh_965[k];

        t_1287[k] = f_1 * mg0_689[k]
                    - f_2 * mg1_689[k]
                    + pb_z[k] * mh_965[k];

        t_1288[k] = pa_z[k] * li_1008[k];
    }

#pragma omp simd aligned(t_1289, t_1290, t_1291, t_1292, t_1293, pa_z, pb_y, pb_z, lh_756, \
                         lh_758, lh_779, li_1009, li_1011, li_1013, mh_966, \
                         mh_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1289[k] = pa_z[k] * li_1009[k];

        t_1290[k] = f_9 * lh_756[k]
                    + pb_z[k] * mh_966[k];

        t_1291[k] = pa_z[k] * li_1011[k];

        t_1292[k] = f_13 * lh_779[k]
                    + pb_y[k] * mh_968[k];

        t_1293[k] = f_10 * lh_758[k]
                    + pa_z[k] * li_1013[k];
    }

#pragma omp simd aligned(t_1294, t_1295, t_1296, t_1297, t_1298, pa_z, pb_y, pb_z, lh_759, \
                         lh_761, lh_782, li_1014, li_1017, li_1018, mh_969, \
                         mh_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1294[k] = pa_z[k] * li_1014[k];

        t_1295[k] = f_9 * lh_759[k]
                    + pb_z[k] * mh_969[k];

        t_1296[k] = f_13 * lh_782[k]
                    + pb_y[k] * mh_971[k];

        t_1297[k] = f_11 * lh_761[k]
                    + pa_z[k] * li_1017[k];

        t_1298[k] = pa_z[k] * li_1018[k];
    }

#pragma omp simd aligned(t_1299, t_1300, t_1301, t_1302, pa_z, pb_y, pb_z, lh_762, lh_763, \
                         lh_765, lh_786, li_1020, li_1022, mh_972, \
                         mh_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1299[k] = f_9 * lh_762[k]
                    + pb_z[k] * mh_972[k];

        t_1300[k] = f_10 * lh_763[k]
                    + pa_z[k] * li_1020[k];

        t_1301[k] = f_13 * lh_786[k]
                    + pb_y[k] * mh_975[k];

        t_1302[k] = f_12 * lh_765[k]
                    + pa_z[k] * li_1022[k];
    }

#pragma omp simd aligned(t_1303, t_1304, t_1305, t_1306, t_1307, t_1308, t_1309, pa_z, pb_x, \
                         li_1029, mh_981, mh_982, mh_983, mh_984, mh_985, \
                         mh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1303[k] = pb_x[k] * mh_981[k];

        t_1304[k] = pb_x[k] * mh_982[k];

        t_1305[k] = pb_x[k] * mh_983[k];

        t_1306[k] = pb_x[k] * mh_984[k];

        t_1307[k] = pb_x[k] * mh_985[k];

        t_1308[k] = pb_x[k] * mh_986[k];

        t_1309[k] = pa_z[k] * li_1029[k];
    }

#pragma omp simd aligned(t_1310, t_1311, t_1312, t_1313, pa_z, pb_z, lh_771, lh_772, lh_773, \
                         lh_774, li_1031, li_1032, li_1033, mh_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1310[k] = f_9 * lh_771[k]
                    + pb_z[k] * mh_981[k];

        t_1311[k] = f_10 * lh_772[k]
                    + pa_z[k] * li_1031[k];

        t_1312[k] = f_11 * lh_773[k]
                    + pa_z[k] * li_1032[k];

        t_1313[k] = f_12 * lh_774[k]
                    + pa_z[k] * li_1033[k];
    }

#pragma omp simd aligned(t_1314, t_1315, t_1316, t_1317, pa_z, pb_x, pb_y, lh_776, lh_797, \
                         lh_798, li_1035, mg0_705, mg1_705, mh_986, \
                         mh_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1314[k] = f_13 * lh_797[k]
                    + pb_y[k] * mh_986[k];

        t_1315[k] = f_14 * lh_776[k]
                    + pa_z[k] * li_1035[k];

        t_1316[k] = f_1 * mg0_705[k]
                    - f_2 * mg1_705[k]
                    + pb_x[k] * mh_987[k];

        t_1317[k] = f_17 * lh_798[k]
                    + pb_y[k] * mh_987[k];
    }

#pragma omp simd aligned(t_1318, t_1319, t_1320, pb_x, pb_y, pb_z, lh_777, lh_800, mg0_708, \
                         mg1_708, mh_987, mh_989, mh_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1318[k] = f_10 * lh_777[k]
                    + pb_z[k] * mh_987[k];

        t_1319[k] = f_7 * mg0_708[k]
                    - f_8 * mg1_708[k]
                    + pb_x[k] * mh_990[k];

        t_1320[k] = f_17 * lh_800[k]
                    + pb_y[k] * mh_989[k];
    }

#pragma omp simd aligned(t_1321, t_1322, t_1323, t_1324, pb_x, pb_y, pb_z, lh_780, lh_803, \
                         mg0_710, mg0_711, mg1_710, mg1_711, mh_990, mh_992, \
                         mh_993 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1321[k] = f_7 * mg0_710[k]
                    - f_8 * mg1_710[k]
                    + pb_x[k] * mh_992[k];

        t_1322[k] = f_5 * mg0_711[k]
                    - f_6 * mg1_711[k]
                    + pb_x[k] * mh_993[k];

        t_1323[k] = f_10 * lh_780[k]
                    + pb_z[k] * mh_990[k];

        t_1324[k] = f_17 * lh_803[k]
                    + pb_y[k] * mh_992[k];
    }

#pragma omp simd aligned(t_1325, t_1326, t_1327, pb_x, pb_z, lh_783, mg0_714, mg0_715, \
                         mg1_714, mg1_715, mh_993, mh_996, mh_997 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1325[k] = f_5 * mg0_714[k]
                    - f_6 * mg1_714[k]
                    + pb_x[k] * mh_996[k];

        t_1326[k] = f_3 * mg0_715[k]
                    - f_4 * mg1_715[k]
                    + pb_x[k] * mh_997[k];

        t_1327[k] = f_10 * lh_783[k]
                    + pb_z[k] * mh_993[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, t_1331, pb_x, pb_y, lh_807, mg0_717, mg0_719, \
                         mg1_717, mg1_719, mh_996, mh_999, mh_1001, \
                         mh_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = f_3 * mg0_717[k]
                    - f_4 * mg1_717[k]
                    + pb_x[k] * mh_999[k];

        t_1329[k] = f_17 * lh_807[k]
                    + pb_y[k] * mh_996[k];

        t_1330[k] = f_3 * mg0_719[k]
                    - f_4 * mg1_719[k]
                    + pb_x[k] * mh_1001[k];

        t_1331[k] = pb_x[k] * mh_1002[k];
    }
}

static auto
compute_prim_mi_electron_repulsion_0_piece10(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t ki0,
                                             const size_t ki1, const size_t lh, const size_t li,
                                             const size_t mg0, const size_t mg1, const size_t mh,
                                             const size_t ncols, const double alpha,
                                             const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_14 = 3.0 / p;
    const auto f_15 = 0.5 / alpha;
    const auto f_16 = 0.5 * beta / (alpha * p);
    const auto f_17 = 3.5 / p;
    const auto f_18 = 3.0 / alpha;
    const auto f_19 = 3.0 * beta / (alpha * p);
    const auto f_20 = 1.0 / alpha;
    const auto f_21 = beta / (alpha * p);
    const auto f_22 = 2.5 / alpha;
    const auto f_23 = 2.5 * beta / (alpha * p);
    const auto f_24 = 1.5 / alpha;
    const auto f_25 = 1.5 * beta / (alpha * p);
    const auto f_26 = 2.5 / p;
    const auto f_27 = 2.0 / alpha;
    const auto f_28 = 2.0 * beta / (alpha * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ki0_805 = buffer.data(ki0 + 805);
    const auto *ki0_833 = buffer.data(ki0 + 833);
    const auto *ki0_861 = buffer.data(ki0 + 861);
    const auto *ki0_867 = buffer.data(ki0 + 867);
    const auto *ki0_889 = buffer.data(ki0 + 889);
    const auto *ki0_895 = buffer.data(ki0 + 895);
    const auto *ki0_917 = buffer.data(ki0 + 917);
    const auto *ki0_923 = buffer.data(ki0 + 923);
    const auto *ki0_951 = buffer.data(ki0 + 951);
    const auto *ki0_979 = buffer.data(ki0 + 979);

    const auto *ki1_805 = buffer.data(ki1 + 805);
    const auto *ki1_833 = buffer.data(ki1 + 833);
    const auto *ki1_861 = buffer.data(ki1 + 861);
    const auto *ki1_867 = buffer.data(ki1 + 867);
    const auto *ki1_889 = buffer.data(ki1 + 889);
    const auto *ki1_895 = buffer.data(ki1 + 895);
    const auto *ki1_917 = buffer.data(ki1 + 917);
    const auto *ki1_923 = buffer.data(ki1 + 923);
    const auto *ki1_951 = buffer.data(ki1 + 951);
    const auto *ki1_979 = buffer.data(ki1 + 979);

    const auto *lh_792 = buffer.data(lh + 792);
    const auto *lh_798 = buffer.data(lh + 798);
    const auto *lh_801 = buffer.data(lh + 801);
    const auto *lh_804 = buffer.data(lh + 804);
    const auto *lh_813 = buffer.data(lh + 813);
    const auto *lh_815 = buffer.data(lh + 815);
    const auto *lh_816 = buffer.data(lh + 816);
    const auto *lh_817 = buffer.data(lh + 817);
    const auto *lh_818 = buffer.data(lh + 818);
    const auto *lh_819 = buffer.data(lh + 819);
    const auto *lh_821 = buffer.data(lh + 821);
    const auto *lh_822 = buffer.data(lh + 822);
    const auto *lh_824 = buffer.data(lh + 824);
    const auto *lh_825 = buffer.data(lh + 825);
    const auto *lh_828 = buffer.data(lh + 828);
    const auto *lh_834 = buffer.data(lh + 834);
    const auto *lh_836 = buffer.data(lh + 836);
    const auto *lh_837 = buffer.data(lh + 837);
    const auto *lh_838 = buffer.data(lh + 838);
    const auto *lh_839 = buffer.data(lh + 839);
    const auto *lh_840 = buffer.data(lh + 840);
    const auto *lh_842 = buffer.data(lh + 842);
    const auto *lh_843 = buffer.data(lh + 843);
    const auto *lh_845 = buffer.data(lh + 845);
    const auto *lh_846 = buffer.data(lh + 846);
    const auto *lh_849 = buffer.data(lh + 849);
    const auto *lh_855 = buffer.data(lh + 855);
    const auto *lh_857 = buffer.data(lh + 857);
    const auto *lh_858 = buffer.data(lh + 858);
    const auto *lh_859 = buffer.data(lh + 859);
    const auto *lh_860 = buffer.data(lh + 860);
    const auto *lh_861 = buffer.data(lh + 861);
    const auto *lh_863 = buffer.data(lh + 863);
    const auto *lh_864 = buffer.data(lh + 864);
    const auto *lh_866 = buffer.data(lh + 866);
    const auto *lh_867 = buffer.data(lh + 867);
    const auto *lh_870 = buffer.data(lh + 870);
    const auto *lh_876 = buffer.data(lh + 876);
    const auto *lh_878 = buffer.data(lh + 878);
    const auto *lh_879 = buffer.data(lh + 879);
    const auto *lh_880 = buffer.data(lh + 880);
    const auto *lh_881 = buffer.data(lh + 881);
    const auto *lh_882 = buffer.data(lh + 882);
    const auto *lh_884 = buffer.data(lh + 884);
    const auto *lh_887 = buffer.data(lh + 887);
    const auto *lh_891 = buffer.data(lh + 891);
    const auto *lh_899 = buffer.data(lh + 899);
    const auto *lh_900 = buffer.data(lh + 900);
    const auto *lh_901 = buffer.data(lh + 901);
    const auto *lh_902 = buffer.data(lh + 902);

    const auto *li_1057 = buffer.data(li + 1057);
    const auto *li_1085 = buffer.data(li + 1085);
    const auto *li_1091 = buffer.data(li + 1091);
    const auto *li_1113 = buffer.data(li + 1113);
    const auto *li_1119 = buffer.data(li + 1119);
    const auto *li_1141 = buffer.data(li + 1141);
    const auto *li_1147 = buffer.data(li + 1147);
    const auto *li_1169 = buffer.data(li + 1169);
    const auto *li_1175 = buffer.data(li + 1175);
    const auto *li_1203 = buffer.data(li + 1203);

    const auto *mg0_717 = buffer.data(mg0 + 717);
    const auto *mg0_718 = buffer.data(mg0 + 718);
    const auto *mg0_719 = buffer.data(mg0 + 719);
    const auto *mg0_720 = buffer.data(mg0 + 720);
    const auto *mg0_723 = buffer.data(mg0 + 723);
    const auto *mg0_725 = buffer.data(mg0 + 725);
    const auto *mg0_726 = buffer.data(mg0 + 726);
    const auto *mg0_729 = buffer.data(mg0 + 729);
    const auto *mg0_730 = buffer.data(mg0 + 730);
    const auto *mg0_732 = buffer.data(mg0 + 732);
    const auto *mg0_733 = buffer.data(mg0 + 733);
    const auto *mg0_734 = buffer.data(mg0 + 734);
    const auto *mg0_735 = buffer.data(mg0 + 735);
    const auto *mg0_738 = buffer.data(mg0 + 738);
    const auto *mg0_740 = buffer.data(mg0 + 740);
    const auto *mg0_741 = buffer.data(mg0 + 741);
    const auto *mg0_744 = buffer.data(mg0 + 744);
    const auto *mg0_745 = buffer.data(mg0 + 745);
    const auto *mg0_747 = buffer.data(mg0 + 747);
    const auto *mg0_748 = buffer.data(mg0 + 748);
    const auto *mg0_749 = buffer.data(mg0 + 749);
    const auto *mg0_750 = buffer.data(mg0 + 750);
    const auto *mg0_753 = buffer.data(mg0 + 753);
    const auto *mg0_755 = buffer.data(mg0 + 755);
    const auto *mg0_756 = buffer.data(mg0 + 756);
    const auto *mg0_759 = buffer.data(mg0 + 759);
    const auto *mg0_760 = buffer.data(mg0 + 760);
    const auto *mg0_762 = buffer.data(mg0 + 762);
    const auto *mg0_763 = buffer.data(mg0 + 763);
    const auto *mg0_764 = buffer.data(mg0 + 764);
    const auto *mg0_765 = buffer.data(mg0 + 765);
    const auto *mg0_768 = buffer.data(mg0 + 768);
    const auto *mg0_770 = buffer.data(mg0 + 770);
    const auto *mg0_771 = buffer.data(mg0 + 771);
    const auto *mg0_774 = buffer.data(mg0 + 774);
    const auto *mg0_775 = buffer.data(mg0 + 775);
    const auto *mg0_777 = buffer.data(mg0 + 777);
    const auto *mg0_778 = buffer.data(mg0 + 778);
    const auto *mg0_779 = buffer.data(mg0 + 779);

    const auto *mg1_717 = buffer.data(mg1 + 717);
    const auto *mg1_718 = buffer.data(mg1 + 718);
    const auto *mg1_719 = buffer.data(mg1 + 719);
    const auto *mg1_720 = buffer.data(mg1 + 720);
    const auto *mg1_723 = buffer.data(mg1 + 723);
    const auto *mg1_725 = buffer.data(mg1 + 725);
    const auto *mg1_726 = buffer.data(mg1 + 726);
    const auto *mg1_729 = buffer.data(mg1 + 729);
    const auto *mg1_730 = buffer.data(mg1 + 730);
    const auto *mg1_732 = buffer.data(mg1 + 732);
    const auto *mg1_733 = buffer.data(mg1 + 733);
    const auto *mg1_734 = buffer.data(mg1 + 734);
    const auto *mg1_735 = buffer.data(mg1 + 735);
    const auto *mg1_738 = buffer.data(mg1 + 738);
    const auto *mg1_740 = buffer.data(mg1 + 740);
    const auto *mg1_741 = buffer.data(mg1 + 741);
    const auto *mg1_744 = buffer.data(mg1 + 744);
    const auto *mg1_745 = buffer.data(mg1 + 745);
    const auto *mg1_747 = buffer.data(mg1 + 747);
    const auto *mg1_748 = buffer.data(mg1 + 748);
    const auto *mg1_749 = buffer.data(mg1 + 749);
    const auto *mg1_750 = buffer.data(mg1 + 750);
    const auto *mg1_753 = buffer.data(mg1 + 753);
    const auto *mg1_755 = buffer.data(mg1 + 755);
    const auto *mg1_756 = buffer.data(mg1 + 756);
    const auto *mg1_759 = buffer.data(mg1 + 759);
    const auto *mg1_760 = buffer.data(mg1 + 760);
    const auto *mg1_762 = buffer.data(mg1 + 762);
    const auto *mg1_763 = buffer.data(mg1 + 763);
    const auto *mg1_764 = buffer.data(mg1 + 764);
    const auto *mg1_765 = buffer.data(mg1 + 765);
    const auto *mg1_768 = buffer.data(mg1 + 768);
    const auto *mg1_770 = buffer.data(mg1 + 770);
    const auto *mg1_771 = buffer.data(mg1 + 771);
    const auto *mg1_774 = buffer.data(mg1 + 774);
    const auto *mg1_775 = buffer.data(mg1 + 775);
    const auto *mg1_777 = buffer.data(mg1 + 777);
    const auto *mg1_778 = buffer.data(mg1 + 778);
    const auto *mg1_779 = buffer.data(mg1 + 779);

    const auto *mh_1002 = buffer.data(mh + 1002);
    const auto *mh_1003 = buffer.data(mh + 1003);
    const auto *mh_1004 = buffer.data(mh + 1004);
    const auto *mh_1005 = buffer.data(mh + 1005);
    const auto *mh_1006 = buffer.data(mh + 1006);
    const auto *mh_1007 = buffer.data(mh + 1007);
    const auto *mh_1008 = buffer.data(mh + 1008);
    const auto *mh_1010 = buffer.data(mh + 1010);
    const auto *mh_1011 = buffer.data(mh + 1011);
    const auto *mh_1013 = buffer.data(mh + 1013);
    const auto *mh_1014 = buffer.data(mh + 1014);
    const auto *mh_1017 = buffer.data(mh + 1017);
    const auto *mh_1018 = buffer.data(mh + 1018);
    const auto *mh_1020 = buffer.data(mh + 1020);
    const auto *mh_1022 = buffer.data(mh + 1022);
    const auto *mh_1023 = buffer.data(mh + 1023);
    const auto *mh_1024 = buffer.data(mh + 1024);
    const auto *mh_1025 = buffer.data(mh + 1025);
    const auto *mh_1026 = buffer.data(mh + 1026);
    const auto *mh_1027 = buffer.data(mh + 1027);
    const auto *mh_1028 = buffer.data(mh + 1028);
    const auto *mh_1029 = buffer.data(mh + 1029);
    const auto *mh_1031 = buffer.data(mh + 1031);
    const auto *mh_1032 = buffer.data(mh + 1032);
    const auto *mh_1034 = buffer.data(mh + 1034);
    const auto *mh_1035 = buffer.data(mh + 1035);
    const auto *mh_1038 = buffer.data(mh + 1038);
    const auto *mh_1039 = buffer.data(mh + 1039);
    const auto *mh_1041 = buffer.data(mh + 1041);
    const auto *mh_1043 = buffer.data(mh + 1043);
    const auto *mh_1044 = buffer.data(mh + 1044);
    const auto *mh_1045 = buffer.data(mh + 1045);
    const auto *mh_1046 = buffer.data(mh + 1046);
    const auto *mh_1047 = buffer.data(mh + 1047);
    const auto *mh_1048 = buffer.data(mh + 1048);
    const auto *mh_1049 = buffer.data(mh + 1049);
    const auto *mh_1050 = buffer.data(mh + 1050);
    const auto *mh_1052 = buffer.data(mh + 1052);
    const auto *mh_1053 = buffer.data(mh + 1053);
    const auto *mh_1055 = buffer.data(mh + 1055);
    const auto *mh_1056 = buffer.data(mh + 1056);
    const auto *mh_1059 = buffer.data(mh + 1059);
    const auto *mh_1060 = buffer.data(mh + 1060);
    const auto *mh_1062 = buffer.data(mh + 1062);
    const auto *mh_1064 = buffer.data(mh + 1064);
    const auto *mh_1065 = buffer.data(mh + 1065);
    const auto *mh_1066 = buffer.data(mh + 1066);
    const auto *mh_1067 = buffer.data(mh + 1067);
    const auto *mh_1068 = buffer.data(mh + 1068);
    const auto *mh_1069 = buffer.data(mh + 1069);
    const auto *mh_1070 = buffer.data(mh + 1070);
    const auto *mh_1071 = buffer.data(mh + 1071);
    const auto *mh_1073 = buffer.data(mh + 1073);
    const auto *mh_1074 = buffer.data(mh + 1074);
    const auto *mh_1076 = buffer.data(mh + 1076);
    const auto *mh_1077 = buffer.data(mh + 1077);
    const auto *mh_1080 = buffer.data(mh + 1080);
    const auto *mh_1081 = buffer.data(mh + 1081);
    const auto *mh_1083 = buffer.data(mh + 1083);
    const auto *mh_1085 = buffer.data(mh + 1085);
    const auto *mh_1086 = buffer.data(mh + 1086);
    const auto *mh_1087 = buffer.data(mh + 1087);
    const auto *mh_1088 = buffer.data(mh + 1088);
    const auto *mh_1089 = buffer.data(mh + 1089);
    const auto *mh_1090 = buffer.data(mh + 1090);
    const auto *mh_1091 = buffer.data(mh + 1091);

#pragma omp simd aligned(t_1332, t_1333, t_1334, t_1335, t_1336, t_1337, pa_z, pb_x, ki0_805, \
                         ki1_805, li_1057, mh_1003, mh_1004, mh_1005, mh_1006, \
                         mh_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1332[k] = pb_x[k] * mh_1003[k];

        t_1333[k] = pb_x[k] * mh_1004[k];

        t_1334[k] = pb_x[k] * mh_1005[k];

        t_1335[k] = pb_x[k] * mh_1006[k];

        t_1336[k] = pb_x[k] * mh_1007[k];

        t_1337[k] = f_15 * ki0_805[k]
                    - f_16 * ki1_805[k]
                    + pa_z[k] * li_1057[k];
    }

#pragma omp simd aligned(t_1338, t_1339, t_1340, pb_y, pb_z, lh_792, lh_815, lh_816, mg0_717, \
                         mg0_718, mg1_717, mg1_718, mh_1002, mh_1004, \
                         mh_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1338[k] = f_10 * lh_792[k]
                    + pb_z[k] * mh_1002[k];

        t_1339[k] = f_17 * lh_815[k]
                    + f_7 * mg0_717[k]
                    - f_8 * mg1_717[k]
                    + pb_y[k] * mh_1004[k];

        t_1340[k] = f_17 * lh_816[k]
                    + f_5 * mg0_718[k]
                    - f_6 * mg1_718[k]
                    + pb_y[k] * mh_1005[k];
    }

#pragma omp simd aligned(t_1341, t_1342, t_1343, pa_y, pb_y, ki0_867, ki1_867, lh_817, lh_818, \
                         li_1091, mg0_719, mg1_719, mh_1006, mh_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1341[k] = f_17 * lh_817[k]
                    + f_3 * mg0_719[k]
                    - f_4 * mg1_719[k]
                    + pb_y[k] * mh_1006[k];

        t_1342[k] = f_17 * lh_818[k]
                    + pb_y[k] * mh_1007[k];

        t_1343[k] = f_18 * ki0_867[k]
                    - f_19 * ki1_867[k]
                    + pa_y[k] * li_1091[k];
    }

#pragma omp simd aligned(t_1344, t_1345, t_1346, t_1347, pb_x, pb_y, pb_z, lh_798, lh_819, \
                         mg0_720, mg0_723, mg1_720, mg1_723, mh_1008, \
                         mh_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1344[k] = f_1 * mg0_720[k]
                    - f_2 * mg1_720[k]
                    + pb_x[k] * mh_1008[k];

        t_1345[k] = f_14 * lh_819[k]
                    + pb_y[k] * mh_1008[k];

        t_1346[k] = f_11 * lh_798[k]
                    + pb_z[k] * mh_1008[k];

        t_1347[k] = f_7 * mg0_723[k]
                    - f_8 * mg1_723[k]
                    + pb_x[k] * mh_1011[k];
    }

#pragma omp simd aligned(t_1348, t_1349, t_1350, pb_x, pb_y, lh_821, mg0_725, mg0_726, \
                         mg1_725, mg1_726, mh_1010, mh_1013, mh_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1348[k] = f_14 * lh_821[k]
                    + pb_y[k] * mh_1010[k];

        t_1349[k] = f_7 * mg0_725[k]
                    - f_8 * mg1_725[k]
                    + pb_x[k] * mh_1013[k];

        t_1350[k] = f_5 * mg0_726[k]
                    - f_6 * mg1_726[k]
                    + pb_x[k] * mh_1014[k];
    }

#pragma omp simd aligned(t_1351, t_1352, t_1353, pb_x, pb_y, pb_z, lh_801, lh_824, mg0_729, \
                         mg1_729, mh_1011, mh_1013, mh_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1351[k] = f_11 * lh_801[k]
                    + pb_z[k] * mh_1011[k];

        t_1352[k] = f_14 * lh_824[k]
                    + pb_y[k] * mh_1013[k];

        t_1353[k] = f_5 * mg0_729[k]
                    - f_6 * mg1_729[k]
                    + pb_x[k] * mh_1017[k];
    }

#pragma omp simd aligned(t_1354, t_1355, t_1356, pb_x, pb_z, lh_804, mg0_730, mg0_732, \
                         mg1_730, mg1_732, mh_1014, mh_1018, mh_1020 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1354[k] = f_3 * mg0_730[k]
                    - f_4 * mg1_730[k]
                    + pb_x[k] * mh_1018[k];

        t_1355[k] = f_11 * lh_804[k]
                    + pb_z[k] * mh_1014[k];

        t_1356[k] = f_3 * mg0_732[k]
                    - f_4 * mg1_732[k]
                    + pb_x[k] * mh_1020[k];
    }

#pragma omp simd aligned(t_1357, t_1358, t_1359, t_1360, t_1361, pb_x, pb_y, lh_828, mg0_734, \
                         mg1_734, mh_1017, mh_1022, mh_1023, mh_1024, \
                         mh_1025 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1357[k] = f_14 * lh_828[k]
                    + pb_y[k] * mh_1017[k];

        t_1358[k] = f_3 * mg0_734[k]
                    - f_4 * mg1_734[k]
                    + pb_x[k] * mh_1022[k];

        t_1359[k] = pb_x[k] * mh_1023[k];

        t_1360[k] = pb_x[k] * mh_1024[k];

        t_1361[k] = pb_x[k] * mh_1025[k];
    }

#pragma omp simd aligned(t_1362, t_1363, t_1364, t_1365, t_1366, pa_z, pb_x, pb_z, ki0_833, \
                         ki1_833, lh_813, li_1085, mh_1023, mh_1026, mh_1027, \
                         mh_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1362[k] = pb_x[k] * mh_1026[k];

        t_1363[k] = pb_x[k] * mh_1027[k];

        t_1364[k] = pb_x[k] * mh_1028[k];

        t_1365[k] = f_20 * ki0_833[k]
                    - f_21 * ki1_833[k]
                    + pa_z[k] * li_1085[k];

        t_1366[k] = f_11 * lh_813[k]
                    + pb_z[k] * mh_1023[k];
    }

#pragma omp simd aligned(t_1367, t_1368, t_1369, pb_y, lh_836, lh_837, lh_838, mg0_732, \
                         mg0_733, mg0_734, mg1_732, mg1_733, mg1_734, mh_1025, mh_1026, \
                         mh_1027 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1367[k] = f_14 * lh_836[k]
                    + f_7 * mg0_732[k]
                    - f_8 * mg1_732[k]
                    + pb_y[k] * mh_1025[k];

        t_1368[k] = f_14 * lh_837[k]
                    + f_5 * mg0_733[k]
                    - f_6 * mg1_733[k]
                    + pb_y[k] * mh_1026[k];

        t_1369[k] = f_14 * lh_838[k]
                    + f_3 * mg0_734[k]
                    - f_4 * mg1_734[k]
                    + pb_y[k] * mh_1027[k];
    }

#pragma omp simd aligned(t_1370, t_1371, t_1372, t_1373, pa_y, pb_x, pb_y, ki0_895, ki1_895, \
                         lh_839, lh_840, li_1119, mg0_735, mg1_735, mh_1028, \
                         mh_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1370[k] = f_14 * lh_839[k]
                    + pb_y[k] * mh_1028[k];

        t_1371[k] = f_22 * ki0_895[k]
                    - f_23 * ki1_895[k]
                    + pa_y[k] * li_1119[k];

        t_1372[k] = f_1 * mg0_735[k]
                    - f_2 * mg1_735[k]
                    + pb_x[k] * mh_1029[k];

        t_1373[k] = f_26 * lh_840[k]
                    + pb_y[k] * mh_1029[k];
    }

#pragma omp simd aligned(t_1374, t_1375, t_1376, pb_x, pb_y, pb_z, lh_819, lh_842, mg0_738, \
                         mg1_738, mh_1029, mh_1031, mh_1032 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1374[k] = f_12 * lh_819[k]
                    + pb_z[k] * mh_1029[k];

        t_1375[k] = f_7 * mg0_738[k]
                    - f_8 * mg1_738[k]
                    + pb_x[k] * mh_1032[k];

        t_1376[k] = f_26 * lh_842[k]
                    + pb_y[k] * mh_1031[k];
    }

#pragma omp simd aligned(t_1377, t_1378, t_1379, t_1380, pb_x, pb_y, pb_z, lh_822, lh_845, \
                         mg0_740, mg0_741, mg1_740, mg1_741, mh_1032, mh_1034, \
                         mh_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1377[k] = f_7 * mg0_740[k]
                    - f_8 * mg1_740[k]
                    + pb_x[k] * mh_1034[k];

        t_1378[k] = f_5 * mg0_741[k]
                    - f_6 * mg1_741[k]
                    + pb_x[k] * mh_1035[k];

        t_1379[k] = f_12 * lh_822[k]
                    + pb_z[k] * mh_1032[k];

        t_1380[k] = f_26 * lh_845[k]
                    + pb_y[k] * mh_1034[k];
    }

#pragma omp simd aligned(t_1381, t_1382, t_1383, pb_x, pb_z, lh_825, mg0_744, mg0_745, \
                         mg1_744, mg1_745, mh_1035, mh_1038, mh_1039 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1381[k] = f_5 * mg0_744[k]
                    - f_6 * mg1_744[k]
                    + pb_x[k] * mh_1038[k];

        t_1382[k] = f_3 * mg0_745[k]
                    - f_4 * mg1_745[k]
                    + pb_x[k] * mh_1039[k];

        t_1383[k] = f_12 * lh_825[k]
                    + pb_z[k] * mh_1035[k];
    }

#pragma omp simd aligned(t_1384, t_1385, t_1386, t_1387, pb_x, pb_y, lh_849, mg0_747, mg0_749, \
                         mg1_747, mg1_749, mh_1038, mh_1041, mh_1043, \
                         mh_1044 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1384[k] = f_3 * mg0_747[k]
                    - f_4 * mg1_747[k]
                    + pb_x[k] * mh_1041[k];

        t_1385[k] = f_26 * lh_849[k]
                    + pb_y[k] * mh_1038[k];

        t_1386[k] = f_3 * mg0_749[k]
                    - f_4 * mg1_749[k]
                    + pb_x[k] * mh_1043[k];

        t_1387[k] = pb_x[k] * mh_1044[k];
    }

#pragma omp simd aligned(t_1388, t_1389, t_1390, t_1391, t_1392, t_1393, pa_z, pb_x, ki0_861, \
                         ki1_861, li_1113, mh_1045, mh_1046, mh_1047, mh_1048, \
                         mh_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1388[k] = pb_x[k] * mh_1045[k];

        t_1389[k] = pb_x[k] * mh_1046[k];

        t_1390[k] = pb_x[k] * mh_1047[k];

        t_1391[k] = pb_x[k] * mh_1048[k];

        t_1392[k] = pb_x[k] * mh_1049[k];

        t_1393[k] = f_24 * ki0_861[k]
                    - f_25 * ki1_861[k]
                    + pa_z[k] * li_1113[k];
    }

#pragma omp simd aligned(t_1394, t_1395, t_1396, pb_y, pb_z, lh_834, lh_857, lh_858, mg0_747, \
                         mg0_748, mg1_747, mg1_748, mh_1044, mh_1046, \
                         mh_1047 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1394[k] = f_12 * lh_834[k]
                    + pb_z[k] * mh_1044[k];

        t_1395[k] = f_26 * lh_857[k]
                    + f_7 * mg0_747[k]
                    - f_8 * mg1_747[k]
                    + pb_y[k] * mh_1046[k];

        t_1396[k] = f_26 * lh_858[k]
                    + f_5 * mg0_748[k]
                    - f_6 * mg1_748[k]
                    + pb_y[k] * mh_1047[k];
    }

#pragma omp simd aligned(t_1397, t_1398, t_1399, pa_y, pb_y, ki0_923, ki1_923, lh_859, lh_860, \
                         li_1147, mg0_749, mg1_749, mh_1048, mh_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1397[k] = f_26 * lh_859[k]
                    + f_3 * mg0_749[k]
                    - f_4 * mg1_749[k]
                    + pb_y[k] * mh_1048[k];

        t_1398[k] = f_26 * lh_860[k]
                    + pb_y[k] * mh_1049[k];

        t_1399[k] = f_27 * ki0_923[k]
                    - f_28 * ki1_923[k]
                    + pa_y[k] * li_1147[k];
    }

#pragma omp simd aligned(t_1400, t_1401, t_1402, t_1403, pb_x, pb_y, pb_z, lh_840, lh_861, \
                         mg0_750, mg0_753, mg1_750, mg1_753, mh_1050, \
                         mh_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1400[k] = f_1 * mg0_750[k]
                    - f_2 * mg1_750[k]
                    + pb_x[k] * mh_1050[k];

        t_1401[k] = f_12 * lh_861[k]
                    + pb_y[k] * mh_1050[k];

        t_1402[k] = f_26 * lh_840[k]
                    + pb_z[k] * mh_1050[k];

        t_1403[k] = f_7 * mg0_753[k]
                    - f_8 * mg1_753[k]
                    + pb_x[k] * mh_1053[k];
    }

#pragma omp simd aligned(t_1404, t_1405, t_1406, pb_x, pb_y, lh_863, mg0_755, mg0_756, \
                         mg1_755, mg1_756, mh_1052, mh_1055, mh_1056 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1404[k] = f_12 * lh_863[k]
                    + pb_y[k] * mh_1052[k];

        t_1405[k] = f_7 * mg0_755[k]
                    - f_8 * mg1_755[k]
                    + pb_x[k] * mh_1055[k];

        t_1406[k] = f_5 * mg0_756[k]
                    - f_6 * mg1_756[k]
                    + pb_x[k] * mh_1056[k];
    }

#pragma omp simd aligned(t_1407, t_1408, t_1409, pb_x, pb_y, pb_z, lh_843, lh_866, mg0_759, \
                         mg1_759, mh_1053, mh_1055, mh_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1407[k] = f_26 * lh_843[k]
                    + pb_z[k] * mh_1053[k];

        t_1408[k] = f_12 * lh_866[k]
                    + pb_y[k] * mh_1055[k];

        t_1409[k] = f_5 * mg0_759[k]
                    - f_6 * mg1_759[k]
                    + pb_x[k] * mh_1059[k];
    }

#pragma omp simd aligned(t_1410, t_1411, t_1412, pb_x, pb_z, lh_846, mg0_760, mg0_762, \
                         mg1_760, mg1_762, mh_1056, mh_1060, mh_1062 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1410[k] = f_3 * mg0_760[k]
                    - f_4 * mg1_760[k]
                    + pb_x[k] * mh_1060[k];

        t_1411[k] = f_26 * lh_846[k]
                    + pb_z[k] * mh_1056[k];

        t_1412[k] = f_3 * mg0_762[k]
                    - f_4 * mg1_762[k]
                    + pb_x[k] * mh_1062[k];
    }

#pragma omp simd aligned(t_1413, t_1414, t_1415, t_1416, t_1417, pb_x, pb_y, lh_870, mg0_764, \
                         mg1_764, mh_1059, mh_1064, mh_1065, mh_1066, \
                         mh_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1413[k] = f_12 * lh_870[k]
                    + pb_y[k] * mh_1059[k];

        t_1414[k] = f_3 * mg0_764[k]
                    - f_4 * mg1_764[k]
                    + pb_x[k] * mh_1064[k];

        t_1415[k] = pb_x[k] * mh_1065[k];

        t_1416[k] = pb_x[k] * mh_1066[k];

        t_1417[k] = pb_x[k] * mh_1067[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, t_1421, t_1422, pa_z, pb_x, pb_z, ki0_889, \
                         ki1_889, lh_855, li_1141, mh_1065, mh_1068, mh_1069, \
                         mh_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = pb_x[k] * mh_1068[k];

        t_1419[k] = pb_x[k] * mh_1069[k];

        t_1420[k] = pb_x[k] * mh_1070[k];

        t_1421[k] = f_27 * ki0_889[k]
                    - f_28 * ki1_889[k]
                    + pa_z[k] * li_1141[k];

        t_1422[k] = f_26 * lh_855[k]
                    + pb_z[k] * mh_1065[k];
    }

#pragma omp simd aligned(t_1423, t_1424, t_1425, pb_y, lh_878, lh_879, lh_880, mg0_762, \
                         mg0_763, mg0_764, mg1_762, mg1_763, mg1_764, mh_1067, mh_1068, \
                         mh_1069 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1423[k] = f_12 * lh_878[k]
                    + f_7 * mg0_762[k]
                    - f_8 * mg1_762[k]
                    + pb_y[k] * mh_1067[k];

        t_1424[k] = f_12 * lh_879[k]
                    + f_5 * mg0_763[k]
                    - f_6 * mg1_763[k]
                    + pb_y[k] * mh_1068[k];

        t_1425[k] = f_12 * lh_880[k]
                    + f_3 * mg0_764[k]
                    - f_4 * mg1_764[k]
                    + pb_y[k] * mh_1069[k];
    }

#pragma omp simd aligned(t_1426, t_1427, t_1428, t_1429, pa_y, pb_x, pb_y, ki0_951, ki1_951, \
                         lh_881, lh_882, li_1175, mg0_765, mg1_765, mh_1070, \
                         mh_1071 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1426[k] = f_12 * lh_881[k]
                    + pb_y[k] * mh_1070[k];

        t_1427[k] = f_24 * ki0_951[k]
                    - f_25 * ki1_951[k]
                    + pa_y[k] * li_1175[k];

        t_1428[k] = f_1 * mg0_765[k]
                    - f_2 * mg1_765[k]
                    + pb_x[k] * mh_1071[k];

        t_1429[k] = f_11 * lh_882[k]
                    + pb_y[k] * mh_1071[k];
    }

#pragma omp simd aligned(t_1430, t_1431, t_1432, pb_x, pb_y, pb_z, lh_861, lh_884, mg0_768, \
                         mg1_768, mh_1071, mh_1073, mh_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1430[k] = f_14 * lh_861[k]
                    + pb_z[k] * mh_1071[k];

        t_1431[k] = f_7 * mg0_768[k]
                    - f_8 * mg1_768[k]
                    + pb_x[k] * mh_1074[k];

        t_1432[k] = f_11 * lh_884[k]
                    + pb_y[k] * mh_1073[k];
    }

#pragma omp simd aligned(t_1433, t_1434, t_1435, t_1436, pb_x, pb_y, pb_z, lh_864, lh_887, \
                         mg0_770, mg0_771, mg1_770, mg1_771, mh_1074, mh_1076, \
                         mh_1077 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1433[k] = f_7 * mg0_770[k]
                    - f_8 * mg1_770[k]
                    + pb_x[k] * mh_1076[k];

        t_1434[k] = f_5 * mg0_771[k]
                    - f_6 * mg1_771[k]
                    + pb_x[k] * mh_1077[k];

        t_1435[k] = f_14 * lh_864[k]
                    + pb_z[k] * mh_1074[k];

        t_1436[k] = f_11 * lh_887[k]
                    + pb_y[k] * mh_1076[k];
    }

#pragma omp simd aligned(t_1437, t_1438, t_1439, pb_x, pb_z, lh_867, mg0_774, mg0_775, \
                         mg1_774, mg1_775, mh_1077, mh_1080, mh_1081 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1437[k] = f_5 * mg0_774[k]
                    - f_6 * mg1_774[k]
                    + pb_x[k] * mh_1080[k];

        t_1438[k] = f_3 * mg0_775[k]
                    - f_4 * mg1_775[k]
                    + pb_x[k] * mh_1081[k];

        t_1439[k] = f_14 * lh_867[k]
                    + pb_z[k] * mh_1077[k];
    }

#pragma omp simd aligned(t_1440, t_1441, t_1442, t_1443, pb_x, pb_y, lh_891, mg0_777, mg0_779, \
                         mg1_777, mg1_779, mh_1080, mh_1083, mh_1085, \
                         mh_1086 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1440[k] = f_3 * mg0_777[k]
                    - f_4 * mg1_777[k]
                    + pb_x[k] * mh_1083[k];

        t_1441[k] = f_11 * lh_891[k]
                    + pb_y[k] * mh_1080[k];

        t_1442[k] = f_3 * mg0_779[k]
                    - f_4 * mg1_779[k]
                    + pb_x[k] * mh_1085[k];

        t_1443[k] = pb_x[k] * mh_1086[k];
    }

#pragma omp simd aligned(t_1444, t_1445, t_1446, t_1447, t_1448, t_1449, pa_z, pb_x, ki0_917, \
                         ki1_917, li_1169, mh_1087, mh_1088, mh_1089, mh_1090, \
                         mh_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1444[k] = pb_x[k] * mh_1087[k];

        t_1445[k] = pb_x[k] * mh_1088[k];

        t_1446[k] = pb_x[k] * mh_1089[k];

        t_1447[k] = pb_x[k] * mh_1090[k];

        t_1448[k] = pb_x[k] * mh_1091[k];

        t_1449[k] = f_22 * ki0_917[k]
                    - f_23 * ki1_917[k]
                    + pa_z[k] * li_1169[k];
    }

#pragma omp simd aligned(t_1450, t_1451, t_1452, pb_y, pb_z, lh_876, lh_899, lh_900, mg0_777, \
                         mg0_778, mg1_777, mg1_778, mh_1086, mh_1088, \
                         mh_1089 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1450[k] = f_14 * lh_876[k]
                    + pb_z[k] * mh_1086[k];

        t_1451[k] = f_11 * lh_899[k]
                    + f_7 * mg0_777[k]
                    - f_8 * mg1_777[k]
                    + pb_y[k] * mh_1088[k];

        t_1452[k] = f_11 * lh_900[k]
                    + f_5 * mg0_778[k]
                    - f_6 * mg1_778[k]
                    + pb_y[k] * mh_1089[k];
    }

#pragma omp simd aligned(t_1453, t_1454, t_1455, pa_y, pb_y, ki0_979, ki1_979, lh_901, lh_902, \
                         li_1203, mg0_779, mg1_779, mh_1090, mh_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1453[k] = f_11 * lh_901[k]
                    + f_3 * mg0_779[k]
                    - f_4 * mg1_779[k]
                    + pb_y[k] * mh_1090[k];

        t_1454[k] = f_11 * lh_902[k]
                    + pb_y[k] * mh_1091[k];

        t_1455[k] = f_20 * ki0_979[k]
                    - f_21 * ki1_979[k]
                    + pa_y[k] * li_1203[k];
    }
}

static auto
compute_prim_mi_electron_repulsion_0_piece11(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t ki0,
                                             const size_t ki1, const size_t lh, const size_t li,
                                             const size_t mg0, const size_t mg1, const size_t mh,
                                             const size_t ncols, const double alpha,
                                             const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_13 = 4.0 / p;
    const auto f_14 = 3.0 / p;
    const auto f_15 = 0.5 / alpha;
    const auto f_16 = 0.5 * beta / (alpha * p);
    const auto f_17 = 3.5 / p;
    const auto f_18 = 3.0 / alpha;
    const auto f_19 = 3.0 * beta / (alpha * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ki0_945 = buffer.data(ki0 + 945);
    const auto *ki0_1007 = buffer.data(ki0 + 1007);

    const auto *ki1_945 = buffer.data(ki1 + 945);
    const auto *ki1_1007 = buffer.data(ki1 + 1007);

    const auto *lh_882 = buffer.data(lh + 882);
    const auto *lh_885 = buffer.data(lh + 885);
    const auto *lh_888 = buffer.data(lh + 888);
    const auto *lh_897 = buffer.data(lh + 897);
    const auto *lh_903 = buffer.data(lh + 903);
    const auto *lh_905 = buffer.data(lh + 905);
    const auto *lh_906 = buffer.data(lh + 906);
    const auto *lh_908 = buffer.data(lh + 908);
    const auto *lh_909 = buffer.data(lh + 909);
    const auto *lh_912 = buffer.data(lh + 912);
    const auto *lh_918 = buffer.data(lh + 918);
    const auto *lh_920 = buffer.data(lh + 920);
    const auto *lh_921 = buffer.data(lh + 921);
    const auto *lh_922 = buffer.data(lh + 922);
    const auto *lh_923 = buffer.data(lh + 923);
    const auto *lh_924 = buffer.data(lh + 924);
    const auto *lh_925 = buffer.data(lh + 925);
    const auto *lh_926 = buffer.data(lh + 926);
    const auto *lh_927 = buffer.data(lh + 927);
    const auto *lh_929 = buffer.data(lh + 929);
    const auto *lh_930 = buffer.data(lh + 930);
    const auto *lh_932 = buffer.data(lh + 932);
    const auto *lh_933 = buffer.data(lh + 933);
    const auto *lh_939 = buffer.data(lh + 939);
    const auto *lh_941 = buffer.data(lh + 941);
    const auto *lh_942 = buffer.data(lh + 942);
    const auto *lh_943 = buffer.data(lh + 943);
    const auto *lh_944 = buffer.data(lh + 944);

    const auto *li_1197 = buffer.data(li + 1197);
    const auto *li_1231 = buffer.data(li + 1231);
    const auto *li_1232 = buffer.data(li + 1232);
    const auto *li_1234 = buffer.data(li + 1234);
    const auto *li_1235 = buffer.data(li + 1235);
    const auto *li_1237 = buffer.data(li + 1237);
    const auto *li_1238 = buffer.data(li + 1238);
    const auto *li_1241 = buffer.data(li + 1241);
    const auto *li_1242 = buffer.data(li + 1242);
    const auto *li_1244 = buffer.data(li + 1244);
    const auto *li_1246 = buffer.data(li + 1246);
    const auto *li_1253 = buffer.data(li + 1253);
    const auto *li_1255 = buffer.data(li + 1255);
    const auto *li_1256 = buffer.data(li + 1256);
    const auto *li_1257 = buffer.data(li + 1257);
    const auto *li_1259 = buffer.data(li + 1259);

    const auto *mg0_780 = buffer.data(mg0 + 780);
    const auto *mg0_783 = buffer.data(mg0 + 783);
    const auto *mg0_785 = buffer.data(mg0 + 785);
    const auto *mg0_786 = buffer.data(mg0 + 786);
    const auto *mg0_789 = buffer.data(mg0 + 789);
    const auto *mg0_790 = buffer.data(mg0 + 790);
    const auto *mg0_792 = buffer.data(mg0 + 792);
    const auto *mg0_793 = buffer.data(mg0 + 793);
    const auto *mg0_794 = buffer.data(mg0 + 794);
    const auto *mg0_810 = buffer.data(mg0 + 810);
    const auto *mg0_813 = buffer.data(mg0 + 813);
    const auto *mg0_815 = buffer.data(mg0 + 815);
    const auto *mg0_816 = buffer.data(mg0 + 816);
    const auto *mg0_819 = buffer.data(mg0 + 819);
    const auto *mg0_820 = buffer.data(mg0 + 820);
    const auto *mg0_822 = buffer.data(mg0 + 822);
    const auto *mg0_823 = buffer.data(mg0 + 823);
    const auto *mg0_824 = buffer.data(mg0 + 824);

    const auto *mg1_780 = buffer.data(mg1 + 780);
    const auto *mg1_783 = buffer.data(mg1 + 783);
    const auto *mg1_785 = buffer.data(mg1 + 785);
    const auto *mg1_786 = buffer.data(mg1 + 786);
    const auto *mg1_789 = buffer.data(mg1 + 789);
    const auto *mg1_790 = buffer.data(mg1 + 790);
    const auto *mg1_792 = buffer.data(mg1 + 792);
    const auto *mg1_793 = buffer.data(mg1 + 793);
    const auto *mg1_794 = buffer.data(mg1 + 794);
    const auto *mg1_810 = buffer.data(mg1 + 810);
    const auto *mg1_813 = buffer.data(mg1 + 813);
    const auto *mg1_815 = buffer.data(mg1 + 815);
    const auto *mg1_816 = buffer.data(mg1 + 816);
    const auto *mg1_819 = buffer.data(mg1 + 819);
    const auto *mg1_820 = buffer.data(mg1 + 820);
    const auto *mg1_822 = buffer.data(mg1 + 822);
    const auto *mg1_823 = buffer.data(mg1 + 823);
    const auto *mg1_824 = buffer.data(mg1 + 824);

    const auto *mh_1092 = buffer.data(mh + 1092);
    const auto *mh_1094 = buffer.data(mh + 1094);
    const auto *mh_1095 = buffer.data(mh + 1095);
    const auto *mh_1097 = buffer.data(mh + 1097);
    const auto *mh_1098 = buffer.data(mh + 1098);
    const auto *mh_1101 = buffer.data(mh + 1101);
    const auto *mh_1102 = buffer.data(mh + 1102);
    const auto *mh_1104 = buffer.data(mh + 1104);
    const auto *mh_1106 = buffer.data(mh + 1106);
    const auto *mh_1107 = buffer.data(mh + 1107);
    const auto *mh_1108 = buffer.data(mh + 1108);
    const auto *mh_1109 = buffer.data(mh + 1109);
    const auto *mh_1110 = buffer.data(mh + 1110);
    const auto *mh_1111 = buffer.data(mh + 1111);
    const auto *mh_1112 = buffer.data(mh + 1112);
    const auto *mh_1113 = buffer.data(mh + 1113);
    const auto *mh_1115 = buffer.data(mh + 1115);
    const auto *mh_1116 = buffer.data(mh + 1116);
    const auto *mh_1118 = buffer.data(mh + 1118);
    const auto *mh_1119 = buffer.data(mh + 1119);
    const auto *mh_1122 = buffer.data(mh + 1122);
    const auto *mh_1128 = buffer.data(mh + 1128);
    const auto *mh_1129 = buffer.data(mh + 1129);
    const auto *mh_1130 = buffer.data(mh + 1130);
    const auto *mh_1131 = buffer.data(mh + 1131);
    const auto *mh_1132 = buffer.data(mh + 1132);
    const auto *mh_1133 = buffer.data(mh + 1133);
    const auto *mh_1134 = buffer.data(mh + 1134);
    const auto *mh_1136 = buffer.data(mh + 1136);
    const auto *mh_1137 = buffer.data(mh + 1137);
    const auto *mh_1139 = buffer.data(mh + 1139);
    const auto *mh_1140 = buffer.data(mh + 1140);
    const auto *mh_1143 = buffer.data(mh + 1143);
    const auto *mh_1144 = buffer.data(mh + 1144);
    const auto *mh_1146 = buffer.data(mh + 1146);
    const auto *mh_1148 = buffer.data(mh + 1148);
    const auto *mh_1149 = buffer.data(mh + 1149);
    const auto *mh_1150 = buffer.data(mh + 1150);
    const auto *mh_1151 = buffer.data(mh + 1151);
    const auto *mh_1152 = buffer.data(mh + 1152);
    const auto *mh_1153 = buffer.data(mh + 1153);
    const auto *mh_1154 = buffer.data(mh + 1154);

#pragma omp simd aligned(t_1456, t_1457, t_1458, t_1459, pb_x, pb_y, pb_z, lh_882, lh_903, \
                         mg0_780, mg0_783, mg1_780, mg1_783, mh_1092, \
                         mh_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1456[k] = f_1 * mg0_780[k]
                    - f_2 * mg1_780[k]
                    + pb_x[k] * mh_1092[k];

        t_1457[k] = f_10 * lh_903[k]
                    + pb_y[k] * mh_1092[k];

        t_1458[k] = f_17 * lh_882[k]
                    + pb_z[k] * mh_1092[k];

        t_1459[k] = f_7 * mg0_783[k]
                    - f_8 * mg1_783[k]
                    + pb_x[k] * mh_1095[k];
    }

#pragma omp simd aligned(t_1460, t_1461, t_1462, pb_x, pb_y, lh_905, mg0_785, mg0_786, \
                         mg1_785, mg1_786, mh_1094, mh_1097, mh_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1460[k] = f_10 * lh_905[k]
                    + pb_y[k] * mh_1094[k];

        t_1461[k] = f_7 * mg0_785[k]
                    - f_8 * mg1_785[k]
                    + pb_x[k] * mh_1097[k];

        t_1462[k] = f_5 * mg0_786[k]
                    - f_6 * mg1_786[k]
                    + pb_x[k] * mh_1098[k];
    }

#pragma omp simd aligned(t_1463, t_1464, t_1465, pb_x, pb_y, pb_z, lh_885, lh_908, mg0_789, \
                         mg1_789, mh_1095, mh_1097, mh_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1463[k] = f_17 * lh_885[k]
                    + pb_z[k] * mh_1095[k];

        t_1464[k] = f_10 * lh_908[k]
                    + pb_y[k] * mh_1097[k];

        t_1465[k] = f_5 * mg0_789[k]
                    - f_6 * mg1_789[k]
                    + pb_x[k] * mh_1101[k];
    }

#pragma omp simd aligned(t_1466, t_1467, t_1468, pb_x, pb_z, lh_888, mg0_790, mg0_792, \
                         mg1_790, mg1_792, mh_1098, mh_1102, mh_1104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1466[k] = f_3 * mg0_790[k]
                    - f_4 * mg1_790[k]
                    + pb_x[k] * mh_1102[k];

        t_1467[k] = f_17 * lh_888[k]
                    + pb_z[k] * mh_1098[k];

        t_1468[k] = f_3 * mg0_792[k]
                    - f_4 * mg1_792[k]
                    + pb_x[k] * mh_1104[k];
    }

#pragma omp simd aligned(t_1469, t_1470, t_1471, t_1472, t_1473, pb_x, pb_y, lh_912, mg0_794, \
                         mg1_794, mh_1101, mh_1106, mh_1107, mh_1108, \
                         mh_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1469[k] = f_10 * lh_912[k]
                    + pb_y[k] * mh_1101[k];

        t_1470[k] = f_3 * mg0_794[k]
                    - f_4 * mg1_794[k]
                    + pb_x[k] * mh_1106[k];

        t_1471[k] = pb_x[k] * mh_1107[k];

        t_1472[k] = pb_x[k] * mh_1108[k];

        t_1473[k] = pb_x[k] * mh_1109[k];
    }

#pragma omp simd aligned(t_1474, t_1475, t_1476, t_1477, t_1478, pa_z, pb_x, pb_z, ki0_945, \
                         ki1_945, lh_897, li_1197, mh_1107, mh_1110, mh_1111, \
                         mh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1474[k] = pb_x[k] * mh_1110[k];

        t_1475[k] = pb_x[k] * mh_1111[k];

        t_1476[k] = pb_x[k] * mh_1112[k];

        t_1477[k] = f_18 * ki0_945[k]
                    - f_19 * ki1_945[k]
                    + pa_z[k] * li_1197[k];

        t_1478[k] = f_17 * lh_897[k]
                    + pb_z[k] * mh_1107[k];
    }

#pragma omp simd aligned(t_1479, t_1480, t_1481, pb_y, lh_920, lh_921, lh_922, mg0_792, \
                         mg0_793, mg0_794, mg1_792, mg1_793, mg1_794, mh_1109, mh_1110, \
                         mh_1111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1479[k] = f_10 * lh_920[k]
                    + f_7 * mg0_792[k]
                    - f_8 * mg1_792[k]
                    + pb_y[k] * mh_1109[k];

        t_1480[k] = f_10 * lh_921[k]
                    + f_5 * mg0_793[k]
                    - f_6 * mg1_793[k]
                    + pb_y[k] * mh_1110[k];

        t_1481[k] = f_10 * lh_922[k]
                    + f_3 * mg0_794[k]
                    - f_4 * mg1_794[k]
                    + pb_y[k] * mh_1111[k];
    }

#pragma omp simd aligned(t_1482, t_1483, t_1484, t_1485, t_1486, pa_y, pb_y, ki0_1007, \
                         ki1_1007, lh_923, lh_924, li_1231, li_1232, li_1234, mh_1112, \
                         mh_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1482[k] = f_10 * lh_923[k]
                    + pb_y[k] * mh_1112[k];

        t_1483[k] = f_15 * ki0_1007[k]
                    - f_16 * ki1_1007[k]
                    + pa_y[k] * li_1231[k];

        t_1484[k] = pa_y[k] * li_1232[k];

        t_1485[k] = f_9 * lh_924[k]
                    + pb_y[k] * mh_1113[k];

        t_1486[k] = pa_y[k] * li_1234[k];
    }

#pragma omp simd aligned(t_1487, t_1488, t_1489, t_1490, pa_y, pb_y, lh_925, lh_926, lh_927, \
                         li_1235, li_1237, li_1238, mh_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1487[k] = f_10 * lh_925[k]
                    + pa_y[k] * li_1235[k];

        t_1488[k] = f_9 * lh_926[k]
                    + pb_y[k] * mh_1115[k];

        t_1489[k] = pa_y[k] * li_1237[k];

        t_1490[k] = f_11 * lh_927[k]
                    + pa_y[k] * li_1238[k];
    }

#pragma omp simd aligned(t_1491, t_1492, t_1493, t_1494, pa_y, pb_y, pb_z, lh_906, lh_929, \
                         lh_930, li_1241, li_1242, mh_1116, mh_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1491[k] = f_13 * lh_906[k]
                    + pb_z[k] * mh_1116[k];

        t_1492[k] = f_9 * lh_929[k]
                    + pb_y[k] * mh_1118[k];

        t_1493[k] = pa_y[k] * li_1241[k];

        t_1494[k] = f_12 * lh_930[k]
                    + pa_y[k] * li_1242[k];
    }

#pragma omp simd aligned(t_1495, t_1496, t_1497, t_1498, pa_y, pb_y, pb_z, lh_909, lh_932, \
                         lh_933, li_1244, li_1246, mh_1119, mh_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1495[k] = f_13 * lh_909[k]
                    + pb_z[k] * mh_1119[k];

        t_1496[k] = f_10 * lh_932[k]
                    + pa_y[k] * li_1244[k];

        t_1497[k] = f_9 * lh_933[k]
                    + pb_y[k] * mh_1122[k];

        t_1498[k] = pa_y[k] * li_1246[k];
    }

#pragma omp simd aligned(t_1499, t_1500, t_1501, t_1502, t_1503, t_1504, pb_x, mh_1128, \
                         mh_1129, mh_1130, mh_1131, mh_1132, mh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1499[k] = pb_x[k] * mh_1128[k];

        t_1500[k] = pb_x[k] * mh_1129[k];

        t_1501[k] = pb_x[k] * mh_1130[k];

        t_1502[k] = pb_x[k] * mh_1131[k];

        t_1503[k] = pb_x[k] * mh_1132[k];

        t_1504[k] = pb_x[k] * mh_1133[k];
    }

#pragma omp simd aligned(t_1505, t_1506, t_1507, t_1508, pa_y, pb_z, lh_918, lh_939, lh_941, \
                         lh_942, li_1253, li_1255, li_1256, mh_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1505[k] = f_14 * lh_939[k]
                    + pa_y[k] * li_1253[k];

        t_1506[k] = f_13 * lh_918[k]
                    + pb_z[k] * mh_1128[k];

        t_1507[k] = f_12 * lh_941[k]
                    + pa_y[k] * li_1255[k];

        t_1508[k] = f_11 * lh_942[k]
                    + pa_y[k] * li_1256[k];
    }

#pragma omp simd aligned(t_1509, t_1510, t_1511, t_1512, t_1513, pa_y, pb_x, pb_y, lh_943, \
                         lh_944, li_1257, li_1259, mg0_810, mg1_810, mh_1133, \
                         mh_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1509[k] = f_10 * lh_943[k]
                    + pa_y[k] * li_1257[k];

        t_1510[k] = f_9 * lh_944[k]
                    + pb_y[k] * mh_1133[k];

        t_1511[k] = pa_y[k] * li_1259[k];

        t_1512[k] = f_1 * mg0_810[k]
                    - f_2 * mg1_810[k]
                    + pb_x[k] * mh_1134[k];

        t_1513[k] = pb_y[k] * mh_1134[k];
    }

#pragma omp simd aligned(t_1514, t_1515, t_1516, t_1517, pb_x, pb_y, pb_z, lh_924, mg0_813, \
                         mg0_815, mg1_813, mg1_815, mh_1134, mh_1136, mh_1137, \
                         mh_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1514[k] = f_0 * lh_924[k]
                    + pb_z[k] * mh_1134[k];

        t_1515[k] = f_7 * mg0_813[k]
                    - f_8 * mg1_813[k]
                    + pb_x[k] * mh_1137[k];

        t_1516[k] = pb_y[k] * mh_1136[k];

        t_1517[k] = f_7 * mg0_815[k]
                    - f_8 * mg1_815[k]
                    + pb_x[k] * mh_1139[k];
    }

#pragma omp simd aligned(t_1518, t_1519, t_1520, t_1521, pb_x, pb_y, pb_z, lh_927, mg0_816, \
                         mg0_819, mg1_816, mg1_819, mh_1137, mh_1139, mh_1140, \
                         mh_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1518[k] = f_5 * mg0_816[k]
                    - f_6 * mg1_816[k]
                    + pb_x[k] * mh_1140[k];

        t_1519[k] = f_0 * lh_927[k]
                    + pb_z[k] * mh_1137[k];

        t_1520[k] = pb_y[k] * mh_1139[k];

        t_1521[k] = f_5 * mg0_819[k]
                    - f_6 * mg1_819[k]
                    + pb_x[k] * mh_1143[k];
    }

#pragma omp simd aligned(t_1522, t_1523, t_1524, t_1525, pb_x, pb_y, pb_z, lh_930, mg0_820, \
                         mg0_822, mg1_820, mg1_822, mh_1140, mh_1143, mh_1144, \
                         mh_1146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1522[k] = f_3 * mg0_820[k]
                    - f_4 * mg1_820[k]
                    + pb_x[k] * mh_1144[k];

        t_1523[k] = f_0 * lh_930[k]
                    + pb_z[k] * mh_1140[k];

        t_1524[k] = f_3 * mg0_822[k]
                    - f_4 * mg1_822[k]
                    + pb_x[k] * mh_1146[k];

        t_1525[k] = pb_y[k] * mh_1143[k];
    }

#pragma omp simd aligned(t_1526, t_1527, t_1528, t_1529, t_1530, t_1531, pb_x, mg0_824, \
                         mg1_824, mh_1148, mh_1149, mh_1150, mh_1151, mh_1152, \
                         mh_1153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1526[k] = f_3 * mg0_824[k]
                    - f_4 * mg1_824[k]
                    + pb_x[k] * mh_1148[k];

        t_1527[k] = pb_x[k] * mh_1149[k];

        t_1528[k] = pb_x[k] * mh_1150[k];

        t_1529[k] = pb_x[k] * mh_1151[k];

        t_1530[k] = pb_x[k] * mh_1152[k];

        t_1531[k] = pb_x[k] * mh_1153[k];
    }

#pragma omp simd aligned(t_1532, t_1533, t_1534, t_1535, pb_x, pb_y, pb_z, lh_939, mg0_820, \
                         mg0_822, mg1_820, mg1_822, mh_1149, mh_1151, \
                         mh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1532[k] = pb_x[k] * mh_1154[k];

        t_1533[k] = f_1 * mg0_820[k]
                    - f_2 * mg1_820[k]
                    + pb_y[k] * mh_1149[k];

        t_1534[k] = f_0 * lh_939[k]
                    + pb_z[k] * mh_1149[k];

        t_1535[k] = f_7 * mg0_822[k]
                    - f_8 * mg1_822[k]
                    + pb_y[k] * mh_1151[k];
    }

#pragma omp simd aligned(t_1536, t_1537, t_1538, t_1539, pb_y, pb_z, lh_944, mg0_823, mg0_824, \
                         mg1_823, mg1_824, mh_1152, mh_1153, mh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1536[k] = f_5 * mg0_823[k]
                    - f_6 * mg1_823[k]
                    + pb_y[k] * mh_1152[k];

        t_1537[k] = f_3 * mg0_824[k]
                    - f_4 * mg1_824[k]
                    + pb_y[k] * mh_1153[k];

        t_1538[k] = pb_y[k] * mh_1154[k];

        t_1539[k] = f_0 * lh_944[k]
                    + f_1 * mg0_824[k]
                    - f_2 * mg1_824[k]
                    + pb_z[k] * mh_1154[k];
    }
}

auto
compute_prim_mi_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ki0, const size_t ki1,
                                     const size_t lh, const size_t li, const size_t mg0,
                                     const size_t mg1, const size_t mh, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    compute_prim_mi_electron_repulsion_0_piece0(buffer, target, pa, pb, ki0, ki1, lh, li, mg0,
                                                mg1, mh, ncols, alpha, beta, p);

    compute_prim_mi_electron_repulsion_0_piece1(buffer, target, pa, pb, ki0, ki1, lh, li, mg0,
                                                mg1, mh, ncols, alpha, beta, p);

    compute_prim_mi_electron_repulsion_0_piece2(buffer, target, pa, pb, ki0, ki1, lh, li, mg0,
                                                mg1, mh, ncols, alpha, beta, p);

    compute_prim_mi_electron_repulsion_0_piece3(buffer, target, pa, pb, ki0, ki1, lh, li, mg0,
                                                mg1, mh, ncols, alpha, beta, p);

    compute_prim_mi_electron_repulsion_0_piece4(buffer, target, pa, pb, ki0, ki1, lh, li, mg0,
                                                mg1, mh, ncols, alpha, beta, p);

    compute_prim_mi_electron_repulsion_0_piece5(buffer, target, pa, pb, ki0, ki1, lh, li, mg0,
                                                mg1, mh, ncols, alpha, beta, p);

    compute_prim_mi_electron_repulsion_0_piece6(buffer, target, pa, pb, ki0, ki1, lh, li, mg0,
                                                mg1, mh, ncols, alpha, beta, p);

    compute_prim_mi_electron_repulsion_0_piece7(buffer, target, pa, pb, ki0, ki1, lh, li, mg0,
                                                mg1, mh, ncols, alpha, beta, p);

    compute_prim_mi_electron_repulsion_0_piece8(buffer, target, pa, pb, lh, li, mh, ncols, p);

    compute_prim_mi_electron_repulsion_0_piece9(buffer, target, pa, pb, lh, li, mg0, mg1, mh,
                                                ncols, alpha, beta, p);

    compute_prim_mi_electron_repulsion_0_piece10(buffer, target, pa, pb, ki0, ki1, lh, li, mg0,
                                                 mg1, mh, ncols, alpha, beta, p);

    compute_prim_mi_electron_repulsion_0_piece11(buffer, target, pa, pb, ki0, ki1, lh, li, mg0,
                                                 mg1, mh, ncols, alpha, beta, p);
}

}  // namespace simdt2ceri
