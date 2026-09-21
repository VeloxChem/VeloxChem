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


#include "SimdElectronRepulsionVrrRecMG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_mg_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kg0,
                                            const size_t kg1, const size_t lf, const size_t lg,
                                            const size_t md0, const size_t md1, const size_t mf,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 4.0 / p;
    const auto f_8 = 2.0 / p;
    const auto f_9 = 0.5 / alpha;
    const auto f_10 = 0.5 * beta / (alpha * p);
    const auto f_11 = 3.5 / p;
    const auto f_12 = 3.0 / alpha;
    const auto f_13 = 3.0 * beta / (alpha * p);
    const auto f_14 = 1.0 / alpha;
    const auto f_15 = beta / (alpha * p);
    const auto f_16 = 1.5 / p;
    const auto f_17 = 3.0 / p;
    const auto f_18 = 2.5 / alpha;
    const auto f_19 = 2.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kg0_0 = buffer.data(kg0 + 0);
    const auto *kg0_15 = buffer.data(kg0 + 15);
    const auto *kg0_30 = buffer.data(kg0 + 30);
    const auto *kg0_55 = buffer.data(kg0 + 55);
    const auto *kg0_89 = buffer.data(kg0 + 89);
    const auto *kg0_100 = buffer.data(kg0 + 100);

    const auto *kg1_0 = buffer.data(kg1 + 0);
    const auto *kg1_15 = buffer.data(kg1 + 15);
    const auto *kg1_30 = buffer.data(kg1 + 30);
    const auto *kg1_55 = buffer.data(kg1 + 55);
    const auto *kg1_89 = buffer.data(kg1 + 89);
    const auto *kg1_100 = buffer.data(kg1 + 100);

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_7 = buffer.data(lf + 7);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);

    const auto *lg_0 = buffer.data(lg + 0);
    const auto *lg_3 = buffer.data(lg + 3);
    const auto *lg_5 = buffer.data(lg + 5);
    const auto *lg_6 = buffer.data(lg + 6);
    const auto *lg_9 = buffer.data(lg + 9);
    const auto *lg_10 = buffer.data(lg + 10);
    const auto *lg_12 = buffer.data(lg + 12);
    const auto *lg_14 = buffer.data(lg + 14);
    const auto *lg_15 = buffer.data(lg + 15);
    const auto *lg_16 = buffer.data(lg + 16);
    const auto *lg_18 = buffer.data(lg + 18);
    const auto *lg_21 = buffer.data(lg + 21);
    const auto *lg_25 = buffer.data(lg + 25);
    const auto *lg_30 = buffer.data(lg + 30);
    const auto *lg_32 = buffer.data(lg + 32);
    const auto *lg_35 = buffer.data(lg + 35);
    const auto *lg_39 = buffer.data(lg + 39);
    const auto *lg_42 = buffer.data(lg + 42);
    const auto *lg_44 = buffer.data(lg + 44);
    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_46 = buffer.data(lg + 46);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_51 = buffer.data(lg + 51);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_75 = buffer.data(lg + 75);
    const auto *lg_77 = buffer.data(lg + 77);
    const auto *lg_78 = buffer.data(lg + 78);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_84 = buffer.data(lg + 84);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_100 = buffer.data(lg + 100);

    const auto *md0_0 = buffer.data(md0 + 0);
    const auto *md0_3 = buffer.data(md0 + 3);
    const auto *md0_5 = buffer.data(md0 + 5);
    const auto *md0_18 = buffer.data(md0 + 18);
    const auto *md0_21 = buffer.data(md0 + 21);
    const auto *md0_23 = buffer.data(md0 + 23);
    const auto *md0_30 = buffer.data(md0 + 30);
    const auto *md0_33 = buffer.data(md0 + 33);
    const auto *md0_35 = buffer.data(md0 + 35);
    const auto *md0_36 = buffer.data(md0 + 36);
    const auto *md0_39 = buffer.data(md0 + 39);
    const auto *md0_41 = buffer.data(md0 + 41);
    const auto *md0_54 = buffer.data(md0 + 54);
    const auto *md0_59 = buffer.data(md0 + 59);

    const auto *md1_0 = buffer.data(md1 + 0);
    const auto *md1_3 = buffer.data(md1 + 3);
    const auto *md1_5 = buffer.data(md1 + 5);
    const auto *md1_18 = buffer.data(md1 + 18);
    const auto *md1_21 = buffer.data(md1 + 21);
    const auto *md1_23 = buffer.data(md1 + 23);
    const auto *md1_30 = buffer.data(md1 + 30);
    const auto *md1_33 = buffer.data(md1 + 33);
    const auto *md1_35 = buffer.data(md1 + 35);
    const auto *md1_36 = buffer.data(md1 + 36);
    const auto *md1_39 = buffer.data(md1 + 39);
    const auto *md1_41 = buffer.data(md1 + 41);
    const auto *md1_54 = buffer.data(md1 + 54);
    const auto *md1_59 = buffer.data(md1 + 59);

    const auto *mf_0 = buffer.data(mf + 0);
    const auto *mf_1 = buffer.data(mf + 1);
    const auto *mf_2 = buffer.data(mf + 2);
    const auto *mf_3 = buffer.data(mf + 3);
    const auto *mf_5 = buffer.data(mf + 5);
    const auto *mf_6 = buffer.data(mf + 6);
    const auto *mf_8 = buffer.data(mf + 8);
    const auto *mf_9 = buffer.data(mf + 9);
    const auto *mf_10 = buffer.data(mf + 10);
    const auto *mf_11 = buffer.data(mf + 11);
    const auto *mf_13 = buffer.data(mf + 13);
    const auto *mf_16 = buffer.data(mf + 16);
    const auto *mf_18 = buffer.data(mf + 18);
    const auto *mf_19 = buffer.data(mf + 19);
    const auto *mf_20 = buffer.data(mf + 20);
    const auto *mf_22 = buffer.data(mf + 22);
    const auto *mf_25 = buffer.data(mf + 25);
    const auto *mf_26 = buffer.data(mf + 26);
    const auto *mf_27 = buffer.data(mf + 27);
    const auto *mf_29 = buffer.data(mf + 29);
    const auto *mf_30 = buffer.data(mf + 30);
    const auto *mf_31 = buffer.data(mf + 31);
    const auto *mf_32 = buffer.data(mf + 32);
    const auto *mf_33 = buffer.data(mf + 33);
    const auto *mf_36 = buffer.data(mf + 36);
    const auto *mf_37 = buffer.data(mf + 37);
    const auto *mf_38 = buffer.data(mf + 38);
    const auto *mf_39 = buffer.data(mf + 39);
    const auto *mf_42 = buffer.data(mf + 42);
    const auto *mf_46 = buffer.data(mf + 46);
    const auto *mf_47 = buffer.data(mf + 47);
    const auto *mf_48 = buffer.data(mf + 48);
    const auto *mf_49 = buffer.data(mf + 49);
    const auto *mf_50 = buffer.data(mf + 50);
    const auto *mf_51 = buffer.data(mf + 51);
    const auto *mf_52 = buffer.data(mf + 52);
    const auto *mf_55 = buffer.data(mf + 55);
    const auto *mf_56 = buffer.data(mf + 56);
    const auto *mf_57 = buffer.data(mf + 57);
    const auto *mf_58 = buffer.data(mf + 58);
    const auto *mf_59 = buffer.data(mf + 59);
    const auto *mf_60 = buffer.data(mf + 60);
    const auto *mf_61 = buffer.data(mf + 61);
    const auto *mf_62 = buffer.data(mf + 62);
    const auto *mf_63 = buffer.data(mf + 63);
    const auto *mf_66 = buffer.data(mf + 66);
    const auto *mf_67 = buffer.data(mf + 67);
    const auto *mf_68 = buffer.data(mf + 68);
    const auto *mf_69 = buffer.data(mf + 69);
    const auto *mf_70 = buffer.data(mf + 70);
    const auto *mf_72 = buffer.data(mf + 72);
    const auto *mf_76 = buffer.data(mf + 76);
    const auto *mf_77 = buffer.data(mf + 77);
    const auto *mf_78 = buffer.data(mf + 78);
    const auto *mf_79 = buffer.data(mf + 79);
    const auto *mf_80 = buffer.data(mf + 80);
    const auto *mf_82 = buffer.data(mf + 82);
    const auto *mf_86 = buffer.data(mf + 86);
    const auto *mf_87 = buffer.data(mf + 87);
    const auto *mf_88 = buffer.data(mf + 88);
    const auto *mf_89 = buffer.data(mf + 89);
    const auto *mf_90 = buffer.data(mf + 90);
    const auto *mf_91 = buffer.data(mf + 91);
    const auto *mf_92 = buffer.data(mf + 92);
    const auto *mf_95 = buffer.data(mf + 95);
    const auto *mf_96 = buffer.data(mf + 96);
    const auto *mf_97 = buffer.data(mf + 97);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, lf_0, md0_0, md1_0, \
                         mf_0, mf_1, mf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lf_0[k]
                 + f_1 * md0_0[k]
                 - f_2 * md1_0[k]
                 + pb_x[k] * mf_0[k];

        t_1[k] = pb_y[k] * mf_0[k];

        t_2[k] = pb_z[k] * mf_0[k];

        t_3[k] = f_3 * md0_0[k]
                 - f_4 * md1_0[k]
                 + pb_y[k] * mf_1[k];

        t_4[k] = pb_y[k] * mf_2[k];

        t_5[k] = f_3 * md0_0[k]
                 - f_4 * md1_0[k]
                 + pb_z[k] * mf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, lf_6, lf_9, md0_3, md1_3, \
                         mf_3, mf_5, mf_6, mf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * lf_6[k]
                 + pb_x[k] * mf_6[k];

        t_7[k] = pb_z[k] * mf_3[k];

        t_8[k] = pb_y[k] * mf_5[k];

        t_9[k] = f_0 * lf_9[k]
                 + pb_x[k] * mf_9[k];

        t_10[k] = f_1 * md0_3[k]
                  - f_2 * md1_3[k]
                  + pb_y[k] * mf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, lg_0, md0_5, md1_5, \
                         mf_6, mf_8, mf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * mf_6[k];

        t_12[k] = f_3 * md0_5[k]
                  - f_4 * md1_5[k]
                  + pb_y[k] * mf_8[k];

        t_13[k] = pb_y[k] * mf_9[k];

        t_14[k] = f_1 * md0_5[k]
                  - f_2 * md1_5[k]
                  + pb_z[k] * mf_9[k];

        t_15[k] = pa_y[k] * lg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pb_y, pb_z, lf_0, lf_1, lg_3, \
                         lg_5, mf_10, mf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * lf_0[k]
                  + pb_y[k] * mf_10[k];

        t_17[k] = pb_z[k] * mf_10[k];

        t_18[k] = f_6 * lf_1[k]
                  + pa_y[k] * lg_3[k];

        t_19[k] = pb_z[k] * mf_11[k];

        t_20[k] = pa_y[k] * lg_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_z, lf_6, lf_16, lf_18, \
                         lg_9, lg_10, mf_13, mf_16, mf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * lf_16[k]
                  + pb_x[k] * mf_16[k];

        t_22[k] = pb_z[k] * mf_13[k];

        t_23[k] = f_7 * lf_18[k]
                  + pb_x[k] * mf_18[k];

        t_24[k] = pa_y[k] * lg_9[k];

        t_25[k] = f_8 * lf_6[k]
                  + pa_y[k] * lg_10[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, pb_z, lf_8, lf_9, \
                         lg_0, lg_12, lg_14, mf_16, mf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * mf_16[k];

        t_27[k] = f_6 * lf_8[k]
                  + pa_y[k] * lg_12[k];

        t_28[k] = f_5 * lf_9[k]
                  + pb_y[k] * mf_19[k];

        t_29[k] = pa_y[k] * lg_14[k];

        t_30[k] = pa_z[k] * lg_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, lf_0, lf_2, \
                         lg_3, lg_5, lg_6, mf_20, mf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * mf_20[k];

        t_32[k] = f_5 * lf_0[k]
                  + pb_z[k] * mf_20[k];

        t_33[k] = pa_z[k] * lg_3[k];

        t_34[k] = pb_y[k] * mf_22[k];

        t_35[k] = f_6 * lf_2[k]
                  + pa_z[k] * lg_5[k];

        t_36[k] = pa_z[k] * lg_6[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, pb_y, lf_27, lf_29, lg_10, mf_25, \
                         mf_27, mf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_7 * lf_27[k]
                  + pb_x[k] * mf_27[k];

        t_38[k] = pb_y[k] * mf_25[k];

        t_39[k] = f_7 * lf_29[k]
                  + pb_x[k] * mf_29[k];

        t_40[k] = pa_z[k] * lg_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_z, pb_y, pb_z, lf_6, lf_7, lf_9, lg_12, \
                         lg_14, mf_26, mf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * lf_6[k]
                  + pb_z[k] * mf_26[k];

        t_42[k] = f_6 * lf_7[k]
                  + pa_z[k] * lg_12[k];

        t_43[k] = pb_y[k] * mf_29[k];

        t_44[k] = f_8 * lf_9[k]
                  + pa_z[k] * lg_14[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_y, pb_z, kg0_0, kg1_0, lf_10, lg_15, \
                         mf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_9 * kg0_0[k]
                  - f_10 * kg1_0[k]
                  + pa_y[k] * lg_15[k];

        t_46[k] = f_6 * lf_10[k]
                  + pb_y[k] * mf_30[k];

        t_47[k] = pb_z[k] * mf_30[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_x, pb_z, lf_33, lf_36, md0_18, md0_21, \
                         md1_18, md1_21, mf_31, mf_32, mf_33, mf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * lf_33[k]
                  + f_3 * md0_21[k]
                  - f_4 * md1_21[k]
                  + pb_x[k] * mf_33[k];

        t_49[k] = pb_z[k] * mf_31[k];

        t_50[k] = f_3 * md0_18[k]
                  - f_4 * md1_18[k]
                  + pb_z[k] * mf_32[k];

        t_51[k] = f_11 * lf_36[k]
                  + pb_x[k] * mf_36[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_z, kg0_55, kg1_55, lf_38, \
                         lf_39, lg_55, mf_33, mf_38, mf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_z[k] * mf_33[k];

        t_53[k] = f_11 * lf_38[k]
                  + pb_x[k] * mf_38[k];

        t_54[k] = f_11 * lf_39[k]
                  + pb_x[k] * mf_39[k];

        t_55[k] = f_12 * kg0_55[k]
                  - f_13 * kg1_55[k]
                  + pa_x[k] * lg_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pb_z, lf_19, md0_21, md0_23, md1_21, \
                         md1_23, mf_36, mf_37, mf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * mf_36[k];

        t_57[k] = f_3 * md0_21[k]
                  - f_4 * md1_21[k]
                  + pb_z[k] * mf_37[k];

        t_58[k] = f_6 * lf_19[k]
                  + pb_y[k] * mf_39[k];

        t_59[k] = f_1 * md0_23[k]
                  - f_2 * md1_23[k]
                  + pb_z[k] * mf_39[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, lf_22, lg_16, \
                         lg_18, lg_30, lg_32, lg_35, mf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * lg_30[k];

        t_61[k] = pa_z[k] * lg_16[k];

        t_62[k] = pa_y[k] * lg_32[k];

        t_63[k] = pa_z[k] * lg_18[k];

        t_64[k] = f_5 * lf_22[k]
                  + pb_y[k] * mf_42[k];

        t_65[k] = pa_y[k] * lg_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, lf_47, lf_48, lg_21, \
                         lg_25, lg_39, mf_47, mf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * lg_21[k];

        t_67[k] = f_11 * lf_47[k]
                  + pb_x[k] * mf_47[k];

        t_68[k] = f_11 * lf_48[k]
                  + pb_x[k] * mf_48[k];

        t_69[k] = pa_y[k] * lg_39[k];

        t_70[k] = pa_z[k] * lg_25[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, pb_z, lf_16, lf_28, lf_29, lg_42, \
                         lg_44, mf_46, mf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * lf_16[k]
                  + pb_z[k] * mf_46[k];

        t_72[k] = f_6 * lf_28[k]
                  + pa_y[k] * lg_42[k];

        t_73[k] = f_5 * lf_29[k]
                  + pb_y[k] * mf_49[k];

        t_74[k] = pa_y[k] * lg_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, kg0_0, kg1_0, lf_20, lg_30, \
                         md0_30, md1_30, mf_50, mf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_9 * kg0_0[k]
                  - f_10 * kg1_0[k]
                  + pa_z[k] * lg_30[k];

        t_76[k] = pb_y[k] * mf_50[k];

        t_77[k] = f_6 * lf_20[k]
                  + pb_z[k] * mf_50[k];

        t_78[k] = f_3 * md0_30[k]
                  - f_4 * md1_30[k]
                  + pb_y[k] * mf_51[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_x, pb_y, lf_55, lf_56, lf_57, \
                         md0_35, md1_35, mf_52, mf_55, mf_56, mf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_y[k] * mf_52[k];

        t_80[k] = f_11 * lf_55[k]
                  + f_3 * md0_35[k]
                  - f_4 * md1_35[k]
                  + pb_x[k] * mf_55[k];

        t_81[k] = f_11 * lf_56[k]
                  + pb_x[k] * mf_56[k];

        t_82[k] = f_11 * lf_57[k]
                  + pb_x[k] * mf_57[k];

        t_83[k] = pb_y[k] * mf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, pb_y, pb_z, lf_26, lf_59, md0_33, \
                         md0_35, md1_33, md1_35, mf_56, mf_58, mf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_11 * lf_59[k]
                  + pb_x[k] * mf_59[k];

        t_85[k] = f_1 * md0_33[k]
                  - f_2 * md1_33[k]
                  + pb_y[k] * mf_56[k];

        t_86[k] = f_6 * lf_26[k]
                  + pb_z[k] * mf_56[k];

        t_87[k] = f_3 * md0_35[k]
                  - f_4 * md1_35[k]
                  + pb_y[k] * mf_58[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pa_y, pb_y, kg0_15, kg0_89, kg1_15, \
                         kg1_89, lf_30, lg_45, lg_89, mf_59, mf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pb_y[k] * mf_59[k];

        t_89[k] = f_12 * kg0_89[k]
                  - f_13 * kg1_89[k]
                  + pa_x[k] * lg_89[k];

        t_90[k] = f_14 * kg0_15[k]
                  - f_15 * kg1_15[k]
                  + pa_y[k] * lg_45[k];

        t_91[k] = f_16 * lf_30[k]
                  + pb_y[k] * mf_60[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_z, lf_63, md0_36, md0_39, md1_36, \
                         md1_39, mf_60, mf_61, mf_62, mf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_z[k] * mf_60[k];

        t_93[k] = f_17 * lf_63[k]
                  + f_3 * md0_39[k]
                  - f_4 * md1_39[k]
                  + pb_x[k] * mf_63[k];

        t_94[k] = pb_z[k] * mf_61[k];

        t_95[k] = f_3 * md0_36[k]
                  - f_4 * md1_36[k]
                  + pb_z[k] * mf_62[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_x, pb_z, lf_66, lf_68, lf_69, mf_63, \
                         mf_66, mf_68, mf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_17 * lf_66[k]
                  + pb_x[k] * mf_66[k];

        t_97[k] = pb_z[k] * mf_63[k];

        t_98[k] = f_17 * lf_68[k]
                  + pb_x[k] * mf_68[k];

        t_99[k] = f_17 * lf_69[k]
                  + pb_x[k] * mf_69[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_y, pb_z, kg0_100, kg1_100, \
                         lf_39, lg_100, md0_39, md1_39, mf_66, mf_67, \
                         mf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_18 * kg0_100[k]
                   - f_19 * kg1_100[k]
                   + pa_x[k] * lg_100[k];

        t_101[k] = pb_z[k] * mf_66[k];

        t_102[k] = f_3 * md0_39[k]
                   - f_4 * md1_39[k]
                   + pb_z[k] * mf_67[k];

        t_103[k] = f_16 * lf_39[k]
                   + pb_y[k] * mf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_z, pb_z, lf_30, lg_45, lg_46, \
                         lg_48, md0_41, md1_41, mf_69, mf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * md0_41[k]
                   - f_2 * md1_41[k]
                   + pb_z[k] * mf_69[k];

        t_105[k] = pa_z[k] * lg_45[k];

        t_106[k] = pa_z[k] * lg_46[k];

        t_107[k] = f_5 * lf_30[k]
                   + pb_z[k] * mf_70[k];

        t_108[k] = pa_z[k] * lg_48[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_z, pb_x, pb_y, lf_32, lf_42, lf_77, \
                         lg_50, lg_51, mf_72, mf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_6 * lf_42[k]
                   + pb_y[k] * mf_72[k];

        t_110[k] = f_6 * lf_32[k]
                   + pa_z[k] * lg_50[k];

        t_111[k] = pa_z[k] * lg_51[k];

        t_112[k] = f_17 * lf_77[k]
                   + pb_x[k] * mf_77[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pb_x, pb_z, lf_36, lf_78, lf_79, \
                         lg_55, mf_76, mf_78, mf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_17 * lf_78[k]
                   + pb_x[k] * mf_78[k];

        t_114[k] = f_17 * lf_79[k]
                   + pb_x[k] * mf_79[k];

        t_115[k] = pa_z[k] * lg_55[k];

        t_116[k] = f_5 * lf_36[k]
                   + pb_z[k] * mf_76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pb_y, lf_37, lf_39, lf_49, \
                         lg_57, lg_59, lg_75, mf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_6 * lf_37[k]
                   + pa_z[k] * lg_57[k];

        t_118[k] = f_6 * lf_49[k]
                   + pb_y[k] * mf_79[k];

        t_119[k] = f_8 * lf_39[k]
                   + pa_z[k] * lg_59[k];

        t_120[k] = pa_y[k] * lg_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pa_y, pb_y, lf_50, lf_51, lf_52, \
                         lg_77, lg_78, lg_80, mf_80, mf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * lf_50[k]
                   + pb_y[k] * mf_80[k];

        t_122[k] = pa_y[k] * lg_77[k];

        t_123[k] = f_6 * lf_51[k]
                   + pa_y[k] * lg_78[k];

        t_124[k] = f_5 * lf_52[k]
                   + pb_y[k] * mf_82[k];

        t_125[k] = pa_y[k] * lg_80[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, pa_y, pb_x, lf_56, lf_86, lf_87, \
                         lf_88, lg_84, lg_85, mf_86, mf_87, mf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_17 * lf_86[k]
                   + pb_x[k] * mf_86[k];

        t_127[k] = f_17 * lf_87[k]
                   + pb_x[k] * mf_87[k];

        t_128[k] = f_17 * lf_88[k]
                   + pb_x[k] * mf_88[k];

        t_129[k] = pa_y[k] * lg_84[k];

        t_130[k] = f_8 * lf_56[k]
                   + pa_y[k] * lg_85[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_y, pb_y, pb_z, lf_46, lf_58, lf_59, \
                         lg_87, lg_89, mf_86, mf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_6 * lf_46[k]
                   + pb_z[k] * mf_86[k];

        t_132[k] = f_6 * lf_58[k]
                   + pa_y[k] * lg_87[k];

        t_133[k] = f_5 * lf_59[k]
                   + pb_y[k] * mf_89[k];

        t_134[k] = pa_y[k] * lg_89[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pa_z, pb_y, pb_z, kg0_30, kg1_30, lf_50, \
                         lg_75, md0_54, md1_54, mf_90, mf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_14 * kg0_30[k]
                   - f_15 * kg1_30[k]
                   + pa_z[k] * lg_75[k];

        t_136[k] = pb_y[k] * mf_90[k];

        t_137[k] = f_16 * lf_50[k]
                   + pb_z[k] * mf_90[k];

        t_138[k] = f_3 * md0_54[k]
                   - f_4 * md1_54[k]
                   + pb_y[k] * mf_91[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pb_x, pb_y, lf_95, lf_96, lf_97, \
                         md0_59, md1_59, mf_92, mf_95, mf_96, mf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * mf_92[k];

        t_140[k] = f_17 * lf_95[k]
                   + f_3 * md0_59[k]
                   - f_4 * md1_59[k]
                   + pb_x[k] * mf_95[k];

        t_141[k] = f_17 * lf_96[k]
                   + pb_x[k] * mf_96[k];

        t_142[k] = f_17 * lf_97[k]
                   + pb_x[k] * mf_97[k];

        t_143[k] = pb_y[k] * mf_95[k];
    }
}

static auto
compute_prim_mg_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kg0,
                                            const size_t kg1, const size_t lf, const size_t lg,
                                            const size_t md0, const size_t md1, const size_t mf,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_8 = 2.0 / p;
    const auto f_9 = 0.5 / alpha;
    const auto f_10 = 0.5 * beta / (alpha * p);
    const auto f_14 = 1.0 / alpha;
    const auto f_15 = beta / (alpha * p);
    const auto f_16 = 1.5 / p;
    const auto f_17 = 3.0 / p;
    const auto f_18 = 2.5 / alpha;
    const auto f_19 = 2.5 * beta / (alpha * p);
    const auto f_20 = 1.5 / alpha;
    const auto f_21 = 1.5 * beta / (alpha * p);
    const auto f_22 = 2.5 / p;
    const auto f_23 = 2.0 / alpha;
    const auto f_24 = 2.0 * beta / (alpha * p);

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

    const auto *kg0_45 = buffer.data(kg0 + 45);
    const auto *kg0_48 = buffer.data(kg0 + 48);
    const auto *kg0_75 = buffer.data(kg0 + 75);
    const auto *kg0_80 = buffer.data(kg0 + 80);
    const auto *kg0_90 = buffer.data(kg0 + 90);
    const auto *kg0_93 = buffer.data(kg0 + 93);
    const auto *kg0_108 = buffer.data(kg0 + 108);
    const auto *kg0_120 = buffer.data(kg0 + 120);
    const auto *kg0_125 = buffer.data(kg0 + 125);
    const auto *kg0_135 = buffer.data(kg0 + 135);
    const auto *kg0_140 = buffer.data(kg0 + 140);
    const auto *kg0_149 = buffer.data(kg0 + 149);
    const auto *kg0_160 = buffer.data(kg0 + 160);
    const auto *kg0_190 = buffer.data(kg0 + 190);
    const auto *kg0_192 = buffer.data(kg0 + 192);
    const auto *kg0_194 = buffer.data(kg0 + 194);
    const auto *kg0_224 = buffer.data(kg0 + 224);
    const auto *kg0_235 = buffer.data(kg0 + 235);
    const auto *kg0_265 = buffer.data(kg0 + 265);
    const auto *kg0_267 = buffer.data(kg0 + 267);
    const auto *kg0_269 = buffer.data(kg0 + 269);

    const auto *kg1_45 = buffer.data(kg1 + 45);
    const auto *kg1_48 = buffer.data(kg1 + 48);
    const auto *kg1_75 = buffer.data(kg1 + 75);
    const auto *kg1_80 = buffer.data(kg1 + 80);
    const auto *kg1_90 = buffer.data(kg1 + 90);
    const auto *kg1_93 = buffer.data(kg1 + 93);
    const auto *kg1_108 = buffer.data(kg1 + 108);
    const auto *kg1_120 = buffer.data(kg1 + 120);
    const auto *kg1_125 = buffer.data(kg1 + 125);
    const auto *kg1_135 = buffer.data(kg1 + 135);
    const auto *kg1_140 = buffer.data(kg1 + 140);
    const auto *kg1_149 = buffer.data(kg1 + 149);
    const auto *kg1_160 = buffer.data(kg1 + 160);
    const auto *kg1_190 = buffer.data(kg1 + 190);
    const auto *kg1_192 = buffer.data(kg1 + 192);
    const auto *kg1_194 = buffer.data(kg1 + 194);
    const auto *kg1_224 = buffer.data(kg1 + 224);
    const auto *kg1_235 = buffer.data(kg1 + 235);
    const auto *kg1_265 = buffer.data(kg1 + 265);
    const auto *kg1_267 = buffer.data(kg1 + 267);
    const auto *kg1_269 = buffer.data(kg1 + 269);

    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);
    const auto *lf_108 = buffer.data(lf + 108);
    const auto *lf_109 = buffer.data(lf + 109);
    const auto *lf_110 = buffer.data(lf + 110);
    const auto *lf_112 = buffer.data(lf + 112);
    const auto *lf_116 = buffer.data(lf + 116);
    const auto *lf_117 = buffer.data(lf + 117);
    const auto *lf_118 = buffer.data(lf + 118);
    const auto *lf_119 = buffer.data(lf + 119);
    const auto *lf_120 = buffer.data(lf + 120);
    const auto *lf_122 = buffer.data(lf + 122);
    const auto *lf_126 = buffer.data(lf + 126);
    const auto *lf_127 = buffer.data(lf + 127);
    const auto *lf_128 = buffer.data(lf + 128);
    const auto *lf_129 = buffer.data(lf + 129);
    const auto *lf_130 = buffer.data(lf + 130);
    const auto *lf_132 = buffer.data(lf + 132);
    const auto *lf_136 = buffer.data(lf + 136);
    const auto *lf_137 = buffer.data(lf + 137);
    const auto *lf_138 = buffer.data(lf + 138);
    const auto *lf_145 = buffer.data(lf + 145);
    const auto *lf_146 = buffer.data(lf + 146);
    const auto *lf_147 = buffer.data(lf + 147);
    const auto *lf_149 = buffer.data(lf + 149);
    const auto *lf_153 = buffer.data(lf + 153);
    const auto *lf_156 = buffer.data(lf + 156);
    const auto *lf_158 = buffer.data(lf + 158);
    const auto *lf_159 = buffer.data(lf + 159);
    const auto *lf_167 = buffer.data(lf + 167);
    const auto *lf_168 = buffer.data(lf + 168);
    const auto *lf_169 = buffer.data(lf + 169);
    const auto *lf_176 = buffer.data(lf + 176);
    const auto *lf_177 = buffer.data(lf + 177);
    const auto *lf_178 = buffer.data(lf + 178);
    const auto *lf_179 = buffer.data(lf + 179);

    const auto *lg_90 = buffer.data(lg + 90);
    const auto *lg_91 = buffer.data(lg + 91);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_95 = buffer.data(lg + 95);
    const auto *lg_96 = buffer.data(lg + 96);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_104 = buffer.data(lg + 104);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_137 = buffer.data(lg + 137);
    const auto *lg_138 = buffer.data(lg + 138);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_144 = buffer.data(lg + 144);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_149 = buffer.data(lg + 149);
    const auto *lg_150 = buffer.data(lg + 150);
    const auto *lg_151 = buffer.data(lg + 151);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_155 = buffer.data(lg + 155);
    const auto *lg_156 = buffer.data(lg + 156);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_168 = buffer.data(lg + 168);
    const auto *lg_180 = buffer.data(lg + 180);
    const auto *lg_183 = buffer.data(lg + 183);
    const auto *lg_185 = buffer.data(lg + 185);
    const auto *lg_190 = buffer.data(lg + 190);
    const auto *lg_192 = buffer.data(lg + 192);
    const auto *lg_194 = buffer.data(lg + 194);
    const auto *lg_195 = buffer.data(lg + 195);
    const auto *lg_200 = buffer.data(lg + 200);
    const auto *lg_224 = buffer.data(lg + 224);
    const auto *lg_235 = buffer.data(lg + 235);
    const auto *lg_265 = buffer.data(lg + 265);
    const auto *lg_267 = buffer.data(lg + 267);
    const auto *lg_269 = buffer.data(lg + 269);

    const auto *md0_57 = buffer.data(md0 + 57);
    const auto *md0_59 = buffer.data(md0 + 59);
    const auto *md0_60 = buffer.data(md0 + 60);
    const auto *md0_63 = buffer.data(md0 + 63);
    const auto *md0_65 = buffer.data(md0 + 65);
    const auto *md0_84 = buffer.data(md0 + 84);
    const auto *md0_87 = buffer.data(md0 + 87);
    const auto *md0_89 = buffer.data(md0 + 89);
    const auto *md0_90 = buffer.data(md0 + 90);
    const auto *md0_93 = buffer.data(md0 + 93);
    const auto *md0_95 = buffer.data(md0 + 95);

    const auto *md1_57 = buffer.data(md1 + 57);
    const auto *md1_59 = buffer.data(md1 + 59);
    const auto *md1_60 = buffer.data(md1 + 60);
    const auto *md1_63 = buffer.data(md1 + 63);
    const auto *md1_65 = buffer.data(md1 + 65);
    const auto *md1_84 = buffer.data(md1 + 84);
    const auto *md1_87 = buffer.data(md1 + 87);
    const auto *md1_89 = buffer.data(md1 + 89);
    const auto *md1_90 = buffer.data(md1 + 90);
    const auto *md1_93 = buffer.data(md1 + 93);
    const auto *md1_95 = buffer.data(md1 + 95);

    const auto *mf_96 = buffer.data(mf + 96);
    const auto *mf_98 = buffer.data(mf + 98);
    const auto *mf_99 = buffer.data(mf + 99);
    const auto *mf_100 = buffer.data(mf + 100);
    const auto *mf_101 = buffer.data(mf + 101);
    const auto *mf_102 = buffer.data(mf + 102);
    const auto *mf_103 = buffer.data(mf + 103);
    const auto *mf_106 = buffer.data(mf + 106);
    const auto *mf_107 = buffer.data(mf + 107);
    const auto *mf_108 = buffer.data(mf + 108);
    const auto *mf_109 = buffer.data(mf + 109);
    const auto *mf_110 = buffer.data(mf + 110);
    const auto *mf_112 = buffer.data(mf + 112);
    const auto *mf_116 = buffer.data(mf + 116);
    const auto *mf_117 = buffer.data(mf + 117);
    const auto *mf_118 = buffer.data(mf + 118);
    const auto *mf_119 = buffer.data(mf + 119);
    const auto *mf_120 = buffer.data(mf + 120);
    const auto *mf_122 = buffer.data(mf + 122);
    const auto *mf_126 = buffer.data(mf + 126);
    const auto *mf_127 = buffer.data(mf + 127);
    const auto *mf_128 = buffer.data(mf + 128);
    const auto *mf_129 = buffer.data(mf + 129);
    const auto *mf_130 = buffer.data(mf + 130);
    const auto *mf_132 = buffer.data(mf + 132);
    const auto *mf_136 = buffer.data(mf + 136);
    const auto *mf_137 = buffer.data(mf + 137);
    const auto *mf_138 = buffer.data(mf + 138);
    const auto *mf_139 = buffer.data(mf + 139);
    const auto *mf_140 = buffer.data(mf + 140);
    const auto *mf_141 = buffer.data(mf + 141);
    const auto *mf_142 = buffer.data(mf + 142);
    const auto *mf_145 = buffer.data(mf + 145);
    const auto *mf_146 = buffer.data(mf + 146);
    const auto *mf_147 = buffer.data(mf + 147);
    const auto *mf_148 = buffer.data(mf + 148);
    const auto *mf_149 = buffer.data(mf + 149);
    const auto *mf_150 = buffer.data(mf + 150);
    const auto *mf_151 = buffer.data(mf + 151);
    const auto *mf_152 = buffer.data(mf + 152);
    const auto *mf_153 = buffer.data(mf + 153);
    const auto *mf_156 = buffer.data(mf + 156);
    const auto *mf_157 = buffer.data(mf + 157);
    const auto *mf_158 = buffer.data(mf + 158);
    const auto *mf_159 = buffer.data(mf + 159);
    const auto *mf_160 = buffer.data(mf + 160);
    const auto *mf_162 = buffer.data(mf + 162);
    const auto *mf_166 = buffer.data(mf + 166);
    const auto *mf_167 = buffer.data(mf + 167);
    const auto *mf_168 = buffer.data(mf + 168);
    const auto *mf_169 = buffer.data(mf + 169);
    const auto *mf_170 = buffer.data(mf + 170);
    const auto *mf_172 = buffer.data(mf + 172);
    const auto *mf_176 = buffer.data(mf + 176);
    const auto *mf_177 = buffer.data(mf + 177);
    const auto *mf_178 = buffer.data(mf + 178);
    const auto *mf_179 = buffer.data(mf + 179);
    const auto *mf_180 = buffer.data(mf + 180);
    const auto *mf_182 = buffer.data(mf + 182);

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_x, pb_y, pb_z, lf_56, lf_99, md0_57, \
                         md0_59, md1_57, md1_59, mf_96, mf_98, mf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_17 * lf_99[k]
                   + pb_x[k] * mf_99[k];

        t_145[k] = f_1 * md0_57[k]
                   - f_2 * md1_57[k]
                   + pb_y[k] * mf_96[k];

        t_146[k] = f_16 * lf_56[k]
                   + pb_z[k] * mf_96[k];

        t_147[k] = f_3 * md0_59[k]
                   - f_4 * md1_59[k]
                   + pb_y[k] * mf_98[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pa_x, pa_y, pb_y, kg0_45, kg0_149, \
                         kg1_45, kg1_149, lf_60, lg_90, lg_149, mf_99, \
                         mf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = pb_y[k] * mf_99[k];

        t_149[k] = f_18 * kg0_149[k]
                   - f_19 * kg1_149[k]
                   + pa_x[k] * lg_149[k];

        t_150[k] = f_20 * kg0_45[k]
                   - f_21 * kg1_45[k]
                   + pa_y[k] * lg_90[k];

        t_151[k] = f_8 * lf_60[k]
                   + pb_y[k] * mf_100[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pb_z, lf_103, md0_60, md0_63, \
                         md1_60, md1_63, mf_100, mf_101, mf_102, \
                         mf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = pb_z[k] * mf_100[k];

        t_153[k] = f_22 * lf_103[k]
                   + f_3 * md0_63[k]
                   - f_4 * md1_63[k]
                   + pb_x[k] * mf_103[k];

        t_154[k] = pb_z[k] * mf_101[k];

        t_155[k] = f_3 * md0_60[k]
                   - f_4 * md1_60[k]
                   + pb_z[k] * mf_102[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_x, pb_z, lf_106, lf_108, lf_109, \
                         mf_103, mf_106, mf_108, mf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_22 * lf_106[k]
                   + pb_x[k] * mf_106[k];

        t_157[k] = pb_z[k] * mf_103[k];

        t_158[k] = f_22 * lf_108[k]
                   + pb_x[k] * mf_108[k];

        t_159[k] = f_22 * lf_109[k]
                   + pb_x[k] * mf_109[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_x, pb_y, pb_z, kg0_160, kg1_160, \
                         lf_69, lg_160, md0_63, md1_63, mf_106, mf_107, \
                         mf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_23 * kg0_160[k]
                   - f_24 * kg1_160[k]
                   + pa_x[k] * lg_160[k];

        t_161[k] = pb_z[k] * mf_106[k];

        t_162[k] = f_3 * md0_63[k]
                   - f_4 * md1_63[k]
                   + pb_z[k] * mf_107[k];

        t_163[k] = f_8 * lf_69[k]
                   + pb_y[k] * mf_109[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, pa_z, pb_z, lf_60, lg_90, lg_91, \
                         lg_93, md0_65, md1_65, mf_109, mf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_1 * md0_65[k]
                   - f_2 * md1_65[k]
                   + pb_z[k] * mf_109[k];

        t_165[k] = pa_z[k] * lg_90[k];

        t_166[k] = pa_z[k] * lg_91[k];

        t_167[k] = f_5 * lf_60[k]
                   + pb_z[k] * mf_110[k];

        t_168[k] = pa_z[k] * lg_93[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pa_z, pb_x, pb_y, lf_62, lf_72, lf_117, \
                         lg_95, lg_96, mf_112, mf_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_16 * lf_72[k]
                   + pb_y[k] * mf_112[k];

        t_170[k] = f_6 * lf_62[k]
                   + pa_z[k] * lg_95[k];

        t_171[k] = pa_z[k] * lg_96[k];

        t_172[k] = f_22 * lf_117[k]
                   + pb_x[k] * mf_117[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pa_z, pb_x, pb_z, lf_66, lf_118, lf_119, \
                         lg_100, mf_116, mf_118, mf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_22 * lf_118[k]
                   + pb_x[k] * mf_118[k];

        t_174[k] = f_22 * lf_119[k]
                   + pb_x[k] * mf_119[k];

        t_175[k] = pa_z[k] * lg_100[k];

        t_176[k] = f_5 * lf_66[k]
                   + pb_z[k] * mf_116[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pa_y, pa_z, pb_y, kg0_75, kg1_75, lf_67, \
                         lf_69, lf_79, lg_102, lg_104, lg_120, mf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_6 * lf_67[k]
                   + pa_z[k] * lg_102[k];

        t_178[k] = f_16 * lf_79[k]
                   + pb_y[k] * mf_119[k];

        t_179[k] = f_8 * lf_69[k]
                   + pa_z[k] * lg_104[k];

        t_180[k] = f_9 * kg0_75[k]
                   - f_10 * kg1_75[k]
                   + pa_y[k] * lg_120[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pa_z, pb_y, pb_z, kg0_48, kg1_48, lf_70, \
                         lf_80, lf_82, lg_108, mf_120, mf_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_6 * lf_80[k]
                   + pb_y[k] * mf_120[k];

        t_182[k] = f_6 * lf_70[k]
                   + pb_z[k] * mf_120[k];

        t_183[k] = f_9 * kg0_48[k]
                   - f_10 * kg1_48[k]
                   + pa_z[k] * lg_108[k];

        t_184[k] = f_6 * lf_82[k]
                   + pb_y[k] * mf_122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pa_y, pb_x, kg0_80, kg1_80, lf_126, \
                         lf_127, lf_128, lg_125, mf_126, mf_127, \
                         mf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_9 * kg0_80[k]
                   - f_10 * kg1_80[k]
                   + pa_y[k] * lg_125[k];

        t_186[k] = f_22 * lf_126[k]
                   + pb_x[k] * mf_126[k];

        t_187[k] = f_22 * lf_127[k]
                   + pb_x[k] * mf_127[k];

        t_188[k] = f_22 * lf_128[k]
                   + pb_x[k] * mf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pa_x, pb_x, pb_z, kg0_190, kg1_190, lf_76, \
                         lf_129, lg_190, mf_126, mf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_22 * lf_129[k]
                   + pb_x[k] * mf_129[k];

        t_190[k] = f_23 * kg0_190[k]
                   - f_24 * kg1_190[k]
                   + pa_x[k] * lg_190[k];

        t_191[k] = f_6 * lf_76[k]
                   + pb_z[k] * mf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_x, pa_y, pb_y, kg0_192, kg0_194, \
                         kg1_192, kg1_194, lf_89, lg_135, lg_192, lg_194, \
                         mf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_23 * kg0_192[k]
                   - f_24 * kg1_192[k]
                   + pa_x[k] * lg_192[k];

        t_193[k] = f_6 * lf_89[k]
                   + pb_y[k] * mf_129[k];

        t_194[k] = f_23 * kg0_194[k]
                   - f_24 * kg1_194[k]
                   + pa_x[k] * lg_194[k];

        t_195[k] = pa_y[k] * lg_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pa_y, pb_y, lf_90, lf_91, lf_92, \
                         lg_137, lg_138, lg_140, mf_130, mf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_5 * lf_90[k]
                   + pb_y[k] * mf_130[k];

        t_197[k] = pa_y[k] * lg_137[k];

        t_198[k] = f_6 * lf_91[k]
                   + pa_y[k] * lg_138[k];

        t_199[k] = f_5 * lf_92[k]
                   + pb_y[k] * mf_132[k];

        t_200[k] = pa_y[k] * lg_140[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, pa_y, pb_x, lf_96, lf_136, lf_137, \
                         lf_138, lg_144, lg_145, mf_136, mf_137, \
                         mf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_22 * lf_136[k]
                   + pb_x[k] * mf_136[k];

        t_202[k] = f_22 * lf_137[k]
                   + pb_x[k] * mf_137[k];

        t_203[k] = f_22 * lf_138[k]
                   + pb_x[k] * mf_138[k];

        t_204[k] = pa_y[k] * lg_144[k];

        t_205[k] = f_8 * lf_96[k]
                   + pa_y[k] * lg_145[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_y, pb_y, pb_z, lf_86, lf_98, lf_99, \
                         lg_147, lg_149, mf_136, mf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_16 * lf_86[k]
                   + pb_z[k] * mf_136[k];

        t_207[k] = f_6 * lf_98[k]
                   + pa_y[k] * lg_147[k];

        t_208[k] = f_5 * lf_99[k]
                   + pb_y[k] * mf_139[k];

        t_209[k] = pa_y[k] * lg_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pa_z, pb_y, pb_z, kg0_75, kg1_75, lf_90, \
                         lg_135, md0_84, md1_84, mf_140, mf_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_20 * kg0_75[k]
                   - f_21 * kg1_75[k]
                   + pa_z[k] * lg_135[k];

        t_211[k] = pb_y[k] * mf_140[k];

        t_212[k] = f_8 * lf_90[k]
                   + pb_z[k] * mf_140[k];

        t_213[k] = f_3 * md0_84[k]
                   - f_4 * md1_84[k]
                   + pb_y[k] * mf_141[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, pb_x, pb_y, lf_145, lf_146, \
                         lf_147, md0_89, md1_89, mf_142, mf_145, mf_146, \
                         mf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = pb_y[k] * mf_142[k];

        t_215[k] = f_22 * lf_145[k]
                   + f_3 * md0_89[k]
                   - f_4 * md1_89[k]
                   + pb_x[k] * mf_145[k];

        t_216[k] = f_22 * lf_146[k]
                   + pb_x[k] * mf_146[k];

        t_217[k] = f_22 * lf_147[k]
                   + pb_x[k] * mf_147[k];

        t_218[k] = pb_y[k] * mf_145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pb_x, pb_y, pb_z, lf_96, lf_149, md0_87, \
                         md0_89, md1_87, md1_89, mf_146, mf_148, \
                         mf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_22 * lf_149[k]
                   + pb_x[k] * mf_149[k];

        t_220[k] = f_1 * md0_87[k]
                   - f_2 * md1_87[k]
                   + pb_y[k] * mf_146[k];

        t_221[k] = f_8 * lf_96[k]
                   + pb_z[k] * mf_146[k];

        t_222[k] = f_3 * md0_89[k]
                   - f_4 * md1_89[k]
                   + pb_y[k] * mf_148[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_x, pa_y, pb_y, kg0_90, kg0_224, \
                         kg1_90, kg1_224, lf_100, lg_150, lg_224, mf_149, \
                         mf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = pb_y[k] * mf_149[k];

        t_224[k] = f_23 * kg0_224[k]
                   - f_24 * kg1_224[k]
                   + pa_x[k] * lg_224[k];

        t_225[k] = f_23 * kg0_90[k]
                   - f_24 * kg1_90[k]
                   + pa_y[k] * lg_150[k];

        t_226[k] = f_22 * lf_100[k]
                   + pb_y[k] * mf_150[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, pb_x, pb_z, lf_153, md0_90, md0_93, \
                         md1_90, md1_93, mf_150, mf_151, mf_152, \
                         mf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = pb_z[k] * mf_150[k];

        t_228[k] = f_8 * lf_153[k]
                   + f_3 * md0_93[k]
                   - f_4 * md1_93[k]
                   + pb_x[k] * mf_153[k];

        t_229[k] = pb_z[k] * mf_151[k];

        t_230[k] = f_3 * md0_90[k]
                   - f_4 * md1_90[k]
                   + pb_z[k] * mf_152[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pb_x, pb_z, lf_156, lf_158, lf_159, \
                         mf_153, mf_156, mf_158, mf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_8 * lf_156[k]
                   + pb_x[k] * mf_156[k];

        t_232[k] = pb_z[k] * mf_153[k];

        t_233[k] = f_8 * lf_158[k]
                   + pb_x[k] * mf_158[k];

        t_234[k] = f_8 * lf_159[k]
                   + pb_x[k] * mf_159[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pa_x, pb_y, pb_z, kg0_235, kg1_235, \
                         lf_109, lg_235, md0_93, md1_93, mf_156, mf_157, \
                         mf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_20 * kg0_235[k]
                   - f_21 * kg1_235[k]
                   + pa_x[k] * lg_235[k];

        t_236[k] = pb_z[k] * mf_156[k];

        t_237[k] = f_3 * md0_93[k]
                   - f_4 * md1_93[k]
                   + pb_z[k] * mf_157[k];

        t_238[k] = f_22 * lf_109[k]
                   + pb_y[k] * mf_159[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pa_z, pb_z, lf_100, lg_150, \
                         lg_151, lg_153, md0_95, md1_95, mf_159, \
                         mf_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * md0_95[k]
                   - f_2 * md1_95[k]
                   + pb_z[k] * mf_159[k];

        t_240[k] = pa_z[k] * lg_150[k];

        t_241[k] = pa_z[k] * lg_151[k];

        t_242[k] = f_5 * lf_100[k]
                   + pb_z[k] * mf_160[k];

        t_243[k] = pa_z[k] * lg_153[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_z, pb_x, pb_y, lf_102, lf_112, lf_167, \
                         lg_155, lg_156, mf_162, mf_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_8 * lf_112[k]
                   + pb_y[k] * mf_162[k];

        t_245[k] = f_6 * lf_102[k]
                   + pa_z[k] * lg_155[k];

        t_246[k] = pa_z[k] * lg_156[k];

        t_247[k] = f_8 * lf_167[k]
                   + pb_x[k] * mf_167[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_z, pb_x, pb_z, lf_106, lf_168, lf_169, \
                         lg_160, mf_166, mf_168, mf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_8 * lf_168[k]
                   + pb_x[k] * mf_168[k];

        t_249[k] = f_8 * lf_169[k]
                   + pb_x[k] * mf_169[k];

        t_250[k] = pa_z[k] * lg_160[k];

        t_251[k] = f_5 * lf_106[k]
                   + pb_z[k] * mf_166[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_y, pa_z, pb_y, kg0_120, kg1_120, \
                         lf_107, lf_109, lf_119, lg_162, lg_164, lg_180, \
                         mf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_6 * lf_107[k]
                   + pa_z[k] * lg_162[k];

        t_253[k] = f_8 * lf_119[k]
                   + pb_y[k] * mf_169[k];

        t_254[k] = f_8 * lf_109[k]
                   + pa_z[k] * lg_164[k];

        t_255[k] = f_14 * kg0_120[k]
                   - f_15 * kg1_120[k]
                   + pa_y[k] * lg_180[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pa_z, pb_y, pb_z, kg0_93, kg1_93, lf_110, \
                         lf_120, lf_122, lg_168, mf_170, mf_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_16 * lf_120[k]
                   + pb_y[k] * mf_170[k];

        t_257[k] = f_6 * lf_110[k]
                   + pb_z[k] * mf_170[k];

        t_258[k] = f_9 * kg0_93[k]
                   - f_10 * kg1_93[k]
                   + pa_z[k] * lg_168[k];

        t_259[k] = f_16 * lf_122[k]
                   + pb_y[k] * mf_172[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pb_x, kg0_125, kg1_125, lf_176, \
                         lf_177, lf_178, lg_185, mf_176, mf_177, \
                         mf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * kg0_125[k]
                   - f_15 * kg1_125[k]
                   + pa_y[k] * lg_185[k];

        t_261[k] = f_8 * lf_176[k]
                   + pb_x[k] * mf_176[k];

        t_262[k] = f_8 * lf_177[k]
                   + pb_x[k] * mf_177[k];

        t_263[k] = f_8 * lf_178[k]
                   + pb_x[k] * mf_178[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pa_x, pb_x, pb_z, kg0_265, kg1_265, lf_116, \
                         lf_179, lg_265, mf_176, mf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_8 * lf_179[k]
                   + pb_x[k] * mf_179[k];

        t_265[k] = f_20 * kg0_265[k]
                   - f_21 * kg1_265[k]
                   + pa_x[k] * lg_265[k];

        t_266[k] = f_6 * lf_116[k]
                   + pb_z[k] * mf_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pa_x, pb_y, kg0_267, kg0_269, kg1_267, kg1_269, \
                         lf_129, lg_267, lg_269, mf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_20 * kg0_267[k]
                   - f_21 * kg1_267[k]
                   + pa_x[k] * lg_267[k];

        t_268[k] = f_16 * lf_129[k]
                   + pb_y[k] * mf_179[k];

        t_269[k] = f_20 * kg0_269[k]
                   - f_21 * kg1_269[k]
                   + pa_x[k] * lg_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pa_y, pb_y, pb_z, kg0_135, kg1_135, lf_120, \
                         lf_130, lg_195, mf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_9 * kg0_135[k]
                   - f_10 * kg1_135[k]
                   + pa_y[k] * lg_195[k];

        t_271[k] = f_6 * lf_130[k]
                   + pb_y[k] * mf_180[k];

        t_272[k] = f_16 * lf_120[k]
                   + pb_z[k] * mf_180[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pa_y, pa_z, pb_y, kg0_108, kg0_140, kg1_108, \
                         kg1_140, lf_132, lg_183, lg_200, mf_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_14 * kg0_108[k]
                   - f_15 * kg1_108[k]
                   + pa_z[k] * lg_183[k];

        t_274[k] = f_6 * lf_132[k]
                   + pb_y[k] * mf_182[k];

        t_275[k] = f_9 * kg0_140[k]
                   - f_10 * kg1_140[k]
                   + pa_y[k] * lg_200[k];
    }
}

static auto
compute_prim_mg_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kg0,
                                            const size_t kg1, const size_t lf, const size_t lg,
                                            const size_t md0, const size_t md1, const size_t mf,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_8 = 2.0 / p;
    const auto f_9 = 0.5 / alpha;
    const auto f_10 = 0.5 * beta / (alpha * p);
    const auto f_14 = 1.0 / alpha;
    const auto f_15 = beta / (alpha * p);
    const auto f_16 = 1.5 / p;
    const auto f_17 = 3.0 / p;
    const auto f_18 = 2.5 / alpha;
    const auto f_19 = 2.5 * beta / (alpha * p);
    const auto f_20 = 1.5 / alpha;
    const auto f_21 = 1.5 * beta / (alpha * p);
    const auto f_22 = 2.5 / p;
    const auto f_23 = 2.0 / alpha;
    const auto f_24 = 2.0 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kg0_135 = buffer.data(kg0 + 135);
    const auto *kg0_150 = buffer.data(kg0 + 150);
    const auto *kg0_153 = buffer.data(kg0 + 153);
    const auto *kg0_168 = buffer.data(kg0 + 168);
    const auto *kg0_180 = buffer.data(kg0 + 180);
    const auto *kg0_183 = buffer.data(kg0 + 183);
    const auto *kg0_185 = buffer.data(kg0 + 185);
    const auto *kg0_195 = buffer.data(kg0 + 195);
    const auto *kg0_200 = buffer.data(kg0 + 200);
    const auto *kg0_210 = buffer.data(kg0 + 210);
    const auto *kg0_215 = buffer.data(kg0 + 215);
    const auto *kg0_280 = buffer.data(kg0 + 280);
    const auto *kg0_282 = buffer.data(kg0 + 282);
    const auto *kg0_284 = buffer.data(kg0 + 284);
    const auto *kg0_314 = buffer.data(kg0 + 314);
    const auto *kg0_325 = buffer.data(kg0 + 325);
    const auto *kg0_355 = buffer.data(kg0 + 355);
    const auto *kg0_357 = buffer.data(kg0 + 357);
    const auto *kg0_359 = buffer.data(kg0 + 359);
    const auto *kg0_370 = buffer.data(kg0 + 370);
    const auto *kg0_372 = buffer.data(kg0 + 372);
    const auto *kg0_374 = buffer.data(kg0 + 374);
    const auto *kg0_385 = buffer.data(kg0 + 385);
    const auto *kg0_387 = buffer.data(kg0 + 387);
    const auto *kg0_389 = buffer.data(kg0 + 389);

    const auto *kg1_135 = buffer.data(kg1 + 135);
    const auto *kg1_150 = buffer.data(kg1 + 150);
    const auto *kg1_153 = buffer.data(kg1 + 153);
    const auto *kg1_168 = buffer.data(kg1 + 168);
    const auto *kg1_180 = buffer.data(kg1 + 180);
    const auto *kg1_183 = buffer.data(kg1 + 183);
    const auto *kg1_185 = buffer.data(kg1 + 185);
    const auto *kg1_195 = buffer.data(kg1 + 195);
    const auto *kg1_200 = buffer.data(kg1 + 200);
    const auto *kg1_210 = buffer.data(kg1 + 210);
    const auto *kg1_215 = buffer.data(kg1 + 215);
    const auto *kg1_280 = buffer.data(kg1 + 280);
    const auto *kg1_282 = buffer.data(kg1 + 282);
    const auto *kg1_284 = buffer.data(kg1 + 284);
    const auto *kg1_314 = buffer.data(kg1 + 314);
    const auto *kg1_325 = buffer.data(kg1 + 325);
    const auto *kg1_355 = buffer.data(kg1 + 355);
    const auto *kg1_357 = buffer.data(kg1 + 357);
    const auto *kg1_359 = buffer.data(kg1 + 359);
    const auto *kg1_370 = buffer.data(kg1 + 370);
    const auto *kg1_372 = buffer.data(kg1 + 372);
    const auto *kg1_374 = buffer.data(kg1 + 374);
    const auto *kg1_385 = buffer.data(kg1 + 385);
    const auto *kg1_387 = buffer.data(kg1 + 387);
    const auto *kg1_389 = buffer.data(kg1 + 389);

    const auto *lf_126 = buffer.data(lf + 126);
    const auto *lf_136 = buffer.data(lf + 136);
    const auto *lf_139 = buffer.data(lf + 139);
    const auto *lf_140 = buffer.data(lf + 140);
    const auto *lf_141 = buffer.data(lf + 141);
    const auto *lf_142 = buffer.data(lf + 142);
    const auto *lf_146 = buffer.data(lf + 146);
    const auto *lf_148 = buffer.data(lf + 148);
    const auto *lf_149 = buffer.data(lf + 149);
    const auto *lf_150 = buffer.data(lf + 150);
    const auto *lf_152 = buffer.data(lf + 152);
    const auto *lf_156 = buffer.data(lf + 156);
    const auto *lf_157 = buffer.data(lf + 157);
    const auto *lf_159 = buffer.data(lf + 159);
    const auto *lf_160 = buffer.data(lf + 160);
    const auto *lf_162 = buffer.data(lf + 162);
    const auto *lf_166 = buffer.data(lf + 166);
    const auto *lf_169 = buffer.data(lf + 169);
    const auto *lf_170 = buffer.data(lf + 170);
    const auto *lf_172 = buffer.data(lf + 172);
    const auto *lf_176 = buffer.data(lf + 176);
    const auto *lf_179 = buffer.data(lf + 179);
    const auto *lf_180 = buffer.data(lf + 180);
    const auto *lf_182 = buffer.data(lf + 182);
    const auto *lf_186 = buffer.data(lf + 186);
    const auto *lf_187 = buffer.data(lf + 187);
    const auto *lf_188 = buffer.data(lf + 188);
    const auto *lf_189 = buffer.data(lf + 189);
    const auto *lf_190 = buffer.data(lf + 190);
    const auto *lf_192 = buffer.data(lf + 192);
    const auto *lf_196 = buffer.data(lf + 196);
    const auto *lf_197 = buffer.data(lf + 197);
    const auto *lf_198 = buffer.data(lf + 198);
    const auto *lf_199 = buffer.data(lf + 199);
    const auto *lf_200 = buffer.data(lf + 200);
    const auto *lf_201 = buffer.data(lf + 201);
    const auto *lf_202 = buffer.data(lf + 202);
    const auto *lf_205 = buffer.data(lf + 205);
    const auto *lf_206 = buffer.data(lf + 206);
    const auto *lf_207 = buffer.data(lf + 207);
    const auto *lf_208 = buffer.data(lf + 208);
    const auto *lf_209 = buffer.data(lf + 209);
    const auto *lf_213 = buffer.data(lf + 213);
    const auto *lf_216 = buffer.data(lf + 216);
    const auto *lf_218 = buffer.data(lf + 218);
    const auto *lf_219 = buffer.data(lf + 219);
    const auto *lf_227 = buffer.data(lf + 227);
    const auto *lf_228 = buffer.data(lf + 228);
    const auto *lf_229 = buffer.data(lf + 229);
    const auto *lf_236 = buffer.data(lf + 236);
    const auto *lf_237 = buffer.data(lf + 237);
    const auto *lf_238 = buffer.data(lf + 238);
    const auto *lf_239 = buffer.data(lf + 239);
    const auto *lf_246 = buffer.data(lf + 246);
    const auto *lf_247 = buffer.data(lf + 247);
    const auto *lf_248 = buffer.data(lf + 248);
    const auto *lf_249 = buffer.data(lf + 249);
    const auto *lf_256 = buffer.data(lf + 256);
    const auto *lf_257 = buffer.data(lf + 257);
    const auto *lf_258 = buffer.data(lf + 258);
    const auto *lf_259 = buffer.data(lf + 259);
    const auto *lf_266 = buffer.data(lf + 266);
    const auto *lf_267 = buffer.data(lf + 267);
    const auto *lf_268 = buffer.data(lf + 268);

    const auto *lg_210 = buffer.data(lg + 210);
    const auto *lg_212 = buffer.data(lg + 212);
    const auto *lg_213 = buffer.data(lg + 213);
    const auto *lg_215 = buffer.data(lg + 215);
    const auto *lg_219 = buffer.data(lg + 219);
    const auto *lg_220 = buffer.data(lg + 220);
    const auto *lg_222 = buffer.data(lg + 222);
    const auto *lg_224 = buffer.data(lg + 224);
    const auto *lg_225 = buffer.data(lg + 225);
    const auto *lg_226 = buffer.data(lg + 226);
    const auto *lg_228 = buffer.data(lg + 228);
    const auto *lg_230 = buffer.data(lg + 230);
    const auto *lg_231 = buffer.data(lg + 231);
    const auto *lg_235 = buffer.data(lg + 235);
    const auto *lg_237 = buffer.data(lg + 237);
    const auto *lg_239 = buffer.data(lg + 239);
    const auto *lg_243 = buffer.data(lg + 243);
    const auto *lg_255 = buffer.data(lg + 255);
    const auto *lg_258 = buffer.data(lg + 258);
    const auto *lg_260 = buffer.data(lg + 260);
    const auto *lg_270 = buffer.data(lg + 270);
    const auto *lg_273 = buffer.data(lg + 273);
    const auto *lg_275 = buffer.data(lg + 275);
    const auto *lg_280 = buffer.data(lg + 280);
    const auto *lg_282 = buffer.data(lg + 282);
    const auto *lg_284 = buffer.data(lg + 284);
    const auto *lg_285 = buffer.data(lg + 285);
    const auto *lg_290 = buffer.data(lg + 290);
    const auto *lg_300 = buffer.data(lg + 300);
    const auto *lg_302 = buffer.data(lg + 302);
    const auto *lg_303 = buffer.data(lg + 303);
    const auto *lg_305 = buffer.data(lg + 305);
    const auto *lg_309 = buffer.data(lg + 309);
    const auto *lg_310 = buffer.data(lg + 310);
    const auto *lg_312 = buffer.data(lg + 312);
    const auto *lg_314 = buffer.data(lg + 314);
    const auto *lg_325 = buffer.data(lg + 325);
    const auto *lg_355 = buffer.data(lg + 355);
    const auto *lg_357 = buffer.data(lg + 357);
    const auto *lg_359 = buffer.data(lg + 359);
    const auto *lg_370 = buffer.data(lg + 370);
    const auto *lg_372 = buffer.data(lg + 372);
    const auto *lg_374 = buffer.data(lg + 374);
    const auto *lg_385 = buffer.data(lg + 385);
    const auto *lg_387 = buffer.data(lg + 387);
    const auto *lg_389 = buffer.data(lg + 389);

    const auto *md0_120 = buffer.data(md0 + 120);
    const auto *md0_123 = buffer.data(md0 + 123);
    const auto *md0_125 = buffer.data(md0 + 125);
    const auto *md0_126 = buffer.data(md0 + 126);
    const auto *md0_129 = buffer.data(md0 + 129);
    const auto *md0_131 = buffer.data(md0 + 131);

    const auto *md1_120 = buffer.data(md1 + 120);
    const auto *md1_123 = buffer.data(md1 + 123);
    const auto *md1_125 = buffer.data(md1 + 125);
    const auto *md1_126 = buffer.data(md1 + 126);
    const auto *md1_129 = buffer.data(md1 + 129);
    const auto *md1_131 = buffer.data(md1 + 131);

    const auto *mf_186 = buffer.data(mf + 186);
    const auto *mf_187 = buffer.data(mf + 187);
    const auto *mf_188 = buffer.data(mf + 188);
    const auto *mf_189 = buffer.data(mf + 189);
    const auto *mf_190 = buffer.data(mf + 190);
    const auto *mf_192 = buffer.data(mf + 192);
    const auto *mf_196 = buffer.data(mf + 196);
    const auto *mf_197 = buffer.data(mf + 197);
    const auto *mf_198 = buffer.data(mf + 198);
    const auto *mf_199 = buffer.data(mf + 199);
    const auto *mf_200 = buffer.data(mf + 200);
    const auto *mf_201 = buffer.data(mf + 201);
    const auto *mf_202 = buffer.data(mf + 202);
    const auto *mf_205 = buffer.data(mf + 205);
    const auto *mf_206 = buffer.data(mf + 206);
    const auto *mf_207 = buffer.data(mf + 207);
    const auto *mf_208 = buffer.data(mf + 208);
    const auto *mf_209 = buffer.data(mf + 209);
    const auto *mf_210 = buffer.data(mf + 210);
    const auto *mf_211 = buffer.data(mf + 211);
    const auto *mf_212 = buffer.data(mf + 212);
    const auto *mf_213 = buffer.data(mf + 213);
    const auto *mf_216 = buffer.data(mf + 216);
    const auto *mf_217 = buffer.data(mf + 217);
    const auto *mf_218 = buffer.data(mf + 218);
    const auto *mf_219 = buffer.data(mf + 219);
    const auto *mf_220 = buffer.data(mf + 220);
    const auto *mf_222 = buffer.data(mf + 222);
    const auto *mf_226 = buffer.data(mf + 226);
    const auto *mf_227 = buffer.data(mf + 227);
    const auto *mf_228 = buffer.data(mf + 228);
    const auto *mf_229 = buffer.data(mf + 229);
    const auto *mf_230 = buffer.data(mf + 230);
    const auto *mf_232 = buffer.data(mf + 232);
    const auto *mf_236 = buffer.data(mf + 236);
    const auto *mf_237 = buffer.data(mf + 237);
    const auto *mf_238 = buffer.data(mf + 238);
    const auto *mf_239 = buffer.data(mf + 239);
    const auto *mf_240 = buffer.data(mf + 240);
    const auto *mf_242 = buffer.data(mf + 242);
    const auto *mf_246 = buffer.data(mf + 246);
    const auto *mf_247 = buffer.data(mf + 247);
    const auto *mf_248 = buffer.data(mf + 248);
    const auto *mf_249 = buffer.data(mf + 249);
    const auto *mf_250 = buffer.data(mf + 250);
    const auto *mf_252 = buffer.data(mf + 252);
    const auto *mf_256 = buffer.data(mf + 256);
    const auto *mf_257 = buffer.data(mf + 257);
    const auto *mf_258 = buffer.data(mf + 258);
    const auto *mf_259 = buffer.data(mf + 259);
    const auto *mf_260 = buffer.data(mf + 260);
    const auto *mf_262 = buffer.data(mf + 262);
    const auto *mf_266 = buffer.data(mf + 266);
    const auto *mf_267 = buffer.data(mf + 267);
    const auto *mf_268 = buffer.data(mf + 268);
    const auto *mf_269 = buffer.data(mf + 269);

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pb_x, lf_186, lf_187, lf_188, lf_189, \
                         mf_186, mf_187, mf_188, mf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_8 * lf_186[k]
                   + pb_x[k] * mf_186[k];

        t_277[k] = f_8 * lf_187[k]
                   + pb_x[k] * mf_187[k];

        t_278[k] = f_8 * lf_188[k]
                   + pb_x[k] * mf_188[k];

        t_279[k] = f_8 * lf_189[k]
                   + pb_x[k] * mf_189[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pa_x, pb_z, kg0_280, kg0_282, kg1_280, kg1_282, \
                         lf_126, lg_280, lg_282, mf_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_20 * kg0_280[k]
                   - f_21 * kg1_280[k]
                   + pa_x[k] * lg_280[k];

        t_281[k] = f_16 * lf_126[k]
                   + pb_z[k] * mf_186[k];

        t_282[k] = f_20 * kg0_282[k]
                   - f_21 * kg1_282[k]
                   + pa_x[k] * lg_282[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_x, pa_y, pb_y, kg0_284, kg1_284, \
                         lf_139, lf_140, lg_210, lg_284, mf_189, \
                         mf_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_6 * lf_139[k]
                   + pb_y[k] * mf_189[k];

        t_284[k] = f_20 * kg0_284[k]
                   - f_21 * kg1_284[k]
                   + pa_x[k] * lg_284[k];

        t_285[k] = pa_y[k] * lg_210[k];

        t_286[k] = f_5 * lf_140[k]
                   + pb_y[k] * mf_190[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, pa_y, pb_x, pb_y, lf_141, lf_142, \
                         lf_196, lg_212, lg_213, lg_215, mf_192, \
                         mf_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = pa_y[k] * lg_212[k];

        t_288[k] = f_6 * lf_141[k]
                   + pa_y[k] * lg_213[k];

        t_289[k] = f_5 * lf_142[k]
                   + pb_y[k] * mf_192[k];

        t_290[k] = pa_y[k] * lg_215[k];

        t_291[k] = f_8 * lf_196[k]
                   + pb_x[k] * mf_196[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pa_y, pb_x, lf_146, lf_197, lf_198, \
                         lg_219, lg_220, mf_197, mf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_8 * lf_197[k]
                   + pb_x[k] * mf_197[k];

        t_293[k] = f_8 * lf_198[k]
                   + pb_x[k] * mf_198[k];

        t_294[k] = pa_y[k] * lg_219[k];

        t_295[k] = f_8 * lf_146[k]
                   + pa_y[k] * lg_220[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pa_y, pb_y, pb_z, lf_136, lf_148, lf_149, \
                         lg_222, lg_224, mf_196, mf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_8 * lf_136[k]
                   + pb_z[k] * mf_196[k];

        t_297[k] = f_6 * lf_148[k]
                   + pa_y[k] * lg_222[k];

        t_298[k] = f_5 * lf_149[k]
                   + pb_y[k] * mf_199[k];

        t_299[k] = pa_y[k] * lg_224[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_z, pb_y, pb_z, kg0_135, kg1_135, \
                         lf_140, lg_210, md0_120, md1_120, mf_200, \
                         mf_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_23 * kg0_135[k]
                   - f_24 * kg1_135[k]
                   + pa_z[k] * lg_210[k];

        t_301[k] = pb_y[k] * mf_200[k];

        t_302[k] = f_22 * lf_140[k]
                   + pb_z[k] * mf_200[k];

        t_303[k] = f_3 * md0_120[k]
                   - f_4 * md1_120[k]
                   + pb_y[k] * mf_201[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pb_x, pb_y, lf_205, lf_206, \
                         lf_207, md0_125, md1_125, mf_202, mf_205, mf_206, \
                         mf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = pb_y[k] * mf_202[k];

        t_305[k] = f_8 * lf_205[k]
                   + f_3 * md0_125[k]
                   - f_4 * md1_125[k]
                   + pb_x[k] * mf_205[k];

        t_306[k] = f_8 * lf_206[k]
                   + pb_x[k] * mf_206[k];

        t_307[k] = f_8 * lf_207[k]
                   + pb_x[k] * mf_207[k];

        t_308[k] = pb_y[k] * mf_205[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pb_x, pb_y, pb_z, lf_146, lf_209, \
                         md0_123, md0_125, md1_123, md1_125, mf_206, mf_208, \
                         mf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_8 * lf_209[k]
                   + pb_x[k] * mf_209[k];

        t_310[k] = f_1 * md0_123[k]
                   - f_2 * md1_123[k]
                   + pb_y[k] * mf_206[k];

        t_311[k] = f_22 * lf_146[k]
                   + pb_z[k] * mf_206[k];

        t_312[k] = f_3 * md0_125[k]
                   - f_4 * md1_125[k]
                   + pb_y[k] * mf_208[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pa_x, pa_y, pb_y, kg0_150, kg0_314, \
                         kg1_150, kg1_314, lf_150, lg_225, lg_314, mf_209, \
                         mf_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = pb_y[k] * mf_209[k];

        t_314[k] = f_20 * kg0_314[k]
                   - f_21 * kg1_314[k]
                   + pa_x[k] * lg_314[k];

        t_315[k] = f_18 * kg0_150[k]
                   - f_19 * kg1_150[k]
                   + pa_y[k] * lg_225[k];

        t_316[k] = f_17 * lf_150[k]
                   + pb_y[k] * mf_210[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, pb_x, pb_z, lf_213, md0_126, md0_129, \
                         md1_126, md1_129, mf_210, mf_211, mf_212, \
                         mf_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = pb_z[k] * mf_210[k];

        t_318[k] = f_16 * lf_213[k]
                   + f_3 * md0_129[k]
                   - f_4 * md1_129[k]
                   + pb_x[k] * mf_213[k];

        t_319[k] = pb_z[k] * mf_211[k];

        t_320[k] = f_3 * md0_126[k]
                   - f_4 * md1_126[k]
                   + pb_z[k] * mf_212[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pb_x, pb_z, lf_216, lf_218, lf_219, \
                         mf_213, mf_216, mf_218, mf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_16 * lf_216[k]
                   + pb_x[k] * mf_216[k];

        t_322[k] = pb_z[k] * mf_213[k];

        t_323[k] = f_16 * lf_218[k]
                   + pb_x[k] * mf_218[k];

        t_324[k] = f_16 * lf_219[k]
                   + pb_x[k] * mf_219[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pa_x, pb_y, pb_z, kg0_325, kg1_325, \
                         lf_159, lg_325, md0_129, md1_129, mf_216, mf_217, \
                         mf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_14 * kg0_325[k]
                   - f_15 * kg1_325[k]
                   + pa_x[k] * lg_325[k];

        t_326[k] = pb_z[k] * mf_216[k];

        t_327[k] = f_3 * md0_129[k]
                   - f_4 * md1_129[k]
                   + pb_z[k] * mf_217[k];

        t_328[k] = f_17 * lf_159[k]
                   + pb_y[k] * mf_219[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, t_333, pa_z, pb_z, lf_150, lg_225, \
                         lg_226, lg_228, md0_131, md1_131, mf_219, \
                         mf_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_1 * md0_131[k]
                   - f_2 * md1_131[k]
                   + pb_z[k] * mf_219[k];

        t_330[k] = pa_z[k] * lg_225[k];

        t_331[k] = pa_z[k] * lg_226[k];

        t_332[k] = f_5 * lf_150[k]
                   + pb_z[k] * mf_220[k];

        t_333[k] = pa_z[k] * lg_228[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pa_z, pb_x, pb_y, lf_152, lf_162, lf_227, \
                         lg_230, lg_231, mf_222, mf_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_22 * lf_162[k]
                   + pb_y[k] * mf_222[k];

        t_335[k] = f_6 * lf_152[k]
                   + pa_z[k] * lg_230[k];

        t_336[k] = pa_z[k] * lg_231[k];

        t_337[k] = f_16 * lf_227[k]
                   + pb_x[k] * mf_227[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, pa_z, pb_x, pb_z, lf_156, lf_228, lf_229, \
                         lg_235, mf_226, mf_228, mf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_16 * lf_228[k]
                   + pb_x[k] * mf_228[k];

        t_339[k] = f_16 * lf_229[k]
                   + pb_x[k] * mf_229[k];

        t_340[k] = pa_z[k] * lg_235[k];

        t_341[k] = f_5 * lf_156[k]
                   + pb_z[k] * mf_226[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pa_y, pa_z, pb_y, kg0_180, kg1_180, \
                         lf_157, lf_159, lf_169, lg_237, lg_239, lg_255, \
                         mf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_6 * lf_157[k]
                   + pa_z[k] * lg_237[k];

        t_343[k] = f_22 * lf_169[k]
                   + pb_y[k] * mf_229[k];

        t_344[k] = f_8 * lf_159[k]
                   + pa_z[k] * lg_239[k];

        t_345[k] = f_20 * kg0_180[k]
                   - f_21 * kg1_180[k]
                   + pa_y[k] * lg_255[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pa_z, pb_y, pb_z, kg0_153, kg1_153, \
                         lf_160, lf_170, lf_172, lg_243, mf_230, \
                         mf_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_8 * lf_170[k]
                   + pb_y[k] * mf_230[k];

        t_347[k] = f_6 * lf_160[k]
                   + pb_z[k] * mf_230[k];

        t_348[k] = f_9 * kg0_153[k]
                   - f_10 * kg1_153[k]
                   + pa_z[k] * lg_243[k];

        t_349[k] = f_8 * lf_172[k]
                   + pb_y[k] * mf_232[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pa_y, pb_x, kg0_185, kg1_185, lf_236, \
                         lf_237, lf_238, lg_260, mf_236, mf_237, \
                         mf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_20 * kg0_185[k]
                   - f_21 * kg1_185[k]
                   + pa_y[k] * lg_260[k];

        t_351[k] = f_16 * lf_236[k]
                   + pb_x[k] * mf_236[k];

        t_352[k] = f_16 * lf_237[k]
                   + pb_x[k] * mf_237[k];

        t_353[k] = f_16 * lf_238[k]
                   + pb_x[k] * mf_238[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pa_x, pb_x, pb_z, kg0_355, kg1_355, lf_166, \
                         lf_239, lg_355, mf_236, mf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_16 * lf_239[k]
                   + pb_x[k] * mf_239[k];

        t_355[k] = f_14 * kg0_355[k]
                   - f_15 * kg1_355[k]
                   + pa_x[k] * lg_355[k];

        t_356[k] = f_6 * lf_166[k]
                   + pb_z[k] * mf_236[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pa_x, pb_y, kg0_357, kg0_359, kg1_357, kg1_359, \
                         lf_179, lg_357, lg_359, mf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_14 * kg0_357[k]
                   - f_15 * kg1_357[k]
                   + pa_x[k] * lg_357[k];

        t_358[k] = f_8 * lf_179[k]
                   + pb_y[k] * mf_239[k];

        t_359[k] = f_14 * kg0_359[k]
                   - f_15 * kg1_359[k]
                   + pa_x[k] * lg_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pa_y, pb_y, pb_z, kg0_195, kg1_195, lf_170, \
                         lf_180, lg_270, mf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_14 * kg0_195[k]
                   - f_15 * kg1_195[k]
                   + pa_y[k] * lg_270[k];

        t_361[k] = f_16 * lf_180[k]
                   + pb_y[k] * mf_240[k];

        t_362[k] = f_16 * lf_170[k]
                   + pb_z[k] * mf_240[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pa_y, pa_z, pb_y, kg0_168, kg0_200, kg1_168, \
                         kg1_200, lf_182, lg_258, lg_275, mf_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_14 * kg0_168[k]
                   - f_15 * kg1_168[k]
                   + pa_z[k] * lg_258[k];

        t_364[k] = f_16 * lf_182[k]
                   + pb_y[k] * mf_242[k];

        t_365[k] = f_14 * kg0_200[k]
                   - f_15 * kg1_200[k]
                   + pa_y[k] * lg_275[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_x, lf_246, lf_247, lf_248, lf_249, \
                         mf_246, mf_247, mf_248, mf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_16 * lf_246[k]
                   + pb_x[k] * mf_246[k];

        t_367[k] = f_16 * lf_247[k]
                   + pb_x[k] * mf_247[k];

        t_368[k] = f_16 * lf_248[k]
                   + pb_x[k] * mf_248[k];

        t_369[k] = f_16 * lf_249[k]
                   + pb_x[k] * mf_249[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pa_x, pb_z, kg0_370, kg0_372, kg1_370, kg1_372, \
                         lf_176, lg_370, lg_372, mf_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_14 * kg0_370[k]
                   - f_15 * kg1_370[k]
                   + pa_x[k] * lg_370[k];

        t_371[k] = f_16 * lf_176[k]
                   + pb_z[k] * mf_246[k];

        t_372[k] = f_14 * kg0_372[k]
                   - f_15 * kg1_372[k]
                   + pa_x[k] * lg_372[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pa_x, pa_y, pb_y, kg0_210, kg0_374, kg1_210, \
                         kg1_374, lf_189, lg_285, lg_374, mf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_16 * lf_189[k]
                   + pb_y[k] * mf_249[k];

        t_374[k] = f_14 * kg0_374[k]
                   - f_15 * kg1_374[k]
                   + pa_x[k] * lg_374[k];

        t_375[k] = f_9 * kg0_210[k]
                   - f_10 * kg1_210[k]
                   + pa_y[k] * lg_285[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pa_z, pb_y, pb_z, kg0_183, kg1_183, \
                         lf_180, lf_190, lf_192, lg_273, mf_250, \
                         mf_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_6 * lf_190[k]
                   + pb_y[k] * mf_250[k];

        t_377[k] = f_8 * lf_180[k]
                   + pb_z[k] * mf_250[k];

        t_378[k] = f_20 * kg0_183[k]
                   - f_21 * kg1_183[k]
                   + pa_z[k] * lg_273[k];

        t_379[k] = f_6 * lf_192[k]
                   + pb_y[k] * mf_252[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pa_y, pb_x, kg0_215, kg1_215, lf_256, \
                         lf_257, lf_258, lg_290, mf_256, mf_257, \
                         mf_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_9 * kg0_215[k]
                   - f_10 * kg1_215[k]
                   + pa_y[k] * lg_290[k];

        t_381[k] = f_16 * lf_256[k]
                   + pb_x[k] * mf_256[k];

        t_382[k] = f_16 * lf_257[k]
                   + pb_x[k] * mf_257[k];

        t_383[k] = f_16 * lf_258[k]
                   + pb_x[k] * mf_258[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pa_x, pb_x, pb_z, kg0_385, kg1_385, lf_186, \
                         lf_259, lg_385, mf_256, mf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_16 * lf_259[k]
                   + pb_x[k] * mf_259[k];

        t_385[k] = f_14 * kg0_385[k]
                   - f_15 * kg1_385[k]
                   + pa_x[k] * lg_385[k];

        t_386[k] = f_8 * lf_186[k]
                   + pb_z[k] * mf_256[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_x, pa_y, pb_y, kg0_387, kg0_389, \
                         kg1_387, kg1_389, lf_199, lg_300, lg_387, lg_389, \
                         mf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_14 * kg0_387[k]
                   - f_15 * kg1_387[k]
                   + pa_x[k] * lg_387[k];

        t_388[k] = f_6 * lf_199[k]
                   + pb_y[k] * mf_259[k];

        t_389[k] = f_14 * kg0_389[k]
                   - f_15 * kg1_389[k]
                   + pa_x[k] * lg_389[k];

        t_390[k] = pa_y[k] * lg_300[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, t_395, pa_y, pb_y, lf_200, lf_201, \
                         lf_202, lg_302, lg_303, lg_305, mf_260, \
                         mf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_5 * lf_200[k]
                   + pb_y[k] * mf_260[k];

        t_392[k] = pa_y[k] * lg_302[k];

        t_393[k] = f_6 * lf_201[k]
                   + pa_y[k] * lg_303[k];

        t_394[k] = f_5 * lf_202[k]
                   + pb_y[k] * mf_262[k];

        t_395[k] = pa_y[k] * lg_305[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, t_400, pa_y, pb_x, lf_206, lf_266, \
                         lf_267, lf_268, lg_309, lg_310, mf_266, mf_267, \
                         mf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_16 * lf_266[k]
                   + pb_x[k] * mf_266[k];

        t_397[k] = f_16 * lf_267[k]
                   + pb_x[k] * mf_267[k];

        t_398[k] = f_16 * lf_268[k]
                   + pb_x[k] * mf_268[k];

        t_399[k] = pa_y[k] * lg_309[k];

        t_400[k] = f_8 * lf_206[k]
                   + pa_y[k] * lg_310[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_y, pb_y, pb_z, lf_196, lf_208, lf_209, \
                         lg_312, lg_314, mf_266, mf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_22 * lf_196[k]
                   + pb_z[k] * mf_266[k];

        t_402[k] = f_6 * lf_208[k]
                   + pa_y[k] * lg_312[k];

        t_403[k] = f_5 * lf_209[k]
                   + pb_y[k] * mf_269[k];

        t_404[k] = pa_y[k] * lg_314[k];
    }
}

static auto
compute_prim_mg_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kg0,
                                            const size_t kg1, const size_t lf, const size_t lg,
                                            const size_t md0, const size_t md1, const size_t mf,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_8 = 2.0 / p;
    const auto f_9 = 0.5 / alpha;
    const auto f_10 = 0.5 * beta / (alpha * p);
    const auto f_11 = 3.5 / p;
    const auto f_12 = 3.0 / alpha;
    const auto f_13 = 3.0 * beta / (alpha * p);
    const auto f_14 = 1.0 / alpha;
    const auto f_15 = beta / (alpha * p);
    const auto f_16 = 1.5 / p;
    const auto f_17 = 3.0 / p;
    const auto f_18 = 2.5 / alpha;
    const auto f_19 = 2.5 * beta / (alpha * p);
    const auto f_20 = 1.5 / alpha;
    const auto f_21 = 1.5 * beta / (alpha * p);
    const auto f_22 = 2.5 / p;
    const auto f_23 = 2.0 / alpha;
    const auto f_24 = 2.0 * beta / (alpha * p);

    auto *t_405 = buffer.data(target + 405);
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
    auto *t_533 = buffer.data(target + 533);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kg0_210 = buffer.data(kg0 + 210);
    const auto *kg0_225 = buffer.data(kg0 + 225);
    const auto *kg0_228 = buffer.data(kg0 + 228);
    const auto *kg0_243 = buffer.data(kg0 + 243);
    const auto *kg0_255 = buffer.data(kg0 + 255);
    const auto *kg0_258 = buffer.data(kg0 + 258);
    const auto *kg0_260 = buffer.data(kg0 + 260);
    const auto *kg0_270 = buffer.data(kg0 + 270);
    const auto *kg0_273 = buffer.data(kg0 + 273);
    const auto *kg0_275 = buffer.data(kg0 + 275);
    const auto *kg0_285 = buffer.data(kg0 + 285);
    const auto *kg0_290 = buffer.data(kg0 + 290);
    const auto *kg0_300 = buffer.data(kg0 + 300);
    const auto *kg0_305 = buffer.data(kg0 + 305);
    const auto *kg0_419 = buffer.data(kg0 + 419);
    const auto *kg0_430 = buffer.data(kg0 + 430);
    const auto *kg0_460 = buffer.data(kg0 + 460);
    const auto *kg0_462 = buffer.data(kg0 + 462);
    const auto *kg0_464 = buffer.data(kg0 + 464);
    const auto *kg0_475 = buffer.data(kg0 + 475);
    const auto *kg0_477 = buffer.data(kg0 + 477);
    const auto *kg0_479 = buffer.data(kg0 + 479);
    const auto *kg0_490 = buffer.data(kg0 + 490);
    const auto *kg0_492 = buffer.data(kg0 + 492);
    const auto *kg0_494 = buffer.data(kg0 + 494);
    const auto *kg0_505 = buffer.data(kg0 + 505);
    const auto *kg0_507 = buffer.data(kg0 + 507);
    const auto *kg0_509 = buffer.data(kg0 + 509);

    const auto *kg1_210 = buffer.data(kg1 + 210);
    const auto *kg1_225 = buffer.data(kg1 + 225);
    const auto *kg1_228 = buffer.data(kg1 + 228);
    const auto *kg1_243 = buffer.data(kg1 + 243);
    const auto *kg1_255 = buffer.data(kg1 + 255);
    const auto *kg1_258 = buffer.data(kg1 + 258);
    const auto *kg1_260 = buffer.data(kg1 + 260);
    const auto *kg1_270 = buffer.data(kg1 + 270);
    const auto *kg1_273 = buffer.data(kg1 + 273);
    const auto *kg1_275 = buffer.data(kg1 + 275);
    const auto *kg1_285 = buffer.data(kg1 + 285);
    const auto *kg1_290 = buffer.data(kg1 + 290);
    const auto *kg1_300 = buffer.data(kg1 + 300);
    const auto *kg1_305 = buffer.data(kg1 + 305);
    const auto *kg1_419 = buffer.data(kg1 + 419);
    const auto *kg1_430 = buffer.data(kg1 + 430);
    const auto *kg1_460 = buffer.data(kg1 + 460);
    const auto *kg1_462 = buffer.data(kg1 + 462);
    const auto *kg1_464 = buffer.data(kg1 + 464);
    const auto *kg1_475 = buffer.data(kg1 + 475);
    const auto *kg1_477 = buffer.data(kg1 + 477);
    const auto *kg1_479 = buffer.data(kg1 + 479);
    const auto *kg1_490 = buffer.data(kg1 + 490);
    const auto *kg1_492 = buffer.data(kg1 + 492);
    const auto *kg1_494 = buffer.data(kg1 + 494);
    const auto *kg1_505 = buffer.data(kg1 + 505);
    const auto *kg1_507 = buffer.data(kg1 + 507);
    const auto *kg1_509 = buffer.data(kg1 + 509);

    const auto *lf_200 = buffer.data(lf + 200);
    const auto *lf_206 = buffer.data(lf + 206);
    const auto *lf_210 = buffer.data(lf + 210);
    const auto *lf_212 = buffer.data(lf + 212);
    const auto *lf_216 = buffer.data(lf + 216);
    const auto *lf_217 = buffer.data(lf + 217);
    const auto *lf_219 = buffer.data(lf + 219);
    const auto *lf_220 = buffer.data(lf + 220);
    const auto *lf_222 = buffer.data(lf + 222);
    const auto *lf_226 = buffer.data(lf + 226);
    const auto *lf_229 = buffer.data(lf + 229);
    const auto *lf_230 = buffer.data(lf + 230);
    const auto *lf_232 = buffer.data(lf + 232);
    const auto *lf_236 = buffer.data(lf + 236);
    const auto *lf_239 = buffer.data(lf + 239);
    const auto *lf_240 = buffer.data(lf + 240);
    const auto *lf_242 = buffer.data(lf + 242);
    const auto *lf_246 = buffer.data(lf + 246);
    const auto *lf_249 = buffer.data(lf + 249);
    const auto *lf_250 = buffer.data(lf + 250);
    const auto *lf_252 = buffer.data(lf + 252);
    const auto *lf_256 = buffer.data(lf + 256);
    const auto *lf_259 = buffer.data(lf + 259);
    const auto *lf_260 = buffer.data(lf + 260);
    const auto *lf_262 = buffer.data(lf + 262);
    const auto *lf_266 = buffer.data(lf + 266);
    const auto *lf_269 = buffer.data(lf + 269);
    const auto *lf_270 = buffer.data(lf + 270);
    const auto *lf_271 = buffer.data(lf + 271);
    const auto *lf_272 = buffer.data(lf + 272);
    const auto *lf_275 = buffer.data(lf + 275);
    const auto *lf_276 = buffer.data(lf + 276);
    const auto *lf_277 = buffer.data(lf + 277);
    const auto *lf_278 = buffer.data(lf + 278);
    const auto *lf_279 = buffer.data(lf + 279);
    const auto *lf_283 = buffer.data(lf + 283);
    const auto *lf_286 = buffer.data(lf + 286);
    const auto *lf_288 = buffer.data(lf + 288);
    const auto *lf_289 = buffer.data(lf + 289);
    const auto *lf_297 = buffer.data(lf + 297);
    const auto *lf_298 = buffer.data(lf + 298);
    const auto *lf_299 = buffer.data(lf + 299);
    const auto *lf_306 = buffer.data(lf + 306);
    const auto *lf_307 = buffer.data(lf + 307);
    const auto *lf_308 = buffer.data(lf + 308);
    const auto *lf_309 = buffer.data(lf + 309);
    const auto *lf_316 = buffer.data(lf + 316);
    const auto *lf_317 = buffer.data(lf + 317);
    const auto *lf_318 = buffer.data(lf + 318);
    const auto *lf_319 = buffer.data(lf + 319);
    const auto *lf_326 = buffer.data(lf + 326);
    const auto *lf_327 = buffer.data(lf + 327);
    const auto *lf_328 = buffer.data(lf + 328);
    const auto *lf_329 = buffer.data(lf + 329);
    const auto *lf_336 = buffer.data(lf + 336);
    const auto *lf_337 = buffer.data(lf + 337);
    const auto *lf_338 = buffer.data(lf + 338);
    const auto *lf_339 = buffer.data(lf + 339);
    const auto *lf_346 = buffer.data(lf + 346);
    const auto *lf_347 = buffer.data(lf + 347);
    const auto *lf_348 = buffer.data(lf + 348);
    const auto *lf_355 = buffer.data(lf + 355);
    const auto *lf_356 = buffer.data(lf + 356);
    const auto *lf_357 = buffer.data(lf + 357);

    const auto *lg_300 = buffer.data(lg + 300);
    const auto *lg_315 = buffer.data(lg + 315);
    const auto *lg_316 = buffer.data(lg + 316);
    const auto *lg_318 = buffer.data(lg + 318);
    const auto *lg_320 = buffer.data(lg + 320);
    const auto *lg_321 = buffer.data(lg + 321);
    const auto *lg_325 = buffer.data(lg + 325);
    const auto *lg_327 = buffer.data(lg + 327);
    const auto *lg_329 = buffer.data(lg + 329);
    const auto *lg_333 = buffer.data(lg + 333);
    const auto *lg_345 = buffer.data(lg + 345);
    const auto *lg_348 = buffer.data(lg + 348);
    const auto *lg_350 = buffer.data(lg + 350);
    const auto *lg_360 = buffer.data(lg + 360);
    const auto *lg_363 = buffer.data(lg + 363);
    const auto *lg_365 = buffer.data(lg + 365);
    const auto *lg_375 = buffer.data(lg + 375);
    const auto *lg_378 = buffer.data(lg + 378);
    const auto *lg_380 = buffer.data(lg + 380);
    const auto *lg_390 = buffer.data(lg + 390);
    const auto *lg_395 = buffer.data(lg + 395);
    const auto *lg_405 = buffer.data(lg + 405);
    const auto *lg_407 = buffer.data(lg + 407);
    const auto *lg_408 = buffer.data(lg + 408);
    const auto *lg_410 = buffer.data(lg + 410);
    const auto *lg_414 = buffer.data(lg + 414);
    const auto *lg_415 = buffer.data(lg + 415);
    const auto *lg_417 = buffer.data(lg + 417);
    const auto *lg_419 = buffer.data(lg + 419);
    const auto *lg_430 = buffer.data(lg + 430);
    const auto *lg_460 = buffer.data(lg + 460);
    const auto *lg_462 = buffer.data(lg + 462);
    const auto *lg_464 = buffer.data(lg + 464);
    const auto *lg_475 = buffer.data(lg + 475);
    const auto *lg_477 = buffer.data(lg + 477);
    const auto *lg_479 = buffer.data(lg + 479);
    const auto *lg_490 = buffer.data(lg + 490);
    const auto *lg_492 = buffer.data(lg + 492);
    const auto *lg_494 = buffer.data(lg + 494);
    const auto *lg_505 = buffer.data(lg + 505);
    const auto *lg_507 = buffer.data(lg + 507);
    const auto *lg_509 = buffer.data(lg + 509);

    const auto *md0_162 = buffer.data(md0 + 162);
    const auto *md0_165 = buffer.data(md0 + 165);
    const auto *md0_167 = buffer.data(md0 + 167);
    const auto *md0_168 = buffer.data(md0 + 168);
    const auto *md0_171 = buffer.data(md0 + 171);
    const auto *md0_173 = buffer.data(md0 + 173);
    const auto *md0_210 = buffer.data(md0 + 210);
    const auto *md0_215 = buffer.data(md0 + 215);

    const auto *md1_162 = buffer.data(md1 + 162);
    const auto *md1_165 = buffer.data(md1 + 165);
    const auto *md1_167 = buffer.data(md1 + 167);
    const auto *md1_168 = buffer.data(md1 + 168);
    const auto *md1_171 = buffer.data(md1 + 171);
    const auto *md1_173 = buffer.data(md1 + 173);
    const auto *md1_210 = buffer.data(md1 + 210);
    const auto *md1_215 = buffer.data(md1 + 215);

    const auto *mf_270 = buffer.data(mf + 270);
    const auto *mf_271 = buffer.data(mf + 271);
    const auto *mf_272 = buffer.data(mf + 272);
    const auto *mf_275 = buffer.data(mf + 275);
    const auto *mf_276 = buffer.data(mf + 276);
    const auto *mf_277 = buffer.data(mf + 277);
    const auto *mf_278 = buffer.data(mf + 278);
    const auto *mf_279 = buffer.data(mf + 279);
    const auto *mf_280 = buffer.data(mf + 280);
    const auto *mf_281 = buffer.data(mf + 281);
    const auto *mf_282 = buffer.data(mf + 282);
    const auto *mf_283 = buffer.data(mf + 283);
    const auto *mf_286 = buffer.data(mf + 286);
    const auto *mf_287 = buffer.data(mf + 287);
    const auto *mf_288 = buffer.data(mf + 288);
    const auto *mf_289 = buffer.data(mf + 289);
    const auto *mf_290 = buffer.data(mf + 290);
    const auto *mf_292 = buffer.data(mf + 292);
    const auto *mf_296 = buffer.data(mf + 296);
    const auto *mf_297 = buffer.data(mf + 297);
    const auto *mf_298 = buffer.data(mf + 298);
    const auto *mf_299 = buffer.data(mf + 299);
    const auto *mf_300 = buffer.data(mf + 300);
    const auto *mf_302 = buffer.data(mf + 302);
    const auto *mf_306 = buffer.data(mf + 306);
    const auto *mf_307 = buffer.data(mf + 307);
    const auto *mf_308 = buffer.data(mf + 308);
    const auto *mf_309 = buffer.data(mf + 309);
    const auto *mf_310 = buffer.data(mf + 310);
    const auto *mf_312 = buffer.data(mf + 312);
    const auto *mf_316 = buffer.data(mf + 316);
    const auto *mf_317 = buffer.data(mf + 317);
    const auto *mf_318 = buffer.data(mf + 318);
    const auto *mf_319 = buffer.data(mf + 319);
    const auto *mf_320 = buffer.data(mf + 320);
    const auto *mf_322 = buffer.data(mf + 322);
    const auto *mf_326 = buffer.data(mf + 326);
    const auto *mf_327 = buffer.data(mf + 327);
    const auto *mf_328 = buffer.data(mf + 328);
    const auto *mf_329 = buffer.data(mf + 329);
    const auto *mf_330 = buffer.data(mf + 330);
    const auto *mf_332 = buffer.data(mf + 332);
    const auto *mf_336 = buffer.data(mf + 336);
    const auto *mf_337 = buffer.data(mf + 337);
    const auto *mf_338 = buffer.data(mf + 338);
    const auto *mf_339 = buffer.data(mf + 339);
    const auto *mf_340 = buffer.data(mf + 340);
    const auto *mf_342 = buffer.data(mf + 342);
    const auto *mf_346 = buffer.data(mf + 346);
    const auto *mf_347 = buffer.data(mf + 347);
    const auto *mf_348 = buffer.data(mf + 348);
    const auto *mf_349 = buffer.data(mf + 349);
    const auto *mf_350 = buffer.data(mf + 350);
    const auto *mf_351 = buffer.data(mf + 351);
    const auto *mf_352 = buffer.data(mf + 352);
    const auto *mf_355 = buffer.data(mf + 355);
    const auto *mf_356 = buffer.data(mf + 356);
    const auto *mf_357 = buffer.data(mf + 357);

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_z, pb_y, pb_z, kg0_210, kg1_210, \
                         lf_200, lg_300, md0_162, md1_162, mf_270, \
                         mf_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_18 * kg0_210[k]
                   - f_19 * kg1_210[k]
                   + pa_z[k] * lg_300[k];

        t_406[k] = pb_y[k] * mf_270[k];

        t_407[k] = f_17 * lf_200[k]
                   + pb_z[k] * mf_270[k];

        t_408[k] = f_3 * md0_162[k]
                   - f_4 * md1_162[k]
                   + pb_y[k] * mf_271[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pb_x, pb_y, lf_275, lf_276, \
                         lf_277, md0_167, md1_167, mf_272, mf_275, mf_276, \
                         mf_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = pb_y[k] * mf_272[k];

        t_410[k] = f_16 * lf_275[k]
                   + f_3 * md0_167[k]
                   - f_4 * md1_167[k]
                   + pb_x[k] * mf_275[k];

        t_411[k] = f_16 * lf_276[k]
                   + pb_x[k] * mf_276[k];

        t_412[k] = f_16 * lf_277[k]
                   + pb_x[k] * mf_277[k];

        t_413[k] = pb_y[k] * mf_275[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pb_x, pb_y, pb_z, lf_206, lf_279, \
                         md0_165, md0_167, md1_165, md1_167, mf_276, mf_278, \
                         mf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_16 * lf_279[k]
                   + pb_x[k] * mf_279[k];

        t_415[k] = f_1 * md0_165[k]
                   - f_2 * md1_165[k]
                   + pb_y[k] * mf_276[k];

        t_416[k] = f_17 * lf_206[k]
                   + pb_z[k] * mf_276[k];

        t_417[k] = f_3 * md0_167[k]
                   - f_4 * md1_167[k]
                   + pb_y[k] * mf_278[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, pa_x, pa_y, pb_y, kg0_225, kg0_419, \
                         kg1_225, kg1_419, lf_210, lg_315, lg_419, mf_279, \
                         mf_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = pb_y[k] * mf_279[k];

        t_419[k] = f_14 * kg0_419[k]
                   - f_15 * kg1_419[k]
                   + pa_x[k] * lg_419[k];

        t_420[k] = f_12 * kg0_225[k]
                   - f_13 * kg1_225[k]
                   + pa_y[k] * lg_315[k];

        t_421[k] = f_11 * lf_210[k]
                   + pb_y[k] * mf_280[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, pb_x, pb_z, lf_283, md0_168, md0_171, \
                         md1_168, md1_171, mf_280, mf_281, mf_282, \
                         mf_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = pb_z[k] * mf_280[k];

        t_423[k] = f_6 * lf_283[k]
                   + f_3 * md0_171[k]
                   - f_4 * md1_171[k]
                   + pb_x[k] * mf_283[k];

        t_424[k] = pb_z[k] * mf_281[k];

        t_425[k] = f_3 * md0_168[k]
                   - f_4 * md1_168[k]
                   + pb_z[k] * mf_282[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, pb_x, pb_z, lf_286, lf_288, lf_289, \
                         mf_283, mf_286, mf_288, mf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_6 * lf_286[k]
                   + pb_x[k] * mf_286[k];

        t_427[k] = pb_z[k] * mf_283[k];

        t_428[k] = f_6 * lf_288[k]
                   + pb_x[k] * mf_288[k];

        t_429[k] = f_6 * lf_289[k]
                   + pb_x[k] * mf_289[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, pa_x, pb_y, pb_z, kg0_430, kg1_430, \
                         lf_219, lg_430, md0_171, md1_171, mf_286, mf_287, \
                         mf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_9 * kg0_430[k]
                   - f_10 * kg1_430[k]
                   + pa_x[k] * lg_430[k];

        t_431[k] = pb_z[k] * mf_286[k];

        t_432[k] = f_3 * md0_171[k]
                   - f_4 * md1_171[k]
                   + pb_z[k] * mf_287[k];

        t_433[k] = f_11 * lf_219[k]
                   + pb_y[k] * mf_289[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, t_438, pa_z, pb_z, lf_210, lg_315, \
                         lg_316, lg_318, md0_173, md1_173, mf_289, \
                         mf_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_1 * md0_173[k]
                   - f_2 * md1_173[k]
                   + pb_z[k] * mf_289[k];

        t_435[k] = pa_z[k] * lg_315[k];

        t_436[k] = pa_z[k] * lg_316[k];

        t_437[k] = f_5 * lf_210[k]
                   + pb_z[k] * mf_290[k];

        t_438[k] = pa_z[k] * lg_318[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, pa_z, pb_x, pb_y, lf_212, lf_222, lf_297, \
                         lg_320, lg_321, mf_292, mf_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_17 * lf_222[k]
                   + pb_y[k] * mf_292[k];

        t_440[k] = f_6 * lf_212[k]
                   + pa_z[k] * lg_320[k];

        t_441[k] = pa_z[k] * lg_321[k];

        t_442[k] = f_6 * lf_297[k]
                   + pb_x[k] * mf_297[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, pa_z, pb_x, pb_z, lf_216, lf_298, lf_299, \
                         lg_325, mf_296, mf_298, mf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_6 * lf_298[k]
                   + pb_x[k] * mf_298[k];

        t_444[k] = f_6 * lf_299[k]
                   + pb_x[k] * mf_299[k];

        t_445[k] = pa_z[k] * lg_325[k];

        t_446[k] = f_5 * lf_216[k]
                   + pb_z[k] * mf_296[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, pa_y, pa_z, pb_y, kg0_255, kg1_255, \
                         lf_217, lf_219, lf_229, lg_327, lg_329, lg_345, \
                         mf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_6 * lf_217[k]
                   + pa_z[k] * lg_327[k];

        t_448[k] = f_17 * lf_229[k]
                   + pb_y[k] * mf_299[k];

        t_449[k] = f_8 * lf_219[k]
                   + pa_z[k] * lg_329[k];

        t_450[k] = f_23 * kg0_255[k]
                   - f_24 * kg1_255[k]
                   + pa_y[k] * lg_345[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, pa_z, pb_y, pb_z, kg0_228, kg1_228, \
                         lf_220, lf_230, lf_232, lg_333, mf_300, \
                         mf_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_22 * lf_230[k]
                   + pb_y[k] * mf_300[k];

        t_452[k] = f_6 * lf_220[k]
                   + pb_z[k] * mf_300[k];

        t_453[k] = f_9 * kg0_228[k]
                   - f_10 * kg1_228[k]
                   + pa_z[k] * lg_333[k];

        t_454[k] = f_22 * lf_232[k]
                   + pb_y[k] * mf_302[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, pa_y, pb_x, kg0_260, kg1_260, lf_306, \
                         lf_307, lf_308, lg_350, mf_306, mf_307, \
                         mf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_23 * kg0_260[k]
                   - f_24 * kg1_260[k]
                   + pa_y[k] * lg_350[k];

        t_456[k] = f_6 * lf_306[k]
                   + pb_x[k] * mf_306[k];

        t_457[k] = f_6 * lf_307[k]
                   + pb_x[k] * mf_307[k];

        t_458[k] = f_6 * lf_308[k]
                   + pb_x[k] * mf_308[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pa_x, pb_x, pb_z, kg0_460, kg1_460, lf_226, \
                         lf_309, lg_460, mf_306, mf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_6 * lf_309[k]
                   + pb_x[k] * mf_309[k];

        t_460[k] = f_9 * kg0_460[k]
                   - f_10 * kg1_460[k]
                   + pa_x[k] * lg_460[k];

        t_461[k] = f_6 * lf_226[k]
                   + pb_z[k] * mf_306[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pa_x, pb_y, kg0_462, kg0_464, kg1_462, kg1_464, \
                         lf_239, lg_462, lg_464, mf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_9 * kg0_462[k]
                   - f_10 * kg1_462[k]
                   + pa_x[k] * lg_462[k];

        t_463[k] = f_22 * lf_239[k]
                   + pb_y[k] * mf_309[k];

        t_464[k] = f_9 * kg0_464[k]
                   - f_10 * kg1_464[k]
                   + pa_x[k] * lg_464[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pa_y, pb_y, pb_z, kg0_270, kg1_270, lf_230, \
                         lf_240, lg_360, mf_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_20 * kg0_270[k]
                   - f_21 * kg1_270[k]
                   + pa_y[k] * lg_360[k];

        t_466[k] = f_8 * lf_240[k]
                   + pb_y[k] * mf_310[k];

        t_467[k] = f_16 * lf_230[k]
                   + pb_z[k] * mf_310[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pa_y, pa_z, pb_y, kg0_243, kg0_275, kg1_243, \
                         kg1_275, lf_242, lg_348, lg_365, mf_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_14 * kg0_243[k]
                   - f_15 * kg1_243[k]
                   + pa_z[k] * lg_348[k];

        t_469[k] = f_8 * lf_242[k]
                   + pb_y[k] * mf_312[k];

        t_470[k] = f_20 * kg0_275[k]
                   - f_21 * kg1_275[k]
                   + pa_y[k] * lg_365[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pb_x, lf_316, lf_317, lf_318, lf_319, \
                         mf_316, mf_317, mf_318, mf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_6 * lf_316[k]
                   + pb_x[k] * mf_316[k];

        t_472[k] = f_6 * lf_317[k]
                   + pb_x[k] * mf_317[k];

        t_473[k] = f_6 * lf_318[k]
                   + pb_x[k] * mf_318[k];

        t_474[k] = f_6 * lf_319[k]
                   + pb_x[k] * mf_319[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, pa_x, pb_z, kg0_475, kg0_477, kg1_475, kg1_477, \
                         lf_236, lg_475, lg_477, mf_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_9 * kg0_475[k]
                   - f_10 * kg1_475[k]
                   + pa_x[k] * lg_475[k];

        t_476[k] = f_16 * lf_236[k]
                   + pb_z[k] * mf_316[k];

        t_477[k] = f_9 * kg0_477[k]
                   - f_10 * kg1_477[k]
                   + pa_x[k] * lg_477[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pa_x, pa_y, pb_y, kg0_285, kg0_479, kg1_285, \
                         kg1_479, lf_249, lg_375, lg_479, mf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_8 * lf_249[k]
                   + pb_y[k] * mf_319[k];

        t_479[k] = f_9 * kg0_479[k]
                   - f_10 * kg1_479[k]
                   + pa_x[k] * lg_479[k];

        t_480[k] = f_14 * kg0_285[k]
                   - f_15 * kg1_285[k]
                   + pa_y[k] * lg_375[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pa_z, pb_y, pb_z, kg0_258, kg1_258, \
                         lf_240, lf_250, lf_252, lg_363, mf_320, \
                         mf_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_16 * lf_250[k]
                   + pb_y[k] * mf_320[k];

        t_482[k] = f_8 * lf_240[k]
                   + pb_z[k] * mf_320[k];

        t_483[k] = f_20 * kg0_258[k]
                   - f_21 * kg1_258[k]
                   + pa_z[k] * lg_363[k];

        t_484[k] = f_16 * lf_252[k]
                   + pb_y[k] * mf_322[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pa_y, pb_x, kg0_290, kg1_290, lf_326, \
                         lf_327, lf_328, lg_380, mf_326, mf_327, \
                         mf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_14 * kg0_290[k]
                   - f_15 * kg1_290[k]
                   + pa_y[k] * lg_380[k];

        t_486[k] = f_6 * lf_326[k]
                   + pb_x[k] * mf_326[k];

        t_487[k] = f_6 * lf_327[k]
                   + pb_x[k] * mf_327[k];

        t_488[k] = f_6 * lf_328[k]
                   + pb_x[k] * mf_328[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pa_x, pb_x, pb_z, kg0_490, kg1_490, lf_246, \
                         lf_329, lg_490, mf_326, mf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_6 * lf_329[k]
                   + pb_x[k] * mf_329[k];

        t_490[k] = f_9 * kg0_490[k]
                   - f_10 * kg1_490[k]
                   + pa_x[k] * lg_490[k];

        t_491[k] = f_8 * lf_246[k]
                   + pb_z[k] * mf_326[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, pa_x, pb_y, kg0_492, kg0_494, kg1_492, kg1_494, \
                         lf_259, lg_492, lg_494, mf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_9 * kg0_492[k]
                   - f_10 * kg1_492[k]
                   + pa_x[k] * lg_492[k];

        t_493[k] = f_16 * lf_259[k]
                   + pb_y[k] * mf_329[k];

        t_494[k] = f_9 * kg0_494[k]
                   - f_10 * kg1_494[k]
                   + pa_x[k] * lg_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, pa_y, pb_y, pb_z, kg0_300, kg1_300, lf_250, \
                         lf_260, lg_390, mf_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_9 * kg0_300[k]
                   - f_10 * kg1_300[k]
                   + pa_y[k] * lg_390[k];

        t_496[k] = f_6 * lf_260[k]
                   + pb_y[k] * mf_330[k];

        t_497[k] = f_22 * lf_250[k]
                   + pb_z[k] * mf_330[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pa_y, pa_z, pb_y, kg0_273, kg0_305, kg1_273, \
                         kg1_305, lf_262, lg_378, lg_395, mf_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_23 * kg0_273[k]
                   - f_24 * kg1_273[k]
                   + pa_z[k] * lg_378[k];

        t_499[k] = f_6 * lf_262[k]
                   + pb_y[k] * mf_332[k];

        t_500[k] = f_9 * kg0_305[k]
                   - f_10 * kg1_305[k]
                   + pa_y[k] * lg_395[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, pb_x, lf_336, lf_337, lf_338, lf_339, \
                         mf_336, mf_337, mf_338, mf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_6 * lf_336[k]
                   + pb_x[k] * mf_336[k];

        t_502[k] = f_6 * lf_337[k]
                   + pb_x[k] * mf_337[k];

        t_503[k] = f_6 * lf_338[k]
                   + pb_x[k] * mf_338[k];

        t_504[k] = f_6 * lf_339[k]
                   + pb_x[k] * mf_339[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, pa_x, pb_z, kg0_505, kg0_507, kg1_505, kg1_507, \
                         lf_256, lg_505, lg_507, mf_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_9 * kg0_505[k]
                   - f_10 * kg1_505[k]
                   + pa_x[k] * lg_505[k];

        t_506[k] = f_22 * lf_256[k]
                   + pb_z[k] * mf_336[k];

        t_507[k] = f_9 * kg0_507[k]
                   - f_10 * kg1_507[k]
                   + pa_x[k] * lg_507[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pa_x, pa_y, pb_y, kg0_509, kg1_509, \
                         lf_269, lf_270, lg_405, lg_509, mf_339, \
                         mf_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_6 * lf_269[k]
                   + pb_y[k] * mf_339[k];

        t_509[k] = f_9 * kg0_509[k]
                   - f_10 * kg1_509[k]
                   + pa_x[k] * lg_509[k];

        t_510[k] = pa_y[k] * lg_405[k];

        t_511[k] = f_5 * lf_270[k]
                   + pb_y[k] * mf_340[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, pa_y, pb_x, pb_y, lf_271, lf_272, \
                         lf_346, lg_407, lg_408, lg_410, mf_342, \
                         mf_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pa_y[k] * lg_407[k];

        t_513[k] = f_6 * lf_271[k]
                   + pa_y[k] * lg_408[k];

        t_514[k] = f_5 * lf_272[k]
                   + pb_y[k] * mf_342[k];

        t_515[k] = pa_y[k] * lg_410[k];

        t_516[k] = f_6 * lf_346[k]
                   + pb_x[k] * mf_346[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, pa_y, pb_x, lf_276, lf_347, lf_348, \
                         lg_414, lg_415, mf_347, mf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_6 * lf_347[k]
                   + pb_x[k] * mf_347[k];

        t_518[k] = f_6 * lf_348[k]
                   + pb_x[k] * mf_348[k];

        t_519[k] = pa_y[k] * lg_414[k];

        t_520[k] = f_8 * lf_276[k]
                   + pa_y[k] * lg_415[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pa_y, pb_y, pb_z, lf_266, lf_278, lf_279, \
                         lg_417, lg_419, mf_346, mf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_17 * lf_266[k]
                   + pb_z[k] * mf_346[k];

        t_522[k] = f_6 * lf_278[k]
                   + pa_y[k] * lg_417[k];

        t_523[k] = f_5 * lf_279[k]
                   + pb_y[k] * mf_349[k];

        t_524[k] = pa_y[k] * lg_419[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, pa_z, pb_y, pb_z, kg0_300, kg1_300, \
                         lf_270, lg_405, md0_210, md1_210, mf_350, \
                         mf_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_12 * kg0_300[k]
                   - f_13 * kg1_300[k]
                   + pa_z[k] * lg_405[k];

        t_526[k] = pb_y[k] * mf_350[k];

        t_527[k] = f_11 * lf_270[k]
                   + pb_z[k] * mf_350[k];

        t_528[k] = f_3 * md0_210[k]
                   - f_4 * md1_210[k]
                   + pb_y[k] * mf_351[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, t_533, pb_x, pb_y, lf_355, lf_356, \
                         lf_357, md0_215, md1_215, mf_352, mf_355, mf_356, \
                         mf_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = pb_y[k] * mf_352[k];

        t_530[k] = f_6 * lf_355[k]
                   + f_3 * md0_215[k]
                   - f_4 * md1_215[k]
                   + pb_x[k] * mf_355[k];

        t_531[k] = f_6 * lf_356[k]
                   + pb_x[k] * mf_356[k];

        t_532[k] = f_6 * lf_357[k]
                   + pb_x[k] * mf_357[k];

        t_533[k] = pb_y[k] * mf_355[k];
    }
}

static auto
compute_prim_mg_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kg0,
                                            const size_t kg1, const size_t lf, const size_t lg,
                                            const size_t md0, const size_t md1, const size_t mf,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 4.0 / p;
    const auto f_8 = 2.0 / p;
    const auto f_9 = 0.5 / alpha;
    const auto f_10 = 0.5 * beta / (alpha * p);
    const auto f_11 = 3.5 / p;
    const auto f_16 = 1.5 / p;
    const auto f_17 = 3.0 / p;
    const auto f_22 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kg0_539 = buffer.data(kg0 + 539);

    const auto *kg1_539 = buffer.data(kg1 + 539);

    const auto *lf_276 = buffer.data(lf + 276);
    const auto *lf_280 = buffer.data(lf + 280);
    const auto *lf_290 = buffer.data(lf + 290);
    const auto *lf_292 = buffer.data(lf + 292);
    const auto *lf_300 = buffer.data(lf + 300);
    const auto *lf_302 = buffer.data(lf + 302);
    const auto *lf_310 = buffer.data(lf + 310);
    const auto *lf_312 = buffer.data(lf + 312);
    const auto *lf_320 = buffer.data(lf + 320);
    const auto *lf_322 = buffer.data(lf + 322);
    const auto *lf_330 = buffer.data(lf + 330);
    const auto *lf_332 = buffer.data(lf + 332);
    const auto *lf_340 = buffer.data(lf + 340);
    const auto *lf_342 = buffer.data(lf + 342);
    const auto *lf_350 = buffer.data(lf + 350);
    const auto *lf_352 = buffer.data(lf + 352);
    const auto *lf_359 = buffer.data(lf + 359);
    const auto *lf_360 = buffer.data(lf + 360);
    const auto *lf_363 = buffer.data(lf + 363);
    const auto *lf_365 = buffer.data(lf + 365);
    const auto *lf_366 = buffer.data(lf + 366);
    const auto *lf_368 = buffer.data(lf + 368);
    const auto *lf_369 = buffer.data(lf + 369);
    const auto *lf_375 = buffer.data(lf + 375);
    const auto *lf_377 = buffer.data(lf + 377);
    const auto *lf_378 = buffer.data(lf + 378);
    const auto *lf_379 = buffer.data(lf + 379);
    const auto *lf_380 = buffer.data(lf + 380);
    const auto *lf_383 = buffer.data(lf + 383);
    const auto *lf_385 = buffer.data(lf + 385);
    const auto *lf_386 = buffer.data(lf + 386);
    const auto *lf_387 = buffer.data(lf + 387);
    const auto *lf_388 = buffer.data(lf + 388);
    const auto *lf_389 = buffer.data(lf + 389);
    const auto *lf_390 = buffer.data(lf + 390);
    const auto *lf_393 = buffer.data(lf + 393);
    const auto *lf_395 = buffer.data(lf + 395);
    const auto *lf_396 = buffer.data(lf + 396);
    const auto *lf_397 = buffer.data(lf + 397);
    const auto *lf_398 = buffer.data(lf + 398);
    const auto *lf_399 = buffer.data(lf + 399);
    const auto *lf_400 = buffer.data(lf + 400);
    const auto *lf_403 = buffer.data(lf + 403);
    const auto *lf_405 = buffer.data(lf + 405);
    const auto *lf_406 = buffer.data(lf + 406);
    const auto *lf_407 = buffer.data(lf + 407);
    const auto *lf_408 = buffer.data(lf + 408);
    const auto *lf_409 = buffer.data(lf + 409);
    const auto *lf_410 = buffer.data(lf + 410);
    const auto *lf_413 = buffer.data(lf + 413);
    const auto *lf_415 = buffer.data(lf + 415);
    const auto *lf_416 = buffer.data(lf + 416);
    const auto *lf_417 = buffer.data(lf + 417);
    const auto *lf_418 = buffer.data(lf + 418);
    const auto *lf_419 = buffer.data(lf + 419);
    const auto *lf_420 = buffer.data(lf + 420);
    const auto *lf_423 = buffer.data(lf + 423);
    const auto *lf_425 = buffer.data(lf + 425);
    const auto *lf_426 = buffer.data(lf + 426);
    const auto *lf_427 = buffer.data(lf + 427);
    const auto *lf_428 = buffer.data(lf + 428);
    const auto *lf_429 = buffer.data(lf + 429);
    const auto *lf_433 = buffer.data(lf + 433);
    const auto *lf_436 = buffer.data(lf + 436);
    const auto *lf_437 = buffer.data(lf + 437);
    const auto *lf_438 = buffer.data(lf + 438);
    const auto *lf_440 = buffer.data(lf + 440);
    const auto *lf_443 = buffer.data(lf + 443);
    const auto *lf_445 = buffer.data(lf + 445);
    const auto *lf_446 = buffer.data(lf + 446);
    const auto *lf_447 = buffer.data(lf + 447);
    const auto *lf_449 = buffer.data(lf + 449);

    const auto *lg_420 = buffer.data(lg + 420);
    const auto *lg_421 = buffer.data(lg + 421);
    const auto *lg_423 = buffer.data(lg + 423);
    const auto *lg_426 = buffer.data(lg + 426);
    const auto *lg_525 = buffer.data(lg + 525);
    const auto *lg_527 = buffer.data(lg + 527);
    const auto *lg_530 = buffer.data(lg + 530);
    const auto *lg_534 = buffer.data(lg + 534);
    const auto *lg_539 = buffer.data(lg + 539);
    const auto *lg_540 = buffer.data(lg + 540);
    const auto *lg_543 = buffer.data(lg + 543);
    const auto *lg_545 = buffer.data(lg + 545);
    const auto *lg_550 = buffer.data(lg + 550);
    const auto *lg_552 = buffer.data(lg + 552);
    const auto *lg_553 = buffer.data(lg + 553);
    const auto *lg_554 = buffer.data(lg + 554);
    const auto *lg_560 = buffer.data(lg + 560);
    const auto *lg_565 = buffer.data(lg + 565);
    const auto *lg_566 = buffer.data(lg + 566);
    const auto *lg_567 = buffer.data(lg + 567);
    const auto *lg_568 = buffer.data(lg + 568);
    const auto *lg_569 = buffer.data(lg + 569);
    const auto *lg_570 = buffer.data(lg + 570);
    const auto *lg_573 = buffer.data(lg + 573);
    const auto *lg_575 = buffer.data(lg + 575);
    const auto *lg_580 = buffer.data(lg + 580);
    const auto *lg_581 = buffer.data(lg + 581);
    const auto *lg_582 = buffer.data(lg + 582);
    const auto *lg_583 = buffer.data(lg + 583);
    const auto *lg_584 = buffer.data(lg + 584);
    const auto *lg_585 = buffer.data(lg + 585);
    const auto *lg_588 = buffer.data(lg + 588);
    const auto *lg_590 = buffer.data(lg + 590);
    const auto *lg_595 = buffer.data(lg + 595);
    const auto *lg_596 = buffer.data(lg + 596);
    const auto *lg_597 = buffer.data(lg + 597);
    const auto *lg_598 = buffer.data(lg + 598);
    const auto *lg_599 = buffer.data(lg + 599);
    const auto *lg_600 = buffer.data(lg + 600);
    const auto *lg_603 = buffer.data(lg + 603);
    const auto *lg_605 = buffer.data(lg + 605);
    const auto *lg_610 = buffer.data(lg + 610);
    const auto *lg_611 = buffer.data(lg + 611);
    const auto *lg_612 = buffer.data(lg + 612);
    const auto *lg_613 = buffer.data(lg + 613);
    const auto *lg_614 = buffer.data(lg + 614);
    const auto *lg_615 = buffer.data(lg + 615);
    const auto *lg_618 = buffer.data(lg + 618);
    const auto *lg_620 = buffer.data(lg + 620);
    const auto *lg_625 = buffer.data(lg + 625);
    const auto *lg_626 = buffer.data(lg + 626);
    const auto *lg_627 = buffer.data(lg + 627);
    const auto *lg_628 = buffer.data(lg + 628);
    const auto *lg_629 = buffer.data(lg + 629);
    const auto *lg_630 = buffer.data(lg + 630);
    const auto *lg_633 = buffer.data(lg + 633);
    const auto *lg_635 = buffer.data(lg + 635);
    const auto *lg_640 = buffer.data(lg + 640);
    const auto *lg_641 = buffer.data(lg + 641);
    const auto *lg_642 = buffer.data(lg + 642);
    const auto *lg_643 = buffer.data(lg + 643);
    const auto *lg_644 = buffer.data(lg + 644);
    const auto *lg_648 = buffer.data(lg + 648);
    const auto *lg_655 = buffer.data(lg + 655);
    const auto *lg_656 = buffer.data(lg + 656);
    const auto *lg_657 = buffer.data(lg + 657);
    const auto *lg_658 = buffer.data(lg + 658);
    const auto *lg_659 = buffer.data(lg + 659);
    const auto *lg_660 = buffer.data(lg + 660);
    const auto *lg_663 = buffer.data(lg + 663);
    const auto *lg_665 = buffer.data(lg + 665);
    const auto *lg_670 = buffer.data(lg + 670);
    const auto *lg_671 = buffer.data(lg + 671);
    const auto *lg_672 = buffer.data(lg + 672);
    const auto *lg_674 = buffer.data(lg + 674);

    const auto *md0_213 = buffer.data(md0 + 213);
    const auto *md0_215 = buffer.data(md0 + 215);
    const auto *md0_270 = buffer.data(md0 + 270);
    const auto *md0_273 = buffer.data(md0 + 273);
    const auto *md0_275 = buffer.data(md0 + 275);

    const auto *md1_213 = buffer.data(md1 + 213);
    const auto *md1_215 = buffer.data(md1 + 215);
    const auto *md1_270 = buffer.data(md1 + 270);
    const auto *md1_273 = buffer.data(md1 + 273);
    const auto *md1_275 = buffer.data(md1 + 275);

    const auto *mf_356 = buffer.data(mf + 356);
    const auto *mf_358 = buffer.data(mf + 358);
    const auto *mf_359 = buffer.data(mf + 359);
    const auto *mf_360 = buffer.data(mf + 360);
    const auto *mf_361 = buffer.data(mf + 361);
    const auto *mf_363 = buffer.data(mf + 363);
    const auto *mf_366 = buffer.data(mf + 366);
    const auto *mf_368 = buffer.data(mf + 368);
    const auto *mf_369 = buffer.data(mf + 369);
    const auto *mf_370 = buffer.data(mf + 370);
    const auto *mf_372 = buffer.data(mf + 372);
    const auto *mf_377 = buffer.data(mf + 377);
    const auto *mf_378 = buffer.data(mf + 378);
    const auto *mf_379 = buffer.data(mf + 379);
    const auto *mf_380 = buffer.data(mf + 380);
    const auto *mf_382 = buffer.data(mf + 382);
    const auto *mf_386 = buffer.data(mf + 386);
    const auto *mf_387 = buffer.data(mf + 387);
    const auto *mf_388 = buffer.data(mf + 388);
    const auto *mf_389 = buffer.data(mf + 389);
    const auto *mf_390 = buffer.data(mf + 390);
    const auto *mf_392 = buffer.data(mf + 392);
    const auto *mf_396 = buffer.data(mf + 396);
    const auto *mf_397 = buffer.data(mf + 397);
    const auto *mf_398 = buffer.data(mf + 398);
    const auto *mf_399 = buffer.data(mf + 399);
    const auto *mf_400 = buffer.data(mf + 400);
    const auto *mf_402 = buffer.data(mf + 402);
    const auto *mf_406 = buffer.data(mf + 406);
    const auto *mf_407 = buffer.data(mf + 407);
    const auto *mf_408 = buffer.data(mf + 408);
    const auto *mf_409 = buffer.data(mf + 409);
    const auto *mf_410 = buffer.data(mf + 410);
    const auto *mf_412 = buffer.data(mf + 412);
    const auto *mf_416 = buffer.data(mf + 416);
    const auto *mf_417 = buffer.data(mf + 417);
    const auto *mf_418 = buffer.data(mf + 418);
    const auto *mf_419 = buffer.data(mf + 419);
    const auto *mf_420 = buffer.data(mf + 420);
    const auto *mf_422 = buffer.data(mf + 422);
    const auto *mf_426 = buffer.data(mf + 426);
    const auto *mf_427 = buffer.data(mf + 427);
    const auto *mf_428 = buffer.data(mf + 428);
    const auto *mf_429 = buffer.data(mf + 429);
    const auto *mf_430 = buffer.data(mf + 430);
    const auto *mf_432 = buffer.data(mf + 432);
    const auto *mf_436 = buffer.data(mf + 436);
    const auto *mf_437 = buffer.data(mf + 437);
    const auto *mf_438 = buffer.data(mf + 438);
    const auto *mf_440 = buffer.data(mf + 440);
    const auto *mf_442 = buffer.data(mf + 442);
    const auto *mf_445 = buffer.data(mf + 445);
    const auto *mf_446 = buffer.data(mf + 446);
    const auto *mf_447 = buffer.data(mf + 447);
    const auto *mf_449 = buffer.data(mf + 449);
    const auto *mf_450 = buffer.data(mf + 450);
    const auto *mf_451 = buffer.data(mf + 451);
    const auto *mf_453 = buffer.data(mf + 453);
    const auto *mf_455 = buffer.data(mf + 455);
    const auto *mf_456 = buffer.data(mf + 456);
    const auto *mf_457 = buffer.data(mf + 457);
    const auto *mf_458 = buffer.data(mf + 458);
    const auto *mf_459 = buffer.data(mf + 459);

#pragma omp simd aligned(t_534, t_535, t_536, t_537, pb_x, pb_y, pb_z, lf_276, lf_359, \
                         md0_213, md0_215, md1_213, md1_215, mf_356, mf_358, \
                         mf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_6 * lf_359[k]
                   + pb_x[k] * mf_359[k];

        t_535[k] = f_1 * md0_213[k]
                   - f_2 * md1_213[k]
                   + pb_y[k] * mf_356[k];

        t_536[k] = f_11 * lf_276[k]
                   + pb_z[k] * mf_356[k];

        t_537[k] = f_3 * md0_215[k]
                   - f_4 * md1_215[k]
                   + pb_y[k] * mf_358[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, t_542, pa_x, pb_y, pb_z, kg0_539, \
                         kg1_539, lf_280, lf_360, lg_539, lg_540, mf_359, \
                         mf_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = pb_y[k] * mf_359[k];

        t_539[k] = f_9 * kg0_539[k]
                   - f_10 * kg1_539[k]
                   + pa_x[k] * lg_539[k];

        t_540[k] = f_8 * lf_360[k]
                   + pa_x[k] * lg_540[k];

        t_541[k] = f_7 * lf_280[k]
                   + pb_y[k] * mf_360[k];

        t_542[k] = pb_z[k] * mf_360[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, t_547, pa_x, pb_x, pb_z, lf_363, lf_365, \
                         lf_366, lg_543, lg_545, mf_361, mf_363, \
                         mf_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_6 * lf_363[k]
                   + pa_x[k] * lg_543[k];

        t_544[k] = pb_z[k] * mf_361[k];

        t_545[k] = f_6 * lf_365[k]
                   + pa_x[k] * lg_545[k];

        t_546[k] = f_5 * lf_366[k]
                   + pb_x[k] * mf_366[k];

        t_547[k] = pb_z[k] * mf_363[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pa_x, pb_x, pb_z, lf_368, lf_369, \
                         lg_550, lg_552, mf_366, mf_368, mf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_5 * lf_368[k]
                   + pb_x[k] * mf_368[k];

        t_549[k] = f_5 * lf_369[k]
                   + pb_x[k] * mf_369[k];

        t_550[k] = pa_x[k] * lg_550[k];

        t_551[k] = pb_z[k] * mf_366[k];

        t_552[k] = pa_x[k] * lg_552[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, t_557, t_558, pa_x, pa_z, pb_z, lf_280, \
                         lg_420, lg_421, lg_423, lg_553, lg_554, \
                         mf_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = pa_x[k] * lg_553[k];

        t_554[k] = pa_x[k] * lg_554[k];

        t_555[k] = pa_z[k] * lg_420[k];

        t_556[k] = pa_z[k] * lg_421[k];

        t_557[k] = f_5 * lf_280[k]
                   + pb_z[k] * mf_370[k];

        t_558[k] = pa_z[k] * lg_423[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pa_x, pa_z, pb_x, pb_y, lf_292, lf_375, \
                         lf_377, lg_426, lg_560, mf_372, mf_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_11 * lf_292[k]
                   + pb_y[k] * mf_372[k];

        t_560[k] = f_6 * lf_375[k]
                   + pa_x[k] * lg_560[k];

        t_561[k] = pa_z[k] * lg_426[k];

        t_562[k] = f_5 * lf_377[k]
                   + pb_x[k] * mf_377[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, t_567, t_568, pa_x, pb_x, lf_378, lf_379, \
                         lg_565, lg_566, lg_567, lg_568, mf_378, \
                         mf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_5 * lf_378[k]
                   + pb_x[k] * mf_378[k];

        t_564[k] = f_5 * lf_379[k]
                   + pb_x[k] * mf_379[k];

        t_565[k] = pa_x[k] * lg_565[k];

        t_566[k] = pa_x[k] * lg_566[k];

        t_567[k] = pa_x[k] * lg_567[k];

        t_568[k] = pa_x[k] * lg_568[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, t_572, t_573, pa_x, pb_y, pb_z, lf_290, lf_300, \
                         lf_380, lf_383, lg_569, lg_570, lg_573, \
                         mf_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = pa_x[k] * lg_569[k];

        t_570[k] = f_8 * lf_380[k]
                   + pa_x[k] * lg_570[k];

        t_571[k] = f_17 * lf_300[k]
                   + pb_y[k] * mf_380[k];

        t_572[k] = f_6 * lf_290[k]
                   + pb_z[k] * mf_380[k];

        t_573[k] = f_6 * lf_383[k]
                   + pa_x[k] * lg_573[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pa_x, pb_x, pb_y, lf_302, lf_385, lf_386, \
                         lf_387, lg_575, mf_382, mf_386, mf_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_17 * lf_302[k]
                   + pb_y[k] * mf_382[k];

        t_575[k] = f_6 * lf_385[k]
                   + pa_x[k] * lg_575[k];

        t_576[k] = f_5 * lf_386[k]
                   + pb_x[k] * mf_386[k];

        t_577[k] = f_5 * lf_387[k]
                   + pb_x[k] * mf_387[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, t_582, t_583, pa_x, pb_x, lf_388, lf_389, \
                         lg_580, lg_581, lg_582, lg_583, mf_388, \
                         mf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_5 * lf_388[k]
                   + pb_x[k] * mf_388[k];

        t_579[k] = f_5 * lf_389[k]
                   + pb_x[k] * mf_389[k];

        t_580[k] = pa_x[k] * lg_580[k];

        t_581[k] = pa_x[k] * lg_581[k];

        t_582[k] = pa_x[k] * lg_582[k];

        t_583[k] = pa_x[k] * lg_583[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, t_588, pa_x, pb_y, pb_z, lf_300, lf_310, \
                         lf_390, lf_393, lg_584, lg_585, lg_588, \
                         mf_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = pa_x[k] * lg_584[k];

        t_585[k] = f_8 * lf_390[k]
                   + pa_x[k] * lg_585[k];

        t_586[k] = f_22 * lf_310[k]
                   + pb_y[k] * mf_390[k];

        t_587[k] = f_16 * lf_300[k]
                   + pb_z[k] * mf_390[k];

        t_588[k] = f_6 * lf_393[k]
                   + pa_x[k] * lg_588[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, pa_x, pb_x, pb_y, lf_312, lf_395, lf_396, \
                         lf_397, lg_590, mf_392, mf_396, mf_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = f_22 * lf_312[k]
                   + pb_y[k] * mf_392[k];

        t_590[k] = f_6 * lf_395[k]
                   + pa_x[k] * lg_590[k];

        t_591[k] = f_5 * lf_396[k]
                   + pb_x[k] * mf_396[k];

        t_592[k] = f_5 * lf_397[k]
                   + pb_x[k] * mf_397[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, t_596, t_597, t_598, pa_x, pb_x, lf_398, lf_399, \
                         lg_595, lg_596, lg_597, lg_598, mf_398, \
                         mf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = f_5 * lf_398[k]
                   + pb_x[k] * mf_398[k];

        t_594[k] = f_5 * lf_399[k]
                   + pb_x[k] * mf_399[k];

        t_595[k] = pa_x[k] * lg_595[k];

        t_596[k] = pa_x[k] * lg_596[k];

        t_597[k] = pa_x[k] * lg_597[k];

        t_598[k] = pa_x[k] * lg_598[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, t_603, pa_x, pb_y, pb_z, lf_310, lf_320, \
                         lf_400, lf_403, lg_599, lg_600, lg_603, \
                         mf_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = pa_x[k] * lg_599[k];

        t_600[k] = f_8 * lf_400[k]
                   + pa_x[k] * lg_600[k];

        t_601[k] = f_8 * lf_320[k]
                   + pb_y[k] * mf_400[k];

        t_602[k] = f_8 * lf_310[k]
                   + pb_z[k] * mf_400[k];

        t_603[k] = f_6 * lf_403[k]
                   + pa_x[k] * lg_603[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pa_x, pb_x, pb_y, lf_322, lf_405, lf_406, \
                         lf_407, lg_605, mf_402, mf_406, mf_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_8 * lf_322[k]
                   + pb_y[k] * mf_402[k];

        t_605[k] = f_6 * lf_405[k]
                   + pa_x[k] * lg_605[k];

        t_606[k] = f_5 * lf_406[k]
                   + pb_x[k] * mf_406[k];

        t_607[k] = f_5 * lf_407[k]
                   + pb_x[k] * mf_407[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, t_612, t_613, pa_x, pb_x, lf_408, lf_409, \
                         lg_610, lg_611, lg_612, lg_613, mf_408, \
                         mf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_5 * lf_408[k]
                   + pb_x[k] * mf_408[k];

        t_609[k] = f_5 * lf_409[k]
                   + pb_x[k] * mf_409[k];

        t_610[k] = pa_x[k] * lg_610[k];

        t_611[k] = pa_x[k] * lg_611[k];

        t_612[k] = pa_x[k] * lg_612[k];

        t_613[k] = pa_x[k] * lg_613[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, pa_x, pb_y, pb_z, lf_320, lf_330, \
                         lf_410, lf_413, lg_614, lg_615, lg_618, \
                         mf_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = pa_x[k] * lg_614[k];

        t_615[k] = f_8 * lf_410[k]
                   + pa_x[k] * lg_615[k];

        t_616[k] = f_16 * lf_330[k]
                   + pb_y[k] * mf_410[k];

        t_617[k] = f_22 * lf_320[k]
                   + pb_z[k] * mf_410[k];

        t_618[k] = f_6 * lf_413[k]
                   + pa_x[k] * lg_618[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, pa_x, pb_x, pb_y, lf_332, lf_415, lf_416, \
                         lf_417, lg_620, mf_412, mf_416, mf_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_16 * lf_332[k]
                   + pb_y[k] * mf_412[k];

        t_620[k] = f_6 * lf_415[k]
                   + pa_x[k] * lg_620[k];

        t_621[k] = f_5 * lf_416[k]
                   + pb_x[k] * mf_416[k];

        t_622[k] = f_5 * lf_417[k]
                   + pb_x[k] * mf_417[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, t_626, t_627, t_628, pa_x, pb_x, lf_418, lf_419, \
                         lg_625, lg_626, lg_627, lg_628, mf_418, \
                         mf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_5 * lf_418[k]
                   + pb_x[k] * mf_418[k];

        t_624[k] = f_5 * lf_419[k]
                   + pb_x[k] * mf_419[k];

        t_625[k] = pa_x[k] * lg_625[k];

        t_626[k] = pa_x[k] * lg_626[k];

        t_627[k] = pa_x[k] * lg_627[k];

        t_628[k] = pa_x[k] * lg_628[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, t_633, pa_x, pb_y, pb_z, lf_330, lf_340, \
                         lf_420, lf_423, lg_629, lg_630, lg_633, \
                         mf_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = pa_x[k] * lg_629[k];

        t_630[k] = f_8 * lf_420[k]
                   + pa_x[k] * lg_630[k];

        t_631[k] = f_6 * lf_340[k]
                   + pb_y[k] * mf_420[k];

        t_632[k] = f_17 * lf_330[k]
                   + pb_z[k] * mf_420[k];

        t_633[k] = f_6 * lf_423[k]
                   + pa_x[k] * lg_633[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pa_x, pb_x, pb_y, lf_342, lf_425, lf_426, \
                         lf_427, lg_635, mf_422, mf_426, mf_427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_6 * lf_342[k]
                   + pb_y[k] * mf_422[k];

        t_635[k] = f_6 * lf_425[k]
                   + pa_x[k] * lg_635[k];

        t_636[k] = f_5 * lf_426[k]
                   + pb_x[k] * mf_426[k];

        t_637[k] = f_5 * lf_427[k]
                   + pb_x[k] * mf_427[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, t_642, t_643, pa_x, pb_x, lf_428, lf_429, \
                         lg_640, lg_641, lg_642, lg_643, mf_428, \
                         mf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_5 * lf_428[k]
                   + pb_x[k] * mf_428[k];

        t_639[k] = f_5 * lf_429[k]
                   + pb_x[k] * mf_429[k];

        t_640[k] = pa_x[k] * lg_640[k];

        t_641[k] = pa_x[k] * lg_641[k];

        t_642[k] = pa_x[k] * lg_642[k];

        t_643[k] = pa_x[k] * lg_643[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, t_647, t_648, pa_x, pa_y, pb_y, lf_350, lf_433, \
                         lg_525, lg_527, lg_644, lg_648, mf_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = pa_x[k] * lg_644[k];

        t_645[k] = pa_y[k] * lg_525[k];

        t_646[k] = f_5 * lf_350[k]
                   + pb_y[k] * mf_430[k];

        t_647[k] = pa_y[k] * lg_527[k];

        t_648[k] = f_6 * lf_433[k]
                   + pa_x[k] * lg_648[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pa_y, pb_x, pb_y, lf_352, lf_436, lf_437, \
                         lg_530, mf_432, mf_436, mf_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_5 * lf_352[k]
                   + pb_y[k] * mf_432[k];

        t_650[k] = pa_y[k] * lg_530[k];

        t_651[k] = f_5 * lf_436[k]
                   + pb_x[k] * mf_436[k];

        t_652[k] = f_5 * lf_437[k]
                   + pb_x[k] * mf_437[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, t_656, t_657, t_658, pa_x, pa_y, pb_x, lf_438, \
                         lg_534, lg_655, lg_656, lg_657, lg_658, \
                         mf_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_5 * lf_438[k]
                   + pb_x[k] * mf_438[k];

        t_654[k] = pa_y[k] * lg_534[k];

        t_655[k] = pa_x[k] * lg_655[k];

        t_656[k] = pa_x[k] * lg_656[k];

        t_657[k] = pa_x[k] * lg_657[k];

        t_658[k] = pa_x[k] * lg_658[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, pa_x, pb_y, pb_z, lf_350, lf_440, \
                         lf_443, lg_659, lg_660, lg_663, mf_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = pa_x[k] * lg_659[k];

        t_660[k] = f_8 * lf_440[k]
                   + pa_x[k] * lg_660[k];

        t_661[k] = pb_y[k] * mf_440[k];

        t_662[k] = f_7 * lf_350[k]
                   + pb_z[k] * mf_440[k];

        t_663[k] = f_6 * lf_443[k]
                   + pa_x[k] * lg_663[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, t_667, t_668, pa_x, pb_x, pb_y, lf_445, lf_446, \
                         lf_447, lg_665, mf_442, mf_445, mf_446, \
                         mf_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = pb_y[k] * mf_442[k];

        t_665[k] = f_6 * lf_445[k]
                   + pa_x[k] * lg_665[k];

        t_666[k] = f_5 * lf_446[k]
                   + pb_x[k] * mf_446[k];

        t_667[k] = f_5 * lf_447[k]
                   + pb_x[k] * mf_447[k];

        t_668[k] = pb_y[k] * mf_445[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, t_672, t_673, t_674, pa_x, pb_x, pb_y, lf_449, \
                         lg_670, lg_671, lg_672, lg_674, mf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_5 * lf_449[k]
                   + pb_x[k] * mf_449[k];

        t_670[k] = pa_x[k] * lg_670[k];

        t_671[k] = pa_x[k] * lg_671[k];

        t_672[k] = pa_x[k] * lg_672[k];

        t_673[k] = pb_y[k] * mf_449[k];

        t_674[k] = pa_x[k] * lg_674[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, pb_x, pb_y, pb_z, lf_360, md0_270, \
                         md0_273, md1_270, md1_273, mf_450, mf_451, \
                         mf_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = f_1 * md0_270[k]
                   - f_2 * md1_270[k]
                   + pb_x[k] * mf_450[k];

        t_676[k] = f_0 * lf_360[k]
                   + pb_y[k] * mf_450[k];

        t_677[k] = pb_z[k] * mf_450[k];

        t_678[k] = f_3 * md0_273[k]
                   - f_4 * md1_273[k]
                   + pb_x[k] * mf_453[k];

        t_679[k] = pb_z[k] * mf_451[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, pb_x, md0_275, md1_275, mf_455, \
                         mf_456, mf_457, mf_458, mf_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_3 * md0_275[k]
                   - f_4 * md1_275[k]
                   + pb_x[k] * mf_455[k];

        t_681[k] = pb_x[k] * mf_456[k];

        t_682[k] = pb_x[k] * mf_457[k];

        t_683[k] = pb_x[k] * mf_458[k];

        t_684[k] = pb_x[k] * mf_459[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, pb_y, pb_z, lf_366, lf_369, \
                         md0_273, md0_275, md1_273, md1_275, mf_456, mf_457, \
                         mf_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = f_0 * lf_366[k]
                   + f_1 * md0_273[k]
                   - f_2 * md1_273[k]
                   + pb_y[k] * mf_456[k];

        t_686[k] = pb_z[k] * mf_456[k];

        t_687[k] = f_3 * md0_273[k]
                   - f_4 * md1_273[k]
                   + pb_z[k] * mf_457[k];

        t_688[k] = f_0 * lf_369[k]
                   + pb_y[k] * mf_459[k];

        t_689[k] = f_1 * md0_275[k]
                   - f_2 * md1_275[k]
                   + pb_z[k] * mf_459[k];
    }
}

static auto
compute_prim_mg_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kg0,
                                            const size_t kg1, const size_t lf, const size_t lg,
                                            const size_t md0, const size_t md1, const size_t mf,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 4.0 / p;
    const auto f_8 = 2.0 / p;
    const auto f_9 = 0.5 / alpha;
    const auto f_10 = 0.5 * beta / (alpha * p);
    const auto f_11 = 3.5 / p;
    const auto f_12 = 3.0 / alpha;
    const auto f_13 = 3.0 * beta / (alpha * p);
    const auto f_14 = 1.0 / alpha;
    const auto f_15 = beta / (alpha * p);
    const auto f_16 = 1.5 / p;
    const auto f_17 = 3.0 / p;
    const auto f_18 = 2.5 / alpha;
    const auto f_19 = 2.5 * beta / (alpha * p);
    const auto f_20 = 1.5 / alpha;
    const auto f_21 = 1.5 * beta / (alpha * p);
    const auto f_22 = 2.5 / p;
    const auto f_23 = 2.0 / alpha;
    const auto f_24 = 2.0 * beta / (alpha * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kg0_430 = buffer.data(kg0 + 430);
    const auto *kg0_445 = buffer.data(kg0 + 445);
    const auto *kg0_460 = buffer.data(kg0 + 460);
    const auto *kg0_464 = buffer.data(kg0 + 464);
    const auto *kg0_475 = buffer.data(kg0 + 475);
    const auto *kg0_479 = buffer.data(kg0 + 479);
    const auto *kg0_490 = buffer.data(kg0 + 490);
    const auto *kg0_494 = buffer.data(kg0 + 494);
    const auto *kg0_505 = buffer.data(kg0 + 505);
    const auto *kg0_509 = buffer.data(kg0 + 509);
    const auto *kg0_524 = buffer.data(kg0 + 524);
    const auto *kg0_539 = buffer.data(kg0 + 539);

    const auto *kg1_430 = buffer.data(kg1 + 430);
    const auto *kg1_445 = buffer.data(kg1 + 445);
    const auto *kg1_460 = buffer.data(kg1 + 460);
    const auto *kg1_464 = buffer.data(kg1 + 464);
    const auto *kg1_475 = buffer.data(kg1 + 475);
    const auto *kg1_479 = buffer.data(kg1 + 479);
    const auto *kg1_490 = buffer.data(kg1 + 490);
    const auto *kg1_494 = buffer.data(kg1 + 494);
    const auto *kg1_505 = buffer.data(kg1 + 505);
    const auto *kg1_509 = buffer.data(kg1 + 509);
    const auto *kg1_524 = buffer.data(kg1 + 524);
    const auto *kg1_539 = buffer.data(kg1 + 539);

    const auto *lf_360 = buffer.data(lf + 360);
    const auto *lf_362 = buffer.data(lf + 362);
    const auto *lf_366 = buffer.data(lf + 366);
    const auto *lf_367 = buffer.data(lf + 367);
    const auto *lf_369 = buffer.data(lf + 369);
    const auto *lf_370 = buffer.data(lf + 370);
    const auto *lf_372 = buffer.data(lf + 372);
    const auto *lf_376 = buffer.data(lf + 376);
    const auto *lf_379 = buffer.data(lf + 379);
    const auto *lf_380 = buffer.data(lf + 380);
    const auto *lf_382 = buffer.data(lf + 382);
    const auto *lf_386 = buffer.data(lf + 386);
    const auto *lf_388 = buffer.data(lf + 388);
    const auto *lf_389 = buffer.data(lf + 389);
    const auto *lf_390 = buffer.data(lf + 390);
    const auto *lf_392 = buffer.data(lf + 392);
    const auto *lf_396 = buffer.data(lf + 396);
    const auto *lf_398 = buffer.data(lf + 398);
    const auto *lf_399 = buffer.data(lf + 399);
    const auto *lf_400 = buffer.data(lf + 400);
    const auto *lf_402 = buffer.data(lf + 402);
    const auto *lf_406 = buffer.data(lf + 406);
    const auto *lf_408 = buffer.data(lf + 408);
    const auto *lf_409 = buffer.data(lf + 409);
    const auto *lf_410 = buffer.data(lf + 410);
    const auto *lf_412 = buffer.data(lf + 412);
    const auto *lf_416 = buffer.data(lf + 416);
    const auto *lf_418 = buffer.data(lf + 418);
    const auto *lf_419 = buffer.data(lf + 419);
    const auto *lf_420 = buffer.data(lf + 420);
    const auto *lf_422 = buffer.data(lf + 422);
    const auto *lf_426 = buffer.data(lf + 426);
    const auto *lf_428 = buffer.data(lf + 428);
    const auto *lf_429 = buffer.data(lf + 429);
    const auto *lf_430 = buffer.data(lf + 430);
    const auto *lf_432 = buffer.data(lf + 432);
    const auto *lf_436 = buffer.data(lf + 436);
    const auto *lf_438 = buffer.data(lf + 438);
    const auto *lf_439 = buffer.data(lf + 439);
    const auto *lf_440 = buffer.data(lf + 440);
    const auto *lf_441 = buffer.data(lf + 441);
    const auto *lf_442 = buffer.data(lf + 442);
    const auto *lf_446 = buffer.data(lf + 446);
    const auto *lf_448 = buffer.data(lf + 448);
    const auto *lf_449 = buffer.data(lf + 449);

    const auto *lg_540 = buffer.data(lg + 540);
    const auto *lg_541 = buffer.data(lg + 541);
    const auto *lg_543 = buffer.data(lg + 543);
    const auto *lg_545 = buffer.data(lg + 545);
    const auto *lg_550 = buffer.data(lg + 550);
    const auto *lg_552 = buffer.data(lg + 552);
    const auto *lg_554 = buffer.data(lg + 554);
    const auto *lg_565 = buffer.data(lg + 565);
    const auto *lg_580 = buffer.data(lg + 580);
    const auto *lg_584 = buffer.data(lg + 584);
    const auto *lg_595 = buffer.data(lg + 595);
    const auto *lg_599 = buffer.data(lg + 599);
    const auto *lg_610 = buffer.data(lg + 610);
    const auto *lg_614 = buffer.data(lg + 614);
    const auto *lg_625 = buffer.data(lg + 625);
    const auto *lg_629 = buffer.data(lg + 629);
    const auto *lg_640 = buffer.data(lg + 640);
    const auto *lg_644 = buffer.data(lg + 644);
    const auto *lg_659 = buffer.data(lg + 659);
    const auto *lg_660 = buffer.data(lg + 660);
    const auto *lg_662 = buffer.data(lg + 662);
    const auto *lg_663 = buffer.data(lg + 663);
    const auto *lg_665 = buffer.data(lg + 665);
    const auto *lg_670 = buffer.data(lg + 670);
    const auto *lg_672 = buffer.data(lg + 672);
    const auto *lg_674 = buffer.data(lg + 674);

    const auto *md0_282 = buffer.data(md0 + 282);
    const auto *md0_285 = buffer.data(md0 + 285);
    const auto *md0_287 = buffer.data(md0 + 287);
    const auto *md0_288 = buffer.data(md0 + 288);
    const auto *md0_291 = buffer.data(md0 + 291);
    const auto *md0_293 = buffer.data(md0 + 293);
    const auto *md0_294 = buffer.data(md0 + 294);
    const auto *md0_297 = buffer.data(md0 + 297);
    const auto *md0_299 = buffer.data(md0 + 299);
    const auto *md0_300 = buffer.data(md0 + 300);
    const auto *md0_303 = buffer.data(md0 + 303);
    const auto *md0_305 = buffer.data(md0 + 305);
    const auto *md0_306 = buffer.data(md0 + 306);
    const auto *md0_309 = buffer.data(md0 + 309);
    const auto *md0_311 = buffer.data(md0 + 311);
    const auto *md0_312 = buffer.data(md0 + 312);
    const auto *md0_315 = buffer.data(md0 + 315);
    const auto *md0_317 = buffer.data(md0 + 317);
    const auto *md0_324 = buffer.data(md0 + 324);
    const auto *md0_327 = buffer.data(md0 + 327);
    const auto *md0_329 = buffer.data(md0 + 329);

    const auto *md1_282 = buffer.data(md1 + 282);
    const auto *md1_285 = buffer.data(md1 + 285);
    const auto *md1_287 = buffer.data(md1 + 287);
    const auto *md1_288 = buffer.data(md1 + 288);
    const auto *md1_291 = buffer.data(md1 + 291);
    const auto *md1_293 = buffer.data(md1 + 293);
    const auto *md1_294 = buffer.data(md1 + 294);
    const auto *md1_297 = buffer.data(md1 + 297);
    const auto *md1_299 = buffer.data(md1 + 299);
    const auto *md1_300 = buffer.data(md1 + 300);
    const auto *md1_303 = buffer.data(md1 + 303);
    const auto *md1_305 = buffer.data(md1 + 305);
    const auto *md1_306 = buffer.data(md1 + 306);
    const auto *md1_309 = buffer.data(md1 + 309);
    const auto *md1_311 = buffer.data(md1 + 311);
    const auto *md1_312 = buffer.data(md1 + 312);
    const auto *md1_315 = buffer.data(md1 + 315);
    const auto *md1_317 = buffer.data(md1 + 317);
    const auto *md1_324 = buffer.data(md1 + 324);
    const auto *md1_327 = buffer.data(md1 + 327);
    const auto *md1_329 = buffer.data(md1 + 329);

    const auto *mf_460 = buffer.data(mf + 460);
    const auto *mf_462 = buffer.data(mf + 462);
    const auto *mf_466 = buffer.data(mf + 466);
    const auto *mf_467 = buffer.data(mf + 467);
    const auto *mf_468 = buffer.data(mf + 468);
    const auto *mf_469 = buffer.data(mf + 469);
    const auto *mf_470 = buffer.data(mf + 470);
    const auto *mf_472 = buffer.data(mf + 472);
    const auto *mf_473 = buffer.data(mf + 473);
    const auto *mf_475 = buffer.data(mf + 475);
    const auto *mf_476 = buffer.data(mf + 476);
    const auto *mf_477 = buffer.data(mf + 477);
    const auto *mf_478 = buffer.data(mf + 478);
    const auto *mf_479 = buffer.data(mf + 479);
    const auto *mf_480 = buffer.data(mf + 480);
    const auto *mf_482 = buffer.data(mf + 482);
    const auto *mf_483 = buffer.data(mf + 483);
    const auto *mf_485 = buffer.data(mf + 485);
    const auto *mf_486 = buffer.data(mf + 486);
    const auto *mf_487 = buffer.data(mf + 487);
    const auto *mf_488 = buffer.data(mf + 488);
    const auto *mf_489 = buffer.data(mf + 489);
    const auto *mf_490 = buffer.data(mf + 490);
    const auto *mf_492 = buffer.data(mf + 492);
    const auto *mf_493 = buffer.data(mf + 493);
    const auto *mf_495 = buffer.data(mf + 495);
    const auto *mf_496 = buffer.data(mf + 496);
    const auto *mf_497 = buffer.data(mf + 497);
    const auto *mf_498 = buffer.data(mf + 498);
    const auto *mf_499 = buffer.data(mf + 499);
    const auto *mf_500 = buffer.data(mf + 500);
    const auto *mf_502 = buffer.data(mf + 502);
    const auto *mf_503 = buffer.data(mf + 503);
    const auto *mf_505 = buffer.data(mf + 505);
    const auto *mf_506 = buffer.data(mf + 506);
    const auto *mf_507 = buffer.data(mf + 507);
    const auto *mf_508 = buffer.data(mf + 508);
    const auto *mf_509 = buffer.data(mf + 509);
    const auto *mf_510 = buffer.data(mf + 510);
    const auto *mf_512 = buffer.data(mf + 512);
    const auto *mf_513 = buffer.data(mf + 513);
    const auto *mf_515 = buffer.data(mf + 515);
    const auto *mf_516 = buffer.data(mf + 516);
    const auto *mf_517 = buffer.data(mf + 517);
    const auto *mf_518 = buffer.data(mf + 518);
    const auto *mf_519 = buffer.data(mf + 519);
    const auto *mf_520 = buffer.data(mf + 520);
    const auto *mf_522 = buffer.data(mf + 522);
    const auto *mf_523 = buffer.data(mf + 523);
    const auto *mf_525 = buffer.data(mf + 525);
    const auto *mf_526 = buffer.data(mf + 526);
    const auto *mf_527 = buffer.data(mf + 527);
    const auto *mf_528 = buffer.data(mf + 528);
    const auto *mf_529 = buffer.data(mf + 529);
    const auto *mf_530 = buffer.data(mf + 530);
    const auto *mf_532 = buffer.data(mf + 532);
    const auto *mf_536 = buffer.data(mf + 536);
    const auto *mf_537 = buffer.data(mf + 537);
    const auto *mf_538 = buffer.data(mf + 538);
    const auto *mf_539 = buffer.data(mf + 539);
    const auto *mf_540 = buffer.data(mf + 540);
    const auto *mf_542 = buffer.data(mf + 542);
    const auto *mf_543 = buffer.data(mf + 543);
    const auto *mf_545 = buffer.data(mf + 545);
    const auto *mf_546 = buffer.data(mf + 546);
    const auto *mf_547 = buffer.data(mf + 547);
    const auto *mf_548 = buffer.data(mf + 548);
    const auto *mf_549 = buffer.data(mf + 549);

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, pa_z, pb_y, pb_z, lf_360, lf_372, \
                         lg_540, lg_541, lg_543, mf_460, mf_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = pa_z[k] * lg_540[k];

        t_691[k] = pa_z[k] * lg_541[k];

        t_692[k] = f_5 * lf_360[k]
                   + pb_z[k] * mf_460[k];

        t_693[k] = pa_z[k] * lg_543[k];

        t_694[k] = f_7 * lf_372[k]
                   + pb_y[k] * mf_462[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, t_700, pa_z, pb_x, lf_362, lg_545, \
                         lg_550, mf_466, mf_467, mf_468, mf_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_6 * lf_362[k]
                   + pa_z[k] * lg_545[k];

        t_696[k] = pb_x[k] * mf_466[k];

        t_697[k] = pb_x[k] * mf_467[k];

        t_698[k] = pb_x[k] * mf_468[k];

        t_699[k] = pb_x[k] * mf_469[k];

        t_700[k] = pa_z[k] * lg_550[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pa_z, pb_y, pb_z, lf_366, lf_367, lf_369, \
                         lf_379, lg_552, lg_554, mf_466, mf_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_5 * lf_366[k]
                   + pb_z[k] * mf_466[k];

        t_702[k] = f_6 * lf_367[k]
                   + pa_z[k] * lg_552[k];

        t_703[k] = f_7 * lf_379[k]
                   + pb_y[k] * mf_469[k];

        t_704[k] = f_8 * lf_369[k]
                   + pa_z[k] * lg_554[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, pb_x, pb_y, pb_z, lf_370, lf_380, \
                         md0_282, md0_285, md1_282, md1_285, mf_470, \
                         mf_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_1 * md0_282[k]
                   - f_2 * md1_282[k]
                   + pb_x[k] * mf_470[k];

        t_706[k] = f_11 * lf_380[k]
                   + pb_y[k] * mf_470[k];

        t_707[k] = f_6 * lf_370[k]
                   + pb_z[k] * mf_470[k];

        t_708[k] = f_3 * md0_285[k]
                   - f_4 * md1_285[k]
                   + pb_x[k] * mf_473[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, t_712, t_713, pb_x, pb_y, lf_382, md0_287, \
                         md1_287, mf_472, mf_475, mf_476, mf_477, \
                         mf_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_11 * lf_382[k]
                   + pb_y[k] * mf_472[k];

        t_710[k] = f_3 * md0_287[k]
                   - f_4 * md1_287[k]
                   + pb_x[k] * mf_475[k];

        t_711[k] = pb_x[k] * mf_476[k];

        t_712[k] = pb_x[k] * mf_477[k];

        t_713[k] = pb_x[k] * mf_478[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pa_z, pb_x, pb_z, kg0_430, kg1_430, lf_376, \
                         lg_565, mf_476, mf_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = pb_x[k] * mf_479[k];

        t_715[k] = f_9 * kg0_430[k]
                   - f_10 * kg1_430[k]
                   + pa_z[k] * lg_565[k];

        t_716[k] = f_6 * lf_376[k]
                   + pb_z[k] * mf_476[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pa_y, pb_y, kg0_464, kg1_464, lf_388, lf_389, \
                         lg_584, md0_287, md1_287, mf_478, mf_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_11 * lf_388[k]
                   + f_3 * md0_287[k]
                   - f_4 * md1_287[k]
                   + pb_y[k] * mf_478[k];

        t_718[k] = f_11 * lf_389[k]
                   + pb_y[k] * mf_479[k];

        t_719[k] = f_12 * kg0_464[k]
                   - f_13 * kg1_464[k]
                   + pa_y[k] * lg_584[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pb_x, pb_y, pb_z, lf_380, lf_390, \
                         md0_288, md0_291, md1_288, md1_291, mf_480, \
                         mf_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_1 * md0_288[k]
                   - f_2 * md1_288[k]
                   + pb_x[k] * mf_480[k];

        t_721[k] = f_17 * lf_390[k]
                   + pb_y[k] * mf_480[k];

        t_722[k] = f_16 * lf_380[k]
                   + pb_z[k] * mf_480[k];

        t_723[k] = f_3 * md0_291[k]
                   - f_4 * md1_291[k]
                   + pb_x[k] * mf_483[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, pb_x, pb_y, lf_392, md0_293, \
                         md1_293, mf_482, mf_485, mf_486, mf_487, \
                         mf_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_17 * lf_392[k]
                   + pb_y[k] * mf_482[k];

        t_725[k] = f_3 * md0_293[k]
                   - f_4 * md1_293[k]
                   + pb_x[k] * mf_485[k];

        t_726[k] = pb_x[k] * mf_486[k];

        t_727[k] = pb_x[k] * mf_487[k];

        t_728[k] = pb_x[k] * mf_488[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pa_z, pb_x, pb_z, kg0_445, kg1_445, lf_386, \
                         lg_580, mf_486, mf_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = pb_x[k] * mf_489[k];

        t_730[k] = f_14 * kg0_445[k]
                   - f_15 * kg1_445[k]
                   + pa_z[k] * lg_580[k];

        t_731[k] = f_16 * lf_386[k]
                   + pb_z[k] * mf_486[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pa_y, pb_y, kg0_479, kg1_479, lf_398, lf_399, \
                         lg_599, md0_293, md1_293, mf_488, mf_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_17 * lf_398[k]
                   + f_3 * md0_293[k]
                   - f_4 * md1_293[k]
                   + pb_y[k] * mf_488[k];

        t_733[k] = f_17 * lf_399[k]
                   + pb_y[k] * mf_489[k];

        t_734[k] = f_18 * kg0_479[k]
                   - f_19 * kg1_479[k]
                   + pa_y[k] * lg_599[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, pb_x, pb_y, pb_z, lf_390, lf_400, \
                         md0_294, md0_297, md1_294, md1_297, mf_490, \
                         mf_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_1 * md0_294[k]
                   - f_2 * md1_294[k]
                   + pb_x[k] * mf_490[k];

        t_736[k] = f_22 * lf_400[k]
                   + pb_y[k] * mf_490[k];

        t_737[k] = f_8 * lf_390[k]
                   + pb_z[k] * mf_490[k];

        t_738[k] = f_3 * md0_297[k]
                   - f_4 * md1_297[k]
                   + pb_x[k] * mf_493[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, t_743, pb_x, pb_y, lf_402, md0_299, \
                         md1_299, mf_492, mf_495, mf_496, mf_497, \
                         mf_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_22 * lf_402[k]
                   + pb_y[k] * mf_492[k];

        t_740[k] = f_3 * md0_299[k]
                   - f_4 * md1_299[k]
                   + pb_x[k] * mf_495[k];

        t_741[k] = pb_x[k] * mf_496[k];

        t_742[k] = pb_x[k] * mf_497[k];

        t_743[k] = pb_x[k] * mf_498[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, pa_z, pb_x, pb_z, kg0_460, kg1_460, lf_396, \
                         lg_595, mf_496, mf_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = pb_x[k] * mf_499[k];

        t_745[k] = f_20 * kg0_460[k]
                   - f_21 * kg1_460[k]
                   + pa_z[k] * lg_595[k];

        t_746[k] = f_8 * lf_396[k]
                   + pb_z[k] * mf_496[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, pa_y, pb_y, kg0_494, kg1_494, lf_408, lf_409, \
                         lg_614, md0_299, md1_299, mf_498, mf_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_22 * lf_408[k]
                   + f_3 * md0_299[k]
                   - f_4 * md1_299[k]
                   + pb_y[k] * mf_498[k];

        t_748[k] = f_22 * lf_409[k]
                   + pb_y[k] * mf_499[k];

        t_749[k] = f_23 * kg0_494[k]
                   - f_24 * kg1_494[k]
                   + pa_y[k] * lg_614[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, pb_x, pb_y, pb_z, lf_400, lf_410, \
                         md0_300, md0_303, md1_300, md1_303, mf_500, \
                         mf_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_1 * md0_300[k]
                   - f_2 * md1_300[k]
                   + pb_x[k] * mf_500[k];

        t_751[k] = f_8 * lf_410[k]
                   + pb_y[k] * mf_500[k];

        t_752[k] = f_22 * lf_400[k]
                   + pb_z[k] * mf_500[k];

        t_753[k] = f_3 * md0_303[k]
                   - f_4 * md1_303[k]
                   + pb_x[k] * mf_503[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, t_758, pb_x, pb_y, lf_412, md0_305, \
                         md1_305, mf_502, mf_505, mf_506, mf_507, \
                         mf_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = f_8 * lf_412[k]
                   + pb_y[k] * mf_502[k];

        t_755[k] = f_3 * md0_305[k]
                   - f_4 * md1_305[k]
                   + pb_x[k] * mf_505[k];

        t_756[k] = pb_x[k] * mf_506[k];

        t_757[k] = pb_x[k] * mf_507[k];

        t_758[k] = pb_x[k] * mf_508[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pa_z, pb_x, pb_z, kg0_475, kg1_475, lf_406, \
                         lg_610, mf_506, mf_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = pb_x[k] * mf_509[k];

        t_760[k] = f_23 * kg0_475[k]
                   - f_24 * kg1_475[k]
                   + pa_z[k] * lg_610[k];

        t_761[k] = f_22 * lf_406[k]
                   + pb_z[k] * mf_506[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pa_y, pb_y, kg0_509, kg1_509, lf_418, lf_419, \
                         lg_629, md0_305, md1_305, mf_508, mf_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_8 * lf_418[k]
                   + f_3 * md0_305[k]
                   - f_4 * md1_305[k]
                   + pb_y[k] * mf_508[k];

        t_763[k] = f_8 * lf_419[k]
                   + pb_y[k] * mf_509[k];

        t_764[k] = f_20 * kg0_509[k]
                   - f_21 * kg1_509[k]
                   + pa_y[k] * lg_629[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, pb_x, pb_y, pb_z, lf_410, lf_420, \
                         md0_306, md0_309, md1_306, md1_309, mf_510, \
                         mf_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_1 * md0_306[k]
                   - f_2 * md1_306[k]
                   + pb_x[k] * mf_510[k];

        t_766[k] = f_16 * lf_420[k]
                   + pb_y[k] * mf_510[k];

        t_767[k] = f_17 * lf_410[k]
                   + pb_z[k] * mf_510[k];

        t_768[k] = f_3 * md0_309[k]
                   - f_4 * md1_309[k]
                   + pb_x[k] * mf_513[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, t_772, t_773, pb_x, pb_y, lf_422, md0_311, \
                         md1_311, mf_512, mf_515, mf_516, mf_517, \
                         mf_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = f_16 * lf_422[k]
                   + pb_y[k] * mf_512[k];

        t_770[k] = f_3 * md0_311[k]
                   - f_4 * md1_311[k]
                   + pb_x[k] * mf_515[k];

        t_771[k] = pb_x[k] * mf_516[k];

        t_772[k] = pb_x[k] * mf_517[k];

        t_773[k] = pb_x[k] * mf_518[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, pa_z, pb_x, pb_z, kg0_490, kg1_490, lf_416, \
                         lg_625, mf_516, mf_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = pb_x[k] * mf_519[k];

        t_775[k] = f_18 * kg0_490[k]
                   - f_19 * kg1_490[k]
                   + pa_z[k] * lg_625[k];

        t_776[k] = f_17 * lf_416[k]
                   + pb_z[k] * mf_516[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pa_y, pb_y, kg0_524, kg1_524, lf_428, lf_429, \
                         lg_644, md0_311, md1_311, mf_518, mf_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_16 * lf_428[k]
                   + f_3 * md0_311[k]
                   - f_4 * md1_311[k]
                   + pb_y[k] * mf_518[k];

        t_778[k] = f_16 * lf_429[k]
                   + pb_y[k] * mf_519[k];

        t_779[k] = f_14 * kg0_524[k]
                   - f_15 * kg1_524[k]
                   + pa_y[k] * lg_644[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, pb_x, pb_y, pb_z, lf_420, lf_430, \
                         md0_312, md0_315, md1_312, md1_315, mf_520, \
                         mf_523 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_1 * md0_312[k]
                   - f_2 * md1_312[k]
                   + pb_x[k] * mf_520[k];

        t_781[k] = f_6 * lf_430[k]
                   + pb_y[k] * mf_520[k];

        t_782[k] = f_11 * lf_420[k]
                   + pb_z[k] * mf_520[k];

        t_783[k] = f_3 * md0_315[k]
                   - f_4 * md1_315[k]
                   + pb_x[k] * mf_523[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, t_788, pb_x, pb_y, lf_432, md0_317, \
                         md1_317, mf_522, mf_525, mf_526, mf_527, \
                         mf_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_6 * lf_432[k]
                   + pb_y[k] * mf_522[k];

        t_785[k] = f_3 * md0_317[k]
                   - f_4 * md1_317[k]
                   + pb_x[k] * mf_525[k];

        t_786[k] = pb_x[k] * mf_526[k];

        t_787[k] = pb_x[k] * mf_527[k];

        t_788[k] = pb_x[k] * mf_528[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, pa_z, pb_x, pb_z, kg0_505, kg1_505, lf_426, \
                         lg_640, mf_526, mf_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = pb_x[k] * mf_529[k];

        t_790[k] = f_12 * kg0_505[k]
                   - f_13 * kg1_505[k]
                   + pa_z[k] * lg_640[k];

        t_791[k] = f_11 * lf_426[k]
                   + pb_z[k] * mf_526[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, pa_y, pb_y, kg0_539, kg1_539, lf_438, \
                         lf_439, lg_659, lg_660, md0_317, md1_317, mf_528, \
                         mf_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = f_6 * lf_438[k]
                   + f_3 * md0_317[k]
                   - f_4 * md1_317[k]
                   + pb_y[k] * mf_528[k];

        t_793[k] = f_6 * lf_439[k]
                   + pb_y[k] * mf_529[k];

        t_794[k] = f_9 * kg0_539[k]
                   - f_10 * kg1_539[k]
                   + pa_y[k] * lg_659[k];

        t_795[k] = pa_y[k] * lg_660[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, t_799, t_800, pa_y, pb_y, lf_440, lf_441, \
                         lf_442, lg_662, lg_663, lg_665, mf_530, \
                         mf_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_5 * lf_440[k]
                   + pb_y[k] * mf_530[k];

        t_797[k] = pa_y[k] * lg_662[k];

        t_798[k] = f_6 * lf_441[k]
                   + pa_y[k] * lg_663[k];

        t_799[k] = f_5 * lf_442[k]
                   + pb_y[k] * mf_532[k];

        t_800[k] = pa_y[k] * lg_665[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, t_804, t_805, t_806, pa_y, pb_x, pb_z, lf_436, \
                         lf_446, lg_670, mf_536, mf_537, mf_538, \
                         mf_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = pb_x[k] * mf_536[k];

        t_802[k] = pb_x[k] * mf_537[k];

        t_803[k] = pb_x[k] * mf_538[k];

        t_804[k] = pb_x[k] * mf_539[k];

        t_805[k] = f_8 * lf_446[k]
                   + pa_y[k] * lg_670[k];

        t_806[k] = f_7 * lf_436[k]
                   + pb_z[k] * mf_536[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, t_810, t_811, pa_y, pb_x, pb_y, lf_448, lf_449, \
                         lg_672, lg_674, md0_324, md1_324, mf_539, \
                         mf_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = f_6 * lf_448[k]
                   + pa_y[k] * lg_672[k];

        t_808[k] = f_5 * lf_449[k]
                   + pb_y[k] * mf_539[k];

        t_809[k] = pa_y[k] * lg_674[k];

        t_810[k] = f_1 * md0_324[k]
                   - f_2 * md1_324[k]
                   + pb_x[k] * mf_540[k];

        t_811[k] = pb_y[k] * mf_540[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pb_x, pb_y, pb_z, lf_440, md0_327, \
                         md0_329, md1_327, md1_329, mf_540, mf_542, mf_543, \
                         mf_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_0 * lf_440[k]
                   + pb_z[k] * mf_540[k];

        t_813[k] = f_3 * md0_327[k]
                   - f_4 * md1_327[k]
                   + pb_x[k] * mf_543[k];

        t_814[k] = pb_y[k] * mf_542[k];

        t_815[k] = f_3 * md0_329[k]
                   - f_4 * md1_329[k]
                   + pb_x[k] * mf_545[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, t_819, t_820, t_821, pb_x, pb_y, pb_z, lf_446, \
                         md0_327, md1_327, mf_546, mf_547, mf_548, \
                         mf_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = pb_x[k] * mf_546[k];

        t_817[k] = pb_x[k] * mf_547[k];

        t_818[k] = pb_x[k] * mf_548[k];

        t_819[k] = pb_x[k] * mf_549[k];

        t_820[k] = f_1 * md0_327[k]
                   - f_2 * md1_327[k]
                   + pb_y[k] * mf_546[k];

        t_821[k] = f_0 * lf_446[k]
                   + pb_z[k] * mf_546[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, pb_y, pb_z, lf_449, md0_329, md1_329, mf_548, \
                         mf_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = f_3 * md0_329[k]
                   - f_4 * md1_329[k]
                   + pb_y[k] * mf_548[k];

        t_823[k] = pb_y[k] * mf_549[k];

        t_824[k] = f_0 * lf_449[k]
                   + f_1 * md0_329[k]
                   - f_2 * md1_329[k]
                   + pb_z[k] * mf_549[k];
    }
}

auto
compute_prim_mg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t kg0, const size_t kg1,
                                     const size_t lf, const size_t lg, const size_t md0,
                                     const size_t md1, const size_t mf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    compute_prim_mg_electron_repulsion_0_piece0(buffer, target, pa, pb, kg0, kg1, lf, lg, md0,
                                                md1, mf, ncols, alpha, beta, p);

    compute_prim_mg_electron_repulsion_0_piece1(buffer, target, pa, pb, kg0, kg1, lf, lg, md0,
                                                md1, mf, ncols, alpha, beta, p);

    compute_prim_mg_electron_repulsion_0_piece2(buffer, target, pa, pb, kg0, kg1, lf, lg, md0,
                                                md1, mf, ncols, alpha, beta, p);

    compute_prim_mg_electron_repulsion_0_piece3(buffer, target, pa, pb, kg0, kg1, lf, lg, md0,
                                                md1, mf, ncols, alpha, beta, p);

    compute_prim_mg_electron_repulsion_0_piece4(buffer, target, pa, pb, kg0, kg1, lf, lg, md0,
                                                md1, mf, ncols, alpha, beta, p);

    compute_prim_mg_electron_repulsion_0_piece5(buffer, target, pa, pb, kg0, kg1, lf, lg, md0,
                                                md1, mf, ncols, alpha, beta, p);
}

}  // namespace simdt2ceri
