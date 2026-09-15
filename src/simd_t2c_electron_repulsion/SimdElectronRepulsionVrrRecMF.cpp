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


#include "SimdElectronRepulsionVrrRecMF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_mf_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kf0,
                                            const size_t kf1, const size_t ld, const size_t lf,
                                            const size_t mp0, const size_t mp1, const size_t md,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 4.0 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 3.5 / p;
    const auto f_10 = 3.0 / alpha;
    const auto f_11 = 3.0 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 3.0 / p;
    const auto f_15 = 2.5 / alpha;
    const auto f_16 = 2.5 * beta / (alpha * p);
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / p;
    const auto f_20 = 2.5 / p;
    const auto f_21 = 2.0 / alpha;
    const auto f_22 = 2.0 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kf0_0 = buffer.data(kf0 + 0);
    const auto *kf0_10 = buffer.data(kf0 + 10);
    const auto *kf0_20 = buffer.data(kf0 + 20);
    const auto *kf0_30 = buffer.data(kf0 + 30);
    const auto *kf0_36 = buffer.data(kf0 + 36);
    const auto *kf0_50 = buffer.data(kf0 + 50);
    const auto *kf0_59 = buffer.data(kf0 + 59);
    const auto *kf0_66 = buffer.data(kf0 + 66);
    const auto *kf0_99 = buffer.data(kf0 + 99);
    const auto *kf0_106 = buffer.data(kf0 + 106);
    const auto *kf0_126 = buffer.data(kf0 + 126);
    const auto *kf0_129 = buffer.data(kf0 + 129);

    const auto *kf1_0 = buffer.data(kf1 + 0);
    const auto *kf1_10 = buffer.data(kf1 + 10);
    const auto *kf1_20 = buffer.data(kf1 + 20);
    const auto *kf1_30 = buffer.data(kf1 + 30);
    const auto *kf1_36 = buffer.data(kf1 + 36);
    const auto *kf1_50 = buffer.data(kf1 + 50);
    const auto *kf1_59 = buffer.data(kf1 + 59);
    const auto *kf1_66 = buffer.data(kf1 + 66);
    const auto *kf1_99 = buffer.data(kf1 + 99);
    const auto *kf1_106 = buffer.data(kf1 + 106);
    const auto *kf1_126 = buffer.data(kf1 + 126);
    const auto *kf1_129 = buffer.data(kf1 + 129);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);
    const auto *ld_54 = buffer.data(ld + 54);
    const auto *ld_57 = buffer.data(ld + 57);
    const auto *ld_59 = buffer.data(ld + 59);
    const auto *ld_63 = buffer.data(ld + 63);
    const auto *ld_65 = buffer.data(ld + 65);
    const auto *ld_70 = buffer.data(ld + 70);
    const auto *ld_71 = buffer.data(ld + 71);
    const auto *ld_75 = buffer.data(ld + 75);
    const auto *ld_76 = buffer.data(ld + 76);
    const auto *ld_77 = buffer.data(ld + 77);
    const auto *ld_81 = buffer.data(ld + 81);
    const auto *ld_82 = buffer.data(ld + 82);
    const auto *ld_87 = buffer.data(ld + 87);
    const auto *ld_89 = buffer.data(ld + 89);

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_126 = buffer.data(lf + 126);
    const auto *lf_129 = buffer.data(lf + 129);

    const auto *mp0_0 = buffer.data(mp0 + 0);
    const auto *mp0_1 = buffer.data(mp0 + 1);
    const auto *mp0_2 = buffer.data(mp0 + 2);
    const auto *mp0_11 = buffer.data(mp0 + 11);
    const auto *mp0_16 = buffer.data(mp0 + 16);
    const auto *mp0_20 = buffer.data(mp0 + 20);
    const auto *mp0_28 = buffer.data(mp0 + 28);
    const auto *mp0_32 = buffer.data(mp0 + 32);

    const auto *mp1_0 = buffer.data(mp1 + 0);
    const auto *mp1_1 = buffer.data(mp1 + 1);
    const auto *mp1_2 = buffer.data(mp1 + 2);
    const auto *mp1_11 = buffer.data(mp1 + 11);
    const auto *mp1_16 = buffer.data(mp1 + 16);
    const auto *mp1_20 = buffer.data(mp1 + 20);
    const auto *mp1_28 = buffer.data(mp1 + 28);
    const auto *mp1_32 = buffer.data(mp1 + 32);

    const auto *md_0 = buffer.data(md + 0);
    const auto *md_2 = buffer.data(md + 2);
    const auto *md_3 = buffer.data(md + 3);
    const auto *md_5 = buffer.data(md + 5);
    const auto *md_6 = buffer.data(md + 6);
    const auto *md_7 = buffer.data(md + 7);
    const auto *md_9 = buffer.data(md + 9);
    const auto *md_11 = buffer.data(md + 11);
    const auto *md_12 = buffer.data(md + 12);
    const auto *md_14 = buffer.data(md + 14);
    const auto *md_15 = buffer.data(md + 15);
    const auto *md_17 = buffer.data(md + 17);
    const auto *md_18 = buffer.data(md + 18);
    const auto *md_19 = buffer.data(md + 19);
    const auto *md_21 = buffer.data(md + 21);
    const auto *md_23 = buffer.data(md + 23);
    const auto *md_27 = buffer.data(md + 27);
    const auto *md_28 = buffer.data(md + 28);
    const auto *md_29 = buffer.data(md + 29);
    const auto *md_30 = buffer.data(md + 30);
    const auto *md_32 = buffer.data(md + 32);
    const auto *md_33 = buffer.data(md + 33);
    const auto *md_35 = buffer.data(md + 35);
    const auto *md_36 = buffer.data(md + 36);
    const auto *md_37 = buffer.data(md + 37);
    const auto *md_39 = buffer.data(md + 39);
    const auto *md_41 = buffer.data(md + 41);
    const auto *md_42 = buffer.data(md + 42);
    const auto *md_45 = buffer.data(md + 45);
    const auto *md_46 = buffer.data(md + 46);
    const auto *md_47 = buffer.data(md + 47);
    const auto *md_48 = buffer.data(md + 48);
    const auto *md_51 = buffer.data(md + 51);
    const auto *md_52 = buffer.data(md + 52);
    const auto *md_53 = buffer.data(md + 53);
    const auto *md_54 = buffer.data(md + 54);
    const auto *md_56 = buffer.data(md + 56);
    const auto *md_57 = buffer.data(md + 57);
    const auto *md_59 = buffer.data(md + 59);
    const auto *md_60 = buffer.data(md + 60);
    const auto *md_61 = buffer.data(md + 61);
    const auto *md_63 = buffer.data(md + 63);
    const auto *md_65 = buffer.data(md + 65);
    const auto *md_66 = buffer.data(md + 66);
    const auto *md_69 = buffer.data(md + 69);
    const auto *md_70 = buffer.data(md + 70);
    const auto *md_71 = buffer.data(md + 71);
    const auto *md_72 = buffer.data(md + 72);
    const auto *md_75 = buffer.data(md + 75);
    const auto *md_76 = buffer.data(md + 76);
    const auto *md_77 = buffer.data(md + 77);
    const auto *md_78 = buffer.data(md + 78);
    const auto *md_81 = buffer.data(md + 81);
    const auto *md_82 = buffer.data(md + 82);
    const auto *md_83 = buffer.data(md + 83);
    const auto *md_84 = buffer.data(md + 84);
    const auto *md_86 = buffer.data(md + 86);
    const auto *md_87 = buffer.data(md + 87);
    const auto *md_89 = buffer.data(md + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ld_0, ld_3, mp0_0, mp1_0, \
                         md_0, md_2, md_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ld_0[k]
                 + f_1 * mp0_0[k]
                 - f_2 * mp1_0[k]
                 + pb_x[k] * md_0[k];

        t_1[k] = pb_y[k] * md_0[k];

        t_2[k] = pb_z[k] * md_0[k];

        t_3[k] = f_0 * ld_3[k]
                 + pb_x[k] * md_3[k];

        t_4[k] = pb_y[k] * md_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, ld_5, mp0_1, mp0_2, mp1_1, \
                         mp1_2, md_3, md_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * ld_5[k]
                 + pb_x[k] * md_5[k];

        t_6[k] = f_1 * mp0_1[k]
                 - f_2 * mp1_1[k]
                 + pb_y[k] * md_3[k];

        t_7[k] = pb_z[k] * md_3[k];

        t_8[k] = pb_y[k] * md_5[k];

        t_9[k] = f_1 * mp0_2[k]
                 - f_2 * mp1_2[k]
                 + pb_z[k] * md_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, ld_0, ld_9, \
                         lf_0, md_6, md_7, md_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * lf_0[k];

        t_11[k] = f_3 * ld_0[k]
                  + pb_y[k] * md_6[k];

        t_12[k] = pb_z[k] * md_6[k];

        t_13[k] = f_4 * ld_9[k]
                  + pb_x[k] * md_9[k];

        t_14[k] = pb_z[k] * md_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, ld_3, ld_5, lf_5, \
                         lf_6, lf_9, md_9, md_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * lf_5[k];

        t_16[k] = f_5 * ld_3[k]
                  + pa_y[k] * lf_6[k];

        t_17[k] = pb_z[k] * md_9[k];

        t_18[k] = f_3 * ld_5[k]
                  + pb_y[k] * md_11[k];

        t_19[k] = pa_y[k] * lf_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, ld_0, lf_0, lf_3, \
                         md_12, md_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * lf_0[k];

        t_21[k] = pb_y[k] * md_12[k];

        t_22[k] = f_3 * ld_0[k]
                  + pb_z[k] * md_12[k];

        t_23[k] = pa_z[k] * lf_3[k];

        t_24[k] = pb_y[k] * md_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, pb_z, ld_3, ld_5, \
                         ld_17, lf_6, lf_9, md_15, md_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * ld_17[k]
                  + pb_x[k] * md_17[k];

        t_26[k] = pa_z[k] * lf_6[k];

        t_27[k] = f_3 * ld_3[k]
                  + pb_z[k] * md_15[k];

        t_28[k] = pb_y[k] * md_17[k];

        t_29[k] = f_5 * ld_5[k]
                  + pa_z[k] * lf_9[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, pb_z, kf0_0, kf1_0, ld_6, \
                         ld_21, lf_10, md_18, md_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_6 * kf0_0[k]
                  - f_7 * kf1_0[k]
                  + pa_y[k] * lf_10[k];

        t_31[k] = f_8 * ld_6[k]
                  + pb_y[k] * md_18[k];

        t_32[k] = pb_z[k] * md_18[k];

        t_33[k] = f_9 * ld_21[k]
                  + pb_x[k] * md_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, kf0_36, kf1_36, ld_23, \
                         lf_36, md_19, md_21, md_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_z[k] * md_19[k];

        t_35[k] = f_9 * ld_23[k]
                  + pb_x[k] * md_23[k];

        t_36[k] = f_10 * kf0_36[k]
                  - f_11 * kf1_36[k]
                  + pa_x[k] * lf_36[k];

        t_37[k] = pb_z[k] * md_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, pb_z, ld_11, lf_11, \
                         lf_20, lf_22, mp0_11, mp1_11, md_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * ld_11[k]
                  + pb_y[k] * md_23[k];

        t_39[k] = f_1 * mp0_11[k]
                  - f_2 * mp1_11[k]
                  + pb_z[k] * md_23[k];

        t_40[k] = pa_y[k] * lf_20[k];

        t_41[k] = pa_z[k] * lf_11[k];

        t_42[k] = pa_y[k] * lf_22[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, pb_z, ld_9, ld_28, \
                         lf_13, lf_16, lf_25, md_27, md_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pa_z[k] * lf_13[k];

        t_44[k] = f_9 * ld_28[k]
                  + pb_x[k] * md_28[k];

        t_45[k] = pa_y[k] * lf_25[k];

        t_46[k] = pa_z[k] * lf_16[k];

        t_47[k] = f_3 * ld_9[k]
                  + pb_z[k] * md_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pa_z, pb_y, kf0_0, kf1_0, ld_17, lf_20, \
                         lf_29, md_29, md_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * ld_17[k]
                  + pb_y[k] * md_29[k];

        t_49[k] = pa_y[k] * lf_29[k];

        t_50[k] = f_6 * kf0_0[k]
                  - f_7 * kf1_0[k]
                  + pa_z[k] * lf_20[k];

        t_51[k] = pb_y[k] * md_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, ld_12, ld_33, ld_35, md_30, \
                         md_32, md_33, md_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_8 * ld_12[k]
                  + pb_z[k] * md_30[k];

        t_53[k] = f_9 * ld_33[k]
                  + pb_x[k] * md_33[k];

        t_54[k] = pb_y[k] * md_32[k];

        t_55[k] = f_9 * ld_35[k]
                  + pb_x[k] * md_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_y, pb_z, kf0_59, kf1_59, ld_15, \
                         lf_59, mp0_16, mp1_16, md_33, md_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * mp0_16[k]
                  - f_2 * mp1_16[k]
                  + pb_y[k] * md_33[k];

        t_57[k] = f_8 * ld_15[k]
                  + pb_z[k] * md_33[k];

        t_58[k] = pb_y[k] * md_35[k];

        t_59[k] = f_10 * kf0_59[k]
                  - f_11 * kf1_59[k]
                  + pa_x[k] * lf_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pb_x, pb_y, pb_z, kf0_10, kf1_10, \
                         ld_18, ld_39, lf_30, md_36, md_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_12 * kf0_10[k]
                  - f_13 * kf1_10[k]
                  + pa_y[k] * lf_30[k];

        t_61[k] = f_5 * ld_18[k]
                  + pb_y[k] * md_36[k];

        t_62[k] = pb_z[k] * md_36[k];

        t_63[k] = f_14 * ld_39[k]
                  + pb_x[k] * md_39[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_z, kf0_66, kf1_66, ld_41, \
                         lf_66, md_37, md_39, md_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_z[k] * md_37[k];

        t_65[k] = f_14 * ld_41[k]
                  + pb_x[k] * md_41[k];

        t_66[k] = f_15 * kf0_66[k]
                  - f_16 * kf1_66[k]
                  + pa_x[k] * lf_66[k];

        t_67[k] = pb_z[k] * md_39[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_z, pb_y, pb_z, ld_18, ld_23, lf_30, \
                         lf_31, mp0_20, mp1_20, md_41, md_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_5 * ld_23[k]
                  + pb_y[k] * md_41[k];

        t_69[k] = f_1 * mp0_20[k]
                  - f_2 * mp1_20[k]
                  + pb_z[k] * md_41[k];

        t_70[k] = pa_z[k] * lf_30[k];

        t_71[k] = pa_z[k] * lf_31[k];

        t_72[k] = f_3 * ld_18[k]
                  + pb_z[k] * md_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, pb_z, ld_21, ld_46, ld_47, \
                         lf_33, lf_36, md_45, md_46, md_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_z[k] * lf_33[k];

        t_74[k] = f_14 * ld_46[k]
                  + pb_x[k] * md_46[k];

        t_75[k] = f_14 * ld_47[k]
                  + pb_x[k] * md_47[k];

        t_76[k] = pa_z[k] * lf_36[k];

        t_77[k] = f_3 * ld_21[k]
                  + pb_z[k] * md_45[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pa_z, pb_y, ld_23, ld_29, ld_30, \
                         lf_39, lf_50, lf_52, md_47, md_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_8 * ld_29[k]
                  + pb_y[k] * md_47[k];

        t_79[k] = f_5 * ld_23[k]
                  + pa_z[k] * lf_39[k];

        t_80[k] = pa_y[k] * lf_50[k];

        t_81[k] = f_3 * ld_30[k]
                  + pb_y[k] * md_48[k];

        t_82[k] = pa_y[k] * lf_52[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, pb_x, pb_z, ld_27, ld_33, ld_51, \
                         ld_52, lf_55, lf_56, md_51, md_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_14 * ld_51[k]
                  + pb_x[k] * md_51[k];

        t_84[k] = f_14 * ld_52[k]
                  + pb_x[k] * md_52[k];

        t_85[k] = pa_y[k] * lf_55[k];

        t_86[k] = f_5 * ld_33[k]
                  + pa_y[k] * lf_56[k];

        t_87[k] = f_8 * ld_27[k]
                  + pb_z[k] * md_51[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, kf0_20, kf1_20, ld_35, \
                         lf_50, lf_59, md_53, md_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_3 * ld_35[k]
                  + pb_y[k] * md_53[k];

        t_89[k] = pa_y[k] * lf_59[k];

        t_90[k] = f_12 * kf0_20[k]
                  - f_13 * kf1_20[k]
                  + pa_z[k] * lf_50[k];

        t_91[k] = pb_y[k] * md_54[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, ld_30, ld_57, ld_59, md_54, \
                         md_56, md_57, md_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_5 * ld_30[k]
                  + pb_z[k] * md_54[k];

        t_93[k] = f_14 * ld_57[k]
                  + pb_x[k] * md_57[k];

        t_94[k] = pb_y[k] * md_56[k];

        t_95[k] = f_14 * ld_59[k]
                  + pb_x[k] * md_59[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_x, pb_y, pb_z, kf0_99, kf1_99, ld_33, \
                         lf_99, mp0_28, mp1_28, md_57, md_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * mp0_28[k]
                  - f_2 * mp1_28[k]
                  + pb_y[k] * md_57[k];

        t_97[k] = f_5 * ld_33[k]
                  + pb_z[k] * md_57[k];

        t_98[k] = pb_y[k] * md_59[k];

        t_99[k] = f_15 * kf0_99[k]
                  - f_16 * kf1_99[k]
                  + pa_x[k] * lf_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_x, pb_y, pb_z, kf0_30, kf1_30, \
                         ld_36, ld_63, lf_60, md_60, md_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_17 * kf0_30[k]
                   - f_18 * kf1_30[k]
                   + pa_y[k] * lf_60[k];

        t_101[k] = f_19 * ld_36[k]
                   + pb_y[k] * md_60[k];

        t_102[k] = pb_z[k] * md_60[k];

        t_103[k] = f_20 * ld_63[k]
                   + pb_x[k] * md_63[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_x, pb_x, pb_z, kf0_106, kf1_106, \
                         ld_65, lf_106, md_61, md_63, md_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_z[k] * md_61[k];

        t_105[k] = f_20 * ld_65[k]
                   + pb_x[k] * md_65[k];

        t_106[k] = f_21 * kf0_106[k]
                   - f_22 * kf1_106[k]
                   + pa_x[k] * lf_106[k];

        t_107[k] = pb_z[k] * md_63[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, ld_36, ld_41, \
                         lf_60, lf_61, mp0_32, mp1_32, md_65, md_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_19 * ld_41[k]
                   + pb_y[k] * md_65[k];

        t_109[k] = f_1 * mp0_32[k]
                   - f_2 * mp1_32[k]
                   + pb_z[k] * md_65[k];

        t_110[k] = pa_z[k] * lf_60[k];

        t_111[k] = pa_z[k] * lf_61[k];

        t_112[k] = f_3 * ld_36[k]
                   + pb_z[k] * md_66[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_z, pb_x, pb_z, ld_39, ld_70, \
                         ld_71, lf_63, lf_66, md_69, md_70, md_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_z[k] * lf_63[k];

        t_114[k] = f_20 * ld_70[k]
                   + pb_x[k] * md_70[k];

        t_115[k] = f_20 * ld_71[k]
                   + pb_x[k] * md_71[k];

        t_116[k] = pa_z[k] * lf_66[k];

        t_117[k] = f_3 * ld_39[k]
                   + pb_z[k] * md_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pa_z, pb_y, kf0_50, kf1_50, ld_41, \
                         ld_47, ld_48, lf_69, lf_80, md_71, md_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_5 * ld_47[k]
                   + pb_y[k] * md_71[k];

        t_119[k] = f_5 * ld_41[k]
                   + pa_z[k] * lf_69[k];

        t_120[k] = f_6 * kf0_50[k]
                   - f_7 * kf1_50[k]
                   + pa_y[k] * lf_80[k];

        t_121[k] = f_8 * ld_48[k]
                   + pb_y[k] * md_72[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pb_x, pb_z, ld_42, ld_75, ld_76, ld_77, \
                         md_72, md_75, md_76, md_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * ld_42[k]
                   + pb_z[k] * md_72[k];

        t_123[k] = f_20 * ld_75[k]
                   + pb_x[k] * md_75[k];

        t_124[k] = f_20 * ld_76[k]
                   + pb_x[k] * md_76[k];

        t_125[k] = f_20 * ld_77[k]
                   + pb_x[k] * md_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_x, pb_y, pb_z, kf0_126, kf1_126, ld_45, \
                         ld_53, lf_126, md_75, md_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_21 * kf0_126[k]
                   - f_22 * kf1_126[k]
                   + pa_x[k] * lf_126[k];

        t_127[k] = f_8 * ld_45[k]
                   + pb_z[k] * md_75[k];

        t_128[k] = f_8 * ld_53[k]
                   + pb_y[k] * md_77[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_x, pa_y, pb_y, kf0_129, kf1_129, \
                         ld_54, lf_90, lf_92, lf_129, md_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_21 * kf0_129[k]
                   - f_22 * kf1_129[k]
                   + pa_x[k] * lf_129[k];

        t_130[k] = pa_y[k] * lf_90[k];

        t_131[k] = f_3 * ld_54[k]
                   + pb_y[k] * md_78[k];

        t_132[k] = pa_y[k] * lf_92[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, pa_y, pb_x, pb_z, ld_51, ld_57, \
                         ld_81, ld_82, lf_95, lf_96, md_81, md_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_20 * ld_81[k]
                   + pb_x[k] * md_81[k];

        t_134[k] = f_20 * ld_82[k]
                   + pb_x[k] * md_82[k];

        t_135[k] = pa_y[k] * lf_95[k];

        t_136[k] = f_5 * ld_57[k]
                   + pa_y[k] * lf_96[k];

        t_137[k] = f_5 * ld_51[k]
                   + pb_z[k] * md_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_y, pa_z, pb_y, kf0_50, kf1_50, ld_59, \
                         lf_90, lf_99, md_83, md_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_3 * ld_59[k]
                   + pb_y[k] * md_83[k];

        t_139[k] = pa_y[k] * lf_99[k];

        t_140[k] = f_17 * kf0_50[k]
                   - f_18 * kf1_50[k]
                   + pa_z[k] * lf_90[k];

        t_141[k] = pb_y[k] * md_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_x, pb_y, pb_z, ld_54, ld_87, ld_89, \
                         md_84, md_86, md_87, md_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_19 * ld_54[k]
                   + pb_z[k] * md_84[k];

        t_143[k] = f_20 * ld_87[k]
                   + pb_x[k] * md_87[k];

        t_144[k] = pb_y[k] * md_86[k];

        t_145[k] = f_20 * ld_89[k]
                   + pb_x[k] * md_89[k];
    }
}

static auto
compute_prim_mf_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kf0,
                                            const size_t kf1, const size_t ld, const size_t lf,
                                            const size_t mp0, const size_t mp1, const size_t md,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 3.0 / p;
    const auto f_15 = 2.5 / alpha;
    const auto f_16 = 2.5 * beta / (alpha * p);
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / p;
    const auto f_20 = 2.5 / p;
    const auto f_21 = 2.0 / alpha;
    const auto f_22 = 2.0 * beta / (alpha * p);

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
    auto *t_277 = buffer.data(target + 277);
    auto *t_278 = buffer.data(target + 278);
    auto *t_279 = buffer.data(target + 279);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kf0_60 = buffer.data(kf0 + 60);
    const auto *kf0_80 = buffer.data(kf0 + 80);
    const auto *kf0_90 = buffer.data(kf0 + 90);
    const auto *kf0_100 = buffer.data(kf0 + 100);
    const auto *kf0_120 = buffer.data(kf0 + 120);
    const auto *kf0_130 = buffer.data(kf0 + 130);
    const auto *kf0_140 = buffer.data(kf0 + 140);
    const auto *kf0_149 = buffer.data(kf0 + 149);
    const auto *kf0_156 = buffer.data(kf0 + 156);
    const auto *kf0_176 = buffer.data(kf0 + 176);
    const auto *kf0_179 = buffer.data(kf0 + 179);
    const auto *kf0_186 = buffer.data(kf0 + 186);
    const auto *kf0_189 = buffer.data(kf0 + 189);
    const auto *kf0_209 = buffer.data(kf0 + 209);
    const auto *kf0_216 = buffer.data(kf0 + 216);
    const auto *kf0_236 = buffer.data(kf0 + 236);
    const auto *kf0_239 = buffer.data(kf0 + 239);
    const auto *kf0_246 = buffer.data(kf0 + 246);
    const auto *kf0_249 = buffer.data(kf0 + 249);
    const auto *kf0_256 = buffer.data(kf0 + 256);
    const auto *kf0_259 = buffer.data(kf0 + 259);
    const auto *kf0_279 = buffer.data(kf0 + 279);

    const auto *kf1_60 = buffer.data(kf1 + 60);
    const auto *kf1_80 = buffer.data(kf1 + 80);
    const auto *kf1_90 = buffer.data(kf1 + 90);
    const auto *kf1_100 = buffer.data(kf1 + 100);
    const auto *kf1_120 = buffer.data(kf1 + 120);
    const auto *kf1_130 = buffer.data(kf1 + 130);
    const auto *kf1_140 = buffer.data(kf1 + 140);
    const auto *kf1_149 = buffer.data(kf1 + 149);
    const auto *kf1_156 = buffer.data(kf1 + 156);
    const auto *kf1_176 = buffer.data(kf1 + 176);
    const auto *kf1_179 = buffer.data(kf1 + 179);
    const auto *kf1_186 = buffer.data(kf1 + 186);
    const auto *kf1_189 = buffer.data(kf1 + 189);
    const auto *kf1_209 = buffer.data(kf1 + 209);
    const auto *kf1_216 = buffer.data(kf1 + 216);
    const auto *kf1_236 = buffer.data(kf1 + 236);
    const auto *kf1_239 = buffer.data(kf1 + 239);
    const auto *kf1_246 = buffer.data(kf1 + 246);
    const auto *kf1_249 = buffer.data(kf1 + 249);
    const auto *kf1_256 = buffer.data(kf1 + 256);
    const auto *kf1_259 = buffer.data(kf1 + 259);
    const auto *kf1_279 = buffer.data(kf1 + 279);

    const auto *ld_57 = buffer.data(ld + 57);
    const auto *ld_60 = buffer.data(ld + 60);
    const auto *ld_63 = buffer.data(ld + 63);
    const auto *ld_65 = buffer.data(ld + 65);
    const auto *ld_66 = buffer.data(ld + 66);
    const auto *ld_69 = buffer.data(ld + 69);
    const auto *ld_71 = buffer.data(ld + 71);
    const auto *ld_72 = buffer.data(ld + 72);
    const auto *ld_75 = buffer.data(ld + 75);
    const auto *ld_77 = buffer.data(ld + 77);
    const auto *ld_78 = buffer.data(ld + 78);
    const auto *ld_81 = buffer.data(ld + 81);
    const auto *ld_83 = buffer.data(ld + 83);
    const auto *ld_84 = buffer.data(ld + 84);
    const auto *ld_87 = buffer.data(ld + 87);
    const auto *ld_89 = buffer.data(ld + 89);
    const auto *ld_90 = buffer.data(ld + 90);
    const auto *ld_93 = buffer.data(ld + 93);
    const auto *ld_95 = buffer.data(ld + 95);
    const auto *ld_96 = buffer.data(ld + 96);
    const auto *ld_99 = buffer.data(ld + 99);
    const auto *ld_100 = buffer.data(ld + 100);
    const auto *ld_101 = buffer.data(ld + 101);
    const auto *ld_102 = buffer.data(ld + 102);
    const auto *ld_105 = buffer.data(ld + 105);
    const auto *ld_106 = buffer.data(ld + 106);
    const auto *ld_107 = buffer.data(ld + 107);
    const auto *ld_108 = buffer.data(ld + 108);
    const auto *ld_111 = buffer.data(ld + 111);
    const auto *ld_112 = buffer.data(ld + 112);
    const auto *ld_113 = buffer.data(ld + 113);
    const auto *ld_114 = buffer.data(ld + 114);
    const auto *ld_117 = buffer.data(ld + 117);
    const auto *ld_118 = buffer.data(ld + 118);
    const auto *ld_119 = buffer.data(ld + 119);
    const auto *ld_120 = buffer.data(ld + 120);
    const auto *ld_123 = buffer.data(ld + 123);
    const auto *ld_125 = buffer.data(ld + 125);
    const auto *ld_129 = buffer.data(ld + 129);
    const auto *ld_131 = buffer.data(ld + 131);
    const auto *ld_136 = buffer.data(ld + 136);
    const auto *ld_137 = buffer.data(ld + 137);
    const auto *ld_141 = buffer.data(ld + 141);
    const auto *ld_142 = buffer.data(ld + 142);
    const auto *ld_143 = buffer.data(ld + 143);
    const auto *ld_147 = buffer.data(ld + 147);
    const auto *ld_148 = buffer.data(ld + 148);
    const auto *ld_149 = buffer.data(ld + 149);
    const auto *ld_153 = buffer.data(ld + 153);
    const auto *ld_154 = buffer.data(ld + 154);
    const auto *ld_155 = buffer.data(ld + 155);
    const auto *ld_159 = buffer.data(ld + 159);
    const auto *ld_160 = buffer.data(ld + 160);
    const auto *ld_165 = buffer.data(ld + 165);
    const auto *ld_167 = buffer.data(ld + 167);

    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_109 = buffer.data(lf + 109);
    const auto *lf_120 = buffer.data(lf + 120);
    const auto *lf_130 = buffer.data(lf + 130);
    const auto *lf_140 = buffer.data(lf + 140);
    const auto *lf_142 = buffer.data(lf + 142);
    const auto *lf_145 = buffer.data(lf + 145);
    const auto *lf_146 = buffer.data(lf + 146);
    const auto *lf_149 = buffer.data(lf + 149);
    const auto *lf_150 = buffer.data(lf + 150);
    const auto *lf_151 = buffer.data(lf + 151);
    const auto *lf_153 = buffer.data(lf + 153);
    const auto *lf_156 = buffer.data(lf + 156);
    const auto *lf_159 = buffer.data(lf + 159);
    const auto *lf_170 = buffer.data(lf + 170);
    const auto *lf_176 = buffer.data(lf + 176);
    const auto *lf_179 = buffer.data(lf + 179);
    const auto *lf_180 = buffer.data(lf + 180);
    const auto *lf_186 = buffer.data(lf + 186);
    const auto *lf_189 = buffer.data(lf + 189);
    const auto *lf_190 = buffer.data(lf + 190);
    const auto *lf_200 = buffer.data(lf + 200);
    const auto *lf_202 = buffer.data(lf + 202);
    const auto *lf_205 = buffer.data(lf + 205);
    const auto *lf_206 = buffer.data(lf + 206);
    const auto *lf_209 = buffer.data(lf + 209);
    const auto *lf_216 = buffer.data(lf + 216);
    const auto *lf_236 = buffer.data(lf + 236);
    const auto *lf_239 = buffer.data(lf + 239);
    const auto *lf_246 = buffer.data(lf + 246);
    const auto *lf_249 = buffer.data(lf + 249);
    const auto *lf_256 = buffer.data(lf + 256);
    const auto *lf_259 = buffer.data(lf + 259);
    const auto *lf_279 = buffer.data(lf + 279);

    const auto *mp0_43 = buffer.data(mp0 + 43);
    const auto *mp0_47 = buffer.data(mp0 + 47);
    const auto *mp0_61 = buffer.data(mp0 + 61);
    const auto *mp0_65 = buffer.data(mp0 + 65);
    const auto *mp0_82 = buffer.data(mp0 + 82);

    const auto *mp1_43 = buffer.data(mp1 + 43);
    const auto *mp1_47 = buffer.data(mp1 + 47);
    const auto *mp1_61 = buffer.data(mp1 + 61);
    const auto *mp1_65 = buffer.data(mp1 + 65);
    const auto *mp1_82 = buffer.data(mp1 + 82);

    const auto *md_87 = buffer.data(md + 87);
    const auto *md_89 = buffer.data(md + 89);
    const auto *md_90 = buffer.data(md + 90);
    const auto *md_91 = buffer.data(md + 91);
    const auto *md_93 = buffer.data(md + 93);
    const auto *md_95 = buffer.data(md + 95);
    const auto *md_96 = buffer.data(md + 96);
    const auto *md_99 = buffer.data(md + 99);
    const auto *md_100 = buffer.data(md + 100);
    const auto *md_101 = buffer.data(md + 101);
    const auto *md_102 = buffer.data(md + 102);
    const auto *md_105 = buffer.data(md + 105);
    const auto *md_106 = buffer.data(md + 106);
    const auto *md_107 = buffer.data(md + 107);
    const auto *md_108 = buffer.data(md + 108);
    const auto *md_111 = buffer.data(md + 111);
    const auto *md_112 = buffer.data(md + 112);
    const auto *md_113 = buffer.data(md + 113);
    const auto *md_114 = buffer.data(md + 114);
    const auto *md_117 = buffer.data(md + 117);
    const auto *md_118 = buffer.data(md + 118);
    const auto *md_119 = buffer.data(md + 119);
    const auto *md_120 = buffer.data(md + 120);
    const auto *md_122 = buffer.data(md + 122);
    const auto *md_123 = buffer.data(md + 123);
    const auto *md_125 = buffer.data(md + 125);
    const auto *md_126 = buffer.data(md + 126);
    const auto *md_127 = buffer.data(md + 127);
    const auto *md_129 = buffer.data(md + 129);
    const auto *md_131 = buffer.data(md + 131);
    const auto *md_132 = buffer.data(md + 132);
    const auto *md_135 = buffer.data(md + 135);
    const auto *md_136 = buffer.data(md + 136);
    const auto *md_137 = buffer.data(md + 137);
    const auto *md_138 = buffer.data(md + 138);
    const auto *md_141 = buffer.data(md + 141);
    const auto *md_142 = buffer.data(md + 142);
    const auto *md_143 = buffer.data(md + 143);
    const auto *md_144 = buffer.data(md + 144);
    const auto *md_147 = buffer.data(md + 147);
    const auto *md_148 = buffer.data(md + 148);
    const auto *md_149 = buffer.data(md + 149);
    const auto *md_150 = buffer.data(md + 150);
    const auto *md_153 = buffer.data(md + 153);
    const auto *md_154 = buffer.data(md + 154);
    const auto *md_155 = buffer.data(md + 155);
    const auto *md_156 = buffer.data(md + 156);
    const auto *md_159 = buffer.data(md + 159);
    const auto *md_160 = buffer.data(md + 160);
    const auto *md_161 = buffer.data(md + 161);
    const auto *md_162 = buffer.data(md + 162);
    const auto *md_164 = buffer.data(md + 164);
    const auto *md_165 = buffer.data(md + 165);
    const auto *md_167 = buffer.data(md + 167);

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pb_y, pb_z, kf0_149, kf1_149, \
                         ld_57, lf_149, mp0_43, mp1_43, md_87, md_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * mp0_43[k]
                   - f_2 * mp1_43[k]
                   + pb_y[k] * md_87[k];

        t_147[k] = f_19 * ld_57[k]
                   + pb_z[k] * md_87[k];

        t_148[k] = pb_y[k] * md_89[k];

        t_149[k] = f_21 * kf0_149[k]
                   - f_22 * kf1_149[k]
                   + pa_x[k] * lf_149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_y, pb_x, pb_y, pb_z, kf0_60, kf1_60, \
                         ld_60, ld_93, lf_100, md_90, md_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_21 * kf0_60[k]
                   - f_22 * kf1_60[k]
                   + pa_y[k] * lf_100[k];

        t_151[k] = f_20 * ld_60[k]
                   + pb_y[k] * md_90[k];

        t_152[k] = pb_z[k] * md_90[k];

        t_153[k] = f_19 * ld_93[k]
                   + pb_x[k] * md_93[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_x, pb_x, pb_z, kf0_156, kf1_156, \
                         ld_95, lf_156, md_91, md_93, md_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_z[k] * md_91[k];

        t_155[k] = f_19 * ld_95[k]
                   + pb_x[k] * md_95[k];

        t_156[k] = f_17 * kf0_156[k]
                   - f_18 * kf1_156[k]
                   + pa_x[k] * lf_156[k];

        t_157[k] = pb_z[k] * md_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_z, pb_y, pb_z, ld_60, ld_65, \
                         lf_100, lf_101, mp0_47, mp1_47, md_95, md_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_20 * ld_65[k]
                   + pb_y[k] * md_95[k];

        t_159[k] = f_1 * mp0_47[k]
                   - f_2 * mp1_47[k]
                   + pb_z[k] * md_95[k];

        t_160[k] = pa_z[k] * lf_100[k];

        t_161[k] = pa_z[k] * lf_101[k];

        t_162[k] = f_3 * ld_60[k]
                   + pb_z[k] * md_96[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pa_z, pb_x, pb_z, ld_63, ld_100, \
                         ld_101, lf_103, lf_106, md_99, md_100, \
                         md_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_z[k] * lf_103[k];

        t_164[k] = f_19 * ld_100[k]
                   + pb_x[k] * md_100[k];

        t_165[k] = f_19 * ld_101[k]
                   + pb_x[k] * md_101[k];

        t_166[k] = pa_z[k] * lf_106[k];

        t_167[k] = f_3 * ld_63[k]
                   + pb_z[k] * md_99[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_y, pa_z, pb_y, kf0_80, kf1_80, ld_65, \
                         ld_71, ld_72, lf_109, lf_120, md_101, md_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_19 * ld_71[k]
                   + pb_y[k] * md_101[k];

        t_169[k] = f_5 * ld_65[k]
                   + pa_z[k] * lf_109[k];

        t_170[k] = f_12 * kf0_80[k]
                   - f_13 * kf1_80[k]
                   + pa_y[k] * lf_120[k];

        t_171[k] = f_5 * ld_72[k]
                   + pb_y[k] * md_102[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_x, pb_z, ld_66, ld_105, ld_106, \
                         ld_107, md_102, md_105, md_106, md_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_8 * ld_66[k]
                   + pb_z[k] * md_102[k];

        t_173[k] = f_19 * ld_105[k]
                   + pb_x[k] * md_105[k];

        t_174[k] = f_19 * ld_106[k]
                   + pb_x[k] * md_106[k];

        t_175[k] = f_19 * ld_107[k]
                   + pb_x[k] * md_107[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, pa_x, pb_y, pb_z, kf0_176, kf1_176, ld_69, \
                         ld_77, lf_176, md_105, md_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_17 * kf0_176[k]
                   - f_18 * kf1_176[k]
                   + pa_x[k] * lf_176[k];

        t_177[k] = f_8 * ld_69[k]
                   + pb_z[k] * md_105[k];

        t_178[k] = f_5 * ld_77[k]
                   + pb_y[k] * md_107[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pa_x, pa_y, pb_y, kf0_90, kf0_179, kf1_90, \
                         kf1_179, ld_78, lf_130, lf_179, md_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_17 * kf0_179[k]
                   - f_18 * kf1_179[k]
                   + pa_x[k] * lf_179[k];

        t_180[k] = f_6 * kf0_90[k]
                   - f_7 * kf1_90[k]
                   + pa_y[k] * lf_130[k];

        t_181[k] = f_8 * ld_78[k]
                   + pb_y[k] * md_108[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pb_x, pb_z, ld_72, ld_111, ld_112, \
                         ld_113, md_108, md_111, md_112, md_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_5 * ld_72[k]
                   + pb_z[k] * md_108[k];

        t_183[k] = f_19 * ld_111[k]
                   + pb_x[k] * md_111[k];

        t_184[k] = f_19 * ld_112[k]
                   + pb_x[k] * md_112[k];

        t_185[k] = f_19 * ld_113[k]
                   + pb_x[k] * md_113[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_x, pb_y, pb_z, kf0_186, kf1_186, ld_75, \
                         ld_83, lf_186, md_111, md_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_17 * kf0_186[k]
                   - f_18 * kf1_186[k]
                   + pa_x[k] * lf_186[k];

        t_187[k] = f_5 * ld_75[k]
                   + pb_z[k] * md_111[k];

        t_188[k] = f_8 * ld_83[k]
                   + pb_y[k] * md_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_x, pa_y, pb_y, kf0_189, kf1_189, \
                         ld_84, lf_140, lf_142, lf_189, md_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_17 * kf0_189[k]
                   - f_18 * kf1_189[k]
                   + pa_x[k] * lf_189[k];

        t_190[k] = pa_y[k] * lf_140[k];

        t_191[k] = f_3 * ld_84[k]
                   + pb_y[k] * md_114[k];

        t_192[k] = pa_y[k] * lf_142[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, pa_y, pb_x, pb_z, ld_81, ld_87, \
                         ld_117, ld_118, lf_145, lf_146, md_117, \
                         md_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_19 * ld_117[k]
                   + pb_x[k] * md_117[k];

        t_194[k] = f_19 * ld_118[k]
                   + pb_x[k] * md_118[k];

        t_195[k] = pa_y[k] * lf_145[k];

        t_196[k] = f_5 * ld_87[k]
                   + pa_y[k] * lf_146[k];

        t_197[k] = f_19 * ld_81[k]
                   + pb_z[k] * md_117[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pa_y, pa_z, pb_y, kf0_90, kf1_90, ld_89, \
                         lf_140, lf_149, md_119, md_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * ld_89[k]
                   + pb_y[k] * md_119[k];

        t_199[k] = pa_y[k] * lf_149[k];

        t_200[k] = f_21 * kf0_90[k]
                   - f_22 * kf1_90[k]
                   + pa_z[k] * lf_140[k];

        t_201[k] = pb_y[k] * md_120[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pb_x, pb_y, pb_z, ld_84, ld_123, ld_125, \
                         md_120, md_122, md_123, md_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_20 * ld_84[k]
                   + pb_z[k] * md_120[k];

        t_203[k] = f_19 * ld_123[k]
                   + pb_x[k] * md_123[k];

        t_204[k] = pb_y[k] * md_122[k];

        t_205[k] = f_19 * ld_125[k]
                   + pb_x[k] * md_125[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_x, pb_y, pb_z, kf0_209, kf1_209, \
                         ld_87, lf_209, mp0_61, mp1_61, md_123, \
                         md_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_1 * mp0_61[k]
                   - f_2 * mp1_61[k]
                   + pb_y[k] * md_123[k];

        t_207[k] = f_20 * ld_87[k]
                   + pb_z[k] * md_123[k];

        t_208[k] = pb_y[k] * md_125[k];

        t_209[k] = f_17 * kf0_209[k]
                   - f_18 * kf1_209[k]
                   + pa_x[k] * lf_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pa_y, pb_x, pb_y, pb_z, kf0_100, kf1_100, \
                         ld_90, ld_129, lf_150, md_126, md_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_15 * kf0_100[k]
                   - f_16 * kf1_100[k]
                   + pa_y[k] * lf_150[k];

        t_211[k] = f_14 * ld_90[k]
                   + pb_y[k] * md_126[k];

        t_212[k] = pb_z[k] * md_126[k];

        t_213[k] = f_5 * ld_129[k]
                   + pb_x[k] * md_129[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_x, pb_x, pb_z, kf0_216, kf1_216, \
                         ld_131, lf_216, md_127, md_129, md_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = pb_z[k] * md_127[k];

        t_215[k] = f_5 * ld_131[k]
                   + pb_x[k] * md_131[k];

        t_216[k] = f_12 * kf0_216[k]
                   - f_13 * kf1_216[k]
                   + pa_x[k] * lf_216[k];

        t_217[k] = pb_z[k] * md_129[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, pa_z, pb_y, pb_z, ld_90, ld_95, \
                         lf_150, lf_151, mp0_65, mp1_65, md_131, \
                         md_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_14 * ld_95[k]
                   + pb_y[k] * md_131[k];

        t_219[k] = f_1 * mp0_65[k]
                   - f_2 * mp1_65[k]
                   + pb_z[k] * md_131[k];

        t_220[k] = pa_z[k] * lf_150[k];

        t_221[k] = pa_z[k] * lf_151[k];

        t_222[k] = f_3 * ld_90[k]
                   + pb_z[k] * md_132[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, pa_z, pb_x, pb_z, ld_93, ld_136, \
                         ld_137, lf_153, lf_156, md_135, md_136, \
                         md_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = pa_z[k] * lf_153[k];

        t_224[k] = f_5 * ld_136[k]
                   + pb_x[k] * md_136[k];

        t_225[k] = f_5 * ld_137[k]
                   + pb_x[k] * md_137[k];

        t_226[k] = pa_z[k] * lf_156[k];

        t_227[k] = f_3 * ld_93[k]
                   + pb_z[k] * md_135[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pa_y, pa_z, pb_y, kf0_120, kf1_120, \
                         ld_95, ld_101, ld_102, lf_159, lf_170, md_137, \
                         md_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_20 * ld_101[k]
                   + pb_y[k] * md_137[k];

        t_229[k] = f_5 * ld_95[k]
                   + pa_z[k] * lf_159[k];

        t_230[k] = f_17 * kf0_120[k]
                   - f_18 * kf1_120[k]
                   + pa_y[k] * lf_170[k];

        t_231[k] = f_19 * ld_102[k]
                   + pb_y[k] * md_138[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pb_x, pb_z, ld_96, ld_141, ld_142, \
                         ld_143, md_138, md_141, md_142, md_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_8 * ld_96[k]
                   + pb_z[k] * md_138[k];

        t_233[k] = f_5 * ld_141[k]
                   + pb_x[k] * md_141[k];

        t_234[k] = f_5 * ld_142[k]
                   + pb_x[k] * md_142[k];

        t_235[k] = f_5 * ld_143[k]
                   + pb_x[k] * md_143[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pa_x, pb_y, pb_z, kf0_236, kf1_236, ld_99, \
                         ld_107, lf_236, md_141, md_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_12 * kf0_236[k]
                   - f_13 * kf1_236[k]
                   + pa_x[k] * lf_236[k];

        t_237[k] = f_8 * ld_99[k]
                   + pb_z[k] * md_141[k];

        t_238[k] = f_19 * ld_107[k]
                   + pb_y[k] * md_143[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pa_x, pa_y, pb_y, kf0_130, kf0_239, kf1_130, \
                         kf1_239, ld_108, lf_180, lf_239, md_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_12 * kf0_239[k]
                   - f_13 * kf1_239[k]
                   + pa_x[k] * lf_239[k];

        t_240[k] = f_12 * kf0_130[k]
                   - f_13 * kf1_130[k]
                   + pa_y[k] * lf_180[k];

        t_241[k] = f_5 * ld_108[k]
                   + pb_y[k] * md_144[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pb_x, pb_z, ld_102, ld_147, ld_148, \
                         ld_149, md_144, md_147, md_148, md_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_5 * ld_102[k]
                   + pb_z[k] * md_144[k];

        t_243[k] = f_5 * ld_147[k]
                   + pb_x[k] * md_147[k];

        t_244[k] = f_5 * ld_148[k]
                   + pb_x[k] * md_148[k];

        t_245[k] = f_5 * ld_149[k]
                   + pb_x[k] * md_149[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pa_x, pb_y, pb_z, kf0_246, kf1_246, ld_105, \
                         ld_113, lf_246, md_147, md_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * kf0_246[k]
                   - f_13 * kf1_246[k]
                   + pa_x[k] * lf_246[k];

        t_247[k] = f_5 * ld_105[k]
                   + pb_z[k] * md_147[k];

        t_248[k] = f_5 * ld_113[k]
                   + pb_y[k] * md_149[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pa_x, pa_y, pb_y, kf0_140, kf0_249, kf1_140, \
                         kf1_249, ld_114, lf_190, lf_249, md_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_12 * kf0_249[k]
                   - f_13 * kf1_249[k]
                   + pa_x[k] * lf_249[k];

        t_250[k] = f_6 * kf0_140[k]
                   - f_7 * kf1_140[k]
                   + pa_y[k] * lf_190[k];

        t_251[k] = f_8 * ld_114[k]
                   + pb_y[k] * md_150[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pb_x, pb_z, ld_108, ld_153, ld_154, \
                         ld_155, md_150, md_153, md_154, md_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_19 * ld_108[k]
                   + pb_z[k] * md_150[k];

        t_253[k] = f_5 * ld_153[k]
                   + pb_x[k] * md_153[k];

        t_254[k] = f_5 * ld_154[k]
                   + pb_x[k] * md_154[k];

        t_255[k] = f_5 * ld_155[k]
                   + pb_x[k] * md_155[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, pa_x, pb_y, pb_z, kf0_256, kf1_256, ld_111, \
                         ld_119, lf_256, md_153, md_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_12 * kf0_256[k]
                   - f_13 * kf1_256[k]
                   + pa_x[k] * lf_256[k];

        t_257[k] = f_19 * ld_111[k]
                   + pb_z[k] * md_153[k];

        t_258[k] = f_8 * ld_119[k]
                   + pb_y[k] * md_155[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, pa_x, pa_y, pb_y, kf0_259, kf1_259, \
                         ld_120, lf_200, lf_202, lf_259, md_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_12 * kf0_259[k]
                   - f_13 * kf1_259[k]
                   + pa_x[k] * lf_259[k];

        t_260[k] = pa_y[k] * lf_200[k];

        t_261[k] = f_3 * ld_120[k]
                   + pb_y[k] * md_156[k];

        t_262[k] = pa_y[k] * lf_202[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, pa_y, pb_x, pb_z, ld_117, ld_123, \
                         ld_159, ld_160, lf_205, lf_206, md_159, \
                         md_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_5 * ld_159[k]
                   + pb_x[k] * md_159[k];

        t_264[k] = f_5 * ld_160[k]
                   + pb_x[k] * md_160[k];

        t_265[k] = pa_y[k] * lf_205[k];

        t_266[k] = f_5 * ld_123[k]
                   + pa_y[k] * lf_206[k];

        t_267[k] = f_20 * ld_117[k]
                   + pb_z[k] * md_159[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_y, pa_z, pb_y, kf0_140, kf1_140, \
                         ld_125, lf_200, lf_209, md_161, md_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_3 * ld_125[k]
                   + pb_y[k] * md_161[k];

        t_269[k] = pa_y[k] * lf_209[k];

        t_270[k] = f_15 * kf0_140[k]
                   - f_16 * kf1_140[k]
                   + pa_z[k] * lf_200[k];

        t_271[k] = pb_y[k] * md_162[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pb_x, pb_y, pb_z, ld_120, ld_165, ld_167, \
                         md_162, md_164, md_165, md_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_14 * ld_120[k]
                   + pb_z[k] * md_162[k];

        t_273[k] = f_5 * ld_165[k]
                   + pb_x[k] * md_165[k];

        t_274[k] = pb_y[k] * md_164[k];

        t_275[k] = f_5 * ld_167[k]
                   + pb_x[k] * md_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pa_x, pb_y, pb_z, kf0_279, kf1_279, \
                         ld_123, lf_279, mp0_82, mp1_82, md_165, \
                         md_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_1 * mp0_82[k]
                   - f_2 * mp1_82[k]
                   + pb_y[k] * md_165[k];

        t_277[k] = f_14 * ld_123[k]
                   + pb_z[k] * md_165[k];

        t_278[k] = pb_y[k] * md_167[k];

        t_279[k] = f_12 * kf0_279[k]
                   - f_13 * kf1_279[k]
                   + pa_x[k] * lf_279[k];
    }
}

static auto
compute_prim_mf_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kf0,
                                            const size_t kf1, const size_t ld, const size_t lf,
                                            const size_t mp0, const size_t mp1, const size_t md,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 4.0 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 3.5 / p;
    const auto f_10 = 3.0 / alpha;
    const auto f_11 = 3.0 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 3.0 / p;
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / p;
    const auto f_20 = 2.5 / p;
    const auto f_21 = 2.0 / alpha;
    const auto f_22 = 2.0 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kf0_150 = buffer.data(kf0 + 150);
    const auto *kf0_170 = buffer.data(kf0 + 170);
    const auto *kf0_180 = buffer.data(kf0 + 180);
    const auto *kf0_190 = buffer.data(kf0 + 190);
    const auto *kf0_200 = buffer.data(kf0 + 200);
    const auto *kf0_286 = buffer.data(kf0 + 286);
    const auto *kf0_306 = buffer.data(kf0 + 306);
    const auto *kf0_309 = buffer.data(kf0 + 309);
    const auto *kf0_316 = buffer.data(kf0 + 316);
    const auto *kf0_319 = buffer.data(kf0 + 319);
    const auto *kf0_326 = buffer.data(kf0 + 326);
    const auto *kf0_329 = buffer.data(kf0 + 329);
    const auto *kf0_336 = buffer.data(kf0 + 336);
    const auto *kf0_339 = buffer.data(kf0 + 339);
    const auto *kf0_359 = buffer.data(kf0 + 359);

    const auto *kf1_150 = buffer.data(kf1 + 150);
    const auto *kf1_170 = buffer.data(kf1 + 170);
    const auto *kf1_180 = buffer.data(kf1 + 180);
    const auto *kf1_190 = buffer.data(kf1 + 190);
    const auto *kf1_200 = buffer.data(kf1 + 200);
    const auto *kf1_286 = buffer.data(kf1 + 286);
    const auto *kf1_306 = buffer.data(kf1 + 306);
    const auto *kf1_309 = buffer.data(kf1 + 309);
    const auto *kf1_316 = buffer.data(kf1 + 316);
    const auto *kf1_319 = buffer.data(kf1 + 319);
    const auto *kf1_326 = buffer.data(kf1 + 326);
    const auto *kf1_329 = buffer.data(kf1 + 329);
    const auto *kf1_336 = buffer.data(kf1 + 336);
    const auto *kf1_339 = buffer.data(kf1 + 339);
    const auto *kf1_359 = buffer.data(kf1 + 359);

    const auto *ld_126 = buffer.data(ld + 126);
    const auto *ld_129 = buffer.data(ld + 129);
    const auto *ld_131 = buffer.data(ld + 131);
    const auto *ld_132 = buffer.data(ld + 132);
    const auto *ld_135 = buffer.data(ld + 135);
    const auto *ld_137 = buffer.data(ld + 137);
    const auto *ld_138 = buffer.data(ld + 138);
    const auto *ld_141 = buffer.data(ld + 141);
    const auto *ld_143 = buffer.data(ld + 143);
    const auto *ld_144 = buffer.data(ld + 144);
    const auto *ld_147 = buffer.data(ld + 147);
    const auto *ld_149 = buffer.data(ld + 149);
    const auto *ld_150 = buffer.data(ld + 150);
    const auto *ld_153 = buffer.data(ld + 153);
    const auto *ld_155 = buffer.data(ld + 155);
    const auto *ld_156 = buffer.data(ld + 156);
    const auto *ld_159 = buffer.data(ld + 159);
    const auto *ld_161 = buffer.data(ld + 161);
    const auto *ld_162 = buffer.data(ld + 162);
    const auto *ld_165 = buffer.data(ld + 165);
    const auto *ld_167 = buffer.data(ld + 167);
    const auto *ld_168 = buffer.data(ld + 168);
    const auto *ld_171 = buffer.data(ld + 171);
    const auto *ld_173 = buffer.data(ld + 173);
    const auto *ld_174 = buffer.data(ld + 174);
    const auto *ld_178 = buffer.data(ld + 178);
    const auto *ld_179 = buffer.data(ld + 179);
    const auto *ld_180 = buffer.data(ld + 180);
    const auto *ld_183 = buffer.data(ld + 183);
    const auto *ld_184 = buffer.data(ld + 184);
    const auto *ld_185 = buffer.data(ld + 185);
    const auto *ld_186 = buffer.data(ld + 186);
    const auto *ld_189 = buffer.data(ld + 189);
    const auto *ld_190 = buffer.data(ld + 190);
    const auto *ld_191 = buffer.data(ld + 191);
    const auto *ld_192 = buffer.data(ld + 192);
    const auto *ld_195 = buffer.data(ld + 195);
    const auto *ld_196 = buffer.data(ld + 196);
    const auto *ld_197 = buffer.data(ld + 197);
    const auto *ld_198 = buffer.data(ld + 198);
    const auto *ld_201 = buffer.data(ld + 201);
    const auto *ld_202 = buffer.data(ld + 202);
    const auto *ld_203 = buffer.data(ld + 203);
    const auto *ld_207 = buffer.data(ld + 207);
    const auto *ld_208 = buffer.data(ld + 208);
    const auto *ld_213 = buffer.data(ld + 213);
    const auto *ld_215 = buffer.data(ld + 215);
    const auto *ld_216 = buffer.data(ld + 216);
    const auto *ld_219 = buffer.data(ld + 219);
    const auto *ld_221 = buffer.data(ld + 221);
    const auto *ld_226 = buffer.data(ld + 226);
    const auto *ld_227 = buffer.data(ld + 227);
    const auto *ld_228 = buffer.data(ld + 228);
    const auto *ld_231 = buffer.data(ld + 231);
    const auto *ld_232 = buffer.data(ld + 232);
    const auto *ld_233 = buffer.data(ld + 233);
    const auto *ld_234 = buffer.data(ld + 234);
    const auto *ld_237 = buffer.data(ld + 237);
    const auto *ld_238 = buffer.data(ld + 238);
    const auto *ld_239 = buffer.data(ld + 239);
    const auto *ld_240 = buffer.data(ld + 240);
    const auto *ld_243 = buffer.data(ld + 243);
    const auto *ld_244 = buffer.data(ld + 244);
    const auto *ld_245 = buffer.data(ld + 245);
    const auto *ld_246 = buffer.data(ld + 246);
    const auto *ld_249 = buffer.data(ld + 249);
    const auto *ld_250 = buffer.data(ld + 250);
    const auto *ld_251 = buffer.data(ld + 251);
    const auto *ld_252 = buffer.data(ld + 252);

    const auto *lf_210 = buffer.data(lf + 210);
    const auto *lf_211 = buffer.data(lf + 211);
    const auto *lf_213 = buffer.data(lf + 213);
    const auto *lf_216 = buffer.data(lf + 216);
    const auto *lf_219 = buffer.data(lf + 219);
    const auto *lf_230 = buffer.data(lf + 230);
    const auto *lf_240 = buffer.data(lf + 240);
    const auto *lf_250 = buffer.data(lf + 250);
    const auto *lf_260 = buffer.data(lf + 260);
    const auto *lf_270 = buffer.data(lf + 270);
    const auto *lf_272 = buffer.data(lf + 272);
    const auto *lf_275 = buffer.data(lf + 275);
    const auto *lf_276 = buffer.data(lf + 276);
    const auto *lf_279 = buffer.data(lf + 279);
    const auto *lf_280 = buffer.data(lf + 280);
    const auto *lf_281 = buffer.data(lf + 281);
    const auto *lf_283 = buffer.data(lf + 283);
    const auto *lf_286 = buffer.data(lf + 286);
    const auto *lf_306 = buffer.data(lf + 306);
    const auto *lf_309 = buffer.data(lf + 309);
    const auto *lf_316 = buffer.data(lf + 316);
    const auto *lf_319 = buffer.data(lf + 319);
    const auto *lf_326 = buffer.data(lf + 326);
    const auto *lf_329 = buffer.data(lf + 329);
    const auto *lf_336 = buffer.data(lf + 336);
    const auto *lf_339 = buffer.data(lf + 339);
    const auto *lf_359 = buffer.data(lf + 359);
    const auto *lf_360 = buffer.data(lf + 360);
    const auto *lf_366 = buffer.data(lf + 366);
    const auto *lf_368 = buffer.data(lf + 368);
    const auto *lf_369 = buffer.data(lf + 369);
    const auto *lf_376 = buffer.data(lf + 376);
    const auto *lf_377 = buffer.data(lf + 377);
    const auto *lf_378 = buffer.data(lf + 378);
    const auto *lf_379 = buffer.data(lf + 379);
    const auto *lf_380 = buffer.data(lf + 380);
    const auto *lf_386 = buffer.data(lf + 386);
    const auto *lf_387 = buffer.data(lf + 387);
    const auto *lf_388 = buffer.data(lf + 388);
    const auto *lf_389 = buffer.data(lf + 389);
    const auto *lf_390 = buffer.data(lf + 390);
    const auto *lf_396 = buffer.data(lf + 396);
    const auto *lf_397 = buffer.data(lf + 397);
    const auto *lf_398 = buffer.data(lf + 398);
    const auto *lf_399 = buffer.data(lf + 399);
    const auto *lf_400 = buffer.data(lf + 400);
    const auto *lf_406 = buffer.data(lf + 406);
    const auto *lf_407 = buffer.data(lf + 407);
    const auto *lf_408 = buffer.data(lf + 408);
    const auto *lf_409 = buffer.data(lf + 409);
    const auto *lf_410 = buffer.data(lf + 410);
    const auto *lf_416 = buffer.data(lf + 416);
    const auto *lf_417 = buffer.data(lf + 417);
    const auto *lf_418 = buffer.data(lf + 418);
    const auto *lf_419 = buffer.data(lf + 419);
    const auto *lf_420 = buffer.data(lf + 420);

    const auto *mp0_86 = buffer.data(mp0 + 86);
    const auto *mp0_106 = buffer.data(mp0 + 106);

    const auto *mp1_86 = buffer.data(mp1 + 86);
    const auto *mp1_106 = buffer.data(mp1 + 106);

    const auto *md_168 = buffer.data(md + 168);
    const auto *md_169 = buffer.data(md + 169);
    const auto *md_171 = buffer.data(md + 171);
    const auto *md_173 = buffer.data(md + 173);
    const auto *md_174 = buffer.data(md + 174);
    const auto *md_177 = buffer.data(md + 177);
    const auto *md_178 = buffer.data(md + 178);
    const auto *md_179 = buffer.data(md + 179);
    const auto *md_180 = buffer.data(md + 180);
    const auto *md_183 = buffer.data(md + 183);
    const auto *md_184 = buffer.data(md + 184);
    const auto *md_185 = buffer.data(md + 185);
    const auto *md_186 = buffer.data(md + 186);
    const auto *md_189 = buffer.data(md + 189);
    const auto *md_190 = buffer.data(md + 190);
    const auto *md_191 = buffer.data(md + 191);
    const auto *md_192 = buffer.data(md + 192);
    const auto *md_195 = buffer.data(md + 195);
    const auto *md_196 = buffer.data(md + 196);
    const auto *md_197 = buffer.data(md + 197);
    const auto *md_198 = buffer.data(md + 198);
    const auto *md_201 = buffer.data(md + 201);
    const auto *md_202 = buffer.data(md + 202);
    const auto *md_203 = buffer.data(md + 203);
    const auto *md_204 = buffer.data(md + 204);
    const auto *md_207 = buffer.data(md + 207);
    const auto *md_208 = buffer.data(md + 208);
    const auto *md_209 = buffer.data(md + 209);
    const auto *md_210 = buffer.data(md + 210);
    const auto *md_212 = buffer.data(md + 212);
    const auto *md_213 = buffer.data(md + 213);
    const auto *md_215 = buffer.data(md + 215);
    const auto *md_216 = buffer.data(md + 216);
    const auto *md_217 = buffer.data(md + 217);
    const auto *md_219 = buffer.data(md + 219);
    const auto *md_221 = buffer.data(md + 221);
    const auto *md_222 = buffer.data(md + 222);
    const auto *md_226 = buffer.data(md + 226);
    const auto *md_227 = buffer.data(md + 227);
    const auto *md_228 = buffer.data(md + 228);
    const auto *md_231 = buffer.data(md + 231);
    const auto *md_232 = buffer.data(md + 232);
    const auto *md_233 = buffer.data(md + 233);
    const auto *md_234 = buffer.data(md + 234);
    const auto *md_237 = buffer.data(md + 237);
    const auto *md_238 = buffer.data(md + 238);
    const auto *md_239 = buffer.data(md + 239);
    const auto *md_240 = buffer.data(md + 240);
    const auto *md_243 = buffer.data(md + 243);
    const auto *md_244 = buffer.data(md + 244);
    const auto *md_245 = buffer.data(md + 245);
    const auto *md_246 = buffer.data(md + 246);
    const auto *md_249 = buffer.data(md + 249);
    const auto *md_250 = buffer.data(md + 250);
    const auto *md_251 = buffer.data(md + 251);

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pa_y, pb_x, pb_y, pb_z, kf0_150, kf1_150, \
                         ld_126, ld_171, lf_210, md_168, md_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_10 * kf0_150[k]
                   - f_11 * kf1_150[k]
                   + pa_y[k] * lf_210[k];

        t_281[k] = f_9 * ld_126[k]
                   + pb_y[k] * md_168[k];

        t_282[k] = pb_z[k] * md_168[k];

        t_283[k] = f_8 * ld_171[k]
                   + pb_x[k] * md_171[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pa_x, pb_x, pb_z, kf0_286, kf1_286, \
                         ld_173, lf_286, md_169, md_171, md_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = pb_z[k] * md_169[k];

        t_285[k] = f_8 * ld_173[k]
                   + pb_x[k] * md_173[k];

        t_286[k] = f_6 * kf0_286[k]
                   - f_7 * kf1_286[k]
                   + pa_x[k] * lf_286[k];

        t_287[k] = pb_z[k] * md_171[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, pa_z, pb_y, pb_z, ld_126, ld_131, \
                         lf_210, lf_211, mp0_86, mp1_86, md_173, \
                         md_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * ld_131[k]
                   + pb_y[k] * md_173[k];

        t_289[k] = f_1 * mp0_86[k]
                   - f_2 * mp1_86[k]
                   + pb_z[k] * md_173[k];

        t_290[k] = pa_z[k] * lf_210[k];

        t_291[k] = pa_z[k] * lf_211[k];

        t_292[k] = f_3 * ld_126[k]
                   + pb_z[k] * md_174[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, t_297, pa_z, pb_x, pb_z, ld_129, ld_178, \
                         ld_179, lf_213, lf_216, md_177, md_178, \
                         md_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = pa_z[k] * lf_213[k];

        t_294[k] = f_8 * ld_178[k]
                   + pb_x[k] * md_178[k];

        t_295[k] = f_8 * ld_179[k]
                   + pb_x[k] * md_179[k];

        t_296[k] = pa_z[k] * lf_216[k];

        t_297[k] = f_3 * ld_129[k]
                   + pb_z[k] * md_177[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pa_y, pa_z, pb_y, kf0_170, kf1_170, \
                         ld_131, ld_137, ld_138, lf_219, lf_230, md_179, \
                         md_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_14 * ld_137[k]
                   + pb_y[k] * md_179[k];

        t_299[k] = f_5 * ld_131[k]
                   + pa_z[k] * lf_219[k];

        t_300[k] = f_21 * kf0_170[k]
                   - f_22 * kf1_170[k]
                   + pa_y[k] * lf_230[k];

        t_301[k] = f_20 * ld_138[k]
                   + pb_y[k] * md_180[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pb_x, pb_z, ld_132, ld_183, ld_184, \
                         ld_185, md_180, md_183, md_184, md_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_8 * ld_132[k]
                   + pb_z[k] * md_180[k];

        t_303[k] = f_8 * ld_183[k]
                   + pb_x[k] * md_183[k];

        t_304[k] = f_8 * ld_184[k]
                   + pb_x[k] * md_184[k];

        t_305[k] = f_8 * ld_185[k]
                   + pb_x[k] * md_185[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, pa_x, pb_y, pb_z, kf0_306, kf1_306, ld_135, \
                         ld_143, lf_306, md_183, md_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_6 * kf0_306[k]
                   - f_7 * kf1_306[k]
                   + pa_x[k] * lf_306[k];

        t_307[k] = f_8 * ld_135[k]
                   + pb_z[k] * md_183[k];

        t_308[k] = f_20 * ld_143[k]
                   + pb_y[k] * md_185[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pa_x, pa_y, pb_y, kf0_180, kf0_309, kf1_180, \
                         kf1_309, ld_144, lf_240, lf_309, md_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_6 * kf0_309[k]
                   - f_7 * kf1_309[k]
                   + pa_x[k] * lf_309[k];

        t_310[k] = f_17 * kf0_180[k]
                   - f_18 * kf1_180[k]
                   + pa_y[k] * lf_240[k];

        t_311[k] = f_19 * ld_144[k]
                   + pb_y[k] * md_186[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pb_x, pb_z, ld_138, ld_189, ld_190, \
                         ld_191, md_186, md_189, md_190, md_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_5 * ld_138[k]
                   + pb_z[k] * md_186[k];

        t_313[k] = f_8 * ld_189[k]
                   + pb_x[k] * md_189[k];

        t_314[k] = f_8 * ld_190[k]
                   + pb_x[k] * md_190[k];

        t_315[k] = f_8 * ld_191[k]
                   + pb_x[k] * md_191[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pa_x, pb_y, pb_z, kf0_316, kf1_316, ld_141, \
                         ld_149, lf_316, md_189, md_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_6 * kf0_316[k]
                   - f_7 * kf1_316[k]
                   + pa_x[k] * lf_316[k];

        t_317[k] = f_5 * ld_141[k]
                   + pb_z[k] * md_189[k];

        t_318[k] = f_19 * ld_149[k]
                   + pb_y[k] * md_191[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_x, pa_y, pb_y, kf0_190, kf0_319, kf1_190, \
                         kf1_319, ld_150, lf_250, lf_319, md_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_6 * kf0_319[k]
                   - f_7 * kf1_319[k]
                   + pa_x[k] * lf_319[k];

        t_320[k] = f_12 * kf0_190[k]
                   - f_13 * kf1_190[k]
                   + pa_y[k] * lf_250[k];

        t_321[k] = f_5 * ld_150[k]
                   + pb_y[k] * md_192[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pb_x, pb_z, ld_144, ld_195, ld_196, \
                         ld_197, md_192, md_195, md_196, md_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_19 * ld_144[k]
                   + pb_z[k] * md_192[k];

        t_323[k] = f_8 * ld_195[k]
                   + pb_x[k] * md_195[k];

        t_324[k] = f_8 * ld_196[k]
                   + pb_x[k] * md_196[k];

        t_325[k] = f_8 * ld_197[k]
                   + pb_x[k] * md_197[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, pa_x, pb_y, pb_z, kf0_326, kf1_326, ld_147, \
                         ld_155, lf_326, md_195, md_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_6 * kf0_326[k]
                   - f_7 * kf1_326[k]
                   + pa_x[k] * lf_326[k];

        t_327[k] = f_19 * ld_147[k]
                   + pb_z[k] * md_195[k];

        t_328[k] = f_5 * ld_155[k]
                   + pb_y[k] * md_197[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pa_x, pa_y, pb_y, kf0_200, kf0_329, kf1_200, \
                         kf1_329, ld_156, lf_260, lf_329, md_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_6 * kf0_329[k]
                   - f_7 * kf1_329[k]
                   + pa_x[k] * lf_329[k];

        t_330[k] = f_6 * kf0_200[k]
                   - f_7 * kf1_200[k]
                   + pa_y[k] * lf_260[k];

        t_331[k] = f_8 * ld_156[k]
                   + pb_y[k] * md_198[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pb_x, pb_z, ld_150, ld_201, ld_202, \
                         ld_203, md_198, md_201, md_202, md_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_20 * ld_150[k]
                   + pb_z[k] * md_198[k];

        t_333[k] = f_8 * ld_201[k]
                   + pb_x[k] * md_201[k];

        t_334[k] = f_8 * ld_202[k]
                   + pb_x[k] * md_202[k];

        t_335[k] = f_8 * ld_203[k]
                   + pb_x[k] * md_203[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pa_x, pb_y, pb_z, kf0_336, kf1_336, ld_153, \
                         ld_161, lf_336, md_201, md_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_6 * kf0_336[k]
                   - f_7 * kf1_336[k]
                   + pa_x[k] * lf_336[k];

        t_337[k] = f_20 * ld_153[k]
                   + pb_z[k] * md_201[k];

        t_338[k] = f_8 * ld_161[k]
                   + pb_y[k] * md_203[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, pa_x, pa_y, pb_y, kf0_339, kf1_339, \
                         ld_162, lf_270, lf_272, lf_339, md_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_6 * kf0_339[k]
                   - f_7 * kf1_339[k]
                   + pa_x[k] * lf_339[k];

        t_340[k] = pa_y[k] * lf_270[k];

        t_341[k] = f_3 * ld_162[k]
                   + pb_y[k] * md_204[k];

        t_342[k] = pa_y[k] * lf_272[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, pa_y, pb_x, pb_z, ld_159, ld_165, \
                         ld_207, ld_208, lf_275, lf_276, md_207, \
                         md_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_8 * ld_207[k]
                   + pb_x[k] * md_207[k];

        t_344[k] = f_8 * ld_208[k]
                   + pb_x[k] * md_208[k];

        t_345[k] = pa_y[k] * lf_275[k];

        t_346[k] = f_5 * ld_165[k]
                   + pa_y[k] * lf_276[k];

        t_347[k] = f_14 * ld_159[k]
                   + pb_z[k] * md_207[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pa_y, pa_z, pb_y, kf0_200, kf1_200, \
                         ld_167, lf_270, lf_279, md_209, md_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_3 * ld_167[k]
                   + pb_y[k] * md_209[k];

        t_349[k] = pa_y[k] * lf_279[k];

        t_350[k] = f_10 * kf0_200[k]
                   - f_11 * kf1_200[k]
                   + pa_z[k] * lf_270[k];

        t_351[k] = pb_y[k] * md_210[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pb_x, pb_y, pb_z, ld_162, ld_213, ld_215, \
                         md_210, md_212, md_213, md_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_9 * ld_162[k]
                   + pb_z[k] * md_210[k];

        t_353[k] = f_8 * ld_213[k]
                   + pb_x[k] * md_213[k];

        t_354[k] = pb_y[k] * md_212[k];

        t_355[k] = f_8 * ld_215[k]
                   + pb_x[k] * md_215[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pa_x, pb_y, pb_z, kf0_359, kf1_359, \
                         ld_165, lf_359, mp0_106, mp1_106, md_213, \
                         md_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_1 * mp0_106[k]
                   - f_2 * mp1_106[k]
                   + pb_y[k] * md_213[k];

        t_357[k] = f_9 * ld_165[k]
                   + pb_z[k] * md_213[k];

        t_358[k] = pb_y[k] * md_215[k];

        t_359[k] = f_6 * kf0_359[k]
                   - f_7 * kf1_359[k]
                   + pa_x[k] * lf_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, pa_x, pb_x, pb_y, pb_z, ld_168, \
                         ld_216, ld_219, lf_360, md_216, md_217, \
                         md_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_5 * ld_216[k]
                   + pa_x[k] * lf_360[k];

        t_361[k] = f_4 * ld_168[k]
                   + pb_y[k] * md_216[k];

        t_362[k] = pb_z[k] * md_216[k];

        t_363[k] = f_3 * ld_219[k]
                   + pb_x[k] * md_219[k];

        t_364[k] = pb_z[k] * md_217[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, pa_x, pb_x, pb_z, ld_221, lf_366, \
                         lf_368, lf_369, md_219, md_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_3 * ld_221[k]
                   + pb_x[k] * md_221[k];

        t_366[k] = pa_x[k] * lf_366[k];

        t_367[k] = pb_z[k] * md_219[k];

        t_368[k] = pa_x[k] * lf_368[k];

        t_369[k] = pa_x[k] * lf_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, pa_z, pb_x, pb_z, ld_168, ld_226, \
                         lf_280, lf_281, lf_283, md_222, md_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pa_z[k] * lf_280[k];

        t_371[k] = pa_z[k] * lf_281[k];

        t_372[k] = f_3 * ld_168[k]
                   + pb_z[k] * md_222[k];

        t_373[k] = pa_z[k] * lf_283[k];

        t_374[k] = f_3 * ld_226[k]
                   + pb_x[k] * md_226[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, t_380, pa_x, pb_x, ld_227, ld_228, \
                         lf_376, lf_377, lf_378, lf_379, lf_380, \
                         md_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_3 * ld_227[k]
                   + pb_x[k] * md_227[k];

        t_376[k] = pa_x[k] * lf_376[k];

        t_377[k] = pa_x[k] * lf_377[k];

        t_378[k] = pa_x[k] * lf_378[k];

        t_379[k] = pa_x[k] * lf_379[k];

        t_380[k] = f_5 * ld_228[k]
                   + pa_x[k] * lf_380[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pb_x, pb_y, pb_z, ld_174, ld_180, ld_231, \
                         ld_232, md_228, md_231, md_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_14 * ld_180[k]
                   + pb_y[k] * md_228[k];

        t_382[k] = f_8 * ld_174[k]
                   + pb_z[k] * md_228[k];

        t_383[k] = f_3 * ld_231[k]
                   + pb_x[k] * md_231[k];

        t_384[k] = f_3 * ld_232[k]
                   + pb_x[k] * md_232[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, t_390, pa_x, pb_x, ld_233, ld_234, \
                         lf_386, lf_387, lf_388, lf_389, lf_390, \
                         md_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_3 * ld_233[k]
                   + pb_x[k] * md_233[k];

        t_386[k] = pa_x[k] * lf_386[k];

        t_387[k] = pa_x[k] * lf_387[k];

        t_388[k] = pa_x[k] * lf_388[k];

        t_389[k] = pa_x[k] * lf_389[k];

        t_390[k] = f_5 * ld_234[k]
                   + pa_x[k] * lf_390[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pb_x, pb_y, pb_z, ld_180, ld_186, ld_237, \
                         ld_238, md_234, md_237, md_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_20 * ld_186[k]
                   + pb_y[k] * md_234[k];

        t_392[k] = f_5 * ld_180[k]
                   + pb_z[k] * md_234[k];

        t_393[k] = f_3 * ld_237[k]
                   + pb_x[k] * md_237[k];

        t_394[k] = f_3 * ld_238[k]
                   + pb_x[k] * md_238[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, t_400, pa_x, pb_x, ld_239, ld_240, \
                         lf_396, lf_397, lf_398, lf_399, lf_400, \
                         md_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_3 * ld_239[k]
                   + pb_x[k] * md_239[k];

        t_396[k] = pa_x[k] * lf_396[k];

        t_397[k] = pa_x[k] * lf_397[k];

        t_398[k] = pa_x[k] * lf_398[k];

        t_399[k] = pa_x[k] * lf_399[k];

        t_400[k] = f_5 * ld_240[k]
                   + pa_x[k] * lf_400[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pb_x, pb_y, pb_z, ld_186, ld_192, ld_243, \
                         ld_244, md_240, md_243, md_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_19 * ld_192[k]
                   + pb_y[k] * md_240[k];

        t_402[k] = f_19 * ld_186[k]
                   + pb_z[k] * md_240[k];

        t_403[k] = f_3 * ld_243[k]
                   + pb_x[k] * md_243[k];

        t_404[k] = f_3 * ld_244[k]
                   + pb_x[k] * md_244[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, t_410, pa_x, pb_x, ld_245, ld_246, \
                         lf_406, lf_407, lf_408, lf_409, lf_410, \
                         md_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_3 * ld_245[k]
                   + pb_x[k] * md_245[k];

        t_406[k] = pa_x[k] * lf_406[k];

        t_407[k] = pa_x[k] * lf_407[k];

        t_408[k] = pa_x[k] * lf_408[k];

        t_409[k] = pa_x[k] * lf_409[k];

        t_410[k] = f_5 * ld_246[k]
                   + pa_x[k] * lf_410[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, pb_x, pb_y, pb_z, ld_192, ld_198, ld_249, \
                         ld_250, md_246, md_249, md_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_5 * ld_198[k]
                   + pb_y[k] * md_246[k];

        t_412[k] = f_20 * ld_192[k]
                   + pb_z[k] * md_246[k];

        t_413[k] = f_3 * ld_249[k]
                   + pb_x[k] * md_249[k];

        t_414[k] = f_3 * ld_250[k]
                   + pb_x[k] * md_250[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, t_420, pa_x, pb_x, ld_251, ld_252, \
                         lf_416, lf_417, lf_418, lf_419, lf_420, \
                         md_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_3 * ld_251[k]
                   + pb_x[k] * md_251[k];

        t_416[k] = pa_x[k] * lf_416[k];

        t_417[k] = pa_x[k] * lf_417[k];

        t_418[k] = pa_x[k] * lf_418[k];

        t_419[k] = pa_x[k] * lf_419[k];

        t_420[k] = f_5 * ld_252[k]
                   + pa_x[k] * lf_420[k];
    }
}

static auto
compute_prim_mf_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kf0,
                                            const size_t kf1, const size_t ld, const size_t lf,
                                            const size_t mp0, const size_t mp1, const size_t md,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 4.0 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 3.5 / p;
    const auto f_10 = 3.0 / alpha;
    const auto f_11 = 3.0 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 3.0 / p;
    const auto f_15 = 2.5 / alpha;
    const auto f_16 = 2.5 * beta / (alpha * p);
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / p;
    const auto f_20 = 2.5 / p;
    const auto f_21 = 2.0 / alpha;
    const auto f_22 = 2.0 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kf0_286 = buffer.data(kf0 + 286);
    const auto *kf0_296 = buffer.data(kf0 + 296);
    const auto *kf0_306 = buffer.data(kf0 + 306);
    const auto *kf0_309 = buffer.data(kf0 + 309);
    const auto *kf0_316 = buffer.data(kf0 + 316);
    const auto *kf0_319 = buffer.data(kf0 + 319);
    const auto *kf0_326 = buffer.data(kf0 + 326);
    const auto *kf0_329 = buffer.data(kf0 + 329);
    const auto *kf0_336 = buffer.data(kf0 + 336);
    const auto *kf0_339 = buffer.data(kf0 + 339);
    const auto *kf0_349 = buffer.data(kf0 + 349);
    const auto *kf0_359 = buffer.data(kf0 + 359);

    const auto *kf1_286 = buffer.data(kf1 + 286);
    const auto *kf1_296 = buffer.data(kf1 + 296);
    const auto *kf1_306 = buffer.data(kf1 + 306);
    const auto *kf1_309 = buffer.data(kf1 + 309);
    const auto *kf1_316 = buffer.data(kf1 + 316);
    const auto *kf1_319 = buffer.data(kf1 + 319);
    const auto *kf1_326 = buffer.data(kf1 + 326);
    const auto *kf1_329 = buffer.data(kf1 + 329);
    const auto *kf1_336 = buffer.data(kf1 + 336);
    const auto *kf1_339 = buffer.data(kf1 + 339);
    const auto *kf1_349 = buffer.data(kf1 + 349);
    const auto *kf1_359 = buffer.data(kf1 + 359);

    const auto *ld_198 = buffer.data(ld + 198);
    const auto *ld_204 = buffer.data(ld + 204);
    const auto *ld_210 = buffer.data(ld + 210);
    const auto *ld_216 = buffer.data(ld + 216);
    const auto *ld_219 = buffer.data(ld + 219);
    const auto *ld_221 = buffer.data(ld + 221);
    const auto *ld_222 = buffer.data(ld + 222);
    const auto *ld_225 = buffer.data(ld + 225);
    const auto *ld_227 = buffer.data(ld + 227);
    const auto *ld_228 = buffer.data(ld + 228);
    const auto *ld_231 = buffer.data(ld + 231);
    const auto *ld_233 = buffer.data(ld + 233);
    const auto *ld_234 = buffer.data(ld + 234);
    const auto *ld_237 = buffer.data(ld + 237);
    const auto *ld_239 = buffer.data(ld + 239);
    const auto *ld_240 = buffer.data(ld + 240);
    const auto *ld_243 = buffer.data(ld + 243);
    const auto *ld_245 = buffer.data(ld + 245);
    const auto *ld_246 = buffer.data(ld + 246);
    const auto *ld_249 = buffer.data(ld + 249);
    const auto *ld_251 = buffer.data(ld + 251);
    const auto *ld_252 = buffer.data(ld + 252);
    const auto *ld_255 = buffer.data(ld + 255);
    const auto *ld_256 = buffer.data(ld + 256);
    const auto *ld_257 = buffer.data(ld + 257);
    const auto *ld_258 = buffer.data(ld + 258);
    const auto *ld_261 = buffer.data(ld + 261);
    const auto *ld_262 = buffer.data(ld + 262);
    const auto *ld_263 = buffer.data(ld + 263);
    const auto *ld_264 = buffer.data(ld + 264);
    const auto *ld_267 = buffer.data(ld + 267);
    const auto *ld_269 = buffer.data(ld + 269);

    const auto *lf_350 = buffer.data(lf + 350);
    const auto *lf_352 = buffer.data(lf + 352);
    const auto *lf_355 = buffer.data(lf + 355);
    const auto *lf_360 = buffer.data(lf + 360);
    const auto *lf_361 = buffer.data(lf + 361);
    const auto *lf_366 = buffer.data(lf + 366);
    const auto *lf_369 = buffer.data(lf + 369);
    const auto *lf_376 = buffer.data(lf + 376);
    const auto *lf_386 = buffer.data(lf + 386);
    const auto *lf_389 = buffer.data(lf + 389);
    const auto *lf_396 = buffer.data(lf + 396);
    const auto *lf_399 = buffer.data(lf + 399);
    const auto *lf_406 = buffer.data(lf + 406);
    const auto *lf_409 = buffer.data(lf + 409);
    const auto *lf_416 = buffer.data(lf + 416);
    const auto *lf_419 = buffer.data(lf + 419);
    const auto *lf_426 = buffer.data(lf + 426);
    const auto *lf_427 = buffer.data(lf + 427);
    const auto *lf_428 = buffer.data(lf + 428);
    const auto *lf_429 = buffer.data(lf + 429);
    const auto *lf_436 = buffer.data(lf + 436);
    const auto *lf_437 = buffer.data(lf + 437);
    const auto *lf_438 = buffer.data(lf + 438);
    const auto *lf_439 = buffer.data(lf + 439);
    const auto *lf_440 = buffer.data(lf + 440);
    const auto *lf_442 = buffer.data(lf + 442);
    const auto *lf_446 = buffer.data(lf + 446);
    const auto *lf_447 = buffer.data(lf + 447);
    const auto *lf_449 = buffer.data(lf + 449);

    const auto *mp0_135 = buffer.data(mp0 + 135);
    const auto *mp0_136 = buffer.data(mp0 + 136);
    const auto *mp0_137 = buffer.data(mp0 + 137);
    const auto *mp0_141 = buffer.data(mp0 + 141);
    const auto *mp0_144 = buffer.data(mp0 + 144);
    const auto *mp0_147 = buffer.data(mp0 + 147);
    const auto *mp0_150 = buffer.data(mp0 + 150);
    const auto *mp0_153 = buffer.data(mp0 + 153);
    const auto *mp0_156 = buffer.data(mp0 + 156);
    const auto *mp0_162 = buffer.data(mp0 + 162);
    const auto *mp0_163 = buffer.data(mp0 + 163);
    const auto *mp0_164 = buffer.data(mp0 + 164);

    const auto *mp1_135 = buffer.data(mp1 + 135);
    const auto *mp1_136 = buffer.data(mp1 + 136);
    const auto *mp1_137 = buffer.data(mp1 + 137);
    const auto *mp1_141 = buffer.data(mp1 + 141);
    const auto *mp1_144 = buffer.data(mp1 + 144);
    const auto *mp1_147 = buffer.data(mp1 + 147);
    const auto *mp1_150 = buffer.data(mp1 + 150);
    const auto *mp1_153 = buffer.data(mp1 + 153);
    const auto *mp1_156 = buffer.data(mp1 + 156);
    const auto *mp1_162 = buffer.data(mp1 + 162);
    const auto *mp1_163 = buffer.data(mp1 + 163);
    const auto *mp1_164 = buffer.data(mp1 + 164);

    const auto *md_252 = buffer.data(md + 252);
    const auto *md_255 = buffer.data(md + 255);
    const auto *md_256 = buffer.data(md + 256);
    const auto *md_257 = buffer.data(md + 257);
    const auto *md_258 = buffer.data(md + 258);
    const auto *md_261 = buffer.data(md + 261);
    const auto *md_262 = buffer.data(md + 262);
    const auto *md_264 = buffer.data(md + 264);
    const auto *md_266 = buffer.data(md + 266);
    const auto *md_267 = buffer.data(md + 267);
    const auto *md_269 = buffer.data(md + 269);
    const auto *md_270 = buffer.data(md + 270);
    const auto *md_273 = buffer.data(md + 273);
    const auto *md_274 = buffer.data(md + 274);
    const auto *md_275 = buffer.data(md + 275);
    const auto *md_276 = buffer.data(md + 276);
    const auto *md_279 = buffer.data(md + 279);
    const auto *md_280 = buffer.data(md + 280);
    const auto *md_281 = buffer.data(md + 281);
    const auto *md_282 = buffer.data(md + 282);
    const auto *md_285 = buffer.data(md + 285);
    const auto *md_286 = buffer.data(md + 286);
    const auto *md_287 = buffer.data(md + 287);
    const auto *md_288 = buffer.data(md + 288);
    const auto *md_291 = buffer.data(md + 291);
    const auto *md_292 = buffer.data(md + 292);
    const auto *md_293 = buffer.data(md + 293);
    const auto *md_294 = buffer.data(md + 294);
    const auto *md_297 = buffer.data(md + 297);
    const auto *md_298 = buffer.data(md + 298);
    const auto *md_299 = buffer.data(md + 299);
    const auto *md_300 = buffer.data(md + 300);
    const auto *md_303 = buffer.data(md + 303);
    const auto *md_304 = buffer.data(md + 304);
    const auto *md_305 = buffer.data(md + 305);
    const auto *md_306 = buffer.data(md + 306);
    const auto *md_309 = buffer.data(md + 309);
    const auto *md_310 = buffer.data(md + 310);
    const auto *md_311 = buffer.data(md + 311);
    const auto *md_312 = buffer.data(md + 312);
    const auto *md_315 = buffer.data(md + 315);
    const auto *md_316 = buffer.data(md + 316);
    const auto *md_317 = buffer.data(md + 317);
    const auto *md_318 = buffer.data(md + 318);
    const auto *md_321 = buffer.data(md + 321);
    const auto *md_322 = buffer.data(md + 322);
    const auto *md_323 = buffer.data(md + 323);
    const auto *md_324 = buffer.data(md + 324);
    const auto *md_327 = buffer.data(md + 327);
    const auto *md_328 = buffer.data(md + 328);
    const auto *md_329 = buffer.data(md + 329);

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pb_x, pb_y, pb_z, ld_198, ld_204, ld_255, \
                         ld_256, md_252, md_255, md_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_8 * ld_204[k]
                   + pb_y[k] * md_252[k];

        t_422[k] = f_14 * ld_198[k]
                   + pb_z[k] * md_252[k];

        t_423[k] = f_3 * ld_255[k]
                   + pb_x[k] * md_255[k];

        t_424[k] = f_3 * ld_256[k]
                   + pb_x[k] * md_256[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, t_430, pa_x, pa_y, pb_x, ld_257, \
                         lf_350, lf_426, lf_427, lf_428, lf_429, \
                         md_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_3 * ld_257[k]
                   + pb_x[k] * md_257[k];

        t_426[k] = pa_x[k] * lf_426[k];

        t_427[k] = pa_x[k] * lf_427[k];

        t_428[k] = pa_x[k] * lf_428[k];

        t_429[k] = pa_x[k] * lf_429[k];

        t_430[k] = pa_y[k] * lf_350[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, t_435, pa_y, pb_x, pb_y, ld_210, ld_261, \
                         ld_262, lf_352, lf_355, md_258, md_261, \
                         md_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_3 * ld_210[k]
                   + pb_y[k] * md_258[k];

        t_432[k] = pa_y[k] * lf_352[k];

        t_433[k] = f_3 * ld_261[k]
                   + pb_x[k] * md_261[k];

        t_434[k] = f_3 * ld_262[k]
                   + pb_x[k] * md_262[k];

        t_435[k] = pa_y[k] * lf_355[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, t_441, pa_x, pb_y, ld_264, lf_436, \
                         lf_437, lf_438, lf_439, lf_440, md_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = pa_x[k] * lf_436[k];

        t_437[k] = pa_x[k] * lf_437[k];

        t_438[k] = pa_x[k] * lf_438[k];

        t_439[k] = pa_x[k] * lf_439[k];

        t_440[k] = f_5 * ld_264[k]
                   + pa_x[k] * lf_440[k];

        t_441[k] = pb_y[k] * md_264[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, pb_x, pb_y, pb_z, ld_210, ld_267, ld_269, \
                         md_264, md_266, md_267, md_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = f_4 * ld_210[k]
                   + pb_z[k] * md_264[k];

        t_443[k] = f_3 * ld_267[k]
                   + pb_x[k] * md_267[k];

        t_444[k] = pb_y[k] * md_266[k];

        t_445[k] = f_3 * ld_269[k]
                   + pb_x[k] * md_269[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, t_449, t_450, pa_x, pb_x, pb_y, lf_446, lf_447, \
                         lf_449, mp0_135, mp1_135, md_269, md_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = pa_x[k] * lf_446[k];

        t_447[k] = pa_x[k] * lf_447[k];

        t_448[k] = pb_y[k] * md_269[k];

        t_449[k] = pa_x[k] * lf_449[k];

        t_450[k] = f_1 * mp0_135[k]
                   - f_2 * mp1_135[k]
                   + pb_x[k] * md_270[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, t_455, pb_x, pb_y, pb_z, ld_216, md_270, \
                         md_273, md_274, md_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_0 * ld_216[k]
                   + pb_y[k] * md_270[k];

        t_452[k] = pb_z[k] * md_270[k];

        t_453[k] = pb_x[k] * md_273[k];

        t_454[k] = pb_x[k] * md_274[k];

        t_455[k] = pb_x[k] * md_275[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pb_y, pb_z, ld_219, ld_221, mp0_136, \
                         mp0_137, mp1_136, mp1_137, md_273, md_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_0 * ld_219[k]
                   + f_1 * mp0_136[k]
                   - f_2 * mp1_136[k]
                   + pb_y[k] * md_273[k];

        t_457[k] = pb_z[k] * md_273[k];

        t_458[k] = f_0 * ld_221[k]
                   + pb_y[k] * md_275[k];

        t_459[k] = f_1 * mp0_137[k]
                   - f_2 * mp1_137[k]
                   + pb_z[k] * md_275[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, t_465, pa_z, pb_x, pb_z, ld_216, \
                         lf_360, lf_361, md_276, md_279, md_280, \
                         md_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = pa_z[k] * lf_360[k];

        t_461[k] = pa_z[k] * lf_361[k];

        t_462[k] = f_3 * ld_216[k]
                   + pb_z[k] * md_276[k];

        t_463[k] = pb_x[k] * md_279[k];

        t_464[k] = pb_x[k] * md_280[k];

        t_465[k] = pb_x[k] * md_281[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_z, pb_y, pb_z, ld_219, ld_221, ld_227, \
                         lf_366, lf_369, md_279, md_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = pa_z[k] * lf_366[k];

        t_467[k] = f_3 * ld_219[k]
                   + pb_z[k] * md_279[k];

        t_468[k] = f_4 * ld_227[k]
                   + pb_y[k] * md_281[k];

        t_469[k] = f_5 * ld_221[k]
                   + pa_z[k] * lf_369[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, pb_x, pb_y, pb_z, ld_222, ld_228, \
                         mp0_141, mp1_141, md_282, md_285, md_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_1 * mp0_141[k]
                   - f_2 * mp1_141[k]
                   + pb_x[k] * md_282[k];

        t_471[k] = f_9 * ld_228[k]
                   + pb_y[k] * md_282[k];

        t_472[k] = f_8 * ld_222[k]
                   + pb_z[k] * md_282[k];

        t_473[k] = pb_x[k] * md_285[k];

        t_474[k] = pb_x[k] * md_286[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_z, pb_x, pb_y, pb_z, kf0_286, kf1_286, \
                         ld_225, ld_233, lf_376, md_285, md_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = pb_x[k] * md_287[k];

        t_476[k] = f_6 * kf0_286[k]
                   - f_7 * kf1_286[k]
                   + pa_z[k] * lf_376[k];

        t_477[k] = f_8 * ld_225[k]
                   + pb_z[k] * md_285[k];

        t_478[k] = f_9 * ld_233[k]
                   + pb_y[k] * md_287[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, pa_y, pb_x, pb_y, pb_z, kf0_309, kf1_309, \
                         ld_228, ld_234, lf_389, mp0_144, mp1_144, \
                         md_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_10 * kf0_309[k]
                   - f_11 * kf1_309[k]
                   + pa_y[k] * lf_389[k];

        t_480[k] = f_1 * mp0_144[k]
                   - f_2 * mp1_144[k]
                   + pb_x[k] * md_288[k];

        t_481[k] = f_14 * ld_234[k]
                   + pb_y[k] * md_288[k];

        t_482[k] = f_5 * ld_228[k]
                   + pb_z[k] * md_288[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, t_487, pa_z, pb_x, pb_z, kf0_296, \
                         kf1_296, ld_231, lf_386, md_291, md_292, \
                         md_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = pb_x[k] * md_291[k];

        t_484[k] = pb_x[k] * md_292[k];

        t_485[k] = pb_x[k] * md_293[k];

        t_486[k] = f_12 * kf0_296[k]
                   - f_13 * kf1_296[k]
                   + pa_z[k] * lf_386[k];

        t_487[k] = f_5 * ld_231[k]
                   + pb_z[k] * md_291[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, pa_y, pb_x, pb_y, kf0_319, kf1_319, \
                         ld_239, ld_240, lf_399, mp0_147, mp1_147, md_293, \
                         md_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_14 * ld_239[k]
                   + pb_y[k] * md_293[k];

        t_489[k] = f_15 * kf0_319[k]
                   - f_16 * kf1_319[k]
                   + pa_y[k] * lf_399[k];

        t_490[k] = f_1 * mp0_147[k]
                   - f_2 * mp1_147[k]
                   + pb_x[k] * md_294[k];

        t_491[k] = f_20 * ld_240[k]
                   + pb_y[k] * md_294[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, pa_z, pb_x, pb_z, kf0_306, \
                         kf1_306, ld_234, lf_396, md_294, md_297, md_298, \
                         md_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_19 * ld_234[k]
                   + pb_z[k] * md_294[k];

        t_493[k] = pb_x[k] * md_297[k];

        t_494[k] = pb_x[k] * md_298[k];

        t_495[k] = pb_x[k] * md_299[k];

        t_496[k] = f_17 * kf0_306[k]
                   - f_18 * kf1_306[k]
                   + pa_z[k] * lf_396[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, pa_y, pb_y, pb_z, kf0_329, kf1_329, ld_237, \
                         ld_245, lf_409, md_297, md_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = f_19 * ld_237[k]
                   + pb_z[k] * md_297[k];

        t_498[k] = f_20 * ld_245[k]
                   + pb_y[k] * md_299[k];

        t_499[k] = f_21 * kf0_329[k]
                   - f_22 * kf1_329[k]
                   + pa_y[k] * lf_409[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, pb_x, pb_y, pb_z, ld_240, ld_246, \
                         mp0_150, mp1_150, md_300, md_303, md_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_1 * mp0_150[k]
                   - f_2 * mp1_150[k]
                   + pb_x[k] * md_300[k];

        t_501[k] = f_19 * ld_246[k]
                   + pb_y[k] * md_300[k];

        t_502[k] = f_20 * ld_240[k]
                   + pb_z[k] * md_300[k];

        t_503[k] = pb_x[k] * md_303[k];

        t_504[k] = pb_x[k] * md_304[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pa_z, pb_x, pb_y, pb_z, kf0_316, kf1_316, \
                         ld_243, ld_251, lf_406, md_303, md_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = pb_x[k] * md_305[k];

        t_506[k] = f_21 * kf0_316[k]
                   - f_22 * kf1_316[k]
                   + pa_z[k] * lf_406[k];

        t_507[k] = f_20 * ld_243[k]
                   + pb_z[k] * md_303[k];

        t_508[k] = f_19 * ld_251[k]
                   + pb_y[k] * md_305[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, pa_y, pb_x, pb_y, pb_z, kf0_339, kf1_339, \
                         ld_246, ld_252, lf_419, mp0_153, mp1_153, \
                         md_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_17 * kf0_339[k]
                   - f_18 * kf1_339[k]
                   + pa_y[k] * lf_419[k];

        t_510[k] = f_1 * mp0_153[k]
                   - f_2 * mp1_153[k]
                   + pb_x[k] * md_306[k];

        t_511[k] = f_5 * ld_252[k]
                   + pb_y[k] * md_306[k];

        t_512[k] = f_14 * ld_246[k]
                   + pb_z[k] * md_306[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, t_516, t_517, pa_z, pb_x, pb_z, kf0_326, \
                         kf1_326, ld_249, lf_416, md_309, md_310, \
                         md_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = pb_x[k] * md_309[k];

        t_514[k] = pb_x[k] * md_310[k];

        t_515[k] = pb_x[k] * md_311[k];

        t_516[k] = f_15 * kf0_326[k]
                   - f_16 * kf1_326[k]
                   + pa_z[k] * lf_416[k];

        t_517[k] = f_14 * ld_249[k]
                   + pb_z[k] * md_309[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pa_y, pb_x, pb_y, kf0_349, kf1_349, \
                         ld_257, ld_258, lf_429, mp0_156, mp1_156, md_311, \
                         md_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_5 * ld_257[k]
                   + pb_y[k] * md_311[k];

        t_519[k] = f_12 * kf0_349[k]
                   - f_13 * kf1_349[k]
                   + pa_y[k] * lf_429[k];

        t_520[k] = f_1 * mp0_156[k]
                   - f_2 * mp1_156[k]
                   + pb_x[k] * md_312[k];

        t_521[k] = f_8 * ld_258[k]
                   + pb_y[k] * md_312[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, pa_z, pb_x, pb_z, kf0_336, \
                         kf1_336, ld_252, lf_426, md_312, md_315, md_316, \
                         md_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_9 * ld_252[k]
                   + pb_z[k] * md_312[k];

        t_523[k] = pb_x[k] * md_315[k];

        t_524[k] = pb_x[k] * md_316[k];

        t_525[k] = pb_x[k] * md_317[k];

        t_526[k] = f_10 * kf0_336[k]
                   - f_11 * kf1_336[k]
                   + pa_z[k] * lf_426[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, pa_y, pb_y, pb_z, kf0_359, kf1_359, \
                         ld_255, ld_263, lf_439, lf_440, md_315, \
                         md_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_9 * ld_255[k]
                   + pb_z[k] * md_315[k];

        t_528[k] = f_8 * ld_263[k]
                   + pb_y[k] * md_317[k];

        t_529[k] = f_6 * kf0_359[k]
                   - f_7 * kf1_359[k]
                   + pa_y[k] * lf_439[k];

        t_530[k] = pa_y[k] * lf_440[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, t_535, pa_y, pb_x, pb_y, ld_264, lf_442, \
                         md_318, md_321, md_322, md_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = f_3 * ld_264[k]
                   + pb_y[k] * md_318[k];

        t_532[k] = pa_y[k] * lf_442[k];

        t_533[k] = pb_x[k] * md_321[k];

        t_534[k] = pb_x[k] * md_322[k];

        t_535[k] = pb_x[k] * md_323[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pa_y, pb_y, pb_z, ld_261, ld_267, ld_269, \
                         lf_446, lf_449, md_321, md_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_5 * ld_267[k]
                   + pa_y[k] * lf_446[k];

        t_537[k] = f_4 * ld_261[k]
                   + pb_z[k] * md_321[k];

        t_538[k] = f_3 * ld_269[k]
                   + pb_y[k] * md_323[k];

        t_539[k] = pa_y[k] * lf_449[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, t_545, pb_x, pb_y, pb_z, ld_264, \
                         mp0_162, mp1_162, md_324, md_327, md_328, \
                         md_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_1 * mp0_162[k]
                   - f_2 * mp1_162[k]
                   + pb_x[k] * md_324[k];

        t_541[k] = pb_y[k] * md_324[k];

        t_542[k] = f_0 * ld_264[k]
                   + pb_z[k] * md_324[k];

        t_543[k] = pb_x[k] * md_327[k];

        t_544[k] = pb_x[k] * md_328[k];

        t_545[k] = pb_x[k] * md_329[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, pb_y, pb_z, ld_267, ld_269, mp0_163, \
                         mp0_164, mp1_163, mp1_164, md_327, md_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = f_1 * mp0_163[k]
                   - f_2 * mp1_163[k]
                   + pb_y[k] * md_327[k];

        t_547[k] = f_0 * ld_267[k]
                   + pb_z[k] * md_327[k];

        t_548[k] = pb_y[k] * md_329[k];

        t_549[k] = f_0 * ld_269[k]
                   + f_1 * mp0_164[k]
                   - f_2 * mp1_164[k]
                   + pb_z[k] * md_329[k];
    }
}

auto
compute_prim_mf_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t kf0, const size_t kf1,
                                     const size_t ld, const size_t lf, const size_t mp0,
                                     const size_t mp1, const size_t md, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    compute_prim_mf_electron_repulsion_0_piece0(buffer, target, pa, pb, kf0, kf1, ld, lf, mp0,
                                                mp1, md, ncols, alpha, beta, p);

    compute_prim_mf_electron_repulsion_0_piece1(buffer, target, pa, pb, kf0, kf1, ld, lf, mp0,
                                                mp1, md, ncols, alpha, beta, p);

    compute_prim_mf_electron_repulsion_0_piece2(buffer, target, pa, pb, kf0, kf1, ld, lf, mp0,
                                                mp1, md, ncols, alpha, beta, p);

    compute_prim_mf_electron_repulsion_0_piece3(buffer, target, pa, pb, kf0, kf1, ld, lf, mp0,
                                                mp1, md, ncols, alpha, beta, p);
}

}  // namespace simdt2ceri
