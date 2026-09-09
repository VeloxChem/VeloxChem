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


#include "SimdElectronRepulsionVrrRecKG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_kg_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t hg0,
                                            const size_t hg1, const size_t if_, const size_t ig,
                                            const size_t kd0, const size_t kd1, const size_t kf,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 3.0 / p;
    const auto f_8 = 2.0 / p;
    const auto f_9 = 0.5 / alpha;
    const auto f_10 = 0.5 * beta / (alpha * p);
    const auto f_11 = 2.5 / p;
    const auto f_12 = 2.0 / alpha;
    const auto f_13 = 2.0 * beta / (alpha * p);
    const auto f_14 = 1.0 / alpha;
    const auto f_15 = beta / (alpha * p);
    const auto f_16 = 1.5 / p;
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);

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

    const auto *hg0_0 = buffer.data(hg0 + 0);
    const auto *hg0_15 = buffer.data(hg0 + 15);
    const auto *hg0_30 = buffer.data(hg0 + 30);
    const auto *hg0_55 = buffer.data(hg0 + 55);
    const auto *hg0_89 = buffer.data(hg0 + 89);
    const auto *hg0_100 = buffer.data(hg0 + 100);

    const auto *hg1_0 = buffer.data(hg1 + 0);
    const auto *hg1_15 = buffer.data(hg1 + 15);
    const auto *hg1_30 = buffer.data(hg1 + 30);
    const auto *hg1_55 = buffer.data(hg1 + 55);
    const auto *hg1_89 = buffer.data(hg1 + 89);
    const auto *hg1_100 = buffer.data(hg1 + 100);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__1 = buffer.data(if_ + 1);
    const auto *if__2 = buffer.data(if_ + 2);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__18 = buffer.data(if_ + 18);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__63 = buffer.data(if_ + 63);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__77 = buffer.data(if_ + 77);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__87 = buffer.data(if_ + 87);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__95 = buffer.data(if_ + 95);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__97 = buffer.data(if_ + 97);

    const auto *ig_0 = buffer.data(ig + 0);
    const auto *ig_3 = buffer.data(ig + 3);
    const auto *ig_5 = buffer.data(ig + 5);
    const auto *ig_6 = buffer.data(ig + 6);
    const auto *ig_9 = buffer.data(ig + 9);
    const auto *ig_10 = buffer.data(ig + 10);
    const auto *ig_12 = buffer.data(ig + 12);
    const auto *ig_14 = buffer.data(ig + 14);
    const auto *ig_15 = buffer.data(ig + 15);
    const auto *ig_16 = buffer.data(ig + 16);
    const auto *ig_18 = buffer.data(ig + 18);
    const auto *ig_21 = buffer.data(ig + 21);
    const auto *ig_25 = buffer.data(ig + 25);
    const auto *ig_30 = buffer.data(ig + 30);
    const auto *ig_32 = buffer.data(ig + 32);
    const auto *ig_35 = buffer.data(ig + 35);
    const auto *ig_39 = buffer.data(ig + 39);
    const auto *ig_42 = buffer.data(ig + 42);
    const auto *ig_44 = buffer.data(ig + 44);
    const auto *ig_45 = buffer.data(ig + 45);
    const auto *ig_46 = buffer.data(ig + 46);
    const auto *ig_48 = buffer.data(ig + 48);
    const auto *ig_50 = buffer.data(ig + 50);
    const auto *ig_51 = buffer.data(ig + 51);
    const auto *ig_55 = buffer.data(ig + 55);
    const auto *ig_57 = buffer.data(ig + 57);
    const auto *ig_59 = buffer.data(ig + 59);
    const auto *ig_75 = buffer.data(ig + 75);
    const auto *ig_77 = buffer.data(ig + 77);
    const auto *ig_78 = buffer.data(ig + 78);
    const auto *ig_80 = buffer.data(ig + 80);
    const auto *ig_84 = buffer.data(ig + 84);
    const auto *ig_85 = buffer.data(ig + 85);
    const auto *ig_87 = buffer.data(ig + 87);
    const auto *ig_89 = buffer.data(ig + 89);
    const auto *ig_100 = buffer.data(ig + 100);

    const auto *kd0_0 = buffer.data(kd0 + 0);
    const auto *kd0_3 = buffer.data(kd0 + 3);
    const auto *kd0_5 = buffer.data(kd0 + 5);
    const auto *kd0_18 = buffer.data(kd0 + 18);
    const auto *kd0_21 = buffer.data(kd0 + 21);
    const auto *kd0_23 = buffer.data(kd0 + 23);
    const auto *kd0_30 = buffer.data(kd0 + 30);
    const auto *kd0_33 = buffer.data(kd0 + 33);
    const auto *kd0_35 = buffer.data(kd0 + 35);
    const auto *kd0_36 = buffer.data(kd0 + 36);
    const auto *kd0_39 = buffer.data(kd0 + 39);
    const auto *kd0_41 = buffer.data(kd0 + 41);
    const auto *kd0_54 = buffer.data(kd0 + 54);
    const auto *kd0_59 = buffer.data(kd0 + 59);

    const auto *kd1_0 = buffer.data(kd1 + 0);
    const auto *kd1_3 = buffer.data(kd1 + 3);
    const auto *kd1_5 = buffer.data(kd1 + 5);
    const auto *kd1_18 = buffer.data(kd1 + 18);
    const auto *kd1_21 = buffer.data(kd1 + 21);
    const auto *kd1_23 = buffer.data(kd1 + 23);
    const auto *kd1_30 = buffer.data(kd1 + 30);
    const auto *kd1_33 = buffer.data(kd1 + 33);
    const auto *kd1_35 = buffer.data(kd1 + 35);
    const auto *kd1_36 = buffer.data(kd1 + 36);
    const auto *kd1_39 = buffer.data(kd1 + 39);
    const auto *kd1_41 = buffer.data(kd1 + 41);
    const auto *kd1_54 = buffer.data(kd1 + 54);
    const auto *kd1_59 = buffer.data(kd1 + 59);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_1 = buffer.data(kf + 1);
    const auto *kf_2 = buffer.data(kf + 2);
    const auto *kf_3 = buffer.data(kf + 3);
    const auto *kf_5 = buffer.data(kf + 5);
    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_9 = buffer.data(kf + 9);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_13 = buffer.data(kf + 13);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_18 = buffer.data(kf + 18);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_25 = buffer.data(kf + 25);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_27 = buffer.data(kf + 27);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_30 = buffer.data(kf + 30);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_33 = buffer.data(kf + 33);
    const auto *kf_36 = buffer.data(kf + 36);
    const auto *kf_37 = buffer.data(kf + 37);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_39 = buffer.data(kf + 39);
    const auto *kf_42 = buffer.data(kf + 42);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_49 = buffer.data(kf + 49);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_51 = buffer.data(kf + 51);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_61 = buffer.data(kf + 61);
    const auto *kf_62 = buffer.data(kf + 62);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_67 = buffer.data(kf + 67);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_69 = buffer.data(kf + 69);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_79 = buffer.data(kf + 79);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_87 = buffer.data(kf + 87);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_95 = buffer.data(kf + 95);
    const auto *kf_96 = buffer.data(kf + 96);
    const auto *kf_97 = buffer.data(kf + 97);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, if__0, kd0_0, kd1_0, \
                         kf_0, kf_1, kf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * if__0[k]
                 + f_1 * kd0_0[k]
                 - f_2 * kd1_0[k]
                 + pb_x[k] * kf_0[k];

        t_1[k] = pb_y[k] * kf_0[k];

        t_2[k] = pb_z[k] * kf_0[k];

        t_3[k] = f_3 * kd0_0[k]
                 - f_4 * kd1_0[k]
                 + pb_y[k] * kf_1[k];

        t_4[k] = pb_y[k] * kf_2[k];

        t_5[k] = f_3 * kd0_0[k]
                 - f_4 * kd1_0[k]
                 + pb_z[k] * kf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, if__6, if__9, kd0_3, \
                         kd1_3, kf_3, kf_5, kf_6, kf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * if__6[k]
                 + pb_x[k] * kf_6[k];

        t_7[k] = pb_z[k] * kf_3[k];

        t_8[k] = pb_y[k] * kf_5[k];

        t_9[k] = f_0 * if__9[k]
                 + pb_x[k] * kf_9[k];

        t_10[k] = f_1 * kd0_3[k]
                  - f_2 * kd1_3[k]
                  + pb_y[k] * kf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, ig_0, kd0_5, kd1_5, \
                         kf_6, kf_8, kf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * kf_6[k];

        t_12[k] = f_3 * kd0_5[k]
                  - f_4 * kd1_5[k]
                  + pb_y[k] * kf_8[k];

        t_13[k] = pb_y[k] * kf_9[k];

        t_14[k] = f_1 * kd0_5[k]
                  - f_2 * kd1_5[k]
                  + pb_z[k] * kf_9[k];

        t_15[k] = pa_y[k] * ig_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pb_y, pb_z, if__0, if__1, ig_3, \
                         ig_5, kf_10, kf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * if__0[k]
                  + pb_y[k] * kf_10[k];

        t_17[k] = pb_z[k] * kf_10[k];

        t_18[k] = f_6 * if__1[k]
                  + pa_y[k] * ig_3[k];

        t_19[k] = pb_z[k] * kf_11[k];

        t_20[k] = pa_y[k] * ig_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_z, if__6, if__16, \
                         if__18, ig_9, ig_10, kf_13, kf_16, kf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * if__16[k]
                  + pb_x[k] * kf_16[k];

        t_22[k] = pb_z[k] * kf_13[k];

        t_23[k] = f_7 * if__18[k]
                  + pb_x[k] * kf_18[k];

        t_24[k] = pa_y[k] * ig_9[k];

        t_25[k] = f_8 * if__6[k]
                  + pa_y[k] * ig_10[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, pb_z, if__8, if__9, \
                         ig_0, ig_12, ig_14, kf_16, kf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * kf_16[k];

        t_27[k] = f_6 * if__8[k]
                  + pa_y[k] * ig_12[k];

        t_28[k] = f_5 * if__9[k]
                  + pb_y[k] * kf_19[k];

        t_29[k] = pa_y[k] * ig_14[k];

        t_30[k] = pa_z[k] * ig_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, if__0, if__2, \
                         ig_3, ig_5, ig_6, kf_20, kf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * kf_20[k];

        t_32[k] = f_5 * if__0[k]
                  + pb_z[k] * kf_20[k];

        t_33[k] = pa_z[k] * ig_3[k];

        t_34[k] = pb_y[k] * kf_22[k];

        t_35[k] = f_6 * if__2[k]
                  + pa_z[k] * ig_5[k];

        t_36[k] = pa_z[k] * ig_6[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, pb_y, if__27, if__29, ig_10, \
                         kf_25, kf_27, kf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_7 * if__27[k]
                  + pb_x[k] * kf_27[k];

        t_38[k] = pb_y[k] * kf_25[k];

        t_39[k] = f_7 * if__29[k]
                  + pb_x[k] * kf_29[k];

        t_40[k] = pa_z[k] * ig_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_z, pb_y, pb_z, if__6, if__7, if__9, ig_12, \
                         ig_14, kf_26, kf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * if__6[k]
                  + pb_z[k] * kf_26[k];

        t_42[k] = f_6 * if__7[k]
                  + pa_z[k] * ig_12[k];

        t_43[k] = pb_y[k] * kf_29[k];

        t_44[k] = f_8 * if__9[k]
                  + pa_z[k] * ig_14[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_y, pb_z, hg0_0, hg1_0, if__10, ig_15, \
                         kf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_9 * hg0_0[k]
                  - f_10 * hg1_0[k]
                  + pa_y[k] * ig_15[k];

        t_46[k] = f_6 * if__10[k]
                  + pb_y[k] * kf_30[k];

        t_47[k] = pb_z[k] * kf_30[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_x, pb_z, if__33, if__36, kd0_18, kd0_21, \
                         kd1_18, kd1_21, kf_31, kf_32, kf_33, kf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_11 * if__33[k]
                  + f_3 * kd0_21[k]
                  - f_4 * kd1_21[k]
                  + pb_x[k] * kf_33[k];

        t_49[k] = pb_z[k] * kf_31[k];

        t_50[k] = f_3 * kd0_18[k]
                  - f_4 * kd1_18[k]
                  + pb_z[k] * kf_32[k];

        t_51[k] = f_11 * if__36[k]
                  + pb_x[k] * kf_36[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_z, hg0_55, hg1_55, if__38, \
                         if__39, ig_55, kf_33, kf_38, kf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_z[k] * kf_33[k];

        t_53[k] = f_11 * if__38[k]
                  + pb_x[k] * kf_38[k];

        t_54[k] = f_11 * if__39[k]
                  + pb_x[k] * kf_39[k];

        t_55[k] = f_12 * hg0_55[k]
                  - f_13 * hg1_55[k]
                  + pa_x[k] * ig_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pb_z, if__19, kd0_21, kd0_23, kd1_21, \
                         kd1_23, kf_36, kf_37, kf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * kf_36[k];

        t_57[k] = f_3 * kd0_21[k]
                  - f_4 * kd1_21[k]
                  + pb_z[k] * kf_37[k];

        t_58[k] = f_6 * if__19[k]
                  + pb_y[k] * kf_39[k];

        t_59[k] = f_1 * kd0_23[k]
                  - f_2 * kd1_23[k]
                  + pb_z[k] * kf_39[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, if__22, ig_16, \
                         ig_18, ig_30, ig_32, ig_35, kf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * ig_30[k];

        t_61[k] = pa_z[k] * ig_16[k];

        t_62[k] = pa_y[k] * ig_32[k];

        t_63[k] = pa_z[k] * ig_18[k];

        t_64[k] = f_5 * if__22[k]
                  + pb_y[k] * kf_42[k];

        t_65[k] = pa_y[k] * ig_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, if__47, if__48, \
                         ig_21, ig_25, ig_39, kf_47, kf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * ig_21[k];

        t_67[k] = f_11 * if__47[k]
                  + pb_x[k] * kf_47[k];

        t_68[k] = f_11 * if__48[k]
                  + pb_x[k] * kf_48[k];

        t_69[k] = pa_y[k] * ig_39[k];

        t_70[k] = pa_z[k] * ig_25[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, pb_z, if__16, if__28, if__29, \
                         ig_42, ig_44, kf_46, kf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * if__16[k]
                  + pb_z[k] * kf_46[k];

        t_72[k] = f_6 * if__28[k]
                  + pa_y[k] * ig_42[k];

        t_73[k] = f_5 * if__29[k]
                  + pb_y[k] * kf_49[k];

        t_74[k] = pa_y[k] * ig_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, hg0_0, hg1_0, if__20, \
                         ig_30, kd0_30, kd1_30, kf_50, kf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_9 * hg0_0[k]
                  - f_10 * hg1_0[k]
                  + pa_z[k] * ig_30[k];

        t_76[k] = pb_y[k] * kf_50[k];

        t_77[k] = f_6 * if__20[k]
                  + pb_z[k] * kf_50[k];

        t_78[k] = f_3 * kd0_30[k]
                  - f_4 * kd1_30[k]
                  + pb_y[k] * kf_51[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_x, pb_y, if__55, if__56, if__57, \
                         kd0_35, kd1_35, kf_52, kf_55, kf_56, kf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_y[k] * kf_52[k];

        t_80[k] = f_11 * if__55[k]
                  + f_3 * kd0_35[k]
                  - f_4 * kd1_35[k]
                  + pb_x[k] * kf_55[k];

        t_81[k] = f_11 * if__56[k]
                  + pb_x[k] * kf_56[k];

        t_82[k] = f_11 * if__57[k]
                  + pb_x[k] * kf_57[k];

        t_83[k] = pb_y[k] * kf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, pb_y, pb_z, if__26, if__59, kd0_33, \
                         kd0_35, kd1_33, kd1_35, kf_56, kf_58, kf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_11 * if__59[k]
                  + pb_x[k] * kf_59[k];

        t_85[k] = f_1 * kd0_33[k]
                  - f_2 * kd1_33[k]
                  + pb_y[k] * kf_56[k];

        t_86[k] = f_6 * if__26[k]
                  + pb_z[k] * kf_56[k];

        t_87[k] = f_3 * kd0_35[k]
                  - f_4 * kd1_35[k]
                  + pb_y[k] * kf_58[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pa_y, pb_y, hg0_15, hg0_89, hg1_15, \
                         hg1_89, if__30, ig_45, ig_89, kf_59, kf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pb_y[k] * kf_59[k];

        t_89[k] = f_12 * hg0_89[k]
                  - f_13 * hg1_89[k]
                  + pa_x[k] * ig_89[k];

        t_90[k] = f_14 * hg0_15[k]
                  - f_15 * hg1_15[k]
                  + pa_y[k] * ig_45[k];

        t_91[k] = f_16 * if__30[k]
                  + pb_y[k] * kf_60[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_z, if__63, kd0_36, kd0_39, kd1_36, \
                         kd1_39, kf_60, kf_61, kf_62, kf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_z[k] * kf_60[k];

        t_93[k] = f_8 * if__63[k]
                  + f_3 * kd0_39[k]
                  - f_4 * kd1_39[k]
                  + pb_x[k] * kf_63[k];

        t_94[k] = pb_z[k] * kf_61[k];

        t_95[k] = f_3 * kd0_36[k]
                  - f_4 * kd1_36[k]
                  + pb_z[k] * kf_62[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_x, pb_z, if__66, if__68, if__69, kf_63, \
                         kf_66, kf_68, kf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_8 * if__66[k]
                  + pb_x[k] * kf_66[k];

        t_97[k] = pb_z[k] * kf_63[k];

        t_98[k] = f_8 * if__68[k]
                  + pb_x[k] * kf_68[k];

        t_99[k] = f_8 * if__69[k]
                  + pb_x[k] * kf_69[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_y, pb_z, hg0_100, hg1_100, \
                         if__39, ig_100, kd0_39, kd1_39, kf_66, kf_67, \
                         kf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_17 * hg0_100[k]
                   - f_18 * hg1_100[k]
                   + pa_x[k] * ig_100[k];

        t_101[k] = pb_z[k] * kf_66[k];

        t_102[k] = f_3 * kd0_39[k]
                   - f_4 * kd1_39[k]
                   + pb_z[k] * kf_67[k];

        t_103[k] = f_16 * if__39[k]
                   + pb_y[k] * kf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_z, pb_z, if__30, ig_45, ig_46, \
                         ig_48, kd0_41, kd1_41, kf_69, kf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * kd0_41[k]
                   - f_2 * kd1_41[k]
                   + pb_z[k] * kf_69[k];

        t_105[k] = pa_z[k] * ig_45[k];

        t_106[k] = pa_z[k] * ig_46[k];

        t_107[k] = f_5 * if__30[k]
                   + pb_z[k] * kf_70[k];

        t_108[k] = pa_z[k] * ig_48[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_z, pb_x, pb_y, if__32, if__42, if__77, \
                         ig_50, ig_51, kf_72, kf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_6 * if__42[k]
                   + pb_y[k] * kf_72[k];

        t_110[k] = f_6 * if__32[k]
                   + pa_z[k] * ig_50[k];

        t_111[k] = pa_z[k] * ig_51[k];

        t_112[k] = f_8 * if__77[k]
                   + pb_x[k] * kf_77[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pb_x, pb_z, if__36, if__78, if__79, \
                         ig_55, kf_76, kf_78, kf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_8 * if__78[k]
                   + pb_x[k] * kf_78[k];

        t_114[k] = f_8 * if__79[k]
                   + pb_x[k] * kf_79[k];

        t_115[k] = pa_z[k] * ig_55[k];

        t_116[k] = f_5 * if__36[k]
                   + pb_z[k] * kf_76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pb_y, if__37, if__39, if__49, \
                         ig_57, ig_59, ig_75, kf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_6 * if__37[k]
                   + pa_z[k] * ig_57[k];

        t_118[k] = f_6 * if__49[k]
                   + pb_y[k] * kf_79[k];

        t_119[k] = f_8 * if__39[k]
                   + pa_z[k] * ig_59[k];

        t_120[k] = pa_y[k] * ig_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pa_y, pb_y, if__50, if__51, \
                         if__52, ig_77, ig_78, ig_80, kf_80, kf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * if__50[k]
                   + pb_y[k] * kf_80[k];

        t_122[k] = pa_y[k] * ig_77[k];

        t_123[k] = f_6 * if__51[k]
                   + pa_y[k] * ig_78[k];

        t_124[k] = f_5 * if__52[k]
                   + pb_y[k] * kf_82[k];

        t_125[k] = pa_y[k] * ig_80[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, pa_y, pb_x, if__56, if__86, \
                         if__87, if__88, ig_84, ig_85, kf_86, kf_87, \
                         kf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_8 * if__86[k]
                   + pb_x[k] * kf_86[k];

        t_127[k] = f_8 * if__87[k]
                   + pb_x[k] * kf_87[k];

        t_128[k] = f_8 * if__88[k]
                   + pb_x[k] * kf_88[k];

        t_129[k] = pa_y[k] * ig_84[k];

        t_130[k] = f_8 * if__56[k]
                   + pa_y[k] * ig_85[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_y, pb_y, pb_z, if__46, if__58, if__59, \
                         ig_87, ig_89, kf_86, kf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_6 * if__46[k]
                   + pb_z[k] * kf_86[k];

        t_132[k] = f_6 * if__58[k]
                   + pa_y[k] * ig_87[k];

        t_133[k] = f_5 * if__59[k]
                   + pb_y[k] * kf_89[k];

        t_134[k] = pa_y[k] * ig_89[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pa_z, pb_y, pb_z, hg0_30, hg1_30, if__50, \
                         ig_75, kd0_54, kd1_54, kf_90, kf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_14 * hg0_30[k]
                   - f_15 * hg1_30[k]
                   + pa_z[k] * ig_75[k];

        t_136[k] = pb_y[k] * kf_90[k];

        t_137[k] = f_16 * if__50[k]
                   + pb_z[k] * kf_90[k];

        t_138[k] = f_3 * kd0_54[k]
                   - f_4 * kd1_54[k]
                   + pb_y[k] * kf_91[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pb_x, pb_y, if__95, if__96, \
                         if__97, kd0_59, kd1_59, kf_92, kf_95, kf_96, \
                         kf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * kf_92[k];

        t_140[k] = f_8 * if__95[k]
                   + f_3 * kd0_59[k]
                   - f_4 * kd1_59[k]
                   + pb_x[k] * kf_95[k];

        t_141[k] = f_8 * if__96[k]
                   + pb_x[k] * kf_96[k];

        t_142[k] = f_8 * if__97[k]
                   + pb_x[k] * kf_97[k];

        t_143[k] = pb_y[k] * kf_95[k];
    }
}

static auto
compute_prim_kg_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t hg0,
                                            const size_t hg1, const size_t if_, const size_t ig,
                                            const size_t kd0, const size_t kd1, const size_t kf,
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
    const auto f_11 = 2.5 / p;
    const auto f_12 = 2.0 / alpha;
    const auto f_13 = 2.0 * beta / (alpha * p);
    const auto f_14 = 1.0 / alpha;
    const auto f_15 = beta / (alpha * p);
    const auto f_16 = 1.5 / p;
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);

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

    const auto *hg0_45 = buffer.data(hg0 + 45);
    const auto *hg0_48 = buffer.data(hg0 + 48);
    const auto *hg0_75 = buffer.data(hg0 + 75);
    const auto *hg0_80 = buffer.data(hg0 + 80);
    const auto *hg0_90 = buffer.data(hg0 + 90);
    const auto *hg0_93 = buffer.data(hg0 + 93);
    const auto *hg0_108 = buffer.data(hg0 + 108);
    const auto *hg0_120 = buffer.data(hg0 + 120);
    const auto *hg0_125 = buffer.data(hg0 + 125);
    const auto *hg0_135 = buffer.data(hg0 + 135);
    const auto *hg0_140 = buffer.data(hg0 + 140);
    const auto *hg0_149 = buffer.data(hg0 + 149);
    const auto *hg0_160 = buffer.data(hg0 + 160);
    const auto *hg0_190 = buffer.data(hg0 + 190);
    const auto *hg0_192 = buffer.data(hg0 + 192);
    const auto *hg0_194 = buffer.data(hg0 + 194);
    const auto *hg0_224 = buffer.data(hg0 + 224);
    const auto *hg0_235 = buffer.data(hg0 + 235);
    const auto *hg0_265 = buffer.data(hg0 + 265);
    const auto *hg0_267 = buffer.data(hg0 + 267);
    const auto *hg0_269 = buffer.data(hg0 + 269);

    const auto *hg1_45 = buffer.data(hg1 + 45);
    const auto *hg1_48 = buffer.data(hg1 + 48);
    const auto *hg1_75 = buffer.data(hg1 + 75);
    const auto *hg1_80 = buffer.data(hg1 + 80);
    const auto *hg1_90 = buffer.data(hg1 + 90);
    const auto *hg1_93 = buffer.data(hg1 + 93);
    const auto *hg1_108 = buffer.data(hg1 + 108);
    const auto *hg1_120 = buffer.data(hg1 + 120);
    const auto *hg1_125 = buffer.data(hg1 + 125);
    const auto *hg1_135 = buffer.data(hg1 + 135);
    const auto *hg1_140 = buffer.data(hg1 + 140);
    const auto *hg1_149 = buffer.data(hg1 + 149);
    const auto *hg1_160 = buffer.data(hg1 + 160);
    const auto *hg1_190 = buffer.data(hg1 + 190);
    const auto *hg1_192 = buffer.data(hg1 + 192);
    const auto *hg1_194 = buffer.data(hg1 + 194);
    const auto *hg1_224 = buffer.data(hg1 + 224);
    const auto *hg1_235 = buffer.data(hg1 + 235);
    const auto *hg1_265 = buffer.data(hg1 + 265);
    const auto *hg1_267 = buffer.data(hg1 + 267);
    const auto *hg1_269 = buffer.data(hg1 + 269);

    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__67 = buffer.data(if_ + 67);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__89 = buffer.data(if_ + 89);
    const auto *if__90 = buffer.data(if_ + 90);
    const auto *if__91 = buffer.data(if_ + 91);
    const auto *if__92 = buffer.data(if_ + 92);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__98 = buffer.data(if_ + 98);
    const auto *if__99 = buffer.data(if_ + 99);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__102 = buffer.data(if_ + 102);
    const auto *if__103 = buffer.data(if_ + 103);
    const auto *if__106 = buffer.data(if_ + 106);
    const auto *if__107 = buffer.data(if_ + 107);
    const auto *if__108 = buffer.data(if_ + 108);
    const auto *if__109 = buffer.data(if_ + 109);
    const auto *if__110 = buffer.data(if_ + 110);
    const auto *if__112 = buffer.data(if_ + 112);
    const auto *if__116 = buffer.data(if_ + 116);
    const auto *if__117 = buffer.data(if_ + 117);
    const auto *if__118 = buffer.data(if_ + 118);
    const auto *if__119 = buffer.data(if_ + 119);
    const auto *if__120 = buffer.data(if_ + 120);
    const auto *if__122 = buffer.data(if_ + 122);
    const auto *if__126 = buffer.data(if_ + 126);
    const auto *if__127 = buffer.data(if_ + 127);
    const auto *if__128 = buffer.data(if_ + 128);
    const auto *if__129 = buffer.data(if_ + 129);
    const auto *if__130 = buffer.data(if_ + 130);
    const auto *if__132 = buffer.data(if_ + 132);
    const auto *if__136 = buffer.data(if_ + 136);
    const auto *if__137 = buffer.data(if_ + 137);
    const auto *if__138 = buffer.data(if_ + 138);
    const auto *if__145 = buffer.data(if_ + 145);
    const auto *if__146 = buffer.data(if_ + 146);
    const auto *if__147 = buffer.data(if_ + 147);
    const auto *if__149 = buffer.data(if_ + 149);
    const auto *if__153 = buffer.data(if_ + 153);
    const auto *if__156 = buffer.data(if_ + 156);
    const auto *if__158 = buffer.data(if_ + 158);
    const auto *if__159 = buffer.data(if_ + 159);
    const auto *if__167 = buffer.data(if_ + 167);
    const auto *if__168 = buffer.data(if_ + 168);
    const auto *if__169 = buffer.data(if_ + 169);
    const auto *if__176 = buffer.data(if_ + 176);
    const auto *if__177 = buffer.data(if_ + 177);
    const auto *if__178 = buffer.data(if_ + 178);
    const auto *if__179 = buffer.data(if_ + 179);

    const auto *ig_90 = buffer.data(ig + 90);
    const auto *ig_91 = buffer.data(ig + 91);
    const auto *ig_93 = buffer.data(ig + 93);
    const auto *ig_95 = buffer.data(ig + 95);
    const auto *ig_96 = buffer.data(ig + 96);
    const auto *ig_100 = buffer.data(ig + 100);
    const auto *ig_102 = buffer.data(ig + 102);
    const auto *ig_104 = buffer.data(ig + 104);
    const auto *ig_108 = buffer.data(ig + 108);
    const auto *ig_120 = buffer.data(ig + 120);
    const auto *ig_125 = buffer.data(ig + 125);
    const auto *ig_135 = buffer.data(ig + 135);
    const auto *ig_137 = buffer.data(ig + 137);
    const auto *ig_138 = buffer.data(ig + 138);
    const auto *ig_140 = buffer.data(ig + 140);
    const auto *ig_144 = buffer.data(ig + 144);
    const auto *ig_145 = buffer.data(ig + 145);
    const auto *ig_147 = buffer.data(ig + 147);
    const auto *ig_149 = buffer.data(ig + 149);
    const auto *ig_150 = buffer.data(ig + 150);
    const auto *ig_151 = buffer.data(ig + 151);
    const auto *ig_153 = buffer.data(ig + 153);
    const auto *ig_155 = buffer.data(ig + 155);
    const auto *ig_156 = buffer.data(ig + 156);
    const auto *ig_160 = buffer.data(ig + 160);
    const auto *ig_162 = buffer.data(ig + 162);
    const auto *ig_164 = buffer.data(ig + 164);
    const auto *ig_168 = buffer.data(ig + 168);
    const auto *ig_180 = buffer.data(ig + 180);
    const auto *ig_183 = buffer.data(ig + 183);
    const auto *ig_185 = buffer.data(ig + 185);
    const auto *ig_190 = buffer.data(ig + 190);
    const auto *ig_192 = buffer.data(ig + 192);
    const auto *ig_194 = buffer.data(ig + 194);
    const auto *ig_195 = buffer.data(ig + 195);
    const auto *ig_200 = buffer.data(ig + 200);
    const auto *ig_224 = buffer.data(ig + 224);
    const auto *ig_235 = buffer.data(ig + 235);
    const auto *ig_265 = buffer.data(ig + 265);
    const auto *ig_267 = buffer.data(ig + 267);
    const auto *ig_269 = buffer.data(ig + 269);

    const auto *kd0_57 = buffer.data(kd0 + 57);
    const auto *kd0_59 = buffer.data(kd0 + 59);
    const auto *kd0_60 = buffer.data(kd0 + 60);
    const auto *kd0_63 = buffer.data(kd0 + 63);
    const auto *kd0_65 = buffer.data(kd0 + 65);
    const auto *kd0_84 = buffer.data(kd0 + 84);
    const auto *kd0_87 = buffer.data(kd0 + 87);
    const auto *kd0_89 = buffer.data(kd0 + 89);
    const auto *kd0_90 = buffer.data(kd0 + 90);
    const auto *kd0_93 = buffer.data(kd0 + 93);
    const auto *kd0_95 = buffer.data(kd0 + 95);

    const auto *kd1_57 = buffer.data(kd1 + 57);
    const auto *kd1_59 = buffer.data(kd1 + 59);
    const auto *kd1_60 = buffer.data(kd1 + 60);
    const auto *kd1_63 = buffer.data(kd1 + 63);
    const auto *kd1_65 = buffer.data(kd1 + 65);
    const auto *kd1_84 = buffer.data(kd1 + 84);
    const auto *kd1_87 = buffer.data(kd1 + 87);
    const auto *kd1_89 = buffer.data(kd1 + 89);
    const auto *kd1_90 = buffer.data(kd1 + 90);
    const auto *kd1_93 = buffer.data(kd1 + 93);
    const auto *kd1_95 = buffer.data(kd1 + 95);

    const auto *kf_96 = buffer.data(kf + 96);
    const auto *kf_98 = buffer.data(kf + 98);
    const auto *kf_99 = buffer.data(kf + 99);
    const auto *kf_100 = buffer.data(kf + 100);
    const auto *kf_101 = buffer.data(kf + 101);
    const auto *kf_102 = buffer.data(kf + 102);
    const auto *kf_103 = buffer.data(kf + 103);
    const auto *kf_106 = buffer.data(kf + 106);
    const auto *kf_107 = buffer.data(kf + 107);
    const auto *kf_108 = buffer.data(kf + 108);
    const auto *kf_109 = buffer.data(kf + 109);
    const auto *kf_110 = buffer.data(kf + 110);
    const auto *kf_112 = buffer.data(kf + 112);
    const auto *kf_116 = buffer.data(kf + 116);
    const auto *kf_117 = buffer.data(kf + 117);
    const auto *kf_118 = buffer.data(kf + 118);
    const auto *kf_119 = buffer.data(kf + 119);
    const auto *kf_120 = buffer.data(kf + 120);
    const auto *kf_122 = buffer.data(kf + 122);
    const auto *kf_126 = buffer.data(kf + 126);
    const auto *kf_127 = buffer.data(kf + 127);
    const auto *kf_128 = buffer.data(kf + 128);
    const auto *kf_129 = buffer.data(kf + 129);
    const auto *kf_130 = buffer.data(kf + 130);
    const auto *kf_132 = buffer.data(kf + 132);
    const auto *kf_136 = buffer.data(kf + 136);
    const auto *kf_137 = buffer.data(kf + 137);
    const auto *kf_138 = buffer.data(kf + 138);
    const auto *kf_139 = buffer.data(kf + 139);
    const auto *kf_140 = buffer.data(kf + 140);
    const auto *kf_141 = buffer.data(kf + 141);
    const auto *kf_142 = buffer.data(kf + 142);
    const auto *kf_145 = buffer.data(kf + 145);
    const auto *kf_146 = buffer.data(kf + 146);
    const auto *kf_147 = buffer.data(kf + 147);
    const auto *kf_148 = buffer.data(kf + 148);
    const auto *kf_149 = buffer.data(kf + 149);
    const auto *kf_150 = buffer.data(kf + 150);
    const auto *kf_151 = buffer.data(kf + 151);
    const auto *kf_152 = buffer.data(kf + 152);
    const auto *kf_153 = buffer.data(kf + 153);
    const auto *kf_156 = buffer.data(kf + 156);
    const auto *kf_157 = buffer.data(kf + 157);
    const auto *kf_158 = buffer.data(kf + 158);
    const auto *kf_159 = buffer.data(kf + 159);
    const auto *kf_160 = buffer.data(kf + 160);
    const auto *kf_162 = buffer.data(kf + 162);
    const auto *kf_166 = buffer.data(kf + 166);
    const auto *kf_167 = buffer.data(kf + 167);
    const auto *kf_168 = buffer.data(kf + 168);
    const auto *kf_169 = buffer.data(kf + 169);
    const auto *kf_170 = buffer.data(kf + 170);
    const auto *kf_172 = buffer.data(kf + 172);
    const auto *kf_176 = buffer.data(kf + 176);
    const auto *kf_177 = buffer.data(kf + 177);
    const auto *kf_178 = buffer.data(kf + 178);
    const auto *kf_179 = buffer.data(kf + 179);
    const auto *kf_180 = buffer.data(kf + 180);
    const auto *kf_182 = buffer.data(kf + 182);

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_x, pb_y, pb_z, if__56, if__99, kd0_57, \
                         kd0_59, kd1_57, kd1_59, kf_96, kf_98, kf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_8 * if__99[k]
                   + pb_x[k] * kf_99[k];

        t_145[k] = f_1 * kd0_57[k]
                   - f_2 * kd1_57[k]
                   + pb_y[k] * kf_96[k];

        t_146[k] = f_16 * if__56[k]
                   + pb_z[k] * kf_96[k];

        t_147[k] = f_3 * kd0_59[k]
                   - f_4 * kd1_59[k]
                   + pb_y[k] * kf_98[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pa_x, pa_y, pb_y, hg0_45, hg0_149, \
                         hg1_45, hg1_149, if__60, ig_90, ig_149, kf_99, \
                         kf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = pb_y[k] * kf_99[k];

        t_149[k] = f_17 * hg0_149[k]
                   - f_18 * hg1_149[k]
                   + pa_x[k] * ig_149[k];

        t_150[k] = f_17 * hg0_45[k]
                   - f_18 * hg1_45[k]
                   + pa_y[k] * ig_90[k];

        t_151[k] = f_8 * if__60[k]
                   + pb_y[k] * kf_100[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pb_z, if__103, kd0_60, kd0_63, \
                         kd1_60, kd1_63, kf_100, kf_101, kf_102, \
                         kf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = pb_z[k] * kf_100[k];

        t_153[k] = f_16 * if__103[k]
                   + f_3 * kd0_63[k]
                   - f_4 * kd1_63[k]
                   + pb_x[k] * kf_103[k];

        t_154[k] = pb_z[k] * kf_101[k];

        t_155[k] = f_3 * kd0_60[k]
                   - f_4 * kd1_60[k]
                   + pb_z[k] * kf_102[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_x, pb_z, if__106, if__108, if__109, \
                         kf_103, kf_106, kf_108, kf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_16 * if__106[k]
                   + pb_x[k] * kf_106[k];

        t_157[k] = pb_z[k] * kf_103[k];

        t_158[k] = f_16 * if__108[k]
                   + pb_x[k] * kf_108[k];

        t_159[k] = f_16 * if__109[k]
                   + pb_x[k] * kf_109[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_x, pb_y, pb_z, hg0_160, hg1_160, \
                         if__69, ig_160, kd0_63, kd1_63, kf_106, kf_107, \
                         kf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_14 * hg0_160[k]
                   - f_15 * hg1_160[k]
                   + pa_x[k] * ig_160[k];

        t_161[k] = pb_z[k] * kf_106[k];

        t_162[k] = f_3 * kd0_63[k]
                   - f_4 * kd1_63[k]
                   + pb_z[k] * kf_107[k];

        t_163[k] = f_8 * if__69[k]
                   + pb_y[k] * kf_109[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, pa_z, pb_z, if__60, ig_90, ig_91, \
                         ig_93, kd0_65, kd1_65, kf_109, kf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_1 * kd0_65[k]
                   - f_2 * kd1_65[k]
                   + pb_z[k] * kf_109[k];

        t_165[k] = pa_z[k] * ig_90[k];

        t_166[k] = pa_z[k] * ig_91[k];

        t_167[k] = f_5 * if__60[k]
                   + pb_z[k] * kf_110[k];

        t_168[k] = pa_z[k] * ig_93[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pa_z, pb_x, pb_y, if__62, if__72, \
                         if__117, ig_95, ig_96, kf_112, kf_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_16 * if__72[k]
                   + pb_y[k] * kf_112[k];

        t_170[k] = f_6 * if__62[k]
                   + pa_z[k] * ig_95[k];

        t_171[k] = pa_z[k] * ig_96[k];

        t_172[k] = f_16 * if__117[k]
                   + pb_x[k] * kf_117[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pa_z, pb_x, pb_z, if__66, if__118, \
                         if__119, ig_100, kf_116, kf_118, kf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_16 * if__118[k]
                   + pb_x[k] * kf_118[k];

        t_174[k] = f_16 * if__119[k]
                   + pb_x[k] * kf_119[k];

        t_175[k] = pa_z[k] * ig_100[k];

        t_176[k] = f_5 * if__66[k]
                   + pb_z[k] * kf_116[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pa_y, pa_z, pb_y, hg0_75, hg1_75, if__67, \
                         if__69, if__79, ig_102, ig_104, ig_120, \
                         kf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_6 * if__67[k]
                   + pa_z[k] * ig_102[k];

        t_178[k] = f_16 * if__79[k]
                   + pb_y[k] * kf_119[k];

        t_179[k] = f_8 * if__69[k]
                   + pa_z[k] * ig_104[k];

        t_180[k] = f_9 * hg0_75[k]
                   - f_10 * hg1_75[k]
                   + pa_y[k] * ig_120[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pa_z, pb_y, pb_z, hg0_48, hg1_48, if__70, \
                         if__80, if__82, ig_108, kf_120, kf_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_6 * if__80[k]
                   + pb_y[k] * kf_120[k];

        t_182[k] = f_6 * if__70[k]
                   + pb_z[k] * kf_120[k];

        t_183[k] = f_9 * hg0_48[k]
                   - f_10 * hg1_48[k]
                   + pa_z[k] * ig_108[k];

        t_184[k] = f_6 * if__82[k]
                   + pb_y[k] * kf_122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pa_y, pb_x, hg0_80, hg1_80, if__126, \
                         if__127, if__128, ig_125, kf_126, kf_127, \
                         kf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_9 * hg0_80[k]
                   - f_10 * hg1_80[k]
                   + pa_y[k] * ig_125[k];

        t_186[k] = f_16 * if__126[k]
                   + pb_x[k] * kf_126[k];

        t_187[k] = f_16 * if__127[k]
                   + pb_x[k] * kf_127[k];

        t_188[k] = f_16 * if__128[k]
                   + pb_x[k] * kf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pa_x, pb_x, pb_z, hg0_190, hg1_190, if__76, \
                         if__129, ig_190, kf_126, kf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_16 * if__129[k]
                   + pb_x[k] * kf_129[k];

        t_190[k] = f_14 * hg0_190[k]
                   - f_15 * hg1_190[k]
                   + pa_x[k] * ig_190[k];

        t_191[k] = f_6 * if__76[k]
                   + pb_z[k] * kf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_x, pa_y, pb_y, hg0_192, hg0_194, \
                         hg1_192, hg1_194, if__89, ig_135, ig_192, ig_194, \
                         kf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_14 * hg0_192[k]
                   - f_15 * hg1_192[k]
                   + pa_x[k] * ig_192[k];

        t_193[k] = f_6 * if__89[k]
                   + pb_y[k] * kf_129[k];

        t_194[k] = f_14 * hg0_194[k]
                   - f_15 * hg1_194[k]
                   + pa_x[k] * ig_194[k];

        t_195[k] = pa_y[k] * ig_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pa_y, pb_y, if__90, if__91, \
                         if__92, ig_137, ig_138, ig_140, kf_130, \
                         kf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_5 * if__90[k]
                   + pb_y[k] * kf_130[k];

        t_197[k] = pa_y[k] * ig_137[k];

        t_198[k] = f_6 * if__91[k]
                   + pa_y[k] * ig_138[k];

        t_199[k] = f_5 * if__92[k]
                   + pb_y[k] * kf_132[k];

        t_200[k] = pa_y[k] * ig_140[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, pa_y, pb_x, if__96, if__136, \
                         if__137, if__138, ig_144, ig_145, kf_136, kf_137, \
                         kf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_16 * if__136[k]
                   + pb_x[k] * kf_136[k];

        t_202[k] = f_16 * if__137[k]
                   + pb_x[k] * kf_137[k];

        t_203[k] = f_16 * if__138[k]
                   + pb_x[k] * kf_138[k];

        t_204[k] = pa_y[k] * ig_144[k];

        t_205[k] = f_8 * if__96[k]
                   + pa_y[k] * ig_145[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_y, pb_y, pb_z, if__86, if__98, if__99, \
                         ig_147, ig_149, kf_136, kf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_16 * if__86[k]
                   + pb_z[k] * kf_136[k];

        t_207[k] = f_6 * if__98[k]
                   + pa_y[k] * ig_147[k];

        t_208[k] = f_5 * if__99[k]
                   + pb_y[k] * kf_139[k];

        t_209[k] = pa_y[k] * ig_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pa_z, pb_y, pb_z, hg0_75, hg1_75, if__90, \
                         ig_135, kd0_84, kd1_84, kf_140, kf_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_17 * hg0_75[k]
                   - f_18 * hg1_75[k]
                   + pa_z[k] * ig_135[k];

        t_211[k] = pb_y[k] * kf_140[k];

        t_212[k] = f_8 * if__90[k]
                   + pb_z[k] * kf_140[k];

        t_213[k] = f_3 * kd0_84[k]
                   - f_4 * kd1_84[k]
                   + pb_y[k] * kf_141[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, pb_x, pb_y, if__145, if__146, \
                         if__147, kd0_89, kd1_89, kf_142, kf_145, kf_146, \
                         kf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = pb_y[k] * kf_142[k];

        t_215[k] = f_16 * if__145[k]
                   + f_3 * kd0_89[k]
                   - f_4 * kd1_89[k]
                   + pb_x[k] * kf_145[k];

        t_216[k] = f_16 * if__146[k]
                   + pb_x[k] * kf_146[k];

        t_217[k] = f_16 * if__147[k]
                   + pb_x[k] * kf_147[k];

        t_218[k] = pb_y[k] * kf_145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pb_x, pb_y, pb_z, if__96, if__149, \
                         kd0_87, kd0_89, kd1_87, kd1_89, kf_146, kf_148, \
                         kf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_16 * if__149[k]
                   + pb_x[k] * kf_149[k];

        t_220[k] = f_1 * kd0_87[k]
                   - f_2 * kd1_87[k]
                   + pb_y[k] * kf_146[k];

        t_221[k] = f_8 * if__96[k]
                   + pb_z[k] * kf_146[k];

        t_222[k] = f_3 * kd0_89[k]
                   - f_4 * kd1_89[k]
                   + pb_y[k] * kf_148[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_x, pa_y, pb_y, hg0_90, hg0_224, \
                         hg1_90, hg1_224, if__100, ig_150, ig_224, kf_149, \
                         kf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = pb_y[k] * kf_149[k];

        t_224[k] = f_14 * hg0_224[k]
                   - f_15 * hg1_224[k]
                   + pa_x[k] * ig_224[k];

        t_225[k] = f_12 * hg0_90[k]
                   - f_13 * hg1_90[k]
                   + pa_y[k] * ig_150[k];

        t_226[k] = f_11 * if__100[k]
                   + pb_y[k] * kf_150[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, pb_x, pb_z, if__153, kd0_90, kd0_93, \
                         kd1_90, kd1_93, kf_150, kf_151, kf_152, \
                         kf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = pb_z[k] * kf_150[k];

        t_228[k] = f_6 * if__153[k]
                   + f_3 * kd0_93[k]
                   - f_4 * kd1_93[k]
                   + pb_x[k] * kf_153[k];

        t_229[k] = pb_z[k] * kf_151[k];

        t_230[k] = f_3 * kd0_90[k]
                   - f_4 * kd1_90[k]
                   + pb_z[k] * kf_152[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pb_x, pb_z, if__156, if__158, if__159, \
                         kf_153, kf_156, kf_158, kf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_6 * if__156[k]
                   + pb_x[k] * kf_156[k];

        t_232[k] = pb_z[k] * kf_153[k];

        t_233[k] = f_6 * if__158[k]
                   + pb_x[k] * kf_158[k];

        t_234[k] = f_6 * if__159[k]
                   + pb_x[k] * kf_159[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pa_x, pb_y, pb_z, hg0_235, hg1_235, \
                         if__109, ig_235, kd0_93, kd1_93, kf_156, kf_157, \
                         kf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_9 * hg0_235[k]
                   - f_10 * hg1_235[k]
                   + pa_x[k] * ig_235[k];

        t_236[k] = pb_z[k] * kf_156[k];

        t_237[k] = f_3 * kd0_93[k]
                   - f_4 * kd1_93[k]
                   + pb_z[k] * kf_157[k];

        t_238[k] = f_11 * if__109[k]
                   + pb_y[k] * kf_159[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pa_z, pb_z, if__100, ig_150, \
                         ig_151, ig_153, kd0_95, kd1_95, kf_159, \
                         kf_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * kd0_95[k]
                   - f_2 * kd1_95[k]
                   + pb_z[k] * kf_159[k];

        t_240[k] = pa_z[k] * ig_150[k];

        t_241[k] = pa_z[k] * ig_151[k];

        t_242[k] = f_5 * if__100[k]
                   + pb_z[k] * kf_160[k];

        t_243[k] = pa_z[k] * ig_153[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_z, pb_x, pb_y, if__102, if__112, \
                         if__167, ig_155, ig_156, kf_162, kf_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_8 * if__112[k]
                   + pb_y[k] * kf_162[k];

        t_245[k] = f_6 * if__102[k]
                   + pa_z[k] * ig_155[k];

        t_246[k] = pa_z[k] * ig_156[k];

        t_247[k] = f_6 * if__167[k]
                   + pb_x[k] * kf_167[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_z, pb_x, pb_z, if__106, if__168, \
                         if__169, ig_160, kf_166, kf_168, kf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_6 * if__168[k]
                   + pb_x[k] * kf_168[k];

        t_249[k] = f_6 * if__169[k]
                   + pb_x[k] * kf_169[k];

        t_250[k] = pa_z[k] * ig_160[k];

        t_251[k] = f_5 * if__106[k]
                   + pb_z[k] * kf_166[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_y, pa_z, pb_y, hg0_120, hg1_120, \
                         if__107, if__109, if__119, ig_162, ig_164, ig_180, \
                         kf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_6 * if__107[k]
                   + pa_z[k] * ig_162[k];

        t_253[k] = f_8 * if__119[k]
                   + pb_y[k] * kf_169[k];

        t_254[k] = f_8 * if__109[k]
                   + pa_z[k] * ig_164[k];

        t_255[k] = f_14 * hg0_120[k]
                   - f_15 * hg1_120[k]
                   + pa_y[k] * ig_180[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pa_z, pb_y, pb_z, hg0_93, hg1_93, \
                         if__110, if__120, if__122, ig_168, kf_170, \
                         kf_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_16 * if__120[k]
                   + pb_y[k] * kf_170[k];

        t_257[k] = f_6 * if__110[k]
                   + pb_z[k] * kf_170[k];

        t_258[k] = f_9 * hg0_93[k]
                   - f_10 * hg1_93[k]
                   + pa_z[k] * ig_168[k];

        t_259[k] = f_16 * if__122[k]
                   + pb_y[k] * kf_172[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pb_x, hg0_125, hg1_125, if__176, \
                         if__177, if__178, ig_185, kf_176, kf_177, \
                         kf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * hg0_125[k]
                   - f_15 * hg1_125[k]
                   + pa_y[k] * ig_185[k];

        t_261[k] = f_6 * if__176[k]
                   + pb_x[k] * kf_176[k];

        t_262[k] = f_6 * if__177[k]
                   + pb_x[k] * kf_177[k];

        t_263[k] = f_6 * if__178[k]
                   + pb_x[k] * kf_178[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pa_x, pb_x, pb_z, hg0_265, hg1_265, if__116, \
                         if__179, ig_265, kf_176, kf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_6 * if__179[k]
                   + pb_x[k] * kf_179[k];

        t_265[k] = f_9 * hg0_265[k]
                   - f_10 * hg1_265[k]
                   + pa_x[k] * ig_265[k];

        t_266[k] = f_6 * if__116[k]
                   + pb_z[k] * kf_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pa_x, pb_y, hg0_267, hg0_269, hg1_267, hg1_269, \
                         if__129, ig_267, ig_269, kf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_9 * hg0_267[k]
                   - f_10 * hg1_267[k]
                   + pa_x[k] * ig_267[k];

        t_268[k] = f_16 * if__129[k]
                   + pb_y[k] * kf_179[k];

        t_269[k] = f_9 * hg0_269[k]
                   - f_10 * hg1_269[k]
                   + pa_x[k] * ig_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pa_y, pb_y, pb_z, hg0_135, hg1_135, if__120, \
                         if__130, ig_195, kf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_9 * hg0_135[k]
                   - f_10 * hg1_135[k]
                   + pa_y[k] * ig_195[k];

        t_271[k] = f_6 * if__130[k]
                   + pb_y[k] * kf_180[k];

        t_272[k] = f_16 * if__120[k]
                   + pb_z[k] * kf_180[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pa_y, pa_z, pb_y, hg0_108, hg0_140, hg1_108, \
                         hg1_140, if__132, ig_183, ig_200, kf_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_14 * hg0_108[k]
                   - f_15 * hg1_108[k]
                   + pa_z[k] * ig_183[k];

        t_274[k] = f_6 * if__132[k]
                   + pb_y[k] * kf_182[k];

        t_275[k] = f_9 * hg0_140[k]
                   - f_10 * hg1_140[k]
                   + pa_y[k] * ig_200[k];
    }
}

static auto
compute_prim_kg_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t hg0,
                                            const size_t hg1, const size_t if_, const size_t ig,
                                            const size_t kd0, const size_t kd1, const size_t kf,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 3.0 / p;
    const auto f_8 = 2.0 / p;
    const auto f_9 = 0.5 / alpha;
    const auto f_10 = 0.5 * beta / (alpha * p);
    const auto f_11 = 2.5 / p;
    const auto f_12 = 2.0 / alpha;
    const auto f_13 = 2.0 * beta / (alpha * p);
    const auto f_16 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hg0_135 = buffer.data(hg0 + 135);
    const auto *hg0_280 = buffer.data(hg0 + 280);
    const auto *hg0_282 = buffer.data(hg0 + 282);
    const auto *hg0_284 = buffer.data(hg0 + 284);
    const auto *hg0_314 = buffer.data(hg0 + 314);

    const auto *hg1_135 = buffer.data(hg1 + 135);
    const auto *hg1_280 = buffer.data(hg1 + 280);
    const auto *hg1_282 = buffer.data(hg1 + 282);
    const auto *hg1_284 = buffer.data(hg1 + 284);
    const auto *hg1_314 = buffer.data(hg1 + 314);

    const auto *if__126 = buffer.data(if_ + 126);
    const auto *if__136 = buffer.data(if_ + 136);
    const auto *if__139 = buffer.data(if_ + 139);
    const auto *if__140 = buffer.data(if_ + 140);
    const auto *if__141 = buffer.data(if_ + 141);
    const auto *if__142 = buffer.data(if_ + 142);
    const auto *if__146 = buffer.data(if_ + 146);
    const auto *if__148 = buffer.data(if_ + 148);
    const auto *if__149 = buffer.data(if_ + 149);
    const auto *if__150 = buffer.data(if_ + 150);
    const auto *if__160 = buffer.data(if_ + 160);
    const auto *if__162 = buffer.data(if_ + 162);
    const auto *if__170 = buffer.data(if_ + 170);
    const auto *if__172 = buffer.data(if_ + 172);
    const auto *if__180 = buffer.data(if_ + 180);
    const auto *if__182 = buffer.data(if_ + 182);
    const auto *if__186 = buffer.data(if_ + 186);
    const auto *if__187 = buffer.data(if_ + 187);
    const auto *if__188 = buffer.data(if_ + 188);
    const auto *if__189 = buffer.data(if_ + 189);
    const auto *if__190 = buffer.data(if_ + 190);
    const auto *if__192 = buffer.data(if_ + 192);
    const auto *if__196 = buffer.data(if_ + 196);
    const auto *if__197 = buffer.data(if_ + 197);
    const auto *if__198 = buffer.data(if_ + 198);
    const auto *if__200 = buffer.data(if_ + 200);
    const auto *if__202 = buffer.data(if_ + 202);
    const auto *if__205 = buffer.data(if_ + 205);
    const auto *if__206 = buffer.data(if_ + 206);
    const auto *if__207 = buffer.data(if_ + 207);
    const auto *if__209 = buffer.data(if_ + 209);
    const auto *if__210 = buffer.data(if_ + 210);
    const auto *if__213 = buffer.data(if_ + 213);
    const auto *if__215 = buffer.data(if_ + 215);
    const auto *if__216 = buffer.data(if_ + 216);
    const auto *if__218 = buffer.data(if_ + 218);
    const auto *if__219 = buffer.data(if_ + 219);
    const auto *if__225 = buffer.data(if_ + 225);
    const auto *if__227 = buffer.data(if_ + 227);
    const auto *if__228 = buffer.data(if_ + 228);
    const auto *if__229 = buffer.data(if_ + 229);
    const auto *if__230 = buffer.data(if_ + 230);
    const auto *if__233 = buffer.data(if_ + 233);
    const auto *if__235 = buffer.data(if_ + 235);
    const auto *if__236 = buffer.data(if_ + 236);
    const auto *if__237 = buffer.data(if_ + 237);
    const auto *if__238 = buffer.data(if_ + 238);
    const auto *if__239 = buffer.data(if_ + 239);
    const auto *if__240 = buffer.data(if_ + 240);
    const auto *if__243 = buffer.data(if_ + 243);
    const auto *if__245 = buffer.data(if_ + 245);
    const auto *if__246 = buffer.data(if_ + 246);
    const auto *if__247 = buffer.data(if_ + 247);
    const auto *if__248 = buffer.data(if_ + 248);
    const auto *if__249 = buffer.data(if_ + 249);
    const auto *if__250 = buffer.data(if_ + 250);
    const auto *if__253 = buffer.data(if_ + 253);
    const auto *if__255 = buffer.data(if_ + 255);
    const auto *if__256 = buffer.data(if_ + 256);
    const auto *if__257 = buffer.data(if_ + 257);
    const auto *if__258 = buffer.data(if_ + 258);
    const auto *if__259 = buffer.data(if_ + 259);
    const auto *if__263 = buffer.data(if_ + 263);
    const auto *if__266 = buffer.data(if_ + 266);
    const auto *if__267 = buffer.data(if_ + 267);
    const auto *if__268 = buffer.data(if_ + 268);
    const auto *if__270 = buffer.data(if_ + 270);
    const auto *if__273 = buffer.data(if_ + 273);
    const auto *if__275 = buffer.data(if_ + 275);
    const auto *if__276 = buffer.data(if_ + 276);
    const auto *if__277 = buffer.data(if_ + 277);
    const auto *if__279 = buffer.data(if_ + 279);

    const auto *ig_210 = buffer.data(ig + 210);
    const auto *ig_212 = buffer.data(ig + 212);
    const auto *ig_213 = buffer.data(ig + 213);
    const auto *ig_215 = buffer.data(ig + 215);
    const auto *ig_219 = buffer.data(ig + 219);
    const auto *ig_220 = buffer.data(ig + 220);
    const auto *ig_222 = buffer.data(ig + 222);
    const auto *ig_224 = buffer.data(ig + 224);
    const auto *ig_225 = buffer.data(ig + 225);
    const auto *ig_226 = buffer.data(ig + 226);
    const auto *ig_228 = buffer.data(ig + 228);
    const auto *ig_231 = buffer.data(ig + 231);
    const auto *ig_280 = buffer.data(ig + 280);
    const auto *ig_282 = buffer.data(ig + 282);
    const auto *ig_284 = buffer.data(ig + 284);
    const auto *ig_300 = buffer.data(ig + 300);
    const auto *ig_302 = buffer.data(ig + 302);
    const auto *ig_305 = buffer.data(ig + 305);
    const auto *ig_309 = buffer.data(ig + 309);
    const auto *ig_314 = buffer.data(ig + 314);
    const auto *ig_315 = buffer.data(ig + 315);
    const auto *ig_318 = buffer.data(ig + 318);
    const auto *ig_320 = buffer.data(ig + 320);
    const auto *ig_325 = buffer.data(ig + 325);
    const auto *ig_327 = buffer.data(ig + 327);
    const auto *ig_328 = buffer.data(ig + 328);
    const auto *ig_329 = buffer.data(ig + 329);
    const auto *ig_335 = buffer.data(ig + 335);
    const auto *ig_340 = buffer.data(ig + 340);
    const auto *ig_341 = buffer.data(ig + 341);
    const auto *ig_342 = buffer.data(ig + 342);
    const auto *ig_343 = buffer.data(ig + 343);
    const auto *ig_344 = buffer.data(ig + 344);
    const auto *ig_345 = buffer.data(ig + 345);
    const auto *ig_348 = buffer.data(ig + 348);
    const auto *ig_350 = buffer.data(ig + 350);
    const auto *ig_355 = buffer.data(ig + 355);
    const auto *ig_356 = buffer.data(ig + 356);
    const auto *ig_357 = buffer.data(ig + 357);
    const auto *ig_358 = buffer.data(ig + 358);
    const auto *ig_359 = buffer.data(ig + 359);
    const auto *ig_360 = buffer.data(ig + 360);
    const auto *ig_363 = buffer.data(ig + 363);
    const auto *ig_365 = buffer.data(ig + 365);
    const auto *ig_370 = buffer.data(ig + 370);
    const auto *ig_371 = buffer.data(ig + 371);
    const auto *ig_372 = buffer.data(ig + 372);
    const auto *ig_373 = buffer.data(ig + 373);
    const auto *ig_374 = buffer.data(ig + 374);
    const auto *ig_375 = buffer.data(ig + 375);
    const auto *ig_378 = buffer.data(ig + 378);
    const auto *ig_380 = buffer.data(ig + 380);
    const auto *ig_385 = buffer.data(ig + 385);
    const auto *ig_386 = buffer.data(ig + 386);
    const auto *ig_387 = buffer.data(ig + 387);
    const auto *ig_388 = buffer.data(ig + 388);
    const auto *ig_389 = buffer.data(ig + 389);
    const auto *ig_393 = buffer.data(ig + 393);
    const auto *ig_400 = buffer.data(ig + 400);
    const auto *ig_401 = buffer.data(ig + 401);
    const auto *ig_402 = buffer.data(ig + 402);
    const auto *ig_403 = buffer.data(ig + 403);
    const auto *ig_404 = buffer.data(ig + 404);
    const auto *ig_405 = buffer.data(ig + 405);
    const auto *ig_408 = buffer.data(ig + 408);
    const auto *ig_410 = buffer.data(ig + 410);
    const auto *ig_415 = buffer.data(ig + 415);
    const auto *ig_416 = buffer.data(ig + 416);
    const auto *ig_417 = buffer.data(ig + 417);
    const auto *ig_419 = buffer.data(ig + 419);

    const auto *kd0_120 = buffer.data(kd0 + 120);
    const auto *kd0_123 = buffer.data(kd0 + 123);
    const auto *kd0_125 = buffer.data(kd0 + 125);
    const auto *kd0_168 = buffer.data(kd0 + 168);
    const auto *kd0_171 = buffer.data(kd0 + 171);
    const auto *kd0_173 = buffer.data(kd0 + 173);

    const auto *kd1_120 = buffer.data(kd1 + 120);
    const auto *kd1_123 = buffer.data(kd1 + 123);
    const auto *kd1_125 = buffer.data(kd1 + 125);
    const auto *kd1_168 = buffer.data(kd1 + 168);
    const auto *kd1_171 = buffer.data(kd1 + 171);
    const auto *kd1_173 = buffer.data(kd1 + 173);

    const auto *kf_186 = buffer.data(kf + 186);
    const auto *kf_187 = buffer.data(kf + 187);
    const auto *kf_188 = buffer.data(kf + 188);
    const auto *kf_189 = buffer.data(kf + 189);
    const auto *kf_190 = buffer.data(kf + 190);
    const auto *kf_192 = buffer.data(kf + 192);
    const auto *kf_196 = buffer.data(kf + 196);
    const auto *kf_197 = buffer.data(kf + 197);
    const auto *kf_198 = buffer.data(kf + 198);
    const auto *kf_199 = buffer.data(kf + 199);
    const auto *kf_200 = buffer.data(kf + 200);
    const auto *kf_201 = buffer.data(kf + 201);
    const auto *kf_202 = buffer.data(kf + 202);
    const auto *kf_205 = buffer.data(kf + 205);
    const auto *kf_206 = buffer.data(kf + 206);
    const auto *kf_207 = buffer.data(kf + 207);
    const auto *kf_208 = buffer.data(kf + 208);
    const auto *kf_209 = buffer.data(kf + 209);
    const auto *kf_210 = buffer.data(kf + 210);
    const auto *kf_211 = buffer.data(kf + 211);
    const auto *kf_213 = buffer.data(kf + 213);
    const auto *kf_216 = buffer.data(kf + 216);
    const auto *kf_218 = buffer.data(kf + 218);
    const auto *kf_219 = buffer.data(kf + 219);
    const auto *kf_220 = buffer.data(kf + 220);
    const auto *kf_222 = buffer.data(kf + 222);
    const auto *kf_227 = buffer.data(kf + 227);
    const auto *kf_228 = buffer.data(kf + 228);
    const auto *kf_229 = buffer.data(kf + 229);
    const auto *kf_230 = buffer.data(kf + 230);
    const auto *kf_232 = buffer.data(kf + 232);
    const auto *kf_236 = buffer.data(kf + 236);
    const auto *kf_237 = buffer.data(kf + 237);
    const auto *kf_238 = buffer.data(kf + 238);
    const auto *kf_239 = buffer.data(kf + 239);
    const auto *kf_240 = buffer.data(kf + 240);
    const auto *kf_242 = buffer.data(kf + 242);
    const auto *kf_246 = buffer.data(kf + 246);
    const auto *kf_247 = buffer.data(kf + 247);
    const auto *kf_248 = buffer.data(kf + 248);
    const auto *kf_249 = buffer.data(kf + 249);
    const auto *kf_250 = buffer.data(kf + 250);
    const auto *kf_252 = buffer.data(kf + 252);
    const auto *kf_256 = buffer.data(kf + 256);
    const auto *kf_257 = buffer.data(kf + 257);
    const auto *kf_258 = buffer.data(kf + 258);
    const auto *kf_259 = buffer.data(kf + 259);
    const auto *kf_260 = buffer.data(kf + 260);
    const auto *kf_262 = buffer.data(kf + 262);
    const auto *kf_266 = buffer.data(kf + 266);
    const auto *kf_267 = buffer.data(kf + 267);
    const auto *kf_268 = buffer.data(kf + 268);
    const auto *kf_270 = buffer.data(kf + 270);
    const auto *kf_272 = buffer.data(kf + 272);
    const auto *kf_275 = buffer.data(kf + 275);
    const auto *kf_276 = buffer.data(kf + 276);
    const auto *kf_277 = buffer.data(kf + 277);
    const auto *kf_279 = buffer.data(kf + 279);
    const auto *kf_280 = buffer.data(kf + 280);
    const auto *kf_281 = buffer.data(kf + 281);
    const auto *kf_283 = buffer.data(kf + 283);
    const auto *kf_285 = buffer.data(kf + 285);
    const auto *kf_286 = buffer.data(kf + 286);
    const auto *kf_287 = buffer.data(kf + 287);
    const auto *kf_288 = buffer.data(kf + 288);
    const auto *kf_289 = buffer.data(kf + 289);

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pb_x, if__186, if__187, if__188, if__189, \
                         kf_186, kf_187, kf_188, kf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * if__186[k]
                   + pb_x[k] * kf_186[k];

        t_277[k] = f_6 * if__187[k]
                   + pb_x[k] * kf_187[k];

        t_278[k] = f_6 * if__188[k]
                   + pb_x[k] * kf_188[k];

        t_279[k] = f_6 * if__189[k]
                   + pb_x[k] * kf_189[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pa_x, pb_z, hg0_280, hg0_282, hg1_280, hg1_282, \
                         if__126, ig_280, ig_282, kf_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_9 * hg0_280[k]
                   - f_10 * hg1_280[k]
                   + pa_x[k] * ig_280[k];

        t_281[k] = f_16 * if__126[k]
                   + pb_z[k] * kf_186[k];

        t_282[k] = f_9 * hg0_282[k]
                   - f_10 * hg1_282[k]
                   + pa_x[k] * ig_282[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_x, pa_y, pb_y, hg0_284, hg1_284, \
                         if__139, if__140, ig_210, ig_284, kf_189, \
                         kf_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_6 * if__139[k]
                   + pb_y[k] * kf_189[k];

        t_284[k] = f_9 * hg0_284[k]
                   - f_10 * hg1_284[k]
                   + pa_x[k] * ig_284[k];

        t_285[k] = pa_y[k] * ig_210[k];

        t_286[k] = f_5 * if__140[k]
                   + pb_y[k] * kf_190[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, pa_y, pb_x, pb_y, if__141, \
                         if__142, if__196, ig_212, ig_213, ig_215, kf_192, \
                         kf_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = pa_y[k] * ig_212[k];

        t_288[k] = f_6 * if__141[k]
                   + pa_y[k] * ig_213[k];

        t_289[k] = f_5 * if__142[k]
                   + pb_y[k] * kf_192[k];

        t_290[k] = pa_y[k] * ig_215[k];

        t_291[k] = f_6 * if__196[k]
                   + pb_x[k] * kf_196[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pa_y, pb_x, if__146, if__197, if__198, \
                         ig_219, ig_220, kf_197, kf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_6 * if__197[k]
                   + pb_x[k] * kf_197[k];

        t_293[k] = f_6 * if__198[k]
                   + pb_x[k] * kf_198[k];

        t_294[k] = pa_y[k] * ig_219[k];

        t_295[k] = f_8 * if__146[k]
                   + pa_y[k] * ig_220[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pa_y, pb_y, pb_z, if__136, if__148, \
                         if__149, ig_222, ig_224, kf_196, kf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_8 * if__136[k]
                   + pb_z[k] * kf_196[k];

        t_297[k] = f_6 * if__148[k]
                   + pa_y[k] * ig_222[k];

        t_298[k] = f_5 * if__149[k]
                   + pb_y[k] * kf_199[k];

        t_299[k] = pa_y[k] * ig_224[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_z, pb_y, pb_z, hg0_135, hg1_135, \
                         if__140, ig_210, kd0_120, kd1_120, kf_200, \
                         kf_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_12 * hg0_135[k]
                   - f_13 * hg1_135[k]
                   + pa_z[k] * ig_210[k];

        t_301[k] = pb_y[k] * kf_200[k];

        t_302[k] = f_11 * if__140[k]
                   + pb_z[k] * kf_200[k];

        t_303[k] = f_3 * kd0_120[k]
                   - f_4 * kd1_120[k]
                   + pb_y[k] * kf_201[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pb_x, pb_y, if__205, if__206, \
                         if__207, kd0_125, kd1_125, kf_202, kf_205, kf_206, \
                         kf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = pb_y[k] * kf_202[k];

        t_305[k] = f_6 * if__205[k]
                   + f_3 * kd0_125[k]
                   - f_4 * kd1_125[k]
                   + pb_x[k] * kf_205[k];

        t_306[k] = f_6 * if__206[k]
                   + pb_x[k] * kf_206[k];

        t_307[k] = f_6 * if__207[k]
                   + pb_x[k] * kf_207[k];

        t_308[k] = pb_y[k] * kf_205[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pb_x, pb_y, pb_z, if__146, if__209, \
                         kd0_123, kd0_125, kd1_123, kd1_125, kf_206, kf_208, \
                         kf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_6 * if__209[k]
                   + pb_x[k] * kf_209[k];

        t_310[k] = f_1 * kd0_123[k]
                   - f_2 * kd1_123[k]
                   + pb_y[k] * kf_206[k];

        t_311[k] = f_11 * if__146[k]
                   + pb_z[k] * kf_206[k];

        t_312[k] = f_3 * kd0_125[k]
                   - f_4 * kd1_125[k]
                   + pb_y[k] * kf_208[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, t_317, pa_x, pb_y, pb_z, hg0_314, \
                         hg1_314, if__150, if__210, ig_314, ig_315, kf_209, \
                         kf_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = pb_y[k] * kf_209[k];

        t_314[k] = f_9 * hg0_314[k]
                   - f_10 * hg1_314[k]
                   + pa_x[k] * ig_314[k];

        t_315[k] = f_8 * if__210[k]
                   + pa_x[k] * ig_315[k];

        t_316[k] = f_7 * if__150[k]
                   + pb_y[k] * kf_210[k];

        t_317[k] = pb_z[k] * kf_210[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, pa_x, pb_x, pb_z, if__213, \
                         if__215, if__216, ig_318, ig_320, kf_211, kf_213, \
                         kf_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_6 * if__213[k]
                   + pa_x[k] * ig_318[k];

        t_319[k] = pb_z[k] * kf_211[k];

        t_320[k] = f_6 * if__215[k]
                   + pa_x[k] * ig_320[k];

        t_321[k] = f_5 * if__216[k]
                   + pb_x[k] * kf_216[k];

        t_322[k] = pb_z[k] * kf_213[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, pa_x, pb_x, pb_z, if__218, \
                         if__219, ig_325, ig_327, kf_216, kf_218, \
                         kf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_5 * if__218[k]
                   + pb_x[k] * kf_218[k];

        t_324[k] = f_5 * if__219[k]
                   + pb_x[k] * kf_219[k];

        t_325[k] = pa_x[k] * ig_325[k];

        t_326[k] = pb_z[k] * kf_216[k];

        t_327[k] = pa_x[k] * ig_327[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, t_332, t_333, pa_x, pa_z, pb_z, if__150, \
                         ig_225, ig_226, ig_228, ig_328, ig_329, \
                         kf_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = pa_x[k] * ig_328[k];

        t_329[k] = pa_x[k] * ig_329[k];

        t_330[k] = pa_z[k] * ig_225[k];

        t_331[k] = pa_z[k] * ig_226[k];

        t_332[k] = f_5 * if__150[k]
                   + pb_z[k] * kf_220[k];

        t_333[k] = pa_z[k] * ig_228[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pa_x, pa_z, pb_x, pb_y, if__162, if__225, \
                         if__227, ig_231, ig_335, kf_222, kf_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_11 * if__162[k]
                   + pb_y[k] * kf_222[k];

        t_335[k] = f_6 * if__225[k]
                   + pa_x[k] * ig_335[k];

        t_336[k] = pa_z[k] * ig_231[k];

        t_337[k] = f_5 * if__227[k]
                   + pb_x[k] * kf_227[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, t_342, t_343, pa_x, pb_x, if__228, \
                         if__229, ig_340, ig_341, ig_342, ig_343, kf_228, \
                         kf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_5 * if__228[k]
                   + pb_x[k] * kf_228[k];

        t_339[k] = f_5 * if__229[k]
                   + pb_x[k] * kf_229[k];

        t_340[k] = pa_x[k] * ig_340[k];

        t_341[k] = pa_x[k] * ig_341[k];

        t_342[k] = pa_x[k] * ig_342[k];

        t_343[k] = pa_x[k] * ig_343[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, pa_x, pb_y, pb_z, if__160, \
                         if__170, if__230, if__233, ig_344, ig_345, ig_348, \
                         kf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = pa_x[k] * ig_344[k];

        t_345[k] = f_8 * if__230[k]
                   + pa_x[k] * ig_345[k];

        t_346[k] = f_8 * if__170[k]
                   + pb_y[k] * kf_230[k];

        t_347[k] = f_6 * if__160[k]
                   + pb_z[k] * kf_230[k];

        t_348[k] = f_6 * if__233[k]
                   + pa_x[k] * ig_348[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pa_x, pb_x, pb_y, if__172, if__235, \
                         if__236, if__237, ig_350, kf_232, kf_236, \
                         kf_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_8 * if__172[k]
                   + pb_y[k] * kf_232[k];

        t_350[k] = f_6 * if__235[k]
                   + pa_x[k] * ig_350[k];

        t_351[k] = f_5 * if__236[k]
                   + pb_x[k] * kf_236[k];

        t_352[k] = f_5 * if__237[k]
                   + pb_x[k] * kf_237[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, t_357, t_358, pa_x, pb_x, if__238, \
                         if__239, ig_355, ig_356, ig_357, ig_358, kf_238, \
                         kf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_5 * if__238[k]
                   + pb_x[k] * kf_238[k];

        t_354[k] = f_5 * if__239[k]
                   + pb_x[k] * kf_239[k];

        t_355[k] = pa_x[k] * ig_355[k];

        t_356[k] = pa_x[k] * ig_356[k];

        t_357[k] = pa_x[k] * ig_357[k];

        t_358[k] = pa_x[k] * ig_358[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, t_363, pa_x, pb_y, pb_z, if__170, \
                         if__180, if__240, if__243, ig_359, ig_360, ig_363, \
                         kf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = pa_x[k] * ig_359[k];

        t_360[k] = f_8 * if__240[k]
                   + pa_x[k] * ig_360[k];

        t_361[k] = f_16 * if__180[k]
                   + pb_y[k] * kf_240[k];

        t_362[k] = f_16 * if__170[k]
                   + pb_z[k] * kf_240[k];

        t_363[k] = f_6 * if__243[k]
                   + pa_x[k] * ig_363[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, pa_x, pb_x, pb_y, if__182, if__245, \
                         if__246, if__247, ig_365, kf_242, kf_246, \
                         kf_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_16 * if__182[k]
                   + pb_y[k] * kf_242[k];

        t_365[k] = f_6 * if__245[k]
                   + pa_x[k] * ig_365[k];

        t_366[k] = f_5 * if__246[k]
                   + pb_x[k] * kf_246[k];

        t_367[k] = f_5 * if__247[k]
                   + pb_x[k] * kf_247[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, t_372, t_373, pa_x, pb_x, if__248, \
                         if__249, ig_370, ig_371, ig_372, ig_373, kf_248, \
                         kf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_5 * if__248[k]
                   + pb_x[k] * kf_248[k];

        t_369[k] = f_5 * if__249[k]
                   + pb_x[k] * kf_249[k];

        t_370[k] = pa_x[k] * ig_370[k];

        t_371[k] = pa_x[k] * ig_371[k];

        t_372[k] = pa_x[k] * ig_372[k];

        t_373[k] = pa_x[k] * ig_373[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, t_378, pa_x, pb_y, pb_z, if__180, \
                         if__190, if__250, if__253, ig_374, ig_375, ig_378, \
                         kf_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_x[k] * ig_374[k];

        t_375[k] = f_8 * if__250[k]
                   + pa_x[k] * ig_375[k];

        t_376[k] = f_6 * if__190[k]
                   + pb_y[k] * kf_250[k];

        t_377[k] = f_8 * if__180[k]
                   + pb_z[k] * kf_250[k];

        t_378[k] = f_6 * if__253[k]
                   + pa_x[k] * ig_378[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pa_x, pb_x, pb_y, if__192, if__255, \
                         if__256, if__257, ig_380, kf_252, kf_256, \
                         kf_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_6 * if__192[k]
                   + pb_y[k] * kf_252[k];

        t_380[k] = f_6 * if__255[k]
                   + pa_x[k] * ig_380[k];

        t_381[k] = f_5 * if__256[k]
                   + pb_x[k] * kf_256[k];

        t_382[k] = f_5 * if__257[k]
                   + pb_x[k] * kf_257[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, t_388, pa_x, pb_x, if__258, \
                         if__259, ig_385, ig_386, ig_387, ig_388, kf_258, \
                         kf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_5 * if__258[k]
                   + pb_x[k] * kf_258[k];

        t_384[k] = f_5 * if__259[k]
                   + pb_x[k] * kf_259[k];

        t_385[k] = pa_x[k] * ig_385[k];

        t_386[k] = pa_x[k] * ig_386[k];

        t_387[k] = pa_x[k] * ig_387[k];

        t_388[k] = pa_x[k] * ig_388[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pa_x, pa_y, pb_y, if__200, \
                         if__263, ig_300, ig_302, ig_389, ig_393, \
                         kf_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = pa_x[k] * ig_389[k];

        t_390[k] = pa_y[k] * ig_300[k];

        t_391[k] = f_5 * if__200[k]
                   + pb_y[k] * kf_260[k];

        t_392[k] = pa_y[k] * ig_302[k];

        t_393[k] = f_6 * if__263[k]
                   + pa_x[k] * ig_393[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pa_y, pb_x, pb_y, if__202, if__266, \
                         if__267, ig_305, kf_262, kf_266, kf_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_5 * if__202[k]
                   + pb_y[k] * kf_262[k];

        t_395[k] = pa_y[k] * ig_305[k];

        t_396[k] = f_5 * if__266[k]
                   + pb_x[k] * kf_266[k];

        t_397[k] = f_5 * if__267[k]
                   + pb_x[k] * kf_267[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, t_402, t_403, pa_x, pa_y, pb_x, if__268, \
                         ig_309, ig_400, ig_401, ig_402, ig_403, \
                         kf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_5 * if__268[k]
                   + pb_x[k] * kf_268[k];

        t_399[k] = pa_y[k] * ig_309[k];

        t_400[k] = pa_x[k] * ig_400[k];

        t_401[k] = pa_x[k] * ig_401[k];

        t_402[k] = pa_x[k] * ig_402[k];

        t_403[k] = pa_x[k] * ig_403[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, t_408, pa_x, pb_y, pb_z, if__200, \
                         if__270, if__273, ig_404, ig_405, ig_408, \
                         kf_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_x[k] * ig_404[k];

        t_405[k] = f_8 * if__270[k]
                   + pa_x[k] * ig_405[k];

        t_406[k] = pb_y[k] * kf_270[k];

        t_407[k] = f_7 * if__200[k]
                   + pb_z[k] * kf_270[k];

        t_408[k] = f_6 * if__273[k]
                   + pa_x[k] * ig_408[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pa_x, pb_x, pb_y, if__275, \
                         if__276, if__277, ig_410, kf_272, kf_275, kf_276, \
                         kf_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = pb_y[k] * kf_272[k];

        t_410[k] = f_6 * if__275[k]
                   + pa_x[k] * ig_410[k];

        t_411[k] = f_5 * if__276[k]
                   + pb_x[k] * kf_276[k];

        t_412[k] = f_5 * if__277[k]
                   + pb_x[k] * kf_277[k];

        t_413[k] = pb_y[k] * kf_275[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, t_419, pa_x, pb_x, pb_y, if__279, \
                         ig_415, ig_416, ig_417, ig_419, kf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_5 * if__279[k]
                   + pb_x[k] * kf_279[k];

        t_415[k] = pa_x[k] * ig_415[k];

        t_416[k] = pa_x[k] * ig_416[k];

        t_417[k] = pa_x[k] * ig_417[k];

        t_418[k] = pb_y[k] * kf_279[k];

        t_419[k] = pa_x[k] * ig_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, pb_x, pb_y, pb_z, if__210, \
                         kd0_168, kd0_171, kd1_168, kd1_171, kf_280, kf_281, \
                         kf_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_1 * kd0_168[k]
                   - f_2 * kd1_168[k]
                   + pb_x[k] * kf_280[k];

        t_421[k] = f_0 * if__210[k]
                   + pb_y[k] * kf_280[k];

        t_422[k] = pb_z[k] * kf_280[k];

        t_423[k] = f_3 * kd0_171[k]
                   - f_4 * kd1_171[k]
                   + pb_x[k] * kf_283[k];

        t_424[k] = pb_z[k] * kf_281[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, pb_x, kd0_173, kd1_173, kf_285, \
                         kf_286, kf_287, kf_288, kf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_3 * kd0_173[k]
                   - f_4 * kd1_173[k]
                   + pb_x[k] * kf_285[k];

        t_426[k] = pb_x[k] * kf_286[k];

        t_427[k] = pb_x[k] * kf_287[k];

        t_428[k] = pb_x[k] * kf_288[k];

        t_429[k] = pb_x[k] * kf_289[k];
    }
}

static auto
compute_prim_kg_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t hg0,
                                            const size_t hg1, const size_t if_, const size_t ig,
                                            const size_t kd0, const size_t kd1, const size_t kf,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 3.0 / p;
    const auto f_8 = 2.0 / p;
    const auto f_9 = 0.5 / alpha;
    const auto f_10 = 0.5 * beta / (alpha * p);
    const auto f_11 = 2.5 / p;
    const auto f_12 = 2.0 / alpha;
    const auto f_13 = 2.0 * beta / (alpha * p);
    const auto f_14 = 1.0 / alpha;
    const auto f_15 = beta / (alpha * p);
    const auto f_16 = 1.5 / p;
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hg0_235 = buffer.data(hg0 + 235);
    const auto *hg0_250 = buffer.data(hg0 + 250);
    const auto *hg0_265 = buffer.data(hg0 + 265);
    const auto *hg0_269 = buffer.data(hg0 + 269);
    const auto *hg0_280 = buffer.data(hg0 + 280);
    const auto *hg0_284 = buffer.data(hg0 + 284);
    const auto *hg0_299 = buffer.data(hg0 + 299);
    const auto *hg0_314 = buffer.data(hg0 + 314);

    const auto *hg1_235 = buffer.data(hg1 + 235);
    const auto *hg1_250 = buffer.data(hg1 + 250);
    const auto *hg1_265 = buffer.data(hg1 + 265);
    const auto *hg1_269 = buffer.data(hg1 + 269);
    const auto *hg1_280 = buffer.data(hg1 + 280);
    const auto *hg1_284 = buffer.data(hg1 + 284);
    const auto *hg1_299 = buffer.data(hg1 + 299);
    const auto *hg1_314 = buffer.data(hg1 + 314);

    const auto *if__210 = buffer.data(if_ + 210);
    const auto *if__212 = buffer.data(if_ + 212);
    const auto *if__216 = buffer.data(if_ + 216);
    const auto *if__217 = buffer.data(if_ + 217);
    const auto *if__219 = buffer.data(if_ + 219);
    const auto *if__220 = buffer.data(if_ + 220);
    const auto *if__222 = buffer.data(if_ + 222);
    const auto *if__226 = buffer.data(if_ + 226);
    const auto *if__229 = buffer.data(if_ + 229);
    const auto *if__230 = buffer.data(if_ + 230);
    const auto *if__232 = buffer.data(if_ + 232);
    const auto *if__236 = buffer.data(if_ + 236);
    const auto *if__238 = buffer.data(if_ + 238);
    const auto *if__239 = buffer.data(if_ + 239);
    const auto *if__240 = buffer.data(if_ + 240);
    const auto *if__242 = buffer.data(if_ + 242);
    const auto *if__246 = buffer.data(if_ + 246);
    const auto *if__248 = buffer.data(if_ + 248);
    const auto *if__249 = buffer.data(if_ + 249);
    const auto *if__250 = buffer.data(if_ + 250);
    const auto *if__252 = buffer.data(if_ + 252);
    const auto *if__256 = buffer.data(if_ + 256);
    const auto *if__258 = buffer.data(if_ + 258);
    const auto *if__259 = buffer.data(if_ + 259);
    const auto *if__260 = buffer.data(if_ + 260);
    const auto *if__262 = buffer.data(if_ + 262);
    const auto *if__266 = buffer.data(if_ + 266);
    const auto *if__268 = buffer.data(if_ + 268);
    const auto *if__269 = buffer.data(if_ + 269);
    const auto *if__270 = buffer.data(if_ + 270);
    const auto *if__271 = buffer.data(if_ + 271);
    const auto *if__272 = buffer.data(if_ + 272);
    const auto *if__276 = buffer.data(if_ + 276);
    const auto *if__278 = buffer.data(if_ + 278);
    const auto *if__279 = buffer.data(if_ + 279);

    const auto *ig_315 = buffer.data(ig + 315);
    const auto *ig_316 = buffer.data(ig + 316);
    const auto *ig_318 = buffer.data(ig + 318);
    const auto *ig_320 = buffer.data(ig + 320);
    const auto *ig_325 = buffer.data(ig + 325);
    const auto *ig_327 = buffer.data(ig + 327);
    const auto *ig_329 = buffer.data(ig + 329);
    const auto *ig_340 = buffer.data(ig + 340);
    const auto *ig_355 = buffer.data(ig + 355);
    const auto *ig_359 = buffer.data(ig + 359);
    const auto *ig_370 = buffer.data(ig + 370);
    const auto *ig_374 = buffer.data(ig + 374);
    const auto *ig_385 = buffer.data(ig + 385);
    const auto *ig_389 = buffer.data(ig + 389);
    const auto *ig_404 = buffer.data(ig + 404);
    const auto *ig_405 = buffer.data(ig + 405);
    const auto *ig_407 = buffer.data(ig + 407);
    const auto *ig_408 = buffer.data(ig + 408);
    const auto *ig_410 = buffer.data(ig + 410);
    const auto *ig_415 = buffer.data(ig + 415);
    const auto *ig_417 = buffer.data(ig + 417);
    const auto *ig_419 = buffer.data(ig + 419);

    const auto *kd0_171 = buffer.data(kd0 + 171);
    const auto *kd0_173 = buffer.data(kd0 + 173);
    const auto *kd0_180 = buffer.data(kd0 + 180);
    const auto *kd0_183 = buffer.data(kd0 + 183);
    const auto *kd0_185 = buffer.data(kd0 + 185);
    const auto *kd0_186 = buffer.data(kd0 + 186);
    const auto *kd0_189 = buffer.data(kd0 + 189);
    const auto *kd0_191 = buffer.data(kd0 + 191);
    const auto *kd0_192 = buffer.data(kd0 + 192);
    const auto *kd0_195 = buffer.data(kd0 + 195);
    const auto *kd0_197 = buffer.data(kd0 + 197);
    const auto *kd0_198 = buffer.data(kd0 + 198);
    const auto *kd0_201 = buffer.data(kd0 + 201);
    const auto *kd0_203 = buffer.data(kd0 + 203);
    const auto *kd0_210 = buffer.data(kd0 + 210);
    const auto *kd0_213 = buffer.data(kd0 + 213);
    const auto *kd0_215 = buffer.data(kd0 + 215);

    const auto *kd1_171 = buffer.data(kd1 + 171);
    const auto *kd1_173 = buffer.data(kd1 + 173);
    const auto *kd1_180 = buffer.data(kd1 + 180);
    const auto *kd1_183 = buffer.data(kd1 + 183);
    const auto *kd1_185 = buffer.data(kd1 + 185);
    const auto *kd1_186 = buffer.data(kd1 + 186);
    const auto *kd1_189 = buffer.data(kd1 + 189);
    const auto *kd1_191 = buffer.data(kd1 + 191);
    const auto *kd1_192 = buffer.data(kd1 + 192);
    const auto *kd1_195 = buffer.data(kd1 + 195);
    const auto *kd1_197 = buffer.data(kd1 + 197);
    const auto *kd1_198 = buffer.data(kd1 + 198);
    const auto *kd1_201 = buffer.data(kd1 + 201);
    const auto *kd1_203 = buffer.data(kd1 + 203);
    const auto *kd1_210 = buffer.data(kd1 + 210);
    const auto *kd1_213 = buffer.data(kd1 + 213);
    const auto *kd1_215 = buffer.data(kd1 + 215);

    const auto *kf_286 = buffer.data(kf + 286);
    const auto *kf_287 = buffer.data(kf + 287);
    const auto *kf_289 = buffer.data(kf + 289);
    const auto *kf_290 = buffer.data(kf + 290);
    const auto *kf_292 = buffer.data(kf + 292);
    const auto *kf_296 = buffer.data(kf + 296);
    const auto *kf_297 = buffer.data(kf + 297);
    const auto *kf_298 = buffer.data(kf + 298);
    const auto *kf_299 = buffer.data(kf + 299);
    const auto *kf_300 = buffer.data(kf + 300);
    const auto *kf_302 = buffer.data(kf + 302);
    const auto *kf_303 = buffer.data(kf + 303);
    const auto *kf_305 = buffer.data(kf + 305);
    const auto *kf_306 = buffer.data(kf + 306);
    const auto *kf_307 = buffer.data(kf + 307);
    const auto *kf_308 = buffer.data(kf + 308);
    const auto *kf_309 = buffer.data(kf + 309);
    const auto *kf_310 = buffer.data(kf + 310);
    const auto *kf_312 = buffer.data(kf + 312);
    const auto *kf_313 = buffer.data(kf + 313);
    const auto *kf_315 = buffer.data(kf + 315);
    const auto *kf_316 = buffer.data(kf + 316);
    const auto *kf_317 = buffer.data(kf + 317);
    const auto *kf_318 = buffer.data(kf + 318);
    const auto *kf_319 = buffer.data(kf + 319);
    const auto *kf_320 = buffer.data(kf + 320);
    const auto *kf_322 = buffer.data(kf + 322);
    const auto *kf_323 = buffer.data(kf + 323);
    const auto *kf_325 = buffer.data(kf + 325);
    const auto *kf_326 = buffer.data(kf + 326);
    const auto *kf_327 = buffer.data(kf + 327);
    const auto *kf_328 = buffer.data(kf + 328);
    const auto *kf_329 = buffer.data(kf + 329);
    const auto *kf_330 = buffer.data(kf + 330);
    const auto *kf_332 = buffer.data(kf + 332);
    const auto *kf_333 = buffer.data(kf + 333);
    const auto *kf_335 = buffer.data(kf + 335);
    const auto *kf_336 = buffer.data(kf + 336);
    const auto *kf_337 = buffer.data(kf + 337);
    const auto *kf_338 = buffer.data(kf + 338);
    const auto *kf_339 = buffer.data(kf + 339);
    const auto *kf_340 = buffer.data(kf + 340);
    const auto *kf_342 = buffer.data(kf + 342);
    const auto *kf_346 = buffer.data(kf + 346);
    const auto *kf_347 = buffer.data(kf + 347);
    const auto *kf_348 = buffer.data(kf + 348);
    const auto *kf_349 = buffer.data(kf + 349);
    const auto *kf_350 = buffer.data(kf + 350);
    const auto *kf_352 = buffer.data(kf + 352);
    const auto *kf_353 = buffer.data(kf + 353);
    const auto *kf_355 = buffer.data(kf + 355);
    const auto *kf_356 = buffer.data(kf + 356);
    const auto *kf_357 = buffer.data(kf + 357);
    const auto *kf_358 = buffer.data(kf + 358);
    const auto *kf_359 = buffer.data(kf + 359);

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, pb_y, pb_z, if__216, if__219, \
                         kd0_171, kd0_173, kd1_171, kd1_173, kf_286, kf_287, \
                         kf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_0 * if__216[k]
                   + f_1 * kd0_171[k]
                   - f_2 * kd1_171[k]
                   + pb_y[k] * kf_286[k];

        t_431[k] = pb_z[k] * kf_286[k];

        t_432[k] = f_3 * kd0_171[k]
                   - f_4 * kd1_171[k]
                   + pb_z[k] * kf_287[k];

        t_433[k] = f_0 * if__219[k]
                   + pb_y[k] * kf_289[k];

        t_434[k] = f_1 * kd0_173[k]
                   - f_2 * kd1_173[k]
                   + pb_z[k] * kf_289[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, pa_z, pb_y, pb_z, if__210, \
                         if__222, ig_315, ig_316, ig_318, kf_290, \
                         kf_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = pa_z[k] * ig_315[k];

        t_436[k] = pa_z[k] * ig_316[k];

        t_437[k] = f_5 * if__210[k]
                   + pb_z[k] * kf_290[k];

        t_438[k] = pa_z[k] * ig_318[k];

        t_439[k] = f_7 * if__222[k]
                   + pb_y[k] * kf_292[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, t_445, pa_z, pb_x, if__212, \
                         ig_320, ig_325, kf_296, kf_297, kf_298, \
                         kf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_6 * if__212[k]
                   + pa_z[k] * ig_320[k];

        t_441[k] = pb_x[k] * kf_296[k];

        t_442[k] = pb_x[k] * kf_297[k];

        t_443[k] = pb_x[k] * kf_298[k];

        t_444[k] = pb_x[k] * kf_299[k];

        t_445[k] = pa_z[k] * ig_325[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, t_449, pa_z, pb_y, pb_z, if__216, if__217, \
                         if__219, if__229, ig_327, ig_329, kf_296, \
                         kf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_5 * if__216[k]
                   + pb_z[k] * kf_296[k];

        t_447[k] = f_6 * if__217[k]
                   + pa_z[k] * ig_327[k];

        t_448[k] = f_7 * if__229[k]
                   + pb_y[k] * kf_299[k];

        t_449[k] = f_8 * if__219[k]
                   + pa_z[k] * ig_329[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, pb_x, pb_y, pb_z, if__220, if__230, \
                         kd0_180, kd0_183, kd1_180, kd1_183, kf_300, \
                         kf_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_1 * kd0_180[k]
                   - f_2 * kd1_180[k]
                   + pb_x[k] * kf_300[k];

        t_451[k] = f_11 * if__230[k]
                   + pb_y[k] * kf_300[k];

        t_452[k] = f_6 * if__220[k]
                   + pb_z[k] * kf_300[k];

        t_453[k] = f_3 * kd0_183[k]
                   - f_4 * kd1_183[k]
                   + pb_x[k] * kf_303[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, pb_x, pb_y, if__232, kd0_185, \
                         kd1_185, kf_302, kf_305, kf_306, kf_307, \
                         kf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = f_11 * if__232[k]
                   + pb_y[k] * kf_302[k];

        t_455[k] = f_3 * kd0_185[k]
                   - f_4 * kd1_185[k]
                   + pb_x[k] * kf_305[k];

        t_456[k] = pb_x[k] * kf_306[k];

        t_457[k] = pb_x[k] * kf_307[k];

        t_458[k] = pb_x[k] * kf_308[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pa_z, pb_x, pb_z, hg0_235, hg1_235, if__226, \
                         ig_340, kf_306, kf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = pb_x[k] * kf_309[k];

        t_460[k] = f_9 * hg0_235[k]
                   - f_10 * hg1_235[k]
                   + pa_z[k] * ig_340[k];

        t_461[k] = f_6 * if__226[k]
                   + pb_z[k] * kf_306[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pa_y, pb_y, hg0_269, hg1_269, if__238, if__239, \
                         ig_359, kd0_185, kd1_185, kf_308, kf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_11 * if__238[k]
                   + f_3 * kd0_185[k]
                   - f_4 * kd1_185[k]
                   + pb_y[k] * kf_308[k];

        t_463[k] = f_11 * if__239[k]
                   + pb_y[k] * kf_309[k];

        t_464[k] = f_12 * hg0_269[k]
                   - f_13 * hg1_269[k]
                   + pa_y[k] * ig_359[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, pb_x, pb_y, pb_z, if__230, if__240, \
                         kd0_186, kd0_189, kd1_186, kd1_189, kf_310, \
                         kf_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_1 * kd0_186[k]
                   - f_2 * kd1_186[k]
                   + pb_x[k] * kf_310[k];

        t_466[k] = f_8 * if__240[k]
                   + pb_y[k] * kf_310[k];

        t_467[k] = f_16 * if__230[k]
                   + pb_z[k] * kf_310[k];

        t_468[k] = f_3 * kd0_189[k]
                   - f_4 * kd1_189[k]
                   + pb_x[k] * kf_313[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, t_473, pb_x, pb_y, if__242, kd0_191, \
                         kd1_191, kf_312, kf_315, kf_316, kf_317, \
                         kf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_8 * if__242[k]
                   + pb_y[k] * kf_312[k];

        t_470[k] = f_3 * kd0_191[k]
                   - f_4 * kd1_191[k]
                   + pb_x[k] * kf_315[k];

        t_471[k] = pb_x[k] * kf_316[k];

        t_472[k] = pb_x[k] * kf_317[k];

        t_473[k] = pb_x[k] * kf_318[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pa_z, pb_x, pb_z, hg0_250, hg1_250, if__236, \
                         ig_355, kf_316, kf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = pb_x[k] * kf_319[k];

        t_475[k] = f_14 * hg0_250[k]
                   - f_15 * hg1_250[k]
                   + pa_z[k] * ig_355[k];

        t_476[k] = f_16 * if__236[k]
                   + pb_z[k] * kf_316[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, pa_y, pb_y, hg0_284, hg1_284, if__248, if__249, \
                         ig_374, kd0_191, kd1_191, kf_318, kf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_8 * if__248[k]
                   + f_3 * kd0_191[k]
                   - f_4 * kd1_191[k]
                   + pb_y[k] * kf_318[k];

        t_478[k] = f_8 * if__249[k]
                   + pb_y[k] * kf_319[k];

        t_479[k] = f_17 * hg0_284[k]
                   - f_18 * hg1_284[k]
                   + pa_y[k] * ig_374[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, pb_x, pb_y, pb_z, if__240, if__250, \
                         kd0_192, kd0_195, kd1_192, kd1_195, kf_320, \
                         kf_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = f_1 * kd0_192[k]
                   - f_2 * kd1_192[k]
                   + pb_x[k] * kf_320[k];

        t_481[k] = f_16 * if__250[k]
                   + pb_y[k] * kf_320[k];

        t_482[k] = f_8 * if__240[k]
                   + pb_z[k] * kf_320[k];

        t_483[k] = f_3 * kd0_195[k]
                   - f_4 * kd1_195[k]
                   + pb_x[k] * kf_323[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, t_488, pb_x, pb_y, if__252, kd0_197, \
                         kd1_197, kf_322, kf_325, kf_326, kf_327, \
                         kf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_16 * if__252[k]
                   + pb_y[k] * kf_322[k];

        t_485[k] = f_3 * kd0_197[k]
                   - f_4 * kd1_197[k]
                   + pb_x[k] * kf_325[k];

        t_486[k] = pb_x[k] * kf_326[k];

        t_487[k] = pb_x[k] * kf_327[k];

        t_488[k] = pb_x[k] * kf_328[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pa_z, pb_x, pb_z, hg0_265, hg1_265, if__246, \
                         ig_370, kf_326, kf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = pb_x[k] * kf_329[k];

        t_490[k] = f_17 * hg0_265[k]
                   - f_18 * hg1_265[k]
                   + pa_z[k] * ig_370[k];

        t_491[k] = f_8 * if__246[k]
                   + pb_z[k] * kf_326[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, pa_y, pb_y, hg0_299, hg1_299, if__258, if__259, \
                         ig_389, kd0_197, kd1_197, kf_328, kf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_16 * if__258[k]
                   + f_3 * kd0_197[k]
                   - f_4 * kd1_197[k]
                   + pb_y[k] * kf_328[k];

        t_493[k] = f_16 * if__259[k]
                   + pb_y[k] * kf_329[k];

        t_494[k] = f_14 * hg0_299[k]
                   - f_15 * hg1_299[k]
                   + pa_y[k] * ig_389[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pb_x, pb_y, pb_z, if__250, if__260, \
                         kd0_198, kd0_201, kd1_198, kd1_201, kf_330, \
                         kf_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_1 * kd0_198[k]
                   - f_2 * kd1_198[k]
                   + pb_x[k] * kf_330[k];

        t_496[k] = f_6 * if__260[k]
                   + pb_y[k] * kf_330[k];

        t_497[k] = f_11 * if__250[k]
                   + pb_z[k] * kf_330[k];

        t_498[k] = f_3 * kd0_201[k]
                   - f_4 * kd1_201[k]
                   + pb_x[k] * kf_333[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, t_503, pb_x, pb_y, if__262, kd0_203, \
                         kd1_203, kf_332, kf_335, kf_336, kf_337, \
                         kf_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_6 * if__262[k]
                   + pb_y[k] * kf_332[k];

        t_500[k] = f_3 * kd0_203[k]
                   - f_4 * kd1_203[k]
                   + pb_x[k] * kf_335[k];

        t_501[k] = pb_x[k] * kf_336[k];

        t_502[k] = pb_x[k] * kf_337[k];

        t_503[k] = pb_x[k] * kf_338[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, pa_z, pb_x, pb_z, hg0_280, hg1_280, if__256, \
                         ig_385, kf_336, kf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = pb_x[k] * kf_339[k];

        t_505[k] = f_12 * hg0_280[k]
                   - f_13 * hg1_280[k]
                   + pa_z[k] * ig_385[k];

        t_506[k] = f_11 * if__256[k]
                   + pb_z[k] * kf_336[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, t_510, pa_y, pb_y, hg0_314, hg1_314, if__268, \
                         if__269, ig_404, ig_405, kd0_203, kd1_203, kf_338, \
                         kf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_6 * if__268[k]
                   + f_3 * kd0_203[k]
                   - f_4 * kd1_203[k]
                   + pb_y[k] * kf_338[k];

        t_508[k] = f_6 * if__269[k]
                   + pb_y[k] * kf_339[k];

        t_509[k] = f_9 * hg0_314[k]
                   - f_10 * hg1_314[k]
                   + pa_y[k] * ig_404[k];

        t_510[k] = pa_y[k] * ig_405[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, t_515, pa_y, pb_y, if__270, if__271, \
                         if__272, ig_407, ig_408, ig_410, kf_340, \
                         kf_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_5 * if__270[k]
                   + pb_y[k] * kf_340[k];

        t_512[k] = pa_y[k] * ig_407[k];

        t_513[k] = f_6 * if__271[k]
                   + pa_y[k] * ig_408[k];

        t_514[k] = f_5 * if__272[k]
                   + pb_y[k] * kf_342[k];

        t_515[k] = pa_y[k] * ig_410[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, t_520, t_521, pa_y, pb_x, pb_z, if__266, \
                         if__276, ig_415, kf_346, kf_347, kf_348, \
                         kf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = pb_x[k] * kf_346[k];

        t_517[k] = pb_x[k] * kf_347[k];

        t_518[k] = pb_x[k] * kf_348[k];

        t_519[k] = pb_x[k] * kf_349[k];

        t_520[k] = f_8 * if__276[k]
                   + pa_y[k] * ig_415[k];

        t_521[k] = f_7 * if__266[k]
                   + pb_z[k] * kf_346[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, pa_y, pb_x, pb_y, if__278, \
                         if__279, ig_417, ig_419, kd0_210, kd1_210, kf_349, \
                         kf_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_6 * if__278[k]
                   + pa_y[k] * ig_417[k];

        t_523[k] = f_5 * if__279[k]
                   + pb_y[k] * kf_349[k];

        t_524[k] = pa_y[k] * ig_419[k];

        t_525[k] = f_1 * kd0_210[k]
                   - f_2 * kd1_210[k]
                   + pb_x[k] * kf_350[k];

        t_526[k] = pb_y[k] * kf_350[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, pb_x, pb_y, pb_z, if__270, kd0_213, \
                         kd0_215, kd1_213, kd1_215, kf_350, kf_352, kf_353, \
                         kf_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_0 * if__270[k]
                   + pb_z[k] * kf_350[k];

        t_528[k] = f_3 * kd0_213[k]
                   - f_4 * kd1_213[k]
                   + pb_x[k] * kf_353[k];

        t_529[k] = pb_y[k] * kf_352[k];

        t_530[k] = f_3 * kd0_215[k]
                   - f_4 * kd1_215[k]
                   + pb_x[k] * kf_355[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, t_535, t_536, pb_x, pb_y, pb_z, if__276, \
                         kd0_213, kd1_213, kf_356, kf_357, kf_358, \
                         kf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = pb_x[k] * kf_356[k];

        t_532[k] = pb_x[k] * kf_357[k];

        t_533[k] = pb_x[k] * kf_358[k];

        t_534[k] = pb_x[k] * kf_359[k];

        t_535[k] = f_1 * kd0_213[k]
                   - f_2 * kd1_213[k]
                   + pb_y[k] * kf_356[k];

        t_536[k] = f_0 * if__276[k]
                   + pb_z[k] * kf_356[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, pb_y, pb_z, if__279, kd0_215, kd1_215, kf_358, \
                         kf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_3 * kd0_215[k]
                   - f_4 * kd1_215[k]
                   + pb_y[k] * kf_358[k];

        t_538[k] = pb_y[k] * kf_359[k];

        t_539[k] = f_0 * if__279[k]
                   + f_1 * kd0_215[k]
                   - f_2 * kd1_215[k]
                   + pb_z[k] * kf_359[k];
    }
}

auto
compute_prim_kg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hg0, const size_t hg1,
                                     const size_t if_, const size_t ig, const size_t kd0,
                                     const size_t kd1, const size_t kf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    compute_prim_kg_electron_repulsion_0_piece0(buffer, target, pa, pb, hg0, hg1, if_, ig, kd0,
                                                kd1, kf, ncols, alpha, beta, p);

    compute_prim_kg_electron_repulsion_0_piece1(buffer, target, pa, pb, hg0, hg1, if_, ig, kd0,
                                                kd1, kf, ncols, alpha, beta, p);

    compute_prim_kg_electron_repulsion_0_piece2(buffer, target, pa, pb, hg0, hg1, if_, ig, kd0,
                                                kd1, kf, ncols, alpha, beta, p);

    compute_prim_kg_electron_repulsion_0_piece3(buffer, target, pa, pb, hg0, hg1, if_, ig, kd0,
                                                kd1, kf, ncols, alpha, beta, p);
}

}  // namespace simdt2ceri
