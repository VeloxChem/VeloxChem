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


#include "SimdElectronRepulsionVrrRecKK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_kk_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t hk0,
                                            const size_t hk1, const size_t ii, const size_t ik,
                                            const size_t kh0, const size_t kh1, const size_t ki,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
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
    const auto f_16 = 3.0 / p;
    const auto f_17 = 0.5 / alpha;
    const auto f_18 = 0.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / alpha;
    const auto f_20 = 2.0 * beta / (alpha * p);

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

    const auto *hk0_0 = buffer.data(hk0 + 0);
    const auto *hk0_136 = buffer.data(hk0 + 136);

    const auto *hk1_0 = buffer.data(hk1 + 0);
    const auto *hk1_136 = buffer.data(hk1 + 136);

    const auto *ii_0 = buffer.data(ii + 0);
    const auto *ii_1 = buffer.data(ii + 1);
    const auto *ii_2 = buffer.data(ii + 2);
    const auto *ii_3 = buffer.data(ii + 3);
    const auto *ii_5 = buffer.data(ii + 5);
    const auto *ii_6 = buffer.data(ii + 6);
    const auto *ii_7 = buffer.data(ii + 7);
    const auto *ii_8 = buffer.data(ii + 8);
    const auto *ii_9 = buffer.data(ii + 9);
    const auto *ii_10 = buffer.data(ii + 10);
    const auto *ii_11 = buffer.data(ii + 11);
    const auto *ii_12 = buffer.data(ii + 12);
    const auto *ii_13 = buffer.data(ii + 13);
    const auto *ii_14 = buffer.data(ii + 14);
    const auto *ii_21 = buffer.data(ii + 21);
    const auto *ii_22 = buffer.data(ii + 22);
    const auto *ii_23 = buffer.data(ii + 23);
    const auto *ii_24 = buffer.data(ii + 24);
    const auto *ii_25 = buffer.data(ii + 25);
    const auto *ii_26 = buffer.data(ii + 26);
    const auto *ii_27 = buffer.data(ii + 27);
    const auto *ii_28 = buffer.data(ii + 28);
    const auto *ii_33 = buffer.data(ii + 33);
    const auto *ii_37 = buffer.data(ii + 37);
    const auto *ii_42 = buffer.data(ii + 42);
    const auto *ii_49 = buffer.data(ii + 49);
    const auto *ii_51 = buffer.data(ii + 51);
    const auto *ii_52 = buffer.data(ii + 52);
    const auto *ii_53 = buffer.data(ii + 53);
    const auto *ii_54 = buffer.data(ii + 54);
    const auto *ii_78 = buffer.data(ii + 78);
    const auto *ii_79 = buffer.data(ii + 79);
    const auto *ii_80 = buffer.data(ii + 80);
    const auto *ii_81 = buffer.data(ii + 81);
    const auto *ii_83 = buffer.data(ii + 83);
    const auto *ii_87 = buffer.data(ii + 87);
    const auto *ii_90 = buffer.data(ii + 90);
    const auto *ii_94 = buffer.data(ii + 94);
    const auto *ii_99 = buffer.data(ii + 99);
    const auto *ii_105 = buffer.data(ii + 105);
    const auto *ii_107 = buffer.data(ii + 107);
    const auto *ii_108 = buffer.data(ii + 108);
    const auto *ii_109 = buffer.data(ii + 109);
    const auto *ii_110 = buffer.data(ii + 110);
    const auto *ii_111 = buffer.data(ii + 111);

    const auto *ik_0 = buffer.data(ik + 0);
    const auto *ik_3 = buffer.data(ik + 3);
    const auto *ik_5 = buffer.data(ik + 5);
    const auto *ik_6 = buffer.data(ik + 6);
    const auto *ik_9 = buffer.data(ik + 9);
    const auto *ik_10 = buffer.data(ik + 10);
    const auto *ik_12 = buffer.data(ik + 12);
    const auto *ik_14 = buffer.data(ik + 14);
    const auto *ik_15 = buffer.data(ik + 15);
    const auto *ik_17 = buffer.data(ik + 17);
    const auto *ik_18 = buffer.data(ik + 18);
    const auto *ik_20 = buffer.data(ik + 20);
    const auto *ik_21 = buffer.data(ik + 21);
    const auto *ik_27 = buffer.data(ik + 27);
    const auto *ik_28 = buffer.data(ik + 28);
    const auto *ik_30 = buffer.data(ik + 30);
    const auto *ik_31 = buffer.data(ik + 31);
    const auto *ik_32 = buffer.data(ik + 32);
    const auto *ik_33 = buffer.data(ik + 33);
    const auto *ik_35 = buffer.data(ik + 35);
    const auto *ik_36 = buffer.data(ik + 36);
    const auto *ik_136 = buffer.data(ik + 136);

    const auto *kh0_0 = buffer.data(kh0 + 0);
    const auto *kh0_1 = buffer.data(kh0 + 1);
    const auto *kh0_2 = buffer.data(kh0 + 2);
    const auto *kh0_3 = buffer.data(kh0 + 3);
    const auto *kh0_5 = buffer.data(kh0 + 5);
    const auto *kh0_6 = buffer.data(kh0 + 6);
    const auto *kh0_8 = buffer.data(kh0 + 8);
    const auto *kh0_9 = buffer.data(kh0 + 9);
    const auto *kh0_15 = buffer.data(kh0 + 15);
    const auto *kh0_17 = buffer.data(kh0 + 17);
    const auto *kh0_18 = buffer.data(kh0 + 18);
    const auto *kh0_19 = buffer.data(kh0 + 19);
    const auto *kh0_20 = buffer.data(kh0 + 20);
    const auto *kh0_63 = buffer.data(kh0 + 63);
    const auto *kh0_65 = buffer.data(kh0 + 65);
    const auto *kh0_66 = buffer.data(kh0 + 66);
    const auto *kh0_68 = buffer.data(kh0 + 68);
    const auto *kh0_69 = buffer.data(kh0 + 69);
    const auto *kh0_70 = buffer.data(kh0 + 70);
    const auto *kh0_72 = buffer.data(kh0 + 72);
    const auto *kh0_73 = buffer.data(kh0 + 73);
    const auto *kh0_78 = buffer.data(kh0 + 78);
    const auto *kh0_79 = buffer.data(kh0 + 79);
    const auto *kh0_80 = buffer.data(kh0 + 80);

    const auto *kh1_0 = buffer.data(kh1 + 0);
    const auto *kh1_1 = buffer.data(kh1 + 1);
    const auto *kh1_2 = buffer.data(kh1 + 2);
    const auto *kh1_3 = buffer.data(kh1 + 3);
    const auto *kh1_5 = buffer.data(kh1 + 5);
    const auto *kh1_6 = buffer.data(kh1 + 6);
    const auto *kh1_8 = buffer.data(kh1 + 8);
    const auto *kh1_9 = buffer.data(kh1 + 9);
    const auto *kh1_15 = buffer.data(kh1 + 15);
    const auto *kh1_17 = buffer.data(kh1 + 17);
    const auto *kh1_18 = buffer.data(kh1 + 18);
    const auto *kh1_19 = buffer.data(kh1 + 19);
    const auto *kh1_20 = buffer.data(kh1 + 20);
    const auto *kh1_63 = buffer.data(kh1 + 63);
    const auto *kh1_65 = buffer.data(kh1 + 65);
    const auto *kh1_66 = buffer.data(kh1 + 66);
    const auto *kh1_68 = buffer.data(kh1 + 68);
    const auto *kh1_69 = buffer.data(kh1 + 69);
    const auto *kh1_70 = buffer.data(kh1 + 70);
    const auto *kh1_72 = buffer.data(kh1 + 72);
    const auto *kh1_73 = buffer.data(kh1 + 73);
    const auto *kh1_78 = buffer.data(kh1 + 78);
    const auto *kh1_79 = buffer.data(kh1 + 79);
    const auto *kh1_80 = buffer.data(kh1 + 80);

    const auto *ki_0 = buffer.data(ki + 0);
    const auto *ki_1 = buffer.data(ki + 1);
    const auto *ki_2 = buffer.data(ki + 2);
    const auto *ki_3 = buffer.data(ki + 3);
    const auto *ki_5 = buffer.data(ki + 5);
    const auto *ki_6 = buffer.data(ki + 6);
    const auto *ki_8 = buffer.data(ki + 8);
    const auto *ki_9 = buffer.data(ki + 9);
    const auto *ki_10 = buffer.data(ki + 10);
    const auto *ki_12 = buffer.data(ki + 12);
    const auto *ki_13 = buffer.data(ki + 13);
    const auto *ki_14 = buffer.data(ki + 14);
    const auto *ki_15 = buffer.data(ki + 15);
    const auto *ki_20 = buffer.data(ki + 20);
    const auto *ki_21 = buffer.data(ki + 21);
    const auto *ki_23 = buffer.data(ki + 23);
    const auto *ki_24 = buffer.data(ki + 24);
    const auto *ki_25 = buffer.data(ki + 25);
    const auto *ki_26 = buffer.data(ki + 26);
    const auto *ki_27 = buffer.data(ki + 27);
    const auto *ki_28 = buffer.data(ki + 28);
    const auto *ki_29 = buffer.data(ki + 29);
    const auto *ki_31 = buffer.data(ki + 31);
    const auto *ki_33 = buffer.data(ki + 33);
    const auto *ki_34 = buffer.data(ki + 34);
    const auto *ki_37 = buffer.data(ki + 37);
    const auto *ki_38 = buffer.data(ki + 38);
    const auto *ki_42 = buffer.data(ki + 42);
    const auto *ki_43 = buffer.data(ki + 43);
    const auto *ki_49 = buffer.data(ki + 49);
    const auto *ki_51 = buffer.data(ki + 51);
    const auto *ki_52 = buffer.data(ki + 52);
    const auto *ki_53 = buffer.data(ki + 53);
    const auto *ki_54 = buffer.data(ki + 54);
    const auto *ki_55 = buffer.data(ki + 55);
    const auto *ki_56 = buffer.data(ki + 56);
    const auto *ki_58 = buffer.data(ki + 58);
    const auto *ki_59 = buffer.data(ki + 59);
    const auto *ki_61 = buffer.data(ki + 61);
    const auto *ki_62 = buffer.data(ki + 62);
    const auto *ki_65 = buffer.data(ki + 65);
    const auto *ki_66 = buffer.data(ki + 66);
    const auto *ki_70 = buffer.data(ki + 70);
    const auto *ki_76 = buffer.data(ki + 76);
    const auto *ki_77 = buffer.data(ki + 77);
    const auto *ki_78 = buffer.data(ki + 78);
    const auto *ki_79 = buffer.data(ki + 79);
    const auto *ki_80 = buffer.data(ki + 80);
    const auto *ki_81 = buffer.data(ki + 81);
    const auto *ki_83 = buffer.data(ki + 83);
    const auto *ki_84 = buffer.data(ki + 84);
    const auto *ki_85 = buffer.data(ki + 85);
    const auto *ki_86 = buffer.data(ki + 86);
    const auto *ki_87 = buffer.data(ki + 87);
    const auto *ki_89 = buffer.data(ki + 89);
    const auto *ki_90 = buffer.data(ki + 90);
    const auto *ki_91 = buffer.data(ki + 91);
    const auto *ki_93 = buffer.data(ki + 93);
    const auto *ki_94 = buffer.data(ki + 94);
    const auto *ki_95 = buffer.data(ki + 95);
    const auto *ki_96 = buffer.data(ki + 96);
    const auto *ki_98 = buffer.data(ki + 98);
    const auto *ki_99 = buffer.data(ki + 99);
    const auto *ki_105 = buffer.data(ki + 105);
    const auto *ki_106 = buffer.data(ki + 106);
    const auto *ki_107 = buffer.data(ki + 107);
    const auto *ki_108 = buffer.data(ki + 108);
    const auto *ki_109 = buffer.data(ki + 109);
    const auto *ki_110 = buffer.data(ki + 110);
    const auto *ki_111 = buffer.data(ki + 111);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, ii_0, kh0_0, kh1_0, \
                         ki_0, ki_1, ki_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ii_0[k]
                 + f_1 * kh0_0[k]
                 - f_2 * kh1_0[k]
                 + pb_x[k] * ki_0[k];

        t_1[k] = pb_y[k] * ki_0[k];

        t_2[k] = pb_z[k] * ki_0[k];

        t_3[k] = f_3 * kh0_0[k]
                 - f_4 * kh1_0[k]
                 + pb_y[k] * ki_1[k];

        t_4[k] = pb_y[k] * ki_2[k];

        t_5[k] = f_3 * kh0_0[k]
                 - f_4 * kh1_0[k]
                 + pb_z[k] * ki_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, kh0_1, kh0_2, kh0_3, kh1_1, \
                         kh1_2, kh1_3, ki_3, ki_5, ki_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * kh0_1[k]
                 - f_6 * kh1_1[k]
                 + pb_y[k] * ki_3[k];

        t_7[k] = pb_z[k] * ki_3[k];

        t_8[k] = pb_y[k] * ki_5[k];

        t_9[k] = f_5 * kh0_2[k]
                 - f_6 * kh1_2[k]
                 + pb_z[k] * ki_5[k];

        t_10[k] = f_7 * kh0_3[k]
                  - f_8 * kh1_3[k]
                  + pb_y[k] * ki_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, kh0_5, kh0_6, kh1_5, \
                         kh1_6, ki_6, ki_8, ki_9, ki_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * ki_6[k];

        t_12[k] = f_3 * kh0_5[k]
                  - f_4 * kh1_5[k]
                  + pb_y[k] * ki_8[k];

        t_13[k] = pb_y[k] * ki_9[k];

        t_14[k] = f_7 * kh0_5[k]
                  - f_8 * kh1_5[k]
                  + pb_z[k] * ki_9[k];

        t_15[k] = f_9 * kh0_6[k]
                  - f_10 * kh1_6[k]
                  + pb_y[k] * ki_10[k];

        t_16[k] = pb_z[k] * ki_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, kh0_8, kh0_9, kh1_8, kh1_9, \
                         ki_12, ki_13, ki_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * kh0_8[k]
                  - f_6 * kh1_8[k]
                  + pb_y[k] * ki_12[k];

        t_18[k] = f_3 * kh0_9[k]
                  - f_4 * kh1_9[k]
                  + pb_y[k] * ki_13[k];

        t_19[k] = pb_y[k] * ki_14[k];

        t_20[k] = f_9 * kh0_9[k]
                  - f_10 * kh1_9[k]
                  + pb_z[k] * ki_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_z, ii_21, ii_23, ii_24, ii_25, \
                         ki_15, ki_21, ki_23, ki_24, ki_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * ii_21[k]
                  + pb_x[k] * ki_21[k];

        t_22[k] = pb_z[k] * ki_15[k];

        t_23[k] = f_0 * ii_23[k]
                  + pb_x[k] * ki_23[k];

        t_24[k] = f_0 * ii_24[k]
                  + pb_x[k] * ki_24[k];

        t_25[k] = f_0 * ii_25[k]
                  + pb_x[k] * ki_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, ii_27, kh0_15, kh1_15, \
                         ki_20, ki_21, ki_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * ki_20[k];

        t_27[k] = f_0 * ii_27[k]
                  + pb_x[k] * ki_27[k];

        t_28[k] = f_1 * kh0_15[k]
                  - f_2 * kh1_15[k]
                  + pb_y[k] * ki_21[k];

        t_29[k] = pb_z[k] * ki_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, kh0_17, kh0_18, kh0_19, kh1_17, kh1_18, \
                         kh1_19, ki_23, ki_24, ki_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_9 * kh0_17[k]
                  - f_10 * kh1_17[k]
                  + pb_y[k] * ki_23[k];

        t_31[k] = f_7 * kh0_18[k]
                  - f_8 * kh1_18[k]
                  + pb_y[k] * ki_24[k];

        t_32[k] = f_5 * kh0_19[k]
                  - f_6 * kh1_19[k]
                  + pb_y[k] * ki_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, ii_0, ik_0, \
                         kh0_20, kh1_20, ki_26, ki_27, ki_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * kh0_20[k]
                  - f_4 * kh1_20[k]
                  + pb_y[k] * ki_26[k];

        t_34[k] = pb_y[k] * ki_27[k];

        t_35[k] = f_1 * kh0_20[k]
                  - f_2 * kh1_20[k]
                  + pb_z[k] * ki_27[k];

        t_36[k] = pa_y[k] * ik_0[k];

        t_37[k] = f_11 * ii_0[k]
                  + pb_y[k] * ki_28[k];

        t_38[k] = pb_z[k] * ki_28[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pb_z, ii_1, ii_3, ik_3, ik_5, \
                         ik_6, ki_29, ki_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * ii_1[k]
                  + pa_y[k] * ik_3[k];

        t_40[k] = pb_z[k] * ki_29[k];

        t_41[k] = pa_y[k] * ik_5[k];

        t_42[k] = f_13 * ii_3[k]
                  + pa_y[k] * ik_6[k];

        t_43[k] = pb_z[k] * ki_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_y, pb_y, pb_z, ii_5, ii_6, ii_8, \
                         ik_9, ik_10, ik_12, ki_33, ki_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_11 * ii_5[k]
                  + pb_y[k] * ki_33[k];

        t_45[k] = pa_y[k] * ik_9[k];

        t_46[k] = f_14 * ii_6[k]
                  + pa_y[k] * ik_10[k];

        t_47[k] = pb_z[k] * ki_34[k];

        t_48[k] = f_12 * ii_8[k]
                  + pa_y[k] * ik_12[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_y, pb_y, pb_z, ii_9, ii_10, ii_12, \
                         ik_14, ik_15, ik_17, ki_37, ki_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_11 * ii_9[k]
                  + pb_y[k] * ki_37[k];

        t_50[k] = pa_y[k] * ik_14[k];

        t_51[k] = f_15 * ii_10[k]
                  + pa_y[k] * ik_15[k];

        t_52[k] = pb_z[k] * ki_38[k];

        t_53[k] = f_13 * ii_12[k]
                  + pa_y[k] * ik_17[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pb_x, pb_y, ii_13, ii_14, ii_49, ik_18, \
                         ik_20, ki_42, ki_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * ii_13[k]
                  + pa_y[k] * ik_18[k];

        t_55[k] = f_11 * ii_14[k]
                  + pb_y[k] * ki_42[k];

        t_56[k] = pa_y[k] * ik_20[k];

        t_57[k] = f_16 * ii_49[k]
                  + pb_x[k] * ki_49[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pb_z, ii_51, ii_52, ii_53, ii_54, \
                         ki_43, ki_51, ki_52, ki_53, ki_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_z[k] * ki_43[k];

        t_59[k] = f_16 * ii_51[k]
                  + pb_x[k] * ki_51[k];

        t_60[k] = f_16 * ii_52[k]
                  + pb_x[k] * ki_52[k];

        t_61[k] = f_16 * ii_53[k]
                  + pb_x[k] * ki_53[k];

        t_62[k] = f_16 * ii_54[k]
                  + pb_x[k] * ki_54[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_z, ii_21, ii_23, ii_24, ik_27, \
                         ik_28, ik_30, ik_31, ki_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_y[k] * ik_27[k];

        t_64[k] = f_0 * ii_21[k]
                  + pa_y[k] * ik_28[k];

        t_65[k] = pb_z[k] * ki_49[k];

        t_66[k] = f_15 * ii_23[k]
                  + pa_y[k] * ik_30[k];

        t_67[k] = f_14 * ii_24[k]
                  + pa_y[k] * ik_31[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pa_z, pb_y, ii_25, ii_26, ii_27, \
                         ik_0, ik_32, ik_33, ik_35, ki_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_13 * ii_25[k]
                  + pa_y[k] * ik_32[k];

        t_69[k] = f_12 * ii_26[k]
                  + pa_y[k] * ik_33[k];

        t_70[k] = f_11 * ii_27[k]
                  + pb_y[k] * ki_55[k];

        t_71[k] = pa_y[k] * ik_35[k];

        t_72[k] = pa_z[k] * ik_0[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, ii_0, ii_2, \
                         ik_3, ik_5, ik_6, ki_56, ki_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_y[k] * ki_56[k];

        t_74[k] = f_11 * ii_0[k]
                  + pb_z[k] * ki_56[k];

        t_75[k] = pa_z[k] * ik_3[k];

        t_76[k] = pb_y[k] * ki_58[k];

        t_77[k] = f_12 * ii_2[k]
                  + pa_z[k] * ik_5[k];

        t_78[k] = pa_z[k] * ik_6[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_z, pb_y, pb_z, ii_3, ii_5, ii_6, \
                         ik_9, ik_10, ki_59, ki_61, ki_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_11 * ii_3[k]
                  + pb_z[k] * ki_59[k];

        t_80[k] = pb_y[k] * ki_61[k];

        t_81[k] = f_13 * ii_5[k]
                  + pa_z[k] * ik_9[k];

        t_82[k] = pa_z[k] * ik_10[k];

        t_83[k] = f_11 * ii_6[k]
                  + pb_z[k] * ki_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pa_z, pb_y, pb_z, ii_7, ii_9, ii_10, \
                         ik_12, ik_14, ik_15, ki_65, ki_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * ii_7[k]
                  + pa_z[k] * ik_12[k];

        t_85[k] = pb_y[k] * ki_65[k];

        t_86[k] = f_14 * ii_9[k]
                  + pa_z[k] * ik_14[k];

        t_87[k] = pa_z[k] * ik_15[k];

        t_88[k] = f_11 * ii_10[k]
                  + pb_z[k] * ki_66[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, pa_z, pb_y, ii_11, ii_12, ii_14, ik_17, \
                         ik_18, ik_20, ik_21, ki_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_12 * ii_11[k]
                  + pa_z[k] * ik_17[k];

        t_90[k] = f_13 * ii_12[k]
                  + pa_z[k] * ik_18[k];

        t_91[k] = pb_y[k] * ki_70[k];

        t_92[k] = f_15 * ii_14[k]
                  + pa_z[k] * ik_20[k];

        t_93[k] = pa_z[k] * ik_21[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pb_x, pb_y, ii_78, ii_79, ii_80, ii_81, \
                         ki_76, ki_78, ki_79, ki_80, ki_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_16 * ii_78[k]
                  + pb_x[k] * ki_78[k];

        t_95[k] = f_16 * ii_79[k]
                  + pb_x[k] * ki_79[k];

        t_96[k] = f_16 * ii_80[k]
                  + pb_x[k] * ki_80[k];

        t_97[k] = f_16 * ii_81[k]
                  + pb_x[k] * ki_81[k];

        t_98[k] = pb_y[k] * ki_76[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_z, pb_x, pb_z, ii_21, ii_22, ii_83, \
                         ik_28, ik_30, ki_77, ki_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_16 * ii_83[k]
                  + pb_x[k] * ki_83[k];

        t_100[k] = pa_z[k] * ik_28[k];

        t_101[k] = f_11 * ii_21[k]
                   + pb_z[k] * ki_77[k];

        t_102[k] = f_12 * ii_22[k]
                   + pa_z[k] * ik_30[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, ii_23, ii_24, ii_25, \
                         ii_27, ik_31, ik_32, ik_33, ik_35, ki_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_13 * ii_23[k]
                   + pa_z[k] * ik_31[k];

        t_104[k] = f_14 * ii_24[k]
                   + pa_z[k] * ik_32[k];

        t_105[k] = f_15 * ii_25[k]
                   + pa_z[k] * ik_33[k];

        t_106[k] = pb_y[k] * ki_83[k];

        t_107[k] = f_0 * ii_27[k]
                   + pa_z[k] * ik_35[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_y, pb_y, pb_z, hk0_0, hk1_0, ii_28, ik_36, \
                         ki_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_17 * hk0_0[k]
                   - f_18 * hk1_0[k]
                   + pa_y[k] * ik_36[k];

        t_109[k] = f_12 * ii_28[k]
                   + pb_y[k] * ki_84[k];

        t_110[k] = pb_z[k] * ki_84[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pb_x, pb_z, ii_87, kh0_63, kh0_66, kh1_63, \
                         kh1_66, ki_85, ki_86, ki_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_15 * ii_87[k]
                   + f_9 * kh0_66[k]
                   - f_10 * kh1_66[k]
                   + pb_x[k] * ki_87[k];

        t_112[k] = pb_z[k] * ki_85[k];

        t_113[k] = f_3 * kh0_63[k]
                   - f_4 * kh1_63[k]
                   + pb_z[k] * ki_86[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, pb_y, pb_z, ii_33, ii_90, kh0_65, \
                         kh0_69, kh1_65, kh1_69, ki_87, ki_89, ki_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_15 * ii_90[k]
                   + f_7 * kh0_69[k]
                   - f_8 * kh1_69[k]
                   + pb_x[k] * ki_90[k];

        t_115[k] = pb_z[k] * ki_87[k];

        t_116[k] = f_12 * ii_33[k]
                   + pb_y[k] * ki_89[k];

        t_117[k] = f_5 * kh0_65[k]
                   - f_6 * kh1_65[k]
                   + pb_z[k] * ki_89[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pb_x, pb_z, ii_94, kh0_66, kh0_73, kh1_66, \
                         kh1_73, ki_90, ki_91, ki_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_15 * ii_94[k]
                   + f_5 * kh0_73[k]
                   - f_6 * kh1_73[k]
                   + pb_x[k] * ki_94[k];

        t_119[k] = pb_z[k] * ki_90[k];

        t_120[k] = f_3 * kh0_66[k]
                   - f_4 * kh1_66[k]
                   + pb_z[k] * ki_91[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_x, pb_y, pb_z, ii_37, ii_99, kh0_68, \
                         kh0_78, kh1_68, kh1_78, ki_93, ki_94, ki_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_12 * ii_37[k]
                   + pb_y[k] * ki_93[k];

        t_122[k] = f_7 * kh0_68[k]
                   - f_8 * kh1_68[k]
                   + pb_z[k] * ki_93[k];

        t_123[k] = f_15 * ii_99[k]
                   + f_3 * kh0_78[k]
                   - f_4 * kh1_78[k]
                   + pb_x[k] * ki_99[k];

        t_124[k] = pb_z[k] * ki_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pb_z, ii_42, kh0_69, kh0_70, \
                         kh0_72, kh1_69, kh1_70, kh1_72, ki_95, ki_96, \
                         ki_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * kh0_69[k]
                   - f_4 * kh1_69[k]
                   + pb_z[k] * ki_95[k];

        t_126[k] = f_5 * kh0_70[k]
                   - f_6 * kh1_70[k]
                   + pb_z[k] * ki_96[k];

        t_127[k] = f_12 * ii_42[k]
                   + pb_y[k] * ki_98[k];

        t_128[k] = f_9 * kh0_72[k]
                   - f_10 * kh1_72[k]
                   + pb_z[k] * ki_98[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pb_x, pb_z, ii_105, ii_107, \
                         ii_108, ii_109, ki_99, ki_105, ki_107, ki_108, \
                         ki_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_15 * ii_105[k]
                   + pb_x[k] * ki_105[k];

        t_130[k] = pb_z[k] * ki_99[k];

        t_131[k] = f_15 * ii_107[k]
                   + pb_x[k] * ki_107[k];

        t_132[k] = f_15 * ii_108[k]
                   + pb_x[k] * ki_108[k];

        t_133[k] = f_15 * ii_109[k]
                   + pb_x[k] * ki_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_x, pb_x, pb_z, hk0_136, hk1_136, \
                         ii_110, ii_111, ik_136, ki_105, ki_110, \
                         ki_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_15 * ii_110[k]
                   + pb_x[k] * ki_110[k];

        t_135[k] = f_15 * ii_111[k]
                   + pb_x[k] * ki_111[k];

        t_136[k] = f_19 * hk0_136[k]
                   - f_20 * hk1_136[k]
                   + pa_x[k] * ik_136[k];

        t_137[k] = pb_z[k] * ki_105[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pb_z, kh0_78, kh0_79, kh0_80, kh1_78, kh1_79, \
                         kh1_80, ki_106, ki_107, ki_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_3 * kh0_78[k]
                   - f_4 * kh1_78[k]
                   + pb_z[k] * ki_106[k];

        t_139[k] = f_5 * kh0_79[k]
                   - f_6 * kh1_79[k]
                   + pb_z[k] * ki_107[k];

        t_140[k] = f_7 * kh0_80[k]
                   - f_8 * kh1_80[k]
                   + pb_z[k] * ki_108[k];
    }
}

static auto
compute_prim_kk_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t hk0,
                                            const size_t hk1, const size_t ii, const size_t ik,
                                            const size_t kh0, const size_t kh1, const size_t ki,
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
    const auto f_17 = 0.5 / alpha;
    const auto f_18 = 0.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / alpha;
    const auto f_20 = 2.0 * beta / (alpha * p);
    const auto f_21 = 1.0 / alpha;
    const auto f_22 = beta / (alpha * p);
    const auto f_23 = 1.5 / alpha;
    const auto f_24 = 1.5 * beta / (alpha * p);

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

    const auto *hk0_0 = buffer.data(hk0 + 0);
    const auto *hk0_36 = buffer.data(hk0 + 36);
    const auto *hk0_215 = buffer.data(hk0 + 215);
    const auto *hk0_244 = buffer.data(hk0 + 244);

    const auto *hk1_0 = buffer.data(hk1 + 0);
    const auto *hk1_36 = buffer.data(hk1 + 36);
    const auto *hk1_215 = buffer.data(hk1 + 215);
    const auto *hk1_244 = buffer.data(hk1 + 244);

    const auto *ii_31 = buffer.data(ii + 31);
    const auto *ii_34 = buffer.data(ii + 34);
    const auto *ii_38 = buffer.data(ii + 38);
    const auto *ii_49 = buffer.data(ii + 49);
    const auto *ii_55 = buffer.data(ii + 55);
    const auto *ii_56 = buffer.data(ii + 56);
    const auto *ii_58 = buffer.data(ii + 58);
    const auto *ii_59 = buffer.data(ii + 59);
    const auto *ii_61 = buffer.data(ii + 61);
    const auto *ii_62 = buffer.data(ii + 62);
    const auto *ii_64 = buffer.data(ii + 64);
    const auto *ii_65 = buffer.data(ii + 65);
    const auto *ii_66 = buffer.data(ii + 66);
    const auto *ii_68 = buffer.data(ii + 68);
    const auto *ii_69 = buffer.data(ii + 69);
    const auto *ii_70 = buffer.data(ii + 70);
    const auto *ii_77 = buffer.data(ii + 77);
    const auto *ii_79 = buffer.data(ii + 79);
    const auto *ii_80 = buffer.data(ii + 80);
    const auto *ii_81 = buffer.data(ii + 81);
    const auto *ii_82 = buffer.data(ii + 82);
    const auto *ii_83 = buffer.data(ii + 83);
    const auto *ii_84 = buffer.data(ii + 84);
    const auto *ii_86 = buffer.data(ii + 86);
    const auto *ii_87 = buffer.data(ii + 87);
    const auto *ii_89 = buffer.data(ii + 89);
    const auto *ii_90 = buffer.data(ii + 90);
    const auto *ii_91 = buffer.data(ii + 91);
    const auto *ii_93 = buffer.data(ii + 93);
    const auto *ii_94 = buffer.data(ii + 94);
    const auto *ii_95 = buffer.data(ii + 95);
    const auto *ii_96 = buffer.data(ii + 96);
    const auto *ii_98 = buffer.data(ii + 98);
    const auto *ii_111 = buffer.data(ii + 111);
    const auto *ii_114 = buffer.data(ii + 114);
    const auto *ii_117 = buffer.data(ii + 117);
    const auto *ii_121 = buffer.data(ii + 121);
    const auto *ii_126 = buffer.data(ii + 126);
    const auto *ii_134 = buffer.data(ii + 134);
    const auto *ii_135 = buffer.data(ii + 135);
    const auto *ii_136 = buffer.data(ii + 136);
    const auto *ii_137 = buffer.data(ii + 137);
    const auto *ii_138 = buffer.data(ii + 138);
    const auto *ii_145 = buffer.data(ii + 145);
    const auto *ii_149 = buffer.data(ii + 149);
    const auto *ii_154 = buffer.data(ii + 154);
    const auto *ii_160 = buffer.data(ii + 160);
    const auto *ii_161 = buffer.data(ii + 161);
    const auto *ii_162 = buffer.data(ii + 162);
    const auto *ii_163 = buffer.data(ii + 163);
    const auto *ii_164 = buffer.data(ii + 164);
    const auto *ii_165 = buffer.data(ii + 165);
    const auto *ii_167 = buffer.data(ii + 167);
    const auto *ii_171 = buffer.data(ii + 171);
    const auto *ii_174 = buffer.data(ii + 174);
    const auto *ii_178 = buffer.data(ii + 178);
    const auto *ii_183 = buffer.data(ii + 183);
    const auto *ii_189 = buffer.data(ii + 189);
    const auto *ii_191 = buffer.data(ii + 191);
    const auto *ii_192 = buffer.data(ii + 192);
    const auto *ii_193 = buffer.data(ii + 193);
    const auto *ii_194 = buffer.data(ii + 194);
    const auto *ii_195 = buffer.data(ii + 195);
    const auto *ii_218 = buffer.data(ii + 218);

    const auto *ik_37 = buffer.data(ik + 37);
    const auto *ik_39 = buffer.data(ik + 39);
    const auto *ik_42 = buffer.data(ik + 42);
    const auto *ik_46 = buffer.data(ik + 46);
    const auto *ik_51 = buffer.data(ik + 51);
    const auto *ik_57 = buffer.data(ik + 57);
    const auto *ik_64 = buffer.data(ik + 64);
    const auto *ik_72 = buffer.data(ik + 72);
    const auto *ik_74 = buffer.data(ik + 74);
    const auto *ik_77 = buffer.data(ik + 77);
    const auto *ik_81 = buffer.data(ik + 81);
    const auto *ik_84 = buffer.data(ik + 84);
    const auto *ik_86 = buffer.data(ik + 86);
    const auto *ik_89 = buffer.data(ik + 89);
    const auto *ik_90 = buffer.data(ik + 90);
    const auto *ik_92 = buffer.data(ik + 92);
    const auto *ik_99 = buffer.data(ik + 99);
    const auto *ik_102 = buffer.data(ik + 102);
    const auto *ik_103 = buffer.data(ik + 103);
    const auto *ik_104 = buffer.data(ik + 104);
    const auto *ik_105 = buffer.data(ik + 105);
    const auto *ik_107 = buffer.data(ik + 107);
    const auto *ik_108 = buffer.data(ik + 108);
    const auto *ik_109 = buffer.data(ik + 109);
    const auto *ik_111 = buffer.data(ik + 111);
    const auto *ik_113 = buffer.data(ik + 113);
    const auto *ik_114 = buffer.data(ik + 114);
    const auto *ik_117 = buffer.data(ik + 117);
    const auto *ik_118 = buffer.data(ik + 118);
    const auto *ik_120 = buffer.data(ik + 120);
    const auto *ik_122 = buffer.data(ik + 122);
    const auto *ik_123 = buffer.data(ik + 123);
    const auto *ik_125 = buffer.data(ik + 125);
    const auto *ik_126 = buffer.data(ik + 126);
    const auto *ik_128 = buffer.data(ik + 128);
    const auto *ik_129 = buffer.data(ik + 129);
    const auto *ik_215 = buffer.data(ik + 215);
    const auto *ik_244 = buffer.data(ik + 244);

    const auto *kh0_81 = buffer.data(kh0 + 81);
    const auto *kh0_83 = buffer.data(kh0 + 83);
    const auto *kh0_105 = buffer.data(kh0 + 105);
    const auto *kh0_106 = buffer.data(kh0 + 106);
    const auto *kh0_108 = buffer.data(kh0 + 108);
    const auto *kh0_110 = buffer.data(kh0 + 110);
    const auto *kh0_111 = buffer.data(kh0 + 111);
    const auto *kh0_113 = buffer.data(kh0 + 113);
    const auto *kh0_114 = buffer.data(kh0 + 114);
    const auto *kh0_119 = buffer.data(kh0 + 119);
    const auto *kh0_120 = buffer.data(kh0 + 120);
    const auto *kh0_122 = buffer.data(kh0 + 122);
    const auto *kh0_123 = buffer.data(kh0 + 123);
    const auto *kh0_124 = buffer.data(kh0 + 124);
    const auto *kh0_125 = buffer.data(kh0 + 125);
    const auto *kh0_126 = buffer.data(kh0 + 126);
    const auto *kh0_128 = buffer.data(kh0 + 128);
    const auto *kh0_129 = buffer.data(kh0 + 129);
    const auto *kh0_131 = buffer.data(kh0 + 131);
    const auto *kh0_132 = buffer.data(kh0 + 132);
    const auto *kh0_133 = buffer.data(kh0 + 133);
    const auto *kh0_135 = buffer.data(kh0 + 135);
    const auto *kh0_136 = buffer.data(kh0 + 136);
    const auto *kh0_141 = buffer.data(kh0 + 141);
    const auto *kh0_142 = buffer.data(kh0 + 142);
    const auto *kh0_143 = buffer.data(kh0 + 143);
    const auto *kh0_144 = buffer.data(kh0 + 144);
    const auto *kh0_146 = buffer.data(kh0 + 146);

    const auto *kh1_81 = buffer.data(kh1 + 81);
    const auto *kh1_83 = buffer.data(kh1 + 83);
    const auto *kh1_105 = buffer.data(kh1 + 105);
    const auto *kh1_106 = buffer.data(kh1 + 106);
    const auto *kh1_108 = buffer.data(kh1 + 108);
    const auto *kh1_110 = buffer.data(kh1 + 110);
    const auto *kh1_111 = buffer.data(kh1 + 111);
    const auto *kh1_113 = buffer.data(kh1 + 113);
    const auto *kh1_114 = buffer.data(kh1 + 114);
    const auto *kh1_119 = buffer.data(kh1 + 119);
    const auto *kh1_120 = buffer.data(kh1 + 120);
    const auto *kh1_122 = buffer.data(kh1 + 122);
    const auto *kh1_123 = buffer.data(kh1 + 123);
    const auto *kh1_124 = buffer.data(kh1 + 124);
    const auto *kh1_125 = buffer.data(kh1 + 125);
    const auto *kh1_126 = buffer.data(kh1 + 126);
    const auto *kh1_128 = buffer.data(kh1 + 128);
    const auto *kh1_129 = buffer.data(kh1 + 129);
    const auto *kh1_131 = buffer.data(kh1 + 131);
    const auto *kh1_132 = buffer.data(kh1 + 132);
    const auto *kh1_133 = buffer.data(kh1 + 133);
    const auto *kh1_135 = buffer.data(kh1 + 135);
    const auto *kh1_136 = buffer.data(kh1 + 136);
    const auto *kh1_141 = buffer.data(kh1 + 141);
    const auto *kh1_142 = buffer.data(kh1 + 142);
    const auto *kh1_143 = buffer.data(kh1 + 143);
    const auto *kh1_144 = buffer.data(kh1 + 144);
    const auto *kh1_146 = buffer.data(kh1 + 146);

    const auto *ki_109 = buffer.data(ki + 109);
    const auto *ki_111 = buffer.data(ki + 111);
    const auto *ki_114 = buffer.data(ki + 114);
    const auto *ki_115 = buffer.data(ki + 115);
    const auto *ki_117 = buffer.data(ki + 117);
    const auto *ki_118 = buffer.data(ki + 118);
    const auto *ki_121 = buffer.data(ki + 121);
    const auto *ki_122 = buffer.data(ki + 122);
    const auto *ki_126 = buffer.data(ki + 126);
    const auto *ki_133 = buffer.data(ki + 133);
    const auto *ki_134 = buffer.data(ki + 134);
    const auto *ki_135 = buffer.data(ki + 135);
    const auto *ki_136 = buffer.data(ki + 136);
    const auto *ki_137 = buffer.data(ki + 137);
    const auto *ki_138 = buffer.data(ki + 138);
    const auto *ki_139 = buffer.data(ki + 139);
    const auto *ki_140 = buffer.data(ki + 140);
    const auto *ki_141 = buffer.data(ki + 141);
    const auto *ki_142 = buffer.data(ki + 142);
    const auto *ki_143 = buffer.data(ki + 143);
    const auto *ki_145 = buffer.data(ki + 145);
    const auto *ki_146 = buffer.data(ki + 146);
    const auto *ki_148 = buffer.data(ki + 148);
    const auto *ki_149 = buffer.data(ki + 149);
    const auto *ki_150 = buffer.data(ki + 150);
    const auto *ki_152 = buffer.data(ki + 152);
    const auto *ki_153 = buffer.data(ki + 153);
    const auto *ki_154 = buffer.data(ki + 154);
    const auto *ki_160 = buffer.data(ki + 160);
    const auto *ki_161 = buffer.data(ki + 161);
    const auto *ki_162 = buffer.data(ki + 162);
    const auto *ki_163 = buffer.data(ki + 163);
    const auto *ki_164 = buffer.data(ki + 164);
    const auto *ki_165 = buffer.data(ki + 165);
    const auto *ki_166 = buffer.data(ki + 166);
    const auto *ki_167 = buffer.data(ki + 167);
    const auto *ki_168 = buffer.data(ki + 168);
    const auto *ki_169 = buffer.data(ki + 169);
    const auto *ki_170 = buffer.data(ki + 170);
    const auto *ki_171 = buffer.data(ki + 171);
    const auto *ki_173 = buffer.data(ki + 173);
    const auto *ki_174 = buffer.data(ki + 174);
    const auto *ki_175 = buffer.data(ki + 175);
    const auto *ki_177 = buffer.data(ki + 177);
    const auto *ki_178 = buffer.data(ki + 178);
    const auto *ki_179 = buffer.data(ki + 179);
    const auto *ki_180 = buffer.data(ki + 180);
    const auto *ki_182 = buffer.data(ki + 182);
    const auto *ki_183 = buffer.data(ki + 183);
    const auto *ki_189 = buffer.data(ki + 189);
    const auto *ki_190 = buffer.data(ki + 190);
    const auto *ki_191 = buffer.data(ki + 191);
    const auto *ki_192 = buffer.data(ki + 192);
    const auto *ki_193 = buffer.data(ki + 193);
    const auto *ki_194 = buffer.data(ki + 194);
    const auto *ki_195 = buffer.data(ki + 195);
    const auto *ki_196 = buffer.data(ki + 196);
    const auto *ki_198 = buffer.data(ki + 198);
    const auto *ki_199 = buffer.data(ki + 199);
    const auto *ki_201 = buffer.data(ki + 201);
    const auto *ki_202 = buffer.data(ki + 202);
    const auto *ki_205 = buffer.data(ki + 205);
    const auto *ki_206 = buffer.data(ki + 206);
    const auto *ki_210 = buffer.data(ki + 210);
    const auto *ki_218 = buffer.data(ki + 218);

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pb_y, pb_z, ii_55, ik_72, kh0_81, \
                         kh0_83, kh1_81, kh1_83, ki_109, ki_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_9 * kh0_81[k]
                   - f_10 * kh1_81[k]
                   + pb_z[k] * ki_109[k];

        t_142[k] = f_12 * ii_55[k]
                   + pb_y[k] * ki_111[k];

        t_143[k] = f_1 * kh0_83[k]
                   - f_2 * kh1_83[k]
                   + pb_z[k] * ki_111[k];

        t_144[k] = pa_y[k] * ik_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, t_150, pa_y, pa_z, pb_y, ii_58, \
                         ik_37, ik_39, ik_42, ik_74, ik_77, ki_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = pa_z[k] * ik_37[k];

        t_146[k] = pa_y[k] * ik_74[k];

        t_147[k] = pa_z[k] * ik_39[k];

        t_148[k] = f_11 * ii_58[k]
                   + pb_y[k] * ki_114[k];

        t_149[k] = pa_y[k] * ik_77[k];

        t_150[k] = pa_z[k] * ik_42[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pa_y, pa_z, pb_y, pb_z, ii_31, ii_61, \
                         ik_46, ik_81, ki_115, ki_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_11 * ii_31[k]
                   + pb_z[k] * ki_115[k];

        t_152[k] = f_11 * ii_61[k]
                   + pb_y[k] * ki_117[k];

        t_153[k] = pa_y[k] * ik_81[k];

        t_154[k] = pa_z[k] * ik_46[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_y, pb_y, pb_z, ii_34, ii_64, ii_65, \
                         ik_84, ik_86, ki_118, ki_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_11 * ii_34[k]
                   + pb_z[k] * ki_118[k];

        t_156[k] = f_12 * ii_64[k]
                   + pa_y[k] * ik_84[k];

        t_157[k] = f_11 * ii_65[k]
                   + pb_y[k] * ki_121[k];

        t_158[k] = pa_y[k] * ik_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_y, pa_z, pb_z, ii_38, ii_68, ii_69, \
                         ik_51, ik_89, ik_90, ki_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * ik_51[k];

        t_160[k] = f_11 * ii_38[k]
                   + pb_z[k] * ki_122[k];

        t_161[k] = f_13 * ii_68[k]
                   + pa_y[k] * ik_89[k];

        t_162[k] = f_12 * ii_69[k]
                   + pa_y[k] * ik_90[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_y, pa_z, pb_x, pb_y, ii_70, ii_134, \
                         ik_57, ik_92, ki_126, ki_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_11 * ii_70[k]
                   + pb_y[k] * ki_126[k];

        t_164[k] = pa_y[k] * ik_92[k];

        t_165[k] = pa_z[k] * ik_57[k];

        t_166[k] = f_15 * ii_134[k]
                   + pb_x[k] * ki_134[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pa_y, pb_x, ii_135, ii_136, \
                         ii_137, ii_138, ik_99, ki_135, ki_136, ki_137, \
                         ki_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_15 * ii_135[k]
                   + pb_x[k] * ki_135[k];

        t_168[k] = f_15 * ii_136[k]
                   + pb_x[k] * ki_136[k];

        t_169[k] = f_15 * ii_137[k]
                   + pb_x[k] * ki_137[k];

        t_170[k] = f_15 * ii_138[k]
                   + pb_x[k] * ki_138[k];

        t_171[k] = pa_y[k] * ik_99[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pa_z, pb_z, ii_49, ii_79, ii_80, \
                         ik_64, ik_102, ik_103, ki_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pa_z[k] * ik_64[k];

        t_173[k] = f_11 * ii_49[k]
                   + pb_z[k] * ki_133[k];

        t_174[k] = f_15 * ii_79[k]
                   + pa_y[k] * ik_102[k];

        t_175[k] = f_14 * ii_80[k]
                   + pa_y[k] * ik_103[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_y, ii_81, ii_82, ii_83, ik_104, \
                         ik_105, ik_107, ki_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * ii_81[k]
                   + pa_y[k] * ik_104[k];

        t_177[k] = f_12 * ii_82[k]
                   + pa_y[k] * ik_105[k];

        t_178[k] = f_11 * ii_83[k]
                   + pb_y[k] * ki_139[k];

        t_179[k] = pa_y[k] * ik_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_z, pb_y, pb_z, hk0_0, hk1_0, ii_56, \
                         ik_72, kh0_105, kh1_105, ki_140, ki_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_17 * hk0_0[k]
                   - f_18 * hk1_0[k]
                   + pa_z[k] * ik_72[k];

        t_181[k] = pb_y[k] * ki_140[k];

        t_182[k] = f_12 * ii_56[k]
                   + pb_z[k] * ki_140[k];

        t_183[k] = f_3 * kh0_105[k]
                   - f_4 * kh1_105[k]
                   + pb_y[k] * ki_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pb_x, pb_y, pb_z, ii_59, ii_145, kh0_106, \
                         kh0_110, kh1_106, kh1_110, ki_142, ki_143, \
                         ki_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = pb_y[k] * ki_142[k];

        t_185[k] = f_15 * ii_145[k]
                   + f_9 * kh0_110[k]
                   - f_10 * kh1_110[k]
                   + pb_x[k] * ki_145[k];

        t_186[k] = f_5 * kh0_106[k]
                   - f_6 * kh1_106[k]
                   + pb_y[k] * ki_143[k];

        t_187[k] = f_12 * ii_59[k]
                   + pb_z[k] * ki_143[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pb_x, pb_y, pb_z, ii_62, ii_149, kh0_108, \
                         kh0_114, kh1_108, kh1_114, ki_145, ki_146, \
                         ki_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pb_y[k] * ki_145[k];

        t_189[k] = f_15 * ii_149[k]
                   + f_7 * kh0_114[k]
                   - f_8 * kh1_114[k]
                   + pb_x[k] * ki_149[k];

        t_190[k] = f_7 * kh0_108[k]
                   - f_8 * kh1_108[k]
                   + pb_y[k] * ki_146[k];

        t_191[k] = f_12 * ii_62[k]
                   + pb_z[k] * ki_146[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_x, pb_y, ii_154, kh0_110, kh0_119, kh1_110, \
                         kh1_119, ki_148, ki_149, ki_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_3 * kh0_110[k]
                   - f_4 * kh1_110[k]
                   + pb_y[k] * ki_148[k];

        t_193[k] = pb_y[k] * ki_149[k];

        t_194[k] = f_15 * ii_154[k]
                   + f_5 * kh0_119[k]
                   - f_6 * kh1_119[k]
                   + pb_x[k] * ki_154[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pb_y, pb_z, ii_66, kh0_111, kh0_113, \
                         kh0_114, kh1_111, kh1_113, kh1_114, ki_150, ki_152, \
                         ki_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_9 * kh0_111[k]
                   - f_10 * kh1_111[k]
                   + pb_y[k] * ki_150[k];

        t_196[k] = f_12 * ii_66[k]
                   + pb_z[k] * ki_150[k];

        t_197[k] = f_5 * kh0_113[k]
                   - f_6 * kh1_113[k]
                   + pb_y[k] * ki_152[k];

        t_198[k] = f_3 * kh0_114[k]
                   - f_4 * kh1_114[k]
                   + pb_y[k] * ki_153[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pb_x, pb_y, ii_160, ii_161, ii_162, \
                         kh0_125, kh1_125, ki_154, ki_160, ki_161, \
                         ki_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = pb_y[k] * ki_154[k];

        t_200[k] = f_15 * ii_160[k]
                   + f_3 * kh0_125[k]
                   - f_4 * kh1_125[k]
                   + pb_x[k] * ki_160[k];

        t_201[k] = f_15 * ii_161[k]
                   + pb_x[k] * ki_161[k];

        t_202[k] = f_15 * ii_162[k]
                   + pb_x[k] * ki_162[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pb_x, pb_y, ii_163, ii_164, \
                         ii_165, ii_167, ki_160, ki_163, ki_164, ki_165, \
                         ki_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_15 * ii_163[k]
                   + pb_x[k] * ki_163[k];

        t_204[k] = f_15 * ii_164[k]
                   + pb_x[k] * ki_164[k];

        t_205[k] = f_15 * ii_165[k]
                   + pb_x[k] * ki_165[k];

        t_206[k] = pb_y[k] * ki_160[k];

        t_207[k] = f_15 * ii_167[k]
                   + pb_x[k] * ki_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_y, pb_z, ii_77, kh0_120, kh0_122, \
                         kh0_123, kh1_120, kh1_122, kh1_123, ki_161, ki_163, \
                         ki_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * kh0_120[k]
                   - f_2 * kh1_120[k]
                   + pb_y[k] * ki_161[k];

        t_209[k] = f_12 * ii_77[k]
                   + pb_z[k] * ki_161[k];

        t_210[k] = f_9 * kh0_122[k]
                   - f_10 * kh1_122[k]
                   + pb_y[k] * ki_163[k];

        t_211[k] = f_7 * kh0_123[k]
                   - f_8 * kh1_123[k]
                   + pb_y[k] * ki_164[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pa_x, pb_y, hk0_215, hk1_215, ik_215, \
                         kh0_124, kh0_125, kh1_124, kh1_125, ki_165, ki_166, \
                         ki_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_5 * kh0_124[k]
                   - f_6 * kh1_124[k]
                   + pb_y[k] * ki_165[k];

        t_213[k] = f_3 * kh0_125[k]
                   - f_4 * kh1_125[k]
                   + pb_y[k] * ki_166[k];

        t_214[k] = pb_y[k] * ki_167[k];

        t_215[k] = f_19 * hk0_215[k]
                   - f_20 * hk1_215[k]
                   + pa_x[k] * ik_215[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, pa_y, pb_y, pb_z, hk0_36, hk1_36, ii_84, ik_108, \
                         ki_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_21 * hk0_36[k]
                   - f_22 * hk1_36[k]
                   + pa_y[k] * ik_108[k];

        t_217[k] = f_13 * ii_84[k]
                   + pb_y[k] * ki_168[k];

        t_218[k] = pb_z[k] * ki_168[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pb_x, pb_z, ii_171, kh0_126, kh0_129, kh1_126, \
                         kh1_129, ki_169, ki_170, ki_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_14 * ii_171[k]
                   + f_9 * kh0_129[k]
                   - f_10 * kh1_129[k]
                   + pb_x[k] * ki_171[k];

        t_220[k] = pb_z[k] * ki_169[k];

        t_221[k] = f_3 * kh0_126[k]
                   - f_4 * kh1_126[k]
                   + pb_z[k] * ki_170[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pb_x, pb_y, pb_z, ii_89, ii_174, kh0_128, \
                         kh0_132, kh1_128, kh1_132, ki_171, ki_173, \
                         ki_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_14 * ii_174[k]
                   + f_7 * kh0_132[k]
                   - f_8 * kh1_132[k]
                   + pb_x[k] * ki_174[k];

        t_223[k] = pb_z[k] * ki_171[k];

        t_224[k] = f_13 * ii_89[k]
                   + pb_y[k] * ki_173[k];

        t_225[k] = f_5 * kh0_128[k]
                   - f_6 * kh1_128[k]
                   + pb_z[k] * ki_173[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pb_x, pb_z, ii_178, kh0_129, kh0_136, kh1_129, \
                         kh1_136, ki_174, ki_175, ki_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_14 * ii_178[k]
                   + f_5 * kh0_136[k]
                   - f_6 * kh1_136[k]
                   + pb_x[k] * ki_178[k];

        t_227[k] = pb_z[k] * ki_174[k];

        t_228[k] = f_3 * kh0_129[k]
                   - f_4 * kh1_129[k]
                   + pb_z[k] * ki_175[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_x, pb_y, pb_z, ii_93, ii_183, kh0_131, \
                         kh0_141, kh1_131, kh1_141, ki_177, ki_178, \
                         ki_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_13 * ii_93[k]
                   + pb_y[k] * ki_177[k];

        t_230[k] = f_7 * kh0_131[k]
                   - f_8 * kh1_131[k]
                   + pb_z[k] * ki_177[k];

        t_231[k] = f_14 * ii_183[k]
                   + f_3 * kh0_141[k]
                   - f_4 * kh1_141[k]
                   + pb_x[k] * ki_183[k];

        t_232[k] = pb_z[k] * ki_178[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pb_y, pb_z, ii_98, kh0_132, kh0_133, \
                         kh0_135, kh1_132, kh1_133, kh1_135, ki_179, ki_180, \
                         ki_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * kh0_132[k]
                   - f_4 * kh1_132[k]
                   + pb_z[k] * ki_179[k];

        t_234[k] = f_5 * kh0_133[k]
                   - f_6 * kh1_133[k]
                   + pb_z[k] * ki_180[k];

        t_235[k] = f_13 * ii_98[k]
                   + pb_y[k] * ki_182[k];

        t_236[k] = f_9 * kh0_135[k]
                   - f_10 * kh1_135[k]
                   + pb_z[k] * ki_182[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, pb_x, pb_z, ii_189, ii_191, \
                         ii_192, ii_193, ki_183, ki_189, ki_191, ki_192, \
                         ki_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_14 * ii_189[k]
                   + pb_x[k] * ki_189[k];

        t_238[k] = pb_z[k] * ki_183[k];

        t_239[k] = f_14 * ii_191[k]
                   + pb_x[k] * ki_191[k];

        t_240[k] = f_14 * ii_192[k]
                   + pb_x[k] * ki_192[k];

        t_241[k] = f_14 * ii_193[k]
                   + pb_x[k] * ki_193[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pa_x, pb_x, pb_z, hk0_244, hk1_244, \
                         ii_194, ii_195, ik_244, ki_189, ki_194, \
                         ki_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_14 * ii_194[k]
                   + pb_x[k] * ki_194[k];

        t_243[k] = f_14 * ii_195[k]
                   + pb_x[k] * ki_195[k];

        t_244[k] = f_23 * hk0_244[k]
                   - f_24 * hk1_244[k]
                   + pa_x[k] * ik_244[k];

        t_245[k] = pb_z[k] * ki_189[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pb_z, kh0_141, kh0_142, kh0_143, kh1_141, \
                         kh1_142, kh1_143, ki_190, ki_191, ki_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_3 * kh0_141[k]
                   - f_4 * kh1_141[k]
                   + pb_z[k] * ki_190[k];

        t_247[k] = f_5 * kh0_142[k]
                   - f_6 * kh1_142[k]
                   + pb_z[k] * ki_191[k];

        t_248[k] = f_7 * kh0_143[k]
                   - f_8 * kh1_143[k]
                   + pb_z[k] * ki_192[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, pa_z, pb_y, pb_z, ii_111, ik_108, \
                         kh0_144, kh0_146, kh1_144, kh1_146, ki_193, \
                         ki_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_9 * kh0_144[k]
                   - f_10 * kh1_144[k]
                   + pb_z[k] * ki_193[k];

        t_250[k] = f_13 * ii_111[k]
                   + pb_y[k] * ki_195[k];

        t_251[k] = f_1 * kh0_146[k]
                   - f_2 * kh1_146[k]
                   + pb_z[k] * ki_195[k];

        t_252[k] = pa_z[k] * ik_108[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, pa_z, pb_y, pb_z, ii_84, ii_86, \
                         ii_114, ik_109, ik_111, ik_113, ki_196, \
                         ki_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = pa_z[k] * ik_109[k];

        t_254[k] = f_11 * ii_84[k]
                   + pb_z[k] * ki_196[k];

        t_255[k] = pa_z[k] * ik_111[k];

        t_256[k] = f_12 * ii_114[k]
                   + pb_y[k] * ki_198[k];

        t_257[k] = f_12 * ii_86[k]
                   + pa_z[k] * ik_113[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, pa_z, pb_y, pb_z, ii_87, ii_89, \
                         ii_117, ik_114, ik_117, ik_118, ki_199, \
                         ki_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = pa_z[k] * ik_114[k];

        t_259[k] = f_11 * ii_87[k]
                   + pb_z[k] * ki_199[k];

        t_260[k] = f_12 * ii_117[k]
                   + pb_y[k] * ki_201[k];

        t_261[k] = f_13 * ii_89[k]
                   + pa_z[k] * ik_117[k];

        t_262[k] = pa_z[k] * ik_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, pa_z, pb_y, pb_z, ii_90, ii_91, ii_93, \
                         ii_121, ik_120, ik_122, ki_202, ki_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_11 * ii_90[k]
                   + pb_z[k] * ki_202[k];

        t_264[k] = f_12 * ii_91[k]
                   + pa_z[k] * ik_120[k];

        t_265[k] = f_12 * ii_121[k]
                   + pb_y[k] * ki_205[k];

        t_266[k] = f_14 * ii_93[k]
                   + pa_z[k] * ik_122[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pa_z, pb_z, ii_94, ii_95, ii_96, ik_123, \
                         ik_125, ik_126, ki_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = pa_z[k] * ik_123[k];

        t_268[k] = f_11 * ii_94[k]
                   + pb_z[k] * ki_206[k];

        t_269[k] = f_12 * ii_95[k]
                   + pa_z[k] * ik_125[k];

        t_270[k] = f_13 * ii_96[k]
                   + pa_z[k] * ik_126[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pa_z, pb_x, pb_y, ii_98, ii_126, ii_218, \
                         ik_128, ik_129, ki_210, ki_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_12 * ii_126[k]
                   + pb_y[k] * ki_210[k];

        t_272[k] = f_15 * ii_98[k]
                   + pa_z[k] * ik_128[k];

        t_273[k] = pa_z[k] * ik_129[k];

        t_274[k] = f_14 * ii_218[k]
                   + pb_x[k] * ki_218[k];
    }
}

static auto
compute_prim_kk_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t hk0,
                                            const size_t hk1, const size_t ii, const size_t ik,
                                            const size_t kh0, const size_t kh1, const size_t ki,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
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
    const auto f_21 = 1.0 / alpha;
    const auto f_22 = beta / (alpha * p);
    const auto f_23 = 1.5 / alpha;
    const auto f_24 = 1.5 * beta / (alpha * p);

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

    const auto *hk0_72 = buffer.data(hk0 + 72);
    const auto *hk0_108 = buffer.data(hk0 + 108);
    const auto *hk0_359 = buffer.data(hk0 + 359);
    const auto *hk0_388 = buffer.data(hk0 + 388);

    const auto *hk1_72 = buffer.data(hk1 + 72);
    const auto *hk1_108 = buffer.data(hk1 + 108);
    const auto *hk1_359 = buffer.data(hk1 + 359);
    const auto *hk1_388 = buffer.data(hk1 + 388);

    const auto *ii_105 = buffer.data(ii + 105);
    const auto *ii_106 = buffer.data(ii + 106);
    const auto *ii_107 = buffer.data(ii + 107);
    const auto *ii_108 = buffer.data(ii + 108);
    const auto *ii_109 = buffer.data(ii + 109);
    const auto *ii_111 = buffer.data(ii + 111);
    const auto *ii_115 = buffer.data(ii + 115);
    const auto *ii_118 = buffer.data(ii + 118);
    const auto *ii_122 = buffer.data(ii + 122);
    const auto *ii_133 = buffer.data(ii + 133);
    const auto *ii_139 = buffer.data(ii + 139);
    const auto *ii_140 = buffer.data(ii + 140);
    const auto *ii_141 = buffer.data(ii + 141);
    const auto *ii_142 = buffer.data(ii + 142);
    const auto *ii_143 = buffer.data(ii + 143);
    const auto *ii_145 = buffer.data(ii + 145);
    const auto *ii_146 = buffer.data(ii + 146);
    const auto *ii_148 = buffer.data(ii + 148);
    const auto *ii_149 = buffer.data(ii + 149);
    const auto *ii_150 = buffer.data(ii + 150);
    const auto *ii_152 = buffer.data(ii + 152);
    const auto *ii_153 = buffer.data(ii + 153);
    const auto *ii_154 = buffer.data(ii + 154);
    const auto *ii_161 = buffer.data(ii + 161);
    const auto *ii_163 = buffer.data(ii + 163);
    const auto *ii_164 = buffer.data(ii + 164);
    const auto *ii_165 = buffer.data(ii + 165);
    const auto *ii_166 = buffer.data(ii + 166);
    const auto *ii_167 = buffer.data(ii + 167);
    const auto *ii_168 = buffer.data(ii + 168);
    const auto *ii_170 = buffer.data(ii + 170);
    const auto *ii_171 = buffer.data(ii + 171);
    const auto *ii_173 = buffer.data(ii + 173);
    const auto *ii_177 = buffer.data(ii + 177);
    const auto *ii_182 = buffer.data(ii + 182);
    const auto *ii_195 = buffer.data(ii + 195);
    const auto *ii_198 = buffer.data(ii + 198);
    const auto *ii_201 = buffer.data(ii + 201);
    const auto *ii_219 = buffer.data(ii + 219);
    const auto *ii_220 = buffer.data(ii + 220);
    const auto *ii_221 = buffer.data(ii + 221);
    const auto *ii_222 = buffer.data(ii + 222);
    const auto *ii_223 = buffer.data(ii + 223);
    const auto *ii_245 = buffer.data(ii + 245);
    const auto *ii_246 = buffer.data(ii + 246);
    const auto *ii_247 = buffer.data(ii + 247);
    const auto *ii_248 = buffer.data(ii + 248);
    const auto *ii_249 = buffer.data(ii + 249);
    const auto *ii_250 = buffer.data(ii + 250);
    const auto *ii_257 = buffer.data(ii + 257);
    const auto *ii_261 = buffer.data(ii + 261);
    const auto *ii_266 = buffer.data(ii + 266);
    const auto *ii_272 = buffer.data(ii + 272);
    const auto *ii_273 = buffer.data(ii + 273);
    const auto *ii_274 = buffer.data(ii + 274);
    const auto *ii_275 = buffer.data(ii + 275);
    const auto *ii_276 = buffer.data(ii + 276);
    const auto *ii_277 = buffer.data(ii + 277);
    const auto *ii_279 = buffer.data(ii + 279);
    const auto *ii_283 = buffer.data(ii + 283);
    const auto *ii_286 = buffer.data(ii + 286);
    const auto *ii_290 = buffer.data(ii + 290);
    const auto *ii_295 = buffer.data(ii + 295);
    const auto *ii_301 = buffer.data(ii + 301);
    const auto *ii_303 = buffer.data(ii + 303);
    const auto *ii_304 = buffer.data(ii + 304);
    const auto *ii_305 = buffer.data(ii + 305);
    const auto *ii_306 = buffer.data(ii + 306);
    const auto *ii_307 = buffer.data(ii + 307);

    const auto *ik_136 = buffer.data(ik + 136);
    const auto *ik_138 = buffer.data(ik + 138);
    const auto *ik_139 = buffer.data(ik + 139);
    const auto *ik_140 = buffer.data(ik + 140);
    const auto *ik_141 = buffer.data(ik + 141);
    const auto *ik_143 = buffer.data(ik + 143);
    const auto *ik_180 = buffer.data(ik + 180);
    const auto *ik_182 = buffer.data(ik + 182);
    const auto *ik_183 = buffer.data(ik + 183);
    const auto *ik_185 = buffer.data(ik + 185);
    const auto *ik_186 = buffer.data(ik + 186);
    const auto *ik_189 = buffer.data(ik + 189);
    const auto *ik_190 = buffer.data(ik + 190);
    const auto *ik_192 = buffer.data(ik + 192);
    const auto *ik_194 = buffer.data(ik + 194);
    const auto *ik_195 = buffer.data(ik + 195);
    const auto *ik_197 = buffer.data(ik + 197);
    const auto *ik_198 = buffer.data(ik + 198);
    const auto *ik_200 = buffer.data(ik + 200);
    const auto *ik_207 = buffer.data(ik + 207);
    const auto *ik_208 = buffer.data(ik + 208);
    const auto *ik_210 = buffer.data(ik + 210);
    const auto *ik_211 = buffer.data(ik + 211);
    const auto *ik_212 = buffer.data(ik + 212);
    const auto *ik_213 = buffer.data(ik + 213);
    const auto *ik_215 = buffer.data(ik + 215);
    const auto *ik_216 = buffer.data(ik + 216);
    const auto *ik_217 = buffer.data(ik + 217);
    const auto *ik_219 = buffer.data(ik + 219);
    const auto *ik_221 = buffer.data(ik + 221);
    const auto *ik_222 = buffer.data(ik + 222);
    const auto *ik_225 = buffer.data(ik + 225);
    const auto *ik_226 = buffer.data(ik + 226);
    const auto *ik_359 = buffer.data(ik + 359);
    const auto *ik_388 = buffer.data(ik + 388);

    const auto *kh0_189 = buffer.data(kh0 + 189);
    const auto *kh0_190 = buffer.data(kh0 + 190);
    const auto *kh0_192 = buffer.data(kh0 + 192);
    const auto *kh0_194 = buffer.data(kh0 + 194);
    const auto *kh0_195 = buffer.data(kh0 + 195);
    const auto *kh0_197 = buffer.data(kh0 + 197);
    const auto *kh0_198 = buffer.data(kh0 + 198);
    const auto *kh0_203 = buffer.data(kh0 + 203);
    const auto *kh0_204 = buffer.data(kh0 + 204);
    const auto *kh0_206 = buffer.data(kh0 + 206);
    const auto *kh0_207 = buffer.data(kh0 + 207);
    const auto *kh0_208 = buffer.data(kh0 + 208);
    const auto *kh0_209 = buffer.data(kh0 + 209);
    const auto *kh0_210 = buffer.data(kh0 + 210);
    const auto *kh0_212 = buffer.data(kh0 + 212);
    const auto *kh0_213 = buffer.data(kh0 + 213);
    const auto *kh0_215 = buffer.data(kh0 + 215);
    const auto *kh0_216 = buffer.data(kh0 + 216);
    const auto *kh0_217 = buffer.data(kh0 + 217);
    const auto *kh0_219 = buffer.data(kh0 + 219);
    const auto *kh0_220 = buffer.data(kh0 + 220);
    const auto *kh0_225 = buffer.data(kh0 + 225);
    const auto *kh0_226 = buffer.data(kh0 + 226);
    const auto *kh0_227 = buffer.data(kh0 + 227);
    const auto *kh0_228 = buffer.data(kh0 + 228);
    const auto *kh0_230 = buffer.data(kh0 + 230);

    const auto *kh1_189 = buffer.data(kh1 + 189);
    const auto *kh1_190 = buffer.data(kh1 + 190);
    const auto *kh1_192 = buffer.data(kh1 + 192);
    const auto *kh1_194 = buffer.data(kh1 + 194);
    const auto *kh1_195 = buffer.data(kh1 + 195);
    const auto *kh1_197 = buffer.data(kh1 + 197);
    const auto *kh1_198 = buffer.data(kh1 + 198);
    const auto *kh1_203 = buffer.data(kh1 + 203);
    const auto *kh1_204 = buffer.data(kh1 + 204);
    const auto *kh1_206 = buffer.data(kh1 + 206);
    const auto *kh1_207 = buffer.data(kh1 + 207);
    const auto *kh1_208 = buffer.data(kh1 + 208);
    const auto *kh1_209 = buffer.data(kh1 + 209);
    const auto *kh1_210 = buffer.data(kh1 + 210);
    const auto *kh1_212 = buffer.data(kh1 + 212);
    const auto *kh1_213 = buffer.data(kh1 + 213);
    const auto *kh1_215 = buffer.data(kh1 + 215);
    const auto *kh1_216 = buffer.data(kh1 + 216);
    const auto *kh1_217 = buffer.data(kh1 + 217);
    const auto *kh1_219 = buffer.data(kh1 + 219);
    const auto *kh1_220 = buffer.data(kh1 + 220);
    const auto *kh1_225 = buffer.data(kh1 + 225);
    const auto *kh1_226 = buffer.data(kh1 + 226);
    const auto *kh1_227 = buffer.data(kh1 + 227);
    const auto *kh1_228 = buffer.data(kh1 + 228);
    const auto *kh1_230 = buffer.data(kh1 + 230);

    const auto *ki_217 = buffer.data(ki + 217);
    const auto *ki_219 = buffer.data(ki + 219);
    const auto *ki_220 = buffer.data(ki + 220);
    const auto *ki_221 = buffer.data(ki + 221);
    const auto *ki_222 = buffer.data(ki + 222);
    const auto *ki_223 = buffer.data(ki + 223);
    const auto *ki_224 = buffer.data(ki + 224);
    const auto *ki_226 = buffer.data(ki + 226);
    const auto *ki_227 = buffer.data(ki + 227);
    const auto *ki_229 = buffer.data(ki + 229);
    const auto *ki_230 = buffer.data(ki + 230);
    const auto *ki_233 = buffer.data(ki + 233);
    const auto *ki_234 = buffer.data(ki + 234);
    const auto *ki_238 = buffer.data(ki + 238);
    const auto *ki_245 = buffer.data(ki + 245);
    const auto *ki_246 = buffer.data(ki + 246);
    const auto *ki_247 = buffer.data(ki + 247);
    const auto *ki_248 = buffer.data(ki + 248);
    const auto *ki_249 = buffer.data(ki + 249);
    const auto *ki_250 = buffer.data(ki + 250);
    const auto *ki_251 = buffer.data(ki + 251);
    const auto *ki_252 = buffer.data(ki + 252);
    const auto *ki_253 = buffer.data(ki + 253);
    const auto *ki_254 = buffer.data(ki + 254);
    const auto *ki_255 = buffer.data(ki + 255);
    const auto *ki_257 = buffer.data(ki + 257);
    const auto *ki_258 = buffer.data(ki + 258);
    const auto *ki_260 = buffer.data(ki + 260);
    const auto *ki_261 = buffer.data(ki + 261);
    const auto *ki_262 = buffer.data(ki + 262);
    const auto *ki_264 = buffer.data(ki + 264);
    const auto *ki_265 = buffer.data(ki + 265);
    const auto *ki_266 = buffer.data(ki + 266);
    const auto *ki_272 = buffer.data(ki + 272);
    const auto *ki_273 = buffer.data(ki + 273);
    const auto *ki_274 = buffer.data(ki + 274);
    const auto *ki_275 = buffer.data(ki + 275);
    const auto *ki_276 = buffer.data(ki + 276);
    const auto *ki_277 = buffer.data(ki + 277);
    const auto *ki_278 = buffer.data(ki + 278);
    const auto *ki_279 = buffer.data(ki + 279);
    const auto *ki_280 = buffer.data(ki + 280);
    const auto *ki_281 = buffer.data(ki + 281);
    const auto *ki_282 = buffer.data(ki + 282);
    const auto *ki_283 = buffer.data(ki + 283);
    const auto *ki_285 = buffer.data(ki + 285);
    const auto *ki_286 = buffer.data(ki + 286);
    const auto *ki_287 = buffer.data(ki + 287);
    const auto *ki_289 = buffer.data(ki + 289);
    const auto *ki_290 = buffer.data(ki + 290);
    const auto *ki_291 = buffer.data(ki + 291);
    const auto *ki_292 = buffer.data(ki + 292);
    const auto *ki_294 = buffer.data(ki + 294);
    const auto *ki_295 = buffer.data(ki + 295);
    const auto *ki_301 = buffer.data(ki + 301);
    const auto *ki_302 = buffer.data(ki + 302);
    const auto *ki_303 = buffer.data(ki + 303);
    const auto *ki_304 = buffer.data(ki + 304);
    const auto *ki_305 = buffer.data(ki + 305);
    const auto *ki_306 = buffer.data(ki + 306);
    const auto *ki_307 = buffer.data(ki + 307);
    const auto *ki_308 = buffer.data(ki + 308);
    const auto *ki_310 = buffer.data(ki + 310);
    const auto *ki_311 = buffer.data(ki + 311);
    const auto *ki_313 = buffer.data(ki + 313);

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, pb_x, ii_219, ii_220, ii_221, \
                         ii_222, ii_223, ki_219, ki_220, ki_221, ki_222, \
                         ki_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_14 * ii_219[k]
                   + pb_x[k] * ki_219[k];

        t_276[k] = f_14 * ii_220[k]
                   + pb_x[k] * ki_220[k];

        t_277[k] = f_14 * ii_221[k]
                   + pb_x[k] * ki_221[k];

        t_278[k] = f_14 * ii_222[k]
                   + pb_x[k] * ki_222[k];

        t_279[k] = f_14 * ii_223[k]
                   + pb_x[k] * ki_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, pa_z, pb_z, ii_105, ii_106, \
                         ii_107, ii_108, ik_136, ik_138, ik_139, ik_140, \
                         ki_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pa_z[k] * ik_136[k];

        t_281[k] = f_11 * ii_105[k]
                   + pb_z[k] * ki_217[k];

        t_282[k] = f_12 * ii_106[k]
                   + pa_z[k] * ik_138[k];

        t_283[k] = f_13 * ii_107[k]
                   + pa_z[k] * ik_139[k];

        t_284[k] = f_14 * ii_108[k]
                   + pa_z[k] * ik_140[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, pa_y, pa_z, pb_y, ii_109, ii_111, ii_139, \
                         ik_141, ik_143, ik_180, ki_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_15 * ii_109[k]
                   + pa_z[k] * ik_141[k];

        t_286[k] = f_12 * ii_139[k]
                   + pb_y[k] * ki_223[k];

        t_287[k] = f_0 * ii_111[k]
                   + pa_z[k] * ik_143[k];

        t_288[k] = pa_y[k] * ik_180[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, t_293, pa_y, pb_y, ii_140, ii_141, \
                         ii_142, ik_182, ik_183, ik_185, ki_224, \
                         ki_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_11 * ii_140[k]
                   + pb_y[k] * ki_224[k];

        t_290[k] = pa_y[k] * ik_182[k];

        t_291[k] = f_12 * ii_141[k]
                   + pa_y[k] * ik_183[k];

        t_292[k] = f_11 * ii_142[k]
                   + pb_y[k] * ki_226[k];

        t_293[k] = pa_y[k] * ik_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_y, pb_y, pb_z, ii_115, ii_143, ii_145, \
                         ik_186, ik_189, ki_227, ki_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_13 * ii_143[k]
                   + pa_y[k] * ik_186[k];

        t_295[k] = f_12 * ii_115[k]
                   + pb_z[k] * ki_227[k];

        t_296[k] = f_11 * ii_145[k]
                   + pb_y[k] * ki_229[k];

        t_297[k] = pa_y[k] * ik_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pa_y, pb_y, pb_z, ii_118, ii_146, ii_148, \
                         ii_149, ik_190, ik_192, ki_230, ki_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_14 * ii_146[k]
                   + pa_y[k] * ik_190[k];

        t_299[k] = f_12 * ii_118[k]
                   + pb_z[k] * ki_230[k];

        t_300[k] = f_12 * ii_148[k]
                   + pa_y[k] * ik_192[k];

        t_301[k] = f_11 * ii_149[k]
                   + pb_y[k] * ki_233[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, pa_y, pb_z, ii_122, ii_150, \
                         ii_152, ii_153, ik_194, ik_195, ik_197, ik_198, \
                         ki_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = pa_y[k] * ik_194[k];

        t_303[k] = f_15 * ii_150[k]
                   + pa_y[k] * ik_195[k];

        t_304[k] = f_12 * ii_122[k]
                   + pb_z[k] * ki_234[k];

        t_305[k] = f_13 * ii_152[k]
                   + pa_y[k] * ik_197[k];

        t_306[k] = f_12 * ii_153[k]
                   + pa_y[k] * ik_198[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pa_y, pb_x, pb_y, ii_154, ii_245, ii_246, \
                         ik_200, ki_238, ki_245, ki_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_11 * ii_154[k]
                   + pb_y[k] * ki_238[k];

        t_308[k] = pa_y[k] * ik_200[k];

        t_309[k] = f_14 * ii_245[k]
                   + pb_x[k] * ki_245[k];

        t_310[k] = f_14 * ii_246[k]
                   + pb_x[k] * ki_246[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, t_315, pa_y, pb_x, ii_247, ii_248, \
                         ii_249, ii_250, ik_207, ki_247, ki_248, ki_249, \
                         ki_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_14 * ii_247[k]
                   + pb_x[k] * ki_247[k];

        t_312[k] = f_14 * ii_248[k]
                   + pb_x[k] * ki_248[k];

        t_313[k] = f_14 * ii_249[k]
                   + pb_x[k] * ki_249[k];

        t_314[k] = f_14 * ii_250[k]
                   + pb_x[k] * ki_250[k];

        t_315[k] = pa_y[k] * ik_207[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pa_y, pb_z, ii_133, ii_161, ii_163, \
                         ii_164, ik_208, ik_210, ik_211, ki_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_0 * ii_161[k]
                   + pa_y[k] * ik_208[k];

        t_317[k] = f_12 * ii_133[k]
                   + pb_z[k] * ki_245[k];

        t_318[k] = f_15 * ii_163[k]
                   + pa_y[k] * ik_210[k];

        t_319[k] = f_14 * ii_164[k]
                   + pa_y[k] * ik_211[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_y, pb_y, ii_165, ii_166, ii_167, \
                         ik_212, ik_213, ik_215, ki_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_13 * ii_165[k]
                   + pa_y[k] * ik_212[k];

        t_321[k] = f_12 * ii_166[k]
                   + pa_y[k] * ik_213[k];

        t_322[k] = f_11 * ii_167[k]
                   + pb_y[k] * ki_251[k];

        t_323[k] = pa_y[k] * ik_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pa_z, pb_y, pb_z, hk0_72, hk1_72, ii_140, \
                         ik_180, kh0_189, kh1_189, ki_252, ki_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_21 * hk0_72[k]
                   - f_22 * hk1_72[k]
                   + pa_z[k] * ik_180[k];

        t_325[k] = pb_y[k] * ki_252[k];

        t_326[k] = f_13 * ii_140[k]
                   + pb_z[k] * ki_252[k];

        t_327[k] = f_3 * kh0_189[k]
                   - f_4 * kh1_189[k]
                   + pb_y[k] * ki_253[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pb_x, pb_y, pb_z, ii_143, ii_257, \
                         kh0_190, kh0_194, kh1_190, kh1_194, ki_254, ki_255, \
                         ki_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = pb_y[k] * ki_254[k];

        t_329[k] = f_14 * ii_257[k]
                   + f_9 * kh0_194[k]
                   - f_10 * kh1_194[k]
                   + pb_x[k] * ki_257[k];

        t_330[k] = f_5 * kh0_190[k]
                   - f_6 * kh1_190[k]
                   + pb_y[k] * ki_255[k];

        t_331[k] = f_13 * ii_143[k]
                   + pb_z[k] * ki_255[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pb_x, pb_y, pb_z, ii_146, ii_261, \
                         kh0_192, kh0_198, kh1_192, kh1_198, ki_257, ki_258, \
                         ki_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pb_y[k] * ki_257[k];

        t_333[k] = f_14 * ii_261[k]
                   + f_7 * kh0_198[k]
                   - f_8 * kh1_198[k]
                   + pb_x[k] * ki_261[k];

        t_334[k] = f_7 * kh0_192[k]
                   - f_8 * kh1_192[k]
                   + pb_y[k] * ki_258[k];

        t_335[k] = f_13 * ii_146[k]
                   + pb_z[k] * ki_258[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pb_x, pb_y, ii_266, kh0_194, kh0_203, kh1_194, \
                         kh1_203, ki_260, ki_261, ki_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_3 * kh0_194[k]
                   - f_4 * kh1_194[k]
                   + pb_y[k] * ki_260[k];

        t_337[k] = pb_y[k] * ki_261[k];

        t_338[k] = f_14 * ii_266[k]
                   + f_5 * kh0_203[k]
                   - f_6 * kh1_203[k]
                   + pb_x[k] * ki_266[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, pb_y, pb_z, ii_150, kh0_195, kh0_197, \
                         kh0_198, kh1_195, kh1_197, kh1_198, ki_262, ki_264, \
                         ki_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_9 * kh0_195[k]
                   - f_10 * kh1_195[k]
                   + pb_y[k] * ki_262[k];

        t_340[k] = f_13 * ii_150[k]
                   + pb_z[k] * ki_262[k];

        t_341[k] = f_5 * kh0_197[k]
                   - f_6 * kh1_197[k]
                   + pb_y[k] * ki_264[k];

        t_342[k] = f_3 * kh0_198[k]
                   - f_4 * kh1_198[k]
                   + pb_y[k] * ki_265[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, pb_x, pb_y, ii_272, ii_273, ii_274, \
                         kh0_209, kh1_209, ki_266, ki_272, ki_273, \
                         ki_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = pb_y[k] * ki_266[k];

        t_344[k] = f_14 * ii_272[k]
                   + f_3 * kh0_209[k]
                   - f_4 * kh1_209[k]
                   + pb_x[k] * ki_272[k];

        t_345[k] = f_14 * ii_273[k]
                   + pb_x[k] * ki_273[k];

        t_346[k] = f_14 * ii_274[k]
                   + pb_x[k] * ki_274[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, pb_x, pb_y, ii_275, ii_276, \
                         ii_277, ii_279, ki_272, ki_275, ki_276, ki_277, \
                         ki_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_14 * ii_275[k]
                   + pb_x[k] * ki_275[k];

        t_348[k] = f_14 * ii_276[k]
                   + pb_x[k] * ki_276[k];

        t_349[k] = f_14 * ii_277[k]
                   + pb_x[k] * ki_277[k];

        t_350[k] = pb_y[k] * ki_272[k];

        t_351[k] = f_14 * ii_279[k]
                   + pb_x[k] * ki_279[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pb_y, pb_z, ii_161, kh0_204, kh0_206, \
                         kh0_207, kh1_204, kh1_206, kh1_207, ki_273, ki_275, \
                         ki_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_1 * kh0_204[k]
                   - f_2 * kh1_204[k]
                   + pb_y[k] * ki_273[k];

        t_353[k] = f_13 * ii_161[k]
                   + pb_z[k] * ki_273[k];

        t_354[k] = f_9 * kh0_206[k]
                   - f_10 * kh1_206[k]
                   + pb_y[k] * ki_275[k];

        t_355[k] = f_7 * kh0_207[k]
                   - f_8 * kh1_207[k]
                   + pb_y[k] * ki_276[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pa_x, pb_y, hk0_359, hk1_359, ik_359, \
                         kh0_208, kh0_209, kh1_208, kh1_209, ki_277, ki_278, \
                         ki_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_5 * kh0_208[k]
                   - f_6 * kh1_208[k]
                   + pb_y[k] * ki_277[k];

        t_357[k] = f_3 * kh0_209[k]
                   - f_4 * kh1_209[k]
                   + pb_y[k] * ki_278[k];

        t_358[k] = pb_y[k] * ki_279[k];

        t_359[k] = f_23 * hk0_359[k]
                   - f_24 * hk1_359[k]
                   + pa_x[k] * ik_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pa_y, pb_y, pb_z, hk0_108, hk1_108, ii_168, \
                         ik_216, ki_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_23 * hk0_108[k]
                   - f_24 * hk1_108[k]
                   + pa_y[k] * ik_216[k];

        t_361[k] = f_14 * ii_168[k]
                   + pb_y[k] * ki_280[k];

        t_362[k] = pb_z[k] * ki_280[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pb_x, pb_z, ii_283, kh0_210, kh0_213, kh1_210, \
                         kh1_213, ki_281, ki_282, ki_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_13 * ii_283[k]
                   + f_9 * kh0_213[k]
                   - f_10 * kh1_213[k]
                   + pb_x[k] * ki_283[k];

        t_364[k] = pb_z[k] * ki_281[k];

        t_365[k] = f_3 * kh0_210[k]
                   - f_4 * kh1_210[k]
                   + pb_z[k] * ki_282[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_x, pb_y, pb_z, ii_173, ii_286, \
                         kh0_212, kh0_216, kh1_212, kh1_216, ki_283, ki_285, \
                         ki_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * ii_286[k]
                   + f_7 * kh0_216[k]
                   - f_8 * kh1_216[k]
                   + pb_x[k] * ki_286[k];

        t_367[k] = pb_z[k] * ki_283[k];

        t_368[k] = f_14 * ii_173[k]
                   + pb_y[k] * ki_285[k];

        t_369[k] = f_5 * kh0_212[k]
                   - f_6 * kh1_212[k]
                   + pb_z[k] * ki_285[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pb_x, pb_z, ii_290, kh0_213, kh0_220, kh1_213, \
                         kh1_220, ki_286, ki_287, ki_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_13 * ii_290[k]
                   + f_5 * kh0_220[k]
                   - f_6 * kh1_220[k]
                   + pb_x[k] * ki_290[k];

        t_371[k] = pb_z[k] * ki_286[k];

        t_372[k] = f_3 * kh0_213[k]
                   - f_4 * kh1_213[k]
                   + pb_z[k] * ki_287[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pb_x, pb_y, pb_z, ii_177, ii_295, \
                         kh0_215, kh0_225, kh1_215, kh1_225, ki_289, ki_290, \
                         ki_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_14 * ii_177[k]
                   + pb_y[k] * ki_289[k];

        t_374[k] = f_7 * kh0_215[k]
                   - f_8 * kh1_215[k]
                   + pb_z[k] * ki_289[k];

        t_375[k] = f_13 * ii_295[k]
                   + f_3 * kh0_225[k]
                   - f_4 * kh1_225[k]
                   + pb_x[k] * ki_295[k];

        t_376[k] = pb_z[k] * ki_290[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pb_y, pb_z, ii_182, kh0_216, kh0_217, \
                         kh0_219, kh1_216, kh1_217, kh1_219, ki_291, ki_292, \
                         ki_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_3 * kh0_216[k]
                   - f_4 * kh1_216[k]
                   + pb_z[k] * ki_291[k];

        t_378[k] = f_5 * kh0_217[k]
                   - f_6 * kh1_217[k]
                   + pb_z[k] * ki_292[k];

        t_379[k] = f_14 * ii_182[k]
                   + pb_y[k] * ki_294[k];

        t_380[k] = f_9 * kh0_219[k]
                   - f_10 * kh1_219[k]
                   + pb_z[k] * ki_294[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, t_385, pb_x, pb_z, ii_301, ii_303, \
                         ii_304, ii_305, ki_295, ki_301, ki_303, ki_304, \
                         ki_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_13 * ii_301[k]
                   + pb_x[k] * ki_301[k];

        t_382[k] = pb_z[k] * ki_295[k];

        t_383[k] = f_13 * ii_303[k]
                   + pb_x[k] * ki_303[k];

        t_384[k] = f_13 * ii_304[k]
                   + pb_x[k] * ki_304[k];

        t_385[k] = f_13 * ii_305[k]
                   + pb_x[k] * ki_305[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pa_x, pb_x, pb_z, hk0_388, hk1_388, \
                         ii_306, ii_307, ik_388, ki_301, ki_306, \
                         ki_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_13 * ii_306[k]
                   + pb_x[k] * ki_306[k];

        t_387[k] = f_13 * ii_307[k]
                   + pb_x[k] * ki_307[k];

        t_388[k] = f_21 * hk0_388[k]
                   - f_22 * hk1_388[k]
                   + pa_x[k] * ik_388[k];

        t_389[k] = pb_z[k] * ki_301[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, pb_z, kh0_225, kh0_226, kh0_227, kh1_225, \
                         kh1_226, kh1_227, ki_302, ki_303, ki_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_3 * kh0_225[k]
                   - f_4 * kh1_225[k]
                   + pb_z[k] * ki_302[k];

        t_391[k] = f_5 * kh0_226[k]
                   - f_6 * kh1_226[k]
                   + pb_z[k] * ki_303[k];

        t_392[k] = f_7 * kh0_227[k]
                   - f_8 * kh1_227[k]
                   + pb_z[k] * ki_304[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, pa_z, pb_y, pb_z, ii_195, ik_216, \
                         kh0_228, kh0_230, kh1_228, kh1_230, ki_305, \
                         ki_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = f_9 * kh0_228[k]
                   - f_10 * kh1_228[k]
                   + pb_z[k] * ki_305[k];

        t_394[k] = f_14 * ii_195[k]
                   + pb_y[k] * ki_307[k];

        t_395[k] = f_1 * kh0_230[k]
                   - f_2 * kh1_230[k]
                   + pb_z[k] * ki_307[k];

        t_396[k] = pa_z[k] * ik_216[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, pa_z, pb_y, pb_z, ii_168, ii_170, \
                         ii_198, ik_217, ik_219, ik_221, ki_308, \
                         ki_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = pa_z[k] * ik_217[k];

        t_398[k] = f_11 * ii_168[k]
                   + pb_z[k] * ki_308[k];

        t_399[k] = pa_z[k] * ik_219[k];

        t_400[k] = f_13 * ii_198[k]
                   + pb_y[k] * ki_310[k];

        t_401[k] = f_12 * ii_170[k]
                   + pa_z[k] * ik_221[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, pa_z, pb_y, pb_z, ii_171, ii_173, \
                         ii_201, ik_222, ik_225, ik_226, ki_311, \
                         ki_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pa_z[k] * ik_222[k];

        t_403[k] = f_11 * ii_171[k]
                   + pb_z[k] * ki_311[k];

        t_404[k] = f_13 * ii_201[k]
                   + pb_y[k] * ki_313[k];

        t_405[k] = f_13 * ii_173[k]
                   + pa_z[k] * ik_225[k];

        t_406[k] = pa_z[k] * ik_226[k];
    }
}

static auto
compute_prim_kk_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t hk0,
                                            const size_t hk1, const size_t ii, const size_t ik,
                                            const size_t kh0, const size_t kh1, const size_t ki,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
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
    const auto f_17 = 0.5 / alpha;
    const auto f_18 = 0.5 * beta / (alpha * p);
    const auto f_21 = 1.0 / alpha;
    const auto f_22 = beta / (alpha * p);
    const auto f_23 = 1.5 / alpha;
    const auto f_24 = 1.5 * beta / (alpha * p);

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

    const auto *hk0_111 = buffer.data(hk0 + 111);
    const auto *hk0_114 = buffer.data(hk0 + 114);
    const auto *hk0_118 = buffer.data(hk0 + 118);
    const auto *hk0_123 = buffer.data(hk0 + 123);
    const auto *hk0_180 = buffer.data(hk0 + 180);
    const auto *hk0_185 = buffer.data(hk0 + 185);
    const auto *hk0_189 = buffer.data(hk0 + 189);
    const auto *hk0_194 = buffer.data(hk0 + 194);
    const auto *hk0_200 = buffer.data(hk0 + 200);
    const auto *hk0_460 = buffer.data(hk0 + 460);
    const auto *hk0_462 = buffer.data(hk0 + 462);
    const auto *hk0_463 = buffer.data(hk0 + 463);
    const auto *hk0_464 = buffer.data(hk0 + 464);
    const auto *hk0_465 = buffer.data(hk0 + 465);
    const auto *hk0_467 = buffer.data(hk0 + 467);

    const auto *hk1_111 = buffer.data(hk1 + 111);
    const auto *hk1_114 = buffer.data(hk1 + 114);
    const auto *hk1_118 = buffer.data(hk1 + 118);
    const auto *hk1_123 = buffer.data(hk1 + 123);
    const auto *hk1_180 = buffer.data(hk1 + 180);
    const auto *hk1_185 = buffer.data(hk1 + 185);
    const auto *hk1_189 = buffer.data(hk1 + 189);
    const auto *hk1_194 = buffer.data(hk1 + 194);
    const auto *hk1_200 = buffer.data(hk1 + 200);
    const auto *hk1_460 = buffer.data(hk1 + 460);
    const auto *hk1_462 = buffer.data(hk1 + 462);
    const auto *hk1_463 = buffer.data(hk1 + 463);
    const auto *hk1_464 = buffer.data(hk1 + 464);
    const auto *hk1_465 = buffer.data(hk1 + 465);
    const auto *hk1_467 = buffer.data(hk1 + 467);

    const auto *ii_174 = buffer.data(ii + 174);
    const auto *ii_175 = buffer.data(ii + 175);
    const auto *ii_177 = buffer.data(ii + 177);
    const auto *ii_178 = buffer.data(ii + 178);
    const auto *ii_179 = buffer.data(ii + 179);
    const auto *ii_180 = buffer.data(ii + 180);
    const auto *ii_182 = buffer.data(ii + 182);
    const auto *ii_189 = buffer.data(ii + 189);
    const auto *ii_190 = buffer.data(ii + 190);
    const auto *ii_191 = buffer.data(ii + 191);
    const auto *ii_192 = buffer.data(ii + 192);
    const auto *ii_193 = buffer.data(ii + 193);
    const auto *ii_195 = buffer.data(ii + 195);
    const auto *ii_196 = buffer.data(ii + 196);
    const auto *ii_199 = buffer.data(ii + 199);
    const auto *ii_202 = buffer.data(ii + 202);
    const auto *ii_205 = buffer.data(ii + 205);
    const auto *ii_206 = buffer.data(ii + 206);
    const auto *ii_210 = buffer.data(ii + 210);
    const auto *ii_217 = buffer.data(ii + 217);
    const auto *ii_223 = buffer.data(ii + 223);
    const auto *ii_224 = buffer.data(ii + 224);
    const auto *ii_226 = buffer.data(ii + 226);
    const auto *ii_227 = buffer.data(ii + 227);
    const auto *ii_229 = buffer.data(ii + 229);
    const auto *ii_230 = buffer.data(ii + 230);
    const auto *ii_233 = buffer.data(ii + 233);
    const auto *ii_234 = buffer.data(ii + 234);
    const auto *ii_238 = buffer.data(ii + 238);
    const auto *ii_245 = buffer.data(ii + 245);
    const auto *ii_251 = buffer.data(ii + 251);
    const auto *ii_252 = buffer.data(ii + 252);
    const auto *ii_253 = buffer.data(ii + 253);
    const auto *ii_254 = buffer.data(ii + 254);
    const auto *ii_255 = buffer.data(ii + 255);
    const auto *ii_257 = buffer.data(ii + 257);
    const auto *ii_258 = buffer.data(ii + 258);
    const auto *ii_260 = buffer.data(ii + 260);
    const auto *ii_261 = buffer.data(ii + 261);
    const auto *ii_262 = buffer.data(ii + 262);
    const auto *ii_264 = buffer.data(ii + 264);
    const auto *ii_265 = buffer.data(ii + 265);
    const auto *ii_266 = buffer.data(ii + 266);
    const auto *ii_273 = buffer.data(ii + 273);
    const auto *ii_275 = buffer.data(ii + 275);
    const auto *ii_276 = buffer.data(ii + 276);
    const auto *ii_277 = buffer.data(ii + 277);
    const auto *ii_278 = buffer.data(ii + 278);
    const auto *ii_279 = buffer.data(ii + 279);
    const auto *ii_330 = buffer.data(ii + 330);
    const auto *ii_331 = buffer.data(ii + 331);
    const auto *ii_332 = buffer.data(ii + 332);
    const auto *ii_333 = buffer.data(ii + 333);
    const auto *ii_334 = buffer.data(ii + 334);
    const auto *ii_335 = buffer.data(ii + 335);
    const auto *ii_348 = buffer.data(ii + 348);
    const auto *ii_353 = buffer.data(ii + 353);
    const auto *ii_354 = buffer.data(ii + 354);
    const auto *ii_357 = buffer.data(ii + 357);
    const auto *ii_358 = buffer.data(ii + 358);
    const auto *ii_359 = buffer.data(ii + 359);
    const auto *ii_360 = buffer.data(ii + 360);
    const auto *ii_361 = buffer.data(ii + 361);
    const auto *ii_362 = buffer.data(ii + 362);
    const auto *ii_363 = buffer.data(ii + 363);
    const auto *ii_385 = buffer.data(ii + 385);
    const auto *ii_386 = buffer.data(ii + 386);
    const auto *ii_387 = buffer.data(ii + 387);
    const auto *ii_388 = buffer.data(ii + 388);
    const auto *ii_389 = buffer.data(ii + 389);
    const auto *ii_390 = buffer.data(ii + 390);
    const auto *ii_397 = buffer.data(ii + 397);
    const auto *ii_401 = buffer.data(ii + 401);
    const auto *ii_406 = buffer.data(ii + 406);
    const auto *ii_412 = buffer.data(ii + 412);
    const auto *ii_413 = buffer.data(ii + 413);
    const auto *ii_414 = buffer.data(ii + 414);
    const auto *ii_415 = buffer.data(ii + 415);
    const auto *ii_416 = buffer.data(ii + 416);
    const auto *ii_417 = buffer.data(ii + 417);
    const auto *ii_419 = buffer.data(ii + 419);

    const auto *ik_228 = buffer.data(ik + 228);
    const auto *ik_230 = buffer.data(ik + 230);
    const auto *ik_231 = buffer.data(ik + 231);
    const auto *ik_233 = buffer.data(ik + 233);
    const auto *ik_234 = buffer.data(ik + 234);
    const auto *ik_236 = buffer.data(ik + 236);
    const auto *ik_237 = buffer.data(ik + 237);
    const auto *ik_244 = buffer.data(ik + 244);
    const auto *ik_246 = buffer.data(ik + 246);
    const auto *ik_247 = buffer.data(ik + 247);
    const auto *ik_248 = buffer.data(ik + 248);
    const auto *ik_249 = buffer.data(ik + 249);
    const auto *ik_251 = buffer.data(ik + 251);
    const auto *ik_255 = buffer.data(ik + 255);
    const auto *ik_258 = buffer.data(ik + 258);
    const auto *ik_262 = buffer.data(ik + 262);
    const auto *ik_267 = buffer.data(ik + 267);
    const auto *ik_288 = buffer.data(ik + 288);
    const auto *ik_293 = buffer.data(ik + 293);
    const auto *ik_297 = buffer.data(ik + 297);
    const auto *ik_302 = buffer.data(ik + 302);
    const auto *ik_308 = buffer.data(ik + 308);
    const auto *ik_324 = buffer.data(ik + 324);
    const auto *ik_326 = buffer.data(ik + 326);
    const auto *ik_327 = buffer.data(ik + 327);
    const auto *ik_329 = buffer.data(ik + 329);
    const auto *ik_330 = buffer.data(ik + 330);
    const auto *ik_333 = buffer.data(ik + 333);
    const auto *ik_334 = buffer.data(ik + 334);
    const auto *ik_336 = buffer.data(ik + 336);
    const auto *ik_338 = buffer.data(ik + 338);
    const auto *ik_339 = buffer.data(ik + 339);
    const auto *ik_341 = buffer.data(ik + 341);
    const auto *ik_342 = buffer.data(ik + 342);
    const auto *ik_344 = buffer.data(ik + 344);
    const auto *ik_351 = buffer.data(ik + 351);
    const auto *ik_352 = buffer.data(ik + 352);
    const auto *ik_354 = buffer.data(ik + 354);
    const auto *ik_355 = buffer.data(ik + 355);
    const auto *ik_356 = buffer.data(ik + 356);
    const auto *ik_357 = buffer.data(ik + 357);
    const auto *ik_359 = buffer.data(ik + 359);
    const auto *ik_460 = buffer.data(ik + 460);
    const auto *ik_462 = buffer.data(ik + 462);
    const auto *ik_463 = buffer.data(ik + 463);
    const auto *ik_464 = buffer.data(ik + 464);
    const auto *ik_465 = buffer.data(ik + 465);
    const auto *ik_467 = buffer.data(ik + 467);

    const auto *kh0_264 = buffer.data(kh0 + 264);
    const auto *kh0_269 = buffer.data(kh0 + 269);
    const auto *kh0_270 = buffer.data(kh0 + 270);
    const auto *kh0_294 = buffer.data(kh0 + 294);
    const auto *kh0_295 = buffer.data(kh0 + 295);
    const auto *kh0_297 = buffer.data(kh0 + 297);
    const auto *kh0_299 = buffer.data(kh0 + 299);
    const auto *kh0_300 = buffer.data(kh0 + 300);
    const auto *kh0_302 = buffer.data(kh0 + 302);
    const auto *kh0_303 = buffer.data(kh0 + 303);
    const auto *kh0_308 = buffer.data(kh0 + 308);
    const auto *kh0_309 = buffer.data(kh0 + 309);
    const auto *kh0_311 = buffer.data(kh0 + 311);
    const auto *kh0_312 = buffer.data(kh0 + 312);
    const auto *kh0_314 = buffer.data(kh0 + 314);

    const auto *kh1_264 = buffer.data(kh1 + 264);
    const auto *kh1_269 = buffer.data(kh1 + 269);
    const auto *kh1_270 = buffer.data(kh1 + 270);
    const auto *kh1_294 = buffer.data(kh1 + 294);
    const auto *kh1_295 = buffer.data(kh1 + 295);
    const auto *kh1_297 = buffer.data(kh1 + 297);
    const auto *kh1_299 = buffer.data(kh1 + 299);
    const auto *kh1_300 = buffer.data(kh1 + 300);
    const auto *kh1_302 = buffer.data(kh1 + 302);
    const auto *kh1_303 = buffer.data(kh1 + 303);
    const auto *kh1_308 = buffer.data(kh1 + 308);
    const auto *kh1_309 = buffer.data(kh1 + 309);
    const auto *kh1_311 = buffer.data(kh1 + 311);
    const auto *kh1_312 = buffer.data(kh1 + 312);
    const auto *kh1_314 = buffer.data(kh1 + 314);

    const auto *ki_314 = buffer.data(ki + 314);
    const auto *ki_317 = buffer.data(ki + 317);
    const auto *ki_318 = buffer.data(ki + 318);
    const auto *ki_322 = buffer.data(ki + 322);
    const auto *ki_329 = buffer.data(ki + 329);
    const auto *ki_330 = buffer.data(ki + 330);
    const auto *ki_331 = buffer.data(ki + 331);
    const auto *ki_332 = buffer.data(ki + 332);
    const auto *ki_333 = buffer.data(ki + 333);
    const auto *ki_334 = buffer.data(ki + 334);
    const auto *ki_335 = buffer.data(ki + 335);
    const auto *ki_336 = buffer.data(ki + 336);
    const auto *ki_338 = buffer.data(ki + 338);
    const auto *ki_339 = buffer.data(ki + 339);
    const auto *ki_341 = buffer.data(ki + 341);
    const auto *ki_342 = buffer.data(ki + 342);
    const auto *ki_345 = buffer.data(ki + 345);
    const auto *ki_346 = buffer.data(ki + 346);
    const auto *ki_348 = buffer.data(ki + 348);
    const auto *ki_350 = buffer.data(ki + 350);
    const auto *ki_353 = buffer.data(ki + 353);
    const auto *ki_354 = buffer.data(ki + 354);
    const auto *ki_357 = buffer.data(ki + 357);
    const auto *ki_358 = buffer.data(ki + 358);
    const auto *ki_359 = buffer.data(ki + 359);
    const auto *ki_360 = buffer.data(ki + 360);
    const auto *ki_361 = buffer.data(ki + 361);
    const auto *ki_362 = buffer.data(ki + 362);
    const auto *ki_363 = buffer.data(ki + 363);
    const auto *ki_364 = buffer.data(ki + 364);
    const auto *ki_366 = buffer.data(ki + 366);
    const auto *ki_367 = buffer.data(ki + 367);
    const auto *ki_369 = buffer.data(ki + 369);
    const auto *ki_370 = buffer.data(ki + 370);
    const auto *ki_373 = buffer.data(ki + 373);
    const auto *ki_374 = buffer.data(ki + 374);
    const auto *ki_378 = buffer.data(ki + 378);
    const auto *ki_385 = buffer.data(ki + 385);
    const auto *ki_386 = buffer.data(ki + 386);
    const auto *ki_387 = buffer.data(ki + 387);
    const auto *ki_388 = buffer.data(ki + 388);
    const auto *ki_389 = buffer.data(ki + 389);
    const auto *ki_390 = buffer.data(ki + 390);
    const auto *ki_391 = buffer.data(ki + 391);
    const auto *ki_392 = buffer.data(ki + 392);
    const auto *ki_393 = buffer.data(ki + 393);
    const auto *ki_394 = buffer.data(ki + 394);
    const auto *ki_395 = buffer.data(ki + 395);
    const auto *ki_397 = buffer.data(ki + 397);
    const auto *ki_398 = buffer.data(ki + 398);
    const auto *ki_400 = buffer.data(ki + 400);
    const auto *ki_401 = buffer.data(ki + 401);
    const auto *ki_402 = buffer.data(ki + 402);
    const auto *ki_404 = buffer.data(ki + 404);
    const auto *ki_405 = buffer.data(ki + 405);
    const auto *ki_406 = buffer.data(ki + 406);
    const auto *ki_412 = buffer.data(ki + 412);
    const auto *ki_413 = buffer.data(ki + 413);
    const auto *ki_414 = buffer.data(ki + 414);
    const auto *ki_415 = buffer.data(ki + 415);
    const auto *ki_416 = buffer.data(ki + 416);
    const auto *ki_417 = buffer.data(ki + 417);
    const auto *ki_419 = buffer.data(ki + 419);

#pragma omp simd aligned(t_407, t_408, t_409, t_410, pa_z, pb_y, pb_z, ii_174, ii_175, ii_177, \
                         ii_205, ik_228, ik_230, ki_314, ki_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_11 * ii_174[k]
                   + pb_z[k] * ki_314[k];

        t_408[k] = f_12 * ii_175[k]
                   + pa_z[k] * ik_228[k];

        t_409[k] = f_13 * ii_205[k]
                   + pb_y[k] * ki_317[k];

        t_410[k] = f_14 * ii_177[k]
                   + pa_z[k] * ik_230[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, pa_z, pb_z, ii_178, ii_179, ii_180, \
                         ik_231, ik_233, ik_234, ki_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = pa_z[k] * ik_231[k];

        t_412[k] = f_11 * ii_178[k]
                   + pb_z[k] * ki_318[k];

        t_413[k] = f_12 * ii_179[k]
                   + pa_z[k] * ik_233[k];

        t_414[k] = f_13 * ii_180[k]
                   + pa_z[k] * ik_234[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, pa_z, pb_x, pb_y, ii_182, ii_210, ii_330, \
                         ik_236, ik_237, ki_322, ki_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_13 * ii_210[k]
                   + pb_y[k] * ki_322[k];

        t_416[k] = f_15 * ii_182[k]
                   + pa_z[k] * ik_236[k];

        t_417[k] = pa_z[k] * ik_237[k];

        t_418[k] = f_13 * ii_330[k]
                   + pb_x[k] * ki_330[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, t_423, pb_x, ii_331, ii_332, ii_333, \
                         ii_334, ii_335, ki_331, ki_332, ki_333, ki_334, \
                         ki_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_13 * ii_331[k]
                   + pb_x[k] * ki_331[k];

        t_420[k] = f_13 * ii_332[k]
                   + pb_x[k] * ki_332[k];

        t_421[k] = f_13 * ii_333[k]
                   + pb_x[k] * ki_333[k];

        t_422[k] = f_13 * ii_334[k]
                   + pb_x[k] * ki_334[k];

        t_423[k] = f_13 * ii_335[k]
                   + pb_x[k] * ki_335[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, t_428, pa_z, pb_z, ii_189, ii_190, \
                         ii_191, ii_192, ik_244, ik_246, ik_247, ik_248, \
                         ki_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pa_z[k] * ik_244[k];

        t_425[k] = f_11 * ii_189[k]
                   + pb_z[k] * ki_329[k];

        t_426[k] = f_12 * ii_190[k]
                   + pa_z[k] * ik_246[k];

        t_427[k] = f_13 * ii_191[k]
                   + pa_z[k] * ik_247[k];

        t_428[k] = f_14 * ii_192[k]
                   + pa_z[k] * ik_248[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pa_y, pa_z, pb_y, hk0_180, hk1_180, \
                         ii_193, ii_195, ii_223, ik_249, ik_251, ik_288, \
                         ki_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_15 * ii_193[k]
                   + pa_z[k] * ik_249[k];

        t_430[k] = f_13 * ii_223[k]
                   + pb_y[k] * ki_335[k];

        t_431[k] = f_0 * ii_195[k]
                   + pa_z[k] * ik_251[k];

        t_432[k] = f_17 * hk0_180[k]
                   - f_18 * hk1_180[k]
                   + pa_y[k] * ik_288[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pa_z, pb_y, pb_z, hk0_111, hk1_111, \
                         ii_196, ii_224, ii_226, ik_255, ki_336, \
                         ki_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_12 * ii_224[k]
                   + pb_y[k] * ki_336[k];

        t_434[k] = f_12 * ii_196[k]
                   + pb_z[k] * ki_336[k];

        t_435[k] = f_17 * hk0_111[k]
                   - f_18 * hk1_111[k]
                   + pa_z[k] * ik_255[k];

        t_436[k] = f_12 * ii_226[k]
                   + pb_y[k] * ki_338[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pa_y, pa_z, pb_z, hk0_114, hk0_185, hk1_114, \
                         hk1_185, ii_199, ik_258, ik_293, ki_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_17 * hk0_185[k]
                   - f_18 * hk1_185[k]
                   + pa_y[k] * ik_293[k];

        t_438[k] = f_17 * hk0_114[k]
                   - f_18 * hk1_114[k]
                   + pa_z[k] * ik_258[k];

        t_439[k] = f_12 * ii_199[k]
                   + pb_z[k] * ki_339[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pa_y, pa_z, pb_y, hk0_118, hk0_189, hk1_118, \
                         hk1_189, ii_229, ik_262, ik_297, ki_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_12 * ii_229[k]
                   + pb_y[k] * ki_341[k];

        t_441[k] = f_17 * hk0_189[k]
                   - f_18 * hk1_189[k]
                   + pa_y[k] * ik_297[k];

        t_442[k] = f_17 * hk0_118[k]
                   - f_18 * hk1_118[k]
                   + pa_z[k] * ik_262[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pb_x, pb_y, pb_z, ii_202, ii_233, ii_348, \
                         kh0_264, kh1_264, ki_342, ki_345, ki_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_12 * ii_202[k]
                   + pb_z[k] * ki_342[k];

        t_444[k] = f_13 * ii_348[k]
                   + f_5 * kh0_264[k]
                   - f_6 * kh1_264[k]
                   + pb_x[k] * ki_348[k];

        t_445[k] = f_12 * ii_233[k]
                   + pb_y[k] * ki_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pa_y, pa_z, pb_z, hk0_123, hk0_194, hk1_123, \
                         hk1_194, ii_206, ik_267, ik_302, ki_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_17 * hk0_194[k]
                   - f_18 * hk1_194[k]
                   + pa_y[k] * ik_302[k];

        t_447[k] = f_17 * hk0_123[k]
                   - f_18 * hk1_123[k]
                   + pa_z[k] * ik_267[k];

        t_448[k] = f_12 * ii_206[k]
                   + pb_z[k] * ki_346[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pb_x, pb_y, ii_238, ii_353, ii_354, kh0_269, \
                         kh0_270, kh1_269, kh1_270, ki_350, ki_353, \
                         ki_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_13 * ii_353[k]
                   + f_3 * kh0_269[k]
                   - f_4 * kh1_269[k]
                   + pb_x[k] * ki_353[k];

        t_450[k] = f_13 * ii_354[k]
                   + f_3 * kh0_270[k]
                   - f_4 * kh1_270[k]
                   + pb_x[k] * ki_354[k];

        t_451[k] = f_12 * ii_238[k]
                   + pb_y[k] * ki_350[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pa_y, pb_x, hk0_200, hk1_200, ii_357, \
                         ii_358, ii_359, ik_308, ki_357, ki_358, \
                         ki_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_17 * hk0_200[k]
                   - f_18 * hk1_200[k]
                   + pa_y[k] * ik_308[k];

        t_453[k] = f_13 * ii_357[k]
                   + pb_x[k] * ki_357[k];

        t_454[k] = f_13 * ii_358[k]
                   + pb_x[k] * ki_358[k];

        t_455[k] = f_13 * ii_359[k]
                   + pb_x[k] * ki_359[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pb_x, ii_360, ii_361, ii_362, ii_363, \
                         ki_360, ki_361, ki_362, ki_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_13 * ii_360[k]
                   + pb_x[k] * ki_360[k];

        t_457[k] = f_13 * ii_361[k]
                   + pb_x[k] * ki_361[k];

        t_458[k] = f_13 * ii_362[k]
                   + pb_x[k] * ki_362[k];

        t_459[k] = f_13 * ii_363[k]
                   + pb_x[k] * ki_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pa_x, pb_z, hk0_460, hk0_462, hk1_460, hk1_462, \
                         ii_217, ik_460, ik_462, ki_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_21 * hk0_460[k]
                   - f_22 * hk1_460[k]
                   + pa_x[k] * ik_460[k];

        t_461[k] = f_12 * ii_217[k]
                   + pb_z[k] * ki_357[k];

        t_462[k] = f_21 * hk0_462[k]
                   - f_22 * hk1_462[k]
                   + pa_x[k] * ik_462[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pa_x, hk0_463, hk0_464, hk0_465, hk1_463, \
                         hk1_464, hk1_465, ik_463, ik_464, ik_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_21 * hk0_463[k]
                   - f_22 * hk1_463[k]
                   + pa_x[k] * ik_463[k];

        t_464[k] = f_21 * hk0_464[k]
                   - f_22 * hk1_464[k]
                   + pa_x[k] * ik_464[k];

        t_465[k] = f_21 * hk0_465[k]
                   - f_22 * hk1_465[k]
                   + pa_x[k] * ik_465[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_x, pa_y, pb_y, hk0_467, hk1_467, \
                         ii_251, ii_252, ik_324, ik_467, ki_363, \
                         ki_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_12 * ii_251[k]
                   + pb_y[k] * ki_363[k];

        t_467[k] = f_21 * hk0_467[k]
                   - f_22 * hk1_467[k]
                   + pa_x[k] * ik_467[k];

        t_468[k] = pa_y[k] * ik_324[k];

        t_469[k] = f_11 * ii_252[k]
                   + pb_y[k] * ki_364[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, pa_y, pb_y, ii_253, ii_254, \
                         ii_255, ik_326, ik_327, ik_329, ik_330, \
                         ki_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = pa_y[k] * ik_326[k];

        t_471[k] = f_12 * ii_253[k]
                   + pa_y[k] * ik_327[k];

        t_472[k] = f_11 * ii_254[k]
                   + pb_y[k] * ki_366[k];

        t_473[k] = pa_y[k] * ik_329[k];

        t_474[k] = f_13 * ii_255[k]
                   + pa_y[k] * ik_330[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_y, pb_y, pb_z, ii_227, ii_257, ii_258, \
                         ik_333, ik_334, ki_367, ki_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_13 * ii_227[k]
                   + pb_z[k] * ki_367[k];

        t_476[k] = f_11 * ii_257[k]
                   + pb_y[k] * ki_369[k];

        t_477[k] = pa_y[k] * ik_333[k];

        t_478[k] = f_14 * ii_258[k]
                   + pa_y[k] * ik_334[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, pa_y, pb_y, pb_z, ii_230, ii_260, ii_261, \
                         ik_336, ik_338, ki_370, ki_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_13 * ii_230[k]
                   + pb_z[k] * ki_370[k];

        t_480[k] = f_12 * ii_260[k]
                   + pa_y[k] * ik_336[k];

        t_481[k] = f_11 * ii_261[k]
                   + pb_y[k] * ki_373[k];

        t_482[k] = pa_y[k] * ik_338[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pa_y, pb_z, ii_234, ii_262, ii_264, \
                         ii_265, ik_339, ik_341, ik_342, ki_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_15 * ii_262[k]
                   + pa_y[k] * ik_339[k];

        t_484[k] = f_13 * ii_234[k]
                   + pb_z[k] * ki_374[k];

        t_485[k] = f_13 * ii_264[k]
                   + pa_y[k] * ik_341[k];

        t_486[k] = f_12 * ii_265[k]
                   + pa_y[k] * ik_342[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, pa_y, pb_x, pb_y, ii_266, ii_385, ii_386, \
                         ik_344, ki_378, ki_385, ki_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_11 * ii_266[k]
                   + pb_y[k] * ki_378[k];

        t_488[k] = pa_y[k] * ik_344[k];

        t_489[k] = f_13 * ii_385[k]
                   + pb_x[k] * ki_385[k];

        t_490[k] = f_13 * ii_386[k]
                   + pb_x[k] * ki_386[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pa_y, pb_x, ii_387, ii_388, \
                         ii_389, ii_390, ik_351, ki_387, ki_388, ki_389, \
                         ki_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_13 * ii_387[k]
                   + pb_x[k] * ki_387[k];

        t_492[k] = f_13 * ii_388[k]
                   + pb_x[k] * ki_388[k];

        t_493[k] = f_13 * ii_389[k]
                   + pb_x[k] * ki_389[k];

        t_494[k] = f_13 * ii_390[k]
                   + pb_x[k] * ki_390[k];

        t_495[k] = pa_y[k] * ik_351[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pa_y, pb_z, ii_245, ii_273, ii_275, \
                         ii_276, ik_352, ik_354, ik_355, ki_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_0 * ii_273[k]
                   + pa_y[k] * ik_352[k];

        t_497[k] = f_13 * ii_245[k]
                   + pb_z[k] * ki_385[k];

        t_498[k] = f_15 * ii_275[k]
                   + pa_y[k] * ik_354[k];

        t_499[k] = f_14 * ii_276[k]
                   + pa_y[k] * ik_355[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pa_y, pb_y, ii_277, ii_278, ii_279, \
                         ik_356, ik_357, ik_359, ki_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_13 * ii_277[k]
                   + pa_y[k] * ik_356[k];

        t_501[k] = f_12 * ii_278[k]
                   + pa_y[k] * ik_357[k];

        t_502[k] = f_11 * ii_279[k]
                   + pb_y[k] * ki_391[k];

        t_503[k] = pa_y[k] * ik_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pa_z, pb_y, pb_z, hk0_180, hk1_180, \
                         ii_252, ik_324, kh0_294, kh1_294, ki_392, \
                         ki_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_23 * hk0_180[k]
                   - f_24 * hk1_180[k]
                   + pa_z[k] * ik_324[k];

        t_505[k] = pb_y[k] * ki_392[k];

        t_506[k] = f_14 * ii_252[k]
                   + pb_z[k] * ki_392[k];

        t_507[k] = f_3 * kh0_294[k]
                   - f_4 * kh1_294[k]
                   + pb_y[k] * ki_393[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pb_x, pb_y, pb_z, ii_255, ii_397, \
                         kh0_295, kh0_299, kh1_295, kh1_299, ki_394, ki_395, \
                         ki_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = pb_y[k] * ki_394[k];

        t_509[k] = f_13 * ii_397[k]
                   + f_9 * kh0_299[k]
                   - f_10 * kh1_299[k]
                   + pb_x[k] * ki_397[k];

        t_510[k] = f_5 * kh0_295[k]
                   - f_6 * kh1_295[k]
                   + pb_y[k] * ki_395[k];

        t_511[k] = f_14 * ii_255[k]
                   + pb_z[k] * ki_395[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pb_x, pb_y, pb_z, ii_258, ii_401, \
                         kh0_297, kh0_303, kh1_297, kh1_303, ki_397, ki_398, \
                         ki_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pb_y[k] * ki_397[k];

        t_513[k] = f_13 * ii_401[k]
                   + f_7 * kh0_303[k]
                   - f_8 * kh1_303[k]
                   + pb_x[k] * ki_401[k];

        t_514[k] = f_7 * kh0_297[k]
                   - f_8 * kh1_297[k]
                   + pb_y[k] * ki_398[k];

        t_515[k] = f_14 * ii_258[k]
                   + pb_z[k] * ki_398[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, pb_x, pb_y, ii_406, kh0_299, kh0_308, kh1_299, \
                         kh1_308, ki_400, ki_401, ki_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_3 * kh0_299[k]
                   - f_4 * kh1_299[k]
                   + pb_y[k] * ki_400[k];

        t_517[k] = pb_y[k] * ki_401[k];

        t_518[k] = f_13 * ii_406[k]
                   + f_5 * kh0_308[k]
                   - f_6 * kh1_308[k]
                   + pb_x[k] * ki_406[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, pb_y, pb_z, ii_262, kh0_300, kh0_302, \
                         kh0_303, kh1_300, kh1_302, kh1_303, ki_402, ki_404, \
                         ki_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_9 * kh0_300[k]
                   - f_10 * kh1_300[k]
                   + pb_y[k] * ki_402[k];

        t_520[k] = f_14 * ii_262[k]
                   + pb_z[k] * ki_402[k];

        t_521[k] = f_5 * kh0_302[k]
                   - f_6 * kh1_302[k]
                   + pb_y[k] * ki_404[k];

        t_522[k] = f_3 * kh0_303[k]
                   - f_4 * kh1_303[k]
                   + pb_y[k] * ki_405[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pb_x, pb_y, ii_412, ii_413, ii_414, \
                         kh0_314, kh1_314, ki_406, ki_412, ki_413, \
                         ki_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = pb_y[k] * ki_406[k];

        t_524[k] = f_13 * ii_412[k]
                   + f_3 * kh0_314[k]
                   - f_4 * kh1_314[k]
                   + pb_x[k] * ki_412[k];

        t_525[k] = f_13 * ii_413[k]
                   + pb_x[k] * ki_413[k];

        t_526[k] = f_13 * ii_414[k]
                   + pb_x[k] * ki_414[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, pb_x, pb_y, ii_415, ii_416, \
                         ii_417, ii_419, ki_412, ki_415, ki_416, ki_417, \
                         ki_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_13 * ii_415[k]
                   + pb_x[k] * ki_415[k];

        t_528[k] = f_13 * ii_416[k]
                   + pb_x[k] * ki_416[k];

        t_529[k] = f_13 * ii_417[k]
                   + pb_x[k] * ki_417[k];

        t_530[k] = pb_y[k] * ki_412[k];

        t_531[k] = f_13 * ii_419[k]
                   + pb_x[k] * ki_419[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pb_y, pb_z, ii_273, kh0_309, kh0_311, \
                         kh0_312, kh1_309, kh1_311, kh1_312, ki_413, ki_415, \
                         ki_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_1 * kh0_309[k]
                   - f_2 * kh1_309[k]
                   + pb_y[k] * ki_413[k];

        t_533[k] = f_14 * ii_273[k]
                   + pb_z[k] * ki_413[k];

        t_534[k] = f_9 * kh0_311[k]
                   - f_10 * kh1_311[k]
                   + pb_y[k] * ki_415[k];

        t_535[k] = f_7 * kh0_312[k]
                   - f_8 * kh1_312[k]
                   + pb_y[k] * ki_416[k];
    }
}

static auto
compute_prim_kk_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t hk0,
                                            const size_t hk1, const size_t ii, const size_t ik,
                                            const size_t kh0, const size_t kh1, const size_t ki,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
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
    const auto f_17 = 0.5 / alpha;
    const auto f_18 = 0.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / alpha;
    const auto f_20 = 2.0 * beta / (alpha * p);
    const auto f_21 = 1.0 / alpha;
    const auto f_22 = beta / (alpha * p);

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

    const auto *hk0_216 = buffer.data(hk0 + 216);
    const auto *hk0_219 = buffer.data(hk0 + 219);
    const auto *hk0_222 = buffer.data(hk0 + 222);
    const auto *hk0_226 = buffer.data(hk0 + 226);
    const auto *hk0_231 = buffer.data(hk0 + 231);
    const auto *hk0_255 = buffer.data(hk0 + 255);
    const auto *hk0_258 = buffer.data(hk0 + 258);
    const auto *hk0_262 = buffer.data(hk0 + 262);
    const auto *hk0_288 = buffer.data(hk0 + 288);
    const auto *hk0_293 = buffer.data(hk0 + 293);
    const auto *hk0_297 = buffer.data(hk0 + 297);
    const auto *hk0_302 = buffer.data(hk0 + 302);
    const auto *hk0_308 = buffer.data(hk0 + 308);
    const auto *hk0_324 = buffer.data(hk0 + 324);
    const auto *hk0_329 = buffer.data(hk0 + 329);
    const auto *hk0_333 = buffer.data(hk0 + 333);
    const auto *hk0_539 = buffer.data(hk0 + 539);
    const auto *hk0_568 = buffer.data(hk0 + 568);
    const auto *hk0_640 = buffer.data(hk0 + 640);
    const auto *hk0_642 = buffer.data(hk0 + 642);
    const auto *hk0_643 = buffer.data(hk0 + 643);
    const auto *hk0_644 = buffer.data(hk0 + 644);
    const auto *hk0_645 = buffer.data(hk0 + 645);
    const auto *hk0_647 = buffer.data(hk0 + 647);

    const auto *hk1_216 = buffer.data(hk1 + 216);
    const auto *hk1_219 = buffer.data(hk1 + 219);
    const auto *hk1_222 = buffer.data(hk1 + 222);
    const auto *hk1_226 = buffer.data(hk1 + 226);
    const auto *hk1_231 = buffer.data(hk1 + 231);
    const auto *hk1_255 = buffer.data(hk1 + 255);
    const auto *hk1_258 = buffer.data(hk1 + 258);
    const auto *hk1_262 = buffer.data(hk1 + 262);
    const auto *hk1_288 = buffer.data(hk1 + 288);
    const auto *hk1_293 = buffer.data(hk1 + 293);
    const auto *hk1_297 = buffer.data(hk1 + 297);
    const auto *hk1_302 = buffer.data(hk1 + 302);
    const auto *hk1_308 = buffer.data(hk1 + 308);
    const auto *hk1_324 = buffer.data(hk1 + 324);
    const auto *hk1_329 = buffer.data(hk1 + 329);
    const auto *hk1_333 = buffer.data(hk1 + 333);
    const auto *hk1_539 = buffer.data(hk1 + 539);
    const auto *hk1_568 = buffer.data(hk1 + 568);
    const auto *hk1_640 = buffer.data(hk1 + 640);
    const auto *hk1_642 = buffer.data(hk1 + 642);
    const auto *hk1_643 = buffer.data(hk1 + 643);
    const auto *hk1_644 = buffer.data(hk1 + 644);
    const auto *hk1_645 = buffer.data(hk1 + 645);
    const auto *hk1_647 = buffer.data(hk1 + 647);

    const auto *ii_280 = buffer.data(ii + 280);
    const auto *ii_282 = buffer.data(ii + 282);
    const auto *ii_283 = buffer.data(ii + 283);
    const auto *ii_285 = buffer.data(ii + 285);
    const auto *ii_286 = buffer.data(ii + 286);
    const auto *ii_287 = buffer.data(ii + 287);
    const auto *ii_289 = buffer.data(ii + 289);
    const auto *ii_290 = buffer.data(ii + 290);
    const auto *ii_291 = buffer.data(ii + 291);
    const auto *ii_292 = buffer.data(ii + 292);
    const auto *ii_294 = buffer.data(ii + 294);
    const auto *ii_301 = buffer.data(ii + 301);
    const auto *ii_302 = buffer.data(ii + 302);
    const auto *ii_303 = buffer.data(ii + 303);
    const auto *ii_304 = buffer.data(ii + 304);
    const auto *ii_305 = buffer.data(ii + 305);
    const auto *ii_307 = buffer.data(ii + 307);
    const auto *ii_308 = buffer.data(ii + 308);
    const auto *ii_310 = buffer.data(ii + 310);
    const auto *ii_311 = buffer.data(ii + 311);
    const auto *ii_313 = buffer.data(ii + 313);
    const auto *ii_314 = buffer.data(ii + 314);
    const auto *ii_317 = buffer.data(ii + 317);
    const auto *ii_318 = buffer.data(ii + 318);
    const auto *ii_322 = buffer.data(ii + 322);
    const auto *ii_329 = buffer.data(ii + 329);
    const auto *ii_335 = buffer.data(ii + 335);
    const auto *ii_336 = buffer.data(ii + 336);
    const auto *ii_338 = buffer.data(ii + 338);
    const auto *ii_339 = buffer.data(ii + 339);
    const auto *ii_341 = buffer.data(ii + 341);
    const auto *ii_345 = buffer.data(ii + 345);
    const auto *ii_350 = buffer.data(ii + 350);
    const auto *ii_363 = buffer.data(ii + 363);
    const auto *ii_364 = buffer.data(ii + 364);
    const auto *ii_366 = buffer.data(ii + 366);
    const auto *ii_369 = buffer.data(ii + 369);
    const auto *ii_423 = buffer.data(ii + 423);
    const auto *ii_426 = buffer.data(ii + 426);
    const auto *ii_430 = buffer.data(ii + 430);
    const auto *ii_435 = buffer.data(ii + 435);
    const auto *ii_441 = buffer.data(ii + 441);
    const auto *ii_443 = buffer.data(ii + 443);
    const auto *ii_444 = buffer.data(ii + 444);
    const auto *ii_445 = buffer.data(ii + 445);
    const auto *ii_446 = buffer.data(ii + 446);
    const auto *ii_447 = buffer.data(ii + 447);
    const auto *ii_470 = buffer.data(ii + 470);
    const auto *ii_471 = buffer.data(ii + 471);
    const auto *ii_472 = buffer.data(ii + 472);
    const auto *ii_473 = buffer.data(ii + 473);
    const auto *ii_474 = buffer.data(ii + 474);
    const auto *ii_475 = buffer.data(ii + 475);
    const auto *ii_488 = buffer.data(ii + 488);
    const auto *ii_493 = buffer.data(ii + 493);
    const auto *ii_494 = buffer.data(ii + 494);
    const auto *ii_497 = buffer.data(ii + 497);
    const auto *ii_498 = buffer.data(ii + 498);
    const auto *ii_499 = buffer.data(ii + 499);
    const auto *ii_500 = buffer.data(ii + 500);
    const auto *ii_501 = buffer.data(ii + 501);
    const auto *ii_502 = buffer.data(ii + 502);
    const auto *ii_503 = buffer.data(ii + 503);

    const auto *ik_360 = buffer.data(ik + 360);
    const auto *ik_361 = buffer.data(ik + 361);
    const auto *ik_363 = buffer.data(ik + 363);
    const auto *ik_365 = buffer.data(ik + 365);
    const auto *ik_366 = buffer.data(ik + 366);
    const auto *ik_369 = buffer.data(ik + 369);
    const auto *ik_370 = buffer.data(ik + 370);
    const auto *ik_372 = buffer.data(ik + 372);
    const auto *ik_374 = buffer.data(ik + 374);
    const auto *ik_375 = buffer.data(ik + 375);
    const auto *ik_377 = buffer.data(ik + 377);
    const auto *ik_378 = buffer.data(ik + 378);
    const auto *ik_380 = buffer.data(ik + 380);
    const auto *ik_381 = buffer.data(ik + 381);
    const auto *ik_388 = buffer.data(ik + 388);
    const auto *ik_390 = buffer.data(ik + 390);
    const auto *ik_391 = buffer.data(ik + 391);
    const auto *ik_392 = buffer.data(ik + 392);
    const auto *ik_393 = buffer.data(ik + 393);
    const auto *ik_395 = buffer.data(ik + 395);
    const auto *ik_399 = buffer.data(ik + 399);
    const auto *ik_402 = buffer.data(ik + 402);
    const auto *ik_406 = buffer.data(ik + 406);
    const auto *ik_411 = buffer.data(ik + 411);
    const auto *ik_432 = buffer.data(ik + 432);
    const auto *ik_435 = buffer.data(ik + 435);
    const auto *ik_437 = buffer.data(ik + 437);
    const auto *ik_438 = buffer.data(ik + 438);
    const auto *ik_441 = buffer.data(ik + 441);
    const auto *ik_442 = buffer.data(ik + 442);
    const auto *ik_446 = buffer.data(ik + 446);
    const auto *ik_452 = buffer.data(ik + 452);
    const auto *ik_468 = buffer.data(ik + 468);
    const auto *ik_473 = buffer.data(ik + 473);
    const auto *ik_477 = buffer.data(ik + 477);
    const auto *ik_539 = buffer.data(ik + 539);
    const auto *ik_568 = buffer.data(ik + 568);
    const auto *ik_640 = buffer.data(ik + 640);
    const auto *ik_642 = buffer.data(ik + 642);
    const auto *ik_643 = buffer.data(ik + 643);
    const auto *ik_644 = buffer.data(ik + 644);
    const auto *ik_645 = buffer.data(ik + 645);
    const auto *ik_647 = buffer.data(ik + 647);

    const auto *kh0_313 = buffer.data(kh0 + 313);
    const auto *kh0_314 = buffer.data(kh0 + 314);
    const auto *kh0_315 = buffer.data(kh0 + 315);
    const auto *kh0_317 = buffer.data(kh0 + 317);
    const auto *kh0_318 = buffer.data(kh0 + 318);
    const auto *kh0_320 = buffer.data(kh0 + 320);
    const auto *kh0_321 = buffer.data(kh0 + 321);
    const auto *kh0_322 = buffer.data(kh0 + 322);
    const auto *kh0_324 = buffer.data(kh0 + 324);
    const auto *kh0_325 = buffer.data(kh0 + 325);
    const auto *kh0_330 = buffer.data(kh0 + 330);
    const auto *kh0_331 = buffer.data(kh0 + 331);
    const auto *kh0_332 = buffer.data(kh0 + 332);
    const auto *kh0_333 = buffer.data(kh0 + 333);
    const auto *kh0_335 = buffer.data(kh0 + 335);
    const auto *kh0_369 = buffer.data(kh0 + 369);
    const auto *kh0_374 = buffer.data(kh0 + 374);
    const auto *kh0_375 = buffer.data(kh0 + 375);

    const auto *kh1_313 = buffer.data(kh1 + 313);
    const auto *kh1_314 = buffer.data(kh1 + 314);
    const auto *kh1_315 = buffer.data(kh1 + 315);
    const auto *kh1_317 = buffer.data(kh1 + 317);
    const auto *kh1_318 = buffer.data(kh1 + 318);
    const auto *kh1_320 = buffer.data(kh1 + 320);
    const auto *kh1_321 = buffer.data(kh1 + 321);
    const auto *kh1_322 = buffer.data(kh1 + 322);
    const auto *kh1_324 = buffer.data(kh1 + 324);
    const auto *kh1_325 = buffer.data(kh1 + 325);
    const auto *kh1_330 = buffer.data(kh1 + 330);
    const auto *kh1_331 = buffer.data(kh1 + 331);
    const auto *kh1_332 = buffer.data(kh1 + 332);
    const auto *kh1_333 = buffer.data(kh1 + 333);
    const auto *kh1_335 = buffer.data(kh1 + 335);
    const auto *kh1_369 = buffer.data(kh1 + 369);
    const auto *kh1_374 = buffer.data(kh1 + 374);
    const auto *kh1_375 = buffer.data(kh1 + 375);

    const auto *ki_417 = buffer.data(ki + 417);
    const auto *ki_418 = buffer.data(ki + 418);
    const auto *ki_419 = buffer.data(ki + 419);
    const auto *ki_420 = buffer.data(ki + 420);
    const auto *ki_421 = buffer.data(ki + 421);
    const auto *ki_422 = buffer.data(ki + 422);
    const auto *ki_423 = buffer.data(ki + 423);
    const auto *ki_425 = buffer.data(ki + 425);
    const auto *ki_426 = buffer.data(ki + 426);
    const auto *ki_427 = buffer.data(ki + 427);
    const auto *ki_429 = buffer.data(ki + 429);
    const auto *ki_430 = buffer.data(ki + 430);
    const auto *ki_431 = buffer.data(ki + 431);
    const auto *ki_432 = buffer.data(ki + 432);
    const auto *ki_434 = buffer.data(ki + 434);
    const auto *ki_435 = buffer.data(ki + 435);
    const auto *ki_441 = buffer.data(ki + 441);
    const auto *ki_442 = buffer.data(ki + 442);
    const auto *ki_443 = buffer.data(ki + 443);
    const auto *ki_444 = buffer.data(ki + 444);
    const auto *ki_445 = buffer.data(ki + 445);
    const auto *ki_446 = buffer.data(ki + 446);
    const auto *ki_447 = buffer.data(ki + 447);
    const auto *ki_448 = buffer.data(ki + 448);
    const auto *ki_450 = buffer.data(ki + 450);
    const auto *ki_451 = buffer.data(ki + 451);
    const auto *ki_453 = buffer.data(ki + 453);
    const auto *ki_454 = buffer.data(ki + 454);
    const auto *ki_457 = buffer.data(ki + 457);
    const auto *ki_458 = buffer.data(ki + 458);
    const auto *ki_462 = buffer.data(ki + 462);
    const auto *ki_469 = buffer.data(ki + 469);
    const auto *ki_470 = buffer.data(ki + 470);
    const auto *ki_471 = buffer.data(ki + 471);
    const auto *ki_472 = buffer.data(ki + 472);
    const auto *ki_473 = buffer.data(ki + 473);
    const auto *ki_474 = buffer.data(ki + 474);
    const auto *ki_475 = buffer.data(ki + 475);
    const auto *ki_476 = buffer.data(ki + 476);
    const auto *ki_478 = buffer.data(ki + 478);
    const auto *ki_479 = buffer.data(ki + 479);
    const auto *ki_481 = buffer.data(ki + 481);
    const auto *ki_482 = buffer.data(ki + 482);
    const auto *ki_485 = buffer.data(ki + 485);
    const auto *ki_486 = buffer.data(ki + 486);
    const auto *ki_488 = buffer.data(ki + 488);
    const auto *ki_490 = buffer.data(ki + 490);
    const auto *ki_493 = buffer.data(ki + 493);
    const auto *ki_494 = buffer.data(ki + 494);
    const auto *ki_497 = buffer.data(ki + 497);
    const auto *ki_498 = buffer.data(ki + 498);
    const auto *ki_499 = buffer.data(ki + 499);
    const auto *ki_500 = buffer.data(ki + 500);
    const auto *ki_501 = buffer.data(ki + 501);
    const auto *ki_502 = buffer.data(ki + 502);
    const auto *ki_503 = buffer.data(ki + 503);
    const auto *ki_504 = buffer.data(ki + 504);
    const auto *ki_506 = buffer.data(ki + 506);
    const auto *ki_507 = buffer.data(ki + 507);
    const auto *ki_509 = buffer.data(ki + 509);

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pa_x, pb_y, hk0_539, hk1_539, ik_539, \
                         kh0_313, kh0_314, kh1_313, kh1_314, ki_417, ki_418, \
                         ki_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_5 * kh0_313[k]
                   - f_6 * kh1_313[k]
                   + pb_y[k] * ki_417[k];

        t_537[k] = f_3 * kh0_314[k]
                   - f_4 * kh1_314[k]
                   + pb_y[k] * ki_418[k];

        t_538[k] = pb_y[k] * ki_419[k];

        t_539[k] = f_21 * hk0_539[k]
                   - f_22 * hk1_539[k]
                   + pa_x[k] * ik_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, pa_y, pb_y, pb_z, hk0_216, hk1_216, ii_280, \
                         ik_360, ki_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_19 * hk0_216[k]
                   - f_20 * hk1_216[k]
                   + pa_y[k] * ik_360[k];

        t_541[k] = f_15 * ii_280[k]
                   + pb_y[k] * ki_420[k];

        t_542[k] = pb_z[k] * ki_420[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, pb_x, pb_z, ii_423, kh0_315, kh0_318, kh1_315, \
                         kh1_318, ki_421, ki_422, ki_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_12 * ii_423[k]
                   + f_9 * kh0_318[k]
                   - f_10 * kh1_318[k]
                   + pb_x[k] * ki_423[k];

        t_544[k] = pb_z[k] * ki_421[k];

        t_545[k] = f_3 * kh0_315[k]
                   - f_4 * kh1_315[k]
                   + pb_z[k] * ki_422[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, pb_x, pb_y, pb_z, ii_285, ii_426, \
                         kh0_317, kh0_321, kh1_317, kh1_321, ki_423, ki_425, \
                         ki_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = f_12 * ii_426[k]
                   + f_7 * kh0_321[k]
                   - f_8 * kh1_321[k]
                   + pb_x[k] * ki_426[k];

        t_547[k] = pb_z[k] * ki_423[k];

        t_548[k] = f_15 * ii_285[k]
                   + pb_y[k] * ki_425[k];

        t_549[k] = f_5 * kh0_317[k]
                   - f_6 * kh1_317[k]
                   + pb_z[k] * ki_425[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, pb_x, pb_z, ii_430, kh0_318, kh0_325, kh1_318, \
                         kh1_325, ki_426, ki_427, ki_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_12 * ii_430[k]
                   + f_5 * kh0_325[k]
                   - f_6 * kh1_325[k]
                   + pb_x[k] * ki_430[k];

        t_551[k] = pb_z[k] * ki_426[k];

        t_552[k] = f_3 * kh0_318[k]
                   - f_4 * kh1_318[k]
                   + pb_z[k] * ki_427[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, pb_x, pb_y, pb_z, ii_289, ii_435, \
                         kh0_320, kh0_330, kh1_320, kh1_330, ki_429, ki_430, \
                         ki_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_15 * ii_289[k]
                   + pb_y[k] * ki_429[k];

        t_554[k] = f_7 * kh0_320[k]
                   - f_8 * kh1_320[k]
                   + pb_z[k] * ki_429[k];

        t_555[k] = f_12 * ii_435[k]
                   + f_3 * kh0_330[k]
                   - f_4 * kh1_330[k]
                   + pb_x[k] * ki_435[k];

        t_556[k] = pb_z[k] * ki_430[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, pb_y, pb_z, ii_294, kh0_321, kh0_322, \
                         kh0_324, kh1_321, kh1_322, kh1_324, ki_431, ki_432, \
                         ki_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_3 * kh0_321[k]
                   - f_4 * kh1_321[k]
                   + pb_z[k] * ki_431[k];

        t_558[k] = f_5 * kh0_322[k]
                   - f_6 * kh1_322[k]
                   + pb_z[k] * ki_432[k];

        t_559[k] = f_15 * ii_294[k]
                   + pb_y[k] * ki_434[k];

        t_560[k] = f_9 * kh0_324[k]
                   - f_10 * kh1_324[k]
                   + pb_z[k] * ki_434[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, t_564, t_565, pb_x, pb_z, ii_441, ii_443, \
                         ii_444, ii_445, ki_435, ki_441, ki_443, ki_444, \
                         ki_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_12 * ii_441[k]
                   + pb_x[k] * ki_441[k];

        t_562[k] = pb_z[k] * ki_435[k];

        t_563[k] = f_12 * ii_443[k]
                   + pb_x[k] * ki_443[k];

        t_564[k] = f_12 * ii_444[k]
                   + pb_x[k] * ki_444[k];

        t_565[k] = f_12 * ii_445[k]
                   + pb_x[k] * ki_445[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pa_x, pb_x, pb_z, hk0_568, hk1_568, \
                         ii_446, ii_447, ik_568, ki_441, ki_446, \
                         ki_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_12 * ii_446[k]
                   + pb_x[k] * ki_446[k];

        t_567[k] = f_12 * ii_447[k]
                   + pb_x[k] * ki_447[k];

        t_568[k] = f_17 * hk0_568[k]
                   - f_18 * hk1_568[k]
                   + pa_x[k] * ik_568[k];

        t_569[k] = pb_z[k] * ki_441[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, pb_z, kh0_330, kh0_331, kh0_332, kh1_330, \
                         kh1_331, kh1_332, ki_442, ki_443, ki_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_3 * kh0_330[k]
                   - f_4 * kh1_330[k]
                   + pb_z[k] * ki_442[k];

        t_571[k] = f_5 * kh0_331[k]
                   - f_6 * kh1_331[k]
                   + pb_z[k] * ki_443[k];

        t_572[k] = f_7 * kh0_332[k]
                   - f_8 * kh1_332[k]
                   + pb_z[k] * ki_444[k];
    }

#pragma omp simd aligned(t_573, t_574, t_575, t_576, pa_z, pb_y, pb_z, ii_307, ik_360, \
                         kh0_333, kh0_335, kh1_333, kh1_335, ki_445, \
                         ki_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_573[k] = f_9 * kh0_333[k]
                   - f_10 * kh1_333[k]
                   + pb_z[k] * ki_445[k];

        t_574[k] = f_15 * ii_307[k]
                   + pb_y[k] * ki_447[k];

        t_575[k] = f_1 * kh0_335[k]
                   - f_2 * kh1_335[k]
                   + pb_z[k] * ki_447[k];

        t_576[k] = pa_z[k] * ik_360[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, pa_z, pb_y, pb_z, ii_280, ii_282, \
                         ii_310, ik_361, ik_363, ik_365, ki_448, \
                         ki_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = pa_z[k] * ik_361[k];

        t_578[k] = f_11 * ii_280[k]
                   + pb_z[k] * ki_448[k];

        t_579[k] = pa_z[k] * ik_363[k];

        t_580[k] = f_14 * ii_310[k]
                   + pb_y[k] * ki_450[k];

        t_581[k] = f_12 * ii_282[k]
                   + pa_z[k] * ik_365[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, pa_z, pb_y, pb_z, ii_283, ii_285, \
                         ii_313, ik_366, ik_369, ik_370, ki_451, \
                         ki_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = pa_z[k] * ik_366[k];

        t_583[k] = f_11 * ii_283[k]
                   + pb_z[k] * ki_451[k];

        t_584[k] = f_14 * ii_313[k]
                   + pb_y[k] * ki_453[k];

        t_585[k] = f_13 * ii_285[k]
                   + pa_z[k] * ik_369[k];

        t_586[k] = pa_z[k] * ik_370[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, t_590, pa_z, pb_y, pb_z, ii_286, ii_287, ii_289, \
                         ii_317, ik_372, ik_374, ki_454, ki_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_11 * ii_286[k]
                   + pb_z[k] * ki_454[k];

        t_588[k] = f_12 * ii_287[k]
                   + pa_z[k] * ik_372[k];

        t_589[k] = f_14 * ii_317[k]
                   + pb_y[k] * ki_457[k];

        t_590[k] = f_14 * ii_289[k]
                   + pa_z[k] * ik_374[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, pa_z, pb_z, ii_290, ii_291, ii_292, \
                         ik_375, ik_377, ik_378, ki_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = pa_z[k] * ik_375[k];

        t_592[k] = f_11 * ii_290[k]
                   + pb_z[k] * ki_458[k];

        t_593[k] = f_12 * ii_291[k]
                   + pa_z[k] * ik_377[k];

        t_594[k] = f_13 * ii_292[k]
                   + pa_z[k] * ik_378[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, pa_z, pb_x, pb_y, ii_294, ii_322, ii_470, \
                         ik_380, ik_381, ki_462, ki_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_14 * ii_322[k]
                   + pb_y[k] * ki_462[k];

        t_596[k] = f_15 * ii_294[k]
                   + pa_z[k] * ik_380[k];

        t_597[k] = pa_z[k] * ik_381[k];

        t_598[k] = f_12 * ii_470[k]
                   + pb_x[k] * ki_470[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, t_603, pb_x, ii_471, ii_472, ii_473, \
                         ii_474, ii_475, ki_471, ki_472, ki_473, ki_474, \
                         ki_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = f_12 * ii_471[k]
                   + pb_x[k] * ki_471[k];

        t_600[k] = f_12 * ii_472[k]
                   + pb_x[k] * ki_472[k];

        t_601[k] = f_12 * ii_473[k]
                   + pb_x[k] * ki_473[k];

        t_602[k] = f_12 * ii_474[k]
                   + pb_x[k] * ki_474[k];

        t_603[k] = f_12 * ii_475[k]
                   + pb_x[k] * ki_475[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, pa_z, pb_z, ii_301, ii_302, \
                         ii_303, ii_304, ik_388, ik_390, ik_391, ik_392, \
                         ki_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = pa_z[k] * ik_388[k];

        t_605[k] = f_11 * ii_301[k]
                   + pb_z[k] * ki_469[k];

        t_606[k] = f_12 * ii_302[k]
                   + pa_z[k] * ik_390[k];

        t_607[k] = f_13 * ii_303[k]
                   + pa_z[k] * ik_391[k];

        t_608[k] = f_14 * ii_304[k]
                   + pa_z[k] * ik_392[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, pa_y, pa_z, pb_y, hk0_288, hk1_288, \
                         ii_305, ii_307, ii_335, ik_393, ik_395, ik_432, \
                         ki_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_15 * ii_305[k]
                   + pa_z[k] * ik_393[k];

        t_610[k] = f_14 * ii_335[k]
                   + pb_y[k] * ki_475[k];

        t_611[k] = f_0 * ii_307[k]
                   + pa_z[k] * ik_395[k];

        t_612[k] = f_21 * hk0_288[k]
                   - f_22 * hk1_288[k]
                   + pa_y[k] * ik_432[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pa_z, pb_y, pb_z, hk0_219, hk1_219, \
                         ii_308, ii_336, ii_338, ik_399, ki_476, \
                         ki_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_13 * ii_336[k]
                   + pb_y[k] * ki_476[k];

        t_614[k] = f_12 * ii_308[k]
                   + pb_z[k] * ki_476[k];

        t_615[k] = f_17 * hk0_219[k]
                   - f_18 * hk1_219[k]
                   + pa_z[k] * ik_399[k];

        t_616[k] = f_13 * ii_338[k]
                   + pb_y[k] * ki_478[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, pa_y, pa_z, pb_z, hk0_222, hk0_293, hk1_222, \
                         hk1_293, ii_311, ik_402, ik_437, ki_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_21 * hk0_293[k]
                   - f_22 * hk1_293[k]
                   + pa_y[k] * ik_437[k];

        t_618[k] = f_17 * hk0_222[k]
                   - f_18 * hk1_222[k]
                   + pa_z[k] * ik_402[k];

        t_619[k] = f_12 * ii_311[k]
                   + pb_z[k] * ki_479[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pa_y, pa_z, pb_y, hk0_226, hk0_297, hk1_226, \
                         hk1_297, ii_341, ik_406, ik_441, ki_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_13 * ii_341[k]
                   + pb_y[k] * ki_481[k];

        t_621[k] = f_21 * hk0_297[k]
                   - f_22 * hk1_297[k]
                   + pa_y[k] * ik_441[k];

        t_622[k] = f_17 * hk0_226[k]
                   - f_18 * hk1_226[k]
                   + pa_z[k] * ik_406[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pb_x, pb_y, pb_z, ii_314, ii_345, ii_488, \
                         kh0_369, kh1_369, ki_482, ki_485, ki_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_12 * ii_314[k]
                   + pb_z[k] * ki_482[k];

        t_624[k] = f_12 * ii_488[k]
                   + f_5 * kh0_369[k]
                   - f_6 * kh1_369[k]
                   + pb_x[k] * ki_488[k];

        t_625[k] = f_13 * ii_345[k]
                   + pb_y[k] * ki_485[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pa_y, pa_z, pb_z, hk0_231, hk0_302, hk1_231, \
                         hk1_302, ii_318, ik_411, ik_446, ki_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_21 * hk0_302[k]
                   - f_22 * hk1_302[k]
                   + pa_y[k] * ik_446[k];

        t_627[k] = f_17 * hk0_231[k]
                   - f_18 * hk1_231[k]
                   + pa_z[k] * ik_411[k];

        t_628[k] = f_12 * ii_318[k]
                   + pb_z[k] * ki_486[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pb_x, pb_y, ii_350, ii_493, ii_494, kh0_374, \
                         kh0_375, kh1_374, kh1_375, ki_490, ki_493, \
                         ki_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_12 * ii_493[k]
                   + f_3 * kh0_374[k]
                   - f_4 * kh1_374[k]
                   + pb_x[k] * ki_493[k];

        t_630[k] = f_12 * ii_494[k]
                   + f_3 * kh0_375[k]
                   - f_4 * kh1_375[k]
                   + pb_x[k] * ki_494[k];

        t_631[k] = f_13 * ii_350[k]
                   + pb_y[k] * ki_490[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, pa_y, pb_x, hk0_308, hk1_308, ii_497, \
                         ii_498, ii_499, ik_452, ki_497, ki_498, \
                         ki_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_21 * hk0_308[k]
                   - f_22 * hk1_308[k]
                   + pa_y[k] * ik_452[k];

        t_633[k] = f_12 * ii_497[k]
                   + pb_x[k] * ki_497[k];

        t_634[k] = f_12 * ii_498[k]
                   + pb_x[k] * ki_498[k];

        t_635[k] = f_12 * ii_499[k]
                   + pb_x[k] * ki_499[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pb_x, ii_500, ii_501, ii_502, ii_503, \
                         ki_500, ki_501, ki_502, ki_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_12 * ii_500[k]
                   + pb_x[k] * ki_500[k];

        t_637[k] = f_12 * ii_501[k]
                   + pb_x[k] * ki_501[k];

        t_638[k] = f_12 * ii_502[k]
                   + pb_x[k] * ki_502[k];

        t_639[k] = f_12 * ii_503[k]
                   + pb_x[k] * ki_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pa_x, pb_z, hk0_640, hk0_642, hk1_640, hk1_642, \
                         ii_329, ik_640, ik_642, ki_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_17 * hk0_640[k]
                   - f_18 * hk1_640[k]
                   + pa_x[k] * ik_640[k];

        t_641[k] = f_12 * ii_329[k]
                   + pb_z[k] * ki_497[k];

        t_642[k] = f_17 * hk0_642[k]
                   - f_18 * hk1_642[k]
                   + pa_x[k] * ik_642[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pa_x, hk0_643, hk0_644, hk0_645, hk1_643, \
                         hk1_644, hk1_645, ik_643, ik_644, ik_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_17 * hk0_643[k]
                   - f_18 * hk1_643[k]
                   + pa_x[k] * ik_643[k];

        t_644[k] = f_17 * hk0_644[k]
                   - f_18 * hk1_644[k]
                   + pa_x[k] * ik_644[k];

        t_645[k] = f_17 * hk0_645[k]
                   - f_18 * hk1_645[k]
                   + pa_x[k] * ik_645[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pa_x, pa_y, pb_y, hk0_324, hk0_647, hk1_324, \
                         hk1_647, ii_363, ik_468, ik_647, ki_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_13 * ii_363[k]
                   + pb_y[k] * ki_503[k];

        t_647[k] = f_17 * hk0_647[k]
                   - f_18 * hk1_647[k]
                   + pa_x[k] * ik_647[k];

        t_648[k] = f_17 * hk0_324[k]
                   - f_18 * hk1_324[k]
                   + pa_y[k] * ik_468[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pa_z, pb_y, pb_z, hk0_255, hk1_255, \
                         ii_336, ii_364, ii_366, ik_435, ki_504, \
                         ki_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_12 * ii_364[k]
                   + pb_y[k] * ki_504[k];

        t_650[k] = f_13 * ii_336[k]
                   + pb_z[k] * ki_504[k];

        t_651[k] = f_21 * hk0_255[k]
                   - f_22 * hk1_255[k]
                   + pa_z[k] * ik_435[k];

        t_652[k] = f_12 * ii_366[k]
                   + pb_y[k] * ki_506[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pa_y, pa_z, pb_z, hk0_258, hk0_329, hk1_258, \
                         hk1_329, ii_339, ik_438, ik_473, ki_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_17 * hk0_329[k]
                   - f_18 * hk1_329[k]
                   + pa_y[k] * ik_473[k];

        t_654[k] = f_21 * hk0_258[k]
                   - f_22 * hk1_258[k]
                   + pa_z[k] * ik_438[k];

        t_655[k] = f_13 * ii_339[k]
                   + pb_z[k] * ki_507[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pa_y, pa_z, pb_y, hk0_262, hk0_333, hk1_262, \
                         hk1_333, ii_369, ik_442, ik_477, ki_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_12 * ii_369[k]
                   + pb_y[k] * ki_509[k];

        t_657[k] = f_17 * hk0_333[k]
                   - f_18 * hk1_333[k]
                   + pa_y[k] * ik_477[k];

        t_658[k] = f_21 * hk0_262[k]
                   - f_22 * hk1_262[k]
                   + pa_z[k] * ik_442[k];
    }
}

static auto
compute_prim_kk_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t hk0,
                                            const size_t hk1, const size_t ii, const size_t ik,
                                            const size_t kh0, const size_t kh1, const size_t ki,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
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
    const auto f_16 = 3.0 / p;
    const auto f_17 = 0.5 / alpha;
    const auto f_18 = 0.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / alpha;
    const auto f_20 = 2.0 * beta / (alpha * p);
    const auto f_21 = 1.0 / alpha;
    const auto f_22 = beta / (alpha * p);

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
    auto *t_786 = buffer.data(target + 786);
    auto *t_787 = buffer.data(target + 787);
    auto *t_788 = buffer.data(target + 788);
    auto *t_789 = buffer.data(target + 789);
    auto *t_790 = buffer.data(target + 790);
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hk0_267 = buffer.data(hk0 + 267);
    const auto *hk0_324 = buffer.data(hk0 + 324);
    const auto *hk0_338 = buffer.data(hk0 + 338);
    const auto *hk0_344 = buffer.data(hk0 + 344);
    const auto *hk0_676 = buffer.data(hk0 + 676);
    const auto *hk0_678 = buffer.data(hk0 + 678);
    const auto *hk0_679 = buffer.data(hk0 + 679);
    const auto *hk0_680 = buffer.data(hk0 + 680);
    const auto *hk0_681 = buffer.data(hk0 + 681);
    const auto *hk0_683 = buffer.data(hk0 + 683);
    const auto *hk0_755 = buffer.data(hk0 + 755);

    const auto *hk1_267 = buffer.data(hk1 + 267);
    const auto *hk1_324 = buffer.data(hk1 + 324);
    const auto *hk1_338 = buffer.data(hk1 + 338);
    const auto *hk1_344 = buffer.data(hk1 + 344);
    const auto *hk1_676 = buffer.data(hk1 + 676);
    const auto *hk1_678 = buffer.data(hk1 + 678);
    const auto *hk1_679 = buffer.data(hk1 + 679);
    const auto *hk1_680 = buffer.data(hk1 + 680);
    const auto *hk1_681 = buffer.data(hk1 + 681);
    const auto *hk1_683 = buffer.data(hk1 + 683);
    const auto *hk1_755 = buffer.data(hk1 + 755);

    const auto *ii_342 = buffer.data(ii + 342);
    const auto *ii_346 = buffer.data(ii + 346);
    const auto *ii_357 = buffer.data(ii + 357);
    const auto *ii_367 = buffer.data(ii + 367);
    const auto *ii_370 = buffer.data(ii + 370);
    const auto *ii_373 = buffer.data(ii + 373);
    const auto *ii_374 = buffer.data(ii + 374);
    const auto *ii_378 = buffer.data(ii + 378);
    const auto *ii_385 = buffer.data(ii + 385);
    const auto *ii_391 = buffer.data(ii + 391);
    const auto *ii_392 = buffer.data(ii + 392);
    const auto *ii_393 = buffer.data(ii + 393);
    const auto *ii_394 = buffer.data(ii + 394);
    const auto *ii_395 = buffer.data(ii + 395);
    const auto *ii_397 = buffer.data(ii + 397);
    const auto *ii_398 = buffer.data(ii + 398);
    const auto *ii_400 = buffer.data(ii + 400);
    const auto *ii_401 = buffer.data(ii + 401);
    const auto *ii_402 = buffer.data(ii + 402);
    const auto *ii_404 = buffer.data(ii + 404);
    const auto *ii_405 = buffer.data(ii + 405);
    const auto *ii_406 = buffer.data(ii + 406);
    const auto *ii_413 = buffer.data(ii + 413);
    const auto *ii_415 = buffer.data(ii + 415);
    const auto *ii_416 = buffer.data(ii + 416);
    const auto *ii_417 = buffer.data(ii + 417);
    const auto *ii_418 = buffer.data(ii + 418);
    const auto *ii_419 = buffer.data(ii + 419);
    const auto *ii_420 = buffer.data(ii + 420);
    const auto *ii_425 = buffer.data(ii + 425);
    const auto *ii_429 = buffer.data(ii + 429);
    const auto *ii_434 = buffer.data(ii + 434);
    const auto *ii_516 = buffer.data(ii + 516);
    const auto *ii_521 = buffer.data(ii + 521);
    const auto *ii_522 = buffer.data(ii + 522);
    const auto *ii_525 = buffer.data(ii + 525);
    const auto *ii_526 = buffer.data(ii + 526);
    const auto *ii_527 = buffer.data(ii + 527);
    const auto *ii_528 = buffer.data(ii + 528);
    const auto *ii_529 = buffer.data(ii + 529);
    const auto *ii_530 = buffer.data(ii + 530);
    const auto *ii_531 = buffer.data(ii + 531);
    const auto *ii_553 = buffer.data(ii + 553);
    const auto *ii_554 = buffer.data(ii + 554);
    const auto *ii_555 = buffer.data(ii + 555);
    const auto *ii_556 = buffer.data(ii + 556);
    const auto *ii_557 = buffer.data(ii + 557);
    const auto *ii_558 = buffer.data(ii + 558);
    const auto *ii_565 = buffer.data(ii + 565);
    const auto *ii_569 = buffer.data(ii + 569);
    const auto *ii_574 = buffer.data(ii + 574);
    const auto *ii_580 = buffer.data(ii + 580);
    const auto *ii_581 = buffer.data(ii + 581);
    const auto *ii_582 = buffer.data(ii + 582);
    const auto *ii_583 = buffer.data(ii + 583);
    const auto *ii_584 = buffer.data(ii + 584);
    const auto *ii_585 = buffer.data(ii + 585);
    const auto *ii_587 = buffer.data(ii + 587);
    const auto *ii_588 = buffer.data(ii + 588);
    const auto *ii_591 = buffer.data(ii + 591);
    const auto *ii_593 = buffer.data(ii + 593);
    const auto *ii_594 = buffer.data(ii + 594);
    const auto *ii_597 = buffer.data(ii + 597);
    const auto *ii_598 = buffer.data(ii + 598);
    const auto *ii_600 = buffer.data(ii + 600);
    const auto *ii_602 = buffer.data(ii + 602);
    const auto *ii_603 = buffer.data(ii + 603);
    const auto *ii_605 = buffer.data(ii + 605);
    const auto *ii_606 = buffer.data(ii + 606);
    const auto *ii_608 = buffer.data(ii + 608);
    const auto *ii_609 = buffer.data(ii + 609);
    const auto *ii_611 = buffer.data(ii + 611);
    const auto *ii_612 = buffer.data(ii + 612);
    const auto *ii_613 = buffer.data(ii + 613);
    const auto *ii_614 = buffer.data(ii + 614);
    const auto *ii_615 = buffer.data(ii + 615);

    const auto *ik_447 = buffer.data(ik + 447);
    const auto *ik_482 = buffer.data(ik + 482);
    const auto *ik_488 = buffer.data(ik + 488);
    const auto *ik_504 = buffer.data(ik + 504);
    const auto *ik_506 = buffer.data(ik + 506);
    const auto *ik_507 = buffer.data(ik + 507);
    const auto *ik_509 = buffer.data(ik + 509);
    const auto *ik_510 = buffer.data(ik + 510);
    const auto *ik_513 = buffer.data(ik + 513);
    const auto *ik_514 = buffer.data(ik + 514);
    const auto *ik_516 = buffer.data(ik + 516);
    const auto *ik_518 = buffer.data(ik + 518);
    const auto *ik_519 = buffer.data(ik + 519);
    const auto *ik_521 = buffer.data(ik + 521);
    const auto *ik_522 = buffer.data(ik + 522);
    const auto *ik_524 = buffer.data(ik + 524);
    const auto *ik_531 = buffer.data(ik + 531);
    const auto *ik_532 = buffer.data(ik + 532);
    const auto *ik_534 = buffer.data(ik + 534);
    const auto *ik_535 = buffer.data(ik + 535);
    const auto *ik_536 = buffer.data(ik + 536);
    const auto *ik_537 = buffer.data(ik + 537);
    const auto *ik_539 = buffer.data(ik + 539);
    const auto *ik_540 = buffer.data(ik + 540);
    const auto *ik_541 = buffer.data(ik + 541);
    const auto *ik_676 = buffer.data(ik + 676);
    const auto *ik_678 = buffer.data(ik + 678);
    const auto *ik_679 = buffer.data(ik + 679);
    const auto *ik_680 = buffer.data(ik + 680);
    const auto *ik_681 = buffer.data(ik + 681);
    const auto *ik_683 = buffer.data(ik + 683);
    const auto *ik_755 = buffer.data(ik + 755);
    const auto *ik_756 = buffer.data(ik + 756);
    const auto *ik_759 = buffer.data(ik + 759);
    const auto *ik_761 = buffer.data(ik + 761);
    const auto *ik_762 = buffer.data(ik + 762);
    const auto *ik_765 = buffer.data(ik + 765);
    const auto *ik_766 = buffer.data(ik + 766);
    const auto *ik_768 = buffer.data(ik + 768);
    const auto *ik_770 = buffer.data(ik + 770);
    const auto *ik_771 = buffer.data(ik + 771);
    const auto *ik_773 = buffer.data(ik + 773);
    const auto *ik_774 = buffer.data(ik + 774);
    const auto *ik_776 = buffer.data(ik + 776);
    const auto *ik_784 = buffer.data(ik + 784);
    const auto *ik_786 = buffer.data(ik + 786);
    const auto *ik_787 = buffer.data(ik + 787);
    const auto *ik_788 = buffer.data(ik + 788);
    const auto *ik_789 = buffer.data(ik + 789);
    const auto *ik_790 = buffer.data(ik + 790);
    const auto *ik_791 = buffer.data(ik + 791);

    const auto *kh0_390 = buffer.data(kh0 + 390);
    const auto *kh0_395 = buffer.data(kh0 + 395);
    const auto *kh0_396 = buffer.data(kh0 + 396);
    const auto *kh0_420 = buffer.data(kh0 + 420);
    const auto *kh0_421 = buffer.data(kh0 + 421);
    const auto *kh0_423 = buffer.data(kh0 + 423);
    const auto *kh0_425 = buffer.data(kh0 + 425);
    const auto *kh0_426 = buffer.data(kh0 + 426);
    const auto *kh0_428 = buffer.data(kh0 + 428);
    const auto *kh0_429 = buffer.data(kh0 + 429);
    const auto *kh0_434 = buffer.data(kh0 + 434);
    const auto *kh0_435 = buffer.data(kh0 + 435);
    const auto *kh0_437 = buffer.data(kh0 + 437);
    const auto *kh0_438 = buffer.data(kh0 + 438);
    const auto *kh0_439 = buffer.data(kh0 + 439);
    const auto *kh0_440 = buffer.data(kh0 + 440);

    const auto *kh1_390 = buffer.data(kh1 + 390);
    const auto *kh1_395 = buffer.data(kh1 + 395);
    const auto *kh1_396 = buffer.data(kh1 + 396);
    const auto *kh1_420 = buffer.data(kh1 + 420);
    const auto *kh1_421 = buffer.data(kh1 + 421);
    const auto *kh1_423 = buffer.data(kh1 + 423);
    const auto *kh1_425 = buffer.data(kh1 + 425);
    const auto *kh1_426 = buffer.data(kh1 + 426);
    const auto *kh1_428 = buffer.data(kh1 + 428);
    const auto *kh1_429 = buffer.data(kh1 + 429);
    const auto *kh1_434 = buffer.data(kh1 + 434);
    const auto *kh1_435 = buffer.data(kh1 + 435);
    const auto *kh1_437 = buffer.data(kh1 + 437);
    const auto *kh1_438 = buffer.data(kh1 + 438);
    const auto *kh1_439 = buffer.data(kh1 + 439);
    const auto *kh1_440 = buffer.data(kh1 + 440);

    const auto *ki_510 = buffer.data(ki + 510);
    const auto *ki_513 = buffer.data(ki + 513);
    const auto *ki_514 = buffer.data(ki + 514);
    const auto *ki_516 = buffer.data(ki + 516);
    const auto *ki_518 = buffer.data(ki + 518);
    const auto *ki_521 = buffer.data(ki + 521);
    const auto *ki_522 = buffer.data(ki + 522);
    const auto *ki_525 = buffer.data(ki + 525);
    const auto *ki_526 = buffer.data(ki + 526);
    const auto *ki_527 = buffer.data(ki + 527);
    const auto *ki_528 = buffer.data(ki + 528);
    const auto *ki_529 = buffer.data(ki + 529);
    const auto *ki_530 = buffer.data(ki + 530);
    const auto *ki_531 = buffer.data(ki + 531);
    const auto *ki_532 = buffer.data(ki + 532);
    const auto *ki_534 = buffer.data(ki + 534);
    const auto *ki_535 = buffer.data(ki + 535);
    const auto *ki_537 = buffer.data(ki + 537);
    const auto *ki_538 = buffer.data(ki + 538);
    const auto *ki_541 = buffer.data(ki + 541);
    const auto *ki_542 = buffer.data(ki + 542);
    const auto *ki_546 = buffer.data(ki + 546);
    const auto *ki_553 = buffer.data(ki + 553);
    const auto *ki_554 = buffer.data(ki + 554);
    const auto *ki_555 = buffer.data(ki + 555);
    const auto *ki_556 = buffer.data(ki + 556);
    const auto *ki_557 = buffer.data(ki + 557);
    const auto *ki_558 = buffer.data(ki + 558);
    const auto *ki_559 = buffer.data(ki + 559);
    const auto *ki_560 = buffer.data(ki + 560);
    const auto *ki_561 = buffer.data(ki + 561);
    const auto *ki_562 = buffer.data(ki + 562);
    const auto *ki_563 = buffer.data(ki + 563);
    const auto *ki_565 = buffer.data(ki + 565);
    const auto *ki_566 = buffer.data(ki + 566);
    const auto *ki_568 = buffer.data(ki + 568);
    const auto *ki_569 = buffer.data(ki + 569);
    const auto *ki_570 = buffer.data(ki + 570);
    const auto *ki_572 = buffer.data(ki + 572);
    const auto *ki_573 = buffer.data(ki + 573);
    const auto *ki_574 = buffer.data(ki + 574);
    const auto *ki_580 = buffer.data(ki + 580);
    const auto *ki_581 = buffer.data(ki + 581);
    const auto *ki_582 = buffer.data(ki + 582);
    const auto *ki_583 = buffer.data(ki + 583);
    const auto *ki_584 = buffer.data(ki + 584);
    const auto *ki_585 = buffer.data(ki + 585);
    const auto *ki_586 = buffer.data(ki + 586);
    const auto *ki_587 = buffer.data(ki + 587);
    const auto *ki_588 = buffer.data(ki + 588);
    const auto *ki_589 = buffer.data(ki + 589);
    const auto *ki_591 = buffer.data(ki + 591);
    const auto *ki_593 = buffer.data(ki + 593);
    const auto *ki_594 = buffer.data(ki + 594);
    const auto *ki_597 = buffer.data(ki + 597);
    const auto *ki_598 = buffer.data(ki + 598);
    const auto *ki_602 = buffer.data(ki + 602);
    const auto *ki_603 = buffer.data(ki + 603);
    const auto *ki_609 = buffer.data(ki + 609);
    const auto *ki_611 = buffer.data(ki + 611);
    const auto *ki_612 = buffer.data(ki + 612);
    const auto *ki_613 = buffer.data(ki + 613);
    const auto *ki_614 = buffer.data(ki + 614);
    const auto *ki_615 = buffer.data(ki + 615);

#pragma omp simd aligned(t_659, t_660, t_661, pb_x, pb_y, pb_z, ii_342, ii_373, ii_516, \
                         kh0_390, kh1_390, ki_510, ki_513, ki_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_13 * ii_342[k]
                   + pb_z[k] * ki_510[k];

        t_660[k] = f_12 * ii_516[k]
                   + f_5 * kh0_390[k]
                   - f_6 * kh1_390[k]
                   + pb_x[k] * ki_516[k];

        t_661[k] = f_12 * ii_373[k]
                   + pb_y[k] * ki_513[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pa_y, pa_z, pb_z, hk0_267, hk0_338, hk1_267, \
                         hk1_338, ii_346, ik_447, ik_482, ki_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_17 * hk0_338[k]
                   - f_18 * hk1_338[k]
                   + pa_y[k] * ik_482[k];

        t_663[k] = f_21 * hk0_267[k]
                   - f_22 * hk1_267[k]
                   + pa_z[k] * ik_447[k];

        t_664[k] = f_13 * ii_346[k]
                   + pb_z[k] * ki_514[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, pb_x, pb_y, ii_378, ii_521, ii_522, kh0_395, \
                         kh0_396, kh1_395, kh1_396, ki_518, ki_521, \
                         ki_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_12 * ii_521[k]
                   + f_3 * kh0_395[k]
                   - f_4 * kh1_395[k]
                   + pb_x[k] * ki_521[k];

        t_666[k] = f_12 * ii_522[k]
                   + f_3 * kh0_396[k]
                   - f_4 * kh1_396[k]
                   + pb_x[k] * ki_522[k];

        t_667[k] = f_12 * ii_378[k]
                   + pb_y[k] * ki_518[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, pa_y, pb_x, hk0_344, hk1_344, ii_525, \
                         ii_526, ii_527, ik_488, ki_525, ki_526, \
                         ki_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_17 * hk0_344[k]
                   - f_18 * hk1_344[k]
                   + pa_y[k] * ik_488[k];

        t_669[k] = f_12 * ii_525[k]
                   + pb_x[k] * ki_525[k];

        t_670[k] = f_12 * ii_526[k]
                   + pb_x[k] * ki_526[k];

        t_671[k] = f_12 * ii_527[k]
                   + pb_x[k] * ki_527[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pb_x, ii_528, ii_529, ii_530, ii_531, \
                         ki_528, ki_529, ki_530, ki_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_12 * ii_528[k]
                   + pb_x[k] * ki_528[k];

        t_673[k] = f_12 * ii_529[k]
                   + pb_x[k] * ki_529[k];

        t_674[k] = f_12 * ii_530[k]
                   + pb_x[k] * ki_530[k];

        t_675[k] = f_12 * ii_531[k]
                   + pb_x[k] * ki_531[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pa_x, pb_z, hk0_676, hk0_678, hk1_676, hk1_678, \
                         ii_357, ik_676, ik_678, ki_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_17 * hk0_676[k]
                   - f_18 * hk1_676[k]
                   + pa_x[k] * ik_676[k];

        t_677[k] = f_13 * ii_357[k]
                   + pb_z[k] * ki_525[k];

        t_678[k] = f_17 * hk0_678[k]
                   - f_18 * hk1_678[k]
                   + pa_x[k] * ik_678[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, pa_x, hk0_679, hk0_680, hk0_681, hk1_679, \
                         hk1_680, hk1_681, ik_679, ik_680, ik_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_17 * hk0_679[k]
                   - f_18 * hk1_679[k]
                   + pa_x[k] * ik_679[k];

        t_680[k] = f_17 * hk0_680[k]
                   - f_18 * hk1_680[k]
                   + pa_x[k] * ik_680[k];

        t_681[k] = f_17 * hk0_681[k]
                   - f_18 * hk1_681[k]
                   + pa_x[k] * ik_681[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, pa_x, pa_y, pb_y, hk0_683, hk1_683, \
                         ii_391, ii_392, ik_504, ik_683, ki_531, \
                         ki_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_12 * ii_391[k]
                   + pb_y[k] * ki_531[k];

        t_683[k] = f_17 * hk0_683[k]
                   - f_18 * hk1_683[k]
                   + pa_x[k] * ik_683[k];

        t_684[k] = pa_y[k] * ik_504[k];

        t_685[k] = f_11 * ii_392[k]
                   + pb_y[k] * ki_532[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, t_690, pa_y, pb_y, ii_393, ii_394, \
                         ii_395, ik_506, ik_507, ik_509, ik_510, \
                         ki_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = pa_y[k] * ik_506[k];

        t_687[k] = f_12 * ii_393[k]
                   + pa_y[k] * ik_507[k];

        t_688[k] = f_11 * ii_394[k]
                   + pb_y[k] * ki_534[k];

        t_689[k] = pa_y[k] * ik_509[k];

        t_690[k] = f_13 * ii_395[k]
                   + pa_y[k] * ik_510[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, t_694, pa_y, pb_y, pb_z, ii_367, ii_397, ii_398, \
                         ik_513, ik_514, ki_535, ki_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_14 * ii_367[k]
                   + pb_z[k] * ki_535[k];

        t_692[k] = f_11 * ii_397[k]
                   + pb_y[k] * ki_537[k];

        t_693[k] = pa_y[k] * ik_513[k];

        t_694[k] = f_14 * ii_398[k]
                   + pa_y[k] * ik_514[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, pa_y, pb_y, pb_z, ii_370, ii_400, ii_401, \
                         ik_516, ik_518, ki_538, ki_541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_14 * ii_370[k]
                   + pb_z[k] * ki_538[k];

        t_696[k] = f_12 * ii_400[k]
                   + pa_y[k] * ik_516[k];

        t_697[k] = f_11 * ii_401[k]
                   + pb_y[k] * ki_541[k];

        t_698[k] = pa_y[k] * ik_518[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, pa_y, pb_z, ii_374, ii_402, ii_404, \
                         ii_405, ik_519, ik_521, ik_522, ki_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_15 * ii_402[k]
                   + pa_y[k] * ik_519[k];

        t_700[k] = f_14 * ii_374[k]
                   + pb_z[k] * ki_542[k];

        t_701[k] = f_13 * ii_404[k]
                   + pa_y[k] * ik_521[k];

        t_702[k] = f_12 * ii_405[k]
                   + pa_y[k] * ik_522[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, t_706, pa_y, pb_x, pb_y, ii_406, ii_553, ii_554, \
                         ik_524, ki_546, ki_553, ki_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_11 * ii_406[k]
                   + pb_y[k] * ki_546[k];

        t_704[k] = pa_y[k] * ik_524[k];

        t_705[k] = f_12 * ii_553[k]
                   + pb_x[k] * ki_553[k];

        t_706[k] = f_12 * ii_554[k]
                   + pb_x[k] * ki_554[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, t_710, t_711, pa_y, pb_x, ii_555, ii_556, \
                         ii_557, ii_558, ik_531, ki_555, ki_556, ki_557, \
                         ki_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = f_12 * ii_555[k]
                   + pb_x[k] * ki_555[k];

        t_708[k] = f_12 * ii_556[k]
                   + pb_x[k] * ki_556[k];

        t_709[k] = f_12 * ii_557[k]
                   + pb_x[k] * ki_557[k];

        t_710[k] = f_12 * ii_558[k]
                   + pb_x[k] * ki_558[k];

        t_711[k] = pa_y[k] * ik_531[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pa_y, pb_z, ii_385, ii_413, ii_415, \
                         ii_416, ik_532, ik_534, ik_535, ki_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_0 * ii_413[k]
                   + pa_y[k] * ik_532[k];

        t_713[k] = f_14 * ii_385[k]
                   + pb_z[k] * ki_553[k];

        t_714[k] = f_15 * ii_415[k]
                   + pa_y[k] * ik_534[k];

        t_715[k] = f_14 * ii_416[k]
                   + pa_y[k] * ik_535[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pa_y, pb_y, ii_417, ii_418, ii_419, \
                         ik_536, ik_537, ik_539, ki_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_13 * ii_417[k]
                   + pa_y[k] * ik_536[k];

        t_717[k] = f_12 * ii_418[k]
                   + pa_y[k] * ik_537[k];

        t_718[k] = f_11 * ii_419[k]
                   + pb_y[k] * ki_559[k];

        t_719[k] = pa_y[k] * ik_539[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pa_z, pb_y, pb_z, hk0_324, hk1_324, \
                         ii_392, ik_504, kh0_420, kh1_420, ki_560, \
                         ki_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_19 * hk0_324[k]
                   - f_20 * hk1_324[k]
                   + pa_z[k] * ik_504[k];

        t_721[k] = pb_y[k] * ki_560[k];

        t_722[k] = f_15 * ii_392[k]
                   + pb_z[k] * ki_560[k];

        t_723[k] = f_3 * kh0_420[k]
                   - f_4 * kh1_420[k]
                   + pb_y[k] * ki_561[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, pb_x, pb_y, pb_z, ii_395, ii_565, \
                         kh0_421, kh0_425, kh1_421, kh1_425, ki_562, ki_563, \
                         ki_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = pb_y[k] * ki_562[k];

        t_725[k] = f_12 * ii_565[k]
                   + f_9 * kh0_425[k]
                   - f_10 * kh1_425[k]
                   + pb_x[k] * ki_565[k];

        t_726[k] = f_5 * kh0_421[k]
                   - f_6 * kh1_421[k]
                   + pb_y[k] * ki_563[k];

        t_727[k] = f_15 * ii_395[k]
                   + pb_z[k] * ki_563[k];
    }

#pragma omp simd aligned(t_728, t_729, t_730, t_731, pb_x, pb_y, pb_z, ii_398, ii_569, \
                         kh0_423, kh0_429, kh1_423, kh1_429, ki_565, ki_566, \
                         ki_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_728[k] = pb_y[k] * ki_565[k];

        t_729[k] = f_12 * ii_569[k]
                   + f_7 * kh0_429[k]
                   - f_8 * kh1_429[k]
                   + pb_x[k] * ki_569[k];

        t_730[k] = f_7 * kh0_423[k]
                   - f_8 * kh1_423[k]
                   + pb_y[k] * ki_566[k];

        t_731[k] = f_15 * ii_398[k]
                   + pb_z[k] * ki_566[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pb_x, pb_y, ii_574, kh0_425, kh0_434, kh1_425, \
                         kh1_434, ki_568, ki_569, ki_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_3 * kh0_425[k]
                   - f_4 * kh1_425[k]
                   + pb_y[k] * ki_568[k];

        t_733[k] = pb_y[k] * ki_569[k];

        t_734[k] = f_12 * ii_574[k]
                   + f_5 * kh0_434[k]
                   - f_6 * kh1_434[k]
                   + pb_x[k] * ki_574[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, pb_y, pb_z, ii_402, kh0_426, kh0_428, \
                         kh0_429, kh1_426, kh1_428, kh1_429, ki_570, ki_572, \
                         ki_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_9 * kh0_426[k]
                   - f_10 * kh1_426[k]
                   + pb_y[k] * ki_570[k];

        t_736[k] = f_15 * ii_402[k]
                   + pb_z[k] * ki_570[k];

        t_737[k] = f_5 * kh0_428[k]
                   - f_6 * kh1_428[k]
                   + pb_y[k] * ki_572[k];

        t_738[k] = f_3 * kh0_429[k]
                   - f_4 * kh1_429[k]
                   + pb_y[k] * ki_573[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, pb_x, pb_y, ii_580, ii_581, ii_582, \
                         kh0_440, kh1_440, ki_574, ki_580, ki_581, \
                         ki_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = pb_y[k] * ki_574[k];

        t_740[k] = f_12 * ii_580[k]
                   + f_3 * kh0_440[k]
                   - f_4 * kh1_440[k]
                   + pb_x[k] * ki_580[k];

        t_741[k] = f_12 * ii_581[k]
                   + pb_x[k] * ki_581[k];

        t_742[k] = f_12 * ii_582[k]
                   + pb_x[k] * ki_582[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, t_747, pb_x, pb_y, ii_583, ii_584, \
                         ii_585, ii_587, ki_580, ki_583, ki_584, ki_585, \
                         ki_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = f_12 * ii_583[k]
                   + pb_x[k] * ki_583[k];

        t_744[k] = f_12 * ii_584[k]
                   + pb_x[k] * ki_584[k];

        t_745[k] = f_12 * ii_585[k]
                   + pb_x[k] * ki_585[k];

        t_746[k] = pb_y[k] * ki_580[k];

        t_747[k] = f_12 * ii_587[k]
                   + pb_x[k] * ki_587[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pb_y, pb_z, ii_413, kh0_435, kh0_437, \
                         kh0_438, kh1_435, kh1_437, kh1_438, ki_581, ki_583, \
                         ki_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_1 * kh0_435[k]
                   - f_2 * kh1_435[k]
                   + pb_y[k] * ki_581[k];

        t_749[k] = f_15 * ii_413[k]
                   + pb_z[k] * ki_581[k];

        t_750[k] = f_9 * kh0_437[k]
                   - f_10 * kh1_437[k]
                   + pb_y[k] * ki_583[k];

        t_751[k] = f_7 * kh0_438[k]
                   - f_8 * kh1_438[k]
                   + pb_y[k] * ki_584[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pa_x, pb_y, hk0_755, hk1_755, ik_755, \
                         kh0_439, kh0_440, kh1_439, kh1_440, ki_585, ki_586, \
                         ki_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_5 * kh0_439[k]
                   - f_6 * kh1_439[k]
                   + pb_y[k] * ki_585[k];

        t_753[k] = f_3 * kh0_440[k]
                   - f_4 * kh1_440[k]
                   + pb_y[k] * ki_586[k];

        t_754[k] = pb_y[k] * ki_587[k];

        t_755[k] = f_17 * hk0_755[k]
                   - f_18 * hk1_755[k]
                   + pa_x[k] * ik_755[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, t_760, pa_x, pb_y, pb_z, ii_420, ii_588, \
                         ii_591, ik_756, ik_759, ki_588, ki_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_0 * ii_588[k]
                   + pa_x[k] * ik_756[k];

        t_757[k] = f_16 * ii_420[k]
                   + pb_y[k] * ki_588[k];

        t_758[k] = pb_z[k] * ki_588[k];

        t_759[k] = f_15 * ii_591[k]
                   + pa_x[k] * ik_759[k];

        t_760[k] = pb_z[k] * ki_589[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, t_764, pa_x, pb_y, pb_z, ii_425, ii_593, ii_594, \
                         ik_761, ik_762, ki_591, ki_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_15 * ii_593[k]
                   + pa_x[k] * ik_761[k];

        t_762[k] = f_14 * ii_594[k]
                   + pa_x[k] * ik_762[k];

        t_763[k] = pb_z[k] * ki_591[k];

        t_764[k] = f_16 * ii_425[k]
                   + pb_y[k] * ki_593[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, pa_x, pb_z, ii_597, ii_598, ii_600, \
                         ik_765, ik_766, ik_768, ki_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_14 * ii_597[k]
                   + pa_x[k] * ik_765[k];

        t_766[k] = f_13 * ii_598[k]
                   + pa_x[k] * ik_766[k];

        t_767[k] = pb_z[k] * ki_594[k];

        t_768[k] = f_13 * ii_600[k]
                   + pa_x[k] * ik_768[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, t_772, pa_x, pb_y, pb_z, ii_429, ii_602, ii_603, \
                         ik_770, ik_771, ki_597, ki_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = f_16 * ii_429[k]
                   + pb_y[k] * ki_597[k];

        t_770[k] = f_13 * ii_602[k]
                   + pa_x[k] * ik_770[k];

        t_771[k] = f_12 * ii_603[k]
                   + pa_x[k] * ik_771[k];

        t_772[k] = pb_z[k] * ki_598[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, t_776, pa_x, pb_y, ii_434, ii_605, ii_606, \
                         ii_608, ik_773, ik_774, ik_776, ki_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_12 * ii_605[k]
                   + pa_x[k] * ik_773[k];

        t_774[k] = f_12 * ii_606[k]
                   + pa_x[k] * ik_774[k];

        t_775[k] = f_16 * ii_434[k]
                   + pb_y[k] * ki_602[k];

        t_776[k] = f_12 * ii_608[k]
                   + pa_x[k] * ik_776[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, t_780, t_781, pb_x, pb_z, ii_609, ii_611, \
                         ii_612, ii_613, ki_603, ki_609, ki_611, ki_612, \
                         ki_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_11 * ii_609[k]
                   + pb_x[k] * ki_609[k];

        t_778[k] = pb_z[k] * ki_603[k];

        t_779[k] = f_11 * ii_611[k]
                   + pb_x[k] * ki_611[k];

        t_780[k] = f_11 * ii_612[k]
                   + pb_x[k] * ki_612[k];

        t_781[k] = f_11 * ii_613[k]
                   + pb_x[k] * ki_613[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, t_785, t_786, pa_x, pb_x, pb_z, ii_614, ii_615, \
                         ik_784, ik_786, ki_609, ki_614, ki_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_11 * ii_614[k]
                   + pb_x[k] * ki_614[k];

        t_783[k] = f_11 * ii_615[k]
                   + pb_x[k] * ki_615[k];

        t_784[k] = pa_x[k] * ik_784[k];

        t_785[k] = pb_z[k] * ki_609[k];

        t_786[k] = pa_x[k] * ik_786[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, t_790, t_791, t_792, t_793, pa_x, pa_z, ik_540, \
                         ik_541, ik_787, ik_788, ik_789, ik_790, \
                         ik_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = pa_x[k] * ik_787[k];

        t_788[k] = pa_x[k] * ik_788[k];

        t_789[k] = pa_x[k] * ik_789[k];

        t_790[k] = pa_x[k] * ik_790[k];

        t_791[k] = pa_x[k] * ik_791[k];

        t_792[k] = pa_z[k] * ik_540[k];

        t_793[k] = pa_z[k] * ik_541[k];
    }
}

static auto
compute_prim_kk_electron_repulsion_0_piece6(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t ii,
                                            const size_t ik, const size_t ki, const size_t ncols,
                                            const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ii_420 = buffer.data(ii + 420);
    const auto *ii_423 = buffer.data(ii + 423);
    const auto *ii_426 = buffer.data(ii + 426);
    const auto *ii_430 = buffer.data(ii + 430);
    const auto *ii_448 = buffer.data(ii + 448);
    const auto *ii_450 = buffer.data(ii + 450);
    const auto *ii_451 = buffer.data(ii + 451);
    const auto *ii_453 = buffer.data(ii + 453);
    const auto *ii_454 = buffer.data(ii + 454);
    const auto *ii_457 = buffer.data(ii + 457);
    const auto *ii_458 = buffer.data(ii + 458);
    const auto *ii_462 = buffer.data(ii + 462);
    const auto *ii_476 = buffer.data(ii + 476);
    const auto *ii_478 = buffer.data(ii + 478);
    const auto *ii_479 = buffer.data(ii + 479);
    const auto *ii_481 = buffer.data(ii + 481);
    const auto *ii_482 = buffer.data(ii + 482);
    const auto *ii_485 = buffer.data(ii + 485);
    const auto *ii_486 = buffer.data(ii + 486);
    const auto *ii_490 = buffer.data(ii + 490);
    const auto *ii_504 = buffer.data(ii + 504);
    const auto *ii_506 = buffer.data(ii + 506);
    const auto *ii_507 = buffer.data(ii + 507);
    const auto *ii_509 = buffer.data(ii + 509);
    const auto *ii_510 = buffer.data(ii + 510);
    const auto *ii_513 = buffer.data(ii + 513);
    const auto *ii_514 = buffer.data(ii + 514);
    const auto *ii_518 = buffer.data(ii + 518);
    const auto *ii_532 = buffer.data(ii + 532);
    const auto *ii_534 = buffer.data(ii + 534);
    const auto *ii_537 = buffer.data(ii + 537);
    const auto *ii_541 = buffer.data(ii + 541);
    const auto *ii_546 = buffer.data(ii + 546);
    const auto *ii_560 = buffer.data(ii + 560);
    const auto *ii_562 = buffer.data(ii + 562);
    const auto *ii_621 = buffer.data(ii + 621);
    const auto *ii_625 = buffer.data(ii + 625);
    const auto *ii_628 = buffer.data(ii + 628);
    const auto *ii_630 = buffer.data(ii + 630);
    const auto *ii_633 = buffer.data(ii + 633);
    const auto *ii_634 = buffer.data(ii + 634);
    const auto *ii_636 = buffer.data(ii + 636);
    const auto *ii_638 = buffer.data(ii + 638);
    const auto *ii_639 = buffer.data(ii + 639);
    const auto *ii_640 = buffer.data(ii + 640);
    const auto *ii_641 = buffer.data(ii + 641);
    const auto *ii_642 = buffer.data(ii + 642);
    const auto *ii_643 = buffer.data(ii + 643);
    const auto *ii_644 = buffer.data(ii + 644);
    const auto *ii_647 = buffer.data(ii + 647);
    const auto *ii_649 = buffer.data(ii + 649);
    const auto *ii_650 = buffer.data(ii + 650);
    const auto *ii_653 = buffer.data(ii + 653);
    const auto *ii_654 = buffer.data(ii + 654);
    const auto *ii_656 = buffer.data(ii + 656);
    const auto *ii_658 = buffer.data(ii + 658);
    const auto *ii_659 = buffer.data(ii + 659);
    const auto *ii_661 = buffer.data(ii + 661);
    const auto *ii_662 = buffer.data(ii + 662);
    const auto *ii_664 = buffer.data(ii + 664);
    const auto *ii_665 = buffer.data(ii + 665);
    const auto *ii_666 = buffer.data(ii + 666);
    const auto *ii_667 = buffer.data(ii + 667);
    const auto *ii_668 = buffer.data(ii + 668);
    const auto *ii_669 = buffer.data(ii + 669);
    const auto *ii_670 = buffer.data(ii + 670);
    const auto *ii_671 = buffer.data(ii + 671);
    const auto *ii_672 = buffer.data(ii + 672);
    const auto *ii_675 = buffer.data(ii + 675);
    const auto *ii_677 = buffer.data(ii + 677);
    const auto *ii_678 = buffer.data(ii + 678);
    const auto *ii_681 = buffer.data(ii + 681);
    const auto *ii_682 = buffer.data(ii + 682);
    const auto *ii_684 = buffer.data(ii + 684);
    const auto *ii_686 = buffer.data(ii + 686);
    const auto *ii_687 = buffer.data(ii + 687);
    const auto *ii_689 = buffer.data(ii + 689);
    const auto *ii_690 = buffer.data(ii + 690);
    const auto *ii_692 = buffer.data(ii + 692);
    const auto *ii_693 = buffer.data(ii + 693);
    const auto *ii_694 = buffer.data(ii + 694);
    const auto *ii_695 = buffer.data(ii + 695);
    const auto *ii_696 = buffer.data(ii + 696);
    const auto *ii_697 = buffer.data(ii + 697);
    const auto *ii_698 = buffer.data(ii + 698);
    const auto *ii_699 = buffer.data(ii + 699);
    const auto *ii_700 = buffer.data(ii + 700);
    const auto *ii_703 = buffer.data(ii + 703);
    const auto *ii_705 = buffer.data(ii + 705);
    const auto *ii_706 = buffer.data(ii + 706);
    const auto *ii_709 = buffer.data(ii + 709);
    const auto *ii_710 = buffer.data(ii + 710);
    const auto *ii_712 = buffer.data(ii + 712);
    const auto *ii_714 = buffer.data(ii + 714);
    const auto *ii_715 = buffer.data(ii + 715);
    const auto *ii_717 = buffer.data(ii + 717);
    const auto *ii_718 = buffer.data(ii + 718);
    const auto *ii_720 = buffer.data(ii + 720);
    const auto *ii_721 = buffer.data(ii + 721);
    const auto *ii_722 = buffer.data(ii + 722);
    const auto *ii_723 = buffer.data(ii + 723);
    const auto *ii_724 = buffer.data(ii + 724);
    const auto *ii_725 = buffer.data(ii + 725);
    const auto *ii_726 = buffer.data(ii + 726);
    const auto *ii_727 = buffer.data(ii + 727);
    const auto *ii_731 = buffer.data(ii + 731);
    const auto *ii_734 = buffer.data(ii + 734);

    const auto *ik_543 = buffer.data(ik + 543);
    const auto *ik_546 = buffer.data(ik + 546);
    const auto *ik_550 = buffer.data(ik + 550);
    const auto *ik_555 = buffer.data(ik + 555);
    const auto *ik_561 = buffer.data(ik + 561);
    const auto *ik_720 = buffer.data(ik + 720);
    const auto *ik_722 = buffer.data(ik + 722);
    const auto *ik_725 = buffer.data(ik + 725);
    const auto *ik_797 = buffer.data(ik + 797);
    const auto *ik_801 = buffer.data(ik + 801);
    const auto *ik_804 = buffer.data(ik + 804);
    const auto *ik_806 = buffer.data(ik + 806);
    const auto *ik_809 = buffer.data(ik + 809);
    const auto *ik_810 = buffer.data(ik + 810);
    const auto *ik_812 = buffer.data(ik + 812);
    const auto *ik_820 = buffer.data(ik + 820);
    const auto *ik_821 = buffer.data(ik + 821);
    const auto *ik_822 = buffer.data(ik + 822);
    const auto *ik_823 = buffer.data(ik + 823);
    const auto *ik_824 = buffer.data(ik + 824);
    const auto *ik_825 = buffer.data(ik + 825);
    const auto *ik_826 = buffer.data(ik + 826);
    const auto *ik_827 = buffer.data(ik + 827);
    const auto *ik_828 = buffer.data(ik + 828);
    const auto *ik_831 = buffer.data(ik + 831);
    const auto *ik_833 = buffer.data(ik + 833);
    const auto *ik_834 = buffer.data(ik + 834);
    const auto *ik_837 = buffer.data(ik + 837);
    const auto *ik_838 = buffer.data(ik + 838);
    const auto *ik_840 = buffer.data(ik + 840);
    const auto *ik_842 = buffer.data(ik + 842);
    const auto *ik_843 = buffer.data(ik + 843);
    const auto *ik_845 = buffer.data(ik + 845);
    const auto *ik_846 = buffer.data(ik + 846);
    const auto *ik_848 = buffer.data(ik + 848);
    const auto *ik_856 = buffer.data(ik + 856);
    const auto *ik_857 = buffer.data(ik + 857);
    const auto *ik_858 = buffer.data(ik + 858);
    const auto *ik_859 = buffer.data(ik + 859);
    const auto *ik_860 = buffer.data(ik + 860);
    const auto *ik_861 = buffer.data(ik + 861);
    const auto *ik_862 = buffer.data(ik + 862);
    const auto *ik_863 = buffer.data(ik + 863);
    const auto *ik_864 = buffer.data(ik + 864);
    const auto *ik_867 = buffer.data(ik + 867);
    const auto *ik_869 = buffer.data(ik + 869);
    const auto *ik_870 = buffer.data(ik + 870);
    const auto *ik_873 = buffer.data(ik + 873);
    const auto *ik_874 = buffer.data(ik + 874);
    const auto *ik_876 = buffer.data(ik + 876);
    const auto *ik_878 = buffer.data(ik + 878);
    const auto *ik_879 = buffer.data(ik + 879);
    const auto *ik_881 = buffer.data(ik + 881);
    const auto *ik_882 = buffer.data(ik + 882);
    const auto *ik_884 = buffer.data(ik + 884);
    const auto *ik_892 = buffer.data(ik + 892);
    const auto *ik_893 = buffer.data(ik + 893);
    const auto *ik_894 = buffer.data(ik + 894);
    const auto *ik_895 = buffer.data(ik + 895);
    const auto *ik_896 = buffer.data(ik + 896);
    const auto *ik_897 = buffer.data(ik + 897);
    const auto *ik_898 = buffer.data(ik + 898);
    const auto *ik_899 = buffer.data(ik + 899);
    const auto *ik_900 = buffer.data(ik + 900);
    const auto *ik_903 = buffer.data(ik + 903);
    const auto *ik_905 = buffer.data(ik + 905);
    const auto *ik_906 = buffer.data(ik + 906);
    const auto *ik_909 = buffer.data(ik + 909);
    const auto *ik_910 = buffer.data(ik + 910);
    const auto *ik_912 = buffer.data(ik + 912);
    const auto *ik_914 = buffer.data(ik + 914);
    const auto *ik_915 = buffer.data(ik + 915);
    const auto *ik_917 = buffer.data(ik + 917);
    const auto *ik_918 = buffer.data(ik + 918);
    const auto *ik_920 = buffer.data(ik + 920);
    const auto *ik_928 = buffer.data(ik + 928);
    const auto *ik_929 = buffer.data(ik + 929);
    const auto *ik_930 = buffer.data(ik + 930);
    const auto *ik_931 = buffer.data(ik + 931);
    const auto *ik_932 = buffer.data(ik + 932);
    const auto *ik_933 = buffer.data(ik + 933);
    const auto *ik_934 = buffer.data(ik + 934);
    const auto *ik_935 = buffer.data(ik + 935);
    const auto *ik_939 = buffer.data(ik + 939);
    const auto *ik_942 = buffer.data(ik + 942);

    const auto *ki_616 = buffer.data(ki + 616);
    const auto *ki_618 = buffer.data(ki + 618);
    const auto *ki_619 = buffer.data(ki + 619);
    const auto *ki_621 = buffer.data(ki + 621);
    const auto *ki_622 = buffer.data(ki + 622);
    const auto *ki_625 = buffer.data(ki + 625);
    const auto *ki_626 = buffer.data(ki + 626);
    const auto *ki_630 = buffer.data(ki + 630);
    const auto *ki_638 = buffer.data(ki + 638);
    const auto *ki_639 = buffer.data(ki + 639);
    const auto *ki_640 = buffer.data(ki + 640);
    const auto *ki_641 = buffer.data(ki + 641);
    const auto *ki_642 = buffer.data(ki + 642);
    const auto *ki_643 = buffer.data(ki + 643);
    const auto *ki_644 = buffer.data(ki + 644);
    const auto *ki_646 = buffer.data(ki + 646);
    const auto *ki_647 = buffer.data(ki + 647);
    const auto *ki_649 = buffer.data(ki + 649);
    const auto *ki_650 = buffer.data(ki + 650);
    const auto *ki_653 = buffer.data(ki + 653);
    const auto *ki_654 = buffer.data(ki + 654);
    const auto *ki_658 = buffer.data(ki + 658);
    const auto *ki_665 = buffer.data(ki + 665);
    const auto *ki_666 = buffer.data(ki + 666);
    const auto *ki_667 = buffer.data(ki + 667);
    const auto *ki_668 = buffer.data(ki + 668);
    const auto *ki_669 = buffer.data(ki + 669);
    const auto *ki_670 = buffer.data(ki + 670);
    const auto *ki_671 = buffer.data(ki + 671);
    const auto *ki_672 = buffer.data(ki + 672);
    const auto *ki_674 = buffer.data(ki + 674);
    const auto *ki_675 = buffer.data(ki + 675);
    const auto *ki_677 = buffer.data(ki + 677);
    const auto *ki_678 = buffer.data(ki + 678);
    const auto *ki_681 = buffer.data(ki + 681);
    const auto *ki_682 = buffer.data(ki + 682);
    const auto *ki_686 = buffer.data(ki + 686);
    const auto *ki_693 = buffer.data(ki + 693);
    const auto *ki_694 = buffer.data(ki + 694);
    const auto *ki_695 = buffer.data(ki + 695);
    const auto *ki_696 = buffer.data(ki + 696);
    const auto *ki_697 = buffer.data(ki + 697);
    const auto *ki_698 = buffer.data(ki + 698);
    const auto *ki_699 = buffer.data(ki + 699);
    const auto *ki_700 = buffer.data(ki + 700);
    const auto *ki_702 = buffer.data(ki + 702);
    const auto *ki_703 = buffer.data(ki + 703);
    const auto *ki_705 = buffer.data(ki + 705);
    const auto *ki_706 = buffer.data(ki + 706);
    const auto *ki_709 = buffer.data(ki + 709);
    const auto *ki_710 = buffer.data(ki + 710);
    const auto *ki_714 = buffer.data(ki + 714);
    const auto *ki_721 = buffer.data(ki + 721);
    const auto *ki_722 = buffer.data(ki + 722);
    const auto *ki_723 = buffer.data(ki + 723);
    const auto *ki_724 = buffer.data(ki + 724);
    const auto *ki_725 = buffer.data(ki + 725);
    const auto *ki_726 = buffer.data(ki + 726);
    const auto *ki_727 = buffer.data(ki + 727);
    const auto *ki_728 = buffer.data(ki + 728);
    const auto *ki_730 = buffer.data(ki + 730);

#pragma omp simd aligned(t_794, t_795, t_796, t_797, pa_x, pa_z, pb_y, pb_z, ii_420, ii_450, \
                         ii_621, ik_543, ik_797, ki_616, ki_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_11 * ii_420[k]
                   + pb_z[k] * ki_616[k];

        t_795[k] = pa_z[k] * ik_543[k];

        t_796[k] = f_15 * ii_450[k]
                   + pb_y[k] * ki_618[k];

        t_797[k] = f_15 * ii_621[k]
                   + pa_x[k] * ik_797[k];
    }

#pragma omp simd aligned(t_798, t_799, t_800, t_801, pa_x, pa_z, pb_y, pb_z, ii_423, ii_453, \
                         ii_625, ik_546, ik_801, ki_619, ki_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_798[k] = pa_z[k] * ik_546[k];

        t_799[k] = f_11 * ii_423[k]
                   + pb_z[k] * ki_619[k];

        t_800[k] = f_15 * ii_453[k]
                   + pb_y[k] * ki_621[k];

        t_801[k] = f_14 * ii_625[k]
                   + pa_x[k] * ik_801[k];
    }

#pragma omp simd aligned(t_802, t_803, t_804, t_805, pa_x, pa_z, pb_y, pb_z, ii_426, ii_457, \
                         ii_628, ik_550, ik_804, ki_622, ki_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_802[k] = pa_z[k] * ik_550[k];

        t_803[k] = f_11 * ii_426[k]
                   + pb_z[k] * ki_622[k];

        t_804[k] = f_13 * ii_628[k]
                   + pa_x[k] * ik_804[k];

        t_805[k] = f_15 * ii_457[k]
                   + pb_y[k] * ki_625[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, t_809, pa_x, pa_z, pb_z, ii_430, ii_630, ii_633, \
                         ik_555, ik_806, ik_809, ki_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_13 * ii_630[k]
                   + pa_x[k] * ik_806[k];

        t_807[k] = pa_z[k] * ik_555[k];

        t_808[k] = f_11 * ii_430[k]
                   + pb_z[k] * ki_626[k];

        t_809[k] = f_12 * ii_633[k]
                   + pa_x[k] * ik_809[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, pa_x, pa_z, pb_y, ii_462, ii_634, ii_636, \
                         ik_561, ik_810, ik_812, ki_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_12 * ii_634[k]
                   + pa_x[k] * ik_810[k];

        t_811[k] = f_15 * ii_462[k]
                   + pb_y[k] * ki_630[k];

        t_812[k] = f_12 * ii_636[k]
                   + pa_x[k] * ik_812[k];

        t_813[k] = pa_z[k] * ik_561[k];
    }

#pragma omp simd aligned(t_814, t_815, t_816, t_817, t_818, pb_x, ii_638, ii_639, ii_640, \
                         ii_641, ii_642, ki_638, ki_639, ki_640, ki_641, \
                         ki_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_814[k] = f_11 * ii_638[k]
                   + pb_x[k] * ki_638[k];

        t_815[k] = f_11 * ii_639[k]
                   + pb_x[k] * ki_639[k];

        t_816[k] = f_11 * ii_640[k]
                   + pb_x[k] * ki_640[k];

        t_817[k] = f_11 * ii_641[k]
                   + pb_x[k] * ki_641[k];

        t_818[k] = f_11 * ii_642[k]
                   + pb_x[k] * ki_642[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, t_822, t_823, t_824, pa_x, pb_x, ii_643, ik_820, \
                         ik_821, ik_822, ik_823, ik_824, ki_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_11 * ii_643[k]
                   + pb_x[k] * ki_643[k];

        t_820[k] = pa_x[k] * ik_820[k];

        t_821[k] = pa_x[k] * ik_821[k];

        t_822[k] = pa_x[k] * ik_822[k];

        t_823[k] = pa_x[k] * ik_823[k];

        t_824[k] = pa_x[k] * ik_824[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, t_829, pa_x, pb_y, ii_476, ii_644, \
                         ik_825, ik_826, ik_827, ik_828, ki_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = pa_x[k] * ik_825[k];

        t_826[k] = pa_x[k] * ik_826[k];

        t_827[k] = pa_x[k] * ik_827[k];

        t_828[k] = f_0 * ii_644[k]
                   + pa_x[k] * ik_828[k];

        t_829[k] = f_14 * ii_476[k]
                   + pb_y[k] * ki_644[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, pa_x, pb_y, pb_z, ii_448, ii_478, ii_647, \
                         ii_649, ik_831, ik_833, ki_644, ki_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_12 * ii_448[k]
                   + pb_z[k] * ki_644[k];

        t_831[k] = f_15 * ii_647[k]
                   + pa_x[k] * ik_831[k];

        t_832[k] = f_14 * ii_478[k]
                   + pb_y[k] * ki_646[k];

        t_833[k] = f_15 * ii_649[k]
                   + pa_x[k] * ik_833[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, t_837, pa_x, pb_y, pb_z, ii_451, ii_481, ii_650, \
                         ii_653, ik_834, ik_837, ki_647, ki_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_14 * ii_650[k]
                   + pa_x[k] * ik_834[k];

        t_835[k] = f_12 * ii_451[k]
                   + pb_z[k] * ki_647[k];

        t_836[k] = f_14 * ii_481[k]
                   + pb_y[k] * ki_649[k];

        t_837[k] = f_14 * ii_653[k]
                   + pa_x[k] * ik_837[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, t_841, pa_x, pb_y, pb_z, ii_454, ii_485, ii_654, \
                         ii_656, ik_838, ik_840, ki_650, ki_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_13 * ii_654[k]
                   + pa_x[k] * ik_838[k];

        t_839[k] = f_12 * ii_454[k]
                   + pb_z[k] * ki_650[k];

        t_840[k] = f_13 * ii_656[k]
                   + pa_x[k] * ik_840[k];

        t_841[k] = f_14 * ii_485[k]
                   + pb_y[k] * ki_653[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, t_845, pa_x, pb_z, ii_458, ii_658, ii_659, \
                         ii_661, ik_842, ik_843, ik_845, ki_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_13 * ii_658[k]
                   + pa_x[k] * ik_842[k];

        t_843[k] = f_12 * ii_659[k]
                   + pa_x[k] * ik_843[k];

        t_844[k] = f_12 * ii_458[k]
                   + pb_z[k] * ki_654[k];

        t_845[k] = f_12 * ii_661[k]
                   + pa_x[k] * ik_845[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, t_849, pa_x, pb_x, pb_y, ii_490, ii_662, ii_664, \
                         ii_665, ik_846, ik_848, ki_658, ki_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = f_12 * ii_662[k]
                   + pa_x[k] * ik_846[k];

        t_847[k] = f_14 * ii_490[k]
                   + pb_y[k] * ki_658[k];

        t_848[k] = f_12 * ii_664[k]
                   + pa_x[k] * ik_848[k];

        t_849[k] = f_11 * ii_665[k]
                   + pb_x[k] * ki_665[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, pb_x, ii_666, ii_667, ii_668, \
                         ii_669, ii_670, ki_666, ki_667, ki_668, ki_669, \
                         ki_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = f_11 * ii_666[k]
                   + pb_x[k] * ki_666[k];

        t_851[k] = f_11 * ii_667[k]
                   + pb_x[k] * ki_667[k];

        t_852[k] = f_11 * ii_668[k]
                   + pb_x[k] * ki_668[k];

        t_853[k] = f_11 * ii_669[k]
                   + pb_x[k] * ki_669[k];

        t_854[k] = f_11 * ii_670[k]
                   + pb_x[k] * ki_670[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, t_860, pa_x, pb_x, ii_671, ik_856, \
                         ik_857, ik_858, ik_859, ik_860, ki_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = f_11 * ii_671[k]
                   + pb_x[k] * ki_671[k];

        t_856[k] = pa_x[k] * ik_856[k];

        t_857[k] = pa_x[k] * ik_857[k];

        t_858[k] = pa_x[k] * ik_858[k];

        t_859[k] = pa_x[k] * ik_859[k];

        t_860[k] = pa_x[k] * ik_860[k];
    }

#pragma omp simd aligned(t_861, t_862, t_863, t_864, t_865, pa_x, pb_y, ii_504, ii_672, \
                         ik_861, ik_862, ik_863, ik_864, ki_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_861[k] = pa_x[k] * ik_861[k];

        t_862[k] = pa_x[k] * ik_862[k];

        t_863[k] = pa_x[k] * ik_863[k];

        t_864[k] = f_0 * ii_672[k]
                   + pa_x[k] * ik_864[k];

        t_865[k] = f_13 * ii_504[k]
                   + pb_y[k] * ki_672[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, t_869, pa_x, pb_y, pb_z, ii_476, ii_506, ii_675, \
                         ii_677, ik_867, ik_869, ki_672, ki_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_13 * ii_476[k]
                   + pb_z[k] * ki_672[k];

        t_867[k] = f_15 * ii_675[k]
                   + pa_x[k] * ik_867[k];

        t_868[k] = f_13 * ii_506[k]
                   + pb_y[k] * ki_674[k];

        t_869[k] = f_15 * ii_677[k]
                   + pa_x[k] * ik_869[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, pa_x, pb_y, pb_z, ii_479, ii_509, ii_678, \
                         ii_681, ik_870, ik_873, ki_675, ki_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = f_14 * ii_678[k]
                   + pa_x[k] * ik_870[k];

        t_871[k] = f_13 * ii_479[k]
                   + pb_z[k] * ki_675[k];

        t_872[k] = f_13 * ii_509[k]
                   + pb_y[k] * ki_677[k];

        t_873[k] = f_14 * ii_681[k]
                   + pa_x[k] * ik_873[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, pa_x, pb_y, pb_z, ii_482, ii_513, ii_682, \
                         ii_684, ik_874, ik_876, ki_678, ki_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = f_13 * ii_682[k]
                   + pa_x[k] * ik_874[k];

        t_875[k] = f_13 * ii_482[k]
                   + pb_z[k] * ki_678[k];

        t_876[k] = f_13 * ii_684[k]
                   + pa_x[k] * ik_876[k];

        t_877[k] = f_13 * ii_513[k]
                   + pb_y[k] * ki_681[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, t_881, pa_x, pb_z, ii_486, ii_686, ii_687, \
                         ii_689, ik_878, ik_879, ik_881, ki_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = f_13 * ii_686[k]
                   + pa_x[k] * ik_878[k];

        t_879[k] = f_12 * ii_687[k]
                   + pa_x[k] * ik_879[k];

        t_880[k] = f_13 * ii_486[k]
                   + pb_z[k] * ki_682[k];

        t_881[k] = f_12 * ii_689[k]
                   + pa_x[k] * ik_881[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pa_x, pb_x, pb_y, ii_518, ii_690, ii_692, \
                         ii_693, ik_882, ik_884, ki_686, ki_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_12 * ii_690[k]
                   + pa_x[k] * ik_882[k];

        t_883[k] = f_13 * ii_518[k]
                   + pb_y[k] * ki_686[k];

        t_884[k] = f_12 * ii_692[k]
                   + pa_x[k] * ik_884[k];

        t_885[k] = f_11 * ii_693[k]
                   + pb_x[k] * ki_693[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, t_890, pb_x, ii_694, ii_695, ii_696, \
                         ii_697, ii_698, ki_694, ki_695, ki_696, ki_697, \
                         ki_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_11 * ii_694[k]
                   + pb_x[k] * ki_694[k];

        t_887[k] = f_11 * ii_695[k]
                   + pb_x[k] * ki_695[k];

        t_888[k] = f_11 * ii_696[k]
                   + pb_x[k] * ki_696[k];

        t_889[k] = f_11 * ii_697[k]
                   + pb_x[k] * ki_697[k];

        t_890[k] = f_11 * ii_698[k]
                   + pb_x[k] * ki_698[k];
    }

#pragma omp simd aligned(t_891, t_892, t_893, t_894, t_895, t_896, pa_x, pb_x, ii_699, ik_892, \
                         ik_893, ik_894, ik_895, ik_896, ki_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_891[k] = f_11 * ii_699[k]
                   + pb_x[k] * ki_699[k];

        t_892[k] = pa_x[k] * ik_892[k];

        t_893[k] = pa_x[k] * ik_893[k];

        t_894[k] = pa_x[k] * ik_894[k];

        t_895[k] = pa_x[k] * ik_895[k];

        t_896[k] = pa_x[k] * ik_896[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, t_900, t_901, pa_x, pb_y, ii_532, ii_700, \
                         ik_897, ik_898, ik_899, ik_900, ki_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = pa_x[k] * ik_897[k];

        t_898[k] = pa_x[k] * ik_898[k];

        t_899[k] = pa_x[k] * ik_899[k];

        t_900[k] = f_0 * ii_700[k]
                   + pa_x[k] * ik_900[k];

        t_901[k] = f_12 * ii_532[k]
                   + pb_y[k] * ki_700[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, t_905, pa_x, pb_y, pb_z, ii_504, ii_534, ii_703, \
                         ii_705, ik_903, ik_905, ki_700, ki_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_14 * ii_504[k]
                   + pb_z[k] * ki_700[k];

        t_903[k] = f_15 * ii_703[k]
                   + pa_x[k] * ik_903[k];

        t_904[k] = f_12 * ii_534[k]
                   + pb_y[k] * ki_702[k];

        t_905[k] = f_15 * ii_705[k]
                   + pa_x[k] * ik_905[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, t_909, pa_x, pb_y, pb_z, ii_507, ii_537, ii_706, \
                         ii_709, ik_906, ik_909, ki_703, ki_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_14 * ii_706[k]
                   + pa_x[k] * ik_906[k];

        t_907[k] = f_14 * ii_507[k]
                   + pb_z[k] * ki_703[k];

        t_908[k] = f_12 * ii_537[k]
                   + pb_y[k] * ki_705[k];

        t_909[k] = f_14 * ii_709[k]
                   + pa_x[k] * ik_909[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, pa_x, pb_y, pb_z, ii_510, ii_541, ii_710, \
                         ii_712, ik_910, ik_912, ki_706, ki_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = f_13 * ii_710[k]
                   + pa_x[k] * ik_910[k];

        t_911[k] = f_14 * ii_510[k]
                   + pb_z[k] * ki_706[k];

        t_912[k] = f_13 * ii_712[k]
                   + pa_x[k] * ik_912[k];

        t_913[k] = f_12 * ii_541[k]
                   + pb_y[k] * ki_709[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, t_917, pa_x, pb_z, ii_514, ii_714, ii_715, \
                         ii_717, ik_914, ik_915, ik_917, ki_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_13 * ii_714[k]
                   + pa_x[k] * ik_914[k];

        t_915[k] = f_12 * ii_715[k]
                   + pa_x[k] * ik_915[k];

        t_916[k] = f_14 * ii_514[k]
                   + pb_z[k] * ki_710[k];

        t_917[k] = f_12 * ii_717[k]
                   + pa_x[k] * ik_917[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, t_921, pa_x, pb_x, pb_y, ii_546, ii_718, ii_720, \
                         ii_721, ik_918, ik_920, ki_714, ki_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_12 * ii_718[k]
                   + pa_x[k] * ik_918[k];

        t_919[k] = f_12 * ii_546[k]
                   + pb_y[k] * ki_714[k];

        t_920[k] = f_12 * ii_720[k]
                   + pa_x[k] * ik_920[k];

        t_921[k] = f_11 * ii_721[k]
                   + pb_x[k] * ki_721[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, t_926, pb_x, ii_722, ii_723, ii_724, \
                         ii_725, ii_726, ki_722, ki_723, ki_724, ki_725, \
                         ki_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_11 * ii_722[k]
                   + pb_x[k] * ki_722[k];

        t_923[k] = f_11 * ii_723[k]
                   + pb_x[k] * ki_723[k];

        t_924[k] = f_11 * ii_724[k]
                   + pb_x[k] * ki_724[k];

        t_925[k] = f_11 * ii_725[k]
                   + pb_x[k] * ki_725[k];

        t_926[k] = f_11 * ii_726[k]
                   + pb_x[k] * ki_726[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, t_931, t_932, pa_x, pb_x, ii_727, ik_928, \
                         ik_929, ik_930, ik_931, ik_932, ki_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_11 * ii_727[k]
                   + pb_x[k] * ki_727[k];

        t_928[k] = pa_x[k] * ik_928[k];

        t_929[k] = pa_x[k] * ik_929[k];

        t_930[k] = pa_x[k] * ik_930[k];

        t_931[k] = pa_x[k] * ik_931[k];

        t_932[k] = pa_x[k] * ik_932[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, t_936, t_937, t_938, pa_x, pa_y, pb_y, ii_560, \
                         ik_720, ik_722, ik_933, ik_934, ik_935, \
                         ki_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = pa_x[k] * ik_933[k];

        t_934[k] = pa_x[k] * ik_934[k];

        t_935[k] = pa_x[k] * ik_935[k];

        t_936[k] = pa_y[k] * ik_720[k];

        t_937[k] = f_11 * ii_560[k]
                   + pb_y[k] * ki_728[k];

        t_938[k] = pa_y[k] * ik_722[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, t_942, pa_x, pa_y, pb_y, ii_562, ii_731, ii_734, \
                         ik_725, ik_939, ik_942, ki_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_15 * ii_731[k]
                   + pa_x[k] * ik_939[k];

        t_940[k] = f_11 * ii_562[k]
                   + pb_y[k] * ki_730[k];

        t_941[k] = pa_y[k] * ik_725[k];

        t_942[k] = f_14 * ii_734[k]
                   + pa_x[k] * ik_942[k];
    }
}

static auto
compute_prim_kk_electron_repulsion_0_piece7(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t ii,
                                            const size_t ik, const size_t kh0, const size_t kh1,
                                            const size_t ki, const size_t ncols,
                                            const double alpha, const double beta,
                                            const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
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
    const auto f_16 = 3.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ii_535 = buffer.data(ii + 535);
    const auto *ii_538 = buffer.data(ii + 538);
    const auto *ii_542 = buffer.data(ii + 542);
    const auto *ii_560 = buffer.data(ii + 560);
    const auto *ii_563 = buffer.data(ii + 563);
    const auto *ii_565 = buffer.data(ii + 565);
    const auto *ii_566 = buffer.data(ii + 566);
    const auto *ii_569 = buffer.data(ii + 569);
    const auto *ii_570 = buffer.data(ii + 570);
    const auto *ii_574 = buffer.data(ii + 574);
    const auto *ii_588 = buffer.data(ii + 588);
    const auto *ii_590 = buffer.data(ii + 590);
    const auto *ii_591 = buffer.data(ii + 591);
    const auto *ii_593 = buffer.data(ii + 593);
    const auto *ii_594 = buffer.data(ii + 594);
    const auto *ii_595 = buffer.data(ii + 595);
    const auto *ii_597 = buffer.data(ii + 597);
    const auto *ii_598 = buffer.data(ii + 598);
    const auto *ii_599 = buffer.data(ii + 599);
    const auto *ii_600 = buffer.data(ii + 600);
    const auto *ii_602 = buffer.data(ii + 602);
    const auto *ii_609 = buffer.data(ii + 609);
    const auto *ii_610 = buffer.data(ii + 610);
    const auto *ii_611 = buffer.data(ii + 611);
    const auto *ii_612 = buffer.data(ii + 612);
    const auto *ii_613 = buffer.data(ii + 613);
    const auto *ii_615 = buffer.data(ii + 615);
    const auto *ii_616 = buffer.data(ii + 616);
    const auto *ii_618 = buffer.data(ii + 618);
    const auto *ii_619 = buffer.data(ii + 619);
    const auto *ii_621 = buffer.data(ii + 621);
    const auto *ii_625 = buffer.data(ii + 625);
    const auto *ii_630 = buffer.data(ii + 630);
    const auto *ii_643 = buffer.data(ii + 643);
    const auto *ii_644 = buffer.data(ii + 644);
    const auto *ii_646 = buffer.data(ii + 646);
    const auto *ii_649 = buffer.data(ii + 649);
    const auto *ii_738 = buffer.data(ii + 738);
    const auto *ii_740 = buffer.data(ii + 740);
    const auto *ii_743 = buffer.data(ii + 743);
    const auto *ii_745 = buffer.data(ii + 745);
    const auto *ii_746 = buffer.data(ii + 746);
    const auto *ii_749 = buffer.data(ii + 749);
    const auto *ii_750 = buffer.data(ii + 750);
    const auto *ii_751 = buffer.data(ii + 751);
    const auto *ii_752 = buffer.data(ii + 752);
    const auto *ii_753 = buffer.data(ii + 753);
    const auto *ii_754 = buffer.data(ii + 754);
    const auto *ii_756 = buffer.data(ii + 756);
    const auto *ii_759 = buffer.data(ii + 759);
    const auto *ii_761 = buffer.data(ii + 761);
    const auto *ii_762 = buffer.data(ii + 762);
    const auto *ii_765 = buffer.data(ii + 765);
    const auto *ii_766 = buffer.data(ii + 766);
    const auto *ii_768 = buffer.data(ii + 768);
    const auto *ii_770 = buffer.data(ii + 770);
    const auto *ii_771 = buffer.data(ii + 771);
    const auto *ii_773 = buffer.data(ii + 773);
    const auto *ii_774 = buffer.data(ii + 774);
    const auto *ii_776 = buffer.data(ii + 776);
    const auto *ii_777 = buffer.data(ii + 777);
    const auto *ii_778 = buffer.data(ii + 778);
    const auto *ii_779 = buffer.data(ii + 779);
    const auto *ii_780 = buffer.data(ii + 780);
    const auto *ii_781 = buffer.data(ii + 781);
    const auto *ii_783 = buffer.data(ii + 783);

    const auto *ik_729 = buffer.data(ik + 729);
    const auto *ik_734 = buffer.data(ik + 734);
    const auto *ik_740 = buffer.data(ik + 740);
    const auto *ik_747 = buffer.data(ik + 747);
    const auto *ik_756 = buffer.data(ik + 756);
    const auto *ik_757 = buffer.data(ik + 757);
    const auto *ik_759 = buffer.data(ik + 759);
    const auto *ik_761 = buffer.data(ik + 761);
    const auto *ik_762 = buffer.data(ik + 762);
    const auto *ik_765 = buffer.data(ik + 765);
    const auto *ik_766 = buffer.data(ik + 766);
    const auto *ik_768 = buffer.data(ik + 768);
    const auto *ik_770 = buffer.data(ik + 770);
    const auto *ik_771 = buffer.data(ik + 771);
    const auto *ik_773 = buffer.data(ik + 773);
    const auto *ik_774 = buffer.data(ik + 774);
    const auto *ik_776 = buffer.data(ik + 776);
    const auto *ik_784 = buffer.data(ik + 784);
    const auto *ik_786 = buffer.data(ik + 786);
    const auto *ik_787 = buffer.data(ik + 787);
    const auto *ik_788 = buffer.data(ik + 788);
    const auto *ik_789 = buffer.data(ik + 789);
    const auto *ik_791 = buffer.data(ik + 791);
    const auto *ik_946 = buffer.data(ik + 946);
    const auto *ik_948 = buffer.data(ik + 948);
    const auto *ik_951 = buffer.data(ik + 951);
    const auto *ik_953 = buffer.data(ik + 953);
    const auto *ik_954 = buffer.data(ik + 954);
    const auto *ik_964 = buffer.data(ik + 964);
    const auto *ik_965 = buffer.data(ik + 965);
    const auto *ik_966 = buffer.data(ik + 966);
    const auto *ik_967 = buffer.data(ik + 967);
    const auto *ik_968 = buffer.data(ik + 968);
    const auto *ik_969 = buffer.data(ik + 969);
    const auto *ik_970 = buffer.data(ik + 970);
    const auto *ik_971 = buffer.data(ik + 971);
    const auto *ik_972 = buffer.data(ik + 972);
    const auto *ik_975 = buffer.data(ik + 975);
    const auto *ik_977 = buffer.data(ik + 977);
    const auto *ik_978 = buffer.data(ik + 978);
    const auto *ik_981 = buffer.data(ik + 981);
    const auto *ik_982 = buffer.data(ik + 982);
    const auto *ik_984 = buffer.data(ik + 984);
    const auto *ik_986 = buffer.data(ik + 986);
    const auto *ik_987 = buffer.data(ik + 987);
    const auto *ik_989 = buffer.data(ik + 989);
    const auto *ik_990 = buffer.data(ik + 990);
    const auto *ik_992 = buffer.data(ik + 992);
    const auto *ik_1000 = buffer.data(ik + 1000);
    const auto *ik_1001 = buffer.data(ik + 1001);
    const auto *ik_1002 = buffer.data(ik + 1002);
    const auto *ik_1003 = buffer.data(ik + 1003);
    const auto *ik_1004 = buffer.data(ik + 1004);
    const auto *ik_1005 = buffer.data(ik + 1005);
    const auto *ik_1007 = buffer.data(ik + 1007);

    const auto *kh0_588 = buffer.data(kh0 + 588);
    const auto *kh0_591 = buffer.data(kh0 + 591);
    const auto *kh0_593 = buffer.data(kh0 + 593);
    const auto *kh0_594 = buffer.data(kh0 + 594);
    const auto *kh0_597 = buffer.data(kh0 + 597);
    const auto *kh0_598 = buffer.data(kh0 + 598);
    const auto *kh0_600 = buffer.data(kh0 + 600);
    const auto *kh0_602 = buffer.data(kh0 + 602);
    const auto *kh0_603 = buffer.data(kh0 + 603);
    const auto *kh0_604 = buffer.data(kh0 + 604);
    const auto *kh0_605 = buffer.data(kh0 + 605);
    const auto *kh0_606 = buffer.data(kh0 + 606);
    const auto *kh0_608 = buffer.data(kh0 + 608);
    const auto *kh0_630 = buffer.data(kh0 + 630);
    const auto *kh0_633 = buffer.data(kh0 + 633);
    const auto *kh0_635 = buffer.data(kh0 + 635);
    const auto *kh0_636 = buffer.data(kh0 + 636);

    const auto *kh1_588 = buffer.data(kh1 + 588);
    const auto *kh1_591 = buffer.data(kh1 + 591);
    const auto *kh1_593 = buffer.data(kh1 + 593);
    const auto *kh1_594 = buffer.data(kh1 + 594);
    const auto *kh1_597 = buffer.data(kh1 + 597);
    const auto *kh1_598 = buffer.data(kh1 + 598);
    const auto *kh1_600 = buffer.data(kh1 + 600);
    const auto *kh1_602 = buffer.data(kh1 + 602);
    const auto *kh1_603 = buffer.data(kh1 + 603);
    const auto *kh1_604 = buffer.data(kh1 + 604);
    const auto *kh1_605 = buffer.data(kh1 + 605);
    const auto *kh1_606 = buffer.data(kh1 + 606);
    const auto *kh1_608 = buffer.data(kh1 + 608);
    const auto *kh1_630 = buffer.data(kh1 + 630);
    const auto *kh1_633 = buffer.data(kh1 + 633);
    const auto *kh1_635 = buffer.data(kh1 + 635);
    const auto *kh1_636 = buffer.data(kh1 + 636);

    const auto *ki_731 = buffer.data(ki + 731);
    const auto *ki_733 = buffer.data(ki + 733);
    const auto *ki_734 = buffer.data(ki + 734);
    const auto *ki_737 = buffer.data(ki + 737);
    const auto *ki_738 = buffer.data(ki + 738);
    const auto *ki_742 = buffer.data(ki + 742);
    const auto *ki_749 = buffer.data(ki + 749);
    const auto *ki_750 = buffer.data(ki + 750);
    const auto *ki_751 = buffer.data(ki + 751);
    const auto *ki_752 = buffer.data(ki + 752);
    const auto *ki_753 = buffer.data(ki + 753);
    const auto *ki_754 = buffer.data(ki + 754);
    const auto *ki_756 = buffer.data(ki + 756);
    const auto *ki_758 = buffer.data(ki + 758);
    const auto *ki_759 = buffer.data(ki + 759);
    const auto *ki_761 = buffer.data(ki + 761);
    const auto *ki_762 = buffer.data(ki + 762);
    const auto *ki_765 = buffer.data(ki + 765);
    const auto *ki_766 = buffer.data(ki + 766);
    const auto *ki_770 = buffer.data(ki + 770);
    const auto *ki_776 = buffer.data(ki + 776);
    const auto *ki_777 = buffer.data(ki + 777);
    const auto *ki_778 = buffer.data(ki + 778);
    const auto *ki_779 = buffer.data(ki + 779);
    const auto *ki_780 = buffer.data(ki + 780);
    const auto *ki_781 = buffer.data(ki + 781);
    const auto *ki_783 = buffer.data(ki + 783);
    const auto *ki_784 = buffer.data(ki + 784);
    const auto *ki_785 = buffer.data(ki + 785);
    const auto *ki_787 = buffer.data(ki + 787);
    const auto *ki_789 = buffer.data(ki + 789);
    const auto *ki_790 = buffer.data(ki + 790);
    const auto *ki_793 = buffer.data(ki + 793);
    const auto *ki_794 = buffer.data(ki + 794);
    const auto *ki_796 = buffer.data(ki + 796);
    const auto *ki_798 = buffer.data(ki + 798);
    const auto *ki_799 = buffer.data(ki + 799);
    const auto *ki_801 = buffer.data(ki + 801);
    const auto *ki_802 = buffer.data(ki + 802);
    const auto *ki_804 = buffer.data(ki + 804);
    const auto *ki_805 = buffer.data(ki + 805);
    const auto *ki_806 = buffer.data(ki + 806);
    const auto *ki_807 = buffer.data(ki + 807);
    const auto *ki_808 = buffer.data(ki + 808);
    const auto *ki_809 = buffer.data(ki + 809);
    const auto *ki_810 = buffer.data(ki + 810);
    const auto *ki_811 = buffer.data(ki + 811);
    const auto *ki_812 = buffer.data(ki + 812);
    const auto *ki_814 = buffer.data(ki + 814);
    const auto *ki_815 = buffer.data(ki + 815);
    const auto *ki_817 = buffer.data(ki + 817);
    const auto *ki_818 = buffer.data(ki + 818);
    const auto *ki_821 = buffer.data(ki + 821);
    const auto *ki_822 = buffer.data(ki + 822);
    const auto *ki_826 = buffer.data(ki + 826);
    const auto *ki_833 = buffer.data(ki + 833);
    const auto *ki_834 = buffer.data(ki + 834);
    const auto *ki_835 = buffer.data(ki + 835);
    const auto *ki_836 = buffer.data(ki + 836);
    const auto *ki_837 = buffer.data(ki + 837);
    const auto *ki_838 = buffer.data(ki + 838);
    const auto *ki_839 = buffer.data(ki + 839);
    const auto *ki_840 = buffer.data(ki + 840);
    const auto *ki_842 = buffer.data(ki + 842);
    const auto *ki_843 = buffer.data(ki + 843);
    const auto *ki_845 = buffer.data(ki + 845);
    const auto *ki_846 = buffer.data(ki + 846);

#pragma omp simd aligned(t_943, t_944, t_945, t_946, pa_x, pa_y, pb_y, pb_z, ii_535, ii_565, \
                         ii_738, ik_729, ik_946, ki_731, ki_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = f_15 * ii_535[k]
                   + pb_z[k] * ki_731[k];

        t_944[k] = f_11 * ii_565[k]
                   + pb_y[k] * ki_733[k];

        t_945[k] = pa_y[k] * ik_729[k];

        t_946[k] = f_13 * ii_738[k]
                   + pa_x[k] * ik_946[k];
    }

#pragma omp simd aligned(t_947, t_948, t_949, t_950, pa_x, pa_y, pb_y, pb_z, ii_538, ii_569, \
                         ii_740, ik_734, ik_948, ki_734, ki_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_947[k] = f_15 * ii_538[k]
                   + pb_z[k] * ki_734[k];

        t_948[k] = f_13 * ii_740[k]
                   + pa_x[k] * ik_948[k];

        t_949[k] = f_11 * ii_569[k]
                   + pb_y[k] * ki_737[k];

        t_950[k] = pa_y[k] * ik_734[k];
    }

#pragma omp simd aligned(t_951, t_952, t_953, t_954, pa_x, pb_z, ii_542, ii_743, ii_745, \
                         ii_746, ik_951, ik_953, ik_954, ki_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_951[k] = f_12 * ii_743[k]
                   + pa_x[k] * ik_951[k];

        t_952[k] = f_15 * ii_542[k]
                   + pb_z[k] * ki_738[k];

        t_953[k] = f_12 * ii_745[k]
                   + pa_x[k] * ik_953[k];

        t_954[k] = f_12 * ii_746[k]
                   + pa_x[k] * ik_954[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, t_958, pa_y, pb_x, pb_y, ii_574, ii_749, ii_750, \
                         ik_740, ki_742, ki_749, ki_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = f_11 * ii_574[k]
                   + pb_y[k] * ki_742[k];

        t_956[k] = pa_y[k] * ik_740[k];

        t_957[k] = f_11 * ii_749[k]
                   + pb_x[k] * ki_749[k];

        t_958[k] = f_11 * ii_750[k]
                   + pb_x[k] * ki_750[k];
    }

#pragma omp simd aligned(t_959, t_960, t_961, t_962, t_963, pa_y, pb_x, ii_751, ii_752, \
                         ii_753, ii_754, ik_747, ki_751, ki_752, ki_753, \
                         ki_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_959[k] = f_11 * ii_751[k]
                   + pb_x[k] * ki_751[k];

        t_960[k] = f_11 * ii_752[k]
                   + pb_x[k] * ki_752[k];

        t_961[k] = f_11 * ii_753[k]
                   + pb_x[k] * ki_753[k];

        t_962[k] = f_11 * ii_754[k]
                   + pb_x[k] * ki_754[k];

        t_963[k] = pa_y[k] * ik_747[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, t_968, t_969, t_970, pa_x, ik_964, \
                         ik_965, ik_966, ik_967, ik_968, ik_969, \
                         ik_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = pa_x[k] * ik_964[k];

        t_965[k] = pa_x[k] * ik_965[k];

        t_966[k] = pa_x[k] * ik_966[k];

        t_967[k] = pa_x[k] * ik_967[k];

        t_968[k] = pa_x[k] * ik_968[k];

        t_969[k] = pa_x[k] * ik_969[k];

        t_970[k] = pa_x[k] * ik_970[k];
    }

#pragma omp simd aligned(t_971, t_972, t_973, t_974, t_975, pa_x, pb_y, pb_z, ii_560, ii_756, \
                         ii_759, ik_971, ik_972, ik_975, ki_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_971[k] = pa_x[k] * ik_971[k];

        t_972[k] = f_0 * ii_756[k]
                   + pa_x[k] * ik_972[k];

        t_973[k] = pb_y[k] * ki_756[k];

        t_974[k] = f_16 * ii_560[k]
                   + pb_z[k] * ki_756[k];

        t_975[k] = f_15 * ii_759[k]
                   + pa_x[k] * ik_975[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, t_979, t_980, pa_x, pb_y, pb_z, ii_563, ii_761, \
                         ii_762, ik_977, ik_978, ki_758, ki_759, \
                         ki_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = pb_y[k] * ki_758[k];

        t_977[k] = f_15 * ii_761[k]
                   + pa_x[k] * ik_977[k];

        t_978[k] = f_14 * ii_762[k]
                   + pa_x[k] * ik_978[k];

        t_979[k] = f_16 * ii_563[k]
                   + pb_z[k] * ki_759[k];

        t_980[k] = pb_y[k] * ki_761[k];
    }

#pragma omp simd aligned(t_981, t_982, t_983, t_984, pa_x, pb_z, ii_566, ii_765, ii_766, \
                         ii_768, ik_981, ik_982, ik_984, ki_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_981[k] = f_14 * ii_765[k]
                   + pa_x[k] * ik_981[k];

        t_982[k] = f_13 * ii_766[k]
                   + pa_x[k] * ik_982[k];

        t_983[k] = f_16 * ii_566[k]
                   + pb_z[k] * ki_762[k];

        t_984[k] = f_13 * ii_768[k]
                   + pa_x[k] * ik_984[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, t_988, pa_x, pb_y, pb_z, ii_570, ii_770, ii_771, \
                         ik_986, ik_987, ki_765, ki_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = pb_y[k] * ki_765[k];

        t_986[k] = f_13 * ii_770[k]
                   + pa_x[k] * ik_986[k];

        t_987[k] = f_12 * ii_771[k]
                   + pa_x[k] * ik_987[k];

        t_988[k] = f_16 * ii_570[k]
                   + pb_z[k] * ki_766[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, t_992, pa_x, pb_y, ii_773, ii_774, ii_776, \
                         ik_989, ik_990, ik_992, ki_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_12 * ii_773[k]
                   + pa_x[k] * ik_989[k];

        t_990[k] = f_12 * ii_774[k]
                   + pa_x[k] * ik_990[k];

        t_991[k] = pb_y[k] * ki_770[k];

        t_992[k] = f_12 * ii_776[k]
                   + pa_x[k] * ik_992[k];
    }

#pragma omp simd aligned(t_993, t_994, t_995, t_996, t_997, pb_x, ii_777, ii_778, ii_779, \
                         ii_780, ii_781, ki_777, ki_778, ki_779, ki_780, \
                         ki_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_993[k] = f_11 * ii_777[k]
                   + pb_x[k] * ki_777[k];

        t_994[k] = f_11 * ii_778[k]
                   + pb_x[k] * ki_778[k];

        t_995[k] = f_11 * ii_779[k]
                   + pb_x[k] * ki_779[k];

        t_996[k] = f_11 * ii_780[k]
                   + pb_x[k] * ki_780[k];

        t_997[k] = f_11 * ii_781[k]
                   + pb_x[k] * ki_781[k];
    }

#pragma omp simd aligned(t_998, t_999, t_1000, t_1001, t_1002, t_1003, pa_x, pb_x, pb_y, \
                         ii_783, ik_1000, ik_1001, ik_1002, ik_1003, ki_776, \
                         ki_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = pb_y[k] * ki_776[k];

        t_999[k] = f_11 * ii_783[k]
                   + pb_x[k] * ki_783[k];

        t_1000[k] = pa_x[k] * ik_1000[k];

        t_1001[k] = pa_x[k] * ik_1001[k];

        t_1002[k] = pa_x[k] * ik_1002[k];

        t_1003[k] = pa_x[k] * ik_1003[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, t_1008, pa_x, pb_x, pb_y, ik_1004, \
                         ik_1005, ik_1007, kh0_588, kh1_588, ki_783, \
                         ki_784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = pa_x[k] * ik_1004[k];

        t_1005[k] = pa_x[k] * ik_1005[k];

        t_1006[k] = pb_y[k] * ki_783[k];

        t_1007[k] = pa_x[k] * ik_1007[k];

        t_1008[k] = f_1 * kh0_588[k]
                    - f_2 * kh1_588[k]
                    + pb_x[k] * ki_784[k];
    }

#pragma omp simd aligned(t_1009, t_1010, t_1011, t_1012, pb_x, pb_y, pb_z, ii_588, kh0_591, \
                         kh1_591, ki_784, ki_785, ki_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1009[k] = f_0 * ii_588[k]
                    + pb_y[k] * ki_784[k];

        t_1010[k] = pb_z[k] * ki_784[k];

        t_1011[k] = f_9 * kh0_591[k]
                    - f_10 * kh1_591[k]
                    + pb_x[k] * ki_787[k];

        t_1012[k] = pb_z[k] * ki_785[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, t_1016, pb_x, pb_y, pb_z, ii_593, kh0_593, \
                         kh0_594, kh1_593, kh1_594, ki_787, ki_789, \
                         ki_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = f_9 * kh0_593[k]
                    - f_10 * kh1_593[k]
                    + pb_x[k] * ki_789[k];

        t_1014[k] = f_7 * kh0_594[k]
                    - f_8 * kh1_594[k]
                    + pb_x[k] * ki_790[k];

        t_1015[k] = pb_z[k] * ki_787[k];

        t_1016[k] = f_0 * ii_593[k]
                    + pb_y[k] * ki_789[k];
    }

#pragma omp simd aligned(t_1017, t_1018, t_1019, t_1020, pb_x, pb_z, kh0_597, kh0_598, \
                         kh0_600, kh1_597, kh1_598, kh1_600, ki_790, ki_793, ki_794, \
                         ki_796 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1017[k] = f_7 * kh0_597[k]
                    - f_8 * kh1_597[k]
                    + pb_x[k] * ki_793[k];

        t_1018[k] = f_5 * kh0_598[k]
                    - f_6 * kh1_598[k]
                    + pb_x[k] * ki_794[k];

        t_1019[k] = pb_z[k] * ki_790[k];

        t_1020[k] = f_5 * kh0_600[k]
                    - f_6 * kh1_600[k]
                    + pb_x[k] * ki_796[k];
    }

#pragma omp simd aligned(t_1021, t_1022, t_1023, t_1024, pb_x, pb_y, pb_z, ii_597, kh0_602, \
                         kh0_603, kh1_602, kh1_603, ki_793, ki_794, ki_798, \
                         ki_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_0 * ii_597[k]
                    + pb_y[k] * ki_793[k];

        t_1022[k] = f_5 * kh0_602[k]
                    - f_6 * kh1_602[k]
                    + pb_x[k] * ki_798[k];

        t_1023[k] = f_3 * kh0_603[k]
                    - f_4 * kh1_603[k]
                    + pb_x[k] * ki_799[k];

        t_1024[k] = pb_z[k] * ki_794[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, pb_x, pb_y, ii_602, kh0_605, kh0_606, \
                         kh1_605, kh1_606, ki_798, ki_801, ki_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_3 * kh0_605[k]
                    - f_4 * kh1_605[k]
                    + pb_x[k] * ki_801[k];

        t_1026[k] = f_3 * kh0_606[k]
                    - f_4 * kh1_606[k]
                    + pb_x[k] * ki_802[k];

        t_1027[k] = f_0 * ii_602[k]
                    + pb_y[k] * ki_798[k];
    }

#pragma omp simd aligned(t_1028, t_1029, t_1030, t_1031, t_1032, t_1033, pb_x, kh0_608, \
                         kh1_608, ki_804, ki_805, ki_806, ki_807, ki_808, \
                         ki_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1028[k] = f_3 * kh0_608[k]
                    - f_4 * kh1_608[k]
                    + pb_x[k] * ki_804[k];

        t_1029[k] = pb_x[k] * ki_805[k];

        t_1030[k] = pb_x[k] * ki_806[k];

        t_1031[k] = pb_x[k] * ki_807[k];

        t_1032[k] = pb_x[k] * ki_808[k];

        t_1033[k] = pb_x[k] * ki_809[k];
    }

#pragma omp simd aligned(t_1034, t_1035, t_1036, t_1037, t_1038, pb_x, pb_y, pb_z, ii_609, \
                         kh0_603, kh1_603, ki_805, ki_806, ki_810, \
                         ki_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = pb_x[k] * ki_810[k];

        t_1035[k] = pb_x[k] * ki_811[k];

        t_1036[k] = f_0 * ii_609[k]
                    + f_1 * kh0_603[k]
                    - f_2 * kh1_603[k]
                    + pb_y[k] * ki_805[k];

        t_1037[k] = pb_z[k] * ki_805[k];

        t_1038[k] = f_3 * kh0_603[k]
                    - f_4 * kh1_603[k]
                    + pb_z[k] * ki_806[k];
    }

#pragma omp simd aligned(t_1039, t_1040, t_1041, pb_z, kh0_604, kh0_605, kh0_606, kh1_604, \
                         kh1_605, kh1_606, ki_807, ki_808, ki_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1039[k] = f_5 * kh0_604[k]
                    - f_6 * kh1_604[k]
                    + pb_z[k] * ki_807[k];

        t_1040[k] = f_7 * kh0_605[k]
                    - f_8 * kh1_605[k]
                    + pb_z[k] * ki_808[k];

        t_1041[k] = f_9 * kh0_606[k]
                    - f_10 * kh1_606[k]
                    + pb_z[k] * ki_809[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, t_1045, t_1046, pa_z, pb_y, pb_z, ii_588, \
                         ii_615, ik_756, ik_757, kh0_608, kh1_608, ki_811, \
                         ki_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = f_0 * ii_615[k]
                    + pb_y[k] * ki_811[k];

        t_1043[k] = f_1 * kh0_608[k]
                    - f_2 * kh1_608[k]
                    + pb_z[k] * ki_811[k];

        t_1044[k] = pa_z[k] * ik_756[k];

        t_1045[k] = pa_z[k] * ik_757[k];

        t_1046[k] = f_11 * ii_588[k]
                    + pb_z[k] * ki_812[k];
    }

#pragma omp simd aligned(t_1047, t_1048, t_1049, t_1050, t_1051, pa_z, pb_y, pb_z, ii_590, \
                         ii_591, ii_618, ik_759, ik_761, ik_762, ki_814, \
                         ki_815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1047[k] = pa_z[k] * ik_759[k];

        t_1048[k] = f_16 * ii_618[k]
                    + pb_y[k] * ki_814[k];

        t_1049[k] = f_12 * ii_590[k]
                    + pa_z[k] * ik_761[k];

        t_1050[k] = pa_z[k] * ik_762[k];

        t_1051[k] = f_11 * ii_591[k]
                    + pb_z[k] * ki_815[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, t_1055, pa_z, pb_y, pb_z, ii_593, ii_594, \
                         ii_621, ik_765, ik_766, ki_817, ki_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_16 * ii_621[k]
                    + pb_y[k] * ki_817[k];

        t_1053[k] = f_13 * ii_593[k]
                    + pa_z[k] * ik_765[k];

        t_1054[k] = pa_z[k] * ik_766[k];

        t_1055[k] = f_11 * ii_594[k]
                    + pb_z[k] * ki_818[k];
    }

#pragma omp simd aligned(t_1056, t_1057, t_1058, t_1059, pa_z, pb_y, ii_595, ii_597, ii_625, \
                         ik_768, ik_770, ik_771, ki_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1056[k] = f_12 * ii_595[k]
                    + pa_z[k] * ik_768[k];

        t_1057[k] = f_16 * ii_625[k]
                    + pb_y[k] * ki_821[k];

        t_1058[k] = f_14 * ii_597[k]
                    + pa_z[k] * ik_770[k];

        t_1059[k] = pa_z[k] * ik_771[k];
    }

#pragma omp simd aligned(t_1060, t_1061, t_1062, t_1063, pa_z, pb_y, pb_z, ii_598, ii_599, \
                         ii_600, ii_630, ik_773, ik_774, ki_822, \
                         ki_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1060[k] = f_11 * ii_598[k]
                    + pb_z[k] * ki_822[k];

        t_1061[k] = f_12 * ii_599[k]
                    + pa_z[k] * ik_773[k];

        t_1062[k] = f_13 * ii_600[k]
                    + pa_z[k] * ik_774[k];

        t_1063[k] = f_16 * ii_630[k]
                    + pb_y[k] * ki_826[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, t_1067, t_1068, t_1069, pa_z, pb_x, ii_602, \
                         ik_776, ki_833, ki_834, ki_835, ki_836, \
                         ki_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = f_15 * ii_602[k]
                    + pa_z[k] * ik_776[k];

        t_1065[k] = pb_x[k] * ki_833[k];

        t_1066[k] = pb_x[k] * ki_834[k];

        t_1067[k] = pb_x[k] * ki_835[k];

        t_1068[k] = pb_x[k] * ki_836[k];

        t_1069[k] = pb_x[k] * ki_837[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, t_1073, t_1074, pa_z, pb_x, pb_z, ii_609, \
                         ii_610, ik_784, ik_786, ki_833, ki_838, \
                         ki_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = pb_x[k] * ki_838[k];

        t_1071[k] = pb_x[k] * ki_839[k];

        t_1072[k] = pa_z[k] * ik_784[k];

        t_1073[k] = f_11 * ii_609[k]
                    + pb_z[k] * ki_833[k];

        t_1074[k] = f_12 * ii_610[k]
                    + pa_z[k] * ik_786[k];
    }

#pragma omp simd aligned(t_1075, t_1076, t_1077, t_1078, pa_z, pb_y, ii_611, ii_612, ii_613, \
                         ii_643, ik_787, ik_788, ik_789, ki_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1075[k] = f_13 * ii_611[k]
                    + pa_z[k] * ik_787[k];

        t_1076[k] = f_14 * ii_612[k]
                    + pa_z[k] * ik_788[k];

        t_1077[k] = f_15 * ii_613[k]
                    + pa_z[k] * ik_789[k];

        t_1078[k] = f_16 * ii_643[k]
                    + pb_y[k] * ki_839[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, pa_z, pb_x, pb_y, pb_z, ii_615, \
                         ii_616, ii_644, ik_791, kh0_630, kh1_630, \
                         ki_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_0 * ii_615[k]
                    + pa_z[k] * ik_791[k];

        t_1080[k] = f_1 * kh0_630[k]
                    - f_2 * kh1_630[k]
                    + pb_x[k] * ki_840[k];

        t_1081[k] = f_15 * ii_644[k]
                    + pb_y[k] * ki_840[k];

        t_1082[k] = f_12 * ii_616[k]
                    + pb_z[k] * ki_840[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, pb_x, pb_y, ii_646, kh0_633, kh0_635, \
                         kh1_633, kh1_635, ki_842, ki_843, ki_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = f_9 * kh0_633[k]
                    - f_10 * kh1_633[k]
                    + pb_x[k] * ki_843[k];

        t_1084[k] = f_15 * ii_646[k]
                    + pb_y[k] * ki_842[k];

        t_1085[k] = f_9 * kh0_635[k]
                    - f_10 * kh1_635[k]
                    + pb_x[k] * ki_845[k];
    }

#pragma omp simd aligned(t_1086, t_1087, t_1088, pb_x, pb_y, pb_z, ii_619, ii_649, kh0_636, \
                         kh1_636, ki_843, ki_845, ki_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1086[k] = f_7 * kh0_636[k]
                    - f_8 * kh1_636[k]
                    + pb_x[k] * ki_846[k];

        t_1087[k] = f_12 * ii_619[k]
                    + pb_z[k] * ki_843[k];

        t_1088[k] = f_15 * ii_649[k]
                    + pb_y[k] * ki_845[k];
    }
}

static auto
compute_prim_kk_electron_repulsion_0_piece8(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t hk0,
                                            const size_t hk1, const size_t ii, const size_t ik,
                                            const size_t kh0, const size_t kh1, const size_t ki,
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
    const auto f_17 = 0.5 / alpha;
    const auto f_18 = 0.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / alpha;
    const auto f_20 = 2.0 * beta / (alpha * p);
    const auto f_21 = 1.0 / alpha;
    const auto f_22 = beta / (alpha * p);
    const auto f_23 = 1.5 / alpha;
    const auto f_24 = 1.5 * beta / (alpha * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hk0_568 = buffer.data(hk0 + 568);
    const auto *hk0_604 = buffer.data(hk0 + 604);
    const auto *hk0_640 = buffer.data(hk0 + 640);
    const auto *hk0_647 = buffer.data(hk0 + 647);
    const auto *hk0_683 = buffer.data(hk0 + 683);
    const auto *hk0_719 = buffer.data(hk0 + 719);

    const auto *hk1_568 = buffer.data(hk1 + 568);
    const auto *hk1_604 = buffer.data(hk1 + 604);
    const auto *hk1_640 = buffer.data(hk1 + 640);
    const auto *hk1_647 = buffer.data(hk1 + 647);
    const auto *hk1_683 = buffer.data(hk1 + 683);
    const auto *hk1_719 = buffer.data(hk1 + 719);

    const auto *ii_622 = buffer.data(ii + 622);
    const auto *ii_626 = buffer.data(ii + 626);
    const auto *ii_637 = buffer.data(ii + 637);
    const auto *ii_644 = buffer.data(ii + 644);
    const auto *ii_647 = buffer.data(ii + 647);
    const auto *ii_650 = buffer.data(ii + 650);
    const auto *ii_653 = buffer.data(ii + 653);
    const auto *ii_654 = buffer.data(ii + 654);
    const auto *ii_658 = buffer.data(ii + 658);
    const auto *ii_665 = buffer.data(ii + 665);
    const auto *ii_667 = buffer.data(ii + 667);
    const auto *ii_668 = buffer.data(ii + 668);
    const auto *ii_669 = buffer.data(ii + 669);
    const auto *ii_670 = buffer.data(ii + 670);
    const auto *ii_671 = buffer.data(ii + 671);
    const auto *ii_672 = buffer.data(ii + 672);
    const auto *ii_674 = buffer.data(ii + 674);
    const auto *ii_675 = buffer.data(ii + 675);
    const auto *ii_677 = buffer.data(ii + 677);
    const auto *ii_678 = buffer.data(ii + 678);
    const auto *ii_681 = buffer.data(ii + 681);
    const auto *ii_682 = buffer.data(ii + 682);
    const auto *ii_686 = buffer.data(ii + 686);
    const auto *ii_693 = buffer.data(ii + 693);
    const auto *ii_695 = buffer.data(ii + 695);
    const auto *ii_696 = buffer.data(ii + 696);
    const auto *ii_697 = buffer.data(ii + 697);
    const auto *ii_698 = buffer.data(ii + 698);
    const auto *ii_699 = buffer.data(ii + 699);
    const auto *ii_700 = buffer.data(ii + 700);
    const auto *ii_702 = buffer.data(ii + 702);
    const auto *ii_703 = buffer.data(ii + 703);
    const auto *ii_705 = buffer.data(ii + 705);
    const auto *ii_706 = buffer.data(ii + 706);
    const auto *ii_709 = buffer.data(ii + 709);
    const auto *ii_710 = buffer.data(ii + 710);
    const auto *ii_714 = buffer.data(ii + 714);
    const auto *ii_723 = buffer.data(ii + 723);
    const auto *ii_724 = buffer.data(ii + 724);
    const auto *ii_725 = buffer.data(ii + 725);
    const auto *ii_726 = buffer.data(ii + 726);
    const auto *ii_727 = buffer.data(ii + 727);
    const auto *ii_728 = buffer.data(ii + 728);
    const auto *ii_730 = buffer.data(ii + 730);
    const auto *ii_733 = buffer.data(ii + 733);
    const auto *ii_737 = buffer.data(ii + 737);
    const auto *ii_742 = buffer.data(ii + 742);

    const auto *ik_820 = buffer.data(ik + 820);
    const auto *ik_856 = buffer.data(ik + 856);
    const auto *ik_863 = buffer.data(ik + 863);
    const auto *ik_892 = buffer.data(ik + 892);
    const auto *ik_899 = buffer.data(ik + 899);
    const auto *ik_935 = buffer.data(ik + 935);

    const auto *kh0_639 = buffer.data(kh0 + 639);
    const auto *kh0_640 = buffer.data(kh0 + 640);
    const auto *kh0_642 = buffer.data(kh0 + 642);
    const auto *kh0_644 = buffer.data(kh0 + 644);
    const auto *kh0_645 = buffer.data(kh0 + 645);
    const auto *kh0_647 = buffer.data(kh0 + 647);
    const auto *kh0_648 = buffer.data(kh0 + 648);
    const auto *kh0_649 = buffer.data(kh0 + 649);
    const auto *kh0_650 = buffer.data(kh0 + 650);
    const auto *kh0_651 = buffer.data(kh0 + 651);
    const auto *kh0_654 = buffer.data(kh0 + 654);
    const auto *kh0_656 = buffer.data(kh0 + 656);
    const auto *kh0_657 = buffer.data(kh0 + 657);
    const auto *kh0_660 = buffer.data(kh0 + 660);
    const auto *kh0_661 = buffer.data(kh0 + 661);
    const auto *kh0_663 = buffer.data(kh0 + 663);
    const auto *kh0_665 = buffer.data(kh0 + 665);
    const auto *kh0_666 = buffer.data(kh0 + 666);
    const auto *kh0_668 = buffer.data(kh0 + 668);
    const auto *kh0_669 = buffer.data(kh0 + 669);
    const auto *kh0_670 = buffer.data(kh0 + 670);
    const auto *kh0_671 = buffer.data(kh0 + 671);
    const auto *kh0_672 = buffer.data(kh0 + 672);
    const auto *kh0_675 = buffer.data(kh0 + 675);
    const auto *kh0_677 = buffer.data(kh0 + 677);
    const auto *kh0_678 = buffer.data(kh0 + 678);
    const auto *kh0_681 = buffer.data(kh0 + 681);
    const auto *kh0_682 = buffer.data(kh0 + 682);
    const auto *kh0_684 = buffer.data(kh0 + 684);
    const auto *kh0_686 = buffer.data(kh0 + 686);
    const auto *kh0_687 = buffer.data(kh0 + 687);
    const auto *kh0_689 = buffer.data(kh0 + 689);
    const auto *kh0_690 = buffer.data(kh0 + 690);
    const auto *kh0_691 = buffer.data(kh0 + 691);
    const auto *kh0_692 = buffer.data(kh0 + 692);
    const auto *kh0_693 = buffer.data(kh0 + 693);
    const auto *kh0_696 = buffer.data(kh0 + 696);
    const auto *kh0_698 = buffer.data(kh0 + 698);
    const auto *kh0_699 = buffer.data(kh0 + 699);
    const auto *kh0_702 = buffer.data(kh0 + 702);
    const auto *kh0_703 = buffer.data(kh0 + 703);
    const auto *kh0_705 = buffer.data(kh0 + 705);
    const auto *kh0_707 = buffer.data(kh0 + 707);
    const auto *kh0_708 = buffer.data(kh0 + 708);
    const auto *kh0_710 = buffer.data(kh0 + 710);
    const auto *kh0_711 = buffer.data(kh0 + 711);
    const auto *kh0_713 = buffer.data(kh0 + 713);

    const auto *kh1_639 = buffer.data(kh1 + 639);
    const auto *kh1_640 = buffer.data(kh1 + 640);
    const auto *kh1_642 = buffer.data(kh1 + 642);
    const auto *kh1_644 = buffer.data(kh1 + 644);
    const auto *kh1_645 = buffer.data(kh1 + 645);
    const auto *kh1_647 = buffer.data(kh1 + 647);
    const auto *kh1_648 = buffer.data(kh1 + 648);
    const auto *kh1_649 = buffer.data(kh1 + 649);
    const auto *kh1_650 = buffer.data(kh1 + 650);
    const auto *kh1_651 = buffer.data(kh1 + 651);
    const auto *kh1_654 = buffer.data(kh1 + 654);
    const auto *kh1_656 = buffer.data(kh1 + 656);
    const auto *kh1_657 = buffer.data(kh1 + 657);
    const auto *kh1_660 = buffer.data(kh1 + 660);
    const auto *kh1_661 = buffer.data(kh1 + 661);
    const auto *kh1_663 = buffer.data(kh1 + 663);
    const auto *kh1_665 = buffer.data(kh1 + 665);
    const auto *kh1_666 = buffer.data(kh1 + 666);
    const auto *kh1_668 = buffer.data(kh1 + 668);
    const auto *kh1_669 = buffer.data(kh1 + 669);
    const auto *kh1_670 = buffer.data(kh1 + 670);
    const auto *kh1_671 = buffer.data(kh1 + 671);
    const auto *kh1_672 = buffer.data(kh1 + 672);
    const auto *kh1_675 = buffer.data(kh1 + 675);
    const auto *kh1_677 = buffer.data(kh1 + 677);
    const auto *kh1_678 = buffer.data(kh1 + 678);
    const auto *kh1_681 = buffer.data(kh1 + 681);
    const auto *kh1_682 = buffer.data(kh1 + 682);
    const auto *kh1_684 = buffer.data(kh1 + 684);
    const auto *kh1_686 = buffer.data(kh1 + 686);
    const auto *kh1_687 = buffer.data(kh1 + 687);
    const auto *kh1_689 = buffer.data(kh1 + 689);
    const auto *kh1_690 = buffer.data(kh1 + 690);
    const auto *kh1_691 = buffer.data(kh1 + 691);
    const auto *kh1_692 = buffer.data(kh1 + 692);
    const auto *kh1_693 = buffer.data(kh1 + 693);
    const auto *kh1_696 = buffer.data(kh1 + 696);
    const auto *kh1_698 = buffer.data(kh1 + 698);
    const auto *kh1_699 = buffer.data(kh1 + 699);
    const auto *kh1_702 = buffer.data(kh1 + 702);
    const auto *kh1_703 = buffer.data(kh1 + 703);
    const auto *kh1_705 = buffer.data(kh1 + 705);
    const auto *kh1_707 = buffer.data(kh1 + 707);
    const auto *kh1_708 = buffer.data(kh1 + 708);
    const auto *kh1_710 = buffer.data(kh1 + 710);
    const auto *kh1_711 = buffer.data(kh1 + 711);
    const auto *kh1_713 = buffer.data(kh1 + 713);

    const auto *ki_846 = buffer.data(ki + 846);
    const auto *ki_849 = buffer.data(ki + 849);
    const auto *ki_850 = buffer.data(ki + 850);
    const auto *ki_852 = buffer.data(ki + 852);
    const auto *ki_854 = buffer.data(ki + 854);
    const auto *ki_855 = buffer.data(ki + 855);
    const auto *ki_857 = buffer.data(ki + 857);
    const auto *ki_858 = buffer.data(ki + 858);
    const auto *ki_860 = buffer.data(ki + 860);
    const auto *ki_861 = buffer.data(ki + 861);
    const auto *ki_862 = buffer.data(ki + 862);
    const auto *ki_863 = buffer.data(ki + 863);
    const auto *ki_864 = buffer.data(ki + 864);
    const auto *ki_865 = buffer.data(ki + 865);
    const auto *ki_866 = buffer.data(ki + 866);
    const auto *ki_867 = buffer.data(ki + 867);
    const auto *ki_868 = buffer.data(ki + 868);
    const auto *ki_870 = buffer.data(ki + 870);
    const auto *ki_871 = buffer.data(ki + 871);
    const auto *ki_873 = buffer.data(ki + 873);
    const auto *ki_874 = buffer.data(ki + 874);
    const auto *ki_877 = buffer.data(ki + 877);
    const auto *ki_878 = buffer.data(ki + 878);
    const auto *ki_880 = buffer.data(ki + 880);
    const auto *ki_882 = buffer.data(ki + 882);
    const auto *ki_883 = buffer.data(ki + 883);
    const auto *ki_885 = buffer.data(ki + 885);
    const auto *ki_886 = buffer.data(ki + 886);
    const auto *ki_888 = buffer.data(ki + 888);
    const auto *ki_889 = buffer.data(ki + 889);
    const auto *ki_890 = buffer.data(ki + 890);
    const auto *ki_891 = buffer.data(ki + 891);
    const auto *ki_892 = buffer.data(ki + 892);
    const auto *ki_893 = buffer.data(ki + 893);
    const auto *ki_894 = buffer.data(ki + 894);
    const auto *ki_895 = buffer.data(ki + 895);
    const auto *ki_896 = buffer.data(ki + 896);
    const auto *ki_898 = buffer.data(ki + 898);
    const auto *ki_899 = buffer.data(ki + 899);
    const auto *ki_901 = buffer.data(ki + 901);
    const auto *ki_902 = buffer.data(ki + 902);
    const auto *ki_905 = buffer.data(ki + 905);
    const auto *ki_906 = buffer.data(ki + 906);
    const auto *ki_908 = buffer.data(ki + 908);
    const auto *ki_910 = buffer.data(ki + 910);
    const auto *ki_911 = buffer.data(ki + 911);
    const auto *ki_913 = buffer.data(ki + 913);
    const auto *ki_914 = buffer.data(ki + 914);
    const auto *ki_916 = buffer.data(ki + 916);
    const auto *ki_917 = buffer.data(ki + 917);
    const auto *ki_918 = buffer.data(ki + 918);
    const auto *ki_919 = buffer.data(ki + 919);
    const auto *ki_920 = buffer.data(ki + 920);
    const auto *ki_921 = buffer.data(ki + 921);
    const auto *ki_922 = buffer.data(ki + 922);
    const auto *ki_923 = buffer.data(ki + 923);
    const auto *ki_924 = buffer.data(ki + 924);
    const auto *ki_926 = buffer.data(ki + 926);
    const auto *ki_927 = buffer.data(ki + 927);
    const auto *ki_929 = buffer.data(ki + 929);
    const auto *ki_930 = buffer.data(ki + 930);
    const auto *ki_933 = buffer.data(ki + 933);
    const auto *ki_934 = buffer.data(ki + 934);
    const auto *ki_936 = buffer.data(ki + 936);
    const auto *ki_938 = buffer.data(ki + 938);
    const auto *ki_939 = buffer.data(ki + 939);
    const auto *ki_941 = buffer.data(ki + 941);
    const auto *ki_942 = buffer.data(ki + 942);
    const auto *ki_944 = buffer.data(ki + 944);
    const auto *ki_945 = buffer.data(ki + 945);

#pragma omp simd aligned(t_1089, t_1090, t_1091, pb_x, pb_z, ii_622, kh0_639, kh0_640, \
                         kh1_639, kh1_640, ki_846, ki_849, ki_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1089[k] = f_7 * kh0_639[k]
                    - f_8 * kh1_639[k]
                    + pb_x[k] * ki_849[k];

        t_1090[k] = f_5 * kh0_640[k]
                    - f_6 * kh1_640[k]
                    + pb_x[k] * ki_850[k];

        t_1091[k] = f_12 * ii_622[k]
                    + pb_z[k] * ki_846[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, pb_x, pb_y, ii_653, kh0_642, kh0_644, \
                         kh1_642, kh1_644, ki_849, ki_852, ki_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = f_5 * kh0_642[k]
                    - f_6 * kh1_642[k]
                    + pb_x[k] * ki_852[k];

        t_1093[k] = f_15 * ii_653[k]
                    + pb_y[k] * ki_849[k];

        t_1094[k] = f_5 * kh0_644[k]
                    - f_6 * kh1_644[k]
                    + pb_x[k] * ki_854[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, pb_x, pb_z, ii_626, kh0_645, kh0_647, \
                         kh1_645, kh1_647, ki_850, ki_855, ki_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = f_3 * kh0_645[k]
                    - f_4 * kh1_645[k]
                    + pb_x[k] * ki_855[k];

        t_1096[k] = f_12 * ii_626[k]
                    + pb_z[k] * ki_850[k];

        t_1097[k] = f_3 * kh0_647[k]
                    - f_4 * kh1_647[k]
                    + pb_x[k] * ki_857[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, t_1101, pb_x, pb_y, ii_658, kh0_648, kh0_650, \
                         kh1_648, kh1_650, ki_854, ki_858, ki_860, \
                         ki_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = f_3 * kh0_648[k]
                    - f_4 * kh1_648[k]
                    + pb_x[k] * ki_858[k];

        t_1099[k] = f_15 * ii_658[k]
                    + pb_y[k] * ki_854[k];

        t_1100[k] = f_3 * kh0_650[k]
                    - f_4 * kh1_650[k]
                    + pb_x[k] * ki_860[k];

        t_1101[k] = pb_x[k] * ki_861[k];
    }

#pragma omp simd aligned(t_1102, t_1103, t_1104, t_1105, t_1106, t_1107, pb_x, ki_862, ki_863, \
                         ki_864, ki_865, ki_866, ki_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1102[k] = pb_x[k] * ki_862[k];

        t_1103[k] = pb_x[k] * ki_863[k];

        t_1104[k] = pb_x[k] * ki_864[k];

        t_1105[k] = pb_x[k] * ki_865[k];

        t_1106[k] = pb_x[k] * ki_866[k];

        t_1107[k] = pb_x[k] * ki_867[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, pa_z, pb_y, pb_z, hk0_568, hk1_568, ii_637, \
                         ii_667, ik_820, kh0_647, kh1_647, ki_861, \
                         ki_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_17 * hk0_568[k]
                    - f_18 * hk1_568[k]
                    + pa_z[k] * ik_820[k];

        t_1109[k] = f_12 * ii_637[k]
                    + pb_z[k] * ki_861[k];

        t_1110[k] = f_15 * ii_667[k]
                    + f_9 * kh0_647[k]
                    - f_10 * kh1_647[k]
                    + pb_y[k] * ki_863[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, pb_y, ii_668, ii_669, ii_670, kh0_648, \
                         kh0_649, kh0_650, kh1_648, kh1_649, kh1_650, ki_864, ki_865, \
                         ki_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_15 * ii_668[k]
                    + f_7 * kh0_648[k]
                    - f_8 * kh1_648[k]
                    + pb_y[k] * ki_864[k];

        t_1112[k] = f_15 * ii_669[k]
                    + f_5 * kh0_649[k]
                    - f_6 * kh1_649[k]
                    + pb_y[k] * ki_865[k];

        t_1113[k] = f_15 * ii_670[k]
                    + f_3 * kh0_650[k]
                    - f_4 * kh1_650[k]
                    + pb_y[k] * ki_866[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, t_1117, pa_y, pb_x, pb_y, hk0_647, hk1_647, \
                         ii_671, ii_672, ik_863, kh0_651, kh1_651, ki_867, \
                         ki_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_15 * ii_671[k]
                    + pb_y[k] * ki_867[k];

        t_1115[k] = f_19 * hk0_647[k]
                    - f_20 * hk1_647[k]
                    + pa_y[k] * ik_863[k];

        t_1116[k] = f_1 * kh0_651[k]
                    - f_2 * kh1_651[k]
                    + pb_x[k] * ki_868[k];

        t_1117[k] = f_14 * ii_672[k]
                    + pb_y[k] * ki_868[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, pb_x, pb_y, pb_z, ii_644, ii_674, kh0_654, \
                         kh1_654, ki_868, ki_870, ki_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = f_13 * ii_644[k]
                    + pb_z[k] * ki_868[k];

        t_1119[k] = f_9 * kh0_654[k]
                    - f_10 * kh1_654[k]
                    + pb_x[k] * ki_871[k];

        t_1120[k] = f_14 * ii_674[k]
                    + pb_y[k] * ki_870[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, t_1124, pb_x, pb_y, pb_z, ii_647, ii_677, \
                         kh0_656, kh0_657, kh1_656, kh1_657, ki_871, ki_873, \
                         ki_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_9 * kh0_656[k]
                    - f_10 * kh1_656[k]
                    + pb_x[k] * ki_873[k];

        t_1122[k] = f_7 * kh0_657[k]
                    - f_8 * kh1_657[k]
                    + pb_x[k] * ki_874[k];

        t_1123[k] = f_13 * ii_647[k]
                    + pb_z[k] * ki_871[k];

        t_1124[k] = f_14 * ii_677[k]
                    + pb_y[k] * ki_873[k];
    }

#pragma omp simd aligned(t_1125, t_1126, t_1127, pb_x, pb_z, ii_650, kh0_660, kh0_661, \
                         kh1_660, kh1_661, ki_874, ki_877, ki_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1125[k] = f_7 * kh0_660[k]
                    - f_8 * kh1_660[k]
                    + pb_x[k] * ki_877[k];

        t_1126[k] = f_5 * kh0_661[k]
                    - f_6 * kh1_661[k]
                    + pb_x[k] * ki_878[k];

        t_1127[k] = f_13 * ii_650[k]
                    + pb_z[k] * ki_874[k];
    }

#pragma omp simd aligned(t_1128, t_1129, t_1130, pb_x, pb_y, ii_681, kh0_663, kh0_665, \
                         kh1_663, kh1_665, ki_877, ki_880, ki_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1128[k] = f_5 * kh0_663[k]
                    - f_6 * kh1_663[k]
                    + pb_x[k] * ki_880[k];

        t_1129[k] = f_14 * ii_681[k]
                    + pb_y[k] * ki_877[k];

        t_1130[k] = f_5 * kh0_665[k]
                    - f_6 * kh1_665[k]
                    + pb_x[k] * ki_882[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, pb_x, pb_z, ii_654, kh0_666, kh0_668, \
                         kh1_666, kh1_668, ki_878, ki_883, ki_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = f_3 * kh0_666[k]
                    - f_4 * kh1_666[k]
                    + pb_x[k] * ki_883[k];

        t_1132[k] = f_13 * ii_654[k]
                    + pb_z[k] * ki_878[k];

        t_1133[k] = f_3 * kh0_668[k]
                    - f_4 * kh1_668[k]
                    + pb_x[k] * ki_885[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, pb_x, pb_y, ii_686, kh0_669, kh0_671, \
                         kh1_669, kh1_671, ki_882, ki_886, ki_888, \
                         ki_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = f_3 * kh0_669[k]
                    - f_4 * kh1_669[k]
                    + pb_x[k] * ki_886[k];

        t_1135[k] = f_14 * ii_686[k]
                    + pb_y[k] * ki_882[k];

        t_1136[k] = f_3 * kh0_671[k]
                    - f_4 * kh1_671[k]
                    + pb_x[k] * ki_888[k];

        t_1137[k] = pb_x[k] * ki_889[k];
    }

#pragma omp simd aligned(t_1138, t_1139, t_1140, t_1141, t_1142, t_1143, pb_x, ki_890, ki_891, \
                         ki_892, ki_893, ki_894, ki_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1138[k] = pb_x[k] * ki_890[k];

        t_1139[k] = pb_x[k] * ki_891[k];

        t_1140[k] = pb_x[k] * ki_892[k];

        t_1141[k] = pb_x[k] * ki_893[k];

        t_1142[k] = pb_x[k] * ki_894[k];

        t_1143[k] = pb_x[k] * ki_895[k];
    }

#pragma omp simd aligned(t_1144, t_1145, t_1146, pa_z, pb_y, pb_z, hk0_604, hk1_604, ii_665, \
                         ii_695, ik_856, kh0_668, kh1_668, ki_889, \
                         ki_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1144[k] = f_21 * hk0_604[k]
                    - f_22 * hk1_604[k]
                    + pa_z[k] * ik_856[k];

        t_1145[k] = f_13 * ii_665[k]
                    + pb_z[k] * ki_889[k];

        t_1146[k] = f_14 * ii_695[k]
                    + f_9 * kh0_668[k]
                    - f_10 * kh1_668[k]
                    + pb_y[k] * ki_891[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, pb_y, ii_696, ii_697, ii_698, kh0_669, \
                         kh0_670, kh0_671, kh1_669, kh1_670, kh1_671, ki_892, ki_893, \
                         ki_894 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = f_14 * ii_696[k]
                    + f_7 * kh0_669[k]
                    - f_8 * kh1_669[k]
                    + pb_y[k] * ki_892[k];

        t_1148[k] = f_14 * ii_697[k]
                    + f_5 * kh0_670[k]
                    - f_6 * kh1_670[k]
                    + pb_y[k] * ki_893[k];

        t_1149[k] = f_14 * ii_698[k]
                    + f_3 * kh0_671[k]
                    - f_4 * kh1_671[k]
                    + pb_y[k] * ki_894[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, t_1153, pa_y, pb_x, pb_y, hk0_683, hk1_683, \
                         ii_699, ii_700, ik_899, kh0_672, kh1_672, ki_895, \
                         ki_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = f_14 * ii_699[k]
                    + pb_y[k] * ki_895[k];

        t_1151[k] = f_23 * hk0_683[k]
                    - f_24 * hk1_683[k]
                    + pa_y[k] * ik_899[k];

        t_1152[k] = f_1 * kh0_672[k]
                    - f_2 * kh1_672[k]
                    + pb_x[k] * ki_896[k];

        t_1153[k] = f_13 * ii_700[k]
                    + pb_y[k] * ki_896[k];
    }

#pragma omp simd aligned(t_1154, t_1155, t_1156, pb_x, pb_y, pb_z, ii_672, ii_702, kh0_675, \
                         kh1_675, ki_896, ki_898, ki_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1154[k] = f_14 * ii_672[k]
                    + pb_z[k] * ki_896[k];

        t_1155[k] = f_9 * kh0_675[k]
                    - f_10 * kh1_675[k]
                    + pb_x[k] * ki_899[k];

        t_1156[k] = f_13 * ii_702[k]
                    + pb_y[k] * ki_898[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, t_1160, pb_x, pb_y, pb_z, ii_675, ii_705, \
                         kh0_677, kh0_678, kh1_677, kh1_678, ki_899, ki_901, \
                         ki_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_9 * kh0_677[k]
                    - f_10 * kh1_677[k]
                    + pb_x[k] * ki_901[k];

        t_1158[k] = f_7 * kh0_678[k]
                    - f_8 * kh1_678[k]
                    + pb_x[k] * ki_902[k];

        t_1159[k] = f_14 * ii_675[k]
                    + pb_z[k] * ki_899[k];

        t_1160[k] = f_13 * ii_705[k]
                    + pb_y[k] * ki_901[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, pb_x, pb_z, ii_678, kh0_681, kh0_682, \
                         kh1_681, kh1_682, ki_902, ki_905, ki_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = f_7 * kh0_681[k]
                    - f_8 * kh1_681[k]
                    + pb_x[k] * ki_905[k];

        t_1162[k] = f_5 * kh0_682[k]
                    - f_6 * kh1_682[k]
                    + pb_x[k] * ki_906[k];

        t_1163[k] = f_14 * ii_678[k]
                    + pb_z[k] * ki_902[k];
    }

#pragma omp simd aligned(t_1164, t_1165, t_1166, pb_x, pb_y, ii_709, kh0_684, kh0_686, \
                         kh1_684, kh1_686, ki_905, ki_908, ki_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1164[k] = f_5 * kh0_684[k]
                    - f_6 * kh1_684[k]
                    + pb_x[k] * ki_908[k];

        t_1165[k] = f_13 * ii_709[k]
                    + pb_y[k] * ki_905[k];

        t_1166[k] = f_5 * kh0_686[k]
                    - f_6 * kh1_686[k]
                    + pb_x[k] * ki_910[k];
    }

#pragma omp simd aligned(t_1167, t_1168, t_1169, pb_x, pb_z, ii_682, kh0_687, kh0_689, \
                         kh1_687, kh1_689, ki_906, ki_911, ki_913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1167[k] = f_3 * kh0_687[k]
                    - f_4 * kh1_687[k]
                    + pb_x[k] * ki_911[k];

        t_1168[k] = f_14 * ii_682[k]
                    + pb_z[k] * ki_906[k];

        t_1169[k] = f_3 * kh0_689[k]
                    - f_4 * kh1_689[k]
                    + pb_x[k] * ki_913[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, t_1173, pb_x, pb_y, ii_714, kh0_690, kh0_692, \
                         kh1_690, kh1_692, ki_910, ki_914, ki_916, \
                         ki_917 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = f_3 * kh0_690[k]
                    - f_4 * kh1_690[k]
                    + pb_x[k] * ki_914[k];

        t_1171[k] = f_13 * ii_714[k]
                    + pb_y[k] * ki_910[k];

        t_1172[k] = f_3 * kh0_692[k]
                    - f_4 * kh1_692[k]
                    + pb_x[k] * ki_916[k];

        t_1173[k] = pb_x[k] * ki_917[k];
    }

#pragma omp simd aligned(t_1174, t_1175, t_1176, t_1177, t_1178, t_1179, pb_x, ki_918, ki_919, \
                         ki_920, ki_921, ki_922, ki_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1174[k] = pb_x[k] * ki_918[k];

        t_1175[k] = pb_x[k] * ki_919[k];

        t_1176[k] = pb_x[k] * ki_920[k];

        t_1177[k] = pb_x[k] * ki_921[k];

        t_1178[k] = pb_x[k] * ki_922[k];

        t_1179[k] = pb_x[k] * ki_923[k];
    }

#pragma omp simd aligned(t_1180, t_1181, t_1182, pa_z, pb_y, pb_z, hk0_640, hk1_640, ii_693, \
                         ii_723, ik_892, kh0_689, kh1_689, ki_917, \
                         ki_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1180[k] = f_23 * hk0_640[k]
                    - f_24 * hk1_640[k]
                    + pa_z[k] * ik_892[k];

        t_1181[k] = f_14 * ii_693[k]
                    + pb_z[k] * ki_917[k];

        t_1182[k] = f_13 * ii_723[k]
                    + f_9 * kh0_689[k]
                    - f_10 * kh1_689[k]
                    + pb_y[k] * ki_919[k];
    }

#pragma omp simd aligned(t_1183, t_1184, t_1185, pb_y, ii_724, ii_725, ii_726, kh0_690, \
                         kh0_691, kh0_692, kh1_690, kh1_691, kh1_692, ki_920, ki_921, \
                         ki_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1183[k] = f_13 * ii_724[k]
                    + f_7 * kh0_690[k]
                    - f_8 * kh1_690[k]
                    + pb_y[k] * ki_920[k];

        t_1184[k] = f_13 * ii_725[k]
                    + f_5 * kh0_691[k]
                    - f_6 * kh1_691[k]
                    + pb_y[k] * ki_921[k];

        t_1185[k] = f_13 * ii_726[k]
                    + f_3 * kh0_692[k]
                    - f_4 * kh1_692[k]
                    + pb_y[k] * ki_922[k];
    }

#pragma omp simd aligned(t_1186, t_1187, t_1188, t_1189, pa_y, pb_x, pb_y, hk0_719, hk1_719, \
                         ii_727, ii_728, ik_935, kh0_693, kh1_693, ki_923, \
                         ki_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1186[k] = f_13 * ii_727[k]
                    + pb_y[k] * ki_923[k];

        t_1187[k] = f_21 * hk0_719[k]
                    - f_22 * hk1_719[k]
                    + pa_y[k] * ik_935[k];

        t_1188[k] = f_1 * kh0_693[k]
                    - f_2 * kh1_693[k]
                    + pb_x[k] * ki_924[k];

        t_1189[k] = f_12 * ii_728[k]
                    + pb_y[k] * ki_924[k];
    }

#pragma omp simd aligned(t_1190, t_1191, t_1192, pb_x, pb_y, pb_z, ii_700, ii_730, kh0_696, \
                         kh1_696, ki_924, ki_926, ki_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1190[k] = f_15 * ii_700[k]
                    + pb_z[k] * ki_924[k];

        t_1191[k] = f_9 * kh0_696[k]
                    - f_10 * kh1_696[k]
                    + pb_x[k] * ki_927[k];

        t_1192[k] = f_12 * ii_730[k]
                    + pb_y[k] * ki_926[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, t_1196, pb_x, pb_y, pb_z, ii_703, ii_733, \
                         kh0_698, kh0_699, kh1_698, kh1_699, ki_927, ki_929, \
                         ki_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = f_9 * kh0_698[k]
                    - f_10 * kh1_698[k]
                    + pb_x[k] * ki_929[k];

        t_1194[k] = f_7 * kh0_699[k]
                    - f_8 * kh1_699[k]
                    + pb_x[k] * ki_930[k];

        t_1195[k] = f_15 * ii_703[k]
                    + pb_z[k] * ki_927[k];

        t_1196[k] = f_12 * ii_733[k]
                    + pb_y[k] * ki_929[k];
    }

#pragma omp simd aligned(t_1197, t_1198, t_1199, pb_x, pb_z, ii_706, kh0_702, kh0_703, \
                         kh1_702, kh1_703, ki_930, ki_933, ki_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1197[k] = f_7 * kh0_702[k]
                    - f_8 * kh1_702[k]
                    + pb_x[k] * ki_933[k];

        t_1198[k] = f_5 * kh0_703[k]
                    - f_6 * kh1_703[k]
                    + pb_x[k] * ki_934[k];

        t_1199[k] = f_15 * ii_706[k]
                    + pb_z[k] * ki_930[k];
    }

#pragma omp simd aligned(t_1200, t_1201, t_1202, pb_x, pb_y, ii_737, kh0_705, kh0_707, \
                         kh1_705, kh1_707, ki_933, ki_936, ki_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1200[k] = f_5 * kh0_705[k]
                    - f_6 * kh1_705[k]
                    + pb_x[k] * ki_936[k];

        t_1201[k] = f_12 * ii_737[k]
                    + pb_y[k] * ki_933[k];

        t_1202[k] = f_5 * kh0_707[k]
                    - f_6 * kh1_707[k]
                    + pb_x[k] * ki_938[k];
    }

#pragma omp simd aligned(t_1203, t_1204, t_1205, pb_x, pb_z, ii_710, kh0_708, kh0_710, \
                         kh1_708, kh1_710, ki_934, ki_939, ki_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1203[k] = f_3 * kh0_708[k]
                    - f_4 * kh1_708[k]
                    + pb_x[k] * ki_939[k];

        t_1204[k] = f_15 * ii_710[k]
                    + pb_z[k] * ki_934[k];

        t_1205[k] = f_3 * kh0_710[k]
                    - f_4 * kh1_710[k]
                    + pb_x[k] * ki_941[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, t_1209, pb_x, pb_y, ii_742, kh0_711, kh0_713, \
                         kh1_711, kh1_713, ki_938, ki_942, ki_944, \
                         ki_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = f_3 * kh0_711[k]
                    - f_4 * kh1_711[k]
                    + pb_x[k] * ki_942[k];

        t_1207[k] = f_12 * ii_742[k]
                    + pb_y[k] * ki_938[k];

        t_1208[k] = f_3 * kh0_713[k]
                    - f_4 * kh1_713[k]
                    + pb_x[k] * ki_944[k];

        t_1209[k] = pb_x[k] * ki_945[k];
    }
}

static auto
compute_prim_kk_electron_repulsion_0_piece9(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t hk0,
                                            const size_t hk1, const size_t ii, const size_t ik,
                                            const size_t kh0, const size_t kh1, const size_t ki,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
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
    const auto f_16 = 3.0 / p;
    const auto f_17 = 0.5 / alpha;
    const auto f_18 = 0.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / alpha;
    const auto f_20 = 2.0 * beta / (alpha * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hk0_676 = buffer.data(hk0 + 676);
    const auto *hk0_755 = buffer.data(hk0 + 755);

    const auto *hk1_676 = buffer.data(hk1 + 676);
    const auto *hk1_755 = buffer.data(hk1 + 755);

    const auto *ii_721 = buffer.data(ii + 721);
    const auto *ii_731 = buffer.data(ii + 731);
    const auto *ii_734 = buffer.data(ii + 734);
    const auto *ii_738 = buffer.data(ii + 738);
    const auto *ii_749 = buffer.data(ii + 749);
    const auto *ii_751 = buffer.data(ii + 751);
    const auto *ii_752 = buffer.data(ii + 752);
    const auto *ii_753 = buffer.data(ii + 753);
    const auto *ii_754 = buffer.data(ii + 754);
    const auto *ii_755 = buffer.data(ii + 755);
    const auto *ii_756 = buffer.data(ii + 756);
    const auto *ii_757 = buffer.data(ii + 757);
    const auto *ii_758 = buffer.data(ii + 758);
    const auto *ii_759 = buffer.data(ii + 759);
    const auto *ii_761 = buffer.data(ii + 761);
    const auto *ii_762 = buffer.data(ii + 762);
    const auto *ii_764 = buffer.data(ii + 764);
    const auto *ii_765 = buffer.data(ii + 765);
    const auto *ii_766 = buffer.data(ii + 766);
    const auto *ii_768 = buffer.data(ii + 768);
    const auto *ii_769 = buffer.data(ii + 769);
    const auto *ii_770 = buffer.data(ii + 770);
    const auto *ii_777 = buffer.data(ii + 777);
    const auto *ii_779 = buffer.data(ii + 779);
    const auto *ii_780 = buffer.data(ii + 780);
    const auto *ii_781 = buffer.data(ii + 781);
    const auto *ii_782 = buffer.data(ii + 782);
    const auto *ii_783 = buffer.data(ii + 783);

    const auto *ik_928 = buffer.data(ik + 928);
    const auto *ik_971 = buffer.data(ik + 971);
    const auto *ik_972 = buffer.data(ik + 972);
    const auto *ik_974 = buffer.data(ik + 974);
    const auto *ik_975 = buffer.data(ik + 975);
    const auto *ik_977 = buffer.data(ik + 977);
    const auto *ik_978 = buffer.data(ik + 978);
    const auto *ik_981 = buffer.data(ik + 981);
    const auto *ik_982 = buffer.data(ik + 982);
    const auto *ik_984 = buffer.data(ik + 984);
    const auto *ik_986 = buffer.data(ik + 986);
    const auto *ik_987 = buffer.data(ik + 987);
    const auto *ik_989 = buffer.data(ik + 989);
    const auto *ik_990 = buffer.data(ik + 990);
    const auto *ik_992 = buffer.data(ik + 992);
    const auto *ik_1000 = buffer.data(ik + 1000);
    const auto *ik_1002 = buffer.data(ik + 1002);
    const auto *ik_1003 = buffer.data(ik + 1003);
    const auto *ik_1004 = buffer.data(ik + 1004);
    const auto *ik_1005 = buffer.data(ik + 1005);
    const auto *ik_1007 = buffer.data(ik + 1007);

    const auto *kh0_710 = buffer.data(kh0 + 710);
    const auto *kh0_711 = buffer.data(kh0 + 711);
    const auto *kh0_712 = buffer.data(kh0 + 712);
    const auto *kh0_713 = buffer.data(kh0 + 713);
    const auto *kh0_735 = buffer.data(kh0 + 735);
    const auto *kh0_738 = buffer.data(kh0 + 738);
    const auto *kh0_740 = buffer.data(kh0 + 740);
    const auto *kh0_741 = buffer.data(kh0 + 741);
    const auto *kh0_744 = buffer.data(kh0 + 744);
    const auto *kh0_745 = buffer.data(kh0 + 745);
    const auto *kh0_747 = buffer.data(kh0 + 747);
    const auto *kh0_749 = buffer.data(kh0 + 749);
    const auto *kh0_750 = buffer.data(kh0 + 750);
    const auto *kh0_752 = buffer.data(kh0 + 752);
    const auto *kh0_753 = buffer.data(kh0 + 753);
    const auto *kh0_754 = buffer.data(kh0 + 754);
    const auto *kh0_755 = buffer.data(kh0 + 755);

    const auto *kh1_710 = buffer.data(kh1 + 710);
    const auto *kh1_711 = buffer.data(kh1 + 711);
    const auto *kh1_712 = buffer.data(kh1 + 712);
    const auto *kh1_713 = buffer.data(kh1 + 713);
    const auto *kh1_735 = buffer.data(kh1 + 735);
    const auto *kh1_738 = buffer.data(kh1 + 738);
    const auto *kh1_740 = buffer.data(kh1 + 740);
    const auto *kh1_741 = buffer.data(kh1 + 741);
    const auto *kh1_744 = buffer.data(kh1 + 744);
    const auto *kh1_745 = buffer.data(kh1 + 745);
    const auto *kh1_747 = buffer.data(kh1 + 747);
    const auto *kh1_749 = buffer.data(kh1 + 749);
    const auto *kh1_750 = buffer.data(kh1 + 750);
    const auto *kh1_752 = buffer.data(kh1 + 752);
    const auto *kh1_753 = buffer.data(kh1 + 753);
    const auto *kh1_754 = buffer.data(kh1 + 754);
    const auto *kh1_755 = buffer.data(kh1 + 755);

    const auto *ki_945 = buffer.data(ki + 945);
    const auto *ki_946 = buffer.data(ki + 946);
    const auto *ki_947 = buffer.data(ki + 947);
    const auto *ki_948 = buffer.data(ki + 948);
    const auto *ki_949 = buffer.data(ki + 949);
    const auto *ki_950 = buffer.data(ki + 950);
    const auto *ki_951 = buffer.data(ki + 951);
    const auto *ki_952 = buffer.data(ki + 952);
    const auto *ki_954 = buffer.data(ki + 954);
    const auto *ki_955 = buffer.data(ki + 955);
    const auto *ki_957 = buffer.data(ki + 957);
    const auto *ki_958 = buffer.data(ki + 958);
    const auto *ki_961 = buffer.data(ki + 961);
    const auto *ki_962 = buffer.data(ki + 962);
    const auto *ki_966 = buffer.data(ki + 966);
    const auto *ki_973 = buffer.data(ki + 973);
    const auto *ki_974 = buffer.data(ki + 974);
    const auto *ki_975 = buffer.data(ki + 975);
    const auto *ki_976 = buffer.data(ki + 976);
    const auto *ki_977 = buffer.data(ki + 977);
    const auto *ki_978 = buffer.data(ki + 978);
    const auto *ki_979 = buffer.data(ki + 979);
    const auto *ki_980 = buffer.data(ki + 980);
    const auto *ki_982 = buffer.data(ki + 982);
    const auto *ki_983 = buffer.data(ki + 983);
    const auto *ki_985 = buffer.data(ki + 985);
    const auto *ki_986 = buffer.data(ki + 986);
    const auto *ki_989 = buffer.data(ki + 989);
    const auto *ki_990 = buffer.data(ki + 990);
    const auto *ki_992 = buffer.data(ki + 992);
    const auto *ki_994 = buffer.data(ki + 994);
    const auto *ki_995 = buffer.data(ki + 995);
    const auto *ki_997 = buffer.data(ki + 997);
    const auto *ki_998 = buffer.data(ki + 998);
    const auto *ki_1000 = buffer.data(ki + 1000);
    const auto *ki_1001 = buffer.data(ki + 1001);
    const auto *ki_1002 = buffer.data(ki + 1002);
    const auto *ki_1003 = buffer.data(ki + 1003);
    const auto *ki_1004 = buffer.data(ki + 1004);
    const auto *ki_1005 = buffer.data(ki + 1005);
    const auto *ki_1006 = buffer.data(ki + 1006);
    const auto *ki_1007 = buffer.data(ki + 1007);

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, t_1214, t_1215, pb_x, ki_946, ki_947, \
                         ki_948, ki_949, ki_950, ki_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = pb_x[k] * ki_946[k];

        t_1211[k] = pb_x[k] * ki_947[k];

        t_1212[k] = pb_x[k] * ki_948[k];

        t_1213[k] = pb_x[k] * ki_949[k];

        t_1214[k] = pb_x[k] * ki_950[k];

        t_1215[k] = pb_x[k] * ki_951[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, pa_z, pb_y, pb_z, hk0_676, hk1_676, ii_721, \
                         ii_751, ik_928, kh0_710, kh1_710, ki_945, \
                         ki_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = f_19 * hk0_676[k]
                    - f_20 * hk1_676[k]
                    + pa_z[k] * ik_928[k];

        t_1217[k] = f_15 * ii_721[k]
                    + pb_z[k] * ki_945[k];

        t_1218[k] = f_12 * ii_751[k]
                    + f_9 * kh0_710[k]
                    - f_10 * kh1_710[k]
                    + pb_y[k] * ki_947[k];
    }

#pragma omp simd aligned(t_1219, t_1220, t_1221, pb_y, ii_752, ii_753, ii_754, kh0_711, \
                         kh0_712, kh0_713, kh1_711, kh1_712, kh1_713, ki_948, ki_949, \
                         ki_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1219[k] = f_12 * ii_752[k]
                    + f_7 * kh0_711[k]
                    - f_8 * kh1_711[k]
                    + pb_y[k] * ki_948[k];

        t_1220[k] = f_12 * ii_753[k]
                    + f_5 * kh0_712[k]
                    - f_6 * kh1_712[k]
                    + pb_y[k] * ki_949[k];

        t_1221[k] = f_12 * ii_754[k]
                    + f_3 * kh0_713[k]
                    - f_4 * kh1_713[k]
                    + pb_y[k] * ki_950[k];
    }

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, t_1226, pa_y, pb_y, hk0_755, hk1_755, \
                         ii_755, ii_756, ik_971, ik_972, ik_974, ki_951, \
                         ki_952 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_12 * ii_755[k]
                    + pb_y[k] * ki_951[k];

        t_1223[k] = f_17 * hk0_755[k]
                    - f_18 * hk1_755[k]
                    + pa_y[k] * ik_971[k];

        t_1224[k] = pa_y[k] * ik_972[k];

        t_1225[k] = f_11 * ii_756[k]
                    + pb_y[k] * ki_952[k];

        t_1226[k] = pa_y[k] * ik_974[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, t_1230, pa_y, pb_y, ii_757, ii_758, ii_759, \
                         ik_975, ik_977, ik_978, ki_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = f_12 * ii_757[k]
                    + pa_y[k] * ik_975[k];

        t_1228[k] = f_11 * ii_758[k]
                    + pb_y[k] * ki_954[k];

        t_1229[k] = pa_y[k] * ik_977[k];

        t_1230[k] = f_13 * ii_759[k]
                    + pa_y[k] * ik_978[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, t_1234, pa_y, pb_y, pb_z, ii_731, ii_761, \
                         ii_762, ik_981, ik_982, ki_955, ki_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = f_16 * ii_731[k]
                    + pb_z[k] * ki_955[k];

        t_1232[k] = f_11 * ii_761[k]
                    + pb_y[k] * ki_957[k];

        t_1233[k] = pa_y[k] * ik_981[k];

        t_1234[k] = f_14 * ii_762[k]
                    + pa_y[k] * ik_982[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, t_1238, pa_y, pb_y, pb_z, ii_734, ii_764, \
                         ii_765, ik_984, ik_986, ki_958, ki_961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_16 * ii_734[k]
                    + pb_z[k] * ki_958[k];

        t_1236[k] = f_12 * ii_764[k]
                    + pa_y[k] * ik_984[k];

        t_1237[k] = f_11 * ii_765[k]
                    + pb_y[k] * ki_961[k];

        t_1238[k] = pa_y[k] * ik_986[k];
    }

#pragma omp simd aligned(t_1239, t_1240, t_1241, t_1242, pa_y, pb_z, ii_738, ii_766, ii_768, \
                         ii_769, ik_987, ik_989, ik_990, ki_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1239[k] = f_15 * ii_766[k]
                    + pa_y[k] * ik_987[k];

        t_1240[k] = f_16 * ii_738[k]
                    + pb_z[k] * ki_962[k];

        t_1241[k] = f_13 * ii_768[k]
                    + pa_y[k] * ik_989[k];

        t_1242[k] = f_12 * ii_769[k]
                    + pa_y[k] * ik_990[k];
    }

#pragma omp simd aligned(t_1243, t_1244, t_1245, t_1246, t_1247, t_1248, pa_y, pb_x, pb_y, \
                         ii_770, ik_992, ki_966, ki_973, ki_974, ki_975, \
                         ki_976 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1243[k] = f_11 * ii_770[k]
                    + pb_y[k] * ki_966[k];

        t_1244[k] = pa_y[k] * ik_992[k];

        t_1245[k] = pb_x[k] * ki_973[k];

        t_1246[k] = pb_x[k] * ki_974[k];

        t_1247[k] = pb_x[k] * ki_975[k];

        t_1248[k] = pb_x[k] * ki_976[k];
    }

#pragma omp simd aligned(t_1249, t_1250, t_1251, t_1252, t_1253, pa_y, pb_x, pb_z, ii_749, \
                         ii_777, ik_1000, ki_973, ki_977, ki_978, \
                         ki_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1249[k] = pb_x[k] * ki_977[k];

        t_1250[k] = pb_x[k] * ki_978[k];

        t_1251[k] = pb_x[k] * ki_979[k];

        t_1252[k] = f_0 * ii_777[k]
                    + pa_y[k] * ik_1000[k];

        t_1253[k] = f_16 * ii_749[k]
                    + pb_z[k] * ki_973[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, t_1257, pa_y, ii_779, ii_780, ii_781, ii_782, \
                         ik_1002, ik_1003, ik_1004, ik_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = f_15 * ii_779[k]
                    + pa_y[k] * ik_1002[k];

        t_1255[k] = f_14 * ii_780[k]
                    + pa_y[k] * ik_1003[k];

        t_1256[k] = f_13 * ii_781[k]
                    + pa_y[k] * ik_1004[k];

        t_1257[k] = f_12 * ii_782[k]
                    + pa_y[k] * ik_1005[k];
    }

#pragma omp simd aligned(t_1258, t_1259, t_1260, t_1261, t_1262, pa_y, pb_x, pb_y, pb_z, \
                         ii_756, ii_783, ik_1007, kh0_735, kh1_735, ki_979, \
                         ki_980 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1258[k] = f_11 * ii_783[k]
                    + pb_y[k] * ki_979[k];

        t_1259[k] = pa_y[k] * ik_1007[k];

        t_1260[k] = f_1 * kh0_735[k]
                    - f_2 * kh1_735[k]
                    + pb_x[k] * ki_980[k];

        t_1261[k] = pb_y[k] * ki_980[k];

        t_1262[k] = f_0 * ii_756[k]
                    + pb_z[k] * ki_980[k];
    }

#pragma omp simd aligned(t_1263, t_1264, t_1265, t_1266, pb_x, pb_y, kh0_738, kh0_740, \
                         kh0_741, kh1_738, kh1_740, kh1_741, ki_982, ki_983, ki_985, \
                         ki_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1263[k] = f_9 * kh0_738[k]
                    - f_10 * kh1_738[k]
                    + pb_x[k] * ki_983[k];

        t_1264[k] = pb_y[k] * ki_982[k];

        t_1265[k] = f_9 * kh0_740[k]
                    - f_10 * kh1_740[k]
                    + pb_x[k] * ki_985[k];

        t_1266[k] = f_7 * kh0_741[k]
                    - f_8 * kh1_741[k]
                    + pb_x[k] * ki_986[k];
    }

#pragma omp simd aligned(t_1267, t_1268, t_1269, t_1270, pb_x, pb_y, pb_z, ii_759, kh0_744, \
                         kh0_745, kh1_744, kh1_745, ki_983, ki_985, ki_989, \
                         ki_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1267[k] = f_0 * ii_759[k]
                    + pb_z[k] * ki_983[k];

        t_1268[k] = pb_y[k] * ki_985[k];

        t_1269[k] = f_7 * kh0_744[k]
                    - f_8 * kh1_744[k]
                    + pb_x[k] * ki_989[k];

        t_1270[k] = f_5 * kh0_745[k]
                    - f_6 * kh1_745[k]
                    + pb_x[k] * ki_990[k];
    }

#pragma omp simd aligned(t_1271, t_1272, t_1273, t_1274, pb_x, pb_y, pb_z, ii_762, kh0_747, \
                         kh0_749, kh1_747, kh1_749, ki_986, ki_989, ki_992, \
                         ki_994 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1271[k] = f_0 * ii_762[k]
                    + pb_z[k] * ki_986[k];

        t_1272[k] = f_5 * kh0_747[k]
                    - f_6 * kh1_747[k]
                    + pb_x[k] * ki_992[k];

        t_1273[k] = pb_y[k] * ki_989[k];

        t_1274[k] = f_5 * kh0_749[k]
                    - f_6 * kh1_749[k]
                    + pb_x[k] * ki_994[k];
    }

#pragma omp simd aligned(t_1275, t_1276, t_1277, pb_x, pb_z, ii_766, kh0_750, kh0_752, \
                         kh1_750, kh1_752, ki_990, ki_995, ki_997 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1275[k] = f_3 * kh0_750[k]
                    - f_4 * kh1_750[k]
                    + pb_x[k] * ki_995[k];

        t_1276[k] = f_0 * ii_766[k]
                    + pb_z[k] * ki_990[k];

        t_1277[k] = f_3 * kh0_752[k]
                    - f_4 * kh1_752[k]
                    + pb_x[k] * ki_997[k];
    }

#pragma omp simd aligned(t_1278, t_1279, t_1280, t_1281, t_1282, pb_x, pb_y, kh0_753, kh0_755, \
                         kh1_753, kh1_755, ki_994, ki_998, ki_1000, ki_1001, \
                         ki_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1278[k] = f_3 * kh0_753[k]
                    - f_4 * kh1_753[k]
                    + pb_x[k] * ki_998[k];

        t_1279[k] = pb_y[k] * ki_994[k];

        t_1280[k] = f_3 * kh0_755[k]
                    - f_4 * kh1_755[k]
                    + pb_x[k] * ki_1000[k];

        t_1281[k] = pb_x[k] * ki_1001[k];

        t_1282[k] = pb_x[k] * ki_1002[k];
    }

#pragma omp simd aligned(t_1283, t_1284, t_1285, t_1286, t_1287, t_1288, pb_x, pb_y, kh0_750, \
                         kh1_750, ki_1001, ki_1003, ki_1004, ki_1005, ki_1006, \
                         ki_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1283[k] = pb_x[k] * ki_1003[k];

        t_1284[k] = pb_x[k] * ki_1004[k];

        t_1285[k] = pb_x[k] * ki_1005[k];

        t_1286[k] = pb_x[k] * ki_1006[k];

        t_1287[k] = pb_x[k] * ki_1007[k];

        t_1288[k] = f_1 * kh0_750[k]
                    - f_2 * kh1_750[k]
                    + pb_y[k] * ki_1001[k];
    }

#pragma omp simd aligned(t_1289, t_1290, t_1291, pb_y, pb_z, ii_777, kh0_752, kh0_753, \
                         kh1_752, kh1_753, ki_1001, ki_1003, ki_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1289[k] = f_0 * ii_777[k]
                    + pb_z[k] * ki_1001[k];

        t_1290[k] = f_9 * kh0_752[k]
                    - f_10 * kh1_752[k]
                    + pb_y[k] * ki_1003[k];

        t_1291[k] = f_7 * kh0_753[k]
                    - f_8 * kh1_753[k]
                    + pb_y[k] * ki_1004[k];
    }

#pragma omp simd aligned(t_1292, t_1293, t_1294, t_1295, pb_y, pb_z, ii_783, kh0_754, kh0_755, \
                         kh1_754, kh1_755, ki_1005, ki_1006, ki_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1292[k] = f_5 * kh0_754[k]
                    - f_6 * kh1_754[k]
                    + pb_y[k] * ki_1005[k];

        t_1293[k] = f_3 * kh0_755[k]
                    - f_4 * kh1_755[k]
                    + pb_y[k] * ki_1006[k];

        t_1294[k] = pb_y[k] * ki_1007[k];

        t_1295[k] = f_0 * ii_783[k]
                    + f_1 * kh0_755[k]
                    - f_2 * kh1_755[k]
                    + pb_z[k] * ki_1007[k];
    }
}

auto
compute_prim_kk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hk0, const size_t hk1,
                                     const size_t ii, const size_t ik, const size_t kh0,
                                     const size_t kh1, const size_t ki, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    compute_prim_kk_electron_repulsion_0_piece0(buffer, target, pa, pb, hk0, hk1, ii, ik, kh0,
                                                kh1, ki, ncols, alpha, beta, p);

    compute_prim_kk_electron_repulsion_0_piece1(buffer, target, pa, pb, hk0, hk1, ii, ik, kh0,
                                                kh1, ki, ncols, alpha, beta, p);

    compute_prim_kk_electron_repulsion_0_piece2(buffer, target, pa, pb, hk0, hk1, ii, ik, kh0,
                                                kh1, ki, ncols, alpha, beta, p);

    compute_prim_kk_electron_repulsion_0_piece3(buffer, target, pa, pb, hk0, hk1, ii, ik, kh0,
                                                kh1, ki, ncols, alpha, beta, p);

    compute_prim_kk_electron_repulsion_0_piece4(buffer, target, pa, pb, hk0, hk1, ii, ik, kh0,
                                                kh1, ki, ncols, alpha, beta, p);

    compute_prim_kk_electron_repulsion_0_piece5(buffer, target, pa, pb, hk0, hk1, ii, ik, kh0,
                                                kh1, ki, ncols, alpha, beta, p);

    compute_prim_kk_electron_repulsion_0_piece6(buffer, target, pa, pb, ii, ik, ki, ncols, p);

    compute_prim_kk_electron_repulsion_0_piece7(buffer, target, pa, pb, ii, ik, kh0, kh1, ki,
                                                ncols, alpha, beta, p);

    compute_prim_kk_electron_repulsion_0_piece8(buffer, target, pa, pb, hk0, hk1, ii, ik, kh0,
                                                kh1, ki, ncols, alpha, beta, p);

    compute_prim_kk_electron_repulsion_0_piece9(buffer, target, pa, pb, hk0, hk1, ii, ik, kh0,
                                                kh1, ki, ncols, alpha, beta, p);
}

}  // namespace simdt2ceri
