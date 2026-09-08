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


#include "SimdElectronRepulsionVrrRecLP.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_lp_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ks, const size_t kp,
                                     const size_t ls, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 3.0 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.5 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ks_0 = buffer.data(ks + 0);
    const auto *ks_1 = buffer.data(ks + 1);
    const auto *ks_2 = buffer.data(ks + 2);
    const auto *ks_3 = buffer.data(ks + 3);
    const auto *ks_5 = buffer.data(ks + 5);
    const auto *ks_6 = buffer.data(ks + 6);
    const auto *ks_7 = buffer.data(ks + 7);
    const auto *ks_8 = buffer.data(ks + 8);
    const auto *ks_9 = buffer.data(ks + 9);
    const auto *ks_10 = buffer.data(ks + 10);
    const auto *ks_11 = buffer.data(ks + 11);
    const auto *ks_12 = buffer.data(ks + 12);
    const auto *ks_13 = buffer.data(ks + 13);
    const auto *ks_14 = buffer.data(ks + 14);
    const auto *ks_15 = buffer.data(ks + 15);
    const auto *ks_16 = buffer.data(ks + 16);
    const auto *ks_17 = buffer.data(ks + 17);
    const auto *ks_18 = buffer.data(ks + 18);
    const auto *ks_19 = buffer.data(ks + 19);
    const auto *ks_20 = buffer.data(ks + 20);
    const auto *ks_21 = buffer.data(ks + 21);
    const auto *ks_23 = buffer.data(ks + 23);
    const auto *ks_24 = buffer.data(ks + 24);
    const auto *ks_25 = buffer.data(ks + 25);
    const auto *ks_27 = buffer.data(ks + 27);
    const auto *ks_28 = buffer.data(ks + 28);
    const auto *ks_29 = buffer.data(ks + 29);
    const auto *ks_30 = buffer.data(ks + 30);
    const auto *ks_31 = buffer.data(ks + 31);
    const auto *ks_32 = buffer.data(ks + 32);
    const auto *ks_33 = buffer.data(ks + 33);
    const auto *ks_34 = buffer.data(ks + 34);
    const auto *ks_35 = buffer.data(ks + 35);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_60 = buffer.data(kp + 60);
    const auto *kp_62 = buffer.data(kp + 62);
    const auto *kp_63 = buffer.data(kp + 63);
    const auto *kp_81 = buffer.data(kp + 81);
    const auto *kp_85 = buffer.data(kp + 85);
    const auto *kp_88 = buffer.data(kp + 88);
    const auto *kp_89 = buffer.data(kp + 89);
    const auto *kp_91 = buffer.data(kp + 91);
    const auto *kp_92 = buffer.data(kp + 92);
    const auto *kp_94 = buffer.data(kp + 94);
    const auto *kp_95 = buffer.data(kp + 95);
    const auto *kp_97 = buffer.data(kp + 97);
    const auto *kp_98 = buffer.data(kp + 98);
    const auto *kp_100 = buffer.data(kp + 100);
    const auto *kp_101 = buffer.data(kp + 101);
    const auto *kp_103 = buffer.data(kp + 103);
    const auto *kp_104 = buffer.data(kp + 104);
    const auto *kp_107 = buffer.data(kp + 107);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pa_y, pa_z, pb_x, pb_y, pb_z, \
                         ks_0, kp_0, ls_0, ls_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ks_0[k]
                 + pb_x[k] * ls_0[k];

        t_1[k] = pb_y[k] * ls_0[k];

        t_2[k] = pb_z[k] * ls_0[k];

        t_3[k] = pa_y[k] * kp_0[k];

        t_4[k] = f_1 * ks_0[k]
                 + pb_y[k] * ls_1[k];

        t_5[k] = pb_z[k] * ls_1[k];

        t_6[k] = pa_z[k] * kp_0[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pa_y, pb_x, pb_y, pb_z, ks_0, ks_1, \
                         ks_3, kp_6, ls_2, ls_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_y[k] * ls_2[k];

        t_8[k] = f_1 * ks_0[k]
                 + pb_z[k] * ls_2[k];

        t_9[k] = f_2 * ks_3[k]
                 + pb_x[k] * ls_3[k];

        t_10[k] = f_3 * ks_1[k]
                  + pb_y[k] * ls_3[k];

        t_11[k] = pb_z[k] * ls_3[k];

        t_12[k] = pa_y[k] * kp_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_y, pa_z, pb_x, pb_y, pb_z, ks_2, \
                         ks_5, kp_4, kp_8, ls_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_z[k] * kp_4[k];

        t_14[k] = pa_y[k] * kp_8[k];

        t_15[k] = f_2 * ks_5[k]
                  + pb_x[k] * ls_5[k];

        t_16[k] = pb_y[k] * ls_5[k];

        t_17[k] = f_3 * ks_2[k]
                  + pb_z[k] * ls_5[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, pa_z, pb_x, pb_y, pb_z, ks_3, \
                         ks_6, kp_9, kp_10, ls_6, ls_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_4 * ks_6[k]
                  + pb_x[k] * ls_6[k];

        t_19[k] = f_5 * ks_3[k]
                  + pb_y[k] * ls_6[k];

        t_20[k] = pb_z[k] * ls_6[k];

        t_21[k] = pa_z[k] * kp_9[k];

        t_22[k] = pa_z[k] * kp_10[k];

        t_23[k] = f_1 * ks_3[k]
                  + pb_z[k] * ls_7[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, pa_y, pb_x, pb_y, pb_z, ks_5, \
                         ks_9, kp_15, kp_17, ls_8, ls_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * kp_15[k];

        t_25[k] = f_1 * ks_5[k]
                  + pb_y[k] * ls_8[k];

        t_26[k] = pa_y[k] * kp_17[k];

        t_27[k] = f_4 * ks_9[k]
                  + pb_x[k] * ls_9[k];

        t_28[k] = pb_y[k] * ls_9[k];

        t_29[k] = f_5 * ks_5[k]
                  + pb_z[k] * ls_9[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pa_z, pb_x, pb_y, pb_z, ks_6, \
                         ks_10, kp_18, kp_19, ls_10, ls_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_6 * ks_10[k]
                  + pb_x[k] * ls_10[k];

        t_31[k] = f_6 * ks_6[k]
                  + pb_y[k] * ls_10[k];

        t_32[k] = pb_z[k] * ls_10[k];

        t_33[k] = pa_z[k] * kp_18[k];

        t_34[k] = pa_z[k] * kp_19[k];

        t_35[k] = f_1 * ks_6[k]
                  + pb_z[k] * ls_11[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_y, pb_x, pb_y, pb_z, ks_7, ks_8, \
                         ks_9, ks_12, kp_27, ls_12, ls_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_6 * ks_12[k]
                  + pb_x[k] * ls_12[k];

        t_37[k] = f_3 * ks_8[k]
                  + pb_y[k] * ls_12[k];

        t_38[k] = f_3 * ks_7[k]
                  + pb_z[k] * ls_12[k];

        t_39[k] = pa_y[k] * kp_27[k];

        t_40[k] = f_1 * ks_9[k]
                  + pb_y[k] * ls_13[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_y, pb_x, pb_y, pb_z, ks_9, ks_14, \
                         ks_15, kp_29, ls_14, ls_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pa_y[k] * kp_29[k];

        t_42[k] = f_6 * ks_14[k]
                  + pb_x[k] * ls_14[k];

        t_43[k] = pb_y[k] * ls_14[k];

        t_44[k] = f_6 * ks_9[k]
                  + pb_z[k] * ls_14[k];

        t_45[k] = f_5 * ks_15[k]
                  + pb_x[k] * ls_15[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_z, pb_y, pb_z, ks_10, kp_30, kp_31, \
                         ls_15, ls_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * ks_10[k]
                  + pb_y[k] * ls_15[k];

        t_47[k] = pb_z[k] * ls_15[k];

        t_48[k] = pa_z[k] * kp_30[k];

        t_49[k] = pa_z[k] * kp_31[k];

        t_50[k] = f_1 * ks_10[k]
                  + pb_z[k] * ls_16[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, pb_x, pb_y, pb_z, ks_11, ks_12, \
                         ks_13, ks_17, ks_18, ls_17, ls_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_5 * ks_17[k]
                  + pb_x[k] * ls_17[k];

        t_52[k] = f_5 * ks_12[k]
                  + pb_y[k] * ls_17[k];

        t_53[k] = f_3 * ks_11[k]
                  + pb_z[k] * ls_17[k];

        t_54[k] = f_5 * ks_18[k]
                  + pb_x[k] * ls_18[k];

        t_55[k] = f_3 * ks_13[k]
                  + pb_y[k] * ls_18[k];

        t_56[k] = f_5 * ks_12[k]
                  + pb_z[k] * ls_18[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, pb_z, ks_14, \
                         ks_20, kp_42, kp_44, ls_19, ls_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_y[k] * kp_42[k];

        t_58[k] = f_1 * ks_14[k]
                  + pb_y[k] * ls_19[k];

        t_59[k] = pa_y[k] * kp_44[k];

        t_60[k] = f_5 * ks_20[k]
                  + pb_x[k] * ls_20[k];

        t_61[k] = pb_y[k] * ls_20[k];

        t_62[k] = f_4 * ks_14[k]
                  + pb_z[k] * ls_20[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, t_68, pa_z, pb_x, pb_y, pb_z, ks_15, \
                         ks_21, kp_45, kp_46, ls_21, ls_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_3 * ks_21[k]
                  + pb_x[k] * ls_21[k];

        t_64[k] = f_2 * ks_15[k]
                  + pb_y[k] * ls_21[k];

        t_65[k] = pb_z[k] * ls_21[k];

        t_66[k] = pa_z[k] * kp_45[k];

        t_67[k] = pa_z[k] * kp_46[k];

        t_68[k] = f_1 * ks_15[k]
                  + pb_z[k] * ls_22[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, t_74, pb_x, pb_y, pb_z, ks_16, ks_17, \
                         ks_18, ks_23, ks_24, ls_23, ls_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * ks_23[k]
                  + pb_x[k] * ls_23[k];

        t_70[k] = f_6 * ks_17[k]
                  + pb_y[k] * ls_23[k];

        t_71[k] = f_3 * ks_16[k]
                  + pb_z[k] * ls_23[k];

        t_72[k] = f_3 * ks_24[k]
                  + pb_x[k] * ls_24[k];

        t_73[k] = f_5 * ks_18[k]
                  + pb_y[k] * ls_24[k];

        t_74[k] = f_5 * ks_17[k]
                  + pb_z[k] * ls_24[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, pb_z, ks_18, ks_19, \
                         ks_20, ks_25, kp_60, ls_25, ls_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_3 * ks_25[k]
                  + pb_x[k] * ls_25[k];

        t_76[k] = f_3 * ks_19[k]
                  + pb_y[k] * ls_25[k];

        t_77[k] = f_6 * ks_18[k]
                  + pb_z[k] * ls_25[k];

        t_78[k] = pa_y[k] * kp_60[k];

        t_79[k] = f_1 * ks_20[k]
                  + pb_y[k] * ls_26[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pa_y, pb_x, pb_y, pb_z, ks_20, ks_27, \
                         ks_28, kp_62, ls_27, ls_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_y[k] * kp_62[k];

        t_81[k] = f_3 * ks_27[k]
                  + pb_x[k] * ls_27[k];

        t_82[k] = pb_y[k] * ls_27[k];

        t_83[k] = f_2 * ks_20[k]
                  + pb_z[k] * ls_27[k];

        t_84[k] = f_1 * ks_28[k]
                  + pb_x[k] * ls_28[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_x, pa_z, pb_z, kp_63, kp_85, kp_88, \
                         kp_89, ls_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pa_x[k] * kp_85[k];

        t_86[k] = pb_z[k] * ls_28[k];

        t_87[k] = pa_z[k] * kp_63[k];

        t_88[k] = pa_x[k] * kp_88[k];

        t_89[k] = pa_x[k] * kp_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, pa_x, pb_x, ks_30, ks_31, kp_91, \
                         kp_92, kp_94, kp_95, ls_30, ls_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_1 * ks_30[k]
                  + pb_x[k] * ls_30[k];

        t_91[k] = pa_x[k] * kp_91[k];

        t_92[k] = pa_x[k] * kp_92[k];

        t_93[k] = f_1 * ks_31[k]
                  + pb_x[k] * ls_31[k];

        t_94[k] = pa_x[k] * kp_94[k];

        t_95[k] = pa_x[k] * kp_95[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, t_101, pa_x, pb_x, ks_32, ks_33, \
                         kp_97, kp_98, kp_100, kp_101, ls_32, ls_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * ks_32[k]
                  + pb_x[k] * ls_32[k];

        t_97[k] = pa_x[k] * kp_97[k];

        t_98[k] = pa_x[k] * kp_98[k];

        t_99[k] = f_1 * ks_33[k]
                  + pb_x[k] * ls_33[k];

        t_100[k] = pa_x[k] * kp_100[k];

        t_101[k] = pa_x[k] * kp_101[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, t_107, pa_x, pa_y, pb_x, pb_y, \
                         ks_35, kp_81, kp_103, kp_104, kp_107, ls_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = pa_y[k] * kp_81[k];

        t_103[k] = pa_x[k] * kp_103[k];

        t_104[k] = pa_x[k] * kp_104[k];

        t_105[k] = f_1 * ks_35[k]
                   + pb_x[k] * ls_35[k];

        t_106[k] = pb_y[k] * ls_35[k];

        t_107[k] = pa_x[k] * kp_107[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, t_113, t_114, pa_z, pb_x, pb_y, \
                         pb_z, ks_28, kp_85, ls_36, ls_37, ls_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_x[k] * ls_36[k];

        t_109[k] = f_0 * ks_28[k]
                   + pb_y[k] * ls_36[k];

        t_110[k] = pb_z[k] * ls_36[k];

        t_111[k] = pb_x[k] * ls_37[k];

        t_112[k] = pa_z[k] * kp_85[k];

        t_113[k] = f_1 * ks_28[k]
                   + pb_z[k] * ls_37[k];

        t_114[k] = pb_x[k] * ls_38[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, t_120, pb_x, pb_y, pb_z, ks_29, \
                         ks_30, ks_31, ls_38, ls_39, ls_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_2 * ks_30[k]
                   + pb_y[k] * ls_38[k];

        t_116[k] = f_3 * ks_29[k]
                   + pb_z[k] * ls_38[k];

        t_117[k] = pb_x[k] * ls_39[k];

        t_118[k] = f_4 * ks_31[k]
                   + pb_y[k] * ls_39[k];

        t_119[k] = f_5 * ks_30[k]
                   + pb_z[k] * ls_39[k];

        t_120[k] = pb_x[k] * ls_40[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, t_126, pb_x, pb_y, pb_z, ks_31, \
                         ks_32, ks_33, ls_40, ls_41, ls_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_6 * ks_32[k]
                   + pb_y[k] * ls_40[k];

        t_122[k] = f_6 * ks_31[k]
                   + pb_z[k] * ls_40[k];

        t_123[k] = pb_x[k] * ls_41[k];

        t_124[k] = f_5 * ks_33[k]
                   + pb_y[k] * ls_41[k];

        t_125[k] = f_4 * ks_32[k]
                   + pb_z[k] * ls_41[k];

        t_126[k] = pb_x[k] * ls_42[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pa_y, pb_x, pb_y, pb_z, ks_33, \
                         ks_34, ks_35, kp_107, ls_42, ls_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_3 * ks_34[k]
                   + pb_y[k] * ls_42[k];

        t_128[k] = f_2 * ks_33[k]
                   + pb_z[k] * ls_42[k];

        t_129[k] = pb_x[k] * ls_43[k];

        t_130[k] = f_1 * ks_35[k]
                   + pb_y[k] * ls_43[k];

        t_131[k] = pa_y[k] * kp_107[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pb_x, pb_y, pb_z, ks_35, \
                         ls_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_x[k] * ls_44[k];

        t_133[k] = pb_y[k] * ls_44[k];

        t_134[k] = f_0 * ks_35[k]
                   + pb_z[k] * ls_44[k];
    }
}

}  // namespace simdt2ceri
