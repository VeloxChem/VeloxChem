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


#include "SimdElectronRepulsionVrrRecPL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_pl_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sk, const size_t sl,
                                     const size_t pk, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 3.0 / p;
    const auto f_2 = 2.5 / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = 1.5 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;

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

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_1 = buffer.data(sk + 1);
    const auto *sk_2 = buffer.data(sk + 2);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_5 = buffer.data(sk + 5);
    const auto *sk_6 = buffer.data(sk + 6);
    const auto *sk_7 = buffer.data(sk + 7);
    const auto *sk_8 = buffer.data(sk + 8);
    const auto *sk_9 = buffer.data(sk + 9);
    const auto *sk_10 = buffer.data(sk + 10);
    const auto *sk_11 = buffer.data(sk + 11);
    const auto *sk_12 = buffer.data(sk + 12);
    const auto *sk_13 = buffer.data(sk + 13);
    const auto *sk_14 = buffer.data(sk + 14);
    const auto *sk_15 = buffer.data(sk + 15);
    const auto *sk_16 = buffer.data(sk + 16);
    const auto *sk_17 = buffer.data(sk + 17);
    const auto *sk_18 = buffer.data(sk + 18);
    const auto *sk_19 = buffer.data(sk + 19);
    const auto *sk_20 = buffer.data(sk + 20);
    const auto *sk_21 = buffer.data(sk + 21);
    const auto *sk_23 = buffer.data(sk + 23);
    const auto *sk_24 = buffer.data(sk + 24);
    const auto *sk_25 = buffer.data(sk + 25);
    const auto *sk_27 = buffer.data(sk + 27);
    const auto *sk_28 = buffer.data(sk + 28);
    const auto *sk_29 = buffer.data(sk + 29);
    const auto *sk_30 = buffer.data(sk + 30);
    const auto *sk_31 = buffer.data(sk + 31);
    const auto *sk_32 = buffer.data(sk + 32);
    const auto *sk_33 = buffer.data(sk + 33);
    const auto *sk_34 = buffer.data(sk + 34);
    const auto *sk_35 = buffer.data(sk + 35);

    const auto *sl_0 = buffer.data(sl + 0);
    const auto *sl_3 = buffer.data(sl + 3);
    const auto *sl_5 = buffer.data(sl + 5);
    const auto *sl_6 = buffer.data(sl + 6);
    const auto *sl_9 = buffer.data(sl + 9);
    const auto *sl_10 = buffer.data(sl + 10);
    const auto *sl_12 = buffer.data(sl + 12);
    const auto *sl_14 = buffer.data(sl + 14);
    const auto *sl_15 = buffer.data(sl + 15);
    const auto *sl_17 = buffer.data(sl + 17);
    const auto *sl_18 = buffer.data(sl + 18);
    const auto *sl_20 = buffer.data(sl + 20);
    const auto *sl_21 = buffer.data(sl + 21);
    const auto *sl_23 = buffer.data(sl + 23);
    const auto *sl_24 = buffer.data(sl + 24);
    const auto *sl_25 = buffer.data(sl + 25);
    const auto *sl_27 = buffer.data(sl + 27);
    const auto *sl_36 = buffer.data(sl + 36);
    const auto *sl_38 = buffer.data(sl + 38);
    const auto *sl_39 = buffer.data(sl + 39);
    const auto *sl_40 = buffer.data(sl + 40);
    const auto *sl_41 = buffer.data(sl + 41);
    const auto *sl_42 = buffer.data(sl + 42);
    const auto *sl_44 = buffer.data(sl + 44);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_2 = buffer.data(pk + 2);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_9 = buffer.data(pk + 9);
    const auto *pk_10 = buffer.data(pk + 10);
    const auto *pk_14 = buffer.data(pk + 14);
    const auto *pk_15 = buffer.data(pk + 15);
    const auto *pk_20 = buffer.data(pk + 20);
    const auto *pk_21 = buffer.data(pk + 21);
    const auto *pk_27 = buffer.data(pk + 27);
    const auto *pk_28 = buffer.data(pk + 28);
    const auto *pk_30 = buffer.data(pk + 30);
    const auto *pk_31 = buffer.data(pk + 31);
    const auto *pk_32 = buffer.data(pk + 32);
    const auto *pk_33 = buffer.data(pk + 33);
    const auto *pk_35 = buffer.data(pk + 35);
    const auto *pk_36 = buffer.data(pk + 36);
    const auto *pk_37 = buffer.data(pk + 37);
    const auto *pk_39 = buffer.data(pk + 39);
    const auto *pk_41 = buffer.data(pk + 41);
    const auto *pk_42 = buffer.data(pk + 42);
    const auto *pk_45 = buffer.data(pk + 45);
    const auto *pk_46 = buffer.data(pk + 46);
    const auto *pk_50 = buffer.data(pk + 50);
    const auto *pk_51 = buffer.data(pk + 51);
    const auto *pk_56 = buffer.data(pk + 56);
    const auto *pk_64 = buffer.data(pk + 64);
    const auto *pk_65 = buffer.data(pk + 65);
    const auto *pk_66 = buffer.data(pk + 66);
    const auto *pk_67 = buffer.data(pk + 67);
    const auto *pk_68 = buffer.data(pk + 68);
    const auto *pk_69 = buffer.data(pk + 69);
    const auto *pk_70 = buffer.data(pk + 70);
    const auto *pk_71 = buffer.data(pk + 71);
    const auto *pk_72 = buffer.data(pk + 72);
    const auto *pk_74 = buffer.data(pk + 74);
    const auto *pk_75 = buffer.data(pk + 75);
    const auto *pk_77 = buffer.data(pk + 77);
    const auto *pk_78 = buffer.data(pk + 78);
    const auto *pk_81 = buffer.data(pk + 81);
    const auto *pk_82 = buffer.data(pk + 82);
    const auto *pk_86 = buffer.data(pk + 86);
    const auto *pk_87 = buffer.data(pk + 87);
    const auto *pk_92 = buffer.data(pk + 92);
    const auto *pk_100 = buffer.data(pk + 100);
    const auto *pk_101 = buffer.data(pk + 101);
    const auto *pk_102 = buffer.data(pk + 102);
    const auto *pk_103 = buffer.data(pk + 103);
    const auto *pk_104 = buffer.data(pk + 104);
    const auto *pk_105 = buffer.data(pk + 105);
    const auto *pk_106 = buffer.data(pk + 106);
    const auto *pk_107 = buffer.data(pk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pb_y, pb_z, sk_0, sk_3, sl_0, sl_3, \
                         pk_0, pk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sk_0[k]
                 + pa_x[k] * sl_0[k];

        t_1[k] = pb_y[k] * pk_0[k];

        t_2[k] = pb_z[k] * pk_0[k];

        t_3[k] = f_1 * sk_3[k]
                 + pa_x[k] * sl_3[k];

        t_4[k] = pb_y[k] * pk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pb_y, pb_z, sk_5, sk_6, sk_9, sl_5, \
                         sl_6, sl_9, pk_3, pk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sk_5[k]
                 + pa_x[k] * sl_5[k];

        t_6[k] = f_2 * sk_6[k]
                 + pa_x[k] * sl_6[k];

        t_7[k] = pb_z[k] * pk_3[k];

        t_8[k] = pb_y[k] * pk_5[k];

        t_9[k] = f_2 * sk_9[k]
                 + pa_x[k] * sl_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pb_y, pb_z, sk_10, sk_12, sk_14, \
                         sl_10, sl_12, sl_14, pk_6, pk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * sk_10[k]
                  + pa_x[k] * sl_10[k];

        t_11[k] = pb_z[k] * pk_6[k];

        t_12[k] = f_3 * sk_12[k]
                  + pa_x[k] * sl_12[k];

        t_13[k] = pb_y[k] * pk_9[k];

        t_14[k] = f_3 * sk_14[k]
                  + pa_x[k] * sl_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_x, pb_y, pb_z, sk_15, sk_17, sk_18, \
                         sl_15, sl_17, sl_18, pk_10, pk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * sk_15[k]
                  + pa_x[k] * sl_15[k];

        t_16[k] = pb_z[k] * pk_10[k];

        t_17[k] = f_4 * sk_17[k]
                  + pa_x[k] * sl_17[k];

        t_18[k] = f_4 * sk_18[k]
                  + pa_x[k] * sl_18[k];

        t_19[k] = pb_y[k] * pk_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_x, pb_z, sk_20, sk_21, sk_23, sk_24, \
                         sl_20, sl_21, sl_23, sl_24, pk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_4 * sk_20[k]
                  + pa_x[k] * sl_20[k];

        t_21[k] = f_5 * sk_21[k]
                  + pa_x[k] * sl_21[k];

        t_22[k] = pb_z[k] * pk_15[k];

        t_23[k] = f_5 * sk_23[k]
                  + pa_x[k] * sl_23[k];

        t_24[k] = f_5 * sk_24[k]
                  + pa_x[k] * sl_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pb_x, pb_y, sk_25, sk_27, sk_28, sl_25, \
                         sl_27, pk_20, pk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_5 * sk_25[k]
                  + pa_x[k] * sl_25[k];

        t_26[k] = pb_y[k] * pk_20[k];

        t_27[k] = f_5 * sk_27[k]
                  + pa_x[k] * sl_27[k];

        t_28[k] = f_6 * sk_28[k]
                  + pb_x[k] * pk_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pb_x, pb_z, sk_30, sk_31, sk_32, sk_33, \
                         pk_21, pk_30, pk_31, pk_32, pk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_z[k] * pk_21[k];

        t_30[k] = f_6 * sk_30[k]
                  + pb_x[k] * pk_30[k];

        t_31[k] = f_6 * sk_31[k]
                  + pb_x[k] * pk_31[k];

        t_32[k] = f_6 * sk_32[k]
                  + pb_x[k] * pk_32[k];

        t_33[k] = f_6 * sk_33[k]
                  + pb_x[k] * pk_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_x, pb_x, pb_y, pb_z, sk_35, sl_36, \
                         sl_38, pk_27, pk_28, pk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * pk_27[k];

        t_35[k] = f_6 * sk_35[k]
                  + pb_x[k] * pk_35[k];

        t_36[k] = pa_x[k] * sl_36[k];

        t_37[k] = pb_z[k] * pk_28[k];

        t_38[k] = pa_x[k] * sl_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, t_44, pa_x, pb_y, sl_39, sl_40, sl_41, \
                         sl_42, sl_44, pk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pa_x[k] * sl_39[k];

        t_40[k] = pa_x[k] * sl_40[k];

        t_41[k] = pa_x[k] * sl_41[k];

        t_42[k] = pa_x[k] * sl_42[k];

        t_43[k] = pb_y[k] * pk_35[k];

        t_44[k] = pa_x[k] * sl_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, pa_y, pb_y, pb_z, sk_0, sk_1, \
                         sl_0, sl_3, sl_5, pk_36, pk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_y[k] * sl_0[k];

        t_46[k] = f_6 * sk_0[k]
                  + pb_y[k] * pk_36[k];

        t_47[k] = pb_z[k] * pk_36[k];

        t_48[k] = f_5 * sk_1[k]
                  + pa_y[k] * sl_3[k];

        t_49[k] = pb_z[k] * pk_37[k];

        t_50[k] = pa_y[k] * sl_5[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pa_y, pb_y, pb_z, sk_3, sk_5, sk_6, \
                         sl_6, sl_9, sl_10, pk_39, pk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_4 * sk_3[k]
                  + pa_y[k] * sl_6[k];

        t_52[k] = pb_z[k] * pk_39[k];

        t_53[k] = f_6 * sk_5[k]
                  + pb_y[k] * pk_41[k];

        t_54[k] = pa_y[k] * sl_9[k];

        t_55[k] = f_3 * sk_6[k]
                  + pa_y[k] * sl_10[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_y, pb_y, pb_z, sk_8, sk_9, sk_10, \
                         sl_12, sl_14, sl_15, pk_42, pk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * pk_42[k];

        t_57[k] = f_5 * sk_8[k]
                  + pa_y[k] * sl_12[k];

        t_58[k] = f_6 * sk_9[k]
                  + pb_y[k] * pk_45[k];

        t_59[k] = pa_y[k] * sl_14[k];

        t_60[k] = f_2 * sk_10[k]
                  + pa_y[k] * sl_15[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pa_y, pb_y, pb_z, sk_12, sk_13, sk_14, \
                         sl_17, sl_18, sl_20, pk_46, pk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_z[k] * pk_46[k];

        t_62[k] = f_4 * sk_12[k]
                  + pa_y[k] * sl_17[k];

        t_63[k] = f_5 * sk_13[k]
                  + pa_y[k] * sl_18[k];

        t_64[k] = f_6 * sk_14[k]
                  + pb_y[k] * pk_50[k];

        t_65[k] = pa_y[k] * sl_20[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_y, pb_z, sk_15, sk_17, sk_18, sk_19, \
                         sl_21, sl_23, sl_24, sl_25, pk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_1 * sk_15[k]
                  + pa_y[k] * sl_21[k];

        t_67[k] = pb_z[k] * pk_51[k];

        t_68[k] = f_3 * sk_17[k]
                  + pa_y[k] * sl_23[k];

        t_69[k] = f_4 * sk_18[k]
                  + pa_y[k] * sl_24[k];

        t_70[k] = f_5 * sk_19[k]
                  + pa_y[k] * sl_25[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, pa_y, pb_x, pb_y, sk_20, sl_27, \
                         pk_56, pk_64, pk_65, pk_66, pk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_6 * sk_20[k]
                  + pb_y[k] * pk_56[k];

        t_72[k] = pa_y[k] * sl_27[k];

        t_73[k] = pb_x[k] * pk_64[k];

        t_74[k] = pb_x[k] * pk_65[k];

        t_75[k] = pb_x[k] * pk_66[k];

        t_76[k] = pb_x[k] * pk_67[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, pa_y, pb_x, pb_z, sk_28, sl_36, \
                         pk_64, pk_68, pk_69, pk_70, pk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_x[k] * pk_68[k];

        t_78[k] = pb_x[k] * pk_69[k];

        t_79[k] = pb_x[k] * pk_70[k];

        t_80[k] = pb_x[k] * pk_71[k];

        t_81[k] = f_0 * sk_28[k]
                  + pa_y[k] * sl_36[k];

        t_82[k] = pb_z[k] * pk_64[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, sk_30, sk_31, sk_32, sk_33, \
                         sk_34, sl_38, sl_39, sl_40, sl_41, sl_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * sk_30[k]
                  + pa_y[k] * sl_38[k];

        t_84[k] = f_2 * sk_31[k]
                  + pa_y[k] * sl_39[k];

        t_85[k] = f_3 * sk_32[k]
                  + pa_y[k] * sl_40[k];

        t_86[k] = f_4 * sk_33[k]
                  + pa_y[k] * sl_41[k];

        t_87[k] = f_5 * sk_34[k]
                  + pa_y[k] * sl_42[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_y, pa_z, pb_y, pb_z, sk_0, sk_35, \
                         sl_0, sl_44, pk_71, pk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_6 * sk_35[k]
                  + pb_y[k] * pk_71[k];

        t_89[k] = pa_y[k] * sl_44[k];

        t_90[k] = pa_z[k] * sl_0[k];

        t_91[k] = pb_y[k] * pk_72[k];

        t_92[k] = f_6 * sk_0[k]
                  + pb_z[k] * pk_72[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pa_z, pb_y, pb_z, sk_2, sk_3, sl_3, \
                         sl_5, sl_6, pk_74, pk_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_z[k] * sl_3[k];

        t_94[k] = pb_y[k] * pk_74[k];

        t_95[k] = f_5 * sk_2[k]
                  + pa_z[k] * sl_5[k];

        t_96[k] = pa_z[k] * sl_6[k];

        t_97[k] = f_6 * sk_3[k]
                  + pb_z[k] * pk_75[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_z, pb_y, pb_z, sk_5, sk_6, sk_7, \
                         sl_9, sl_10, sl_12, pk_77, pk_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = pb_y[k] * pk_77[k];

        t_99[k] = f_4 * sk_5[k]
                  + pa_z[k] * sl_9[k];

        t_100[k] = pa_z[k] * sl_10[k];

        t_101[k] = f_6 * sk_6[k]
                   + pb_z[k] * pk_78[k];

        t_102[k] = f_5 * sk_7[k]
                   + pa_z[k] * sl_12[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, pb_z, sk_9, sk_10, \
                         sk_11, sl_14, sl_15, sl_17, pk_81, pk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pb_y[k] * pk_81[k];

        t_104[k] = f_3 * sk_9[k]
                   + pa_z[k] * sl_14[k];

        t_105[k] = pa_z[k] * sl_15[k];

        t_106[k] = f_6 * sk_10[k]
                   + pb_z[k] * pk_82[k];

        t_107[k] = f_5 * sk_11[k]
                   + pa_z[k] * sl_17[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, sk_12, sk_14, \
                         sk_15, sl_18, sl_20, sl_21, pk_86, pk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_4 * sk_12[k]
                   + pa_z[k] * sl_18[k];

        t_109[k] = pb_y[k] * pk_86[k];

        t_110[k] = f_2 * sk_14[k]
                   + pa_z[k] * sl_20[k];

        t_111[k] = pa_z[k] * sl_21[k];

        t_112[k] = f_6 * sk_15[k]
                   + pb_z[k] * pk_87[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_z, pb_y, sk_16, sk_17, sk_18, \
                         sk_20, sl_23, sl_24, sl_25, sl_27, pk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_5 * sk_16[k]
                   + pa_z[k] * sl_23[k];

        t_114[k] = f_4 * sk_17[k]
                   + pa_z[k] * sl_24[k];

        t_115[k] = f_3 * sk_18[k]
                   + pa_z[k] * sl_25[k];

        t_116[k] = pb_y[k] * pk_92[k];

        t_117[k] = f_1 * sk_20[k]
                   + pa_z[k] * sl_27[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, t_123, t_124, pb_x, pk_100, \
                         pk_101, pk_102, pk_103, pk_104, pk_105, \
                         pk_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pb_x[k] * pk_100[k];

        t_119[k] = pb_x[k] * pk_101[k];

        t_120[k] = pb_x[k] * pk_102[k];

        t_121[k] = pb_x[k] * pk_103[k];

        t_122[k] = pb_x[k] * pk_104[k];

        t_123[k] = pb_x[k] * pk_105[k];

        t_124[k] = pb_x[k] * pk_106[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, pa_z, pb_x, pb_z, sk_28, sk_29, \
                         sk_30, sl_36, sl_38, sl_39, pk_100, pk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_x[k] * pk_107[k];

        t_126[k] = pa_z[k] * sl_36[k];

        t_127[k] = f_6 * sk_28[k]
                   + pb_z[k] * pk_100[k];

        t_128[k] = f_5 * sk_29[k]
                   + pa_z[k] * sl_38[k];

        t_129[k] = f_4 * sk_30[k]
                   + pa_z[k] * sl_39[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, pa_z, pb_y, sk_31, sk_32, sk_33, \
                         sk_35, sl_40, sl_41, sl_42, sl_44, pk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * sk_31[k]
                   + pa_z[k] * sl_40[k];

        t_131[k] = f_2 * sk_32[k]
                   + pa_z[k] * sl_41[k];

        t_132[k] = f_1 * sk_33[k]
                   + pa_z[k] * sl_42[k];

        t_133[k] = pb_y[k] * pk_107[k];

        t_134[k] = f_0 * sk_35[k]
                   + pa_z[k] * sl_44[k];
    }
}

}  // namespace simdt2ceri
