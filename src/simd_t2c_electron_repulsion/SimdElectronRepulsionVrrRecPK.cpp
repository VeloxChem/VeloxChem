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


#include "SimdElectronRepulsionVrrRecPK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_pk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t si, const size_t sk,
                                     const size_t pi, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_11 = buffer.data(si + 11);
    const auto *si_12 = buffer.data(si + 12);
    const auto *si_13 = buffer.data(si + 13);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_15 = buffer.data(si + 15);
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);
    const auto *si_20 = buffer.data(si + 20);
    const auto *si_21 = buffer.data(si + 21);
    const auto *si_22 = buffer.data(si + 22);
    const auto *si_23 = buffer.data(si + 23);
    const auto *si_24 = buffer.data(si + 24);

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_1 = buffer.data(sk + 1);
    const auto *sk_2 = buffer.data(sk + 2);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_4 = buffer.data(sk + 4);
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

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_13 = buffer.data(pi + 13);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_16 = buffer.data(pi + 16);
    const auto *pi_17 = buffer.data(pi + 17);
    const auto *pi_18 = buffer.data(pi + 18);
    const auto *pi_19 = buffer.data(pi + 19);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_21 = buffer.data(pi + 21);
    const auto *pi_22 = buffer.data(pi + 22);
    const auto *pi_23 = buffer.data(pi + 23);
    const auto *pi_24 = buffer.data(pi + 24);
    const auto *pi_25 = buffer.data(pi + 25);
    const auto *pi_26 = buffer.data(pi + 26);
    const auto *pi_27 = buffer.data(pi + 27);
    const auto *pi_28 = buffer.data(pi + 28);
    const auto *pi_29 = buffer.data(pi + 29);
    const auto *pi_30 = buffer.data(pi + 30);
    const auto *pi_31 = buffer.data(pi + 31);
    const auto *pi_32 = buffer.data(pi + 32);
    const auto *pi_33 = buffer.data(pi + 33);
    const auto *pi_34 = buffer.data(pi + 34);
    const auto *pi_35 = buffer.data(pi + 35);
    const auto *pi_36 = buffer.data(pi + 36);
    const auto *pi_37 = buffer.data(pi + 37);
    const auto *pi_38 = buffer.data(pi + 38);
    const auto *pi_39 = buffer.data(pi + 39);
    const auto *pi_40 = buffer.data(pi + 40);
    const auto *pi_41 = buffer.data(pi + 41);
    const auto *pi_42 = buffer.data(pi + 42);
    const auto *pi_43 = buffer.data(pi + 43);
    const auto *pi_44 = buffer.data(pi + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pb_y, pb_z, si_0, si_3, sk_0, sk_1, \
                         pi_0, pi_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si_0[k]
                 + pa_x[k] * sk_0[k];

        t_1[k] = pb_y[k] * pi_0[k];

        t_2[k] = pb_z[k] * pi_0[k];

        t_3[k] = f_1 * si_3[k]
                 + pa_x[k] * sk_1[k];

        t_4[k] = pb_y[k] * pi_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pb_y, pb_z, si_4, si_5, si_8, sk_2, \
                         sk_3, sk_4, pi_2, pi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * si_4[k]
                 + pa_x[k] * sk_2[k];

        t_6[k] = f_2 * si_5[k]
                 + pa_x[k] * sk_3[k];

        t_7[k] = pb_z[k] * pi_2[k];

        t_8[k] = pb_y[k] * pi_3[k];

        t_9[k] = f_2 * si_8[k]
                 + pa_x[k] * sk_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pb_y, pb_z, si_9, si_11, si_13, \
                         sk_5, sk_6, sk_7, pi_4, pi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * si_9[k]
                  + pa_x[k] * sk_5[k];

        t_11[k] = pb_z[k] * pi_4[k];

        t_12[k] = f_3 * si_11[k]
                  + pa_x[k] * sk_6[k];

        t_13[k] = pb_y[k] * pi_5[k];

        t_14[k] = f_3 * si_13[k]
                  + pa_x[k] * sk_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_x, pb_y, pb_z, si_14, si_15, si_16, \
                         sk_8, sk_9, sk_10, pi_6, pi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * si_14[k]
                  + pa_x[k] * sk_8[k];

        t_16[k] = pb_z[k] * pi_6[k];

        t_17[k] = f_4 * si_15[k]
                  + pa_x[k] * sk_9[k];

        t_18[k] = f_4 * si_16[k]
                  + pa_x[k] * sk_10[k];

        t_19[k] = pb_y[k] * pi_7[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, si_17, si_18, si_20, sk_11, \
                         pi_8, pi_10, pi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_4 * si_17[k]
                  + pa_x[k] * sk_11[k];

        t_21[k] = f_5 * si_18[k]
                  + pb_x[k] * pi_10[k];

        t_22[k] = pb_z[k] * pi_8[k];

        t_23[k] = f_5 * si_20[k]
                  + pb_x[k] * pi_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_x, pb_x, pb_y, si_21, si_22, si_24, \
                         sk_12, pi_9, pi_12, pi_13, pi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * si_21[k]
                  + pb_x[k] * pi_12[k];

        t_25[k] = f_5 * si_22[k]
                  + pb_x[k] * pi_13[k];

        t_26[k] = pb_y[k] * pi_9[k];

        t_27[k] = f_5 * si_24[k]
                  + pb_x[k] * pi_14[k];

        t_28[k] = pa_x[k] * sk_12[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, pa_x, pb_y, pb_z, sk_13, sk_14, \
                         sk_15, sk_16, pi_10, pi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_z[k] * pi_10[k];

        t_30[k] = pa_x[k] * sk_13[k];

        t_31[k] = pa_x[k] * sk_14[k];

        t_32[k] = pa_x[k] * sk_15[k];

        t_33[k] = pa_x[k] * sk_16[k];

        t_34[k] = pb_y[k] * pi_14[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_x, pa_y, pb_y, pb_z, si_0, si_1, \
                         sk_0, sk_1, sk_17, pi_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_x[k] * sk_17[k];

        t_36[k] = pa_y[k] * sk_0[k];

        t_37[k] = f_5 * si_0[k]
                  + pb_y[k] * pi_15[k];

        t_38[k] = pb_z[k] * pi_15[k];

        t_39[k] = f_4 * si_1[k]
                  + pa_y[k] * sk_1[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_y, pb_y, pb_z, si_3, si_4, sk_2, \
                         sk_3, pi_16, pi_17, pi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_z[k] * pi_16[k];

        t_41[k] = pa_y[k] * sk_2[k];

        t_42[k] = f_3 * si_3[k]
                  + pa_y[k] * sk_3[k];

        t_43[k] = pb_z[k] * pi_17[k];

        t_44[k] = f_5 * si_4[k]
                  + pb_y[k] * pi_18[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_y, pb_y, pb_z, si_5, si_7, si_8, \
                         sk_4, sk_5, sk_6, pi_19, pi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_y[k] * sk_4[k];

        t_46[k] = f_2 * si_5[k]
                  + pa_y[k] * sk_5[k];

        t_47[k] = pb_z[k] * pi_19[k];

        t_48[k] = f_4 * si_7[k]
                  + pa_y[k] * sk_6[k];

        t_49[k] = f_5 * si_8[k]
                  + pb_y[k] * pi_20[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_y, pb_z, si_9, si_11, si_12, sk_7, \
                         sk_8, sk_9, sk_10, pi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * sk_7[k];

        t_51[k] = f_1 * si_9[k]
                  + pa_y[k] * sk_8[k];

        t_52[k] = pb_z[k] * pi_21[k];

        t_53[k] = f_3 * si_11[k]
                  + pa_y[k] * sk_9[k];

        t_54[k] = f_4 * si_12[k]
                  + pa_y[k] * sk_10[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, pa_y, pb_x, pb_y, si_13, sk_11, \
                         pi_22, pi_23, pi_24, pi_25, pi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_5 * si_13[k]
                  + pb_y[k] * pi_22[k];

        t_56[k] = pa_y[k] * sk_11[k];

        t_57[k] = pb_x[k] * pi_23[k];

        t_58[k] = pb_x[k] * pi_24[k];

        t_59[k] = pb_x[k] * pi_25[k];

        t_60[k] = pb_x[k] * pi_26[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pa_y, pb_x, pb_z, si_18, sk_12, pi_23, \
                         pi_27, pi_28, pi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_x[k] * pi_27[k];

        t_62[k] = pb_x[k] * pi_28[k];

        t_63[k] = pb_x[k] * pi_29[k];

        t_64[k] = f_0 * si_18[k]
                  + pa_y[k] * sk_12[k];

        t_65[k] = pb_z[k] * pi_23[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_y, si_20, si_21, si_22, si_23, sk_13, \
                         sk_14, sk_15, sk_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_1 * si_20[k]
                  + pa_y[k] * sk_13[k];

        t_67[k] = f_2 * si_21[k]
                  + pa_y[k] * sk_14[k];

        t_68[k] = f_3 * si_22[k]
                  + pa_y[k] * sk_15[k];

        t_69[k] = f_4 * si_23[k]
                  + pa_y[k] * sk_16[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pa_y, pa_z, pb_y, pb_z, si_0, si_24, \
                         sk_0, sk_17, pi_29, pi_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_5 * si_24[k]
                  + pb_y[k] * pi_29[k];

        t_71[k] = pa_y[k] * sk_17[k];

        t_72[k] = pa_z[k] * sk_0[k];

        t_73[k] = pb_y[k] * pi_30[k];

        t_74[k] = f_5 * si_0[k]
                  + pb_z[k] * pi_30[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pa_z, pb_y, pb_z, si_2, si_3, sk_1, \
                         sk_2, sk_3, pi_31, pi_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = pa_z[k] * sk_1[k];

        t_76[k] = pb_y[k] * pi_31[k];

        t_77[k] = f_4 * si_2[k]
                  + pa_z[k] * sk_2[k];

        t_78[k] = pa_z[k] * sk_3[k];

        t_79[k] = f_5 * si_3[k]
                  + pb_z[k] * pi_32[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pa_z, pb_y, pb_z, si_4, si_5, si_6, \
                         sk_4, sk_5, sk_6, pi_33, pi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_y[k] * pi_33[k];

        t_81[k] = f_3 * si_4[k]
                  + pa_z[k] * sk_4[k];

        t_82[k] = pa_z[k] * sk_5[k];

        t_83[k] = f_5 * si_5[k]
                  + pb_z[k] * pi_34[k];

        t_84[k] = f_4 * si_6[k]
                  + pa_z[k] * sk_6[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_z, pb_y, pb_z, si_8, si_9, si_10, \
                         sk_7, sk_8, sk_9, pi_35, pi_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pb_y[k] * pi_35[k];

        t_86[k] = f_2 * si_8[k]
                  + pa_z[k] * sk_7[k];

        t_87[k] = pa_z[k] * sk_8[k];

        t_88[k] = f_5 * si_9[k]
                  + pb_z[k] * pi_36[k];

        t_89[k] = f_4 * si_10[k]
                  + pa_z[k] * sk_9[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pa_z, pb_x, pb_y, si_11, si_13, sk_10, \
                         sk_11, pi_37, pi_38, pi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_3 * si_11[k]
                  + pa_z[k] * sk_10[k];

        t_91[k] = pb_y[k] * pi_37[k];

        t_92[k] = f_1 * si_13[k]
                  + pa_z[k] * sk_11[k];

        t_93[k] = pb_x[k] * pi_38[k];

        t_94[k] = pb_x[k] * pi_39[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, t_100, pa_z, pb_x, sk_12, pi_40, pi_41, \
                         pi_42, pi_43, pi_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pb_x[k] * pi_40[k];

        t_96[k] = pb_x[k] * pi_41[k];

        t_97[k] = pb_x[k] * pi_42[k];

        t_98[k] = pb_x[k] * pi_43[k];

        t_99[k] = pb_x[k] * pi_44[k];

        t_100[k] = pa_z[k] * sk_12[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_z, pb_z, si_18, si_19, si_20, si_21, \
                         sk_13, sk_14, sk_15, pi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_5 * si_18[k]
                   + pb_z[k] * pi_38[k];

        t_102[k] = f_4 * si_19[k]
                   + pa_z[k] * sk_13[k];

        t_103[k] = f_3 * si_20[k]
                   + pa_z[k] * sk_14[k];

        t_104[k] = f_2 * si_21[k]
                   + pa_z[k] * sk_15[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_z, pb_y, si_22, si_24, sk_16, sk_17, \
                         pi_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_1 * si_22[k]
                   + pa_z[k] * sk_16[k];

        t_106[k] = pb_y[k] * pi_44[k];

        t_107[k] = f_0 * si_24[k]
                   + pa_z[k] * sk_17[k];
    }
}

auto
compute_prim_pk_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t si, const size_t sk,
                                     const size_t pi, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_13 = buffer.data(si + 13);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_15 = buffer.data(si + 15);
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);
    const auto *si_20 = buffer.data(si + 20);

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_1 = buffer.data(sk + 1);
    const auto *sk_2 = buffer.data(sk + 2);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_4 = buffer.data(sk + 4);
    const auto *sk_5 = buffer.data(sk + 5);
    const auto *sk_7 = buffer.data(sk + 7);
    const auto *sk_8 = buffer.data(sk + 8);
    const auto *sk_11 = buffer.data(sk + 11);
    const auto *sk_17 = buffer.data(sk + 17);
    const auto *sk_18 = buffer.data(sk + 18);
    const auto *sk_19 = buffer.data(sk + 19);
    const auto *sk_20 = buffer.data(sk + 20);
    const auto *sk_21 = buffer.data(sk + 21);
    const auto *sk_23 = buffer.data(sk + 23);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pb_z, si_0, si_3, si_4, si_5, sk_0, \
                         sk_1, sk_2, sk_3, pi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si_0[k]
                 + pa_x[k] * sk_0[k];

        t_1[k] = pb_z[k] * pi_0[k];

        t_2[k] = f_1 * si_3[k]
                 + pa_x[k] * sk_1[k];

        t_3[k] = f_1 * si_4[k]
                 + pa_x[k] * sk_2[k];

        t_4[k] = f_2 * si_5[k]
                 + pa_x[k] * sk_3[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, si_6, si_7, si_9, si_10, si_13, sk_4, \
                         sk_5, sk_7, sk_8, sk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_2 * si_6[k]
                 + pa_x[k] * sk_4[k];

        t_6[k] = f_3 * si_7[k]
                 + pa_x[k] * sk_5[k];

        t_7[k] = f_3 * si_9[k]
                 + pa_x[k] * sk_7[k];

        t_8[k] = f_4 * si_10[k]
                 + pa_x[k] * sk_8[k];

        t_9[k] = f_4 * si_13[k]
                 + pa_x[k] * sk_11[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_x, pb_x, si_14, si_20, sk_17, \
                         sk_18, sk_19, sk_20, pi_1, pi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * si_14[k]
                  + pb_x[k] * pi_1[k];

        t_11[k] = f_5 * si_20[k]
                  + pb_x[k] * pi_2[k];

        t_12[k] = pa_x[k] * sk_17[k];

        t_13[k] = pa_x[k] * sk_18[k];

        t_14[k] = pa_x[k] * sk_19[k];

        t_15[k] = pa_x[k] * sk_20[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_x, pa_y, pb_y, si_0, si_1, sk_0, \
                         sk_1, sk_21, sk_23, pi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_x[k] * sk_21[k];

        t_17[k] = pa_x[k] * sk_23[k];

        t_18[k] = pa_y[k] * sk_0[k];

        t_19[k] = f_5 * si_0[k]
                  + pb_y[k] * pi_3[k];

        t_20[k] = f_4 * si_1[k]
                  + pa_y[k] * sk_1[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, pa_y, si_3, si_5, si_7, sk_2, \
                         sk_3, sk_4, sk_5, sk_7, sk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_y[k] * sk_2[k];

        t_22[k] = f_3 * si_3[k]
                  + pa_y[k] * sk_3[k];

        t_23[k] = pa_y[k] * sk_4[k];

        t_24[k] = f_2 * si_5[k]
                  + pa_y[k] * sk_5[k];

        t_25[k] = pa_y[k] * sk_7[k];

        t_26[k] = f_1 * si_7[k]
                  + pa_y[k] * sk_8[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_y, pb_z, si_14, si_16, si_17, sk_11, \
                         sk_17, sk_18, sk_19, pi_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_y[k] * sk_11[k];

        t_28[k] = f_0 * si_14[k]
                  + pa_y[k] * sk_17[k];

        t_29[k] = pb_z[k] * pi_4[k];

        t_30[k] = f_1 * si_16[k]
                  + pa_y[k] * sk_18[k];

        t_31[k] = f_2 * si_17[k]
                  + pa_y[k] * sk_19[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pa_y, pa_z, pb_y, si_18, si_19, si_20, \
                         sk_0, sk_20, sk_21, sk_23, pi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * si_18[k]
                  + pa_y[k] * sk_20[k];

        t_33[k] = f_4 * si_19[k]
                  + pa_y[k] * sk_21[k];

        t_34[k] = f_5 * si_20[k]
                  + pb_y[k] * pi_5[k];

        t_35[k] = pa_y[k] * sk_23[k];

        t_36[k] = pa_z[k] * sk_0[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pa_z, pb_z, si_0, si_2, si_4, sk_1, \
                         sk_2, sk_3, sk_4, pi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * si_0[k]
                  + pb_z[k] * pi_6[k];

        t_38[k] = pa_z[k] * sk_1[k];

        t_39[k] = f_4 * si_2[k]
                  + pa_z[k] * sk_2[k];

        t_40[k] = pa_z[k] * sk_3[k];

        t_41[k] = f_3 * si_4[k]
                  + pa_z[k] * sk_4[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_z, si_6, si_9, sk_5, sk_7, sk_8, \
                         sk_11, sk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * sk_5[k];

        t_43[k] = f_2 * si_6[k]
                  + pa_z[k] * sk_7[k];

        t_44[k] = pa_z[k] * sk_8[k];

        t_45[k] = f_1 * si_9[k]
                  + pa_z[k] * sk_11[k];

        t_46[k] = pa_z[k] * sk_17[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_z, pb_z, si_14, si_15, si_16, si_17, \
                         sk_18, sk_19, sk_20, pi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * si_14[k]
                  + pb_z[k] * pi_7[k];

        t_48[k] = f_4 * si_15[k]
                  + pa_z[k] * sk_18[k];

        t_49[k] = f_3 * si_16[k]
                  + pa_z[k] * sk_19[k];

        t_50[k] = f_2 * si_17[k]
                  + pa_z[k] * sk_20[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_z, pb_y, si_18, si_20, sk_21, sk_23, \
                         pi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * si_18[k]
                  + pa_z[k] * sk_21[k];

        t_52[k] = pb_y[k] * pi_8[k];

        t_53[k] = f_0 * si_20[k]
                  + pa_z[k] * sk_23[k];
    }
}

auto
compute_prim_pk_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t si, const size_t sk,
                                     const size_t pi, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_11 = buffer.data(si + 11);
    const auto *si_12 = buffer.data(si + 12);
    const auto *si_13 = buffer.data(si + 13);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_15 = buffer.data(si + 15);
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);
    const auto *si_20 = buffer.data(si + 20);
    const auto *si_21 = buffer.data(si + 21);
    const auto *si_22 = buffer.data(si + 22);

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_1 = buffer.data(sk + 1);
    const auto *sk_2 = buffer.data(sk + 2);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_4 = buffer.data(sk + 4);
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

    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_18 = buffer.data(pi + 18);
    const auto *pi_19 = buffer.data(pi + 19);
    const auto *pi_24 = buffer.data(pi + 24);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, si_0, si_3, si_4, si_5, si_8, sk_0, \
                         sk_1, sk_2, sk_3, sk_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si_0[k]
                 + pa_x[k] * sk_0[k];

        t_1[k] = f_1 * si_3[k]
                 + pa_x[k] * sk_1[k];

        t_2[k] = f_1 * si_4[k]
                 + pa_x[k] * sk_2[k];

        t_3[k] = f_2 * si_5[k]
                 + pa_x[k] * sk_3[k];

        t_4[k] = f_2 * si_8[k]
                 + pa_x[k] * sk_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_x, si_9, si_13, si_14, si_15, sk_5, sk_7, \
                         sk_8, sk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * si_9[k]
                 + pa_x[k] * sk_5[k];

        t_6[k] = f_3 * si_13[k]
                 + pa_x[k] * sk_7[k];

        t_7[k] = f_4 * si_14[k]
                 + pa_x[k] * sk_8[k];

        t_8[k] = f_4 * si_15[k]
                 + pa_x[k] * sk_11[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_y, pb_x, pb_y, si_0, si_1, si_16, si_22, \
                         sk_1, pi_7, pi_8, pi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * si_16[k]
                 + pb_x[k] * pi_7[k];

        t_10[k] = f_5 * si_22[k]
                  + pb_x[k] * pi_8[k];

        t_11[k] = f_5 * si_0[k]
                  + pb_y[k] * pi_9[k];

        t_12[k] = f_4 * si_1[k]
                  + pa_y[k] * sk_1[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_y, si_3, si_5, si_7, si_9, si_11, \
                         sk_3, sk_5, sk_6, sk_8, sk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * si_3[k]
                  + pa_y[k] * sk_3[k];

        t_14[k] = f_2 * si_5[k]
                  + pa_y[k] * sk_5[k];

        t_15[k] = f_4 * si_7[k]
                  + pa_y[k] * sk_6[k];

        t_16[k] = f_1 * si_9[k]
                  + pa_y[k] * sk_8[k];

        t_17[k] = f_3 * si_11[k]
                  + pa_y[k] * sk_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_y, si_12, si_16, si_18, si_19, \
                         si_20, sk_10, sk_12, sk_13, sk_14, sk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_4 * si_12[k]
                  + pa_y[k] * sk_10[k];

        t_19[k] = f_0 * si_16[k]
                  + pa_y[k] * sk_12[k];

        t_20[k] = f_1 * si_18[k]
                  + pa_y[k] * sk_13[k];

        t_21[k] = f_2 * si_19[k]
                  + pa_y[k] * sk_14[k];

        t_22[k] = f_3 * si_20[k]
                  + pa_y[k] * sk_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pa_z, pb_y, si_21, si_22, sk_0, sk_16, \
                         sk_17, pi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * si_21[k]
                  + pa_y[k] * sk_16[k];

        t_24[k] = f_5 * si_22[k]
                  + pb_y[k] * pi_18[k];

        t_25[k] = pa_y[k] * sk_17[k];

        t_26[k] = pa_z[k] * sk_0[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_z, pb_z, si_0, si_2, si_4, si_6, sk_2, \
                         sk_4, sk_6, pi_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * si_0[k]
                  + pb_z[k] * pi_19[k];

        t_28[k] = f_4 * si_2[k]
                  + pa_z[k] * sk_2[k];

        t_29[k] = f_3 * si_4[k]
                  + pa_z[k] * sk_4[k];

        t_30[k] = f_4 * si_6[k]
                  + pa_z[k] * sk_6[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_z, si_8, si_10, si_11, si_13, sk_7, \
                         sk_9, sk_10, sk_11, sk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_2 * si_8[k]
                  + pa_z[k] * sk_7[k];

        t_32[k] = f_4 * si_10[k]
                  + pa_z[k] * sk_9[k];

        t_33[k] = f_3 * si_11[k]
                  + pa_z[k] * sk_10[k];

        t_34[k] = f_1 * si_13[k]
                  + pa_z[k] * sk_11[k];

        t_35[k] = pa_z[k] * sk_12[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_z, si_16, si_17, si_18, si_19, \
                         sk_13, sk_14, sk_15, pi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_5 * si_16[k]
                  + pb_z[k] * pi_24[k];

        t_37[k] = f_4 * si_17[k]
                  + pa_z[k] * sk_13[k];

        t_38[k] = f_3 * si_18[k]
                  + pa_z[k] * sk_14[k];

        t_39[k] = f_2 * si_19[k]
                  + pa_z[k] * sk_15[k];
    }

#pragma omp simd aligned(t_40, t_41, pa_z, si_20, si_22, sk_16, sk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_1 * si_20[k]
                  + pa_z[k] * sk_16[k];

        t_41[k] = f_0 * si_22[k]
                  + pa_z[k] * sk_17[k];
    }
}

auto
compute_prim_pk_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t si, const size_t sk,
                                     const size_t pi, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_11 = buffer.data(si + 11);
    const auto *si_12 = buffer.data(si + 12);
    const auto *si_13 = buffer.data(si + 13);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);
    const auto *si_20 = buffer.data(si + 20);
    const auto *si_21 = buffer.data(si + 21);
    const auto *si_22 = buffer.data(si + 22);
    const auto *si_23 = buffer.data(si + 23);
    const auto *si_24 = buffer.data(si + 24);

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_4 = buffer.data(sk + 4);
    const auto *sk_5 = buffer.data(sk + 5);
    const auto *sk_8 = buffer.data(sk + 8);
    const auto *sk_9 = buffer.data(sk + 9);
    const auto *sk_11 = buffer.data(sk + 11);
    const auto *sk_13 = buffer.data(sk + 13);
    const auto *sk_14 = buffer.data(sk + 14);
    const auto *sk_16 = buffer.data(sk + 16);
    const auto *sk_17 = buffer.data(sk + 17);
    const auto *sk_19 = buffer.data(sk + 19);
    const auto *sk_22 = buffer.data(sk + 22);
    const auto *sk_24 = buffer.data(sk + 24);
    const auto *sk_25 = buffer.data(sk + 25);
    const auto *sk_26 = buffer.data(sk + 26);
    const auto *sk_27 = buffer.data(sk + 27);
    const auto *sk_29 = buffer.data(sk + 29);

    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_13 = buffer.data(pi + 13);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_16 = buffer.data(pi + 16);
    const auto *pi_18 = buffer.data(pi + 18);
    const auto *pi_19 = buffer.data(pi + 19);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_21 = buffer.data(pi + 21);
    const auto *pi_23 = buffer.data(pi + 23);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, si_0, si_3, si_4, si_5, si_8, sk_0, \
                         sk_3, sk_4, sk_5, sk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si_0[k]
                 + pa_x[k] * sk_0[k];

        t_1[k] = f_1 * si_3[k]
                 + pa_x[k] * sk_3[k];

        t_2[k] = f_1 * si_4[k]
                 + pa_x[k] * sk_4[k];

        t_3[k] = f_2 * si_5[k]
                 + pa_x[k] * sk_5[k];

        t_4[k] = f_2 * si_8[k]
                 + pa_x[k] * sk_8[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, si_9, si_13, si_14, si_17, sk_9, \
                         sk_13, sk_14, sk_19, sk_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * si_9[k]
                 + pa_x[k] * sk_9[k];

        t_6[k] = f_3 * si_13[k]
                 + pa_x[k] * sk_13[k];

        t_7[k] = f_4 * si_14[k]
                 + pa_x[k] * sk_14[k];

        t_8[k] = f_4 * si_17[k]
                 + pa_x[k] * sk_19[k];

        t_9[k] = pa_x[k] * sk_22[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_x, pa_y, sk_0, sk_24, sk_25, \
                         sk_26, sk_27, sk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * sk_24[k];

        t_11[k] = pa_x[k] * sk_25[k];

        t_12[k] = pa_x[k] * sk_26[k];

        t_13[k] = pa_x[k] * sk_27[k];

        t_14[k] = pa_x[k] * sk_29[k];

        t_15[k] = pa_y[k] * sk_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, si_1, si_3, si_5, sk_3, sk_4, \
                         sk_5, sk_8, sk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * si_1[k]
                  + pa_y[k] * sk_3[k];

        t_17[k] = pa_y[k] * sk_4[k];

        t_18[k] = f_3 * si_3[k]
                  + pa_y[k] * sk_5[k];

        t_19[k] = pa_y[k] * sk_8[k];

        t_20[k] = f_2 * si_5[k]
                  + pa_y[k] * sk_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, si_7, si_9, si_11, si_12, sk_11, \
                         sk_13, sk_14, sk_16, sk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_4 * si_7[k]
                  + pa_y[k] * sk_11[k];

        t_22[k] = pa_y[k] * sk_13[k];

        t_23[k] = f_1 * si_9[k]
                  + pa_y[k] * sk_14[k];

        t_24[k] = f_3 * si_11[k]
                  + pa_y[k] * sk_16[k];

        t_25[k] = f_4 * si_12[k]
                  + pa_y[k] * sk_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, pa_y, pb_x, si_18, sk_19, sk_22, \
                         pi_8, pi_9, pi_10, pi_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_y[k] * sk_19[k];

        t_27[k] = pb_x[k] * pi_8[k];

        t_28[k] = pb_x[k] * pi_9[k];

        t_29[k] = pb_x[k] * pi_10[k];

        t_30[k] = pb_x[k] * pi_11[k];

        t_31[k] = f_0 * si_18[k]
                  + pa_y[k] * sk_22[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_y, si_20, si_21, si_22, si_23, sk_24, \
                         sk_25, sk_26, sk_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * si_20[k]
                  + pa_y[k] * sk_24[k];

        t_33[k] = f_2 * si_21[k]
                  + pa_y[k] * sk_25[k];

        t_34[k] = f_3 * si_22[k]
                  + pa_y[k] * sk_26[k];

        t_35[k] = f_4 * si_23[k]
                  + pa_y[k] * sk_27[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_y, pa_z, pb_y, pb_z, si_0, si_24, \
                         sk_0, sk_3, sk_29, pi_12, pi_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_5 * si_24[k]
                  + pb_y[k] * pi_12[k];

        t_37[k] = pa_y[k] * sk_29[k];

        t_38[k] = pa_z[k] * sk_0[k];

        t_39[k] = f_5 * si_0[k]
                  + pb_z[k] * pi_13[k];

        t_40[k] = pa_z[k] * sk_3[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_z, pb_y, si_2, si_4, sk_4, sk_5, \
                         sk_8, sk_9, pi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_4 * si_2[k]
                  + pa_z[k] * sk_4[k];

        t_42[k] = pa_z[k] * sk_5[k];

        t_43[k] = pb_y[k] * pi_14[k];

        t_44[k] = f_3 * si_4[k]
                  + pa_z[k] * sk_8[k];

        t_45[k] = pa_z[k] * sk_9[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_z, pb_y, si_6, si_8, si_10, sk_11, \
                         sk_13, sk_14, sk_16, pi_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * si_6[k]
                  + pa_z[k] * sk_11[k];

        t_47[k] = pb_y[k] * pi_15[k];

        t_48[k] = f_2 * si_8[k]
                  + pa_z[k] * sk_13[k];

        t_49[k] = pa_z[k] * sk_14[k];

        t_50[k] = f_4 * si_10[k]
                  + pa_z[k] * sk_16[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pa_z, pb_x, pb_y, si_11, si_13, sk_17, \
                         sk_19, pi_16, pi_19, pi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * si_11[k]
                  + pa_z[k] * sk_17[k];

        t_52[k] = pb_y[k] * pi_16[k];

        t_53[k] = f_1 * si_13[k]
                  + pa_z[k] * sk_19[k];

        t_54[k] = pb_x[k] * pi_19[k];

        t_55[k] = pb_x[k] * pi_20[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_z, pb_x, pb_z, si_18, si_19, sk_22, \
                         sk_24, pi_18, pi_21, pi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * pi_21[k];

        t_57[k] = pb_x[k] * pi_23[k];

        t_58[k] = pa_z[k] * sk_22[k];

        t_59[k] = f_5 * si_18[k]
                  + pb_z[k] * pi_18[k];

        t_60[k] = f_4 * si_19[k]
                  + pa_z[k] * sk_24[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pa_z, pb_y, si_20, si_21, si_22, si_24, \
                         sk_25, sk_26, sk_27, sk_29, pi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * si_20[k]
                  + pa_z[k] * sk_25[k];

        t_62[k] = f_2 * si_21[k]
                  + pa_z[k] * sk_26[k];

        t_63[k] = f_1 * si_22[k]
                  + pa_z[k] * sk_27[k];

        t_64[k] = pb_y[k] * pi_23[k];

        t_65[k] = f_0 * si_24[k]
                  + pa_z[k] * sk_29[k];
    }
}

auto
compute_prim_pk_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t si, const size_t sk,
                                     const size_t pi, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_13 = buffer.data(si + 13);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_15 = buffer.data(si + 15);
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);
    const auto *si_20 = buffer.data(si + 20);

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_1 = buffer.data(sk + 1);
    const auto *sk_2 = buffer.data(sk + 2);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_4 = buffer.data(sk + 4);
    const auto *sk_5 = buffer.data(sk + 5);
    const auto *sk_7 = buffer.data(sk + 7);
    const auto *sk_8 = buffer.data(sk + 8);
    const auto *sk_11 = buffer.data(sk + 11);
    const auto *sk_17 = buffer.data(sk + 17);
    const auto *sk_18 = buffer.data(sk + 18);
    const auto *sk_19 = buffer.data(sk + 19);
    const auto *sk_20 = buffer.data(sk + 20);
    const auto *sk_21 = buffer.data(sk + 21);
    const auto *sk_23 = buffer.data(sk + 23);

    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, si_0, si_3, si_4, si_5, si_6, sk_0, \
                         sk_1, sk_2, sk_3, sk_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si_0[k]
                 + pa_x[k] * sk_0[k];

        t_1[k] = f_1 * si_3[k]
                 + pa_x[k] * sk_1[k];

        t_2[k] = f_1 * si_4[k]
                 + pa_x[k] * sk_2[k];

        t_3[k] = f_2 * si_5[k]
                 + pa_x[k] * sk_3[k];

        t_4[k] = f_2 * si_6[k]
                 + pa_x[k] * sk_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, si_7, si_9, si_10, si_13, sk_5, sk_7, \
                         sk_8, sk_11, sk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * si_7[k]
                 + pa_x[k] * sk_5[k];

        t_6[k] = f_3 * si_9[k]
                 + pa_x[k] * sk_7[k];

        t_7[k] = f_4 * si_10[k]
                 + pa_x[k] * sk_8[k];

        t_8[k] = f_4 * si_13[k]
                 + pa_x[k] * sk_11[k];

        t_9[k] = pa_x[k] * sk_17[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pa_y, si_1, si_3, si_5, sk_0, \
                         sk_1, sk_3, sk_5, sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * sk_23[k];

        t_11[k] = pa_y[k] * sk_0[k];

        t_12[k] = f_4 * si_1[k]
                  + pa_y[k] * sk_1[k];

        t_13[k] = f_3 * si_3[k]
                  + pa_y[k] * sk_3[k];

        t_14[k] = f_2 * si_5[k]
                  + pa_y[k] * sk_5[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, si_7, si_14, si_16, si_17, si_18, \
                         sk_8, sk_17, sk_18, sk_19, sk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * si_7[k]
                  + pa_y[k] * sk_8[k];

        t_16[k] = f_0 * si_14[k]
                  + pa_y[k] * sk_17[k];

        t_17[k] = f_1 * si_16[k]
                  + pa_y[k] * sk_18[k];

        t_18[k] = f_2 * si_17[k]
                  + pa_y[k] * sk_19[k];

        t_19[k] = f_3 * si_18[k]
                  + pa_y[k] * sk_20[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_y, pa_z, pb_y, si_19, si_20, sk_0, sk_21, \
                         sk_23, pi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_4 * si_19[k]
                  + pa_y[k] * sk_21[k];

        t_21[k] = f_5 * si_20[k]
                  + pb_y[k] * pi_2[k];

        t_22[k] = pa_y[k] * sk_23[k];

        t_23[k] = pa_z[k] * sk_0[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_z, si_0, si_2, si_4, si_6, sk_2, \
                         sk_4, sk_7, pi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * si_0[k]
                  + pb_z[k] * pi_3[k];

        t_25[k] = f_4 * si_2[k]
                  + pa_z[k] * sk_2[k];

        t_26[k] = f_3 * si_4[k]
                  + pa_z[k] * sk_4[k];

        t_27[k] = f_2 * si_6[k]
                  + pa_z[k] * sk_7[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_z, pb_z, si_9, si_14, si_15, si_16, \
                         sk_11, sk_17, sk_18, sk_19, pi_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * si_9[k]
                  + pa_z[k] * sk_11[k];

        t_29[k] = pa_z[k] * sk_17[k];

        t_30[k] = f_5 * si_14[k]
                  + pb_z[k] * pi_4[k];

        t_31[k] = f_4 * si_15[k]
                  + pa_z[k] * sk_18[k];

        t_32[k] = f_3 * si_16[k]
                  + pa_z[k] * sk_19[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_z, si_17, si_18, si_20, sk_20, sk_21, \
                         sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_2 * si_17[k]
                  + pa_z[k] * sk_20[k];

        t_34[k] = f_1 * si_18[k]
                  + pa_z[k] * sk_21[k];

        t_35[k] = f_0 * si_20[k]
                  + pa_z[k] * sk_23[k];
    }
}

auto
compute_prim_pk_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t si, const size_t sk, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_11 = buffer.data(si + 11);
    const auto *si_12 = buffer.data(si + 12);
    const auto *si_13 = buffer.data(si + 13);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_15 = buffer.data(si + 15);
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_1 = buffer.data(sk + 1);
    const auto *sk_2 = buffer.data(sk + 2);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_4 = buffer.data(sk + 4);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, si_0, si_1, si_3, si_5, sk_0, sk_1, \
                         sk_3, sk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si_0[k]
                 + pa_x[k] * sk_0[k];

        t_1[k] = f_1 * si_1[k]
                 + pa_y[k] * sk_1[k];

        t_2[k] = f_2 * si_3[k]
                 + pa_y[k] * sk_3[k];

        t_3[k] = f_3 * si_5[k]
                 + pa_y[k] * sk_5[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, si_2, si_4, si_8, si_12, sk_2, sk_4, \
                         sk_8, sk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * si_8[k]
                 + pa_y[k] * sk_8[k];

        t_5[k] = f_0 * si_12[k]
                 + pa_y[k] * sk_12[k];

        t_6[k] = f_1 * si_2[k]
                 + pa_z[k] * sk_2[k];

        t_7[k] = f_2 * si_4[k]
                 + pa_z[k] * sk_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, pa_z, si_6, si_7, si_9, si_10, si_11, \
                         sk_6, sk_7, sk_9, sk_10, sk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * si_6[k]
                 + pa_z[k] * sk_6[k];

        t_9[k] = f_3 * si_7[k]
                 + pa_z[k] * sk_7[k];

        t_10[k] = f_1 * si_9[k]
                  + pa_z[k] * sk_9[k];

        t_11[k] = f_2 * si_10[k]
                  + pa_z[k] * sk_10[k];

        t_12[k] = f_4 * si_11[k]
                  + pa_z[k] * sk_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_z, si_13, si_14, si_15, si_16, \
                         si_17, sk_13, sk_14, sk_15, sk_16, sk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * si_13[k]
                  + pa_z[k] * sk_13[k];

        t_14[k] = f_2 * si_14[k]
                  + pa_z[k] * sk_14[k];

        t_15[k] = f_3 * si_15[k]
                  + pa_z[k] * sk_15[k];

        t_16[k] = f_4 * si_16[k]
                  + pa_z[k] * sk_16[k];

        t_17[k] = f_0 * si_17[k]
                  + pa_z[k] * sk_17[k];
    }
}

auto
compute_prim_pk_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t si, const size_t sk,
                                     const size_t pi, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_11 = buffer.data(si + 11);
    const auto *si_12 = buffer.data(si + 12);
    const auto *si_15 = buffer.data(si + 15);
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);
    const auto *si_20 = buffer.data(si + 20);
    const auto *si_21 = buffer.data(si + 21);
    const auto *si_22 = buffer.data(si + 22);

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_4 = buffer.data(sk + 4);
    const auto *sk_5 = buffer.data(sk + 5);
    const auto *sk_7 = buffer.data(sk + 7);
    const auto *sk_8 = buffer.data(sk + 8);
    const auto *sk_10 = buffer.data(sk + 10);
    const auto *sk_11 = buffer.data(sk + 11);
    const auto *sk_12 = buffer.data(sk + 12);
    const auto *sk_14 = buffer.data(sk + 14);
    const auto *sk_15 = buffer.data(sk + 15);
    const auto *sk_16 = buffer.data(sk + 16);
    const auto *sk_17 = buffer.data(sk + 17);
    const auto *sk_19 = buffer.data(sk + 19);
    const auto *sk_20 = buffer.data(sk + 20);
    const auto *sk_21 = buffer.data(sk + 21);
    const auto *sk_22 = buffer.data(sk + 22);
    const auto *sk_23 = buffer.data(sk + 23);

    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_17 = buffer.data(pi + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, si_0, si_3, si_4, si_5, si_7, sk_0, \
                         sk_3, sk_4, sk_5, sk_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si_0[k]
                 + pa_x[k] * sk_0[k];

        t_1[k] = f_1 * si_3[k]
                 + pa_x[k] * sk_3[k];

        t_2[k] = f_1 * si_4[k]
                 + pa_x[k] * sk_4[k];

        t_3[k] = f_2 * si_5[k]
                 + pa_x[k] * sk_5[k];

        t_4[k] = f_2 * si_7[k]
                 + pa_x[k] * sk_7[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, si_8, si_11, si_12, si_15, sk_8, \
                         sk_11, sk_12, sk_16, sk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * si_8[k]
                 + pa_x[k] * sk_8[k];

        t_6[k] = f_3 * si_11[k]
                 + pa_x[k] * sk_11[k];

        t_7[k] = f_4 * si_12[k]
                 + pa_x[k] * sk_12[k];

        t_8[k] = f_4 * si_15[k]
                 + pa_x[k] * sk_16[k];

        t_9[k] = pa_x[k] * sk_17[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pa_y, si_1, si_3, si_5, sk_0, \
                         sk_3, sk_5, sk_8, sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * sk_23[k];

        t_11[k] = pa_y[k] * sk_0[k];

        t_12[k] = f_4 * si_1[k]
                  + pa_y[k] * sk_3[k];

        t_13[k] = f_3 * si_3[k]
                  + pa_y[k] * sk_5[k];

        t_14[k] = f_2 * si_5[k]
                  + pa_y[k] * sk_8[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pb_x, si_8, si_16, si_18, si_19, \
                         sk_12, sk_17, sk_19, sk_20, pi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * si_8[k]
                  + pa_y[k] * sk_12[k];

        t_16[k] = pb_x[k] * pi_5[k];

        t_17[k] = f_0 * si_16[k]
                  + pa_y[k] * sk_17[k];

        t_18[k] = f_1 * si_18[k]
                  + pa_y[k] * sk_19[k];

        t_19[k] = f_2 * si_19[k]
                  + pa_y[k] * sk_20[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, pb_y, si_20, si_21, si_22, \
                         sk_0, sk_21, sk_22, sk_23, pi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_3 * si_20[k]
                  + pa_y[k] * sk_21[k];

        t_21[k] = f_4 * si_21[k]
                  + pa_y[k] * sk_22[k];

        t_22[k] = f_5 * si_22[k]
                  + pb_y[k] * pi_6[k];

        t_23[k] = pa_y[k] * sk_23[k];

        t_24[k] = pa_z[k] * sk_0[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_z, pb_y, pb_z, si_0, si_2, si_4, sk_4, \
                         sk_7, pi_7, pi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_5 * si_0[k]
                  + pb_z[k] * pi_7[k];

        t_26[k] = f_4 * si_2[k]
                  + pa_z[k] * sk_4[k];

        t_27[k] = pb_y[k] * pi_8[k];

        t_28[k] = f_3 * si_4[k]
                  + pa_z[k] * sk_7[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_z, pb_y, si_6, si_7, si_9, si_10, \
                         sk_10, sk_11, sk_14, sk_15, pi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_4 * si_6[k]
                  + pa_z[k] * sk_10[k];

        t_30[k] = pb_y[k] * pi_9[k];

        t_31[k] = f_2 * si_7[k]
                  + pa_z[k] * sk_11[k];

        t_32[k] = f_4 * si_9[k]
                  + pa_z[k] * sk_14[k];

        t_33[k] = f_3 * si_10[k]
                  + pa_z[k] * sk_15[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_z, pb_x, pb_y, pb_z, si_11, si_16, \
                         sk_16, sk_17, pi_10, pi_12, pi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * pi_10[k];

        t_35[k] = f_1 * si_11[k]
                  + pa_z[k] * sk_16[k];

        t_36[k] = pb_x[k] * pi_17[k];

        t_37[k] = pa_z[k] * sk_17[k];

        t_38[k] = f_5 * si_16[k]
                  + pb_z[k] * pi_12[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_z, pb_y, si_17, si_18, si_19, si_20, \
                         sk_19, sk_20, sk_21, sk_22, pi_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_4 * si_17[k]
                  + pa_z[k] * sk_19[k];

        t_40[k] = f_3 * si_18[k]
                  + pa_z[k] * sk_20[k];

        t_41[k] = f_2 * si_19[k]
                  + pa_z[k] * sk_21[k];

        t_42[k] = f_1 * si_20[k]
                  + pa_z[k] * sk_22[k];

        t_43[k] = pb_y[k] * pi_17[k];
    }

#pragma omp simd aligned(t_44, pa_z, si_22, sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * si_22[k]
                  + pa_z[k] * sk_23[k];
    }
}

auto
compute_prim_pk_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t si, const size_t sk, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_15 = buffer.data(si + 15);
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_20 = buffer.data(si + 20);

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_1 = buffer.data(sk + 1);
    const auto *sk_2 = buffer.data(sk + 2);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_4 = buffer.data(sk + 4);
    const auto *sk_5 = buffer.data(sk + 5);
    const auto *sk_7 = buffer.data(sk + 7);
    const auto *sk_8 = buffer.data(sk + 8);
    const auto *sk_11 = buffer.data(sk + 11);
    const auto *sk_17 = buffer.data(sk + 17);
    const auto *sk_18 = buffer.data(sk + 18);
    const auto *sk_19 = buffer.data(sk + 19);
    const auto *sk_20 = buffer.data(sk + 20);
    const auto *sk_21 = buffer.data(sk + 21);
    const auto *sk_23 = buffer.data(sk + 23);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, si_0, si_1, si_3, si_5, sk_0, sk_1, \
                         sk_3, sk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si_0[k]
                 + pa_x[k] * sk_0[k];

        t_1[k] = f_1 * si_1[k]
                 + pa_y[k] * sk_1[k];

        t_2[k] = f_2 * si_3[k]
                 + pa_y[k] * sk_3[k];

        t_3[k] = f_3 * si_5[k]
                 + pa_y[k] * sk_5[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, si_2, si_4, si_7, si_14, sk_2, sk_4, \
                         sk_8, sk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * si_7[k]
                 + pa_y[k] * sk_8[k];

        t_5[k] = f_0 * si_14[k]
                 + pa_y[k] * sk_17[k];

        t_6[k] = f_1 * si_2[k]
                 + pa_z[k] * sk_2[k];

        t_7[k] = f_2 * si_4[k]
                 + pa_z[k] * sk_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, pa_z, si_6, si_9, si_15, si_16, si_17, \
                         sk_7, sk_11, sk_18, sk_19, sk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_3 * si_6[k]
                 + pa_z[k] * sk_7[k];

        t_9[k] = f_4 * si_9[k]
                 + pa_z[k] * sk_11[k];

        t_10[k] = f_1 * si_15[k]
                  + pa_z[k] * sk_18[k];

        t_11[k] = f_2 * si_16[k]
                  + pa_z[k] * sk_19[k];

        t_12[k] = f_3 * si_17[k]
                  + pa_z[k] * sk_20[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_z, si_18, si_20, sk_21, sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_4 * si_18[k]
                  + pa_z[k] * sk_21[k];

        t_14[k] = f_0 * si_20[k]
                  + pa_z[k] * sk_23[k];
    }
}

auto
compute_prim_pk_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t si, const size_t sk, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_1 = buffer.data(sk + 1);
    const auto *sk_2 = buffer.data(sk + 2);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, si_0, si_1, si_2, sk_0, sk_1, \
                         sk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si_0[k]
                 + pa_x[k] * sk_0[k];

        t_1[k] = f_0 * si_1[k]
                 + pa_y[k] * sk_1[k];

        t_2[k] = f_0 * si_2[k]
                 + pa_z[k] * sk_2[k];
    }
}

auto
compute_prim_pk_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t si, const size_t sk, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_11 = buffer.data(si + 11);
    const auto *si_12 = buffer.data(si + 12);
    const auto *si_13 = buffer.data(si + 13);
    const auto *si_14 = buffer.data(si + 14);

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_1 = buffer.data(sk + 1);
    const auto *sk_2 = buffer.data(sk + 2);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_4 = buffer.data(sk + 4);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, si_0, si_1, si_3, si_5, sk_0, sk_1, \
                         sk_3, sk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si_0[k]
                 + pa_x[k] * sk_0[k];

        t_1[k] = f_1 * si_1[k]
                 + pa_y[k] * sk_1[k];

        t_2[k] = f_2 * si_3[k]
                 + pa_y[k] * sk_3[k];

        t_3[k] = f_3 * si_5[k]
                 + pa_y[k] * sk_5[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, si_2, si_4, si_7, si_9, sk_2, sk_4, \
                         sk_7, sk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * si_7[k]
                 + pa_y[k] * sk_7[k];

        t_5[k] = f_0 * si_9[k]
                 + pa_y[k] * sk_9[k];

        t_6[k] = f_1 * si_2[k]
                 + pa_z[k] * sk_2[k];

        t_7[k] = f_2 * si_4[k]
                 + pa_z[k] * sk_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, pa_z, si_6, si_8, si_10, si_11, si_12, \
                         sk_6, sk_8, sk_10, sk_11, sk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_3 * si_6[k]
                 + pa_z[k] * sk_6[k];

        t_9[k] = f_4 * si_8[k]
                 + pa_z[k] * sk_8[k];

        t_10[k] = f_1 * si_10[k]
                  + pa_z[k] * sk_10[k];

        t_11[k] = f_2 * si_11[k]
                  + pa_z[k] * sk_11[k];

        t_12[k] = f_3 * si_12[k]
                  + pa_z[k] * sk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_z, si_13, si_14, sk_13, sk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_4 * si_13[k]
                  + pa_z[k] * sk_13[k];

        t_14[k] = f_0 * si_14[k]
                  + pa_z[k] * sk_14[k];
    }
}

}  // namespace simdt2ceri
