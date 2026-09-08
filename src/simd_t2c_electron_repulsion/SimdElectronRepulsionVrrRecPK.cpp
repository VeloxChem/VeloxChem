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
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_20 = buffer.data(si + 20);
    const auto *si_21 = buffer.data(si + 21);
    const auto *si_22 = buffer.data(si + 22);
    const auto *si_23 = buffer.data(si + 23);
    const auto *si_24 = buffer.data(si + 24);
    const auto *si_25 = buffer.data(si + 25);
    const auto *si_26 = buffer.data(si + 26);
    const auto *si_27 = buffer.data(si + 27);

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_5 = buffer.data(sk + 5);
    const auto *sk_6 = buffer.data(sk + 6);
    const auto *sk_9 = buffer.data(sk + 9);
    const auto *sk_10 = buffer.data(sk + 10);
    const auto *sk_12 = buffer.data(sk + 12);
    const auto *sk_14 = buffer.data(sk + 14);
    const auto *sk_15 = buffer.data(sk + 15);
    const auto *sk_17 = buffer.data(sk + 17);
    const auto *sk_18 = buffer.data(sk + 18);
    const auto *sk_20 = buffer.data(sk + 20);
    const auto *sk_28 = buffer.data(sk + 28);
    const auto *sk_30 = buffer.data(sk + 30);
    const auto *sk_31 = buffer.data(sk + 31);
    const auto *sk_32 = buffer.data(sk + 32);
    const auto *sk_33 = buffer.data(sk + 33);
    const auto *sk_35 = buffer.data(sk + 35);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_21 = buffer.data(pi + 21);
    const auto *pi_23 = buffer.data(pi + 23);
    const auto *pi_24 = buffer.data(pi + 24);
    const auto *pi_25 = buffer.data(pi + 25);
    const auto *pi_27 = buffer.data(pi + 27);
    const auto *pi_28 = buffer.data(pi + 28);
    const auto *pi_29 = buffer.data(pi + 29);
    const auto *pi_31 = buffer.data(pi + 31);
    const auto *pi_33 = buffer.data(pi + 33);
    const auto *pi_34 = buffer.data(pi + 34);
    const auto *pi_37 = buffer.data(pi + 37);
    const auto *pi_38 = buffer.data(pi + 38);
    const auto *pi_42 = buffer.data(pi + 42);
    const auto *pi_49 = buffer.data(pi + 49);
    const auto *pi_50 = buffer.data(pi + 50);
    const auto *pi_51 = buffer.data(pi + 51);
    const auto *pi_52 = buffer.data(pi + 52);
    const auto *pi_53 = buffer.data(pi + 53);
    const auto *pi_54 = buffer.data(pi + 54);
    const auto *pi_55 = buffer.data(pi + 55);
    const auto *pi_56 = buffer.data(pi + 56);
    const auto *pi_58 = buffer.data(pi + 58);
    const auto *pi_59 = buffer.data(pi + 59);
    const auto *pi_61 = buffer.data(pi + 61);
    const auto *pi_62 = buffer.data(pi + 62);
    const auto *pi_65 = buffer.data(pi + 65);
    const auto *pi_66 = buffer.data(pi + 66);
    const auto *pi_70 = buffer.data(pi + 70);
    const auto *pi_77 = buffer.data(pi + 77);
    const auto *pi_78 = buffer.data(pi + 78);
    const auto *pi_79 = buffer.data(pi + 79);
    const auto *pi_80 = buffer.data(pi + 80);
    const auto *pi_81 = buffer.data(pi + 81);
    const auto *pi_82 = buffer.data(pi + 82);
    const auto *pi_83 = buffer.data(pi + 83);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pb_y, pb_z, si_0, si_3, sk_0, sk_3, \
                         pi_0, pi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si_0[k]
                 + pa_x[k] * sk_0[k];

        t_1[k] = pb_y[k] * pi_0[k];

        t_2[k] = pb_z[k] * pi_0[k];

        t_3[k] = f_1 * si_3[k]
                 + pa_x[k] * sk_3[k];

        t_4[k] = pb_y[k] * pi_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pb_y, pb_z, si_5, si_6, si_9, sk_5, \
                         sk_6, sk_9, pi_3, pi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * si_5[k]
                 + pa_x[k] * sk_5[k];

        t_6[k] = f_2 * si_6[k]
                 + pa_x[k] * sk_6[k];

        t_7[k] = pb_z[k] * pi_3[k];

        t_8[k] = pb_y[k] * pi_5[k];

        t_9[k] = f_2 * si_9[k]
                 + pa_x[k] * sk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pb_y, pb_z, si_10, si_12, si_14, \
                         sk_10, sk_12, sk_14, pi_6, pi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * si_10[k]
                  + pa_x[k] * sk_10[k];

        t_11[k] = pb_z[k] * pi_6[k];

        t_12[k] = f_3 * si_12[k]
                  + pa_x[k] * sk_12[k];

        t_13[k] = pb_y[k] * pi_9[k];

        t_14[k] = f_3 * si_14[k]
                  + pa_x[k] * sk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_x, pb_y, pb_z, si_15, si_17, si_18, \
                         sk_15, sk_17, sk_18, pi_10, pi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * si_15[k]
                  + pa_x[k] * sk_15[k];

        t_16[k] = pb_z[k] * pi_10[k];

        t_17[k] = f_4 * si_17[k]
                  + pa_x[k] * sk_17[k];

        t_18[k] = f_4 * si_18[k]
                  + pa_x[k] * sk_18[k];

        t_19[k] = pb_y[k] * pi_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, si_20, si_21, si_23, sk_20, \
                         pi_15, pi_21, pi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_4 * si_20[k]
                  + pa_x[k] * sk_20[k];

        t_21[k] = f_5 * si_21[k]
                  + pb_x[k] * pi_21[k];

        t_22[k] = pb_z[k] * pi_15[k];

        t_23[k] = f_5 * si_23[k]
                  + pb_x[k] * pi_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_x, pb_x, pb_y, si_24, si_25, si_27, \
                         sk_28, pi_20, pi_24, pi_25, pi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * si_24[k]
                  + pb_x[k] * pi_24[k];

        t_25[k] = f_5 * si_25[k]
                  + pb_x[k] * pi_25[k];

        t_26[k] = pb_y[k] * pi_20[k];

        t_27[k] = f_5 * si_27[k]
                  + pb_x[k] * pi_27[k];

        t_28[k] = pa_x[k] * sk_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, pa_x, pb_y, pb_z, sk_30, sk_31, \
                         sk_32, sk_33, pi_21, pi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_z[k] * pi_21[k];

        t_30[k] = pa_x[k] * sk_30[k];

        t_31[k] = pa_x[k] * sk_31[k];

        t_32[k] = pa_x[k] * sk_32[k];

        t_33[k] = pa_x[k] * sk_33[k];

        t_34[k] = pb_y[k] * pi_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_x, pa_y, pb_y, pb_z, si_0, si_1, \
                         sk_0, sk_3, sk_35, pi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_x[k] * sk_35[k];

        t_36[k] = pa_y[k] * sk_0[k];

        t_37[k] = f_5 * si_0[k]
                  + pb_y[k] * pi_28[k];

        t_38[k] = pb_z[k] * pi_28[k];

        t_39[k] = f_4 * si_1[k]
                  + pa_y[k] * sk_3[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_y, pb_y, pb_z, si_3, si_5, sk_5, \
                         sk_6, pi_29, pi_31, pi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_z[k] * pi_29[k];

        t_41[k] = pa_y[k] * sk_5[k];

        t_42[k] = f_3 * si_3[k]
                  + pa_y[k] * sk_6[k];

        t_43[k] = pb_z[k] * pi_31[k];

        t_44[k] = f_5 * si_5[k]
                  + pb_y[k] * pi_33[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_y, pb_y, pb_z, si_6, si_8, si_9, \
                         sk_9, sk_10, sk_12, pi_34, pi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_y[k] * sk_9[k];

        t_46[k] = f_2 * si_6[k]
                  + pa_y[k] * sk_10[k];

        t_47[k] = pb_z[k] * pi_34[k];

        t_48[k] = f_4 * si_8[k]
                  + pa_y[k] * sk_12[k];

        t_49[k] = f_5 * si_9[k]
                  + pb_y[k] * pi_37[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_y, pb_z, si_10, si_12, si_13, sk_14, \
                         sk_15, sk_17, sk_18, pi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * sk_14[k];

        t_51[k] = f_1 * si_10[k]
                  + pa_y[k] * sk_15[k];

        t_52[k] = pb_z[k] * pi_38[k];

        t_53[k] = f_3 * si_12[k]
                  + pa_y[k] * sk_17[k];

        t_54[k] = f_4 * si_13[k]
                  + pa_y[k] * sk_18[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, pa_y, pb_x, pb_y, si_14, sk_20, \
                         pi_42, pi_49, pi_50, pi_51, pi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_5 * si_14[k]
                  + pb_y[k] * pi_42[k];

        t_56[k] = pa_y[k] * sk_20[k];

        t_57[k] = pb_x[k] * pi_49[k];

        t_58[k] = pb_x[k] * pi_50[k];

        t_59[k] = pb_x[k] * pi_51[k];

        t_60[k] = pb_x[k] * pi_52[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pa_y, pb_x, pb_z, si_21, sk_28, pi_49, \
                         pi_53, pi_54, pi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_x[k] * pi_53[k];

        t_62[k] = pb_x[k] * pi_54[k];

        t_63[k] = pb_x[k] * pi_55[k];

        t_64[k] = f_0 * si_21[k]
                  + pa_y[k] * sk_28[k];

        t_65[k] = pb_z[k] * pi_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_y, si_23, si_24, si_25, si_26, sk_30, \
                         sk_31, sk_32, sk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_1 * si_23[k]
                  + pa_y[k] * sk_30[k];

        t_67[k] = f_2 * si_24[k]
                  + pa_y[k] * sk_31[k];

        t_68[k] = f_3 * si_25[k]
                  + pa_y[k] * sk_32[k];

        t_69[k] = f_4 * si_26[k]
                  + pa_y[k] * sk_33[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pa_y, pa_z, pb_y, pb_z, si_0, si_27, \
                         sk_0, sk_35, pi_55, pi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_5 * si_27[k]
                  + pb_y[k] * pi_55[k];

        t_71[k] = pa_y[k] * sk_35[k];

        t_72[k] = pa_z[k] * sk_0[k];

        t_73[k] = pb_y[k] * pi_56[k];

        t_74[k] = f_5 * si_0[k]
                  + pb_z[k] * pi_56[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pa_z, pb_y, pb_z, si_2, si_3, sk_3, \
                         sk_5, sk_6, pi_58, pi_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = pa_z[k] * sk_3[k];

        t_76[k] = pb_y[k] * pi_58[k];

        t_77[k] = f_4 * si_2[k]
                  + pa_z[k] * sk_5[k];

        t_78[k] = pa_z[k] * sk_6[k];

        t_79[k] = f_5 * si_3[k]
                  + pb_z[k] * pi_59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pa_z, pb_y, pb_z, si_5, si_6, si_7, \
                         sk_9, sk_10, sk_12, pi_61, pi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_y[k] * pi_61[k];

        t_81[k] = f_3 * si_5[k]
                  + pa_z[k] * sk_9[k];

        t_82[k] = pa_z[k] * sk_10[k];

        t_83[k] = f_5 * si_6[k]
                  + pb_z[k] * pi_62[k];

        t_84[k] = f_4 * si_7[k]
                  + pa_z[k] * sk_12[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_z, pb_y, pb_z, si_9, si_10, si_11, \
                         sk_14, sk_15, sk_17, pi_65, pi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pb_y[k] * pi_65[k];

        t_86[k] = f_2 * si_9[k]
                  + pa_z[k] * sk_14[k];

        t_87[k] = pa_z[k] * sk_15[k];

        t_88[k] = f_5 * si_10[k]
                  + pb_z[k] * pi_66[k];

        t_89[k] = f_4 * si_11[k]
                  + pa_z[k] * sk_17[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pa_z, pb_x, pb_y, si_12, si_14, sk_18, \
                         sk_20, pi_70, pi_77, pi_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_3 * si_12[k]
                  + pa_z[k] * sk_18[k];

        t_91[k] = pb_y[k] * pi_70[k];

        t_92[k] = f_1 * si_14[k]
                  + pa_z[k] * sk_20[k];

        t_93[k] = pb_x[k] * pi_77[k];

        t_94[k] = pb_x[k] * pi_78[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, t_100, pa_z, pb_x, sk_28, pi_79, pi_80, \
                         pi_81, pi_82, pi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pb_x[k] * pi_79[k];

        t_96[k] = pb_x[k] * pi_80[k];

        t_97[k] = pb_x[k] * pi_81[k];

        t_98[k] = pb_x[k] * pi_82[k];

        t_99[k] = pb_x[k] * pi_83[k];

        t_100[k] = pa_z[k] * sk_28[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_z, pb_z, si_21, si_22, si_23, si_24, \
                         sk_30, sk_31, sk_32, pi_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_5 * si_21[k]
                   + pb_z[k] * pi_77[k];

        t_102[k] = f_4 * si_22[k]
                   + pa_z[k] * sk_30[k];

        t_103[k] = f_3 * si_23[k]
                   + pa_z[k] * sk_31[k];

        t_104[k] = f_2 * si_24[k]
                   + pa_z[k] * sk_32[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_z, pb_y, si_25, si_27, sk_33, sk_35, \
                         pi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_1 * si_25[k]
                   + pa_z[k] * sk_33[k];

        t_106[k] = pb_y[k] * pi_83[k];

        t_107[k] = f_0 * si_27[k]
                   + pa_z[k] * sk_35[k];
    }
}

}  // namespace simdt2ceri
