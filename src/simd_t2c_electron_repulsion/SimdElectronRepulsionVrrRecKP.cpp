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


#include "SimdElectronRepulsionVrrRecKP.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_kp_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t is, const size_t ip,
                                     const size_t ks, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 2.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.5 / p;

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

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);
    const auto *is_23 = buffer.data(is + 23);
    const auto *is_24 = buffer.data(is + 24);
    const auto *is_25 = buffer.data(is + 25);
    const auto *is_26 = buffer.data(is + 26);
    const auto *is_27 = buffer.data(is + 27);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_30 = buffer.data(ip + 30);
    const auto *ip_31 = buffer.data(ip + 31);
    const auto *ip_42 = buffer.data(ip + 42);
    const auto *ip_44 = buffer.data(ip + 44);
    const auto *ip_45 = buffer.data(ip + 45);
    const auto *ip_60 = buffer.data(ip + 60);
    const auto *ip_64 = buffer.data(ip + 64);
    const auto *ip_67 = buffer.data(ip + 67);
    const auto *ip_68 = buffer.data(ip + 68);
    const auto *ip_70 = buffer.data(ip + 70);
    const auto *ip_71 = buffer.data(ip + 71);
    const auto *ip_73 = buffer.data(ip + 73);
    const auto *ip_74 = buffer.data(ip + 74);
    const auto *ip_76 = buffer.data(ip + 76);
    const auto *ip_77 = buffer.data(ip + 77);
    const auto *ip_79 = buffer.data(ip + 79);
    const auto *ip_80 = buffer.data(ip + 80);
    const auto *ip_83 = buffer.data(ip + 83);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pa_y, pa_z, pb_x, pb_y, pb_z, \
                         is_0, ip_0, ks_0, ks_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * is_0[k]
                 + pb_x[k] * ks_0[k];

        t_1[k] = pb_y[k] * ks_0[k];

        t_2[k] = pb_z[k] * ks_0[k];

        t_3[k] = pa_y[k] * ip_0[k];

        t_4[k] = f_1 * is_0[k]
                 + pb_y[k] * ks_1[k];

        t_5[k] = pb_z[k] * ks_1[k];

        t_6[k] = pa_z[k] * ip_0[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pa_y, pb_x, pb_y, pb_z, is_0, is_1, \
                         is_3, ip_6, ks_2, ks_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_y[k] * ks_2[k];

        t_8[k] = f_1 * is_0[k]
                 + pb_z[k] * ks_2[k];

        t_9[k] = f_2 * is_3[k]
                 + pb_x[k] * ks_3[k];

        t_10[k] = f_3 * is_1[k]
                  + pb_y[k] * ks_3[k];

        t_11[k] = pb_z[k] * ks_3[k];

        t_12[k] = pa_y[k] * ip_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_y, pa_z, pb_x, pb_y, pb_z, is_2, \
                         is_5, ip_4, ip_8, ks_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_z[k] * ip_4[k];

        t_14[k] = pa_y[k] * ip_8[k];

        t_15[k] = f_2 * is_5[k]
                  + pb_x[k] * ks_5[k];

        t_16[k] = pb_y[k] * ks_5[k];

        t_17[k] = f_3 * is_2[k]
                  + pb_z[k] * ks_5[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, pa_z, pb_x, pb_y, pb_z, is_3, \
                         is_6, ip_9, ip_10, ks_6, ks_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_4 * is_6[k]
                  + pb_x[k] * ks_6[k];

        t_19[k] = f_5 * is_3[k]
                  + pb_y[k] * ks_6[k];

        t_20[k] = pb_z[k] * ks_6[k];

        t_21[k] = pa_z[k] * ip_9[k];

        t_22[k] = pa_z[k] * ip_10[k];

        t_23[k] = f_1 * is_3[k]
                  + pb_z[k] * ks_7[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, pa_y, pb_x, pb_y, pb_z, is_5, \
                         is_9, ip_15, ip_17, ks_8, ks_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * ip_15[k];

        t_25[k] = f_1 * is_5[k]
                  + pb_y[k] * ks_8[k];

        t_26[k] = pa_y[k] * ip_17[k];

        t_27[k] = f_4 * is_9[k]
                  + pb_x[k] * ks_9[k];

        t_28[k] = pb_y[k] * ks_9[k];

        t_29[k] = f_5 * is_5[k]
                  + pb_z[k] * ks_9[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pa_z, pb_x, pb_y, pb_z, is_6, \
                         is_10, ip_18, ip_19, ks_10, ks_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * is_10[k]
                  + pb_x[k] * ks_10[k];

        t_31[k] = f_4 * is_6[k]
                  + pb_y[k] * ks_10[k];

        t_32[k] = pb_z[k] * ks_10[k];

        t_33[k] = pa_z[k] * ip_18[k];

        t_34[k] = pa_z[k] * ip_19[k];

        t_35[k] = f_1 * is_6[k]
                  + pb_z[k] * ks_11[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_y, pb_x, pb_y, pb_z, is_7, is_8, \
                         is_9, is_12, ip_27, ks_12, ks_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_5 * is_12[k]
                  + pb_x[k] * ks_12[k];

        t_37[k] = f_3 * is_8[k]
                  + pb_y[k] * ks_12[k];

        t_38[k] = f_3 * is_7[k]
                  + pb_z[k] * ks_12[k];

        t_39[k] = pa_y[k] * ip_27[k];

        t_40[k] = f_1 * is_9[k]
                  + pb_y[k] * ks_13[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_y, pb_x, pb_y, pb_z, is_9, is_14, \
                         is_15, ip_29, ks_14, ks_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pa_y[k] * ip_29[k];

        t_42[k] = f_5 * is_14[k]
                  + pb_x[k] * ks_14[k];

        t_43[k] = pb_y[k] * ks_14[k];

        t_44[k] = f_4 * is_9[k]
                  + pb_z[k] * ks_14[k];

        t_45[k] = f_3 * is_15[k]
                  + pb_x[k] * ks_15[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_z, pb_y, pb_z, is_10, ip_30, ip_31, \
                         ks_15, ks_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_2 * is_10[k]
                  + pb_y[k] * ks_15[k];

        t_47[k] = pb_z[k] * ks_15[k];

        t_48[k] = pa_z[k] * ip_30[k];

        t_49[k] = pa_z[k] * ip_31[k];

        t_50[k] = f_1 * is_10[k]
                  + pb_z[k] * ks_16[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, pb_x, pb_y, pb_z, is_11, is_12, \
                         is_13, is_17, is_18, ks_17, ks_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * is_17[k]
                  + pb_x[k] * ks_17[k];

        t_52[k] = f_5 * is_12[k]
                  + pb_y[k] * ks_17[k];

        t_53[k] = f_3 * is_11[k]
                  + pb_z[k] * ks_17[k];

        t_54[k] = f_3 * is_18[k]
                  + pb_x[k] * ks_18[k];

        t_55[k] = f_3 * is_13[k]
                  + pb_y[k] * ks_18[k];

        t_56[k] = f_5 * is_12[k]
                  + pb_z[k] * ks_18[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, pb_z, is_14, \
                         is_20, ip_42, ip_44, ks_19, ks_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_y[k] * ip_42[k];

        t_58[k] = f_1 * is_14[k]
                  + pb_y[k] * ks_19[k];

        t_59[k] = pa_y[k] * ip_44[k];

        t_60[k] = f_3 * is_20[k]
                  + pb_x[k] * ks_20[k];

        t_61[k] = pb_y[k] * ks_20[k];

        t_62[k] = f_2 * is_14[k]
                  + pb_z[k] * ks_20[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, t_68, pa_x, pa_z, pb_x, pb_z, is_21, \
                         ip_45, ip_64, ip_67, ip_68, ks_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_1 * is_21[k]
                  + pb_x[k] * ks_21[k];

        t_64[k] = pa_x[k] * ip_64[k];

        t_65[k] = pb_z[k] * ks_21[k];

        t_66[k] = pa_z[k] * ip_45[k];

        t_67[k] = pa_x[k] * ip_67[k];

        t_68[k] = pa_x[k] * ip_68[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, t_74, pa_x, pb_x, is_23, is_24, ip_70, \
                         ip_71, ip_73, ip_74, ks_23, ks_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * is_23[k]
                  + pb_x[k] * ks_23[k];

        t_70[k] = pa_x[k] * ip_70[k];

        t_71[k] = pa_x[k] * ip_71[k];

        t_72[k] = f_1 * is_24[k]
                  + pb_x[k] * ks_24[k];

        t_73[k] = pa_x[k] * ip_73[k];

        t_74[k] = pa_x[k] * ip_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, t_80, pa_x, pa_y, pb_x, is_25, ip_60, \
                         ip_76, ip_77, ip_79, ip_80, ks_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * is_25[k]
                  + pb_x[k] * ks_25[k];

        t_76[k] = pa_x[k] * ip_76[k];

        t_77[k] = pa_x[k] * ip_77[k];

        t_78[k] = pa_y[k] * ip_60[k];

        t_79[k] = pa_x[k] * ip_79[k];

        t_80[k] = pa_x[k] * ip_80[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, pa_x, pb_x, pb_y, pb_z, is_21, \
                         is_27, ip_83, ks_27, ks_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_1 * is_27[k]
                  + pb_x[k] * ks_27[k];

        t_82[k] = pb_y[k] * ks_27[k];

        t_83[k] = pa_x[k] * ip_83[k];

        t_84[k] = pb_x[k] * ks_28[k];

        t_85[k] = f_0 * is_21[k]
                  + pb_y[k] * ks_28[k];

        t_86[k] = pb_z[k] * ks_28[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, pa_z, pb_x, pb_y, pb_z, is_21, \
                         is_22, is_23, ip_64, ks_29, ks_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pb_x[k] * ks_29[k];

        t_88[k] = pa_z[k] * ip_64[k];

        t_89[k] = f_1 * is_21[k]
                  + pb_z[k] * ks_29[k];

        t_90[k] = pb_x[k] * ks_30[k];

        t_91[k] = f_2 * is_23[k]
                  + pb_y[k] * ks_30[k];

        t_92[k] = f_3 * is_22[k]
                  + pb_z[k] * ks_30[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, pb_x, pb_y, pb_z, is_23, \
                         is_24, is_25, ks_31, ks_32, ks_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pb_x[k] * ks_31[k];

        t_94[k] = f_4 * is_24[k]
                  + pb_y[k] * ks_31[k];

        t_95[k] = f_5 * is_23[k]
                  + pb_z[k] * ks_31[k];

        t_96[k] = pb_x[k] * ks_32[k];

        t_97[k] = f_5 * is_25[k]
                  + pb_y[k] * ks_32[k];

        t_98[k] = f_4 * is_24[k]
                  + pb_z[k] * ks_32[k];

        t_99[k] = pb_x[k] * ks_33[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pa_y, pb_x, pb_y, pb_z, is_25, \
                         is_26, is_27, ip_83, ks_33, ks_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_3 * is_26[k]
                   + pb_y[k] * ks_33[k];

        t_101[k] = f_2 * is_25[k]
                   + pb_z[k] * ks_33[k];

        t_102[k] = pb_x[k] * ks_34[k];

        t_103[k] = f_1 * is_27[k]
                   + pb_y[k] * ks_34[k];

        t_104[k] = pa_y[k] * ip_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pb_x, pb_y, pb_z, is_27, \
                         ks_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pb_x[k] * ks_35[k];

        t_106[k] = pb_y[k] * ks_35[k];

        t_107[k] = f_0 * is_27[k]
                   + pb_z[k] * ks_35[k];
    }
}

}  // namespace simdt2ceri
