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


#include "SimdTransferKP.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hrr_kp(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t ks, const size_t ls, const size_t nmax) -> void
{
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *ks_0 = buffer.data(ks + 0);
    const auto *ks_1 = buffer.data(ks + 1);
    const auto *ks_2 = buffer.data(ks + 2);
    const auto *ks_3 = buffer.data(ks + 3);
    const auto *ks_4 = buffer.data(ks + 4);
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
    const auto *ks_22 = buffer.data(ks + 22);
    const auto *ks_23 = buffer.data(ks + 23);
    const auto *ks_24 = buffer.data(ks + 24);
    const auto *ks_25 = buffer.data(ks + 25);
    const auto *ks_26 = buffer.data(ks + 26);
    const auto *ks_27 = buffer.data(ks + 27);
    const auto *ks_28 = buffer.data(ks + 28);
    const auto *ks_29 = buffer.data(ks + 29);
    const auto *ks_30 = buffer.data(ks + 30);
    const auto *ks_31 = buffer.data(ks + 31);
    const auto *ks_32 = buffer.data(ks + 32);
    const auto *ks_33 = buffer.data(ks + 33);
    const auto *ks_34 = buffer.data(ks + 34);
    const auto *ks_35 = buffer.data(ks + 35);

    const auto *ls_0 = buffer.data(ls + 0);
    const auto *ls_1 = buffer.data(ls + 1);
    const auto *ls_2 = buffer.data(ls + 2);
    const auto *ls_3 = buffer.data(ls + 3);
    const auto *ls_4 = buffer.data(ls + 4);
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
    const auto *ls_29 = buffer.data(ls + 29);
    const auto *ls_30 = buffer.data(ls + 30);
    const auto *ls_31 = buffer.data(ls + 31);
    const auto *ls_32 = buffer.data(ls + 32);
    const auto *ls_33 = buffer.data(ls + 33);
    const auto *ls_34 = buffer.data(ls + 34);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, ks_0, ks_1, ls_0, \
                         ls_1, ls_2, ls_3, ls_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * ks_0[k]
                 + ls_0[k];

        t_1[k] = ab_y[k] * ks_0[k]
                 + ls_1[k];

        t_2[k] = ab_z[k] * ks_0[k]
                 + ls_2[k];

        t_3[k] = ab_x[k] * ks_1[k]
                 + ls_1[k];

        t_4[k] = ab_y[k] * ks_1[k]
                 + ls_3[k];

        t_5[k] = ab_z[k] * ks_1[k]
                 + ls_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, ks_2, ks_3, ls_2, ls_3, \
                         ls_4, ls_5, ls_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_6[k] = ab_x[k] * ks_2[k]
                 + ls_2[k];

        t_7[k] = ab_y[k] * ks_2[k]
                 + ls_4[k];

        t_8[k] = ab_z[k] * ks_2[k]
                 + ls_5[k];

        t_9[k] = ab_x[k] * ks_3[k]
                 + ls_3[k];

        t_10[k] = ab_y[k] * ks_3[k]
                  + ls_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, ks_3, ks_4, \
                         ks_5, ls_4, ls_5, ls_7, ls_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_11[k] = ab_z[k] * ks_3[k]
                  + ls_7[k];

        t_12[k] = ab_x[k] * ks_4[k]
                  + ls_4[k];

        t_13[k] = ab_y[k] * ks_4[k]
                  + ls_7[k];

        t_14[k] = ab_z[k] * ks_4[k]
                  + ls_8[k];

        t_15[k] = ab_x[k] * ks_5[k]
                  + ls_5[k];

        t_16[k] = ab_y[k] * ks_5[k]
                  + ls_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, ks_5, ks_6, ks_7, \
                         ls_6, ls_7, ls_9, ls_10, ls_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_17[k] = ab_z[k] * ks_5[k]
                  + ls_9[k];

        t_18[k] = ab_x[k] * ks_6[k]
                  + ls_6[k];

        t_19[k] = ab_y[k] * ks_6[k]
                  + ls_10[k];

        t_20[k] = ab_z[k] * ks_6[k]
                  + ls_11[k];

        t_21[k] = ab_x[k] * ks_7[k]
                  + ls_7[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, ks_7, ks_8, ls_8, \
                         ls_11, ls_12, ls_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_22[k] = ab_y[k] * ks_7[k]
                  + ls_11[k];

        t_23[k] = ab_z[k] * ks_7[k]
                  + ls_12[k];

        t_24[k] = ab_x[k] * ks_8[k]
                  + ls_8[k];

        t_25[k] = ab_y[k] * ks_8[k]
                  + ls_12[k];

        t_26[k] = ab_z[k] * ks_8[k]
                  + ls_13[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, ks_9, ks_10, ls_9, \
                         ls_10, ls_13, ls_14, ls_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_27[k] = ab_x[k] * ks_9[k]
                  + ls_9[k];

        t_28[k] = ab_y[k] * ks_9[k]
                  + ls_13[k];

        t_29[k] = ab_z[k] * ks_9[k]
                  + ls_14[k];

        t_30[k] = ab_x[k] * ks_10[k]
                  + ls_10[k];

        t_31[k] = ab_y[k] * ks_10[k]
                  + ls_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ab_x, ab_y, ab_z, ks_10, ks_11, \
                         ks_12, ls_11, ls_12, ls_16, ls_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_32[k] = ab_z[k] * ks_10[k]
                  + ls_16[k];

        t_33[k] = ab_x[k] * ks_11[k]
                  + ls_11[k];

        t_34[k] = ab_y[k] * ks_11[k]
                  + ls_16[k];

        t_35[k] = ab_z[k] * ks_11[k]
                  + ls_17[k];

        t_36[k] = ab_x[k] * ks_12[k]
                  + ls_12[k];

        t_37[k] = ab_y[k] * ks_12[k]
                  + ls_17[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, ks_12, ks_13, \
                         ks_14, ls_13, ls_14, ls_18, ls_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_38[k] = ab_z[k] * ks_12[k]
                  + ls_18[k];

        t_39[k] = ab_x[k] * ks_13[k]
                  + ls_13[k];

        t_40[k] = ab_y[k] * ks_13[k]
                  + ls_18[k];

        t_41[k] = ab_z[k] * ks_13[k]
                  + ls_19[k];

        t_42[k] = ab_x[k] * ks_14[k]
                  + ls_14[k];

        t_43[k] = ab_y[k] * ks_14[k]
                  + ls_19[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, ks_14, ks_15, ks_16, \
                         ls_15, ls_16, ls_20, ls_21, ls_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_44[k] = ab_z[k] * ks_14[k]
                  + ls_20[k];

        t_45[k] = ab_x[k] * ks_15[k]
                  + ls_15[k];

        t_46[k] = ab_y[k] * ks_15[k]
                  + ls_21[k];

        t_47[k] = ab_z[k] * ks_15[k]
                  + ls_22[k];

        t_48[k] = ab_x[k] * ks_16[k]
                  + ls_16[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, ks_16, ks_17, ls_17, \
                         ls_22, ls_23, ls_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_49[k] = ab_y[k] * ks_16[k]
                  + ls_22[k];

        t_50[k] = ab_z[k] * ks_16[k]
                  + ls_23[k];

        t_51[k] = ab_x[k] * ks_17[k]
                  + ls_17[k];

        t_52[k] = ab_y[k] * ks_17[k]
                  + ls_23[k];

        t_53[k] = ab_z[k] * ks_17[k]
                  + ls_24[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, ks_18, ks_19, \
                         ls_18, ls_19, ls_24, ls_25, ls_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_54[k] = ab_x[k] * ks_18[k]
                  + ls_18[k];

        t_55[k] = ab_y[k] * ks_18[k]
                  + ls_24[k];

        t_56[k] = ab_z[k] * ks_18[k]
                  + ls_25[k];

        t_57[k] = ab_x[k] * ks_19[k]
                  + ls_19[k];

        t_58[k] = ab_y[k] * ks_19[k]
                  + ls_25[k];

        t_59[k] = ab_z[k] * ks_19[k]
                  + ls_26[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, ab_y, ab_z, ks_20, ks_21, ls_20, \
                         ls_21, ls_26, ls_27, ls_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = ab_x[k] * ks_20[k]
                  + ls_20[k];

        t_61[k] = ab_y[k] * ks_20[k]
                  + ls_26[k];

        t_62[k] = ab_z[k] * ks_20[k]
                  + ls_27[k];

        t_63[k] = ab_x[k] * ks_21[k]
                  + ls_21[k];

        t_64[k] = ab_y[k] * ks_21[k]
                  + ls_28[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, t_70, ab_x, ab_y, ab_z, ks_21, ks_22, \
                         ks_23, ls_22, ls_23, ls_29, ls_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = ab_z[k] * ks_21[k]
                  + ls_29[k];

        t_66[k] = ab_x[k] * ks_22[k]
                  + ls_22[k];

        t_67[k] = ab_y[k] * ks_22[k]
                  + ls_29[k];

        t_68[k] = ab_z[k] * ks_22[k]
                  + ls_30[k];

        t_69[k] = ab_x[k] * ks_23[k]
                  + ls_23[k];

        t_70[k] = ab_y[k] * ks_23[k]
                  + ls_30[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, ab_z, ks_23, ks_24, \
                         ks_25, ls_24, ls_25, ls_31, ls_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_71[k] = ab_z[k] * ks_23[k]
                  + ls_31[k];

        t_72[k] = ab_x[k] * ks_24[k]
                  + ls_24[k];

        t_73[k] = ab_y[k] * ks_24[k]
                  + ls_31[k];

        t_74[k] = ab_z[k] * ks_24[k]
                  + ls_32[k];

        t_75[k] = ab_x[k] * ks_25[k]
                  + ls_25[k];

        t_76[k] = ab_y[k] * ks_25[k]
                  + ls_32[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, ab_x, ab_y, ab_z, ks_25, ks_26, \
                         ks_27, ls_26, ls_27, ls_33, ls_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_77[k] = ab_z[k] * ks_25[k]
                  + ls_33[k];

        t_78[k] = ab_x[k] * ks_26[k]
                  + ls_26[k];

        t_79[k] = ab_y[k] * ks_26[k]
                  + ls_33[k];

        t_80[k] = ab_z[k] * ks_26[k]
                  + ls_34[k];

        t_81[k] = ab_x[k] * ks_27[k]
                  + ls_27[k];

        t_82[k] = ab_y[k] * ks_27[k]
                  + ls_34[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, ab_x, ab_y, ab_z, ks_27, ks_28, ks_29, \
                         ls_28, ls_29, ls_35, ls_36, ls_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_83[k] = ab_z[k] * ks_27[k]
                  + ls_35[k];

        t_84[k] = ab_x[k] * ks_28[k]
                  + ls_28[k];

        t_85[k] = ab_y[k] * ks_28[k]
                  + ls_36[k];

        t_86[k] = ab_z[k] * ks_28[k]
                  + ls_37[k];

        t_87[k] = ab_x[k] * ks_29[k]
                  + ls_29[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, ab_x, ab_y, ab_z, ks_29, ks_30, ls_30, \
                         ls_37, ls_38, ls_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_88[k] = ab_y[k] * ks_29[k]
                  + ls_37[k];

        t_89[k] = ab_z[k] * ks_29[k]
                  + ls_38[k];

        t_90[k] = ab_x[k] * ks_30[k]
                  + ls_30[k];

        t_91[k] = ab_y[k] * ks_30[k]
                  + ls_38[k];

        t_92[k] = ab_z[k] * ks_30[k]
                  + ls_39[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, ab_x, ab_y, ab_z, ks_31, ks_32, \
                         ls_31, ls_32, ls_39, ls_40, ls_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_93[k] = ab_x[k] * ks_31[k]
                  + ls_31[k];

        t_94[k] = ab_y[k] * ks_31[k]
                  + ls_39[k];

        t_95[k] = ab_z[k] * ks_31[k]
                  + ls_40[k];

        t_96[k] = ab_x[k] * ks_32[k]
                  + ls_32[k];

        t_97[k] = ab_y[k] * ks_32[k]
                  + ls_40[k];

        t_98[k] = ab_z[k] * ks_32[k]
                  + ls_41[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, ab_x, ab_y, ab_z, ks_33, \
                         ks_34, ls_33, ls_34, ls_41, ls_42, ls_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_99[k] = ab_x[k] * ks_33[k]
                  + ls_33[k];

        t_100[k] = ab_y[k] * ks_33[k]
                   + ls_41[k];

        t_101[k] = ab_z[k] * ks_33[k]
                   + ls_42[k];

        t_102[k] = ab_x[k] * ks_34[k]
                   + ls_34[k];

        t_103[k] = ab_y[k] * ks_34[k]
                   + ls_42[k];

        t_104[k] = ab_z[k] * ks_34[k]
                   + ls_43[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, ab_x, ab_y, ab_z, ks_35, ls_35, ls_43, \
                         ls_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = ab_x[k] * ks_35[k]
                   + ls_35[k];

        t_106[k] = ab_y[k] * ks_35[k]
                   + ls_43[k];

        t_107[k] = ab_z[k] * ks_35[k]
                   + ls_44[k];
    }
}

}  // namespace simdovl
