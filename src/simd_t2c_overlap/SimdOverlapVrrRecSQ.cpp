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


#include "SimdOverlapVrrRecSQ.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_sq_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                          const size_t sn, const size_t so, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / p;
    const auto f_1 = 4.5 / p;
    const auto f_2 = 4.0 / p;
    const auto f_3 = 3.5 / p;
    const auto f_4 = 3.0 / p;
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 0.5 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sn_0 = buffer.data(sn + 0);
    const auto *sn_3 = buffer.data(sn + 3);
    const auto *sn_5 = buffer.data(sn + 5);
    const auto *sn_6 = buffer.data(sn + 6);
    const auto *sn_9 = buffer.data(sn + 9);
    const auto *sn_10 = buffer.data(sn + 10);
    const auto *sn_12 = buffer.data(sn + 12);
    const auto *sn_14 = buffer.data(sn + 14);
    const auto *sn_15 = buffer.data(sn + 15);
    const auto *sn_17 = buffer.data(sn + 17);
    const auto *sn_18 = buffer.data(sn + 18);
    const auto *sn_20 = buffer.data(sn + 20);
    const auto *sn_21 = buffer.data(sn + 21);
    const auto *sn_23 = buffer.data(sn + 23);
    const auto *sn_24 = buffer.data(sn + 24);
    const auto *sn_25 = buffer.data(sn + 25);
    const auto *sn_27 = buffer.data(sn + 27);
    const auto *sn_28 = buffer.data(sn + 28);
    const auto *sn_30 = buffer.data(sn + 30);
    const auto *sn_31 = buffer.data(sn + 31);
    const auto *sn_32 = buffer.data(sn + 32);
    const auto *sn_33 = buffer.data(sn + 33);
    const auto *sn_35 = buffer.data(sn + 35);
    const auto *sn_36 = buffer.data(sn + 36);
    const auto *sn_38 = buffer.data(sn + 38);
    const auto *sn_39 = buffer.data(sn + 39);
    const auto *sn_40 = buffer.data(sn + 40);
    const auto *sn_41 = buffer.data(sn + 41);
    const auto *sn_42 = buffer.data(sn + 42);
    const auto *sn_44 = buffer.data(sn + 44);
    const auto *sn_45 = buffer.data(sn + 45);
    const auto *sn_47 = buffer.data(sn + 47);
    const auto *sn_48 = buffer.data(sn + 48);
    const auto *sn_49 = buffer.data(sn + 49);
    const auto *sn_50 = buffer.data(sn + 50);
    const auto *sn_51 = buffer.data(sn + 51);
    const auto *sn_52 = buffer.data(sn + 52);
    const auto *sn_54 = buffer.data(sn + 54);
    const auto *sn_55 = buffer.data(sn + 55);
    const auto *sn_57 = buffer.data(sn + 57);
    const auto *sn_58 = buffer.data(sn + 58);
    const auto *sn_59 = buffer.data(sn + 59);
    const auto *sn_60 = buffer.data(sn + 60);
    const auto *sn_61 = buffer.data(sn + 61);
    const auto *sn_62 = buffer.data(sn + 62);
    const auto *sn_63 = buffer.data(sn + 63);
    const auto *sn_64 = buffer.data(sn + 64);
    const auto *sn_65 = buffer.data(sn + 65);

    const auto *so_0 = buffer.data(so + 0);
    const auto *so_2 = buffer.data(so + 2);
    const auto *so_3 = buffer.data(so + 3);
    const auto *so_5 = buffer.data(so + 5);
    const auto *so_6 = buffer.data(so + 6);
    const auto *so_9 = buffer.data(so + 9);
    const auto *so_10 = buffer.data(so + 10);
    const auto *so_12 = buffer.data(so + 12);
    const auto *so_14 = buffer.data(so + 14);
    const auto *so_15 = buffer.data(so + 15);
    const auto *so_17 = buffer.data(so + 17);
    const auto *so_18 = buffer.data(so + 18);
    const auto *so_20 = buffer.data(so + 20);
    const auto *so_21 = buffer.data(so + 21);
    const auto *so_23 = buffer.data(so + 23);
    const auto *so_24 = buffer.data(so + 24);
    const auto *so_25 = buffer.data(so + 25);
    const auto *so_27 = buffer.data(so + 27);
    const auto *so_28 = buffer.data(so + 28);
    const auto *so_30 = buffer.data(so + 30);
    const auto *so_31 = buffer.data(so + 31);
    const auto *so_32 = buffer.data(so + 32);
    const auto *so_33 = buffer.data(so + 33);
    const auto *so_35 = buffer.data(so + 35);
    const auto *so_36 = buffer.data(so + 36);
    const auto *so_38 = buffer.data(so + 38);
    const auto *so_39 = buffer.data(so + 39);
    const auto *so_40 = buffer.data(so + 40);
    const auto *so_41 = buffer.data(so + 41);
    const auto *so_42 = buffer.data(so + 42);
    const auto *so_44 = buffer.data(so + 44);
    const auto *so_45 = buffer.data(so + 45);
    const auto *so_47 = buffer.data(so + 47);
    const auto *so_48 = buffer.data(so + 48);
    const auto *so_49 = buffer.data(so + 49);
    const auto *so_50 = buffer.data(so + 50);
    const auto *so_51 = buffer.data(so + 51);
    const auto *so_52 = buffer.data(so + 52);
    const auto *so_54 = buffer.data(so + 54);
    const auto *so_55 = buffer.data(so + 55);
    const auto *so_57 = buffer.data(so + 57);
    const auto *so_58 = buffer.data(so + 58);
    const auto *so_59 = buffer.data(so + 59);
    const auto *so_60 = buffer.data(so + 60);
    const auto *so_61 = buffer.data(so + 61);
    const auto *so_62 = buffer.data(so + 62);
    const auto *so_63 = buffer.data(so + 63);
    const auto *so_65 = buffer.data(so + 65);
    const auto *so_66 = buffer.data(so + 66);
    const auto *so_67 = buffer.data(so + 67);
    const auto *so_68 = buffer.data(so + 68);
    const auto *so_69 = buffer.data(so + 69);
    const auto *so_70 = buffer.data(so + 70);
    const auto *so_71 = buffer.data(so + 71);
    const auto *so_72 = buffer.data(so + 72);
    const auto *so_73 = buffer.data(so + 73);
    const auto *so_74 = buffer.data(so + 74);
    const auto *so_75 = buffer.data(so + 75);
    const auto *so_76 = buffer.data(so + 76);
    const auto *so_77 = buffer.data(so + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, sn_0, sn_3, sn_5, \
                         so_0, so_2, so_3, so_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sn_0[k]
                 + pb_x[k] * so_0[k];

        t_1[k] = pb_y[k] * so_0[k];

        t_2[k] = pb_z[k] * so_0[k];

        t_3[k] = f_1 * sn_3[k]
                 + pb_x[k] * so_3[k];

        t_4[k] = pb_y[k] * so_2[k];

        t_5[k] = f_1 * sn_5[k]
                 + pb_x[k] * so_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, sn_6, sn_9, sn_10, so_3, \
                         so_5, so_6, so_9, so_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * sn_6[k]
                 + pb_x[k] * so_6[k];

        t_7[k] = pb_z[k] * so_3[k];

        t_8[k] = pb_y[k] * so_5[k];

        t_9[k] = f_2 * sn_9[k]
                 + pb_x[k] * so_9[k];

        t_10[k] = f_3 * sn_10[k]
                  + pb_x[k] * so_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pb_x, pb_y, pb_z, sn_12, sn_14, sn_15, \
                         so_6, so_9, so_12, so_14, so_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * so_6[k];

        t_12[k] = f_3 * sn_12[k]
                  + pb_x[k] * so_12[k];

        t_13[k] = pb_y[k] * so_9[k];

        t_14[k] = f_3 * sn_14[k]
                  + pb_x[k] * so_14[k];

        t_15[k] = f_4 * sn_15[k]
                  + pb_x[k] * so_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pb_x, pb_y, pb_z, sn_17, sn_18, sn_20, \
                         so_10, so_14, so_17, so_18, so_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_z[k] * so_10[k];

        t_17[k] = f_4 * sn_17[k]
                  + pb_x[k] * so_17[k];

        t_18[k] = f_4 * sn_18[k]
                  + pb_x[k] * so_18[k];

        t_19[k] = pb_y[k] * so_14[k];

        t_20[k] = f_4 * sn_20[k]
                  + pb_x[k] * so_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_z, sn_21, sn_23, sn_24, sn_25, \
                         so_15, so_21, so_23, so_24, so_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * sn_21[k]
                  + pb_x[k] * so_21[k];

        t_22[k] = pb_z[k] * so_15[k];

        t_23[k] = f_5 * sn_23[k]
                  + pb_x[k] * so_23[k];

        t_24[k] = f_5 * sn_24[k]
                  + pb_x[k] * so_24[k];

        t_25[k] = f_5 * sn_25[k]
                  + pb_x[k] * so_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pb_x, pb_y, pb_z, sn_27, sn_28, sn_30, \
                         so_20, so_21, so_27, so_28, so_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * so_20[k];

        t_27[k] = f_5 * sn_27[k]
                  + pb_x[k] * so_27[k];

        t_28[k] = f_6 * sn_28[k]
                  + pb_x[k] * so_28[k];

        t_29[k] = pb_z[k] * so_21[k];

        t_30[k] = f_6 * sn_30[k]
                  + pb_x[k] * so_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pb_x, pb_y, sn_31, sn_32, sn_33, sn_35, \
                         so_27, so_31, so_32, so_33, so_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * sn_31[k]
                  + pb_x[k] * so_31[k];

        t_32[k] = f_6 * sn_32[k]
                  + pb_x[k] * so_32[k];

        t_33[k] = f_6 * sn_33[k]
                  + pb_x[k] * so_33[k];

        t_34[k] = pb_y[k] * so_27[k];

        t_35[k] = f_6 * sn_35[k]
                  + pb_x[k] * so_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pb_x, pb_z, sn_36, sn_38, sn_39, sn_40, \
                         so_28, so_36, so_38, so_39, so_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_7 * sn_36[k]
                  + pb_x[k] * so_36[k];

        t_37[k] = pb_z[k] * so_28[k];

        t_38[k] = f_7 * sn_38[k]
                  + pb_x[k] * so_38[k];

        t_39[k] = f_7 * sn_39[k]
                  + pb_x[k] * so_39[k];

        t_40[k] = f_7 * sn_40[k]
                  + pb_x[k] * so_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pb_x, pb_y, sn_41, sn_42, sn_44, sn_45, \
                         so_35, so_41, so_42, so_44, so_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * sn_41[k]
                  + pb_x[k] * so_41[k];

        t_42[k] = f_7 * sn_42[k]
                  + pb_x[k] * so_42[k];

        t_43[k] = pb_y[k] * so_35[k];

        t_44[k] = f_7 * sn_44[k]
                  + pb_x[k] * so_44[k];

        t_45[k] = f_8 * sn_45[k]
                  + pb_x[k] * so_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pb_x, pb_z, sn_47, sn_48, sn_49, sn_50, \
                         so_36, so_47, so_48, so_49, so_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_z[k] * so_36[k];

        t_47[k] = f_8 * sn_47[k]
                  + pb_x[k] * so_47[k];

        t_48[k] = f_8 * sn_48[k]
                  + pb_x[k] * so_48[k];

        t_49[k] = f_8 * sn_49[k]
                  + pb_x[k] * so_49[k];

        t_50[k] = f_8 * sn_50[k]
                  + pb_x[k] * so_50[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pb_x, pb_y, sn_51, sn_52, sn_54, sn_55, \
                         so_44, so_51, so_52, so_54, so_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_8 * sn_51[k]
                  + pb_x[k] * so_51[k];

        t_52[k] = f_8 * sn_52[k]
                  + pb_x[k] * so_52[k];

        t_53[k] = pb_y[k] * so_44[k];

        t_54[k] = f_8 * sn_54[k]
                  + pb_x[k] * so_54[k];

        t_55[k] = f_9 * sn_55[k]
                  + pb_x[k] * so_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pb_x, pb_z, sn_57, sn_58, sn_59, sn_60, \
                         so_45, so_57, so_58, so_59, so_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * so_45[k];

        t_57[k] = f_9 * sn_57[k]
                  + pb_x[k] * so_57[k];

        t_58[k] = f_9 * sn_58[k]
                  + pb_x[k] * so_58[k];

        t_59[k] = f_9 * sn_59[k]
                  + pb_x[k] * so_59[k];

        t_60[k] = f_9 * sn_60[k]
                  + pb_x[k] * so_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pb_x, pb_y, sn_61, sn_62, sn_63, sn_65, \
                         so_54, so_61, so_62, so_63, so_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_9 * sn_61[k]
                  + pb_x[k] * so_61[k];

        t_62[k] = f_9 * sn_62[k]
                  + pb_x[k] * so_62[k];

        t_63[k] = f_9 * sn_63[k]
                  + pb_x[k] * so_63[k];

        t_64[k] = pb_y[k] * so_54[k];

        t_65[k] = f_9 * sn_65[k]
                  + pb_x[k] * so_65[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, t_72, pb_x, so_66, so_67, so_68, \
                         so_69, so_70, so_71, so_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_x[k] * so_66[k];

        t_67[k] = pb_x[k] * so_67[k];

        t_68[k] = pb_x[k] * so_68[k];

        t_69[k] = pb_x[k] * so_69[k];

        t_70[k] = pb_x[k] * so_70[k];

        t_71[k] = pb_x[k] * so_71[k];

        t_72[k] = pb_x[k] * so_72[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, t_78, pb_x, pb_y, sn_55, so_66, so_73, \
                         so_74, so_75, so_76, so_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_x[k] * so_73[k];

        t_74[k] = pb_x[k] * so_74[k];

        t_75[k] = pb_x[k] * so_75[k];

        t_76[k] = pb_x[k] * so_76[k];

        t_77[k] = pb_x[k] * so_77[k];

        t_78[k] = f_0 * sn_55[k]
                  + pb_y[k] * so_66[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_y, pb_z, sn_57, sn_58, sn_59, sn_60, \
                         so_66, so_68, so_69, so_70, so_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_z[k] * so_66[k];

        t_80[k] = f_1 * sn_57[k]
                  + pb_y[k] * so_68[k];

        t_81[k] = f_2 * sn_58[k]
                  + pb_y[k] * so_69[k];

        t_82[k] = f_3 * sn_59[k]
                  + pb_y[k] * so_70[k];

        t_83[k] = f_4 * sn_60[k]
                  + pb_y[k] * so_71[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pb_y, sn_61, sn_62, sn_63, sn_64, \
                         sn_65, so_72, so_73, so_74, so_75, so_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_5 * sn_61[k]
                  + pb_y[k] * so_72[k];

        t_85[k] = f_6 * sn_62[k]
                  + pb_y[k] * so_73[k];

        t_86[k] = f_7 * sn_63[k]
                  + pb_y[k] * so_74[k];

        t_87[k] = f_8 * sn_64[k]
                  + pb_y[k] * so_75[k];

        t_88[k] = f_9 * sn_65[k]
                  + pb_y[k] * so_76[k];
    }

#pragma omp simd aligned(t_89, t_90, pb_y, pb_z, sn_65, so_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pb_y[k] * so_77[k];

        t_90[k] = f_0 * sn_65[k]
                  + pb_z[k] * so_77[k];
    }
}

}  // namespace simdovl
