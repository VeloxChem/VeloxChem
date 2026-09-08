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


#include "SimdElectronRepulsionVrrRecIP.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_ip_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hs, const size_t hp,
                                     const size_t is, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_16 = buffer.data(hs + 16);
    const auto *hs_17 = buffer.data(hs + 17);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
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
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);
    const auto *is_23 = buffer.data(is + 23);
    const auto *is_24 = buffer.data(is + 24);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pa_y, pa_z, pb_x, pb_y, pb_z, \
                         hs_0, hp_0, is_0, is_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = pa_y[k] * hp_0[k];

        t_4[k] = f_1 * hs_0[k]
                 + pb_y[k] * is_1[k];

        t_5[k] = pb_z[k] * is_1[k];

        t_6[k] = pa_z[k] * hp_0[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pa_y, pb_x, pb_y, pb_z, hs_0, hs_1, \
                         hs_3, hp_2, is_2, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_y[k] * is_2[k];

        t_8[k] = f_1 * hs_0[k]
                 + pb_z[k] * is_2[k];

        t_9[k] = f_2 * hs_3[k]
                 + pb_x[k] * is_3[k];

        t_10[k] = f_3 * hs_1[k]
                  + pb_y[k] * is_3[k];

        t_11[k] = pb_z[k] * is_3[k];

        t_12[k] = pa_y[k] * hp_2[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_y, pa_z, pb_x, pb_y, pb_z, hs_2, \
                         hs_4, hp_1, hp_3, is_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_z[k] * hp_1[k];

        t_14[k] = pa_y[k] * hp_3[k];

        t_15[k] = f_2 * hs_4[k]
                  + pb_x[k] * is_4[k];

        t_16[k] = pb_y[k] * is_4[k];

        t_17[k] = f_3 * hs_2[k]
                  + pb_z[k] * is_4[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, pa_z, pb_x, pb_y, pb_z, hs_3, \
                         hs_5, hp_4, hp_5, is_5, is_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_4 * hs_5[k]
                  + pb_x[k] * is_5[k];

        t_19[k] = f_4 * hs_3[k]
                  + pb_y[k] * is_5[k];

        t_20[k] = pb_z[k] * is_5[k];

        t_21[k] = pa_z[k] * hp_4[k];

        t_22[k] = pa_z[k] * hp_5[k];

        t_23[k] = f_1 * hs_3[k]
                  + pb_z[k] * is_6[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, pa_y, pb_x, pb_y, pb_z, hs_4, \
                         hs_8, hp_6, hp_7, is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * hp_6[k];

        t_25[k] = f_1 * hs_4[k]
                  + pb_y[k] * is_7[k];

        t_26[k] = pa_y[k] * hp_7[k];

        t_27[k] = f_4 * hs_8[k]
                  + pb_x[k] * is_8[k];

        t_28[k] = pb_y[k] * is_8[k];

        t_29[k] = f_4 * hs_4[k]
                  + pb_z[k] * is_8[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pa_z, pb_x, pb_y, pb_z, hs_5, \
                         hs_9, hp_8, hp_9, is_9, is_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * hs_9[k]
                  + pb_x[k] * is_9[k];

        t_31[k] = f_2 * hs_5[k]
                  + pb_y[k] * is_9[k];

        t_32[k] = pb_z[k] * is_9[k];

        t_33[k] = pa_z[k] * hp_8[k];

        t_34[k] = pa_z[k] * hp_9[k];

        t_35[k] = f_1 * hs_5[k]
                  + pb_z[k] * is_10[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_y, pb_x, pb_y, pb_z, hs_6, hs_7, \
                         hs_8, hs_10, hp_10, is_11, is_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_3 * hs_10[k]
                  + pb_x[k] * is_11[k];

        t_37[k] = f_3 * hs_7[k]
                  + pb_y[k] * is_11[k];

        t_38[k] = f_3 * hs_6[k]
                  + pb_z[k] * is_11[k];

        t_39[k] = pa_y[k] * hp_10[k];

        t_40[k] = f_1 * hs_8[k]
                  + pb_y[k] * is_12[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_y, pb_x, pb_y, pb_z, hs_8, hs_11, \
                         hs_12, hp_11, is_13, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pa_y[k] * hp_11[k];

        t_42[k] = f_3 * hs_11[k]
                  + pb_x[k] * is_13[k];

        t_43[k] = pb_y[k] * is_13[k];

        t_44[k] = f_2 * hs_8[k]
                  + pb_z[k] * is_13[k];

        t_45[k] = f_1 * hs_12[k]
                  + pb_x[k] * is_14[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_x, pa_z, pb_z, hp_12, hp_14, hp_15, \
                         hp_16, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_x[k] * hp_14[k];

        t_47[k] = pb_z[k] * is_14[k];

        t_48[k] = pa_z[k] * hp_12[k];

        t_49[k] = pa_x[k] * hp_15[k];

        t_50[k] = pa_x[k] * hp_16[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, pa_x, pb_x, hs_14, hs_15, hp_17, \
                         hp_18, hp_19, hp_20, is_15, is_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * hs_14[k]
                  + pb_x[k] * is_15[k];

        t_52[k] = pa_x[k] * hp_17[k];

        t_53[k] = pa_x[k] * hp_18[k];

        t_54[k] = f_1 * hs_15[k]
                  + pb_x[k] * is_16[k];

        t_55[k] = pa_x[k] * hp_19[k];

        t_56[k] = pa_x[k] * hp_20[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, t_62, pa_x, pa_y, pb_x, pb_y, hs_17, \
                         hp_13, hp_21, hp_22, hp_23, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_y[k] * hp_13[k];

        t_58[k] = pa_x[k] * hp_21[k];

        t_59[k] = pa_x[k] * hp_22[k];

        t_60[k] = f_1 * hs_17[k]
                  + pb_x[k] * is_17[k];

        t_61[k] = pb_y[k] * is_17[k];

        t_62[k] = pa_x[k] * hp_23[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, t_68, t_69, pa_z, pb_x, pb_y, pb_z, \
                         hs_12, hp_14, is_18, is_19, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * is_18[k];

        t_64[k] = f_0 * hs_12[k]
                  + pb_y[k] * is_18[k];

        t_65[k] = pb_z[k] * is_18[k];

        t_66[k] = pb_x[k] * is_19[k];

        t_67[k] = pa_z[k] * hp_14[k];

        t_68[k] = f_1 * hs_12[k]
                  + pb_z[k] * is_19[k];

        t_69[k] = pb_x[k] * is_20[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, t_75, pb_x, pb_y, pb_z, hs_13, hs_14, \
                         hs_15, is_20, is_21, is_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_2 * hs_14[k]
                  + pb_y[k] * is_20[k];

        t_71[k] = f_3 * hs_13[k]
                  + pb_z[k] * is_20[k];

        t_72[k] = pb_x[k] * is_21[k];

        t_73[k] = f_4 * hs_15[k]
                  + pb_y[k] * is_21[k];

        t_74[k] = f_4 * hs_14[k]
                  + pb_z[k] * is_21[k];

        t_75[k] = pb_x[k] * is_22[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pa_y, pb_x, pb_y, pb_z, hs_15, hs_16, \
                         hs_17, hp_23, is_22, is_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_3 * hs_16[k]
                  + pb_y[k] * is_22[k];

        t_77[k] = f_2 * hs_15[k]
                  + pb_z[k] * is_22[k];

        t_78[k] = pb_x[k] * is_23[k];

        t_79[k] = f_1 * hs_17[k]
                  + pb_y[k] * is_23[k];

        t_80[k] = pa_y[k] * hp_23[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_x, pb_y, pb_z, hs_17, \
                         is_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pb_x[k] * is_24[k];

        t_82[k] = pb_y[k] * is_24[k];

        t_83[k] = f_0 * hs_17[k]
                  + pb_z[k] * is_24[k];
    }
}

auto
compute_prim_ip_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hs, const size_t hp,
                                     const size_t is, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_16 = buffer.data(hs + 16);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_29 = buffer.data(hp + 29);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
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
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_y, pa_z, pb_x, pb_y, pb_z, hs_0, \
                         hp_0, is_0, is_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = pa_y[k] * hp_0[k];

        t_4[k] = f_1 * hs_0[k]
                 + pb_y[k] * is_1[k];

        t_5[k] = pa_z[k] * hp_0[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_y, pb_x, pb_y, pb_z, hs_0, hs_1, hs_3, \
                         hp_4, is_2, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * hs_0[k]
                 + pb_z[k] * is_2[k];

        t_7[k] = f_2 * hs_3[k]
                 + pb_x[k] * is_3[k];

        t_8[k] = f_3 * hs_1[k]
                 + pb_y[k] * is_3[k];

        t_9[k] = pb_z[k] * is_3[k];

        t_10[k] = pa_y[k] * hp_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, pb_z, hs_2, hs_3, \
                         hs_4, hs_5, is_4, is_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * hs_4[k]
                  + pb_x[k] * is_4[k];

        t_12[k] = pb_y[k] * is_4[k];

        t_13[k] = f_3 * hs_2[k]
                  + pb_z[k] * is_4[k];

        t_14[k] = f_4 * hs_5[k]
                  + pb_x[k] * is_5[k];

        t_15[k] = f_4 * hs_3[k]
                  + pb_y[k] * is_5[k];

        t_16[k] = pb_z[k] * is_5[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_y, pb_x, pb_y, pb_z, hs_3, hs_4, \
                         hs_8, hp_8, is_6, is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * hs_3[k]
                  + pb_z[k] * is_6[k];

        t_18[k] = f_1 * hs_4[k]
                  + pb_y[k] * is_7[k];

        t_19[k] = pa_y[k] * hp_8[k];

        t_20[k] = f_4 * hs_8[k]
                  + pb_x[k] * is_8[k];

        t_21[k] = pb_y[k] * is_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pb_x, pb_y, pb_z, hs_4, hs_5, hs_9, \
                         is_8, is_9, is_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_4 * hs_4[k]
                  + pb_z[k] * is_8[k];

        t_23[k] = f_3 * hs_9[k]
                  + pb_x[k] * is_9[k];

        t_24[k] = f_2 * hs_5[k]
                  + pb_y[k] * is_9[k];

        t_25[k] = pb_z[k] * is_9[k];

        t_26[k] = f_1 * hs_5[k]
                  + pb_z[k] * is_10[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pb_y, pb_z, hs_6, hs_7, hs_8, hp_12, \
                         is_11, is_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_3 * hs_7[k]
                  + pb_y[k] * is_11[k];

        t_28[k] = f_3 * hs_6[k]
                  + pb_z[k] * is_11[k];

        t_29[k] = f_1 * hs_8[k]
                  + pb_y[k] * is_12[k];

        t_30[k] = pa_y[k] * hp_12[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pb_x, pb_y, pb_z, hs_8, hs_10, \
                         hs_11, hp_16, is_13, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * hs_10[k]
                  + pb_x[k] * is_13[k];

        t_32[k] = pb_y[k] * is_13[k];

        t_33[k] = f_2 * hs_8[k]
                  + pb_z[k] * is_13[k];

        t_34[k] = f_1 * hs_11[k]
                  + pb_x[k] * is_14[k];

        t_35[k] = pa_x[k] * hp_16[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, pa_x, hp_18, hp_20, hp_21, hp_23, \
                         hp_24, hp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_x[k] * hp_18[k];

        t_37[k] = pa_x[k] * hp_20[k];

        t_38[k] = pa_x[k] * hp_21[k];

        t_39[k] = pa_x[k] * hp_23[k];

        t_40[k] = pa_x[k] * hp_24[k];

        t_41[k] = pa_x[k] * hp_25[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_x, pb_x, pb_y, pb_z, hs_11, hs_16, \
                         hp_29, is_15, is_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_1 * hs_16[k]
                  + pb_x[k] * is_15[k];

        t_43[k] = pa_x[k] * hp_29[k];

        t_44[k] = pb_x[k] * is_16[k];

        t_45[k] = f_0 * hs_11[k]
                  + pb_y[k] * is_16[k];

        t_46[k] = pb_z[k] * is_16[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, pb_z, hs_11, hs_12, \
                         hs_13, hp_16, is_17, is_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_z[k] * hp_16[k];

        t_48[k] = f_1 * hs_11[k]
                  + pb_z[k] * is_17[k];

        t_49[k] = pb_x[k] * is_18[k];

        t_50[k] = f_2 * hs_13[k]
                  + pb_y[k] * is_18[k];

        t_51[k] = f_3 * hs_12[k]
                  + pb_z[k] * is_18[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, pb_x, pb_y, pb_z, hs_13, hs_14, \
                         hs_15, is_19, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_x[k] * is_19[k];

        t_53[k] = f_4 * hs_14[k]
                  + pb_y[k] * is_19[k];

        t_54[k] = f_4 * hs_13[k]
                  + pb_z[k] * is_19[k];

        t_55[k] = pb_x[k] * is_20[k];

        t_56[k] = f_3 * hs_15[k]
                  + pb_y[k] * is_20[k];

        t_57[k] = f_2 * hs_14[k]
                  + pb_z[k] * is_20[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, pb_z, hs_16, hp_29, \
                         is_21, is_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * hs_16[k]
                  + pb_y[k] * is_21[k];

        t_59[k] = pa_y[k] * hp_29[k];

        t_60[k] = pb_x[k] * is_22[k];

        t_61[k] = pb_y[k] * is_22[k];

        t_62[k] = f_0 * hs_16[k]
                  + pb_z[k] * is_22[k];
    }
}

auto
compute_prim_ip_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t hs, const size_t is, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_8 = buffer.data(hs + 8);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pb_x, pb_y, pb_z, hs_0, is_0, \
                         is_1, is_2, is_3, is_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = pb_z[k] * is_1[k];

        t_4[k] = pb_y[k] * is_2[k];

        t_5[k] = pb_z[k] * is_3[k];

        t_6[k] = pb_y[k] * is_4[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, t_13, pb_x, pb_y, pb_z, hs_5, is_5, \
                         is_6, is_7, is_8, is_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_z[k] * is_5[k];

        t_8[k] = pb_y[k] * is_6[k];

        t_9[k] = pb_x[k] * is_7[k];

        t_10[k] = f_0 * hs_5[k]
                  + pb_y[k] * is_7[k];

        t_11[k] = pb_z[k] * is_7[k];

        t_12[k] = pb_x[k] * is_8[k];

        t_13[k] = pb_x[k] * is_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_x, pb_y, pb_z, hs_8, is_10, \
                         is_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pb_x[k] * is_10[k];

        t_15[k] = pb_x[k] * is_11[k];

        t_16[k] = pb_y[k] * is_11[k];

        t_17[k] = f_0 * hs_8[k]
                  + pb_z[k] * is_11[k];
    }
}

auto
compute_prim_ip_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t hs, const size_t is, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_11 = buffer.data(hs + 11);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pb_x, pb_y, pb_z, hs_0, is_0, \
                         is_1, is_2, is_3, is_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = pb_z[k] * is_1[k];

        t_4[k] = pb_y[k] * is_2[k];

        t_5[k] = pb_z[k] * is_3[k];

        t_6[k] = pb_y[k] * is_4[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, t_13, pb_x, pb_y, pb_z, hs_7, is_5, \
                         is_6, is_7, is_8, is_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_z[k] * is_5[k];

        t_8[k] = pb_y[k] * is_6[k];

        t_9[k] = pb_x[k] * is_7[k];

        t_10[k] = f_0 * hs_7[k]
                  + pb_y[k] * is_7[k];

        t_11[k] = pb_z[k] * is_7[k];

        t_12[k] = pb_x[k] * is_8[k];

        t_13[k] = pb_x[k] * is_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_x, pb_y, pb_z, hs_11, is_10, \
                         is_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pb_x[k] * is_10[k];

        t_15[k] = pb_x[k] * is_11[k];

        t_16[k] = pb_y[k] * is_11[k];

        t_17[k] = f_0 * hs_11[k]
                  + pb_z[k] * is_11[k];
    }
}

auto
compute_prim_ip_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hs, const size_t hp,
                                     const size_t is, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_15 = buffer.data(hs + 15);

    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_21 = buffer.data(hp + 21);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
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
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pb_x, pb_y, pb_z, hs_0, hs_1, \
                         is_0, is_1, is_2, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = f_1 * hs_0[k]
                 + pb_y[k] * is_1[k];

        t_4[k] = f_1 * hs_0[k]
                 + pb_z[k] * is_2[k];

        t_5[k] = f_2 * hs_1[k]
                 + pb_y[k] * is_3[k];

        t_6[k] = pb_z[k] * is_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pa_y, pb_y, pb_z, hs_2, hs_3, hp_3, \
                         is_4, is_5, is_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_y[k] * hp_3[k];

        t_8[k] = pb_y[k] * is_4[k];

        t_9[k] = f_2 * hs_2[k]
                 + pb_z[k] * is_4[k];

        t_10[k] = f_3 * hs_3[k]
                  + pb_y[k] * is_5[k];

        t_11[k] = pb_z[k] * is_5[k];

        t_12[k] = f_1 * hs_3[k]
                  + pb_z[k] * is_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pa_y, pb_y, pb_z, hs_4, hs_5, \
                         hp_6, is_7, is_8, is_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * hp_6[k];

        t_14[k] = pb_y[k] * is_7[k];

        t_15[k] = f_3 * hs_4[k]
                  + pb_z[k] * is_7[k];

        t_16[k] = f_4 * hs_5[k]
                  + pb_y[k] * is_8[k];

        t_17[k] = pb_z[k] * is_8[k];

        t_18[k] = f_1 * hs_5[k]
                  + pb_z[k] * is_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_x, pa_y, pb_y, pb_z, hs_6, hs_7, \
                         hp_9, hp_11, is_10, is_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * hs_6[k]
                  + pb_z[k] * is_10[k];

        t_20[k] = pa_y[k] * hp_9[k];

        t_21[k] = pb_y[k] * is_11[k];

        t_22[k] = f_4 * hs_7[k]
                  + pb_z[k] * is_11[k];

        t_23[k] = pa_x[k] * hp_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, pa_x, pb_x, pb_y, pb_z, hs_10, \
                         hp_21, is_12, is_13, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_x[k] * hp_21[k];

        t_25[k] = pb_x[k] * is_12[k];

        t_26[k] = f_0 * hs_10[k]
                  + pb_y[k] * is_12[k];

        t_27[k] = pb_z[k] * is_12[k];

        t_28[k] = f_1 * hs_10[k]
                  + pb_z[k] * is_13[k];

        t_29[k] = pb_x[k] * is_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pa_y, pb_x, pb_z, hs_11, hs_12, \
                         hs_13, hp_21, is_14, is_15, is_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * hs_11[k]
                  + pb_z[k] * is_14[k];

        t_31[k] = pb_x[k] * is_15[k];

        t_32[k] = f_3 * hs_12[k]
                  + pb_z[k] * is_15[k];

        t_33[k] = pb_x[k] * is_16[k];

        t_34[k] = f_4 * hs_13[k]
                  + pb_z[k] * is_16[k];

        t_35[k] = pa_y[k] * hp_21[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_x, pb_y, pb_z, hs_15, \
                         is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_x[k] * is_17[k];

        t_37[k] = pb_y[k] * is_17[k];

        t_38[k] = f_0 * hs_15[k]
                  + pb_z[k] * is_17[k];
    }
}

auto
compute_prim_ip_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hs, const size_t hp,
                                     const size_t is, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_z, pb_x, pb_y, pb_z, hs_0, hs_3, hp_0, \
                         is_0, is_1, is_2, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = f_1 * hs_0[k]
                 + pb_y[k] * is_1[k];

        t_2[k] = pa_z[k] * hp_0[k];

        t_3[k] = f_1 * hs_0[k]
                 + pb_z[k] * is_2[k];

        t_4[k] = f_2 * hs_3[k]
                 + pb_x[k] * is_3[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, hs_1, hs_2, hs_3, hs_4, \
                         hs_5, is_3, is_4, is_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * hs_1[k]
                 + pb_y[k] * is_3[k];

        t_6[k] = f_2 * hs_4[k]
                 + pb_x[k] * is_4[k];

        t_7[k] = f_3 * hs_2[k]
                 + pb_z[k] * is_4[k];

        t_8[k] = f_4 * hs_5[k]
                 + pb_x[k] * is_5[k];

        t_9[k] = f_4 * hs_3[k]
                 + pb_y[k] * is_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_x, pb_y, pb_z, hs_4, hs_5, hs_6, \
                         hs_7, hs_9, is_6, is_7, is_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_4 * hs_6[k]
                  + pb_x[k] * is_6[k];

        t_11[k] = f_4 * hs_4[k]
                  + pb_z[k] * is_6[k];

        t_12[k] = f_3 * hs_7[k]
                  + pb_x[k] * is_7[k];

        t_13[k] = f_2 * hs_5[k]
                  + pb_y[k] * is_7[k];

        t_14[k] = f_3 * hs_9[k]
                  + pb_x[k] * is_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, hs_6, hs_10, hs_15, is_9, \
                         is_10, is_13, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_2 * hs_6[k]
                  + pb_z[k] * is_9[k];

        t_16[k] = f_1 * hs_10[k]
                  + pb_x[k] * is_10[k];

        t_17[k] = f_1 * hs_15[k]
                  + pb_x[k] * is_13[k];

        t_18[k] = f_0 * hs_10[k]
                  + pb_y[k] * is_14[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_z, pb_y, pb_z, hs_10, hs_11, hs_12, \
                         hs_13, hp_1, is_15, is_16, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_z[k] * hp_1[k];

        t_20[k] = f_1 * hs_10[k]
                  + pb_z[k] * is_15[k];

        t_21[k] = f_2 * hs_12[k]
                  + pb_y[k] * is_16[k];

        t_22[k] = f_3 * hs_11[k]
                  + pb_z[k] * is_16[k];

        t_23[k] = f_4 * hs_13[k]
                  + pb_y[k] * is_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pb_y, pb_z, hs_12, hs_13, hs_14, \
                         hs_15, hp_2, is_17, is_18, is_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_4 * hs_12[k]
                  + pb_z[k] * is_17[k];

        t_25[k] = f_3 * hs_14[k]
                  + pb_y[k] * is_18[k];

        t_26[k] = f_2 * hs_13[k]
                  + pb_z[k] * is_18[k];

        t_27[k] = f_1 * hs_15[k]
                  + pb_y[k] * is_19[k];

        t_28[k] = pa_y[k] * hp_2[k];
    }

#pragma omp simd aligned(t_29, pb_z, hs_15, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * hs_15[k]
                  + pb_z[k] * is_20[k];
    }
}

auto
compute_prim_ip_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hs, const size_t hp,
                                     const size_t is, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_16 = buffer.data(hs + 16);
    const auto *hs_17 = buffer.data(hs + 17);

    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_23 = buffer.data(hp + 23);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pb_x, pb_y, pb_z, hs_0, hs_1, \
                         is_0, is_1, is_2, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = f_1 * hs_0[k]
                 + pb_y[k] * is_1[k];

        t_4[k] = f_1 * hs_0[k]
                 + pb_z[k] * is_2[k];

        t_5[k] = f_2 * hs_1[k]
                 + pb_y[k] * is_3[k];

        t_6[k] = pb_z[k] * is_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pa_y, pb_y, pb_z, hs_2, hs_3, hp_4, \
                         is_4, is_5, is_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_y[k] * hp_4[k];

        t_8[k] = pb_y[k] * is_4[k];

        t_9[k] = f_2 * hs_2[k]
                 + pb_z[k] * is_4[k];

        t_10[k] = f_3 * hs_3[k]
                  + pb_y[k] * is_5[k];

        t_11[k] = pb_z[k] * is_5[k];

        t_12[k] = f_1 * hs_3[k]
                  + pb_z[k] * is_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pa_y, pb_y, pb_z, hs_4, hs_5, \
                         hp_6, is_7, is_8, is_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * hs_4[k]
                  + pb_y[k] * is_7[k];

        t_14[k] = pa_y[k] * hp_6[k];

        t_15[k] = pb_y[k] * is_8[k];

        t_16[k] = f_3 * hs_4[k]
                  + pb_z[k] * is_8[k];

        t_17[k] = f_4 * hs_5[k]
                  + pb_y[k] * is_9[k];

        t_18[k] = pb_z[k] * is_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pb_y, pb_z, hs_5, hs_6, hs_7, \
                         hs_8, hp_8, is_10, is_11, is_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * hs_5[k]
                  + pb_z[k] * is_10[k];

        t_20[k] = f_2 * hs_7[k]
                  + pb_y[k] * is_11[k];

        t_21[k] = f_2 * hs_6[k]
                  + pb_z[k] * is_11[k];

        t_22[k] = f_1 * hs_8[k]
                  + pb_y[k] * is_12[k];

        t_23[k] = pa_y[k] * hp_8[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, pa_x, pb_y, pb_z, hs_8, hp_12, \
                         hp_14, hp_15, hp_16, is_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * is_13[k];

        t_25[k] = f_4 * hs_8[k]
                  + pb_z[k] * is_13[k];

        t_26[k] = pa_x[k] * hp_12[k];

        t_27[k] = pa_x[k] * hp_14[k];

        t_28[k] = pa_x[k] * hp_15[k];

        t_29[k] = pa_x[k] * hp_16[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pa_x, pb_x, pb_y, hs_12, hp_17, \
                         hp_18, hp_19, hp_23, is_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_x[k] * hp_17[k];

        t_31[k] = pa_x[k] * hp_18[k];

        t_32[k] = pa_x[k] * hp_19[k];

        t_33[k] = pa_x[k] * hp_23[k];

        t_34[k] = pb_x[k] * is_16[k];

        t_35[k] = f_0 * hs_12[k]
                  + pb_y[k] * is_16[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, pb_x, pb_y, pb_z, hs_12, hs_13, \
                         hs_14, is_16, is_17, is_18, is_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_z[k] * is_16[k];

        t_37[k] = f_1 * hs_12[k]
                  + pb_z[k] * is_17[k];

        t_38[k] = pb_x[k] * is_18[k];

        t_39[k] = f_4 * hs_14[k]
                  + pb_y[k] * is_18[k];

        t_40[k] = f_2 * hs_13[k]
                  + pb_z[k] * is_18[k];

        t_41[k] = pb_x[k] * is_19[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, hs_14, hs_15, \
                         hs_16, hs_17, is_19, is_20, is_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * hs_15[k]
                  + pb_y[k] * is_19[k];

        t_43[k] = f_3 * hs_14[k]
                  + pb_z[k] * is_19[k];

        t_44[k] = pb_x[k] * is_20[k];

        t_45[k] = f_2 * hs_16[k]
                  + pb_y[k] * is_20[k];

        t_46[k] = f_4 * hs_15[k]
                  + pb_z[k] * is_20[k];

        t_47[k] = f_1 * hs_17[k]
                  + pb_y[k] * is_21[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pb_x, pb_y, pb_z, hs_17, hp_23, \
                         is_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_y[k] * hp_23[k];

        t_49[k] = pb_x[k] * is_22[k];

        t_50[k] = pb_y[k] * is_22[k];

        t_51[k] = f_0 * hs_17[k]
                  + pb_z[k] * is_22[k];
    }
}

auto
compute_prim_ip_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hs, const size_t hp,
                                     const size_t is, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);

    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_8 = buffer.data(hp + 8);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pb_x, pb_y, pb_z, hs_0, hs_1, \
                         is_0, is_1, is_2, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = f_1 * hs_0[k]
                 + pb_y[k] * is_1[k];

        t_4[k] = f_1 * hs_0[k]
                 + pb_z[k] * is_2[k];

        t_5[k] = f_2 * hs_1[k]
                 + pb_y[k] * is_3[k];

        t_6[k] = pb_z[k] * is_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pb_y, pb_z, hs_2, hs_3, hs_4, is_4, \
                         is_5, is_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_y[k] * is_4[k];

        t_8[k] = f_2 * hs_2[k]
                 + pb_z[k] * is_4[k];

        t_9[k] = f_3 * hs_3[k]
                 + pb_y[k] * is_5[k];

        t_10[k] = pb_z[k] * is_5[k];

        t_11[k] = pb_y[k] * is_6[k];

        t_12[k] = f_3 * hs_4[k]
                  + pb_z[k] * is_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pa_x, pb_y, pb_z, hs_5, hs_6, \
                         hp_4, hp_8, is_7, is_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_4 * hs_5[k]
                  + pb_y[k] * is_7[k];

        t_14[k] = pb_z[k] * is_7[k];

        t_15[k] = pb_y[k] * is_9[k];

        t_16[k] = f_4 * hs_6[k]
                  + pb_z[k] * is_9[k];

        t_17[k] = pa_x[k] * hp_4[k];

        t_18[k] = pa_x[k] * hp_8[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, t_25, pb_x, pb_y, pb_z, hs_10, \
                         hs_11, hs_12, is_14, is_15, is_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_x[k] * is_14[k];

        t_20[k] = f_0 * hs_10[k]
                  + pb_y[k] * is_14[k];

        t_21[k] = pb_z[k] * is_14[k];

        t_22[k] = f_1 * hs_10[k]
                  + pb_z[k] * is_15[k];

        t_23[k] = pb_x[k] * is_16[k];

        t_24[k] = f_4 * hs_12[k]
                  + pb_y[k] * is_16[k];

        t_25[k] = f_2 * hs_11[k]
                  + pb_z[k] * is_16[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, pb_x, pb_y, pb_z, hs_12, hs_13, \
                         hs_14, is_17, is_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_x[k] * is_17[k];

        t_27[k] = f_3 * hs_13[k]
                  + pb_y[k] * is_17[k];

        t_28[k] = f_3 * hs_12[k]
                  + pb_z[k] * is_17[k];

        t_29[k] = pb_x[k] * is_18[k];

        t_30[k] = f_2 * hs_14[k]
                  + pb_y[k] * is_18[k];

        t_31[k] = f_4 * hs_13[k]
                  + pb_z[k] * is_18[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pa_y, pb_x, pb_y, pb_z, hs_15, hp_8, \
                         is_19, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * hs_15[k]
                  + pb_y[k] * is_19[k];

        t_33[k] = pa_y[k] * hp_8[k];

        t_34[k] = pb_x[k] * is_20[k];

        t_35[k] = pb_y[k] * is_20[k];

        t_36[k] = f_0 * hs_15[k]
                  + pb_z[k] * is_20[k];
    }
}

auto
compute_prim_ip_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t hs, const size_t is, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_11 = buffer.data(hs + 11);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_14 = buffer.data(is + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pb_x, pb_y, pb_z, hs_0, is_0, \
                         is_1, is_2, is_3, is_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = pb_z[k] * is_1[k];

        t_4[k] = pb_y[k] * is_2[k];

        t_5[k] = pb_z[k] * is_3[k];

        t_6[k] = pb_y[k] * is_4[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, t_13, pb_x, pb_y, pb_z, hs_7, is_5, \
                         is_6, is_9, is_10, is_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_z[k] * is_5[k];

        t_8[k] = pb_y[k] * is_6[k];

        t_9[k] = pb_x[k] * is_9[k];

        t_10[k] = f_0 * hs_7[k]
                  + pb_y[k] * is_9[k];

        t_11[k] = pb_z[k] * is_9[k];

        t_12[k] = pb_x[k] * is_10[k];

        t_13[k] = pb_x[k] * is_11[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_x, pb_y, pb_z, hs_11, is_12, \
                         is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pb_x[k] * is_12[k];

        t_15[k] = pb_x[k] * is_14[k];

        t_16[k] = pb_y[k] * is_14[k];

        t_17[k] = f_0 * hs_11[k]
                  + pb_z[k] * is_14[k];
    }
}

auto
compute_prim_ip_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hs, const size_t hp,
                                     const size_t is, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_16 = buffer.data(hs + 16);

    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_21 = buffer.data(hp + 21);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_20 = buffer.data(is + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pb_x, pb_y, pb_z, hs_0, hs_1, \
                         is_0, is_1, is_2, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = f_1 * hs_0[k]
                 + pb_y[k] * is_1[k];

        t_4[k] = f_1 * hs_0[k]
                 + pb_z[k] * is_2[k];

        t_5[k] = f_2 * hs_1[k]
                 + pb_y[k] * is_3[k];

        t_6[k] = pb_z[k] * is_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pa_y, pb_y, pb_z, hs_2, hs_3, hp_3, \
                         is_4, is_5, is_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_y[k] * hp_3[k];

        t_8[k] = pb_y[k] * is_4[k];

        t_9[k] = f_2 * hs_2[k]
                 + pb_z[k] * is_4[k];

        t_10[k] = f_3 * hs_3[k]
                  + pb_y[k] * is_5[k];

        t_11[k] = pb_z[k] * is_5[k];

        t_12[k] = f_1 * hs_3[k]
                  + pb_z[k] * is_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pa_y, pb_y, pb_z, hs_4, hs_5, \
                         hp_6, is_7, is_8, is_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * hp_6[k];

        t_14[k] = pb_y[k] * is_7[k];

        t_15[k] = f_3 * hs_4[k]
                  + pb_z[k] * is_7[k];

        t_16[k] = f_4 * hs_5[k]
                  + pb_y[k] * is_8[k];

        t_17[k] = pb_z[k] * is_8[k];

        t_18[k] = f_1 * hs_5[k]
                  + pb_z[k] * is_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_x, pa_y, pb_y, pb_z, hs_6, hs_7, \
                         hp_9, hp_11, is_10, is_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * hs_6[k]
                  + pb_z[k] * is_10[k];

        t_20[k] = pa_y[k] * hp_9[k];

        t_21[k] = pb_y[k] * is_11[k];

        t_22[k] = f_4 * hs_7[k]
                  + pb_z[k] * is_11[k];

        t_23[k] = pa_x[k] * hp_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, pa_x, pb_x, pb_y, pb_z, hs_11, \
                         hp_21, is_14, is_15, is_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_x[k] * hp_21[k];

        t_25[k] = pb_x[k] * is_14[k];

        t_26[k] = f_0 * hs_11[k]
                  + pb_y[k] * is_14[k];

        t_27[k] = pb_z[k] * is_14[k];

        t_28[k] = f_1 * hs_11[k]
                  + pb_z[k] * is_15[k];

        t_29[k] = pb_x[k] * is_16[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pa_y, pb_x, pb_z, hs_12, hs_13, \
                         hs_14, hp_21, is_16, is_17, is_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * hs_12[k]
                  + pb_z[k] * is_16[k];

        t_31[k] = pb_x[k] * is_17[k];

        t_32[k] = f_3 * hs_13[k]
                  + pb_z[k] * is_17[k];

        t_33[k] = pb_x[k] * is_18[k];

        t_34[k] = f_4 * hs_14[k]
                  + pb_z[k] * is_18[k];

        t_35[k] = pa_y[k] * hp_21[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_x, pb_y, pb_z, hs_16, \
                         is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_x[k] * is_20[k];

        t_37[k] = pb_y[k] * is_20[k];

        t_38[k] = f_0 * hs_16[k]
                  + pb_z[k] * is_20[k];
    }
}

auto
compute_prim_ip_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hs, const size_t hp,
                                      const size_t is, const size_t ncols,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_14 = buffer.data(hs + 14);

    const auto *hp_8 = buffer.data(hp + 8);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_19 = buffer.data(is + 19);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pb_x, pb_y, pb_z, hs_0, hs_1, \
                         is_0, is_1, is_2, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = f_1 * hs_0[k]
                 + pb_z[k] * is_1[k];

        t_4[k] = pb_z[k] * is_2[k];

        t_5[k] = pb_y[k] * is_3[k];

        t_6[k] = f_2 * hs_1[k]
                 + pb_z[k] * is_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pb_y, pb_z, hs_3, hs_5, is_4, is_5, \
                         is_6, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_z[k] * is_4[k];

        t_8[k] = pb_y[k] * is_5[k];

        t_9[k] = f_3 * hs_3[k]
                 + pb_z[k] * is_5[k];

        t_10[k] = pb_z[k] * is_6[k];

        t_11[k] = pb_y[k] * is_8[k];

        t_12[k] = f_4 * hs_5[k]
                  + pb_z[k] * is_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, t_19, pb_x, pb_y, pb_z, hs_9, \
                         hs_10, is_13, is_14, is_15, is_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_x[k] * is_13[k];

        t_14[k] = f_0 * hs_9[k]
                  + pb_y[k] * is_13[k];

        t_15[k] = pb_z[k] * is_13[k];

        t_16[k] = f_1 * hs_9[k]
                  + pb_z[k] * is_14[k];

        t_17[k] = pb_x[k] * is_15[k];

        t_18[k] = f_2 * hs_10[k]
                  + pb_z[k] * is_15[k];

        t_19[k] = pb_x[k] * is_16[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_y, pb_z, hs_11, \
                         hs_12, hp_8, is_16, is_17, is_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_3 * hs_11[k]
                  + pb_z[k] * is_16[k];

        t_21[k] = pb_x[k] * is_17[k];

        t_22[k] = f_4 * hs_12[k]
                  + pb_z[k] * is_17[k];

        t_23[k] = pa_y[k] * hp_8[k];

        t_24[k] = pb_x[k] * is_19[k];

        t_25[k] = pb_y[k] * is_19[k];
    }

#pragma omp simd aligned(t_26, pb_z, hs_14, is_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hs_14[k]
                  + pb_z[k] * is_19[k];
    }
}

auto
compute_prim_ip_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t hs, const size_t is, const size_t ncols,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_12 = buffer.data(hs + 12);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_17 = buffer.data(is + 17);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hs_0, hs_8, hs_12, is_0, is_12, \
                         is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = f_0 * hs_8[k]
                 + pb_y[k] * is_12[k];

        t_2[k] = f_0 * hs_12[k]
                 + pb_z[k] * is_17[k];
    }
}

auto
compute_prim_ip_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t hs, const size_t is, const size_t ncols,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_12 = buffer.data(hs + 12);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_24 = buffer.data(is + 24);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hs_0, hs_8, hs_12, is_0, is_18, \
                         is_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = f_0 * hs_8[k]
                 + pb_y[k] * is_18[k];

        t_2[k] = f_0 * hs_12[k]
                 + pb_z[k] * is_24[k];
    }
}

auto
compute_prim_ip_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hs, const size_t hp,
                                      const size_t is, const size_t ncols,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);

    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, hs_0, hs_1, is_0, \
                         is_1, is_2, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = f_1 * hs_0[k]
                 + pb_y[k] * is_1[k];

        t_4[k] = f_1 * hs_0[k]
                 + pb_z[k] * is_2[k];

        t_5[k] = f_2 * hs_1[k]
                 + pb_y[k] * is_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_y, pb_z, hs_2, hs_3, hs_4, hs_5, is_4, is_5, \
                         is_6, is_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * hs_2[k]
                 + pb_z[k] * is_4[k];

        t_7[k] = f_3 * hs_3[k]
                 + pb_y[k] * is_5[k];

        t_8[k] = f_3 * hs_4[k]
                 + pb_z[k] * is_6[k];

        t_9[k] = f_4 * hs_5[k]
                 + pb_y[k] * is_7[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_x, pb_x, pb_y, pb_z, hs_6, \
                         hs_10, hp_1, hp_2, is_9, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_4 * hs_6[k]
                  + pb_z[k] * is_9[k];

        t_11[k] = pa_x[k] * hp_1[k];

        t_12[k] = pa_x[k] * hp_2[k];

        t_13[k] = pb_x[k] * is_14[k];

        t_14[k] = f_0 * hs_10[k]
                  + pb_y[k] * is_14[k];

        t_15[k] = pb_z[k] * is_14[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pb_y, pb_z, hs_10, hs_11, hs_12, hs_13, \
                         is_15, is_16, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * hs_10[k]
                  + pb_z[k] * is_15[k];

        t_17[k] = f_4 * hs_12[k]
                  + pb_y[k] * is_16[k];

        t_18[k] = f_2 * hs_11[k]
                  + pb_z[k] * is_16[k];

        t_19[k] = f_3 * hs_13[k]
                  + pb_y[k] * is_17[k];

        t_20[k] = f_3 * hs_12[k]
                  + pb_z[k] * is_17[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_y, pb_z, hs_13, hs_14, \
                         hs_15, hp_2, is_18, is_19, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_2 * hs_14[k]
                  + pb_y[k] * is_18[k];

        t_22[k] = f_4 * hs_13[k]
                  + pb_z[k] * is_18[k];

        t_23[k] = f_1 * hs_15[k]
                  + pb_y[k] * is_19[k];

        t_24[k] = pa_y[k] * hp_2[k];

        t_25[k] = pb_x[k] * is_20[k];
    }

#pragma omp simd aligned(t_26, t_27, pb_y, pb_z, hs_15, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * is_20[k];

        t_27[k] = f_0 * hs_15[k]
                  + pb_z[k] * is_20[k];
    }
}

auto
compute_prim_ip_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t hs, const size_t is, const size_t ncols,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_12 = buffer.data(hs + 12);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_17 = buffer.data(is + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, hs_0, hs_8, \
                         is_0, is_12, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = pb_x[k] * is_12[k];

        t_4[k] = f_0 * hs_8[k]
                 + pb_y[k] * is_12[k];

        t_5[k] = pb_z[k] * is_12[k];

        t_6[k] = pb_x[k] * is_17[k];

        t_7[k] = pb_y[k] * is_17[k];
    }

#pragma omp simd aligned(t_8, pb_z, hs_12, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hs_12[k]
                 + pb_z[k] * is_17[k];
    }
}

auto
compute_prim_ip_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t hs, const size_t is, const size_t ncols,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_11 = buffer.data(hs + 11);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_14 = buffer.data(is + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hs_0, hs_7, hs_11, is_0, is_9, \
                         is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = f_0 * hs_7[k]
                 + pb_y[k] * is_9[k];

        t_2[k] = f_0 * hs_11[k]
                 + pb_z[k] * is_14[k];
    }
}

auto
compute_prim_ip_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hs, const size_t hp,
                                      const size_t is, const size_t ncols,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_15 = buffer.data(hs + 15);

    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_15 = buffer.data(hp + 15);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_22 = buffer.data(is + 22);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_y, pb_x, pb_y, pb_z, hs_0, hp_3, \
                         is_0, is_2, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = f_1 * hs_0[k]
                 + pb_z[k] * is_2[k];

        t_4[k] = pb_z[k] * is_3[k];

        t_5[k] = pa_y[k] * hp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pa_y, pb_y, pb_z, hs_1, hs_2, hp_4, \
                         is_4, is_5, is_6, is_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_y[k] * is_4[k];

        t_7[k] = f_2 * hs_1[k]
                 + pb_z[k] * is_4[k];

        t_8[k] = pb_z[k] * is_5[k];

        t_9[k] = f_1 * hs_2[k]
                 + pb_z[k] * is_6[k];

        t_10[k] = pa_y[k] * hp_4[k];

        t_11[k] = pb_y[k] * is_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pa_y, pb_z, hs_3, hs_4, hs_5, hp_5, \
                         is_7, is_8, is_9, is_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * hs_3[k]
                  + pb_z[k] * is_7[k];

        t_13[k] = pb_z[k] * is_8[k];

        t_14[k] = f_1 * hs_4[k]
                  + pb_z[k] * is_9[k];

        t_15[k] = f_2 * hs_5[k]
                  + pb_z[k] * is_10[k];

        t_16[k] = pa_y[k] * hp_5[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, t_22, t_23, pb_x, pb_y, pb_z, hs_6, \
                         hs_10, is_11, is_16, is_17, is_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = pb_y[k] * is_11[k];

        t_18[k] = f_4 * hs_6[k]
                  + pb_z[k] * is_11[k];

        t_19[k] = pb_x[k] * is_16[k];

        t_20[k] = f_0 * hs_10[k]
                  + pb_y[k] * is_16[k];

        t_21[k] = pb_z[k] * is_16[k];

        t_22[k] = f_1 * hs_10[k]
                  + pb_z[k] * is_17[k];

        t_23[k] = pb_x[k] * is_18[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, pa_y, pb_x, pb_z, hs_11, hs_12, \
                         hs_13, hp_15, is_18, is_19, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_2 * hs_11[k]
                  + pb_z[k] * is_18[k];

        t_25[k] = pb_x[k] * is_19[k];

        t_26[k] = f_3 * hs_12[k]
                  + pb_z[k] * is_19[k];

        t_27[k] = pb_x[k] * is_20[k];

        t_28[k] = f_4 * hs_13[k]
                  + pb_z[k] * is_20[k];

        t_29[k] = pa_y[k] * hp_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_x, pb_y, pb_z, hs_15, \
                         is_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_x[k] * is_22[k];

        t_31[k] = pb_y[k] * is_22[k];

        t_32[k] = f_0 * hs_15[k]
                  + pb_z[k] * is_22[k];
    }
}

auto
compute_prim_ip_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t hs, const size_t is, const size_t ncols,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_11 = buffer.data(hs + 11);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_14 = buffer.data(is + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, hs_0, hs_7, \
                         is_0, is_9, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pb_x[k] * is_0[k];

        t_1[k] = pb_y[k] * is_0[k];

        t_2[k] = pb_z[k] * is_0[k];

        t_3[k] = pb_x[k] * is_9[k];

        t_4[k] = f_0 * hs_7[k]
                 + pb_y[k] * is_9[k];

        t_5[k] = pb_z[k] * is_9[k];

        t_6[k] = pb_x[k] * is_14[k];

        t_7[k] = pb_y[k] * is_14[k];
    }

#pragma omp simd aligned(t_8, pb_z, hs_11, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hs_11[k]
                 + pb_z[k] * is_14[k];
    }
}

}  // namespace simdt2ceri
