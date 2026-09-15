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


#include "SimdElectronRepulsionGeom10VrrRecIP.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_ip_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t hp, const size_t kp,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);
    const auto *hp_30 = buffer.data(hp + 30);
    const auto *hp_31 = buffer.data(hp + 31);
    const auto *hp_32 = buffer.data(hp + 32);
    const auto *hp_33 = buffer.data(hp + 33);
    const auto *hp_34 = buffer.data(hp + 34);
    const auto *hp_35 = buffer.data(hp + 35);
    const auto *hp_36 = buffer.data(hp + 36);
    const auto *hp_37 = buffer.data(hp + 37);
    const auto *hp_38 = buffer.data(hp + 38);
    const auto *hp_39 = buffer.data(hp + 39);
    const auto *hp_40 = buffer.data(hp + 40);
    const auto *hp_41 = buffer.data(hp + 41);
    const auto *hp_42 = buffer.data(hp + 42);
    const auto *hp_43 = buffer.data(hp + 43);
    const auto *hp_44 = buffer.data(hp + 44);
    const auto *hp_45 = buffer.data(hp + 45);
    const auto *hp_46 = buffer.data(hp + 46);
    const auto *hp_47 = buffer.data(hp + 47);
    const auto *hp_48 = buffer.data(hp + 48);
    const auto *hp_49 = buffer.data(hp + 49);
    const auto *hp_50 = buffer.data(hp + 50);
    const auto *hp_51 = buffer.data(hp + 51);
    const auto *hp_52 = buffer.data(hp + 52);
    const auto *hp_53 = buffer.data(hp + 53);
    const auto *hp_54 = buffer.data(hp + 54);
    const auto *hp_55 = buffer.data(hp + 55);
    const auto *hp_56 = buffer.data(hp + 56);
    const auto *hp_57 = buffer.data(hp + 57);
    const auto *hp_58 = buffer.data(hp + 58);
    const auto *hp_59 = buffer.data(hp + 59);
    const auto *hp_60 = buffer.data(hp + 60);
    const auto *hp_61 = buffer.data(hp + 61);
    const auto *hp_62 = buffer.data(hp + 62);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_34 = buffer.data(kp + 34);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_39 = buffer.data(kp + 39);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_48 = buffer.data(kp + 48);
    const auto *kp_49 = buffer.data(kp + 49);
    const auto *kp_50 = buffer.data(kp + 50);
    const auto *kp_51 = buffer.data(kp + 51);
    const auto *kp_52 = buffer.data(kp + 52);
    const auto *kp_53 = buffer.data(kp + 53);
    const auto *kp_54 = buffer.data(kp + 54);
    const auto *kp_55 = buffer.data(kp + 55);
    const auto *kp_56 = buffer.data(kp + 56);
    const auto *kp_57 = buffer.data(kp + 57);
    const auto *kp_58 = buffer.data(kp + 58);
    const auto *kp_59 = buffer.data(kp + 59);
    const auto *kp_60 = buffer.data(kp + 60);
    const auto *kp_61 = buffer.data(kp + 61);
    const auto *kp_62 = buffer.data(kp + 62);
    const auto *kp_63 = buffer.data(kp + 63);
    const auto *kp_64 = buffer.data(kp + 64);
    const auto *kp_65 = buffer.data(kp + 65);
    const auto *kp_66 = buffer.data(kp + 66);
    const auto *kp_67 = buffer.data(kp + 67);
    const auto *kp_68 = buffer.data(kp + 68);
    const auto *kp_69 = buffer.data(kp + 69);
    const auto *kp_70 = buffer.data(kp + 70);
    const auto *kp_71 = buffer.data(kp + 71);
    const auto *kp_72 = buffer.data(kp + 72);
    const auto *kp_73 = buffer.data(kp + 73);
    const auto *kp_74 = buffer.data(kp + 74);
    const auto *kp_75 = buffer.data(kp + 75);
    const auto *kp_76 = buffer.data(kp + 76);
    const auto *kp_77 = buffer.data(kp + 77);
    const auto *kp_78 = buffer.data(kp + 78);
    const auto *kp_79 = buffer.data(kp + 79);
    const auto *kp_80 = buffer.data(kp + 80);
    const auto *kp_81 = buffer.data(kp + 81);
    const auto *kp_82 = buffer.data(kp + 82);
    const auto *kp_83 = buffer.data(kp + 83);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, hp_0, hp_1, hp_2, hp_3, hp_4, kp_0, kp_1, \
                         kp_2, kp_3, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -6.0 * hp_0[k]
                 + f_0 * kp_0[k];

        t_1[k] = -6.0 * hp_1[k]
                 + f_0 * kp_1[k];

        t_2[k] = -6.0 * hp_2[k]
                 + f_0 * kp_2[k];

        t_3[k] = -5.0 * hp_3[k]
                 + f_0 * kp_3[k];

        t_4[k] = -5.0 * hp_4[k]
                 + f_0 * kp_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, hp_5, hp_6, hp_7, hp_8, hp_9, kp_5, kp_6, \
                         kp_7, kp_8, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -5.0 * hp_5[k]
                 + f_0 * kp_5[k];

        t_6[k] = -5.0 * hp_6[k]
                 + f_0 * kp_6[k];

        t_7[k] = -5.0 * hp_7[k]
                 + f_0 * kp_7[k];

        t_8[k] = -5.0 * hp_8[k]
                 + f_0 * kp_8[k];

        t_9[k] = -4.0 * hp_9[k]
                 + f_0 * kp_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, hp_10, hp_11, hp_12, hp_13, hp_14, \
                         kp_10, kp_11, kp_12, kp_13, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -4.0 * hp_10[k]
                  + f_0 * kp_10[k];

        t_11[k] = -4.0 * hp_11[k]
                  + f_0 * kp_11[k];

        t_12[k] = -4.0 * hp_12[k]
                  + f_0 * kp_12[k];

        t_13[k] = -4.0 * hp_13[k]
                  + f_0 * kp_13[k];

        t_14[k] = -4.0 * hp_14[k]
                  + f_0 * kp_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, hp_15, hp_16, hp_17, hp_18, hp_19, \
                         kp_15, kp_16, kp_17, kp_18, kp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -4.0 * hp_15[k]
                  + f_0 * kp_15[k];

        t_16[k] = -4.0 * hp_16[k]
                  + f_0 * kp_16[k];

        t_17[k] = -4.0 * hp_17[k]
                  + f_0 * kp_17[k];

        t_18[k] = -3.0 * hp_18[k]
                  + f_0 * kp_18[k];

        t_19[k] = -3.0 * hp_19[k]
                  + f_0 * kp_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, hp_20, hp_21, hp_22, hp_23, hp_24, \
                         kp_20, kp_21, kp_22, kp_23, kp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -3.0 * hp_20[k]
                  + f_0 * kp_20[k];

        t_21[k] = -3.0 * hp_21[k]
                  + f_0 * kp_21[k];

        t_22[k] = -3.0 * hp_22[k]
                  + f_0 * kp_22[k];

        t_23[k] = -3.0 * hp_23[k]
                  + f_0 * kp_23[k];

        t_24[k] = -3.0 * hp_24[k]
                  + f_0 * kp_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, hp_25, hp_26, hp_27, hp_28, hp_29, \
                         kp_25, kp_26, kp_27, kp_28, kp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -3.0 * hp_25[k]
                  + f_0 * kp_25[k];

        t_26[k] = -3.0 * hp_26[k]
                  + f_0 * kp_26[k];

        t_27[k] = -3.0 * hp_27[k]
                  + f_0 * kp_27[k];

        t_28[k] = -3.0 * hp_28[k]
                  + f_0 * kp_28[k];

        t_29[k] = -3.0 * hp_29[k]
                  + f_0 * kp_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, hp_30, hp_31, hp_32, hp_33, hp_34, \
                         kp_30, kp_31, kp_32, kp_33, kp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -2.0 * hp_30[k]
                  + f_0 * kp_30[k];

        t_31[k] = -2.0 * hp_31[k]
                  + f_0 * kp_31[k];

        t_32[k] = -2.0 * hp_32[k]
                  + f_0 * kp_32[k];

        t_33[k] = -2.0 * hp_33[k]
                  + f_0 * kp_33[k];

        t_34[k] = -2.0 * hp_34[k]
                  + f_0 * kp_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, hp_35, hp_36, hp_37, hp_38, hp_39, \
                         kp_35, kp_36, kp_37, kp_38, kp_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -2.0 * hp_35[k]
                  + f_0 * kp_35[k];

        t_36[k] = -2.0 * hp_36[k]
                  + f_0 * kp_36[k];

        t_37[k] = -2.0 * hp_37[k]
                  + f_0 * kp_37[k];

        t_38[k] = -2.0 * hp_38[k]
                  + f_0 * kp_38[k];

        t_39[k] = -2.0 * hp_39[k]
                  + f_0 * kp_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, hp_40, hp_41, hp_42, hp_43, hp_44, \
                         kp_40, kp_41, kp_42, kp_43, kp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -2.0 * hp_40[k]
                  + f_0 * kp_40[k];

        t_41[k] = -2.0 * hp_41[k]
                  + f_0 * kp_41[k];

        t_42[k] = -2.0 * hp_42[k]
                  + f_0 * kp_42[k];

        t_43[k] = -2.0 * hp_43[k]
                  + f_0 * kp_43[k];

        t_44[k] = -2.0 * hp_44[k]
                  + f_0 * kp_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, hp_45, hp_46, hp_47, hp_48, hp_49, \
                         kp_45, kp_46, kp_47, kp_48, kp_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -hp_45[k]
                  + f_0 * kp_45[k];

        t_46[k] = -hp_46[k]
                  + f_0 * kp_46[k];

        t_47[k] = -hp_47[k]
                  + f_0 * kp_47[k];

        t_48[k] = -hp_48[k]
                  + f_0 * kp_48[k];

        t_49[k] = -hp_49[k]
                  + f_0 * kp_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, hp_50, hp_51, hp_52, hp_53, hp_54, \
                         kp_50, kp_51, kp_52, kp_53, kp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -hp_50[k]
                  + f_0 * kp_50[k];

        t_51[k] = -hp_51[k]
                  + f_0 * kp_51[k];

        t_52[k] = -hp_52[k]
                  + f_0 * kp_52[k];

        t_53[k] = -hp_53[k]
                  + f_0 * kp_53[k];

        t_54[k] = -hp_54[k]
                  + f_0 * kp_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, hp_55, hp_56, hp_57, hp_58, hp_59, \
                         kp_55, kp_56, kp_57, kp_58, kp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -hp_55[k]
                  + f_0 * kp_55[k];

        t_56[k] = -hp_56[k]
                  + f_0 * kp_56[k];

        t_57[k] = -hp_57[k]
                  + f_0 * kp_57[k];

        t_58[k] = -hp_58[k]
                  + f_0 * kp_58[k];

        t_59[k] = -hp_59[k]
                  + f_0 * kp_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, hp_60, hp_61, hp_62, kp_60, \
                         kp_61, kp_62, kp_63, kp_64, kp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -hp_60[k]
                  + f_0 * kp_60[k];

        t_61[k] = -hp_61[k]
                  + f_0 * kp_61[k];

        t_62[k] = -hp_62[k]
                  + f_0 * kp_62[k];

        t_63[k] = f_0 * kp_63[k];

        t_64[k] = f_0 * kp_64[k];

        t_65[k] = f_0 * kp_65[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, t_72, t_73, kp_66, kp_67, kp_68, \
                         kp_69, kp_70, kp_71, kp_72, kp_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_0 * kp_66[k];

        t_67[k] = f_0 * kp_67[k];

        t_68[k] = f_0 * kp_68[k];

        t_69[k] = f_0 * kp_69[k];

        t_70[k] = f_0 * kp_70[k];

        t_71[k] = f_0 * kp_71[k];

        t_72[k] = f_0 * kp_72[k];

        t_73[k] = f_0 * kp_73[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, t_80, t_81, kp_74, kp_75, kp_76, \
                         kp_77, kp_78, kp_79, kp_80, kp_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * kp_74[k];

        t_75[k] = f_0 * kp_75[k];

        t_76[k] = f_0 * kp_76[k];

        t_77[k] = f_0 * kp_77[k];

        t_78[k] = f_0 * kp_78[k];

        t_79[k] = f_0 * kp_79[k];

        t_80[k] = f_0 * kp_80[k];

        t_81[k] = f_0 * kp_81[k];
    }

#pragma omp simd aligned(t_82, t_83, kp_82, kp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * kp_82[k];

        t_83[k] = f_0 * kp_83[k];
    }
}

auto
compute_prim_geom_10_ip_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t hp, const size_t kp,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);
    const auto *hp_30 = buffer.data(hp + 30);
    const auto *hp_31 = buffer.data(hp + 31);
    const auto *hp_32 = buffer.data(hp + 32);
    const auto *hp_33 = buffer.data(hp + 33);
    const auto *hp_34 = buffer.data(hp + 34);
    const auto *hp_35 = buffer.data(hp + 35);
    const auto *hp_36 = buffer.data(hp + 36);
    const auto *hp_37 = buffer.data(hp + 37);
    const auto *hp_38 = buffer.data(hp + 38);
    const auto *hp_39 = buffer.data(hp + 39);
    const auto *hp_40 = buffer.data(hp + 40);
    const auto *hp_41 = buffer.data(hp + 41);
    const auto *hp_42 = buffer.data(hp + 42);
    const auto *hp_43 = buffer.data(hp + 43);
    const auto *hp_44 = buffer.data(hp + 44);
    const auto *hp_45 = buffer.data(hp + 45);
    const auto *hp_46 = buffer.data(hp + 46);
    const auto *hp_47 = buffer.data(hp + 47);
    const auto *hp_48 = buffer.data(hp + 48);
    const auto *hp_49 = buffer.data(hp + 49);
    const auto *hp_50 = buffer.data(hp + 50);
    const auto *hp_51 = buffer.data(hp + 51);
    const auto *hp_52 = buffer.data(hp + 52);
    const auto *hp_53 = buffer.data(hp + 53);
    const auto *hp_54 = buffer.data(hp + 54);
    const auto *hp_55 = buffer.data(hp + 55);
    const auto *hp_56 = buffer.data(hp + 56);
    const auto *hp_57 = buffer.data(hp + 57);
    const auto *hp_58 = buffer.data(hp + 58);
    const auto *hp_59 = buffer.data(hp + 59);
    const auto *hp_60 = buffer.data(hp + 60);
    const auto *hp_61 = buffer.data(hp + 61);
    const auto *hp_62 = buffer.data(hp + 62);

    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_34 = buffer.data(kp + 34);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_39 = buffer.data(kp + 39);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_48 = buffer.data(kp + 48);
    const auto *kp_49 = buffer.data(kp + 49);
    const auto *kp_50 = buffer.data(kp + 50);
    const auto *kp_51 = buffer.data(kp + 51);
    const auto *kp_52 = buffer.data(kp + 52);
    const auto *kp_53 = buffer.data(kp + 53);
    const auto *kp_54 = buffer.data(kp + 54);
    const auto *kp_55 = buffer.data(kp + 55);
    const auto *kp_56 = buffer.data(kp + 56);
    const auto *kp_57 = buffer.data(kp + 57);
    const auto *kp_58 = buffer.data(kp + 58);
    const auto *kp_59 = buffer.data(kp + 59);
    const auto *kp_63 = buffer.data(kp + 63);
    const auto *kp_64 = buffer.data(kp + 64);
    const auto *kp_65 = buffer.data(kp + 65);
    const auto *kp_66 = buffer.data(kp + 66);
    const auto *kp_67 = buffer.data(kp + 67);
    const auto *kp_68 = buffer.data(kp + 68);
    const auto *kp_69 = buffer.data(kp + 69);
    const auto *kp_70 = buffer.data(kp + 70);
    const auto *kp_71 = buffer.data(kp + 71);
    const auto *kp_72 = buffer.data(kp + 72);
    const auto *kp_73 = buffer.data(kp + 73);
    const auto *kp_74 = buffer.data(kp + 74);
    const auto *kp_75 = buffer.data(kp + 75);
    const auto *kp_76 = buffer.data(kp + 76);
    const auto *kp_77 = buffer.data(kp + 77);
    const auto *kp_78 = buffer.data(kp + 78);
    const auto *kp_79 = buffer.data(kp + 79);
    const auto *kp_80 = buffer.data(kp + 80);
    const auto *kp_84 = buffer.data(kp + 84);
    const auto *kp_85 = buffer.data(kp + 85);
    const auto *kp_86 = buffer.data(kp + 86);
    const auto *kp_87 = buffer.data(kp + 87);
    const auto *kp_88 = buffer.data(kp + 88);
    const auto *kp_89 = buffer.data(kp + 89);
    const auto *kp_90 = buffer.data(kp + 90);
    const auto *kp_91 = buffer.data(kp + 91);
    const auto *kp_92 = buffer.data(kp + 92);
    const auto *kp_93 = buffer.data(kp + 93);
    const auto *kp_94 = buffer.data(kp + 94);
    const auto *kp_95 = buffer.data(kp + 95);
    const auto *kp_96 = buffer.data(kp + 96);
    const auto *kp_97 = buffer.data(kp + 97);
    const auto *kp_98 = buffer.data(kp + 98);
    const auto *kp_99 = buffer.data(kp + 99);
    const auto *kp_100 = buffer.data(kp + 100);
    const auto *kp_101 = buffer.data(kp + 101);
    const auto *kp_102 = buffer.data(kp + 102);
    const auto *kp_103 = buffer.data(kp + 103);
    const auto *kp_104 = buffer.data(kp + 104);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, hp_0, hp_1, hp_2, kp_3, kp_4, kp_5, \
                         kp_9, kp_10, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kp_3[k];

        t_1[k] = f_0 * kp_4[k];

        t_2[k] = f_0 * kp_5[k];

        t_3[k] = -hp_0[k]
                 + f_0 * kp_9[k];

        t_4[k] = -hp_1[k]
                 + f_0 * kp_10[k];

        t_5[k] = -hp_2[k]
                 + f_0 * kp_11[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, hp_3, hp_4, hp_5, kp_12, kp_13, \
                         kp_14, kp_18, kp_19, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * kp_12[k];

        t_7[k] = f_0 * kp_13[k];

        t_8[k] = f_0 * kp_14[k];

        t_9[k] = -2.0 * hp_3[k]
                 + f_0 * kp_18[k];

        t_10[k] = -2.0 * hp_4[k]
                  + f_0 * kp_19[k];

        t_11[k] = -2.0 * hp_5[k]
                  + f_0 * kp_20[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, hp_6, hp_7, hp_8, kp_21, kp_22, \
                         kp_23, kp_24, kp_25, kp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -hp_6[k]
                  + f_0 * kp_21[k];

        t_13[k] = -hp_7[k]
                  + f_0 * kp_22[k];

        t_14[k] = -hp_8[k]
                  + f_0 * kp_23[k];

        t_15[k] = f_0 * kp_24[k];

        t_16[k] = f_0 * kp_25[k];

        t_17[k] = f_0 * kp_26[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, hp_9, hp_10, hp_11, hp_12, hp_13, \
                         kp_30, kp_31, kp_32, kp_33, kp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -3.0 * hp_9[k]
                  + f_0 * kp_30[k];

        t_19[k] = -3.0 * hp_10[k]
                  + f_0 * kp_31[k];

        t_20[k] = -3.0 * hp_11[k]
                  + f_0 * kp_32[k];

        t_21[k] = -2.0 * hp_12[k]
                  + f_0 * kp_33[k];

        t_22[k] = -2.0 * hp_13[k]
                  + f_0 * kp_34[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, t_28, hp_14, hp_15, hp_16, hp_17, \
                         kp_35, kp_36, kp_37, kp_38, kp_39, kp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -2.0 * hp_14[k]
                  + f_0 * kp_35[k];

        t_24[k] = -hp_15[k]
                  + f_0 * kp_36[k];

        t_25[k] = -hp_16[k]
                  + f_0 * kp_37[k];

        t_26[k] = -hp_17[k]
                  + f_0 * kp_38[k];

        t_27[k] = f_0 * kp_39[k];

        t_28[k] = f_0 * kp_40[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, hp_18, hp_19, hp_20, hp_21, kp_41, \
                         kp_45, kp_46, kp_47, kp_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * kp_41[k];

        t_30[k] = -4.0 * hp_18[k]
                  + f_0 * kp_45[k];

        t_31[k] = -4.0 * hp_19[k]
                  + f_0 * kp_46[k];

        t_32[k] = -4.0 * hp_20[k]
                  + f_0 * kp_47[k];

        t_33[k] = -3.0 * hp_21[k]
                  + f_0 * kp_48[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, hp_22, hp_23, hp_24, hp_25, hp_26, \
                         kp_49, kp_50, kp_51, kp_52, kp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = -3.0 * hp_22[k]
                  + f_0 * kp_49[k];

        t_35[k] = -3.0 * hp_23[k]
                  + f_0 * kp_50[k];

        t_36[k] = -2.0 * hp_24[k]
                  + f_0 * kp_51[k];

        t_37[k] = -2.0 * hp_25[k]
                  + f_0 * kp_52[k];

        t_38[k] = -2.0 * hp_26[k]
                  + f_0 * kp_53[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, t_44, hp_27, hp_28, hp_29, kp_54, \
                         kp_55, kp_56, kp_57, kp_58, kp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -hp_27[k]
                  + f_0 * kp_54[k];

        t_40[k] = -hp_28[k]
                  + f_0 * kp_55[k];

        t_41[k] = -hp_29[k]
                  + f_0 * kp_56[k];

        t_42[k] = f_0 * kp_57[k];

        t_43[k] = f_0 * kp_58[k];

        t_44[k] = f_0 * kp_59[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, hp_30, hp_31, hp_32, hp_33, hp_34, \
                         kp_63, kp_64, kp_65, kp_66, kp_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -5.0 * hp_30[k]
                  + f_0 * kp_63[k];

        t_46[k] = -5.0 * hp_31[k]
                  + f_0 * kp_64[k];

        t_47[k] = -5.0 * hp_32[k]
                  + f_0 * kp_65[k];

        t_48[k] = -4.0 * hp_33[k]
                  + f_0 * kp_66[k];

        t_49[k] = -4.0 * hp_34[k]
                  + f_0 * kp_67[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, hp_35, hp_36, hp_37, hp_38, hp_39, \
                         kp_68, kp_69, kp_70, kp_71, kp_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -4.0 * hp_35[k]
                  + f_0 * kp_68[k];

        t_51[k] = -3.0 * hp_36[k]
                  + f_0 * kp_69[k];

        t_52[k] = -3.0 * hp_37[k]
                  + f_0 * kp_70[k];

        t_53[k] = -3.0 * hp_38[k]
                  + f_0 * kp_71[k];

        t_54[k] = -2.0 * hp_39[k]
                  + f_0 * kp_72[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, hp_40, hp_41, hp_42, hp_43, hp_44, \
                         kp_73, kp_74, kp_75, kp_76, kp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * hp_40[k]
                  + f_0 * kp_73[k];

        t_56[k] = -2.0 * hp_41[k]
                  + f_0 * kp_74[k];

        t_57[k] = -hp_42[k]
                  + f_0 * kp_75[k];

        t_58[k] = -hp_43[k]
                  + f_0 * kp_76[k];

        t_59[k] = -hp_44[k]
                  + f_0 * kp_77[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, hp_45, hp_46, hp_47, kp_78, \
                         kp_79, kp_80, kp_84, kp_85, kp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * kp_78[k];

        t_61[k] = f_0 * kp_79[k];

        t_62[k] = f_0 * kp_80[k];

        t_63[k] = -6.0 * hp_45[k]
                  + f_0 * kp_84[k];

        t_64[k] = -6.0 * hp_46[k]
                  + f_0 * kp_85[k];

        t_65[k] = -6.0 * hp_47[k]
                  + f_0 * kp_86[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, hp_48, hp_49, hp_50, hp_51, hp_52, \
                         kp_87, kp_88, kp_89, kp_90, kp_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -5.0 * hp_48[k]
                  + f_0 * kp_87[k];

        t_67[k] = -5.0 * hp_49[k]
                  + f_0 * kp_88[k];

        t_68[k] = -5.0 * hp_50[k]
                  + f_0 * kp_89[k];

        t_69[k] = -4.0 * hp_51[k]
                  + f_0 * kp_90[k];

        t_70[k] = -4.0 * hp_52[k]
                  + f_0 * kp_91[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, hp_53, hp_54, hp_55, hp_56, hp_57, \
                         kp_92, kp_93, kp_94, kp_95, kp_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -4.0 * hp_53[k]
                  + f_0 * kp_92[k];

        t_72[k] = -3.0 * hp_54[k]
                  + f_0 * kp_93[k];

        t_73[k] = -3.0 * hp_55[k]
                  + f_0 * kp_94[k];

        t_74[k] = -3.0 * hp_56[k]
                  + f_0 * kp_95[k];

        t_75[k] = -2.0 * hp_57[k]
                  + f_0 * kp_96[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, hp_58, hp_59, hp_60, hp_61, hp_62, \
                         kp_97, kp_98, kp_99, kp_100, kp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -2.0 * hp_58[k]
                  + f_0 * kp_97[k];

        t_77[k] = -2.0 * hp_59[k]
                  + f_0 * kp_98[k];

        t_78[k] = -hp_60[k]
                  + f_0 * kp_99[k];

        t_79[k] = -hp_61[k]
                  + f_0 * kp_100[k];

        t_80[k] = -hp_62[k]
                  + f_0 * kp_101[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, kp_102, kp_103, kp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_0 * kp_102[k];

        t_82[k] = f_0 * kp_103[k];

        t_83[k] = f_0 * kp_104[k];
    }
}

auto
compute_prim_geom_10_ip_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t hp, const size_t kp,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);
    const auto *hp_30 = buffer.data(hp + 30);
    const auto *hp_31 = buffer.data(hp + 31);
    const auto *hp_32 = buffer.data(hp + 32);
    const auto *hp_33 = buffer.data(hp + 33);
    const auto *hp_34 = buffer.data(hp + 34);
    const auto *hp_35 = buffer.data(hp + 35);
    const auto *hp_36 = buffer.data(hp + 36);
    const auto *hp_37 = buffer.data(hp + 37);
    const auto *hp_38 = buffer.data(hp + 38);
    const auto *hp_39 = buffer.data(hp + 39);
    const auto *hp_40 = buffer.data(hp + 40);
    const auto *hp_41 = buffer.data(hp + 41);
    const auto *hp_42 = buffer.data(hp + 42);
    const auto *hp_43 = buffer.data(hp + 43);
    const auto *hp_44 = buffer.data(hp + 44);
    const auto *hp_45 = buffer.data(hp + 45);
    const auto *hp_46 = buffer.data(hp + 46);
    const auto *hp_47 = buffer.data(hp + 47);
    const auto *hp_48 = buffer.data(hp + 48);
    const auto *hp_49 = buffer.data(hp + 49);
    const auto *hp_50 = buffer.data(hp + 50);
    const auto *hp_51 = buffer.data(hp + 51);
    const auto *hp_52 = buffer.data(hp + 52);
    const auto *hp_53 = buffer.data(hp + 53);
    const auto *hp_54 = buffer.data(hp + 54);
    const auto *hp_55 = buffer.data(hp + 55);
    const auto *hp_56 = buffer.data(hp + 56);
    const auto *hp_57 = buffer.data(hp + 57);
    const auto *hp_58 = buffer.data(hp + 58);
    const auto *hp_59 = buffer.data(hp + 59);
    const auto *hp_60 = buffer.data(hp + 60);
    const auto *hp_61 = buffer.data(hp + 61);
    const auto *hp_62 = buffer.data(hp + 62);

    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_34 = buffer.data(kp + 34);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_39 = buffer.data(kp + 39);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_48 = buffer.data(kp + 48);
    const auto *kp_49 = buffer.data(kp + 49);
    const auto *kp_50 = buffer.data(kp + 50);
    const auto *kp_51 = buffer.data(kp + 51);
    const auto *kp_52 = buffer.data(kp + 52);
    const auto *kp_53 = buffer.data(kp + 53);
    const auto *kp_54 = buffer.data(kp + 54);
    const auto *kp_55 = buffer.data(kp + 55);
    const auto *kp_56 = buffer.data(kp + 56);
    const auto *kp_57 = buffer.data(kp + 57);
    const auto *kp_58 = buffer.data(kp + 58);
    const auto *kp_59 = buffer.data(kp + 59);
    const auto *kp_60 = buffer.data(kp + 60);
    const auto *kp_61 = buffer.data(kp + 61);
    const auto *kp_62 = buffer.data(kp + 62);
    const auto *kp_66 = buffer.data(kp + 66);
    const auto *kp_67 = buffer.data(kp + 67);
    const auto *kp_68 = buffer.data(kp + 68);
    const auto *kp_69 = buffer.data(kp + 69);
    const auto *kp_70 = buffer.data(kp + 70);
    const auto *kp_71 = buffer.data(kp + 71);
    const auto *kp_72 = buffer.data(kp + 72);
    const auto *kp_73 = buffer.data(kp + 73);
    const auto *kp_74 = buffer.data(kp + 74);
    const auto *kp_75 = buffer.data(kp + 75);
    const auto *kp_76 = buffer.data(kp + 76);
    const auto *kp_77 = buffer.data(kp + 77);
    const auto *kp_78 = buffer.data(kp + 78);
    const auto *kp_79 = buffer.data(kp + 79);
    const auto *kp_80 = buffer.data(kp + 80);
    const auto *kp_81 = buffer.data(kp + 81);
    const auto *kp_82 = buffer.data(kp + 82);
    const auto *kp_83 = buffer.data(kp + 83);
    const auto *kp_87 = buffer.data(kp + 87);
    const auto *kp_88 = buffer.data(kp + 88);
    const auto *kp_89 = buffer.data(kp + 89);
    const auto *kp_90 = buffer.data(kp + 90);
    const auto *kp_91 = buffer.data(kp + 91);
    const auto *kp_92 = buffer.data(kp + 92);
    const auto *kp_93 = buffer.data(kp + 93);
    const auto *kp_94 = buffer.data(kp + 94);
    const auto *kp_95 = buffer.data(kp + 95);
    const auto *kp_96 = buffer.data(kp + 96);
    const auto *kp_97 = buffer.data(kp + 97);
    const auto *kp_98 = buffer.data(kp + 98);
    const auto *kp_99 = buffer.data(kp + 99);
    const auto *kp_100 = buffer.data(kp + 100);
    const auto *kp_101 = buffer.data(kp + 101);
    const auto *kp_102 = buffer.data(kp + 102);
    const auto *kp_103 = buffer.data(kp + 103);
    const auto *kp_104 = buffer.data(kp + 104);
    const auto *kp_105 = buffer.data(kp + 105);
    const auto *kp_106 = buffer.data(kp + 106);
    const auto *kp_107 = buffer.data(kp + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, hp_0, kp_6, kp_7, kp_8, kp_12, \
                         kp_13, kp_14, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kp_6[k];

        t_1[k] = f_0 * kp_7[k];

        t_2[k] = f_0 * kp_8[k];

        t_3[k] = f_0 * kp_12[k];

        t_4[k] = f_0 * kp_13[k];

        t_5[k] = f_0 * kp_14[k];

        t_6[k] = -hp_0[k]
                 + f_0 * kp_15[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, hp_1, hp_2, hp_3, kp_16, kp_17, \
                         kp_21, kp_22, kp_23, kp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -hp_1[k]
                 + f_0 * kp_16[k];

        t_8[k] = -hp_2[k]
                 + f_0 * kp_17[k];

        t_9[k] = f_0 * kp_21[k];

        t_10[k] = f_0 * kp_22[k];

        t_11[k] = f_0 * kp_23[k];

        t_12[k] = -hp_3[k]
                  + f_0 * kp_24[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, hp_4, hp_5, hp_6, hp_7, hp_8, kp_25, \
                         kp_26, kp_27, kp_28, kp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -hp_4[k]
                  + f_0 * kp_25[k];

        t_14[k] = -hp_5[k]
                  + f_0 * kp_26[k];

        t_15[k] = -2.0 * hp_6[k]
                  + f_0 * kp_27[k];

        t_16[k] = -2.0 * hp_7[k]
                  + f_0 * kp_28[k];

        t_17[k] = -2.0 * hp_8[k]
                  + f_0 * kp_29[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, hp_9, hp_10, hp_11, kp_33, kp_34, \
                         kp_35, kp_36, kp_37, kp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * kp_33[k];

        t_19[k] = f_0 * kp_34[k];

        t_20[k] = f_0 * kp_35[k];

        t_21[k] = -hp_9[k]
                  + f_0 * kp_36[k];

        t_22[k] = -hp_10[k]
                  + f_0 * kp_37[k];

        t_23[k] = -hp_11[k]
                  + f_0 * kp_38[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, hp_12, hp_13, hp_14, hp_15, hp_16, \
                         kp_39, kp_40, kp_41, kp_42, kp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -2.0 * hp_12[k]
                  + f_0 * kp_39[k];

        t_25[k] = -2.0 * hp_13[k]
                  + f_0 * kp_40[k];

        t_26[k] = -2.0 * hp_14[k]
                  + f_0 * kp_41[k];

        t_27[k] = -3.0 * hp_15[k]
                  + f_0 * kp_42[k];

        t_28[k] = -3.0 * hp_16[k]
                  + f_0 * kp_43[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, hp_17, hp_18, hp_19, kp_44, \
                         kp_48, kp_49, kp_50, kp_51, kp_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -3.0 * hp_17[k]
                  + f_0 * kp_44[k];

        t_30[k] = f_0 * kp_48[k];

        t_31[k] = f_0 * kp_49[k];

        t_32[k] = f_0 * kp_50[k];

        t_33[k] = -hp_18[k]
                  + f_0 * kp_51[k];

        t_34[k] = -hp_19[k]
                  + f_0 * kp_52[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, hp_20, hp_21, hp_22, hp_23, hp_24, \
                         kp_53, kp_54, kp_55, kp_56, kp_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -hp_20[k]
                  + f_0 * kp_53[k];

        t_36[k] = -2.0 * hp_21[k]
                  + f_0 * kp_54[k];

        t_37[k] = -2.0 * hp_22[k]
                  + f_0 * kp_55[k];

        t_38[k] = -2.0 * hp_23[k]
                  + f_0 * kp_56[k];

        t_39[k] = -3.0 * hp_24[k]
                  + f_0 * kp_57[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, hp_25, hp_26, hp_27, hp_28, hp_29, \
                         kp_58, kp_59, kp_60, kp_61, kp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -3.0 * hp_25[k]
                  + f_0 * kp_58[k];

        t_41[k] = -3.0 * hp_26[k]
                  + f_0 * kp_59[k];

        t_42[k] = -4.0 * hp_27[k]
                  + f_0 * kp_60[k];

        t_43[k] = -4.0 * hp_28[k]
                  + f_0 * kp_61[k];

        t_44[k] = -4.0 * hp_29[k]
                  + f_0 * kp_62[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, hp_30, hp_31, hp_32, kp_66, \
                         kp_67, kp_68, kp_69, kp_70, kp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_0 * kp_66[k];

        t_46[k] = f_0 * kp_67[k];

        t_47[k] = f_0 * kp_68[k];

        t_48[k] = -hp_30[k]
                  + f_0 * kp_69[k];

        t_49[k] = -hp_31[k]
                  + f_0 * kp_70[k];

        t_50[k] = -hp_32[k]
                  + f_0 * kp_71[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, hp_33, hp_34, hp_35, hp_36, hp_37, \
                         kp_72, kp_73, kp_74, kp_75, kp_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -2.0 * hp_33[k]
                  + f_0 * kp_72[k];

        t_52[k] = -2.0 * hp_34[k]
                  + f_0 * kp_73[k];

        t_53[k] = -2.0 * hp_35[k]
                  + f_0 * kp_74[k];

        t_54[k] = -3.0 * hp_36[k]
                  + f_0 * kp_75[k];

        t_55[k] = -3.0 * hp_37[k]
                  + f_0 * kp_76[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, hp_38, hp_39, hp_40, hp_41, hp_42, \
                         kp_77, kp_78, kp_79, kp_80, kp_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -3.0 * hp_38[k]
                  + f_0 * kp_77[k];

        t_57[k] = -4.0 * hp_39[k]
                  + f_0 * kp_78[k];

        t_58[k] = -4.0 * hp_40[k]
                  + f_0 * kp_79[k];

        t_59[k] = -4.0 * hp_41[k]
                  + f_0 * kp_80[k];

        t_60[k] = -5.0 * hp_42[k]
                  + f_0 * kp_81[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, hp_43, hp_44, hp_45, kp_82, \
                         kp_83, kp_87, kp_88, kp_89, kp_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -5.0 * hp_43[k]
                  + f_0 * kp_82[k];

        t_62[k] = -5.0 * hp_44[k]
                  + f_0 * kp_83[k];

        t_63[k] = f_0 * kp_87[k];

        t_64[k] = f_0 * kp_88[k];

        t_65[k] = f_0 * kp_89[k];

        t_66[k] = -hp_45[k]
                  + f_0 * kp_90[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, hp_46, hp_47, hp_48, hp_49, hp_50, \
                         kp_91, kp_92, kp_93, kp_94, kp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -hp_46[k]
                  + f_0 * kp_91[k];

        t_68[k] = -hp_47[k]
                  + f_0 * kp_92[k];

        t_69[k] = -2.0 * hp_48[k]
                  + f_0 * kp_93[k];

        t_70[k] = -2.0 * hp_49[k]
                  + f_0 * kp_94[k];

        t_71[k] = -2.0 * hp_50[k]
                  + f_0 * kp_95[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, hp_51, hp_52, hp_53, hp_54, hp_55, \
                         kp_96, kp_97, kp_98, kp_99, kp_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -3.0 * hp_51[k]
                  + f_0 * kp_96[k];

        t_73[k] = -3.0 * hp_52[k]
                  + f_0 * kp_97[k];

        t_74[k] = -3.0 * hp_53[k]
                  + f_0 * kp_98[k];

        t_75[k] = -4.0 * hp_54[k]
                  + f_0 * kp_99[k];

        t_76[k] = -4.0 * hp_55[k]
                  + f_0 * kp_100[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, hp_56, hp_57, hp_58, hp_59, hp_60, \
                         kp_101, kp_102, kp_103, kp_104, kp_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -4.0 * hp_56[k]
                  + f_0 * kp_101[k];

        t_78[k] = -5.0 * hp_57[k]
                  + f_0 * kp_102[k];

        t_79[k] = -5.0 * hp_58[k]
                  + f_0 * kp_103[k];

        t_80[k] = -5.0 * hp_59[k]
                  + f_0 * kp_104[k];

        t_81[k] = -6.0 * hp_60[k]
                  + f_0 * kp_105[k];
    }

#pragma omp simd aligned(t_82, t_83, hp_61, hp_62, kp_106, kp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -6.0 * hp_61[k]
                  + f_0 * kp_106[k];

        t_83[k] = -6.0 * hp_62[k]
                  + f_0 * kp_107[k];
    }
}

}  // namespace simdt2ceri
