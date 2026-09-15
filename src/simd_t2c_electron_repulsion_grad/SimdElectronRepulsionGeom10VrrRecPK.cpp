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


#include "SimdElectronRepulsionGeom10VrrRecPK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_pk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t sk, const size_t dk,
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
    const auto *sk_18 = buffer.data(sk + 18);
    const auto *sk_19 = buffer.data(sk + 19);
    const auto *sk_20 = buffer.data(sk + 20);
    const auto *sk_21 = buffer.data(sk + 21);
    const auto *sk_22 = buffer.data(sk + 22);
    const auto *sk_23 = buffer.data(sk + 23);
    const auto *sk_24 = buffer.data(sk + 24);
    const auto *sk_25 = buffer.data(sk + 25);
    const auto *sk_26 = buffer.data(sk + 26);
    const auto *sk_27 = buffer.data(sk + 27);
    const auto *sk_28 = buffer.data(sk + 28);
    const auto *sk_29 = buffer.data(sk + 29);
    const auto *sk_30 = buffer.data(sk + 30);
    const auto *sk_31 = buffer.data(sk + 31);
    const auto *sk_32 = buffer.data(sk + 32);
    const auto *sk_33 = buffer.data(sk + 33);
    const auto *sk_34 = buffer.data(sk + 34);
    const auto *sk_35 = buffer.data(sk + 35);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_4 = buffer.data(dk + 4);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_6 = buffer.data(dk + 6);
    const auto *dk_7 = buffer.data(dk + 7);
    const auto *dk_8 = buffer.data(dk + 8);
    const auto *dk_9 = buffer.data(dk + 9);
    const auto *dk_10 = buffer.data(dk + 10);
    const auto *dk_11 = buffer.data(dk + 11);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_13 = buffer.data(dk + 13);
    const auto *dk_14 = buffer.data(dk + 14);
    const auto *dk_15 = buffer.data(dk + 15);
    const auto *dk_16 = buffer.data(dk + 16);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_18 = buffer.data(dk + 18);
    const auto *dk_19 = buffer.data(dk + 19);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_22 = buffer.data(dk + 22);
    const auto *dk_23 = buffer.data(dk + 23);
    const auto *dk_24 = buffer.data(dk + 24);
    const auto *dk_25 = buffer.data(dk + 25);
    const auto *dk_26 = buffer.data(dk + 26);
    const auto *dk_27 = buffer.data(dk + 27);
    const auto *dk_28 = buffer.data(dk + 28);
    const auto *dk_29 = buffer.data(dk + 29);
    const auto *dk_30 = buffer.data(dk + 30);
    const auto *dk_31 = buffer.data(dk + 31);
    const auto *dk_32 = buffer.data(dk + 32);
    const auto *dk_33 = buffer.data(dk + 33);
    const auto *dk_34 = buffer.data(dk + 34);
    const auto *dk_35 = buffer.data(dk + 35);
    const auto *dk_36 = buffer.data(dk + 36);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_38 = buffer.data(dk + 38);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_40 = buffer.data(dk + 40);
    const auto *dk_41 = buffer.data(dk + 41);
    const auto *dk_42 = buffer.data(dk + 42);
    const auto *dk_43 = buffer.data(dk + 43);
    const auto *dk_44 = buffer.data(dk + 44);
    const auto *dk_45 = buffer.data(dk + 45);
    const auto *dk_46 = buffer.data(dk + 46);
    const auto *dk_47 = buffer.data(dk + 47);
    const auto *dk_48 = buffer.data(dk + 48);
    const auto *dk_49 = buffer.data(dk + 49);
    const auto *dk_50 = buffer.data(dk + 50);
    const auto *dk_51 = buffer.data(dk + 51);
    const auto *dk_52 = buffer.data(dk + 52);
    const auto *dk_53 = buffer.data(dk + 53);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_55 = buffer.data(dk + 55);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_58 = buffer.data(dk + 58);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_60 = buffer.data(dk + 60);
    const auto *dk_61 = buffer.data(dk + 61);
    const auto *dk_62 = buffer.data(dk + 62);
    const auto *dk_63 = buffer.data(dk + 63);
    const auto *dk_64 = buffer.data(dk + 64);
    const auto *dk_65 = buffer.data(dk + 65);
    const auto *dk_66 = buffer.data(dk + 66);
    const auto *dk_67 = buffer.data(dk + 67);
    const auto *dk_68 = buffer.data(dk + 68);
    const auto *dk_69 = buffer.data(dk + 69);
    const auto *dk_70 = buffer.data(dk + 70);
    const auto *dk_71 = buffer.data(dk + 71);
    const auto *dk_72 = buffer.data(dk + 72);
    const auto *dk_73 = buffer.data(dk + 73);
    const auto *dk_74 = buffer.data(dk + 74);
    const auto *dk_75 = buffer.data(dk + 75);
    const auto *dk_76 = buffer.data(dk + 76);
    const auto *dk_77 = buffer.data(dk + 77);
    const auto *dk_78 = buffer.data(dk + 78);
    const auto *dk_79 = buffer.data(dk + 79);
    const auto *dk_80 = buffer.data(dk + 80);
    const auto *dk_81 = buffer.data(dk + 81);
    const auto *dk_82 = buffer.data(dk + 82);
    const auto *dk_83 = buffer.data(dk + 83);
    const auto *dk_84 = buffer.data(dk + 84);
    const auto *dk_85 = buffer.data(dk + 85);
    const auto *dk_86 = buffer.data(dk + 86);
    const auto *dk_87 = buffer.data(dk + 87);
    const auto *dk_88 = buffer.data(dk + 88);
    const auto *dk_89 = buffer.data(dk + 89);
    const auto *dk_90 = buffer.data(dk + 90);
    const auto *dk_91 = buffer.data(dk + 91);
    const auto *dk_92 = buffer.data(dk + 92);
    const auto *dk_93 = buffer.data(dk + 93);
    const auto *dk_94 = buffer.data(dk + 94);
    const auto *dk_95 = buffer.data(dk + 95);
    const auto *dk_96 = buffer.data(dk + 96);
    const auto *dk_97 = buffer.data(dk + 97);
    const auto *dk_98 = buffer.data(dk + 98);
    const auto *dk_99 = buffer.data(dk + 99);
    const auto *dk_100 = buffer.data(dk + 100);
    const auto *dk_101 = buffer.data(dk + 101);
    const auto *dk_102 = buffer.data(dk + 102);
    const auto *dk_103 = buffer.data(dk + 103);
    const auto *dk_104 = buffer.data(dk + 104);
    const auto *dk_105 = buffer.data(dk + 105);
    const auto *dk_106 = buffer.data(dk + 106);
    const auto *dk_107 = buffer.data(dk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, sk_0, sk_1, sk_2, sk_3, sk_4, dk_0, dk_1, \
                         dk_2, dk_3, dk_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -sk_0[k]
                 + f_0 * dk_0[k];

        t_1[k] = -sk_1[k]
                 + f_0 * dk_1[k];

        t_2[k] = -sk_2[k]
                 + f_0 * dk_2[k];

        t_3[k] = -sk_3[k]
                 + f_0 * dk_3[k];

        t_4[k] = -sk_4[k]
                 + f_0 * dk_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, sk_5, sk_6, sk_7, sk_8, sk_9, dk_5, dk_6, \
                         dk_7, dk_8, dk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -sk_5[k]
                 + f_0 * dk_5[k];

        t_6[k] = -sk_6[k]
                 + f_0 * dk_6[k];

        t_7[k] = -sk_7[k]
                 + f_0 * dk_7[k];

        t_8[k] = -sk_8[k]
                 + f_0 * dk_8[k];

        t_9[k] = -sk_9[k]
                 + f_0 * dk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, sk_10, sk_11, sk_12, sk_13, sk_14, \
                         dk_10, dk_11, dk_12, dk_13, dk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -sk_10[k]
                  + f_0 * dk_10[k];

        t_11[k] = -sk_11[k]
                  + f_0 * dk_11[k];

        t_12[k] = -sk_12[k]
                  + f_0 * dk_12[k];

        t_13[k] = -sk_13[k]
                  + f_0 * dk_13[k];

        t_14[k] = -sk_14[k]
                  + f_0 * dk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, sk_15, sk_16, sk_17, sk_18, sk_19, \
                         dk_15, dk_16, dk_17, dk_18, dk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -sk_15[k]
                  + f_0 * dk_15[k];

        t_16[k] = -sk_16[k]
                  + f_0 * dk_16[k];

        t_17[k] = -sk_17[k]
                  + f_0 * dk_17[k];

        t_18[k] = -sk_18[k]
                  + f_0 * dk_18[k];

        t_19[k] = -sk_19[k]
                  + f_0 * dk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, sk_20, sk_21, sk_22, sk_23, sk_24, \
                         dk_20, dk_21, dk_22, dk_23, dk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -sk_20[k]
                  + f_0 * dk_20[k];

        t_21[k] = -sk_21[k]
                  + f_0 * dk_21[k];

        t_22[k] = -sk_22[k]
                  + f_0 * dk_22[k];

        t_23[k] = -sk_23[k]
                  + f_0 * dk_23[k];

        t_24[k] = -sk_24[k]
                  + f_0 * dk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, sk_25, sk_26, sk_27, sk_28, sk_29, \
                         dk_25, dk_26, dk_27, dk_28, dk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -sk_25[k]
                  + f_0 * dk_25[k];

        t_26[k] = -sk_26[k]
                  + f_0 * dk_26[k];

        t_27[k] = -sk_27[k]
                  + f_0 * dk_27[k];

        t_28[k] = -sk_28[k]
                  + f_0 * dk_28[k];

        t_29[k] = -sk_29[k]
                  + f_0 * dk_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, sk_30, sk_31, sk_32, sk_33, sk_34, \
                         dk_30, dk_31, dk_32, dk_33, dk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -sk_30[k]
                  + f_0 * dk_30[k];

        t_31[k] = -sk_31[k]
                  + f_0 * dk_31[k];

        t_32[k] = -sk_32[k]
                  + f_0 * dk_32[k];

        t_33[k] = -sk_33[k]
                  + f_0 * dk_33[k];

        t_34[k] = -sk_34[k]
                  + f_0 * dk_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, t_40, t_41, sk_35, dk_35, dk_36, dk_37, \
                         dk_38, dk_39, dk_40, dk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -sk_35[k]
                  + f_0 * dk_35[k];

        t_36[k] = f_0 * dk_36[k];

        t_37[k] = f_0 * dk_37[k];

        t_38[k] = f_0 * dk_38[k];

        t_39[k] = f_0 * dk_39[k];

        t_40[k] = f_0 * dk_40[k];

        t_41[k] = f_0 * dk_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, t_48, t_49, dk_42, dk_43, dk_44, \
                         dk_45, dk_46, dk_47, dk_48, dk_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_0 * dk_42[k];

        t_43[k] = f_0 * dk_43[k];

        t_44[k] = f_0 * dk_44[k];

        t_45[k] = f_0 * dk_45[k];

        t_46[k] = f_0 * dk_46[k];

        t_47[k] = f_0 * dk_47[k];

        t_48[k] = f_0 * dk_48[k];

        t_49[k] = f_0 * dk_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, t_55, t_56, t_57, dk_50, dk_51, dk_52, \
                         dk_53, dk_54, dk_55, dk_56, dk_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_0 * dk_50[k];

        t_51[k] = f_0 * dk_51[k];

        t_52[k] = f_0 * dk_52[k];

        t_53[k] = f_0 * dk_53[k];

        t_54[k] = f_0 * dk_54[k];

        t_55[k] = f_0 * dk_55[k];

        t_56[k] = f_0 * dk_56[k];

        t_57[k] = f_0 * dk_57[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, t_63, t_64, t_65, dk_58, dk_59, dk_60, \
                         dk_61, dk_62, dk_63, dk_64, dk_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_0 * dk_58[k];

        t_59[k] = f_0 * dk_59[k];

        t_60[k] = f_0 * dk_60[k];

        t_61[k] = f_0 * dk_61[k];

        t_62[k] = f_0 * dk_62[k];

        t_63[k] = f_0 * dk_63[k];

        t_64[k] = f_0 * dk_64[k];

        t_65[k] = f_0 * dk_65[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, t_72, t_73, dk_66, dk_67, dk_68, \
                         dk_69, dk_70, dk_71, dk_72, dk_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_0 * dk_66[k];

        t_67[k] = f_0 * dk_67[k];

        t_68[k] = f_0 * dk_68[k];

        t_69[k] = f_0 * dk_69[k];

        t_70[k] = f_0 * dk_70[k];

        t_71[k] = f_0 * dk_71[k];

        t_72[k] = f_0 * dk_72[k];

        t_73[k] = f_0 * dk_73[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, t_80, t_81, dk_74, dk_75, dk_76, \
                         dk_77, dk_78, dk_79, dk_80, dk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * dk_74[k];

        t_75[k] = f_0 * dk_75[k];

        t_76[k] = f_0 * dk_76[k];

        t_77[k] = f_0 * dk_77[k];

        t_78[k] = f_0 * dk_78[k];

        t_79[k] = f_0 * dk_79[k];

        t_80[k] = f_0 * dk_80[k];

        t_81[k] = f_0 * dk_81[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, t_87, t_88, t_89, dk_82, dk_83, dk_84, \
                         dk_85, dk_86, dk_87, dk_88, dk_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * dk_82[k];

        t_83[k] = f_0 * dk_83[k];

        t_84[k] = f_0 * dk_84[k];

        t_85[k] = f_0 * dk_85[k];

        t_86[k] = f_0 * dk_86[k];

        t_87[k] = f_0 * dk_87[k];

        t_88[k] = f_0 * dk_88[k];

        t_89[k] = f_0 * dk_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, dk_90, dk_91, dk_92, \
                         dk_93, dk_94, dk_95, dk_96, dk_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * dk_90[k];

        t_91[k] = f_0 * dk_91[k];

        t_92[k] = f_0 * dk_92[k];

        t_93[k] = f_0 * dk_93[k];

        t_94[k] = f_0 * dk_94[k];

        t_95[k] = f_0 * dk_95[k];

        t_96[k] = f_0 * dk_96[k];

        t_97[k] = f_0 * dk_97[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, t_104, t_105, dk_98, dk_99, \
                         dk_100, dk_101, dk_102, dk_103, dk_104, \
                         dk_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * dk_98[k];

        t_99[k] = f_0 * dk_99[k];

        t_100[k] = f_0 * dk_100[k];

        t_101[k] = f_0 * dk_101[k];

        t_102[k] = f_0 * dk_102[k];

        t_103[k] = f_0 * dk_103[k];

        t_104[k] = f_0 * dk_104[k];

        t_105[k] = f_0 * dk_105[k];
    }

#pragma omp simd aligned(t_106, t_107, dk_106, dk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * dk_106[k];

        t_107[k] = f_0 * dk_107[k];
    }
}

auto
compute_prim_geom_10_pk_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t sk, const size_t dk,
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
    const auto *sk_18 = buffer.data(sk + 18);
    const auto *sk_19 = buffer.data(sk + 19);
    const auto *sk_20 = buffer.data(sk + 20);
    const auto *sk_21 = buffer.data(sk + 21);
    const auto *sk_22 = buffer.data(sk + 22);
    const auto *sk_23 = buffer.data(sk + 23);
    const auto *sk_24 = buffer.data(sk + 24);
    const auto *sk_25 = buffer.data(sk + 25);
    const auto *sk_26 = buffer.data(sk + 26);
    const auto *sk_27 = buffer.data(sk + 27);
    const auto *sk_28 = buffer.data(sk + 28);
    const auto *sk_29 = buffer.data(sk + 29);
    const auto *sk_30 = buffer.data(sk + 30);
    const auto *sk_31 = buffer.data(sk + 31);
    const auto *sk_32 = buffer.data(sk + 32);
    const auto *sk_33 = buffer.data(sk + 33);
    const auto *sk_34 = buffer.data(sk + 34);
    const auto *sk_35 = buffer.data(sk + 35);

    const auto *dk_36 = buffer.data(dk + 36);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_38 = buffer.data(dk + 38);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_40 = buffer.data(dk + 40);
    const auto *dk_41 = buffer.data(dk + 41);
    const auto *dk_42 = buffer.data(dk + 42);
    const auto *dk_43 = buffer.data(dk + 43);
    const auto *dk_44 = buffer.data(dk + 44);
    const auto *dk_45 = buffer.data(dk + 45);
    const auto *dk_46 = buffer.data(dk + 46);
    const auto *dk_47 = buffer.data(dk + 47);
    const auto *dk_48 = buffer.data(dk + 48);
    const auto *dk_49 = buffer.data(dk + 49);
    const auto *dk_50 = buffer.data(dk + 50);
    const auto *dk_51 = buffer.data(dk + 51);
    const auto *dk_52 = buffer.data(dk + 52);
    const auto *dk_53 = buffer.data(dk + 53);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_55 = buffer.data(dk + 55);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_58 = buffer.data(dk + 58);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_60 = buffer.data(dk + 60);
    const auto *dk_61 = buffer.data(dk + 61);
    const auto *dk_62 = buffer.data(dk + 62);
    const auto *dk_63 = buffer.data(dk + 63);
    const auto *dk_64 = buffer.data(dk + 64);
    const auto *dk_65 = buffer.data(dk + 65);
    const auto *dk_66 = buffer.data(dk + 66);
    const auto *dk_67 = buffer.data(dk + 67);
    const auto *dk_68 = buffer.data(dk + 68);
    const auto *dk_69 = buffer.data(dk + 69);
    const auto *dk_70 = buffer.data(dk + 70);
    const auto *dk_71 = buffer.data(dk + 71);
    const auto *dk_108 = buffer.data(dk + 108);
    const auto *dk_109 = buffer.data(dk + 109);
    const auto *dk_110 = buffer.data(dk + 110);
    const auto *dk_111 = buffer.data(dk + 111);
    const auto *dk_112 = buffer.data(dk + 112);
    const auto *dk_113 = buffer.data(dk + 113);
    const auto *dk_114 = buffer.data(dk + 114);
    const auto *dk_115 = buffer.data(dk + 115);
    const auto *dk_116 = buffer.data(dk + 116);
    const auto *dk_117 = buffer.data(dk + 117);
    const auto *dk_118 = buffer.data(dk + 118);
    const auto *dk_119 = buffer.data(dk + 119);
    const auto *dk_120 = buffer.data(dk + 120);
    const auto *dk_121 = buffer.data(dk + 121);
    const auto *dk_122 = buffer.data(dk + 122);
    const auto *dk_123 = buffer.data(dk + 123);
    const auto *dk_124 = buffer.data(dk + 124);
    const auto *dk_125 = buffer.data(dk + 125);
    const auto *dk_126 = buffer.data(dk + 126);
    const auto *dk_127 = buffer.data(dk + 127);
    const auto *dk_128 = buffer.data(dk + 128);
    const auto *dk_129 = buffer.data(dk + 129);
    const auto *dk_130 = buffer.data(dk + 130);
    const auto *dk_131 = buffer.data(dk + 131);
    const auto *dk_132 = buffer.data(dk + 132);
    const auto *dk_133 = buffer.data(dk + 133);
    const auto *dk_134 = buffer.data(dk + 134);
    const auto *dk_135 = buffer.data(dk + 135);
    const auto *dk_136 = buffer.data(dk + 136);
    const auto *dk_137 = buffer.data(dk + 137);
    const auto *dk_138 = buffer.data(dk + 138);
    const auto *dk_139 = buffer.data(dk + 139);
    const auto *dk_140 = buffer.data(dk + 140);
    const auto *dk_141 = buffer.data(dk + 141);
    const auto *dk_142 = buffer.data(dk + 142);
    const auto *dk_143 = buffer.data(dk + 143);
    const auto *dk_144 = buffer.data(dk + 144);
    const auto *dk_145 = buffer.data(dk + 145);
    const auto *dk_146 = buffer.data(dk + 146);
    const auto *dk_147 = buffer.data(dk + 147);
    const auto *dk_148 = buffer.data(dk + 148);
    const auto *dk_149 = buffer.data(dk + 149);
    const auto *dk_150 = buffer.data(dk + 150);
    const auto *dk_151 = buffer.data(dk + 151);
    const auto *dk_152 = buffer.data(dk + 152);
    const auto *dk_153 = buffer.data(dk + 153);
    const auto *dk_154 = buffer.data(dk + 154);
    const auto *dk_155 = buffer.data(dk + 155);
    const auto *dk_156 = buffer.data(dk + 156);
    const auto *dk_157 = buffer.data(dk + 157);
    const auto *dk_158 = buffer.data(dk + 158);
    const auto *dk_159 = buffer.data(dk + 159);
    const auto *dk_160 = buffer.data(dk + 160);
    const auto *dk_161 = buffer.data(dk + 161);
    const auto *dk_162 = buffer.data(dk + 162);
    const auto *dk_163 = buffer.data(dk + 163);
    const auto *dk_164 = buffer.data(dk + 164);
    const auto *dk_165 = buffer.data(dk + 165);
    const auto *dk_166 = buffer.data(dk + 166);
    const auto *dk_167 = buffer.data(dk + 167);
    const auto *dk_168 = buffer.data(dk + 168);
    const auto *dk_169 = buffer.data(dk + 169);
    const auto *dk_170 = buffer.data(dk + 170);
    const auto *dk_171 = buffer.data(dk + 171);
    const auto *dk_172 = buffer.data(dk + 172);
    const auto *dk_173 = buffer.data(dk + 173);
    const auto *dk_174 = buffer.data(dk + 174);
    const auto *dk_175 = buffer.data(dk + 175);
    const auto *dk_176 = buffer.data(dk + 176);
    const auto *dk_177 = buffer.data(dk + 177);
    const auto *dk_178 = buffer.data(dk + 178);
    const auto *dk_179 = buffer.data(dk + 179);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, dk_36, dk_37, dk_38, dk_39, \
                         dk_40, dk_41, dk_42, dk_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_36[k];

        t_1[k] = f_0 * dk_37[k];

        t_2[k] = f_0 * dk_38[k];

        t_3[k] = f_0 * dk_39[k];

        t_4[k] = f_0 * dk_40[k];

        t_5[k] = f_0 * dk_41[k];

        t_6[k] = f_0 * dk_42[k];

        t_7[k] = f_0 * dk_43[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, dk_44, dk_45, dk_46, \
                         dk_47, dk_48, dk_49, dk_50, dk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dk_44[k];

        t_9[k] = f_0 * dk_45[k];

        t_10[k] = f_0 * dk_46[k];

        t_11[k] = f_0 * dk_47[k];

        t_12[k] = f_0 * dk_48[k];

        t_13[k] = f_0 * dk_49[k];

        t_14[k] = f_0 * dk_50[k];

        t_15[k] = f_0 * dk_51[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, dk_52, dk_53, dk_54, \
                         dk_55, dk_56, dk_57, dk_58, dk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * dk_52[k];

        t_17[k] = f_0 * dk_53[k];

        t_18[k] = f_0 * dk_54[k];

        t_19[k] = f_0 * dk_55[k];

        t_20[k] = f_0 * dk_56[k];

        t_21[k] = f_0 * dk_57[k];

        t_22[k] = f_0 * dk_58[k];

        t_23[k] = f_0 * dk_59[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, dk_60, dk_61, dk_62, \
                         dk_63, dk_64, dk_65, dk_66, dk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * dk_60[k];

        t_25[k] = f_0 * dk_61[k];

        t_26[k] = f_0 * dk_62[k];

        t_27[k] = f_0 * dk_63[k];

        t_28[k] = f_0 * dk_64[k];

        t_29[k] = f_0 * dk_65[k];

        t_30[k] = f_0 * dk_66[k];

        t_31[k] = f_0 * dk_67[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, sk_0, sk_1, dk_68, dk_69, dk_70, \
                         dk_71, dk_108, dk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * dk_68[k];

        t_33[k] = f_0 * dk_69[k];

        t_34[k] = f_0 * dk_70[k];

        t_35[k] = f_0 * dk_71[k];

        t_36[k] = -sk_0[k]
                  + f_0 * dk_108[k];

        t_37[k] = -sk_1[k]
                  + f_0 * dk_109[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, sk_2, sk_3, sk_4, sk_5, sk_6, dk_110, \
                         dk_111, dk_112, dk_113, dk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -sk_2[k]
                  + f_0 * dk_110[k];

        t_39[k] = -sk_3[k]
                  + f_0 * dk_111[k];

        t_40[k] = -sk_4[k]
                  + f_0 * dk_112[k];

        t_41[k] = -sk_5[k]
                  + f_0 * dk_113[k];

        t_42[k] = -sk_6[k]
                  + f_0 * dk_114[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, sk_7, sk_8, sk_9, sk_10, sk_11, dk_115, \
                         dk_116, dk_117, dk_118, dk_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -sk_7[k]
                  + f_0 * dk_115[k];

        t_44[k] = -sk_8[k]
                  + f_0 * dk_116[k];

        t_45[k] = -sk_9[k]
                  + f_0 * dk_117[k];

        t_46[k] = -sk_10[k]
                  + f_0 * dk_118[k];

        t_47[k] = -sk_11[k]
                  + f_0 * dk_119[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, sk_12, sk_13, sk_14, sk_15, sk_16, \
                         dk_120, dk_121, dk_122, dk_123, dk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -sk_12[k]
                  + f_0 * dk_120[k];

        t_49[k] = -sk_13[k]
                  + f_0 * dk_121[k];

        t_50[k] = -sk_14[k]
                  + f_0 * dk_122[k];

        t_51[k] = -sk_15[k]
                  + f_0 * dk_123[k];

        t_52[k] = -sk_16[k]
                  + f_0 * dk_124[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, sk_17, sk_18, sk_19, sk_20, sk_21, \
                         dk_125, dk_126, dk_127, dk_128, dk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -sk_17[k]
                  + f_0 * dk_125[k];

        t_54[k] = -sk_18[k]
                  + f_0 * dk_126[k];

        t_55[k] = -sk_19[k]
                  + f_0 * dk_127[k];

        t_56[k] = -sk_20[k]
                  + f_0 * dk_128[k];

        t_57[k] = -sk_21[k]
                  + f_0 * dk_129[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, sk_22, sk_23, sk_24, sk_25, sk_26, \
                         dk_130, dk_131, dk_132, dk_133, dk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -sk_22[k]
                  + f_0 * dk_130[k];

        t_59[k] = -sk_23[k]
                  + f_0 * dk_131[k];

        t_60[k] = -sk_24[k]
                  + f_0 * dk_132[k];

        t_61[k] = -sk_25[k]
                  + f_0 * dk_133[k];

        t_62[k] = -sk_26[k]
                  + f_0 * dk_134[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, sk_27, sk_28, sk_29, sk_30, sk_31, \
                         dk_135, dk_136, dk_137, dk_138, dk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -sk_27[k]
                  + f_0 * dk_135[k];

        t_64[k] = -sk_28[k]
                  + f_0 * dk_136[k];

        t_65[k] = -sk_29[k]
                  + f_0 * dk_137[k];

        t_66[k] = -sk_30[k]
                  + f_0 * dk_138[k];

        t_67[k] = -sk_31[k]
                  + f_0 * dk_139[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, sk_32, sk_33, sk_34, sk_35, \
                         dk_140, dk_141, dk_142, dk_143, dk_144, \
                         dk_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -sk_32[k]
                  + f_0 * dk_140[k];

        t_69[k] = -sk_33[k]
                  + f_0 * dk_141[k];

        t_70[k] = -sk_34[k]
                  + f_0 * dk_142[k];

        t_71[k] = -sk_35[k]
                  + f_0 * dk_143[k];

        t_72[k] = f_0 * dk_144[k];

        t_73[k] = f_0 * dk_145[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, t_80, t_81, dk_146, dk_147, \
                         dk_148, dk_149, dk_150, dk_151, dk_152, \
                         dk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * dk_146[k];

        t_75[k] = f_0 * dk_147[k];

        t_76[k] = f_0 * dk_148[k];

        t_77[k] = f_0 * dk_149[k];

        t_78[k] = f_0 * dk_150[k];

        t_79[k] = f_0 * dk_151[k];

        t_80[k] = f_0 * dk_152[k];

        t_81[k] = f_0 * dk_153[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, t_87, t_88, t_89, dk_154, dk_155, \
                         dk_156, dk_157, dk_158, dk_159, dk_160, \
                         dk_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * dk_154[k];

        t_83[k] = f_0 * dk_155[k];

        t_84[k] = f_0 * dk_156[k];

        t_85[k] = f_0 * dk_157[k];

        t_86[k] = f_0 * dk_158[k];

        t_87[k] = f_0 * dk_159[k];

        t_88[k] = f_0 * dk_160[k];

        t_89[k] = f_0 * dk_161[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, dk_162, dk_163, \
                         dk_164, dk_165, dk_166, dk_167, dk_168, \
                         dk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * dk_162[k];

        t_91[k] = f_0 * dk_163[k];

        t_92[k] = f_0 * dk_164[k];

        t_93[k] = f_0 * dk_165[k];

        t_94[k] = f_0 * dk_166[k];

        t_95[k] = f_0 * dk_167[k];

        t_96[k] = f_0 * dk_168[k];

        t_97[k] = f_0 * dk_169[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, t_104, t_105, dk_170, dk_171, \
                         dk_172, dk_173, dk_174, dk_175, dk_176, \
                         dk_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * dk_170[k];

        t_99[k] = f_0 * dk_171[k];

        t_100[k] = f_0 * dk_172[k];

        t_101[k] = f_0 * dk_173[k];

        t_102[k] = f_0 * dk_174[k];

        t_103[k] = f_0 * dk_175[k];

        t_104[k] = f_0 * dk_176[k];

        t_105[k] = f_0 * dk_177[k];
    }

#pragma omp simd aligned(t_106, t_107, dk_178, dk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * dk_178[k];

        t_107[k] = f_0 * dk_179[k];
    }
}

auto
compute_prim_geom_10_pk_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t sk, const size_t dk,
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
    const auto *sk_18 = buffer.data(sk + 18);
    const auto *sk_19 = buffer.data(sk + 19);
    const auto *sk_20 = buffer.data(sk + 20);
    const auto *sk_21 = buffer.data(sk + 21);
    const auto *sk_22 = buffer.data(sk + 22);
    const auto *sk_23 = buffer.data(sk + 23);
    const auto *sk_24 = buffer.data(sk + 24);
    const auto *sk_25 = buffer.data(sk + 25);
    const auto *sk_26 = buffer.data(sk + 26);
    const auto *sk_27 = buffer.data(sk + 27);
    const auto *sk_28 = buffer.data(sk + 28);
    const auto *sk_29 = buffer.data(sk + 29);
    const auto *sk_30 = buffer.data(sk + 30);
    const auto *sk_31 = buffer.data(sk + 31);
    const auto *sk_32 = buffer.data(sk + 32);
    const auto *sk_33 = buffer.data(sk + 33);
    const auto *sk_34 = buffer.data(sk + 34);
    const auto *sk_35 = buffer.data(sk + 35);

    const auto *dk_72 = buffer.data(dk + 72);
    const auto *dk_73 = buffer.data(dk + 73);
    const auto *dk_74 = buffer.data(dk + 74);
    const auto *dk_75 = buffer.data(dk + 75);
    const auto *dk_76 = buffer.data(dk + 76);
    const auto *dk_77 = buffer.data(dk + 77);
    const auto *dk_78 = buffer.data(dk + 78);
    const auto *dk_79 = buffer.data(dk + 79);
    const auto *dk_80 = buffer.data(dk + 80);
    const auto *dk_81 = buffer.data(dk + 81);
    const auto *dk_82 = buffer.data(dk + 82);
    const auto *dk_83 = buffer.data(dk + 83);
    const auto *dk_84 = buffer.data(dk + 84);
    const auto *dk_85 = buffer.data(dk + 85);
    const auto *dk_86 = buffer.data(dk + 86);
    const auto *dk_87 = buffer.data(dk + 87);
    const auto *dk_88 = buffer.data(dk + 88);
    const auto *dk_89 = buffer.data(dk + 89);
    const auto *dk_90 = buffer.data(dk + 90);
    const auto *dk_91 = buffer.data(dk + 91);
    const auto *dk_92 = buffer.data(dk + 92);
    const auto *dk_93 = buffer.data(dk + 93);
    const auto *dk_94 = buffer.data(dk + 94);
    const auto *dk_95 = buffer.data(dk + 95);
    const auto *dk_96 = buffer.data(dk + 96);
    const auto *dk_97 = buffer.data(dk + 97);
    const auto *dk_98 = buffer.data(dk + 98);
    const auto *dk_99 = buffer.data(dk + 99);
    const auto *dk_100 = buffer.data(dk + 100);
    const auto *dk_101 = buffer.data(dk + 101);
    const auto *dk_102 = buffer.data(dk + 102);
    const auto *dk_103 = buffer.data(dk + 103);
    const auto *dk_104 = buffer.data(dk + 104);
    const auto *dk_105 = buffer.data(dk + 105);
    const auto *dk_106 = buffer.data(dk + 106);
    const auto *dk_107 = buffer.data(dk + 107);
    const auto *dk_144 = buffer.data(dk + 144);
    const auto *dk_145 = buffer.data(dk + 145);
    const auto *dk_146 = buffer.data(dk + 146);
    const auto *dk_147 = buffer.data(dk + 147);
    const auto *dk_148 = buffer.data(dk + 148);
    const auto *dk_149 = buffer.data(dk + 149);
    const auto *dk_150 = buffer.data(dk + 150);
    const auto *dk_151 = buffer.data(dk + 151);
    const auto *dk_152 = buffer.data(dk + 152);
    const auto *dk_153 = buffer.data(dk + 153);
    const auto *dk_154 = buffer.data(dk + 154);
    const auto *dk_155 = buffer.data(dk + 155);
    const auto *dk_156 = buffer.data(dk + 156);
    const auto *dk_157 = buffer.data(dk + 157);
    const auto *dk_158 = buffer.data(dk + 158);
    const auto *dk_159 = buffer.data(dk + 159);
    const auto *dk_160 = buffer.data(dk + 160);
    const auto *dk_161 = buffer.data(dk + 161);
    const auto *dk_162 = buffer.data(dk + 162);
    const auto *dk_163 = buffer.data(dk + 163);
    const auto *dk_164 = buffer.data(dk + 164);
    const auto *dk_165 = buffer.data(dk + 165);
    const auto *dk_166 = buffer.data(dk + 166);
    const auto *dk_167 = buffer.data(dk + 167);
    const auto *dk_168 = buffer.data(dk + 168);
    const auto *dk_169 = buffer.data(dk + 169);
    const auto *dk_170 = buffer.data(dk + 170);
    const auto *dk_171 = buffer.data(dk + 171);
    const auto *dk_172 = buffer.data(dk + 172);
    const auto *dk_173 = buffer.data(dk + 173);
    const auto *dk_174 = buffer.data(dk + 174);
    const auto *dk_175 = buffer.data(dk + 175);
    const auto *dk_176 = buffer.data(dk + 176);
    const auto *dk_177 = buffer.data(dk + 177);
    const auto *dk_178 = buffer.data(dk + 178);
    const auto *dk_179 = buffer.data(dk + 179);
    const auto *dk_180 = buffer.data(dk + 180);
    const auto *dk_181 = buffer.data(dk + 181);
    const auto *dk_182 = buffer.data(dk + 182);
    const auto *dk_183 = buffer.data(dk + 183);
    const auto *dk_184 = buffer.data(dk + 184);
    const auto *dk_185 = buffer.data(dk + 185);
    const auto *dk_186 = buffer.data(dk + 186);
    const auto *dk_187 = buffer.data(dk + 187);
    const auto *dk_188 = buffer.data(dk + 188);
    const auto *dk_189 = buffer.data(dk + 189);
    const auto *dk_190 = buffer.data(dk + 190);
    const auto *dk_191 = buffer.data(dk + 191);
    const auto *dk_192 = buffer.data(dk + 192);
    const auto *dk_193 = buffer.data(dk + 193);
    const auto *dk_194 = buffer.data(dk + 194);
    const auto *dk_195 = buffer.data(dk + 195);
    const auto *dk_196 = buffer.data(dk + 196);
    const auto *dk_197 = buffer.data(dk + 197);
    const auto *dk_198 = buffer.data(dk + 198);
    const auto *dk_199 = buffer.data(dk + 199);
    const auto *dk_200 = buffer.data(dk + 200);
    const auto *dk_201 = buffer.data(dk + 201);
    const auto *dk_202 = buffer.data(dk + 202);
    const auto *dk_203 = buffer.data(dk + 203);
    const auto *dk_204 = buffer.data(dk + 204);
    const auto *dk_205 = buffer.data(dk + 205);
    const auto *dk_206 = buffer.data(dk + 206);
    const auto *dk_207 = buffer.data(dk + 207);
    const auto *dk_208 = buffer.data(dk + 208);
    const auto *dk_209 = buffer.data(dk + 209);
    const auto *dk_210 = buffer.data(dk + 210);
    const auto *dk_211 = buffer.data(dk + 211);
    const auto *dk_212 = buffer.data(dk + 212);
    const auto *dk_213 = buffer.data(dk + 213);
    const auto *dk_214 = buffer.data(dk + 214);
    const auto *dk_215 = buffer.data(dk + 215);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, dk_72, dk_73, dk_74, dk_75, \
                         dk_76, dk_77, dk_78, dk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dk_72[k];

        t_1[k] = f_0 * dk_73[k];

        t_2[k] = f_0 * dk_74[k];

        t_3[k] = f_0 * dk_75[k];

        t_4[k] = f_0 * dk_76[k];

        t_5[k] = f_0 * dk_77[k];

        t_6[k] = f_0 * dk_78[k];

        t_7[k] = f_0 * dk_79[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, dk_80, dk_81, dk_82, \
                         dk_83, dk_84, dk_85, dk_86, dk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dk_80[k];

        t_9[k] = f_0 * dk_81[k];

        t_10[k] = f_0 * dk_82[k];

        t_11[k] = f_0 * dk_83[k];

        t_12[k] = f_0 * dk_84[k];

        t_13[k] = f_0 * dk_85[k];

        t_14[k] = f_0 * dk_86[k];

        t_15[k] = f_0 * dk_87[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, dk_88, dk_89, dk_90, \
                         dk_91, dk_92, dk_93, dk_94, dk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * dk_88[k];

        t_17[k] = f_0 * dk_89[k];

        t_18[k] = f_0 * dk_90[k];

        t_19[k] = f_0 * dk_91[k];

        t_20[k] = f_0 * dk_92[k];

        t_21[k] = f_0 * dk_93[k];

        t_22[k] = f_0 * dk_94[k];

        t_23[k] = f_0 * dk_95[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, dk_96, dk_97, dk_98, \
                         dk_99, dk_100, dk_101, dk_102, dk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * dk_96[k];

        t_25[k] = f_0 * dk_97[k];

        t_26[k] = f_0 * dk_98[k];

        t_27[k] = f_0 * dk_99[k];

        t_28[k] = f_0 * dk_100[k];

        t_29[k] = f_0 * dk_101[k];

        t_30[k] = f_0 * dk_102[k];

        t_31[k] = f_0 * dk_103[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, dk_104, dk_105, \
                         dk_106, dk_107, dk_144, dk_145, dk_146, \
                         dk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * dk_104[k];

        t_33[k] = f_0 * dk_105[k];

        t_34[k] = f_0 * dk_106[k];

        t_35[k] = f_0 * dk_107[k];

        t_36[k] = f_0 * dk_144[k];

        t_37[k] = f_0 * dk_145[k];

        t_38[k] = f_0 * dk_146[k];

        t_39[k] = f_0 * dk_147[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, dk_148, dk_149, \
                         dk_150, dk_151, dk_152, dk_153, dk_154, \
                         dk_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * dk_148[k];

        t_41[k] = f_0 * dk_149[k];

        t_42[k] = f_0 * dk_150[k];

        t_43[k] = f_0 * dk_151[k];

        t_44[k] = f_0 * dk_152[k];

        t_45[k] = f_0 * dk_153[k];

        t_46[k] = f_0 * dk_154[k];

        t_47[k] = f_0 * dk_155[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, dk_156, dk_157, \
                         dk_158, dk_159, dk_160, dk_161, dk_162, \
                         dk_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * dk_156[k];

        t_49[k] = f_0 * dk_157[k];

        t_50[k] = f_0 * dk_158[k];

        t_51[k] = f_0 * dk_159[k];

        t_52[k] = f_0 * dk_160[k];

        t_53[k] = f_0 * dk_161[k];

        t_54[k] = f_0 * dk_162[k];

        t_55[k] = f_0 * dk_163[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, t_63, dk_164, dk_165, \
                         dk_166, dk_167, dk_168, dk_169, dk_170, \
                         dk_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * dk_164[k];

        t_57[k] = f_0 * dk_165[k];

        t_58[k] = f_0 * dk_166[k];

        t_59[k] = f_0 * dk_167[k];

        t_60[k] = f_0 * dk_168[k];

        t_61[k] = f_0 * dk_169[k];

        t_62[k] = f_0 * dk_170[k];

        t_63[k] = f_0 * dk_171[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, t_70, t_71, dk_172, dk_173, \
                         dk_174, dk_175, dk_176, dk_177, dk_178, \
                         dk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * dk_172[k];

        t_65[k] = f_0 * dk_173[k];

        t_66[k] = f_0 * dk_174[k];

        t_67[k] = f_0 * dk_175[k];

        t_68[k] = f_0 * dk_176[k];

        t_69[k] = f_0 * dk_177[k];

        t_70[k] = f_0 * dk_178[k];

        t_71[k] = f_0 * dk_179[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, sk_0, sk_1, sk_2, sk_3, sk_4, dk_180, \
                         dk_181, dk_182, dk_183, dk_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -sk_0[k]
                  + f_0 * dk_180[k];

        t_73[k] = -sk_1[k]
                  + f_0 * dk_181[k];

        t_74[k] = -sk_2[k]
                  + f_0 * dk_182[k];

        t_75[k] = -sk_3[k]
                  + f_0 * dk_183[k];

        t_76[k] = -sk_4[k]
                  + f_0 * dk_184[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, sk_5, sk_6, sk_7, sk_8, sk_9, dk_185, \
                         dk_186, dk_187, dk_188, dk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -sk_5[k]
                  + f_0 * dk_185[k];

        t_78[k] = -sk_6[k]
                  + f_0 * dk_186[k];

        t_79[k] = -sk_7[k]
                  + f_0 * dk_187[k];

        t_80[k] = -sk_8[k]
                  + f_0 * dk_188[k];

        t_81[k] = -sk_9[k]
                  + f_0 * dk_189[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, sk_10, sk_11, sk_12, sk_13, sk_14, \
                         dk_190, dk_191, dk_192, dk_193, dk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -sk_10[k]
                  + f_0 * dk_190[k];

        t_83[k] = -sk_11[k]
                  + f_0 * dk_191[k];

        t_84[k] = -sk_12[k]
                  + f_0 * dk_192[k];

        t_85[k] = -sk_13[k]
                  + f_0 * dk_193[k];

        t_86[k] = -sk_14[k]
                  + f_0 * dk_194[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, sk_15, sk_16, sk_17, sk_18, sk_19, \
                         dk_195, dk_196, dk_197, dk_198, dk_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -sk_15[k]
                  + f_0 * dk_195[k];

        t_88[k] = -sk_16[k]
                  + f_0 * dk_196[k];

        t_89[k] = -sk_17[k]
                  + f_0 * dk_197[k];

        t_90[k] = -sk_18[k]
                  + f_0 * dk_198[k];

        t_91[k] = -sk_19[k]
                  + f_0 * dk_199[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, sk_20, sk_21, sk_22, sk_23, sk_24, \
                         dk_200, dk_201, dk_202, dk_203, dk_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = -sk_20[k]
                  + f_0 * dk_200[k];

        t_93[k] = -sk_21[k]
                  + f_0 * dk_201[k];

        t_94[k] = -sk_22[k]
                  + f_0 * dk_202[k];

        t_95[k] = -sk_23[k]
                  + f_0 * dk_203[k];

        t_96[k] = -sk_24[k]
                  + f_0 * dk_204[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, sk_25, sk_26, sk_27, sk_28, sk_29, \
                         dk_205, dk_206, dk_207, dk_208, dk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -sk_25[k]
                  + f_0 * dk_205[k];

        t_98[k] = -sk_26[k]
                  + f_0 * dk_206[k];

        t_99[k] = -sk_27[k]
                  + f_0 * dk_207[k];

        t_100[k] = -sk_28[k]
                   + f_0 * dk_208[k];

        t_101[k] = -sk_29[k]
                   + f_0 * dk_209[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, sk_30, sk_31, sk_32, sk_33, sk_34, \
                         dk_210, dk_211, dk_212, dk_213, dk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -sk_30[k]
                   + f_0 * dk_210[k];

        t_103[k] = -sk_31[k]
                   + f_0 * dk_211[k];

        t_104[k] = -sk_32[k]
                   + f_0 * dk_212[k];

        t_105[k] = -sk_33[k]
                   + f_0 * dk_213[k];

        t_106[k] = -sk_34[k]
                   + f_0 * dk_214[k];
    }

#pragma omp simd aligned(t_107, sk_35, dk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -sk_35[k]
                   + f_0 * dk_215[k];
    }
}

}  // namespace simdt2ceri
