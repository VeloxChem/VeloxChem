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


#include "SimdElectronRepulsionGeom10VrrRecSK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_sk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t pk, const size_t ncols,
                                             const double alpha) -> void
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

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_4 = buffer.data(pk + 4);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_7 = buffer.data(pk + 7);
    const auto *pk_8 = buffer.data(pk + 8);
    const auto *pk_9 = buffer.data(pk + 9);
    const auto *pk_10 = buffer.data(pk + 10);
    const auto *pk_11 = buffer.data(pk + 11);
    const auto *pk_12 = buffer.data(pk + 12);
    const auto *pk_13 = buffer.data(pk + 13);
    const auto *pk_14 = buffer.data(pk + 14);
    const auto *pk_15 = buffer.data(pk + 15);
    const auto *pk_16 = buffer.data(pk + 16);
    const auto *pk_17 = buffer.data(pk + 17);
    const auto *pk_18 = buffer.data(pk + 18);
    const auto *pk_19 = buffer.data(pk + 19);
    const auto *pk_20 = buffer.data(pk + 20);
    const auto *pk_21 = buffer.data(pk + 21);
    const auto *pk_22 = buffer.data(pk + 22);
    const auto *pk_23 = buffer.data(pk + 23);
    const auto *pk_24 = buffer.data(pk + 24);
    const auto *pk_25 = buffer.data(pk + 25);
    const auto *pk_26 = buffer.data(pk + 26);
    const auto *pk_27 = buffer.data(pk + 27);
    const auto *pk_28 = buffer.data(pk + 28);
    const auto *pk_29 = buffer.data(pk + 29);
    const auto *pk_30 = buffer.data(pk + 30);
    const auto *pk_31 = buffer.data(pk + 31);
    const auto *pk_32 = buffer.data(pk + 32);
    const auto *pk_33 = buffer.data(pk + 33);
    const auto *pk_34 = buffer.data(pk + 34);
    const auto *pk_35 = buffer.data(pk + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, pk_0, pk_1, pk_2, pk_3, pk_4, \
                         pk_5, pk_6, pk_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_0[k];

        t_1[k] = f_0 * pk_1[k];

        t_2[k] = f_0 * pk_2[k];

        t_3[k] = f_0 * pk_3[k];

        t_4[k] = f_0 * pk_4[k];

        t_5[k] = f_0 * pk_5[k];

        t_6[k] = f_0 * pk_6[k];

        t_7[k] = f_0 * pk_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, pk_8, pk_9, pk_10, \
                         pk_11, pk_12, pk_13, pk_14, pk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * pk_8[k];

        t_9[k] = f_0 * pk_9[k];

        t_10[k] = f_0 * pk_10[k];

        t_11[k] = f_0 * pk_11[k];

        t_12[k] = f_0 * pk_12[k];

        t_13[k] = f_0 * pk_13[k];

        t_14[k] = f_0 * pk_14[k];

        t_15[k] = f_0 * pk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, pk_16, pk_17, pk_18, \
                         pk_19, pk_20, pk_21, pk_22, pk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * pk_16[k];

        t_17[k] = f_0 * pk_17[k];

        t_18[k] = f_0 * pk_18[k];

        t_19[k] = f_0 * pk_19[k];

        t_20[k] = f_0 * pk_20[k];

        t_21[k] = f_0 * pk_21[k];

        t_22[k] = f_0 * pk_22[k];

        t_23[k] = f_0 * pk_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, pk_24, pk_25, pk_26, \
                         pk_27, pk_28, pk_29, pk_30, pk_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * pk_24[k];

        t_25[k] = f_0 * pk_25[k];

        t_26[k] = f_0 * pk_26[k];

        t_27[k] = f_0 * pk_27[k];

        t_28[k] = f_0 * pk_28[k];

        t_29[k] = f_0 * pk_29[k];

        t_30[k] = f_0 * pk_30[k];

        t_31[k] = f_0 * pk_31[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pk_32, pk_33, pk_34, \
                         pk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * pk_32[k];

        t_33[k] = f_0 * pk_33[k];

        t_34[k] = f_0 * pk_34[k];

        t_35[k] = f_0 * pk_35[k];
    }
}

auto
compute_prim_geom_10_sk_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t pk, const size_t ncols,
                                             const double alpha) -> void
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

    const auto *pk_36 = buffer.data(pk + 36);
    const auto *pk_37 = buffer.data(pk + 37);
    const auto *pk_38 = buffer.data(pk + 38);
    const auto *pk_39 = buffer.data(pk + 39);
    const auto *pk_40 = buffer.data(pk + 40);
    const auto *pk_41 = buffer.data(pk + 41);
    const auto *pk_42 = buffer.data(pk + 42);
    const auto *pk_43 = buffer.data(pk + 43);
    const auto *pk_44 = buffer.data(pk + 44);
    const auto *pk_45 = buffer.data(pk + 45);
    const auto *pk_46 = buffer.data(pk + 46);
    const auto *pk_47 = buffer.data(pk + 47);
    const auto *pk_48 = buffer.data(pk + 48);
    const auto *pk_49 = buffer.data(pk + 49);
    const auto *pk_50 = buffer.data(pk + 50);
    const auto *pk_51 = buffer.data(pk + 51);
    const auto *pk_52 = buffer.data(pk + 52);
    const auto *pk_53 = buffer.data(pk + 53);
    const auto *pk_54 = buffer.data(pk + 54);
    const auto *pk_55 = buffer.data(pk + 55);
    const auto *pk_56 = buffer.data(pk + 56);
    const auto *pk_57 = buffer.data(pk + 57);
    const auto *pk_58 = buffer.data(pk + 58);
    const auto *pk_59 = buffer.data(pk + 59);
    const auto *pk_60 = buffer.data(pk + 60);
    const auto *pk_61 = buffer.data(pk + 61);
    const auto *pk_62 = buffer.data(pk + 62);
    const auto *pk_63 = buffer.data(pk + 63);
    const auto *pk_64 = buffer.data(pk + 64);
    const auto *pk_65 = buffer.data(pk + 65);
    const auto *pk_66 = buffer.data(pk + 66);
    const auto *pk_67 = buffer.data(pk + 67);
    const auto *pk_68 = buffer.data(pk + 68);
    const auto *pk_69 = buffer.data(pk + 69);
    const auto *pk_70 = buffer.data(pk + 70);
    const auto *pk_71 = buffer.data(pk + 71);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, pk_36, pk_37, pk_38, pk_39, \
                         pk_40, pk_41, pk_42, pk_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_36[k];

        t_1[k] = f_0 * pk_37[k];

        t_2[k] = f_0 * pk_38[k];

        t_3[k] = f_0 * pk_39[k];

        t_4[k] = f_0 * pk_40[k];

        t_5[k] = f_0 * pk_41[k];

        t_6[k] = f_0 * pk_42[k];

        t_7[k] = f_0 * pk_43[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, pk_44, pk_45, pk_46, \
                         pk_47, pk_48, pk_49, pk_50, pk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * pk_44[k];

        t_9[k] = f_0 * pk_45[k];

        t_10[k] = f_0 * pk_46[k];

        t_11[k] = f_0 * pk_47[k];

        t_12[k] = f_0 * pk_48[k];

        t_13[k] = f_0 * pk_49[k];

        t_14[k] = f_0 * pk_50[k];

        t_15[k] = f_0 * pk_51[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, pk_52, pk_53, pk_54, \
                         pk_55, pk_56, pk_57, pk_58, pk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * pk_52[k];

        t_17[k] = f_0 * pk_53[k];

        t_18[k] = f_0 * pk_54[k];

        t_19[k] = f_0 * pk_55[k];

        t_20[k] = f_0 * pk_56[k];

        t_21[k] = f_0 * pk_57[k];

        t_22[k] = f_0 * pk_58[k];

        t_23[k] = f_0 * pk_59[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, pk_60, pk_61, pk_62, \
                         pk_63, pk_64, pk_65, pk_66, pk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * pk_60[k];

        t_25[k] = f_0 * pk_61[k];

        t_26[k] = f_0 * pk_62[k];

        t_27[k] = f_0 * pk_63[k];

        t_28[k] = f_0 * pk_64[k];

        t_29[k] = f_0 * pk_65[k];

        t_30[k] = f_0 * pk_66[k];

        t_31[k] = f_0 * pk_67[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pk_68, pk_69, pk_70, \
                         pk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * pk_68[k];

        t_33[k] = f_0 * pk_69[k];

        t_34[k] = f_0 * pk_70[k];

        t_35[k] = f_0 * pk_71[k];
    }
}

auto
compute_prim_geom_10_sk_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t pk, const size_t ncols,
                                             const double alpha) -> void
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

    const auto *pk_72 = buffer.data(pk + 72);
    const auto *pk_73 = buffer.data(pk + 73);
    const auto *pk_74 = buffer.data(pk + 74);
    const auto *pk_75 = buffer.data(pk + 75);
    const auto *pk_76 = buffer.data(pk + 76);
    const auto *pk_77 = buffer.data(pk + 77);
    const auto *pk_78 = buffer.data(pk + 78);
    const auto *pk_79 = buffer.data(pk + 79);
    const auto *pk_80 = buffer.data(pk + 80);
    const auto *pk_81 = buffer.data(pk + 81);
    const auto *pk_82 = buffer.data(pk + 82);
    const auto *pk_83 = buffer.data(pk + 83);
    const auto *pk_84 = buffer.data(pk + 84);
    const auto *pk_85 = buffer.data(pk + 85);
    const auto *pk_86 = buffer.data(pk + 86);
    const auto *pk_87 = buffer.data(pk + 87);
    const auto *pk_88 = buffer.data(pk + 88);
    const auto *pk_89 = buffer.data(pk + 89);
    const auto *pk_90 = buffer.data(pk + 90);
    const auto *pk_91 = buffer.data(pk + 91);
    const auto *pk_92 = buffer.data(pk + 92);
    const auto *pk_93 = buffer.data(pk + 93);
    const auto *pk_94 = buffer.data(pk + 94);
    const auto *pk_95 = buffer.data(pk + 95);
    const auto *pk_96 = buffer.data(pk + 96);
    const auto *pk_97 = buffer.data(pk + 97);
    const auto *pk_98 = buffer.data(pk + 98);
    const auto *pk_99 = buffer.data(pk + 99);
    const auto *pk_100 = buffer.data(pk + 100);
    const auto *pk_101 = buffer.data(pk + 101);
    const auto *pk_102 = buffer.data(pk + 102);
    const auto *pk_103 = buffer.data(pk + 103);
    const auto *pk_104 = buffer.data(pk + 104);
    const auto *pk_105 = buffer.data(pk + 105);
    const auto *pk_106 = buffer.data(pk + 106);
    const auto *pk_107 = buffer.data(pk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, pk_72, pk_73, pk_74, pk_75, \
                         pk_76, pk_77, pk_78, pk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pk_72[k];

        t_1[k] = f_0 * pk_73[k];

        t_2[k] = f_0 * pk_74[k];

        t_3[k] = f_0 * pk_75[k];

        t_4[k] = f_0 * pk_76[k];

        t_5[k] = f_0 * pk_77[k];

        t_6[k] = f_0 * pk_78[k];

        t_7[k] = f_0 * pk_79[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, pk_80, pk_81, pk_82, \
                         pk_83, pk_84, pk_85, pk_86, pk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * pk_80[k];

        t_9[k] = f_0 * pk_81[k];

        t_10[k] = f_0 * pk_82[k];

        t_11[k] = f_0 * pk_83[k];

        t_12[k] = f_0 * pk_84[k];

        t_13[k] = f_0 * pk_85[k];

        t_14[k] = f_0 * pk_86[k];

        t_15[k] = f_0 * pk_87[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, pk_88, pk_89, pk_90, \
                         pk_91, pk_92, pk_93, pk_94, pk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * pk_88[k];

        t_17[k] = f_0 * pk_89[k];

        t_18[k] = f_0 * pk_90[k];

        t_19[k] = f_0 * pk_91[k];

        t_20[k] = f_0 * pk_92[k];

        t_21[k] = f_0 * pk_93[k];

        t_22[k] = f_0 * pk_94[k];

        t_23[k] = f_0 * pk_95[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, pk_96, pk_97, pk_98, \
                         pk_99, pk_100, pk_101, pk_102, pk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * pk_96[k];

        t_25[k] = f_0 * pk_97[k];

        t_26[k] = f_0 * pk_98[k];

        t_27[k] = f_0 * pk_99[k];

        t_28[k] = f_0 * pk_100[k];

        t_29[k] = f_0 * pk_101[k];

        t_30[k] = f_0 * pk_102[k];

        t_31[k] = f_0 * pk_103[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pk_104, pk_105, pk_106, \
                         pk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * pk_104[k];

        t_33[k] = f_0 * pk_105[k];

        t_34[k] = f_0 * pk_106[k];

        t_35[k] = f_0 * pk_107[k];
    }
}

}  // namespace simdt2ceri
