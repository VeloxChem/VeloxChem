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


#include "SimdElectronRepulsionGeom10VrrRecSI.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_si_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t pi, const size_t ncols,
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, pi_0, pi_1, pi_2, pi_3, pi_4, \
                         pi_5, pi_6, pi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k];

        t_1[k] = f_0 * pi_1[k];

        t_2[k] = f_0 * pi_2[k];

        t_3[k] = f_0 * pi_3[k];

        t_4[k] = f_0 * pi_4[k];

        t_5[k] = f_0 * pi_5[k];

        t_6[k] = f_0 * pi_6[k];

        t_7[k] = f_0 * pi_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, pi_8, pi_9, pi_10, \
                         pi_11, pi_12, pi_13, pi_14, pi_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * pi_8[k];

        t_9[k] = f_0 * pi_9[k];

        t_10[k] = f_0 * pi_10[k];

        t_11[k] = f_0 * pi_11[k];

        t_12[k] = f_0 * pi_12[k];

        t_13[k] = f_0 * pi_13[k];

        t_14[k] = f_0 * pi_14[k];

        t_15[k] = f_0 * pi_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, pi_16, pi_17, pi_18, \
                         pi_19, pi_20, pi_21, pi_22, pi_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * pi_16[k];

        t_17[k] = f_0 * pi_17[k];

        t_18[k] = f_0 * pi_18[k];

        t_19[k] = f_0 * pi_19[k];

        t_20[k] = f_0 * pi_20[k];

        t_21[k] = f_0 * pi_21[k];

        t_22[k] = f_0 * pi_22[k];

        t_23[k] = f_0 * pi_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pi_24, pi_25, pi_26, \
                         pi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * pi_24[k];

        t_25[k] = f_0 * pi_25[k];

        t_26[k] = f_0 * pi_26[k];

        t_27[k] = f_0 * pi_27[k];
    }
}

auto
compute_prim_geom_10_si_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t pi, const size_t ncols,
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
    const auto *pi_45 = buffer.data(pi + 45);
    const auto *pi_46 = buffer.data(pi + 46);
    const auto *pi_47 = buffer.data(pi + 47);
    const auto *pi_48 = buffer.data(pi + 48);
    const auto *pi_49 = buffer.data(pi + 49);
    const auto *pi_50 = buffer.data(pi + 50);
    const auto *pi_51 = buffer.data(pi + 51);
    const auto *pi_52 = buffer.data(pi + 52);
    const auto *pi_53 = buffer.data(pi + 53);
    const auto *pi_54 = buffer.data(pi + 54);
    const auto *pi_55 = buffer.data(pi + 55);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, pi_28, pi_29, pi_30, pi_31, \
                         pi_32, pi_33, pi_34, pi_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_28[k];

        t_1[k] = f_0 * pi_29[k];

        t_2[k] = f_0 * pi_30[k];

        t_3[k] = f_0 * pi_31[k];

        t_4[k] = f_0 * pi_32[k];

        t_5[k] = f_0 * pi_33[k];

        t_6[k] = f_0 * pi_34[k];

        t_7[k] = f_0 * pi_35[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, pi_36, pi_37, pi_38, \
                         pi_39, pi_40, pi_41, pi_42, pi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * pi_36[k];

        t_9[k] = f_0 * pi_37[k];

        t_10[k] = f_0 * pi_38[k];

        t_11[k] = f_0 * pi_39[k];

        t_12[k] = f_0 * pi_40[k];

        t_13[k] = f_0 * pi_41[k];

        t_14[k] = f_0 * pi_42[k];

        t_15[k] = f_0 * pi_43[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, pi_44, pi_45, pi_46, \
                         pi_47, pi_48, pi_49, pi_50, pi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * pi_44[k];

        t_17[k] = f_0 * pi_45[k];

        t_18[k] = f_0 * pi_46[k];

        t_19[k] = f_0 * pi_47[k];

        t_20[k] = f_0 * pi_48[k];

        t_21[k] = f_0 * pi_49[k];

        t_22[k] = f_0 * pi_50[k];

        t_23[k] = f_0 * pi_51[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pi_52, pi_53, pi_54, \
                         pi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * pi_52[k];

        t_25[k] = f_0 * pi_53[k];

        t_26[k] = f_0 * pi_54[k];

        t_27[k] = f_0 * pi_55[k];
    }
}

auto
compute_prim_geom_10_si_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t pi, const size_t ncols,
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

    const auto *pi_56 = buffer.data(pi + 56);
    const auto *pi_57 = buffer.data(pi + 57);
    const auto *pi_58 = buffer.data(pi + 58);
    const auto *pi_59 = buffer.data(pi + 59);
    const auto *pi_60 = buffer.data(pi + 60);
    const auto *pi_61 = buffer.data(pi + 61);
    const auto *pi_62 = buffer.data(pi + 62);
    const auto *pi_63 = buffer.data(pi + 63);
    const auto *pi_64 = buffer.data(pi + 64);
    const auto *pi_65 = buffer.data(pi + 65);
    const auto *pi_66 = buffer.data(pi + 66);
    const auto *pi_67 = buffer.data(pi + 67);
    const auto *pi_68 = buffer.data(pi + 68);
    const auto *pi_69 = buffer.data(pi + 69);
    const auto *pi_70 = buffer.data(pi + 70);
    const auto *pi_71 = buffer.data(pi + 71);
    const auto *pi_72 = buffer.data(pi + 72);
    const auto *pi_73 = buffer.data(pi + 73);
    const auto *pi_74 = buffer.data(pi + 74);
    const auto *pi_75 = buffer.data(pi + 75);
    const auto *pi_76 = buffer.data(pi + 76);
    const auto *pi_77 = buffer.data(pi + 77);
    const auto *pi_78 = buffer.data(pi + 78);
    const auto *pi_79 = buffer.data(pi + 79);
    const auto *pi_80 = buffer.data(pi + 80);
    const auto *pi_81 = buffer.data(pi + 81);
    const auto *pi_82 = buffer.data(pi + 82);
    const auto *pi_83 = buffer.data(pi + 83);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, pi_56, pi_57, pi_58, pi_59, \
                         pi_60, pi_61, pi_62, pi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_56[k];

        t_1[k] = f_0 * pi_57[k];

        t_2[k] = f_0 * pi_58[k];

        t_3[k] = f_0 * pi_59[k];

        t_4[k] = f_0 * pi_60[k];

        t_5[k] = f_0 * pi_61[k];

        t_6[k] = f_0 * pi_62[k];

        t_7[k] = f_0 * pi_63[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, pi_64, pi_65, pi_66, \
                         pi_67, pi_68, pi_69, pi_70, pi_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * pi_64[k];

        t_9[k] = f_0 * pi_65[k];

        t_10[k] = f_0 * pi_66[k];

        t_11[k] = f_0 * pi_67[k];

        t_12[k] = f_0 * pi_68[k];

        t_13[k] = f_0 * pi_69[k];

        t_14[k] = f_0 * pi_70[k];

        t_15[k] = f_0 * pi_71[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, pi_72, pi_73, pi_74, \
                         pi_75, pi_76, pi_77, pi_78, pi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * pi_72[k];

        t_17[k] = f_0 * pi_73[k];

        t_18[k] = f_0 * pi_74[k];

        t_19[k] = f_0 * pi_75[k];

        t_20[k] = f_0 * pi_76[k];

        t_21[k] = f_0 * pi_77[k];

        t_22[k] = f_0 * pi_78[k];

        t_23[k] = f_0 * pi_79[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pi_80, pi_81, pi_82, \
                         pi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * pi_80[k];

        t_25[k] = f_0 * pi_81[k];

        t_26[k] = f_0 * pi_82[k];

        t_27[k] = f_0 * pi_83[k];
    }
}

}  // namespace simdt2ceri
