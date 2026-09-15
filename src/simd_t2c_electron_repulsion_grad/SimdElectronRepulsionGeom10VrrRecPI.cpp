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


#include "SimdElectronRepulsionGeom10VrrRecPI.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_pi_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t si, const size_t di,
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
    const auto *si_25 = buffer.data(si + 25);
    const auto *si_26 = buffer.data(si + 26);
    const auto *si_27 = buffer.data(si + 27);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_26 = buffer.data(di + 26);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_50 = buffer.data(di + 50);
    const auto *di_51 = buffer.data(di + 51);
    const auto *di_52 = buffer.data(di + 52);
    const auto *di_53 = buffer.data(di + 53);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_55 = buffer.data(di + 55);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_57 = buffer.data(di + 57);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_63 = buffer.data(di + 63);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_66 = buffer.data(di + 66);
    const auto *di_67 = buffer.data(di + 67);
    const auto *di_68 = buffer.data(di + 68);
    const auto *di_69 = buffer.data(di + 69);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_71 = buffer.data(di + 71);
    const auto *di_72 = buffer.data(di + 72);
    const auto *di_73 = buffer.data(di + 73);
    const auto *di_74 = buffer.data(di + 74);
    const auto *di_75 = buffer.data(di + 75);
    const auto *di_76 = buffer.data(di + 76);
    const auto *di_77 = buffer.data(di + 77);
    const auto *di_78 = buffer.data(di + 78);
    const auto *di_79 = buffer.data(di + 79);
    const auto *di_80 = buffer.data(di + 80);
    const auto *di_81 = buffer.data(di + 81);
    const auto *di_82 = buffer.data(di + 82);
    const auto *di_83 = buffer.data(di + 83);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, si_0, si_1, si_2, si_3, si_4, di_0, di_1, \
                         di_2, di_3, di_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -si_0[k]
                 + f_0 * di_0[k];

        t_1[k] = -si_1[k]
                 + f_0 * di_1[k];

        t_2[k] = -si_2[k]
                 + f_0 * di_2[k];

        t_3[k] = -si_3[k]
                 + f_0 * di_3[k];

        t_4[k] = -si_4[k]
                 + f_0 * di_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, si_5, si_6, si_7, si_8, si_9, di_5, di_6, \
                         di_7, di_8, di_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -si_5[k]
                 + f_0 * di_5[k];

        t_6[k] = -si_6[k]
                 + f_0 * di_6[k];

        t_7[k] = -si_7[k]
                 + f_0 * di_7[k];

        t_8[k] = -si_8[k]
                 + f_0 * di_8[k];

        t_9[k] = -si_9[k]
                 + f_0 * di_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, si_10, si_11, si_12, si_13, si_14, \
                         di_10, di_11, di_12, di_13, di_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -si_10[k]
                  + f_0 * di_10[k];

        t_11[k] = -si_11[k]
                  + f_0 * di_11[k];

        t_12[k] = -si_12[k]
                  + f_0 * di_12[k];

        t_13[k] = -si_13[k]
                  + f_0 * di_13[k];

        t_14[k] = -si_14[k]
                  + f_0 * di_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, si_15, si_16, si_17, si_18, si_19, \
                         di_15, di_16, di_17, di_18, di_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -si_15[k]
                  + f_0 * di_15[k];

        t_16[k] = -si_16[k]
                  + f_0 * di_16[k];

        t_17[k] = -si_17[k]
                  + f_0 * di_17[k];

        t_18[k] = -si_18[k]
                  + f_0 * di_18[k];

        t_19[k] = -si_19[k]
                  + f_0 * di_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, si_20, si_21, si_22, si_23, si_24, \
                         di_20, di_21, di_22, di_23, di_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -si_20[k]
                  + f_0 * di_20[k];

        t_21[k] = -si_21[k]
                  + f_0 * di_21[k];

        t_22[k] = -si_22[k]
                  + f_0 * di_22[k];

        t_23[k] = -si_23[k]
                  + f_0 * di_23[k];

        t_24[k] = -si_24[k]
                  + f_0 * di_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, si_25, si_26, si_27, di_25, \
                         di_26, di_27, di_28, di_29, di_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -si_25[k]
                  + f_0 * di_25[k];

        t_26[k] = -si_26[k]
                  + f_0 * di_26[k];

        t_27[k] = -si_27[k]
                  + f_0 * di_27[k];

        t_28[k] = f_0 * di_28[k];

        t_29[k] = f_0 * di_29[k];

        t_30[k] = f_0 * di_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, t_37, t_38, di_31, di_32, di_33, \
                         di_34, di_35, di_36, di_37, di_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * di_31[k];

        t_32[k] = f_0 * di_32[k];

        t_33[k] = f_0 * di_33[k];

        t_34[k] = f_0 * di_34[k];

        t_35[k] = f_0 * di_35[k];

        t_36[k] = f_0 * di_36[k];

        t_37[k] = f_0 * di_37[k];

        t_38[k] = f_0 * di_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, t_44, t_45, t_46, di_39, di_40, di_41, \
                         di_42, di_43, di_44, di_45, di_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * di_39[k];

        t_40[k] = f_0 * di_40[k];

        t_41[k] = f_0 * di_41[k];

        t_42[k] = f_0 * di_42[k];

        t_43[k] = f_0 * di_43[k];

        t_44[k] = f_0 * di_44[k];

        t_45[k] = f_0 * di_45[k];

        t_46[k] = f_0 * di_46[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, t_53, t_54, di_47, di_48, di_49, \
                         di_50, di_51, di_52, di_53, di_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * di_47[k];

        t_48[k] = f_0 * di_48[k];

        t_49[k] = f_0 * di_49[k];

        t_50[k] = f_0 * di_50[k];

        t_51[k] = f_0 * di_51[k];

        t_52[k] = f_0 * di_52[k];

        t_53[k] = f_0 * di_53[k];

        t_54[k] = f_0 * di_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, t_62, di_55, di_56, di_57, \
                         di_58, di_59, di_60, di_61, di_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * di_55[k];

        t_56[k] = f_0 * di_56[k];

        t_57[k] = f_0 * di_57[k];

        t_58[k] = f_0 * di_58[k];

        t_59[k] = f_0 * di_59[k];

        t_60[k] = f_0 * di_60[k];

        t_61[k] = f_0 * di_61[k];

        t_62[k] = f_0 * di_62[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, t_68, t_69, t_70, di_63, di_64, di_65, \
                         di_66, di_67, di_68, di_69, di_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_0 * di_63[k];

        t_64[k] = f_0 * di_64[k];

        t_65[k] = f_0 * di_65[k];

        t_66[k] = f_0 * di_66[k];

        t_67[k] = f_0 * di_67[k];

        t_68[k] = f_0 * di_68[k];

        t_69[k] = f_0 * di_69[k];

        t_70[k] = f_0 * di_70[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, t_77, t_78, di_71, di_72, di_73, \
                         di_74, di_75, di_76, di_77, di_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_0 * di_71[k];

        t_72[k] = f_0 * di_72[k];

        t_73[k] = f_0 * di_73[k];

        t_74[k] = f_0 * di_74[k];

        t_75[k] = f_0 * di_75[k];

        t_76[k] = f_0 * di_76[k];

        t_77[k] = f_0 * di_77[k];

        t_78[k] = f_0 * di_78[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, di_79, di_80, di_81, di_82, \
                         di_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_0 * di_79[k];

        t_80[k] = f_0 * di_80[k];

        t_81[k] = f_0 * di_81[k];

        t_82[k] = f_0 * di_82[k];

        t_83[k] = f_0 * di_83[k];
    }
}

auto
compute_prim_geom_10_pi_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t si, const size_t di,
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
    const auto *si_25 = buffer.data(si + 25);
    const auto *si_26 = buffer.data(si + 26);
    const auto *si_27 = buffer.data(si + 27);

    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_50 = buffer.data(di + 50);
    const auto *di_51 = buffer.data(di + 51);
    const auto *di_52 = buffer.data(di + 52);
    const auto *di_53 = buffer.data(di + 53);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_55 = buffer.data(di + 55);
    const auto *di_84 = buffer.data(di + 84);
    const auto *di_85 = buffer.data(di + 85);
    const auto *di_86 = buffer.data(di + 86);
    const auto *di_87 = buffer.data(di + 87);
    const auto *di_88 = buffer.data(di + 88);
    const auto *di_89 = buffer.data(di + 89);
    const auto *di_90 = buffer.data(di + 90);
    const auto *di_91 = buffer.data(di + 91);
    const auto *di_92 = buffer.data(di + 92);
    const auto *di_93 = buffer.data(di + 93);
    const auto *di_94 = buffer.data(di + 94);
    const auto *di_95 = buffer.data(di + 95);
    const auto *di_96 = buffer.data(di + 96);
    const auto *di_97 = buffer.data(di + 97);
    const auto *di_98 = buffer.data(di + 98);
    const auto *di_99 = buffer.data(di + 99);
    const auto *di_100 = buffer.data(di + 100);
    const auto *di_101 = buffer.data(di + 101);
    const auto *di_102 = buffer.data(di + 102);
    const auto *di_103 = buffer.data(di + 103);
    const auto *di_104 = buffer.data(di + 104);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_106 = buffer.data(di + 106);
    const auto *di_107 = buffer.data(di + 107);
    const auto *di_108 = buffer.data(di + 108);
    const auto *di_109 = buffer.data(di + 109);
    const auto *di_110 = buffer.data(di + 110);
    const auto *di_111 = buffer.data(di + 111);
    const auto *di_112 = buffer.data(di + 112);
    const auto *di_113 = buffer.data(di + 113);
    const auto *di_114 = buffer.data(di + 114);
    const auto *di_115 = buffer.data(di + 115);
    const auto *di_116 = buffer.data(di + 116);
    const auto *di_117 = buffer.data(di + 117);
    const auto *di_118 = buffer.data(di + 118);
    const auto *di_119 = buffer.data(di + 119);
    const auto *di_120 = buffer.data(di + 120);
    const auto *di_121 = buffer.data(di + 121);
    const auto *di_122 = buffer.data(di + 122);
    const auto *di_123 = buffer.data(di + 123);
    const auto *di_124 = buffer.data(di + 124);
    const auto *di_125 = buffer.data(di + 125);
    const auto *di_126 = buffer.data(di + 126);
    const auto *di_127 = buffer.data(di + 127);
    const auto *di_128 = buffer.data(di + 128);
    const auto *di_129 = buffer.data(di + 129);
    const auto *di_130 = buffer.data(di + 130);
    const auto *di_131 = buffer.data(di + 131);
    const auto *di_132 = buffer.data(di + 132);
    const auto *di_133 = buffer.data(di + 133);
    const auto *di_134 = buffer.data(di + 134);
    const auto *di_135 = buffer.data(di + 135);
    const auto *di_136 = buffer.data(di + 136);
    const auto *di_137 = buffer.data(di + 137);
    const auto *di_138 = buffer.data(di + 138);
    const auto *di_139 = buffer.data(di + 139);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, di_28, di_29, di_30, di_31, \
                         di_32, di_33, di_34, di_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_28[k];

        t_1[k] = f_0 * di_29[k];

        t_2[k] = f_0 * di_30[k];

        t_3[k] = f_0 * di_31[k];

        t_4[k] = f_0 * di_32[k];

        t_5[k] = f_0 * di_33[k];

        t_6[k] = f_0 * di_34[k];

        t_7[k] = f_0 * di_35[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, di_36, di_37, di_38, \
                         di_39, di_40, di_41, di_42, di_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_36[k];

        t_9[k] = f_0 * di_37[k];

        t_10[k] = f_0 * di_38[k];

        t_11[k] = f_0 * di_39[k];

        t_12[k] = f_0 * di_40[k];

        t_13[k] = f_0 * di_41[k];

        t_14[k] = f_0 * di_42[k];

        t_15[k] = f_0 * di_43[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, di_44, di_45, di_46, \
                         di_47, di_48, di_49, di_50, di_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * di_44[k];

        t_17[k] = f_0 * di_45[k];

        t_18[k] = f_0 * di_46[k];

        t_19[k] = f_0 * di_47[k];

        t_20[k] = f_0 * di_48[k];

        t_21[k] = f_0 * di_49[k];

        t_22[k] = f_0 * di_50[k];

        t_23[k] = f_0 * di_51[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, si_0, si_1, di_52, di_53, di_54, \
                         di_55, di_84, di_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * di_52[k];

        t_25[k] = f_0 * di_53[k];

        t_26[k] = f_0 * di_54[k];

        t_27[k] = f_0 * di_55[k];

        t_28[k] = -si_0[k]
                  + f_0 * di_84[k];

        t_29[k] = -si_1[k]
                  + f_0 * di_85[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, si_2, si_3, si_4, si_5, si_6, di_86, \
                         di_87, di_88, di_89, di_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -si_2[k]
                  + f_0 * di_86[k];

        t_31[k] = -si_3[k]
                  + f_0 * di_87[k];

        t_32[k] = -si_4[k]
                  + f_0 * di_88[k];

        t_33[k] = -si_5[k]
                  + f_0 * di_89[k];

        t_34[k] = -si_6[k]
                  + f_0 * di_90[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, si_7, si_8, si_9, si_10, si_11, di_91, \
                         di_92, di_93, di_94, di_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -si_7[k]
                  + f_0 * di_91[k];

        t_36[k] = -si_8[k]
                  + f_0 * di_92[k];

        t_37[k] = -si_9[k]
                  + f_0 * di_93[k];

        t_38[k] = -si_10[k]
                  + f_0 * di_94[k];

        t_39[k] = -si_11[k]
                  + f_0 * di_95[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, si_12, si_13, si_14, si_15, si_16, \
                         di_96, di_97, di_98, di_99, di_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -si_12[k]
                  + f_0 * di_96[k];

        t_41[k] = -si_13[k]
                  + f_0 * di_97[k];

        t_42[k] = -si_14[k]
                  + f_0 * di_98[k];

        t_43[k] = -si_15[k]
                  + f_0 * di_99[k];

        t_44[k] = -si_16[k]
                  + f_0 * di_100[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, si_17, si_18, si_19, si_20, si_21, \
                         di_101, di_102, di_103, di_104, di_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -si_17[k]
                  + f_0 * di_101[k];

        t_46[k] = -si_18[k]
                  + f_0 * di_102[k];

        t_47[k] = -si_19[k]
                  + f_0 * di_103[k];

        t_48[k] = -si_20[k]
                  + f_0 * di_104[k];

        t_49[k] = -si_21[k]
                  + f_0 * di_105[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, si_22, si_23, si_24, si_25, si_26, \
                         di_106, di_107, di_108, di_109, di_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -si_22[k]
                  + f_0 * di_106[k];

        t_51[k] = -si_23[k]
                  + f_0 * di_107[k];

        t_52[k] = -si_24[k]
                  + f_0 * di_108[k];

        t_53[k] = -si_25[k]
                  + f_0 * di_109[k];

        t_54[k] = -si_26[k]
                  + f_0 * di_110[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, si_27, di_111, di_112, \
                         di_113, di_114, di_115, di_116, di_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -si_27[k]
                  + f_0 * di_111[k];

        t_56[k] = f_0 * di_112[k];

        t_57[k] = f_0 * di_113[k];

        t_58[k] = f_0 * di_114[k];

        t_59[k] = f_0 * di_115[k];

        t_60[k] = f_0 * di_116[k];

        t_61[k] = f_0 * di_117[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, t_67, t_68, t_69, di_118, di_119, \
                         di_120, di_121, di_122, di_123, di_124, \
                         di_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_0 * di_118[k];

        t_63[k] = f_0 * di_119[k];

        t_64[k] = f_0 * di_120[k];

        t_65[k] = f_0 * di_121[k];

        t_66[k] = f_0 * di_122[k];

        t_67[k] = f_0 * di_123[k];

        t_68[k] = f_0 * di_124[k];

        t_69[k] = f_0 * di_125[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, t_75, t_76, t_77, di_126, di_127, \
                         di_128, di_129, di_130, di_131, di_132, \
                         di_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_0 * di_126[k];

        t_71[k] = f_0 * di_127[k];

        t_72[k] = f_0 * di_128[k];

        t_73[k] = f_0 * di_129[k];

        t_74[k] = f_0 * di_130[k];

        t_75[k] = f_0 * di_131[k];

        t_76[k] = f_0 * di_132[k];

        t_77[k] = f_0 * di_133[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, di_134, di_135, di_136, di_137, \
                         di_138, di_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * di_134[k];

        t_79[k] = f_0 * di_135[k];

        t_80[k] = f_0 * di_136[k];

        t_81[k] = f_0 * di_137[k];

        t_82[k] = f_0 * di_138[k];

        t_83[k] = f_0 * di_139[k];
    }
}

auto
compute_prim_geom_10_pi_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t si, const size_t di,
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
    const auto *si_25 = buffer.data(si + 25);
    const auto *si_26 = buffer.data(si + 26);
    const auto *si_27 = buffer.data(si + 27);

    const auto *di_56 = buffer.data(di + 56);
    const auto *di_57 = buffer.data(di + 57);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_63 = buffer.data(di + 63);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_66 = buffer.data(di + 66);
    const auto *di_67 = buffer.data(di + 67);
    const auto *di_68 = buffer.data(di + 68);
    const auto *di_69 = buffer.data(di + 69);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_71 = buffer.data(di + 71);
    const auto *di_72 = buffer.data(di + 72);
    const auto *di_73 = buffer.data(di + 73);
    const auto *di_74 = buffer.data(di + 74);
    const auto *di_75 = buffer.data(di + 75);
    const auto *di_76 = buffer.data(di + 76);
    const auto *di_77 = buffer.data(di + 77);
    const auto *di_78 = buffer.data(di + 78);
    const auto *di_79 = buffer.data(di + 79);
    const auto *di_80 = buffer.data(di + 80);
    const auto *di_81 = buffer.data(di + 81);
    const auto *di_82 = buffer.data(di + 82);
    const auto *di_83 = buffer.data(di + 83);
    const auto *di_112 = buffer.data(di + 112);
    const auto *di_113 = buffer.data(di + 113);
    const auto *di_114 = buffer.data(di + 114);
    const auto *di_115 = buffer.data(di + 115);
    const auto *di_116 = buffer.data(di + 116);
    const auto *di_117 = buffer.data(di + 117);
    const auto *di_118 = buffer.data(di + 118);
    const auto *di_119 = buffer.data(di + 119);
    const auto *di_120 = buffer.data(di + 120);
    const auto *di_121 = buffer.data(di + 121);
    const auto *di_122 = buffer.data(di + 122);
    const auto *di_123 = buffer.data(di + 123);
    const auto *di_124 = buffer.data(di + 124);
    const auto *di_125 = buffer.data(di + 125);
    const auto *di_126 = buffer.data(di + 126);
    const auto *di_127 = buffer.data(di + 127);
    const auto *di_128 = buffer.data(di + 128);
    const auto *di_129 = buffer.data(di + 129);
    const auto *di_130 = buffer.data(di + 130);
    const auto *di_131 = buffer.data(di + 131);
    const auto *di_132 = buffer.data(di + 132);
    const auto *di_133 = buffer.data(di + 133);
    const auto *di_134 = buffer.data(di + 134);
    const auto *di_135 = buffer.data(di + 135);
    const auto *di_136 = buffer.data(di + 136);
    const auto *di_137 = buffer.data(di + 137);
    const auto *di_138 = buffer.data(di + 138);
    const auto *di_139 = buffer.data(di + 139);
    const auto *di_140 = buffer.data(di + 140);
    const auto *di_141 = buffer.data(di + 141);
    const auto *di_142 = buffer.data(di + 142);
    const auto *di_143 = buffer.data(di + 143);
    const auto *di_144 = buffer.data(di + 144);
    const auto *di_145 = buffer.data(di + 145);
    const auto *di_146 = buffer.data(di + 146);
    const auto *di_147 = buffer.data(di + 147);
    const auto *di_148 = buffer.data(di + 148);
    const auto *di_149 = buffer.data(di + 149);
    const auto *di_150 = buffer.data(di + 150);
    const auto *di_151 = buffer.data(di + 151);
    const auto *di_152 = buffer.data(di + 152);
    const auto *di_153 = buffer.data(di + 153);
    const auto *di_154 = buffer.data(di + 154);
    const auto *di_155 = buffer.data(di + 155);
    const auto *di_156 = buffer.data(di + 156);
    const auto *di_157 = buffer.data(di + 157);
    const auto *di_158 = buffer.data(di + 158);
    const auto *di_159 = buffer.data(di + 159);
    const auto *di_160 = buffer.data(di + 160);
    const auto *di_161 = buffer.data(di + 161);
    const auto *di_162 = buffer.data(di + 162);
    const auto *di_163 = buffer.data(di + 163);
    const auto *di_164 = buffer.data(di + 164);
    const auto *di_165 = buffer.data(di + 165);
    const auto *di_166 = buffer.data(di + 166);
    const auto *di_167 = buffer.data(di + 167);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, di_56, di_57, di_58, di_59, \
                         di_60, di_61, di_62, di_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * di_56[k];

        t_1[k] = f_0 * di_57[k];

        t_2[k] = f_0 * di_58[k];

        t_3[k] = f_0 * di_59[k];

        t_4[k] = f_0 * di_60[k];

        t_5[k] = f_0 * di_61[k];

        t_6[k] = f_0 * di_62[k];

        t_7[k] = f_0 * di_63[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, di_64, di_65, di_66, \
                         di_67, di_68, di_69, di_70, di_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * di_64[k];

        t_9[k] = f_0 * di_65[k];

        t_10[k] = f_0 * di_66[k];

        t_11[k] = f_0 * di_67[k];

        t_12[k] = f_0 * di_68[k];

        t_13[k] = f_0 * di_69[k];

        t_14[k] = f_0 * di_70[k];

        t_15[k] = f_0 * di_71[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, di_72, di_73, di_74, \
                         di_75, di_76, di_77, di_78, di_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * di_72[k];

        t_17[k] = f_0 * di_73[k];

        t_18[k] = f_0 * di_74[k];

        t_19[k] = f_0 * di_75[k];

        t_20[k] = f_0 * di_76[k];

        t_21[k] = f_0 * di_77[k];

        t_22[k] = f_0 * di_78[k];

        t_23[k] = f_0 * di_79[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, di_80, di_81, di_82, \
                         di_83, di_112, di_113, di_114, di_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * di_80[k];

        t_25[k] = f_0 * di_81[k];

        t_26[k] = f_0 * di_82[k];

        t_27[k] = f_0 * di_83[k];

        t_28[k] = f_0 * di_112[k];

        t_29[k] = f_0 * di_113[k];

        t_30[k] = f_0 * di_114[k];

        t_31[k] = f_0 * di_115[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, di_116, di_117, \
                         di_118, di_119, di_120, di_121, di_122, \
                         di_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * di_116[k];

        t_33[k] = f_0 * di_117[k];

        t_34[k] = f_0 * di_118[k];

        t_35[k] = f_0 * di_119[k];

        t_36[k] = f_0 * di_120[k];

        t_37[k] = f_0 * di_121[k];

        t_38[k] = f_0 * di_122[k];

        t_39[k] = f_0 * di_123[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, di_124, di_125, \
                         di_126, di_127, di_128, di_129, di_130, \
                         di_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * di_124[k];

        t_41[k] = f_0 * di_125[k];

        t_42[k] = f_0 * di_126[k];

        t_43[k] = f_0 * di_127[k];

        t_44[k] = f_0 * di_128[k];

        t_45[k] = f_0 * di_129[k];

        t_46[k] = f_0 * di_130[k];

        t_47[k] = f_0 * di_131[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, di_132, di_133, \
                         di_134, di_135, di_136, di_137, di_138, \
                         di_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * di_132[k];

        t_49[k] = f_0 * di_133[k];

        t_50[k] = f_0 * di_134[k];

        t_51[k] = f_0 * di_135[k];

        t_52[k] = f_0 * di_136[k];

        t_53[k] = f_0 * di_137[k];

        t_54[k] = f_0 * di_138[k];

        t_55[k] = f_0 * di_139[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, si_0, si_1, si_2, si_3, si_4, di_140, \
                         di_141, di_142, di_143, di_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -si_0[k]
                  + f_0 * di_140[k];

        t_57[k] = -si_1[k]
                  + f_0 * di_141[k];

        t_58[k] = -si_2[k]
                  + f_0 * di_142[k];

        t_59[k] = -si_3[k]
                  + f_0 * di_143[k];

        t_60[k] = -si_4[k]
                  + f_0 * di_144[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, si_5, si_6, si_7, si_8, si_9, di_145, \
                         di_146, di_147, di_148, di_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -si_5[k]
                  + f_0 * di_145[k];

        t_62[k] = -si_6[k]
                  + f_0 * di_146[k];

        t_63[k] = -si_7[k]
                  + f_0 * di_147[k];

        t_64[k] = -si_8[k]
                  + f_0 * di_148[k];

        t_65[k] = -si_9[k]
                  + f_0 * di_149[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, si_10, si_11, si_12, si_13, si_14, \
                         di_150, di_151, di_152, di_153, di_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -si_10[k]
                  + f_0 * di_150[k];

        t_67[k] = -si_11[k]
                  + f_0 * di_151[k];

        t_68[k] = -si_12[k]
                  + f_0 * di_152[k];

        t_69[k] = -si_13[k]
                  + f_0 * di_153[k];

        t_70[k] = -si_14[k]
                  + f_0 * di_154[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, si_15, si_16, si_17, si_18, si_19, \
                         di_155, di_156, di_157, di_158, di_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -si_15[k]
                  + f_0 * di_155[k];

        t_72[k] = -si_16[k]
                  + f_0 * di_156[k];

        t_73[k] = -si_17[k]
                  + f_0 * di_157[k];

        t_74[k] = -si_18[k]
                  + f_0 * di_158[k];

        t_75[k] = -si_19[k]
                  + f_0 * di_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, si_20, si_21, si_22, si_23, si_24, \
                         di_160, di_161, di_162, di_163, di_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -si_20[k]
                  + f_0 * di_160[k];

        t_77[k] = -si_21[k]
                  + f_0 * di_161[k];

        t_78[k] = -si_22[k]
                  + f_0 * di_162[k];

        t_79[k] = -si_23[k]
                  + f_0 * di_163[k];

        t_80[k] = -si_24[k]
                  + f_0 * di_164[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, si_25, si_26, si_27, di_165, di_166, \
                         di_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -si_25[k]
                  + f_0 * di_165[k];

        t_82[k] = -si_26[k]
                  + f_0 * di_166[k];

        t_83[k] = -si_27[k]
                  + f_0 * di_167[k];
    }
}

}  // namespace simdt2ceri
