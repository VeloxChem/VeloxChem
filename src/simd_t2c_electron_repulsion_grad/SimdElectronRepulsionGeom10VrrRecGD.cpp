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


#include "SimdElectronRepulsionGeom10VrrRecGD.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_gd_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t fd, const size_t hd,
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

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_24 = buffer.data(fd + 24);
    const auto *fd_25 = buffer.data(fd + 25);
    const auto *fd_26 = buffer.data(fd + 26);
    const auto *fd_27 = buffer.data(fd + 27);
    const auto *fd_28 = buffer.data(fd + 28);
    const auto *fd_29 = buffer.data(fd + 29);
    const auto *fd_30 = buffer.data(fd + 30);
    const auto *fd_31 = buffer.data(fd + 31);
    const auto *fd_32 = buffer.data(fd + 32);
    const auto *fd_33 = buffer.data(fd + 33);
    const auto *fd_34 = buffer.data(fd + 34);
    const auto *fd_35 = buffer.data(fd + 35);
    const auto *fd_36 = buffer.data(fd + 36);
    const auto *fd_37 = buffer.data(fd + 37);
    const auto *fd_38 = buffer.data(fd + 38);
    const auto *fd_39 = buffer.data(fd + 39);
    const auto *fd_40 = buffer.data(fd + 40);
    const auto *fd_41 = buffer.data(fd + 41);
    const auto *fd_42 = buffer.data(fd + 42);
    const auto *fd_43 = buffer.data(fd + 43);
    const auto *fd_44 = buffer.data(fd + 44);
    const auto *fd_45 = buffer.data(fd + 45);
    const auto *fd_46 = buffer.data(fd + 46);
    const auto *fd_47 = buffer.data(fd + 47);
    const auto *fd_48 = buffer.data(fd + 48);
    const auto *fd_49 = buffer.data(fd + 49);
    const auto *fd_50 = buffer.data(fd + 50);
    const auto *fd_51 = buffer.data(fd + 51);
    const auto *fd_52 = buffer.data(fd + 52);
    const auto *fd_53 = buffer.data(fd + 53);
    const auto *fd_54 = buffer.data(fd + 54);
    const auto *fd_55 = buffer.data(fd + 55);
    const auto *fd_56 = buffer.data(fd + 56);
    const auto *fd_57 = buffer.data(fd + 57);
    const auto *fd_58 = buffer.data(fd + 58);
    const auto *fd_59 = buffer.data(fd + 59);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_62 = buffer.data(hd + 62);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_64 = buffer.data(hd + 64);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_67 = buffer.data(hd + 67);
    const auto *hd_68 = buffer.data(hd + 68);
    const auto *hd_69 = buffer.data(hd + 69);
    const auto *hd_70 = buffer.data(hd + 70);
    const auto *hd_71 = buffer.data(hd + 71);
    const auto *hd_72 = buffer.data(hd + 72);
    const auto *hd_73 = buffer.data(hd + 73);
    const auto *hd_74 = buffer.data(hd + 74);
    const auto *hd_75 = buffer.data(hd + 75);
    const auto *hd_76 = buffer.data(hd + 76);
    const auto *hd_77 = buffer.data(hd + 77);
    const auto *hd_78 = buffer.data(hd + 78);
    const auto *hd_79 = buffer.data(hd + 79);
    const auto *hd_80 = buffer.data(hd + 80);
    const auto *hd_81 = buffer.data(hd + 81);
    const auto *hd_82 = buffer.data(hd + 82);
    const auto *hd_83 = buffer.data(hd + 83);
    const auto *hd_84 = buffer.data(hd + 84);
    const auto *hd_85 = buffer.data(hd + 85);
    const auto *hd_86 = buffer.data(hd + 86);
    const auto *hd_87 = buffer.data(hd + 87);
    const auto *hd_88 = buffer.data(hd + 88);
    const auto *hd_89 = buffer.data(hd + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, fd_0, fd_1, fd_2, fd_3, fd_4, hd_0, hd_1, \
                         hd_2, hd_3, hd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -4.0 * fd_0[k]
                 + f_0 * hd_0[k];

        t_1[k] = -4.0 * fd_1[k]
                 + f_0 * hd_1[k];

        t_2[k] = -4.0 * fd_2[k]
                 + f_0 * hd_2[k];

        t_3[k] = -4.0 * fd_3[k]
                 + f_0 * hd_3[k];

        t_4[k] = -4.0 * fd_4[k]
                 + f_0 * hd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, fd_5, fd_6, fd_7, fd_8, fd_9, hd_5, hd_6, \
                         hd_7, hd_8, hd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -4.0 * fd_5[k]
                 + f_0 * hd_5[k];

        t_6[k] = -3.0 * fd_6[k]
                 + f_0 * hd_6[k];

        t_7[k] = -3.0 * fd_7[k]
                 + f_0 * hd_7[k];

        t_8[k] = -3.0 * fd_8[k]
                 + f_0 * hd_8[k];

        t_9[k] = -3.0 * fd_9[k]
                 + f_0 * hd_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, fd_10, fd_11, fd_12, fd_13, fd_14, \
                         hd_10, hd_11, hd_12, hd_13, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -3.0 * fd_10[k]
                  + f_0 * hd_10[k];

        t_11[k] = -3.0 * fd_11[k]
                  + f_0 * hd_11[k];

        t_12[k] = -3.0 * fd_12[k]
                  + f_0 * hd_12[k];

        t_13[k] = -3.0 * fd_13[k]
                  + f_0 * hd_13[k];

        t_14[k] = -3.0 * fd_14[k]
                  + f_0 * hd_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, fd_15, fd_16, fd_17, fd_18, fd_19, \
                         hd_15, hd_16, hd_17, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -3.0 * fd_15[k]
                  + f_0 * hd_15[k];

        t_16[k] = -3.0 * fd_16[k]
                  + f_0 * hd_16[k];

        t_17[k] = -3.0 * fd_17[k]
                  + f_0 * hd_17[k];

        t_18[k] = -2.0 * fd_18[k]
                  + f_0 * hd_18[k];

        t_19[k] = -2.0 * fd_19[k]
                  + f_0 * hd_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, fd_20, fd_21, fd_22, fd_23, fd_24, \
                         hd_20, hd_21, hd_22, hd_23, hd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -2.0 * fd_20[k]
                  + f_0 * hd_20[k];

        t_21[k] = -2.0 * fd_21[k]
                  + f_0 * hd_21[k];

        t_22[k] = -2.0 * fd_22[k]
                  + f_0 * hd_22[k];

        t_23[k] = -2.0 * fd_23[k]
                  + f_0 * hd_23[k];

        t_24[k] = -2.0 * fd_24[k]
                  + f_0 * hd_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, fd_25, fd_26, fd_27, fd_28, fd_29, \
                         hd_25, hd_26, hd_27, hd_28, hd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -2.0 * fd_25[k]
                  + f_0 * hd_25[k];

        t_26[k] = -2.0 * fd_26[k]
                  + f_0 * hd_26[k];

        t_27[k] = -2.0 * fd_27[k]
                  + f_0 * hd_27[k];

        t_28[k] = -2.0 * fd_28[k]
                  + f_0 * hd_28[k];

        t_29[k] = -2.0 * fd_29[k]
                  + f_0 * hd_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, fd_30, fd_31, fd_32, fd_33, fd_34, \
                         hd_30, hd_31, hd_32, hd_33, hd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -2.0 * fd_30[k]
                  + f_0 * hd_30[k];

        t_31[k] = -2.0 * fd_31[k]
                  + f_0 * hd_31[k];

        t_32[k] = -2.0 * fd_32[k]
                  + f_0 * hd_32[k];

        t_33[k] = -2.0 * fd_33[k]
                  + f_0 * hd_33[k];

        t_34[k] = -2.0 * fd_34[k]
                  + f_0 * hd_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, fd_35, fd_36, fd_37, fd_38, fd_39, \
                         hd_35, hd_36, hd_37, hd_38, hd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -2.0 * fd_35[k]
                  + f_0 * hd_35[k];

        t_36[k] = -fd_36[k]
                  + f_0 * hd_36[k];

        t_37[k] = -fd_37[k]
                  + f_0 * hd_37[k];

        t_38[k] = -fd_38[k]
                  + f_0 * hd_38[k];

        t_39[k] = -fd_39[k]
                  + f_0 * hd_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, fd_40, fd_41, fd_42, fd_43, fd_44, \
                         hd_40, hd_41, hd_42, hd_43, hd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -fd_40[k]
                  + f_0 * hd_40[k];

        t_41[k] = -fd_41[k]
                  + f_0 * hd_41[k];

        t_42[k] = -fd_42[k]
                  + f_0 * hd_42[k];

        t_43[k] = -fd_43[k]
                  + f_0 * hd_43[k];

        t_44[k] = -fd_44[k]
                  + f_0 * hd_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, fd_45, fd_46, fd_47, fd_48, fd_49, \
                         hd_45, hd_46, hd_47, hd_48, hd_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -fd_45[k]
                  + f_0 * hd_45[k];

        t_46[k] = -fd_46[k]
                  + f_0 * hd_46[k];

        t_47[k] = -fd_47[k]
                  + f_0 * hd_47[k];

        t_48[k] = -fd_48[k]
                  + f_0 * hd_48[k];

        t_49[k] = -fd_49[k]
                  + f_0 * hd_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, fd_50, fd_51, fd_52, fd_53, fd_54, \
                         hd_50, hd_51, hd_52, hd_53, hd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -fd_50[k]
                  + f_0 * hd_50[k];

        t_51[k] = -fd_51[k]
                  + f_0 * hd_51[k];

        t_52[k] = -fd_52[k]
                  + f_0 * hd_52[k];

        t_53[k] = -fd_53[k]
                  + f_0 * hd_53[k];

        t_54[k] = -fd_54[k]
                  + f_0 * hd_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, fd_55, fd_56, fd_57, fd_58, fd_59, \
                         hd_55, hd_56, hd_57, hd_58, hd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -fd_55[k]
                  + f_0 * hd_55[k];

        t_56[k] = -fd_56[k]
                  + f_0 * hd_56[k];

        t_57[k] = -fd_57[k]
                  + f_0 * hd_57[k];

        t_58[k] = -fd_58[k]
                  + f_0 * hd_58[k];

        t_59[k] = -fd_59[k]
                  + f_0 * hd_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, t_66, t_67, hd_60, hd_61, hd_62, \
                         hd_63, hd_64, hd_65, hd_66, hd_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * hd_60[k];

        t_61[k] = f_0 * hd_61[k];

        t_62[k] = f_0 * hd_62[k];

        t_63[k] = f_0 * hd_63[k];

        t_64[k] = f_0 * hd_64[k];

        t_65[k] = f_0 * hd_65[k];

        t_66[k] = f_0 * hd_66[k];

        t_67[k] = f_0 * hd_67[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, t_74, t_75, hd_68, hd_69, hd_70, \
                         hd_71, hd_72, hd_73, hd_74, hd_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * hd_68[k];

        t_69[k] = f_0 * hd_69[k];

        t_70[k] = f_0 * hd_70[k];

        t_71[k] = f_0 * hd_71[k];

        t_72[k] = f_0 * hd_72[k];

        t_73[k] = f_0 * hd_73[k];

        t_74[k] = f_0 * hd_74[k];

        t_75[k] = f_0 * hd_75[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, t_81, t_82, t_83, hd_76, hd_77, hd_78, \
                         hd_79, hd_80, hd_81, hd_82, hd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_0 * hd_76[k];

        t_77[k] = f_0 * hd_77[k];

        t_78[k] = f_0 * hd_78[k];

        t_79[k] = f_0 * hd_79[k];

        t_80[k] = f_0 * hd_80[k];

        t_81[k] = f_0 * hd_81[k];

        t_82[k] = f_0 * hd_82[k];

        t_83[k] = f_0 * hd_83[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, hd_84, hd_85, hd_86, hd_87, \
                         hd_88, hd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_0 * hd_84[k];

        t_85[k] = f_0 * hd_85[k];

        t_86[k] = f_0 * hd_86[k];

        t_87[k] = f_0 * hd_87[k];

        t_88[k] = f_0 * hd_88[k];

        t_89[k] = f_0 * hd_89[k];
    }
}

auto
compute_prim_geom_10_gd_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t fd, const size_t hd,
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

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_24 = buffer.data(fd + 24);
    const auto *fd_25 = buffer.data(fd + 25);
    const auto *fd_26 = buffer.data(fd + 26);
    const auto *fd_27 = buffer.data(fd + 27);
    const auto *fd_28 = buffer.data(fd + 28);
    const auto *fd_29 = buffer.data(fd + 29);
    const auto *fd_30 = buffer.data(fd + 30);
    const auto *fd_31 = buffer.data(fd + 31);
    const auto *fd_32 = buffer.data(fd + 32);
    const auto *fd_33 = buffer.data(fd + 33);
    const auto *fd_34 = buffer.data(fd + 34);
    const auto *fd_35 = buffer.data(fd + 35);
    const auto *fd_36 = buffer.data(fd + 36);
    const auto *fd_37 = buffer.data(fd + 37);
    const auto *fd_38 = buffer.data(fd + 38);
    const auto *fd_39 = buffer.data(fd + 39);
    const auto *fd_40 = buffer.data(fd + 40);
    const auto *fd_41 = buffer.data(fd + 41);
    const auto *fd_42 = buffer.data(fd + 42);
    const auto *fd_43 = buffer.data(fd + 43);
    const auto *fd_44 = buffer.data(fd + 44);
    const auto *fd_45 = buffer.data(fd + 45);
    const auto *fd_46 = buffer.data(fd + 46);
    const auto *fd_47 = buffer.data(fd + 47);
    const auto *fd_48 = buffer.data(fd + 48);
    const auto *fd_49 = buffer.data(fd + 49);
    const auto *fd_50 = buffer.data(fd + 50);
    const auto *fd_51 = buffer.data(fd + 51);
    const auto *fd_52 = buffer.data(fd + 52);
    const auto *fd_53 = buffer.data(fd + 53);
    const auto *fd_54 = buffer.data(fd + 54);
    const auto *fd_55 = buffer.data(fd + 55);
    const auto *fd_56 = buffer.data(fd + 56);
    const auto *fd_57 = buffer.data(fd + 57);
    const auto *fd_58 = buffer.data(fd + 58);
    const auto *fd_59 = buffer.data(fd + 59);

    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_62 = buffer.data(hd + 62);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_64 = buffer.data(hd + 64);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_67 = buffer.data(hd + 67);
    const auto *hd_68 = buffer.data(hd + 68);
    const auto *hd_69 = buffer.data(hd + 69);
    const auto *hd_70 = buffer.data(hd + 70);
    const auto *hd_71 = buffer.data(hd + 71);
    const auto *hd_72 = buffer.data(hd + 72);
    const auto *hd_73 = buffer.data(hd + 73);
    const auto *hd_74 = buffer.data(hd + 74);
    const auto *hd_75 = buffer.data(hd + 75);
    const auto *hd_76 = buffer.data(hd + 76);
    const auto *hd_77 = buffer.data(hd + 77);
    const auto *hd_78 = buffer.data(hd + 78);
    const auto *hd_79 = buffer.data(hd + 79);
    const auto *hd_80 = buffer.data(hd + 80);
    const auto *hd_81 = buffer.data(hd + 81);
    const auto *hd_82 = buffer.data(hd + 82);
    const auto *hd_83 = buffer.data(hd + 83);
    const auto *hd_90 = buffer.data(hd + 90);
    const auto *hd_91 = buffer.data(hd + 91);
    const auto *hd_92 = buffer.data(hd + 92);
    const auto *hd_93 = buffer.data(hd + 93);
    const auto *hd_94 = buffer.data(hd + 94);
    const auto *hd_95 = buffer.data(hd + 95);
    const auto *hd_96 = buffer.data(hd + 96);
    const auto *hd_97 = buffer.data(hd + 97);
    const auto *hd_98 = buffer.data(hd + 98);
    const auto *hd_99 = buffer.data(hd + 99);
    const auto *hd_100 = buffer.data(hd + 100);
    const auto *hd_101 = buffer.data(hd + 101);
    const auto *hd_102 = buffer.data(hd + 102);
    const auto *hd_103 = buffer.data(hd + 103);
    const auto *hd_104 = buffer.data(hd + 104);
    const auto *hd_105 = buffer.data(hd + 105);
    const auto *hd_106 = buffer.data(hd + 106);
    const auto *hd_107 = buffer.data(hd + 107);
    const auto *hd_108 = buffer.data(hd + 108);
    const auto *hd_109 = buffer.data(hd + 109);
    const auto *hd_110 = buffer.data(hd + 110);
    const auto *hd_111 = buffer.data(hd + 111);
    const auto *hd_112 = buffer.data(hd + 112);
    const auto *hd_113 = buffer.data(hd + 113);
    const auto *hd_114 = buffer.data(hd + 114);
    const auto *hd_115 = buffer.data(hd + 115);
    const auto *hd_116 = buffer.data(hd + 116);
    const auto *hd_117 = buffer.data(hd + 117);
    const auto *hd_118 = buffer.data(hd + 118);
    const auto *hd_119 = buffer.data(hd + 119);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, fd_0, hd_6, hd_7, hd_8, hd_9, \
                         hd_10, hd_11, hd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_6[k];

        t_1[k] = f_0 * hd_7[k];

        t_2[k] = f_0 * hd_8[k];

        t_3[k] = f_0 * hd_9[k];

        t_4[k] = f_0 * hd_10[k];

        t_5[k] = f_0 * hd_11[k];

        t_6[k] = -fd_0[k]
                 + f_0 * hd_18[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, fd_1, fd_2, fd_3, fd_4, fd_5, hd_19, \
                         hd_20, hd_21, hd_22, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -fd_1[k]
                 + f_0 * hd_19[k];

        t_8[k] = -fd_2[k]
                 + f_0 * hd_20[k];

        t_9[k] = -fd_3[k]
                 + f_0 * hd_21[k];

        t_10[k] = -fd_4[k]
                  + f_0 * hd_22[k];

        t_11[k] = -fd_5[k]
                  + f_0 * hd_23[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, t_18, fd_6, hd_24, hd_25, hd_26, \
                         hd_27, hd_28, hd_29, hd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * hd_24[k];

        t_13[k] = f_0 * hd_25[k];

        t_14[k] = f_0 * hd_26[k];

        t_15[k] = f_0 * hd_27[k];

        t_16[k] = f_0 * hd_28[k];

        t_17[k] = f_0 * hd_29[k];

        t_18[k] = -2.0 * fd_6[k]
                  + f_0 * hd_36[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, fd_7, fd_8, fd_9, fd_10, fd_11, hd_37, \
                         hd_38, hd_39, hd_40, hd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -2.0 * fd_7[k]
                  + f_0 * hd_37[k];

        t_20[k] = -2.0 * fd_8[k]
                  + f_0 * hd_38[k];

        t_21[k] = -2.0 * fd_9[k]
                  + f_0 * hd_39[k];

        t_22[k] = -2.0 * fd_10[k]
                  + f_0 * hd_40[k];

        t_23[k] = -2.0 * fd_11[k]
                  + f_0 * hd_41[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, fd_12, fd_13, fd_14, fd_15, fd_16, \
                         hd_42, hd_43, hd_44, hd_45, hd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -fd_12[k]
                  + f_0 * hd_42[k];

        t_25[k] = -fd_13[k]
                  + f_0 * hd_43[k];

        t_26[k] = -fd_14[k]
                  + f_0 * hd_44[k];

        t_27[k] = -fd_15[k]
                  + f_0 * hd_45[k];

        t_28[k] = -fd_16[k]
                  + f_0 * hd_46[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, t_35, fd_17, hd_47, hd_48, hd_49, \
                         hd_50, hd_51, hd_52, hd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -fd_17[k]
                  + f_0 * hd_47[k];

        t_30[k] = f_0 * hd_48[k];

        t_31[k] = f_0 * hd_49[k];

        t_32[k] = f_0 * hd_50[k];

        t_33[k] = f_0 * hd_51[k];

        t_34[k] = f_0 * hd_52[k];

        t_35[k] = f_0 * hd_53[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, fd_18, fd_19, fd_20, fd_21, fd_22, \
                         hd_60, hd_61, hd_62, hd_63, hd_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -3.0 * fd_18[k]
                  + f_0 * hd_60[k];

        t_37[k] = -3.0 * fd_19[k]
                  + f_0 * hd_61[k];

        t_38[k] = -3.0 * fd_20[k]
                  + f_0 * hd_62[k];

        t_39[k] = -3.0 * fd_21[k]
                  + f_0 * hd_63[k];

        t_40[k] = -3.0 * fd_22[k]
                  + f_0 * hd_64[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, fd_23, fd_24, fd_25, fd_26, fd_27, \
                         hd_65, hd_66, hd_67, hd_68, hd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -3.0 * fd_23[k]
                  + f_0 * hd_65[k];

        t_42[k] = -2.0 * fd_24[k]
                  + f_0 * hd_66[k];

        t_43[k] = -2.0 * fd_25[k]
                  + f_0 * hd_67[k];

        t_44[k] = -2.0 * fd_26[k]
                  + f_0 * hd_68[k];

        t_45[k] = -2.0 * fd_27[k]
                  + f_0 * hd_69[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, fd_28, fd_29, fd_30, fd_31, fd_32, \
                         hd_70, hd_71, hd_72, hd_73, hd_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -2.0 * fd_28[k]
                  + f_0 * hd_70[k];

        t_47[k] = -2.0 * fd_29[k]
                  + f_0 * hd_71[k];

        t_48[k] = -fd_30[k]
                  + f_0 * hd_72[k];

        t_49[k] = -fd_31[k]
                  + f_0 * hd_73[k];

        t_50[k] = -fd_32[k]
                  + f_0 * hd_74[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, fd_33, fd_34, fd_35, hd_75, \
                         hd_76, hd_77, hd_78, hd_79, hd_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -fd_33[k]
                  + f_0 * hd_75[k];

        t_52[k] = -fd_34[k]
                  + f_0 * hd_76[k];

        t_53[k] = -fd_35[k]
                  + f_0 * hd_77[k];

        t_54[k] = f_0 * hd_78[k];

        t_55[k] = f_0 * hd_79[k];

        t_56[k] = f_0 * hd_80[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, t_62, fd_36, fd_37, fd_38, hd_81, \
                         hd_82, hd_83, hd_90, hd_91, hd_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_0 * hd_81[k];

        t_58[k] = f_0 * hd_82[k];

        t_59[k] = f_0 * hd_83[k];

        t_60[k] = -4.0 * fd_36[k]
                  + f_0 * hd_90[k];

        t_61[k] = -4.0 * fd_37[k]
                  + f_0 * hd_91[k];

        t_62[k] = -4.0 * fd_38[k]
                  + f_0 * hd_92[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, fd_39, fd_40, fd_41, fd_42, fd_43, \
                         hd_93, hd_94, hd_95, hd_96, hd_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -4.0 * fd_39[k]
                  + f_0 * hd_93[k];

        t_64[k] = -4.0 * fd_40[k]
                  + f_0 * hd_94[k];

        t_65[k] = -4.0 * fd_41[k]
                  + f_0 * hd_95[k];

        t_66[k] = -3.0 * fd_42[k]
                  + f_0 * hd_96[k];

        t_67[k] = -3.0 * fd_43[k]
                  + f_0 * hd_97[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, fd_44, fd_45, fd_46, fd_47, fd_48, \
                         hd_98, hd_99, hd_100, hd_101, hd_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -3.0 * fd_44[k]
                  + f_0 * hd_98[k];

        t_69[k] = -3.0 * fd_45[k]
                  + f_0 * hd_99[k];

        t_70[k] = -3.0 * fd_46[k]
                  + f_0 * hd_100[k];

        t_71[k] = -3.0 * fd_47[k]
                  + f_0 * hd_101[k];

        t_72[k] = -2.0 * fd_48[k]
                  + f_0 * hd_102[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, fd_49, fd_50, fd_51, fd_52, fd_53, \
                         hd_103, hd_104, hd_105, hd_106, hd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = -2.0 * fd_49[k]
                  + f_0 * hd_103[k];

        t_74[k] = -2.0 * fd_50[k]
                  + f_0 * hd_104[k];

        t_75[k] = -2.0 * fd_51[k]
                  + f_0 * hd_105[k];

        t_76[k] = -2.0 * fd_52[k]
                  + f_0 * hd_106[k];

        t_77[k] = -2.0 * fd_53[k]
                  + f_0 * hd_107[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, fd_54, fd_55, fd_56, fd_57, fd_58, \
                         hd_108, hd_109, hd_110, hd_111, hd_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = -fd_54[k]
                  + f_0 * hd_108[k];

        t_79[k] = -fd_55[k]
                  + f_0 * hd_109[k];

        t_80[k] = -fd_56[k]
                  + f_0 * hd_110[k];

        t_81[k] = -fd_57[k]
                  + f_0 * hd_111[k];

        t_82[k] = -fd_58[k]
                  + f_0 * hd_112[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, t_89, fd_59, hd_113, hd_114, \
                         hd_115, hd_116, hd_117, hd_118, hd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = -fd_59[k]
                  + f_0 * hd_113[k];

        t_84[k] = f_0 * hd_114[k];

        t_85[k] = f_0 * hd_115[k];

        t_86[k] = f_0 * hd_116[k];

        t_87[k] = f_0 * hd_117[k];

        t_88[k] = f_0 * hd_118[k];

        t_89[k] = f_0 * hd_119[k];
    }
}

auto
compute_prim_geom_10_gd_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t fd, const size_t hd,
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

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_24 = buffer.data(fd + 24);
    const auto *fd_25 = buffer.data(fd + 25);
    const auto *fd_26 = buffer.data(fd + 26);
    const auto *fd_27 = buffer.data(fd + 27);
    const auto *fd_28 = buffer.data(fd + 28);
    const auto *fd_29 = buffer.data(fd + 29);
    const auto *fd_30 = buffer.data(fd + 30);
    const auto *fd_31 = buffer.data(fd + 31);
    const auto *fd_32 = buffer.data(fd + 32);
    const auto *fd_33 = buffer.data(fd + 33);
    const auto *fd_34 = buffer.data(fd + 34);
    const auto *fd_35 = buffer.data(fd + 35);
    const auto *fd_36 = buffer.data(fd + 36);
    const auto *fd_37 = buffer.data(fd + 37);
    const auto *fd_38 = buffer.data(fd + 38);
    const auto *fd_39 = buffer.data(fd + 39);
    const auto *fd_40 = buffer.data(fd + 40);
    const auto *fd_41 = buffer.data(fd + 41);
    const auto *fd_42 = buffer.data(fd + 42);
    const auto *fd_43 = buffer.data(fd + 43);
    const auto *fd_44 = buffer.data(fd + 44);
    const auto *fd_45 = buffer.data(fd + 45);
    const auto *fd_46 = buffer.data(fd + 46);
    const auto *fd_47 = buffer.data(fd + 47);
    const auto *fd_48 = buffer.data(fd + 48);
    const auto *fd_49 = buffer.data(fd + 49);
    const auto *fd_50 = buffer.data(fd + 50);
    const auto *fd_51 = buffer.data(fd + 51);
    const auto *fd_52 = buffer.data(fd + 52);
    const auto *fd_53 = buffer.data(fd + 53);
    const auto *fd_54 = buffer.data(fd + 54);
    const auto *fd_55 = buffer.data(fd + 55);
    const auto *fd_56 = buffer.data(fd + 56);
    const auto *fd_57 = buffer.data(fd + 57);
    const auto *fd_58 = buffer.data(fd + 58);
    const auto *fd_59 = buffer.data(fd + 59);

    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_67 = buffer.data(hd + 67);
    const auto *hd_68 = buffer.data(hd + 68);
    const auto *hd_69 = buffer.data(hd + 69);
    const auto *hd_70 = buffer.data(hd + 70);
    const auto *hd_71 = buffer.data(hd + 71);
    const auto *hd_72 = buffer.data(hd + 72);
    const auto *hd_73 = buffer.data(hd + 73);
    const auto *hd_74 = buffer.data(hd + 74);
    const auto *hd_75 = buffer.data(hd + 75);
    const auto *hd_76 = buffer.data(hd + 76);
    const auto *hd_77 = buffer.data(hd + 77);
    const auto *hd_78 = buffer.data(hd + 78);
    const auto *hd_79 = buffer.data(hd + 79);
    const auto *hd_80 = buffer.data(hd + 80);
    const auto *hd_81 = buffer.data(hd + 81);
    const auto *hd_82 = buffer.data(hd + 82);
    const auto *hd_83 = buffer.data(hd + 83);
    const auto *hd_84 = buffer.data(hd + 84);
    const auto *hd_85 = buffer.data(hd + 85);
    const auto *hd_86 = buffer.data(hd + 86);
    const auto *hd_87 = buffer.data(hd + 87);
    const auto *hd_88 = buffer.data(hd + 88);
    const auto *hd_89 = buffer.data(hd + 89);
    const auto *hd_96 = buffer.data(hd + 96);
    const auto *hd_97 = buffer.data(hd + 97);
    const auto *hd_98 = buffer.data(hd + 98);
    const auto *hd_99 = buffer.data(hd + 99);
    const auto *hd_100 = buffer.data(hd + 100);
    const auto *hd_101 = buffer.data(hd + 101);
    const auto *hd_102 = buffer.data(hd + 102);
    const auto *hd_103 = buffer.data(hd + 103);
    const auto *hd_104 = buffer.data(hd + 104);
    const auto *hd_105 = buffer.data(hd + 105);
    const auto *hd_106 = buffer.data(hd + 106);
    const auto *hd_107 = buffer.data(hd + 107);
    const auto *hd_108 = buffer.data(hd + 108);
    const auto *hd_109 = buffer.data(hd + 109);
    const auto *hd_110 = buffer.data(hd + 110);
    const auto *hd_111 = buffer.data(hd + 111);
    const auto *hd_112 = buffer.data(hd + 112);
    const auto *hd_113 = buffer.data(hd + 113);
    const auto *hd_114 = buffer.data(hd + 114);
    const auto *hd_115 = buffer.data(hd + 115);
    const auto *hd_116 = buffer.data(hd + 116);
    const auto *hd_117 = buffer.data(hd + 117);
    const auto *hd_118 = buffer.data(hd + 118);
    const auto *hd_119 = buffer.data(hd + 119);
    const auto *hd_120 = buffer.data(hd + 120);
    const auto *hd_121 = buffer.data(hd + 121);
    const auto *hd_122 = buffer.data(hd + 122);
    const auto *hd_123 = buffer.data(hd + 123);
    const auto *hd_124 = buffer.data(hd + 124);
    const auto *hd_125 = buffer.data(hd + 125);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, hd_12, hd_13, hd_14, hd_15, \
                         hd_16, hd_17, hd_24, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_12[k];

        t_1[k] = f_0 * hd_13[k];

        t_2[k] = f_0 * hd_14[k];

        t_3[k] = f_0 * hd_15[k];

        t_4[k] = f_0 * hd_16[k];

        t_5[k] = f_0 * hd_17[k];

        t_6[k] = f_0 * hd_24[k];

        t_7[k] = f_0 * hd_25[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, fd_0, fd_1, hd_26, hd_27, hd_28, \
                         hd_29, hd_30, hd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hd_26[k];

        t_9[k] = f_0 * hd_27[k];

        t_10[k] = f_0 * hd_28[k];

        t_11[k] = f_0 * hd_29[k];

        t_12[k] = -fd_0[k]
                  + f_0 * hd_30[k];

        t_13[k] = -fd_1[k]
                  + f_0 * hd_31[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, fd_2, fd_3, fd_4, fd_5, hd_32, \
                         hd_33, hd_34, hd_35, hd_42, hd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -fd_2[k]
                  + f_0 * hd_32[k];

        t_15[k] = -fd_3[k]
                  + f_0 * hd_33[k];

        t_16[k] = -fd_4[k]
                  + f_0 * hd_34[k];

        t_17[k] = -fd_5[k]
                  + f_0 * hd_35[k];

        t_18[k] = f_0 * hd_42[k];

        t_19[k] = f_0 * hd_43[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, fd_6, fd_7, hd_44, hd_45, hd_46, \
                         hd_47, hd_48, hd_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * hd_44[k];

        t_21[k] = f_0 * hd_45[k];

        t_22[k] = f_0 * hd_46[k];

        t_23[k] = f_0 * hd_47[k];

        t_24[k] = -fd_6[k]
                  + f_0 * hd_48[k];

        t_25[k] = -fd_7[k]
                  + f_0 * hd_49[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, fd_8, fd_9, fd_10, fd_11, fd_12, hd_50, \
                         hd_51, hd_52, hd_53, hd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -fd_8[k]
                  + f_0 * hd_50[k];

        t_27[k] = -fd_9[k]
                  + f_0 * hd_51[k];

        t_28[k] = -fd_10[k]
                  + f_0 * hd_52[k];

        t_29[k] = -fd_11[k]
                  + f_0 * hd_53[k];

        t_30[k] = -2.0 * fd_12[k]
                  + f_0 * hd_54[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, fd_13, fd_14, fd_15, fd_16, fd_17, \
                         hd_55, hd_56, hd_57, hd_58, hd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -2.0 * fd_13[k]
                  + f_0 * hd_55[k];

        t_32[k] = -2.0 * fd_14[k]
                  + f_0 * hd_56[k];

        t_33[k] = -2.0 * fd_15[k]
                  + f_0 * hd_57[k];

        t_34[k] = -2.0 * fd_16[k]
                  + f_0 * hd_58[k];

        t_35[k] = -2.0 * fd_17[k]
                  + f_0 * hd_59[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, t_42, fd_18, hd_66, hd_67, hd_68, \
                         hd_69, hd_70, hd_71, hd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * hd_66[k];

        t_37[k] = f_0 * hd_67[k];

        t_38[k] = f_0 * hd_68[k];

        t_39[k] = f_0 * hd_69[k];

        t_40[k] = f_0 * hd_70[k];

        t_41[k] = f_0 * hd_71[k];

        t_42[k] = -fd_18[k]
                  + f_0 * hd_72[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, fd_19, fd_20, fd_21, fd_22, fd_23, \
                         hd_73, hd_74, hd_75, hd_76, hd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -fd_19[k]
                  + f_0 * hd_73[k];

        t_44[k] = -fd_20[k]
                  + f_0 * hd_74[k];

        t_45[k] = -fd_21[k]
                  + f_0 * hd_75[k];

        t_46[k] = -fd_22[k]
                  + f_0 * hd_76[k];

        t_47[k] = -fd_23[k]
                  + f_0 * hd_77[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, fd_24, fd_25, fd_26, fd_27, fd_28, \
                         hd_78, hd_79, hd_80, hd_81, hd_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -2.0 * fd_24[k]
                  + f_0 * hd_78[k];

        t_49[k] = -2.0 * fd_25[k]
                  + f_0 * hd_79[k];

        t_50[k] = -2.0 * fd_26[k]
                  + f_0 * hd_80[k];

        t_51[k] = -2.0 * fd_27[k]
                  + f_0 * hd_81[k];

        t_52[k] = -2.0 * fd_28[k]
                  + f_0 * hd_82[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, fd_29, fd_30, fd_31, fd_32, fd_33, \
                         hd_83, hd_84, hd_85, hd_86, hd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -2.0 * fd_29[k]
                  + f_0 * hd_83[k];

        t_54[k] = -3.0 * fd_30[k]
                  + f_0 * hd_84[k];

        t_55[k] = -3.0 * fd_31[k]
                  + f_0 * hd_85[k];

        t_56[k] = -3.0 * fd_32[k]
                  + f_0 * hd_86[k];

        t_57[k] = -3.0 * fd_33[k]
                  + f_0 * hd_87[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, t_63, t_64, fd_34, fd_35, hd_88, hd_89, \
                         hd_96, hd_97, hd_98, hd_99, hd_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -3.0 * fd_34[k]
                  + f_0 * hd_88[k];

        t_59[k] = -3.0 * fd_35[k]
                  + f_0 * hd_89[k];

        t_60[k] = f_0 * hd_96[k];

        t_61[k] = f_0 * hd_97[k];

        t_62[k] = f_0 * hd_98[k];

        t_63[k] = f_0 * hd_99[k];

        t_64[k] = f_0 * hd_100[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, fd_36, fd_37, fd_38, fd_39, hd_101, \
                         hd_102, hd_103, hd_104, hd_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_0 * hd_101[k];

        t_66[k] = -fd_36[k]
                  + f_0 * hd_102[k];

        t_67[k] = -fd_37[k]
                  + f_0 * hd_103[k];

        t_68[k] = -fd_38[k]
                  + f_0 * hd_104[k];

        t_69[k] = -fd_39[k]
                  + f_0 * hd_105[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, fd_40, fd_41, fd_42, fd_43, fd_44, \
                         hd_106, hd_107, hd_108, hd_109, hd_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -fd_40[k]
                  + f_0 * hd_106[k];

        t_71[k] = -fd_41[k]
                  + f_0 * hd_107[k];

        t_72[k] = -2.0 * fd_42[k]
                  + f_0 * hd_108[k];

        t_73[k] = -2.0 * fd_43[k]
                  + f_0 * hd_109[k];

        t_74[k] = -2.0 * fd_44[k]
                  + f_0 * hd_110[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, fd_45, fd_46, fd_47, fd_48, fd_49, \
                         hd_111, hd_112, hd_113, hd_114, hd_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * fd_45[k]
                  + f_0 * hd_111[k];

        t_76[k] = -2.0 * fd_46[k]
                  + f_0 * hd_112[k];

        t_77[k] = -2.0 * fd_47[k]
                  + f_0 * hd_113[k];

        t_78[k] = -3.0 * fd_48[k]
                  + f_0 * hd_114[k];

        t_79[k] = -3.0 * fd_49[k]
                  + f_0 * hd_115[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, fd_50, fd_51, fd_52, fd_53, fd_54, \
                         hd_116, hd_117, hd_118, hd_119, hd_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -3.0 * fd_50[k]
                  + f_0 * hd_116[k];

        t_81[k] = -3.0 * fd_51[k]
                  + f_0 * hd_117[k];

        t_82[k] = -3.0 * fd_52[k]
                  + f_0 * hd_118[k];

        t_83[k] = -3.0 * fd_53[k]
                  + f_0 * hd_119[k];

        t_84[k] = -4.0 * fd_54[k]
                  + f_0 * hd_120[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, fd_55, fd_56, fd_57, fd_58, fd_59, \
                         hd_121, hd_122, hd_123, hd_124, hd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -4.0 * fd_55[k]
                  + f_0 * hd_121[k];

        t_86[k] = -4.0 * fd_56[k]
                  + f_0 * hd_122[k];

        t_87[k] = -4.0 * fd_57[k]
                  + f_0 * hd_123[k];

        t_88[k] = -4.0 * fd_58[k]
                  + f_0 * hd_124[k];

        t_89[k] = -4.0 * fd_59[k]
                  + f_0 * hd_125[k];
    }
}

}  // namespace simdt2ceri
