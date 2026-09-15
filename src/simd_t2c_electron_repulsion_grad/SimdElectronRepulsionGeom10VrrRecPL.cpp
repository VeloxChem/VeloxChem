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


#include "SimdElectronRepulsionGeom10VrrRecPL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_pl_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t sl, const size_t dl,
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
    auto *t_108 = buffer.data(target + 108);
    auto *t_109 = buffer.data(target + 109);
    auto *t_110 = buffer.data(target + 110);
    auto *t_111 = buffer.data(target + 111);
    auto *t_112 = buffer.data(target + 112);
    auto *t_113 = buffer.data(target + 113);
    auto *t_114 = buffer.data(target + 114);
    auto *t_115 = buffer.data(target + 115);
    auto *t_116 = buffer.data(target + 116);
    auto *t_117 = buffer.data(target + 117);
    auto *t_118 = buffer.data(target + 118);
    auto *t_119 = buffer.data(target + 119);
    auto *t_120 = buffer.data(target + 120);
    auto *t_121 = buffer.data(target + 121);
    auto *t_122 = buffer.data(target + 122);
    auto *t_123 = buffer.data(target + 123);
    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);

    const auto *sl_0 = buffer.data(sl + 0);
    const auto *sl_1 = buffer.data(sl + 1);
    const auto *sl_2 = buffer.data(sl + 2);
    const auto *sl_3 = buffer.data(sl + 3);
    const auto *sl_4 = buffer.data(sl + 4);
    const auto *sl_5 = buffer.data(sl + 5);
    const auto *sl_6 = buffer.data(sl + 6);
    const auto *sl_7 = buffer.data(sl + 7);
    const auto *sl_8 = buffer.data(sl + 8);
    const auto *sl_9 = buffer.data(sl + 9);
    const auto *sl_10 = buffer.data(sl + 10);
    const auto *sl_11 = buffer.data(sl + 11);
    const auto *sl_12 = buffer.data(sl + 12);
    const auto *sl_13 = buffer.data(sl + 13);
    const auto *sl_14 = buffer.data(sl + 14);
    const auto *sl_15 = buffer.data(sl + 15);
    const auto *sl_16 = buffer.data(sl + 16);
    const auto *sl_17 = buffer.data(sl + 17);
    const auto *sl_18 = buffer.data(sl + 18);
    const auto *sl_19 = buffer.data(sl + 19);
    const auto *sl_20 = buffer.data(sl + 20);
    const auto *sl_21 = buffer.data(sl + 21);
    const auto *sl_22 = buffer.data(sl + 22);
    const auto *sl_23 = buffer.data(sl + 23);
    const auto *sl_24 = buffer.data(sl + 24);
    const auto *sl_25 = buffer.data(sl + 25);
    const auto *sl_26 = buffer.data(sl + 26);
    const auto *sl_27 = buffer.data(sl + 27);
    const auto *sl_28 = buffer.data(sl + 28);
    const auto *sl_29 = buffer.data(sl + 29);
    const auto *sl_30 = buffer.data(sl + 30);
    const auto *sl_31 = buffer.data(sl + 31);
    const auto *sl_32 = buffer.data(sl + 32);
    const auto *sl_33 = buffer.data(sl + 33);
    const auto *sl_34 = buffer.data(sl + 34);
    const auto *sl_35 = buffer.data(sl + 35);
    const auto *sl_36 = buffer.data(sl + 36);
    const auto *sl_37 = buffer.data(sl + 37);
    const auto *sl_38 = buffer.data(sl + 38);
    const auto *sl_39 = buffer.data(sl + 39);
    const auto *sl_40 = buffer.data(sl + 40);
    const auto *sl_41 = buffer.data(sl + 41);
    const auto *sl_42 = buffer.data(sl + 42);
    const auto *sl_43 = buffer.data(sl + 43);
    const auto *sl_44 = buffer.data(sl + 44);

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);
    const auto *dl_3 = buffer.data(dl + 3);
    const auto *dl_4 = buffer.data(dl + 4);
    const auto *dl_5 = buffer.data(dl + 5);
    const auto *dl_6 = buffer.data(dl + 6);
    const auto *dl_7 = buffer.data(dl + 7);
    const auto *dl_8 = buffer.data(dl + 8);
    const auto *dl_9 = buffer.data(dl + 9);
    const auto *dl_10 = buffer.data(dl + 10);
    const auto *dl_11 = buffer.data(dl + 11);
    const auto *dl_12 = buffer.data(dl + 12);
    const auto *dl_13 = buffer.data(dl + 13);
    const auto *dl_14 = buffer.data(dl + 14);
    const auto *dl_15 = buffer.data(dl + 15);
    const auto *dl_16 = buffer.data(dl + 16);
    const auto *dl_17 = buffer.data(dl + 17);
    const auto *dl_18 = buffer.data(dl + 18);
    const auto *dl_19 = buffer.data(dl + 19);
    const auto *dl_20 = buffer.data(dl + 20);
    const auto *dl_21 = buffer.data(dl + 21);
    const auto *dl_22 = buffer.data(dl + 22);
    const auto *dl_23 = buffer.data(dl + 23);
    const auto *dl_24 = buffer.data(dl + 24);
    const auto *dl_25 = buffer.data(dl + 25);
    const auto *dl_26 = buffer.data(dl + 26);
    const auto *dl_27 = buffer.data(dl + 27);
    const auto *dl_28 = buffer.data(dl + 28);
    const auto *dl_29 = buffer.data(dl + 29);
    const auto *dl_30 = buffer.data(dl + 30);
    const auto *dl_31 = buffer.data(dl + 31);
    const auto *dl_32 = buffer.data(dl + 32);
    const auto *dl_33 = buffer.data(dl + 33);
    const auto *dl_34 = buffer.data(dl + 34);
    const auto *dl_35 = buffer.data(dl + 35);
    const auto *dl_36 = buffer.data(dl + 36);
    const auto *dl_37 = buffer.data(dl + 37);
    const auto *dl_38 = buffer.data(dl + 38);
    const auto *dl_39 = buffer.data(dl + 39);
    const auto *dl_40 = buffer.data(dl + 40);
    const auto *dl_41 = buffer.data(dl + 41);
    const auto *dl_42 = buffer.data(dl + 42);
    const auto *dl_43 = buffer.data(dl + 43);
    const auto *dl_44 = buffer.data(dl + 44);
    const auto *dl_45 = buffer.data(dl + 45);
    const auto *dl_46 = buffer.data(dl + 46);
    const auto *dl_47 = buffer.data(dl + 47);
    const auto *dl_48 = buffer.data(dl + 48);
    const auto *dl_49 = buffer.data(dl + 49);
    const auto *dl_50 = buffer.data(dl + 50);
    const auto *dl_51 = buffer.data(dl + 51);
    const auto *dl_52 = buffer.data(dl + 52);
    const auto *dl_53 = buffer.data(dl + 53);
    const auto *dl_54 = buffer.data(dl + 54);
    const auto *dl_55 = buffer.data(dl + 55);
    const auto *dl_56 = buffer.data(dl + 56);
    const auto *dl_57 = buffer.data(dl + 57);
    const auto *dl_58 = buffer.data(dl + 58);
    const auto *dl_59 = buffer.data(dl + 59);
    const auto *dl_60 = buffer.data(dl + 60);
    const auto *dl_61 = buffer.data(dl + 61);
    const auto *dl_62 = buffer.data(dl + 62);
    const auto *dl_63 = buffer.data(dl + 63);
    const auto *dl_64 = buffer.data(dl + 64);
    const auto *dl_65 = buffer.data(dl + 65);
    const auto *dl_66 = buffer.data(dl + 66);
    const auto *dl_67 = buffer.data(dl + 67);
    const auto *dl_68 = buffer.data(dl + 68);
    const auto *dl_69 = buffer.data(dl + 69);
    const auto *dl_70 = buffer.data(dl + 70);
    const auto *dl_71 = buffer.data(dl + 71);
    const auto *dl_72 = buffer.data(dl + 72);
    const auto *dl_73 = buffer.data(dl + 73);
    const auto *dl_74 = buffer.data(dl + 74);
    const auto *dl_75 = buffer.data(dl + 75);
    const auto *dl_76 = buffer.data(dl + 76);
    const auto *dl_77 = buffer.data(dl + 77);
    const auto *dl_78 = buffer.data(dl + 78);
    const auto *dl_79 = buffer.data(dl + 79);
    const auto *dl_80 = buffer.data(dl + 80);
    const auto *dl_81 = buffer.data(dl + 81);
    const auto *dl_82 = buffer.data(dl + 82);
    const auto *dl_83 = buffer.data(dl + 83);
    const auto *dl_84 = buffer.data(dl + 84);
    const auto *dl_85 = buffer.data(dl + 85);
    const auto *dl_86 = buffer.data(dl + 86);
    const auto *dl_87 = buffer.data(dl + 87);
    const auto *dl_88 = buffer.data(dl + 88);
    const auto *dl_89 = buffer.data(dl + 89);
    const auto *dl_90 = buffer.data(dl + 90);
    const auto *dl_91 = buffer.data(dl + 91);
    const auto *dl_92 = buffer.data(dl + 92);
    const auto *dl_93 = buffer.data(dl + 93);
    const auto *dl_94 = buffer.data(dl + 94);
    const auto *dl_95 = buffer.data(dl + 95);
    const auto *dl_96 = buffer.data(dl + 96);
    const auto *dl_97 = buffer.data(dl + 97);
    const auto *dl_98 = buffer.data(dl + 98);
    const auto *dl_99 = buffer.data(dl + 99);
    const auto *dl_100 = buffer.data(dl + 100);
    const auto *dl_101 = buffer.data(dl + 101);
    const auto *dl_102 = buffer.data(dl + 102);
    const auto *dl_103 = buffer.data(dl + 103);
    const auto *dl_104 = buffer.data(dl + 104);
    const auto *dl_105 = buffer.data(dl + 105);
    const auto *dl_106 = buffer.data(dl + 106);
    const auto *dl_107 = buffer.data(dl + 107);
    const auto *dl_108 = buffer.data(dl + 108);
    const auto *dl_109 = buffer.data(dl + 109);
    const auto *dl_110 = buffer.data(dl + 110);
    const auto *dl_111 = buffer.data(dl + 111);
    const auto *dl_112 = buffer.data(dl + 112);
    const auto *dl_113 = buffer.data(dl + 113);
    const auto *dl_114 = buffer.data(dl + 114);
    const auto *dl_115 = buffer.data(dl + 115);
    const auto *dl_116 = buffer.data(dl + 116);
    const auto *dl_117 = buffer.data(dl + 117);
    const auto *dl_118 = buffer.data(dl + 118);
    const auto *dl_119 = buffer.data(dl + 119);
    const auto *dl_120 = buffer.data(dl + 120);
    const auto *dl_121 = buffer.data(dl + 121);
    const auto *dl_122 = buffer.data(dl + 122);
    const auto *dl_123 = buffer.data(dl + 123);
    const auto *dl_124 = buffer.data(dl + 124);
    const auto *dl_125 = buffer.data(dl + 125);
    const auto *dl_126 = buffer.data(dl + 126);
    const auto *dl_127 = buffer.data(dl + 127);
    const auto *dl_128 = buffer.data(dl + 128);
    const auto *dl_129 = buffer.data(dl + 129);
    const auto *dl_130 = buffer.data(dl + 130);
    const auto *dl_131 = buffer.data(dl + 131);
    const auto *dl_132 = buffer.data(dl + 132);
    const auto *dl_133 = buffer.data(dl + 133);
    const auto *dl_134 = buffer.data(dl + 134);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, sl_0, sl_1, sl_2, sl_3, sl_4, dl_0, dl_1, \
                         dl_2, dl_3, dl_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -sl_0[k]
                 + f_0 * dl_0[k];

        t_1[k] = -sl_1[k]
                 + f_0 * dl_1[k];

        t_2[k] = -sl_2[k]
                 + f_0 * dl_2[k];

        t_3[k] = -sl_3[k]
                 + f_0 * dl_3[k];

        t_4[k] = -sl_4[k]
                 + f_0 * dl_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, sl_5, sl_6, sl_7, sl_8, sl_9, dl_5, dl_6, \
                         dl_7, dl_8, dl_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -sl_5[k]
                 + f_0 * dl_5[k];

        t_6[k] = -sl_6[k]
                 + f_0 * dl_6[k];

        t_7[k] = -sl_7[k]
                 + f_0 * dl_7[k];

        t_8[k] = -sl_8[k]
                 + f_0 * dl_8[k];

        t_9[k] = -sl_9[k]
                 + f_0 * dl_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, sl_10, sl_11, sl_12, sl_13, sl_14, \
                         dl_10, dl_11, dl_12, dl_13, dl_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -sl_10[k]
                  + f_0 * dl_10[k];

        t_11[k] = -sl_11[k]
                  + f_0 * dl_11[k];

        t_12[k] = -sl_12[k]
                  + f_0 * dl_12[k];

        t_13[k] = -sl_13[k]
                  + f_0 * dl_13[k];

        t_14[k] = -sl_14[k]
                  + f_0 * dl_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, sl_15, sl_16, sl_17, sl_18, sl_19, \
                         dl_15, dl_16, dl_17, dl_18, dl_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -sl_15[k]
                  + f_0 * dl_15[k];

        t_16[k] = -sl_16[k]
                  + f_0 * dl_16[k];

        t_17[k] = -sl_17[k]
                  + f_0 * dl_17[k];

        t_18[k] = -sl_18[k]
                  + f_0 * dl_18[k];

        t_19[k] = -sl_19[k]
                  + f_0 * dl_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, sl_20, sl_21, sl_22, sl_23, sl_24, \
                         dl_20, dl_21, dl_22, dl_23, dl_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -sl_20[k]
                  + f_0 * dl_20[k];

        t_21[k] = -sl_21[k]
                  + f_0 * dl_21[k];

        t_22[k] = -sl_22[k]
                  + f_0 * dl_22[k];

        t_23[k] = -sl_23[k]
                  + f_0 * dl_23[k];

        t_24[k] = -sl_24[k]
                  + f_0 * dl_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, sl_25, sl_26, sl_27, sl_28, sl_29, \
                         dl_25, dl_26, dl_27, dl_28, dl_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -sl_25[k]
                  + f_0 * dl_25[k];

        t_26[k] = -sl_26[k]
                  + f_0 * dl_26[k];

        t_27[k] = -sl_27[k]
                  + f_0 * dl_27[k];

        t_28[k] = -sl_28[k]
                  + f_0 * dl_28[k];

        t_29[k] = -sl_29[k]
                  + f_0 * dl_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, sl_30, sl_31, sl_32, sl_33, sl_34, \
                         dl_30, dl_31, dl_32, dl_33, dl_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -sl_30[k]
                  + f_0 * dl_30[k];

        t_31[k] = -sl_31[k]
                  + f_0 * dl_31[k];

        t_32[k] = -sl_32[k]
                  + f_0 * dl_32[k];

        t_33[k] = -sl_33[k]
                  + f_0 * dl_33[k];

        t_34[k] = -sl_34[k]
                  + f_0 * dl_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, sl_35, sl_36, sl_37, sl_38, sl_39, \
                         dl_35, dl_36, dl_37, dl_38, dl_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -sl_35[k]
                  + f_0 * dl_35[k];

        t_36[k] = -sl_36[k]
                  + f_0 * dl_36[k];

        t_37[k] = -sl_37[k]
                  + f_0 * dl_37[k];

        t_38[k] = -sl_38[k]
                  + f_0 * dl_38[k];

        t_39[k] = -sl_39[k]
                  + f_0 * dl_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, sl_40, sl_41, sl_42, sl_43, sl_44, \
                         dl_40, dl_41, dl_42, dl_43, dl_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -sl_40[k]
                  + f_0 * dl_40[k];

        t_41[k] = -sl_41[k]
                  + f_0 * dl_41[k];

        t_42[k] = -sl_42[k]
                  + f_0 * dl_42[k];

        t_43[k] = -sl_43[k]
                  + f_0 * dl_43[k];

        t_44[k] = -sl_44[k]
                  + f_0 * dl_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, t_51, t_52, dl_45, dl_46, dl_47, \
                         dl_48, dl_49, dl_50, dl_51, dl_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_0 * dl_45[k];

        t_46[k] = f_0 * dl_46[k];

        t_47[k] = f_0 * dl_47[k];

        t_48[k] = f_0 * dl_48[k];

        t_49[k] = f_0 * dl_49[k];

        t_50[k] = f_0 * dl_50[k];

        t_51[k] = f_0 * dl_51[k];

        t_52[k] = f_0 * dl_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, t_59, t_60, dl_53, dl_54, dl_55, \
                         dl_56, dl_57, dl_58, dl_59, dl_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * dl_53[k];

        t_54[k] = f_0 * dl_54[k];

        t_55[k] = f_0 * dl_55[k];

        t_56[k] = f_0 * dl_56[k];

        t_57[k] = f_0 * dl_57[k];

        t_58[k] = f_0 * dl_58[k];

        t_59[k] = f_0 * dl_59[k];

        t_60[k] = f_0 * dl_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, t_67, t_68, dl_61, dl_62, dl_63, \
                         dl_64, dl_65, dl_66, dl_67, dl_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * dl_61[k];

        t_62[k] = f_0 * dl_62[k];

        t_63[k] = f_0 * dl_63[k];

        t_64[k] = f_0 * dl_64[k];

        t_65[k] = f_0 * dl_65[k];

        t_66[k] = f_0 * dl_66[k];

        t_67[k] = f_0 * dl_67[k];

        t_68[k] = f_0 * dl_68[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, t_74, t_75, t_76, dl_69, dl_70, dl_71, \
                         dl_72, dl_73, dl_74, dl_75, dl_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_0 * dl_69[k];

        t_70[k] = f_0 * dl_70[k];

        t_71[k] = f_0 * dl_71[k];

        t_72[k] = f_0 * dl_72[k];

        t_73[k] = f_0 * dl_73[k];

        t_74[k] = f_0 * dl_74[k];

        t_75[k] = f_0 * dl_75[k];

        t_76[k] = f_0 * dl_76[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, t_83, t_84, dl_77, dl_78, dl_79, \
                         dl_80, dl_81, dl_82, dl_83, dl_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_0 * dl_77[k];

        t_78[k] = f_0 * dl_78[k];

        t_79[k] = f_0 * dl_79[k];

        t_80[k] = f_0 * dl_80[k];

        t_81[k] = f_0 * dl_81[k];

        t_82[k] = f_0 * dl_82[k];

        t_83[k] = f_0 * dl_83[k];

        t_84[k] = f_0 * dl_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, t_90, t_91, t_92, dl_85, dl_86, dl_87, \
                         dl_88, dl_89, dl_90, dl_91, dl_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_0 * dl_85[k];

        t_86[k] = f_0 * dl_86[k];

        t_87[k] = f_0 * dl_87[k];

        t_88[k] = f_0 * dl_88[k];

        t_89[k] = f_0 * dl_89[k];

        t_90[k] = f_0 * dl_90[k];

        t_91[k] = f_0 * dl_91[k];

        t_92[k] = f_0 * dl_92[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, t_100, dl_93, dl_94, dl_95, \
                         dl_96, dl_97, dl_98, dl_99, dl_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * dl_93[k];

        t_94[k] = f_0 * dl_94[k];

        t_95[k] = f_0 * dl_95[k];

        t_96[k] = f_0 * dl_96[k];

        t_97[k] = f_0 * dl_97[k];

        t_98[k] = f_0 * dl_98[k];

        t_99[k] = f_0 * dl_99[k];

        t_100[k] = f_0 * dl_100[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, t_107, t_108, dl_101, \
                         dl_102, dl_103, dl_104, dl_105, dl_106, dl_107, \
                         dl_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * dl_101[k];

        t_102[k] = f_0 * dl_102[k];

        t_103[k] = f_0 * dl_103[k];

        t_104[k] = f_0 * dl_104[k];

        t_105[k] = f_0 * dl_105[k];

        t_106[k] = f_0 * dl_106[k];

        t_107[k] = f_0 * dl_107[k];

        t_108[k] = f_0 * dl_108[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, t_114, t_115, t_116, dl_109, \
                         dl_110, dl_111, dl_112, dl_113, dl_114, dl_115, \
                         dl_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_0 * dl_109[k];

        t_110[k] = f_0 * dl_110[k];

        t_111[k] = f_0 * dl_111[k];

        t_112[k] = f_0 * dl_112[k];

        t_113[k] = f_0 * dl_113[k];

        t_114[k] = f_0 * dl_114[k];

        t_115[k] = f_0 * dl_115[k];

        t_116[k] = f_0 * dl_116[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, t_122, t_123, t_124, dl_117, \
                         dl_118, dl_119, dl_120, dl_121, dl_122, dl_123, \
                         dl_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_0 * dl_117[k];

        t_118[k] = f_0 * dl_118[k];

        t_119[k] = f_0 * dl_119[k];

        t_120[k] = f_0 * dl_120[k];

        t_121[k] = f_0 * dl_121[k];

        t_122[k] = f_0 * dl_122[k];

        t_123[k] = f_0 * dl_123[k];

        t_124[k] = f_0 * dl_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, t_130, t_131, t_132, dl_125, \
                         dl_126, dl_127, dl_128, dl_129, dl_130, dl_131, \
                         dl_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * dl_125[k];

        t_126[k] = f_0 * dl_126[k];

        t_127[k] = f_0 * dl_127[k];

        t_128[k] = f_0 * dl_128[k];

        t_129[k] = f_0 * dl_129[k];

        t_130[k] = f_0 * dl_130[k];

        t_131[k] = f_0 * dl_131[k];

        t_132[k] = f_0 * dl_132[k];
    }

#pragma omp simd aligned(t_133, t_134, dl_133, dl_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_0 * dl_133[k];

        t_134[k] = f_0 * dl_134[k];
    }
}

auto
compute_prim_geom_10_pl_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t sl, const size_t dl,
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
    auto *t_108 = buffer.data(target + 108);
    auto *t_109 = buffer.data(target + 109);
    auto *t_110 = buffer.data(target + 110);
    auto *t_111 = buffer.data(target + 111);
    auto *t_112 = buffer.data(target + 112);
    auto *t_113 = buffer.data(target + 113);
    auto *t_114 = buffer.data(target + 114);
    auto *t_115 = buffer.data(target + 115);
    auto *t_116 = buffer.data(target + 116);
    auto *t_117 = buffer.data(target + 117);
    auto *t_118 = buffer.data(target + 118);
    auto *t_119 = buffer.data(target + 119);
    auto *t_120 = buffer.data(target + 120);
    auto *t_121 = buffer.data(target + 121);
    auto *t_122 = buffer.data(target + 122);
    auto *t_123 = buffer.data(target + 123);
    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);

    const auto *sl_0 = buffer.data(sl + 0);
    const auto *sl_1 = buffer.data(sl + 1);
    const auto *sl_2 = buffer.data(sl + 2);
    const auto *sl_3 = buffer.data(sl + 3);
    const auto *sl_4 = buffer.data(sl + 4);
    const auto *sl_5 = buffer.data(sl + 5);
    const auto *sl_6 = buffer.data(sl + 6);
    const auto *sl_7 = buffer.data(sl + 7);
    const auto *sl_8 = buffer.data(sl + 8);
    const auto *sl_9 = buffer.data(sl + 9);
    const auto *sl_10 = buffer.data(sl + 10);
    const auto *sl_11 = buffer.data(sl + 11);
    const auto *sl_12 = buffer.data(sl + 12);
    const auto *sl_13 = buffer.data(sl + 13);
    const auto *sl_14 = buffer.data(sl + 14);
    const auto *sl_15 = buffer.data(sl + 15);
    const auto *sl_16 = buffer.data(sl + 16);
    const auto *sl_17 = buffer.data(sl + 17);
    const auto *sl_18 = buffer.data(sl + 18);
    const auto *sl_19 = buffer.data(sl + 19);
    const auto *sl_20 = buffer.data(sl + 20);
    const auto *sl_21 = buffer.data(sl + 21);
    const auto *sl_22 = buffer.data(sl + 22);
    const auto *sl_23 = buffer.data(sl + 23);
    const auto *sl_24 = buffer.data(sl + 24);
    const auto *sl_25 = buffer.data(sl + 25);
    const auto *sl_26 = buffer.data(sl + 26);
    const auto *sl_27 = buffer.data(sl + 27);
    const auto *sl_28 = buffer.data(sl + 28);
    const auto *sl_29 = buffer.data(sl + 29);
    const auto *sl_30 = buffer.data(sl + 30);
    const auto *sl_31 = buffer.data(sl + 31);
    const auto *sl_32 = buffer.data(sl + 32);
    const auto *sl_33 = buffer.data(sl + 33);
    const auto *sl_34 = buffer.data(sl + 34);
    const auto *sl_35 = buffer.data(sl + 35);
    const auto *sl_36 = buffer.data(sl + 36);
    const auto *sl_37 = buffer.data(sl + 37);
    const auto *sl_38 = buffer.data(sl + 38);
    const auto *sl_39 = buffer.data(sl + 39);
    const auto *sl_40 = buffer.data(sl + 40);
    const auto *sl_41 = buffer.data(sl + 41);
    const auto *sl_42 = buffer.data(sl + 42);
    const auto *sl_43 = buffer.data(sl + 43);
    const auto *sl_44 = buffer.data(sl + 44);

    const auto *dl_45 = buffer.data(dl + 45);
    const auto *dl_46 = buffer.data(dl + 46);
    const auto *dl_47 = buffer.data(dl + 47);
    const auto *dl_48 = buffer.data(dl + 48);
    const auto *dl_49 = buffer.data(dl + 49);
    const auto *dl_50 = buffer.data(dl + 50);
    const auto *dl_51 = buffer.data(dl + 51);
    const auto *dl_52 = buffer.data(dl + 52);
    const auto *dl_53 = buffer.data(dl + 53);
    const auto *dl_54 = buffer.data(dl + 54);
    const auto *dl_55 = buffer.data(dl + 55);
    const auto *dl_56 = buffer.data(dl + 56);
    const auto *dl_57 = buffer.data(dl + 57);
    const auto *dl_58 = buffer.data(dl + 58);
    const auto *dl_59 = buffer.data(dl + 59);
    const auto *dl_60 = buffer.data(dl + 60);
    const auto *dl_61 = buffer.data(dl + 61);
    const auto *dl_62 = buffer.data(dl + 62);
    const auto *dl_63 = buffer.data(dl + 63);
    const auto *dl_64 = buffer.data(dl + 64);
    const auto *dl_65 = buffer.data(dl + 65);
    const auto *dl_66 = buffer.data(dl + 66);
    const auto *dl_67 = buffer.data(dl + 67);
    const auto *dl_68 = buffer.data(dl + 68);
    const auto *dl_69 = buffer.data(dl + 69);
    const auto *dl_70 = buffer.data(dl + 70);
    const auto *dl_71 = buffer.data(dl + 71);
    const auto *dl_72 = buffer.data(dl + 72);
    const auto *dl_73 = buffer.data(dl + 73);
    const auto *dl_74 = buffer.data(dl + 74);
    const auto *dl_75 = buffer.data(dl + 75);
    const auto *dl_76 = buffer.data(dl + 76);
    const auto *dl_77 = buffer.data(dl + 77);
    const auto *dl_78 = buffer.data(dl + 78);
    const auto *dl_79 = buffer.data(dl + 79);
    const auto *dl_80 = buffer.data(dl + 80);
    const auto *dl_81 = buffer.data(dl + 81);
    const auto *dl_82 = buffer.data(dl + 82);
    const auto *dl_83 = buffer.data(dl + 83);
    const auto *dl_84 = buffer.data(dl + 84);
    const auto *dl_85 = buffer.data(dl + 85);
    const auto *dl_86 = buffer.data(dl + 86);
    const auto *dl_87 = buffer.data(dl + 87);
    const auto *dl_88 = buffer.data(dl + 88);
    const auto *dl_89 = buffer.data(dl + 89);
    const auto *dl_135 = buffer.data(dl + 135);
    const auto *dl_136 = buffer.data(dl + 136);
    const auto *dl_137 = buffer.data(dl + 137);
    const auto *dl_138 = buffer.data(dl + 138);
    const auto *dl_139 = buffer.data(dl + 139);
    const auto *dl_140 = buffer.data(dl + 140);
    const auto *dl_141 = buffer.data(dl + 141);
    const auto *dl_142 = buffer.data(dl + 142);
    const auto *dl_143 = buffer.data(dl + 143);
    const auto *dl_144 = buffer.data(dl + 144);
    const auto *dl_145 = buffer.data(dl + 145);
    const auto *dl_146 = buffer.data(dl + 146);
    const auto *dl_147 = buffer.data(dl + 147);
    const auto *dl_148 = buffer.data(dl + 148);
    const auto *dl_149 = buffer.data(dl + 149);
    const auto *dl_150 = buffer.data(dl + 150);
    const auto *dl_151 = buffer.data(dl + 151);
    const auto *dl_152 = buffer.data(dl + 152);
    const auto *dl_153 = buffer.data(dl + 153);
    const auto *dl_154 = buffer.data(dl + 154);
    const auto *dl_155 = buffer.data(dl + 155);
    const auto *dl_156 = buffer.data(dl + 156);
    const auto *dl_157 = buffer.data(dl + 157);
    const auto *dl_158 = buffer.data(dl + 158);
    const auto *dl_159 = buffer.data(dl + 159);
    const auto *dl_160 = buffer.data(dl + 160);
    const auto *dl_161 = buffer.data(dl + 161);
    const auto *dl_162 = buffer.data(dl + 162);
    const auto *dl_163 = buffer.data(dl + 163);
    const auto *dl_164 = buffer.data(dl + 164);
    const auto *dl_165 = buffer.data(dl + 165);
    const auto *dl_166 = buffer.data(dl + 166);
    const auto *dl_167 = buffer.data(dl + 167);
    const auto *dl_168 = buffer.data(dl + 168);
    const auto *dl_169 = buffer.data(dl + 169);
    const auto *dl_170 = buffer.data(dl + 170);
    const auto *dl_171 = buffer.data(dl + 171);
    const auto *dl_172 = buffer.data(dl + 172);
    const auto *dl_173 = buffer.data(dl + 173);
    const auto *dl_174 = buffer.data(dl + 174);
    const auto *dl_175 = buffer.data(dl + 175);
    const auto *dl_176 = buffer.data(dl + 176);
    const auto *dl_177 = buffer.data(dl + 177);
    const auto *dl_178 = buffer.data(dl + 178);
    const auto *dl_179 = buffer.data(dl + 179);
    const auto *dl_180 = buffer.data(dl + 180);
    const auto *dl_181 = buffer.data(dl + 181);
    const auto *dl_182 = buffer.data(dl + 182);
    const auto *dl_183 = buffer.data(dl + 183);
    const auto *dl_184 = buffer.data(dl + 184);
    const auto *dl_185 = buffer.data(dl + 185);
    const auto *dl_186 = buffer.data(dl + 186);
    const auto *dl_187 = buffer.data(dl + 187);
    const auto *dl_188 = buffer.data(dl + 188);
    const auto *dl_189 = buffer.data(dl + 189);
    const auto *dl_190 = buffer.data(dl + 190);
    const auto *dl_191 = buffer.data(dl + 191);
    const auto *dl_192 = buffer.data(dl + 192);
    const auto *dl_193 = buffer.data(dl + 193);
    const auto *dl_194 = buffer.data(dl + 194);
    const auto *dl_195 = buffer.data(dl + 195);
    const auto *dl_196 = buffer.data(dl + 196);
    const auto *dl_197 = buffer.data(dl + 197);
    const auto *dl_198 = buffer.data(dl + 198);
    const auto *dl_199 = buffer.data(dl + 199);
    const auto *dl_200 = buffer.data(dl + 200);
    const auto *dl_201 = buffer.data(dl + 201);
    const auto *dl_202 = buffer.data(dl + 202);
    const auto *dl_203 = buffer.data(dl + 203);
    const auto *dl_204 = buffer.data(dl + 204);
    const auto *dl_205 = buffer.data(dl + 205);
    const auto *dl_206 = buffer.data(dl + 206);
    const auto *dl_207 = buffer.data(dl + 207);
    const auto *dl_208 = buffer.data(dl + 208);
    const auto *dl_209 = buffer.data(dl + 209);
    const auto *dl_210 = buffer.data(dl + 210);
    const auto *dl_211 = buffer.data(dl + 211);
    const auto *dl_212 = buffer.data(dl + 212);
    const auto *dl_213 = buffer.data(dl + 213);
    const auto *dl_214 = buffer.data(dl + 214);
    const auto *dl_215 = buffer.data(dl + 215);
    const auto *dl_216 = buffer.data(dl + 216);
    const auto *dl_217 = buffer.data(dl + 217);
    const auto *dl_218 = buffer.data(dl + 218);
    const auto *dl_219 = buffer.data(dl + 219);
    const auto *dl_220 = buffer.data(dl + 220);
    const auto *dl_221 = buffer.data(dl + 221);
    const auto *dl_222 = buffer.data(dl + 222);
    const auto *dl_223 = buffer.data(dl + 223);
    const auto *dl_224 = buffer.data(dl + 224);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, dl_45, dl_46, dl_47, dl_48, \
                         dl_49, dl_50, dl_51, dl_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dl_45[k];

        t_1[k] = f_0 * dl_46[k];

        t_2[k] = f_0 * dl_47[k];

        t_3[k] = f_0 * dl_48[k];

        t_4[k] = f_0 * dl_49[k];

        t_5[k] = f_0 * dl_50[k];

        t_6[k] = f_0 * dl_51[k];

        t_7[k] = f_0 * dl_52[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, dl_53, dl_54, dl_55, \
                         dl_56, dl_57, dl_58, dl_59, dl_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dl_53[k];

        t_9[k] = f_0 * dl_54[k];

        t_10[k] = f_0 * dl_55[k];

        t_11[k] = f_0 * dl_56[k];

        t_12[k] = f_0 * dl_57[k];

        t_13[k] = f_0 * dl_58[k];

        t_14[k] = f_0 * dl_59[k];

        t_15[k] = f_0 * dl_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, dl_61, dl_62, dl_63, \
                         dl_64, dl_65, dl_66, dl_67, dl_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * dl_61[k];

        t_17[k] = f_0 * dl_62[k];

        t_18[k] = f_0 * dl_63[k];

        t_19[k] = f_0 * dl_64[k];

        t_20[k] = f_0 * dl_65[k];

        t_21[k] = f_0 * dl_66[k];

        t_22[k] = f_0 * dl_67[k];

        t_23[k] = f_0 * dl_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, dl_69, dl_70, dl_71, \
                         dl_72, dl_73, dl_74, dl_75, dl_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * dl_69[k];

        t_25[k] = f_0 * dl_70[k];

        t_26[k] = f_0 * dl_71[k];

        t_27[k] = f_0 * dl_72[k];

        t_28[k] = f_0 * dl_73[k];

        t_29[k] = f_0 * dl_74[k];

        t_30[k] = f_0 * dl_75[k];

        t_31[k] = f_0 * dl_76[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, dl_77, dl_78, dl_79, \
                         dl_80, dl_81, dl_82, dl_83, dl_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * dl_77[k];

        t_33[k] = f_0 * dl_78[k];

        t_34[k] = f_0 * dl_79[k];

        t_35[k] = f_0 * dl_80[k];

        t_36[k] = f_0 * dl_81[k];

        t_37[k] = f_0 * dl_82[k];

        t_38[k] = f_0 * dl_83[k];

        t_39[k] = f_0 * dl_84[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, sl_0, sl_1, dl_85, dl_86, \
                         dl_87, dl_88, dl_89, dl_135, dl_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * dl_85[k];

        t_41[k] = f_0 * dl_86[k];

        t_42[k] = f_0 * dl_87[k];

        t_43[k] = f_0 * dl_88[k];

        t_44[k] = f_0 * dl_89[k];

        t_45[k] = -sl_0[k]
                  + f_0 * dl_135[k];

        t_46[k] = -sl_1[k]
                  + f_0 * dl_136[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, sl_2, sl_3, sl_4, sl_5, sl_6, dl_137, \
                         dl_138, dl_139, dl_140, dl_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -sl_2[k]
                  + f_0 * dl_137[k];

        t_48[k] = -sl_3[k]
                  + f_0 * dl_138[k];

        t_49[k] = -sl_4[k]
                  + f_0 * dl_139[k];

        t_50[k] = -sl_5[k]
                  + f_0 * dl_140[k];

        t_51[k] = -sl_6[k]
                  + f_0 * dl_141[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, sl_7, sl_8, sl_9, sl_10, sl_11, dl_142, \
                         dl_143, dl_144, dl_145, dl_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = -sl_7[k]
                  + f_0 * dl_142[k];

        t_53[k] = -sl_8[k]
                  + f_0 * dl_143[k];

        t_54[k] = -sl_9[k]
                  + f_0 * dl_144[k];

        t_55[k] = -sl_10[k]
                  + f_0 * dl_145[k];

        t_56[k] = -sl_11[k]
                  + f_0 * dl_146[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, sl_12, sl_13, sl_14, sl_15, sl_16, \
                         dl_147, dl_148, dl_149, dl_150, dl_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = -sl_12[k]
                  + f_0 * dl_147[k];

        t_58[k] = -sl_13[k]
                  + f_0 * dl_148[k];

        t_59[k] = -sl_14[k]
                  + f_0 * dl_149[k];

        t_60[k] = -sl_15[k]
                  + f_0 * dl_150[k];

        t_61[k] = -sl_16[k]
                  + f_0 * dl_151[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, sl_17, sl_18, sl_19, sl_20, sl_21, \
                         dl_152, dl_153, dl_154, dl_155, dl_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -sl_17[k]
                  + f_0 * dl_152[k];

        t_63[k] = -sl_18[k]
                  + f_0 * dl_153[k];

        t_64[k] = -sl_19[k]
                  + f_0 * dl_154[k];

        t_65[k] = -sl_20[k]
                  + f_0 * dl_155[k];

        t_66[k] = -sl_21[k]
                  + f_0 * dl_156[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, sl_22, sl_23, sl_24, sl_25, sl_26, \
                         dl_157, dl_158, dl_159, dl_160, dl_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -sl_22[k]
                  + f_0 * dl_157[k];

        t_68[k] = -sl_23[k]
                  + f_0 * dl_158[k];

        t_69[k] = -sl_24[k]
                  + f_0 * dl_159[k];

        t_70[k] = -sl_25[k]
                  + f_0 * dl_160[k];

        t_71[k] = -sl_26[k]
                  + f_0 * dl_161[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, sl_27, sl_28, sl_29, sl_30, sl_31, \
                         dl_162, dl_163, dl_164, dl_165, dl_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -sl_27[k]
                  + f_0 * dl_162[k];

        t_73[k] = -sl_28[k]
                  + f_0 * dl_163[k];

        t_74[k] = -sl_29[k]
                  + f_0 * dl_164[k];

        t_75[k] = -sl_30[k]
                  + f_0 * dl_165[k];

        t_76[k] = -sl_31[k]
                  + f_0 * dl_166[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, sl_32, sl_33, sl_34, sl_35, sl_36, \
                         dl_167, dl_168, dl_169, dl_170, dl_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -sl_32[k]
                  + f_0 * dl_167[k];

        t_78[k] = -sl_33[k]
                  + f_0 * dl_168[k];

        t_79[k] = -sl_34[k]
                  + f_0 * dl_169[k];

        t_80[k] = -sl_35[k]
                  + f_0 * dl_170[k];

        t_81[k] = -sl_36[k]
                  + f_0 * dl_171[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, sl_37, sl_38, sl_39, sl_40, sl_41, \
                         dl_172, dl_173, dl_174, dl_175, dl_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -sl_37[k]
                  + f_0 * dl_172[k];

        t_83[k] = -sl_38[k]
                  + f_0 * dl_173[k];

        t_84[k] = -sl_39[k]
                  + f_0 * dl_174[k];

        t_85[k] = -sl_40[k]
                  + f_0 * dl_175[k];

        t_86[k] = -sl_41[k]
                  + f_0 * dl_176[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, sl_42, sl_43, sl_44, dl_177, \
                         dl_178, dl_179, dl_180, dl_181, dl_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -sl_42[k]
                  + f_0 * dl_177[k];

        t_88[k] = -sl_43[k]
                  + f_0 * dl_178[k];

        t_89[k] = -sl_44[k]
                  + f_0 * dl_179[k];

        t_90[k] = f_0 * dl_180[k];

        t_91[k] = f_0 * dl_181[k];

        t_92[k] = f_0 * dl_182[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, t_100, dl_183, dl_184, \
                         dl_185, dl_186, dl_187, dl_188, dl_189, \
                         dl_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * dl_183[k];

        t_94[k] = f_0 * dl_184[k];

        t_95[k] = f_0 * dl_185[k];

        t_96[k] = f_0 * dl_186[k];

        t_97[k] = f_0 * dl_187[k];

        t_98[k] = f_0 * dl_188[k];

        t_99[k] = f_0 * dl_189[k];

        t_100[k] = f_0 * dl_190[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, t_107, t_108, dl_191, \
                         dl_192, dl_193, dl_194, dl_195, dl_196, dl_197, \
                         dl_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * dl_191[k];

        t_102[k] = f_0 * dl_192[k];

        t_103[k] = f_0 * dl_193[k];

        t_104[k] = f_0 * dl_194[k];

        t_105[k] = f_0 * dl_195[k];

        t_106[k] = f_0 * dl_196[k];

        t_107[k] = f_0 * dl_197[k];

        t_108[k] = f_0 * dl_198[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, t_114, t_115, t_116, dl_199, \
                         dl_200, dl_201, dl_202, dl_203, dl_204, dl_205, \
                         dl_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_0 * dl_199[k];

        t_110[k] = f_0 * dl_200[k];

        t_111[k] = f_0 * dl_201[k];

        t_112[k] = f_0 * dl_202[k];

        t_113[k] = f_0 * dl_203[k];

        t_114[k] = f_0 * dl_204[k];

        t_115[k] = f_0 * dl_205[k];

        t_116[k] = f_0 * dl_206[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, t_122, t_123, t_124, dl_207, \
                         dl_208, dl_209, dl_210, dl_211, dl_212, dl_213, \
                         dl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_0 * dl_207[k];

        t_118[k] = f_0 * dl_208[k];

        t_119[k] = f_0 * dl_209[k];

        t_120[k] = f_0 * dl_210[k];

        t_121[k] = f_0 * dl_211[k];

        t_122[k] = f_0 * dl_212[k];

        t_123[k] = f_0 * dl_213[k];

        t_124[k] = f_0 * dl_214[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, t_130, t_131, t_132, dl_215, \
                         dl_216, dl_217, dl_218, dl_219, dl_220, dl_221, \
                         dl_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * dl_215[k];

        t_126[k] = f_0 * dl_216[k];

        t_127[k] = f_0 * dl_217[k];

        t_128[k] = f_0 * dl_218[k];

        t_129[k] = f_0 * dl_219[k];

        t_130[k] = f_0 * dl_220[k];

        t_131[k] = f_0 * dl_221[k];

        t_132[k] = f_0 * dl_222[k];
    }

#pragma omp simd aligned(t_133, t_134, dl_223, dl_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_0 * dl_223[k];

        t_134[k] = f_0 * dl_224[k];
    }
}

auto
compute_prim_geom_10_pl_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t sl, const size_t dl,
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
    auto *t_108 = buffer.data(target + 108);
    auto *t_109 = buffer.data(target + 109);
    auto *t_110 = buffer.data(target + 110);
    auto *t_111 = buffer.data(target + 111);
    auto *t_112 = buffer.data(target + 112);
    auto *t_113 = buffer.data(target + 113);
    auto *t_114 = buffer.data(target + 114);
    auto *t_115 = buffer.data(target + 115);
    auto *t_116 = buffer.data(target + 116);
    auto *t_117 = buffer.data(target + 117);
    auto *t_118 = buffer.data(target + 118);
    auto *t_119 = buffer.data(target + 119);
    auto *t_120 = buffer.data(target + 120);
    auto *t_121 = buffer.data(target + 121);
    auto *t_122 = buffer.data(target + 122);
    auto *t_123 = buffer.data(target + 123);
    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);

    const auto *sl_0 = buffer.data(sl + 0);
    const auto *sl_1 = buffer.data(sl + 1);
    const auto *sl_2 = buffer.data(sl + 2);
    const auto *sl_3 = buffer.data(sl + 3);
    const auto *sl_4 = buffer.data(sl + 4);
    const auto *sl_5 = buffer.data(sl + 5);
    const auto *sl_6 = buffer.data(sl + 6);
    const auto *sl_7 = buffer.data(sl + 7);
    const auto *sl_8 = buffer.data(sl + 8);
    const auto *sl_9 = buffer.data(sl + 9);
    const auto *sl_10 = buffer.data(sl + 10);
    const auto *sl_11 = buffer.data(sl + 11);
    const auto *sl_12 = buffer.data(sl + 12);
    const auto *sl_13 = buffer.data(sl + 13);
    const auto *sl_14 = buffer.data(sl + 14);
    const auto *sl_15 = buffer.data(sl + 15);
    const auto *sl_16 = buffer.data(sl + 16);
    const auto *sl_17 = buffer.data(sl + 17);
    const auto *sl_18 = buffer.data(sl + 18);
    const auto *sl_19 = buffer.data(sl + 19);
    const auto *sl_20 = buffer.data(sl + 20);
    const auto *sl_21 = buffer.data(sl + 21);
    const auto *sl_22 = buffer.data(sl + 22);
    const auto *sl_23 = buffer.data(sl + 23);
    const auto *sl_24 = buffer.data(sl + 24);
    const auto *sl_25 = buffer.data(sl + 25);
    const auto *sl_26 = buffer.data(sl + 26);
    const auto *sl_27 = buffer.data(sl + 27);
    const auto *sl_28 = buffer.data(sl + 28);
    const auto *sl_29 = buffer.data(sl + 29);
    const auto *sl_30 = buffer.data(sl + 30);
    const auto *sl_31 = buffer.data(sl + 31);
    const auto *sl_32 = buffer.data(sl + 32);
    const auto *sl_33 = buffer.data(sl + 33);
    const auto *sl_34 = buffer.data(sl + 34);
    const auto *sl_35 = buffer.data(sl + 35);
    const auto *sl_36 = buffer.data(sl + 36);
    const auto *sl_37 = buffer.data(sl + 37);
    const auto *sl_38 = buffer.data(sl + 38);
    const auto *sl_39 = buffer.data(sl + 39);
    const auto *sl_40 = buffer.data(sl + 40);
    const auto *sl_41 = buffer.data(sl + 41);
    const auto *sl_42 = buffer.data(sl + 42);
    const auto *sl_43 = buffer.data(sl + 43);
    const auto *sl_44 = buffer.data(sl + 44);

    const auto *dl_90 = buffer.data(dl + 90);
    const auto *dl_91 = buffer.data(dl + 91);
    const auto *dl_92 = buffer.data(dl + 92);
    const auto *dl_93 = buffer.data(dl + 93);
    const auto *dl_94 = buffer.data(dl + 94);
    const auto *dl_95 = buffer.data(dl + 95);
    const auto *dl_96 = buffer.data(dl + 96);
    const auto *dl_97 = buffer.data(dl + 97);
    const auto *dl_98 = buffer.data(dl + 98);
    const auto *dl_99 = buffer.data(dl + 99);
    const auto *dl_100 = buffer.data(dl + 100);
    const auto *dl_101 = buffer.data(dl + 101);
    const auto *dl_102 = buffer.data(dl + 102);
    const auto *dl_103 = buffer.data(dl + 103);
    const auto *dl_104 = buffer.data(dl + 104);
    const auto *dl_105 = buffer.data(dl + 105);
    const auto *dl_106 = buffer.data(dl + 106);
    const auto *dl_107 = buffer.data(dl + 107);
    const auto *dl_108 = buffer.data(dl + 108);
    const auto *dl_109 = buffer.data(dl + 109);
    const auto *dl_110 = buffer.data(dl + 110);
    const auto *dl_111 = buffer.data(dl + 111);
    const auto *dl_112 = buffer.data(dl + 112);
    const auto *dl_113 = buffer.data(dl + 113);
    const auto *dl_114 = buffer.data(dl + 114);
    const auto *dl_115 = buffer.data(dl + 115);
    const auto *dl_116 = buffer.data(dl + 116);
    const auto *dl_117 = buffer.data(dl + 117);
    const auto *dl_118 = buffer.data(dl + 118);
    const auto *dl_119 = buffer.data(dl + 119);
    const auto *dl_120 = buffer.data(dl + 120);
    const auto *dl_121 = buffer.data(dl + 121);
    const auto *dl_122 = buffer.data(dl + 122);
    const auto *dl_123 = buffer.data(dl + 123);
    const auto *dl_124 = buffer.data(dl + 124);
    const auto *dl_125 = buffer.data(dl + 125);
    const auto *dl_126 = buffer.data(dl + 126);
    const auto *dl_127 = buffer.data(dl + 127);
    const auto *dl_128 = buffer.data(dl + 128);
    const auto *dl_129 = buffer.data(dl + 129);
    const auto *dl_130 = buffer.data(dl + 130);
    const auto *dl_131 = buffer.data(dl + 131);
    const auto *dl_132 = buffer.data(dl + 132);
    const auto *dl_133 = buffer.data(dl + 133);
    const auto *dl_134 = buffer.data(dl + 134);
    const auto *dl_180 = buffer.data(dl + 180);
    const auto *dl_181 = buffer.data(dl + 181);
    const auto *dl_182 = buffer.data(dl + 182);
    const auto *dl_183 = buffer.data(dl + 183);
    const auto *dl_184 = buffer.data(dl + 184);
    const auto *dl_185 = buffer.data(dl + 185);
    const auto *dl_186 = buffer.data(dl + 186);
    const auto *dl_187 = buffer.data(dl + 187);
    const auto *dl_188 = buffer.data(dl + 188);
    const auto *dl_189 = buffer.data(dl + 189);
    const auto *dl_190 = buffer.data(dl + 190);
    const auto *dl_191 = buffer.data(dl + 191);
    const auto *dl_192 = buffer.data(dl + 192);
    const auto *dl_193 = buffer.data(dl + 193);
    const auto *dl_194 = buffer.data(dl + 194);
    const auto *dl_195 = buffer.data(dl + 195);
    const auto *dl_196 = buffer.data(dl + 196);
    const auto *dl_197 = buffer.data(dl + 197);
    const auto *dl_198 = buffer.data(dl + 198);
    const auto *dl_199 = buffer.data(dl + 199);
    const auto *dl_200 = buffer.data(dl + 200);
    const auto *dl_201 = buffer.data(dl + 201);
    const auto *dl_202 = buffer.data(dl + 202);
    const auto *dl_203 = buffer.data(dl + 203);
    const auto *dl_204 = buffer.data(dl + 204);
    const auto *dl_205 = buffer.data(dl + 205);
    const auto *dl_206 = buffer.data(dl + 206);
    const auto *dl_207 = buffer.data(dl + 207);
    const auto *dl_208 = buffer.data(dl + 208);
    const auto *dl_209 = buffer.data(dl + 209);
    const auto *dl_210 = buffer.data(dl + 210);
    const auto *dl_211 = buffer.data(dl + 211);
    const auto *dl_212 = buffer.data(dl + 212);
    const auto *dl_213 = buffer.data(dl + 213);
    const auto *dl_214 = buffer.data(dl + 214);
    const auto *dl_215 = buffer.data(dl + 215);
    const auto *dl_216 = buffer.data(dl + 216);
    const auto *dl_217 = buffer.data(dl + 217);
    const auto *dl_218 = buffer.data(dl + 218);
    const auto *dl_219 = buffer.data(dl + 219);
    const auto *dl_220 = buffer.data(dl + 220);
    const auto *dl_221 = buffer.data(dl + 221);
    const auto *dl_222 = buffer.data(dl + 222);
    const auto *dl_223 = buffer.data(dl + 223);
    const auto *dl_224 = buffer.data(dl + 224);
    const auto *dl_225 = buffer.data(dl + 225);
    const auto *dl_226 = buffer.data(dl + 226);
    const auto *dl_227 = buffer.data(dl + 227);
    const auto *dl_228 = buffer.data(dl + 228);
    const auto *dl_229 = buffer.data(dl + 229);
    const auto *dl_230 = buffer.data(dl + 230);
    const auto *dl_231 = buffer.data(dl + 231);
    const auto *dl_232 = buffer.data(dl + 232);
    const auto *dl_233 = buffer.data(dl + 233);
    const auto *dl_234 = buffer.data(dl + 234);
    const auto *dl_235 = buffer.data(dl + 235);
    const auto *dl_236 = buffer.data(dl + 236);
    const auto *dl_237 = buffer.data(dl + 237);
    const auto *dl_238 = buffer.data(dl + 238);
    const auto *dl_239 = buffer.data(dl + 239);
    const auto *dl_240 = buffer.data(dl + 240);
    const auto *dl_241 = buffer.data(dl + 241);
    const auto *dl_242 = buffer.data(dl + 242);
    const auto *dl_243 = buffer.data(dl + 243);
    const auto *dl_244 = buffer.data(dl + 244);
    const auto *dl_245 = buffer.data(dl + 245);
    const auto *dl_246 = buffer.data(dl + 246);
    const auto *dl_247 = buffer.data(dl + 247);
    const auto *dl_248 = buffer.data(dl + 248);
    const auto *dl_249 = buffer.data(dl + 249);
    const auto *dl_250 = buffer.data(dl + 250);
    const auto *dl_251 = buffer.data(dl + 251);
    const auto *dl_252 = buffer.data(dl + 252);
    const auto *dl_253 = buffer.data(dl + 253);
    const auto *dl_254 = buffer.data(dl + 254);
    const auto *dl_255 = buffer.data(dl + 255);
    const auto *dl_256 = buffer.data(dl + 256);
    const auto *dl_257 = buffer.data(dl + 257);
    const auto *dl_258 = buffer.data(dl + 258);
    const auto *dl_259 = buffer.data(dl + 259);
    const auto *dl_260 = buffer.data(dl + 260);
    const auto *dl_261 = buffer.data(dl + 261);
    const auto *dl_262 = buffer.data(dl + 262);
    const auto *dl_263 = buffer.data(dl + 263);
    const auto *dl_264 = buffer.data(dl + 264);
    const auto *dl_265 = buffer.data(dl + 265);
    const auto *dl_266 = buffer.data(dl + 266);
    const auto *dl_267 = buffer.data(dl + 267);
    const auto *dl_268 = buffer.data(dl + 268);
    const auto *dl_269 = buffer.data(dl + 269);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, dl_90, dl_91, dl_92, dl_93, \
                         dl_94, dl_95, dl_96, dl_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dl_90[k];

        t_1[k] = f_0 * dl_91[k];

        t_2[k] = f_0 * dl_92[k];

        t_3[k] = f_0 * dl_93[k];

        t_4[k] = f_0 * dl_94[k];

        t_5[k] = f_0 * dl_95[k];

        t_6[k] = f_0 * dl_96[k];

        t_7[k] = f_0 * dl_97[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, dl_98, dl_99, dl_100, \
                         dl_101, dl_102, dl_103, dl_104, dl_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dl_98[k];

        t_9[k] = f_0 * dl_99[k];

        t_10[k] = f_0 * dl_100[k];

        t_11[k] = f_0 * dl_101[k];

        t_12[k] = f_0 * dl_102[k];

        t_13[k] = f_0 * dl_103[k];

        t_14[k] = f_0 * dl_104[k];

        t_15[k] = f_0 * dl_105[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, dl_106, dl_107, \
                         dl_108, dl_109, dl_110, dl_111, dl_112, \
                         dl_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * dl_106[k];

        t_17[k] = f_0 * dl_107[k];

        t_18[k] = f_0 * dl_108[k];

        t_19[k] = f_0 * dl_109[k];

        t_20[k] = f_0 * dl_110[k];

        t_21[k] = f_0 * dl_111[k];

        t_22[k] = f_0 * dl_112[k];

        t_23[k] = f_0 * dl_113[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, dl_114, dl_115, \
                         dl_116, dl_117, dl_118, dl_119, dl_120, \
                         dl_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * dl_114[k];

        t_25[k] = f_0 * dl_115[k];

        t_26[k] = f_0 * dl_116[k];

        t_27[k] = f_0 * dl_117[k];

        t_28[k] = f_0 * dl_118[k];

        t_29[k] = f_0 * dl_119[k];

        t_30[k] = f_0 * dl_120[k];

        t_31[k] = f_0 * dl_121[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, dl_122, dl_123, \
                         dl_124, dl_125, dl_126, dl_127, dl_128, \
                         dl_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * dl_122[k];

        t_33[k] = f_0 * dl_123[k];

        t_34[k] = f_0 * dl_124[k];

        t_35[k] = f_0 * dl_125[k];

        t_36[k] = f_0 * dl_126[k];

        t_37[k] = f_0 * dl_127[k];

        t_38[k] = f_0 * dl_128[k];

        t_39[k] = f_0 * dl_129[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, dl_130, dl_131, \
                         dl_132, dl_133, dl_134, dl_180, dl_181, \
                         dl_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * dl_130[k];

        t_41[k] = f_0 * dl_131[k];

        t_42[k] = f_0 * dl_132[k];

        t_43[k] = f_0 * dl_133[k];

        t_44[k] = f_0 * dl_134[k];

        t_45[k] = f_0 * dl_180[k];

        t_46[k] = f_0 * dl_181[k];

        t_47[k] = f_0 * dl_182[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, dl_183, dl_184, \
                         dl_185, dl_186, dl_187, dl_188, dl_189, \
                         dl_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * dl_183[k];

        t_49[k] = f_0 * dl_184[k];

        t_50[k] = f_0 * dl_185[k];

        t_51[k] = f_0 * dl_186[k];

        t_52[k] = f_0 * dl_187[k];

        t_53[k] = f_0 * dl_188[k];

        t_54[k] = f_0 * dl_189[k];

        t_55[k] = f_0 * dl_190[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, t_63, dl_191, dl_192, \
                         dl_193, dl_194, dl_195, dl_196, dl_197, \
                         dl_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * dl_191[k];

        t_57[k] = f_0 * dl_192[k];

        t_58[k] = f_0 * dl_193[k];

        t_59[k] = f_0 * dl_194[k];

        t_60[k] = f_0 * dl_195[k];

        t_61[k] = f_0 * dl_196[k];

        t_62[k] = f_0 * dl_197[k];

        t_63[k] = f_0 * dl_198[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, t_70, t_71, dl_199, dl_200, \
                         dl_201, dl_202, dl_203, dl_204, dl_205, \
                         dl_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * dl_199[k];

        t_65[k] = f_0 * dl_200[k];

        t_66[k] = f_0 * dl_201[k];

        t_67[k] = f_0 * dl_202[k];

        t_68[k] = f_0 * dl_203[k];

        t_69[k] = f_0 * dl_204[k];

        t_70[k] = f_0 * dl_205[k];

        t_71[k] = f_0 * dl_206[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, t_78, t_79, dl_207, dl_208, \
                         dl_209, dl_210, dl_211, dl_212, dl_213, \
                         dl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_0 * dl_207[k];

        t_73[k] = f_0 * dl_208[k];

        t_74[k] = f_0 * dl_209[k];

        t_75[k] = f_0 * dl_210[k];

        t_76[k] = f_0 * dl_211[k];

        t_77[k] = f_0 * dl_212[k];

        t_78[k] = f_0 * dl_213[k];

        t_79[k] = f_0 * dl_214[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, t_86, t_87, dl_215, dl_216, \
                         dl_217, dl_218, dl_219, dl_220, dl_221, \
                         dl_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * dl_215[k];

        t_81[k] = f_0 * dl_216[k];

        t_82[k] = f_0 * dl_217[k];

        t_83[k] = f_0 * dl_218[k];

        t_84[k] = f_0 * dl_219[k];

        t_85[k] = f_0 * dl_220[k];

        t_86[k] = f_0 * dl_221[k];

        t_87[k] = f_0 * dl_222[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, t_93, sl_0, sl_1, sl_2, sl_3, dl_223, \
                         dl_224, dl_225, dl_226, dl_227, dl_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_0 * dl_223[k];

        t_89[k] = f_0 * dl_224[k];

        t_90[k] = -sl_0[k]
                  + f_0 * dl_225[k];

        t_91[k] = -sl_1[k]
                  + f_0 * dl_226[k];

        t_92[k] = -sl_2[k]
                  + f_0 * dl_227[k];

        t_93[k] = -sl_3[k]
                  + f_0 * dl_228[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, sl_4, sl_5, sl_6, sl_7, sl_8, dl_229, \
                         dl_230, dl_231, dl_232, dl_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -sl_4[k]
                  + f_0 * dl_229[k];

        t_95[k] = -sl_5[k]
                  + f_0 * dl_230[k];

        t_96[k] = -sl_6[k]
                  + f_0 * dl_231[k];

        t_97[k] = -sl_7[k]
                  + f_0 * dl_232[k];

        t_98[k] = -sl_8[k]
                  + f_0 * dl_233[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, sl_9, sl_10, sl_11, sl_12, sl_13, \
                         dl_234, dl_235, dl_236, dl_237, dl_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -sl_9[k]
                  + f_0 * dl_234[k];

        t_100[k] = -sl_10[k]
                   + f_0 * dl_235[k];

        t_101[k] = -sl_11[k]
                   + f_0 * dl_236[k];

        t_102[k] = -sl_12[k]
                   + f_0 * dl_237[k];

        t_103[k] = -sl_13[k]
                   + f_0 * dl_238[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, sl_14, sl_15, sl_16, sl_17, sl_18, \
                         dl_239, dl_240, dl_241, dl_242, dl_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -sl_14[k]
                   + f_0 * dl_239[k];

        t_105[k] = -sl_15[k]
                   + f_0 * dl_240[k];

        t_106[k] = -sl_16[k]
                   + f_0 * dl_241[k];

        t_107[k] = -sl_17[k]
                   + f_0 * dl_242[k];

        t_108[k] = -sl_18[k]
                   + f_0 * dl_243[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, sl_19, sl_20, sl_21, sl_22, sl_23, \
                         dl_244, dl_245, dl_246, dl_247, dl_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -sl_19[k]
                   + f_0 * dl_244[k];

        t_110[k] = -sl_20[k]
                   + f_0 * dl_245[k];

        t_111[k] = -sl_21[k]
                   + f_0 * dl_246[k];

        t_112[k] = -sl_22[k]
                   + f_0 * dl_247[k];

        t_113[k] = -sl_23[k]
                   + f_0 * dl_248[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, sl_24, sl_25, sl_26, sl_27, sl_28, \
                         dl_249, dl_250, dl_251, dl_252, dl_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -sl_24[k]
                   + f_0 * dl_249[k];

        t_115[k] = -sl_25[k]
                   + f_0 * dl_250[k];

        t_116[k] = -sl_26[k]
                   + f_0 * dl_251[k];

        t_117[k] = -sl_27[k]
                   + f_0 * dl_252[k];

        t_118[k] = -sl_28[k]
                   + f_0 * dl_253[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, sl_29, sl_30, sl_31, sl_32, sl_33, \
                         dl_254, dl_255, dl_256, dl_257, dl_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -sl_29[k]
                   + f_0 * dl_254[k];

        t_120[k] = -sl_30[k]
                   + f_0 * dl_255[k];

        t_121[k] = -sl_31[k]
                   + f_0 * dl_256[k];

        t_122[k] = -sl_32[k]
                   + f_0 * dl_257[k];

        t_123[k] = -sl_33[k]
                   + f_0 * dl_258[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, sl_34, sl_35, sl_36, sl_37, sl_38, \
                         dl_259, dl_260, dl_261, dl_262, dl_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -sl_34[k]
                   + f_0 * dl_259[k];

        t_125[k] = -sl_35[k]
                   + f_0 * dl_260[k];

        t_126[k] = -sl_36[k]
                   + f_0 * dl_261[k];

        t_127[k] = -sl_37[k]
                   + f_0 * dl_262[k];

        t_128[k] = -sl_38[k]
                   + f_0 * dl_263[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, sl_39, sl_40, sl_41, sl_42, sl_43, \
                         dl_264, dl_265, dl_266, dl_267, dl_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -sl_39[k]
                   + f_0 * dl_264[k];

        t_130[k] = -sl_40[k]
                   + f_0 * dl_265[k];

        t_131[k] = -sl_41[k]
                   + f_0 * dl_266[k];

        t_132[k] = -sl_42[k]
                   + f_0 * dl_267[k];

        t_133[k] = -sl_43[k]
                   + f_0 * dl_268[k];
    }

#pragma omp simd aligned(t_134, sl_44, dl_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -sl_44[k]
                   + f_0 * dl_269[k];
    }
}

}  // namespace simdt2ceri
