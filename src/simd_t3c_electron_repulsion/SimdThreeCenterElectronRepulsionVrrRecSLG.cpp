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


#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_slg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skg0,
                                                          const size_t skf, const size_t skg1,
                                                          const size_t sld0, const size_t sld1,
                                                          const size_t slf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 3.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skg0_0 = buffer.data(skg0 + 0);
    const auto *skg0_3 = buffer.data(skg0 + 3);
    const auto *skg0_5 = buffer.data(skg0 + 5);
    const auto *skg0_10 = buffer.data(skg0 + 10);
    const auto *skg0_14 = buffer.data(skg0 + 14);
    const auto *skg0_18 = buffer.data(skg0 + 18);
    const auto *skg0_25 = buffer.data(skg0 + 25);
    const auto *skg0_30 = buffer.data(skg0 + 30);
    const auto *skg0_35 = buffer.data(skg0 + 35);
    const auto *skg0_44 = buffer.data(skg0 + 44);
    const auto *skg0_45 = buffer.data(skg0 + 45);
    const auto *skg0_48 = buffer.data(skg0 + 48);
    const auto *skg0_55 = buffer.data(skg0 + 55);
    const auto *skg0_75 = buffer.data(skg0 + 75);
    const auto *skg0_78 = buffer.data(skg0 + 78);

    const auto *skf_0 = buffer.data(skf + 0);
    const auto *skf_1 = buffer.data(skf + 1);
    const auto *skf_2 = buffer.data(skf + 2);
    const auto *skf_3 = buffer.data(skf + 3);
    const auto *skf_5 = buffer.data(skf + 5);
    const auto *skf_6 = buffer.data(skf + 6);
    const auto *skf_7 = buffer.data(skf + 7);
    const auto *skf_8 = buffer.data(skf + 8);
    const auto *skf_9 = buffer.data(skf + 9);
    const auto *skf_10 = buffer.data(skf + 10);
    const auto *skf_12 = buffer.data(skf + 12);
    const auto *skf_16 = buffer.data(skf + 16);
    const auto *skf_17 = buffer.data(skf + 17);
    const auto *skf_18 = buffer.data(skf + 18);
    const auto *skf_19 = buffer.data(skf + 19);
    const auto *skf_20 = buffer.data(skf + 20);
    const auto *skf_22 = buffer.data(skf + 22);
    const auto *skf_26 = buffer.data(skf + 26);
    const auto *skf_27 = buffer.data(skf + 27);
    const auto *skf_28 = buffer.data(skf + 28);
    const auto *skf_29 = buffer.data(skf + 29);
    const auto *skf_30 = buffer.data(skf + 30);
    const auto *skf_32 = buffer.data(skf + 32);
    const auto *skf_33 = buffer.data(skf + 33);
    const auto *skf_35 = buffer.data(skf + 35);
    const auto *skf_36 = buffer.data(skf + 36);
    const auto *skf_37 = buffer.data(skf + 37);
    const auto *skf_38 = buffer.data(skf + 38);
    const auto *skf_39 = buffer.data(skf + 39);
    const auto *skf_40 = buffer.data(skf + 40);
    const auto *skf_42 = buffer.data(skf + 42);
    const auto *skf_46 = buffer.data(skf + 46);
    const auto *skf_47 = buffer.data(skf + 47);
    const auto *skf_48 = buffer.data(skf + 48);
    const auto *skf_49 = buffer.data(skf + 49);
    const auto *skf_50 = buffer.data(skf + 50);
    const auto *skf_51 = buffer.data(skf + 51);
    const auto *skf_52 = buffer.data(skf + 52);
    const auto *skf_53 = buffer.data(skf + 53);
    const auto *skf_55 = buffer.data(skf + 55);
    const auto *skf_56 = buffer.data(skf + 56);
    const auto *skf_57 = buffer.data(skf + 57);
    const auto *skf_58 = buffer.data(skf + 58);
    const auto *skf_59 = buffer.data(skf + 59);
    const auto *skf_60 = buffer.data(skf + 60);
    const auto *skf_63 = buffer.data(skf + 63);
    const auto *skf_65 = buffer.data(skf + 65);
    const auto *skf_66 = buffer.data(skf + 66);
    const auto *skf_67 = buffer.data(skf + 67);
    const auto *skf_68 = buffer.data(skf + 68);
    const auto *skf_69 = buffer.data(skf + 69);
    const auto *skf_75 = buffer.data(skf + 75);
    const auto *skf_76 = buffer.data(skf + 76);
    const auto *skf_77 = buffer.data(skf + 77);
    const auto *skf_78 = buffer.data(skf + 78);
    const auto *skf_79 = buffer.data(skf + 79);

    const auto *skg1_0 = buffer.data(skg1 + 0);
    const auto *skg1_3 = buffer.data(skg1 + 3);
    const auto *skg1_5 = buffer.data(skg1 + 5);
    const auto *skg1_10 = buffer.data(skg1 + 10);
    const auto *skg1_14 = buffer.data(skg1 + 14);
    const auto *skg1_18 = buffer.data(skg1 + 18);
    const auto *skg1_25 = buffer.data(skg1 + 25);
    const auto *skg1_30 = buffer.data(skg1 + 30);
    const auto *skg1_35 = buffer.data(skg1 + 35);
    const auto *skg1_44 = buffer.data(skg1 + 44);
    const auto *skg1_45 = buffer.data(skg1 + 45);
    const auto *skg1_48 = buffer.data(skg1 + 48);
    const auto *skg1_55 = buffer.data(skg1 + 55);
    const auto *skg1_75 = buffer.data(skg1 + 75);
    const auto *skg1_78 = buffer.data(skg1 + 78);

    const auto *sld0_0 = buffer.data(sld0 + 0);
    const auto *sld0_3 = buffer.data(sld0 + 3);
    const auto *sld0_5 = buffer.data(sld0 + 5);
    const auto *sld0_9 = buffer.data(sld0 + 9);
    const auto *sld0_11 = buffer.data(sld0 + 11);
    const auto *sld0_17 = buffer.data(sld0 + 17);
    const auto *sld0_18 = buffer.data(sld0 + 18);
    const auto *sld0_21 = buffer.data(sld0 + 21);
    const auto *sld0_23 = buffer.data(sld0 + 23);
    const auto *sld0_29 = buffer.data(sld0 + 29);
    const auto *sld0_30 = buffer.data(sld0 + 30);
    const auto *sld0_33 = buffer.data(sld0 + 33);
    const auto *sld0_35 = buffer.data(sld0 + 35);
    const auto *sld0_36 = buffer.data(sld0 + 36);
    const auto *sld0_39 = buffer.data(sld0 + 39);
    const auto *sld0_41 = buffer.data(sld0 + 41);
    const auto *sld0_47 = buffer.data(sld0 + 47);

    const auto *sld1_0 = buffer.data(sld1 + 0);
    const auto *sld1_3 = buffer.data(sld1 + 3);
    const auto *sld1_5 = buffer.data(sld1 + 5);
    const auto *sld1_9 = buffer.data(sld1 + 9);
    const auto *sld1_11 = buffer.data(sld1 + 11);
    const auto *sld1_17 = buffer.data(sld1 + 17);
    const auto *sld1_18 = buffer.data(sld1 + 18);
    const auto *sld1_21 = buffer.data(sld1 + 21);
    const auto *sld1_23 = buffer.data(sld1 + 23);
    const auto *sld1_29 = buffer.data(sld1 + 29);
    const auto *sld1_30 = buffer.data(sld1 + 30);
    const auto *sld1_33 = buffer.data(sld1 + 33);
    const auto *sld1_35 = buffer.data(sld1 + 35);
    const auto *sld1_36 = buffer.data(sld1 + 36);
    const auto *sld1_39 = buffer.data(sld1 + 39);
    const auto *sld1_41 = buffer.data(sld1 + 41);
    const auto *sld1_47 = buffer.data(sld1 + 47);

    const auto *slf_0 = buffer.data(slf + 0);
    const auto *slf_2 = buffer.data(slf + 2);
    const auto *slf_3 = buffer.data(slf + 3);
    const auto *slf_5 = buffer.data(slf + 5);
    const auto *slf_6 = buffer.data(slf + 6);
    const auto *slf_7 = buffer.data(slf + 7);
    const auto *slf_8 = buffer.data(slf + 8);
    const auto *slf_9 = buffer.data(slf + 9);
    const auto *slf_10 = buffer.data(slf + 10);
    const auto *slf_12 = buffer.data(slf + 12);
    const auto *slf_16 = buffer.data(slf + 16);
    const auto *slf_17 = buffer.data(slf + 17);
    const auto *slf_18 = buffer.data(slf + 18);
    const auto *slf_19 = buffer.data(slf + 19);
    const auto *slf_20 = buffer.data(slf + 20);
    const auto *slf_22 = buffer.data(slf + 22);
    const auto *slf_26 = buffer.data(slf + 26);
    const auto *slf_27 = buffer.data(slf + 27);
    const auto *slf_28 = buffer.data(slf + 28);
    const auto *slf_29 = buffer.data(slf + 29);
    const auto *slf_30 = buffer.data(slf + 30);
    const auto *slf_32 = buffer.data(slf + 32);
    const auto *slf_33 = buffer.data(slf + 33);
    const auto *slf_35 = buffer.data(slf + 35);
    const auto *slf_36 = buffer.data(slf + 36);
    const auto *slf_37 = buffer.data(slf + 37);
    const auto *slf_38 = buffer.data(slf + 38);
    const auto *slf_39 = buffer.data(slf + 39);
    const auto *slf_40 = buffer.data(slf + 40);
    const auto *slf_42 = buffer.data(slf + 42);
    const auto *slf_46 = buffer.data(slf + 46);
    const auto *slf_47 = buffer.data(slf + 47);
    const auto *slf_48 = buffer.data(slf + 48);
    const auto *slf_49 = buffer.data(slf + 49);
    const auto *slf_50 = buffer.data(slf + 50);
    const auto *slf_52 = buffer.data(slf + 52);
    const auto *slf_53 = buffer.data(slf + 53);
    const auto *slf_55 = buffer.data(slf + 55);
    const auto *slf_56 = buffer.data(slf + 56);
    const auto *slf_57 = buffer.data(slf + 57);
    const auto *slf_58 = buffer.data(slf + 58);
    const auto *slf_59 = buffer.data(slf + 59);
    const auto *slf_60 = buffer.data(slf + 60);
    const auto *slf_62 = buffer.data(slf + 62);
    const auto *slf_63 = buffer.data(slf + 63);
    const auto *slf_65 = buffer.data(slf + 65);
    const auto *slf_66 = buffer.data(slf + 66);
    const auto *slf_67 = buffer.data(slf + 67);
    const auto *slf_68 = buffer.data(slf + 68);
    const auto *slf_69 = buffer.data(slf + 69);
    const auto *slf_70 = buffer.data(slf + 70);
    const auto *slf_72 = buffer.data(slf + 72);
    const auto *slf_75 = buffer.data(slf + 75);
    const auto *slf_76 = buffer.data(slf + 76);
    const auto *slf_77 = buffer.data(slf + 77);
    const auto *slf_78 = buffer.data(slf + 78);
    const auto *slf_79 = buffer.data(slf + 79);
    const auto *slf_80 = buffer.data(slf + 80);
    const auto *slf_82 = buffer.data(slf + 82);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, skf_0, skf_3, sld0_0, sld0_3, \
                         sld1_0, sld1_3, slf_0, slf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * skf_0[k]
                 + f_1 * sld0_0[k]
                 - f_2 * sld1_0[k]
                 + f_3 * pc_x[k] * slf_0[k];

        t_1[k] = f_3 * pc_y[k] * slf_0[k];

        t_2[k] = f_3 * pc_z[k] * slf_0[k];

        t_3[k] = f_0 * skf_3[k]
                 + f_4 * sld0_3[k]
                 - f_5 * sld1_3[k]
                 + f_3 * pc_x[k] * slf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pc_x, pc_y, skf_5, skf_6, skf_7, sld0_5, sld1_5, \
                         slf_2, slf_5, slf_6, slf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * slf_2[k];

        t_5[k] = f_0 * skf_5[k]
                 + f_4 * sld0_5[k]
                 - f_5 * sld1_5[k]
                 + f_3 * pc_x[k] * slf_5[k];

        t_6[k] = f_0 * skf_6[k]
                 + f_3 * pc_x[k] * slf_6[k];

        t_7[k] = f_0 * skf_7[k]
                 + f_3 * pc_x[k] * slf_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pc_x, pc_y, pc_z, skf_8, skf_9, sld0_3, sld1_3, \
                         slf_6, slf_8, slf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * skf_8[k]
                 + f_3 * pc_x[k] * slf_8[k];

        t_9[k] = f_0 * skf_9[k]
                 + f_3 * pc_x[k] * slf_9[k];

        t_10[k] = f_1 * sld0_3[k]
                  - f_2 * sld1_3[k]
                  + f_3 * pc_y[k] * slf_6[k];

        t_11[k] = f_3 * pc_z[k] * slf_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_y, pc_y, pc_z, skg0_0, skf_0, \
                         skg1_0, sld0_5, sld1_5, slf_8, slf_9, slf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_4 * sld0_5[k]
                  - f_5 * sld1_5[k]
                  + f_3 * pc_y[k] * slf_8[k];

        t_13[k] = f_3 * pc_y[k] * slf_9[k];

        t_14[k] = f_1 * sld0_5[k]
                  - f_2 * sld1_5[k]
                  + f_3 * pc_z[k] * slf_9[k];

        t_15[k] = pb_y[k] * skg0_0[k]
                  - f_6 * pc_y[k] * skg1_0[k];

        t_16[k] = f_7 * skf_0[k]
                  + f_3 * pc_y[k] * slf_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pc_y, pc_z, skg0_3, skg0_5, skf_1, \
                         skf_2, skg1_3, skg1_5, slf_10, slf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * pc_z[k] * slf_10[k];

        t_18[k] = pb_y[k] * skg0_3[k]
                  + f_8 * skf_1[k]
                  - f_6 * pc_y[k] * skg1_3[k];

        t_19[k] = f_7 * skf_2[k]
                  + f_3 * pc_y[k] * slf_12[k];

        t_20[k] = pb_y[k] * skg0_5[k]
                  - f_6 * pc_y[k] * skg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, skf_16, skf_17, skf_18, skf_19, slf_16, \
                         slf_17, slf_18, slf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * skf_16[k]
                  + f_3 * pc_x[k] * slf_16[k];

        t_22[k] = f_9 * skf_17[k]
                  + f_3 * pc_x[k] * slf_17[k];

        t_23[k] = f_9 * skf_18[k]
                  + f_3 * pc_x[k] * slf_18[k];

        t_24[k] = f_9 * skf_19[k]
                  + f_3 * pc_x[k] * slf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, skf_6, skf_8, skf_9, sld0_9, \
                         sld0_11, sld1_9, sld1_11, slf_16, slf_18, \
                         slf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * skf_6[k]
                  + f_1 * sld0_9[k]
                  - f_2 * sld1_9[k]
                  + f_3 * pc_y[k] * slf_16[k];

        t_26[k] = f_3 * pc_z[k] * slf_16[k];

        t_27[k] = f_7 * skf_8[k]
                  + f_4 * sld0_11[k]
                  - f_5 * sld1_11[k]
                  + f_3 * pc_y[k] * slf_18[k];

        t_28[k] = f_7 * skf_9[k]
                  + f_3 * pc_y[k] * slf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pb_y, pb_z, pc_y, pc_z, skg0_0, skg0_14, \
                         skf_0, skg1_0, skg1_14, slf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * skg0_14[k]
                  - f_6 * pc_y[k] * skg1_14[k];

        t_30[k] = pb_z[k] * skg0_0[k]
                  - f_6 * pc_z[k] * skg1_0[k];

        t_31[k] = f_3 * pc_y[k] * slf_20[k];

        t_32[k] = f_7 * skf_0[k]
                  + f_3 * pc_z[k] * slf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_z, pc_x, pc_y, pc_z, skg0_3, skg0_5, \
                         skf_2, skf_26, skg1_3, skg1_5, slf_22, \
                         slf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_z[k] * skg0_3[k]
                  - f_6 * pc_z[k] * skg1_3[k];

        t_34[k] = f_3 * pc_y[k] * slf_22[k];

        t_35[k] = pb_z[k] * skg0_5[k]
                  + f_8 * skf_2[k]
                  - f_6 * pc_z[k] * skg1_5[k];

        t_36[k] = f_9 * skf_26[k]
                  + f_3 * pc_x[k] * slf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_z, pc_x, pc_z, skg0_10, skf_27, skf_28, \
                         skf_29, skg1_10, slf_27, slf_28, slf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * skf_27[k]
                  + f_3 * pc_x[k] * slf_27[k];

        t_38[k] = f_9 * skf_28[k]
                  + f_3 * pc_x[k] * slf_28[k];

        t_39[k] = f_9 * skf_29[k]
                  + f_3 * pc_x[k] * slf_29[k];

        t_40[k] = pb_z[k] * skg0_10[k]
                  - f_6 * pc_z[k] * skg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, skf_6, skf_9, sld0_17, sld1_17, \
                         slf_26, slf_28, slf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * skf_6[k]
                  + f_3 * pc_z[k] * slf_26[k];

        t_42[k] = f_4 * sld0_17[k]
                  - f_5 * sld1_17[k]
                  + f_3 * pc_y[k] * slf_28[k];

        t_43[k] = f_3 * pc_y[k] * slf_29[k];

        t_44[k] = f_7 * skf_9[k]
                  + f_1 * sld0_17[k]
                  - f_2 * sld1_17[k]
                  + f_3 * pc_z[k] * slf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, skf_10, skf_30, skf_33, \
                         sld0_18, sld0_21, sld1_18, sld1_21, slf_30, \
                         slf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_10 * skf_30[k]
                  + f_1 * sld0_18[k]
                  - f_2 * sld1_18[k]
                  + f_3 * pc_x[k] * slf_30[k];

        t_46[k] = f_8 * skf_10[k]
                  + f_3 * pc_y[k] * slf_30[k];

        t_47[k] = f_3 * pc_z[k] * slf_30[k];

        t_48[k] = f_10 * skf_33[k]
                  + f_4 * sld0_21[k]
                  - f_5 * sld1_21[k]
                  + f_3 * pc_x[k] * slf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, skf_12, skf_35, skf_36, skf_37, \
                         sld0_23, sld1_23, slf_32, slf_35, slf_36, \
                         slf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * skf_12[k]
                  + f_3 * pc_y[k] * slf_32[k];

        t_50[k] = f_10 * skf_35[k]
                  + f_4 * sld0_23[k]
                  - f_5 * sld1_23[k]
                  + f_3 * pc_x[k] * slf_35[k];

        t_51[k] = f_10 * skf_36[k]
                  + f_3 * pc_x[k] * slf_36[k];

        t_52[k] = f_10 * skf_37[k]
                  + f_3 * pc_x[k] * slf_37[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pc_x, pc_y, pc_z, skf_16, skf_38, skf_39, \
                         sld0_21, sld1_21, slf_36, slf_38, slf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * skf_38[k]
                  + f_3 * pc_x[k] * slf_38[k];

        t_54[k] = f_10 * skf_39[k]
                  + f_3 * pc_x[k] * slf_39[k];

        t_55[k] = f_8 * skf_16[k]
                  + f_1 * sld0_21[k]
                  - f_2 * sld1_21[k]
                  + f_3 * pc_y[k] * slf_36[k];

        t_56[k] = f_3 * pc_z[k] * slf_36[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, pc_z, skg0_30, skf_18, skf_19, \
                         skg1_30, sld0_23, sld1_23, slf_38, slf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_8 * skf_18[k]
                  + f_4 * sld0_23[k]
                  - f_5 * sld1_23[k]
                  + f_3 * pc_y[k] * slf_38[k];

        t_58[k] = f_8 * skf_19[k]
                  + f_3 * pc_y[k] * slf_39[k];

        t_59[k] = f_1 * sld0_23[k]
                  - f_2 * sld1_23[k]
                  + f_3 * pc_z[k] * slf_39[k];

        t_60[k] = pb_y[k] * skg0_30[k]
                  - f_6 * pc_y[k] * skg1_30[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, skg0_18, skf_10, skf_20, \
                         skf_22, skg1_18, slf_40, slf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_7 * skf_20[k]
                  + f_3 * pc_y[k] * slf_40[k];

        t_62[k] = f_7 * skf_10[k]
                  + f_3 * pc_z[k] * slf_40[k];

        t_63[k] = pb_z[k] * skg0_18[k]
                  - f_6 * pc_z[k] * skg1_18[k];

        t_64[k] = f_7 * skf_22[k]
                  + f_3 * pc_y[k] * slf_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_x, pc_y, skg0_35, skf_46, skf_47, \
                         skf_48, skg1_35, slf_46, slf_47, slf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * skg0_35[k]
                  - f_6 * pc_y[k] * skg1_35[k];

        t_66[k] = f_10 * skf_46[k]
                  + f_3 * pc_x[k] * slf_46[k];

        t_67[k] = f_10 * skf_47[k]
                  + f_3 * pc_x[k] * slf_47[k];

        t_68[k] = f_10 * skf_48[k]
                  + f_3 * pc_x[k] * slf_48[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_z, pc_x, pc_z, skg0_25, skf_16, skf_49, skg1_25, \
                         slf_46, slf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * skf_49[k]
                  + f_3 * pc_x[k] * slf_49[k];

        t_70[k] = pb_z[k] * skg0_25[k]
                  - f_6 * pc_z[k] * skg1_25[k];

        t_71[k] = f_7 * skf_16[k]
                  + f_3 * pc_z[k] * slf_46[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_y, pc_y, skg0_44, skf_28, skf_29, skg1_44, \
                         sld0_29, sld1_29, slf_48, slf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * skf_28[k]
                  + f_4 * sld0_29[k]
                  - f_5 * sld1_29[k]
                  + f_3 * pc_y[k] * slf_48[k];

        t_73[k] = f_7 * skf_29[k]
                  + f_3 * pc_y[k] * slf_49[k];

        t_74[k] = pb_y[k] * skg0_44[k]
                  - f_6 * pc_y[k] * skg1_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, pc_z, skf_20, skf_50, skf_53, \
                         sld0_30, sld0_33, sld1_30, sld1_33, slf_50, \
                         slf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_10 * skf_50[k]
                  + f_1 * sld0_30[k]
                  - f_2 * sld1_30[k]
                  + f_3 * pc_x[k] * slf_50[k];

        t_76[k] = f_3 * pc_y[k] * slf_50[k];

        t_77[k] = f_8 * skf_20[k]
                  + f_3 * pc_z[k] * slf_50[k];

        t_78[k] = f_10 * skf_53[k]
                  + f_4 * sld0_33[k]
                  - f_5 * sld1_33[k]
                  + f_3 * pc_x[k] * slf_53[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, skf_55, skf_56, skf_57, sld0_35, \
                         sld1_35, slf_52, slf_55, slf_56, slf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_y[k] * slf_52[k];

        t_80[k] = f_10 * skf_55[k]
                  + f_4 * sld0_35[k]
                  - f_5 * sld1_35[k]
                  + f_3 * pc_x[k] * slf_55[k];

        t_81[k] = f_10 * skf_56[k]
                  + f_3 * pc_x[k] * slf_56[k];

        t_82[k] = f_10 * skf_57[k]
                  + f_3 * pc_x[k] * slf_57[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, pc_z, skf_26, skf_58, skf_59, \
                         sld0_33, sld1_33, slf_56, slf_58, slf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_10 * skf_58[k]
                  + f_3 * pc_x[k] * slf_58[k];

        t_84[k] = f_10 * skf_59[k]
                  + f_3 * pc_x[k] * slf_59[k];

        t_85[k] = f_1 * sld0_33[k]
                  - f_2 * sld1_33[k]
                  + f_3 * pc_y[k] * slf_56[k];

        t_86[k] = f_8 * skf_26[k]
                  + f_3 * pc_z[k] * slf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, skf_29, skf_60, sld0_35, \
                         sld0_36, sld1_35, sld1_36, slf_58, slf_59, \
                         slf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * sld0_35[k]
                  - f_5 * sld1_35[k]
                  + f_3 * pc_y[k] * slf_58[k];

        t_88[k] = f_3 * pc_y[k] * slf_59[k];

        t_89[k] = f_8 * skf_29[k]
                  + f_1 * sld0_35[k]
                  - f_2 * sld1_35[k]
                  + f_3 * pc_z[k] * slf_59[k];

        t_90[k] = f_11 * skf_60[k]
                  + f_1 * sld0_36[k]
                  - f_2 * sld1_36[k]
                  + f_3 * pc_x[k] * slf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, skf_30, skf_32, skf_63, \
                         sld0_39, sld1_39, slf_60, slf_62, slf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_12 * skf_30[k]
                  + f_3 * pc_y[k] * slf_60[k];

        t_92[k] = f_3 * pc_z[k] * slf_60[k];

        t_93[k] = f_11 * skf_63[k]
                  + f_4 * sld0_39[k]
                  - f_5 * sld1_39[k]
                  + f_3 * pc_x[k] * slf_63[k];

        t_94[k] = f_12 * skf_32[k]
                  + f_3 * pc_y[k] * slf_62[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, skf_65, skf_66, skf_67, skf_68, \
                         sld0_41, sld1_41, slf_65, slf_66, slf_67, \
                         slf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_11 * skf_65[k]
                  + f_4 * sld0_41[k]
                  - f_5 * sld1_41[k]
                  + f_3 * pc_x[k] * slf_65[k];

        t_96[k] = f_11 * skf_66[k]
                  + f_3 * pc_x[k] * slf_66[k];

        t_97[k] = f_11 * skf_67[k]
                  + f_3 * pc_x[k] * slf_67[k];

        t_98[k] = f_11 * skf_68[k]
                  + f_3 * pc_x[k] * slf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pc_x, pc_y, pc_z, skf_36, skf_69, sld0_39, \
                         sld1_39, slf_66, slf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_11 * skf_69[k]
                  + f_3 * pc_x[k] * slf_69[k];

        t_100[k] = f_12 * skf_36[k]
                   + f_1 * sld0_39[k]
                   - f_2 * sld1_39[k]
                   + f_3 * pc_y[k] * slf_66[k];

        t_101[k] = f_3 * pc_z[k] * slf_66[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pb_z, pc_y, pc_z, skg0_45, skf_38, \
                         skf_39, skg1_45, sld0_41, sld1_41, slf_68, \
                         slf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_12 * skf_38[k]
                   + f_4 * sld0_41[k]
                   - f_5 * sld1_41[k]
                   + f_3 * pc_y[k] * slf_68[k];

        t_103[k] = f_12 * skf_39[k]
                   + f_3 * pc_y[k] * slf_69[k];

        t_104[k] = f_1 * sld0_41[k]
                   - f_2 * sld1_41[k]
                   + f_3 * pc_z[k] * slf_69[k];

        t_105[k] = pb_z[k] * skg0_45[k]
                   - f_6 * pc_z[k] * skg1_45[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_z, pc_y, pc_z, skg0_48, skf_30, \
                         skf_40, skf_42, skg1_48, slf_70, slf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_8 * skf_40[k]
                   + f_3 * pc_y[k] * slf_70[k];

        t_107[k] = f_7 * skf_30[k]
                   + f_3 * pc_z[k] * slf_70[k];

        t_108[k] = pb_z[k] * skg0_48[k]
                   - f_6 * pc_z[k] * skg1_48[k];

        t_109[k] = f_8 * skf_42[k]
                   + f_3 * pc_y[k] * slf_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, skf_75, skf_76, skf_77, skf_78, \
                         sld0_47, sld1_47, slf_75, slf_76, slf_77, \
                         slf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_11 * skf_75[k]
                   + f_4 * sld0_47[k]
                   - f_5 * sld1_47[k]
                   + f_3 * pc_x[k] * slf_75[k];

        t_111[k] = f_11 * skf_76[k]
                   + f_3 * pc_x[k] * slf_76[k];

        t_112[k] = f_11 * skf_77[k]
                   + f_3 * pc_x[k] * slf_77[k];

        t_113[k] = f_11 * skf_78[k]
                   + f_3 * pc_x[k] * slf_78[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_z, pc_x, pc_z, skg0_55, skf_36, skf_79, \
                         skg1_55, slf_76, slf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_11 * skf_79[k]
                   + f_3 * pc_x[k] * slf_79[k];

        t_115[k] = pb_z[k] * skg0_55[k]
                   - f_6 * pc_z[k] * skg1_55[k];

        t_116[k] = f_7 * skf_36[k]
                   + f_3 * pc_z[k] * slf_76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_y, pc_y, pc_z, skg0_75, skf_39, \
                         skf_48, skf_49, skg1_75, sld0_47, sld1_47, slf_78, \
                         slf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_8 * skf_48[k]
                   + f_4 * sld0_47[k]
                   - f_5 * sld1_47[k]
                   + f_3 * pc_y[k] * slf_78[k];

        t_118[k] = f_8 * skf_49[k]
                   + f_3 * pc_y[k] * slf_79[k];

        t_119[k] = f_7 * skf_39[k]
                   + f_1 * sld0_47[k]
                   - f_2 * sld1_47[k]
                   + f_3 * pc_z[k] * slf_79[k];

        t_120[k] = pb_y[k] * skg0_75[k]
                   - f_6 * pc_y[k] * skg1_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, pc_y, pc_z, skg0_78, skf_40, \
                         skf_50, skf_51, skf_52, skg1_78, slf_80, \
                         slf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_7 * skf_50[k]
                   + f_3 * pc_y[k] * slf_80[k];

        t_122[k] = f_8 * skf_40[k]
                   + f_3 * pc_z[k] * slf_80[k];

        t_123[k] = pb_y[k] * skg0_78[k]
                   + f_8 * skf_51[k]
                   - f_6 * pc_y[k] * skg1_78[k];

        t_124[k] = f_7 * skf_52[k]
                   + f_3 * pc_y[k] * slf_82[k];
    }
}

static auto
compute_prim_slg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skg0,
                                                          const size_t skf, const size_t skg1,
                                                          const size_t sld0, const size_t sld1,
                                                          const size_t slf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.0 / q;

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
    auto *t_135 = buffer.data(target + 135);
    auto *t_136 = buffer.data(target + 136);
    auto *t_137 = buffer.data(target + 137);
    auto *t_138 = buffer.data(target + 138);
    auto *t_139 = buffer.data(target + 139);
    auto *t_140 = buffer.data(target + 140);
    auto *t_141 = buffer.data(target + 141);
    auto *t_142 = buffer.data(target + 142);
    auto *t_143 = buffer.data(target + 143);
    auto *t_144 = buffer.data(target + 144);
    auto *t_145 = buffer.data(target + 145);
    auto *t_146 = buffer.data(target + 146);
    auto *t_147 = buffer.data(target + 147);
    auto *t_148 = buffer.data(target + 148);
    auto *t_149 = buffer.data(target + 149);
    auto *t_150 = buffer.data(target + 150);
    auto *t_151 = buffer.data(target + 151);
    auto *t_152 = buffer.data(target + 152);
    auto *t_153 = buffer.data(target + 153);
    auto *t_154 = buffer.data(target + 154);
    auto *t_155 = buffer.data(target + 155);
    auto *t_156 = buffer.data(target + 156);
    auto *t_157 = buffer.data(target + 157);
    auto *t_158 = buffer.data(target + 158);
    auto *t_159 = buffer.data(target + 159);
    auto *t_160 = buffer.data(target + 160);
    auto *t_161 = buffer.data(target + 161);
    auto *t_162 = buffer.data(target + 162);
    auto *t_163 = buffer.data(target + 163);
    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);
    auto *t_168 = buffer.data(target + 168);
    auto *t_169 = buffer.data(target + 169);
    auto *t_170 = buffer.data(target + 170);
    auto *t_171 = buffer.data(target + 171);
    auto *t_172 = buffer.data(target + 172);
    auto *t_173 = buffer.data(target + 173);
    auto *t_174 = buffer.data(target + 174);
    auto *t_175 = buffer.data(target + 175);
    auto *t_176 = buffer.data(target + 176);
    auto *t_177 = buffer.data(target + 177);
    auto *t_178 = buffer.data(target + 178);
    auto *t_179 = buffer.data(target + 179);
    auto *t_180 = buffer.data(target + 180);
    auto *t_181 = buffer.data(target + 181);
    auto *t_182 = buffer.data(target + 182);
    auto *t_183 = buffer.data(target + 183);
    auto *t_184 = buffer.data(target + 184);
    auto *t_185 = buffer.data(target + 185);
    auto *t_186 = buffer.data(target + 186);
    auto *t_187 = buffer.data(target + 187);
    auto *t_188 = buffer.data(target + 188);
    auto *t_189 = buffer.data(target + 189);
    auto *t_190 = buffer.data(target + 190);
    auto *t_191 = buffer.data(target + 191);
    auto *t_192 = buffer.data(target + 192);
    auto *t_193 = buffer.data(target + 193);
    auto *t_194 = buffer.data(target + 194);
    auto *t_195 = buffer.data(target + 195);
    auto *t_196 = buffer.data(target + 196);
    auto *t_197 = buffer.data(target + 197);
    auto *t_198 = buffer.data(target + 198);
    auto *t_199 = buffer.data(target + 199);
    auto *t_200 = buffer.data(target + 200);
    auto *t_201 = buffer.data(target + 201);
    auto *t_202 = buffer.data(target + 202);
    auto *t_203 = buffer.data(target + 203);
    auto *t_204 = buffer.data(target + 204);
    auto *t_205 = buffer.data(target + 205);
    auto *t_206 = buffer.data(target + 206);
    auto *t_207 = buffer.data(target + 207);
    auto *t_208 = buffer.data(target + 208);
    auto *t_209 = buffer.data(target + 209);
    auto *t_210 = buffer.data(target + 210);
    auto *t_211 = buffer.data(target + 211);
    auto *t_212 = buffer.data(target + 212);
    auto *t_213 = buffer.data(target + 213);
    auto *t_214 = buffer.data(target + 214);
    auto *t_215 = buffer.data(target + 215);
    auto *t_216 = buffer.data(target + 216);
    auto *t_217 = buffer.data(target + 217);
    auto *t_218 = buffer.data(target + 218);
    auto *t_219 = buffer.data(target + 219);
    auto *t_220 = buffer.data(target + 220);
    auto *t_221 = buffer.data(target + 221);
    auto *t_222 = buffer.data(target + 222);
    auto *t_223 = buffer.data(target + 223);
    auto *t_224 = buffer.data(target + 224);
    auto *t_225 = buffer.data(target + 225);
    auto *t_226 = buffer.data(target + 226);
    auto *t_227 = buffer.data(target + 227);
    auto *t_228 = buffer.data(target + 228);
    auto *t_229 = buffer.data(target + 229);
    auto *t_230 = buffer.data(target + 230);
    auto *t_231 = buffer.data(target + 231);
    auto *t_232 = buffer.data(target + 232);
    auto *t_233 = buffer.data(target + 233);
    auto *t_234 = buffer.data(target + 234);
    auto *t_235 = buffer.data(target + 235);
    auto *t_236 = buffer.data(target + 236);
    auto *t_237 = buffer.data(target + 237);
    auto *t_238 = buffer.data(target + 238);
    auto *t_239 = buffer.data(target + 239);
    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skg0_80 = buffer.data(skg0 + 80);
    const auto *skg0_89 = buffer.data(skg0 + 89);
    const auto *skg0_90 = buffer.data(skg0 + 90);
    const auto *skg0_93 = buffer.data(skg0 + 93);
    const auto *skg0_100 = buffer.data(skg0 + 100);
    const auto *skg0_135 = buffer.data(skg0 + 135);
    const auto *skg0_138 = buffer.data(skg0 + 138);
    const auto *skg0_140 = buffer.data(skg0 + 140);
    const auto *skg0_149 = buffer.data(skg0 + 149);
    const auto *skg0_150 = buffer.data(skg0 + 150);
    const auto *skg0_153 = buffer.data(skg0 + 153);

    const auto *skf_46 = buffer.data(skf + 46);
    const auto *skf_50 = buffer.data(skf + 50);
    const auto *skf_56 = buffer.data(skf + 56);
    const auto *skf_58 = buffer.data(skf + 58);
    const auto *skf_59 = buffer.data(skf + 59);
    const auto *skf_60 = buffer.data(skf + 60);
    const auto *skf_62 = buffer.data(skf + 62);
    const auto *skf_66 = buffer.data(skf + 66);
    const auto *skf_68 = buffer.data(skf + 68);
    const auto *skf_69 = buffer.data(skf + 69);
    const auto *skf_70 = buffer.data(skf + 70);
    const auto *skf_72 = buffer.data(skf + 72);
    const auto *skf_76 = buffer.data(skf + 76);
    const auto *skf_78 = buffer.data(skf + 78);
    const auto *skf_79 = buffer.data(skf + 79);
    const auto *skf_80 = buffer.data(skf + 80);
    const auto *skf_82 = buffer.data(skf + 82);
    const auto *skf_86 = buffer.data(skf + 86);
    const auto *skf_87 = buffer.data(skf + 87);
    const auto *skf_88 = buffer.data(skf + 88);
    const auto *skf_89 = buffer.data(skf + 89);
    const auto *skf_90 = buffer.data(skf + 90);
    const auto *skf_91 = buffer.data(skf + 91);
    const auto *skf_92 = buffer.data(skf + 92);
    const auto *skf_93 = buffer.data(skf + 93);
    const auto *skf_95 = buffer.data(skf + 95);
    const auto *skf_96 = buffer.data(skf + 96);
    const auto *skf_97 = buffer.data(skf + 97);
    const auto *skf_98 = buffer.data(skf + 98);
    const auto *skf_99 = buffer.data(skf + 99);
    const auto *skf_100 = buffer.data(skf + 100);
    const auto *skf_102 = buffer.data(skf + 102);
    const auto *skf_103 = buffer.data(skf + 103);
    const auto *skf_105 = buffer.data(skf + 105);
    const auto *skf_106 = buffer.data(skf + 106);
    const auto *skf_107 = buffer.data(skf + 107);
    const auto *skf_108 = buffer.data(skf + 108);
    const auto *skf_109 = buffer.data(skf + 109);
    const auto *skf_110 = buffer.data(skf + 110);
    const auto *skf_112 = buffer.data(skf + 112);
    const auto *skf_115 = buffer.data(skf + 115);
    const auto *skf_116 = buffer.data(skf + 116);
    const auto *skf_117 = buffer.data(skf + 117);
    const auto *skf_118 = buffer.data(skf + 118);
    const auto *skf_119 = buffer.data(skf + 119);
    const auto *skf_120 = buffer.data(skf + 120);
    const auto *skf_123 = buffer.data(skf + 123);
    const auto *skf_125 = buffer.data(skf + 125);
    const auto *skf_126 = buffer.data(skf + 126);
    const auto *skf_127 = buffer.data(skf + 127);
    const auto *skf_128 = buffer.data(skf + 128);
    const auto *skf_129 = buffer.data(skf + 129);
    const auto *skf_136 = buffer.data(skf + 136);
    const auto *skf_137 = buffer.data(skf + 137);
    const auto *skf_138 = buffer.data(skf + 138);
    const auto *skf_139 = buffer.data(skf + 139);
    const auto *skf_140 = buffer.data(skf + 140);
    const auto *skf_143 = buffer.data(skf + 143);
    const auto *skf_145 = buffer.data(skf + 145);
    const auto *skf_146 = buffer.data(skf + 146);
    const auto *skf_147 = buffer.data(skf + 147);
    const auto *skf_148 = buffer.data(skf + 148);
    const auto *skf_149 = buffer.data(skf + 149);
    const auto *skf_150 = buffer.data(skf + 150);
    const auto *skf_153 = buffer.data(skf + 153);
    const auto *skf_155 = buffer.data(skf + 155);
    const auto *skf_156 = buffer.data(skf + 156);
    const auto *skf_157 = buffer.data(skf + 157);
    const auto *skf_158 = buffer.data(skf + 158);
    const auto *skf_159 = buffer.data(skf + 159);

    const auto *skg1_80 = buffer.data(skg1 + 80);
    const auto *skg1_89 = buffer.data(skg1 + 89);
    const auto *skg1_90 = buffer.data(skg1 + 90);
    const auto *skg1_93 = buffer.data(skg1 + 93);
    const auto *skg1_100 = buffer.data(skg1 + 100);
    const auto *skg1_135 = buffer.data(skg1 + 135);
    const auto *skg1_138 = buffer.data(skg1 + 138);
    const auto *skg1_140 = buffer.data(skg1 + 140);
    const auto *skg1_149 = buffer.data(skg1 + 149);
    const auto *skg1_150 = buffer.data(skg1 + 150);
    const auto *skg1_153 = buffer.data(skg1 + 153);

    const auto *sld0_51 = buffer.data(sld0 + 51);
    const auto *sld0_53 = buffer.data(sld0 + 53);
    const auto *sld0_54 = buffer.data(sld0 + 54);
    const auto *sld0_57 = buffer.data(sld0 + 57);
    const auto *sld0_59 = buffer.data(sld0 + 59);
    const auto *sld0_60 = buffer.data(sld0 + 60);
    const auto *sld0_63 = buffer.data(sld0 + 63);
    const auto *sld0_65 = buffer.data(sld0 + 65);
    const auto *sld0_71 = buffer.data(sld0 + 71);
    const auto *sld0_72 = buffer.data(sld0 + 72);
    const auto *sld0_75 = buffer.data(sld0 + 75);
    const auto *sld0_77 = buffer.data(sld0 + 77);
    const auto *sld0_81 = buffer.data(sld0 + 81);
    const auto *sld0_83 = buffer.data(sld0 + 83);
    const auto *sld0_84 = buffer.data(sld0 + 84);
    const auto *sld0_87 = buffer.data(sld0 + 87);
    const auto *sld0_89 = buffer.data(sld0 + 89);
    const auto *sld0_90 = buffer.data(sld0 + 90);
    const auto *sld0_93 = buffer.data(sld0 + 93);
    const auto *sld0_95 = buffer.data(sld0 + 95);

    const auto *sld1_51 = buffer.data(sld1 + 51);
    const auto *sld1_53 = buffer.data(sld1 + 53);
    const auto *sld1_54 = buffer.data(sld1 + 54);
    const auto *sld1_57 = buffer.data(sld1 + 57);
    const auto *sld1_59 = buffer.data(sld1 + 59);
    const auto *sld1_60 = buffer.data(sld1 + 60);
    const auto *sld1_63 = buffer.data(sld1 + 63);
    const auto *sld1_65 = buffer.data(sld1 + 65);
    const auto *sld1_71 = buffer.data(sld1 + 71);
    const auto *sld1_72 = buffer.data(sld1 + 72);
    const auto *sld1_75 = buffer.data(sld1 + 75);
    const auto *sld1_77 = buffer.data(sld1 + 77);
    const auto *sld1_81 = buffer.data(sld1 + 81);
    const auto *sld1_83 = buffer.data(sld1 + 83);
    const auto *sld1_84 = buffer.data(sld1 + 84);
    const auto *sld1_87 = buffer.data(sld1 + 87);
    const auto *sld1_89 = buffer.data(sld1 + 89);
    const auto *sld1_90 = buffer.data(sld1 + 90);
    const auto *sld1_93 = buffer.data(sld1 + 93);
    const auto *sld1_95 = buffer.data(sld1 + 95);

    const auto *slf_86 = buffer.data(slf + 86);
    const auto *slf_87 = buffer.data(slf + 87);
    const auto *slf_88 = buffer.data(slf + 88);
    const auto *slf_89 = buffer.data(slf + 89);
    const auto *slf_90 = buffer.data(slf + 90);
    const auto *slf_92 = buffer.data(slf + 92);
    const auto *slf_93 = buffer.data(slf + 93);
    const auto *slf_95 = buffer.data(slf + 95);
    const auto *slf_96 = buffer.data(slf + 96);
    const auto *slf_97 = buffer.data(slf + 97);
    const auto *slf_98 = buffer.data(slf + 98);
    const auto *slf_99 = buffer.data(slf + 99);
    const auto *slf_100 = buffer.data(slf + 100);
    const auto *slf_102 = buffer.data(slf + 102);
    const auto *slf_103 = buffer.data(slf + 103);
    const auto *slf_105 = buffer.data(slf + 105);
    const auto *slf_106 = buffer.data(slf + 106);
    const auto *slf_107 = buffer.data(slf + 107);
    const auto *slf_108 = buffer.data(slf + 108);
    const auto *slf_109 = buffer.data(slf + 109);
    const auto *slf_110 = buffer.data(slf + 110);
    const auto *slf_112 = buffer.data(slf + 112);
    const auto *slf_115 = buffer.data(slf + 115);
    const auto *slf_116 = buffer.data(slf + 116);
    const auto *slf_117 = buffer.data(slf + 117);
    const auto *slf_118 = buffer.data(slf + 118);
    const auto *slf_119 = buffer.data(slf + 119);
    const auto *slf_120 = buffer.data(slf + 120);
    const auto *slf_122 = buffer.data(slf + 122);
    const auto *slf_123 = buffer.data(slf + 123);
    const auto *slf_125 = buffer.data(slf + 125);
    const auto *slf_126 = buffer.data(slf + 126);
    const auto *slf_127 = buffer.data(slf + 127);
    const auto *slf_128 = buffer.data(slf + 128);
    const auto *slf_129 = buffer.data(slf + 129);
    const auto *slf_130 = buffer.data(slf + 130);
    const auto *slf_132 = buffer.data(slf + 132);
    const auto *slf_136 = buffer.data(slf + 136);
    const auto *slf_137 = buffer.data(slf + 137);
    const auto *slf_138 = buffer.data(slf + 138);
    const auto *slf_139 = buffer.data(slf + 139);
    const auto *slf_140 = buffer.data(slf + 140);
    const auto *slf_142 = buffer.data(slf + 142);
    const auto *slf_143 = buffer.data(slf + 143);
    const auto *slf_145 = buffer.data(slf + 145);
    const auto *slf_146 = buffer.data(slf + 146);
    const auto *slf_147 = buffer.data(slf + 147);
    const auto *slf_148 = buffer.data(slf + 148);
    const auto *slf_149 = buffer.data(slf + 149);
    const auto *slf_150 = buffer.data(slf + 150);
    const auto *slf_152 = buffer.data(slf + 152);
    const auto *slf_153 = buffer.data(slf + 153);
    const auto *slf_155 = buffer.data(slf + 155);
    const auto *slf_156 = buffer.data(slf + 156);
    const auto *slf_157 = buffer.data(slf + 157);
    const auto *slf_158 = buffer.data(slf + 158);
    const auto *slf_159 = buffer.data(slf + 159);
    const auto *slf_160 = buffer.data(slf + 160);
    const auto *slf_162 = buffer.data(slf + 162);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pc_x, pc_y, skg0_80, skf_86, \
                         skf_87, skf_88, skg1_80, slf_86, slf_87, \
                         slf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_y[k] * skg0_80[k]
                   - f_6 * pc_y[k] * skg1_80[k];

        t_126[k] = f_11 * skf_86[k]
                   + f_3 * pc_x[k] * slf_86[k];

        t_127[k] = f_11 * skf_87[k]
                   + f_3 * pc_x[k] * slf_87[k];

        t_128[k] = f_11 * skf_88[k]
                   + f_3 * pc_x[k] * slf_88[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pc_x, pc_y, pc_z, skf_46, skf_56, skf_89, \
                         sld0_51, sld1_51, slf_86, slf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_11 * skf_89[k]
                   + f_3 * pc_x[k] * slf_89[k];

        t_130[k] = f_7 * skf_56[k]
                   + f_1 * sld0_51[k]
                   - f_2 * sld1_51[k]
                   + f_3 * pc_y[k] * slf_86[k];

        t_131[k] = f_8 * skf_46[k]
                   + f_3 * pc_z[k] * slf_86[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pb_y, pc_y, skg0_89, skf_58, skf_59, skg1_89, \
                         sld0_53, sld1_53, slf_88, slf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_7 * skf_58[k]
                   + f_4 * sld0_53[k]
                   - f_5 * sld1_53[k]
                   + f_3 * pc_y[k] * slf_88[k];

        t_133[k] = f_7 * skf_59[k]
                   + f_3 * pc_y[k] * slf_89[k];

        t_134[k] = pb_y[k] * skg0_89[k]
                   - f_6 * pc_y[k] * skg1_89[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pc_x, pc_y, pc_z, skf_50, skf_90, skf_93, \
                         sld0_54, sld0_57, sld1_54, sld1_57, slf_90, \
                         slf_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_11 * skf_90[k]
                   + f_1 * sld0_54[k]
                   - f_2 * sld1_54[k]
                   + f_3 * pc_x[k] * slf_90[k];

        t_136[k] = f_3 * pc_y[k] * slf_90[k];

        t_137[k] = f_12 * skf_50[k]
                   + f_3 * pc_z[k] * slf_90[k];

        t_138[k] = f_11 * skf_93[k]
                   + f_4 * sld0_57[k]
                   - f_5 * sld1_57[k]
                   + f_3 * pc_x[k] * slf_93[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pc_x, pc_y, skf_95, skf_96, skf_97, \
                         sld0_59, sld1_59, slf_92, slf_95, slf_96, \
                         slf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_3 * pc_y[k] * slf_92[k];

        t_140[k] = f_11 * skf_95[k]
                   + f_4 * sld0_59[k]
                   - f_5 * sld1_59[k]
                   + f_3 * pc_x[k] * slf_95[k];

        t_141[k] = f_11 * skf_96[k]
                   + f_3 * pc_x[k] * slf_96[k];

        t_142[k] = f_11 * skf_97[k]
                   + f_3 * pc_x[k] * slf_97[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_x, pc_y, pc_z, skf_56, skf_98, skf_99, \
                         sld0_57, sld1_57, slf_96, slf_98, slf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_11 * skf_98[k]
                   + f_3 * pc_x[k] * slf_98[k];

        t_144[k] = f_11 * skf_99[k]
                   + f_3 * pc_x[k] * slf_99[k];

        t_145[k] = f_1 * sld0_57[k]
                   - f_2 * sld1_57[k]
                   + f_3 * pc_y[k] * slf_96[k];

        t_146[k] = f_12 * skf_56[k]
                   + f_3 * pc_z[k] * slf_96[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, skf_59, skf_100, \
                         sld0_59, sld0_60, sld1_59, sld1_60, slf_98, slf_99, \
                         slf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * sld0_59[k]
                   - f_5 * sld1_59[k]
                   + f_3 * pc_y[k] * slf_98[k];

        t_148[k] = f_3 * pc_y[k] * slf_99[k];

        t_149[k] = f_12 * skf_59[k]
                   + f_1 * sld0_59[k]
                   - f_2 * sld1_59[k]
                   + f_3 * pc_z[k] * slf_99[k];

        t_150[k] = f_13 * skf_100[k]
                   + f_1 * sld0_60[k]
                   - f_2 * sld1_60[k]
                   + f_3 * pc_x[k] * slf_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, pc_z, skf_60, skf_62, \
                         skf_103, sld0_63, sld1_63, slf_100, slf_102, \
                         slf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_13 * skf_60[k]
                   + f_3 * pc_y[k] * slf_100[k];

        t_152[k] = f_3 * pc_z[k] * slf_100[k];

        t_153[k] = f_13 * skf_103[k]
                   + f_4 * sld0_63[k]
                   - f_5 * sld1_63[k]
                   + f_3 * pc_x[k] * slf_103[k];

        t_154[k] = f_13 * skf_62[k]
                   + f_3 * pc_y[k] * slf_102[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, skf_105, skf_106, skf_107, skf_108, \
                         sld0_65, sld1_65, slf_105, slf_106, slf_107, \
                         slf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_13 * skf_105[k]
                   + f_4 * sld0_65[k]
                   - f_5 * sld1_65[k]
                   + f_3 * pc_x[k] * slf_105[k];

        t_156[k] = f_13 * skf_106[k]
                   + f_3 * pc_x[k] * slf_106[k];

        t_157[k] = f_13 * skf_107[k]
                   + f_3 * pc_x[k] * slf_107[k];

        t_158[k] = f_13 * skf_108[k]
                   + f_3 * pc_x[k] * slf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pc_x, pc_y, pc_z, skf_66, skf_109, sld0_63, \
                         sld1_63, slf_106, slf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_13 * skf_109[k]
                   + f_3 * pc_x[k] * slf_109[k];

        t_160[k] = f_13 * skf_66[k]
                   + f_1 * sld0_63[k]
                   - f_2 * sld1_63[k]
                   + f_3 * pc_y[k] * slf_106[k];

        t_161[k] = f_3 * pc_z[k] * slf_106[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pb_z, pc_y, pc_z, skg0_90, skf_68, \
                         skf_69, skg1_90, sld0_65, sld1_65, slf_108, \
                         slf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_13 * skf_68[k]
                   + f_4 * sld0_65[k]
                   - f_5 * sld1_65[k]
                   + f_3 * pc_y[k] * slf_108[k];

        t_163[k] = f_13 * skf_69[k]
                   + f_3 * pc_y[k] * slf_109[k];

        t_164[k] = f_1 * sld0_65[k]
                   - f_2 * sld1_65[k]
                   + f_3 * pc_z[k] * slf_109[k];

        t_165[k] = pb_z[k] * skg0_90[k]
                   - f_6 * pc_z[k] * skg1_90[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, skg0_93, skf_60, \
                         skf_70, skf_72, skg1_93, slf_110, slf_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_12 * skf_70[k]
                   + f_3 * pc_y[k] * slf_110[k];

        t_167[k] = f_7 * skf_60[k]
                   + f_3 * pc_z[k] * slf_110[k];

        t_168[k] = pb_z[k] * skg0_93[k]
                   - f_6 * pc_z[k] * skg1_93[k];

        t_169[k] = f_12 * skf_72[k]
                   + f_3 * pc_y[k] * slf_112[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, skf_115, skf_116, skf_117, skf_118, \
                         sld0_71, sld1_71, slf_115, slf_116, slf_117, \
                         slf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_13 * skf_115[k]
                   + f_4 * sld0_71[k]
                   - f_5 * sld1_71[k]
                   + f_3 * pc_x[k] * slf_115[k];

        t_171[k] = f_13 * skf_116[k]
                   + f_3 * pc_x[k] * slf_116[k];

        t_172[k] = f_13 * skf_117[k]
                   + f_3 * pc_x[k] * slf_117[k];

        t_173[k] = f_13 * skf_118[k]
                   + f_3 * pc_x[k] * slf_118[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pb_z, pc_x, pc_z, skg0_100, skf_66, skf_119, \
                         skg1_100, slf_116, slf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * skf_119[k]
                   + f_3 * pc_x[k] * slf_119[k];

        t_175[k] = pb_z[k] * skg0_100[k]
                   - f_6 * pc_z[k] * skg1_100[k];

        t_176[k] = f_7 * skf_66[k]
                   + f_3 * pc_z[k] * slf_116[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pc_y, pc_z, skf_69, skf_78, skf_79, sld0_71, \
                         sld1_71, slf_118, slf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_12 * skf_78[k]
                   + f_4 * sld0_71[k]
                   - f_5 * sld1_71[k]
                   + f_3 * pc_y[k] * slf_118[k];

        t_178[k] = f_12 * skf_79[k]
                   + f_3 * pc_y[k] * slf_119[k];

        t_179[k] = f_7 * skf_69[k]
                   + f_1 * sld0_71[k]
                   - f_2 * sld1_71[k]
                   + f_3 * pc_z[k] * slf_119[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pc_x, pc_y, pc_z, skf_70, skf_80, skf_120, \
                         sld0_72, sld1_72, slf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_13 * skf_120[k]
                   + f_1 * sld0_72[k]
                   - f_2 * sld1_72[k]
                   + f_3 * pc_x[k] * slf_120[k];

        t_181[k] = f_8 * skf_80[k]
                   + f_3 * pc_y[k] * slf_120[k];

        t_182[k] = f_8 * skf_70[k]
                   + f_3 * pc_z[k] * slf_120[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_x, pc_y, skf_82, skf_123, skf_125, sld0_75, \
                         sld0_77, sld1_75, sld1_77, slf_122, slf_123, \
                         slf_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_13 * skf_123[k]
                   + f_4 * sld0_75[k]
                   - f_5 * sld1_75[k]
                   + f_3 * pc_x[k] * slf_123[k];

        t_184[k] = f_8 * skf_82[k]
                   + f_3 * pc_y[k] * slf_122[k];

        t_185[k] = f_13 * skf_125[k]
                   + f_4 * sld0_77[k]
                   - f_5 * sld1_77[k]
                   + f_3 * pc_x[k] * slf_125[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pc_x, skf_126, skf_127, skf_128, skf_129, \
                         slf_126, slf_127, slf_128, slf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_13 * skf_126[k]
                   + f_3 * pc_x[k] * slf_126[k];

        t_187[k] = f_13 * skf_127[k]
                   + f_3 * pc_x[k] * slf_127[k];

        t_188[k] = f_13 * skf_128[k]
                   + f_3 * pc_x[k] * slf_128[k];

        t_189[k] = f_13 * skf_129[k]
                   + f_3 * pc_x[k] * slf_129[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_y, pc_z, skf_76, skf_86, skf_88, sld0_75, \
                         sld0_77, sld1_75, sld1_77, slf_126, slf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_8 * skf_86[k]
                   + f_1 * sld0_75[k]
                   - f_2 * sld1_75[k]
                   + f_3 * pc_y[k] * slf_126[k];

        t_191[k] = f_8 * skf_76[k]
                   + f_3 * pc_z[k] * slf_126[k];

        t_192[k] = f_8 * skf_88[k]
                   + f_4 * sld0_77[k]
                   - f_5 * sld1_77[k]
                   + f_3 * pc_y[k] * slf_128[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pb_y, pc_y, pc_z, skg0_135, skf_79, \
                         skf_89, skf_90, skg1_135, sld0_77, sld1_77, slf_129, \
                         slf_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * skf_89[k]
                   + f_3 * pc_y[k] * slf_129[k];

        t_194[k] = f_8 * skf_79[k]
                   + f_1 * sld0_77[k]
                   - f_2 * sld1_77[k]
                   + f_3 * pc_z[k] * slf_129[k];

        t_195[k] = pb_y[k] * skg0_135[k]
                   - f_6 * pc_y[k] * skg1_135[k];

        t_196[k] = f_7 * skf_90[k]
                   + f_3 * pc_y[k] * slf_130[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pb_y, pc_y, pc_z, skg0_138, skg0_140, \
                         skf_80, skf_91, skf_92, skg1_138, skg1_140, slf_130, \
                         slf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * skf_80[k]
                   + f_3 * pc_z[k] * slf_130[k];

        t_198[k] = pb_y[k] * skg0_138[k]
                   + f_8 * skf_91[k]
                   - f_6 * pc_y[k] * skg1_138[k];

        t_199[k] = f_7 * skf_92[k]
                   + f_3 * pc_y[k] * slf_132[k];

        t_200[k] = pb_y[k] * skg0_140[k]
                   - f_6 * pc_y[k] * skg1_140[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pc_x, skf_136, skf_137, skf_138, skf_139, \
                         slf_136, slf_137, slf_138, slf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_13 * skf_136[k]
                   + f_3 * pc_x[k] * slf_136[k];

        t_202[k] = f_13 * skf_137[k]
                   + f_3 * pc_x[k] * slf_137[k];

        t_203[k] = f_13 * skf_138[k]
                   + f_3 * pc_x[k] * slf_138[k];

        t_204[k] = f_13 * skf_139[k]
                   + f_3 * pc_x[k] * slf_139[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pc_y, pc_z, skf_86, skf_96, skf_98, sld0_81, \
                         sld0_83, sld1_81, sld1_83, slf_136, slf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_7 * skf_96[k]
                   + f_1 * sld0_81[k]
                   - f_2 * sld1_81[k]
                   + f_3 * pc_y[k] * slf_136[k];

        t_206[k] = f_12 * skf_86[k]
                   + f_3 * pc_z[k] * slf_136[k];

        t_207[k] = f_7 * skf_98[k]
                   + f_4 * sld0_83[k]
                   - f_5 * sld1_83[k]
                   + f_3 * pc_y[k] * slf_138[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_y, pc_x, pc_y, skg0_149, skf_99, \
                         skf_140, skg1_149, sld0_84, sld1_84, slf_139, \
                         slf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_7 * skf_99[k]
                   + f_3 * pc_y[k] * slf_139[k];

        t_209[k] = pb_y[k] * skg0_149[k]
                   - f_6 * pc_y[k] * skg1_149[k];

        t_210[k] = f_13 * skf_140[k]
                   + f_1 * sld0_84[k]
                   - f_2 * sld1_84[k]
                   + f_3 * pc_x[k] * slf_140[k];

        t_211[k] = f_3 * pc_y[k] * slf_140[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pc_x, pc_y, pc_z, skf_90, skf_143, sld0_87, \
                         sld1_87, slf_140, slf_142, slf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_13 * skf_90[k]
                   + f_3 * pc_z[k] * slf_140[k];

        t_213[k] = f_13 * skf_143[k]
                   + f_4 * sld0_87[k]
                   - f_5 * sld1_87[k]
                   + f_3 * pc_x[k] * slf_143[k];

        t_214[k] = f_3 * pc_y[k] * slf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, skf_145, skf_146, skf_147, skf_148, \
                         sld0_89, sld1_89, slf_145, slf_146, slf_147, \
                         slf_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_13 * skf_145[k]
                   + f_4 * sld0_89[k]
                   - f_5 * sld1_89[k]
                   + f_3 * pc_x[k] * slf_145[k];

        t_216[k] = f_13 * skf_146[k]
                   + f_3 * pc_x[k] * slf_146[k];

        t_217[k] = f_13 * skf_147[k]
                   + f_3 * pc_x[k] * slf_147[k];

        t_218[k] = f_13 * skf_148[k]
                   + f_3 * pc_x[k] * slf_148[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pc_x, pc_y, pc_z, skf_96, skf_149, \
                         sld0_87, sld0_89, sld1_87, sld1_89, slf_146, slf_148, \
                         slf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_13 * skf_149[k]
                   + f_3 * pc_x[k] * slf_149[k];

        t_220[k] = f_1 * sld0_87[k]
                   - f_2 * sld1_87[k]
                   + f_3 * pc_y[k] * slf_146[k];

        t_221[k] = f_13 * skf_96[k]
                   + f_3 * pc_z[k] * slf_146[k];

        t_222[k] = f_4 * sld0_89[k]
                   - f_5 * sld1_89[k]
                   + f_3 * pc_y[k] * slf_148[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pc_x, pc_y, pc_z, skf_99, skf_100, \
                         skf_150, sld0_89, sld0_90, sld1_89, sld1_90, slf_149, \
                         slf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_y[k] * slf_149[k];

        t_224[k] = f_13 * skf_99[k]
                   + f_1 * sld0_89[k]
                   - f_2 * sld1_89[k]
                   + f_3 * pc_z[k] * slf_149[k];

        t_225[k] = f_12 * skf_150[k]
                   + f_1 * sld0_90[k]
                   - f_2 * sld1_90[k]
                   + f_3 * pc_x[k] * slf_150[k];

        t_226[k] = f_11 * skf_100[k]
                   + f_3 * pc_y[k] * slf_150[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_x, pc_y, pc_z, skf_102, skf_153, sld0_93, \
                         sld1_93, slf_150, slf_152, slf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_3 * pc_z[k] * slf_150[k];

        t_228[k] = f_12 * skf_153[k]
                   + f_4 * sld0_93[k]
                   - f_5 * sld1_93[k]
                   + f_3 * pc_x[k] * slf_153[k];

        t_229[k] = f_11 * skf_102[k]
                   + f_3 * pc_y[k] * slf_152[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, skf_155, skf_156, skf_157, skf_158, \
                         sld0_95, sld1_95, slf_155, slf_156, slf_157, \
                         slf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_12 * skf_155[k]
                   + f_4 * sld0_95[k]
                   - f_5 * sld1_95[k]
                   + f_3 * pc_x[k] * slf_155[k];

        t_231[k] = f_12 * skf_156[k]
                   + f_3 * pc_x[k] * slf_156[k];

        t_232[k] = f_12 * skf_157[k]
                   + f_3 * pc_x[k] * slf_157[k];

        t_233[k] = f_12 * skf_158[k]
                   + f_3 * pc_x[k] * slf_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, pc_y, pc_z, skf_106, skf_159, sld0_93, \
                         sld1_93, slf_156, slf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_12 * skf_159[k]
                   + f_3 * pc_x[k] * slf_159[k];

        t_235[k] = f_11 * skf_106[k]
                   + f_1 * sld0_93[k]
                   - f_2 * sld1_93[k]
                   + f_3 * pc_y[k] * slf_156[k];

        t_236[k] = f_3 * pc_z[k] * slf_156[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pb_z, pc_y, pc_z, skg0_150, skf_108, \
                         skf_109, skg1_150, sld0_95, sld1_95, slf_158, \
                         slf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_11 * skf_108[k]
                   + f_4 * sld0_95[k]
                   - f_5 * sld1_95[k]
                   + f_3 * pc_y[k] * slf_158[k];

        t_238[k] = f_11 * skf_109[k]
                   + f_3 * pc_y[k] * slf_159[k];

        t_239[k] = f_1 * sld0_95[k]
                   - f_2 * sld1_95[k]
                   + f_3 * pc_z[k] * slf_159[k];

        t_240[k] = pb_z[k] * skg0_150[k]
                   - f_6 * pc_z[k] * skg1_150[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pb_z, pc_y, pc_z, skg0_153, skf_100, \
                         skf_110, skf_112, skg1_153, slf_160, slf_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_13 * skf_110[k]
                   + f_3 * pc_y[k] * slf_160[k];

        t_242[k] = f_7 * skf_100[k]
                   + f_3 * pc_z[k] * slf_160[k];

        t_243[k] = pb_z[k] * skg0_153[k]
                   - f_6 * pc_z[k] * skg1_153[k];

        t_244[k] = f_13 * skf_112[k]
                   + f_3 * pc_y[k] * slf_162[k];
    }
}

static auto
compute_prim_slg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skg0,
                                                          const size_t skf, const size_t skg1,
                                                          const size_t sld0, const size_t sld1,
                                                          const size_t slf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 3.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.0 / q;

    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);
    auto *t_248 = buffer.data(target + 248);
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);
    auto *t_256 = buffer.data(target + 256);
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);
    auto *t_259 = buffer.data(target + 259);
    auto *t_260 = buffer.data(target + 260);
    auto *t_261 = buffer.data(target + 261);
    auto *t_262 = buffer.data(target + 262);
    auto *t_263 = buffer.data(target + 263);
    auto *t_264 = buffer.data(target + 264);
    auto *t_265 = buffer.data(target + 265);
    auto *t_266 = buffer.data(target + 266);
    auto *t_267 = buffer.data(target + 267);
    auto *t_268 = buffer.data(target + 268);
    auto *t_269 = buffer.data(target + 269);
    auto *t_270 = buffer.data(target + 270);
    auto *t_271 = buffer.data(target + 271);
    auto *t_272 = buffer.data(target + 272);
    auto *t_273 = buffer.data(target + 273);
    auto *t_274 = buffer.data(target + 274);
    auto *t_275 = buffer.data(target + 275);
    auto *t_276 = buffer.data(target + 276);
    auto *t_277 = buffer.data(target + 277);
    auto *t_278 = buffer.data(target + 278);
    auto *t_279 = buffer.data(target + 279);
    auto *t_280 = buffer.data(target + 280);
    auto *t_281 = buffer.data(target + 281);
    auto *t_282 = buffer.data(target + 282);
    auto *t_283 = buffer.data(target + 283);
    auto *t_284 = buffer.data(target + 284);
    auto *t_285 = buffer.data(target + 285);
    auto *t_286 = buffer.data(target + 286);
    auto *t_287 = buffer.data(target + 287);
    auto *t_288 = buffer.data(target + 288);
    auto *t_289 = buffer.data(target + 289);
    auto *t_290 = buffer.data(target + 290);
    auto *t_291 = buffer.data(target + 291);
    auto *t_292 = buffer.data(target + 292);
    auto *t_293 = buffer.data(target + 293);
    auto *t_294 = buffer.data(target + 294);
    auto *t_295 = buffer.data(target + 295);
    auto *t_296 = buffer.data(target + 296);
    auto *t_297 = buffer.data(target + 297);
    auto *t_298 = buffer.data(target + 298);
    auto *t_299 = buffer.data(target + 299);
    auto *t_300 = buffer.data(target + 300);
    auto *t_301 = buffer.data(target + 301);
    auto *t_302 = buffer.data(target + 302);
    auto *t_303 = buffer.data(target + 303);
    auto *t_304 = buffer.data(target + 304);
    auto *t_305 = buffer.data(target + 305);
    auto *t_306 = buffer.data(target + 306);
    auto *t_307 = buffer.data(target + 307);
    auto *t_308 = buffer.data(target + 308);
    auto *t_309 = buffer.data(target + 309);
    auto *t_310 = buffer.data(target + 310);
    auto *t_311 = buffer.data(target + 311);
    auto *t_312 = buffer.data(target + 312);
    auto *t_313 = buffer.data(target + 313);
    auto *t_314 = buffer.data(target + 314);
    auto *t_315 = buffer.data(target + 315);
    auto *t_316 = buffer.data(target + 316);
    auto *t_317 = buffer.data(target + 317);
    auto *t_318 = buffer.data(target + 318);
    auto *t_319 = buffer.data(target + 319);
    auto *t_320 = buffer.data(target + 320);
    auto *t_321 = buffer.data(target + 321);
    auto *t_322 = buffer.data(target + 322);
    auto *t_323 = buffer.data(target + 323);
    auto *t_324 = buffer.data(target + 324);
    auto *t_325 = buffer.data(target + 325);
    auto *t_326 = buffer.data(target + 326);
    auto *t_327 = buffer.data(target + 327);
    auto *t_328 = buffer.data(target + 328);
    auto *t_329 = buffer.data(target + 329);
    auto *t_330 = buffer.data(target + 330);
    auto *t_331 = buffer.data(target + 331);
    auto *t_332 = buffer.data(target + 332);
    auto *t_333 = buffer.data(target + 333);
    auto *t_334 = buffer.data(target + 334);
    auto *t_335 = buffer.data(target + 335);
    auto *t_336 = buffer.data(target + 336);
    auto *t_337 = buffer.data(target + 337);
    auto *t_338 = buffer.data(target + 338);
    auto *t_339 = buffer.data(target + 339);
    auto *t_340 = buffer.data(target + 340);
    auto *t_341 = buffer.data(target + 341);
    auto *t_342 = buffer.data(target + 342);
    auto *t_343 = buffer.data(target + 343);
    auto *t_344 = buffer.data(target + 344);
    auto *t_345 = buffer.data(target + 345);
    auto *t_346 = buffer.data(target + 346);
    auto *t_347 = buffer.data(target + 347);
    auto *t_348 = buffer.data(target + 348);
    auto *t_349 = buffer.data(target + 349);
    auto *t_350 = buffer.data(target + 350);
    auto *t_351 = buffer.data(target + 351);
    auto *t_352 = buffer.data(target + 352);
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);
    auto *t_357 = buffer.data(target + 357);
    auto *t_358 = buffer.data(target + 358);
    auto *t_359 = buffer.data(target + 359);
    auto *t_360 = buffer.data(target + 360);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skg0_160 = buffer.data(skg0 + 160);
    const auto *skg0_210 = buffer.data(skg0 + 210);
    const auto *skg0_213 = buffer.data(skg0 + 213);
    const auto *skg0_215 = buffer.data(skg0 + 215);
    const auto *skg0_224 = buffer.data(skg0 + 224);
    const auto *skg0_225 = buffer.data(skg0 + 225);
    const auto *skg0_228 = buffer.data(skg0 + 228);
    const auto *skg0_235 = buffer.data(skg0 + 235);

    const auto *skf_106 = buffer.data(skf + 106);
    const auto *skf_109 = buffer.data(skf + 109);
    const auto *skf_110 = buffer.data(skf + 110);
    const auto *skf_116 = buffer.data(skf + 116);
    const auto *skf_118 = buffer.data(skf + 118);
    const auto *skf_119 = buffer.data(skf + 119);
    const auto *skf_120 = buffer.data(skf + 120);
    const auto *skf_122 = buffer.data(skf + 122);
    const auto *skf_126 = buffer.data(skf + 126);
    const auto *skf_128 = buffer.data(skf + 128);
    const auto *skf_129 = buffer.data(skf + 129);
    const auto *skf_130 = buffer.data(skf + 130);
    const auto *skf_132 = buffer.data(skf + 132);
    const auto *skf_136 = buffer.data(skf + 136);
    const auto *skf_138 = buffer.data(skf + 138);
    const auto *skf_139 = buffer.data(skf + 139);
    const auto *skf_140 = buffer.data(skf + 140);
    const auto *skf_141 = buffer.data(skf + 141);
    const auto *skf_142 = buffer.data(skf + 142);
    const auto *skf_146 = buffer.data(skf + 146);
    const auto *skf_148 = buffer.data(skf + 148);
    const auto *skf_149 = buffer.data(skf + 149);
    const auto *skf_150 = buffer.data(skf + 150);
    const auto *skf_152 = buffer.data(skf + 152);
    const auto *skf_156 = buffer.data(skf + 156);
    const auto *skf_158 = buffer.data(skf + 158);
    const auto *skf_159 = buffer.data(skf + 159);
    const auto *skf_160 = buffer.data(skf + 160);
    const auto *skf_162 = buffer.data(skf + 162);
    const auto *skf_165 = buffer.data(skf + 165);
    const auto *skf_166 = buffer.data(skf + 166);
    const auto *skf_167 = buffer.data(skf + 167);
    const auto *skf_168 = buffer.data(skf + 168);
    const auto *skf_169 = buffer.data(skf + 169);
    const auto *skf_170 = buffer.data(skf + 170);
    const auto *skf_172 = buffer.data(skf + 172);
    const auto *skf_173 = buffer.data(skf + 173);
    const auto *skf_175 = buffer.data(skf + 175);
    const auto *skf_176 = buffer.data(skf + 176);
    const auto *skf_177 = buffer.data(skf + 177);
    const auto *skf_178 = buffer.data(skf + 178);
    const auto *skf_179 = buffer.data(skf + 179);
    const auto *skf_180 = buffer.data(skf + 180);
    const auto *skf_183 = buffer.data(skf + 183);
    const auto *skf_185 = buffer.data(skf + 185);
    const auto *skf_186 = buffer.data(skf + 186);
    const auto *skf_187 = buffer.data(skf + 187);
    const auto *skf_188 = buffer.data(skf + 188);
    const auto *skf_189 = buffer.data(skf + 189);
    const auto *skf_196 = buffer.data(skf + 196);
    const auto *skf_197 = buffer.data(skf + 197);
    const auto *skf_198 = buffer.data(skf + 198);
    const auto *skf_199 = buffer.data(skf + 199);
    const auto *skf_200 = buffer.data(skf + 200);
    const auto *skf_203 = buffer.data(skf + 203);
    const auto *skf_205 = buffer.data(skf + 205);
    const auto *skf_206 = buffer.data(skf + 206);
    const auto *skf_207 = buffer.data(skf + 207);
    const auto *skf_208 = buffer.data(skf + 208);
    const auto *skf_209 = buffer.data(skf + 209);
    const auto *skf_210 = buffer.data(skf + 210);
    const auto *skf_213 = buffer.data(skf + 213);
    const auto *skf_215 = buffer.data(skf + 215);
    const auto *skf_216 = buffer.data(skf + 216);
    const auto *skf_217 = buffer.data(skf + 217);
    const auto *skf_218 = buffer.data(skf + 218);
    const auto *skf_219 = buffer.data(skf + 219);
    const auto *skf_225 = buffer.data(skf + 225);
    const auto *skf_226 = buffer.data(skf + 226);
    const auto *skf_227 = buffer.data(skf + 227);
    const auto *skf_228 = buffer.data(skf + 228);
    const auto *skf_229 = buffer.data(skf + 229);
    const auto *skf_230 = buffer.data(skf + 230);
    const auto *skf_233 = buffer.data(skf + 233);
    const auto *skf_235 = buffer.data(skf + 235);
    const auto *skf_236 = buffer.data(skf + 236);
    const auto *skf_237 = buffer.data(skf + 237);
    const auto *skf_238 = buffer.data(skf + 238);
    const auto *skf_239 = buffer.data(skf + 239);
    const auto *skf_240 = buffer.data(skf + 240);

    const auto *skg1_160 = buffer.data(skg1 + 160);
    const auto *skg1_210 = buffer.data(skg1 + 210);
    const auto *skg1_213 = buffer.data(skg1 + 213);
    const auto *skg1_215 = buffer.data(skg1 + 215);
    const auto *skg1_224 = buffer.data(skg1 + 224);
    const auto *skg1_225 = buffer.data(skg1 + 225);
    const auto *skg1_228 = buffer.data(skg1 + 228);
    const auto *skg1_235 = buffer.data(skg1 + 235);

    const auto *sld0_101 = buffer.data(sld0 + 101);
    const auto *sld0_102 = buffer.data(sld0 + 102);
    const auto *sld0_105 = buffer.data(sld0 + 105);
    const auto *sld0_107 = buffer.data(sld0 + 107);
    const auto *sld0_108 = buffer.data(sld0 + 108);
    const auto *sld0_111 = buffer.data(sld0 + 111);
    const auto *sld0_113 = buffer.data(sld0 + 113);
    const auto *sld0_117 = buffer.data(sld0 + 117);
    const auto *sld0_119 = buffer.data(sld0 + 119);
    const auto *sld0_120 = buffer.data(sld0 + 120);
    const auto *sld0_123 = buffer.data(sld0 + 123);
    const auto *sld0_125 = buffer.data(sld0 + 125);
    const auto *sld0_126 = buffer.data(sld0 + 126);
    const auto *sld0_129 = buffer.data(sld0 + 129);
    const auto *sld0_131 = buffer.data(sld0 + 131);
    const auto *sld0_137 = buffer.data(sld0 + 137);
    const auto *sld0_138 = buffer.data(sld0 + 138);
    const auto *sld0_141 = buffer.data(sld0 + 141);
    const auto *sld0_143 = buffer.data(sld0 + 143);
    const auto *sld0_144 = buffer.data(sld0 + 144);

    const auto *sld1_101 = buffer.data(sld1 + 101);
    const auto *sld1_102 = buffer.data(sld1 + 102);
    const auto *sld1_105 = buffer.data(sld1 + 105);
    const auto *sld1_107 = buffer.data(sld1 + 107);
    const auto *sld1_108 = buffer.data(sld1 + 108);
    const auto *sld1_111 = buffer.data(sld1 + 111);
    const auto *sld1_113 = buffer.data(sld1 + 113);
    const auto *sld1_117 = buffer.data(sld1 + 117);
    const auto *sld1_119 = buffer.data(sld1 + 119);
    const auto *sld1_120 = buffer.data(sld1 + 120);
    const auto *sld1_123 = buffer.data(sld1 + 123);
    const auto *sld1_125 = buffer.data(sld1 + 125);
    const auto *sld1_126 = buffer.data(sld1 + 126);
    const auto *sld1_129 = buffer.data(sld1 + 129);
    const auto *sld1_131 = buffer.data(sld1 + 131);
    const auto *sld1_137 = buffer.data(sld1 + 137);
    const auto *sld1_138 = buffer.data(sld1 + 138);
    const auto *sld1_141 = buffer.data(sld1 + 141);
    const auto *sld1_143 = buffer.data(sld1 + 143);
    const auto *sld1_144 = buffer.data(sld1 + 144);

    const auto *slf_165 = buffer.data(slf + 165);
    const auto *slf_166 = buffer.data(slf + 166);
    const auto *slf_167 = buffer.data(slf + 167);
    const auto *slf_168 = buffer.data(slf + 168);
    const auto *slf_169 = buffer.data(slf + 169);
    const auto *slf_170 = buffer.data(slf + 170);
    const auto *slf_172 = buffer.data(slf + 172);
    const auto *slf_173 = buffer.data(slf + 173);
    const auto *slf_175 = buffer.data(slf + 175);
    const auto *slf_176 = buffer.data(slf + 176);
    const auto *slf_177 = buffer.data(slf + 177);
    const auto *slf_178 = buffer.data(slf + 178);
    const auto *slf_179 = buffer.data(slf + 179);
    const auto *slf_180 = buffer.data(slf + 180);
    const auto *slf_182 = buffer.data(slf + 182);
    const auto *slf_183 = buffer.data(slf + 183);
    const auto *slf_185 = buffer.data(slf + 185);
    const auto *slf_186 = buffer.data(slf + 186);
    const auto *slf_187 = buffer.data(slf + 187);
    const auto *slf_188 = buffer.data(slf + 188);
    const auto *slf_189 = buffer.data(slf + 189);
    const auto *slf_190 = buffer.data(slf + 190);
    const auto *slf_192 = buffer.data(slf + 192);
    const auto *slf_196 = buffer.data(slf + 196);
    const auto *slf_197 = buffer.data(slf + 197);
    const auto *slf_198 = buffer.data(slf + 198);
    const auto *slf_199 = buffer.data(slf + 199);
    const auto *slf_200 = buffer.data(slf + 200);
    const auto *slf_202 = buffer.data(slf + 202);
    const auto *slf_203 = buffer.data(slf + 203);
    const auto *slf_205 = buffer.data(slf + 205);
    const auto *slf_206 = buffer.data(slf + 206);
    const auto *slf_207 = buffer.data(slf + 207);
    const auto *slf_208 = buffer.data(slf + 208);
    const auto *slf_209 = buffer.data(slf + 209);
    const auto *slf_210 = buffer.data(slf + 210);
    const auto *slf_212 = buffer.data(slf + 212);
    const auto *slf_213 = buffer.data(slf + 213);
    const auto *slf_215 = buffer.data(slf + 215);
    const auto *slf_216 = buffer.data(slf + 216);
    const auto *slf_217 = buffer.data(slf + 217);
    const auto *slf_218 = buffer.data(slf + 218);
    const auto *slf_219 = buffer.data(slf + 219);
    const auto *slf_220 = buffer.data(slf + 220);
    const auto *slf_222 = buffer.data(slf + 222);
    const auto *slf_225 = buffer.data(slf + 225);
    const auto *slf_226 = buffer.data(slf + 226);
    const auto *slf_227 = buffer.data(slf + 227);
    const auto *slf_228 = buffer.data(slf + 228);
    const auto *slf_229 = buffer.data(slf + 229);
    const auto *slf_230 = buffer.data(slf + 230);
    const auto *slf_232 = buffer.data(slf + 232);
    const auto *slf_233 = buffer.data(slf + 233);
    const auto *slf_235 = buffer.data(slf + 235);
    const auto *slf_236 = buffer.data(slf + 236);
    const auto *slf_237 = buffer.data(slf + 237);
    const auto *slf_238 = buffer.data(slf + 238);
    const auto *slf_239 = buffer.data(slf + 239);
    const auto *slf_240 = buffer.data(slf + 240);

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pc_x, skf_165, skf_166, skf_167, skf_168, \
                         sld0_101, sld1_101, slf_165, slf_166, slf_167, \
                         slf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_12 * skf_165[k]
                   + f_4 * sld0_101[k]
                   - f_5 * sld1_101[k]
                   + f_3 * pc_x[k] * slf_165[k];

        t_246[k] = f_12 * skf_166[k]
                   + f_3 * pc_x[k] * slf_166[k];

        t_247[k] = f_12 * skf_167[k]
                   + f_3 * pc_x[k] * slf_167[k];

        t_248[k] = f_12 * skf_168[k]
                   + f_3 * pc_x[k] * slf_168[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pb_z, pc_x, pc_z, skg0_160, skf_106, skf_169, \
                         skg1_160, slf_166, slf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_12 * skf_169[k]
                   + f_3 * pc_x[k] * slf_169[k];

        t_250[k] = pb_z[k] * skg0_160[k]
                   - f_6 * pc_z[k] * skg1_160[k];

        t_251[k] = f_7 * skf_106[k]
                   + f_3 * pc_z[k] * slf_166[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, pc_y, pc_z, skf_109, skf_118, skf_119, sld0_101, \
                         sld1_101, slf_168, slf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_13 * skf_118[k]
                   + f_4 * sld0_101[k]
                   - f_5 * sld1_101[k]
                   + f_3 * pc_y[k] * slf_168[k];

        t_253[k] = f_13 * skf_119[k]
                   + f_3 * pc_y[k] * slf_169[k];

        t_254[k] = f_7 * skf_109[k]
                   + f_1 * sld0_101[k]
                   - f_2 * sld1_101[k]
                   + f_3 * pc_z[k] * slf_169[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, pc_z, skf_110, skf_120, skf_170, \
                         sld0_102, sld1_102, slf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_12 * skf_170[k]
                   + f_1 * sld0_102[k]
                   - f_2 * sld1_102[k]
                   + f_3 * pc_x[k] * slf_170[k];

        t_256[k] = f_12 * skf_120[k]
                   + f_3 * pc_y[k] * slf_170[k];

        t_257[k] = f_8 * skf_110[k]
                   + f_3 * pc_z[k] * slf_170[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_x, pc_y, skf_122, skf_173, skf_175, sld0_105, \
                         sld0_107, sld1_105, sld1_107, slf_172, slf_173, \
                         slf_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_12 * skf_173[k]
                   + f_4 * sld0_105[k]
                   - f_5 * sld1_105[k]
                   + f_3 * pc_x[k] * slf_173[k];

        t_259[k] = f_12 * skf_122[k]
                   + f_3 * pc_y[k] * slf_172[k];

        t_260[k] = f_12 * skf_175[k]
                   + f_4 * sld0_107[k]
                   - f_5 * sld1_107[k]
                   + f_3 * pc_x[k] * slf_175[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pc_x, skf_176, skf_177, skf_178, skf_179, \
                         slf_176, slf_177, slf_178, slf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_12 * skf_176[k]
                   + f_3 * pc_x[k] * slf_176[k];

        t_262[k] = f_12 * skf_177[k]
                   + f_3 * pc_x[k] * slf_177[k];

        t_263[k] = f_12 * skf_178[k]
                   + f_3 * pc_x[k] * slf_178[k];

        t_264[k] = f_12 * skf_179[k]
                   + f_3 * pc_x[k] * slf_179[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, pc_y, pc_z, skf_116, skf_126, skf_128, sld0_105, \
                         sld0_107, sld1_105, sld1_107, slf_176, \
                         slf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_12 * skf_126[k]
                   + f_1 * sld0_105[k]
                   - f_2 * sld1_105[k]
                   + f_3 * pc_y[k] * slf_176[k];

        t_266[k] = f_8 * skf_116[k]
                   + f_3 * pc_z[k] * slf_176[k];

        t_267[k] = f_12 * skf_128[k]
                   + f_4 * sld0_107[k]
                   - f_5 * sld1_107[k]
                   + f_3 * pc_y[k] * slf_178[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_x, pc_y, pc_z, skf_119, skf_129, skf_180, \
                         sld0_107, sld0_108, sld1_107, sld1_108, slf_179, \
                         slf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_12 * skf_129[k]
                   + f_3 * pc_y[k] * slf_179[k];

        t_269[k] = f_8 * skf_119[k]
                   + f_1 * sld0_107[k]
                   - f_2 * sld1_107[k]
                   + f_3 * pc_z[k] * slf_179[k];

        t_270[k] = f_12 * skf_180[k]
                   + f_1 * sld0_108[k]
                   - f_2 * sld1_108[k]
                   + f_3 * pc_x[k] * slf_180[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, skf_120, skf_130, \
                         skf_132, skf_183, sld0_111, sld1_111, slf_180, slf_182, \
                         slf_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_8 * skf_130[k]
                   + f_3 * pc_y[k] * slf_180[k];

        t_272[k] = f_12 * skf_120[k]
                   + f_3 * pc_z[k] * slf_180[k];

        t_273[k] = f_12 * skf_183[k]
                   + f_4 * sld0_111[k]
                   - f_5 * sld1_111[k]
                   + f_3 * pc_x[k] * slf_183[k];

        t_274[k] = f_8 * skf_132[k]
                   + f_3 * pc_y[k] * slf_182[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pc_x, skf_185, skf_186, skf_187, skf_188, \
                         sld0_113, sld1_113, slf_185, slf_186, slf_187, \
                         slf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_12 * skf_185[k]
                   + f_4 * sld0_113[k]
                   - f_5 * sld1_113[k]
                   + f_3 * pc_x[k] * slf_185[k];

        t_276[k] = f_12 * skf_186[k]
                   + f_3 * pc_x[k] * slf_186[k];

        t_277[k] = f_12 * skf_187[k]
                   + f_3 * pc_x[k] * slf_187[k];

        t_278[k] = f_12 * skf_188[k]
                   + f_3 * pc_x[k] * slf_188[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pc_x, pc_y, pc_z, skf_126, skf_136, skf_189, \
                         sld0_111, sld1_111, slf_186, slf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_12 * skf_189[k]
                   + f_3 * pc_x[k] * slf_189[k];

        t_280[k] = f_8 * skf_136[k]
                   + f_1 * sld0_111[k]
                   - f_2 * sld1_111[k]
                   + f_3 * pc_y[k] * slf_186[k];

        t_281[k] = f_12 * skf_126[k]
                   + f_3 * pc_z[k] * slf_186[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pb_y, pc_y, pc_z, skg0_210, skf_129, \
                         skf_138, skf_139, skg1_210, sld0_113, sld1_113, slf_188, \
                         slf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_8 * skf_138[k]
                   + f_4 * sld0_113[k]
                   - f_5 * sld1_113[k]
                   + f_3 * pc_y[k] * slf_188[k];

        t_283[k] = f_8 * skf_139[k]
                   + f_3 * pc_y[k] * slf_189[k];

        t_284[k] = f_12 * skf_129[k]
                   + f_1 * sld0_113[k]
                   - f_2 * sld1_113[k]
                   + f_3 * pc_z[k] * slf_189[k];

        t_285[k] = pb_y[k] * skg0_210[k]
                   - f_6 * pc_y[k] * skg1_210[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pc_y, pc_z, skg0_213, skf_130, \
                         skf_140, skf_141, skf_142, skg1_213, slf_190, \
                         slf_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_7 * skf_140[k]
                   + f_3 * pc_y[k] * slf_190[k];

        t_287[k] = f_13 * skf_130[k]
                   + f_3 * pc_z[k] * slf_190[k];

        t_288[k] = pb_y[k] * skg0_213[k]
                   + f_8 * skf_141[k]
                   - f_6 * pc_y[k] * skg1_213[k];

        t_289[k] = f_7 * skf_142[k]
                   + f_3 * pc_y[k] * slf_192[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_y, pc_x, pc_y, skg0_215, skf_196, \
                         skf_197, skf_198, skg1_215, slf_196, slf_197, \
                         slf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_y[k] * skg0_215[k]
                   - f_6 * pc_y[k] * skg1_215[k];

        t_291[k] = f_12 * skf_196[k]
                   + f_3 * pc_x[k] * slf_196[k];

        t_292[k] = f_12 * skf_197[k]
                   + f_3 * pc_x[k] * slf_197[k];

        t_293[k] = f_12 * skf_198[k]
                   + f_3 * pc_x[k] * slf_198[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pc_x, pc_y, pc_z, skf_136, skf_146, skf_199, \
                         sld0_117, sld1_117, slf_196, slf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_12 * skf_199[k]
                   + f_3 * pc_x[k] * slf_199[k];

        t_295[k] = f_7 * skf_146[k]
                   + f_1 * sld0_117[k]
                   - f_2 * sld1_117[k]
                   + f_3 * pc_y[k] * slf_196[k];

        t_296[k] = f_13 * skf_136[k]
                   + f_3 * pc_z[k] * slf_196[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_y, pc_y, skg0_224, skf_148, skf_149, \
                         skg1_224, sld0_119, sld1_119, slf_198, \
                         slf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_7 * skf_148[k]
                   + f_4 * sld0_119[k]
                   - f_5 * sld1_119[k]
                   + f_3 * pc_y[k] * slf_198[k];

        t_298[k] = f_7 * skf_149[k]
                   + f_3 * pc_y[k] * slf_199[k];

        t_299[k] = pb_y[k] * skg0_224[k]
                   - f_6 * pc_y[k] * skg1_224[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pc_x, pc_y, pc_z, skf_140, skf_200, \
                         skf_203, sld0_120, sld0_123, sld1_120, sld1_123, slf_200, \
                         slf_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_12 * skf_200[k]
                   + f_1 * sld0_120[k]
                   - f_2 * sld1_120[k]
                   + f_3 * pc_x[k] * slf_200[k];

        t_301[k] = f_3 * pc_y[k] * slf_200[k];

        t_302[k] = f_11 * skf_140[k]
                   + f_3 * pc_z[k] * slf_200[k];

        t_303[k] = f_12 * skf_203[k]
                   + f_4 * sld0_123[k]
                   - f_5 * sld1_123[k]
                   + f_3 * pc_x[k] * slf_203[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_x, pc_y, skf_205, skf_206, skf_207, \
                         sld0_125, sld1_125, slf_202, slf_205, slf_206, \
                         slf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_3 * pc_y[k] * slf_202[k];

        t_305[k] = f_12 * skf_205[k]
                   + f_4 * sld0_125[k]
                   - f_5 * sld1_125[k]
                   + f_3 * pc_x[k] * slf_205[k];

        t_306[k] = f_12 * skf_206[k]
                   + f_3 * pc_x[k] * slf_206[k];

        t_307[k] = f_12 * skf_207[k]
                   + f_3 * pc_x[k] * slf_207[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pc_x, pc_y, pc_z, skf_146, skf_208, \
                         skf_209, sld0_123, sld1_123, slf_206, slf_208, \
                         slf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_12 * skf_208[k]
                   + f_3 * pc_x[k] * slf_208[k];

        t_309[k] = f_12 * skf_209[k]
                   + f_3 * pc_x[k] * slf_209[k];

        t_310[k] = f_1 * sld0_123[k]
                   - f_2 * sld1_123[k]
                   + f_3 * pc_y[k] * slf_206[k];

        t_311[k] = f_11 * skf_146[k]
                   + f_3 * pc_z[k] * slf_206[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pc_x, pc_y, pc_z, skf_149, skf_210, \
                         sld0_125, sld0_126, sld1_125, sld1_126, slf_208, slf_209, \
                         slf_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_4 * sld0_125[k]
                   - f_5 * sld1_125[k]
                   + f_3 * pc_y[k] * slf_208[k];

        t_313[k] = f_3 * pc_y[k] * slf_209[k];

        t_314[k] = f_11 * skf_149[k]
                   + f_1 * sld0_125[k]
                   - f_2 * sld1_125[k]
                   + f_3 * pc_z[k] * slf_209[k];

        t_315[k] = f_8 * skf_210[k]
                   + f_1 * sld0_126[k]
                   - f_2 * sld1_126[k]
                   + f_3 * pc_x[k] * slf_210[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pc_x, pc_y, pc_z, skf_150, skf_152, \
                         skf_213, sld0_129, sld1_129, slf_210, slf_212, \
                         slf_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_10 * skf_150[k]
                   + f_3 * pc_y[k] * slf_210[k];

        t_317[k] = f_3 * pc_z[k] * slf_210[k];

        t_318[k] = f_8 * skf_213[k]
                   + f_4 * sld0_129[k]
                   - f_5 * sld1_129[k]
                   + f_3 * pc_x[k] * slf_213[k];

        t_319[k] = f_10 * skf_152[k]
                   + f_3 * pc_y[k] * slf_212[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pc_x, skf_215, skf_216, skf_217, skf_218, \
                         sld0_131, sld1_131, slf_215, slf_216, slf_217, \
                         slf_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_8 * skf_215[k]
                   + f_4 * sld0_131[k]
                   - f_5 * sld1_131[k]
                   + f_3 * pc_x[k] * slf_215[k];

        t_321[k] = f_8 * skf_216[k]
                   + f_3 * pc_x[k] * slf_216[k];

        t_322[k] = f_8 * skf_217[k]
                   + f_3 * pc_x[k] * slf_217[k];

        t_323[k] = f_8 * skf_218[k]
                   + f_3 * pc_x[k] * slf_218[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, pc_x, pc_y, pc_z, skf_156, skf_219, sld0_129, \
                         sld1_129, slf_216, slf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_8 * skf_219[k]
                   + f_3 * pc_x[k] * slf_219[k];

        t_325[k] = f_10 * skf_156[k]
                   + f_1 * sld0_129[k]
                   - f_2 * sld1_129[k]
                   + f_3 * pc_y[k] * slf_216[k];

        t_326[k] = f_3 * pc_z[k] * slf_216[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pb_z, pc_y, pc_z, skg0_225, skf_158, \
                         skf_159, skg1_225, sld0_131, sld1_131, slf_218, \
                         slf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_10 * skf_158[k]
                   + f_4 * sld0_131[k]
                   - f_5 * sld1_131[k]
                   + f_3 * pc_y[k] * slf_218[k];

        t_328[k] = f_10 * skf_159[k]
                   + f_3 * pc_y[k] * slf_219[k];

        t_329[k] = f_1 * sld0_131[k]
                   - f_2 * sld1_131[k]
                   + f_3 * pc_z[k] * slf_219[k];

        t_330[k] = pb_z[k] * skg0_225[k]
                   - f_6 * pc_z[k] * skg1_225[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pb_z, pc_y, pc_z, skg0_228, skf_150, \
                         skf_160, skf_162, skg1_228, slf_220, slf_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_11 * skf_160[k]
                   + f_3 * pc_y[k] * slf_220[k];

        t_332[k] = f_7 * skf_150[k]
                   + f_3 * pc_z[k] * slf_220[k];

        t_333[k] = pb_z[k] * skg0_228[k]
                   - f_6 * pc_z[k] * skg1_228[k];

        t_334[k] = f_11 * skf_162[k]
                   + f_3 * pc_y[k] * slf_222[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pc_x, skf_225, skf_226, skf_227, skf_228, \
                         sld0_137, sld1_137, slf_225, slf_226, slf_227, \
                         slf_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_8 * skf_225[k]
                   + f_4 * sld0_137[k]
                   - f_5 * sld1_137[k]
                   + f_3 * pc_x[k] * slf_225[k];

        t_336[k] = f_8 * skf_226[k]
                   + f_3 * pc_x[k] * slf_226[k];

        t_337[k] = f_8 * skf_227[k]
                   + f_3 * pc_x[k] * slf_227[k];

        t_338[k] = f_8 * skf_228[k]
                   + f_3 * pc_x[k] * slf_228[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pb_z, pc_x, pc_z, skg0_235, skf_156, skf_229, \
                         skg1_235, slf_226, slf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_8 * skf_229[k]
                   + f_3 * pc_x[k] * slf_229[k];

        t_340[k] = pb_z[k] * skg0_235[k]
                   - f_6 * pc_z[k] * skg1_235[k];

        t_341[k] = f_7 * skf_156[k]
                   + f_3 * pc_z[k] * slf_226[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_y, pc_z, skf_159, skf_168, skf_169, sld0_137, \
                         sld1_137, slf_228, slf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_11 * skf_168[k]
                   + f_4 * sld0_137[k]
                   - f_5 * sld1_137[k]
                   + f_3 * pc_y[k] * slf_228[k];

        t_343[k] = f_11 * skf_169[k]
                   + f_3 * pc_y[k] * slf_229[k];

        t_344[k] = f_7 * skf_159[k]
                   + f_1 * sld0_137[k]
                   - f_2 * sld1_137[k]
                   + f_3 * pc_z[k] * slf_229[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_y, pc_z, skf_160, skf_170, skf_230, \
                         sld0_138, sld1_138, slf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_8 * skf_230[k]
                   + f_1 * sld0_138[k]
                   - f_2 * sld1_138[k]
                   + f_3 * pc_x[k] * slf_230[k];

        t_346[k] = f_13 * skf_170[k]
                   + f_3 * pc_y[k] * slf_230[k];

        t_347[k] = f_8 * skf_160[k]
                   + f_3 * pc_z[k] * slf_230[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, skf_172, skf_233, skf_235, sld0_141, \
                         sld0_143, sld1_141, sld1_143, slf_232, slf_233, \
                         slf_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_8 * skf_233[k]
                   + f_4 * sld0_141[k]
                   - f_5 * sld1_141[k]
                   + f_3 * pc_x[k] * slf_233[k];

        t_349[k] = f_13 * skf_172[k]
                   + f_3 * pc_y[k] * slf_232[k];

        t_350[k] = f_8 * skf_235[k]
                   + f_4 * sld0_143[k]
                   - f_5 * sld1_143[k]
                   + f_3 * pc_x[k] * slf_235[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, pc_x, skf_236, skf_237, skf_238, skf_239, \
                         slf_236, slf_237, slf_238, slf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_8 * skf_236[k]
                   + f_3 * pc_x[k] * slf_236[k];

        t_352[k] = f_8 * skf_237[k]
                   + f_3 * pc_x[k] * slf_237[k];

        t_353[k] = f_8 * skf_238[k]
                   + f_3 * pc_x[k] * slf_238[k];

        t_354[k] = f_8 * skf_239[k]
                   + f_3 * pc_x[k] * slf_239[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, pc_y, pc_z, skf_166, skf_176, skf_178, sld0_141, \
                         sld0_143, sld1_141, sld1_143, slf_236, \
                         slf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_13 * skf_176[k]
                   + f_1 * sld0_141[k]
                   - f_2 * sld1_141[k]
                   + f_3 * pc_y[k] * slf_236[k];

        t_356[k] = f_8 * skf_166[k]
                   + f_3 * pc_z[k] * slf_236[k];

        t_357[k] = f_13 * skf_178[k]
                   + f_4 * sld0_143[k]
                   - f_5 * sld1_143[k]
                   + f_3 * pc_y[k] * slf_238[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, pc_x, pc_y, pc_z, skf_169, skf_179, skf_240, \
                         sld0_143, sld0_144, sld1_143, sld1_144, slf_239, \
                         slf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_13 * skf_179[k]
                   + f_3 * pc_y[k] * slf_239[k];

        t_359[k] = f_8 * skf_169[k]
                   + f_1 * sld0_143[k]
                   - f_2 * sld1_143[k]
                   + f_3 * pc_z[k] * slf_239[k];

        t_360[k] = f_8 * skf_240[k]
                   + f_1 * sld0_144[k]
                   - f_2 * sld1_144[k]
                   + f_3 * pc_x[k] * slf_240[k];
    }
}

static auto
compute_prim_slg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skg0,
                                                          const size_t skf, const size_t skg1,
                                                          const size_t sld0, const size_t sld1,
                                                          const size_t slf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 3.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.0 / q;

    auto *t_361 = buffer.data(target + 361);
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);
    auto *t_378 = buffer.data(target + 378);
    auto *t_379 = buffer.data(target + 379);
    auto *t_380 = buffer.data(target + 380);
    auto *t_381 = buffer.data(target + 381);
    auto *t_382 = buffer.data(target + 382);
    auto *t_383 = buffer.data(target + 383);
    auto *t_384 = buffer.data(target + 384);
    auto *t_385 = buffer.data(target + 385);
    auto *t_386 = buffer.data(target + 386);
    auto *t_387 = buffer.data(target + 387);
    auto *t_388 = buffer.data(target + 388);
    auto *t_389 = buffer.data(target + 389);
    auto *t_390 = buffer.data(target + 390);
    auto *t_391 = buffer.data(target + 391);
    auto *t_392 = buffer.data(target + 392);
    auto *t_393 = buffer.data(target + 393);
    auto *t_394 = buffer.data(target + 394);
    auto *t_395 = buffer.data(target + 395);
    auto *t_396 = buffer.data(target + 396);
    auto *t_397 = buffer.data(target + 397);
    auto *t_398 = buffer.data(target + 398);
    auto *t_399 = buffer.data(target + 399);
    auto *t_400 = buffer.data(target + 400);
    auto *t_401 = buffer.data(target + 401);
    auto *t_402 = buffer.data(target + 402);
    auto *t_403 = buffer.data(target + 403);
    auto *t_404 = buffer.data(target + 404);
    auto *t_405 = buffer.data(target + 405);
    auto *t_406 = buffer.data(target + 406);
    auto *t_407 = buffer.data(target + 407);
    auto *t_408 = buffer.data(target + 408);
    auto *t_409 = buffer.data(target + 409);
    auto *t_410 = buffer.data(target + 410);
    auto *t_411 = buffer.data(target + 411);
    auto *t_412 = buffer.data(target + 412);
    auto *t_413 = buffer.data(target + 413);
    auto *t_414 = buffer.data(target + 414);
    auto *t_415 = buffer.data(target + 415);
    auto *t_416 = buffer.data(target + 416);
    auto *t_417 = buffer.data(target + 417);
    auto *t_418 = buffer.data(target + 418);
    auto *t_419 = buffer.data(target + 419);
    auto *t_420 = buffer.data(target + 420);
    auto *t_421 = buffer.data(target + 421);
    auto *t_422 = buffer.data(target + 422);
    auto *t_423 = buffer.data(target + 423);
    auto *t_424 = buffer.data(target + 424);
    auto *t_425 = buffer.data(target + 425);
    auto *t_426 = buffer.data(target + 426);
    auto *t_427 = buffer.data(target + 427);
    auto *t_428 = buffer.data(target + 428);
    auto *t_429 = buffer.data(target + 429);
    auto *t_430 = buffer.data(target + 430);
    auto *t_431 = buffer.data(target + 431);
    auto *t_432 = buffer.data(target + 432);
    auto *t_433 = buffer.data(target + 433);
    auto *t_434 = buffer.data(target + 434);
    auto *t_435 = buffer.data(target + 435);
    auto *t_436 = buffer.data(target + 436);
    auto *t_437 = buffer.data(target + 437);
    auto *t_438 = buffer.data(target + 438);
    auto *t_439 = buffer.data(target + 439);
    auto *t_440 = buffer.data(target + 440);
    auto *t_441 = buffer.data(target + 441);
    auto *t_442 = buffer.data(target + 442);
    auto *t_443 = buffer.data(target + 443);
    auto *t_444 = buffer.data(target + 444);
    auto *t_445 = buffer.data(target + 445);
    auto *t_446 = buffer.data(target + 446);
    auto *t_447 = buffer.data(target + 447);
    auto *t_448 = buffer.data(target + 448);
    auto *t_449 = buffer.data(target + 449);
    auto *t_450 = buffer.data(target + 450);
    auto *t_451 = buffer.data(target + 451);
    auto *t_452 = buffer.data(target + 452);
    auto *t_453 = buffer.data(target + 453);
    auto *t_454 = buffer.data(target + 454);
    auto *t_455 = buffer.data(target + 455);
    auto *t_456 = buffer.data(target + 456);
    auto *t_457 = buffer.data(target + 457);
    auto *t_458 = buffer.data(target + 458);
    auto *t_459 = buffer.data(target + 459);
    auto *t_460 = buffer.data(target + 460);
    auto *t_461 = buffer.data(target + 461);
    auto *t_462 = buffer.data(target + 462);
    auto *t_463 = buffer.data(target + 463);
    auto *t_464 = buffer.data(target + 464);
    auto *t_465 = buffer.data(target + 465);
    auto *t_466 = buffer.data(target + 466);
    auto *t_467 = buffer.data(target + 467);
    auto *t_468 = buffer.data(target + 468);
    auto *t_469 = buffer.data(target + 469);
    auto *t_470 = buffer.data(target + 470);
    auto *t_471 = buffer.data(target + 471);
    auto *t_472 = buffer.data(target + 472);
    auto *t_473 = buffer.data(target + 473);
    auto *t_474 = buffer.data(target + 474);
    auto *t_475 = buffer.data(target + 475);
    auto *t_476 = buffer.data(target + 476);
    auto *t_477 = buffer.data(target + 477);
    auto *t_478 = buffer.data(target + 478);
    auto *t_479 = buffer.data(target + 479);
    auto *t_480 = buffer.data(target + 480);
    auto *t_481 = buffer.data(target + 481);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skg0_300 = buffer.data(skg0 + 300);
    const auto *skg0_303 = buffer.data(skg0 + 303);
    const auto *skg0_305 = buffer.data(skg0 + 305);
    const auto *skg0_314 = buffer.data(skg0 + 314);
    const auto *skg0_315 = buffer.data(skg0 + 315);
    const auto *skg0_318 = buffer.data(skg0 + 318);
    const auto *skg0_420 = buffer.data(skg0 + 420);
    const auto *skg0_423 = buffer.data(skg0 + 423);
    const auto *skg0_425 = buffer.data(skg0 + 425);
    const auto *skg0_430 = buffer.data(skg0 + 430);
    const auto *skg0_432 = buffer.data(skg0 + 432);
    const auto *skg0_434 = buffer.data(skg0 + 434);
    const auto *skg0_440 = buffer.data(skg0 + 440);
    const auto *skg0_445 = buffer.data(skg0 + 445);
    const auto *skg0_447 = buffer.data(skg0 + 447);
    const auto *skg0_449 = buffer.data(skg0 + 449);
    const auto *skg0_450 = buffer.data(skg0 + 450);
    const auto *skg0_453 = buffer.data(skg0 + 453);
    const auto *skg0_455 = buffer.data(skg0 + 455);
    const auto *skg0_460 = buffer.data(skg0 + 460);
    const auto *skg0_462 = buffer.data(skg0 + 462);
    const auto *skg0_464 = buffer.data(skg0 + 464);
    const auto *skg0_465 = buffer.data(skg0 + 465);
    const auto *skg0_468 = buffer.data(skg0 + 468);
    const auto *skg0_470 = buffer.data(skg0 + 470);
    const auto *skg0_475 = buffer.data(skg0 + 475);
    const auto *skg0_477 = buffer.data(skg0 + 477);
    const auto *skg0_479 = buffer.data(skg0 + 479);
    const auto *skg0_480 = buffer.data(skg0 + 480);

    const auto *skf_170 = buffer.data(skf + 170);
    const auto *skf_176 = buffer.data(skf + 176);
    const auto *skf_179 = buffer.data(skf + 179);
    const auto *skf_180 = buffer.data(skf + 180);
    const auto *skf_182 = buffer.data(skf + 182);
    const auto *skf_186 = buffer.data(skf + 186);
    const auto *skf_188 = buffer.data(skf + 188);
    const auto *skf_189 = buffer.data(skf + 189);
    const auto *skf_190 = buffer.data(skf + 190);
    const auto *skf_192 = buffer.data(skf + 192);
    const auto *skf_196 = buffer.data(skf + 196);
    const auto *skf_198 = buffer.data(skf + 198);
    const auto *skf_199 = buffer.data(skf + 199);
    const auto *skf_200 = buffer.data(skf + 200);
    const auto *skf_201 = buffer.data(skf + 201);
    const auto *skf_202 = buffer.data(skf + 202);
    const auto *skf_206 = buffer.data(skf + 206);
    const auto *skf_208 = buffer.data(skf + 208);
    const auto *skf_209 = buffer.data(skf + 209);
    const auto *skf_210 = buffer.data(skf + 210);
    const auto *skf_212 = buffer.data(skf + 212);
    const auto *skf_216 = buffer.data(skf + 216);
    const auto *skf_219 = buffer.data(skf + 219);
    const auto *skf_220 = buffer.data(skf + 220);
    const auto *skf_222 = buffer.data(skf + 222);
    const auto *skf_226 = buffer.data(skf + 226);
    const auto *skf_229 = buffer.data(skf + 229);
    const auto *skf_230 = buffer.data(skf + 230);
    const auto *skf_232 = buffer.data(skf + 232);
    const auto *skf_236 = buffer.data(skf + 236);
    const auto *skf_239 = buffer.data(skf + 239);
    const auto *skf_240 = buffer.data(skf + 240);
    const auto *skf_242 = buffer.data(skf + 242);
    const auto *skf_243 = buffer.data(skf + 243);
    const auto *skf_245 = buffer.data(skf + 245);
    const auto *skf_246 = buffer.data(skf + 246);
    const auto *skf_247 = buffer.data(skf + 247);
    const auto *skf_248 = buffer.data(skf + 248);
    const auto *skf_249 = buffer.data(skf + 249);
    const auto *skf_250 = buffer.data(skf + 250);
    const auto *skf_253 = buffer.data(skf + 253);
    const auto *skf_255 = buffer.data(skf + 255);
    const auto *skf_256 = buffer.data(skf + 256);
    const auto *skf_257 = buffer.data(skf + 257);
    const auto *skf_258 = buffer.data(skf + 258);
    const auto *skf_259 = buffer.data(skf + 259);
    const auto *skf_266 = buffer.data(skf + 266);
    const auto *skf_267 = buffer.data(skf + 267);
    const auto *skf_268 = buffer.data(skf + 268);
    const auto *skf_269 = buffer.data(skf + 269);
    const auto *skf_270 = buffer.data(skf + 270);
    const auto *skf_273 = buffer.data(skf + 273);
    const auto *skf_275 = buffer.data(skf + 275);
    const auto *skf_276 = buffer.data(skf + 276);
    const auto *skf_277 = buffer.data(skf + 277);
    const auto *skf_278 = buffer.data(skf + 278);
    const auto *skf_279 = buffer.data(skf + 279);
    const auto *skf_280 = buffer.data(skf + 280);
    const auto *skf_283 = buffer.data(skf + 283);
    const auto *skf_285 = buffer.data(skf + 285);
    const auto *skf_286 = buffer.data(skf + 286);
    const auto *skf_287 = buffer.data(skf + 287);
    const auto *skf_288 = buffer.data(skf + 288);
    const auto *skf_289 = buffer.data(skf + 289);
    const auto *skf_295 = buffer.data(skf + 295);
    const auto *skf_296 = buffer.data(skf + 296);
    const auto *skf_297 = buffer.data(skf + 297);
    const auto *skf_298 = buffer.data(skf + 298);
    const auto *skf_299 = buffer.data(skf + 299);
    const auto *skf_300 = buffer.data(skf + 300);
    const auto *skf_303 = buffer.data(skf + 303);
    const auto *skf_305 = buffer.data(skf + 305);
    const auto *skf_306 = buffer.data(skf + 306);
    const auto *skf_307 = buffer.data(skf + 307);
    const auto *skf_308 = buffer.data(skf + 308);
    const auto *skf_309 = buffer.data(skf + 309);
    const auto *skf_310 = buffer.data(skf + 310);
    const auto *skf_313 = buffer.data(skf + 313);
    const auto *skf_315 = buffer.data(skf + 315);
    const auto *skf_316 = buffer.data(skf + 316);
    const auto *skf_317 = buffer.data(skf + 317);
    const auto *skf_318 = buffer.data(skf + 318);
    const auto *skf_319 = buffer.data(skf + 319);
    const auto *skf_320 = buffer.data(skf + 320);

    const auto *skg1_300 = buffer.data(skg1 + 300);
    const auto *skg1_303 = buffer.data(skg1 + 303);
    const auto *skg1_305 = buffer.data(skg1 + 305);
    const auto *skg1_314 = buffer.data(skg1 + 314);
    const auto *skg1_315 = buffer.data(skg1 + 315);
    const auto *skg1_318 = buffer.data(skg1 + 318);
    const auto *skg1_420 = buffer.data(skg1 + 420);
    const auto *skg1_423 = buffer.data(skg1 + 423);
    const auto *skg1_425 = buffer.data(skg1 + 425);
    const auto *skg1_430 = buffer.data(skg1 + 430);
    const auto *skg1_432 = buffer.data(skg1 + 432);
    const auto *skg1_434 = buffer.data(skg1 + 434);
    const auto *skg1_440 = buffer.data(skg1 + 440);
    const auto *skg1_445 = buffer.data(skg1 + 445);
    const auto *skg1_447 = buffer.data(skg1 + 447);
    const auto *skg1_449 = buffer.data(skg1 + 449);
    const auto *skg1_450 = buffer.data(skg1 + 450);
    const auto *skg1_453 = buffer.data(skg1 + 453);
    const auto *skg1_455 = buffer.data(skg1 + 455);
    const auto *skg1_460 = buffer.data(skg1 + 460);
    const auto *skg1_462 = buffer.data(skg1 + 462);
    const auto *skg1_464 = buffer.data(skg1 + 464);
    const auto *skg1_465 = buffer.data(skg1 + 465);
    const auto *skg1_468 = buffer.data(skg1 + 468);
    const auto *skg1_470 = buffer.data(skg1 + 470);
    const auto *skg1_475 = buffer.data(skg1 + 475);
    const auto *skg1_477 = buffer.data(skg1 + 477);
    const auto *skg1_479 = buffer.data(skg1 + 479);
    const auto *skg1_480 = buffer.data(skg1 + 480);

    const auto *sld0_147 = buffer.data(sld0 + 147);
    const auto *sld0_149 = buffer.data(sld0 + 149);
    const auto *sld0_150 = buffer.data(sld0 + 150);
    const auto *sld0_153 = buffer.data(sld0 + 153);
    const auto *sld0_155 = buffer.data(sld0 + 155);
    const auto *sld0_159 = buffer.data(sld0 + 159);
    const auto *sld0_161 = buffer.data(sld0 + 161);
    const auto *sld0_162 = buffer.data(sld0 + 162);
    const auto *sld0_165 = buffer.data(sld0 + 165);
    const auto *sld0_167 = buffer.data(sld0 + 167);

    const auto *sld1_147 = buffer.data(sld1 + 147);
    const auto *sld1_149 = buffer.data(sld1 + 149);
    const auto *sld1_150 = buffer.data(sld1 + 150);
    const auto *sld1_153 = buffer.data(sld1 + 153);
    const auto *sld1_155 = buffer.data(sld1 + 155);
    const auto *sld1_159 = buffer.data(sld1 + 159);
    const auto *sld1_161 = buffer.data(sld1 + 161);
    const auto *sld1_162 = buffer.data(sld1 + 162);
    const auto *sld1_165 = buffer.data(sld1 + 165);
    const auto *sld1_167 = buffer.data(sld1 + 167);

    const auto *slf_240 = buffer.data(slf + 240);
    const auto *slf_242 = buffer.data(slf + 242);
    const auto *slf_243 = buffer.data(slf + 243);
    const auto *slf_245 = buffer.data(slf + 245);
    const auto *slf_246 = buffer.data(slf + 246);
    const auto *slf_247 = buffer.data(slf + 247);
    const auto *slf_248 = buffer.data(slf + 248);
    const auto *slf_249 = buffer.data(slf + 249);
    const auto *slf_250 = buffer.data(slf + 250);
    const auto *slf_252 = buffer.data(slf + 252);
    const auto *slf_253 = buffer.data(slf + 253);
    const auto *slf_255 = buffer.data(slf + 255);
    const auto *slf_256 = buffer.data(slf + 256);
    const auto *slf_257 = buffer.data(slf + 257);
    const auto *slf_258 = buffer.data(slf + 258);
    const auto *slf_259 = buffer.data(slf + 259);
    const auto *slf_260 = buffer.data(slf + 260);
    const auto *slf_262 = buffer.data(slf + 262);
    const auto *slf_266 = buffer.data(slf + 266);
    const auto *slf_267 = buffer.data(slf + 267);
    const auto *slf_268 = buffer.data(slf + 268);
    const auto *slf_269 = buffer.data(slf + 269);
    const auto *slf_270 = buffer.data(slf + 270);
    const auto *slf_272 = buffer.data(slf + 272);
    const auto *slf_273 = buffer.data(slf + 273);
    const auto *slf_275 = buffer.data(slf + 275);
    const auto *slf_276 = buffer.data(slf + 276);
    const auto *slf_277 = buffer.data(slf + 277);
    const auto *slf_278 = buffer.data(slf + 278);
    const auto *slf_279 = buffer.data(slf + 279);
    const auto *slf_280 = buffer.data(slf + 280);
    const auto *slf_282 = buffer.data(slf + 282);
    const auto *slf_286 = buffer.data(slf + 286);
    const auto *slf_287 = buffer.data(slf + 287);
    const auto *slf_288 = buffer.data(slf + 288);
    const auto *slf_289 = buffer.data(slf + 289);
    const auto *slf_290 = buffer.data(slf + 290);
    const auto *slf_292 = buffer.data(slf + 292);
    const auto *slf_296 = buffer.data(slf + 296);
    const auto *slf_297 = buffer.data(slf + 297);
    const auto *slf_298 = buffer.data(slf + 298);
    const auto *slf_299 = buffer.data(slf + 299);
    const auto *slf_300 = buffer.data(slf + 300);
    const auto *slf_302 = buffer.data(slf + 302);
    const auto *slf_306 = buffer.data(slf + 306);
    const auto *slf_307 = buffer.data(slf + 307);
    const auto *slf_308 = buffer.data(slf + 308);
    const auto *slf_309 = buffer.data(slf + 309);
    const auto *slf_310 = buffer.data(slf + 310);
    const auto *slf_312 = buffer.data(slf + 312);
    const auto *slf_316 = buffer.data(slf + 316);
    const auto *slf_317 = buffer.data(slf + 317);
    const auto *slf_318 = buffer.data(slf + 318);
    const auto *slf_319 = buffer.data(slf + 319);
    const auto *slf_320 = buffer.data(slf + 320);

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pc_x, pc_y, pc_z, skf_170, skf_180, \
                         skf_182, skf_243, sld0_147, sld1_147, slf_240, slf_242, \
                         slf_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_12 * skf_180[k]
                   + f_3 * pc_y[k] * slf_240[k];

        t_362[k] = f_12 * skf_170[k]
                   + f_3 * pc_z[k] * slf_240[k];

        t_363[k] = f_8 * skf_243[k]
                   + f_4 * sld0_147[k]
                   - f_5 * sld1_147[k]
                   + f_3 * pc_x[k] * slf_243[k];

        t_364[k] = f_12 * skf_182[k]
                   + f_3 * pc_y[k] * slf_242[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, skf_245, skf_246, skf_247, skf_248, \
                         sld0_149, sld1_149, slf_245, slf_246, slf_247, \
                         slf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_8 * skf_245[k]
                   + f_4 * sld0_149[k]
                   - f_5 * sld1_149[k]
                   + f_3 * pc_x[k] * slf_245[k];

        t_366[k] = f_8 * skf_246[k]
                   + f_3 * pc_x[k] * slf_246[k];

        t_367[k] = f_8 * skf_247[k]
                   + f_3 * pc_x[k] * slf_247[k];

        t_368[k] = f_8 * skf_248[k]
                   + f_3 * pc_x[k] * slf_248[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pc_x, pc_y, pc_z, skf_176, skf_186, skf_249, \
                         sld0_147, sld1_147, slf_246, slf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_8 * skf_249[k]
                   + f_3 * pc_x[k] * slf_249[k];

        t_370[k] = f_12 * skf_186[k]
                   + f_1 * sld0_147[k]
                   - f_2 * sld1_147[k]
                   + f_3 * pc_y[k] * slf_246[k];

        t_371[k] = f_12 * skf_176[k]
                   + f_3 * pc_z[k] * slf_246[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, skf_179, skf_188, skf_189, sld0_149, \
                         sld1_149, slf_248, slf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_12 * skf_188[k]
                   + f_4 * sld0_149[k]
                   - f_5 * sld1_149[k]
                   + f_3 * pc_y[k] * slf_248[k];

        t_373[k] = f_12 * skf_189[k]
                   + f_3 * pc_y[k] * slf_249[k];

        t_374[k] = f_12 * skf_179[k]
                   + f_1 * sld0_149[k]
                   - f_2 * sld1_149[k]
                   + f_3 * pc_z[k] * slf_249[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, pc_y, pc_z, skf_180, skf_190, skf_250, \
                         sld0_150, sld1_150, slf_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_8 * skf_250[k]
                   + f_1 * sld0_150[k]
                   - f_2 * sld1_150[k]
                   + f_3 * pc_x[k] * slf_250[k];

        t_376[k] = f_8 * skf_190[k]
                   + f_3 * pc_y[k] * slf_250[k];

        t_377[k] = f_13 * skf_180[k]
                   + f_3 * pc_z[k] * slf_250[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, pc_x, pc_y, skf_192, skf_253, skf_255, sld0_153, \
                         sld0_155, sld1_153, sld1_155, slf_252, slf_253, \
                         slf_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_8 * skf_253[k]
                   + f_4 * sld0_153[k]
                   - f_5 * sld1_153[k]
                   + f_3 * pc_x[k] * slf_253[k];

        t_379[k] = f_8 * skf_192[k]
                   + f_3 * pc_y[k] * slf_252[k];

        t_380[k] = f_8 * skf_255[k]
                   + f_4 * sld0_155[k]
                   - f_5 * sld1_155[k]
                   + f_3 * pc_x[k] * slf_255[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, skf_256, skf_257, skf_258, skf_259, \
                         slf_256, slf_257, slf_258, slf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_8 * skf_256[k]
                   + f_3 * pc_x[k] * slf_256[k];

        t_382[k] = f_8 * skf_257[k]
                   + f_3 * pc_x[k] * slf_257[k];

        t_383[k] = f_8 * skf_258[k]
                   + f_3 * pc_x[k] * slf_258[k];

        t_384[k] = f_8 * skf_259[k]
                   + f_3 * pc_x[k] * slf_259[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, skf_186, skf_196, skf_198, sld0_153, \
                         sld0_155, sld1_153, sld1_155, slf_256, \
                         slf_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_8 * skf_196[k]
                   + f_1 * sld0_153[k]
                   - f_2 * sld1_153[k]
                   + f_3 * pc_y[k] * slf_256[k];

        t_386[k] = f_13 * skf_186[k]
                   + f_3 * pc_z[k] * slf_256[k];

        t_387[k] = f_8 * skf_198[k]
                   + f_4 * sld0_155[k]
                   - f_5 * sld1_155[k]
                   + f_3 * pc_y[k] * slf_258[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pb_y, pc_y, pc_z, skg0_300, skf_189, \
                         skf_199, skf_200, skg1_300, sld0_155, sld1_155, slf_259, \
                         slf_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_8 * skf_199[k]
                   + f_3 * pc_y[k] * slf_259[k];

        t_389[k] = f_13 * skf_189[k]
                   + f_1 * sld0_155[k]
                   - f_2 * sld1_155[k]
                   + f_3 * pc_z[k] * slf_259[k];

        t_390[k] = pb_y[k] * skg0_300[k]
                   - f_6 * pc_y[k] * skg1_300[k];

        t_391[k] = f_7 * skf_200[k]
                   + f_3 * pc_y[k] * slf_260[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, pb_y, pc_y, pc_z, skg0_303, skg0_305, \
                         skf_190, skf_201, skf_202, skg1_303, skg1_305, slf_260, \
                         slf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_11 * skf_190[k]
                   + f_3 * pc_z[k] * slf_260[k];

        t_393[k] = pb_y[k] * skg0_303[k]
                   + f_8 * skf_201[k]
                   - f_6 * pc_y[k] * skg1_303[k];

        t_394[k] = f_7 * skf_202[k]
                   + f_3 * pc_y[k] * slf_262[k];

        t_395[k] = pb_y[k] * skg0_305[k]
                   - f_6 * pc_y[k] * skg1_305[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, pc_x, skf_266, skf_267, skf_268, skf_269, \
                         slf_266, slf_267, slf_268, slf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_8 * skf_266[k]
                   + f_3 * pc_x[k] * slf_266[k];

        t_397[k] = f_8 * skf_267[k]
                   + f_3 * pc_x[k] * slf_267[k];

        t_398[k] = f_8 * skf_268[k]
                   + f_3 * pc_x[k] * slf_268[k];

        t_399[k] = f_8 * skf_269[k]
                   + f_3 * pc_x[k] * slf_269[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pc_y, pc_z, skf_196, skf_206, skf_208, sld0_159, \
                         sld0_161, sld1_159, sld1_161, slf_266, \
                         slf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_7 * skf_206[k]
                   + f_1 * sld0_159[k]
                   - f_2 * sld1_159[k]
                   + f_3 * pc_y[k] * slf_266[k];

        t_401[k] = f_11 * skf_196[k]
                   + f_3 * pc_z[k] * slf_266[k];

        t_402[k] = f_7 * skf_208[k]
                   + f_4 * sld0_161[k]
                   - f_5 * sld1_161[k]
                   + f_3 * pc_y[k] * slf_268[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, pb_y, pc_x, pc_y, skg0_314, skf_209, \
                         skf_270, skg1_314, sld0_162, sld1_162, slf_269, \
                         slf_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_7 * skf_209[k]
                   + f_3 * pc_y[k] * slf_269[k];

        t_404[k] = pb_y[k] * skg0_314[k]
                   - f_6 * pc_y[k] * skg1_314[k];

        t_405[k] = f_8 * skf_270[k]
                   + f_1 * sld0_162[k]
                   - f_2 * sld1_162[k]
                   + f_3 * pc_x[k] * slf_270[k];

        t_406[k] = f_3 * pc_y[k] * slf_270[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pc_x, pc_y, pc_z, skf_200, skf_273, sld0_165, \
                         sld1_165, slf_270, slf_272, slf_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_10 * skf_200[k]
                   + f_3 * pc_z[k] * slf_270[k];

        t_408[k] = f_8 * skf_273[k]
                   + f_4 * sld0_165[k]
                   - f_5 * sld1_165[k]
                   + f_3 * pc_x[k] * slf_273[k];

        t_409[k] = f_3 * pc_y[k] * slf_272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, skf_275, skf_276, skf_277, skf_278, \
                         sld0_167, sld1_167, slf_275, slf_276, slf_277, \
                         slf_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_8 * skf_275[k]
                   + f_4 * sld0_167[k]
                   - f_5 * sld1_167[k]
                   + f_3 * pc_x[k] * slf_275[k];

        t_411[k] = f_8 * skf_276[k]
                   + f_3 * pc_x[k] * slf_276[k];

        t_412[k] = f_8 * skf_277[k]
                   + f_3 * pc_x[k] * slf_277[k];

        t_413[k] = f_8 * skf_278[k]
                   + f_3 * pc_x[k] * slf_278[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pc_x, pc_y, pc_z, skf_206, skf_279, \
                         sld0_165, sld0_167, sld1_165, sld1_167, slf_276, slf_278, \
                         slf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_8 * skf_279[k]
                   + f_3 * pc_x[k] * slf_279[k];

        t_415[k] = f_1 * sld0_165[k]
                   - f_2 * sld1_165[k]
                   + f_3 * pc_y[k] * slf_276[k];

        t_416[k] = f_10 * skf_206[k]
                   + f_3 * pc_z[k] * slf_276[k];

        t_417[k] = f_4 * sld0_167[k]
                   - f_5 * sld1_167[k]
                   + f_3 * pc_y[k] * slf_278[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, pb_x, pc_x, pc_y, pc_z, skg0_420, skf_209, \
                         skf_280, skg1_420, sld0_167, sld1_167, \
                         slf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_3 * pc_y[k] * slf_279[k];

        t_419[k] = f_10 * skf_209[k]
                   + f_1 * sld0_167[k]
                   - f_2 * sld1_167[k]
                   + f_3 * pc_z[k] * slf_279[k];

        t_420[k] = pb_x[k] * skg0_420[k]
                   + f_13 * skf_280[k]
                   - f_6 * pc_x[k] * skg1_420[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pb_x, pc_x, pc_y, pc_z, skg0_423, \
                         skf_210, skf_212, skf_283, skg1_423, slf_280, \
                         slf_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_9 * skf_210[k]
                   + f_3 * pc_y[k] * slf_280[k];

        t_422[k] = f_3 * pc_z[k] * slf_280[k];

        t_423[k] = pb_x[k] * skg0_423[k]
                   + f_8 * skf_283[k]
                   - f_6 * pc_x[k] * skg1_423[k];

        t_424[k] = f_9 * skf_212[k]
                   + f_3 * pc_y[k] * slf_282[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pb_x, pc_x, skg0_425, skf_285, skf_286, \
                         skf_287, skf_288, skg1_425, slf_286, slf_287, \
                         slf_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = pb_x[k] * skg0_425[k]
                   + f_8 * skf_285[k]
                   - f_6 * pc_x[k] * skg1_425[k];

        t_426[k] = f_7 * skf_286[k]
                   + f_3 * pc_x[k] * slf_286[k];

        t_427[k] = f_7 * skf_287[k]
                   + f_3 * pc_x[k] * slf_287[k];

        t_428[k] = f_7 * skf_288[k]
                   + f_3 * pc_x[k] * slf_288[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pb_x, pc_x, pc_z, skg0_430, skg0_432, \
                         skf_289, skg1_430, skg1_432, slf_286, \
                         slf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_7 * skf_289[k]
                   + f_3 * pc_x[k] * slf_289[k];

        t_430[k] = pb_x[k] * skg0_430[k]
                   - f_6 * pc_x[k] * skg1_430[k];

        t_431[k] = f_3 * pc_z[k] * slf_286[k];

        t_432[k] = pb_x[k] * skg0_432[k]
                   - f_6 * pc_x[k] * skg1_432[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, pb_x, pb_z, pc_x, pc_y, pc_z, skg0_315, \
                         skg0_434, skf_219, skg1_315, skg1_434, \
                         slf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_9 * skf_219[k]
                   + f_3 * pc_y[k] * slf_289[k];

        t_434[k] = pb_x[k] * skg0_434[k]
                   - f_6 * pc_x[k] * skg1_434[k];

        t_435[k] = pb_z[k] * skg0_315[k]
                   - f_6 * pc_z[k] * skg1_315[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, pb_z, pc_y, pc_z, skg0_318, skf_210, \
                         skf_220, skf_222, skg1_318, slf_290, slf_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_10 * skf_220[k]
                   + f_3 * pc_y[k] * slf_290[k];

        t_437[k] = f_7 * skf_210[k]
                   + f_3 * pc_z[k] * slf_290[k];

        t_438[k] = pb_z[k] * skg0_318[k]
                   - f_6 * pc_z[k] * skg1_318[k];

        t_439[k] = f_10 * skf_222[k]
                   + f_3 * pc_y[k] * slf_292[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, pb_x, pc_x, skg0_440, skf_295, skf_296, \
                         skf_297, skf_298, skg1_440, slf_296, slf_297, \
                         slf_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = pb_x[k] * skg0_440[k]
                   + f_8 * skf_295[k]
                   - f_6 * pc_x[k] * skg1_440[k];

        t_441[k] = f_7 * skf_296[k]
                   + f_3 * pc_x[k] * slf_296[k];

        t_442[k] = f_7 * skf_297[k]
                   + f_3 * pc_x[k] * slf_297[k];

        t_443[k] = f_7 * skf_298[k]
                   + f_3 * pc_x[k] * slf_298[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, pb_x, pc_x, pc_z, skg0_445, skg0_447, \
                         skf_216, skf_299, skg1_445, skg1_447, slf_296, \
                         slf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_7 * skf_299[k]
                   + f_3 * pc_x[k] * slf_299[k];

        t_445[k] = pb_x[k] * skg0_445[k]
                   - f_6 * pc_x[k] * skg1_445[k];

        t_446[k] = f_7 * skf_216[k]
                   + f_3 * pc_z[k] * slf_296[k];

        t_447[k] = pb_x[k] * skg0_447[k]
                   - f_6 * pc_x[k] * skg1_447[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, pb_x, pc_x, pc_y, skg0_449, skg0_450, \
                         skf_229, skf_230, skf_300, skg1_449, skg1_450, slf_299, \
                         slf_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = f_10 * skf_229[k]
                   + f_3 * pc_y[k] * slf_299[k];

        t_449[k] = pb_x[k] * skg0_449[k]
                   - f_6 * pc_x[k] * skg1_449[k];

        t_450[k] = pb_x[k] * skg0_450[k]
                   + f_13 * skf_300[k]
                   - f_6 * pc_x[k] * skg1_450[k];

        t_451[k] = f_11 * skf_230[k]
                   + f_3 * pc_y[k] * slf_300[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, pb_x, pc_x, pc_y, pc_z, skg0_453, skf_220, \
                         skf_232, skf_303, skg1_453, slf_300, slf_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_8 * skf_220[k]
                   + f_3 * pc_z[k] * slf_300[k];

        t_453[k] = pb_x[k] * skg0_453[k]
                   + f_8 * skf_303[k]
                   - f_6 * pc_x[k] * skg1_453[k];

        t_454[k] = f_11 * skf_232[k]
                   + f_3 * pc_y[k] * slf_302[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, pb_x, pc_x, skg0_455, skf_305, skf_306, \
                         skf_307, skf_308, skg1_455, slf_306, slf_307, \
                         slf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = pb_x[k] * skg0_455[k]
                   + f_8 * skf_305[k]
                   - f_6 * pc_x[k] * skg1_455[k];

        t_456[k] = f_7 * skf_306[k]
                   + f_3 * pc_x[k] * slf_306[k];

        t_457[k] = f_7 * skf_307[k]
                   + f_3 * pc_x[k] * slf_307[k];

        t_458[k] = f_7 * skf_308[k]
                   + f_3 * pc_x[k] * slf_308[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, t_462, pb_x, pc_x, pc_z, skg0_460, skg0_462, \
                         skf_226, skf_309, skg1_460, skg1_462, slf_306, \
                         slf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_7 * skf_309[k]
                   + f_3 * pc_x[k] * slf_309[k];

        t_460[k] = pb_x[k] * skg0_460[k]
                   - f_6 * pc_x[k] * skg1_460[k];

        t_461[k] = f_8 * skf_226[k]
                   + f_3 * pc_z[k] * slf_306[k];

        t_462[k] = pb_x[k] * skg0_462[k]
                   - f_6 * pc_x[k] * skg1_462[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, pb_x, pc_x, pc_y, skg0_464, skg0_465, \
                         skf_239, skf_240, skf_310, skg1_464, skg1_465, slf_309, \
                         slf_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_11 * skf_239[k]
                   + f_3 * pc_y[k] * slf_309[k];

        t_464[k] = pb_x[k] * skg0_464[k]
                   - f_6 * pc_x[k] * skg1_464[k];

        t_465[k] = pb_x[k] * skg0_465[k]
                   + f_13 * skf_310[k]
                   - f_6 * pc_x[k] * skg1_465[k];

        t_466[k] = f_13 * skf_240[k]
                   + f_3 * pc_y[k] * slf_310[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, pb_x, pc_x, pc_y, pc_z, skg0_468, skf_230, \
                         skf_242, skf_313, skg1_468, slf_310, slf_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_12 * skf_230[k]
                   + f_3 * pc_z[k] * slf_310[k];

        t_468[k] = pb_x[k] * skg0_468[k]
                   + f_8 * skf_313[k]
                   - f_6 * pc_x[k] * skg1_468[k];

        t_469[k] = f_13 * skf_242[k]
                   + f_3 * pc_y[k] * slf_312[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pb_x, pc_x, skg0_470, skf_315, skf_316, \
                         skf_317, skf_318, skg1_470, slf_316, slf_317, \
                         slf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = pb_x[k] * skg0_470[k]
                   + f_8 * skf_315[k]
                   - f_6 * pc_x[k] * skg1_470[k];

        t_471[k] = f_7 * skf_316[k]
                   + f_3 * pc_x[k] * slf_316[k];

        t_472[k] = f_7 * skf_317[k]
                   + f_3 * pc_x[k] * slf_317[k];

        t_473[k] = f_7 * skf_318[k]
                   + f_3 * pc_x[k] * slf_318[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_x, pc_x, pc_z, skg0_475, skg0_477, \
                         skf_236, skf_319, skg1_475, skg1_477, slf_316, \
                         slf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_7 * skf_319[k]
                   + f_3 * pc_x[k] * slf_319[k];

        t_475[k] = pb_x[k] * skg0_475[k]
                   - f_6 * pc_x[k] * skg1_475[k];

        t_476[k] = f_12 * skf_236[k]
                   + f_3 * pc_z[k] * slf_316[k];

        t_477[k] = pb_x[k] * skg0_477[k]
                   - f_6 * pc_x[k] * skg1_477[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, pb_x, pc_x, pc_y, skg0_479, skg0_480, \
                         skf_249, skf_250, skf_320, skg1_479, skg1_480, slf_319, \
                         slf_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_13 * skf_249[k]
                   + f_3 * pc_y[k] * slf_319[k];

        t_479[k] = pb_x[k] * skg0_479[k]
                   - f_6 * pc_x[k] * skg1_479[k];

        t_480[k] = pb_x[k] * skg0_480[k]
                   + f_13 * skf_320[k]
                   - f_6 * pc_x[k] * skg1_480[k];

        t_481[k] = f_12 * skf_250[k]
                   + f_3 * pc_y[k] * slf_320[k];
    }
}

static auto
compute_prim_slg_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skg0,
                                                          const size_t skf, const size_t skg1,
                                                          const size_t sld0, const size_t sld1,
                                                          const size_t slf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 3.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.0 / q;

    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);
    auto *t_495 = buffer.data(target + 495);
    auto *t_496 = buffer.data(target + 496);
    auto *t_497 = buffer.data(target + 497);
    auto *t_498 = buffer.data(target + 498);
    auto *t_499 = buffer.data(target + 499);
    auto *t_500 = buffer.data(target + 500);
    auto *t_501 = buffer.data(target + 501);
    auto *t_502 = buffer.data(target + 502);
    auto *t_503 = buffer.data(target + 503);
    auto *t_504 = buffer.data(target + 504);
    auto *t_505 = buffer.data(target + 505);
    auto *t_506 = buffer.data(target + 506);
    auto *t_507 = buffer.data(target + 507);
    auto *t_508 = buffer.data(target + 508);
    auto *t_509 = buffer.data(target + 509);
    auto *t_510 = buffer.data(target + 510);
    auto *t_511 = buffer.data(target + 511);
    auto *t_512 = buffer.data(target + 512);
    auto *t_513 = buffer.data(target + 513);
    auto *t_514 = buffer.data(target + 514);
    auto *t_515 = buffer.data(target + 515);
    auto *t_516 = buffer.data(target + 516);
    auto *t_517 = buffer.data(target + 517);
    auto *t_518 = buffer.data(target + 518);
    auto *t_519 = buffer.data(target + 519);
    auto *t_520 = buffer.data(target + 520);
    auto *t_521 = buffer.data(target + 521);
    auto *t_522 = buffer.data(target + 522);
    auto *t_523 = buffer.data(target + 523);
    auto *t_524 = buffer.data(target + 524);
    auto *t_525 = buffer.data(target + 525);
    auto *t_526 = buffer.data(target + 526);
    auto *t_527 = buffer.data(target + 527);
    auto *t_528 = buffer.data(target + 528);
    auto *t_529 = buffer.data(target + 529);
    auto *t_530 = buffer.data(target + 530);
    auto *t_531 = buffer.data(target + 531);
    auto *t_532 = buffer.data(target + 532);
    auto *t_533 = buffer.data(target + 533);
    auto *t_534 = buffer.data(target + 534);
    auto *t_535 = buffer.data(target + 535);
    auto *t_536 = buffer.data(target + 536);
    auto *t_537 = buffer.data(target + 537);
    auto *t_538 = buffer.data(target + 538);
    auto *t_539 = buffer.data(target + 539);
    auto *t_540 = buffer.data(target + 540);
    auto *t_541 = buffer.data(target + 541);
    auto *t_542 = buffer.data(target + 542);
    auto *t_543 = buffer.data(target + 543);
    auto *t_544 = buffer.data(target + 544);
    auto *t_545 = buffer.data(target + 545);
    auto *t_546 = buffer.data(target + 546);
    auto *t_547 = buffer.data(target + 547);
    auto *t_548 = buffer.data(target + 548);
    auto *t_549 = buffer.data(target + 549);
    auto *t_550 = buffer.data(target + 550);
    auto *t_551 = buffer.data(target + 551);
    auto *t_552 = buffer.data(target + 552);
    auto *t_553 = buffer.data(target + 553);
    auto *t_554 = buffer.data(target + 554);
    auto *t_555 = buffer.data(target + 555);
    auto *t_556 = buffer.data(target + 556);
    auto *t_557 = buffer.data(target + 557);
    auto *t_558 = buffer.data(target + 558);
    auto *t_559 = buffer.data(target + 559);
    auto *t_560 = buffer.data(target + 560);
    auto *t_561 = buffer.data(target + 561);
    auto *t_562 = buffer.data(target + 562);
    auto *t_563 = buffer.data(target + 563);
    auto *t_564 = buffer.data(target + 564);
    auto *t_565 = buffer.data(target + 565);
    auto *t_566 = buffer.data(target + 566);
    auto *t_567 = buffer.data(target + 567);
    auto *t_568 = buffer.data(target + 568);
    auto *t_569 = buffer.data(target + 569);
    auto *t_570 = buffer.data(target + 570);
    auto *t_571 = buffer.data(target + 571);
    auto *t_572 = buffer.data(target + 572);
    auto *t_573 = buffer.data(target + 573);
    auto *t_574 = buffer.data(target + 574);
    auto *t_575 = buffer.data(target + 575);
    auto *t_576 = buffer.data(target + 576);
    auto *t_577 = buffer.data(target + 577);
    auto *t_578 = buffer.data(target + 578);
    auto *t_579 = buffer.data(target + 579);
    auto *t_580 = buffer.data(target + 580);
    auto *t_581 = buffer.data(target + 581);
    auto *t_582 = buffer.data(target + 582);
    auto *t_583 = buffer.data(target + 583);
    auto *t_584 = buffer.data(target + 584);
    auto *t_585 = buffer.data(target + 585);
    auto *t_586 = buffer.data(target + 586);
    auto *t_587 = buffer.data(target + 587);
    auto *t_588 = buffer.data(target + 588);
    auto *t_589 = buffer.data(target + 589);
    auto *t_590 = buffer.data(target + 590);
    auto *t_591 = buffer.data(target + 591);
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);
    auto *t_604 = buffer.data(target + 604);
    auto *t_605 = buffer.data(target + 605);
    auto *t_606 = buffer.data(target + 606);
    auto *t_607 = buffer.data(target + 607);
    auto *t_608 = buffer.data(target + 608);
    auto *t_609 = buffer.data(target + 609);
    auto *t_610 = buffer.data(target + 610);
    auto *t_611 = buffer.data(target + 611);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skg0_405 = buffer.data(skg0 + 405);
    const auto *skg0_410 = buffer.data(skg0 + 410);
    const auto *skg0_420 = buffer.data(skg0 + 420);
    const auto *skg0_423 = buffer.data(skg0 + 423);
    const auto *skg0_430 = buffer.data(skg0 + 430);
    const auto *skg0_432 = buffer.data(skg0 + 432);
    const auto *skg0_483 = buffer.data(skg0 + 483);
    const auto *skg0_485 = buffer.data(skg0 + 485);
    const auto *skg0_490 = buffer.data(skg0 + 490);
    const auto *skg0_492 = buffer.data(skg0 + 492);
    const auto *skg0_494 = buffer.data(skg0 + 494);
    const auto *skg0_495 = buffer.data(skg0 + 495);
    const auto *skg0_498 = buffer.data(skg0 + 498);
    const auto *skg0_500 = buffer.data(skg0 + 500);
    const auto *skg0_505 = buffer.data(skg0 + 505);
    const auto *skg0_507 = buffer.data(skg0 + 507);
    const auto *skg0_509 = buffer.data(skg0 + 509);
    const auto *skg0_513 = buffer.data(skg0 + 513);
    const auto *skg0_520 = buffer.data(skg0 + 520);
    const auto *skg0_522 = buffer.data(skg0 + 522);
    const auto *skg0_524 = buffer.data(skg0 + 524);
    const auto *skg0_525 = buffer.data(skg0 + 525);
    const auto *skg0_528 = buffer.data(skg0 + 528);
    const auto *skg0_530 = buffer.data(skg0 + 530);
    const auto *skg0_535 = buffer.data(skg0 + 535);
    const auto *skg0_537 = buffer.data(skg0 + 537);
    const auto *skg0_539 = buffer.data(skg0 + 539);

    const auto *skf_240 = buffer.data(skf + 240);
    const auto *skf_246 = buffer.data(skf + 246);
    const auto *skf_250 = buffer.data(skf + 250);
    const auto *skf_252 = buffer.data(skf + 252);
    const auto *skf_256 = buffer.data(skf + 256);
    const auto *skf_259 = buffer.data(skf + 259);
    const auto *skf_260 = buffer.data(skf + 260);
    const auto *skf_262 = buffer.data(skf + 262);
    const auto *skf_266 = buffer.data(skf + 266);
    const auto *skf_269 = buffer.data(skf + 269);
    const auto *skf_270 = buffer.data(skf + 270);
    const auto *skf_272 = buffer.data(skf + 272);
    const auto *skf_276 = buffer.data(skf + 276);
    const auto *skf_279 = buffer.data(skf + 279);
    const auto *skf_280 = buffer.data(skf + 280);
    const auto *skf_282 = buffer.data(skf + 282);
    const auto *skf_286 = buffer.data(skf + 286);
    const auto *skf_287 = buffer.data(skf + 287);
    const auto *skf_288 = buffer.data(skf + 288);
    const auto *skf_289 = buffer.data(skf + 289);
    const auto *skf_290 = buffer.data(skf + 290);
    const auto *skf_292 = buffer.data(skf + 292);
    const auto *skf_296 = buffer.data(skf + 296);
    const auto *skf_299 = buffer.data(skf + 299);
    const auto *skf_300 = buffer.data(skf + 300);
    const auto *skf_302 = buffer.data(skf + 302);
    const auto *skf_306 = buffer.data(skf + 306);
    const auto *skf_308 = buffer.data(skf + 308);
    const auto *skf_309 = buffer.data(skf + 309);
    const auto *skf_310 = buffer.data(skf + 310);
    const auto *skf_312 = buffer.data(skf + 312);
    const auto *skf_316 = buffer.data(skf + 316);
    const auto *skf_318 = buffer.data(skf + 318);
    const auto *skf_319 = buffer.data(skf + 319);
    const auto *skf_320 = buffer.data(skf + 320);
    const auto *skf_322 = buffer.data(skf + 322);
    const auto *skf_323 = buffer.data(skf + 323);
    const auto *skf_325 = buffer.data(skf + 325);
    const auto *skf_326 = buffer.data(skf + 326);
    const auto *skf_327 = buffer.data(skf + 327);
    const auto *skf_328 = buffer.data(skf + 328);
    const auto *skf_329 = buffer.data(skf + 329);
    const auto *skf_330 = buffer.data(skf + 330);
    const auto *skf_333 = buffer.data(skf + 333);
    const auto *skf_335 = buffer.data(skf + 335);
    const auto *skf_336 = buffer.data(skf + 336);
    const auto *skf_337 = buffer.data(skf + 337);
    const auto *skf_338 = buffer.data(skf + 338);
    const auto *skf_339 = buffer.data(skf + 339);
    const auto *skf_343 = buffer.data(skf + 343);
    const auto *skf_346 = buffer.data(skf + 346);
    const auto *skf_347 = buffer.data(skf + 347);
    const auto *skf_348 = buffer.data(skf + 348);
    const auto *skf_349 = buffer.data(skf + 349);
    const auto *skf_350 = buffer.data(skf + 350);
    const auto *skf_353 = buffer.data(skf + 353);
    const auto *skf_355 = buffer.data(skf + 355);
    const auto *skf_356 = buffer.data(skf + 356);
    const auto *skf_357 = buffer.data(skf + 357);
    const auto *skf_358 = buffer.data(skf + 358);
    const auto *skf_359 = buffer.data(skf + 359);

    const auto *skg1_405 = buffer.data(skg1 + 405);
    const auto *skg1_410 = buffer.data(skg1 + 410);
    const auto *skg1_420 = buffer.data(skg1 + 420);
    const auto *skg1_423 = buffer.data(skg1 + 423);
    const auto *skg1_430 = buffer.data(skg1 + 430);
    const auto *skg1_432 = buffer.data(skg1 + 432);
    const auto *skg1_483 = buffer.data(skg1 + 483);
    const auto *skg1_485 = buffer.data(skg1 + 485);
    const auto *skg1_490 = buffer.data(skg1 + 490);
    const auto *skg1_492 = buffer.data(skg1 + 492);
    const auto *skg1_494 = buffer.data(skg1 + 494);
    const auto *skg1_495 = buffer.data(skg1 + 495);
    const auto *skg1_498 = buffer.data(skg1 + 498);
    const auto *skg1_500 = buffer.data(skg1 + 500);
    const auto *skg1_505 = buffer.data(skg1 + 505);
    const auto *skg1_507 = buffer.data(skg1 + 507);
    const auto *skg1_509 = buffer.data(skg1 + 509);
    const auto *skg1_513 = buffer.data(skg1 + 513);
    const auto *skg1_520 = buffer.data(skg1 + 520);
    const auto *skg1_522 = buffer.data(skg1 + 522);
    const auto *skg1_524 = buffer.data(skg1 + 524);
    const auto *skg1_525 = buffer.data(skg1 + 525);
    const auto *skg1_528 = buffer.data(skg1 + 528);
    const auto *skg1_530 = buffer.data(skg1 + 530);
    const auto *skg1_535 = buffer.data(skg1 + 535);
    const auto *skg1_537 = buffer.data(skg1 + 537);
    const auto *skg1_539 = buffer.data(skg1 + 539);

    const auto *sld0_216 = buffer.data(sld0 + 216);
    const auto *sld0_219 = buffer.data(sld0 + 219);
    const auto *sld0_221 = buffer.data(sld0 + 221);
    const auto *sld0_227 = buffer.data(sld0 + 227);
    const auto *sld0_228 = buffer.data(sld0 + 228);
    const auto *sld0_231 = buffer.data(sld0 + 231);
    const auto *sld0_233 = buffer.data(sld0 + 233);
    const auto *sld0_234 = buffer.data(sld0 + 234);
    const auto *sld0_237 = buffer.data(sld0 + 237);
    const auto *sld0_239 = buffer.data(sld0 + 239);
    const auto *sld0_240 = buffer.data(sld0 + 240);
    const auto *sld0_243 = buffer.data(sld0 + 243);
    const auto *sld0_245 = buffer.data(sld0 + 245);

    const auto *sld1_216 = buffer.data(sld1 + 216);
    const auto *sld1_219 = buffer.data(sld1 + 219);
    const auto *sld1_221 = buffer.data(sld1 + 221);
    const auto *sld1_227 = buffer.data(sld1 + 227);
    const auto *sld1_228 = buffer.data(sld1 + 228);
    const auto *sld1_231 = buffer.data(sld1 + 231);
    const auto *sld1_233 = buffer.data(sld1 + 233);
    const auto *sld1_234 = buffer.data(sld1 + 234);
    const auto *sld1_237 = buffer.data(sld1 + 237);
    const auto *sld1_239 = buffer.data(sld1 + 239);
    const auto *sld1_240 = buffer.data(sld1 + 240);
    const auto *sld1_243 = buffer.data(sld1 + 243);
    const auto *sld1_245 = buffer.data(sld1 + 245);

    const auto *slf_320 = buffer.data(slf + 320);
    const auto *slf_322 = buffer.data(slf + 322);
    const auto *slf_326 = buffer.data(slf + 326);
    const auto *slf_327 = buffer.data(slf + 327);
    const auto *slf_328 = buffer.data(slf + 328);
    const auto *slf_329 = buffer.data(slf + 329);
    const auto *slf_330 = buffer.data(slf + 330);
    const auto *slf_332 = buffer.data(slf + 332);
    const auto *slf_336 = buffer.data(slf + 336);
    const auto *slf_337 = buffer.data(slf + 337);
    const auto *slf_338 = buffer.data(slf + 338);
    const auto *slf_339 = buffer.data(slf + 339);
    const auto *slf_340 = buffer.data(slf + 340);
    const auto *slf_342 = buffer.data(slf + 342);
    const auto *slf_346 = buffer.data(slf + 346);
    const auto *slf_347 = buffer.data(slf + 347);
    const auto *slf_348 = buffer.data(slf + 348);
    const auto *slf_349 = buffer.data(slf + 349);
    const auto *slf_350 = buffer.data(slf + 350);
    const auto *slf_352 = buffer.data(slf + 352);
    const auto *slf_356 = buffer.data(slf + 356);
    const auto *slf_357 = buffer.data(slf + 357);
    const auto *slf_358 = buffer.data(slf + 358);
    const auto *slf_359 = buffer.data(slf + 359);
    const auto *slf_360 = buffer.data(slf + 360);
    const auto *slf_362 = buffer.data(slf + 362);
    const auto *slf_363 = buffer.data(slf + 363);
    const auto *slf_365 = buffer.data(slf + 365);
    const auto *slf_366 = buffer.data(slf + 366);
    const auto *slf_367 = buffer.data(slf + 367);
    const auto *slf_368 = buffer.data(slf + 368);
    const auto *slf_369 = buffer.data(slf + 369);
    const auto *slf_370 = buffer.data(slf + 370);
    const auto *slf_372 = buffer.data(slf + 372);
    const auto *slf_375 = buffer.data(slf + 375);
    const auto *slf_376 = buffer.data(slf + 376);
    const auto *slf_377 = buffer.data(slf + 377);
    const auto *slf_378 = buffer.data(slf + 378);
    const auto *slf_379 = buffer.data(slf + 379);
    const auto *slf_380 = buffer.data(slf + 380);
    const auto *slf_382 = buffer.data(slf + 382);
    const auto *slf_383 = buffer.data(slf + 383);
    const auto *slf_385 = buffer.data(slf + 385);
    const auto *slf_386 = buffer.data(slf + 386);
    const auto *slf_387 = buffer.data(slf + 387);
    const auto *slf_388 = buffer.data(slf + 388);
    const auto *slf_389 = buffer.data(slf + 389);
    const auto *slf_390 = buffer.data(slf + 390);
    const auto *slf_392 = buffer.data(slf + 392);
    const auto *slf_393 = buffer.data(slf + 393);
    const auto *slf_395 = buffer.data(slf + 395);
    const auto *slf_396 = buffer.data(slf + 396);
    const auto *slf_397 = buffer.data(slf + 397);
    const auto *slf_398 = buffer.data(slf + 398);
    const auto *slf_399 = buffer.data(slf + 399);
    const auto *slf_400 = buffer.data(slf + 400);
    const auto *slf_402 = buffer.data(slf + 402);
    const auto *slf_403 = buffer.data(slf + 403);
    const auto *slf_405 = buffer.data(slf + 405);
    const auto *slf_406 = buffer.data(slf + 406);
    const auto *slf_407 = buffer.data(slf + 407);
    const auto *slf_408 = buffer.data(slf + 408);
    const auto *slf_409 = buffer.data(slf + 409);

#pragma omp simd aligned(t_482, t_483, t_484, pb_x, pc_x, pc_y, pc_z, skg0_483, skf_240, \
                         skf_252, skf_323, skg1_483, slf_320, slf_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_13 * skf_240[k]
                   + f_3 * pc_z[k] * slf_320[k];

        t_483[k] = pb_x[k] * skg0_483[k]
                   + f_8 * skf_323[k]
                   - f_6 * pc_x[k] * skg1_483[k];

        t_484[k] = f_12 * skf_252[k]
                   + f_3 * pc_y[k] * slf_322[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pb_x, pc_x, skg0_485, skf_325, skf_326, \
                         skf_327, skf_328, skg1_485, slf_326, slf_327, \
                         slf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pb_x[k] * skg0_485[k]
                   + f_8 * skf_325[k]
                   - f_6 * pc_x[k] * skg1_485[k];

        t_486[k] = f_7 * skf_326[k]
                   + f_3 * pc_x[k] * slf_326[k];

        t_487[k] = f_7 * skf_327[k]
                   + f_3 * pc_x[k] * slf_327[k];

        t_488[k] = f_7 * skf_328[k]
                   + f_3 * pc_x[k] * slf_328[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, pb_x, pc_x, pc_z, skg0_490, skg0_492, \
                         skf_246, skf_329, skg1_490, skg1_492, slf_326, \
                         slf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_7 * skf_329[k]
                   + f_3 * pc_x[k] * slf_329[k];

        t_490[k] = pb_x[k] * skg0_490[k]
                   - f_6 * pc_x[k] * skg1_490[k];

        t_491[k] = f_13 * skf_246[k]
                   + f_3 * pc_z[k] * slf_326[k];

        t_492[k] = pb_x[k] * skg0_492[k]
                   - f_6 * pc_x[k] * skg1_492[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, t_496, pb_x, pc_x, pc_y, skg0_494, skg0_495, \
                         skf_259, skf_260, skf_330, skg1_494, skg1_495, slf_329, \
                         slf_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = f_12 * skf_259[k]
                   + f_3 * pc_y[k] * slf_329[k];

        t_494[k] = pb_x[k] * skg0_494[k]
                   - f_6 * pc_x[k] * skg1_494[k];

        t_495[k] = pb_x[k] * skg0_495[k]
                   + f_13 * skf_330[k]
                   - f_6 * pc_x[k] * skg1_495[k];

        t_496[k] = f_8 * skf_260[k]
                   + f_3 * pc_y[k] * slf_330[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, pb_x, pc_x, pc_y, pc_z, skg0_498, skf_250, \
                         skf_262, skf_333, skg1_498, slf_330, slf_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = f_11 * skf_250[k]
                   + f_3 * pc_z[k] * slf_330[k];

        t_498[k] = pb_x[k] * skg0_498[k]
                   + f_8 * skf_333[k]
                   - f_6 * pc_x[k] * skg1_498[k];

        t_499[k] = f_8 * skf_262[k]
                   + f_3 * pc_y[k] * slf_332[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pb_x, pc_x, skg0_500, skf_335, skf_336, \
                         skf_337, skf_338, skg1_500, slf_336, slf_337, \
                         slf_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = pb_x[k] * skg0_500[k]
                   + f_8 * skf_335[k]
                   - f_6 * pc_x[k] * skg1_500[k];

        t_501[k] = f_7 * skf_336[k]
                   + f_3 * pc_x[k] * slf_336[k];

        t_502[k] = f_7 * skf_337[k]
                   + f_3 * pc_x[k] * slf_337[k];

        t_503[k] = f_7 * skf_338[k]
                   + f_3 * pc_x[k] * slf_338[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pb_x, pc_x, pc_z, skg0_505, skg0_507, \
                         skf_256, skf_339, skg1_505, skg1_507, slf_336, \
                         slf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_7 * skf_339[k]
                   + f_3 * pc_x[k] * slf_339[k];

        t_505[k] = pb_x[k] * skg0_505[k]
                   - f_6 * pc_x[k] * skg1_505[k];

        t_506[k] = f_11 * skf_256[k]
                   + f_3 * pc_z[k] * slf_336[k];

        t_507[k] = pb_x[k] * skg0_507[k]
                   - f_6 * pc_x[k] * skg1_507[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pb_x, pb_y, pc_x, pc_y, skg0_405, \
                         skg0_509, skf_269, skf_270, skg1_405, skg1_509, slf_339, \
                         slf_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_8 * skf_269[k]
                   + f_3 * pc_y[k] * slf_339[k];

        t_509[k] = pb_x[k] * skg0_509[k]
                   - f_6 * pc_x[k] * skg1_509[k];

        t_510[k] = pb_y[k] * skg0_405[k]
                   - f_6 * pc_y[k] * skg1_405[k];

        t_511[k] = f_7 * skf_270[k]
                   + f_3 * pc_y[k] * slf_340[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pb_x, pc_x, pc_y, pc_z, skg0_513, skf_260, \
                         skf_272, skf_343, skg1_513, slf_340, slf_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_10 * skf_260[k]
                   + f_3 * pc_z[k] * slf_340[k];

        t_513[k] = pb_x[k] * skg0_513[k]
                   + f_8 * skf_343[k]
                   - f_6 * pc_x[k] * skg1_513[k];

        t_514[k] = f_7 * skf_272[k]
                   + f_3 * pc_y[k] * slf_342[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, pb_y, pc_x, pc_y, skg0_410, skf_346, \
                         skf_347, skf_348, skg1_410, slf_346, slf_347, \
                         slf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = pb_y[k] * skg0_410[k]
                   - f_6 * pc_y[k] * skg1_410[k];

        t_516[k] = f_7 * skf_346[k]
                   + f_3 * pc_x[k] * slf_346[k];

        t_517[k] = f_7 * skf_347[k]
                   + f_3 * pc_x[k] * slf_347[k];

        t_518[k] = f_7 * skf_348[k]
                   + f_3 * pc_x[k] * slf_348[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, pb_x, pc_x, pc_z, skg0_520, skg0_522, \
                         skf_266, skf_349, skg1_520, skg1_522, slf_346, \
                         slf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_7 * skf_349[k]
                   + f_3 * pc_x[k] * slf_349[k];

        t_520[k] = pb_x[k] * skg0_520[k]
                   - f_6 * pc_x[k] * skg1_520[k];

        t_521[k] = f_10 * skf_266[k]
                   + f_3 * pc_z[k] * slf_346[k];

        t_522[k] = pb_x[k] * skg0_522[k]
                   - f_6 * pc_x[k] * skg1_522[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pb_x, pc_x, pc_y, skg0_524, skg0_525, \
                         skf_279, skf_350, skg1_524, skg1_525, slf_349, \
                         slf_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_7 * skf_279[k]
                   + f_3 * pc_y[k] * slf_349[k];

        t_524[k] = pb_x[k] * skg0_524[k]
                   - f_6 * pc_x[k] * skg1_524[k];

        t_525[k] = pb_x[k] * skg0_525[k]
                   + f_13 * skf_350[k]
                   - f_6 * pc_x[k] * skg1_525[k];

        t_526[k] = f_3 * pc_y[k] * slf_350[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, pb_x, pc_x, pc_y, pc_z, skg0_528, skf_270, \
                         skf_353, skg1_528, slf_350, slf_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_9 * skf_270[k]
                   + f_3 * pc_z[k] * slf_350[k];

        t_528[k] = pb_x[k] * skg0_528[k]
                   + f_8 * skf_353[k]
                   - f_6 * pc_x[k] * skg1_528[k];

        t_529[k] = f_3 * pc_y[k] * slf_352[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pb_x, pc_x, skg0_530, skf_355, skf_356, \
                         skf_357, skf_358, skg1_530, slf_356, slf_357, \
                         slf_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = pb_x[k] * skg0_530[k]
                   + f_8 * skf_355[k]
                   - f_6 * pc_x[k] * skg1_530[k];

        t_531[k] = f_7 * skf_356[k]
                   + f_3 * pc_x[k] * slf_356[k];

        t_532[k] = f_7 * skf_357[k]
                   + f_3 * pc_x[k] * slf_357[k];

        t_533[k] = f_7 * skf_358[k]
                   + f_3 * pc_x[k] * slf_358[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, pb_x, pc_x, pc_z, skg0_535, skg0_537, \
                         skf_276, skf_359, skg1_535, skg1_537, slf_356, \
                         slf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_7 * skf_359[k]
                   + f_3 * pc_x[k] * slf_359[k];

        t_535[k] = pb_x[k] * skg0_535[k]
                   - f_6 * pc_x[k] * skg1_535[k];

        t_536[k] = f_9 * skf_276[k]
                   + f_3 * pc_z[k] * slf_356[k];

        t_537[k] = pb_x[k] * skg0_537[k]
                   - f_6 * pc_x[k] * skg1_537[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, t_542, pb_x, pc_x, pc_y, pc_z, skg0_539, \
                         skf_280, skg1_539, sld0_216, sld1_216, slf_359, \
                         slf_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_3 * pc_y[k] * slf_359[k];

        t_539[k] = pb_x[k] * skg0_539[k]
                   - f_6 * pc_x[k] * skg1_539[k];

        t_540[k] = f_1 * sld0_216[k]
                   - f_2 * sld1_216[k]
                   + f_3 * pc_x[k] * slf_360[k];

        t_541[k] = f_0 * skf_280[k]
                   + f_3 * pc_y[k] * slf_360[k];

        t_542[k] = f_3 * pc_z[k] * slf_360[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, pc_x, pc_y, skf_282, sld0_219, sld0_221, \
                         sld1_219, sld1_221, slf_362, slf_363, slf_365, \
                         slf_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_4 * sld0_219[k]
                   - f_5 * sld1_219[k]
                   + f_3 * pc_x[k] * slf_363[k];

        t_544[k] = f_0 * skf_282[k]
                   + f_3 * pc_y[k] * slf_362[k];

        t_545[k] = f_4 * sld0_221[k]
                   - f_5 * sld1_221[k]
                   + f_3 * pc_x[k] * slf_365[k];

        t_546[k] = f_3 * pc_x[k] * slf_366[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, pc_x, pc_y, pc_z, skf_286, \
                         sld0_219, sld1_219, slf_366, slf_367, slf_368, \
                         slf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_3 * pc_x[k] * slf_367[k];

        t_548[k] = f_3 * pc_x[k] * slf_368[k];

        t_549[k] = f_3 * pc_x[k] * slf_369[k];

        t_550[k] = f_0 * skf_286[k]
                   + f_1 * sld0_219[k]
                   - f_2 * sld1_219[k]
                   + f_3 * pc_y[k] * slf_366[k];

        t_551[k] = f_3 * pc_z[k] * slf_366[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pb_z, pc_y, pc_z, skg0_420, skf_288, \
                         skf_289, skg1_420, sld0_221, sld1_221, slf_368, \
                         slf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_0 * skf_288[k]
                   + f_4 * sld0_221[k]
                   - f_5 * sld1_221[k]
                   + f_3 * pc_y[k] * slf_368[k];

        t_553[k] = f_0 * skf_289[k]
                   + f_3 * pc_y[k] * slf_369[k];

        t_554[k] = f_1 * sld0_221[k]
                   - f_2 * sld1_221[k]
                   + f_3 * pc_z[k] * slf_369[k];

        t_555[k] = pb_z[k] * skg0_420[k]
                   - f_6 * pc_z[k] * skg1_420[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, pb_z, pc_y, pc_z, skg0_423, skf_280, \
                         skf_290, skf_292, skg1_423, slf_370, slf_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_9 * skf_290[k]
                   + f_3 * pc_y[k] * slf_370[k];

        t_557[k] = f_7 * skf_280[k]
                   + f_3 * pc_z[k] * slf_370[k];

        t_558[k] = pb_z[k] * skg0_423[k]
                   - f_6 * pc_z[k] * skg1_423[k];

        t_559[k] = f_9 * skf_292[k]
                   + f_3 * pc_y[k] * slf_372[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, pc_x, sld0_227, sld1_227, slf_375, \
                         slf_376, slf_377, slf_378, slf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_4 * sld0_227[k]
                   - f_5 * sld1_227[k]
                   + f_3 * pc_x[k] * slf_375[k];

        t_561[k] = f_3 * pc_x[k] * slf_376[k];

        t_562[k] = f_3 * pc_x[k] * slf_377[k];

        t_563[k] = f_3 * pc_x[k] * slf_378[k];

        t_564[k] = f_3 * pc_x[k] * slf_379[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, pb_z, pc_y, pc_z, skg0_430, skg0_432, \
                         skf_286, skf_287, skf_299, skg1_430, skg1_432, slf_376, \
                         slf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = pb_z[k] * skg0_430[k]
                   - f_6 * pc_z[k] * skg1_430[k];

        t_566[k] = f_7 * skf_286[k]
                   + f_3 * pc_z[k] * slf_376[k];

        t_567[k] = pb_z[k] * skg0_432[k]
                   + f_8 * skf_287[k]
                   - f_6 * pc_z[k] * skg1_432[k];

        t_568[k] = f_9 * skf_299[k]
                   + f_3 * pc_y[k] * slf_379[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, t_572, pc_x, pc_y, pc_z, skf_289, skf_290, \
                         skf_300, sld0_227, sld0_228, sld1_227, sld1_228, slf_379, \
                         slf_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_7 * skf_289[k]
                   + f_1 * sld0_227[k]
                   - f_2 * sld1_227[k]
                   + f_3 * pc_z[k] * slf_379[k];

        t_570[k] = f_1 * sld0_228[k]
                   - f_2 * sld1_228[k]
                   + f_3 * pc_x[k] * slf_380[k];

        t_571[k] = f_10 * skf_300[k]
                   + f_3 * pc_y[k] * slf_380[k];

        t_572[k] = f_8 * skf_290[k]
                   + f_3 * pc_z[k] * slf_380[k];
    }

#pragma omp simd aligned(t_573, t_574, t_575, t_576, pc_x, pc_y, skf_302, sld0_231, sld0_233, \
                         sld1_231, sld1_233, slf_382, slf_383, slf_385, \
                         slf_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_573[k] = f_4 * sld0_231[k]
                   - f_5 * sld1_231[k]
                   + f_3 * pc_x[k] * slf_383[k];

        t_574[k] = f_10 * skf_302[k]
                   + f_3 * pc_y[k] * slf_382[k];

        t_575[k] = f_4 * sld0_233[k]
                   - f_5 * sld1_233[k]
                   + f_3 * pc_x[k] * slf_385[k];

        t_576[k] = f_3 * pc_x[k] * slf_386[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, pc_x, pc_y, pc_z, skf_296, \
                         skf_306, sld0_231, sld1_231, slf_386, slf_387, slf_388, \
                         slf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = f_3 * pc_x[k] * slf_387[k];

        t_578[k] = f_3 * pc_x[k] * slf_388[k];

        t_579[k] = f_3 * pc_x[k] * slf_389[k];

        t_580[k] = f_10 * skf_306[k]
                   + f_1 * sld0_231[k]
                   - f_2 * sld1_231[k]
                   + f_3 * pc_y[k] * slf_386[k];

        t_581[k] = f_8 * skf_296[k]
                   + f_3 * pc_z[k] * slf_386[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, pc_y, pc_z, skf_299, skf_308, skf_309, sld0_233, \
                         sld1_233, slf_388, slf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_10 * skf_308[k]
                   + f_4 * sld0_233[k]
                   - f_5 * sld1_233[k]
                   + f_3 * pc_y[k] * slf_388[k];

        t_583[k] = f_10 * skf_309[k]
                   + f_3 * pc_y[k] * slf_389[k];

        t_584[k] = f_8 * skf_299[k]
                   + f_1 * sld0_233[k]
                   - f_2 * sld1_233[k]
                   + f_3 * pc_z[k] * slf_389[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, pc_x, pc_y, pc_z, skf_300, skf_310, \
                         sld0_234, sld0_237, sld1_234, sld1_237, slf_390, \
                         slf_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_1 * sld0_234[k]
                   - f_2 * sld1_234[k]
                   + f_3 * pc_x[k] * slf_390[k];

        t_586[k] = f_11 * skf_310[k]
                   + f_3 * pc_y[k] * slf_390[k];

        t_587[k] = f_12 * skf_300[k]
                   + f_3 * pc_z[k] * slf_390[k];

        t_588[k] = f_4 * sld0_237[k]
                   - f_5 * sld1_237[k]
                   + f_3 * pc_x[k] * slf_393[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, t_593, pc_x, pc_y, skf_312, sld0_239, \
                         sld1_239, slf_392, slf_395, slf_396, slf_397, \
                         slf_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = f_11 * skf_312[k]
                   + f_3 * pc_y[k] * slf_392[k];

        t_590[k] = f_4 * sld0_239[k]
                   - f_5 * sld1_239[k]
                   + f_3 * pc_x[k] * slf_395[k];

        t_591[k] = f_3 * pc_x[k] * slf_396[k];

        t_592[k] = f_3 * pc_x[k] * slf_397[k];

        t_593[k] = f_3 * pc_x[k] * slf_398[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, pc_x, pc_y, pc_z, skf_306, skf_316, sld0_237, \
                         sld1_237, slf_396, slf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_3 * pc_x[k] * slf_399[k];

        t_595[k] = f_11 * skf_316[k]
                   + f_1 * sld0_237[k]
                   - f_2 * sld1_237[k]
                   + f_3 * pc_y[k] * slf_396[k];

        t_596[k] = f_12 * skf_306[k]
                   + f_3 * pc_z[k] * slf_396[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, pc_y, pc_z, skf_309, skf_318, skf_319, sld0_239, \
                         sld1_239, slf_398, slf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = f_11 * skf_318[k]
                   + f_4 * sld0_239[k]
                   - f_5 * sld1_239[k]
                   + f_3 * pc_y[k] * slf_398[k];

        t_598[k] = f_11 * skf_319[k]
                   + f_3 * pc_y[k] * slf_399[k];

        t_599[k] = f_12 * skf_309[k]
                   + f_1 * sld0_239[k]
                   - f_2 * sld1_239[k]
                   + f_3 * pc_z[k] * slf_399[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, pc_z, skf_310, skf_320, \
                         sld0_240, sld0_243, sld1_240, sld1_243, slf_400, \
                         slf_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_1 * sld0_240[k]
                   - f_2 * sld1_240[k]
                   + f_3 * pc_x[k] * slf_400[k];

        t_601[k] = f_13 * skf_320[k]
                   + f_3 * pc_y[k] * slf_400[k];

        t_602[k] = f_13 * skf_310[k]
                   + f_3 * pc_z[k] * slf_400[k];

        t_603[k] = f_4 * sld0_243[k]
                   - f_5 * sld1_243[k]
                   + f_3 * pc_x[k] * slf_403[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, pc_x, pc_y, skf_322, sld0_245, \
                         sld1_245, slf_402, slf_405, slf_406, slf_407, \
                         slf_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_13 * skf_322[k]
                   + f_3 * pc_y[k] * slf_402[k];

        t_605[k] = f_4 * sld0_245[k]
                   - f_5 * sld1_245[k]
                   + f_3 * pc_x[k] * slf_405[k];

        t_606[k] = f_3 * pc_x[k] * slf_406[k];

        t_607[k] = f_3 * pc_x[k] * slf_407[k];

        t_608[k] = f_3 * pc_x[k] * slf_408[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pc_x, pc_y, pc_z, skf_316, skf_326, sld0_243, \
                         sld1_243, slf_406, slf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_3 * pc_x[k] * slf_409[k];

        t_610[k] = f_13 * skf_326[k]
                   + f_1 * sld0_243[k]
                   - f_2 * sld1_243[k]
                   + f_3 * pc_y[k] * slf_406[k];

        t_611[k] = f_13 * skf_316[k]
                   + f_3 * pc_z[k] * slf_406[k];
    }
}

static auto
compute_prim_slg_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skg0,
                                                          const size_t skf, const size_t skg1,
                                                          const size_t sld0, const size_t sld1,
                                                          const size_t slf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 3.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.0 / q;

    auto *t_612 = buffer.data(target + 612);
    auto *t_613 = buffer.data(target + 613);
    auto *t_614 = buffer.data(target + 614);
    auto *t_615 = buffer.data(target + 615);
    auto *t_616 = buffer.data(target + 616);
    auto *t_617 = buffer.data(target + 617);
    auto *t_618 = buffer.data(target + 618);
    auto *t_619 = buffer.data(target + 619);
    auto *t_620 = buffer.data(target + 620);
    auto *t_621 = buffer.data(target + 621);
    auto *t_622 = buffer.data(target + 622);
    auto *t_623 = buffer.data(target + 623);
    auto *t_624 = buffer.data(target + 624);
    auto *t_625 = buffer.data(target + 625);
    auto *t_626 = buffer.data(target + 626);
    auto *t_627 = buffer.data(target + 627);
    auto *t_628 = buffer.data(target + 628);
    auto *t_629 = buffer.data(target + 629);
    auto *t_630 = buffer.data(target + 630);
    auto *t_631 = buffer.data(target + 631);
    auto *t_632 = buffer.data(target + 632);
    auto *t_633 = buffer.data(target + 633);
    auto *t_634 = buffer.data(target + 634);
    auto *t_635 = buffer.data(target + 635);
    auto *t_636 = buffer.data(target + 636);
    auto *t_637 = buffer.data(target + 637);
    auto *t_638 = buffer.data(target + 638);
    auto *t_639 = buffer.data(target + 639);
    auto *t_640 = buffer.data(target + 640);
    auto *t_641 = buffer.data(target + 641);
    auto *t_642 = buffer.data(target + 642);
    auto *t_643 = buffer.data(target + 643);
    auto *t_644 = buffer.data(target + 644);
    auto *t_645 = buffer.data(target + 645);
    auto *t_646 = buffer.data(target + 646);
    auto *t_647 = buffer.data(target + 647);
    auto *t_648 = buffer.data(target + 648);
    auto *t_649 = buffer.data(target + 649);
    auto *t_650 = buffer.data(target + 650);
    auto *t_651 = buffer.data(target + 651);
    auto *t_652 = buffer.data(target + 652);
    auto *t_653 = buffer.data(target + 653);
    auto *t_654 = buffer.data(target + 654);
    auto *t_655 = buffer.data(target + 655);
    auto *t_656 = buffer.data(target + 656);
    auto *t_657 = buffer.data(target + 657);
    auto *t_658 = buffer.data(target + 658);
    auto *t_659 = buffer.data(target + 659);
    auto *t_660 = buffer.data(target + 660);
    auto *t_661 = buffer.data(target + 661);
    auto *t_662 = buffer.data(target + 662);
    auto *t_663 = buffer.data(target + 663);
    auto *t_664 = buffer.data(target + 664);
    auto *t_665 = buffer.data(target + 665);
    auto *t_666 = buffer.data(target + 666);
    auto *t_667 = buffer.data(target + 667);
    auto *t_668 = buffer.data(target + 668);
    auto *t_669 = buffer.data(target + 669);
    auto *t_670 = buffer.data(target + 670);
    auto *t_671 = buffer.data(target + 671);
    auto *t_672 = buffer.data(target + 672);
    auto *t_673 = buffer.data(target + 673);
    auto *t_674 = buffer.data(target + 674);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skg0_525 = buffer.data(skg0 + 525);
    const auto *skg0_530 = buffer.data(skg0 + 530);
    const auto *skg0_535 = buffer.data(skg0 + 535);
    const auto *skg0_537 = buffer.data(skg0 + 537);
    const auto *skg0_539 = buffer.data(skg0 + 539);

    const auto *skf_319 = buffer.data(skf + 319);
    const auto *skf_320 = buffer.data(skf + 320);
    const auto *skf_326 = buffer.data(skf + 326);
    const auto *skf_328 = buffer.data(skf + 328);
    const auto *skf_329 = buffer.data(skf + 329);
    const auto *skf_330 = buffer.data(skf + 330);
    const auto *skf_332 = buffer.data(skf + 332);
    const auto *skf_336 = buffer.data(skf + 336);
    const auto *skf_338 = buffer.data(skf + 338);
    const auto *skf_339 = buffer.data(skf + 339);
    const auto *skf_340 = buffer.data(skf + 340);
    const auto *skf_342 = buffer.data(skf + 342);
    const auto *skf_346 = buffer.data(skf + 346);
    const auto *skf_348 = buffer.data(skf + 348);
    const auto *skf_349 = buffer.data(skf + 349);
    const auto *skf_350 = buffer.data(skf + 350);
    const auto *skf_352 = buffer.data(skf + 352);
    const auto *skf_356 = buffer.data(skf + 356);
    const auto *skf_358 = buffer.data(skf + 358);
    const auto *skf_359 = buffer.data(skf + 359);

    const auto *skg1_525 = buffer.data(skg1 + 525);
    const auto *skg1_530 = buffer.data(skg1 + 530);
    const auto *skg1_535 = buffer.data(skg1 + 535);
    const auto *skg1_537 = buffer.data(skg1 + 537);
    const auto *skg1_539 = buffer.data(skg1 + 539);

    const auto *sld0_245 = buffer.data(sld0 + 245);
    const auto *sld0_246 = buffer.data(sld0 + 246);
    const auto *sld0_249 = buffer.data(sld0 + 249);
    const auto *sld0_251 = buffer.data(sld0 + 251);
    const auto *sld0_252 = buffer.data(sld0 + 252);
    const auto *sld0_255 = buffer.data(sld0 + 255);
    const auto *sld0_257 = buffer.data(sld0 + 257);
    const auto *sld0_261 = buffer.data(sld0 + 261);
    const auto *sld0_264 = buffer.data(sld0 + 264);
    const auto *sld0_267 = buffer.data(sld0 + 267);
    const auto *sld0_269 = buffer.data(sld0 + 269);

    const auto *sld1_245 = buffer.data(sld1 + 245);
    const auto *sld1_246 = buffer.data(sld1 + 246);
    const auto *sld1_249 = buffer.data(sld1 + 249);
    const auto *sld1_251 = buffer.data(sld1 + 251);
    const auto *sld1_252 = buffer.data(sld1 + 252);
    const auto *sld1_255 = buffer.data(sld1 + 255);
    const auto *sld1_257 = buffer.data(sld1 + 257);
    const auto *sld1_261 = buffer.data(sld1 + 261);
    const auto *sld1_264 = buffer.data(sld1 + 264);
    const auto *sld1_267 = buffer.data(sld1 + 267);
    const auto *sld1_269 = buffer.data(sld1 + 269);

    const auto *slf_408 = buffer.data(slf + 408);
    const auto *slf_409 = buffer.data(slf + 409);
    const auto *slf_410 = buffer.data(slf + 410);
    const auto *slf_412 = buffer.data(slf + 412);
    const auto *slf_413 = buffer.data(slf + 413);
    const auto *slf_415 = buffer.data(slf + 415);
    const auto *slf_416 = buffer.data(slf + 416);
    const auto *slf_417 = buffer.data(slf + 417);
    const auto *slf_418 = buffer.data(slf + 418);
    const auto *slf_419 = buffer.data(slf + 419);
    const auto *slf_420 = buffer.data(slf + 420);
    const auto *slf_422 = buffer.data(slf + 422);
    const auto *slf_423 = buffer.data(slf + 423);
    const auto *slf_425 = buffer.data(slf + 425);
    const auto *slf_426 = buffer.data(slf + 426);
    const auto *slf_427 = buffer.data(slf + 427);
    const auto *slf_428 = buffer.data(slf + 428);
    const auto *slf_429 = buffer.data(slf + 429);
    const auto *slf_430 = buffer.data(slf + 430);
    const auto *slf_432 = buffer.data(slf + 432);
    const auto *slf_433 = buffer.data(slf + 433);
    const auto *slf_436 = buffer.data(slf + 436);
    const auto *slf_437 = buffer.data(slf + 437);
    const auto *slf_438 = buffer.data(slf + 438);
    const auto *slf_439 = buffer.data(slf + 439);
    const auto *slf_440 = buffer.data(slf + 440);
    const auto *slf_442 = buffer.data(slf + 442);
    const auto *slf_443 = buffer.data(slf + 443);
    const auto *slf_445 = buffer.data(slf + 445);
    const auto *slf_446 = buffer.data(slf + 446);
    const auto *slf_447 = buffer.data(slf + 447);
    const auto *slf_448 = buffer.data(slf + 448);
    const auto *slf_449 = buffer.data(slf + 449);

#pragma omp simd aligned(t_612, t_613, t_614, pc_y, pc_z, skf_319, skf_328, skf_329, sld0_245, \
                         sld1_245, slf_408, slf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_13 * skf_328[k]
                   + f_4 * sld0_245[k]
                   - f_5 * sld1_245[k]
                   + f_3 * pc_y[k] * slf_408[k];

        t_613[k] = f_13 * skf_329[k]
                   + f_3 * pc_y[k] * slf_409[k];

        t_614[k] = f_13 * skf_319[k]
                   + f_1 * sld0_245[k]
                   - f_2 * sld1_245[k]
                   + f_3 * pc_z[k] * slf_409[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, pc_x, pc_y, pc_z, skf_320, skf_330, \
                         sld0_246, sld0_249, sld1_246, sld1_249, slf_410, \
                         slf_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_1 * sld0_246[k]
                   - f_2 * sld1_246[k]
                   + f_3 * pc_x[k] * slf_410[k];

        t_616[k] = f_12 * skf_330[k]
                   + f_3 * pc_y[k] * slf_410[k];

        t_617[k] = f_11 * skf_320[k]
                   + f_3 * pc_z[k] * slf_410[k];

        t_618[k] = f_4 * sld0_249[k]
                   - f_5 * sld1_249[k]
                   + f_3 * pc_x[k] * slf_413[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, t_623, pc_x, pc_y, skf_332, sld0_251, \
                         sld1_251, slf_412, slf_415, slf_416, slf_417, \
                         slf_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_12 * skf_332[k]
                   + f_3 * pc_y[k] * slf_412[k];

        t_620[k] = f_4 * sld0_251[k]
                   - f_5 * sld1_251[k]
                   + f_3 * pc_x[k] * slf_415[k];

        t_621[k] = f_3 * pc_x[k] * slf_416[k];

        t_622[k] = f_3 * pc_x[k] * slf_417[k];

        t_623[k] = f_3 * pc_x[k] * slf_418[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pc_x, pc_y, pc_z, skf_326, skf_336, sld0_249, \
                         sld1_249, slf_416, slf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_3 * pc_x[k] * slf_419[k];

        t_625[k] = f_12 * skf_336[k]
                   + f_1 * sld0_249[k]
                   - f_2 * sld1_249[k]
                   + f_3 * pc_y[k] * slf_416[k];

        t_626[k] = f_11 * skf_326[k]
                   + f_3 * pc_z[k] * slf_416[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pc_y, pc_z, skf_329, skf_338, skf_339, sld0_251, \
                         sld1_251, slf_418, slf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_12 * skf_338[k]
                   + f_4 * sld0_251[k]
                   - f_5 * sld1_251[k]
                   + f_3 * pc_y[k] * slf_418[k];

        t_628[k] = f_12 * skf_339[k]
                   + f_3 * pc_y[k] * slf_419[k];

        t_629[k] = f_11 * skf_329[k]
                   + f_1 * sld0_251[k]
                   - f_2 * sld1_251[k]
                   + f_3 * pc_z[k] * slf_419[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, pc_y, pc_z, skf_330, skf_340, \
                         sld0_252, sld0_255, sld1_252, sld1_255, slf_420, \
                         slf_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_1 * sld0_252[k]
                   - f_2 * sld1_252[k]
                   + f_3 * pc_x[k] * slf_420[k];

        t_631[k] = f_8 * skf_340[k]
                   + f_3 * pc_y[k] * slf_420[k];

        t_632[k] = f_10 * skf_330[k]
                   + f_3 * pc_z[k] * slf_420[k];

        t_633[k] = f_4 * sld0_255[k]
                   - f_5 * sld1_255[k]
                   + f_3 * pc_x[k] * slf_423[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, pc_x, pc_y, skf_342, sld0_257, \
                         sld1_257, slf_422, slf_425, slf_426, slf_427, \
                         slf_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_8 * skf_342[k]
                   + f_3 * pc_y[k] * slf_422[k];

        t_635[k] = f_4 * sld0_257[k]
                   - f_5 * sld1_257[k]
                   + f_3 * pc_x[k] * slf_425[k];

        t_636[k] = f_3 * pc_x[k] * slf_426[k];

        t_637[k] = f_3 * pc_x[k] * slf_427[k];

        t_638[k] = f_3 * pc_x[k] * slf_428[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pc_x, pc_y, pc_z, skf_336, skf_346, sld0_255, \
                         sld1_255, slf_426, slf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_3 * pc_x[k] * slf_429[k];

        t_640[k] = f_8 * skf_346[k]
                   + f_1 * sld0_255[k]
                   - f_2 * sld1_255[k]
                   + f_3 * pc_y[k] * slf_426[k];

        t_641[k] = f_10 * skf_336[k]
                   + f_3 * pc_z[k] * slf_426[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pb_y, pc_y, pc_z, skg0_525, skf_339, \
                         skf_348, skf_349, skg1_525, sld0_257, sld1_257, slf_428, \
                         slf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_8 * skf_348[k]
                   + f_4 * sld0_257[k]
                   - f_5 * sld1_257[k]
                   + f_3 * pc_y[k] * slf_428[k];

        t_643[k] = f_8 * skf_349[k]
                   + f_3 * pc_y[k] * slf_429[k];

        t_644[k] = f_10 * skf_339[k]
                   + f_1 * sld0_257[k]
                   - f_2 * sld1_257[k]
                   + f_3 * pc_z[k] * slf_429[k];

        t_645[k] = pb_y[k] * skg0_525[k]
                   - f_6 * pc_y[k] * skg1_525[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, pc_x, pc_y, pc_z, skf_340, skf_350, \
                         skf_352, sld0_261, sld1_261, slf_430, slf_432, \
                         slf_433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_7 * skf_350[k]
                   + f_3 * pc_y[k] * slf_430[k];

        t_647[k] = f_9 * skf_340[k]
                   + f_3 * pc_z[k] * slf_430[k];

        t_648[k] = f_4 * sld0_261[k]
                   - f_5 * sld1_261[k]
                   + f_3 * pc_x[k] * slf_433[k];

        t_649[k] = f_7 * skf_352[k]
                   + f_3 * pc_y[k] * slf_432[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, pb_y, pc_x, pc_y, skg0_530, \
                         skg1_530, slf_436, slf_437, slf_438, slf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = pb_y[k] * skg0_530[k]
                   - f_6 * pc_y[k] * skg1_530[k];

        t_651[k] = f_3 * pc_x[k] * slf_436[k];

        t_652[k] = f_3 * pc_x[k] * slf_437[k];

        t_653[k] = f_3 * pc_x[k] * slf_438[k];

        t_654[k] = f_3 * pc_x[k] * slf_439[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, pb_y, pc_y, pc_z, skg0_535, skg0_537, skf_346, \
                         skf_356, skf_358, skg1_535, skg1_537, \
                         slf_436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = pb_y[k] * skg0_535[k]
                   + f_13 * skf_356[k]
                   - f_6 * pc_y[k] * skg1_535[k];

        t_656[k] = f_9 * skf_346[k]
                   + f_3 * pc_z[k] * slf_436[k];

        t_657[k] = pb_y[k] * skg0_537[k]
                   + f_8 * skf_358[k]
                   - f_6 * pc_y[k] * skg1_537[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, pb_y, pc_x, pc_y, skg0_539, skf_359, \
                         skg1_539, sld0_264, sld1_264, slf_439, \
                         slf_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = f_7 * skf_359[k]
                   + f_3 * pc_y[k] * slf_439[k];

        t_659[k] = pb_y[k] * skg0_539[k]
                   - f_6 * pc_y[k] * skg1_539[k];

        t_660[k] = f_1 * sld0_264[k]
                   - f_2 * sld1_264[k]
                   + f_3 * pc_x[k] * slf_440[k];

        t_661[k] = f_3 * pc_y[k] * slf_440[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, pc_x, pc_y, pc_z, skf_350, sld0_267, \
                         sld0_269, sld1_267, sld1_269, slf_440, slf_442, slf_443, \
                         slf_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_0 * skf_350[k]
                   + f_3 * pc_z[k] * slf_440[k];

        t_663[k] = f_4 * sld0_267[k]
                   - f_5 * sld1_267[k]
                   + f_3 * pc_x[k] * slf_443[k];

        t_664[k] = f_3 * pc_y[k] * slf_442[k];

        t_665[k] = f_4 * sld0_269[k]
                   - f_5 * sld1_269[k]
                   + f_3 * pc_x[k] * slf_445[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, t_670, t_671, pc_x, pc_y, pc_z, skf_356, \
                         sld0_267, sld1_267, slf_446, slf_447, slf_448, \
                         slf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_3 * pc_x[k] * slf_446[k];

        t_667[k] = f_3 * pc_x[k] * slf_447[k];

        t_668[k] = f_3 * pc_x[k] * slf_448[k];

        t_669[k] = f_3 * pc_x[k] * slf_449[k];

        t_670[k] = f_1 * sld0_267[k]
                   - f_2 * sld1_267[k]
                   + f_3 * pc_y[k] * slf_446[k];

        t_671[k] = f_0 * skf_356[k]
                   + f_3 * pc_z[k] * slf_446[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, pc_y, pc_z, skf_359, sld0_269, sld1_269, \
                         slf_448, slf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_4 * sld0_269[k]
                   - f_5 * sld1_269[k]
                   + f_3 * pc_y[k] * slf_448[k];

        t_673[k] = f_3 * pc_y[k] * slf_449[k];

        t_674[k] = f_0 * skf_359[k]
                   + f_1 * sld0_269[k]
                   - f_2 * sld1_269[k]
                   + f_3 * pc_z[k] * slf_449[k];
    }
}

auto
compute_prim_slg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t skg0, const size_t skf,
                                                   const size_t skg1, const size_t sld0,
                                                   const size_t sld1, const size_t slf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_slg_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, skg0, skf,
                                                              skg1, sld0, sld1, slf, ncols,
                                                              gamma, p, q);

    compute_prim_slg_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, skg0, skf,
                                                              skg1, sld0, sld1, slf, ncols,
                                                              gamma, p, q);

    compute_prim_slg_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, skg0, skf,
                                                              skg1, sld0, sld1, slf, ncols,
                                                              gamma, p, q);

    compute_prim_slg_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, skg0, skf,
                                                              skg1, sld0, sld1, slf, ncols,
                                                              gamma, p, q);

    compute_prim_slg_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, skg0, skf,
                                                              skg1, sld0, sld1, slf, ncols,
                                                              gamma, p, q);

    compute_prim_slg_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, skg0, skf,
                                                              skg1, sld0, sld1, slf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
