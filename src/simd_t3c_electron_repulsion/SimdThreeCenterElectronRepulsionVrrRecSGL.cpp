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


#include "SimdThreeCenterElectronRepulsionVrrRecSGL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sgl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfl0,
                                                          const size_t sfk, const size_t sfl1,
                                                          const size_t sgi0, const size_t sgi1,
                                                          const size_t sgk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfl0_0 = buffer.data(sfl0 + 0);
    const auto *sfl0_3 = buffer.data(sfl0 + 3);
    const auto *sfl0_5 = buffer.data(sfl0 + 5);
    const auto *sfl0_6 = buffer.data(sfl0 + 6);
    const auto *sfl0_9 = buffer.data(sfl0 + 9);
    const auto *sfl0_10 = buffer.data(sfl0 + 10);
    const auto *sfl0_12 = buffer.data(sfl0 + 12);
    const auto *sfl0_14 = buffer.data(sfl0 + 14);
    const auto *sfl0_15 = buffer.data(sfl0 + 15);
    const auto *sfl0_17 = buffer.data(sfl0 + 17);
    const auto *sfl0_18 = buffer.data(sfl0 + 18);
    const auto *sfl0_20 = buffer.data(sfl0 + 20);
    const auto *sfl0_21 = buffer.data(sfl0 + 21);
    const auto *sfl0_23 = buffer.data(sfl0 + 23);
    const auto *sfl0_24 = buffer.data(sfl0 + 24);
    const auto *sfl0_25 = buffer.data(sfl0 + 25);
    const auto *sfl0_27 = buffer.data(sfl0 + 27);
    const auto *sfl0_44 = buffer.data(sfl0 + 44);

    const auto *sfk_0 = buffer.data(sfk + 0);
    const auto *sfk_1 = buffer.data(sfk + 1);
    const auto *sfk_2 = buffer.data(sfk + 2);
    const auto *sfk_3 = buffer.data(sfk + 3);
    const auto *sfk_5 = buffer.data(sfk + 5);
    const auto *sfk_6 = buffer.data(sfk + 6);
    const auto *sfk_7 = buffer.data(sfk + 7);
    const auto *sfk_8 = buffer.data(sfk + 8);
    const auto *sfk_9 = buffer.data(sfk + 9);
    const auto *sfk_10 = buffer.data(sfk + 10);
    const auto *sfk_11 = buffer.data(sfk + 11);
    const auto *sfk_12 = buffer.data(sfk + 12);
    const auto *sfk_13 = buffer.data(sfk + 13);
    const auto *sfk_14 = buffer.data(sfk + 14);
    const auto *sfk_15 = buffer.data(sfk + 15);
    const auto *sfk_16 = buffer.data(sfk + 16);
    const auto *sfk_17 = buffer.data(sfk + 17);
    const auto *sfk_18 = buffer.data(sfk + 18);
    const auto *sfk_19 = buffer.data(sfk + 19);
    const auto *sfk_20 = buffer.data(sfk + 20);
    const auto *sfk_21 = buffer.data(sfk + 21);
    const auto *sfk_23 = buffer.data(sfk + 23);
    const auto *sfk_24 = buffer.data(sfk + 24);
    const auto *sfk_25 = buffer.data(sfk + 25);
    const auto *sfk_27 = buffer.data(sfk + 27);
    const auto *sfk_28 = buffer.data(sfk + 28);
    const auto *sfk_29 = buffer.data(sfk + 29);
    const auto *sfk_30 = buffer.data(sfk + 30);
    const auto *sfk_31 = buffer.data(sfk + 31);
    const auto *sfk_32 = buffer.data(sfk + 32);
    const auto *sfk_33 = buffer.data(sfk + 33);
    const auto *sfk_34 = buffer.data(sfk + 34);
    const auto *sfk_35 = buffer.data(sfk + 35);
    const auto *sfk_64 = buffer.data(sfk + 64);
    const auto *sfk_65 = buffer.data(sfk + 65);
    const auto *sfk_66 = buffer.data(sfk + 66);
    const auto *sfk_67 = buffer.data(sfk + 67);
    const auto *sfk_68 = buffer.data(sfk + 68);
    const auto *sfk_69 = buffer.data(sfk + 69);
    const auto *sfk_70 = buffer.data(sfk + 70);
    const auto *sfk_71 = buffer.data(sfk + 71);

    const auto *sfl1_0 = buffer.data(sfl1 + 0);
    const auto *sfl1_3 = buffer.data(sfl1 + 3);
    const auto *sfl1_5 = buffer.data(sfl1 + 5);
    const auto *sfl1_6 = buffer.data(sfl1 + 6);
    const auto *sfl1_9 = buffer.data(sfl1 + 9);
    const auto *sfl1_10 = buffer.data(sfl1 + 10);
    const auto *sfl1_12 = buffer.data(sfl1 + 12);
    const auto *sfl1_14 = buffer.data(sfl1 + 14);
    const auto *sfl1_15 = buffer.data(sfl1 + 15);
    const auto *sfl1_17 = buffer.data(sfl1 + 17);
    const auto *sfl1_18 = buffer.data(sfl1 + 18);
    const auto *sfl1_20 = buffer.data(sfl1 + 20);
    const auto *sfl1_21 = buffer.data(sfl1 + 21);
    const auto *sfl1_23 = buffer.data(sfl1 + 23);
    const auto *sfl1_24 = buffer.data(sfl1 + 24);
    const auto *sfl1_25 = buffer.data(sfl1 + 25);
    const auto *sfl1_27 = buffer.data(sfl1 + 27);
    const auto *sfl1_44 = buffer.data(sfl1 + 44);

    const auto *sgi0_0 = buffer.data(sgi0 + 0);
    const auto *sgi0_3 = buffer.data(sgi0 + 3);
    const auto *sgi0_5 = buffer.data(sgi0 + 5);
    const auto *sgi0_6 = buffer.data(sgi0 + 6);
    const auto *sgi0_9 = buffer.data(sgi0 + 9);
    const auto *sgi0_10 = buffer.data(sgi0 + 10);
    const auto *sgi0_12 = buffer.data(sgi0 + 12);
    const auto *sgi0_14 = buffer.data(sgi0 + 14);
    const auto *sgi0_15 = buffer.data(sgi0 + 15);
    const auto *sgi0_17 = buffer.data(sgi0 + 17);
    const auto *sgi0_18 = buffer.data(sgi0 + 18);
    const auto *sgi0_20 = buffer.data(sgi0 + 20);
    const auto *sgi0_21 = buffer.data(sgi0 + 21);
    const auto *sgi0_23 = buffer.data(sgi0 + 23);
    const auto *sgi0_24 = buffer.data(sgi0 + 24);
    const auto *sgi0_25 = buffer.data(sgi0 + 25);
    const auto *sgi0_26 = buffer.data(sgi0 + 26);
    const auto *sgi0_27 = buffer.data(sgi0 + 27);
    const auto *sgi0_49 = buffer.data(sgi0 + 49);
    const auto *sgi0_51 = buffer.data(sgi0 + 51);
    const auto *sgi0_52 = buffer.data(sgi0 + 52);
    const auto *sgi0_53 = buffer.data(sgi0 + 53);
    const auto *sgi0_54 = buffer.data(sgi0 + 54);
    const auto *sgi0_55 = buffer.data(sgi0 + 55);

    const auto *sgi1_0 = buffer.data(sgi1 + 0);
    const auto *sgi1_3 = buffer.data(sgi1 + 3);
    const auto *sgi1_5 = buffer.data(sgi1 + 5);
    const auto *sgi1_6 = buffer.data(sgi1 + 6);
    const auto *sgi1_9 = buffer.data(sgi1 + 9);
    const auto *sgi1_10 = buffer.data(sgi1 + 10);
    const auto *sgi1_12 = buffer.data(sgi1 + 12);
    const auto *sgi1_14 = buffer.data(sgi1 + 14);
    const auto *sgi1_15 = buffer.data(sgi1 + 15);
    const auto *sgi1_17 = buffer.data(sgi1 + 17);
    const auto *sgi1_18 = buffer.data(sgi1 + 18);
    const auto *sgi1_20 = buffer.data(sgi1 + 20);
    const auto *sgi1_21 = buffer.data(sgi1 + 21);
    const auto *sgi1_23 = buffer.data(sgi1 + 23);
    const auto *sgi1_24 = buffer.data(sgi1 + 24);
    const auto *sgi1_25 = buffer.data(sgi1 + 25);
    const auto *sgi1_26 = buffer.data(sgi1 + 26);
    const auto *sgi1_27 = buffer.data(sgi1 + 27);
    const auto *sgi1_49 = buffer.data(sgi1 + 49);
    const auto *sgi1_51 = buffer.data(sgi1 + 51);
    const auto *sgi1_52 = buffer.data(sgi1 + 52);
    const auto *sgi1_53 = buffer.data(sgi1 + 53);
    const auto *sgi1_54 = buffer.data(sgi1 + 54);
    const auto *sgi1_55 = buffer.data(sgi1 + 55);

    const auto *sgk_0 = buffer.data(sgk + 0);
    const auto *sgk_2 = buffer.data(sgk + 2);
    const auto *sgk_3 = buffer.data(sgk + 3);
    const auto *sgk_5 = buffer.data(sgk + 5);
    const auto *sgk_6 = buffer.data(sgk + 6);
    const auto *sgk_9 = buffer.data(sgk + 9);
    const auto *sgk_10 = buffer.data(sgk + 10);
    const auto *sgk_12 = buffer.data(sgk + 12);
    const auto *sgk_14 = buffer.data(sgk + 14);
    const auto *sgk_15 = buffer.data(sgk + 15);
    const auto *sgk_17 = buffer.data(sgk + 17);
    const auto *sgk_18 = buffer.data(sgk + 18);
    const auto *sgk_20 = buffer.data(sgk + 20);
    const auto *sgk_21 = buffer.data(sgk + 21);
    const auto *sgk_23 = buffer.data(sgk + 23);
    const auto *sgk_24 = buffer.data(sgk + 24);
    const auto *sgk_25 = buffer.data(sgk + 25);
    const auto *sgk_27 = buffer.data(sgk + 27);
    const auto *sgk_28 = buffer.data(sgk + 28);
    const auto *sgk_29 = buffer.data(sgk + 29);
    const auto *sgk_30 = buffer.data(sgk + 30);
    const auto *sgk_31 = buffer.data(sgk + 31);
    const auto *sgk_32 = buffer.data(sgk + 32);
    const auto *sgk_33 = buffer.data(sgk + 33);
    const auto *sgk_34 = buffer.data(sgk + 34);
    const auto *sgk_35 = buffer.data(sgk + 35);
    const auto *sgk_36 = buffer.data(sgk + 36);
    const auto *sgk_38 = buffer.data(sgk + 38);
    const auto *sgk_39 = buffer.data(sgk + 39);
    const auto *sgk_41 = buffer.data(sgk + 41);
    const auto *sgk_42 = buffer.data(sgk + 42);
    const auto *sgk_45 = buffer.data(sgk + 45);
    const auto *sgk_46 = buffer.data(sgk + 46);
    const auto *sgk_50 = buffer.data(sgk + 50);
    const auto *sgk_51 = buffer.data(sgk + 51);
    const auto *sgk_56 = buffer.data(sgk + 56);
    const auto *sgk_64 = buffer.data(sgk + 64);
    const auto *sgk_65 = buffer.data(sgk + 65);
    const auto *sgk_66 = buffer.data(sgk + 66);
    const auto *sgk_67 = buffer.data(sgk + 67);
    const auto *sgk_68 = buffer.data(sgk + 68);
    const auto *sgk_69 = buffer.data(sgk + 69);
    const auto *sgk_70 = buffer.data(sgk + 70);
    const auto *sgk_71 = buffer.data(sgk + 71);
    const auto *sgk_72 = buffer.data(sgk + 72);
    const auto *sgk_74 = buffer.data(sgk + 74);
    const auto *sgk_75 = buffer.data(sgk + 75);
    const auto *sgk_77 = buffer.data(sgk + 77);
    const auto *sgk_78 = buffer.data(sgk + 78);
    const auto *sgk_81 = buffer.data(sgk + 81);
    const auto *sgk_82 = buffer.data(sgk + 82);
    const auto *sgk_86 = buffer.data(sgk + 86);
    const auto *sgk_87 = buffer.data(sgk + 87);
    const auto *sgk_92 = buffer.data(sgk + 92);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sfk_0, sfk_3, sgi0_0, sgi0_3, \
                         sgi1_0, sgi1_3, sgk_0, sgk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfk_0[k]
                 + f_1 * sgi0_0[k]
                 - f_2 * sgi1_0[k]
                 + f_3 * pc_x[k] * sgk_0[k];

        t_1[k] = f_3 * pc_y[k] * sgk_0[k];

        t_2[k] = f_3 * pc_z[k] * sgk_0[k];

        t_3[k] = f_0 * sfk_3[k]
                 + f_4 * sgi0_3[k]
                 - f_5 * sgi1_3[k]
                 + f_3 * pc_x[k] * sgk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sfk_5, sfk_6, sgi0_5, sgi0_6, sgi1_5, \
                         sgi1_6, sgk_2, sgk_5, sgk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sgk_2[k];

        t_5[k] = f_0 * sfk_5[k]
                 + f_4 * sgi0_5[k]
                 - f_5 * sgi1_5[k]
                 + f_3 * pc_x[k] * sgk_5[k];

        t_6[k] = f_0 * sfk_6[k]
                 + f_6 * sgi0_6[k]
                 - f_7 * sgi1_6[k]
                 + f_3 * pc_x[k] * sgk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sfk_9, sgi0_9, sgi1_9, sgk_3, sgk_5, \
                         sgk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sgk_3[k];

        t_8[k] = f_3 * pc_y[k] * sgk_5[k];

        t_9[k] = f_0 * sfk_9[k]
                 + f_6 * sgi0_9[k]
                 - f_7 * sgi1_9[k]
                 + f_3 * pc_x[k] * sgk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sfk_10, sfk_12, sgi0_10, sgi0_12, \
                         sgi1_10, sgi1_12, sgk_6, sgk_10, sgk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sfk_10[k]
                  + f_8 * sgi0_10[k]
                  - f_9 * sgi1_10[k]
                  + f_3 * pc_x[k] * sgk_10[k];

        t_11[k] = f_3 * pc_z[k] * sgk_6[k];

        t_12[k] = f_0 * sfk_12[k]
                  + f_8 * sgi0_12[k]
                  - f_9 * sgi1_12[k]
                  + f_3 * pc_x[k] * sgk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, sfk_14, sfk_15, sgi0_14, sgi0_15, \
                         sgi1_14, sgi1_15, sgk_9, sgk_14, sgk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sgk_9[k];

        t_14[k] = f_0 * sfk_14[k]
                  + f_8 * sgi0_14[k]
                  - f_9 * sgi1_14[k]
                  + f_3 * pc_x[k] * sgk_14[k];

        t_15[k] = f_0 * sfk_15[k]
                  + f_10 * sgi0_15[k]
                  - f_11 * sgi1_15[k]
                  + f_3 * pc_x[k] * sgk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, sfk_17, sfk_18, sgi0_17, sgi0_18, \
                         sgi1_17, sgi1_18, sgk_10, sgk_17, sgk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * sgk_10[k];

        t_17[k] = f_0 * sfk_17[k]
                  + f_10 * sgi0_17[k]
                  - f_11 * sgi1_17[k]
                  + f_3 * pc_x[k] * sgk_17[k];

        t_18[k] = f_0 * sfk_18[k]
                  + f_10 * sgi0_18[k]
                  - f_11 * sgi1_18[k]
                  + f_3 * pc_x[k] * sgk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pc_x, pc_y, sfk_20, sfk_21, sgi0_20, sgi0_21, \
                         sgi1_20, sgi1_21, sgk_14, sgk_20, sgk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sgk_14[k];

        t_20[k] = f_0 * sfk_20[k]
                  + f_10 * sgi0_20[k]
                  - f_11 * sgi1_20[k]
                  + f_3 * pc_x[k] * sgk_20[k];

        t_21[k] = f_0 * sfk_21[k]
                  + f_12 * sgi0_21[k]
                  - f_13 * sgi1_21[k]
                  + f_3 * pc_x[k] * sgk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pc_x, pc_z, sfk_23, sfk_24, sgi0_23, sgi0_24, \
                         sgi1_23, sgi1_24, sgk_15, sgk_23, sgk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * pc_z[k] * sgk_15[k];

        t_23[k] = f_0 * sfk_23[k]
                  + f_12 * sgi0_23[k]
                  - f_13 * sgi1_23[k]
                  + f_3 * pc_x[k] * sgk_23[k];

        t_24[k] = f_0 * sfk_24[k]
                  + f_12 * sgi0_24[k]
                  - f_13 * sgi1_24[k]
                  + f_3 * pc_x[k] * sgk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pc_x, pc_y, sfk_25, sfk_27, sgi0_25, sgi0_27, \
                         sgi1_25, sgi1_27, sgk_20, sgk_25, sgk_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * sfk_25[k]
                  + f_12 * sgi0_25[k]
                  - f_13 * sgi1_25[k]
                  + f_3 * pc_x[k] * sgk_25[k];

        t_26[k] = f_3 * pc_y[k] * sgk_20[k];

        t_27[k] = f_0 * sfk_27[k]
                  + f_12 * sgi0_27[k]
                  - f_13 * sgi1_27[k]
                  + f_3 * pc_x[k] * sgk_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pc_x, sfk_28, sfk_29, sfk_30, sfk_31, \
                         sfk_32, sgk_28, sgk_29, sgk_30, sgk_31, \
                         sgk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * sfk_28[k]
                  + f_3 * pc_x[k] * sgk_28[k];

        t_29[k] = f_0 * sfk_29[k]
                  + f_3 * pc_x[k] * sgk_29[k];

        t_30[k] = f_0 * sfk_30[k]
                  + f_3 * pc_x[k] * sgk_30[k];

        t_31[k] = f_0 * sfk_31[k]
                  + f_3 * pc_x[k] * sgk_31[k];

        t_32[k] = f_0 * sfk_32[k]
                  + f_3 * pc_x[k] * sgk_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pc_x, pc_y, sfk_33, sfk_34, sfk_35, sgi0_21, \
                         sgi1_21, sgk_28, sgk_33, sgk_34, sgk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * sfk_33[k]
                  + f_3 * pc_x[k] * sgk_33[k];

        t_34[k] = f_0 * sfk_34[k]
                  + f_3 * pc_x[k] * sgk_34[k];

        t_35[k] = f_0 * sfk_35[k]
                  + f_3 * pc_x[k] * sgk_35[k];

        t_36[k] = f_1 * sgi0_21[k]
                  - f_2 * sgi1_21[k]
                  + f_3 * pc_y[k] * sgk_28[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pc_y, pc_z, sgi0_23, sgi0_24, sgi0_25, \
                         sgi1_23, sgi1_24, sgi1_25, sgk_28, sgk_30, sgk_31, \
                         sgk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * pc_z[k] * sgk_28[k];

        t_38[k] = f_4 * sgi0_23[k]
                  - f_5 * sgi1_23[k]
                  + f_3 * pc_y[k] * sgk_30[k];

        t_39[k] = f_6 * sgi0_24[k]
                  - f_7 * sgi1_24[k]
                  + f_3 * pc_y[k] * sgk_31[k];

        t_40[k] = f_8 * sgi0_25[k]
                  - f_9 * sgi1_25[k]
                  + f_3 * pc_y[k] * sgk_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, sgi0_26, sgi0_27, sgi1_26, \
                         sgi1_27, sgk_33, sgk_34, sgk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * sgi0_26[k]
                  - f_11 * sgi1_26[k]
                  + f_3 * pc_y[k] * sgk_33[k];

        t_42[k] = f_12 * sgi0_27[k]
                  - f_13 * sgi1_27[k]
                  + f_3 * pc_y[k] * sgk_34[k];

        t_43[k] = f_3 * pc_y[k] * sgk_35[k];

        t_44[k] = f_1 * sgi0_27[k]
                  - f_2 * sgi1_27[k]
                  + f_3 * pc_z[k] * sgk_35[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_y, pc_y, pc_z, sfl0_0, sfl0_3, sfk_0, \
                         sfk_1, sfl1_0, sfl1_3, sgk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pb_y[k] * sfl0_0[k]
                  - f_14 * pc_y[k] * sfl1_0[k];

        t_46[k] = f_15 * sfk_0[k]
                  + f_3 * pc_y[k] * sgk_36[k];

        t_47[k] = f_3 * pc_z[k] * sgk_36[k];

        t_48[k] = pb_y[k] * sfl0_3[k]
                  + f_16 * sfk_1[k]
                  - f_14 * pc_y[k] * sfl1_3[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_y, pc_y, pc_z, sfl0_5, sfl0_6, sfk_2, \
                         sfk_3, sfl1_5, sfl1_6, sgk_38, sgk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_15 * sfk_2[k]
                  + f_3 * pc_y[k] * sgk_38[k];

        t_50[k] = pb_y[k] * sfl0_5[k]
                  - f_14 * pc_y[k] * sfl1_5[k];

        t_51[k] = pb_y[k] * sfl0_6[k]
                  + f_17 * sfk_3[k]
                  - f_14 * pc_y[k] * sfl1_6[k];

        t_52[k] = f_3 * pc_z[k] * sgk_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_y, pc_y, pc_z, sfl0_9, sfl0_10, sfk_5, \
                         sfk_6, sfl1_9, sfl1_10, sgk_41, sgk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_15 * sfk_5[k]
                  + f_3 * pc_y[k] * sgk_41[k];

        t_54[k] = pb_y[k] * sfl0_9[k]
                  - f_14 * pc_y[k] * sfl1_9[k];

        t_55[k] = pb_y[k] * sfl0_10[k]
                  + f_0 * sfk_6[k]
                  - f_14 * pc_y[k] * sfl1_10[k];

        t_56[k] = f_3 * pc_z[k] * sgk_42[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, sfl0_12, sfl0_14, sfl0_15, sfk_8, \
                         sfk_9, sfk_10, sfl1_12, sfl1_14, sfl1_15, \
                         sgk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_y[k] * sfl0_12[k]
                  + f_16 * sfk_8[k]
                  - f_14 * pc_y[k] * sfl1_12[k];

        t_58[k] = f_15 * sfk_9[k]
                  + f_3 * pc_y[k] * sgk_45[k];

        t_59[k] = pb_y[k] * sfl0_14[k]
                  - f_14 * pc_y[k] * sfl1_14[k];

        t_60[k] = pb_y[k] * sfl0_15[k]
                  + f_18 * sfk_10[k]
                  - f_14 * pc_y[k] * sfl1_15[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_y, pc_y, pc_z, sfl0_17, sfl0_18, sfk_12, \
                         sfk_13, sfk_14, sfl1_17, sfl1_18, sgk_46, \
                         sgk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * pc_z[k] * sgk_46[k];

        t_62[k] = pb_y[k] * sfl0_17[k]
                  + f_17 * sfk_12[k]
                  - f_14 * pc_y[k] * sfl1_17[k];

        t_63[k] = pb_y[k] * sfl0_18[k]
                  + f_16 * sfk_13[k]
                  - f_14 * pc_y[k] * sfl1_18[k];

        t_64[k] = f_15 * sfk_14[k]
                  + f_3 * pc_y[k] * sgk_50[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_y, pc_z, sfl0_20, sfl0_21, sfl0_23, \
                         sfk_15, sfk_17, sfl1_20, sfl1_21, sfl1_23, \
                         sgk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * sfl0_20[k]
                  - f_14 * pc_y[k] * sfl1_20[k];

        t_66[k] = pb_y[k] * sfl0_21[k]
                  + f_19 * sfk_15[k]
                  - f_14 * pc_y[k] * sfl1_21[k];

        t_67[k] = f_3 * pc_z[k] * sgk_51[k];

        t_68[k] = pb_y[k] * sfl0_23[k]
                  + f_0 * sfk_17[k]
                  - f_14 * pc_y[k] * sfl1_23[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_y, pc_y, sfl0_24, sfl0_25, sfl0_27, \
                         sfk_18, sfk_19, sfk_20, sfl1_24, sfl1_25, sfl1_27, \
                         sgk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_y[k] * sfl0_24[k]
                  + f_17 * sfk_18[k]
                  - f_14 * pc_y[k] * sfl1_24[k];

        t_70[k] = pb_y[k] * sfl0_25[k]
                  + f_16 * sfk_19[k]
                  - f_14 * pc_y[k] * sfl1_25[k];

        t_71[k] = f_15 * sfk_20[k]
                  + f_3 * pc_y[k] * sgk_56[k];

        t_72[k] = pb_y[k] * sfl0_27[k]
                  - f_14 * pc_y[k] * sfl1_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pc_x, sfk_64, sfk_65, sfk_66, sfk_67, \
                         sfk_68, sgk_64, sgk_65, sgk_66, sgk_67, \
                         sgk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_17 * sfk_64[k]
                  + f_3 * pc_x[k] * sgk_64[k];

        t_74[k] = f_17 * sfk_65[k]
                  + f_3 * pc_x[k] * sgk_65[k];

        t_75[k] = f_17 * sfk_66[k]
                  + f_3 * pc_x[k] * sgk_66[k];

        t_76[k] = f_17 * sfk_67[k]
                  + f_3 * pc_x[k] * sgk_67[k];

        t_77[k] = f_17 * sfk_68[k]
                  + f_3 * pc_x[k] * sgk_68[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pc_x, pc_y, sfk_28, sfk_69, sfk_70, sfk_71, \
                         sgi0_49, sgi1_49, sgk_64, sgk_69, sgk_70, \
                         sgk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_17 * sfk_69[k]
                  + f_3 * pc_x[k] * sgk_69[k];

        t_79[k] = f_17 * sfk_70[k]
                  + f_3 * pc_x[k] * sgk_70[k];

        t_80[k] = f_17 * sfk_71[k]
                  + f_3 * pc_x[k] * sgk_71[k];

        t_81[k] = f_15 * sfk_28[k]
                  + f_1 * sgi0_49[k]
                  - f_2 * sgi1_49[k]
                  + f_3 * pc_y[k] * sgk_64[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pc_y, pc_z, sfk_30, sfk_31, sgi0_51, sgi0_52, \
                         sgi1_51, sgi1_52, sgk_64, sgk_66, sgk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * pc_z[k] * sgk_64[k];

        t_83[k] = f_15 * sfk_30[k]
                  + f_4 * sgi0_51[k]
                  - f_5 * sgi1_51[k]
                  + f_3 * pc_y[k] * sgk_66[k];

        t_84[k] = f_15 * sfk_31[k]
                  + f_6 * sgi0_52[k]
                  - f_7 * sgi1_52[k]
                  + f_3 * pc_y[k] * sgk_67[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_y, sfk_32, sfk_33, sfk_34, sgi0_53, sgi0_54, \
                         sgi0_55, sgi1_53, sgi1_54, sgi1_55, sgk_68, sgk_69, \
                         sgk_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_15 * sfk_32[k]
                  + f_8 * sgi0_53[k]
                  - f_9 * sgi1_53[k]
                  + f_3 * pc_y[k] * sgk_68[k];

        t_86[k] = f_15 * sfk_33[k]
                  + f_10 * sgi0_54[k]
                  - f_11 * sgi1_54[k]
                  + f_3 * pc_y[k] * sgk_69[k];

        t_87[k] = f_15 * sfk_34[k]
                  + f_12 * sgi0_55[k]
                  - f_13 * sgi1_55[k]
                  + f_3 * pc_y[k] * sgk_70[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pb_y, pb_z, pc_y, pc_z, sfl0_0, sfl0_44, \
                         sfk_35, sfl1_0, sfl1_44, sgk_71, sgk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * sfk_35[k]
                  + f_3 * pc_y[k] * sgk_71[k];

        t_89[k] = pb_y[k] * sfl0_44[k]
                  - f_14 * pc_y[k] * sfl1_44[k];

        t_90[k] = pb_z[k] * sfl0_0[k]
                  - f_14 * pc_z[k] * sfl1_0[k];

        t_91[k] = f_3 * pc_y[k] * sgk_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_z, pc_y, pc_z, sfl0_3, sfl0_5, sfk_0, \
                         sfk_2, sfl1_3, sfl1_5, sgk_72, sgk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * sfk_0[k]
                  + f_3 * pc_z[k] * sgk_72[k];

        t_93[k] = pb_z[k] * sfl0_3[k]
                  - f_14 * pc_z[k] * sfl1_3[k];

        t_94[k] = f_3 * pc_y[k] * sgk_74[k];

        t_95[k] = pb_z[k] * sfl0_5[k]
                  + f_16 * sfk_2[k]
                  - f_14 * pc_z[k] * sfl1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_z, pc_y, pc_z, sfl0_6, sfl0_9, sfk_3, \
                         sfk_5, sfl1_6, sfl1_9, sgk_75, sgk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pb_z[k] * sfl0_6[k]
                  - f_14 * pc_z[k] * sfl1_6[k];

        t_97[k] = f_15 * sfk_3[k]
                  + f_3 * pc_z[k] * sgk_75[k];

        t_98[k] = f_3 * pc_y[k] * sgk_77[k];

        t_99[k] = pb_z[k] * sfl0_9[k]
                  + f_17 * sfk_5[k]
                  - f_14 * pc_z[k] * sfl1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pb_z, pc_y, pc_z, sfl0_10, sfl0_12, \
                         sfk_6, sfk_7, sfl1_10, sfl1_12, sgk_78, \
                         sgk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pb_z[k] * sfl0_10[k]
                   - f_14 * pc_z[k] * sfl1_10[k];

        t_101[k] = f_15 * sfk_6[k]
                   + f_3 * pc_z[k] * sgk_78[k];

        t_102[k] = pb_z[k] * sfl0_12[k]
                   + f_16 * sfk_7[k]
                   - f_14 * pc_z[k] * sfl1_12[k];

        t_103[k] = f_3 * pc_y[k] * sgk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_z, pc_z, sfl0_14, sfl0_15, sfl0_17, \
                         sfk_9, sfk_10, sfk_11, sfl1_14, sfl1_15, sfl1_17, \
                         sgk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_z[k] * sfl0_14[k]
                   + f_0 * sfk_9[k]
                   - f_14 * pc_z[k] * sfl1_14[k];

        t_105[k] = pb_z[k] * sfl0_15[k]
                   - f_14 * pc_z[k] * sfl1_15[k];

        t_106[k] = f_15 * sfk_10[k]
                   + f_3 * pc_z[k] * sgk_82[k];

        t_107[k] = pb_z[k] * sfl0_17[k]
                   + f_16 * sfk_11[k]
                   - f_14 * pc_z[k] * sfl1_17[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pb_z, pc_y, pc_z, sfl0_18, sfl0_20, \
                         sfl0_21, sfk_12, sfk_14, sfl1_18, sfl1_20, sfl1_21, \
                         sgk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_z[k] * sfl0_18[k]
                   + f_17 * sfk_12[k]
                   - f_14 * pc_z[k] * sfl1_18[k];

        t_109[k] = f_3 * pc_y[k] * sgk_86[k];

        t_110[k] = pb_z[k] * sfl0_20[k]
                   + f_18 * sfk_14[k]
                   - f_14 * pc_z[k] * sfl1_20[k];

        t_111[k] = pb_z[k] * sfl0_21[k]
                   - f_14 * pc_z[k] * sfl1_21[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_z, pc_z, sfl0_23, sfl0_24, sfk_15, sfk_16, \
                         sfk_17, sfl1_23, sfl1_24, sgk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_15 * sfk_15[k]
                   + f_3 * pc_z[k] * sgk_87[k];

        t_113[k] = pb_z[k] * sfl0_23[k]
                   + f_16 * sfk_16[k]
                   - f_14 * pc_z[k] * sfl1_23[k];

        t_114[k] = pb_z[k] * sfl0_24[k]
                   + f_17 * sfk_17[k]
                   - f_14 * pc_z[k] * sfl1_24[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pb_z, pc_y, pc_z, sfl0_25, sfl0_27, sfk_18, \
                         sfk_20, sfl1_25, sfl1_27, sgk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_z[k] * sfl0_25[k]
                   + f_0 * sfk_18[k]
                   - f_14 * pc_z[k] * sfl1_25[k];

        t_116[k] = f_3 * pc_y[k] * sgk_92[k];

        t_117[k] = pb_z[k] * sfl0_27[k]
                   + f_19 * sfk_20[k]
                   - f_14 * pc_z[k] * sfl1_27[k];
    }
}

static auto
compute_prim_sgl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfl0,
                                                          const size_t sfk, const size_t sfl1,
                                                          const size_t sgi0, const size_t sgi1,
                                                          const size_t sgk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfl0_36 = buffer.data(sfl0 + 36);
    const auto *sfl0_48 = buffer.data(sfl0 + 48);
    const auto *sfl0_51 = buffer.data(sfl0 + 51);
    const auto *sfl0_55 = buffer.data(sfl0 + 55);
    const auto *sfl0_60 = buffer.data(sfl0 + 60);
    const auto *sfl0_66 = buffer.data(sfl0 + 66);
    const auto *sfl0_81 = buffer.data(sfl0 + 81);
    const auto *sfl0_90 = buffer.data(sfl0 + 90);
    const auto *sfl0_95 = buffer.data(sfl0 + 95);
    const auto *sfl0_99 = buffer.data(sfl0 + 99);
    const auto *sfl0_102 = buffer.data(sfl0 + 102);
    const auto *sfl0_104 = buffer.data(sfl0 + 104);
    const auto *sfl0_107 = buffer.data(sfl0 + 107);
    const auto *sfl0_108 = buffer.data(sfl0 + 108);
    const auto *sfl0_110 = buffer.data(sfl0 + 110);
    const auto *sfl0_113 = buffer.data(sfl0 + 113);
    const auto *sfl0_114 = buffer.data(sfl0 + 114);
    const auto *sfl0_115 = buffer.data(sfl0 + 115);
    const auto *sfl0_117 = buffer.data(sfl0 + 117);
    const auto *sfl0_134 = buffer.data(sfl0 + 134);

    const auto *sfk_28 = buffer.data(sfk + 28);
    const auto *sfk_35 = buffer.data(sfk + 35);
    const auto *sfk_36 = buffer.data(sfk + 36);
    const auto *sfk_38 = buffer.data(sfk + 38);
    const auto *sfk_39 = buffer.data(sfk + 39);
    const auto *sfk_41 = buffer.data(sfk + 41);
    const auto *sfk_42 = buffer.data(sfk + 42);
    const auto *sfk_45 = buffer.data(sfk + 45);
    const auto *sfk_46 = buffer.data(sfk + 46);
    const auto *sfk_50 = buffer.data(sfk + 50);
    const auto *sfk_51 = buffer.data(sfk + 51);
    const auto *sfk_56 = buffer.data(sfk + 56);
    const auto *sfk_64 = buffer.data(sfk + 64);
    const auto *sfk_66 = buffer.data(sfk + 66);
    const auto *sfk_67 = buffer.data(sfk + 67);
    const auto *sfk_68 = buffer.data(sfk + 68);
    const auto *sfk_69 = buffer.data(sfk + 69);
    const auto *sfk_70 = buffer.data(sfk + 70);
    const auto *sfk_71 = buffer.data(sfk + 71);
    const auto *sfk_72 = buffer.data(sfk + 72);
    const auto *sfk_74 = buffer.data(sfk + 74);
    const auto *sfk_75 = buffer.data(sfk + 75);
    const auto *sfk_77 = buffer.data(sfk + 77);
    const auto *sfk_80 = buffer.data(sfk + 80);
    const auto *sfk_81 = buffer.data(sfk + 81);
    const auto *sfk_84 = buffer.data(sfk + 84);
    const auto *sfk_85 = buffer.data(sfk + 85);
    const auto *sfk_86 = buffer.data(sfk + 86);
    const auto *sfk_89 = buffer.data(sfk + 89);
    const auto *sfk_90 = buffer.data(sfk + 90);
    const auto *sfk_91 = buffer.data(sfk + 91);
    const auto *sfk_92 = buffer.data(sfk + 92);
    const auto *sfk_100 = buffer.data(sfk + 100);
    const auto *sfk_101 = buffer.data(sfk + 101);
    const auto *sfk_102 = buffer.data(sfk + 102);
    const auto *sfk_103 = buffer.data(sfk + 103);
    const auto *sfk_104 = buffer.data(sfk + 104);
    const auto *sfk_105 = buffer.data(sfk + 105);
    const auto *sfk_106 = buffer.data(sfk + 106);
    const auto *sfk_107 = buffer.data(sfk + 107);
    const auto *sfk_108 = buffer.data(sfk + 108);
    const auto *sfk_111 = buffer.data(sfk + 111);
    const auto *sfk_113 = buffer.data(sfk + 113);
    const auto *sfk_114 = buffer.data(sfk + 114);
    const auto *sfk_117 = buffer.data(sfk + 117);
    const auto *sfk_118 = buffer.data(sfk + 118);
    const auto *sfk_120 = buffer.data(sfk + 120);
    const auto *sfk_122 = buffer.data(sfk + 122);
    const auto *sfk_123 = buffer.data(sfk + 123);
    const auto *sfk_125 = buffer.data(sfk + 125);
    const auto *sfk_126 = buffer.data(sfk + 126);
    const auto *sfk_128 = buffer.data(sfk + 128);
    const auto *sfk_129 = buffer.data(sfk + 129);
    const auto *sfk_131 = buffer.data(sfk + 131);
    const auto *sfk_132 = buffer.data(sfk + 132);
    const auto *sfk_133 = buffer.data(sfk + 133);
    const auto *sfk_135 = buffer.data(sfk + 135);
    const auto *sfk_136 = buffer.data(sfk + 136);
    const auto *sfk_137 = buffer.data(sfk + 137);
    const auto *sfk_138 = buffer.data(sfk + 138);
    const auto *sfk_139 = buffer.data(sfk + 139);
    const auto *sfk_140 = buffer.data(sfk + 140);
    const auto *sfk_141 = buffer.data(sfk + 141);
    const auto *sfk_142 = buffer.data(sfk + 142);
    const auto *sfk_143 = buffer.data(sfk + 143);
    const auto *sfk_172 = buffer.data(sfk + 172);
    const auto *sfk_173 = buffer.data(sfk + 173);
    const auto *sfk_174 = buffer.data(sfk + 174);
    const auto *sfk_175 = buffer.data(sfk + 175);
    const auto *sfk_176 = buffer.data(sfk + 176);
    const auto *sfk_177 = buffer.data(sfk + 177);
    const auto *sfk_178 = buffer.data(sfk + 178);
    const auto *sfk_179 = buffer.data(sfk + 179);
    const auto *sfk_180 = buffer.data(sfk + 180);
    const auto *sfk_183 = buffer.data(sfk + 183);
    const auto *sfk_185 = buffer.data(sfk + 185);
    const auto *sfk_186 = buffer.data(sfk + 186);

    const auto *sfl1_36 = buffer.data(sfl1 + 36);
    const auto *sfl1_48 = buffer.data(sfl1 + 48);
    const auto *sfl1_51 = buffer.data(sfl1 + 51);
    const auto *sfl1_55 = buffer.data(sfl1 + 55);
    const auto *sfl1_60 = buffer.data(sfl1 + 60);
    const auto *sfl1_66 = buffer.data(sfl1 + 66);
    const auto *sfl1_81 = buffer.data(sfl1 + 81);
    const auto *sfl1_90 = buffer.data(sfl1 + 90);
    const auto *sfl1_95 = buffer.data(sfl1 + 95);
    const auto *sfl1_99 = buffer.data(sfl1 + 99);
    const auto *sfl1_102 = buffer.data(sfl1 + 102);
    const auto *sfl1_104 = buffer.data(sfl1 + 104);
    const auto *sfl1_107 = buffer.data(sfl1 + 107);
    const auto *sfl1_108 = buffer.data(sfl1 + 108);
    const auto *sfl1_110 = buffer.data(sfl1 + 110);
    const auto *sfl1_113 = buffer.data(sfl1 + 113);
    const auto *sfl1_114 = buffer.data(sfl1 + 114);
    const auto *sfl1_115 = buffer.data(sfl1 + 115);
    const auto *sfl1_117 = buffer.data(sfl1 + 117);
    const auto *sfl1_134 = buffer.data(sfl1 + 134);

    const auto *sgi0_79 = buffer.data(sgi0 + 79);
    const auto *sgi0_80 = buffer.data(sgi0 + 80);
    const auto *sgi0_81 = buffer.data(sgi0 + 81);
    const auto *sgi0_82 = buffer.data(sgi0 + 82);
    const auto *sgi0_83 = buffer.data(sgi0 + 83);
    const auto *sgi0_84 = buffer.data(sgi0 + 84);
    const auto *sgi0_87 = buffer.data(sgi0 + 87);
    const auto *sgi0_89 = buffer.data(sgi0 + 89);
    const auto *sgi0_90 = buffer.data(sgi0 + 90);
    const auto *sgi0_93 = buffer.data(sgi0 + 93);
    const auto *sgi0_94 = buffer.data(sgi0 + 94);
    const auto *sgi0_96 = buffer.data(sgi0 + 96);
    const auto *sgi0_98 = buffer.data(sgi0 + 98);
    const auto *sgi0_99 = buffer.data(sgi0 + 99);
    const auto *sgi0_101 = buffer.data(sgi0 + 101);
    const auto *sgi0_102 = buffer.data(sgi0 + 102);
    const auto *sgi0_104 = buffer.data(sgi0 + 104);
    const auto *sgi0_105 = buffer.data(sgi0 + 105);
    const auto *sgi0_107 = buffer.data(sgi0 + 107);
    const auto *sgi0_108 = buffer.data(sgi0 + 108);
    const auto *sgi0_109 = buffer.data(sgi0 + 109);
    const auto *sgi0_110 = buffer.data(sgi0 + 110);
    const auto *sgi0_111 = buffer.data(sgi0 + 111);
    const auto *sgi0_135 = buffer.data(sgi0 + 135);
    const auto *sgi0_136 = buffer.data(sgi0 + 136);
    const auto *sgi0_137 = buffer.data(sgi0 + 137);
    const auto *sgi0_138 = buffer.data(sgi0 + 138);
    const auto *sgi0_139 = buffer.data(sgi0 + 139);
    const auto *sgi0_140 = buffer.data(sgi0 + 140);
    const auto *sgi0_143 = buffer.data(sgi0 + 143);
    const auto *sgi0_145 = buffer.data(sgi0 + 145);
    const auto *sgi0_146 = buffer.data(sgi0 + 146);

    const auto *sgi1_79 = buffer.data(sgi1 + 79);
    const auto *sgi1_80 = buffer.data(sgi1 + 80);
    const auto *sgi1_81 = buffer.data(sgi1 + 81);
    const auto *sgi1_82 = buffer.data(sgi1 + 82);
    const auto *sgi1_83 = buffer.data(sgi1 + 83);
    const auto *sgi1_84 = buffer.data(sgi1 + 84);
    const auto *sgi1_87 = buffer.data(sgi1 + 87);
    const auto *sgi1_89 = buffer.data(sgi1 + 89);
    const auto *sgi1_90 = buffer.data(sgi1 + 90);
    const auto *sgi1_93 = buffer.data(sgi1 + 93);
    const auto *sgi1_94 = buffer.data(sgi1 + 94);
    const auto *sgi1_96 = buffer.data(sgi1 + 96);
    const auto *sgi1_98 = buffer.data(sgi1 + 98);
    const auto *sgi1_99 = buffer.data(sgi1 + 99);
    const auto *sgi1_101 = buffer.data(sgi1 + 101);
    const auto *sgi1_102 = buffer.data(sgi1 + 102);
    const auto *sgi1_104 = buffer.data(sgi1 + 104);
    const auto *sgi1_105 = buffer.data(sgi1 + 105);
    const auto *sgi1_107 = buffer.data(sgi1 + 107);
    const auto *sgi1_108 = buffer.data(sgi1 + 108);
    const auto *sgi1_109 = buffer.data(sgi1 + 109);
    const auto *sgi1_110 = buffer.data(sgi1 + 110);
    const auto *sgi1_111 = buffer.data(sgi1 + 111);
    const auto *sgi1_135 = buffer.data(sgi1 + 135);
    const auto *sgi1_136 = buffer.data(sgi1 + 136);
    const auto *sgi1_137 = buffer.data(sgi1 + 137);
    const auto *sgi1_138 = buffer.data(sgi1 + 138);
    const auto *sgi1_139 = buffer.data(sgi1 + 139);
    const auto *sgi1_140 = buffer.data(sgi1 + 140);
    const auto *sgi1_143 = buffer.data(sgi1 + 143);
    const auto *sgi1_145 = buffer.data(sgi1 + 145);
    const auto *sgi1_146 = buffer.data(sgi1 + 146);

    const auto *sgk_100 = buffer.data(sgk + 100);
    const auto *sgk_101 = buffer.data(sgk + 101);
    const auto *sgk_102 = buffer.data(sgk + 102);
    const auto *sgk_103 = buffer.data(sgk + 103);
    const auto *sgk_104 = buffer.data(sgk + 104);
    const auto *sgk_105 = buffer.data(sgk + 105);
    const auto *sgk_106 = buffer.data(sgk + 106);
    const auto *sgk_107 = buffer.data(sgk + 107);
    const auto *sgk_108 = buffer.data(sgk + 108);
    const auto *sgk_110 = buffer.data(sgk + 110);
    const auto *sgk_111 = buffer.data(sgk + 111);
    const auto *sgk_113 = buffer.data(sgk + 113);
    const auto *sgk_114 = buffer.data(sgk + 114);
    const auto *sgk_117 = buffer.data(sgk + 117);
    const auto *sgk_118 = buffer.data(sgk + 118);
    const auto *sgk_120 = buffer.data(sgk + 120);
    const auto *sgk_122 = buffer.data(sgk + 122);
    const auto *sgk_123 = buffer.data(sgk + 123);
    const auto *sgk_125 = buffer.data(sgk + 125);
    const auto *sgk_126 = buffer.data(sgk + 126);
    const auto *sgk_128 = buffer.data(sgk + 128);
    const auto *sgk_129 = buffer.data(sgk + 129);
    const auto *sgk_131 = buffer.data(sgk + 131);
    const auto *sgk_132 = buffer.data(sgk + 132);
    const auto *sgk_133 = buffer.data(sgk + 133);
    const auto *sgk_135 = buffer.data(sgk + 135);
    const auto *sgk_136 = buffer.data(sgk + 136);
    const auto *sgk_137 = buffer.data(sgk + 137);
    const auto *sgk_138 = buffer.data(sgk + 138);
    const auto *sgk_139 = buffer.data(sgk + 139);
    const auto *sgk_140 = buffer.data(sgk + 140);
    const auto *sgk_141 = buffer.data(sgk + 141);
    const auto *sgk_142 = buffer.data(sgk + 142);
    const auto *sgk_143 = buffer.data(sgk + 143);
    const auto *sgk_144 = buffer.data(sgk + 144);
    const auto *sgk_146 = buffer.data(sgk + 146);
    const auto *sgk_147 = buffer.data(sgk + 147);
    const auto *sgk_149 = buffer.data(sgk + 149);
    const auto *sgk_150 = buffer.data(sgk + 150);
    const auto *sgk_153 = buffer.data(sgk + 153);
    const auto *sgk_154 = buffer.data(sgk + 154);
    const auto *sgk_158 = buffer.data(sgk + 158);
    const auto *sgk_159 = buffer.data(sgk + 159);
    const auto *sgk_164 = buffer.data(sgk + 164);
    const auto *sgk_172 = buffer.data(sgk + 172);
    const auto *sgk_173 = buffer.data(sgk + 173);
    const auto *sgk_174 = buffer.data(sgk + 174);
    const auto *sgk_175 = buffer.data(sgk + 175);
    const auto *sgk_176 = buffer.data(sgk + 176);
    const auto *sgk_177 = buffer.data(sgk + 177);
    const auto *sgk_178 = buffer.data(sgk + 178);
    const auto *sgk_179 = buffer.data(sgk + 179);
    const auto *sgk_180 = buffer.data(sgk + 180);
    const auto *sgk_182 = buffer.data(sgk + 182);
    const auto *sgk_183 = buffer.data(sgk + 183);
    const auto *sgk_185 = buffer.data(sgk + 185);
    const auto *sgk_186 = buffer.data(sgk + 186);

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pc_x, sfk_100, sfk_101, sfk_102, \
                         sfk_103, sfk_104, sgk_100, sgk_101, sgk_102, sgk_103, \
                         sgk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_17 * sfk_100[k]
                   + f_3 * pc_x[k] * sgk_100[k];

        t_119[k] = f_17 * sfk_101[k]
                   + f_3 * pc_x[k] * sgk_101[k];

        t_120[k] = f_17 * sfk_102[k]
                   + f_3 * pc_x[k] * sgk_102[k];

        t_121[k] = f_17 * sfk_103[k]
                   + f_3 * pc_x[k] * sgk_103[k];

        t_122[k] = f_17 * sfk_104[k]
                   + f_3 * pc_x[k] * sgk_104[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pb_z, pc_x, pc_z, sfl0_36, sfk_105, \
                         sfk_106, sfk_107, sfl1_36, sgk_105, sgk_106, \
                         sgk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_17 * sfk_105[k]
                   + f_3 * pc_x[k] * sgk_105[k];

        t_124[k] = f_17 * sfk_106[k]
                   + f_3 * pc_x[k] * sgk_106[k];

        t_125[k] = f_17 * sfk_107[k]
                   + f_3 * pc_x[k] * sgk_107[k];

        t_126[k] = pb_z[k] * sfl0_36[k]
                   - f_14 * pc_z[k] * sfl1_36[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pc_y, pc_z, sfk_28, sgi0_79, sgi0_80, sgi1_79, \
                         sgi1_80, sgk_100, sgk_102, sgk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_15 * sfk_28[k]
                   + f_3 * pc_z[k] * sgk_100[k];

        t_128[k] = f_4 * sgi0_79[k]
                   - f_5 * sgi1_79[k]
                   + f_3 * pc_y[k] * sgk_102[k];

        t_129[k] = f_6 * sgi0_80[k]
                   - f_7 * sgi1_80[k]
                   + f_3 * pc_y[k] * sgk_103[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_y, sgi0_81, sgi0_82, sgi0_83, sgi1_81, \
                         sgi1_82, sgi1_83, sgk_104, sgk_105, sgk_106, \
                         sgk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_8 * sgi0_81[k]
                   - f_9 * sgi1_81[k]
                   + f_3 * pc_y[k] * sgk_104[k];

        t_131[k] = f_10 * sgi0_82[k]
                   - f_11 * sgi1_82[k]
                   + f_3 * pc_y[k] * sgk_105[k];

        t_132[k] = f_12 * sgi0_83[k]
                   - f_13 * sgi1_83[k]
                   + f_3 * pc_y[k] * sgk_106[k];

        t_133[k] = f_3 * pc_y[k] * sgk_107[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, sfk_35, sfk_36, \
                         sfk_108, sgi0_83, sgi0_84, sgi1_83, sgi1_84, sgk_107, \
                         sgk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_15 * sfk_35[k]
                   + f_1 * sgi0_83[k]
                   - f_2 * sgi1_83[k]
                   + f_3 * pc_z[k] * sgk_107[k];

        t_135[k] = f_16 * sfk_108[k]
                   + f_1 * sgi0_84[k]
                   - f_2 * sgi1_84[k]
                   + f_3 * pc_x[k] * sgk_108[k];

        t_136[k] = f_16 * sfk_36[k]
                   + f_3 * pc_y[k] * sgk_108[k];

        t_137[k] = f_3 * pc_z[k] * sgk_108[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_x, pc_y, sfk_38, sfk_111, sfk_113, sgi0_87, \
                         sgi0_89, sgi1_87, sgi1_89, sgk_110, sgk_111, \
                         sgk_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_16 * sfk_111[k]
                   + f_4 * sgi0_87[k]
                   - f_5 * sgi1_87[k]
                   + f_3 * pc_x[k] * sgk_111[k];

        t_139[k] = f_16 * sfk_38[k]
                   + f_3 * pc_y[k] * sgk_110[k];

        t_140[k] = f_16 * sfk_113[k]
                   + f_4 * sgi0_89[k]
                   - f_5 * sgi1_89[k]
                   + f_3 * pc_x[k] * sgk_113[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_x, pc_y, pc_z, sfk_41, sfk_114, sgi0_90, \
                         sgi1_90, sgk_111, sgk_113, sgk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_16 * sfk_114[k]
                   + f_6 * sgi0_90[k]
                   - f_7 * sgi1_90[k]
                   + f_3 * pc_x[k] * sgk_114[k];

        t_142[k] = f_3 * pc_z[k] * sgk_111[k];

        t_143[k] = f_16 * sfk_41[k]
                   + f_3 * pc_y[k] * sgk_113[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_z, sfk_117, sfk_118, sgi0_93, sgi0_94, \
                         sgi1_93, sgi1_94, sgk_114, sgk_117, sgk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_16 * sfk_117[k]
                   + f_6 * sgi0_93[k]
                   - f_7 * sgi1_93[k]
                   + f_3 * pc_x[k] * sgk_117[k];

        t_145[k] = f_16 * sfk_118[k]
                   + f_8 * sgi0_94[k]
                   - f_9 * sgi1_94[k]
                   + f_3 * pc_x[k] * sgk_118[k];

        t_146[k] = f_3 * pc_z[k] * sgk_114[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_x, pc_y, sfk_45, sfk_120, sfk_122, sgi0_96, \
                         sgi0_98, sgi1_96, sgi1_98, sgk_117, sgk_120, \
                         sgk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_16 * sfk_120[k]
                   + f_8 * sgi0_96[k]
                   - f_9 * sgi1_96[k]
                   + f_3 * pc_x[k] * sgk_120[k];

        t_148[k] = f_16 * sfk_45[k]
                   + f_3 * pc_y[k] * sgk_117[k];

        t_149[k] = f_16 * sfk_122[k]
                   + f_8 * sgi0_98[k]
                   - f_9 * sgi1_98[k]
                   + f_3 * pc_x[k] * sgk_122[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pc_x, pc_z, sfk_123, sfk_125, sgi0_99, sgi0_101, \
                         sgi1_99, sgi1_101, sgk_118, sgk_123, sgk_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_16 * sfk_123[k]
                   + f_10 * sgi0_99[k]
                   - f_11 * sgi1_99[k]
                   + f_3 * pc_x[k] * sgk_123[k];

        t_151[k] = f_3 * pc_z[k] * sgk_118[k];

        t_152[k] = f_16 * sfk_125[k]
                   + f_10 * sgi0_101[k]
                   - f_11 * sgi1_101[k]
                   + f_3 * pc_x[k] * sgk_125[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_x, pc_y, sfk_50, sfk_126, sfk_128, sgi0_102, \
                         sgi0_104, sgi1_102, sgi1_104, sgk_122, sgk_126, \
                         sgk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_16 * sfk_126[k]
                   + f_10 * sgi0_102[k]
                   - f_11 * sgi1_102[k]
                   + f_3 * pc_x[k] * sgk_126[k];

        t_154[k] = f_16 * sfk_50[k]
                   + f_3 * pc_y[k] * sgk_122[k];

        t_155[k] = f_16 * sfk_128[k]
                   + f_10 * sgi0_104[k]
                   - f_11 * sgi1_104[k]
                   + f_3 * pc_x[k] * sgk_128[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pc_x, pc_z, sfk_129, sfk_131, sgi0_105, \
                         sgi0_107, sgi1_105, sgi1_107, sgk_123, sgk_129, \
                         sgk_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_16 * sfk_129[k]
                   + f_12 * sgi0_105[k]
                   - f_13 * sgi1_105[k]
                   + f_3 * pc_x[k] * sgk_129[k];

        t_157[k] = f_3 * pc_z[k] * sgk_123[k];

        t_158[k] = f_16 * sfk_131[k]
                   + f_12 * sgi0_107[k]
                   - f_13 * sgi1_107[k]
                   + f_3 * pc_x[k] * sgk_131[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pc_x, pc_y, sfk_56, sfk_132, sfk_133, sgi0_108, \
                         sgi0_109, sgi1_108, sgi1_109, sgk_128, sgk_132, \
                         sgk_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_16 * sfk_132[k]
                   + f_12 * sgi0_108[k]
                   - f_13 * sgi1_108[k]
                   + f_3 * pc_x[k] * sgk_132[k];

        t_160[k] = f_16 * sfk_133[k]
                   + f_12 * sgi0_109[k]
                   - f_13 * sgi1_109[k]
                   + f_3 * pc_x[k] * sgk_133[k];

        t_161[k] = f_16 * sfk_56[k]
                   + f_3 * pc_y[k] * sgk_128[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, sfk_135, sfk_136, sfk_137, sfk_138, \
                         sgi0_111, sgi1_111, sgk_135, sgk_136, sgk_137, \
                         sgk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_16 * sfk_135[k]
                   + f_12 * sgi0_111[k]
                   - f_13 * sgi1_111[k]
                   + f_3 * pc_x[k] * sgk_135[k];

        t_163[k] = f_16 * sfk_136[k]
                   + f_3 * pc_x[k] * sgk_136[k];

        t_164[k] = f_16 * sfk_137[k]
                   + f_3 * pc_x[k] * sgk_137[k];

        t_165[k] = f_16 * sfk_138[k]
                   + f_3 * pc_x[k] * sgk_138[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, sfk_139, sfk_140, sfk_141, \
                         sfk_142, sfk_143, sgk_139, sgk_140, sgk_141, sgk_142, \
                         sgk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_16 * sfk_139[k]
                   + f_3 * pc_x[k] * sgk_139[k];

        t_167[k] = f_16 * sfk_140[k]
                   + f_3 * pc_x[k] * sgk_140[k];

        t_168[k] = f_16 * sfk_141[k]
                   + f_3 * pc_x[k] * sgk_141[k];

        t_169[k] = f_16 * sfk_142[k]
                   + f_3 * pc_x[k] * sgk_142[k];

        t_170[k] = f_16 * sfk_143[k]
                   + f_3 * pc_x[k] * sgk_143[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pc_y, pc_z, sfk_64, sfk_66, sgi0_105, sgi0_107, \
                         sgi1_105, sgi1_107, sgk_136, sgk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_16 * sfk_64[k]
                   + f_1 * sgi0_105[k]
                   - f_2 * sgi1_105[k]
                   + f_3 * pc_y[k] * sgk_136[k];

        t_172[k] = f_3 * pc_z[k] * sgk_136[k];

        t_173[k] = f_16 * sfk_66[k]
                   + f_4 * sgi0_107[k]
                   - f_5 * sgi1_107[k]
                   + f_3 * pc_y[k] * sgk_138[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, sfk_67, sfk_68, sfk_69, sgi0_108, \
                         sgi0_109, sgi0_110, sgi1_108, sgi1_109, sgi1_110, sgk_139, sgk_140, \
                         sgk_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_16 * sfk_67[k]
                   + f_6 * sgi0_108[k]
                   - f_7 * sgi1_108[k]
                   + f_3 * pc_y[k] * sgk_139[k];

        t_175[k] = f_16 * sfk_68[k]
                   + f_8 * sgi0_109[k]
                   - f_9 * sgi1_109[k]
                   + f_3 * pc_y[k] * sgk_140[k];

        t_176[k] = f_16 * sfk_69[k]
                   + f_10 * sgi0_110[k]
                   - f_11 * sgi1_110[k]
                   + f_3 * pc_y[k] * sgk_141[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pb_y, pc_y, pc_z, sfl0_90, sfk_70, \
                         sfk_71, sfl1_90, sgi0_111, sgi1_111, sgk_142, \
                         sgk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_16 * sfk_70[k]
                   + f_12 * sgi0_111[k]
                   - f_13 * sgi1_111[k]
                   + f_3 * pc_y[k] * sgk_142[k];

        t_178[k] = f_16 * sfk_71[k]
                   + f_3 * pc_y[k] * sgk_143[k];

        t_179[k] = f_1 * sgi0_111[k]
                   - f_2 * sgi1_111[k]
                   + f_3 * pc_z[k] * sgk_143[k];

        t_180[k] = pb_y[k] * sfl0_90[k]
                   - f_14 * pc_y[k] * sfl1_90[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_z, pc_y, pc_z, sfl0_48, sfk_36, \
                         sfk_72, sfk_74, sfl1_48, sgk_144, sgk_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_15 * sfk_72[k]
                   + f_3 * pc_y[k] * sgk_144[k];

        t_182[k] = f_15 * sfk_36[k]
                   + f_3 * pc_z[k] * sgk_144[k];

        t_183[k] = pb_z[k] * sfl0_48[k]
                   - f_14 * pc_z[k] * sfl1_48[k];

        t_184[k] = f_15 * sfk_74[k]
                   + f_3 * pc_y[k] * sgk_146[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_y, pb_z, pc_y, pc_z, sfl0_51, sfl0_95, \
                         sfk_39, sfk_77, sfl1_51, sfl1_95, sgk_147, \
                         sgk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pb_y[k] * sfl0_95[k]
                   - f_14 * pc_y[k] * sfl1_95[k];

        t_186[k] = pb_z[k] * sfl0_51[k]
                   - f_14 * pc_z[k] * sfl1_51[k];

        t_187[k] = f_15 * sfk_39[k]
                   + f_3 * pc_z[k] * sgk_147[k];

        t_188[k] = f_15 * sfk_77[k]
                   + f_3 * pc_y[k] * sgk_149[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pb_y, pb_z, pc_y, pc_z, sfl0_55, sfl0_99, \
                         sfk_42, sfl1_55, sfl1_99, sgk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pb_y[k] * sfl0_99[k]
                   - f_14 * pc_y[k] * sfl1_99[k];

        t_190[k] = pb_z[k] * sfl0_55[k]
                   - f_14 * pc_z[k] * sfl1_55[k];

        t_191[k] = f_15 * sfk_42[k]
                   + f_3 * pc_z[k] * sgk_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_y, pc_y, sfl0_102, sfl0_104, sfk_80, sfk_81, \
                         sfl1_102, sfl1_104, sgk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pb_y[k] * sfl0_102[k]
                   + f_16 * sfk_80[k]
                   - f_14 * pc_y[k] * sfl1_102[k];

        t_193[k] = f_15 * sfk_81[k]
                   + f_3 * pc_y[k] * sgk_153[k];

        t_194[k] = pb_y[k] * sfl0_104[k]
                   - f_14 * pc_y[k] * sfl1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pb_y, pb_z, pc_y, pc_z, sfl0_60, sfl0_107, \
                         sfk_46, sfk_84, sfl1_60, sfl1_107, sgk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pb_z[k] * sfl0_60[k]
                   - f_14 * pc_z[k] * sfl1_60[k];

        t_196[k] = f_15 * sfk_46[k]
                   + f_3 * pc_z[k] * sgk_154[k];

        t_197[k] = pb_y[k] * sfl0_107[k]
                   + f_17 * sfk_84[k]
                   - f_14 * pc_y[k] * sfl1_107[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pb_y, pc_y, sfl0_108, sfl0_110, sfk_85, sfk_86, \
                         sfl1_108, sfl1_110, sgk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pb_y[k] * sfl0_108[k]
                   + f_16 * sfk_85[k]
                   - f_14 * pc_y[k] * sfl1_108[k];

        t_199[k] = f_15 * sfk_86[k]
                   + f_3 * pc_y[k] * sgk_158[k];

        t_200[k] = pb_y[k] * sfl0_110[k]
                   - f_14 * pc_y[k] * sfl1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pb_y, pb_z, pc_y, pc_z, sfl0_66, sfl0_113, \
                         sfk_51, sfk_89, sfl1_66, sfl1_113, sgk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pb_z[k] * sfl0_66[k]
                   - f_14 * pc_z[k] * sfl1_66[k];

        t_202[k] = f_15 * sfk_51[k]
                   + f_3 * pc_z[k] * sgk_159[k];

        t_203[k] = pb_y[k] * sfl0_113[k]
                   + f_0 * sfk_89[k]
                   - f_14 * pc_y[k] * sfl1_113[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_y, pc_y, sfl0_114, sfl0_115, sfl0_117, \
                         sfk_90, sfk_91, sfk_92, sfl1_114, sfl1_115, sfl1_117, \
                         sgk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pb_y[k] * sfl0_114[k]
                   + f_17 * sfk_90[k]
                   - f_14 * pc_y[k] * sfl1_114[k];

        t_205[k] = pb_y[k] * sfl0_115[k]
                   + f_16 * sfk_91[k]
                   - f_14 * pc_y[k] * sfl1_115[k];

        t_206[k] = f_15 * sfk_92[k]
                   + f_3 * pc_y[k] * sgk_164[k];

        t_207[k] = pb_y[k] * sfl0_117[k]
                   - f_14 * pc_y[k] * sfl1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pc_x, sfk_172, sfk_173, sfk_174, \
                         sfk_175, sfk_176, sgk_172, sgk_173, sgk_174, sgk_175, \
                         sgk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_16 * sfk_172[k]
                   + f_3 * pc_x[k] * sgk_172[k];

        t_209[k] = f_16 * sfk_173[k]
                   + f_3 * pc_x[k] * sgk_173[k];

        t_210[k] = f_16 * sfk_174[k]
                   + f_3 * pc_x[k] * sgk_174[k];

        t_211[k] = f_16 * sfk_175[k]
                   + f_3 * pc_x[k] * sgk_175[k];

        t_212[k] = f_16 * sfk_176[k]
                   + f_3 * pc_x[k] * sgk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pb_z, pc_x, pc_z, sfl0_81, sfk_177, \
                         sfk_178, sfk_179, sfl1_81, sgk_177, sgk_178, \
                         sgk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_16 * sfk_177[k]
                   + f_3 * pc_x[k] * sgk_177[k];

        t_214[k] = f_16 * sfk_178[k]
                   + f_3 * pc_x[k] * sgk_178[k];

        t_215[k] = f_16 * sfk_179[k]
                   + f_3 * pc_x[k] * sgk_179[k];

        t_216[k] = pb_z[k] * sfl0_81[k]
                   - f_14 * pc_z[k] * sfl1_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pc_y, pc_z, sfk_64, sfk_102, sfk_103, sgi0_135, \
                         sgi0_136, sgi1_135, sgi1_136, sgk_172, sgk_174, \
                         sgk_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * sfk_64[k]
                   + f_3 * pc_z[k] * sgk_172[k];

        t_218[k] = f_15 * sfk_102[k]
                   + f_4 * sgi0_135[k]
                   - f_5 * sgi1_135[k]
                   + f_3 * pc_y[k] * sgk_174[k];

        t_219[k] = f_15 * sfk_103[k]
                   + f_6 * sgi0_136[k]
                   - f_7 * sgi1_136[k]
                   + f_3 * pc_y[k] * sgk_175[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_y, sfk_104, sfk_105, sfk_106, sgi0_137, \
                         sgi0_138, sgi0_139, sgi1_137, sgi1_138, sgi1_139, sgk_176, sgk_177, \
                         sgk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * sfk_104[k]
                   + f_8 * sgi0_137[k]
                   - f_9 * sgi1_137[k]
                   + f_3 * pc_y[k] * sgk_176[k];

        t_221[k] = f_15 * sfk_105[k]
                   + f_10 * sgi0_138[k]
                   - f_11 * sgi1_138[k]
                   + f_3 * pc_y[k] * sgk_177[k];

        t_222[k] = f_15 * sfk_106[k]
                   + f_12 * sgi0_139[k]
                   - f_13 * sgi1_139[k]
                   + f_3 * pc_y[k] * sgk_178[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pb_y, pc_x, pc_y, sfl0_134, sfk_107, \
                         sfk_180, sfl1_134, sgi0_140, sgi1_140, sgk_179, \
                         sgk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * sfk_107[k]
                   + f_3 * pc_y[k] * sgk_179[k];

        t_224[k] = pb_y[k] * sfl0_134[k]
                   - f_14 * pc_y[k] * sfl1_134[k];

        t_225[k] = f_16 * sfk_180[k]
                   + f_1 * sgi0_140[k]
                   - f_2 * sgi1_140[k]
                   + f_3 * pc_x[k] * sgk_180[k];

        t_226[k] = f_3 * pc_y[k] * sgk_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_x, pc_y, pc_z, sfk_72, sfk_183, sgi0_143, \
                         sgi1_143, sgk_180, sgk_182, sgk_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_16 * sfk_72[k]
                   + f_3 * pc_z[k] * sgk_180[k];

        t_228[k] = f_16 * sfk_183[k]
                   + f_4 * sgi0_143[k]
                   - f_5 * sgi1_143[k]
                   + f_3 * pc_x[k] * sgk_183[k];

        t_229[k] = f_3 * pc_y[k] * sgk_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, pc_x, pc_z, sfk_75, sfk_185, sfk_186, sgi0_145, \
                         sgi0_146, sgi1_145, sgi1_146, sgk_183, sgk_185, \
                         sgk_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_16 * sfk_185[k]
                   + f_4 * sgi0_145[k]
                   - f_5 * sgi1_145[k]
                   + f_3 * pc_x[k] * sgk_185[k];

        t_231[k] = f_16 * sfk_186[k]
                   + f_6 * sgi0_146[k]
                   - f_7 * sgi1_146[k]
                   + f_3 * pc_x[k] * sgk_186[k];

        t_232[k] = f_16 * sfk_75[k]
                   + f_3 * pc_z[k] * sgk_183[k];
    }
}

static auto
compute_prim_sgl_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfl0,
                                                          const size_t sfk, const size_t sfl1,
                                                          const size_t sgi0, const size_t sgi1,
                                                          const size_t sgk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 4.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfl0_135 = buffer.data(sfl0 + 135);
    const auto *sfl0_138 = buffer.data(sfl0 + 138);
    const auto *sfl0_141 = buffer.data(sfl0 + 141);
    const auto *sfl0_145 = buffer.data(sfl0 + 145);
    const auto *sfl0_150 = buffer.data(sfl0 + 150);
    const auto *sfl0_156 = buffer.data(sfl0 + 156);
    const auto *sfl0_270 = buffer.data(sfl0 + 270);
    const auto *sfl0_273 = buffer.data(sfl0 + 273);
    const auto *sfl0_275 = buffer.data(sfl0 + 275);
    const auto *sfl0_276 = buffer.data(sfl0 + 276);
    const auto *sfl0_279 = buffer.data(sfl0 + 279);
    const auto *sfl0_280 = buffer.data(sfl0 + 280);
    const auto *sfl0_282 = buffer.data(sfl0 + 282);
    const auto *sfl0_284 = buffer.data(sfl0 + 284);
    const auto *sfl0_285 = buffer.data(sfl0 + 285);
    const auto *sfl0_287 = buffer.data(sfl0 + 287);
    const auto *sfl0_288 = buffer.data(sfl0 + 288);
    const auto *sfl0_290 = buffer.data(sfl0 + 290);
    const auto *sfl0_291 = buffer.data(sfl0 + 291);
    const auto *sfl0_293 = buffer.data(sfl0 + 293);
    const auto *sfl0_294 = buffer.data(sfl0 + 294);
    const auto *sfl0_295 = buffer.data(sfl0 + 295);
    const auto *sfl0_297 = buffer.data(sfl0 + 297);
    const auto *sfl0_306 = buffer.data(sfl0 + 306);
    const auto *sfl0_308 = buffer.data(sfl0 + 308);
    const auto *sfl0_309 = buffer.data(sfl0 + 309);
    const auto *sfl0_310 = buffer.data(sfl0 + 310);
    const auto *sfl0_311 = buffer.data(sfl0 + 311);
    const auto *sfl0_312 = buffer.data(sfl0 + 312);
    const auto *sfl0_314 = buffer.data(sfl0 + 314);
    const auto *sfl0_320 = buffer.data(sfl0 + 320);
    const auto *sfl0_324 = buffer.data(sfl0 + 324);
    const auto *sfl0_327 = buffer.data(sfl0 + 327);
    const auto *sfl0_329 = buffer.data(sfl0 + 329);
    const auto *sfl0_332 = buffer.data(sfl0 + 332);
    const auto *sfl0_333 = buffer.data(sfl0 + 333);
    const auto *sfl0_335 = buffer.data(sfl0 + 335);
    const auto *sfl0_338 = buffer.data(sfl0 + 338);
    const auto *sfl0_339 = buffer.data(sfl0 + 339);
    const auto *sfl0_340 = buffer.data(sfl0 + 340);
    const auto *sfl0_342 = buffer.data(sfl0 + 342);

    const auto *sfk_78 = buffer.data(sfk + 78);
    const auto *sfk_82 = buffer.data(sfk + 82);
    const auto *sfk_87 = buffer.data(sfk + 87);
    const auto *sfk_100 = buffer.data(sfk + 100);
    const auto *sfk_107 = buffer.data(sfk + 107);
    const auto *sfk_108 = buffer.data(sfk + 108);
    const auto *sfk_110 = buffer.data(sfk + 110);
    const auto *sfk_111 = buffer.data(sfk + 111);
    const auto *sfk_113 = buffer.data(sfk + 113);
    const auto *sfk_114 = buffer.data(sfk + 114);
    const auto *sfk_117 = buffer.data(sfk + 117);
    const auto *sfk_118 = buffer.data(sfk + 118);
    const auto *sfk_122 = buffer.data(sfk + 122);
    const auto *sfk_123 = buffer.data(sfk + 123);
    const auto *sfk_128 = buffer.data(sfk + 128);
    const auto *sfk_143 = buffer.data(sfk + 143);
    const auto *sfk_144 = buffer.data(sfk + 144);
    const auto *sfk_146 = buffer.data(sfk + 146);
    const auto *sfk_149 = buffer.data(sfk + 149);
    const auto *sfk_153 = buffer.data(sfk + 153);
    const auto *sfk_158 = buffer.data(sfk + 158);
    const auto *sfk_164 = buffer.data(sfk + 164);
    const auto *sfk_189 = buffer.data(sfk + 189);
    const auto *sfk_190 = buffer.data(sfk + 190);
    const auto *sfk_192 = buffer.data(sfk + 192);
    const auto *sfk_194 = buffer.data(sfk + 194);
    const auto *sfk_195 = buffer.data(sfk + 195);
    const auto *sfk_197 = buffer.data(sfk + 197);
    const auto *sfk_198 = buffer.data(sfk + 198);
    const auto *sfk_200 = buffer.data(sfk + 200);
    const auto *sfk_201 = buffer.data(sfk + 201);
    const auto *sfk_203 = buffer.data(sfk + 203);
    const auto *sfk_204 = buffer.data(sfk + 204);
    const auto *sfk_205 = buffer.data(sfk + 205);
    const auto *sfk_207 = buffer.data(sfk + 207);
    const auto *sfk_208 = buffer.data(sfk + 208);
    const auto *sfk_209 = buffer.data(sfk + 209);
    const auto *sfk_210 = buffer.data(sfk + 210);
    const auto *sfk_211 = buffer.data(sfk + 211);
    const auto *sfk_212 = buffer.data(sfk + 212);
    const auto *sfk_213 = buffer.data(sfk + 213);
    const auto *sfk_214 = buffer.data(sfk + 214);
    const auto *sfk_215 = buffer.data(sfk + 215);
    const auto *sfk_216 = buffer.data(sfk + 216);
    const auto *sfk_219 = buffer.data(sfk + 219);
    const auto *sfk_221 = buffer.data(sfk + 221);
    const auto *sfk_222 = buffer.data(sfk + 222);
    const auto *sfk_225 = buffer.data(sfk + 225);
    const auto *sfk_226 = buffer.data(sfk + 226);
    const auto *sfk_228 = buffer.data(sfk + 228);
    const auto *sfk_230 = buffer.data(sfk + 230);
    const auto *sfk_231 = buffer.data(sfk + 231);
    const auto *sfk_233 = buffer.data(sfk + 233);
    const auto *sfk_234 = buffer.data(sfk + 234);
    const auto *sfk_236 = buffer.data(sfk + 236);
    const auto *sfk_237 = buffer.data(sfk + 237);
    const auto *sfk_239 = buffer.data(sfk + 239);
    const auto *sfk_240 = buffer.data(sfk + 240);
    const auto *sfk_241 = buffer.data(sfk + 241);
    const auto *sfk_243 = buffer.data(sfk + 243);
    const auto *sfk_244 = buffer.data(sfk + 244);
    const auto *sfk_245 = buffer.data(sfk + 245);
    const auto *sfk_246 = buffer.data(sfk + 246);
    const auto *sfk_247 = buffer.data(sfk + 247);
    const auto *sfk_248 = buffer.data(sfk + 248);
    const auto *sfk_249 = buffer.data(sfk + 249);
    const auto *sfk_250 = buffer.data(sfk + 250);
    const auto *sfk_251 = buffer.data(sfk + 251);
    const auto *sfk_257 = buffer.data(sfk + 257);
    const auto *sfk_261 = buffer.data(sfk + 261);
    const auto *sfk_264 = buffer.data(sfk + 264);
    const auto *sfk_266 = buffer.data(sfk + 266);
    const auto *sfk_269 = buffer.data(sfk + 269);
    const auto *sfk_270 = buffer.data(sfk + 270);
    const auto *sfk_272 = buffer.data(sfk + 272);
    const auto *sfk_275 = buffer.data(sfk + 275);
    const auto *sfk_276 = buffer.data(sfk + 276);
    const auto *sfk_277 = buffer.data(sfk + 277);
    const auto *sfk_279 = buffer.data(sfk + 279);
    const auto *sfk_280 = buffer.data(sfk + 280);
    const auto *sfk_281 = buffer.data(sfk + 281);
    const auto *sfk_282 = buffer.data(sfk + 282);
    const auto *sfk_283 = buffer.data(sfk + 283);
    const auto *sfk_284 = buffer.data(sfk + 284);
    const auto *sfk_285 = buffer.data(sfk + 285);
    const auto *sfk_286 = buffer.data(sfk + 286);

    const auto *sfl1_135 = buffer.data(sfl1 + 135);
    const auto *sfl1_138 = buffer.data(sfl1 + 138);
    const auto *sfl1_141 = buffer.data(sfl1 + 141);
    const auto *sfl1_145 = buffer.data(sfl1 + 145);
    const auto *sfl1_150 = buffer.data(sfl1 + 150);
    const auto *sfl1_156 = buffer.data(sfl1 + 156);
    const auto *sfl1_270 = buffer.data(sfl1 + 270);
    const auto *sfl1_273 = buffer.data(sfl1 + 273);
    const auto *sfl1_275 = buffer.data(sfl1 + 275);
    const auto *sfl1_276 = buffer.data(sfl1 + 276);
    const auto *sfl1_279 = buffer.data(sfl1 + 279);
    const auto *sfl1_280 = buffer.data(sfl1 + 280);
    const auto *sfl1_282 = buffer.data(sfl1 + 282);
    const auto *sfl1_284 = buffer.data(sfl1 + 284);
    const auto *sfl1_285 = buffer.data(sfl1 + 285);
    const auto *sfl1_287 = buffer.data(sfl1 + 287);
    const auto *sfl1_288 = buffer.data(sfl1 + 288);
    const auto *sfl1_290 = buffer.data(sfl1 + 290);
    const auto *sfl1_291 = buffer.data(sfl1 + 291);
    const auto *sfl1_293 = buffer.data(sfl1 + 293);
    const auto *sfl1_294 = buffer.data(sfl1 + 294);
    const auto *sfl1_295 = buffer.data(sfl1 + 295);
    const auto *sfl1_297 = buffer.data(sfl1 + 297);
    const auto *sfl1_306 = buffer.data(sfl1 + 306);
    const auto *sfl1_308 = buffer.data(sfl1 + 308);
    const auto *sfl1_309 = buffer.data(sfl1 + 309);
    const auto *sfl1_310 = buffer.data(sfl1 + 310);
    const auto *sfl1_311 = buffer.data(sfl1 + 311);
    const auto *sfl1_312 = buffer.data(sfl1 + 312);
    const auto *sfl1_314 = buffer.data(sfl1 + 314);
    const auto *sfl1_320 = buffer.data(sfl1 + 320);
    const auto *sfl1_324 = buffer.data(sfl1 + 324);
    const auto *sfl1_327 = buffer.data(sfl1 + 327);
    const auto *sfl1_329 = buffer.data(sfl1 + 329);
    const auto *sfl1_332 = buffer.data(sfl1 + 332);
    const auto *sfl1_333 = buffer.data(sfl1 + 333);
    const auto *sfl1_335 = buffer.data(sfl1 + 335);
    const auto *sfl1_338 = buffer.data(sfl1 + 338);
    const auto *sfl1_339 = buffer.data(sfl1 + 339);
    const auto *sfl1_340 = buffer.data(sfl1 + 340);
    const auto *sfl1_342 = buffer.data(sfl1 + 342);

    const auto *sgi0_149 = buffer.data(sgi0 + 149);
    const auto *sgi0_150 = buffer.data(sgi0 + 150);
    const auto *sgi0_152 = buffer.data(sgi0 + 152);
    const auto *sgi0_154 = buffer.data(sgi0 + 154);
    const auto *sgi0_155 = buffer.data(sgi0 + 155);
    const auto *sgi0_157 = buffer.data(sgi0 + 157);
    const auto *sgi0_158 = buffer.data(sgi0 + 158);
    const auto *sgi0_160 = buffer.data(sgi0 + 160);
    const auto *sgi0_161 = buffer.data(sgi0 + 161);
    const auto *sgi0_163 = buffer.data(sgi0 + 163);
    const auto *sgi0_164 = buffer.data(sgi0 + 164);
    const auto *sgi0_165 = buffer.data(sgi0 + 165);
    const auto *sgi0_166 = buffer.data(sgi0 + 166);
    const auto *sgi0_167 = buffer.data(sgi0 + 167);

    const auto *sgi1_149 = buffer.data(sgi1 + 149);
    const auto *sgi1_150 = buffer.data(sgi1 + 150);
    const auto *sgi1_152 = buffer.data(sgi1 + 152);
    const auto *sgi1_154 = buffer.data(sgi1 + 154);
    const auto *sgi1_155 = buffer.data(sgi1 + 155);
    const auto *sgi1_157 = buffer.data(sgi1 + 157);
    const auto *sgi1_158 = buffer.data(sgi1 + 158);
    const auto *sgi1_160 = buffer.data(sgi1 + 160);
    const auto *sgi1_161 = buffer.data(sgi1 + 161);
    const auto *sgi1_163 = buffer.data(sgi1 + 163);
    const auto *sgi1_164 = buffer.data(sgi1 + 164);
    const auto *sgi1_165 = buffer.data(sgi1 + 165);
    const auto *sgi1_166 = buffer.data(sgi1 + 166);
    const auto *sgi1_167 = buffer.data(sgi1 + 167);

    const auto *sgk_185 = buffer.data(sgk + 185);
    const auto *sgk_186 = buffer.data(sgk + 186);
    const auto *sgk_189 = buffer.data(sgk + 189);
    const auto *sgk_190 = buffer.data(sgk + 190);
    const auto *sgk_192 = buffer.data(sgk + 192);
    const auto *sgk_194 = buffer.data(sgk + 194);
    const auto *sgk_195 = buffer.data(sgk + 195);
    const auto *sgk_197 = buffer.data(sgk + 197);
    const auto *sgk_198 = buffer.data(sgk + 198);
    const auto *sgk_200 = buffer.data(sgk + 200);
    const auto *sgk_201 = buffer.data(sgk + 201);
    const auto *sgk_203 = buffer.data(sgk + 203);
    const auto *sgk_204 = buffer.data(sgk + 204);
    const auto *sgk_205 = buffer.data(sgk + 205);
    const auto *sgk_207 = buffer.data(sgk + 207);
    const auto *sgk_208 = buffer.data(sgk + 208);
    const auto *sgk_209 = buffer.data(sgk + 209);
    const auto *sgk_210 = buffer.data(sgk + 210);
    const auto *sgk_211 = buffer.data(sgk + 211);
    const auto *sgk_212 = buffer.data(sgk + 212);
    const auto *sgk_213 = buffer.data(sgk + 213);
    const auto *sgk_214 = buffer.data(sgk + 214);
    const auto *sgk_215 = buffer.data(sgk + 215);
    const auto *sgk_216 = buffer.data(sgk + 216);
    const auto *sgk_218 = buffer.data(sgk + 218);
    const auto *sgk_219 = buffer.data(sgk + 219);
    const auto *sgk_221 = buffer.data(sgk + 221);
    const auto *sgk_222 = buffer.data(sgk + 222);
    const auto *sgk_225 = buffer.data(sgk + 225);
    const auto *sgk_226 = buffer.data(sgk + 226);
    const auto *sgk_230 = buffer.data(sgk + 230);
    const auto *sgk_231 = buffer.data(sgk + 231);
    const auto *sgk_236 = buffer.data(sgk + 236);
    const auto *sgk_244 = buffer.data(sgk + 244);
    const auto *sgk_245 = buffer.data(sgk + 245);
    const auto *sgk_246 = buffer.data(sgk + 246);
    const auto *sgk_247 = buffer.data(sgk + 247);
    const auto *sgk_248 = buffer.data(sgk + 248);
    const auto *sgk_249 = buffer.data(sgk + 249);
    const auto *sgk_250 = buffer.data(sgk + 250);
    const auto *sgk_251 = buffer.data(sgk + 251);
    const auto *sgk_252 = buffer.data(sgk + 252);
    const auto *sgk_254 = buffer.data(sgk + 254);
    const auto *sgk_255 = buffer.data(sgk + 255);
    const auto *sgk_257 = buffer.data(sgk + 257);
    const auto *sgk_258 = buffer.data(sgk + 258);
    const auto *sgk_261 = buffer.data(sgk + 261);
    const auto *sgk_262 = buffer.data(sgk + 262);
    const auto *sgk_266 = buffer.data(sgk + 266);
    const auto *sgk_267 = buffer.data(sgk + 267);
    const auto *sgk_272 = buffer.data(sgk + 272);
    const auto *sgk_280 = buffer.data(sgk + 280);
    const auto *sgk_281 = buffer.data(sgk + 281);
    const auto *sgk_282 = buffer.data(sgk + 282);
    const auto *sgk_283 = buffer.data(sgk + 283);
    const auto *sgk_284 = buffer.data(sgk + 284);
    const auto *sgk_285 = buffer.data(sgk + 285);
    const auto *sgk_286 = buffer.data(sgk + 286);

#pragma omp simd aligned(t_233, t_234, t_235, pc_x, pc_y, sfk_189, sfk_190, sgi0_149, \
                         sgi0_150, sgi1_149, sgi1_150, sgk_185, sgk_189, \
                         sgk_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * pc_y[k] * sgk_185[k];

        t_234[k] = f_16 * sfk_189[k]
                   + f_6 * sgi0_149[k]
                   - f_7 * sgi1_149[k]
                   + f_3 * pc_x[k] * sgk_189[k];

        t_235[k] = f_16 * sfk_190[k]
                   + f_8 * sgi0_150[k]
                   - f_9 * sgi1_150[k]
                   + f_3 * pc_x[k] * sgk_190[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pc_x, pc_y, pc_z, sfk_78, sfk_192, sgi0_152, \
                         sgi1_152, sgk_186, sgk_189, sgk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_16 * sfk_78[k]
                   + f_3 * pc_z[k] * sgk_186[k];

        t_237[k] = f_16 * sfk_192[k]
                   + f_8 * sgi0_152[k]
                   - f_9 * sgi1_152[k]
                   + f_3 * pc_x[k] * sgk_192[k];

        t_238[k] = f_3 * pc_y[k] * sgk_189[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, pc_z, sfk_82, sfk_194, sfk_195, sgi0_154, \
                         sgi0_155, sgi1_154, sgi1_155, sgk_190, sgk_194, \
                         sgk_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_16 * sfk_194[k]
                   + f_8 * sgi0_154[k]
                   - f_9 * sgi1_154[k]
                   + f_3 * pc_x[k] * sgk_194[k];

        t_240[k] = f_16 * sfk_195[k]
                   + f_10 * sgi0_155[k]
                   - f_11 * sgi1_155[k]
                   + f_3 * pc_x[k] * sgk_195[k];

        t_241[k] = f_16 * sfk_82[k]
                   + f_3 * pc_z[k] * sgk_190[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, pc_x, pc_y, sfk_197, sfk_198, sgi0_157, \
                         sgi0_158, sgi1_157, sgi1_158, sgk_194, sgk_197, \
                         sgk_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_16 * sfk_197[k]
                   + f_10 * sgi0_157[k]
                   - f_11 * sgi1_157[k]
                   + f_3 * pc_x[k] * sgk_197[k];

        t_243[k] = f_16 * sfk_198[k]
                   + f_10 * sgi0_158[k]
                   - f_11 * sgi1_158[k]
                   + f_3 * pc_x[k] * sgk_198[k];

        t_244[k] = f_3 * pc_y[k] * sgk_194[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_x, pc_z, sfk_87, sfk_200, sfk_201, sgi0_160, \
                         sgi0_161, sgi1_160, sgi1_161, sgk_195, sgk_200, \
                         sgk_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_16 * sfk_200[k]
                   + f_10 * sgi0_160[k]
                   - f_11 * sgi1_160[k]
                   + f_3 * pc_x[k] * sgk_200[k];

        t_246[k] = f_16 * sfk_201[k]
                   + f_12 * sgi0_161[k]
                   - f_13 * sgi1_161[k]
                   + f_3 * pc_x[k] * sgk_201[k];

        t_247[k] = f_16 * sfk_87[k]
                   + f_3 * pc_z[k] * sgk_195[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_x, sfk_203, sfk_204, sfk_205, sgi0_163, \
                         sgi0_164, sgi0_165, sgi1_163, sgi1_164, sgi1_165, sgk_203, sgk_204, \
                         sgk_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_16 * sfk_203[k]
                   + f_12 * sgi0_163[k]
                   - f_13 * sgi1_163[k]
                   + f_3 * pc_x[k] * sgk_203[k];

        t_249[k] = f_16 * sfk_204[k]
                   + f_12 * sgi0_164[k]
                   - f_13 * sgi1_164[k]
                   + f_3 * pc_x[k] * sgk_204[k];

        t_250[k] = f_16 * sfk_205[k]
                   + f_12 * sgi0_165[k]
                   - f_13 * sgi1_165[k]
                   + f_3 * pc_x[k] * sgk_205[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pc_x, pc_y, sfk_207, sfk_208, sfk_209, \
                         sgi0_167, sgi1_167, sgk_200, sgk_207, sgk_208, \
                         sgk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_3 * pc_y[k] * sgk_200[k];

        t_252[k] = f_16 * sfk_207[k]
                   + f_12 * sgi0_167[k]
                   - f_13 * sgi1_167[k]
                   + f_3 * pc_x[k] * sgk_207[k];

        t_253[k] = f_16 * sfk_208[k]
                   + f_3 * pc_x[k] * sgk_208[k];

        t_254[k] = f_16 * sfk_209[k]
                   + f_3 * pc_x[k] * sgk_209[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pc_x, sfk_210, sfk_211, sfk_212, \
                         sfk_213, sfk_214, sgk_210, sgk_211, sgk_212, sgk_213, \
                         sgk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_16 * sfk_210[k]
                   + f_3 * pc_x[k] * sgk_210[k];

        t_256[k] = f_16 * sfk_211[k]
                   + f_3 * pc_x[k] * sgk_211[k];

        t_257[k] = f_16 * sfk_212[k]
                   + f_3 * pc_x[k] * sgk_212[k];

        t_258[k] = f_16 * sfk_213[k]
                   + f_3 * pc_x[k] * sgk_213[k];

        t_259[k] = f_16 * sfk_214[k]
                   + f_3 * pc_x[k] * sgk_214[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, pc_z, sfk_100, sfk_215, \
                         sgi0_161, sgi0_163, sgi1_161, sgi1_163, sgk_208, sgk_210, \
                         sgk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_16 * sfk_215[k]
                   + f_3 * pc_x[k] * sgk_215[k];

        t_261[k] = f_1 * sgi0_161[k]
                   - f_2 * sgi1_161[k]
                   + f_3 * pc_y[k] * sgk_208[k];

        t_262[k] = f_16 * sfk_100[k]
                   + f_3 * pc_z[k] * sgk_208[k];

        t_263[k] = f_4 * sgi0_163[k]
                   - f_5 * sgi1_163[k]
                   + f_3 * pc_y[k] * sgk_210[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, sgi0_164, sgi0_165, sgi0_166, sgi1_164, \
                         sgi1_165, sgi1_166, sgk_211, sgk_212, \
                         sgk_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_6 * sgi0_164[k]
                   - f_7 * sgi1_164[k]
                   + f_3 * pc_y[k] * sgk_211[k];

        t_265[k] = f_8 * sgi0_165[k]
                   - f_9 * sgi1_165[k]
                   + f_3 * pc_y[k] * sgk_212[k];

        t_266[k] = f_10 * sgi0_166[k]
                   - f_11 * sgi1_166[k]
                   + f_3 * pc_y[k] * sgk_213[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pb_x, pc_x, pc_y, pc_z, sfl0_270, \
                         sfk_107, sfk_216, sfl1_270, sgi0_167, sgi1_167, sgk_214, \
                         sgk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_12 * sgi0_167[k]
                   - f_13 * sgi1_167[k]
                   + f_3 * pc_y[k] * sgk_214[k];

        t_268[k] = f_3 * pc_y[k] * sgk_215[k];

        t_269[k] = f_16 * sfk_107[k]
                   + f_1 * sgi0_167[k]
                   - f_2 * sgi1_167[k]
                   + f_3 * pc_z[k] * sgk_215[k];

        t_270[k] = pb_x[k] * sfl0_270[k]
                   + f_20 * sfk_216[k]
                   - f_14 * pc_x[k] * sfl1_270[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pb_x, pc_x, pc_y, pc_z, sfl0_273, \
                         sfk_108, sfk_110, sfk_219, sfl1_273, sgk_216, \
                         sgk_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_17 * sfk_108[k]
                   + f_3 * pc_y[k] * sgk_216[k];

        t_272[k] = f_3 * pc_z[k] * sgk_216[k];

        t_273[k] = pb_x[k] * sfl0_273[k]
                   + f_19 * sfk_219[k]
                   - f_14 * pc_x[k] * sfl1_273[k];

        t_274[k] = f_17 * sfk_110[k]
                   + f_3 * pc_y[k] * sgk_218[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pb_x, pc_x, pc_z, sfl0_275, sfl0_276, sfk_221, \
                         sfk_222, sfl1_275, sfl1_276, sgk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = pb_x[k] * sfl0_275[k]
                   + f_19 * sfk_221[k]
                   - f_14 * pc_x[k] * sfl1_275[k];

        t_276[k] = pb_x[k] * sfl0_276[k]
                   + f_18 * sfk_222[k]
                   - f_14 * pc_x[k] * sfl1_276[k];

        t_277[k] = f_3 * pc_z[k] * sgk_219[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pb_x, pc_x, pc_y, sfl0_279, sfl0_280, sfk_113, \
                         sfk_225, sfk_226, sfl1_279, sfl1_280, \
                         sgk_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_17 * sfk_113[k]
                   + f_3 * pc_y[k] * sgk_221[k];

        t_279[k] = pb_x[k] * sfl0_279[k]
                   + f_18 * sfk_225[k]
                   - f_14 * pc_x[k] * sfl1_279[k];

        t_280[k] = pb_x[k] * sfl0_280[k]
                   + f_0 * sfk_226[k]
                   - f_14 * pc_x[k] * sfl1_280[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, pb_x, pc_x, pc_y, pc_z, sfl0_282, sfk_117, \
                         sfk_228, sfl1_282, sgk_222, sgk_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_3 * pc_z[k] * sgk_222[k];

        t_282[k] = pb_x[k] * sfl0_282[k]
                   + f_0 * sfk_228[k]
                   - f_14 * pc_x[k] * sfl1_282[k];

        t_283[k] = f_17 * sfk_117[k]
                   + f_3 * pc_y[k] * sgk_225[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pb_x, pc_x, pc_z, sfl0_284, sfl0_285, sfk_230, \
                         sfk_231, sfl1_284, sfl1_285, sgk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = pb_x[k] * sfl0_284[k]
                   + f_0 * sfk_230[k]
                   - f_14 * pc_x[k] * sfl1_284[k];

        t_285[k] = pb_x[k] * sfl0_285[k]
                   + f_17 * sfk_231[k]
                   - f_14 * pc_x[k] * sfl1_285[k];

        t_286[k] = f_3 * pc_z[k] * sgk_226[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pb_x, pc_x, pc_y, sfl0_287, sfl0_288, sfk_122, \
                         sfk_233, sfk_234, sfl1_287, sfl1_288, \
                         sgk_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = pb_x[k] * sfl0_287[k]
                   + f_17 * sfk_233[k]
                   - f_14 * pc_x[k] * sfl1_287[k];

        t_288[k] = pb_x[k] * sfl0_288[k]
                   + f_17 * sfk_234[k]
                   - f_14 * pc_x[k] * sfl1_288[k];

        t_289[k] = f_17 * sfk_122[k]
                   + f_3 * pc_y[k] * sgk_230[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pb_x, pc_x, pc_z, sfl0_290, sfl0_291, sfk_236, \
                         sfk_237, sfl1_290, sfl1_291, sgk_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_x[k] * sfl0_290[k]
                   + f_17 * sfk_236[k]
                   - f_14 * pc_x[k] * sfl1_290[k];

        t_291[k] = pb_x[k] * sfl0_291[k]
                   + f_16 * sfk_237[k]
                   - f_14 * pc_x[k] * sfl1_291[k];

        t_292[k] = f_3 * pc_z[k] * sgk_231[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, pb_x, pc_x, sfl0_293, sfl0_294, sfl0_295, \
                         sfk_239, sfk_240, sfk_241, sfl1_293, sfl1_294, \
                         sfl1_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = pb_x[k] * sfl0_293[k]
                   + f_16 * sfk_239[k]
                   - f_14 * pc_x[k] * sfl1_293[k];

        t_294[k] = pb_x[k] * sfl0_294[k]
                   + f_16 * sfk_240[k]
                   - f_14 * pc_x[k] * sfl1_294[k];

        t_295[k] = pb_x[k] * sfl0_295[k]
                   + f_16 * sfk_241[k]
                   - f_14 * pc_x[k] * sfl1_295[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pb_x, pc_x, pc_y, sfl0_297, sfk_128, \
                         sfk_243, sfk_244, sfk_245, sfl1_297, sgk_236, sgk_244, \
                         sgk_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_17 * sfk_128[k]
                   + f_3 * pc_y[k] * sgk_236[k];

        t_297[k] = pb_x[k] * sfl0_297[k]
                   + f_16 * sfk_243[k]
                   - f_14 * pc_x[k] * sfl1_297[k];

        t_298[k] = f_15 * sfk_244[k]
                   + f_3 * pc_x[k] * sgk_244[k];

        t_299[k] = f_15 * sfk_245[k]
                   + f_3 * pc_x[k] * sgk_245[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, pc_x, sfk_246, sfk_247, sfk_248, \
                         sfk_249, sfk_250, sgk_246, sgk_247, sgk_248, sgk_249, \
                         sgk_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_15 * sfk_246[k]
                   + f_3 * pc_x[k] * sgk_246[k];

        t_301[k] = f_15 * sfk_247[k]
                   + f_3 * pc_x[k] * sgk_247[k];

        t_302[k] = f_15 * sfk_248[k]
                   + f_3 * pc_x[k] * sgk_248[k];

        t_303[k] = f_15 * sfk_249[k]
                   + f_3 * pc_x[k] * sgk_249[k];

        t_304[k] = f_15 * sfk_250[k]
                   + f_3 * pc_x[k] * sgk_250[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pb_x, pc_x, pc_z, sfl0_306, sfl0_308, \
                         sfk_251, sfl1_306, sfl1_308, sgk_244, \
                         sgk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_15 * sfk_251[k]
                   + f_3 * pc_x[k] * sgk_251[k];

        t_306[k] = pb_x[k] * sfl0_306[k]
                   - f_14 * pc_x[k] * sfl1_306[k];

        t_307[k] = f_3 * pc_z[k] * sgk_244[k];

        t_308[k] = pb_x[k] * sfl0_308[k]
                   - f_14 * pc_x[k] * sfl1_308[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pb_x, pc_x, sfl0_309, sfl0_310, sfl0_311, \
                         sfl0_312, sfl1_309, sfl1_310, sfl1_311, \
                         sfl1_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = pb_x[k] * sfl0_309[k]
                   - f_14 * pc_x[k] * sfl1_309[k];

        t_310[k] = pb_x[k] * sfl0_310[k]
                   - f_14 * pc_x[k] * sfl1_310[k];

        t_311[k] = pb_x[k] * sfl0_311[k]
                   - f_14 * pc_x[k] * sfl1_311[k];

        t_312[k] = pb_x[k] * sfl0_312[k]
                   - f_14 * pc_x[k] * sfl1_312[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pb_x, pb_z, pc_x, pc_y, pc_z, sfl0_135, \
                         sfl0_314, sfk_143, sfl1_135, sfl1_314, \
                         sgk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_17 * sfk_143[k]
                   + f_3 * pc_y[k] * sgk_251[k];

        t_314[k] = pb_x[k] * sfl0_314[k]
                   - f_14 * pc_x[k] * sfl1_314[k];

        t_315[k] = pb_z[k] * sfl0_135[k]
                   - f_14 * pc_z[k] * sfl1_135[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pb_z, pc_y, pc_z, sfl0_138, sfk_108, \
                         sfk_144, sfk_146, sfl1_138, sgk_252, sgk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_16 * sfk_144[k]
                   + f_3 * pc_y[k] * sgk_252[k];

        t_317[k] = f_15 * sfk_108[k]
                   + f_3 * pc_z[k] * sgk_252[k];

        t_318[k] = pb_z[k] * sfl0_138[k]
                   - f_14 * pc_z[k] * sfl1_138[k];

        t_319[k] = f_16 * sfk_146[k]
                   + f_3 * pc_y[k] * sgk_254[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, pb_x, pb_z, pc_x, pc_z, sfl0_141, sfl0_320, \
                         sfk_111, sfk_257, sfl1_141, sfl1_320, \
                         sgk_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = pb_x[k] * sfl0_320[k]
                   + f_19 * sfk_257[k]
                   - f_14 * pc_x[k] * sfl1_320[k];

        t_321[k] = pb_z[k] * sfl0_141[k]
                   - f_14 * pc_z[k] * sfl1_141[k];

        t_322[k] = f_15 * sfk_111[k]
                   + f_3 * pc_z[k] * sgk_255[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, pb_x, pb_z, pc_x, pc_y, pc_z, sfl0_145, \
                         sfl0_324, sfk_149, sfk_261, sfl1_145, sfl1_324, \
                         sgk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_16 * sfk_149[k]
                   + f_3 * pc_y[k] * sgk_257[k];

        t_324[k] = pb_x[k] * sfl0_324[k]
                   + f_18 * sfk_261[k]
                   - f_14 * pc_x[k] * sfl1_324[k];

        t_325[k] = pb_z[k] * sfl0_145[k]
                   - f_14 * pc_z[k] * sfl1_145[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, pb_x, pc_x, pc_y, pc_z, sfl0_327, sfk_114, \
                         sfk_153, sfk_264, sfl1_327, sgk_258, sgk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_15 * sfk_114[k]
                   + f_3 * pc_z[k] * sgk_258[k];

        t_327[k] = pb_x[k] * sfl0_327[k]
                   + f_0 * sfk_264[k]
                   - f_14 * pc_x[k] * sfl1_327[k];

        t_328[k] = f_16 * sfk_153[k]
                   + f_3 * pc_y[k] * sgk_261[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pb_x, pb_z, pc_x, pc_z, sfl0_150, sfl0_329, \
                         sfk_118, sfk_266, sfl1_150, sfl1_329, \
                         sgk_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = pb_x[k] * sfl0_329[k]
                   + f_0 * sfk_266[k]
                   - f_14 * pc_x[k] * sfl1_329[k];

        t_330[k] = pb_z[k] * sfl0_150[k]
                   - f_14 * pc_z[k] * sfl1_150[k];

        t_331[k] = f_15 * sfk_118[k]
                   + f_3 * pc_z[k] * sgk_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pb_x, pc_x, pc_y, sfl0_332, sfl0_333, sfk_158, \
                         sfk_269, sfk_270, sfl1_332, sfl1_333, \
                         sgk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pb_x[k] * sfl0_332[k]
                   + f_17 * sfk_269[k]
                   - f_14 * pc_x[k] * sfl1_332[k];

        t_333[k] = pb_x[k] * sfl0_333[k]
                   + f_17 * sfk_270[k]
                   - f_14 * pc_x[k] * sfl1_333[k];

        t_334[k] = f_16 * sfk_158[k]
                   + f_3 * pc_y[k] * sgk_266[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pb_x, pb_z, pc_x, pc_z, sfl0_156, sfl0_335, \
                         sfk_123, sfk_272, sfl1_156, sfl1_335, \
                         sgk_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = pb_x[k] * sfl0_335[k]
                   + f_17 * sfk_272[k]
                   - f_14 * pc_x[k] * sfl1_335[k];

        t_336[k] = pb_z[k] * sfl0_156[k]
                   - f_14 * pc_z[k] * sfl1_156[k];

        t_337[k] = f_15 * sfk_123[k]
                   + f_3 * pc_z[k] * sgk_267[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pb_x, pc_x, sfl0_338, sfl0_339, sfl0_340, \
                         sfk_275, sfk_276, sfk_277, sfl1_338, sfl1_339, \
                         sfl1_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pb_x[k] * sfl0_338[k]
                   + f_16 * sfk_275[k]
                   - f_14 * pc_x[k] * sfl1_338[k];

        t_339[k] = pb_x[k] * sfl0_339[k]
                   + f_16 * sfk_276[k]
                   - f_14 * pc_x[k] * sfl1_339[k];

        t_340[k] = pb_x[k] * sfl0_340[k]
                   + f_16 * sfk_277[k]
                   - f_14 * pc_x[k] * sfl1_340[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pb_x, pc_x, pc_y, sfl0_342, sfk_164, \
                         sfk_279, sfk_280, sfk_281, sfl1_342, sgk_272, sgk_280, \
                         sgk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_16 * sfk_164[k]
                   + f_3 * pc_y[k] * sgk_272[k];

        t_342[k] = pb_x[k] * sfl0_342[k]
                   + f_16 * sfk_279[k]
                   - f_14 * pc_x[k] * sfl1_342[k];

        t_343[k] = f_15 * sfk_280[k]
                   + f_3 * pc_x[k] * sgk_280[k];

        t_344[k] = f_15 * sfk_281[k]
                   + f_3 * pc_x[k] * sgk_281[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pc_x, sfk_282, sfk_283, sfk_284, \
                         sfk_285, sfk_286, sgk_282, sgk_283, sgk_284, sgk_285, \
                         sgk_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_15 * sfk_282[k]
                   + f_3 * pc_x[k] * sgk_282[k];

        t_346[k] = f_15 * sfk_283[k]
                   + f_3 * pc_x[k] * sgk_283[k];

        t_347[k] = f_15 * sfk_284[k]
                   + f_3 * pc_x[k] * sgk_284[k];

        t_348[k] = f_15 * sfk_285[k]
                   + f_3 * pc_x[k] * sgk_285[k];

        t_349[k] = f_15 * sfk_286[k]
                   + f_3 * pc_x[k] * sgk_286[k];
    }
}

static auto
compute_prim_sgl_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfl0,
                                                          const size_t sfk, const size_t sfl1,
                                                          const size_t sgi0, const size_t sgi1,
                                                          const size_t sgk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 4.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfl0_225 = buffer.data(sfl0 + 225);
    const auto *sfl0_230 = buffer.data(sfl0 + 230);
    const auto *sfl0_234 = buffer.data(sfl0 + 234);
    const auto *sfl0_239 = buffer.data(sfl0 + 239);
    const auto *sfl0_245 = buffer.data(sfl0 + 245);
    const auto *sfl0_252 = buffer.data(sfl0 + 252);
    const auto *sfl0_351 = buffer.data(sfl0 + 351);
    const auto *sfl0_353 = buffer.data(sfl0 + 353);
    const auto *sfl0_354 = buffer.data(sfl0 + 354);
    const auto *sfl0_355 = buffer.data(sfl0 + 355);
    const auto *sfl0_356 = buffer.data(sfl0 + 356);
    const auto *sfl0_357 = buffer.data(sfl0 + 357);
    const auto *sfl0_359 = buffer.data(sfl0 + 359);
    const auto *sfl0_363 = buffer.data(sfl0 + 363);
    const auto *sfl0_366 = buffer.data(sfl0 + 366);
    const auto *sfl0_370 = buffer.data(sfl0 + 370);
    const auto *sfl0_372 = buffer.data(sfl0 + 372);
    const auto *sfl0_375 = buffer.data(sfl0 + 375);
    const auto *sfl0_377 = buffer.data(sfl0 + 377);
    const auto *sfl0_378 = buffer.data(sfl0 + 378);
    const auto *sfl0_381 = buffer.data(sfl0 + 381);
    const auto *sfl0_383 = buffer.data(sfl0 + 383);
    const auto *sfl0_384 = buffer.data(sfl0 + 384);
    const auto *sfl0_385 = buffer.data(sfl0 + 385);
    const auto *sfl0_396 = buffer.data(sfl0 + 396);
    const auto *sfl0_398 = buffer.data(sfl0 + 398);
    const auto *sfl0_399 = buffer.data(sfl0 + 399);
    const auto *sfl0_400 = buffer.data(sfl0 + 400);
    const auto *sfl0_401 = buffer.data(sfl0 + 401);
    const auto *sfl0_402 = buffer.data(sfl0 + 402);
    const auto *sfl0_404 = buffer.data(sfl0 + 404);
    const auto *sfl0_405 = buffer.data(sfl0 + 405);
    const auto *sfl0_408 = buffer.data(sfl0 + 408);
    const auto *sfl0_410 = buffer.data(sfl0 + 410);
    const auto *sfl0_411 = buffer.data(sfl0 + 411);
    const auto *sfl0_414 = buffer.data(sfl0 + 414);
    const auto *sfl0_415 = buffer.data(sfl0 + 415);
    const auto *sfl0_417 = buffer.data(sfl0 + 417);
    const auto *sfl0_419 = buffer.data(sfl0 + 419);
    const auto *sfl0_420 = buffer.data(sfl0 + 420);
    const auto *sfl0_422 = buffer.data(sfl0 + 422);
    const auto *sfl0_423 = buffer.data(sfl0 + 423);
    const auto *sfl0_425 = buffer.data(sfl0 + 425);
    const auto *sfl0_426 = buffer.data(sfl0 + 426);
    const auto *sfl0_428 = buffer.data(sfl0 + 428);
    const auto *sfl0_429 = buffer.data(sfl0 + 429);
    const auto *sfl0_430 = buffer.data(sfl0 + 430);
    const auto *sfl0_432 = buffer.data(sfl0 + 432);
    const auto *sfl0_441 = buffer.data(sfl0 + 441);
    const auto *sfl0_443 = buffer.data(sfl0 + 443);
    const auto *sfl0_444 = buffer.data(sfl0 + 444);
    const auto *sfl0_445 = buffer.data(sfl0 + 445);
    const auto *sfl0_446 = buffer.data(sfl0 + 446);
    const auto *sfl0_447 = buffer.data(sfl0 + 447);
    const auto *sfl0_449 = buffer.data(sfl0 + 449);

    const auto *sfk_136 = buffer.data(sfk + 136);
    const auto *sfk_144 = buffer.data(sfk + 144);
    const auto *sfk_147 = buffer.data(sfk + 147);
    const auto *sfk_150 = buffer.data(sfk + 150);
    const auto *sfk_154 = buffer.data(sfk + 154);
    const auto *sfk_159 = buffer.data(sfk + 159);
    const auto *sfk_172 = buffer.data(sfk + 172);
    const auto *sfk_179 = buffer.data(sfk + 179);
    const auto *sfk_180 = buffer.data(sfk + 180);
    const auto *sfk_182 = buffer.data(sfk + 182);
    const auto *sfk_183 = buffer.data(sfk + 183);
    const auto *sfk_185 = buffer.data(sfk + 185);
    const auto *sfk_186 = buffer.data(sfk + 186);
    const auto *sfk_189 = buffer.data(sfk + 189);
    const auto *sfk_190 = buffer.data(sfk + 190);
    const auto *sfk_194 = buffer.data(sfk + 194);
    const auto *sfk_195 = buffer.data(sfk + 195);
    const auto *sfk_200 = buffer.data(sfk + 200);
    const auto *sfk_208 = buffer.data(sfk + 208);
    const auto *sfk_215 = buffer.data(sfk + 215);
    const auto *sfk_216 = buffer.data(sfk + 216);
    const auto *sfk_218 = buffer.data(sfk + 218);
    const auto *sfk_221 = buffer.data(sfk + 221);
    const auto *sfk_225 = buffer.data(sfk + 225);
    const auto *sfk_230 = buffer.data(sfk + 230);
    const auto *sfk_287 = buffer.data(sfk + 287);
    const auto *sfk_291 = buffer.data(sfk + 291);
    const auto *sfk_294 = buffer.data(sfk + 294);
    const auto *sfk_298 = buffer.data(sfk + 298);
    const auto *sfk_300 = buffer.data(sfk + 300);
    const auto *sfk_303 = buffer.data(sfk + 303);
    const auto *sfk_305 = buffer.data(sfk + 305);
    const auto *sfk_306 = buffer.data(sfk + 306);
    const auto *sfk_309 = buffer.data(sfk + 309);
    const auto *sfk_311 = buffer.data(sfk + 311);
    const auto *sfk_312 = buffer.data(sfk + 312);
    const auto *sfk_313 = buffer.data(sfk + 313);
    const auto *sfk_316 = buffer.data(sfk + 316);
    const auto *sfk_317 = buffer.data(sfk + 317);
    const auto *sfk_318 = buffer.data(sfk + 318);
    const auto *sfk_319 = buffer.data(sfk + 319);
    const auto *sfk_320 = buffer.data(sfk + 320);
    const auto *sfk_321 = buffer.data(sfk + 321);
    const auto *sfk_322 = buffer.data(sfk + 322);
    const auto *sfk_323 = buffer.data(sfk + 323);
    const auto *sfk_324 = buffer.data(sfk + 324);
    const auto *sfk_327 = buffer.data(sfk + 327);
    const auto *sfk_329 = buffer.data(sfk + 329);
    const auto *sfk_330 = buffer.data(sfk + 330);
    const auto *sfk_333 = buffer.data(sfk + 333);
    const auto *sfk_334 = buffer.data(sfk + 334);
    const auto *sfk_336 = buffer.data(sfk + 336);
    const auto *sfk_338 = buffer.data(sfk + 338);
    const auto *sfk_339 = buffer.data(sfk + 339);
    const auto *sfk_341 = buffer.data(sfk + 341);
    const auto *sfk_342 = buffer.data(sfk + 342);
    const auto *sfk_344 = buffer.data(sfk + 344);
    const auto *sfk_345 = buffer.data(sfk + 345);
    const auto *sfk_347 = buffer.data(sfk + 347);
    const auto *sfk_348 = buffer.data(sfk + 348);
    const auto *sfk_349 = buffer.data(sfk + 349);
    const auto *sfk_351 = buffer.data(sfk + 351);
    const auto *sfk_352 = buffer.data(sfk + 352);
    const auto *sfk_353 = buffer.data(sfk + 353);
    const auto *sfk_354 = buffer.data(sfk + 354);
    const auto *sfk_355 = buffer.data(sfk + 355);
    const auto *sfk_356 = buffer.data(sfk + 356);
    const auto *sfk_357 = buffer.data(sfk + 357);
    const auto *sfk_358 = buffer.data(sfk + 358);
    const auto *sfk_359 = buffer.data(sfk + 359);

    const auto *sfl1_225 = buffer.data(sfl1 + 225);
    const auto *sfl1_230 = buffer.data(sfl1 + 230);
    const auto *sfl1_234 = buffer.data(sfl1 + 234);
    const auto *sfl1_239 = buffer.data(sfl1 + 239);
    const auto *sfl1_245 = buffer.data(sfl1 + 245);
    const auto *sfl1_252 = buffer.data(sfl1 + 252);
    const auto *sfl1_351 = buffer.data(sfl1 + 351);
    const auto *sfl1_353 = buffer.data(sfl1 + 353);
    const auto *sfl1_354 = buffer.data(sfl1 + 354);
    const auto *sfl1_355 = buffer.data(sfl1 + 355);
    const auto *sfl1_356 = buffer.data(sfl1 + 356);
    const auto *sfl1_357 = buffer.data(sfl1 + 357);
    const auto *sfl1_359 = buffer.data(sfl1 + 359);
    const auto *sfl1_363 = buffer.data(sfl1 + 363);
    const auto *sfl1_366 = buffer.data(sfl1 + 366);
    const auto *sfl1_370 = buffer.data(sfl1 + 370);
    const auto *sfl1_372 = buffer.data(sfl1 + 372);
    const auto *sfl1_375 = buffer.data(sfl1 + 375);
    const auto *sfl1_377 = buffer.data(sfl1 + 377);
    const auto *sfl1_378 = buffer.data(sfl1 + 378);
    const auto *sfl1_381 = buffer.data(sfl1 + 381);
    const auto *sfl1_383 = buffer.data(sfl1 + 383);
    const auto *sfl1_384 = buffer.data(sfl1 + 384);
    const auto *sfl1_385 = buffer.data(sfl1 + 385);
    const auto *sfl1_396 = buffer.data(sfl1 + 396);
    const auto *sfl1_398 = buffer.data(sfl1 + 398);
    const auto *sfl1_399 = buffer.data(sfl1 + 399);
    const auto *sfl1_400 = buffer.data(sfl1 + 400);
    const auto *sfl1_401 = buffer.data(sfl1 + 401);
    const auto *sfl1_402 = buffer.data(sfl1 + 402);
    const auto *sfl1_404 = buffer.data(sfl1 + 404);
    const auto *sfl1_405 = buffer.data(sfl1 + 405);
    const auto *sfl1_408 = buffer.data(sfl1 + 408);
    const auto *sfl1_410 = buffer.data(sfl1 + 410);
    const auto *sfl1_411 = buffer.data(sfl1 + 411);
    const auto *sfl1_414 = buffer.data(sfl1 + 414);
    const auto *sfl1_415 = buffer.data(sfl1 + 415);
    const auto *sfl1_417 = buffer.data(sfl1 + 417);
    const auto *sfl1_419 = buffer.data(sfl1 + 419);
    const auto *sfl1_420 = buffer.data(sfl1 + 420);
    const auto *sfl1_422 = buffer.data(sfl1 + 422);
    const auto *sfl1_423 = buffer.data(sfl1 + 423);
    const auto *sfl1_425 = buffer.data(sfl1 + 425);
    const auto *sfl1_426 = buffer.data(sfl1 + 426);
    const auto *sfl1_428 = buffer.data(sfl1 + 428);
    const auto *sfl1_429 = buffer.data(sfl1 + 429);
    const auto *sfl1_430 = buffer.data(sfl1 + 430);
    const auto *sfl1_432 = buffer.data(sfl1 + 432);
    const auto *sfl1_441 = buffer.data(sfl1 + 441);
    const auto *sfl1_443 = buffer.data(sfl1 + 443);
    const auto *sfl1_444 = buffer.data(sfl1 + 444);
    const auto *sfl1_445 = buffer.data(sfl1 + 445);
    const auto *sfl1_446 = buffer.data(sfl1 + 446);
    const auto *sfl1_447 = buffer.data(sfl1 + 447);
    const auto *sfl1_449 = buffer.data(sfl1 + 449);

    const auto *sgi0_280 = buffer.data(sgi0 + 280);
    const auto *sgi0_283 = buffer.data(sgi0 + 283);
    const auto *sgi0_285 = buffer.data(sgi0 + 285);
    const auto *sgi0_286 = buffer.data(sgi0 + 286);
    const auto *sgi0_289 = buffer.data(sgi0 + 289);
    const auto *sgi0_290 = buffer.data(sgi0 + 290);
    const auto *sgi0_292 = buffer.data(sgi0 + 292);
    const auto *sgi0_294 = buffer.data(sgi0 + 294);
    const auto *sgi0_295 = buffer.data(sgi0 + 295);
    const auto *sgi0_297 = buffer.data(sgi0 + 297);
    const auto *sgi0_298 = buffer.data(sgi0 + 298);
    const auto *sgi0_300 = buffer.data(sgi0 + 300);

    const auto *sgi1_280 = buffer.data(sgi1 + 280);
    const auto *sgi1_283 = buffer.data(sgi1 + 283);
    const auto *sgi1_285 = buffer.data(sgi1 + 285);
    const auto *sgi1_286 = buffer.data(sgi1 + 286);
    const auto *sgi1_289 = buffer.data(sgi1 + 289);
    const auto *sgi1_290 = buffer.data(sgi1 + 290);
    const auto *sgi1_292 = buffer.data(sgi1 + 292);
    const auto *sgi1_294 = buffer.data(sgi1 + 294);
    const auto *sgi1_295 = buffer.data(sgi1 + 295);
    const auto *sgi1_297 = buffer.data(sgi1 + 297);
    const auto *sgi1_298 = buffer.data(sgi1 + 298);
    const auto *sgi1_300 = buffer.data(sgi1 + 300);

    const auto *sgk_280 = buffer.data(sgk + 280);
    const auto *sgk_287 = buffer.data(sgk + 287);
    const auto *sgk_288 = buffer.data(sgk + 288);
    const auto *sgk_290 = buffer.data(sgk + 290);
    const auto *sgk_291 = buffer.data(sgk + 291);
    const auto *sgk_293 = buffer.data(sgk + 293);
    const auto *sgk_294 = buffer.data(sgk + 294);
    const auto *sgk_297 = buffer.data(sgk + 297);
    const auto *sgk_298 = buffer.data(sgk + 298);
    const auto *sgk_302 = buffer.data(sgk + 302);
    const auto *sgk_303 = buffer.data(sgk + 303);
    const auto *sgk_308 = buffer.data(sgk + 308);
    const auto *sgk_316 = buffer.data(sgk + 316);
    const auto *sgk_317 = buffer.data(sgk + 317);
    const auto *sgk_318 = buffer.data(sgk + 318);
    const auto *sgk_319 = buffer.data(sgk + 319);
    const auto *sgk_320 = buffer.data(sgk + 320);
    const auto *sgk_321 = buffer.data(sgk + 321);
    const auto *sgk_322 = buffer.data(sgk + 322);
    const auto *sgk_323 = buffer.data(sgk + 323);
    const auto *sgk_324 = buffer.data(sgk + 324);
    const auto *sgk_326 = buffer.data(sgk + 326);
    const auto *sgk_327 = buffer.data(sgk + 327);
    const auto *sgk_329 = buffer.data(sgk + 329);
    const auto *sgk_330 = buffer.data(sgk + 330);
    const auto *sgk_333 = buffer.data(sgk + 333);
    const auto *sgk_334 = buffer.data(sgk + 334);
    const auto *sgk_338 = buffer.data(sgk + 338);
    const auto *sgk_339 = buffer.data(sgk + 339);
    const auto *sgk_344 = buffer.data(sgk + 344);
    const auto *sgk_352 = buffer.data(sgk + 352);
    const auto *sgk_353 = buffer.data(sgk + 353);
    const auto *sgk_354 = buffer.data(sgk + 354);
    const auto *sgk_355 = buffer.data(sgk + 355);
    const auto *sgk_356 = buffer.data(sgk + 356);
    const auto *sgk_357 = buffer.data(sgk + 357);
    const auto *sgk_358 = buffer.data(sgk + 358);
    const auto *sgk_359 = buffer.data(sgk + 359);
    const auto *sgk_360 = buffer.data(sgk + 360);
    const auto *sgk_362 = buffer.data(sgk + 362);
    const auto *sgk_363 = buffer.data(sgk + 363);
    const auto *sgk_365 = buffer.data(sgk + 365);
    const auto *sgk_366 = buffer.data(sgk + 366);
    const auto *sgk_369 = buffer.data(sgk + 369);
    const auto *sgk_370 = buffer.data(sgk + 370);
    const auto *sgk_372 = buffer.data(sgk + 372);
    const auto *sgk_374 = buffer.data(sgk + 374);
    const auto *sgk_375 = buffer.data(sgk + 375);
    const auto *sgk_377 = buffer.data(sgk + 377);
    const auto *sgk_378 = buffer.data(sgk + 378);
    const auto *sgk_380 = buffer.data(sgk + 380);

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pb_x, pc_x, pc_z, sfl0_351, sfl0_353, \
                         sfk_136, sfk_287, sfl1_351, sfl1_353, sgk_280, \
                         sgk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_15 * sfk_287[k]
                   + f_3 * pc_x[k] * sgk_287[k];

        t_351[k] = pb_x[k] * sfl0_351[k]
                   - f_14 * pc_x[k] * sfl1_351[k];

        t_352[k] = f_15 * sfk_136[k]
                   + f_3 * pc_z[k] * sgk_280[k];

        t_353[k] = pb_x[k] * sfl0_353[k]
                   - f_14 * pc_x[k] * sfl1_353[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pb_x, pc_x, sfl0_354, sfl0_355, sfl0_356, \
                         sfl0_357, sfl1_354, sfl1_355, sfl1_356, \
                         sfl1_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = pb_x[k] * sfl0_354[k]
                   - f_14 * pc_x[k] * sfl1_354[k];

        t_355[k] = pb_x[k] * sfl0_355[k]
                   - f_14 * pc_x[k] * sfl1_355[k];

        t_356[k] = pb_x[k] * sfl0_356[k]
                   - f_14 * pc_x[k] * sfl1_356[k];

        t_357[k] = pb_x[k] * sfl0_357[k]
                   - f_14 * pc_x[k] * sfl1_357[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pb_x, pb_y, pc_x, pc_y, sfl0_225, \
                         sfl0_359, sfk_179, sfk_180, sfl1_225, sfl1_359, sgk_287, \
                         sgk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_16 * sfk_179[k]
                   + f_3 * pc_y[k] * sgk_287[k];

        t_359[k] = pb_x[k] * sfl0_359[k]
                   - f_14 * pc_x[k] * sfl1_359[k];

        t_360[k] = pb_y[k] * sfl0_225[k]
                   - f_14 * pc_y[k] * sfl1_225[k];

        t_361[k] = f_15 * sfk_180[k]
                   + f_3 * pc_y[k] * sgk_288[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pb_x, pc_x, pc_y, pc_z, sfl0_363, sfk_144, \
                         sfk_182, sfk_291, sfl1_363, sgk_288, sgk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_16 * sfk_144[k]
                   + f_3 * pc_z[k] * sgk_288[k];

        t_363[k] = pb_x[k] * sfl0_363[k]
                   + f_19 * sfk_291[k]
                   - f_14 * pc_x[k] * sfl1_363[k];

        t_364[k] = f_15 * sfk_182[k]
                   + f_3 * pc_y[k] * sgk_290[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, pb_x, pb_y, pc_x, pc_y, pc_z, sfl0_230, \
                         sfl0_366, sfk_147, sfk_294, sfl1_230, sfl1_366, \
                         sgk_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = pb_y[k] * sfl0_230[k]
                   - f_14 * pc_y[k] * sfl1_230[k];

        t_366[k] = pb_x[k] * sfl0_366[k]
                   + f_18 * sfk_294[k]
                   - f_14 * pc_x[k] * sfl1_366[k];

        t_367[k] = f_16 * sfk_147[k]
                   + f_3 * pc_z[k] * sgk_291[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, pb_x, pb_y, pc_x, pc_y, sfl0_234, sfl0_370, \
                         sfk_185, sfk_298, sfl1_234, sfl1_370, \
                         sgk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_15 * sfk_185[k]
                   + f_3 * pc_y[k] * sgk_293[k];

        t_369[k] = pb_y[k] * sfl0_234[k]
                   - f_14 * pc_y[k] * sfl1_234[k];

        t_370[k] = pb_x[k] * sfl0_370[k]
                   + f_0 * sfk_298[k]
                   - f_14 * pc_x[k] * sfl1_370[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, pb_x, pc_x, pc_y, pc_z, sfl0_372, sfk_150, \
                         sfk_189, sfk_300, sfl1_372, sgk_294, sgk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_16 * sfk_150[k]
                   + f_3 * pc_z[k] * sgk_294[k];

        t_372[k] = pb_x[k] * sfl0_372[k]
                   + f_0 * sfk_300[k]
                   - f_14 * pc_x[k] * sfl1_372[k];

        t_373[k] = f_15 * sfk_189[k]
                   + f_3 * pc_y[k] * sgk_297[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pb_x, pb_y, pc_x, pc_y, pc_z, sfl0_239, \
                         sfl0_375, sfk_154, sfk_303, sfl1_239, sfl1_375, \
                         sgk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pb_y[k] * sfl0_239[k]
                   - f_14 * pc_y[k] * sfl1_239[k];

        t_375[k] = pb_x[k] * sfl0_375[k]
                   + f_17 * sfk_303[k]
                   - f_14 * pc_x[k] * sfl1_375[k];

        t_376[k] = f_16 * sfk_154[k]
                   + f_3 * pc_z[k] * sgk_298[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, pb_x, pc_x, pc_y, sfl0_377, sfl0_378, sfk_194, \
                         sfk_305, sfk_306, sfl1_377, sfl1_378, \
                         sgk_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = pb_x[k] * sfl0_377[k]
                   + f_17 * sfk_305[k]
                   - f_14 * pc_x[k] * sfl1_377[k];

        t_378[k] = pb_x[k] * sfl0_378[k]
                   + f_17 * sfk_306[k]
                   - f_14 * pc_x[k] * sfl1_378[k];

        t_379[k] = f_15 * sfk_194[k]
                   + f_3 * pc_y[k] * sgk_302[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, pb_x, pb_y, pc_x, pc_y, pc_z, sfl0_245, \
                         sfl0_381, sfk_159, sfk_309, sfl1_245, sfl1_381, \
                         sgk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = pb_y[k] * sfl0_245[k]
                   - f_14 * pc_y[k] * sfl1_245[k];

        t_381[k] = pb_x[k] * sfl0_381[k]
                   + f_16 * sfk_309[k]
                   - f_14 * pc_x[k] * sfl1_381[k];

        t_382[k] = f_16 * sfk_159[k]
                   + f_3 * pc_z[k] * sgk_303[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pb_x, pc_x, sfl0_383, sfl0_384, sfl0_385, \
                         sfk_311, sfk_312, sfk_313, sfl1_383, sfl1_384, \
                         sfl1_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = pb_x[k] * sfl0_383[k]
                   + f_16 * sfk_311[k]
                   - f_14 * pc_x[k] * sfl1_383[k];

        t_384[k] = pb_x[k] * sfl0_384[k]
                   + f_16 * sfk_312[k]
                   - f_14 * pc_x[k] * sfl1_384[k];

        t_385[k] = pb_x[k] * sfl0_385[k]
                   + f_16 * sfk_313[k]
                   - f_14 * pc_x[k] * sfl1_385[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pb_y, pc_x, pc_y, sfl0_252, sfk_200, \
                         sfk_316, sfk_317, sfl1_252, sgk_308, sgk_316, \
                         sgk_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_15 * sfk_200[k]
                   + f_3 * pc_y[k] * sgk_308[k];

        t_387[k] = pb_y[k] * sfl0_252[k]
                   - f_14 * pc_y[k] * sfl1_252[k];

        t_388[k] = f_15 * sfk_316[k]
                   + f_3 * pc_x[k] * sgk_316[k];

        t_389[k] = f_15 * sfk_317[k]
                   + f_3 * pc_x[k] * sgk_317[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, pc_x, sfk_318, sfk_319, sfk_320, \
                         sfk_321, sfk_322, sgk_318, sgk_319, sgk_320, sgk_321, \
                         sgk_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_15 * sfk_318[k]
                   + f_3 * pc_x[k] * sgk_318[k];

        t_391[k] = f_15 * sfk_319[k]
                   + f_3 * pc_x[k] * sgk_319[k];

        t_392[k] = f_15 * sfk_320[k]
                   + f_3 * pc_x[k] * sgk_320[k];

        t_393[k] = f_15 * sfk_321[k]
                   + f_3 * pc_x[k] * sgk_321[k];

        t_394[k] = f_15 * sfk_322[k]
                   + f_3 * pc_x[k] * sgk_322[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pb_x, pc_x, pc_z, sfl0_396, sfl0_398, \
                         sfk_172, sfk_323, sfl1_396, sfl1_398, sgk_316, \
                         sgk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_15 * sfk_323[k]
                   + f_3 * pc_x[k] * sgk_323[k];

        t_396[k] = pb_x[k] * sfl0_396[k]
                   - f_14 * pc_x[k] * sfl1_396[k];

        t_397[k] = f_16 * sfk_172[k]
                   + f_3 * pc_z[k] * sgk_316[k];

        t_398[k] = pb_x[k] * sfl0_398[k]
                   - f_14 * pc_x[k] * sfl1_398[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, pb_x, pc_x, sfl0_399, sfl0_400, sfl0_401, \
                         sfl0_402, sfl1_399, sfl1_400, sfl1_401, \
                         sfl1_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = pb_x[k] * sfl0_399[k]
                   - f_14 * pc_x[k] * sfl1_399[k];

        t_400[k] = pb_x[k] * sfl0_400[k]
                   - f_14 * pc_x[k] * sfl1_400[k];

        t_401[k] = pb_x[k] * sfl0_401[k]
                   - f_14 * pc_x[k] * sfl1_401[k];

        t_402[k] = pb_x[k] * sfl0_402[k]
                   - f_14 * pc_x[k] * sfl1_402[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, pb_x, pc_x, pc_y, sfl0_404, sfl0_405, \
                         sfk_215, sfk_324, sfl1_404, sfl1_405, sgk_323, \
                         sgk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_15 * sfk_215[k]
                   + f_3 * pc_y[k] * sgk_323[k];

        t_404[k] = pb_x[k] * sfl0_404[k]
                   - f_14 * pc_x[k] * sfl1_404[k];

        t_405[k] = pb_x[k] * sfl0_405[k]
                   + f_20 * sfk_324[k]
                   - f_14 * pc_x[k] * sfl1_405[k];

        t_406[k] = f_3 * pc_y[k] * sgk_324[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pb_x, pc_x, pc_y, pc_z, sfl0_408, sfk_180, \
                         sfk_327, sfl1_408, sgk_324, sgk_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_17 * sfk_180[k]
                   + f_3 * pc_z[k] * sgk_324[k];

        t_408[k] = pb_x[k] * sfl0_408[k]
                   + f_19 * sfk_327[k]
                   - f_14 * pc_x[k] * sfl1_408[k];

        t_409[k] = f_3 * pc_y[k] * sgk_326[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pb_x, pc_x, pc_z, sfl0_410, sfl0_411, sfk_183, \
                         sfk_329, sfk_330, sfl1_410, sfl1_411, \
                         sgk_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = pb_x[k] * sfl0_410[k]
                   + f_19 * sfk_329[k]
                   - f_14 * pc_x[k] * sfl1_410[k];

        t_411[k] = pb_x[k] * sfl0_411[k]
                   + f_18 * sfk_330[k]
                   - f_14 * pc_x[k] * sfl1_411[k];

        t_412[k] = f_17 * sfk_183[k]
                   + f_3 * pc_z[k] * sgk_327[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pb_x, pc_x, pc_y, sfl0_414, sfl0_415, sfk_333, \
                         sfk_334, sfl1_414, sfl1_415, sgk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_3 * pc_y[k] * sgk_329[k];

        t_414[k] = pb_x[k] * sfl0_414[k]
                   + f_18 * sfk_333[k]
                   - f_14 * pc_x[k] * sfl1_414[k];

        t_415[k] = pb_x[k] * sfl0_415[k]
                   + f_0 * sfk_334[k]
                   - f_14 * pc_x[k] * sfl1_415[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, pb_x, pc_x, pc_y, pc_z, sfl0_417, sfk_186, \
                         sfk_336, sfl1_417, sgk_330, sgk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_17 * sfk_186[k]
                   + f_3 * pc_z[k] * sgk_330[k];

        t_417[k] = pb_x[k] * sfl0_417[k]
                   + f_0 * sfk_336[k]
                   - f_14 * pc_x[k] * sfl1_417[k];

        t_418[k] = f_3 * pc_y[k] * sgk_333[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, pb_x, pc_x, pc_z, sfl0_419, sfl0_420, sfk_190, \
                         sfk_338, sfk_339, sfl1_419, sfl1_420, \
                         sgk_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = pb_x[k] * sfl0_419[k]
                   + f_0 * sfk_338[k]
                   - f_14 * pc_x[k] * sfl1_419[k];

        t_420[k] = pb_x[k] * sfl0_420[k]
                   + f_17 * sfk_339[k]
                   - f_14 * pc_x[k] * sfl1_420[k];

        t_421[k] = f_17 * sfk_190[k]
                   + f_3 * pc_z[k] * sgk_334[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, pb_x, pc_x, pc_y, sfl0_422, sfl0_423, sfk_341, \
                         sfk_342, sfl1_422, sfl1_423, sgk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = pb_x[k] * sfl0_422[k]
                   + f_17 * sfk_341[k]
                   - f_14 * pc_x[k] * sfl1_422[k];

        t_423[k] = pb_x[k] * sfl0_423[k]
                   + f_17 * sfk_342[k]
                   - f_14 * pc_x[k] * sfl1_423[k];

        t_424[k] = f_3 * pc_y[k] * sgk_338[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, pb_x, pc_x, pc_z, sfl0_425, sfl0_426, sfk_195, \
                         sfk_344, sfk_345, sfl1_425, sfl1_426, \
                         sgk_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = pb_x[k] * sfl0_425[k]
                   + f_17 * sfk_344[k]
                   - f_14 * pc_x[k] * sfl1_425[k];

        t_426[k] = pb_x[k] * sfl0_426[k]
                   + f_16 * sfk_345[k]
                   - f_14 * pc_x[k] * sfl1_426[k];

        t_427[k] = f_17 * sfk_195[k]
                   + f_3 * pc_z[k] * sgk_339[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, pb_x, pc_x, sfl0_428, sfl0_429, sfl0_430, \
                         sfk_347, sfk_348, sfk_349, sfl1_428, sfl1_429, \
                         sfl1_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = pb_x[k] * sfl0_428[k]
                   + f_16 * sfk_347[k]
                   - f_14 * pc_x[k] * sfl1_428[k];

        t_429[k] = pb_x[k] * sfl0_429[k]
                   + f_16 * sfk_348[k]
                   - f_14 * pc_x[k] * sfl1_429[k];

        t_430[k] = pb_x[k] * sfl0_430[k]
                   + f_16 * sfk_349[k]
                   - f_14 * pc_x[k] * sfl1_430[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pb_x, pc_x, pc_y, sfl0_432, sfk_351, \
                         sfk_352, sfk_353, sfl1_432, sgk_344, sgk_352, \
                         sgk_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_3 * pc_y[k] * sgk_344[k];

        t_432[k] = pb_x[k] * sfl0_432[k]
                   + f_16 * sfk_351[k]
                   - f_14 * pc_x[k] * sfl1_432[k];

        t_433[k] = f_15 * sfk_352[k]
                   + f_3 * pc_x[k] * sgk_352[k];

        t_434[k] = f_15 * sfk_353[k]
                   + f_3 * pc_x[k] * sgk_353[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, pc_x, sfk_354, sfk_355, sfk_356, \
                         sfk_357, sfk_358, sgk_354, sgk_355, sgk_356, sgk_357, \
                         sgk_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_15 * sfk_354[k]
                   + f_3 * pc_x[k] * sgk_354[k];

        t_436[k] = f_15 * sfk_355[k]
                   + f_3 * pc_x[k] * sgk_355[k];

        t_437[k] = f_15 * sfk_356[k]
                   + f_3 * pc_x[k] * sgk_356[k];

        t_438[k] = f_15 * sfk_357[k]
                   + f_3 * pc_x[k] * sgk_357[k];

        t_439[k] = f_15 * sfk_358[k]
                   + f_3 * pc_x[k] * sgk_358[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, pb_x, pc_x, pc_z, sfl0_441, sfl0_443, \
                         sfk_208, sfk_359, sfl1_441, sfl1_443, sgk_352, \
                         sgk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_15 * sfk_359[k]
                   + f_3 * pc_x[k] * sgk_359[k];

        t_441[k] = pb_x[k] * sfl0_441[k]
                   - f_14 * pc_x[k] * sfl1_441[k];

        t_442[k] = f_17 * sfk_208[k]
                   + f_3 * pc_z[k] * sgk_352[k];

        t_443[k] = pb_x[k] * sfl0_443[k]
                   - f_14 * pc_x[k] * sfl1_443[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, pb_x, pc_x, sfl0_444, sfl0_445, sfl0_446, \
                         sfl0_447, sfl1_444, sfl1_445, sfl1_446, \
                         sfl1_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = pb_x[k] * sfl0_444[k]
                   - f_14 * pc_x[k] * sfl1_444[k];

        t_445[k] = pb_x[k] * sfl0_445[k]
                   - f_14 * pc_x[k] * sfl1_445[k];

        t_446[k] = pb_x[k] * sfl0_446[k]
                   - f_14 * pc_x[k] * sfl1_446[k];

        t_447[k] = pb_x[k] * sfl0_447[k]
                   - f_14 * pc_x[k] * sfl1_447[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, t_452, pb_x, pc_x, pc_y, pc_z, sfl0_449, \
                         sfk_216, sfl1_449, sgi0_280, sgi1_280, sgk_359, \
                         sgk_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = f_3 * pc_y[k] * sgk_359[k];

        t_449[k] = pb_x[k] * sfl0_449[k]
                   - f_14 * pc_x[k] * sfl1_449[k];

        t_450[k] = f_1 * sgi0_280[k]
                   - f_2 * sgi1_280[k]
                   + f_3 * pc_x[k] * sgk_360[k];

        t_451[k] = f_0 * sfk_216[k]
                   + f_3 * pc_y[k] * sgk_360[k];

        t_452[k] = f_3 * pc_z[k] * sgk_360[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pc_x, pc_y, sfk_218, sgi0_283, sgi0_285, \
                         sgi1_283, sgi1_285, sgk_362, sgk_363, \
                         sgk_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_4 * sgi0_283[k]
                   - f_5 * sgi1_283[k]
                   + f_3 * pc_x[k] * sgk_363[k];

        t_454[k] = f_0 * sfk_218[k]
                   + f_3 * pc_y[k] * sgk_362[k];

        t_455[k] = f_4 * sgi0_285[k]
                   - f_5 * sgi1_285[k]
                   + f_3 * pc_x[k] * sgk_365[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, pc_y, pc_z, sfk_221, sgi0_286, \
                         sgi0_289, sgi1_286, sgi1_289, sgk_363, sgk_365, sgk_366, \
                         sgk_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_6 * sgi0_286[k]
                   - f_7 * sgi1_286[k]
                   + f_3 * pc_x[k] * sgk_366[k];

        t_457[k] = f_3 * pc_z[k] * sgk_363[k];

        t_458[k] = f_0 * sfk_221[k]
                   + f_3 * pc_y[k] * sgk_365[k];

        t_459[k] = f_6 * sgi0_289[k]
                   - f_7 * sgi1_289[k]
                   + f_3 * pc_x[k] * sgk_369[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, pc_x, pc_y, pc_z, sfk_225, sgi0_290, \
                         sgi0_292, sgi1_290, sgi1_292, sgk_366, sgk_369, sgk_370, \
                         sgk_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_8 * sgi0_290[k]
                   - f_9 * sgi1_290[k]
                   + f_3 * pc_x[k] * sgk_370[k];

        t_461[k] = f_3 * pc_z[k] * sgk_366[k];

        t_462[k] = f_8 * sgi0_292[k]
                   - f_9 * sgi1_292[k]
                   + f_3 * pc_x[k] * sgk_372[k];

        t_463[k] = f_0 * sfk_225[k]
                   + f_3 * pc_y[k] * sgk_369[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, pc_x, pc_z, sgi0_294, sgi0_295, sgi0_297, \
                         sgi1_294, sgi1_295, sgi1_297, sgk_370, sgk_374, sgk_375, \
                         sgk_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_8 * sgi0_294[k]
                   - f_9 * sgi1_294[k]
                   + f_3 * pc_x[k] * sgk_374[k];

        t_465[k] = f_10 * sgi0_295[k]
                   - f_11 * sgi1_295[k]
                   + f_3 * pc_x[k] * sgk_375[k];

        t_466[k] = f_3 * pc_z[k] * sgk_370[k];

        t_467[k] = f_10 * sgi0_297[k]
                   - f_11 * sgi1_297[k]
                   + f_3 * pc_x[k] * sgk_377[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_x, pc_y, sfk_230, sgi0_298, sgi0_300, \
                         sgi1_298, sgi1_300, sgk_374, sgk_378, \
                         sgk_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_10 * sgi0_298[k]
                   - f_11 * sgi1_298[k]
                   + f_3 * pc_x[k] * sgk_378[k];

        t_469[k] = f_0 * sfk_230[k]
                   + f_3 * pc_y[k] * sgk_374[k];

        t_470[k] = f_10 * sgi0_300[k]
                   - f_11 * sgi1_300[k]
                   + f_3 * pc_x[k] * sgk_380[k];
    }
}

static auto
compute_prim_sgl_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfl0,
                                                          const size_t sfk, const size_t sfl1,
                                                          const size_t sgi0, const size_t sgi1,
                                                          const size_t sgk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfl0_270 = buffer.data(sfl0 + 270);
    const auto *sfl0_273 = buffer.data(sfl0 + 273);
    const auto *sfl0_276 = buffer.data(sfl0 + 276);
    const auto *sfl0_280 = buffer.data(sfl0 + 280);
    const auto *sfl0_285 = buffer.data(sfl0 + 285);
    const auto *sfl0_291 = buffer.data(sfl0 + 291);
    const auto *sfl0_306 = buffer.data(sfl0 + 306);
    const auto *sfl0_308 = buffer.data(sfl0 + 308);
    const auto *sfl0_309 = buffer.data(sfl0 + 309);
    const auto *sfl0_310 = buffer.data(sfl0 + 310);
    const auto *sfl0_311 = buffer.data(sfl0 + 311);
    const auto *sfl0_312 = buffer.data(sfl0 + 312);
    const auto *sfl0_405 = buffer.data(sfl0 + 405);
    const auto *sfl0_410 = buffer.data(sfl0 + 410);

    const auto *sfk_216 = buffer.data(sfk + 216);
    const auto *sfk_219 = buffer.data(sfk + 219);
    const auto *sfk_222 = buffer.data(sfk + 222);
    const auto *sfk_226 = buffer.data(sfk + 226);
    const auto *sfk_231 = buffer.data(sfk + 231);
    const auto *sfk_236 = buffer.data(sfk + 236);
    const auto *sfk_244 = buffer.data(sfk + 244);
    const auto *sfk_245 = buffer.data(sfk + 245);
    const auto *sfk_246 = buffer.data(sfk + 246);
    const auto *sfk_247 = buffer.data(sfk + 247);
    const auto *sfk_248 = buffer.data(sfk + 248);
    const auto *sfk_249 = buffer.data(sfk + 249);
    const auto *sfk_250 = buffer.data(sfk + 250);
    const auto *sfk_251 = buffer.data(sfk + 251);
    const auto *sfk_252 = buffer.data(sfk + 252);
    const auto *sfk_254 = buffer.data(sfk + 254);
    const auto *sfk_255 = buffer.data(sfk + 255);
    const auto *sfk_257 = buffer.data(sfk + 257);
    const auto *sfk_258 = buffer.data(sfk + 258);
    const auto *sfk_261 = buffer.data(sfk + 261);
    const auto *sfk_262 = buffer.data(sfk + 262);
    const auto *sfk_266 = buffer.data(sfk + 266);
    const auto *sfk_267 = buffer.data(sfk + 267);
    const auto *sfk_272 = buffer.data(sfk + 272);
    const auto *sfk_280 = buffer.data(sfk + 280);
    const auto *sfk_287 = buffer.data(sfk + 287);
    const auto *sfk_288 = buffer.data(sfk + 288);
    const auto *sfk_290 = buffer.data(sfk + 290);
    const auto *sfk_291 = buffer.data(sfk + 291);
    const auto *sfk_293 = buffer.data(sfk + 293);
    const auto *sfk_297 = buffer.data(sfk + 297);
    const auto *sfk_302 = buffer.data(sfk + 302);
    const auto *sfk_308 = buffer.data(sfk + 308);
    const auto *sfk_316 = buffer.data(sfk + 316);
    const auto *sfk_318 = buffer.data(sfk + 318);
    const auto *sfk_319 = buffer.data(sfk + 319);
    const auto *sfk_320 = buffer.data(sfk + 320);
    const auto *sfk_321 = buffer.data(sfk + 321);
    const auto *sfk_322 = buffer.data(sfk + 322);
    const auto *sfk_323 = buffer.data(sfk + 323);
    const auto *sfk_324 = buffer.data(sfk + 324);
    const auto *sfk_326 = buffer.data(sfk + 326);
    const auto *sfk_329 = buffer.data(sfk + 329);

    const auto *sfl1_270 = buffer.data(sfl1 + 270);
    const auto *sfl1_273 = buffer.data(sfl1 + 273);
    const auto *sfl1_276 = buffer.data(sfl1 + 276);
    const auto *sfl1_280 = buffer.data(sfl1 + 280);
    const auto *sfl1_285 = buffer.data(sfl1 + 285);
    const auto *sfl1_291 = buffer.data(sfl1 + 291);
    const auto *sfl1_306 = buffer.data(sfl1 + 306);
    const auto *sfl1_308 = buffer.data(sfl1 + 308);
    const auto *sfl1_309 = buffer.data(sfl1 + 309);
    const auto *sfl1_310 = buffer.data(sfl1 + 310);
    const auto *sfl1_311 = buffer.data(sfl1 + 311);
    const auto *sfl1_312 = buffer.data(sfl1 + 312);
    const auto *sfl1_405 = buffer.data(sfl1 + 405);
    const auto *sfl1_410 = buffer.data(sfl1 + 410);

    const auto *sgi0_301 = buffer.data(sgi0 + 301);
    const auto *sgi0_303 = buffer.data(sgi0 + 303);
    const auto *sgi0_304 = buffer.data(sgi0 + 304);
    const auto *sgi0_305 = buffer.data(sgi0 + 305);
    const auto *sgi0_306 = buffer.data(sgi0 + 306);
    const auto *sgi0_307 = buffer.data(sgi0 + 307);
    const auto *sgi0_313 = buffer.data(sgi0 + 313);
    const auto *sgi0_317 = buffer.data(sgi0 + 317);
    const auto *sgi0_320 = buffer.data(sgi0 + 320);
    const auto *sgi0_322 = buffer.data(sgi0 + 322);
    const auto *sgi0_325 = buffer.data(sgi0 + 325);
    const auto *sgi0_326 = buffer.data(sgi0 + 326);
    const auto *sgi0_328 = buffer.data(sgi0 + 328);
    const auto *sgi0_331 = buffer.data(sgi0 + 331);
    const auto *sgi0_332 = buffer.data(sgi0 + 332);
    const auto *sgi0_333 = buffer.data(sgi0 + 333);
    const auto *sgi0_335 = buffer.data(sgi0 + 335);
    const auto *sgi0_336 = buffer.data(sgi0 + 336);
    const auto *sgi0_339 = buffer.data(sgi0 + 339);
    const auto *sgi0_341 = buffer.data(sgi0 + 341);
    const auto *sgi0_342 = buffer.data(sgi0 + 342);
    const auto *sgi0_345 = buffer.data(sgi0 + 345);
    const auto *sgi0_346 = buffer.data(sgi0 + 346);
    const auto *sgi0_348 = buffer.data(sgi0 + 348);
    const auto *sgi0_350 = buffer.data(sgi0 + 350);
    const auto *sgi0_351 = buffer.data(sgi0 + 351);
    const auto *sgi0_353 = buffer.data(sgi0 + 353);
    const auto *sgi0_354 = buffer.data(sgi0 + 354);
    const auto *sgi0_356 = buffer.data(sgi0 + 356);
    const auto *sgi0_357 = buffer.data(sgi0 + 357);
    const auto *sgi0_359 = buffer.data(sgi0 + 359);
    const auto *sgi0_360 = buffer.data(sgi0 + 360);
    const auto *sgi0_361 = buffer.data(sgi0 + 361);
    const auto *sgi0_362 = buffer.data(sgi0 + 362);
    const auto *sgi0_363 = buffer.data(sgi0 + 363);
    const auto *sgi0_367 = buffer.data(sgi0 + 367);
    const auto *sgi0_370 = buffer.data(sgi0 + 370);

    const auto *sgi1_301 = buffer.data(sgi1 + 301);
    const auto *sgi1_303 = buffer.data(sgi1 + 303);
    const auto *sgi1_304 = buffer.data(sgi1 + 304);
    const auto *sgi1_305 = buffer.data(sgi1 + 305);
    const auto *sgi1_306 = buffer.data(sgi1 + 306);
    const auto *sgi1_307 = buffer.data(sgi1 + 307);
    const auto *sgi1_313 = buffer.data(sgi1 + 313);
    const auto *sgi1_317 = buffer.data(sgi1 + 317);
    const auto *sgi1_320 = buffer.data(sgi1 + 320);
    const auto *sgi1_322 = buffer.data(sgi1 + 322);
    const auto *sgi1_325 = buffer.data(sgi1 + 325);
    const auto *sgi1_326 = buffer.data(sgi1 + 326);
    const auto *sgi1_328 = buffer.data(sgi1 + 328);
    const auto *sgi1_331 = buffer.data(sgi1 + 331);
    const auto *sgi1_332 = buffer.data(sgi1 + 332);
    const auto *sgi1_333 = buffer.data(sgi1 + 333);
    const auto *sgi1_335 = buffer.data(sgi1 + 335);
    const auto *sgi1_336 = buffer.data(sgi1 + 336);
    const auto *sgi1_339 = buffer.data(sgi1 + 339);
    const auto *sgi1_341 = buffer.data(sgi1 + 341);
    const auto *sgi1_342 = buffer.data(sgi1 + 342);
    const auto *sgi1_345 = buffer.data(sgi1 + 345);
    const auto *sgi1_346 = buffer.data(sgi1 + 346);
    const auto *sgi1_348 = buffer.data(sgi1 + 348);
    const auto *sgi1_350 = buffer.data(sgi1 + 350);
    const auto *sgi1_351 = buffer.data(sgi1 + 351);
    const auto *sgi1_353 = buffer.data(sgi1 + 353);
    const auto *sgi1_354 = buffer.data(sgi1 + 354);
    const auto *sgi1_356 = buffer.data(sgi1 + 356);
    const auto *sgi1_357 = buffer.data(sgi1 + 357);
    const auto *sgi1_359 = buffer.data(sgi1 + 359);
    const auto *sgi1_360 = buffer.data(sgi1 + 360);
    const auto *sgi1_361 = buffer.data(sgi1 + 361);
    const auto *sgi1_362 = buffer.data(sgi1 + 362);
    const auto *sgi1_363 = buffer.data(sgi1 + 363);
    const auto *sgi1_367 = buffer.data(sgi1 + 367);
    const auto *sgi1_370 = buffer.data(sgi1 + 370);

    const auto *sgk_375 = buffer.data(sgk + 375);
    const auto *sgk_380 = buffer.data(sgk + 380);
    const auto *sgk_381 = buffer.data(sgk + 381);
    const auto *sgk_383 = buffer.data(sgk + 383);
    const auto *sgk_384 = buffer.data(sgk + 384);
    const auto *sgk_385 = buffer.data(sgk + 385);
    const auto *sgk_387 = buffer.data(sgk + 387);
    const auto *sgk_388 = buffer.data(sgk + 388);
    const auto *sgk_389 = buffer.data(sgk + 389);
    const auto *sgk_390 = buffer.data(sgk + 390);
    const auto *sgk_391 = buffer.data(sgk + 391);
    const auto *sgk_392 = buffer.data(sgk + 392);
    const auto *sgk_393 = buffer.data(sgk + 393);
    const auto *sgk_394 = buffer.data(sgk + 394);
    const auto *sgk_395 = buffer.data(sgk + 395);
    const auto *sgk_396 = buffer.data(sgk + 396);
    const auto *sgk_398 = buffer.data(sgk + 398);
    const auto *sgk_399 = buffer.data(sgk + 399);
    const auto *sgk_401 = buffer.data(sgk + 401);
    const auto *sgk_402 = buffer.data(sgk + 402);
    const auto *sgk_405 = buffer.data(sgk + 405);
    const auto *sgk_406 = buffer.data(sgk + 406);
    const auto *sgk_408 = buffer.data(sgk + 408);
    const auto *sgk_410 = buffer.data(sgk + 410);
    const auto *sgk_411 = buffer.data(sgk + 411);
    const auto *sgk_413 = buffer.data(sgk + 413);
    const auto *sgk_414 = buffer.data(sgk + 414);
    const auto *sgk_416 = buffer.data(sgk + 416);
    const auto *sgk_419 = buffer.data(sgk + 419);
    const auto *sgk_420 = buffer.data(sgk + 420);
    const auto *sgk_421 = buffer.data(sgk + 421);
    const auto *sgk_423 = buffer.data(sgk + 423);
    const auto *sgk_424 = buffer.data(sgk + 424);
    const auto *sgk_425 = buffer.data(sgk + 425);
    const auto *sgk_426 = buffer.data(sgk + 426);
    const auto *sgk_427 = buffer.data(sgk + 427);
    const auto *sgk_428 = buffer.data(sgk + 428);
    const auto *sgk_429 = buffer.data(sgk + 429);
    const auto *sgk_430 = buffer.data(sgk + 430);
    const auto *sgk_431 = buffer.data(sgk + 431);
    const auto *sgk_432 = buffer.data(sgk + 432);
    const auto *sgk_434 = buffer.data(sgk + 434);
    const auto *sgk_435 = buffer.data(sgk + 435);
    const auto *sgk_437 = buffer.data(sgk + 437);
    const auto *sgk_438 = buffer.data(sgk + 438);
    const auto *sgk_441 = buffer.data(sgk + 441);
    const auto *sgk_442 = buffer.data(sgk + 442);
    const auto *sgk_444 = buffer.data(sgk + 444);
    const auto *sgk_446 = buffer.data(sgk + 446);
    const auto *sgk_447 = buffer.data(sgk + 447);
    const auto *sgk_449 = buffer.data(sgk + 449);
    const auto *sgk_450 = buffer.data(sgk + 450);
    const auto *sgk_452 = buffer.data(sgk + 452);
    const auto *sgk_453 = buffer.data(sgk + 453);
    const auto *sgk_455 = buffer.data(sgk + 455);
    const auto *sgk_456 = buffer.data(sgk + 456);
    const auto *sgk_457 = buffer.data(sgk + 457);
    const auto *sgk_459 = buffer.data(sgk + 459);
    const auto *sgk_460 = buffer.data(sgk + 460);
    const auto *sgk_461 = buffer.data(sgk + 461);
    const auto *sgk_462 = buffer.data(sgk + 462);
    const auto *sgk_463 = buffer.data(sgk + 463);
    const auto *sgk_464 = buffer.data(sgk + 464);
    const auto *sgk_465 = buffer.data(sgk + 465);
    const auto *sgk_466 = buffer.data(sgk + 466);
    const auto *sgk_467 = buffer.data(sgk + 467);
    const auto *sgk_468 = buffer.data(sgk + 468);
    const auto *sgk_470 = buffer.data(sgk + 470);
    const auto *sgk_471 = buffer.data(sgk + 471);
    const auto *sgk_473 = buffer.data(sgk + 473);
    const auto *sgk_474 = buffer.data(sgk + 474);

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pc_x, pc_z, sgi0_301, sgi0_303, sgi0_304, \
                         sgi1_301, sgi1_303, sgi1_304, sgk_375, sgk_381, sgk_383, \
                         sgk_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_12 * sgi0_301[k]
                   - f_13 * sgi1_301[k]
                   + f_3 * pc_x[k] * sgk_381[k];

        t_472[k] = f_3 * pc_z[k] * sgk_375[k];

        t_473[k] = f_12 * sgi0_303[k]
                   - f_13 * sgi1_303[k]
                   + f_3 * pc_x[k] * sgk_383[k];

        t_474[k] = f_12 * sgi0_304[k]
                   - f_13 * sgi1_304[k]
                   + f_3 * pc_x[k] * sgk_384[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pc_x, pc_y, sfk_236, sgi0_305, sgi0_307, \
                         sgi1_305, sgi1_307, sgk_380, sgk_385, sgk_387, \
                         sgk_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_12 * sgi0_305[k]
                   - f_13 * sgi1_305[k]
                   + f_3 * pc_x[k] * sgk_385[k];

        t_476[k] = f_0 * sfk_236[k]
                   + f_3 * pc_y[k] * sgk_380[k];

        t_477[k] = f_12 * sgi0_307[k]
                   - f_13 * sgi1_307[k]
                   + f_3 * pc_x[k] * sgk_387[k];

        t_478[k] = f_3 * pc_x[k] * sgk_388[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, t_483, t_484, t_485, pc_x, sgk_389, \
                         sgk_390, sgk_391, sgk_392, sgk_393, sgk_394, \
                         sgk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_3 * pc_x[k] * sgk_389[k];

        t_480[k] = f_3 * pc_x[k] * sgk_390[k];

        t_481[k] = f_3 * pc_x[k] * sgk_391[k];

        t_482[k] = f_3 * pc_x[k] * sgk_392[k];

        t_483[k] = f_3 * pc_x[k] * sgk_393[k];

        t_484[k] = f_3 * pc_x[k] * sgk_394[k];

        t_485[k] = f_3 * pc_x[k] * sgk_395[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, pc_y, pc_z, sfk_244, sfk_246, sgi0_301, \
                         sgi0_303, sgi1_301, sgi1_303, sgk_388, \
                         sgk_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_0 * sfk_244[k]
                   + f_1 * sgi0_301[k]
                   - f_2 * sgi1_301[k]
                   + f_3 * pc_y[k] * sgk_388[k];

        t_487[k] = f_3 * pc_z[k] * sgk_388[k];

        t_488[k] = f_0 * sfk_246[k]
                   + f_4 * sgi0_303[k]
                   - f_5 * sgi1_303[k]
                   + f_3 * pc_y[k] * sgk_390[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_y, sfk_247, sfk_248, sfk_249, sgi0_304, \
                         sgi0_305, sgi0_306, sgi1_304, sgi1_305, sgi1_306, sgk_391, sgk_392, \
                         sgk_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_0 * sfk_247[k]
                   + f_6 * sgi0_304[k]
                   - f_7 * sgi1_304[k]
                   + f_3 * pc_y[k] * sgk_391[k];

        t_490[k] = f_0 * sfk_248[k]
                   + f_8 * sgi0_305[k]
                   - f_9 * sgi1_305[k]
                   + f_3 * pc_y[k] * sgk_392[k];

        t_491[k] = f_0 * sfk_249[k]
                   + f_10 * sgi0_306[k]
                   - f_11 * sgi1_306[k]
                   + f_3 * pc_y[k] * sgk_393[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pb_z, pc_y, pc_z, sfl0_270, sfk_250, \
                         sfk_251, sfl1_270, sgi0_307, sgi1_307, sgk_394, \
                         sgk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_0 * sfk_250[k]
                   + f_12 * sgi0_307[k]
                   - f_13 * sgi1_307[k]
                   + f_3 * pc_y[k] * sgk_394[k];

        t_493[k] = f_0 * sfk_251[k]
                   + f_3 * pc_y[k] * sgk_395[k];

        t_494[k] = f_1 * sgi0_307[k]
                   - f_2 * sgi1_307[k]
                   + f_3 * pc_z[k] * sgk_395[k];

        t_495[k] = pb_z[k] * sfl0_270[k]
                   - f_14 * pc_z[k] * sfl1_270[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pb_z, pc_y, pc_z, sfl0_273, sfk_216, \
                         sfk_252, sfk_254, sfl1_273, sgk_396, sgk_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_17 * sfk_252[k]
                   + f_3 * pc_y[k] * sgk_396[k];

        t_497[k] = f_15 * sfk_216[k]
                   + f_3 * pc_z[k] * sgk_396[k];

        t_498[k] = pb_z[k] * sfl0_273[k]
                   - f_14 * pc_z[k] * sfl1_273[k];

        t_499[k] = f_17 * sfk_254[k]
                   + f_3 * pc_y[k] * sgk_398[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pb_z, pc_x, pc_y, pc_z, sfl0_276, \
                         sfk_219, sfk_257, sfl1_276, sgi0_313, sgi1_313, sgk_399, \
                         sgk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_4 * sgi0_313[k]
                   - f_5 * sgi1_313[k]
                   + f_3 * pc_x[k] * sgk_401[k];

        t_501[k] = pb_z[k] * sfl0_276[k]
                   - f_14 * pc_z[k] * sfl1_276[k];

        t_502[k] = f_15 * sfk_219[k]
                   + f_3 * pc_z[k] * sgk_399[k];

        t_503[k] = f_17 * sfk_257[k]
                   + f_3 * pc_y[k] * sgk_401[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, pb_z, pc_x, pc_z, sfl0_280, sfk_222, sfl1_280, \
                         sgi0_317, sgi1_317, sgk_402, sgk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_6 * sgi0_317[k]
                   - f_7 * sgi1_317[k]
                   + f_3 * pc_x[k] * sgk_405[k];

        t_505[k] = pb_z[k] * sfl0_280[k]
                   - f_14 * pc_z[k] * sfl1_280[k];

        t_506[k] = f_15 * sfk_222[k]
                   + f_3 * pc_z[k] * sgk_402[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, pc_x, pc_y, sfk_261, sgi0_320, sgi0_322, \
                         sgi1_320, sgi1_322, sgk_405, sgk_408, \
                         sgk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_8 * sgi0_320[k]
                   - f_9 * sgi1_320[k]
                   + f_3 * pc_x[k] * sgk_408[k];

        t_508[k] = f_17 * sfk_261[k]
                   + f_3 * pc_y[k] * sgk_405[k];

        t_509[k] = f_8 * sgi0_322[k]
                   - f_9 * sgi1_322[k]
                   + f_3 * pc_x[k] * sgk_410[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, pb_z, pc_x, pc_z, sfl0_285, sfk_226, sfl1_285, \
                         sgi0_325, sgi1_325, sgk_406, sgk_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = pb_z[k] * sfl0_285[k]
                   - f_14 * pc_z[k] * sfl1_285[k];

        t_511[k] = f_15 * sfk_226[k]
                   + f_3 * pc_z[k] * sgk_406[k];

        t_512[k] = f_10 * sgi0_325[k]
                   - f_11 * sgi1_325[k]
                   + f_3 * pc_x[k] * sgk_413[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, pc_x, pc_y, sfk_266, sgi0_326, sgi0_328, \
                         sgi1_326, sgi1_328, sgk_410, sgk_414, \
                         sgk_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = f_10 * sgi0_326[k]
                   - f_11 * sgi1_326[k]
                   + f_3 * pc_x[k] * sgk_414[k];

        t_514[k] = f_17 * sfk_266[k]
                   + f_3 * pc_y[k] * sgk_410[k];

        t_515[k] = f_10 * sgi0_328[k]
                   - f_11 * sgi1_328[k]
                   + f_3 * pc_x[k] * sgk_416[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, pb_z, pc_x, pc_z, sfl0_291, sfk_231, sfl1_291, \
                         sgi0_331, sgi1_331, sgk_411, sgk_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = pb_z[k] * sfl0_291[k]
                   - f_14 * pc_z[k] * sfl1_291[k];

        t_517[k] = f_15 * sfk_231[k]
                   + f_3 * pc_z[k] * sgk_411[k];

        t_518[k] = f_12 * sgi0_331[k]
                   - f_13 * sgi1_331[k]
                   + f_3 * pc_x[k] * sgk_419[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_x, pc_y, sfk_272, sgi0_332, sgi0_333, \
                         sgi1_332, sgi1_333, sgk_416, sgk_420, \
                         sgk_421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_12 * sgi0_332[k]
                   - f_13 * sgi1_332[k]
                   + f_3 * pc_x[k] * sgk_420[k];

        t_520[k] = f_12 * sgi0_333[k]
                   - f_13 * sgi1_333[k]
                   + f_3 * pc_x[k] * sgk_421[k];

        t_521[k] = f_17 * sfk_272[k]
                   + f_3 * pc_y[k] * sgk_416[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, t_527, pc_x, sgi0_335, sgi1_335, \
                         sgk_423, sgk_424, sgk_425, sgk_426, sgk_427, \
                         sgk_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_12 * sgi0_335[k]
                   - f_13 * sgi1_335[k]
                   + f_3 * pc_x[k] * sgk_423[k];

        t_523[k] = f_3 * pc_x[k] * sgk_424[k];

        t_524[k] = f_3 * pc_x[k] * sgk_425[k];

        t_525[k] = f_3 * pc_x[k] * sgk_426[k];

        t_526[k] = f_3 * pc_x[k] * sgk_427[k];

        t_527[k] = f_3 * pc_x[k] * sgk_428[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, t_532, pb_z, pc_x, pc_z, sfl0_306, \
                         sfk_244, sfl1_306, sgk_424, sgk_429, sgk_430, \
                         sgk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = f_3 * pc_x[k] * sgk_429[k];

        t_529[k] = f_3 * pc_x[k] * sgk_430[k];

        t_530[k] = f_3 * pc_x[k] * sgk_431[k];

        t_531[k] = pb_z[k] * sfl0_306[k]
                   - f_14 * pc_z[k] * sfl1_306[k];

        t_532[k] = f_15 * sfk_244[k]
                   + f_3 * pc_z[k] * sgk_424[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, pb_z, pc_z, sfl0_308, sfl0_309, sfl0_310, \
                         sfk_245, sfk_246, sfk_247, sfl1_308, sfl1_309, \
                         sfl1_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = pb_z[k] * sfl0_308[k]
                   + f_16 * sfk_245[k]
                   - f_14 * pc_z[k] * sfl1_308[k];

        t_534[k] = pb_z[k] * sfl0_309[k]
                   + f_17 * sfk_246[k]
                   - f_14 * pc_z[k] * sfl1_309[k];

        t_535[k] = pb_z[k] * sfl0_310[k]
                   + f_0 * sfk_247[k]
                   - f_14 * pc_z[k] * sfl1_310[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, pb_z, pc_y, pc_z, sfl0_311, sfl0_312, sfk_248, \
                         sfk_249, sfk_287, sfl1_311, sfl1_312, \
                         sgk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = pb_z[k] * sfl0_311[k]
                   + f_18 * sfk_248[k]
                   - f_14 * pc_z[k] * sfl1_311[k];

        t_537[k] = pb_z[k] * sfl0_312[k]
                   + f_19 * sfk_249[k]
                   - f_14 * pc_z[k] * sfl1_312[k];

        t_538[k] = f_17 * sfk_287[k]
                   + f_3 * pc_y[k] * sgk_431[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, t_542, pc_x, pc_y, pc_z, sfk_251, sfk_252, \
                         sfk_288, sgi0_335, sgi0_336, sgi1_335, sgi1_336, sgk_431, \
                         sgk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_15 * sfk_251[k]
                   + f_1 * sgi0_335[k]
                   - f_2 * sgi1_335[k]
                   + f_3 * pc_z[k] * sgk_431[k];

        t_540[k] = f_1 * sgi0_336[k]
                   - f_2 * sgi1_336[k]
                   + f_3 * pc_x[k] * sgk_432[k];

        t_541[k] = f_16 * sfk_288[k]
                   + f_3 * pc_y[k] * sgk_432[k];

        t_542[k] = f_16 * sfk_252[k]
                   + f_3 * pc_z[k] * sgk_432[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, pc_x, pc_y, sfk_290, sgi0_339, sgi0_341, \
                         sgi1_339, sgi1_341, sgk_434, sgk_435, \
                         sgk_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_4 * sgi0_339[k]
                   - f_5 * sgi1_339[k]
                   + f_3 * pc_x[k] * sgk_435[k];

        t_544[k] = f_16 * sfk_290[k]
                   + f_3 * pc_y[k] * sgk_434[k];

        t_545[k] = f_4 * sgi0_341[k]
                   - f_5 * sgi1_341[k]
                   + f_3 * pc_x[k] * sgk_437[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, pc_x, pc_y, pc_z, sfk_255, sfk_293, sgi0_342, \
                         sgi1_342, sgk_435, sgk_437, sgk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = f_6 * sgi0_342[k]
                   - f_7 * sgi1_342[k]
                   + f_3 * pc_x[k] * sgk_438[k];

        t_547[k] = f_16 * sfk_255[k]
                   + f_3 * pc_z[k] * sgk_435[k];

        t_548[k] = f_16 * sfk_293[k]
                   + f_3 * pc_y[k] * sgk_437[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, pc_x, pc_z, sfk_258, sgi0_345, sgi0_346, \
                         sgi1_345, sgi1_346, sgk_438, sgk_441, \
                         sgk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_6 * sgi0_345[k]
                   - f_7 * sgi1_345[k]
                   + f_3 * pc_x[k] * sgk_441[k];

        t_550[k] = f_8 * sgi0_346[k]
                   - f_9 * sgi1_346[k]
                   + f_3 * pc_x[k] * sgk_442[k];

        t_551[k] = f_16 * sfk_258[k]
                   + f_3 * pc_z[k] * sgk_438[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pc_x, pc_y, sfk_297, sgi0_348, sgi0_350, \
                         sgi1_348, sgi1_350, sgk_441, sgk_444, \
                         sgk_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_8 * sgi0_348[k]
                   - f_9 * sgi1_348[k]
                   + f_3 * pc_x[k] * sgk_444[k];

        t_553[k] = f_16 * sfk_297[k]
                   + f_3 * pc_y[k] * sgk_441[k];

        t_554[k] = f_8 * sgi0_350[k]
                   - f_9 * sgi1_350[k]
                   + f_3 * pc_x[k] * sgk_446[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, pc_x, pc_z, sfk_262, sgi0_351, sgi0_353, \
                         sgi1_351, sgi1_353, sgk_442, sgk_447, \
                         sgk_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_10 * sgi0_351[k]
                   - f_11 * sgi1_351[k]
                   + f_3 * pc_x[k] * sgk_447[k];

        t_556[k] = f_16 * sfk_262[k]
                   + f_3 * pc_z[k] * sgk_442[k];

        t_557[k] = f_10 * sgi0_353[k]
                   - f_11 * sgi1_353[k]
                   + f_3 * pc_x[k] * sgk_449[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, pc_x, pc_y, sfk_302, sgi0_354, sgi0_356, \
                         sgi1_354, sgi1_356, sgk_446, sgk_450, \
                         sgk_452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_10 * sgi0_354[k]
                   - f_11 * sgi1_354[k]
                   + f_3 * pc_x[k] * sgk_450[k];

        t_559[k] = f_16 * sfk_302[k]
                   + f_3 * pc_y[k] * sgk_446[k];

        t_560[k] = f_10 * sgi0_356[k]
                   - f_11 * sgi1_356[k]
                   + f_3 * pc_x[k] * sgk_452[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pc_x, pc_z, sfk_267, sgi0_357, sgi0_359, \
                         sgi1_357, sgi1_359, sgk_447, sgk_453, \
                         sgk_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_12 * sgi0_357[k]
                   - f_13 * sgi1_357[k]
                   + f_3 * pc_x[k] * sgk_453[k];

        t_562[k] = f_16 * sfk_267[k]
                   + f_3 * pc_z[k] * sgk_447[k];

        t_563[k] = f_12 * sgi0_359[k]
                   - f_13 * sgi1_359[k]
                   + f_3 * pc_x[k] * sgk_455[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pc_x, pc_y, sfk_308, sgi0_360, sgi0_361, \
                         sgi1_360, sgi1_361, sgk_452, sgk_456, \
                         sgk_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_12 * sgi0_360[k]
                   - f_13 * sgi1_360[k]
                   + f_3 * pc_x[k] * sgk_456[k];

        t_565[k] = f_12 * sgi0_361[k]
                   - f_13 * sgi1_361[k]
                   + f_3 * pc_x[k] * sgk_457[k];

        t_566[k] = f_16 * sfk_308[k]
                   + f_3 * pc_y[k] * sgk_452[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, t_572, pc_x, sgi0_363, sgi1_363, \
                         sgk_459, sgk_460, sgk_461, sgk_462, sgk_463, \
                         sgk_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_12 * sgi0_363[k]
                   - f_13 * sgi1_363[k]
                   + f_3 * pc_x[k] * sgk_459[k];

        t_568[k] = f_3 * pc_x[k] * sgk_460[k];

        t_569[k] = f_3 * pc_x[k] * sgk_461[k];

        t_570[k] = f_3 * pc_x[k] * sgk_462[k];

        t_571[k] = f_3 * pc_x[k] * sgk_463[k];

        t_572[k] = f_3 * pc_x[k] * sgk_464[k];
    }

#pragma omp simd aligned(t_573, t_574, t_575, t_576, t_577, pc_x, pc_y, pc_z, sfk_280, \
                         sfk_316, sgi0_357, sgi1_357, sgk_460, sgk_465, sgk_466, \
                         sgk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_573[k] = f_3 * pc_x[k] * sgk_465[k];

        t_574[k] = f_3 * pc_x[k] * sgk_466[k];

        t_575[k] = f_3 * pc_x[k] * sgk_467[k];

        t_576[k] = f_16 * sfk_316[k]
                   + f_1 * sgi0_357[k]
                   - f_2 * sgi1_357[k]
                   + f_3 * pc_y[k] * sgk_460[k];

        t_577[k] = f_16 * sfk_280[k]
                   + f_3 * pc_z[k] * sgk_460[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pc_y, sfk_318, sfk_319, sfk_320, sgi0_359, \
                         sgi0_360, sgi0_361, sgi1_359, sgi1_360, sgi1_361, sgk_462, sgk_463, \
                         sgk_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_16 * sfk_318[k]
                   + f_4 * sgi0_359[k]
                   - f_5 * sgi1_359[k]
                   + f_3 * pc_y[k] * sgk_462[k];

        t_579[k] = f_16 * sfk_319[k]
                   + f_6 * sgi0_360[k]
                   - f_7 * sgi1_360[k]
                   + f_3 * pc_y[k] * sgk_463[k];

        t_580[k] = f_16 * sfk_320[k]
                   + f_8 * sgi0_361[k]
                   - f_9 * sgi1_361[k]
                   + f_3 * pc_y[k] * sgk_464[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, sfk_321, sfk_322, sfk_323, sgi0_362, \
                         sgi0_363, sgi1_362, sgi1_363, sgk_465, sgk_466, \
                         sgk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_16 * sfk_321[k]
                   + f_10 * sgi0_362[k]
                   - f_11 * sgi1_362[k]
                   + f_3 * pc_y[k] * sgk_465[k];

        t_582[k] = f_16 * sfk_322[k]
                   + f_12 * sgi0_363[k]
                   - f_13 * sgi1_363[k]
                   + f_3 * pc_y[k] * sgk_466[k];

        t_583[k] = f_16 * sfk_323[k]
                   + f_3 * pc_y[k] * sgk_467[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pb_y, pc_y, pc_z, sfl0_405, sfk_287, \
                         sfk_288, sfk_324, sfl1_405, sgi0_363, sgi1_363, sgk_467, \
                         sgk_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * sfk_287[k]
                   + f_1 * sgi0_363[k]
                   - f_2 * sgi1_363[k]
                   + f_3 * pc_z[k] * sgk_467[k];

        t_585[k] = pb_y[k] * sfl0_405[k]
                   - f_14 * pc_y[k] * sfl1_405[k];

        t_586[k] = f_15 * sfk_324[k]
                   + f_3 * pc_y[k] * sgk_468[k];

        t_587[k] = f_17 * sfk_288[k]
                   + f_3 * pc_z[k] * sgk_468[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, pb_y, pc_x, pc_y, sfl0_410, sfk_326, sfl1_410, \
                         sgi0_367, sgi1_367, sgk_470, sgk_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_4 * sgi0_367[k]
                   - f_5 * sgi1_367[k]
                   + f_3 * pc_x[k] * sgk_471[k];

        t_589[k] = f_15 * sfk_326[k]
                   + f_3 * pc_y[k] * sgk_470[k];

        t_590[k] = pb_y[k] * sfl0_410[k]
                   - f_14 * pc_y[k] * sfl1_410[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, pc_x, pc_y, pc_z, sfk_291, sfk_329, sgi0_370, \
                         sgi1_370, sgk_471, sgk_473, sgk_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_6 * sgi0_370[k]
                   - f_7 * sgi1_370[k]
                   + f_3 * pc_x[k] * sgk_474[k];

        t_592[k] = f_17 * sfk_291[k]
                   + f_3 * pc_z[k] * sgk_471[k];

        t_593[k] = f_15 * sfk_329[k]
                   + f_3 * pc_y[k] * sgk_473[k];
    }
}

static auto
compute_prim_sgl_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfl0,
                                                          const size_t sfk, const size_t sfl1,
                                                          const size_t sgi0, const size_t sgi1,
                                                          const size_t sgk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 4.0 / q;

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

    const auto *sfl0_414 = buffer.data(sfl0 + 414);
    const auto *sfl0_419 = buffer.data(sfl0 + 419);
    const auto *sfl0_425 = buffer.data(sfl0 + 425);
    const auto *sfl0_432 = buffer.data(sfl0 + 432);
    const auto *sfl0_441 = buffer.data(sfl0 + 441);
    const auto *sfl0_443 = buffer.data(sfl0 + 443);
    const auto *sfl0_444 = buffer.data(sfl0 + 444);
    const auto *sfl0_445 = buffer.data(sfl0 + 445);
    const auto *sfl0_446 = buffer.data(sfl0 + 446);
    const auto *sfl0_447 = buffer.data(sfl0 + 447);
    const auto *sfl0_449 = buffer.data(sfl0 + 449);

    const auto *sfk_294 = buffer.data(sfk + 294);
    const auto *sfk_298 = buffer.data(sfk + 298);
    const auto *sfk_303 = buffer.data(sfk + 303);
    const auto *sfk_316 = buffer.data(sfk + 316);
    const auto *sfk_324 = buffer.data(sfk + 324);
    const auto *sfk_327 = buffer.data(sfk + 327);
    const auto *sfk_330 = buffer.data(sfk + 330);
    const auto *sfk_333 = buffer.data(sfk + 333);
    const auto *sfk_334 = buffer.data(sfk + 334);
    const auto *sfk_338 = buffer.data(sfk + 338);
    const auto *sfk_339 = buffer.data(sfk + 339);
    const auto *sfk_344 = buffer.data(sfk + 344);
    const auto *sfk_352 = buffer.data(sfk + 352);
    const auto *sfk_354 = buffer.data(sfk + 354);
    const auto *sfk_355 = buffer.data(sfk + 355);
    const auto *sfk_356 = buffer.data(sfk + 356);
    const auto *sfk_357 = buffer.data(sfk + 357);
    const auto *sfk_358 = buffer.data(sfk + 358);
    const auto *sfk_359 = buffer.data(sfk + 359);

    const auto *sfl1_414 = buffer.data(sfl1 + 414);
    const auto *sfl1_419 = buffer.data(sfl1 + 419);
    const auto *sfl1_425 = buffer.data(sfl1 + 425);
    const auto *sfl1_432 = buffer.data(sfl1 + 432);
    const auto *sfl1_441 = buffer.data(sfl1 + 441);
    const auto *sfl1_443 = buffer.data(sfl1 + 443);
    const auto *sfl1_444 = buffer.data(sfl1 + 444);
    const auto *sfl1_445 = buffer.data(sfl1 + 445);
    const auto *sfl1_446 = buffer.data(sfl1 + 446);
    const auto *sfl1_447 = buffer.data(sfl1 + 447);
    const auto *sfl1_449 = buffer.data(sfl1 + 449);

    const auto *sgi0_374 = buffer.data(sgi0 + 374);
    const auto *sgi0_376 = buffer.data(sgi0 + 376);
    const auto *sgi0_379 = buffer.data(sgi0 + 379);
    const auto *sgi0_381 = buffer.data(sgi0 + 381);
    const auto *sgi0_382 = buffer.data(sgi0 + 382);
    const auto *sgi0_385 = buffer.data(sgi0 + 385);
    const auto *sgi0_387 = buffer.data(sgi0 + 387);
    const auto *sgi0_388 = buffer.data(sgi0 + 388);
    const auto *sgi0_389 = buffer.data(sgi0 + 389);
    const auto *sgi0_392 = buffer.data(sgi0 + 392);
    const auto *sgi0_395 = buffer.data(sgi0 + 395);
    const auto *sgi0_397 = buffer.data(sgi0 + 397);
    const auto *sgi0_398 = buffer.data(sgi0 + 398);
    const auto *sgi0_401 = buffer.data(sgi0 + 401);
    const auto *sgi0_402 = buffer.data(sgi0 + 402);
    const auto *sgi0_404 = buffer.data(sgi0 + 404);
    const auto *sgi0_406 = buffer.data(sgi0 + 406);
    const auto *sgi0_407 = buffer.data(sgi0 + 407);
    const auto *sgi0_409 = buffer.data(sgi0 + 409);
    const auto *sgi0_410 = buffer.data(sgi0 + 410);
    const auto *sgi0_412 = buffer.data(sgi0 + 412);
    const auto *sgi0_413 = buffer.data(sgi0 + 413);
    const auto *sgi0_415 = buffer.data(sgi0 + 415);
    const auto *sgi0_416 = buffer.data(sgi0 + 416);
    const auto *sgi0_417 = buffer.data(sgi0 + 417);
    const auto *sgi0_418 = buffer.data(sgi0 + 418);
    const auto *sgi0_419 = buffer.data(sgi0 + 419);

    const auto *sgi1_374 = buffer.data(sgi1 + 374);
    const auto *sgi1_376 = buffer.data(sgi1 + 376);
    const auto *sgi1_379 = buffer.data(sgi1 + 379);
    const auto *sgi1_381 = buffer.data(sgi1 + 381);
    const auto *sgi1_382 = buffer.data(sgi1 + 382);
    const auto *sgi1_385 = buffer.data(sgi1 + 385);
    const auto *sgi1_387 = buffer.data(sgi1 + 387);
    const auto *sgi1_388 = buffer.data(sgi1 + 388);
    const auto *sgi1_389 = buffer.data(sgi1 + 389);
    const auto *sgi1_392 = buffer.data(sgi1 + 392);
    const auto *sgi1_395 = buffer.data(sgi1 + 395);
    const auto *sgi1_397 = buffer.data(sgi1 + 397);
    const auto *sgi1_398 = buffer.data(sgi1 + 398);
    const auto *sgi1_401 = buffer.data(sgi1 + 401);
    const auto *sgi1_402 = buffer.data(sgi1 + 402);
    const auto *sgi1_404 = buffer.data(sgi1 + 404);
    const auto *sgi1_406 = buffer.data(sgi1 + 406);
    const auto *sgi1_407 = buffer.data(sgi1 + 407);
    const auto *sgi1_409 = buffer.data(sgi1 + 409);
    const auto *sgi1_410 = buffer.data(sgi1 + 410);
    const auto *sgi1_412 = buffer.data(sgi1 + 412);
    const auto *sgi1_413 = buffer.data(sgi1 + 413);
    const auto *sgi1_415 = buffer.data(sgi1 + 415);
    const auto *sgi1_416 = buffer.data(sgi1 + 416);
    const auto *sgi1_417 = buffer.data(sgi1 + 417);
    const auto *sgi1_418 = buffer.data(sgi1 + 418);
    const auto *sgi1_419 = buffer.data(sgi1 + 419);

    const auto *sgk_474 = buffer.data(sgk + 474);
    const auto *sgk_477 = buffer.data(sgk + 477);
    const auto *sgk_478 = buffer.data(sgk + 478);
    const auto *sgk_480 = buffer.data(sgk + 480);
    const auto *sgk_482 = buffer.data(sgk + 482);
    const auto *sgk_483 = buffer.data(sgk + 483);
    const auto *sgk_485 = buffer.data(sgk + 485);
    const auto *sgk_486 = buffer.data(sgk + 486);
    const auto *sgk_488 = buffer.data(sgk + 488);
    const auto *sgk_489 = buffer.data(sgk + 489);
    const auto *sgk_491 = buffer.data(sgk + 491);
    const auto *sgk_492 = buffer.data(sgk + 492);
    const auto *sgk_493 = buffer.data(sgk + 493);
    const auto *sgk_496 = buffer.data(sgk + 496);
    const auto *sgk_497 = buffer.data(sgk + 497);
    const auto *sgk_498 = buffer.data(sgk + 498);
    const auto *sgk_499 = buffer.data(sgk + 499);
    const auto *sgk_500 = buffer.data(sgk + 500);
    const auto *sgk_501 = buffer.data(sgk + 501);
    const auto *sgk_502 = buffer.data(sgk + 502);
    const auto *sgk_503 = buffer.data(sgk + 503);
    const auto *sgk_504 = buffer.data(sgk + 504);
    const auto *sgk_506 = buffer.data(sgk + 506);
    const auto *sgk_507 = buffer.data(sgk + 507);
    const auto *sgk_509 = buffer.data(sgk + 509);
    const auto *sgk_510 = buffer.data(sgk + 510);
    const auto *sgk_513 = buffer.data(sgk + 513);
    const auto *sgk_514 = buffer.data(sgk + 514);
    const auto *sgk_516 = buffer.data(sgk + 516);
    const auto *sgk_518 = buffer.data(sgk + 518);
    const auto *sgk_519 = buffer.data(sgk + 519);
    const auto *sgk_521 = buffer.data(sgk + 521);
    const auto *sgk_522 = buffer.data(sgk + 522);
    const auto *sgk_524 = buffer.data(sgk + 524);
    const auto *sgk_525 = buffer.data(sgk + 525);
    const auto *sgk_527 = buffer.data(sgk + 527);
    const auto *sgk_528 = buffer.data(sgk + 528);
    const auto *sgk_529 = buffer.data(sgk + 529);
    const auto *sgk_531 = buffer.data(sgk + 531);
    const auto *sgk_532 = buffer.data(sgk + 532);
    const auto *sgk_533 = buffer.data(sgk + 533);
    const auto *sgk_534 = buffer.data(sgk + 534);
    const auto *sgk_535 = buffer.data(sgk + 535);
    const auto *sgk_536 = buffer.data(sgk + 536);
    const auto *sgk_537 = buffer.data(sgk + 537);
    const auto *sgk_538 = buffer.data(sgk + 538);
    const auto *sgk_539 = buffer.data(sgk + 539);

#pragma omp simd aligned(t_594, t_595, t_596, pb_y, pc_x, pc_y, pc_z, sfl0_414, sfk_294, \
                         sfl1_414, sgi0_374, sgi1_374, sgk_474, \
                         sgk_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = pb_y[k] * sfl0_414[k]
                   - f_14 * pc_y[k] * sfl1_414[k];

        t_595[k] = f_8 * sgi0_374[k]
                   - f_9 * sgi1_374[k]
                   + f_3 * pc_x[k] * sgk_478[k];

        t_596[k] = f_17 * sfk_294[k]
                   + f_3 * pc_z[k] * sgk_474[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, pb_y, pc_x, pc_y, sfl0_419, sfk_333, sfl1_419, \
                         sgi0_376, sgi1_376, sgk_477, sgk_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = f_8 * sgi0_376[k]
                   - f_9 * sgi1_376[k]
                   + f_3 * pc_x[k] * sgk_480[k];

        t_598[k] = f_15 * sfk_333[k]
                   + f_3 * pc_y[k] * sgk_477[k];

        t_599[k] = pb_y[k] * sfl0_419[k]
                   - f_14 * pc_y[k] * sfl1_419[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pc_x, pc_z, sfk_298, sgi0_379, sgi0_381, \
                         sgi1_379, sgi1_381, sgk_478, sgk_483, \
                         sgk_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_10 * sgi0_379[k]
                   - f_11 * sgi1_379[k]
                   + f_3 * pc_x[k] * sgk_483[k];

        t_601[k] = f_17 * sfk_298[k]
                   + f_3 * pc_z[k] * sgk_478[k];

        t_602[k] = f_10 * sgi0_381[k]
                   - f_11 * sgi1_381[k]
                   + f_3 * pc_x[k] * sgk_485[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, pb_y, pc_x, pc_y, sfl0_425, sfk_338, sfl1_425, \
                         sgi0_382, sgi1_382, sgk_482, sgk_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = f_10 * sgi0_382[k]
                   - f_11 * sgi1_382[k]
                   + f_3 * pc_x[k] * sgk_486[k];

        t_604[k] = f_15 * sfk_338[k]
                   + f_3 * pc_y[k] * sgk_482[k];

        t_605[k] = pb_y[k] * sfl0_425[k]
                   - f_14 * pc_y[k] * sfl1_425[k];
    }

#pragma omp simd aligned(t_606, t_607, t_608, pc_x, pc_z, sfk_303, sgi0_385, sgi0_387, \
                         sgi1_385, sgi1_387, sgk_483, sgk_489, \
                         sgk_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_606[k] = f_12 * sgi0_385[k]
                   - f_13 * sgi1_385[k]
                   + f_3 * pc_x[k] * sgk_489[k];

        t_607[k] = f_17 * sfk_303[k]
                   + f_3 * pc_z[k] * sgk_483[k];

        t_608[k] = f_12 * sgi0_387[k]
                   - f_13 * sgi1_387[k]
                   + f_3 * pc_x[k] * sgk_491[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pc_x, pc_y, sfk_344, sgi0_388, sgi0_389, \
                         sgi1_388, sgi1_389, sgk_488, sgk_492, \
                         sgk_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_12 * sgi0_388[k]
                   - f_13 * sgi1_388[k]
                   + f_3 * pc_x[k] * sgk_492[k];

        t_610[k] = f_12 * sgi0_389[k]
                   - f_13 * sgi1_389[k]
                   + f_3 * pc_x[k] * sgk_493[k];

        t_611[k] = f_15 * sfk_344[k]
                   + f_3 * pc_y[k] * sgk_488[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, t_616, t_617, pb_y, pc_x, pc_y, sfl0_432, \
                         sfl1_432, sgk_496, sgk_497, sgk_498, sgk_499, \
                         sgk_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = pb_y[k] * sfl0_432[k]
                   - f_14 * pc_y[k] * sfl1_432[k];

        t_613[k] = f_3 * pc_x[k] * sgk_496[k];

        t_614[k] = f_3 * pc_x[k] * sgk_497[k];

        t_615[k] = f_3 * pc_x[k] * sgk_498[k];

        t_616[k] = f_3 * pc_x[k] * sgk_499[k];

        t_617[k] = f_3 * pc_x[k] * sgk_500[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, pb_y, pc_x, pc_y, sfl0_441, sfk_352, \
                         sfl1_441, sgk_501, sgk_502, sgk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_3 * pc_x[k] * sgk_501[k];

        t_619[k] = f_3 * pc_x[k] * sgk_502[k];

        t_620[k] = f_3 * pc_x[k] * sgk_503[k];

        t_621[k] = pb_y[k] * sfl0_441[k]
                   + f_20 * sfk_352[k]
                   - f_14 * pc_y[k] * sfl1_441[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, pb_y, pc_y, pc_z, sfl0_443, sfl0_444, sfk_316, \
                         sfk_354, sfk_355, sfl1_443, sfl1_444, \
                         sgk_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_17 * sfk_316[k]
                   + f_3 * pc_z[k] * sgk_496[k];

        t_623[k] = pb_y[k] * sfl0_443[k]
                   + f_19 * sfk_354[k]
                   - f_14 * pc_y[k] * sfl1_443[k];

        t_624[k] = pb_y[k] * sfl0_444[k]
                   + f_18 * sfk_355[k]
                   - f_14 * pc_y[k] * sfl1_444[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, pb_y, pc_y, sfl0_445, sfl0_446, sfl0_447, \
                         sfk_356, sfk_357, sfk_358, sfl1_445, sfl1_446, \
                         sfl1_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = pb_y[k] * sfl0_445[k]
                   + f_0 * sfk_356[k]
                   - f_14 * pc_y[k] * sfl1_445[k];

        t_626[k] = pb_y[k] * sfl0_446[k]
                   + f_17 * sfk_357[k]
                   - f_14 * pc_y[k] * sfl1_446[k];

        t_627[k] = pb_y[k] * sfl0_447[k]
                   + f_16 * sfk_358[k]
                   - f_14 * pc_y[k] * sfl1_447[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, t_631, pb_y, pc_x, pc_y, sfl0_449, sfk_359, \
                         sfl1_449, sgi0_392, sgi1_392, sgk_503, \
                         sgk_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_15 * sfk_359[k]
                   + f_3 * pc_y[k] * sgk_503[k];

        t_629[k] = pb_y[k] * sfl0_449[k]
                   - f_14 * pc_y[k] * sfl1_449[k];

        t_630[k] = f_1 * sgi0_392[k]
                   - f_2 * sgi1_392[k]
                   + f_3 * pc_x[k] * sgk_504[k];

        t_631[k] = f_3 * pc_y[k] * sgk_504[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, pc_x, pc_y, pc_z, sfk_324, sgi0_395, \
                         sgi0_397, sgi1_395, sgi1_397, sgk_504, sgk_506, sgk_507, \
                         sgk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_0 * sfk_324[k]
                   + f_3 * pc_z[k] * sgk_504[k];

        t_633[k] = f_4 * sgi0_395[k]
                   - f_5 * sgi1_395[k]
                   + f_3 * pc_x[k] * sgk_507[k];

        t_634[k] = f_3 * pc_y[k] * sgk_506[k];

        t_635[k] = f_4 * sgi0_397[k]
                   - f_5 * sgi1_397[k]
                   + f_3 * pc_x[k] * sgk_509[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pc_x, pc_y, pc_z, sfk_327, sgi0_398, \
                         sgi0_401, sgi1_398, sgi1_401, sgk_507, sgk_509, sgk_510, \
                         sgk_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_6 * sgi0_398[k]
                   - f_7 * sgi1_398[k]
                   + f_3 * pc_x[k] * sgk_510[k];

        t_637[k] = f_0 * sfk_327[k]
                   + f_3 * pc_z[k] * sgk_507[k];

        t_638[k] = f_3 * pc_y[k] * sgk_509[k];

        t_639[k] = f_6 * sgi0_401[k]
                   - f_7 * sgi1_401[k]
                   + f_3 * pc_x[k] * sgk_513[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, pc_x, pc_y, pc_z, sfk_330, sgi0_402, \
                         sgi0_404, sgi1_402, sgi1_404, sgk_510, sgk_513, sgk_514, \
                         sgk_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_8 * sgi0_402[k]
                   - f_9 * sgi1_402[k]
                   + f_3 * pc_x[k] * sgk_514[k];

        t_641[k] = f_0 * sfk_330[k]
                   + f_3 * pc_z[k] * sgk_510[k];

        t_642[k] = f_8 * sgi0_404[k]
                   - f_9 * sgi1_404[k]
                   + f_3 * pc_x[k] * sgk_516[k];

        t_643[k] = f_3 * pc_y[k] * sgk_513[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, pc_x, pc_z, sfk_334, sgi0_406, sgi0_407, \
                         sgi1_406, sgi1_407, sgk_514, sgk_518, \
                         sgk_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_8 * sgi0_406[k]
                   - f_9 * sgi1_406[k]
                   + f_3 * pc_x[k] * sgk_518[k];

        t_645[k] = f_10 * sgi0_407[k]
                   - f_11 * sgi1_407[k]
                   + f_3 * pc_x[k] * sgk_519[k];

        t_646[k] = f_0 * sfk_334[k]
                   + f_3 * pc_z[k] * sgk_514[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, pc_x, pc_y, sgi0_409, sgi0_410, sgi0_412, \
                         sgi1_409, sgi1_410, sgi1_412, sgk_518, sgk_521, sgk_522, \
                         sgk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_10 * sgi0_409[k]
                   - f_11 * sgi1_409[k]
                   + f_3 * pc_x[k] * sgk_521[k];

        t_648[k] = f_10 * sgi0_410[k]
                   - f_11 * sgi1_410[k]
                   + f_3 * pc_x[k] * sgk_522[k];

        t_649[k] = f_3 * pc_y[k] * sgk_518[k];

        t_650[k] = f_10 * sgi0_412[k]
                   - f_11 * sgi1_412[k]
                   + f_3 * pc_x[k] * sgk_524[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pc_x, pc_z, sfk_339, sgi0_413, sgi0_415, \
                         sgi1_413, sgi1_415, sgk_519, sgk_525, \
                         sgk_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_12 * sgi0_413[k]
                   - f_13 * sgi1_413[k]
                   + f_3 * pc_x[k] * sgk_525[k];

        t_652[k] = f_0 * sfk_339[k]
                   + f_3 * pc_z[k] * sgk_519[k];

        t_653[k] = f_12 * sgi0_415[k]
                   - f_13 * sgi1_415[k]
                   + f_3 * pc_x[k] * sgk_527[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, t_657, pc_x, pc_y, sgi0_416, sgi0_417, sgi0_419, \
                         sgi1_416, sgi1_417, sgi1_419, sgk_524, sgk_528, sgk_529, \
                         sgk_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_12 * sgi0_416[k]
                   - f_13 * sgi1_416[k]
                   + f_3 * pc_x[k] * sgk_528[k];

        t_655[k] = f_12 * sgi0_417[k]
                   - f_13 * sgi1_417[k]
                   + f_3 * pc_x[k] * sgk_529[k];

        t_656[k] = f_3 * pc_y[k] * sgk_524[k];

        t_657[k] = f_12 * sgi0_419[k]
                   - f_13 * sgi1_419[k]
                   + f_3 * pc_x[k] * sgk_531[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, t_662, t_663, t_664, pc_x, sgk_532, \
                         sgk_533, sgk_534, sgk_535, sgk_536, sgk_537, \
                         sgk_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = f_3 * pc_x[k] * sgk_532[k];

        t_659[k] = f_3 * pc_x[k] * sgk_533[k];

        t_660[k] = f_3 * pc_x[k] * sgk_534[k];

        t_661[k] = f_3 * pc_x[k] * sgk_535[k];

        t_662[k] = f_3 * pc_x[k] * sgk_536[k];

        t_663[k] = f_3 * pc_x[k] * sgk_537[k];

        t_664[k] = f_3 * pc_x[k] * sgk_538[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, pc_x, pc_y, pc_z, sfk_352, sgi0_413, \
                         sgi0_415, sgi1_413, sgi1_415, sgk_532, sgk_534, \
                         sgk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_3 * pc_x[k] * sgk_539[k];

        t_666[k] = f_1 * sgi0_413[k]
                   - f_2 * sgi1_413[k]
                   + f_3 * pc_y[k] * sgk_532[k];

        t_667[k] = f_0 * sfk_352[k]
                   + f_3 * pc_z[k] * sgk_532[k];

        t_668[k] = f_4 * sgi0_415[k]
                   - f_5 * sgi1_415[k]
                   + f_3 * pc_y[k] * sgk_534[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, pc_y, sgi0_416, sgi0_417, sgi0_418, sgi1_416, \
                         sgi1_417, sgi1_418, sgk_535, sgk_536, \
                         sgk_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_6 * sgi0_416[k]
                   - f_7 * sgi1_416[k]
                   + f_3 * pc_y[k] * sgk_535[k];

        t_670[k] = f_8 * sgi0_417[k]
                   - f_9 * sgi1_417[k]
                   + f_3 * pc_y[k] * sgk_536[k];

        t_671[k] = f_10 * sgi0_418[k]
                   - f_11 * sgi1_418[k]
                   + f_3 * pc_y[k] * sgk_537[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, pc_y, pc_z, sfk_359, sgi0_419, sgi1_419, \
                         sgk_538, sgk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_12 * sgi0_419[k]
                   - f_13 * sgi1_419[k]
                   + f_3 * pc_y[k] * sgk_538[k];

        t_673[k] = f_3 * pc_y[k] * sgk_539[k];

        t_674[k] = f_0 * sfk_359[k]
                   + f_1 * sgi0_419[k]
                   - f_2 * sgi1_419[k]
                   + f_3 * pc_z[k] * sgk_539[k];
    }
}

auto
compute_prim_sgl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sfl0, const size_t sfk,
                                                   const size_t sfl1, const size_t sgi0,
                                                   const size_t sgi1, const size_t sgk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sgl_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sfl0, sfk,
                                                              sfl1, sgi0, sgi1, sgk, ncols,
                                                              gamma, p, q);

    compute_prim_sgl_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sfl0, sfk,
                                                              sfl1, sgi0, sgi1, sgk, ncols,
                                                              gamma, p, q);

    compute_prim_sgl_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sfl0, sfk,
                                                              sfl1, sgi0, sgi1, sgk, ncols,
                                                              gamma, p, q);

    compute_prim_sgl_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sfl0, sfk,
                                                              sfl1, sgi0, sgi1, sgk, ncols,
                                                              gamma, p, q);

    compute_prim_sgl_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sfl0, sfk,
                                                              sfl1, sgi0, sgi1, sgk, ncols,
                                                              gamma, p, q);

    compute_prim_sgl_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, sfl0, sfk,
                                                              sfl1, sgi0, sgi1, sgk, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
