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


#include "SimdThreeCenterElectronRepulsionVrrRecSHL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_shl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgl0,
                                                          const size_t sgk, const size_t sgl1,
                                                          const size_t shi0, const size_t shi1,
                                                          const size_t shk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_18 = 2.0 / q;
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

    const auto *sgl0_0 = buffer.data(sgl0 + 0);
    const auto *sgl0_3 = buffer.data(sgl0 + 3);
    const auto *sgl0_5 = buffer.data(sgl0 + 5);
    const auto *sgl0_6 = buffer.data(sgl0 + 6);
    const auto *sgl0_9 = buffer.data(sgl0 + 9);
    const auto *sgl0_10 = buffer.data(sgl0 + 10);
    const auto *sgl0_12 = buffer.data(sgl0 + 12);
    const auto *sgl0_14 = buffer.data(sgl0 + 14);
    const auto *sgl0_15 = buffer.data(sgl0 + 15);
    const auto *sgl0_17 = buffer.data(sgl0 + 17);
    const auto *sgl0_18 = buffer.data(sgl0 + 18);
    const auto *sgl0_20 = buffer.data(sgl0 + 20);
    const auto *sgl0_21 = buffer.data(sgl0 + 21);
    const auto *sgl0_23 = buffer.data(sgl0 + 23);
    const auto *sgl0_24 = buffer.data(sgl0 + 24);
    const auto *sgl0_25 = buffer.data(sgl0 + 25);
    const auto *sgl0_27 = buffer.data(sgl0 + 27);
    const auto *sgl0_44 = buffer.data(sgl0 + 44);

    const auto *sgk_0 = buffer.data(sgk + 0);
    const auto *sgk_1 = buffer.data(sgk + 1);
    const auto *sgk_2 = buffer.data(sgk + 2);
    const auto *sgk_3 = buffer.data(sgk + 3);
    const auto *sgk_5 = buffer.data(sgk + 5);
    const auto *sgk_6 = buffer.data(sgk + 6);
    const auto *sgk_7 = buffer.data(sgk + 7);
    const auto *sgk_8 = buffer.data(sgk + 8);
    const auto *sgk_9 = buffer.data(sgk + 9);
    const auto *sgk_10 = buffer.data(sgk + 10);
    const auto *sgk_11 = buffer.data(sgk + 11);
    const auto *sgk_12 = buffer.data(sgk + 12);
    const auto *sgk_13 = buffer.data(sgk + 13);
    const auto *sgk_14 = buffer.data(sgk + 14);
    const auto *sgk_15 = buffer.data(sgk + 15);
    const auto *sgk_16 = buffer.data(sgk + 16);
    const auto *sgk_17 = buffer.data(sgk + 17);
    const auto *sgk_18 = buffer.data(sgk + 18);
    const auto *sgk_19 = buffer.data(sgk + 19);
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
    const auto *sgk_64 = buffer.data(sgk + 64);
    const auto *sgk_65 = buffer.data(sgk + 65);
    const auto *sgk_66 = buffer.data(sgk + 66);
    const auto *sgk_67 = buffer.data(sgk + 67);
    const auto *sgk_68 = buffer.data(sgk + 68);
    const auto *sgk_69 = buffer.data(sgk + 69);
    const auto *sgk_70 = buffer.data(sgk + 70);
    const auto *sgk_71 = buffer.data(sgk + 71);

    const auto *sgl1_0 = buffer.data(sgl1 + 0);
    const auto *sgl1_3 = buffer.data(sgl1 + 3);
    const auto *sgl1_5 = buffer.data(sgl1 + 5);
    const auto *sgl1_6 = buffer.data(sgl1 + 6);
    const auto *sgl1_9 = buffer.data(sgl1 + 9);
    const auto *sgl1_10 = buffer.data(sgl1 + 10);
    const auto *sgl1_12 = buffer.data(sgl1 + 12);
    const auto *sgl1_14 = buffer.data(sgl1 + 14);
    const auto *sgl1_15 = buffer.data(sgl1 + 15);
    const auto *sgl1_17 = buffer.data(sgl1 + 17);
    const auto *sgl1_18 = buffer.data(sgl1 + 18);
    const auto *sgl1_20 = buffer.data(sgl1 + 20);
    const auto *sgl1_21 = buffer.data(sgl1 + 21);
    const auto *sgl1_23 = buffer.data(sgl1 + 23);
    const auto *sgl1_24 = buffer.data(sgl1 + 24);
    const auto *sgl1_25 = buffer.data(sgl1 + 25);
    const auto *sgl1_27 = buffer.data(sgl1 + 27);
    const auto *sgl1_44 = buffer.data(sgl1 + 44);

    const auto *shi0_0 = buffer.data(shi0 + 0);
    const auto *shi0_3 = buffer.data(shi0 + 3);
    const auto *shi0_5 = buffer.data(shi0 + 5);
    const auto *shi0_6 = buffer.data(shi0 + 6);
    const auto *shi0_9 = buffer.data(shi0 + 9);
    const auto *shi0_10 = buffer.data(shi0 + 10);
    const auto *shi0_12 = buffer.data(shi0 + 12);
    const auto *shi0_14 = buffer.data(shi0 + 14);
    const auto *shi0_15 = buffer.data(shi0 + 15);
    const auto *shi0_17 = buffer.data(shi0 + 17);
    const auto *shi0_18 = buffer.data(shi0 + 18);
    const auto *shi0_20 = buffer.data(shi0 + 20);
    const auto *shi0_21 = buffer.data(shi0 + 21);
    const auto *shi0_23 = buffer.data(shi0 + 23);
    const auto *shi0_24 = buffer.data(shi0 + 24);
    const auto *shi0_25 = buffer.data(shi0 + 25);
    const auto *shi0_26 = buffer.data(shi0 + 26);
    const auto *shi0_27 = buffer.data(shi0 + 27);
    const auto *shi0_49 = buffer.data(shi0 + 49);
    const auto *shi0_51 = buffer.data(shi0 + 51);
    const auto *shi0_52 = buffer.data(shi0 + 52);
    const auto *shi0_53 = buffer.data(shi0 + 53);
    const auto *shi0_54 = buffer.data(shi0 + 54);
    const auto *shi0_55 = buffer.data(shi0 + 55);

    const auto *shi1_0 = buffer.data(shi1 + 0);
    const auto *shi1_3 = buffer.data(shi1 + 3);
    const auto *shi1_5 = buffer.data(shi1 + 5);
    const auto *shi1_6 = buffer.data(shi1 + 6);
    const auto *shi1_9 = buffer.data(shi1 + 9);
    const auto *shi1_10 = buffer.data(shi1 + 10);
    const auto *shi1_12 = buffer.data(shi1 + 12);
    const auto *shi1_14 = buffer.data(shi1 + 14);
    const auto *shi1_15 = buffer.data(shi1 + 15);
    const auto *shi1_17 = buffer.data(shi1 + 17);
    const auto *shi1_18 = buffer.data(shi1 + 18);
    const auto *shi1_20 = buffer.data(shi1 + 20);
    const auto *shi1_21 = buffer.data(shi1 + 21);
    const auto *shi1_23 = buffer.data(shi1 + 23);
    const auto *shi1_24 = buffer.data(shi1 + 24);
    const auto *shi1_25 = buffer.data(shi1 + 25);
    const auto *shi1_26 = buffer.data(shi1 + 26);
    const auto *shi1_27 = buffer.data(shi1 + 27);
    const auto *shi1_49 = buffer.data(shi1 + 49);
    const auto *shi1_51 = buffer.data(shi1 + 51);
    const auto *shi1_52 = buffer.data(shi1 + 52);
    const auto *shi1_53 = buffer.data(shi1 + 53);
    const auto *shi1_54 = buffer.data(shi1 + 54);
    const auto *shi1_55 = buffer.data(shi1 + 55);

    const auto *shk_0 = buffer.data(shk + 0);
    const auto *shk_2 = buffer.data(shk + 2);
    const auto *shk_3 = buffer.data(shk + 3);
    const auto *shk_5 = buffer.data(shk + 5);
    const auto *shk_6 = buffer.data(shk + 6);
    const auto *shk_9 = buffer.data(shk + 9);
    const auto *shk_10 = buffer.data(shk + 10);
    const auto *shk_12 = buffer.data(shk + 12);
    const auto *shk_14 = buffer.data(shk + 14);
    const auto *shk_15 = buffer.data(shk + 15);
    const auto *shk_17 = buffer.data(shk + 17);
    const auto *shk_18 = buffer.data(shk + 18);
    const auto *shk_20 = buffer.data(shk + 20);
    const auto *shk_21 = buffer.data(shk + 21);
    const auto *shk_23 = buffer.data(shk + 23);
    const auto *shk_24 = buffer.data(shk + 24);
    const auto *shk_25 = buffer.data(shk + 25);
    const auto *shk_27 = buffer.data(shk + 27);
    const auto *shk_28 = buffer.data(shk + 28);
    const auto *shk_29 = buffer.data(shk + 29);
    const auto *shk_30 = buffer.data(shk + 30);
    const auto *shk_31 = buffer.data(shk + 31);
    const auto *shk_32 = buffer.data(shk + 32);
    const auto *shk_33 = buffer.data(shk + 33);
    const auto *shk_34 = buffer.data(shk + 34);
    const auto *shk_35 = buffer.data(shk + 35);
    const auto *shk_36 = buffer.data(shk + 36);
    const auto *shk_38 = buffer.data(shk + 38);
    const auto *shk_39 = buffer.data(shk + 39);
    const auto *shk_41 = buffer.data(shk + 41);
    const auto *shk_42 = buffer.data(shk + 42);
    const auto *shk_45 = buffer.data(shk + 45);
    const auto *shk_46 = buffer.data(shk + 46);
    const auto *shk_50 = buffer.data(shk + 50);
    const auto *shk_51 = buffer.data(shk + 51);
    const auto *shk_56 = buffer.data(shk + 56);
    const auto *shk_64 = buffer.data(shk + 64);
    const auto *shk_65 = buffer.data(shk + 65);
    const auto *shk_66 = buffer.data(shk + 66);
    const auto *shk_67 = buffer.data(shk + 67);
    const auto *shk_68 = buffer.data(shk + 68);
    const auto *shk_69 = buffer.data(shk + 69);
    const auto *shk_70 = buffer.data(shk + 70);
    const auto *shk_71 = buffer.data(shk + 71);
    const auto *shk_72 = buffer.data(shk + 72);
    const auto *shk_74 = buffer.data(shk + 74);
    const auto *shk_75 = buffer.data(shk + 75);
    const auto *shk_77 = buffer.data(shk + 77);
    const auto *shk_78 = buffer.data(shk + 78);
    const auto *shk_81 = buffer.data(shk + 81);
    const auto *shk_82 = buffer.data(shk + 82);
    const auto *shk_86 = buffer.data(shk + 86);
    const auto *shk_87 = buffer.data(shk + 87);
    const auto *shk_92 = buffer.data(shk + 92);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sgk_0, sgk_3, shi0_0, shi0_3, \
                         shi1_0, shi1_3, shk_0, shk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sgk_0[k]
                 + f_1 * shi0_0[k]
                 - f_2 * shi1_0[k]
                 + f_3 * pc_x[k] * shk_0[k];

        t_1[k] = f_3 * pc_y[k] * shk_0[k];

        t_2[k] = f_3 * pc_z[k] * shk_0[k];

        t_3[k] = f_0 * sgk_3[k]
                 + f_4 * shi0_3[k]
                 - f_5 * shi1_3[k]
                 + f_3 * pc_x[k] * shk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sgk_5, sgk_6, shi0_5, shi0_6, shi1_5, \
                         shi1_6, shk_2, shk_5, shk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * shk_2[k];

        t_5[k] = f_0 * sgk_5[k]
                 + f_4 * shi0_5[k]
                 - f_5 * shi1_5[k]
                 + f_3 * pc_x[k] * shk_5[k];

        t_6[k] = f_0 * sgk_6[k]
                 + f_6 * shi0_6[k]
                 - f_7 * shi1_6[k]
                 + f_3 * pc_x[k] * shk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sgk_9, shi0_9, shi1_9, shk_3, shk_5, \
                         shk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * shk_3[k];

        t_8[k] = f_3 * pc_y[k] * shk_5[k];

        t_9[k] = f_0 * sgk_9[k]
                 + f_6 * shi0_9[k]
                 - f_7 * shi1_9[k]
                 + f_3 * pc_x[k] * shk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sgk_10, sgk_12, shi0_10, shi0_12, \
                         shi1_10, shi1_12, shk_6, shk_10, shk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sgk_10[k]
                  + f_8 * shi0_10[k]
                  - f_9 * shi1_10[k]
                  + f_3 * pc_x[k] * shk_10[k];

        t_11[k] = f_3 * pc_z[k] * shk_6[k];

        t_12[k] = f_0 * sgk_12[k]
                  + f_8 * shi0_12[k]
                  - f_9 * shi1_12[k]
                  + f_3 * pc_x[k] * shk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, sgk_14, sgk_15, shi0_14, shi0_15, \
                         shi1_14, shi1_15, shk_9, shk_14, shk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * shk_9[k];

        t_14[k] = f_0 * sgk_14[k]
                  + f_8 * shi0_14[k]
                  - f_9 * shi1_14[k]
                  + f_3 * pc_x[k] * shk_14[k];

        t_15[k] = f_0 * sgk_15[k]
                  + f_10 * shi0_15[k]
                  - f_11 * shi1_15[k]
                  + f_3 * pc_x[k] * shk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, sgk_17, sgk_18, shi0_17, shi0_18, \
                         shi1_17, shi1_18, shk_10, shk_17, shk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * shk_10[k];

        t_17[k] = f_0 * sgk_17[k]
                  + f_10 * shi0_17[k]
                  - f_11 * shi1_17[k]
                  + f_3 * pc_x[k] * shk_17[k];

        t_18[k] = f_0 * sgk_18[k]
                  + f_10 * shi0_18[k]
                  - f_11 * shi1_18[k]
                  + f_3 * pc_x[k] * shk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pc_x, pc_y, sgk_20, sgk_21, shi0_20, shi0_21, \
                         shi1_20, shi1_21, shk_14, shk_20, shk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * shk_14[k];

        t_20[k] = f_0 * sgk_20[k]
                  + f_10 * shi0_20[k]
                  - f_11 * shi1_20[k]
                  + f_3 * pc_x[k] * shk_20[k];

        t_21[k] = f_0 * sgk_21[k]
                  + f_12 * shi0_21[k]
                  - f_13 * shi1_21[k]
                  + f_3 * pc_x[k] * shk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pc_x, pc_z, sgk_23, sgk_24, shi0_23, shi0_24, \
                         shi1_23, shi1_24, shk_15, shk_23, shk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * pc_z[k] * shk_15[k];

        t_23[k] = f_0 * sgk_23[k]
                  + f_12 * shi0_23[k]
                  - f_13 * shi1_23[k]
                  + f_3 * pc_x[k] * shk_23[k];

        t_24[k] = f_0 * sgk_24[k]
                  + f_12 * shi0_24[k]
                  - f_13 * shi1_24[k]
                  + f_3 * pc_x[k] * shk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pc_x, pc_y, sgk_25, sgk_27, shi0_25, shi0_27, \
                         shi1_25, shi1_27, shk_20, shk_25, shk_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * sgk_25[k]
                  + f_12 * shi0_25[k]
                  - f_13 * shi1_25[k]
                  + f_3 * pc_x[k] * shk_25[k];

        t_26[k] = f_3 * pc_y[k] * shk_20[k];

        t_27[k] = f_0 * sgk_27[k]
                  + f_12 * shi0_27[k]
                  - f_13 * shi1_27[k]
                  + f_3 * pc_x[k] * shk_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pc_x, sgk_28, sgk_29, sgk_30, sgk_31, \
                         sgk_32, shk_28, shk_29, shk_30, shk_31, \
                         shk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * sgk_28[k]
                  + f_3 * pc_x[k] * shk_28[k];

        t_29[k] = f_0 * sgk_29[k]
                  + f_3 * pc_x[k] * shk_29[k];

        t_30[k] = f_0 * sgk_30[k]
                  + f_3 * pc_x[k] * shk_30[k];

        t_31[k] = f_0 * sgk_31[k]
                  + f_3 * pc_x[k] * shk_31[k];

        t_32[k] = f_0 * sgk_32[k]
                  + f_3 * pc_x[k] * shk_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pc_x, pc_y, sgk_33, sgk_34, sgk_35, shi0_21, \
                         shi1_21, shk_28, shk_33, shk_34, shk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * sgk_33[k]
                  + f_3 * pc_x[k] * shk_33[k];

        t_34[k] = f_0 * sgk_34[k]
                  + f_3 * pc_x[k] * shk_34[k];

        t_35[k] = f_0 * sgk_35[k]
                  + f_3 * pc_x[k] * shk_35[k];

        t_36[k] = f_1 * shi0_21[k]
                  - f_2 * shi1_21[k]
                  + f_3 * pc_y[k] * shk_28[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pc_y, pc_z, shi0_23, shi0_24, shi0_25, \
                         shi1_23, shi1_24, shi1_25, shk_28, shk_30, shk_31, \
                         shk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * pc_z[k] * shk_28[k];

        t_38[k] = f_4 * shi0_23[k]
                  - f_5 * shi1_23[k]
                  + f_3 * pc_y[k] * shk_30[k];

        t_39[k] = f_6 * shi0_24[k]
                  - f_7 * shi1_24[k]
                  + f_3 * pc_y[k] * shk_31[k];

        t_40[k] = f_8 * shi0_25[k]
                  - f_9 * shi1_25[k]
                  + f_3 * pc_y[k] * shk_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, shi0_26, shi0_27, shi1_26, \
                         shi1_27, shk_33, shk_34, shk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * shi0_26[k]
                  - f_11 * shi1_26[k]
                  + f_3 * pc_y[k] * shk_33[k];

        t_42[k] = f_12 * shi0_27[k]
                  - f_13 * shi1_27[k]
                  + f_3 * pc_y[k] * shk_34[k];

        t_43[k] = f_3 * pc_y[k] * shk_35[k];

        t_44[k] = f_1 * shi0_27[k]
                  - f_2 * shi1_27[k]
                  + f_3 * pc_z[k] * shk_35[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_y, pc_y, pc_z, sgl0_0, sgl0_3, sgk_0, \
                         sgk_1, sgl1_0, sgl1_3, shk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pb_y[k] * sgl0_0[k]
                  - f_14 * pc_y[k] * sgl1_0[k];

        t_46[k] = f_15 * sgk_0[k]
                  + f_3 * pc_y[k] * shk_36[k];

        t_47[k] = f_3 * pc_z[k] * shk_36[k];

        t_48[k] = pb_y[k] * sgl0_3[k]
                  + f_16 * sgk_1[k]
                  - f_14 * pc_y[k] * sgl1_3[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_y, pc_y, pc_z, sgl0_5, sgl0_6, sgk_2, \
                         sgk_3, sgl1_5, sgl1_6, shk_38, shk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_15 * sgk_2[k]
                  + f_3 * pc_y[k] * shk_38[k];

        t_50[k] = pb_y[k] * sgl0_5[k]
                  - f_14 * pc_y[k] * sgl1_5[k];

        t_51[k] = pb_y[k] * sgl0_6[k]
                  + f_17 * sgk_3[k]
                  - f_14 * pc_y[k] * sgl1_6[k];

        t_52[k] = f_3 * pc_z[k] * shk_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_y, pc_y, pc_z, sgl0_9, sgl0_10, sgk_5, \
                         sgk_6, sgl1_9, sgl1_10, shk_41, shk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_15 * sgk_5[k]
                  + f_3 * pc_y[k] * shk_41[k];

        t_54[k] = pb_y[k] * sgl0_9[k]
                  - f_14 * pc_y[k] * sgl1_9[k];

        t_55[k] = pb_y[k] * sgl0_10[k]
                  + f_18 * sgk_6[k]
                  - f_14 * pc_y[k] * sgl1_10[k];

        t_56[k] = f_3 * pc_z[k] * shk_42[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, sgl0_12, sgl0_14, sgl0_15, sgk_8, \
                         sgk_9, sgk_10, sgl1_12, sgl1_14, sgl1_15, \
                         shk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_y[k] * sgl0_12[k]
                  + f_16 * sgk_8[k]
                  - f_14 * pc_y[k] * sgl1_12[k];

        t_58[k] = f_15 * sgk_9[k]
                  + f_3 * pc_y[k] * shk_45[k];

        t_59[k] = pb_y[k] * sgl0_14[k]
                  - f_14 * pc_y[k] * sgl1_14[k];

        t_60[k] = pb_y[k] * sgl0_15[k]
                  + f_0 * sgk_10[k]
                  - f_14 * pc_y[k] * sgl1_15[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_y, pc_y, pc_z, sgl0_17, sgl0_18, sgk_12, \
                         sgk_13, sgk_14, sgl1_17, sgl1_18, shk_46, \
                         shk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * pc_z[k] * shk_46[k];

        t_62[k] = pb_y[k] * sgl0_17[k]
                  + f_17 * sgk_12[k]
                  - f_14 * pc_y[k] * sgl1_17[k];

        t_63[k] = pb_y[k] * sgl0_18[k]
                  + f_16 * sgk_13[k]
                  - f_14 * pc_y[k] * sgl1_18[k];

        t_64[k] = f_15 * sgk_14[k]
                  + f_3 * pc_y[k] * shk_50[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_y, pc_z, sgl0_20, sgl0_21, sgl0_23, \
                         sgk_15, sgk_17, sgl1_20, sgl1_21, sgl1_23, \
                         shk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * sgl0_20[k]
                  - f_14 * pc_y[k] * sgl1_20[k];

        t_66[k] = pb_y[k] * sgl0_21[k]
                  + f_19 * sgk_15[k]
                  - f_14 * pc_y[k] * sgl1_21[k];

        t_67[k] = f_3 * pc_z[k] * shk_51[k];

        t_68[k] = pb_y[k] * sgl0_23[k]
                  + f_18 * sgk_17[k]
                  - f_14 * pc_y[k] * sgl1_23[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_y, pc_y, sgl0_24, sgl0_25, sgl0_27, \
                         sgk_18, sgk_19, sgk_20, sgl1_24, sgl1_25, sgl1_27, \
                         shk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_y[k] * sgl0_24[k]
                  + f_17 * sgk_18[k]
                  - f_14 * pc_y[k] * sgl1_24[k];

        t_70[k] = pb_y[k] * sgl0_25[k]
                  + f_16 * sgk_19[k]
                  - f_14 * pc_y[k] * sgl1_25[k];

        t_71[k] = f_15 * sgk_20[k]
                  + f_3 * pc_y[k] * shk_56[k];

        t_72[k] = pb_y[k] * sgl0_27[k]
                  - f_14 * pc_y[k] * sgl1_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pc_x, sgk_64, sgk_65, sgk_66, sgk_67, \
                         sgk_68, shk_64, shk_65, shk_66, shk_67, \
                         shk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_18 * sgk_64[k]
                  + f_3 * pc_x[k] * shk_64[k];

        t_74[k] = f_18 * sgk_65[k]
                  + f_3 * pc_x[k] * shk_65[k];

        t_75[k] = f_18 * sgk_66[k]
                  + f_3 * pc_x[k] * shk_66[k];

        t_76[k] = f_18 * sgk_67[k]
                  + f_3 * pc_x[k] * shk_67[k];

        t_77[k] = f_18 * sgk_68[k]
                  + f_3 * pc_x[k] * shk_68[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pc_x, pc_y, sgk_28, sgk_69, sgk_70, sgk_71, \
                         shi0_49, shi1_49, shk_64, shk_69, shk_70, \
                         shk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_18 * sgk_69[k]
                  + f_3 * pc_x[k] * shk_69[k];

        t_79[k] = f_18 * sgk_70[k]
                  + f_3 * pc_x[k] * shk_70[k];

        t_80[k] = f_18 * sgk_71[k]
                  + f_3 * pc_x[k] * shk_71[k];

        t_81[k] = f_15 * sgk_28[k]
                  + f_1 * shi0_49[k]
                  - f_2 * shi1_49[k]
                  + f_3 * pc_y[k] * shk_64[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pc_y, pc_z, sgk_30, sgk_31, shi0_51, shi0_52, \
                         shi1_51, shi1_52, shk_64, shk_66, shk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * pc_z[k] * shk_64[k];

        t_83[k] = f_15 * sgk_30[k]
                  + f_4 * shi0_51[k]
                  - f_5 * shi1_51[k]
                  + f_3 * pc_y[k] * shk_66[k];

        t_84[k] = f_15 * sgk_31[k]
                  + f_6 * shi0_52[k]
                  - f_7 * shi1_52[k]
                  + f_3 * pc_y[k] * shk_67[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_y, sgk_32, sgk_33, sgk_34, shi0_53, shi0_54, \
                         shi0_55, shi1_53, shi1_54, shi1_55, shk_68, shk_69, \
                         shk_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_15 * sgk_32[k]
                  + f_8 * shi0_53[k]
                  - f_9 * shi1_53[k]
                  + f_3 * pc_y[k] * shk_68[k];

        t_86[k] = f_15 * sgk_33[k]
                  + f_10 * shi0_54[k]
                  - f_11 * shi1_54[k]
                  + f_3 * pc_y[k] * shk_69[k];

        t_87[k] = f_15 * sgk_34[k]
                  + f_12 * shi0_55[k]
                  - f_13 * shi1_55[k]
                  + f_3 * pc_y[k] * shk_70[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pb_y, pb_z, pc_y, pc_z, sgl0_0, sgl0_44, \
                         sgk_35, sgl1_0, sgl1_44, shk_71, shk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * sgk_35[k]
                  + f_3 * pc_y[k] * shk_71[k];

        t_89[k] = pb_y[k] * sgl0_44[k]
                  - f_14 * pc_y[k] * sgl1_44[k];

        t_90[k] = pb_z[k] * sgl0_0[k]
                  - f_14 * pc_z[k] * sgl1_0[k];

        t_91[k] = f_3 * pc_y[k] * shk_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_z, pc_y, pc_z, sgl0_3, sgl0_5, sgk_0, \
                         sgk_2, sgl1_3, sgl1_5, shk_72, shk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * sgk_0[k]
                  + f_3 * pc_z[k] * shk_72[k];

        t_93[k] = pb_z[k] * sgl0_3[k]
                  - f_14 * pc_z[k] * sgl1_3[k];

        t_94[k] = f_3 * pc_y[k] * shk_74[k];

        t_95[k] = pb_z[k] * sgl0_5[k]
                  + f_16 * sgk_2[k]
                  - f_14 * pc_z[k] * sgl1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_z, pc_y, pc_z, sgl0_6, sgl0_9, sgk_3, \
                         sgk_5, sgl1_6, sgl1_9, shk_75, shk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pb_z[k] * sgl0_6[k]
                  - f_14 * pc_z[k] * sgl1_6[k];

        t_97[k] = f_15 * sgk_3[k]
                  + f_3 * pc_z[k] * shk_75[k];

        t_98[k] = f_3 * pc_y[k] * shk_77[k];

        t_99[k] = pb_z[k] * sgl0_9[k]
                  + f_17 * sgk_5[k]
                  - f_14 * pc_z[k] * sgl1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pb_z, pc_y, pc_z, sgl0_10, sgl0_12, \
                         sgk_6, sgk_7, sgl1_10, sgl1_12, shk_78, \
                         shk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pb_z[k] * sgl0_10[k]
                   - f_14 * pc_z[k] * sgl1_10[k];

        t_101[k] = f_15 * sgk_6[k]
                   + f_3 * pc_z[k] * shk_78[k];

        t_102[k] = pb_z[k] * sgl0_12[k]
                   + f_16 * sgk_7[k]
                   - f_14 * pc_z[k] * sgl1_12[k];

        t_103[k] = f_3 * pc_y[k] * shk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_z, pc_z, sgl0_14, sgl0_15, sgl0_17, \
                         sgk_9, sgk_10, sgk_11, sgl1_14, sgl1_15, sgl1_17, \
                         shk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_z[k] * sgl0_14[k]
                   + f_18 * sgk_9[k]
                   - f_14 * pc_z[k] * sgl1_14[k];

        t_105[k] = pb_z[k] * sgl0_15[k]
                   - f_14 * pc_z[k] * sgl1_15[k];

        t_106[k] = f_15 * sgk_10[k]
                   + f_3 * pc_z[k] * shk_82[k];

        t_107[k] = pb_z[k] * sgl0_17[k]
                   + f_16 * sgk_11[k]
                   - f_14 * pc_z[k] * sgl1_17[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pb_z, pc_y, pc_z, sgl0_18, sgl0_20, \
                         sgl0_21, sgk_12, sgk_14, sgl1_18, sgl1_20, sgl1_21, \
                         shk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_z[k] * sgl0_18[k]
                   + f_17 * sgk_12[k]
                   - f_14 * pc_z[k] * sgl1_18[k];

        t_109[k] = f_3 * pc_y[k] * shk_86[k];

        t_110[k] = pb_z[k] * sgl0_20[k]
                   + f_0 * sgk_14[k]
                   - f_14 * pc_z[k] * sgl1_20[k];

        t_111[k] = pb_z[k] * sgl0_21[k]
                   - f_14 * pc_z[k] * sgl1_21[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_z, pc_z, sgl0_23, sgl0_24, sgk_15, sgk_16, \
                         sgk_17, sgl1_23, sgl1_24, shk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_15 * sgk_15[k]
                   + f_3 * pc_z[k] * shk_87[k];

        t_113[k] = pb_z[k] * sgl0_23[k]
                   + f_16 * sgk_16[k]
                   - f_14 * pc_z[k] * sgl1_23[k];

        t_114[k] = pb_z[k] * sgl0_24[k]
                   + f_17 * sgk_17[k]
                   - f_14 * pc_z[k] * sgl1_24[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pb_z, pc_y, pc_z, sgl0_25, sgl0_27, sgk_18, \
                         sgk_20, sgl1_25, sgl1_27, shk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_z[k] * sgl0_25[k]
                   + f_18 * sgk_18[k]
                   - f_14 * pc_z[k] * sgl1_25[k];

        t_116[k] = f_3 * pc_y[k] * shk_92[k];

        t_117[k] = pb_z[k] * sgl0_27[k]
                   + f_19 * sgk_20[k]
                   - f_14 * pc_z[k] * sgl1_27[k];
    }
}

static auto
compute_prim_shl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgl0,
                                                          const size_t sgk, const size_t sgl1,
                                                          const size_t shi0, const size_t shi1,
                                                          const size_t shk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_18 = 2.0 / q;

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

    const auto *sgl0_36 = buffer.data(sgl0 + 36);
    const auto *sgl0_48 = buffer.data(sgl0 + 48);
    const auto *sgl0_51 = buffer.data(sgl0 + 51);
    const auto *sgl0_55 = buffer.data(sgl0 + 55);
    const auto *sgl0_60 = buffer.data(sgl0 + 60);
    const auto *sgl0_66 = buffer.data(sgl0 + 66);
    const auto *sgl0_81 = buffer.data(sgl0 + 81);
    const auto *sgl0_90 = buffer.data(sgl0 + 90);
    const auto *sgl0_95 = buffer.data(sgl0 + 95);
    const auto *sgl0_99 = buffer.data(sgl0 + 99);
    const auto *sgl0_102 = buffer.data(sgl0 + 102);
    const auto *sgl0_104 = buffer.data(sgl0 + 104);
    const auto *sgl0_107 = buffer.data(sgl0 + 107);
    const auto *sgl0_108 = buffer.data(sgl0 + 108);
    const auto *sgl0_110 = buffer.data(sgl0 + 110);
    const auto *sgl0_113 = buffer.data(sgl0 + 113);
    const auto *sgl0_114 = buffer.data(sgl0 + 114);
    const auto *sgl0_115 = buffer.data(sgl0 + 115);
    const auto *sgl0_117 = buffer.data(sgl0 + 117);
    const auto *sgl0_134 = buffer.data(sgl0 + 134);

    const auto *sgk_28 = buffer.data(sgk + 28);
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
    const auto *sgk_80 = buffer.data(sgk + 80);
    const auto *sgk_81 = buffer.data(sgk + 81);
    const auto *sgk_84 = buffer.data(sgk + 84);
    const auto *sgk_85 = buffer.data(sgk + 85);
    const auto *sgk_86 = buffer.data(sgk + 86);
    const auto *sgk_89 = buffer.data(sgk + 89);
    const auto *sgk_90 = buffer.data(sgk + 90);
    const auto *sgk_91 = buffer.data(sgk + 91);
    const auto *sgk_92 = buffer.data(sgk + 92);
    const auto *sgk_100 = buffer.data(sgk + 100);
    const auto *sgk_101 = buffer.data(sgk + 101);
    const auto *sgk_102 = buffer.data(sgk + 102);
    const auto *sgk_103 = buffer.data(sgk + 103);
    const auto *sgk_104 = buffer.data(sgk + 104);
    const auto *sgk_105 = buffer.data(sgk + 105);
    const auto *sgk_106 = buffer.data(sgk + 106);
    const auto *sgk_107 = buffer.data(sgk + 107);
    const auto *sgk_108 = buffer.data(sgk + 108);
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
    const auto *sgk_172 = buffer.data(sgk + 172);
    const auto *sgk_173 = buffer.data(sgk + 173);
    const auto *sgk_174 = buffer.data(sgk + 174);
    const auto *sgk_175 = buffer.data(sgk + 175);
    const auto *sgk_176 = buffer.data(sgk + 176);
    const auto *sgk_177 = buffer.data(sgk + 177);
    const auto *sgk_178 = buffer.data(sgk + 178);
    const auto *sgk_179 = buffer.data(sgk + 179);
    const auto *sgk_180 = buffer.data(sgk + 180);
    const auto *sgk_183 = buffer.data(sgk + 183);
    const auto *sgk_185 = buffer.data(sgk + 185);
    const auto *sgk_186 = buffer.data(sgk + 186);

    const auto *sgl1_36 = buffer.data(sgl1 + 36);
    const auto *sgl1_48 = buffer.data(sgl1 + 48);
    const auto *sgl1_51 = buffer.data(sgl1 + 51);
    const auto *sgl1_55 = buffer.data(sgl1 + 55);
    const auto *sgl1_60 = buffer.data(sgl1 + 60);
    const auto *sgl1_66 = buffer.data(sgl1 + 66);
    const auto *sgl1_81 = buffer.data(sgl1 + 81);
    const auto *sgl1_90 = buffer.data(sgl1 + 90);
    const auto *sgl1_95 = buffer.data(sgl1 + 95);
    const auto *sgl1_99 = buffer.data(sgl1 + 99);
    const auto *sgl1_102 = buffer.data(sgl1 + 102);
    const auto *sgl1_104 = buffer.data(sgl1 + 104);
    const auto *sgl1_107 = buffer.data(sgl1 + 107);
    const auto *sgl1_108 = buffer.data(sgl1 + 108);
    const auto *sgl1_110 = buffer.data(sgl1 + 110);
    const auto *sgl1_113 = buffer.data(sgl1 + 113);
    const auto *sgl1_114 = buffer.data(sgl1 + 114);
    const auto *sgl1_115 = buffer.data(sgl1 + 115);
    const auto *sgl1_117 = buffer.data(sgl1 + 117);
    const auto *sgl1_134 = buffer.data(sgl1 + 134);

    const auto *shi0_79 = buffer.data(shi0 + 79);
    const auto *shi0_80 = buffer.data(shi0 + 80);
    const auto *shi0_81 = buffer.data(shi0 + 81);
    const auto *shi0_82 = buffer.data(shi0 + 82);
    const auto *shi0_83 = buffer.data(shi0 + 83);
    const auto *shi0_84 = buffer.data(shi0 + 84);
    const auto *shi0_87 = buffer.data(shi0 + 87);
    const auto *shi0_89 = buffer.data(shi0 + 89);
    const auto *shi0_90 = buffer.data(shi0 + 90);
    const auto *shi0_93 = buffer.data(shi0 + 93);
    const auto *shi0_94 = buffer.data(shi0 + 94);
    const auto *shi0_96 = buffer.data(shi0 + 96);
    const auto *shi0_98 = buffer.data(shi0 + 98);
    const auto *shi0_99 = buffer.data(shi0 + 99);
    const auto *shi0_101 = buffer.data(shi0 + 101);
    const auto *shi0_102 = buffer.data(shi0 + 102);
    const auto *shi0_104 = buffer.data(shi0 + 104);
    const auto *shi0_105 = buffer.data(shi0 + 105);
    const auto *shi0_107 = buffer.data(shi0 + 107);
    const auto *shi0_108 = buffer.data(shi0 + 108);
    const auto *shi0_109 = buffer.data(shi0 + 109);
    const auto *shi0_110 = buffer.data(shi0 + 110);
    const auto *shi0_111 = buffer.data(shi0 + 111);
    const auto *shi0_135 = buffer.data(shi0 + 135);
    const auto *shi0_136 = buffer.data(shi0 + 136);
    const auto *shi0_137 = buffer.data(shi0 + 137);
    const auto *shi0_138 = buffer.data(shi0 + 138);
    const auto *shi0_139 = buffer.data(shi0 + 139);
    const auto *shi0_140 = buffer.data(shi0 + 140);
    const auto *shi0_143 = buffer.data(shi0 + 143);
    const auto *shi0_145 = buffer.data(shi0 + 145);
    const auto *shi0_146 = buffer.data(shi0 + 146);

    const auto *shi1_79 = buffer.data(shi1 + 79);
    const auto *shi1_80 = buffer.data(shi1 + 80);
    const auto *shi1_81 = buffer.data(shi1 + 81);
    const auto *shi1_82 = buffer.data(shi1 + 82);
    const auto *shi1_83 = buffer.data(shi1 + 83);
    const auto *shi1_84 = buffer.data(shi1 + 84);
    const auto *shi1_87 = buffer.data(shi1 + 87);
    const auto *shi1_89 = buffer.data(shi1 + 89);
    const auto *shi1_90 = buffer.data(shi1 + 90);
    const auto *shi1_93 = buffer.data(shi1 + 93);
    const auto *shi1_94 = buffer.data(shi1 + 94);
    const auto *shi1_96 = buffer.data(shi1 + 96);
    const auto *shi1_98 = buffer.data(shi1 + 98);
    const auto *shi1_99 = buffer.data(shi1 + 99);
    const auto *shi1_101 = buffer.data(shi1 + 101);
    const auto *shi1_102 = buffer.data(shi1 + 102);
    const auto *shi1_104 = buffer.data(shi1 + 104);
    const auto *shi1_105 = buffer.data(shi1 + 105);
    const auto *shi1_107 = buffer.data(shi1 + 107);
    const auto *shi1_108 = buffer.data(shi1 + 108);
    const auto *shi1_109 = buffer.data(shi1 + 109);
    const auto *shi1_110 = buffer.data(shi1 + 110);
    const auto *shi1_111 = buffer.data(shi1 + 111);
    const auto *shi1_135 = buffer.data(shi1 + 135);
    const auto *shi1_136 = buffer.data(shi1 + 136);
    const auto *shi1_137 = buffer.data(shi1 + 137);
    const auto *shi1_138 = buffer.data(shi1 + 138);
    const auto *shi1_139 = buffer.data(shi1 + 139);
    const auto *shi1_140 = buffer.data(shi1 + 140);
    const auto *shi1_143 = buffer.data(shi1 + 143);
    const auto *shi1_145 = buffer.data(shi1 + 145);
    const auto *shi1_146 = buffer.data(shi1 + 146);

    const auto *shk_100 = buffer.data(shk + 100);
    const auto *shk_101 = buffer.data(shk + 101);
    const auto *shk_102 = buffer.data(shk + 102);
    const auto *shk_103 = buffer.data(shk + 103);
    const auto *shk_104 = buffer.data(shk + 104);
    const auto *shk_105 = buffer.data(shk + 105);
    const auto *shk_106 = buffer.data(shk + 106);
    const auto *shk_107 = buffer.data(shk + 107);
    const auto *shk_108 = buffer.data(shk + 108);
    const auto *shk_110 = buffer.data(shk + 110);
    const auto *shk_111 = buffer.data(shk + 111);
    const auto *shk_113 = buffer.data(shk + 113);
    const auto *shk_114 = buffer.data(shk + 114);
    const auto *shk_117 = buffer.data(shk + 117);
    const auto *shk_118 = buffer.data(shk + 118);
    const auto *shk_120 = buffer.data(shk + 120);
    const auto *shk_122 = buffer.data(shk + 122);
    const auto *shk_123 = buffer.data(shk + 123);
    const auto *shk_125 = buffer.data(shk + 125);
    const auto *shk_126 = buffer.data(shk + 126);
    const auto *shk_128 = buffer.data(shk + 128);
    const auto *shk_129 = buffer.data(shk + 129);
    const auto *shk_131 = buffer.data(shk + 131);
    const auto *shk_132 = buffer.data(shk + 132);
    const auto *shk_133 = buffer.data(shk + 133);
    const auto *shk_135 = buffer.data(shk + 135);
    const auto *shk_136 = buffer.data(shk + 136);
    const auto *shk_137 = buffer.data(shk + 137);
    const auto *shk_138 = buffer.data(shk + 138);
    const auto *shk_139 = buffer.data(shk + 139);
    const auto *shk_140 = buffer.data(shk + 140);
    const auto *shk_141 = buffer.data(shk + 141);
    const auto *shk_142 = buffer.data(shk + 142);
    const auto *shk_143 = buffer.data(shk + 143);
    const auto *shk_144 = buffer.data(shk + 144);
    const auto *shk_146 = buffer.data(shk + 146);
    const auto *shk_147 = buffer.data(shk + 147);
    const auto *shk_149 = buffer.data(shk + 149);
    const auto *shk_150 = buffer.data(shk + 150);
    const auto *shk_153 = buffer.data(shk + 153);
    const auto *shk_154 = buffer.data(shk + 154);
    const auto *shk_158 = buffer.data(shk + 158);
    const auto *shk_159 = buffer.data(shk + 159);
    const auto *shk_164 = buffer.data(shk + 164);
    const auto *shk_172 = buffer.data(shk + 172);
    const auto *shk_173 = buffer.data(shk + 173);
    const auto *shk_174 = buffer.data(shk + 174);
    const auto *shk_175 = buffer.data(shk + 175);
    const auto *shk_176 = buffer.data(shk + 176);
    const auto *shk_177 = buffer.data(shk + 177);
    const auto *shk_178 = buffer.data(shk + 178);
    const auto *shk_179 = buffer.data(shk + 179);
    const auto *shk_180 = buffer.data(shk + 180);
    const auto *shk_182 = buffer.data(shk + 182);
    const auto *shk_183 = buffer.data(shk + 183);
    const auto *shk_185 = buffer.data(shk + 185);
    const auto *shk_186 = buffer.data(shk + 186);

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pc_x, sgk_100, sgk_101, sgk_102, \
                         sgk_103, sgk_104, shk_100, shk_101, shk_102, shk_103, \
                         shk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_18 * sgk_100[k]
                   + f_3 * pc_x[k] * shk_100[k];

        t_119[k] = f_18 * sgk_101[k]
                   + f_3 * pc_x[k] * shk_101[k];

        t_120[k] = f_18 * sgk_102[k]
                   + f_3 * pc_x[k] * shk_102[k];

        t_121[k] = f_18 * sgk_103[k]
                   + f_3 * pc_x[k] * shk_103[k];

        t_122[k] = f_18 * sgk_104[k]
                   + f_3 * pc_x[k] * shk_104[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pb_z, pc_x, pc_z, sgl0_36, sgk_105, \
                         sgk_106, sgk_107, sgl1_36, shk_105, shk_106, \
                         shk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_18 * sgk_105[k]
                   + f_3 * pc_x[k] * shk_105[k];

        t_124[k] = f_18 * sgk_106[k]
                   + f_3 * pc_x[k] * shk_106[k];

        t_125[k] = f_18 * sgk_107[k]
                   + f_3 * pc_x[k] * shk_107[k];

        t_126[k] = pb_z[k] * sgl0_36[k]
                   - f_14 * pc_z[k] * sgl1_36[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pc_y, pc_z, sgk_28, shi0_79, shi0_80, shi1_79, \
                         shi1_80, shk_100, shk_102, shk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_15 * sgk_28[k]
                   + f_3 * pc_z[k] * shk_100[k];

        t_128[k] = f_4 * shi0_79[k]
                   - f_5 * shi1_79[k]
                   + f_3 * pc_y[k] * shk_102[k];

        t_129[k] = f_6 * shi0_80[k]
                   - f_7 * shi1_80[k]
                   + f_3 * pc_y[k] * shk_103[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_y, shi0_81, shi0_82, shi0_83, shi1_81, \
                         shi1_82, shi1_83, shk_104, shk_105, shk_106, \
                         shk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_8 * shi0_81[k]
                   - f_9 * shi1_81[k]
                   + f_3 * pc_y[k] * shk_104[k];

        t_131[k] = f_10 * shi0_82[k]
                   - f_11 * shi1_82[k]
                   + f_3 * pc_y[k] * shk_105[k];

        t_132[k] = f_12 * shi0_83[k]
                   - f_13 * shi1_83[k]
                   + f_3 * pc_y[k] * shk_106[k];

        t_133[k] = f_3 * pc_y[k] * shk_107[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, sgk_35, sgk_36, \
                         sgk_108, shi0_83, shi0_84, shi1_83, shi1_84, shk_107, \
                         shk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_15 * sgk_35[k]
                   + f_1 * shi0_83[k]
                   - f_2 * shi1_83[k]
                   + f_3 * pc_z[k] * shk_107[k];

        t_135[k] = f_17 * sgk_108[k]
                   + f_1 * shi0_84[k]
                   - f_2 * shi1_84[k]
                   + f_3 * pc_x[k] * shk_108[k];

        t_136[k] = f_16 * sgk_36[k]
                   + f_3 * pc_y[k] * shk_108[k];

        t_137[k] = f_3 * pc_z[k] * shk_108[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_x, pc_y, sgk_38, sgk_111, sgk_113, shi0_87, \
                         shi0_89, shi1_87, shi1_89, shk_110, shk_111, \
                         shk_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_17 * sgk_111[k]
                   + f_4 * shi0_87[k]
                   - f_5 * shi1_87[k]
                   + f_3 * pc_x[k] * shk_111[k];

        t_139[k] = f_16 * sgk_38[k]
                   + f_3 * pc_y[k] * shk_110[k];

        t_140[k] = f_17 * sgk_113[k]
                   + f_4 * shi0_89[k]
                   - f_5 * shi1_89[k]
                   + f_3 * pc_x[k] * shk_113[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_x, pc_y, pc_z, sgk_41, sgk_114, shi0_90, \
                         shi1_90, shk_111, shk_113, shk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_17 * sgk_114[k]
                   + f_6 * shi0_90[k]
                   - f_7 * shi1_90[k]
                   + f_3 * pc_x[k] * shk_114[k];

        t_142[k] = f_3 * pc_z[k] * shk_111[k];

        t_143[k] = f_16 * sgk_41[k]
                   + f_3 * pc_y[k] * shk_113[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_z, sgk_117, sgk_118, shi0_93, shi0_94, \
                         shi1_93, shi1_94, shk_114, shk_117, shk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_17 * sgk_117[k]
                   + f_6 * shi0_93[k]
                   - f_7 * shi1_93[k]
                   + f_3 * pc_x[k] * shk_117[k];

        t_145[k] = f_17 * sgk_118[k]
                   + f_8 * shi0_94[k]
                   - f_9 * shi1_94[k]
                   + f_3 * pc_x[k] * shk_118[k];

        t_146[k] = f_3 * pc_z[k] * shk_114[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_x, pc_y, sgk_45, sgk_120, sgk_122, shi0_96, \
                         shi0_98, shi1_96, shi1_98, shk_117, shk_120, \
                         shk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_17 * sgk_120[k]
                   + f_8 * shi0_96[k]
                   - f_9 * shi1_96[k]
                   + f_3 * pc_x[k] * shk_120[k];

        t_148[k] = f_16 * sgk_45[k]
                   + f_3 * pc_y[k] * shk_117[k];

        t_149[k] = f_17 * sgk_122[k]
                   + f_8 * shi0_98[k]
                   - f_9 * shi1_98[k]
                   + f_3 * pc_x[k] * shk_122[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pc_x, pc_z, sgk_123, sgk_125, shi0_99, shi0_101, \
                         shi1_99, shi1_101, shk_118, shk_123, shk_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_17 * sgk_123[k]
                   + f_10 * shi0_99[k]
                   - f_11 * shi1_99[k]
                   + f_3 * pc_x[k] * shk_123[k];

        t_151[k] = f_3 * pc_z[k] * shk_118[k];

        t_152[k] = f_17 * sgk_125[k]
                   + f_10 * shi0_101[k]
                   - f_11 * shi1_101[k]
                   + f_3 * pc_x[k] * shk_125[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_x, pc_y, sgk_50, sgk_126, sgk_128, shi0_102, \
                         shi0_104, shi1_102, shi1_104, shk_122, shk_126, \
                         shk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_17 * sgk_126[k]
                   + f_10 * shi0_102[k]
                   - f_11 * shi1_102[k]
                   + f_3 * pc_x[k] * shk_126[k];

        t_154[k] = f_16 * sgk_50[k]
                   + f_3 * pc_y[k] * shk_122[k];

        t_155[k] = f_17 * sgk_128[k]
                   + f_10 * shi0_104[k]
                   - f_11 * shi1_104[k]
                   + f_3 * pc_x[k] * shk_128[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pc_x, pc_z, sgk_129, sgk_131, shi0_105, \
                         shi0_107, shi1_105, shi1_107, shk_123, shk_129, \
                         shk_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_17 * sgk_129[k]
                   + f_12 * shi0_105[k]
                   - f_13 * shi1_105[k]
                   + f_3 * pc_x[k] * shk_129[k];

        t_157[k] = f_3 * pc_z[k] * shk_123[k];

        t_158[k] = f_17 * sgk_131[k]
                   + f_12 * shi0_107[k]
                   - f_13 * shi1_107[k]
                   + f_3 * pc_x[k] * shk_131[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pc_x, pc_y, sgk_56, sgk_132, sgk_133, shi0_108, \
                         shi0_109, shi1_108, shi1_109, shk_128, shk_132, \
                         shk_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_17 * sgk_132[k]
                   + f_12 * shi0_108[k]
                   - f_13 * shi1_108[k]
                   + f_3 * pc_x[k] * shk_132[k];

        t_160[k] = f_17 * sgk_133[k]
                   + f_12 * shi0_109[k]
                   - f_13 * shi1_109[k]
                   + f_3 * pc_x[k] * shk_133[k];

        t_161[k] = f_16 * sgk_56[k]
                   + f_3 * pc_y[k] * shk_128[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, sgk_135, sgk_136, sgk_137, sgk_138, \
                         shi0_111, shi1_111, shk_135, shk_136, shk_137, \
                         shk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_17 * sgk_135[k]
                   + f_12 * shi0_111[k]
                   - f_13 * shi1_111[k]
                   + f_3 * pc_x[k] * shk_135[k];

        t_163[k] = f_17 * sgk_136[k]
                   + f_3 * pc_x[k] * shk_136[k];

        t_164[k] = f_17 * sgk_137[k]
                   + f_3 * pc_x[k] * shk_137[k];

        t_165[k] = f_17 * sgk_138[k]
                   + f_3 * pc_x[k] * shk_138[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, sgk_139, sgk_140, sgk_141, \
                         sgk_142, sgk_143, shk_139, shk_140, shk_141, shk_142, \
                         shk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_17 * sgk_139[k]
                   + f_3 * pc_x[k] * shk_139[k];

        t_167[k] = f_17 * sgk_140[k]
                   + f_3 * pc_x[k] * shk_140[k];

        t_168[k] = f_17 * sgk_141[k]
                   + f_3 * pc_x[k] * shk_141[k];

        t_169[k] = f_17 * sgk_142[k]
                   + f_3 * pc_x[k] * shk_142[k];

        t_170[k] = f_17 * sgk_143[k]
                   + f_3 * pc_x[k] * shk_143[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pc_y, pc_z, sgk_64, sgk_66, shi0_105, shi0_107, \
                         shi1_105, shi1_107, shk_136, shk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_16 * sgk_64[k]
                   + f_1 * shi0_105[k]
                   - f_2 * shi1_105[k]
                   + f_3 * pc_y[k] * shk_136[k];

        t_172[k] = f_3 * pc_z[k] * shk_136[k];

        t_173[k] = f_16 * sgk_66[k]
                   + f_4 * shi0_107[k]
                   - f_5 * shi1_107[k]
                   + f_3 * pc_y[k] * shk_138[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, sgk_67, sgk_68, sgk_69, shi0_108, \
                         shi0_109, shi0_110, shi1_108, shi1_109, shi1_110, shk_139, shk_140, \
                         shk_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_16 * sgk_67[k]
                   + f_6 * shi0_108[k]
                   - f_7 * shi1_108[k]
                   + f_3 * pc_y[k] * shk_139[k];

        t_175[k] = f_16 * sgk_68[k]
                   + f_8 * shi0_109[k]
                   - f_9 * shi1_109[k]
                   + f_3 * pc_y[k] * shk_140[k];

        t_176[k] = f_16 * sgk_69[k]
                   + f_10 * shi0_110[k]
                   - f_11 * shi1_110[k]
                   + f_3 * pc_y[k] * shk_141[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pb_y, pc_y, pc_z, sgl0_90, sgk_70, \
                         sgk_71, sgl1_90, shi0_111, shi1_111, shk_142, \
                         shk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_16 * sgk_70[k]
                   + f_12 * shi0_111[k]
                   - f_13 * shi1_111[k]
                   + f_3 * pc_y[k] * shk_142[k];

        t_178[k] = f_16 * sgk_71[k]
                   + f_3 * pc_y[k] * shk_143[k];

        t_179[k] = f_1 * shi0_111[k]
                   - f_2 * shi1_111[k]
                   + f_3 * pc_z[k] * shk_143[k];

        t_180[k] = pb_y[k] * sgl0_90[k]
                   - f_14 * pc_y[k] * sgl1_90[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_z, pc_y, pc_z, sgl0_48, sgk_36, \
                         sgk_72, sgk_74, sgl1_48, shk_144, shk_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_15 * sgk_72[k]
                   + f_3 * pc_y[k] * shk_144[k];

        t_182[k] = f_15 * sgk_36[k]
                   + f_3 * pc_z[k] * shk_144[k];

        t_183[k] = pb_z[k] * sgl0_48[k]
                   - f_14 * pc_z[k] * sgl1_48[k];

        t_184[k] = f_15 * sgk_74[k]
                   + f_3 * pc_y[k] * shk_146[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_y, pb_z, pc_y, pc_z, sgl0_51, sgl0_95, \
                         sgk_39, sgk_77, sgl1_51, sgl1_95, shk_147, \
                         shk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pb_y[k] * sgl0_95[k]
                   - f_14 * pc_y[k] * sgl1_95[k];

        t_186[k] = pb_z[k] * sgl0_51[k]
                   - f_14 * pc_z[k] * sgl1_51[k];

        t_187[k] = f_15 * sgk_39[k]
                   + f_3 * pc_z[k] * shk_147[k];

        t_188[k] = f_15 * sgk_77[k]
                   + f_3 * pc_y[k] * shk_149[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pb_y, pb_z, pc_y, pc_z, sgl0_55, sgl0_99, \
                         sgk_42, sgl1_55, sgl1_99, shk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pb_y[k] * sgl0_99[k]
                   - f_14 * pc_y[k] * sgl1_99[k];

        t_190[k] = pb_z[k] * sgl0_55[k]
                   - f_14 * pc_z[k] * sgl1_55[k];

        t_191[k] = f_15 * sgk_42[k]
                   + f_3 * pc_z[k] * shk_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_y, pc_y, sgl0_102, sgl0_104, sgk_80, sgk_81, \
                         sgl1_102, sgl1_104, shk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pb_y[k] * sgl0_102[k]
                   + f_16 * sgk_80[k]
                   - f_14 * pc_y[k] * sgl1_102[k];

        t_193[k] = f_15 * sgk_81[k]
                   + f_3 * pc_y[k] * shk_153[k];

        t_194[k] = pb_y[k] * sgl0_104[k]
                   - f_14 * pc_y[k] * sgl1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pb_y, pb_z, pc_y, pc_z, sgl0_60, sgl0_107, \
                         sgk_46, sgk_84, sgl1_60, sgl1_107, shk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pb_z[k] * sgl0_60[k]
                   - f_14 * pc_z[k] * sgl1_60[k];

        t_196[k] = f_15 * sgk_46[k]
                   + f_3 * pc_z[k] * shk_154[k];

        t_197[k] = pb_y[k] * sgl0_107[k]
                   + f_17 * sgk_84[k]
                   - f_14 * pc_y[k] * sgl1_107[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pb_y, pc_y, sgl0_108, sgl0_110, sgk_85, sgk_86, \
                         sgl1_108, sgl1_110, shk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pb_y[k] * sgl0_108[k]
                   + f_16 * sgk_85[k]
                   - f_14 * pc_y[k] * sgl1_108[k];

        t_199[k] = f_15 * sgk_86[k]
                   + f_3 * pc_y[k] * shk_158[k];

        t_200[k] = pb_y[k] * sgl0_110[k]
                   - f_14 * pc_y[k] * sgl1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pb_y, pb_z, pc_y, pc_z, sgl0_66, sgl0_113, \
                         sgk_51, sgk_89, sgl1_66, sgl1_113, shk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pb_z[k] * sgl0_66[k]
                   - f_14 * pc_z[k] * sgl1_66[k];

        t_202[k] = f_15 * sgk_51[k]
                   + f_3 * pc_z[k] * shk_159[k];

        t_203[k] = pb_y[k] * sgl0_113[k]
                   + f_18 * sgk_89[k]
                   - f_14 * pc_y[k] * sgl1_113[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_y, pc_y, sgl0_114, sgl0_115, sgl0_117, \
                         sgk_90, sgk_91, sgk_92, sgl1_114, sgl1_115, sgl1_117, \
                         shk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pb_y[k] * sgl0_114[k]
                   + f_17 * sgk_90[k]
                   - f_14 * pc_y[k] * sgl1_114[k];

        t_205[k] = pb_y[k] * sgl0_115[k]
                   + f_16 * sgk_91[k]
                   - f_14 * pc_y[k] * sgl1_115[k];

        t_206[k] = f_15 * sgk_92[k]
                   + f_3 * pc_y[k] * shk_164[k];

        t_207[k] = pb_y[k] * sgl0_117[k]
                   - f_14 * pc_y[k] * sgl1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pc_x, sgk_172, sgk_173, sgk_174, \
                         sgk_175, sgk_176, shk_172, shk_173, shk_174, shk_175, \
                         shk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_17 * sgk_172[k]
                   + f_3 * pc_x[k] * shk_172[k];

        t_209[k] = f_17 * sgk_173[k]
                   + f_3 * pc_x[k] * shk_173[k];

        t_210[k] = f_17 * sgk_174[k]
                   + f_3 * pc_x[k] * shk_174[k];

        t_211[k] = f_17 * sgk_175[k]
                   + f_3 * pc_x[k] * shk_175[k];

        t_212[k] = f_17 * sgk_176[k]
                   + f_3 * pc_x[k] * shk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pb_z, pc_x, pc_z, sgl0_81, sgk_177, \
                         sgk_178, sgk_179, sgl1_81, shk_177, shk_178, \
                         shk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_17 * sgk_177[k]
                   + f_3 * pc_x[k] * shk_177[k];

        t_214[k] = f_17 * sgk_178[k]
                   + f_3 * pc_x[k] * shk_178[k];

        t_215[k] = f_17 * sgk_179[k]
                   + f_3 * pc_x[k] * shk_179[k];

        t_216[k] = pb_z[k] * sgl0_81[k]
                   - f_14 * pc_z[k] * sgl1_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pc_y, pc_z, sgk_64, sgk_102, sgk_103, shi0_135, \
                         shi0_136, shi1_135, shi1_136, shk_172, shk_174, \
                         shk_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * sgk_64[k]
                   + f_3 * pc_z[k] * shk_172[k];

        t_218[k] = f_15 * sgk_102[k]
                   + f_4 * shi0_135[k]
                   - f_5 * shi1_135[k]
                   + f_3 * pc_y[k] * shk_174[k];

        t_219[k] = f_15 * sgk_103[k]
                   + f_6 * shi0_136[k]
                   - f_7 * shi1_136[k]
                   + f_3 * pc_y[k] * shk_175[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_y, sgk_104, sgk_105, sgk_106, shi0_137, \
                         shi0_138, shi0_139, shi1_137, shi1_138, shi1_139, shk_176, shk_177, \
                         shk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * sgk_104[k]
                   + f_8 * shi0_137[k]
                   - f_9 * shi1_137[k]
                   + f_3 * pc_y[k] * shk_176[k];

        t_221[k] = f_15 * sgk_105[k]
                   + f_10 * shi0_138[k]
                   - f_11 * shi1_138[k]
                   + f_3 * pc_y[k] * shk_177[k];

        t_222[k] = f_15 * sgk_106[k]
                   + f_12 * shi0_139[k]
                   - f_13 * shi1_139[k]
                   + f_3 * pc_y[k] * shk_178[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pb_y, pc_x, pc_y, sgl0_134, sgk_107, \
                         sgk_180, sgl1_134, shi0_140, shi1_140, shk_179, \
                         shk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * sgk_107[k]
                   + f_3 * pc_y[k] * shk_179[k];

        t_224[k] = pb_y[k] * sgl0_134[k]
                   - f_14 * pc_y[k] * sgl1_134[k];

        t_225[k] = f_17 * sgk_180[k]
                   + f_1 * shi0_140[k]
                   - f_2 * shi1_140[k]
                   + f_3 * pc_x[k] * shk_180[k];

        t_226[k] = f_3 * pc_y[k] * shk_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_x, pc_y, pc_z, sgk_72, sgk_183, shi0_143, \
                         shi1_143, shk_180, shk_182, shk_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_16 * sgk_72[k]
                   + f_3 * pc_z[k] * shk_180[k];

        t_228[k] = f_17 * sgk_183[k]
                   + f_4 * shi0_143[k]
                   - f_5 * shi1_143[k]
                   + f_3 * pc_x[k] * shk_183[k];

        t_229[k] = f_3 * pc_y[k] * shk_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, pc_x, pc_z, sgk_75, sgk_185, sgk_186, shi0_145, \
                         shi0_146, shi1_145, shi1_146, shk_183, shk_185, \
                         shk_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_17 * sgk_185[k]
                   + f_4 * shi0_145[k]
                   - f_5 * shi1_145[k]
                   + f_3 * pc_x[k] * shk_185[k];

        t_231[k] = f_17 * sgk_186[k]
                   + f_6 * shi0_146[k]
                   - f_7 * shi1_146[k]
                   + f_3 * pc_x[k] * shk_186[k];

        t_232[k] = f_16 * sgk_75[k]
                   + f_3 * pc_z[k] * shk_183[k];
    }
}

static auto
compute_prim_shl_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgl0,
                                                          const size_t sgk, const size_t sgl1,
                                                          const size_t shi0, const size_t shi1,
                                                          const size_t shk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_18 = 2.0 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgl0_135 = buffer.data(sgl0 + 135);
    const auto *sgl0_138 = buffer.data(sgl0 + 138);
    const auto *sgl0_141 = buffer.data(sgl0 + 141);
    const auto *sgl0_145 = buffer.data(sgl0 + 145);
    const auto *sgl0_147 = buffer.data(sgl0 + 147);
    const auto *sgl0_150 = buffer.data(sgl0 + 150);
    const auto *sgl0_152 = buffer.data(sgl0 + 152);
    const auto *sgl0_153 = buffer.data(sgl0 + 153);
    const auto *sgl0_156 = buffer.data(sgl0 + 156);
    const auto *sgl0_158 = buffer.data(sgl0 + 158);
    const auto *sgl0_159 = buffer.data(sgl0 + 159);
    const auto *sgl0_160 = buffer.data(sgl0 + 160);

    const auto *sgk_78 = buffer.data(sgk + 78);
    const auto *sgk_82 = buffer.data(sgk + 82);
    const auto *sgk_87 = buffer.data(sgk + 87);
    const auto *sgk_100 = buffer.data(sgk + 100);
    const auto *sgk_107 = buffer.data(sgk + 107);
    const auto *sgk_108 = buffer.data(sgk + 108);
    const auto *sgk_110 = buffer.data(sgk + 110);
    const auto *sgk_111 = buffer.data(sgk + 111);
    const auto *sgk_113 = buffer.data(sgk + 113);
    const auto *sgk_114 = buffer.data(sgk + 114);
    const auto *sgk_115 = buffer.data(sgk + 115);
    const auto *sgk_117 = buffer.data(sgk + 117);
    const auto *sgk_118 = buffer.data(sgk + 118);
    const auto *sgk_119 = buffer.data(sgk + 119);
    const auto *sgk_120 = buffer.data(sgk + 120);
    const auto *sgk_122 = buffer.data(sgk + 122);
    const auto *sgk_123 = buffer.data(sgk + 123);
    const auto *sgk_124 = buffer.data(sgk + 124);
    const auto *sgk_125 = buffer.data(sgk + 125);
    const auto *sgk_126 = buffer.data(sgk + 126);
    const auto *sgk_128 = buffer.data(sgk + 128);
    const auto *sgk_136 = buffer.data(sgk + 136);
    const auto *sgk_138 = buffer.data(sgk + 138);
    const auto *sgk_139 = buffer.data(sgk + 139);
    const auto *sgk_140 = buffer.data(sgk + 140);
    const auto *sgk_141 = buffer.data(sgk + 141);
    const auto *sgk_142 = buffer.data(sgk + 142);
    const auto *sgk_143 = buffer.data(sgk + 143);
    const auto *sgk_144 = buffer.data(sgk + 144);
    const auto *sgk_146 = buffer.data(sgk + 146);
    const auto *sgk_149 = buffer.data(sgk + 149);
    const auto *sgk_153 = buffer.data(sgk + 153);
    const auto *sgk_158 = buffer.data(sgk + 158);
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
    const auto *sgk_219 = buffer.data(sgk + 219);
    const auto *sgk_221 = buffer.data(sgk + 221);
    const auto *sgk_222 = buffer.data(sgk + 222);
    const auto *sgk_225 = buffer.data(sgk + 225);
    const auto *sgk_226 = buffer.data(sgk + 226);
    const auto *sgk_228 = buffer.data(sgk + 228);
    const auto *sgk_230 = buffer.data(sgk + 230);
    const auto *sgk_231 = buffer.data(sgk + 231);
    const auto *sgk_233 = buffer.data(sgk + 233);
    const auto *sgk_234 = buffer.data(sgk + 234);
    const auto *sgk_236 = buffer.data(sgk + 236);
    const auto *sgk_237 = buffer.data(sgk + 237);
    const auto *sgk_239 = buffer.data(sgk + 239);
    const auto *sgk_240 = buffer.data(sgk + 240);
    const auto *sgk_241 = buffer.data(sgk + 241);
    const auto *sgk_243 = buffer.data(sgk + 243);
    const auto *sgk_244 = buffer.data(sgk + 244);
    const auto *sgk_245 = buffer.data(sgk + 245);
    const auto *sgk_246 = buffer.data(sgk + 246);
    const auto *sgk_247 = buffer.data(sgk + 247);
    const auto *sgk_248 = buffer.data(sgk + 248);
    const auto *sgk_249 = buffer.data(sgk + 249);
    const auto *sgk_250 = buffer.data(sgk + 250);
    const auto *sgk_251 = buffer.data(sgk + 251);
    const auto *sgk_257 = buffer.data(sgk + 257);
    const auto *sgk_261 = buffer.data(sgk + 261);
    const auto *sgk_266 = buffer.data(sgk + 266);
    const auto *sgk_272 = buffer.data(sgk + 272);

    const auto *sgl1_135 = buffer.data(sgl1 + 135);
    const auto *sgl1_138 = buffer.data(sgl1 + 138);
    const auto *sgl1_141 = buffer.data(sgl1 + 141);
    const auto *sgl1_145 = buffer.data(sgl1 + 145);
    const auto *sgl1_147 = buffer.data(sgl1 + 147);
    const auto *sgl1_150 = buffer.data(sgl1 + 150);
    const auto *sgl1_152 = buffer.data(sgl1 + 152);
    const auto *sgl1_153 = buffer.data(sgl1 + 153);
    const auto *sgl1_156 = buffer.data(sgl1 + 156);
    const auto *sgl1_158 = buffer.data(sgl1 + 158);
    const auto *sgl1_159 = buffer.data(sgl1 + 159);
    const auto *sgl1_160 = buffer.data(sgl1 + 160);

    const auto *shi0_149 = buffer.data(shi0 + 149);
    const auto *shi0_150 = buffer.data(shi0 + 150);
    const auto *shi0_152 = buffer.data(shi0 + 152);
    const auto *shi0_154 = buffer.data(shi0 + 154);
    const auto *shi0_155 = buffer.data(shi0 + 155);
    const auto *shi0_157 = buffer.data(shi0 + 157);
    const auto *shi0_158 = buffer.data(shi0 + 158);
    const auto *shi0_160 = buffer.data(shi0 + 160);
    const auto *shi0_161 = buffer.data(shi0 + 161);
    const auto *shi0_163 = buffer.data(shi0 + 163);
    const auto *shi0_164 = buffer.data(shi0 + 164);
    const auto *shi0_165 = buffer.data(shi0 + 165);
    const auto *shi0_166 = buffer.data(shi0 + 166);
    const auto *shi0_167 = buffer.data(shi0 + 167);
    const auto *shi0_168 = buffer.data(shi0 + 168);
    const auto *shi0_171 = buffer.data(shi0 + 171);
    const auto *shi0_173 = buffer.data(shi0 + 173);
    const auto *shi0_174 = buffer.data(shi0 + 174);
    const auto *shi0_177 = buffer.data(shi0 + 177);
    const auto *shi0_178 = buffer.data(shi0 + 178);
    const auto *shi0_180 = buffer.data(shi0 + 180);
    const auto *shi0_182 = buffer.data(shi0 + 182);
    const auto *shi0_183 = buffer.data(shi0 + 183);
    const auto *shi0_185 = buffer.data(shi0 + 185);
    const auto *shi0_186 = buffer.data(shi0 + 186);
    const auto *shi0_188 = buffer.data(shi0 + 188);
    const auto *shi0_189 = buffer.data(shi0 + 189);
    const auto *shi0_191 = buffer.data(shi0 + 191);
    const auto *shi0_192 = buffer.data(shi0 + 192);
    const auto *shi0_193 = buffer.data(shi0 + 193);
    const auto *shi0_194 = buffer.data(shi0 + 194);
    const auto *shi0_195 = buffer.data(shi0 + 195);
    const auto *shi0_201 = buffer.data(shi0 + 201);
    const auto *shi0_205 = buffer.data(shi0 + 205);
    const auto *shi0_210 = buffer.data(shi0 + 210);
    const auto *shi0_216 = buffer.data(shi0 + 216);

    const auto *shi1_149 = buffer.data(shi1 + 149);
    const auto *shi1_150 = buffer.data(shi1 + 150);
    const auto *shi1_152 = buffer.data(shi1 + 152);
    const auto *shi1_154 = buffer.data(shi1 + 154);
    const auto *shi1_155 = buffer.data(shi1 + 155);
    const auto *shi1_157 = buffer.data(shi1 + 157);
    const auto *shi1_158 = buffer.data(shi1 + 158);
    const auto *shi1_160 = buffer.data(shi1 + 160);
    const auto *shi1_161 = buffer.data(shi1 + 161);
    const auto *shi1_163 = buffer.data(shi1 + 163);
    const auto *shi1_164 = buffer.data(shi1 + 164);
    const auto *shi1_165 = buffer.data(shi1 + 165);
    const auto *shi1_166 = buffer.data(shi1 + 166);
    const auto *shi1_167 = buffer.data(shi1 + 167);
    const auto *shi1_168 = buffer.data(shi1 + 168);
    const auto *shi1_171 = buffer.data(shi1 + 171);
    const auto *shi1_173 = buffer.data(shi1 + 173);
    const auto *shi1_174 = buffer.data(shi1 + 174);
    const auto *shi1_177 = buffer.data(shi1 + 177);
    const auto *shi1_178 = buffer.data(shi1 + 178);
    const auto *shi1_180 = buffer.data(shi1 + 180);
    const auto *shi1_182 = buffer.data(shi1 + 182);
    const auto *shi1_183 = buffer.data(shi1 + 183);
    const auto *shi1_185 = buffer.data(shi1 + 185);
    const auto *shi1_186 = buffer.data(shi1 + 186);
    const auto *shi1_188 = buffer.data(shi1 + 188);
    const auto *shi1_189 = buffer.data(shi1 + 189);
    const auto *shi1_191 = buffer.data(shi1 + 191);
    const auto *shi1_192 = buffer.data(shi1 + 192);
    const auto *shi1_193 = buffer.data(shi1 + 193);
    const auto *shi1_194 = buffer.data(shi1 + 194);
    const auto *shi1_195 = buffer.data(shi1 + 195);
    const auto *shi1_201 = buffer.data(shi1 + 201);
    const auto *shi1_205 = buffer.data(shi1 + 205);
    const auto *shi1_210 = buffer.data(shi1 + 210);
    const auto *shi1_216 = buffer.data(shi1 + 216);

    const auto *shk_185 = buffer.data(shk + 185);
    const auto *shk_186 = buffer.data(shk + 186);
    const auto *shk_189 = buffer.data(shk + 189);
    const auto *shk_190 = buffer.data(shk + 190);
    const auto *shk_192 = buffer.data(shk + 192);
    const auto *shk_194 = buffer.data(shk + 194);
    const auto *shk_195 = buffer.data(shk + 195);
    const auto *shk_197 = buffer.data(shk + 197);
    const auto *shk_198 = buffer.data(shk + 198);
    const auto *shk_200 = buffer.data(shk + 200);
    const auto *shk_201 = buffer.data(shk + 201);
    const auto *shk_203 = buffer.data(shk + 203);
    const auto *shk_204 = buffer.data(shk + 204);
    const auto *shk_205 = buffer.data(shk + 205);
    const auto *shk_207 = buffer.data(shk + 207);
    const auto *shk_208 = buffer.data(shk + 208);
    const auto *shk_209 = buffer.data(shk + 209);
    const auto *shk_210 = buffer.data(shk + 210);
    const auto *shk_211 = buffer.data(shk + 211);
    const auto *shk_212 = buffer.data(shk + 212);
    const auto *shk_213 = buffer.data(shk + 213);
    const auto *shk_214 = buffer.data(shk + 214);
    const auto *shk_215 = buffer.data(shk + 215);
    const auto *shk_216 = buffer.data(shk + 216);
    const auto *shk_218 = buffer.data(shk + 218);
    const auto *shk_219 = buffer.data(shk + 219);
    const auto *shk_221 = buffer.data(shk + 221);
    const auto *shk_222 = buffer.data(shk + 222);
    const auto *shk_225 = buffer.data(shk + 225);
    const auto *shk_226 = buffer.data(shk + 226);
    const auto *shk_228 = buffer.data(shk + 228);
    const auto *shk_230 = buffer.data(shk + 230);
    const auto *shk_231 = buffer.data(shk + 231);
    const auto *shk_233 = buffer.data(shk + 233);
    const auto *shk_234 = buffer.data(shk + 234);
    const auto *shk_236 = buffer.data(shk + 236);
    const auto *shk_237 = buffer.data(shk + 237);
    const auto *shk_239 = buffer.data(shk + 239);
    const auto *shk_240 = buffer.data(shk + 240);
    const auto *shk_241 = buffer.data(shk + 241);
    const auto *shk_243 = buffer.data(shk + 243);
    const auto *shk_244 = buffer.data(shk + 244);
    const auto *shk_245 = buffer.data(shk + 245);
    const auto *shk_246 = buffer.data(shk + 246);
    const auto *shk_247 = buffer.data(shk + 247);
    const auto *shk_248 = buffer.data(shk + 248);
    const auto *shk_249 = buffer.data(shk + 249);
    const auto *shk_250 = buffer.data(shk + 250);
    const auto *shk_251 = buffer.data(shk + 251);
    const auto *shk_252 = buffer.data(shk + 252);
    const auto *shk_254 = buffer.data(shk + 254);
    const auto *shk_255 = buffer.data(shk + 255);
    const auto *shk_257 = buffer.data(shk + 257);
    const auto *shk_258 = buffer.data(shk + 258);
    const auto *shk_261 = buffer.data(shk + 261);
    const auto *shk_262 = buffer.data(shk + 262);
    const auto *shk_266 = buffer.data(shk + 266);
    const auto *shk_267 = buffer.data(shk + 267);
    const auto *shk_272 = buffer.data(shk + 272);

#pragma omp simd aligned(t_233, t_234, t_235, pc_x, pc_y, sgk_189, sgk_190, shi0_149, \
                         shi0_150, shi1_149, shi1_150, shk_185, shk_189, \
                         shk_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * pc_y[k] * shk_185[k];

        t_234[k] = f_17 * sgk_189[k]
                   + f_6 * shi0_149[k]
                   - f_7 * shi1_149[k]
                   + f_3 * pc_x[k] * shk_189[k];

        t_235[k] = f_17 * sgk_190[k]
                   + f_8 * shi0_150[k]
                   - f_9 * shi1_150[k]
                   + f_3 * pc_x[k] * shk_190[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pc_x, pc_y, pc_z, sgk_78, sgk_192, shi0_152, \
                         shi1_152, shk_186, shk_189, shk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_16 * sgk_78[k]
                   + f_3 * pc_z[k] * shk_186[k];

        t_237[k] = f_17 * sgk_192[k]
                   + f_8 * shi0_152[k]
                   - f_9 * shi1_152[k]
                   + f_3 * pc_x[k] * shk_192[k];

        t_238[k] = f_3 * pc_y[k] * shk_189[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, pc_z, sgk_82, sgk_194, sgk_195, shi0_154, \
                         shi0_155, shi1_154, shi1_155, shk_190, shk_194, \
                         shk_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_17 * sgk_194[k]
                   + f_8 * shi0_154[k]
                   - f_9 * shi1_154[k]
                   + f_3 * pc_x[k] * shk_194[k];

        t_240[k] = f_17 * sgk_195[k]
                   + f_10 * shi0_155[k]
                   - f_11 * shi1_155[k]
                   + f_3 * pc_x[k] * shk_195[k];

        t_241[k] = f_16 * sgk_82[k]
                   + f_3 * pc_z[k] * shk_190[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, pc_x, pc_y, sgk_197, sgk_198, shi0_157, \
                         shi0_158, shi1_157, shi1_158, shk_194, shk_197, \
                         shk_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_17 * sgk_197[k]
                   + f_10 * shi0_157[k]
                   - f_11 * shi1_157[k]
                   + f_3 * pc_x[k] * shk_197[k];

        t_243[k] = f_17 * sgk_198[k]
                   + f_10 * shi0_158[k]
                   - f_11 * shi1_158[k]
                   + f_3 * pc_x[k] * shk_198[k];

        t_244[k] = f_3 * pc_y[k] * shk_194[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_x, pc_z, sgk_87, sgk_200, sgk_201, shi0_160, \
                         shi0_161, shi1_160, shi1_161, shk_195, shk_200, \
                         shk_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_17 * sgk_200[k]
                   + f_10 * shi0_160[k]
                   - f_11 * shi1_160[k]
                   + f_3 * pc_x[k] * shk_200[k];

        t_246[k] = f_17 * sgk_201[k]
                   + f_12 * shi0_161[k]
                   - f_13 * shi1_161[k]
                   + f_3 * pc_x[k] * shk_201[k];

        t_247[k] = f_16 * sgk_87[k]
                   + f_3 * pc_z[k] * shk_195[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_x, sgk_203, sgk_204, sgk_205, shi0_163, \
                         shi0_164, shi0_165, shi1_163, shi1_164, shi1_165, shk_203, shk_204, \
                         shk_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_17 * sgk_203[k]
                   + f_12 * shi0_163[k]
                   - f_13 * shi1_163[k]
                   + f_3 * pc_x[k] * shk_203[k];

        t_249[k] = f_17 * sgk_204[k]
                   + f_12 * shi0_164[k]
                   - f_13 * shi1_164[k]
                   + f_3 * pc_x[k] * shk_204[k];

        t_250[k] = f_17 * sgk_205[k]
                   + f_12 * shi0_165[k]
                   - f_13 * shi1_165[k]
                   + f_3 * pc_x[k] * shk_205[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pc_x, pc_y, sgk_207, sgk_208, sgk_209, \
                         shi0_167, shi1_167, shk_200, shk_207, shk_208, \
                         shk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_3 * pc_y[k] * shk_200[k];

        t_252[k] = f_17 * sgk_207[k]
                   + f_12 * shi0_167[k]
                   - f_13 * shi1_167[k]
                   + f_3 * pc_x[k] * shk_207[k];

        t_253[k] = f_17 * sgk_208[k]
                   + f_3 * pc_x[k] * shk_208[k];

        t_254[k] = f_17 * sgk_209[k]
                   + f_3 * pc_x[k] * shk_209[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pc_x, sgk_210, sgk_211, sgk_212, \
                         sgk_213, sgk_214, shk_210, shk_211, shk_212, shk_213, \
                         shk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_17 * sgk_210[k]
                   + f_3 * pc_x[k] * shk_210[k];

        t_256[k] = f_17 * sgk_211[k]
                   + f_3 * pc_x[k] * shk_211[k];

        t_257[k] = f_17 * sgk_212[k]
                   + f_3 * pc_x[k] * shk_212[k];

        t_258[k] = f_17 * sgk_213[k]
                   + f_3 * pc_x[k] * shk_213[k];

        t_259[k] = f_17 * sgk_214[k]
                   + f_3 * pc_x[k] * shk_214[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, pc_z, sgk_100, sgk_215, \
                         shi0_161, shi0_163, shi1_161, shi1_163, shk_208, shk_210, \
                         shk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_17 * sgk_215[k]
                   + f_3 * pc_x[k] * shk_215[k];

        t_261[k] = f_1 * shi0_161[k]
                   - f_2 * shi1_161[k]
                   + f_3 * pc_y[k] * shk_208[k];

        t_262[k] = f_16 * sgk_100[k]
                   + f_3 * pc_z[k] * shk_208[k];

        t_263[k] = f_4 * shi0_163[k]
                   - f_5 * shi1_163[k]
                   + f_3 * pc_y[k] * shk_210[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, shi0_164, shi0_165, shi0_166, shi1_164, \
                         shi1_165, shi1_166, shk_211, shk_212, \
                         shk_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_6 * shi0_164[k]
                   - f_7 * shi1_164[k]
                   + f_3 * pc_y[k] * shk_211[k];

        t_265[k] = f_8 * shi0_165[k]
                   - f_9 * shi1_165[k]
                   + f_3 * pc_y[k] * shk_212[k];

        t_266[k] = f_10 * shi0_166[k]
                   - f_11 * shi1_166[k]
                   + f_3 * pc_y[k] * shk_213[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pc_x, pc_y, pc_z, sgk_107, sgk_216, \
                         shi0_167, shi0_168, shi1_167, shi1_168, shk_214, shk_215, \
                         shk_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_12 * shi0_167[k]
                   - f_13 * shi1_167[k]
                   + f_3 * pc_y[k] * shk_214[k];

        t_268[k] = f_3 * pc_y[k] * shk_215[k];

        t_269[k] = f_16 * sgk_107[k]
                   + f_1 * shi0_167[k]
                   - f_2 * shi1_167[k]
                   + f_3 * pc_z[k] * shk_215[k];

        t_270[k] = f_16 * sgk_216[k]
                   + f_1 * shi0_168[k]
                   - f_2 * shi1_168[k]
                   + f_3 * pc_x[k] * shk_216[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, sgk_108, sgk_110, \
                         sgk_219, shi0_171, shi1_171, shk_216, shk_218, \
                         shk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_17 * sgk_108[k]
                   + f_3 * pc_y[k] * shk_216[k];

        t_272[k] = f_3 * pc_z[k] * shk_216[k];

        t_273[k] = f_16 * sgk_219[k]
                   + f_4 * shi0_171[k]
                   - f_5 * shi1_171[k]
                   + f_3 * pc_x[k] * shk_219[k];

        t_274[k] = f_17 * sgk_110[k]
                   + f_3 * pc_y[k] * shk_218[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pc_x, pc_z, sgk_221, sgk_222, shi0_173, \
                         shi0_174, shi1_173, shi1_174, shk_219, shk_221, \
                         shk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_16 * sgk_221[k]
                   + f_4 * shi0_173[k]
                   - f_5 * shi1_173[k]
                   + f_3 * pc_x[k] * shk_221[k];

        t_276[k] = f_16 * sgk_222[k]
                   + f_6 * shi0_174[k]
                   - f_7 * shi1_174[k]
                   + f_3 * pc_x[k] * shk_222[k];

        t_277[k] = f_3 * pc_z[k] * shk_219[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pc_x, pc_y, sgk_113, sgk_225, sgk_226, shi0_177, \
                         shi0_178, shi1_177, shi1_178, shk_221, shk_225, \
                         shk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_17 * sgk_113[k]
                   + f_3 * pc_y[k] * shk_221[k];

        t_279[k] = f_16 * sgk_225[k]
                   + f_6 * shi0_177[k]
                   - f_7 * shi1_177[k]
                   + f_3 * pc_x[k] * shk_225[k];

        t_280[k] = f_16 * sgk_226[k]
                   + f_8 * shi0_178[k]
                   - f_9 * shi1_178[k]
                   + f_3 * pc_x[k] * shk_226[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, pc_x, pc_y, pc_z, sgk_117, sgk_228, shi0_180, \
                         shi1_180, shk_222, shk_225, shk_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_3 * pc_z[k] * shk_222[k];

        t_282[k] = f_16 * sgk_228[k]
                   + f_8 * shi0_180[k]
                   - f_9 * shi1_180[k]
                   + f_3 * pc_x[k] * shk_228[k];

        t_283[k] = f_17 * sgk_117[k]
                   + f_3 * pc_y[k] * shk_225[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pc_x, pc_z, sgk_230, sgk_231, shi0_182, \
                         shi0_183, shi1_182, shi1_183, shk_226, shk_230, \
                         shk_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_16 * sgk_230[k]
                   + f_8 * shi0_182[k]
                   - f_9 * shi1_182[k]
                   + f_3 * pc_x[k] * shk_230[k];

        t_285[k] = f_16 * sgk_231[k]
                   + f_10 * shi0_183[k]
                   - f_11 * shi1_183[k]
                   + f_3 * pc_x[k] * shk_231[k];

        t_286[k] = f_3 * pc_z[k] * shk_226[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pc_x, pc_y, sgk_122, sgk_233, sgk_234, shi0_185, \
                         shi0_186, shi1_185, shi1_186, shk_230, shk_233, \
                         shk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_16 * sgk_233[k]
                   + f_10 * shi0_185[k]
                   - f_11 * shi1_185[k]
                   + f_3 * pc_x[k] * shk_233[k];

        t_288[k] = f_16 * sgk_234[k]
                   + f_10 * shi0_186[k]
                   - f_11 * shi1_186[k]
                   + f_3 * pc_x[k] * shk_234[k];

        t_289[k] = f_17 * sgk_122[k]
                   + f_3 * pc_y[k] * shk_230[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pc_x, pc_z, sgk_236, sgk_237, shi0_188, \
                         shi0_189, shi1_188, shi1_189, shk_231, shk_236, \
                         shk_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_16 * sgk_236[k]
                   + f_10 * shi0_188[k]
                   - f_11 * shi1_188[k]
                   + f_3 * pc_x[k] * shk_236[k];

        t_291[k] = f_16 * sgk_237[k]
                   + f_12 * shi0_189[k]
                   - f_13 * shi1_189[k]
                   + f_3 * pc_x[k] * shk_237[k];

        t_292[k] = f_3 * pc_z[k] * shk_231[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, pc_x, sgk_239, sgk_240, sgk_241, shi0_191, \
                         shi0_192, shi0_193, shi1_191, shi1_192, shi1_193, shk_239, shk_240, \
                         shk_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_16 * sgk_239[k]
                   + f_12 * shi0_191[k]
                   - f_13 * shi1_191[k]
                   + f_3 * pc_x[k] * shk_239[k];

        t_294[k] = f_16 * sgk_240[k]
                   + f_12 * shi0_192[k]
                   - f_13 * shi1_192[k]
                   + f_3 * pc_x[k] * shk_240[k];

        t_295[k] = f_16 * sgk_241[k]
                   + f_12 * shi0_193[k]
                   - f_13 * shi1_193[k]
                   + f_3 * pc_x[k] * shk_241[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pc_x, pc_y, sgk_128, sgk_243, sgk_244, \
                         sgk_245, shi0_195, shi1_195, shk_236, shk_243, shk_244, \
                         shk_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_17 * sgk_128[k]
                   + f_3 * pc_y[k] * shk_236[k];

        t_297[k] = f_16 * sgk_243[k]
                   + f_12 * shi0_195[k]
                   - f_13 * shi1_195[k]
                   + f_3 * pc_x[k] * shk_243[k];

        t_298[k] = f_16 * sgk_244[k]
                   + f_3 * pc_x[k] * shk_244[k];

        t_299[k] = f_16 * sgk_245[k]
                   + f_3 * pc_x[k] * shk_245[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, pc_x, sgk_246, sgk_247, sgk_248, \
                         sgk_249, sgk_250, shk_246, shk_247, shk_248, shk_249, \
                         shk_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_16 * sgk_246[k]
                   + f_3 * pc_x[k] * shk_246[k];

        t_301[k] = f_16 * sgk_247[k]
                   + f_3 * pc_x[k] * shk_247[k];

        t_302[k] = f_16 * sgk_248[k]
                   + f_3 * pc_x[k] * shk_248[k];

        t_303[k] = f_16 * sgk_249[k]
                   + f_3 * pc_x[k] * shk_249[k];

        t_304[k] = f_16 * sgk_250[k]
                   + f_3 * pc_x[k] * shk_250[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, pc_x, pc_y, pc_z, sgk_136, sgk_251, shi0_189, \
                         shi1_189, shk_244, shk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_16 * sgk_251[k]
                   + f_3 * pc_x[k] * shk_251[k];

        t_306[k] = f_17 * sgk_136[k]
                   + f_1 * shi0_189[k]
                   - f_2 * shi1_189[k]
                   + f_3 * pc_y[k] * shk_244[k];

        t_307[k] = f_3 * pc_z[k] * shk_244[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, pc_y, sgk_138, sgk_139, sgk_140, shi0_191, \
                         shi0_192, shi0_193, shi1_191, shi1_192, shi1_193, shk_246, shk_247, \
                         shk_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_17 * sgk_138[k]
                   + f_4 * shi0_191[k]
                   - f_5 * shi1_191[k]
                   + f_3 * pc_y[k] * shk_246[k];

        t_309[k] = f_17 * sgk_139[k]
                   + f_6 * shi0_192[k]
                   - f_7 * shi1_192[k]
                   + f_3 * pc_y[k] * shk_247[k];

        t_310[k] = f_17 * sgk_140[k]
                   + f_8 * shi0_193[k]
                   - f_9 * shi1_193[k]
                   + f_3 * pc_y[k] * shk_248[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, sgk_141, sgk_142, sgk_143, \
                         shi0_194, shi0_195, shi1_194, shi1_195, shk_249, shk_250, \
                         shk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_17 * sgk_141[k]
                   + f_10 * shi0_194[k]
                   - f_11 * shi1_194[k]
                   + f_3 * pc_y[k] * shk_249[k];

        t_312[k] = f_17 * sgk_142[k]
                   + f_12 * shi0_195[k]
                   - f_13 * shi1_195[k]
                   + f_3 * pc_y[k] * shk_250[k];

        t_313[k] = f_17 * sgk_143[k]
                   + f_3 * pc_y[k] * shk_251[k];

        t_314[k] = f_1 * shi0_195[k]
                   - f_2 * shi1_195[k]
                   + f_3 * pc_z[k] * shk_251[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pb_z, pc_y, pc_z, sgl0_135, sgl0_138, \
                         sgk_108, sgk_144, sgl1_135, sgl1_138, \
                         shk_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pb_z[k] * sgl0_135[k]
                   - f_14 * pc_z[k] * sgl1_135[k];

        t_316[k] = f_16 * sgk_144[k]
                   + f_3 * pc_y[k] * shk_252[k];

        t_317[k] = f_15 * sgk_108[k]
                   + f_3 * pc_z[k] * shk_252[k];

        t_318[k] = pb_z[k] * sgl0_138[k]
                   - f_14 * pc_z[k] * sgl1_138[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pb_z, pc_x, pc_y, pc_z, sgl0_141, sgk_146, \
                         sgk_257, sgl1_141, shi0_201, shi1_201, shk_254, \
                         shk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_16 * sgk_146[k]
                   + f_3 * pc_y[k] * shk_254[k];

        t_320[k] = f_16 * sgk_257[k]
                   + f_4 * shi0_201[k]
                   - f_5 * shi1_201[k]
                   + f_3 * pc_x[k] * shk_257[k];

        t_321[k] = pb_z[k] * sgl0_141[k]
                   - f_14 * pc_z[k] * sgl1_141[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, pc_x, pc_y, pc_z, sgk_111, sgk_149, sgk_261, \
                         shi0_205, shi1_205, shk_255, shk_257, \
                         shk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_15 * sgk_111[k]
                   + f_3 * pc_z[k] * shk_255[k];

        t_323[k] = f_16 * sgk_149[k]
                   + f_3 * pc_y[k] * shk_257[k];

        t_324[k] = f_16 * sgk_261[k]
                   + f_6 * shi0_205[k]
                   - f_7 * shi1_205[k]
                   + f_3 * pc_x[k] * shk_261[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pb_z, pc_y, pc_z, sgl0_145, sgl0_147, \
                         sgk_114, sgk_115, sgk_153, sgl1_145, sgl1_147, shk_258, \
                         shk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pb_z[k] * sgl0_145[k]
                   - f_14 * pc_z[k] * sgl1_145[k];

        t_326[k] = f_15 * sgk_114[k]
                   + f_3 * pc_z[k] * shk_258[k];

        t_327[k] = pb_z[k] * sgl0_147[k]
                   + f_16 * sgk_115[k]
                   - f_14 * pc_z[k] * sgl1_147[k];

        t_328[k] = f_16 * sgk_153[k]
                   + f_3 * pc_y[k] * shk_261[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pb_z, pc_x, pc_z, sgl0_150, sgk_118, sgk_266, \
                         sgl1_150, shi0_210, shi1_210, shk_262, \
                         shk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_16 * sgk_266[k]
                   + f_8 * shi0_210[k]
                   - f_9 * shi1_210[k]
                   + f_3 * pc_x[k] * shk_266[k];

        t_330[k] = pb_z[k] * sgl0_150[k]
                   - f_14 * pc_z[k] * sgl1_150[k];

        t_331[k] = f_15 * sgk_118[k]
                   + f_3 * pc_z[k] * shk_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pb_z, pc_y, pc_z, sgl0_152, sgl0_153, sgk_119, \
                         sgk_120, sgk_158, sgl1_152, sgl1_153, \
                         shk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pb_z[k] * sgl0_152[k]
                   + f_16 * sgk_119[k]
                   - f_14 * pc_z[k] * sgl1_152[k];

        t_333[k] = pb_z[k] * sgl0_153[k]
                   + f_17 * sgk_120[k]
                   - f_14 * pc_z[k] * sgl1_153[k];

        t_334[k] = f_16 * sgk_158[k]
                   + f_3 * pc_y[k] * shk_266[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pb_z, pc_x, pc_z, sgl0_156, sgk_123, sgk_272, \
                         sgl1_156, shi0_216, shi1_216, shk_267, \
                         shk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_16 * sgk_272[k]
                   + f_10 * shi0_216[k]
                   - f_11 * shi1_216[k]
                   + f_3 * pc_x[k] * shk_272[k];

        t_336[k] = pb_z[k] * sgl0_156[k]
                   - f_14 * pc_z[k] * sgl1_156[k];

        t_337[k] = f_15 * sgk_123[k]
                   + f_3 * pc_z[k] * shk_267[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pb_z, pc_z, sgl0_158, sgl0_159, sgl0_160, \
                         sgk_124, sgk_125, sgk_126, sgl1_158, sgl1_159, \
                         sgl1_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pb_z[k] * sgl0_158[k]
                   + f_16 * sgk_124[k]
                   - f_14 * pc_z[k] * sgl1_158[k];

        t_339[k] = pb_z[k] * sgl0_159[k]
                   + f_17 * sgk_125[k]
                   - f_14 * pc_z[k] * sgl1_159[k];

        t_340[k] = pb_z[k] * sgl0_160[k]
                   + f_18 * sgk_126[k]
                   - f_14 * pc_z[k] * sgl1_160[k];
    }
}

static auto
compute_prim_shl_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgl0,
                                                          const size_t sgk, const size_t sgl1,
                                                          const size_t shi0, const size_t shi1,
                                                          const size_t shk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 4.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgl0_171 = buffer.data(sgl0 + 171);
    const auto *sgl0_225 = buffer.data(sgl0 + 225);
    const auto *sgl0_228 = buffer.data(sgl0 + 228);
    const auto *sgl0_230 = buffer.data(sgl0 + 230);
    const auto *sgl0_231 = buffer.data(sgl0 + 231);
    const auto *sgl0_234 = buffer.data(sgl0 + 234);
    const auto *sgl0_235 = buffer.data(sgl0 + 235);
    const auto *sgl0_237 = buffer.data(sgl0 + 237);
    const auto *sgl0_239 = buffer.data(sgl0 + 239);
    const auto *sgl0_240 = buffer.data(sgl0 + 240);
    const auto *sgl0_242 = buffer.data(sgl0 + 242);
    const auto *sgl0_243 = buffer.data(sgl0 + 243);
    const auto *sgl0_245 = buffer.data(sgl0 + 245);
    const auto *sgl0_246 = buffer.data(sgl0 + 246);
    const auto *sgl0_248 = buffer.data(sgl0 + 248);
    const auto *sgl0_249 = buffer.data(sgl0 + 249);
    const auto *sgl0_250 = buffer.data(sgl0 + 250);
    const auto *sgl0_252 = buffer.data(sgl0 + 252);
    const auto *sgl0_269 = buffer.data(sgl0 + 269);
    const auto *sgl0_450 = buffer.data(sgl0 + 450);
    const auto *sgl0_453 = buffer.data(sgl0 + 453);

    const auto *sgk_136 = buffer.data(sgk + 136);
    const auto *sgk_143 = buffer.data(sgk + 143);
    const auto *sgk_144 = buffer.data(sgk + 144);
    const auto *sgk_147 = buffer.data(sgk + 147);
    const auto *sgk_150 = buffer.data(sgk + 150);
    const auto *sgk_154 = buffer.data(sgk + 154);
    const auto *sgk_159 = buffer.data(sgk + 159);
    const auto *sgk_164 = buffer.data(sgk + 164);
    const auto *sgk_172 = buffer.data(sgk + 172);
    const auto *sgk_174 = buffer.data(sgk + 174);
    const auto *sgk_175 = buffer.data(sgk + 175);
    const auto *sgk_176 = buffer.data(sgk + 176);
    const auto *sgk_177 = buffer.data(sgk + 177);
    const auto *sgk_178 = buffer.data(sgk + 178);
    const auto *sgk_179 = buffer.data(sgk + 179);
    const auto *sgk_180 = buffer.data(sgk + 180);
    const auto *sgk_181 = buffer.data(sgk + 181);
    const auto *sgk_182 = buffer.data(sgk + 182);
    const auto *sgk_183 = buffer.data(sgk + 183);
    const auto *sgk_185 = buffer.data(sgk + 185);
    const auto *sgk_186 = buffer.data(sgk + 186);
    const auto *sgk_188 = buffer.data(sgk + 188);
    const auto *sgk_189 = buffer.data(sgk + 189);
    const auto *sgk_190 = buffer.data(sgk + 190);
    const auto *sgk_192 = buffer.data(sgk + 192);
    const auto *sgk_193 = buffer.data(sgk + 193);
    const auto *sgk_194 = buffer.data(sgk + 194);
    const auto *sgk_195 = buffer.data(sgk + 195);
    const auto *sgk_197 = buffer.data(sgk + 197);
    const auto *sgk_198 = buffer.data(sgk + 198);
    const auto *sgk_199 = buffer.data(sgk + 199);
    const auto *sgk_200 = buffer.data(sgk + 200);
    const auto *sgk_208 = buffer.data(sgk + 208);
    const auto *sgk_210 = buffer.data(sgk + 210);
    const auto *sgk_211 = buffer.data(sgk + 211);
    const auto *sgk_212 = buffer.data(sgk + 212);
    const auto *sgk_213 = buffer.data(sgk + 213);
    const auto *sgk_214 = buffer.data(sgk + 214);
    const auto *sgk_215 = buffer.data(sgk + 215);
    const auto *sgk_216 = buffer.data(sgk + 216);
    const auto *sgk_218 = buffer.data(sgk + 218);
    const auto *sgk_279 = buffer.data(sgk + 279);
    const auto *sgk_280 = buffer.data(sgk + 280);
    const auto *sgk_281 = buffer.data(sgk + 281);
    const auto *sgk_282 = buffer.data(sgk + 282);
    const auto *sgk_283 = buffer.data(sgk + 283);
    const auto *sgk_284 = buffer.data(sgk + 284);
    const auto *sgk_285 = buffer.data(sgk + 285);
    const auto *sgk_286 = buffer.data(sgk + 286);
    const auto *sgk_287 = buffer.data(sgk + 287);
    const auto *sgk_316 = buffer.data(sgk + 316);
    const auto *sgk_317 = buffer.data(sgk + 317);
    const auto *sgk_318 = buffer.data(sgk + 318);
    const auto *sgk_319 = buffer.data(sgk + 319);
    const auto *sgk_320 = buffer.data(sgk + 320);
    const auto *sgk_321 = buffer.data(sgk + 321);
    const auto *sgk_322 = buffer.data(sgk + 322);
    const auto *sgk_323 = buffer.data(sgk + 323);
    const auto *sgk_324 = buffer.data(sgk + 324);
    const auto *sgk_327 = buffer.data(sgk + 327);
    const auto *sgk_329 = buffer.data(sgk + 329);
    const auto *sgk_330 = buffer.data(sgk + 330);
    const auto *sgk_333 = buffer.data(sgk + 333);
    const auto *sgk_334 = buffer.data(sgk + 334);
    const auto *sgk_336 = buffer.data(sgk + 336);
    const auto *sgk_338 = buffer.data(sgk + 338);
    const auto *sgk_339 = buffer.data(sgk + 339);
    const auto *sgk_341 = buffer.data(sgk + 341);
    const auto *sgk_342 = buffer.data(sgk + 342);
    const auto *sgk_344 = buffer.data(sgk + 344);
    const auto *sgk_345 = buffer.data(sgk + 345);
    const auto *sgk_347 = buffer.data(sgk + 347);
    const auto *sgk_348 = buffer.data(sgk + 348);
    const auto *sgk_349 = buffer.data(sgk + 349);
    const auto *sgk_351 = buffer.data(sgk + 351);
    const auto *sgk_352 = buffer.data(sgk + 352);
    const auto *sgk_353 = buffer.data(sgk + 353);
    const auto *sgk_354 = buffer.data(sgk + 354);
    const auto *sgk_355 = buffer.data(sgk + 355);
    const auto *sgk_356 = buffer.data(sgk + 356);
    const auto *sgk_357 = buffer.data(sgk + 357);
    const auto *sgk_358 = buffer.data(sgk + 358);
    const auto *sgk_359 = buffer.data(sgk + 359);
    const auto *sgk_360 = buffer.data(sgk + 360);
    const auto *sgk_363 = buffer.data(sgk + 363);

    const auto *sgl1_171 = buffer.data(sgl1 + 171);
    const auto *sgl1_225 = buffer.data(sgl1 + 225);
    const auto *sgl1_228 = buffer.data(sgl1 + 228);
    const auto *sgl1_230 = buffer.data(sgl1 + 230);
    const auto *sgl1_231 = buffer.data(sgl1 + 231);
    const auto *sgl1_234 = buffer.data(sgl1 + 234);
    const auto *sgl1_235 = buffer.data(sgl1 + 235);
    const auto *sgl1_237 = buffer.data(sgl1 + 237);
    const auto *sgl1_239 = buffer.data(sgl1 + 239);
    const auto *sgl1_240 = buffer.data(sgl1 + 240);
    const auto *sgl1_242 = buffer.data(sgl1 + 242);
    const auto *sgl1_243 = buffer.data(sgl1 + 243);
    const auto *sgl1_245 = buffer.data(sgl1 + 245);
    const auto *sgl1_246 = buffer.data(sgl1 + 246);
    const auto *sgl1_248 = buffer.data(sgl1 + 248);
    const auto *sgl1_249 = buffer.data(sgl1 + 249);
    const auto *sgl1_250 = buffer.data(sgl1 + 250);
    const auto *sgl1_252 = buffer.data(sgl1 + 252);
    const auto *sgl1_269 = buffer.data(sgl1 + 269);
    const auto *sgl1_450 = buffer.data(sgl1 + 450);
    const auto *sgl1_453 = buffer.data(sgl1 + 453);

    const auto *shi0_219 = buffer.data(shi0 + 219);
    const auto *shi0_220 = buffer.data(shi0 + 220);
    const auto *shi0_221 = buffer.data(shi0 + 221);
    const auto *shi0_222 = buffer.data(shi0 + 222);
    const auto *shi0_223 = buffer.data(shi0 + 223);
    const auto *shi0_245 = buffer.data(shi0 + 245);
    const auto *shi0_247 = buffer.data(shi0 + 247);
    const auto *shi0_248 = buffer.data(shi0 + 248);
    const auto *shi0_249 = buffer.data(shi0 + 249);
    const auto *shi0_250 = buffer.data(shi0 + 250);
    const auto *shi0_251 = buffer.data(shi0 + 251);
    const auto *shi0_252 = buffer.data(shi0 + 252);
    const auto *shi0_255 = buffer.data(shi0 + 255);
    const auto *shi0_257 = buffer.data(shi0 + 257);
    const auto *shi0_258 = buffer.data(shi0 + 258);
    const auto *shi0_261 = buffer.data(shi0 + 261);
    const auto *shi0_262 = buffer.data(shi0 + 262);
    const auto *shi0_264 = buffer.data(shi0 + 264);
    const auto *shi0_266 = buffer.data(shi0 + 266);
    const auto *shi0_267 = buffer.data(shi0 + 267);
    const auto *shi0_269 = buffer.data(shi0 + 269);
    const auto *shi0_270 = buffer.data(shi0 + 270);
    const auto *shi0_272 = buffer.data(shi0 + 272);
    const auto *shi0_273 = buffer.data(shi0 + 273);
    const auto *shi0_275 = buffer.data(shi0 + 275);
    const auto *shi0_276 = buffer.data(shi0 + 276);
    const auto *shi0_277 = buffer.data(shi0 + 277);
    const auto *shi0_278 = buffer.data(shi0 + 278);
    const auto *shi0_279 = buffer.data(shi0 + 279);

    const auto *shi1_219 = buffer.data(shi1 + 219);
    const auto *shi1_220 = buffer.data(shi1 + 220);
    const auto *shi1_221 = buffer.data(shi1 + 221);
    const auto *shi1_222 = buffer.data(shi1 + 222);
    const auto *shi1_223 = buffer.data(shi1 + 223);
    const auto *shi1_245 = buffer.data(shi1 + 245);
    const auto *shi1_247 = buffer.data(shi1 + 247);
    const auto *shi1_248 = buffer.data(shi1 + 248);
    const auto *shi1_249 = buffer.data(shi1 + 249);
    const auto *shi1_250 = buffer.data(shi1 + 250);
    const auto *shi1_251 = buffer.data(shi1 + 251);
    const auto *shi1_252 = buffer.data(shi1 + 252);
    const auto *shi1_255 = buffer.data(shi1 + 255);
    const auto *shi1_257 = buffer.data(shi1 + 257);
    const auto *shi1_258 = buffer.data(shi1 + 258);
    const auto *shi1_261 = buffer.data(shi1 + 261);
    const auto *shi1_262 = buffer.data(shi1 + 262);
    const auto *shi1_264 = buffer.data(shi1 + 264);
    const auto *shi1_266 = buffer.data(shi1 + 266);
    const auto *shi1_267 = buffer.data(shi1 + 267);
    const auto *shi1_269 = buffer.data(shi1 + 269);
    const auto *shi1_270 = buffer.data(shi1 + 270);
    const auto *shi1_272 = buffer.data(shi1 + 272);
    const auto *shi1_273 = buffer.data(shi1 + 273);
    const auto *shi1_275 = buffer.data(shi1 + 275);
    const auto *shi1_276 = buffer.data(shi1 + 276);
    const auto *shi1_277 = buffer.data(shi1 + 277);
    const auto *shi1_278 = buffer.data(shi1 + 278);
    const auto *shi1_279 = buffer.data(shi1 + 279);

    const auto *shk_272 = buffer.data(shk + 272);
    const auto *shk_279 = buffer.data(shk + 279);
    const auto *shk_280 = buffer.data(shk + 280);
    const auto *shk_281 = buffer.data(shk + 281);
    const auto *shk_282 = buffer.data(shk + 282);
    const auto *shk_283 = buffer.data(shk + 283);
    const auto *shk_284 = buffer.data(shk + 284);
    const auto *shk_285 = buffer.data(shk + 285);
    const auto *shk_286 = buffer.data(shk + 286);
    const auto *shk_287 = buffer.data(shk + 287);
    const auto *shk_288 = buffer.data(shk + 288);
    const auto *shk_290 = buffer.data(shk + 290);
    const auto *shk_291 = buffer.data(shk + 291);
    const auto *shk_293 = buffer.data(shk + 293);
    const auto *shk_294 = buffer.data(shk + 294);
    const auto *shk_297 = buffer.data(shk + 297);
    const auto *shk_298 = buffer.data(shk + 298);
    const auto *shk_302 = buffer.data(shk + 302);
    const auto *shk_303 = buffer.data(shk + 303);
    const auto *shk_308 = buffer.data(shk + 308);
    const auto *shk_316 = buffer.data(shk + 316);
    const auto *shk_317 = buffer.data(shk + 317);
    const auto *shk_318 = buffer.data(shk + 318);
    const auto *shk_319 = buffer.data(shk + 319);
    const auto *shk_320 = buffer.data(shk + 320);
    const auto *shk_321 = buffer.data(shk + 321);
    const auto *shk_322 = buffer.data(shk + 322);
    const auto *shk_323 = buffer.data(shk + 323);
    const auto *shk_324 = buffer.data(shk + 324);
    const auto *shk_326 = buffer.data(shk + 326);
    const auto *shk_327 = buffer.data(shk + 327);
    const auto *shk_329 = buffer.data(shk + 329);
    const auto *shk_330 = buffer.data(shk + 330);
    const auto *shk_333 = buffer.data(shk + 333);
    const auto *shk_334 = buffer.data(shk + 334);
    const auto *shk_336 = buffer.data(shk + 336);
    const auto *shk_338 = buffer.data(shk + 338);
    const auto *shk_339 = buffer.data(shk + 339);
    const auto *shk_341 = buffer.data(shk + 341);
    const auto *shk_342 = buffer.data(shk + 342);
    const auto *shk_344 = buffer.data(shk + 344);
    const auto *shk_345 = buffer.data(shk + 345);
    const auto *shk_347 = buffer.data(shk + 347);
    const auto *shk_348 = buffer.data(shk + 348);
    const auto *shk_349 = buffer.data(shk + 349);
    const auto *shk_351 = buffer.data(shk + 351);
    const auto *shk_352 = buffer.data(shk + 352);
    const auto *shk_353 = buffer.data(shk + 353);
    const auto *shk_354 = buffer.data(shk + 354);
    const auto *shk_355 = buffer.data(shk + 355);
    const auto *shk_356 = buffer.data(shk + 356);
    const auto *shk_357 = buffer.data(shk + 357);
    const auto *shk_358 = buffer.data(shk + 358);
    const auto *shk_359 = buffer.data(shk + 359);
    const auto *shk_360 = buffer.data(shk + 360);
    const auto *shk_362 = buffer.data(shk + 362);

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_y, sgk_164, sgk_279, sgk_280, \
                         sgk_281, shi0_223, shi1_223, shk_272, shk_279, shk_280, \
                         shk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_16 * sgk_164[k]
                   + f_3 * pc_y[k] * shk_272[k];

        t_342[k] = f_16 * sgk_279[k]
                   + f_12 * shi0_223[k]
                   - f_13 * shi1_223[k]
                   + f_3 * pc_x[k] * shk_279[k];

        t_343[k] = f_16 * sgk_280[k]
                   + f_3 * pc_x[k] * shk_280[k];

        t_344[k] = f_16 * sgk_281[k]
                   + f_3 * pc_x[k] * shk_281[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pc_x, sgk_282, sgk_283, sgk_284, \
                         sgk_285, sgk_286, shk_282, shk_283, shk_284, shk_285, \
                         shk_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_16 * sgk_282[k]
                   + f_3 * pc_x[k] * shk_282[k];

        t_346[k] = f_16 * sgk_283[k]
                   + f_3 * pc_x[k] * shk_283[k];

        t_347[k] = f_16 * sgk_284[k]
                   + f_3 * pc_x[k] * shk_284[k];

        t_348[k] = f_16 * sgk_285[k]
                   + f_3 * pc_x[k] * shk_285[k];

        t_349[k] = f_16 * sgk_286[k]
                   + f_3 * pc_x[k] * shk_286[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, pb_z, pc_x, pc_z, sgl0_171, sgk_136, sgk_287, \
                         sgl1_171, shk_280, shk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_16 * sgk_287[k]
                   + f_3 * pc_x[k] * shk_287[k];

        t_351[k] = pb_z[k] * sgl0_171[k]
                   - f_14 * pc_z[k] * sgl1_171[k];

        t_352[k] = f_15 * sgk_136[k]
                   + f_3 * pc_z[k] * shk_280[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, sgk_174, sgk_175, sgk_176, shi0_219, \
                         shi0_220, shi0_221, shi1_219, shi1_220, shi1_221, shk_282, shk_283, \
                         shk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_16 * sgk_174[k]
                   + f_4 * shi0_219[k]
                   - f_5 * shi1_219[k]
                   + f_3 * pc_y[k] * shk_282[k];

        t_354[k] = f_16 * sgk_175[k]
                   + f_6 * shi0_220[k]
                   - f_7 * shi1_220[k]
                   + f_3 * pc_y[k] * shk_283[k];

        t_355[k] = f_16 * sgk_176[k]
                   + f_8 * shi0_221[k]
                   - f_9 * shi1_221[k]
                   + f_3 * pc_y[k] * shk_284[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_y, sgk_177, sgk_178, sgk_179, shi0_222, \
                         shi0_223, shi1_222, shi1_223, shk_285, shk_286, \
                         shk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_16 * sgk_177[k]
                   + f_10 * shi0_222[k]
                   - f_11 * shi1_222[k]
                   + f_3 * pc_y[k] * shk_285[k];

        t_357[k] = f_16 * sgk_178[k]
                   + f_12 * shi0_223[k]
                   - f_13 * shi1_223[k]
                   + f_3 * pc_y[k] * shk_286[k];

        t_358[k] = f_16 * sgk_179[k]
                   + f_3 * pc_y[k] * shk_287[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pb_y, pc_y, pc_z, sgl0_225, sgk_143, \
                         sgk_144, sgk_180, sgl1_225, shi0_223, shi1_223, shk_287, \
                         shk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_15 * sgk_143[k]
                   + f_1 * shi0_223[k]
                   - f_2 * shi1_223[k]
                   + f_3 * pc_z[k] * shk_287[k];

        t_360[k] = pb_y[k] * sgl0_225[k]
                   - f_14 * pc_y[k] * sgl1_225[k];

        t_361[k] = f_15 * sgk_180[k]
                   + f_3 * pc_y[k] * shk_288[k];

        t_362[k] = f_16 * sgk_144[k]
                   + f_3 * pc_z[k] * shk_288[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pb_y, pc_y, sgl0_228, sgl0_230, sgl0_231, \
                         sgk_181, sgk_182, sgk_183, sgl1_228, sgl1_230, sgl1_231, \
                         shk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pb_y[k] * sgl0_228[k]
                   + f_16 * sgk_181[k]
                   - f_14 * pc_y[k] * sgl1_228[k];

        t_364[k] = f_15 * sgk_182[k]
                   + f_3 * pc_y[k] * shk_290[k];

        t_365[k] = pb_y[k] * sgl0_230[k]
                   - f_14 * pc_y[k] * sgl1_230[k];

        t_366[k] = pb_y[k] * sgl0_231[k]
                   + f_17 * sgk_183[k]
                   - f_14 * pc_y[k] * sgl1_231[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pb_y, pc_y, pc_z, sgl0_234, sgl0_235, \
                         sgk_147, sgk_185, sgk_186, sgl1_234, sgl1_235, shk_291, \
                         shk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_16 * sgk_147[k]
                   + f_3 * pc_z[k] * shk_291[k];

        t_368[k] = f_15 * sgk_185[k]
                   + f_3 * pc_y[k] * shk_293[k];

        t_369[k] = pb_y[k] * sgl0_234[k]
                   - f_14 * pc_y[k] * sgl1_234[k];

        t_370[k] = pb_y[k] * sgl0_235[k]
                   + f_18 * sgk_186[k]
                   - f_14 * pc_y[k] * sgl1_235[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pb_y, pc_y, pc_z, sgl0_237, sgl0_239, \
                         sgk_150, sgk_188, sgk_189, sgl1_237, sgl1_239, shk_294, \
                         shk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_16 * sgk_150[k]
                   + f_3 * pc_z[k] * shk_294[k];

        t_372[k] = pb_y[k] * sgl0_237[k]
                   + f_16 * sgk_188[k]
                   - f_14 * pc_y[k] * sgl1_237[k];

        t_373[k] = f_15 * sgk_189[k]
                   + f_3 * pc_y[k] * shk_297[k];

        t_374[k] = pb_y[k] * sgl0_239[k]
                   - f_14 * pc_y[k] * sgl1_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pb_y, pc_y, pc_z, sgl0_240, sgl0_242, sgk_154, \
                         sgk_190, sgk_192, sgl1_240, sgl1_242, \
                         shk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = pb_y[k] * sgl0_240[k]
                   + f_0 * sgk_190[k]
                   - f_14 * pc_y[k] * sgl1_240[k];

        t_376[k] = f_16 * sgk_154[k]
                   + f_3 * pc_z[k] * shk_298[k];

        t_377[k] = pb_y[k] * sgl0_242[k]
                   + f_17 * sgk_192[k]
                   - f_14 * pc_y[k] * sgl1_242[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pb_y, pc_y, sgl0_243, sgl0_245, sgl0_246, \
                         sgk_193, sgk_194, sgk_195, sgl1_243, sgl1_245, sgl1_246, \
                         shk_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pb_y[k] * sgl0_243[k]
                   + f_16 * sgk_193[k]
                   - f_14 * pc_y[k] * sgl1_243[k];

        t_379[k] = f_15 * sgk_194[k]
                   + f_3 * pc_y[k] * shk_302[k];

        t_380[k] = pb_y[k] * sgl0_245[k]
                   - f_14 * pc_y[k] * sgl1_245[k];

        t_381[k] = pb_y[k] * sgl0_246[k]
                   + f_19 * sgk_195[k]
                   - f_14 * pc_y[k] * sgl1_246[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pb_y, pc_y, pc_z, sgl0_248, sgl0_249, sgk_159, \
                         sgk_197, sgk_198, sgl1_248, sgl1_249, \
                         shk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_16 * sgk_159[k]
                   + f_3 * pc_z[k] * shk_303[k];

        t_383[k] = pb_y[k] * sgl0_248[k]
                   + f_18 * sgk_197[k]
                   - f_14 * pc_y[k] * sgl1_248[k];

        t_384[k] = pb_y[k] * sgl0_249[k]
                   + f_17 * sgk_198[k]
                   - f_14 * pc_y[k] * sgl1_249[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pb_y, pc_x, pc_y, sgl0_250, sgl0_252, \
                         sgk_199, sgk_200, sgk_316, sgl1_250, sgl1_252, shk_308, \
                         shk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pb_y[k] * sgl0_250[k]
                   + f_16 * sgk_199[k]
                   - f_14 * pc_y[k] * sgl1_250[k];

        t_386[k] = f_15 * sgk_200[k]
                   + f_3 * pc_y[k] * shk_308[k];

        t_387[k] = pb_y[k] * sgl0_252[k]
                   - f_14 * pc_y[k] * sgl1_252[k];

        t_388[k] = f_16 * sgk_316[k]
                   + f_3 * pc_x[k] * shk_316[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pc_x, sgk_317, sgk_318, sgk_319, \
                         sgk_320, sgk_321, shk_317, shk_318, shk_319, shk_320, \
                         shk_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_16 * sgk_317[k]
                   + f_3 * pc_x[k] * shk_317[k];

        t_390[k] = f_16 * sgk_318[k]
                   + f_3 * pc_x[k] * shk_318[k];

        t_391[k] = f_16 * sgk_319[k]
                   + f_3 * pc_x[k] * shk_319[k];

        t_392[k] = f_16 * sgk_320[k]
                   + f_3 * pc_x[k] * shk_320[k];

        t_393[k] = f_16 * sgk_321[k]
                   + f_3 * pc_x[k] * shk_321[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pc_x, pc_y, pc_z, sgk_172, sgk_208, \
                         sgk_322, sgk_323, shi0_245, shi1_245, shk_316, shk_322, \
                         shk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_16 * sgk_322[k]
                   + f_3 * pc_x[k] * shk_322[k];

        t_395[k] = f_16 * sgk_323[k]
                   + f_3 * pc_x[k] * shk_323[k];

        t_396[k] = f_15 * sgk_208[k]
                   + f_1 * shi0_245[k]
                   - f_2 * shi1_245[k]
                   + f_3 * pc_y[k] * shk_316[k];

        t_397[k] = f_16 * sgk_172[k]
                   + f_3 * pc_z[k] * shk_316[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, sgk_210, sgk_211, sgk_212, shi0_247, \
                         shi0_248, shi0_249, shi1_247, shi1_248, shi1_249, shk_318, shk_319, \
                         shk_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_15 * sgk_210[k]
                   + f_4 * shi0_247[k]
                   - f_5 * shi1_247[k]
                   + f_3 * pc_y[k] * shk_318[k];

        t_399[k] = f_15 * sgk_211[k]
                   + f_6 * shi0_248[k]
                   - f_7 * shi1_248[k]
                   + f_3 * pc_y[k] * shk_319[k];

        t_400[k] = f_15 * sgk_212[k]
                   + f_8 * shi0_249[k]
                   - f_9 * shi1_249[k]
                   + f_3 * pc_y[k] * shk_320[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_y, sgk_213, sgk_214, sgk_215, shi0_250, \
                         shi0_251, shi1_250, shi1_251, shk_321, shk_322, \
                         shk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * sgk_213[k]
                   + f_10 * shi0_250[k]
                   - f_11 * shi1_250[k]
                   + f_3 * pc_y[k] * shk_321[k];

        t_402[k] = f_15 * sgk_214[k]
                   + f_12 * shi0_251[k]
                   - f_13 * shi1_251[k]
                   + f_3 * pc_y[k] * shk_322[k];

        t_403[k] = f_15 * sgk_215[k]
                   + f_3 * pc_y[k] * shk_323[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pb_y, pc_x, pc_y, pc_z, sgl0_269, \
                         sgk_180, sgk_324, sgl1_269, shi0_252, shi1_252, \
                         shk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pb_y[k] * sgl0_269[k]
                   - f_14 * pc_y[k] * sgl1_269[k];

        t_405[k] = f_16 * sgk_324[k]
                   + f_1 * shi0_252[k]
                   - f_2 * shi1_252[k]
                   + f_3 * pc_x[k] * shk_324[k];

        t_406[k] = f_3 * pc_y[k] * shk_324[k];

        t_407[k] = f_17 * sgk_180[k]
                   + f_3 * pc_z[k] * shk_324[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pc_x, pc_y, sgk_327, sgk_329, shi0_255, \
                         shi0_257, shi1_255, shi1_257, shk_326, shk_327, \
                         shk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_16 * sgk_327[k]
                   + f_4 * shi0_255[k]
                   - f_5 * shi1_255[k]
                   + f_3 * pc_x[k] * shk_327[k];

        t_409[k] = f_3 * pc_y[k] * shk_326[k];

        t_410[k] = f_16 * sgk_329[k]
                   + f_4 * shi0_257[k]
                   - f_5 * shi1_257[k]
                   + f_3 * pc_x[k] * shk_329[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pc_x, pc_y, pc_z, sgk_183, sgk_330, shi0_258, \
                         shi1_258, shk_327, shk_329, shk_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_16 * sgk_330[k]
                   + f_6 * shi0_258[k]
                   - f_7 * shi1_258[k]
                   + f_3 * pc_x[k] * shk_330[k];

        t_412[k] = f_17 * sgk_183[k]
                   + f_3 * pc_z[k] * shk_327[k];

        t_413[k] = f_3 * pc_y[k] * shk_329[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_z, sgk_186, sgk_333, sgk_334, shi0_261, \
                         shi0_262, shi1_261, shi1_262, shk_330, shk_333, \
                         shk_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_16 * sgk_333[k]
                   + f_6 * shi0_261[k]
                   - f_7 * shi1_261[k]
                   + f_3 * pc_x[k] * shk_333[k];

        t_415[k] = f_16 * sgk_334[k]
                   + f_8 * shi0_262[k]
                   - f_9 * shi1_262[k]
                   + f_3 * pc_x[k] * shk_334[k];

        t_416[k] = f_17 * sgk_186[k]
                   + f_3 * pc_z[k] * shk_330[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_x, pc_y, sgk_336, sgk_338, shi0_264, \
                         shi0_266, shi1_264, shi1_266, shk_333, shk_336, \
                         shk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_16 * sgk_336[k]
                   + f_8 * shi0_264[k]
                   - f_9 * shi1_264[k]
                   + f_3 * pc_x[k] * shk_336[k];

        t_418[k] = f_3 * pc_y[k] * shk_333[k];

        t_419[k] = f_16 * sgk_338[k]
                   + f_8 * shi0_266[k]
                   - f_9 * shi1_266[k]
                   + f_3 * pc_x[k] * shk_338[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pc_x, pc_z, sgk_190, sgk_339, sgk_341, shi0_267, \
                         shi0_269, shi1_267, shi1_269, shk_334, shk_339, \
                         shk_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_16 * sgk_339[k]
                   + f_10 * shi0_267[k]
                   - f_11 * shi1_267[k]
                   + f_3 * pc_x[k] * shk_339[k];

        t_421[k] = f_17 * sgk_190[k]
                   + f_3 * pc_z[k] * shk_334[k];

        t_422[k] = f_16 * sgk_341[k]
                   + f_10 * shi0_269[k]
                   - f_11 * shi1_269[k]
                   + f_3 * pc_x[k] * shk_341[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pc_x, pc_y, sgk_342, sgk_344, shi0_270, \
                         shi0_272, shi1_270, shi1_272, shk_338, shk_342, \
                         shk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_16 * sgk_342[k]
                   + f_10 * shi0_270[k]
                   - f_11 * shi1_270[k]
                   + f_3 * pc_x[k] * shk_342[k];

        t_424[k] = f_3 * pc_y[k] * shk_338[k];

        t_425[k] = f_16 * sgk_344[k]
                   + f_10 * shi0_272[k]
                   - f_11 * shi1_272[k]
                   + f_3 * pc_x[k] * shk_344[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pc_x, pc_z, sgk_195, sgk_345, sgk_347, shi0_273, \
                         shi0_275, shi1_273, shi1_275, shk_339, shk_345, \
                         shk_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_16 * sgk_345[k]
                   + f_12 * shi0_273[k]
                   - f_13 * shi1_273[k]
                   + f_3 * pc_x[k] * shk_345[k];

        t_427[k] = f_17 * sgk_195[k]
                   + f_3 * pc_z[k] * shk_339[k];

        t_428[k] = f_16 * sgk_347[k]
                   + f_12 * shi0_275[k]
                   - f_13 * shi1_275[k]
                   + f_3 * pc_x[k] * shk_347[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_x, pc_y, sgk_348, sgk_349, shi0_276, \
                         shi0_277, shi1_276, shi1_277, shk_344, shk_348, \
                         shk_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_16 * sgk_348[k]
                   + f_12 * shi0_276[k]
                   - f_13 * shi1_276[k]
                   + f_3 * pc_x[k] * shk_348[k];

        t_430[k] = f_16 * sgk_349[k]
                   + f_12 * shi0_277[k]
                   - f_13 * shi1_277[k]
                   + f_3 * pc_x[k] * shk_349[k];

        t_431[k] = f_3 * pc_y[k] * shk_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, sgk_351, sgk_352, sgk_353, sgk_354, \
                         shi0_279, shi1_279, shk_351, shk_352, shk_353, \
                         shk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_16 * sgk_351[k]
                   + f_12 * shi0_279[k]
                   - f_13 * shi1_279[k]
                   + f_3 * pc_x[k] * shk_351[k];

        t_433[k] = f_16 * sgk_352[k]
                   + f_3 * pc_x[k] * shk_352[k];

        t_434[k] = f_16 * sgk_353[k]
                   + f_3 * pc_x[k] * shk_353[k];

        t_435[k] = f_16 * sgk_354[k]
                   + f_3 * pc_x[k] * shk_354[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, sgk_355, sgk_356, sgk_357, \
                         sgk_358, sgk_359, shk_355, shk_356, shk_357, shk_358, \
                         shk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_16 * sgk_355[k]
                   + f_3 * pc_x[k] * shk_355[k];

        t_437[k] = f_16 * sgk_356[k]
                   + f_3 * pc_x[k] * shk_356[k];

        t_438[k] = f_16 * sgk_357[k]
                   + f_3 * pc_x[k] * shk_357[k];

        t_439[k] = f_16 * sgk_358[k]
                   + f_3 * pc_x[k] * shk_358[k];

        t_440[k] = f_16 * sgk_359[k]
                   + f_3 * pc_x[k] * shk_359[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_y, pc_z, sgk_208, shi0_273, shi0_275, \
                         shi0_276, shi1_273, shi1_275, shi1_276, shk_352, shk_354, \
                         shk_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * shi0_273[k]
                   - f_2 * shi1_273[k]
                   + f_3 * pc_y[k] * shk_352[k];

        t_442[k] = f_17 * sgk_208[k]
                   + f_3 * pc_z[k] * shk_352[k];

        t_443[k] = f_4 * shi0_275[k]
                   - f_5 * shi1_275[k]
                   + f_3 * pc_y[k] * shk_354[k];

        t_444[k] = f_6 * shi0_276[k]
                   - f_7 * shi1_276[k]
                   + f_3 * pc_y[k] * shk_355[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_y, shi0_277, shi0_278, shi0_279, \
                         shi1_277, shi1_278, shi1_279, shk_356, shk_357, shk_358, \
                         shk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_8 * shi0_277[k]
                   - f_9 * shi1_277[k]
                   + f_3 * pc_y[k] * shk_356[k];

        t_446[k] = f_10 * shi0_278[k]
                   - f_11 * shi1_278[k]
                   + f_3 * pc_y[k] * shk_357[k];

        t_447[k] = f_12 * shi0_279[k]
                   - f_13 * shi1_279[k]
                   + f_3 * pc_y[k] * shk_358[k];

        t_448[k] = f_3 * pc_y[k] * shk_359[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pb_x, pc_x, pc_y, pc_z, sgl0_450, sgk_215, \
                         sgk_216, sgk_360, sgl1_450, shi0_279, shi1_279, shk_359, \
                         shk_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_17 * sgk_215[k]
                   + f_1 * shi0_279[k]
                   - f_2 * shi1_279[k]
                   + f_3 * pc_z[k] * shk_359[k];

        t_450[k] = pb_x[k] * sgl0_450[k]
                   + f_20 * sgk_360[k]
                   - f_14 * pc_x[k] * sgl1_450[k];

        t_451[k] = f_18 * sgk_216[k]
                   + f_3 * pc_y[k] * shk_360[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, pb_x, pc_x, pc_y, pc_z, sgl0_453, sgk_218, \
                         sgk_363, sgl1_453, shk_360, shk_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_3 * pc_z[k] * shk_360[k];

        t_453[k] = pb_x[k] * sgl0_453[k]
                   + f_19 * sgk_363[k]
                   - f_14 * pc_x[k] * sgl1_453[k];

        t_454[k] = f_18 * sgk_218[k]
                   + f_3 * pc_y[k] * shk_362[k];
    }
}

static auto
compute_prim_shl_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgl0,
                                                          const size_t sgk, const size_t sgl1,
                                                          const size_t shk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 4.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgl0_270 = buffer.data(sgl0 + 270);
    const auto *sgl0_273 = buffer.data(sgl0 + 273);
    const auto *sgl0_276 = buffer.data(sgl0 + 276);
    const auto *sgl0_280 = buffer.data(sgl0 + 280);
    const auto *sgl0_285 = buffer.data(sgl0 + 285);
    const auto *sgl0_291 = buffer.data(sgl0 + 291);
    const auto *sgl0_455 = buffer.data(sgl0 + 455);
    const auto *sgl0_456 = buffer.data(sgl0 + 456);
    const auto *sgl0_459 = buffer.data(sgl0 + 459);
    const auto *sgl0_460 = buffer.data(sgl0 + 460);
    const auto *sgl0_462 = buffer.data(sgl0 + 462);
    const auto *sgl0_464 = buffer.data(sgl0 + 464);
    const auto *sgl0_465 = buffer.data(sgl0 + 465);
    const auto *sgl0_467 = buffer.data(sgl0 + 467);
    const auto *sgl0_468 = buffer.data(sgl0 + 468);
    const auto *sgl0_470 = buffer.data(sgl0 + 470);
    const auto *sgl0_471 = buffer.data(sgl0 + 471);
    const auto *sgl0_473 = buffer.data(sgl0 + 473);
    const auto *sgl0_474 = buffer.data(sgl0 + 474);
    const auto *sgl0_475 = buffer.data(sgl0 + 475);
    const auto *sgl0_477 = buffer.data(sgl0 + 477);
    const auto *sgl0_486 = buffer.data(sgl0 + 486);
    const auto *sgl0_488 = buffer.data(sgl0 + 488);
    const auto *sgl0_489 = buffer.data(sgl0 + 489);
    const auto *sgl0_490 = buffer.data(sgl0 + 490);
    const auto *sgl0_491 = buffer.data(sgl0 + 491);
    const auto *sgl0_492 = buffer.data(sgl0 + 492);
    const auto *sgl0_494 = buffer.data(sgl0 + 494);
    const auto *sgl0_500 = buffer.data(sgl0 + 500);
    const auto *sgl0_504 = buffer.data(sgl0 + 504);
    const auto *sgl0_507 = buffer.data(sgl0 + 507);
    const auto *sgl0_509 = buffer.data(sgl0 + 509);
    const auto *sgl0_512 = buffer.data(sgl0 + 512);
    const auto *sgl0_513 = buffer.data(sgl0 + 513);
    const auto *sgl0_515 = buffer.data(sgl0 + 515);
    const auto *sgl0_518 = buffer.data(sgl0 + 518);
    const auto *sgl0_519 = buffer.data(sgl0 + 519);
    const auto *sgl0_520 = buffer.data(sgl0 + 520);
    const auto *sgl0_522 = buffer.data(sgl0 + 522);
    const auto *sgl0_531 = buffer.data(sgl0 + 531);
    const auto *sgl0_533 = buffer.data(sgl0 + 533);
    const auto *sgl0_534 = buffer.data(sgl0 + 534);
    const auto *sgl0_535 = buffer.data(sgl0 + 535);
    const auto *sgl0_536 = buffer.data(sgl0 + 536);
    const auto *sgl0_537 = buffer.data(sgl0 + 537);
    const auto *sgl0_539 = buffer.data(sgl0 + 539);
    const auto *sgl0_540 = buffer.data(sgl0 + 540);
    const auto *sgl0_543 = buffer.data(sgl0 + 543);
    const auto *sgl0_545 = buffer.data(sgl0 + 545);
    const auto *sgl0_546 = buffer.data(sgl0 + 546);
    const auto *sgl0_549 = buffer.data(sgl0 + 549);
    const auto *sgl0_550 = buffer.data(sgl0 + 550);
    const auto *sgl0_552 = buffer.data(sgl0 + 552);
    const auto *sgl0_554 = buffer.data(sgl0 + 554);
    const auto *sgl0_555 = buffer.data(sgl0 + 555);
    const auto *sgl0_557 = buffer.data(sgl0 + 557);
    const auto *sgl0_558 = buffer.data(sgl0 + 558);
    const auto *sgl0_560 = buffer.data(sgl0 + 560);
    const auto *sgl0_561 = buffer.data(sgl0 + 561);
    const auto *sgl0_563 = buffer.data(sgl0 + 563);
    const auto *sgl0_564 = buffer.data(sgl0 + 564);
    const auto *sgl0_565 = buffer.data(sgl0 + 565);
    const auto *sgl0_567 = buffer.data(sgl0 + 567);

    const auto *sgk_216 = buffer.data(sgk + 216);
    const auto *sgk_219 = buffer.data(sgk + 219);
    const auto *sgk_221 = buffer.data(sgk + 221);
    const auto *sgk_222 = buffer.data(sgk + 222);
    const auto *sgk_225 = buffer.data(sgk + 225);
    const auto *sgk_226 = buffer.data(sgk + 226);
    const auto *sgk_230 = buffer.data(sgk + 230);
    const auto *sgk_231 = buffer.data(sgk + 231);
    const auto *sgk_236 = buffer.data(sgk + 236);
    const auto *sgk_244 = buffer.data(sgk + 244);
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
    const auto *sgk_287 = buffer.data(sgk + 287);
    const auto *sgk_288 = buffer.data(sgk + 288);
    const auto *sgk_290 = buffer.data(sgk + 290);
    const auto *sgk_293 = buffer.data(sgk + 293);
    const auto *sgk_297 = buffer.data(sgk + 297);
    const auto *sgk_302 = buffer.data(sgk + 302);
    const auto *sgk_308 = buffer.data(sgk + 308);
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
    const auto *sgk_401 = buffer.data(sgk + 401);
    const auto *sgk_405 = buffer.data(sgk + 405);
    const auto *sgk_408 = buffer.data(sgk + 408);
    const auto *sgk_410 = buffer.data(sgk + 410);
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

    const auto *sgl1_270 = buffer.data(sgl1 + 270);
    const auto *sgl1_273 = buffer.data(sgl1 + 273);
    const auto *sgl1_276 = buffer.data(sgl1 + 276);
    const auto *sgl1_280 = buffer.data(sgl1 + 280);
    const auto *sgl1_285 = buffer.data(sgl1 + 285);
    const auto *sgl1_291 = buffer.data(sgl1 + 291);
    const auto *sgl1_455 = buffer.data(sgl1 + 455);
    const auto *sgl1_456 = buffer.data(sgl1 + 456);
    const auto *sgl1_459 = buffer.data(sgl1 + 459);
    const auto *sgl1_460 = buffer.data(sgl1 + 460);
    const auto *sgl1_462 = buffer.data(sgl1 + 462);
    const auto *sgl1_464 = buffer.data(sgl1 + 464);
    const auto *sgl1_465 = buffer.data(sgl1 + 465);
    const auto *sgl1_467 = buffer.data(sgl1 + 467);
    const auto *sgl1_468 = buffer.data(sgl1 + 468);
    const auto *sgl1_470 = buffer.data(sgl1 + 470);
    const auto *sgl1_471 = buffer.data(sgl1 + 471);
    const auto *sgl1_473 = buffer.data(sgl1 + 473);
    const auto *sgl1_474 = buffer.data(sgl1 + 474);
    const auto *sgl1_475 = buffer.data(sgl1 + 475);
    const auto *sgl1_477 = buffer.data(sgl1 + 477);
    const auto *sgl1_486 = buffer.data(sgl1 + 486);
    const auto *sgl1_488 = buffer.data(sgl1 + 488);
    const auto *sgl1_489 = buffer.data(sgl1 + 489);
    const auto *sgl1_490 = buffer.data(sgl1 + 490);
    const auto *sgl1_491 = buffer.data(sgl1 + 491);
    const auto *sgl1_492 = buffer.data(sgl1 + 492);
    const auto *sgl1_494 = buffer.data(sgl1 + 494);
    const auto *sgl1_500 = buffer.data(sgl1 + 500);
    const auto *sgl1_504 = buffer.data(sgl1 + 504);
    const auto *sgl1_507 = buffer.data(sgl1 + 507);
    const auto *sgl1_509 = buffer.data(sgl1 + 509);
    const auto *sgl1_512 = buffer.data(sgl1 + 512);
    const auto *sgl1_513 = buffer.data(sgl1 + 513);
    const auto *sgl1_515 = buffer.data(sgl1 + 515);
    const auto *sgl1_518 = buffer.data(sgl1 + 518);
    const auto *sgl1_519 = buffer.data(sgl1 + 519);
    const auto *sgl1_520 = buffer.data(sgl1 + 520);
    const auto *sgl1_522 = buffer.data(sgl1 + 522);
    const auto *sgl1_531 = buffer.data(sgl1 + 531);
    const auto *sgl1_533 = buffer.data(sgl1 + 533);
    const auto *sgl1_534 = buffer.data(sgl1 + 534);
    const auto *sgl1_535 = buffer.data(sgl1 + 535);
    const auto *sgl1_536 = buffer.data(sgl1 + 536);
    const auto *sgl1_537 = buffer.data(sgl1 + 537);
    const auto *sgl1_539 = buffer.data(sgl1 + 539);
    const auto *sgl1_540 = buffer.data(sgl1 + 540);
    const auto *sgl1_543 = buffer.data(sgl1 + 543);
    const auto *sgl1_545 = buffer.data(sgl1 + 545);
    const auto *sgl1_546 = buffer.data(sgl1 + 546);
    const auto *sgl1_549 = buffer.data(sgl1 + 549);
    const auto *sgl1_550 = buffer.data(sgl1 + 550);
    const auto *sgl1_552 = buffer.data(sgl1 + 552);
    const auto *sgl1_554 = buffer.data(sgl1 + 554);
    const auto *sgl1_555 = buffer.data(sgl1 + 555);
    const auto *sgl1_557 = buffer.data(sgl1 + 557);
    const auto *sgl1_558 = buffer.data(sgl1 + 558);
    const auto *sgl1_560 = buffer.data(sgl1 + 560);
    const auto *sgl1_561 = buffer.data(sgl1 + 561);
    const auto *sgl1_563 = buffer.data(sgl1 + 563);
    const auto *sgl1_564 = buffer.data(sgl1 + 564);
    const auto *sgl1_565 = buffer.data(sgl1 + 565);
    const auto *sgl1_567 = buffer.data(sgl1 + 567);

    const auto *shk_363 = buffer.data(shk + 363);
    const auto *shk_365 = buffer.data(shk + 365);
    const auto *shk_366 = buffer.data(shk + 366);
    const auto *shk_369 = buffer.data(shk + 369);
    const auto *shk_370 = buffer.data(shk + 370);
    const auto *shk_374 = buffer.data(shk + 374);
    const auto *shk_375 = buffer.data(shk + 375);
    const auto *shk_380 = buffer.data(shk + 380);
    const auto *shk_388 = buffer.data(shk + 388);
    const auto *shk_389 = buffer.data(shk + 389);
    const auto *shk_390 = buffer.data(shk + 390);
    const auto *shk_391 = buffer.data(shk + 391);
    const auto *shk_392 = buffer.data(shk + 392);
    const auto *shk_393 = buffer.data(shk + 393);
    const auto *shk_394 = buffer.data(shk + 394);
    const auto *shk_395 = buffer.data(shk + 395);
    const auto *shk_396 = buffer.data(shk + 396);
    const auto *shk_398 = buffer.data(shk + 398);
    const auto *shk_399 = buffer.data(shk + 399);
    const auto *shk_401 = buffer.data(shk + 401);
    const auto *shk_402 = buffer.data(shk + 402);
    const auto *shk_405 = buffer.data(shk + 405);
    const auto *shk_406 = buffer.data(shk + 406);
    const auto *shk_410 = buffer.data(shk + 410);
    const auto *shk_411 = buffer.data(shk + 411);
    const auto *shk_416 = buffer.data(shk + 416);
    const auto *shk_424 = buffer.data(shk + 424);
    const auto *shk_425 = buffer.data(shk + 425);
    const auto *shk_426 = buffer.data(shk + 426);
    const auto *shk_427 = buffer.data(shk + 427);
    const auto *shk_428 = buffer.data(shk + 428);
    const auto *shk_429 = buffer.data(shk + 429);
    const auto *shk_430 = buffer.data(shk + 430);
    const auto *shk_431 = buffer.data(shk + 431);
    const auto *shk_432 = buffer.data(shk + 432);
    const auto *shk_434 = buffer.data(shk + 434);
    const auto *shk_435 = buffer.data(shk + 435);
    const auto *shk_437 = buffer.data(shk + 437);
    const auto *shk_438 = buffer.data(shk + 438);
    const auto *shk_441 = buffer.data(shk + 441);
    const auto *shk_442 = buffer.data(shk + 442);
    const auto *shk_446 = buffer.data(shk + 446);
    const auto *shk_447 = buffer.data(shk + 447);
    const auto *shk_452 = buffer.data(shk + 452);
    const auto *shk_460 = buffer.data(shk + 460);
    const auto *shk_461 = buffer.data(shk + 461);

#pragma omp simd aligned(t_455, t_456, t_457, pb_x, pc_x, pc_z, sgl0_455, sgl0_456, sgk_365, \
                         sgk_366, sgl1_455, sgl1_456, shk_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = pb_x[k] * sgl0_455[k]
                   + f_19 * sgk_365[k]
                   - f_14 * pc_x[k] * sgl1_455[k];

        t_456[k] = pb_x[k] * sgl0_456[k]
                   + f_0 * sgk_366[k]
                   - f_14 * pc_x[k] * sgl1_456[k];

        t_457[k] = f_3 * pc_z[k] * shk_363[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, pb_x, pc_x, pc_y, sgl0_459, sgl0_460, sgk_221, \
                         sgk_369, sgk_370, sgl1_459, sgl1_460, \
                         shk_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = f_18 * sgk_221[k]
                   + f_3 * pc_y[k] * shk_365[k];

        t_459[k] = pb_x[k] * sgl0_459[k]
                   + f_0 * sgk_369[k]
                   - f_14 * pc_x[k] * sgl1_459[k];

        t_460[k] = pb_x[k] * sgl0_460[k]
                   + f_18 * sgk_370[k]
                   - f_14 * pc_x[k] * sgl1_460[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, pb_x, pc_x, pc_y, pc_z, sgl0_462, sgk_225, \
                         sgk_372, sgl1_462, shk_366, shk_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_3 * pc_z[k] * shk_366[k];

        t_462[k] = pb_x[k] * sgl0_462[k]
                   + f_18 * sgk_372[k]
                   - f_14 * pc_x[k] * sgl1_462[k];

        t_463[k] = f_18 * sgk_225[k]
                   + f_3 * pc_y[k] * shk_369[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, pb_x, pc_x, pc_z, sgl0_464, sgl0_465, sgk_374, \
                         sgk_375, sgl1_464, sgl1_465, shk_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = pb_x[k] * sgl0_464[k]
                   + f_18 * sgk_374[k]
                   - f_14 * pc_x[k] * sgl1_464[k];

        t_465[k] = pb_x[k] * sgl0_465[k]
                   + f_17 * sgk_375[k]
                   - f_14 * pc_x[k] * sgl1_465[k];

        t_466[k] = f_3 * pc_z[k] * shk_370[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, pb_x, pc_x, pc_y, sgl0_467, sgl0_468, sgk_230, \
                         sgk_377, sgk_378, sgl1_467, sgl1_468, \
                         shk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = pb_x[k] * sgl0_467[k]
                   + f_17 * sgk_377[k]
                   - f_14 * pc_x[k] * sgl1_467[k];

        t_468[k] = pb_x[k] * sgl0_468[k]
                   + f_17 * sgk_378[k]
                   - f_14 * pc_x[k] * sgl1_468[k];

        t_469[k] = f_18 * sgk_230[k]
                   + f_3 * pc_y[k] * shk_374[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pb_x, pc_x, pc_z, sgl0_470, sgl0_471, sgk_380, \
                         sgk_381, sgl1_470, sgl1_471, shk_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = pb_x[k] * sgl0_470[k]
                   + f_17 * sgk_380[k]
                   - f_14 * pc_x[k] * sgl1_470[k];

        t_471[k] = pb_x[k] * sgl0_471[k]
                   + f_16 * sgk_381[k]
                   - f_14 * pc_x[k] * sgl1_471[k];

        t_472[k] = f_3 * pc_z[k] * shk_375[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pb_x, pc_x, sgl0_473, sgl0_474, sgl0_475, \
                         sgk_383, sgk_384, sgk_385, sgl1_473, sgl1_474, \
                         sgl1_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = pb_x[k] * sgl0_473[k]
                   + f_16 * sgk_383[k]
                   - f_14 * pc_x[k] * sgl1_473[k];

        t_474[k] = pb_x[k] * sgl0_474[k]
                   + f_16 * sgk_384[k]
                   - f_14 * pc_x[k] * sgl1_474[k];

        t_475[k] = pb_x[k] * sgl0_475[k]
                   + f_16 * sgk_385[k]
                   - f_14 * pc_x[k] * sgl1_475[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, pb_x, pc_x, pc_y, sgl0_477, sgk_236, \
                         sgk_387, sgk_388, sgk_389, sgl1_477, shk_380, shk_388, \
                         shk_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_18 * sgk_236[k]
                   + f_3 * pc_y[k] * shk_380[k];

        t_477[k] = pb_x[k] * sgl0_477[k]
                   + f_16 * sgk_387[k]
                   - f_14 * pc_x[k] * sgl1_477[k];

        t_478[k] = f_15 * sgk_388[k]
                   + f_3 * pc_x[k] * shk_388[k];

        t_479[k] = f_15 * sgk_389[k]
                   + f_3 * pc_x[k] * shk_389[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, pc_x, sgk_390, sgk_391, sgk_392, \
                         sgk_393, sgk_394, shk_390, shk_391, shk_392, shk_393, \
                         shk_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = f_15 * sgk_390[k]
                   + f_3 * pc_x[k] * shk_390[k];

        t_481[k] = f_15 * sgk_391[k]
                   + f_3 * pc_x[k] * shk_391[k];

        t_482[k] = f_15 * sgk_392[k]
                   + f_3 * pc_x[k] * shk_392[k];

        t_483[k] = f_15 * sgk_393[k]
                   + f_3 * pc_x[k] * shk_393[k];

        t_484[k] = f_15 * sgk_394[k]
                   + f_3 * pc_x[k] * shk_394[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pb_x, pc_x, pc_z, sgl0_486, sgl0_488, \
                         sgk_395, sgl1_486, sgl1_488, shk_388, \
                         shk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_15 * sgk_395[k]
                   + f_3 * pc_x[k] * shk_395[k];

        t_486[k] = pb_x[k] * sgl0_486[k]
                   - f_14 * pc_x[k] * sgl1_486[k];

        t_487[k] = f_3 * pc_z[k] * shk_388[k];

        t_488[k] = pb_x[k] * sgl0_488[k]
                   - f_14 * pc_x[k] * sgl1_488[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, pb_x, pc_x, sgl0_489, sgl0_490, sgl0_491, \
                         sgl0_492, sgl1_489, sgl1_490, sgl1_491, \
                         sgl1_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = pb_x[k] * sgl0_489[k]
                   - f_14 * pc_x[k] * sgl1_489[k];

        t_490[k] = pb_x[k] * sgl0_490[k]
                   - f_14 * pc_x[k] * sgl1_490[k];

        t_491[k] = pb_x[k] * sgl0_491[k]
                   - f_14 * pc_x[k] * sgl1_491[k];

        t_492[k] = pb_x[k] * sgl0_492[k]
                   - f_14 * pc_x[k] * sgl1_492[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, pb_x, pb_z, pc_x, pc_y, pc_z, sgl0_270, \
                         sgl0_494, sgk_251, sgl1_270, sgl1_494, \
                         shk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = f_18 * sgk_251[k]
                   + f_3 * pc_y[k] * shk_395[k];

        t_494[k] = pb_x[k] * sgl0_494[k]
                   - f_14 * pc_x[k] * sgl1_494[k];

        t_495[k] = pb_z[k] * sgl0_270[k]
                   - f_14 * pc_z[k] * sgl1_270[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pb_z, pc_y, pc_z, sgl0_273, sgk_216, \
                         sgk_252, sgk_254, sgl1_273, shk_396, shk_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_17 * sgk_252[k]
                   + f_3 * pc_y[k] * shk_396[k];

        t_497[k] = f_15 * sgk_216[k]
                   + f_3 * pc_z[k] * shk_396[k];

        t_498[k] = pb_z[k] * sgl0_273[k]
                   - f_14 * pc_z[k] * sgl1_273[k];

        t_499[k] = f_17 * sgk_254[k]
                   + f_3 * pc_y[k] * shk_398[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, pb_x, pb_z, pc_x, pc_z, sgl0_276, sgl0_500, \
                         sgk_219, sgk_401, sgl1_276, sgl1_500, \
                         shk_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = pb_x[k] * sgl0_500[k]
                   + f_19 * sgk_401[k]
                   - f_14 * pc_x[k] * sgl1_500[k];

        t_501[k] = pb_z[k] * sgl0_276[k]
                   - f_14 * pc_z[k] * sgl1_276[k];

        t_502[k] = f_15 * sgk_219[k]
                   + f_3 * pc_z[k] * shk_399[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, pb_x, pb_z, pc_x, pc_y, pc_z, sgl0_280, \
                         sgl0_504, sgk_257, sgk_405, sgl1_280, sgl1_504, \
                         shk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_17 * sgk_257[k]
                   + f_3 * pc_y[k] * shk_401[k];

        t_504[k] = pb_x[k] * sgl0_504[k]
                   + f_0 * sgk_405[k]
                   - f_14 * pc_x[k] * sgl1_504[k];

        t_505[k] = pb_z[k] * sgl0_280[k]
                   - f_14 * pc_z[k] * sgl1_280[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pb_x, pc_x, pc_y, pc_z, sgl0_507, sgk_222, \
                         sgk_261, sgk_408, sgl1_507, shk_402, shk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_15 * sgk_222[k]
                   + f_3 * pc_z[k] * shk_402[k];

        t_507[k] = pb_x[k] * sgl0_507[k]
                   + f_18 * sgk_408[k]
                   - f_14 * pc_x[k] * sgl1_507[k];

        t_508[k] = f_17 * sgk_261[k]
                   + f_3 * pc_y[k] * shk_405[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pb_x, pb_z, pc_x, pc_z, sgl0_285, sgl0_509, \
                         sgk_226, sgk_410, sgl1_285, sgl1_509, \
                         shk_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = pb_x[k] * sgl0_509[k]
                   + f_18 * sgk_410[k]
                   - f_14 * pc_x[k] * sgl1_509[k];

        t_510[k] = pb_z[k] * sgl0_285[k]
                   - f_14 * pc_z[k] * sgl1_285[k];

        t_511[k] = f_15 * sgk_226[k]
                   + f_3 * pc_z[k] * shk_406[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pb_x, pc_x, pc_y, sgl0_512, sgl0_513, sgk_266, \
                         sgk_413, sgk_414, sgl1_512, sgl1_513, \
                         shk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pb_x[k] * sgl0_512[k]
                   + f_17 * sgk_413[k]
                   - f_14 * pc_x[k] * sgl1_512[k];

        t_513[k] = pb_x[k] * sgl0_513[k]
                   + f_17 * sgk_414[k]
                   - f_14 * pc_x[k] * sgl1_513[k];

        t_514[k] = f_17 * sgk_266[k]
                   + f_3 * pc_y[k] * shk_410[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pb_x, pb_z, pc_x, pc_z, sgl0_291, sgl0_515, \
                         sgk_231, sgk_416, sgl1_291, sgl1_515, \
                         shk_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = pb_x[k] * sgl0_515[k]
                   + f_17 * sgk_416[k]
                   - f_14 * pc_x[k] * sgl1_515[k];

        t_516[k] = pb_z[k] * sgl0_291[k]
                   - f_14 * pc_z[k] * sgl1_291[k];

        t_517[k] = f_15 * sgk_231[k]
                   + f_3 * pc_z[k] * shk_411[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, pb_x, pc_x, sgl0_518, sgl0_519, sgl0_520, \
                         sgk_419, sgk_420, sgk_421, sgl1_518, sgl1_519, \
                         sgl1_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = pb_x[k] * sgl0_518[k]
                   + f_16 * sgk_419[k]
                   - f_14 * pc_x[k] * sgl1_518[k];

        t_519[k] = pb_x[k] * sgl0_519[k]
                   + f_16 * sgk_420[k]
                   - f_14 * pc_x[k] * sgl1_519[k];

        t_520[k] = pb_x[k] * sgl0_520[k]
                   + f_16 * sgk_421[k]
                   - f_14 * pc_x[k] * sgl1_520[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pb_x, pc_x, pc_y, sgl0_522, sgk_272, \
                         sgk_423, sgk_424, sgk_425, sgl1_522, shk_416, shk_424, \
                         shk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_17 * sgk_272[k]
                   + f_3 * pc_y[k] * shk_416[k];

        t_522[k] = pb_x[k] * sgl0_522[k]
                   + f_16 * sgk_423[k]
                   - f_14 * pc_x[k] * sgl1_522[k];

        t_523[k] = f_15 * sgk_424[k]
                   + f_3 * pc_x[k] * shk_424[k];

        t_524[k] = f_15 * sgk_425[k]
                   + f_3 * pc_x[k] * shk_425[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pc_x, sgk_426, sgk_427, sgk_428, \
                         sgk_429, sgk_430, shk_426, shk_427, shk_428, shk_429, \
                         shk_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_15 * sgk_426[k]
                   + f_3 * pc_x[k] * shk_426[k];

        t_526[k] = f_15 * sgk_427[k]
                   + f_3 * pc_x[k] * shk_427[k];

        t_527[k] = f_15 * sgk_428[k]
                   + f_3 * pc_x[k] * shk_428[k];

        t_528[k] = f_15 * sgk_429[k]
                   + f_3 * pc_x[k] * shk_429[k];

        t_529[k] = f_15 * sgk_430[k]
                   + f_3 * pc_x[k] * shk_430[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pb_x, pc_x, pc_z, sgl0_531, sgl0_533, \
                         sgk_244, sgk_431, sgl1_531, sgl1_533, shk_424, \
                         shk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_15 * sgk_431[k]
                   + f_3 * pc_x[k] * shk_431[k];

        t_531[k] = pb_x[k] * sgl0_531[k]
                   - f_14 * pc_x[k] * sgl1_531[k];

        t_532[k] = f_15 * sgk_244[k]
                   + f_3 * pc_z[k] * shk_424[k];

        t_533[k] = pb_x[k] * sgl0_533[k]
                   - f_14 * pc_x[k] * sgl1_533[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, pb_x, pc_x, sgl0_534, sgl0_535, sgl0_536, \
                         sgl0_537, sgl1_534, sgl1_535, sgl1_536, \
                         sgl1_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = pb_x[k] * sgl0_534[k]
                   - f_14 * pc_x[k] * sgl1_534[k];

        t_535[k] = pb_x[k] * sgl0_535[k]
                   - f_14 * pc_x[k] * sgl1_535[k];

        t_536[k] = pb_x[k] * sgl0_536[k]
                   - f_14 * pc_x[k] * sgl1_536[k];

        t_537[k] = pb_x[k] * sgl0_537[k]
                   - f_14 * pc_x[k] * sgl1_537[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, pb_x, pc_x, pc_y, sgl0_539, sgl0_540, \
                         sgk_287, sgk_288, sgk_432, sgl1_539, sgl1_540, shk_431, \
                         shk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_17 * sgk_287[k]
                   + f_3 * pc_y[k] * shk_431[k];

        t_539[k] = pb_x[k] * sgl0_539[k]
                   - f_14 * pc_x[k] * sgl1_539[k];

        t_540[k] = pb_x[k] * sgl0_540[k]
                   + f_20 * sgk_432[k]
                   - f_14 * pc_x[k] * sgl1_540[k];

        t_541[k] = f_16 * sgk_288[k]
                   + f_3 * pc_y[k] * shk_432[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, pb_x, pc_x, pc_y, pc_z, sgl0_543, sgk_252, \
                         sgk_290, sgk_435, sgl1_543, shk_432, shk_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_16 * sgk_252[k]
                   + f_3 * pc_z[k] * shk_432[k];

        t_543[k] = pb_x[k] * sgl0_543[k]
                   + f_19 * sgk_435[k]
                   - f_14 * pc_x[k] * sgl1_543[k];

        t_544[k] = f_16 * sgk_290[k]
                   + f_3 * pc_y[k] * shk_434[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pb_x, pc_x, pc_z, sgl0_545, sgl0_546, sgk_255, \
                         sgk_437, sgk_438, sgl1_545, sgl1_546, \
                         shk_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = pb_x[k] * sgl0_545[k]
                   + f_19 * sgk_437[k]
                   - f_14 * pc_x[k] * sgl1_545[k];

        t_546[k] = pb_x[k] * sgl0_546[k]
                   + f_0 * sgk_438[k]
                   - f_14 * pc_x[k] * sgl1_546[k];

        t_547[k] = f_16 * sgk_255[k]
                   + f_3 * pc_z[k] * shk_435[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, pb_x, pc_x, pc_y, sgl0_549, sgl0_550, sgk_293, \
                         sgk_441, sgk_442, sgl1_549, sgl1_550, \
                         shk_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_16 * sgk_293[k]
                   + f_3 * pc_y[k] * shk_437[k];

        t_549[k] = pb_x[k] * sgl0_549[k]
                   + f_0 * sgk_441[k]
                   - f_14 * pc_x[k] * sgl1_549[k];

        t_550[k] = pb_x[k] * sgl0_550[k]
                   + f_18 * sgk_442[k]
                   - f_14 * pc_x[k] * sgl1_550[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, pb_x, pc_x, pc_y, pc_z, sgl0_552, sgk_258, \
                         sgk_297, sgk_444, sgl1_552, shk_438, shk_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_16 * sgk_258[k]
                   + f_3 * pc_z[k] * shk_438[k];

        t_552[k] = pb_x[k] * sgl0_552[k]
                   + f_18 * sgk_444[k]
                   - f_14 * pc_x[k] * sgl1_552[k];

        t_553[k] = f_16 * sgk_297[k]
                   + f_3 * pc_y[k] * shk_441[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pb_x, pc_x, pc_z, sgl0_554, sgl0_555, sgk_262, \
                         sgk_446, sgk_447, sgl1_554, sgl1_555, \
                         shk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = pb_x[k] * sgl0_554[k]
                   + f_18 * sgk_446[k]
                   - f_14 * pc_x[k] * sgl1_554[k];

        t_555[k] = pb_x[k] * sgl0_555[k]
                   + f_17 * sgk_447[k]
                   - f_14 * pc_x[k] * sgl1_555[k];

        t_556[k] = f_16 * sgk_262[k]
                   + f_3 * pc_z[k] * shk_442[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pb_x, pc_x, pc_y, sgl0_557, sgl0_558, sgk_302, \
                         sgk_449, sgk_450, sgl1_557, sgl1_558, \
                         shk_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = pb_x[k] * sgl0_557[k]
                   + f_17 * sgk_449[k]
                   - f_14 * pc_x[k] * sgl1_557[k];

        t_558[k] = pb_x[k] * sgl0_558[k]
                   + f_17 * sgk_450[k]
                   - f_14 * pc_x[k] * sgl1_558[k];

        t_559[k] = f_16 * sgk_302[k]
                   + f_3 * pc_y[k] * shk_446[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, pb_x, pc_x, pc_z, sgl0_560, sgl0_561, sgk_267, \
                         sgk_452, sgk_453, sgl1_560, sgl1_561, \
                         shk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = pb_x[k] * sgl0_560[k]
                   + f_17 * sgk_452[k]
                   - f_14 * pc_x[k] * sgl1_560[k];

        t_561[k] = pb_x[k] * sgl0_561[k]
                   + f_16 * sgk_453[k]
                   - f_14 * pc_x[k] * sgl1_561[k];

        t_562[k] = f_16 * sgk_267[k]
                   + f_3 * pc_z[k] * shk_447[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pb_x, pc_x, sgl0_563, sgl0_564, sgl0_565, \
                         sgk_455, sgk_456, sgk_457, sgl1_563, sgl1_564, \
                         sgl1_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = pb_x[k] * sgl0_563[k]
                   + f_16 * sgk_455[k]
                   - f_14 * pc_x[k] * sgl1_563[k];

        t_564[k] = pb_x[k] * sgl0_564[k]
                   + f_16 * sgk_456[k]
                   - f_14 * pc_x[k] * sgl1_564[k];

        t_565[k] = pb_x[k] * sgl0_565[k]
                   + f_16 * sgk_457[k]
                   - f_14 * pc_x[k] * sgl1_565[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pb_x, pc_x, pc_y, sgl0_567, sgk_308, \
                         sgk_459, sgk_460, sgk_461, sgl1_567, shk_452, shk_460, \
                         shk_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_16 * sgk_308[k]
                   + f_3 * pc_y[k] * shk_452[k];

        t_567[k] = pb_x[k] * sgl0_567[k]
                   + f_16 * sgk_459[k]
                   - f_14 * pc_x[k] * sgl1_567[k];

        t_568[k] = f_15 * sgk_460[k]
                   + f_3 * pc_x[k] * shk_460[k];

        t_569[k] = f_15 * sgk_461[k]
                   + f_3 * pc_x[k] * shk_461[k];
    }
}

static auto
compute_prim_shl_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgl0,
                                                          const size_t sgk, const size_t sgl1,
                                                          const size_t shi0, const size_t shi1,
                                                          const size_t shk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 4.0 / q;

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
    auto *t_675 = buffer.data(target + 675);
    auto *t_676 = buffer.data(target + 676);
    auto *t_677 = buffer.data(target + 677);
    auto *t_678 = buffer.data(target + 678);
    auto *t_679 = buffer.data(target + 679);
    auto *t_680 = buffer.data(target + 680);
    auto *t_681 = buffer.data(target + 681);
    auto *t_682 = buffer.data(target + 682);
    auto *t_683 = buffer.data(target + 683);
    auto *t_684 = buffer.data(target + 684);
    auto *t_685 = buffer.data(target + 685);
    auto *t_686 = buffer.data(target + 686);
    auto *t_687 = buffer.data(target + 687);
    auto *t_688 = buffer.data(target + 688);
    auto *t_689 = buffer.data(target + 689);
    auto *t_690 = buffer.data(target + 690);
    auto *t_691 = buffer.data(target + 691);
    auto *t_692 = buffer.data(target + 692);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgl0_405 = buffer.data(sgl0 + 405);
    const auto *sgl0_410 = buffer.data(sgl0 + 410);
    const auto *sgl0_414 = buffer.data(sgl0 + 414);
    const auto *sgl0_419 = buffer.data(sgl0 + 419);
    const auto *sgl0_425 = buffer.data(sgl0 + 425);
    const auto *sgl0_432 = buffer.data(sgl0 + 432);
    const auto *sgl0_576 = buffer.data(sgl0 + 576);
    const auto *sgl0_578 = buffer.data(sgl0 + 578);
    const auto *sgl0_579 = buffer.data(sgl0 + 579);
    const auto *sgl0_580 = buffer.data(sgl0 + 580);
    const auto *sgl0_581 = buffer.data(sgl0 + 581);
    const auto *sgl0_582 = buffer.data(sgl0 + 582);
    const auto *sgl0_584 = buffer.data(sgl0 + 584);
    const auto *sgl0_588 = buffer.data(sgl0 + 588);
    const auto *sgl0_591 = buffer.data(sgl0 + 591);
    const auto *sgl0_595 = buffer.data(sgl0 + 595);
    const auto *sgl0_597 = buffer.data(sgl0 + 597);
    const auto *sgl0_600 = buffer.data(sgl0 + 600);
    const auto *sgl0_602 = buffer.data(sgl0 + 602);
    const auto *sgl0_603 = buffer.data(sgl0 + 603);
    const auto *sgl0_606 = buffer.data(sgl0 + 606);
    const auto *sgl0_608 = buffer.data(sgl0 + 608);
    const auto *sgl0_609 = buffer.data(sgl0 + 609);
    const auto *sgl0_610 = buffer.data(sgl0 + 610);
    const auto *sgl0_621 = buffer.data(sgl0 + 621);
    const auto *sgl0_623 = buffer.data(sgl0 + 623);
    const auto *sgl0_624 = buffer.data(sgl0 + 624);
    const auto *sgl0_625 = buffer.data(sgl0 + 625);
    const auto *sgl0_626 = buffer.data(sgl0 + 626);
    const auto *sgl0_627 = buffer.data(sgl0 + 627);
    const auto *sgl0_629 = buffer.data(sgl0 + 629);
    const auto *sgl0_630 = buffer.data(sgl0 + 630);
    const auto *sgl0_633 = buffer.data(sgl0 + 633);
    const auto *sgl0_635 = buffer.data(sgl0 + 635);
    const auto *sgl0_636 = buffer.data(sgl0 + 636);
    const auto *sgl0_639 = buffer.data(sgl0 + 639);
    const auto *sgl0_640 = buffer.data(sgl0 + 640);
    const auto *sgl0_642 = buffer.data(sgl0 + 642);
    const auto *sgl0_644 = buffer.data(sgl0 + 644);
    const auto *sgl0_645 = buffer.data(sgl0 + 645);
    const auto *sgl0_647 = buffer.data(sgl0 + 647);
    const auto *sgl0_648 = buffer.data(sgl0 + 648);
    const auto *sgl0_650 = buffer.data(sgl0 + 650);
    const auto *sgl0_651 = buffer.data(sgl0 + 651);
    const auto *sgl0_653 = buffer.data(sgl0 + 653);
    const auto *sgl0_654 = buffer.data(sgl0 + 654);
    const auto *sgl0_655 = buffer.data(sgl0 + 655);
    const auto *sgl0_657 = buffer.data(sgl0 + 657);
    const auto *sgl0_666 = buffer.data(sgl0 + 666);
    const auto *sgl0_668 = buffer.data(sgl0 + 668);
    const auto *sgl0_669 = buffer.data(sgl0 + 669);
    const auto *sgl0_670 = buffer.data(sgl0 + 670);
    const auto *sgl0_671 = buffer.data(sgl0 + 671);
    const auto *sgl0_672 = buffer.data(sgl0 + 672);
    const auto *sgl0_674 = buffer.data(sgl0 + 674);

    const auto *sgk_280 = buffer.data(sgk + 280);
    const auto *sgk_288 = buffer.data(sgk + 288);
    const auto *sgk_291 = buffer.data(sgk + 291);
    const auto *sgk_294 = buffer.data(sgk + 294);
    const auto *sgk_298 = buffer.data(sgk + 298);
    const auto *sgk_303 = buffer.data(sgk + 303);
    const auto *sgk_316 = buffer.data(sgk + 316);
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
    const auto *sgk_359 = buffer.data(sgk + 359);
    const auto *sgk_360 = buffer.data(sgk + 360);
    const auto *sgk_362 = buffer.data(sgk + 362);
    const auto *sgk_365 = buffer.data(sgk + 365);
    const auto *sgk_369 = buffer.data(sgk + 369);
    const auto *sgk_462 = buffer.data(sgk + 462);
    const auto *sgk_463 = buffer.data(sgk + 463);
    const auto *sgk_464 = buffer.data(sgk + 464);
    const auto *sgk_465 = buffer.data(sgk + 465);
    const auto *sgk_466 = buffer.data(sgk + 466);
    const auto *sgk_467 = buffer.data(sgk + 467);
    const auto *sgk_471 = buffer.data(sgk + 471);
    const auto *sgk_474 = buffer.data(sgk + 474);
    const auto *sgk_478 = buffer.data(sgk + 478);
    const auto *sgk_480 = buffer.data(sgk + 480);
    const auto *sgk_483 = buffer.data(sgk + 483);
    const auto *sgk_485 = buffer.data(sgk + 485);
    const auto *sgk_486 = buffer.data(sgk + 486);
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

    const auto *sgl1_405 = buffer.data(sgl1 + 405);
    const auto *sgl1_410 = buffer.data(sgl1 + 410);
    const auto *sgl1_414 = buffer.data(sgl1 + 414);
    const auto *sgl1_419 = buffer.data(sgl1 + 419);
    const auto *sgl1_425 = buffer.data(sgl1 + 425);
    const auto *sgl1_432 = buffer.data(sgl1 + 432);
    const auto *sgl1_576 = buffer.data(sgl1 + 576);
    const auto *sgl1_578 = buffer.data(sgl1 + 578);
    const auto *sgl1_579 = buffer.data(sgl1 + 579);
    const auto *sgl1_580 = buffer.data(sgl1 + 580);
    const auto *sgl1_581 = buffer.data(sgl1 + 581);
    const auto *sgl1_582 = buffer.data(sgl1 + 582);
    const auto *sgl1_584 = buffer.data(sgl1 + 584);
    const auto *sgl1_588 = buffer.data(sgl1 + 588);
    const auto *sgl1_591 = buffer.data(sgl1 + 591);
    const auto *sgl1_595 = buffer.data(sgl1 + 595);
    const auto *sgl1_597 = buffer.data(sgl1 + 597);
    const auto *sgl1_600 = buffer.data(sgl1 + 600);
    const auto *sgl1_602 = buffer.data(sgl1 + 602);
    const auto *sgl1_603 = buffer.data(sgl1 + 603);
    const auto *sgl1_606 = buffer.data(sgl1 + 606);
    const auto *sgl1_608 = buffer.data(sgl1 + 608);
    const auto *sgl1_609 = buffer.data(sgl1 + 609);
    const auto *sgl1_610 = buffer.data(sgl1 + 610);
    const auto *sgl1_621 = buffer.data(sgl1 + 621);
    const auto *sgl1_623 = buffer.data(sgl1 + 623);
    const auto *sgl1_624 = buffer.data(sgl1 + 624);
    const auto *sgl1_625 = buffer.data(sgl1 + 625);
    const auto *sgl1_626 = buffer.data(sgl1 + 626);
    const auto *sgl1_627 = buffer.data(sgl1 + 627);
    const auto *sgl1_629 = buffer.data(sgl1 + 629);
    const auto *sgl1_630 = buffer.data(sgl1 + 630);
    const auto *sgl1_633 = buffer.data(sgl1 + 633);
    const auto *sgl1_635 = buffer.data(sgl1 + 635);
    const auto *sgl1_636 = buffer.data(sgl1 + 636);
    const auto *sgl1_639 = buffer.data(sgl1 + 639);
    const auto *sgl1_640 = buffer.data(sgl1 + 640);
    const auto *sgl1_642 = buffer.data(sgl1 + 642);
    const auto *sgl1_644 = buffer.data(sgl1 + 644);
    const auto *sgl1_645 = buffer.data(sgl1 + 645);
    const auto *sgl1_647 = buffer.data(sgl1 + 647);
    const auto *sgl1_648 = buffer.data(sgl1 + 648);
    const auto *sgl1_650 = buffer.data(sgl1 + 650);
    const auto *sgl1_651 = buffer.data(sgl1 + 651);
    const auto *sgl1_653 = buffer.data(sgl1 + 653);
    const auto *sgl1_654 = buffer.data(sgl1 + 654);
    const auto *sgl1_655 = buffer.data(sgl1 + 655);
    const auto *sgl1_657 = buffer.data(sgl1 + 657);
    const auto *sgl1_666 = buffer.data(sgl1 + 666);
    const auto *sgl1_668 = buffer.data(sgl1 + 668);
    const auto *sgl1_669 = buffer.data(sgl1 + 669);
    const auto *sgl1_670 = buffer.data(sgl1 + 670);
    const auto *sgl1_671 = buffer.data(sgl1 + 671);
    const auto *sgl1_672 = buffer.data(sgl1 + 672);
    const auto *sgl1_674 = buffer.data(sgl1 + 674);

    const auto *shi0_420 = buffer.data(shi0 + 420);
    const auto *shi0_423 = buffer.data(shi0 + 423);
    const auto *shi0_425 = buffer.data(shi0 + 425);
    const auto *shi0_426 = buffer.data(shi0 + 426);
    const auto *shi0_429 = buffer.data(shi0 + 429);
    const auto *shi0_430 = buffer.data(shi0 + 430);
    const auto *shi0_432 = buffer.data(shi0 + 432);
    const auto *shi0_434 = buffer.data(shi0 + 434);
    const auto *shi0_435 = buffer.data(shi0 + 435);
    const auto *shi0_437 = buffer.data(shi0 + 437);

    const auto *shi1_420 = buffer.data(shi1 + 420);
    const auto *shi1_423 = buffer.data(shi1 + 423);
    const auto *shi1_425 = buffer.data(shi1 + 425);
    const auto *shi1_426 = buffer.data(shi1 + 426);
    const auto *shi1_429 = buffer.data(shi1 + 429);
    const auto *shi1_430 = buffer.data(shi1 + 430);
    const auto *shi1_432 = buffer.data(shi1 + 432);
    const auto *shi1_434 = buffer.data(shi1 + 434);
    const auto *shi1_435 = buffer.data(shi1 + 435);
    const auto *shi1_437 = buffer.data(shi1 + 437);

    const auto *shk_460 = buffer.data(shk + 460);
    const auto *shk_462 = buffer.data(shk + 462);
    const auto *shk_463 = buffer.data(shk + 463);
    const auto *shk_464 = buffer.data(shk + 464);
    const auto *shk_465 = buffer.data(shk + 465);
    const auto *shk_466 = buffer.data(shk + 466);
    const auto *shk_467 = buffer.data(shk + 467);
    const auto *shk_468 = buffer.data(shk + 468);
    const auto *shk_470 = buffer.data(shk + 470);
    const auto *shk_471 = buffer.data(shk + 471);
    const auto *shk_473 = buffer.data(shk + 473);
    const auto *shk_474 = buffer.data(shk + 474);
    const auto *shk_477 = buffer.data(shk + 477);
    const auto *shk_478 = buffer.data(shk + 478);
    const auto *shk_482 = buffer.data(shk + 482);
    const auto *shk_483 = buffer.data(shk + 483);
    const auto *shk_488 = buffer.data(shk + 488);
    const auto *shk_496 = buffer.data(shk + 496);
    const auto *shk_497 = buffer.data(shk + 497);
    const auto *shk_498 = buffer.data(shk + 498);
    const auto *shk_499 = buffer.data(shk + 499);
    const auto *shk_500 = buffer.data(shk + 500);
    const auto *shk_501 = buffer.data(shk + 501);
    const auto *shk_502 = buffer.data(shk + 502);
    const auto *shk_503 = buffer.data(shk + 503);
    const auto *shk_504 = buffer.data(shk + 504);
    const auto *shk_506 = buffer.data(shk + 506);
    const auto *shk_507 = buffer.data(shk + 507);
    const auto *shk_509 = buffer.data(shk + 509);
    const auto *shk_510 = buffer.data(shk + 510);
    const auto *shk_513 = buffer.data(shk + 513);
    const auto *shk_514 = buffer.data(shk + 514);
    const auto *shk_518 = buffer.data(shk + 518);
    const auto *shk_519 = buffer.data(shk + 519);
    const auto *shk_524 = buffer.data(shk + 524);
    const auto *shk_532 = buffer.data(shk + 532);
    const auto *shk_533 = buffer.data(shk + 533);
    const auto *shk_534 = buffer.data(shk + 534);
    const auto *shk_535 = buffer.data(shk + 535);
    const auto *shk_536 = buffer.data(shk + 536);
    const auto *shk_537 = buffer.data(shk + 537);
    const auto *shk_538 = buffer.data(shk + 538);
    const auto *shk_539 = buffer.data(shk + 539);
    const auto *shk_540 = buffer.data(shk + 540);
    const auto *shk_542 = buffer.data(shk + 542);
    const auto *shk_543 = buffer.data(shk + 543);
    const auto *shk_545 = buffer.data(shk + 545);
    const auto *shk_546 = buffer.data(shk + 546);
    const auto *shk_549 = buffer.data(shk + 549);
    const auto *shk_550 = buffer.data(shk + 550);
    const auto *shk_552 = buffer.data(shk + 552);
    const auto *shk_554 = buffer.data(shk + 554);
    const auto *shk_555 = buffer.data(shk + 555);
    const auto *shk_557 = buffer.data(shk + 557);

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, pc_x, sgk_462, sgk_463, sgk_464, \
                         sgk_465, sgk_466, shk_462, shk_463, shk_464, shk_465, \
                         shk_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_15 * sgk_462[k]
                   + f_3 * pc_x[k] * shk_462[k];

        t_571[k] = f_15 * sgk_463[k]
                   + f_3 * pc_x[k] * shk_463[k];

        t_572[k] = f_15 * sgk_464[k]
                   + f_3 * pc_x[k] * shk_464[k];

        t_573[k] = f_15 * sgk_465[k]
                   + f_3 * pc_x[k] * shk_465[k];

        t_574[k] = f_15 * sgk_466[k]
                   + f_3 * pc_x[k] * shk_466[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, pb_x, pc_x, pc_z, sgl0_576, sgl0_578, \
                         sgk_280, sgk_467, sgl1_576, sgl1_578, shk_460, \
                         shk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_15 * sgk_467[k]
                   + f_3 * pc_x[k] * shk_467[k];

        t_576[k] = pb_x[k] * sgl0_576[k]
                   - f_14 * pc_x[k] * sgl1_576[k];

        t_577[k] = f_16 * sgk_280[k]
                   + f_3 * pc_z[k] * shk_460[k];

        t_578[k] = pb_x[k] * sgl0_578[k]
                   - f_14 * pc_x[k] * sgl1_578[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, t_582, pb_x, pc_x, sgl0_579, sgl0_580, sgl0_581, \
                         sgl0_582, sgl1_579, sgl1_580, sgl1_581, \
                         sgl1_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = pb_x[k] * sgl0_579[k]
                   - f_14 * pc_x[k] * sgl1_579[k];

        t_580[k] = pb_x[k] * sgl0_580[k]
                   - f_14 * pc_x[k] * sgl1_580[k];

        t_581[k] = pb_x[k] * sgl0_581[k]
                   - f_14 * pc_x[k] * sgl1_581[k];

        t_582[k] = pb_x[k] * sgl0_582[k]
                   - f_14 * pc_x[k] * sgl1_582[k];
    }

#pragma omp simd aligned(t_583, t_584, t_585, t_586, pb_x, pb_y, pc_x, pc_y, sgl0_405, \
                         sgl0_584, sgk_323, sgk_324, sgl1_405, sgl1_584, shk_467, \
                         shk_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_583[k] = f_16 * sgk_323[k]
                   + f_3 * pc_y[k] * shk_467[k];

        t_584[k] = pb_x[k] * sgl0_584[k]
                   - f_14 * pc_x[k] * sgl1_584[k];

        t_585[k] = pb_y[k] * sgl0_405[k]
                   - f_14 * pc_y[k] * sgl1_405[k];

        t_586[k] = f_15 * sgk_324[k]
                   + f_3 * pc_y[k] * shk_468[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pb_x, pc_x, pc_y, pc_z, sgl0_588, sgk_288, \
                         sgk_326, sgk_471, sgl1_588, shk_468, shk_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_17 * sgk_288[k]
                   + f_3 * pc_z[k] * shk_468[k];

        t_588[k] = pb_x[k] * sgl0_588[k]
                   + f_19 * sgk_471[k]
                   - f_14 * pc_x[k] * sgl1_588[k];

        t_589[k] = f_15 * sgk_326[k]
                   + f_3 * pc_y[k] * shk_470[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, pb_x, pb_y, pc_x, pc_y, pc_z, sgl0_410, \
                         sgl0_591, sgk_291, sgk_474, sgl1_410, sgl1_591, \
                         shk_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = pb_y[k] * sgl0_410[k]
                   - f_14 * pc_y[k] * sgl1_410[k];

        t_591[k] = pb_x[k] * sgl0_591[k]
                   + f_0 * sgk_474[k]
                   - f_14 * pc_x[k] * sgl1_591[k];

        t_592[k] = f_17 * sgk_291[k]
                   + f_3 * pc_z[k] * shk_471[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pb_x, pb_y, pc_x, pc_y, sgl0_414, sgl0_595, \
                         sgk_329, sgk_478, sgl1_414, sgl1_595, \
                         shk_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = f_15 * sgk_329[k]
                   + f_3 * pc_y[k] * shk_473[k];

        t_594[k] = pb_y[k] * sgl0_414[k]
                   - f_14 * pc_y[k] * sgl1_414[k];

        t_595[k] = pb_x[k] * sgl0_595[k]
                   + f_18 * sgk_478[k]
                   - f_14 * pc_x[k] * sgl1_595[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, pb_x, pc_x, pc_y, pc_z, sgl0_597, sgk_294, \
                         sgk_333, sgk_480, sgl1_597, shk_474, shk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_17 * sgk_294[k]
                   + f_3 * pc_z[k] * shk_474[k];

        t_597[k] = pb_x[k] * sgl0_597[k]
                   + f_18 * sgk_480[k]
                   - f_14 * pc_x[k] * sgl1_597[k];

        t_598[k] = f_15 * sgk_333[k]
                   + f_3 * pc_y[k] * shk_477[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, pb_x, pb_y, pc_x, pc_y, pc_z, sgl0_419, \
                         sgl0_600, sgk_298, sgk_483, sgl1_419, sgl1_600, \
                         shk_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = pb_y[k] * sgl0_419[k]
                   - f_14 * pc_y[k] * sgl1_419[k];

        t_600[k] = pb_x[k] * sgl0_600[k]
                   + f_17 * sgk_483[k]
                   - f_14 * pc_x[k] * sgl1_600[k];

        t_601[k] = f_17 * sgk_298[k]
                   + f_3 * pc_z[k] * shk_478[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, pb_x, pc_x, pc_y, sgl0_602, sgl0_603, sgk_338, \
                         sgk_485, sgk_486, sgl1_602, sgl1_603, \
                         shk_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = pb_x[k] * sgl0_602[k]
                   + f_17 * sgk_485[k]
                   - f_14 * pc_x[k] * sgl1_602[k];

        t_603[k] = pb_x[k] * sgl0_603[k]
                   + f_17 * sgk_486[k]
                   - f_14 * pc_x[k] * sgl1_603[k];

        t_604[k] = f_15 * sgk_338[k]
                   + f_3 * pc_y[k] * shk_482[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, pb_x, pb_y, pc_x, pc_y, pc_z, sgl0_425, \
                         sgl0_606, sgk_303, sgk_489, sgl1_425, sgl1_606, \
                         shk_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = pb_y[k] * sgl0_425[k]
                   - f_14 * pc_y[k] * sgl1_425[k];

        t_606[k] = pb_x[k] * sgl0_606[k]
                   + f_16 * sgk_489[k]
                   - f_14 * pc_x[k] * sgl1_606[k];

        t_607[k] = f_17 * sgk_303[k]
                   + f_3 * pc_z[k] * shk_483[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, pb_x, pc_x, sgl0_608, sgl0_609, sgl0_610, \
                         sgk_491, sgk_492, sgk_493, sgl1_608, sgl1_609, \
                         sgl1_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = pb_x[k] * sgl0_608[k]
                   + f_16 * sgk_491[k]
                   - f_14 * pc_x[k] * sgl1_608[k];

        t_609[k] = pb_x[k] * sgl0_609[k]
                   + f_16 * sgk_492[k]
                   - f_14 * pc_x[k] * sgl1_609[k];

        t_610[k] = pb_x[k] * sgl0_610[k]
                   + f_16 * sgk_493[k]
                   - f_14 * pc_x[k] * sgl1_610[k];
    }

#pragma omp simd aligned(t_611, t_612, t_613, t_614, pb_y, pc_x, pc_y, sgl0_432, sgk_344, \
                         sgk_496, sgk_497, sgl1_432, shk_488, shk_496, \
                         shk_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = f_15 * sgk_344[k]
                   + f_3 * pc_y[k] * shk_488[k];

        t_612[k] = pb_y[k] * sgl0_432[k]
                   - f_14 * pc_y[k] * sgl1_432[k];

        t_613[k] = f_15 * sgk_496[k]
                   + f_3 * pc_x[k] * shk_496[k];

        t_614[k] = f_15 * sgk_497[k]
                   + f_3 * pc_x[k] * shk_497[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, pc_x, sgk_498, sgk_499, sgk_500, \
                         sgk_501, sgk_502, shk_498, shk_499, shk_500, shk_501, \
                         shk_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_15 * sgk_498[k]
                   + f_3 * pc_x[k] * shk_498[k];

        t_616[k] = f_15 * sgk_499[k]
                   + f_3 * pc_x[k] * shk_499[k];

        t_617[k] = f_15 * sgk_500[k]
                   + f_3 * pc_x[k] * shk_500[k];

        t_618[k] = f_15 * sgk_501[k]
                   + f_3 * pc_x[k] * shk_501[k];

        t_619[k] = f_15 * sgk_502[k]
                   + f_3 * pc_x[k] * shk_502[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, pb_x, pc_x, pc_z, sgl0_621, sgl0_623, \
                         sgk_316, sgk_503, sgl1_621, sgl1_623, shk_496, \
                         shk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_15 * sgk_503[k]
                   + f_3 * pc_x[k] * shk_503[k];

        t_621[k] = pb_x[k] * sgl0_621[k]
                   - f_14 * pc_x[k] * sgl1_621[k];

        t_622[k] = f_17 * sgk_316[k]
                   + f_3 * pc_z[k] * shk_496[k];

        t_623[k] = pb_x[k] * sgl0_623[k]
                   - f_14 * pc_x[k] * sgl1_623[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, pb_x, pc_x, sgl0_624, sgl0_625, sgl0_626, \
                         sgl0_627, sgl1_624, sgl1_625, sgl1_626, \
                         sgl1_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = pb_x[k] * sgl0_624[k]
                   - f_14 * pc_x[k] * sgl1_624[k];

        t_625[k] = pb_x[k] * sgl0_625[k]
                   - f_14 * pc_x[k] * sgl1_625[k];

        t_626[k] = pb_x[k] * sgl0_626[k]
                   - f_14 * pc_x[k] * sgl1_626[k];

        t_627[k] = pb_x[k] * sgl0_627[k]
                   - f_14 * pc_x[k] * sgl1_627[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, t_631, pb_x, pc_x, pc_y, sgl0_629, sgl0_630, \
                         sgk_359, sgk_504, sgl1_629, sgl1_630, shk_503, \
                         shk_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_15 * sgk_359[k]
                   + f_3 * pc_y[k] * shk_503[k];

        t_629[k] = pb_x[k] * sgl0_629[k]
                   - f_14 * pc_x[k] * sgl1_629[k];

        t_630[k] = pb_x[k] * sgl0_630[k]
                   + f_20 * sgk_504[k]
                   - f_14 * pc_x[k] * sgl1_630[k];

        t_631[k] = f_3 * pc_y[k] * shk_504[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, pb_x, pc_x, pc_y, pc_z, sgl0_633, sgk_324, \
                         sgk_507, sgl1_633, shk_504, shk_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_18 * sgk_324[k]
                   + f_3 * pc_z[k] * shk_504[k];

        t_633[k] = pb_x[k] * sgl0_633[k]
                   + f_19 * sgk_507[k]
                   - f_14 * pc_x[k] * sgl1_633[k];

        t_634[k] = f_3 * pc_y[k] * shk_506[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, pb_x, pc_x, pc_z, sgl0_635, sgl0_636, sgk_327, \
                         sgk_509, sgk_510, sgl1_635, sgl1_636, \
                         shk_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = pb_x[k] * sgl0_635[k]
                   + f_19 * sgk_509[k]
                   - f_14 * pc_x[k] * sgl1_635[k];

        t_636[k] = pb_x[k] * sgl0_636[k]
                   + f_0 * sgk_510[k]
                   - f_14 * pc_x[k] * sgl1_636[k];

        t_637[k] = f_18 * sgk_327[k]
                   + f_3 * pc_z[k] * shk_507[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, pb_x, pc_x, pc_y, sgl0_639, sgl0_640, sgk_513, \
                         sgk_514, sgl1_639, sgl1_640, shk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_3 * pc_y[k] * shk_509[k];

        t_639[k] = pb_x[k] * sgl0_639[k]
                   + f_0 * sgk_513[k]
                   - f_14 * pc_x[k] * sgl1_639[k];

        t_640[k] = pb_x[k] * sgl0_640[k]
                   + f_18 * sgk_514[k]
                   - f_14 * pc_x[k] * sgl1_640[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, pb_x, pc_x, pc_y, pc_z, sgl0_642, sgk_330, \
                         sgk_516, sgl1_642, shk_510, shk_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_18 * sgk_330[k]
                   + f_3 * pc_z[k] * shk_510[k];

        t_642[k] = pb_x[k] * sgl0_642[k]
                   + f_18 * sgk_516[k]
                   - f_14 * pc_x[k] * sgl1_642[k];

        t_643[k] = f_3 * pc_y[k] * shk_513[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, pb_x, pc_x, pc_z, sgl0_644, sgl0_645, sgk_334, \
                         sgk_518, sgk_519, sgl1_644, sgl1_645, \
                         shk_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = pb_x[k] * sgl0_644[k]
                   + f_18 * sgk_518[k]
                   - f_14 * pc_x[k] * sgl1_644[k];

        t_645[k] = pb_x[k] * sgl0_645[k]
                   + f_17 * sgk_519[k]
                   - f_14 * pc_x[k] * sgl1_645[k];

        t_646[k] = f_18 * sgk_334[k]
                   + f_3 * pc_z[k] * shk_514[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pb_x, pc_x, pc_y, sgl0_647, sgl0_648, sgk_521, \
                         sgk_522, sgl1_647, sgl1_648, shk_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = pb_x[k] * sgl0_647[k]
                   + f_17 * sgk_521[k]
                   - f_14 * pc_x[k] * sgl1_647[k];

        t_648[k] = pb_x[k] * sgl0_648[k]
                   + f_17 * sgk_522[k]
                   - f_14 * pc_x[k] * sgl1_648[k];

        t_649[k] = f_3 * pc_y[k] * shk_518[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pb_x, pc_x, pc_z, sgl0_650, sgl0_651, sgk_339, \
                         sgk_524, sgk_525, sgl1_650, sgl1_651, \
                         shk_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = pb_x[k] * sgl0_650[k]
                   + f_17 * sgk_524[k]
                   - f_14 * pc_x[k] * sgl1_650[k];

        t_651[k] = pb_x[k] * sgl0_651[k]
                   + f_16 * sgk_525[k]
                   - f_14 * pc_x[k] * sgl1_651[k];

        t_652[k] = f_18 * sgk_339[k]
                   + f_3 * pc_z[k] * shk_519[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pb_x, pc_x, sgl0_653, sgl0_654, sgl0_655, \
                         sgk_527, sgk_528, sgk_529, sgl1_653, sgl1_654, \
                         sgl1_655 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = pb_x[k] * sgl0_653[k]
                   + f_16 * sgk_527[k]
                   - f_14 * pc_x[k] * sgl1_653[k];

        t_654[k] = pb_x[k] * sgl0_654[k]
                   + f_16 * sgk_528[k]
                   - f_14 * pc_x[k] * sgl1_654[k];

        t_655[k] = pb_x[k] * sgl0_655[k]
                   + f_16 * sgk_529[k]
                   - f_14 * pc_x[k] * sgl1_655[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, t_659, pb_x, pc_x, pc_y, sgl0_657, sgk_531, \
                         sgk_532, sgk_533, sgl1_657, shk_524, shk_532, \
                         shk_533 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_3 * pc_y[k] * shk_524[k];

        t_657[k] = pb_x[k] * sgl0_657[k]
                   + f_16 * sgk_531[k]
                   - f_14 * pc_x[k] * sgl1_657[k];

        t_658[k] = f_15 * sgk_532[k]
                   + f_3 * pc_x[k] * shk_532[k];

        t_659[k] = f_15 * sgk_533[k]
                   + f_3 * pc_x[k] * shk_533[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, pc_x, sgk_534, sgk_535, sgk_536, \
                         sgk_537, sgk_538, shk_534, shk_535, shk_536, shk_537, \
                         shk_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_15 * sgk_534[k]
                   + f_3 * pc_x[k] * shk_534[k];

        t_661[k] = f_15 * sgk_535[k]
                   + f_3 * pc_x[k] * shk_535[k];

        t_662[k] = f_15 * sgk_536[k]
                   + f_3 * pc_x[k] * shk_536[k];

        t_663[k] = f_15 * sgk_537[k]
                   + f_3 * pc_x[k] * shk_537[k];

        t_664[k] = f_15 * sgk_538[k]
                   + f_3 * pc_x[k] * shk_538[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, pb_x, pc_x, pc_z, sgl0_666, sgl0_668, \
                         sgk_352, sgk_539, sgl1_666, sgl1_668, shk_532, \
                         shk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_15 * sgk_539[k]
                   + f_3 * pc_x[k] * shk_539[k];

        t_666[k] = pb_x[k] * sgl0_666[k]
                   - f_14 * pc_x[k] * sgl1_666[k];

        t_667[k] = f_18 * sgk_352[k]
                   + f_3 * pc_z[k] * shk_532[k];

        t_668[k] = pb_x[k] * sgl0_668[k]
                   - f_14 * pc_x[k] * sgl1_668[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, t_672, pb_x, pc_x, sgl0_669, sgl0_670, sgl0_671, \
                         sgl0_672, sgl1_669, sgl1_670, sgl1_671, \
                         sgl1_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = pb_x[k] * sgl0_669[k]
                   - f_14 * pc_x[k] * sgl1_669[k];

        t_670[k] = pb_x[k] * sgl0_670[k]
                   - f_14 * pc_x[k] * sgl1_670[k];

        t_671[k] = pb_x[k] * sgl0_671[k]
                   - f_14 * pc_x[k] * sgl1_671[k];

        t_672[k] = pb_x[k] * sgl0_672[k]
                   - f_14 * pc_x[k] * sgl1_672[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, t_677, pb_x, pc_x, pc_y, pc_z, sgl0_674, \
                         sgk_360, sgl1_674, shi0_420, shi1_420, shk_539, \
                         shk_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_3 * pc_y[k] * shk_539[k];

        t_674[k] = pb_x[k] * sgl0_674[k]
                   - f_14 * pc_x[k] * sgl1_674[k];

        t_675[k] = f_1 * shi0_420[k]
                   - f_2 * shi1_420[k]
                   + f_3 * pc_x[k] * shk_540[k];

        t_676[k] = f_0 * sgk_360[k]
                   + f_3 * pc_y[k] * shk_540[k];

        t_677[k] = f_3 * pc_z[k] * shk_540[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, pc_x, pc_y, sgk_362, shi0_423, shi0_425, \
                         shi1_423, shi1_425, shk_542, shk_543, \
                         shk_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = f_4 * shi0_423[k]
                   - f_5 * shi1_423[k]
                   + f_3 * pc_x[k] * shk_543[k];

        t_679[k] = f_0 * sgk_362[k]
                   + f_3 * pc_y[k] * shk_542[k];

        t_680[k] = f_4 * shi0_425[k]
                   - f_5 * shi1_425[k]
                   + f_3 * pc_x[k] * shk_545[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, t_684, pc_x, pc_y, pc_z, sgk_365, shi0_426, \
                         shi0_429, shi1_426, shi1_429, shk_543, shk_545, shk_546, \
                         shk_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = f_6 * shi0_426[k]
                   - f_7 * shi1_426[k]
                   + f_3 * pc_x[k] * shk_546[k];

        t_682[k] = f_3 * pc_z[k] * shk_543[k];

        t_683[k] = f_0 * sgk_365[k]
                   + f_3 * pc_y[k] * shk_545[k];

        t_684[k] = f_6 * shi0_429[k]
                   - f_7 * shi1_429[k]
                   + f_3 * pc_x[k] * shk_549[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, pc_x, pc_y, pc_z, sgk_369, shi0_430, \
                         shi0_432, shi1_430, shi1_432, shk_546, shk_549, shk_550, \
                         shk_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = f_8 * shi0_430[k]
                   - f_9 * shi1_430[k]
                   + f_3 * pc_x[k] * shk_550[k];

        t_686[k] = f_3 * pc_z[k] * shk_546[k];

        t_687[k] = f_8 * shi0_432[k]
                   - f_9 * shi1_432[k]
                   + f_3 * pc_x[k] * shk_552[k];

        t_688[k] = f_0 * sgk_369[k]
                   + f_3 * pc_y[k] * shk_549[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, t_692, pc_x, pc_z, shi0_434, shi0_435, shi0_437, \
                         shi1_434, shi1_435, shi1_437, shk_550, shk_554, shk_555, \
                         shk_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = f_8 * shi0_434[k]
                   - f_9 * shi1_434[k]
                   + f_3 * pc_x[k] * shk_554[k];

        t_690[k] = f_10 * shi0_435[k]
                   - f_11 * shi1_435[k]
                   + f_3 * pc_x[k] * shk_555[k];

        t_691[k] = f_3 * pc_z[k] * shk_550[k];

        t_692[k] = f_10 * shi0_437[k]
                   - f_11 * shi1_437[k]
                   + f_3 * pc_x[k] * shk_557[k];
    }
}

static auto
compute_prim_shl_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgl0,
                                                          const size_t sgk, const size_t sgl1,
                                                          const size_t shi0, const size_t shi1,
                                                          const size_t shk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;

    auto *t_693 = buffer.data(target + 693);
    auto *t_694 = buffer.data(target + 694);
    auto *t_695 = buffer.data(target + 695);
    auto *t_696 = buffer.data(target + 696);
    auto *t_697 = buffer.data(target + 697);
    auto *t_698 = buffer.data(target + 698);
    auto *t_699 = buffer.data(target + 699);
    auto *t_700 = buffer.data(target + 700);
    auto *t_701 = buffer.data(target + 701);
    auto *t_702 = buffer.data(target + 702);
    auto *t_703 = buffer.data(target + 703);
    auto *t_704 = buffer.data(target + 704);
    auto *t_705 = buffer.data(target + 705);
    auto *t_706 = buffer.data(target + 706);
    auto *t_707 = buffer.data(target + 707);
    auto *t_708 = buffer.data(target + 708);
    auto *t_709 = buffer.data(target + 709);
    auto *t_710 = buffer.data(target + 710);
    auto *t_711 = buffer.data(target + 711);
    auto *t_712 = buffer.data(target + 712);
    auto *t_713 = buffer.data(target + 713);
    auto *t_714 = buffer.data(target + 714);
    auto *t_715 = buffer.data(target + 715);
    auto *t_716 = buffer.data(target + 716);
    auto *t_717 = buffer.data(target + 717);
    auto *t_718 = buffer.data(target + 718);
    auto *t_719 = buffer.data(target + 719);
    auto *t_720 = buffer.data(target + 720);
    auto *t_721 = buffer.data(target + 721);
    auto *t_722 = buffer.data(target + 722);
    auto *t_723 = buffer.data(target + 723);
    auto *t_724 = buffer.data(target + 724);
    auto *t_725 = buffer.data(target + 725);
    auto *t_726 = buffer.data(target + 726);
    auto *t_727 = buffer.data(target + 727);
    auto *t_728 = buffer.data(target + 728);
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);
    auto *t_731 = buffer.data(target + 731);
    auto *t_732 = buffer.data(target + 732);
    auto *t_733 = buffer.data(target + 733);
    auto *t_734 = buffer.data(target + 734);
    auto *t_735 = buffer.data(target + 735);
    auto *t_736 = buffer.data(target + 736);
    auto *t_737 = buffer.data(target + 737);
    auto *t_738 = buffer.data(target + 738);
    auto *t_739 = buffer.data(target + 739);
    auto *t_740 = buffer.data(target + 740);
    auto *t_741 = buffer.data(target + 741);
    auto *t_742 = buffer.data(target + 742);
    auto *t_743 = buffer.data(target + 743);
    auto *t_744 = buffer.data(target + 744);
    auto *t_745 = buffer.data(target + 745);
    auto *t_746 = buffer.data(target + 746);
    auto *t_747 = buffer.data(target + 747);
    auto *t_748 = buffer.data(target + 748);
    auto *t_749 = buffer.data(target + 749);
    auto *t_750 = buffer.data(target + 750);
    auto *t_751 = buffer.data(target + 751);
    auto *t_752 = buffer.data(target + 752);
    auto *t_753 = buffer.data(target + 753);
    auto *t_754 = buffer.data(target + 754);
    auto *t_755 = buffer.data(target + 755);
    auto *t_756 = buffer.data(target + 756);
    auto *t_757 = buffer.data(target + 757);
    auto *t_758 = buffer.data(target + 758);
    auto *t_759 = buffer.data(target + 759);
    auto *t_760 = buffer.data(target + 760);
    auto *t_761 = buffer.data(target + 761);
    auto *t_762 = buffer.data(target + 762);
    auto *t_763 = buffer.data(target + 763);
    auto *t_764 = buffer.data(target + 764);
    auto *t_765 = buffer.data(target + 765);
    auto *t_766 = buffer.data(target + 766);
    auto *t_767 = buffer.data(target + 767);
    auto *t_768 = buffer.data(target + 768);
    auto *t_769 = buffer.data(target + 769);
    auto *t_770 = buffer.data(target + 770);
    auto *t_771 = buffer.data(target + 771);
    auto *t_772 = buffer.data(target + 772);
    auto *t_773 = buffer.data(target + 773);
    auto *t_774 = buffer.data(target + 774);
    auto *t_775 = buffer.data(target + 775);
    auto *t_776 = buffer.data(target + 776);
    auto *t_777 = buffer.data(target + 777);
    auto *t_778 = buffer.data(target + 778);
    auto *t_779 = buffer.data(target + 779);
    auto *t_780 = buffer.data(target + 780);
    auto *t_781 = buffer.data(target + 781);
    auto *t_782 = buffer.data(target + 782);
    auto *t_783 = buffer.data(target + 783);
    auto *t_784 = buffer.data(target + 784);
    auto *t_785 = buffer.data(target + 785);
    auto *t_786 = buffer.data(target + 786);
    auto *t_787 = buffer.data(target + 787);
    auto *t_788 = buffer.data(target + 788);
    auto *t_789 = buffer.data(target + 789);
    auto *t_790 = buffer.data(target + 790);
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);
    auto *t_794 = buffer.data(target + 794);
    auto *t_795 = buffer.data(target + 795);
    auto *t_796 = buffer.data(target + 796);
    auto *t_797 = buffer.data(target + 797);
    auto *t_798 = buffer.data(target + 798);
    auto *t_799 = buffer.data(target + 799);
    auto *t_800 = buffer.data(target + 800);
    auto *t_801 = buffer.data(target + 801);
    auto *t_802 = buffer.data(target + 802);
    auto *t_803 = buffer.data(target + 803);
    auto *t_804 = buffer.data(target + 804);
    auto *t_805 = buffer.data(target + 805);
    auto *t_806 = buffer.data(target + 806);
    auto *t_807 = buffer.data(target + 807);
    auto *t_808 = buffer.data(target + 808);
    auto *t_809 = buffer.data(target + 809);
    auto *t_810 = buffer.data(target + 810);
    auto *t_811 = buffer.data(target + 811);
    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);
    auto *t_815 = buffer.data(target + 815);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgl0_450 = buffer.data(sgl0 + 450);
    const auto *sgl0_453 = buffer.data(sgl0 + 453);
    const auto *sgl0_456 = buffer.data(sgl0 + 456);
    const auto *sgl0_460 = buffer.data(sgl0 + 460);
    const auto *sgl0_465 = buffer.data(sgl0 + 465);
    const auto *sgl0_471 = buffer.data(sgl0 + 471);
    const auto *sgl0_486 = buffer.data(sgl0 + 486);
    const auto *sgl0_488 = buffer.data(sgl0 + 488);
    const auto *sgl0_489 = buffer.data(sgl0 + 489);
    const auto *sgl0_490 = buffer.data(sgl0 + 490);
    const auto *sgl0_491 = buffer.data(sgl0 + 491);
    const auto *sgl0_492 = buffer.data(sgl0 + 492);

    const auto *sgk_360 = buffer.data(sgk + 360);
    const auto *sgk_363 = buffer.data(sgk + 363);
    const auto *sgk_366 = buffer.data(sgk + 366);
    const auto *sgk_370 = buffer.data(sgk + 370);
    const auto *sgk_374 = buffer.data(sgk + 374);
    const auto *sgk_375 = buffer.data(sgk + 375);
    const auto *sgk_380 = buffer.data(sgk + 380);
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
    const auto *sgk_410 = buffer.data(sgk + 410);
    const auto *sgk_411 = buffer.data(sgk + 411);
    const auto *sgk_416 = buffer.data(sgk + 416);
    const auto *sgk_424 = buffer.data(sgk + 424);
    const auto *sgk_431 = buffer.data(sgk + 431);
    const auto *sgk_432 = buffer.data(sgk + 432);
    const auto *sgk_434 = buffer.data(sgk + 434);
    const auto *sgk_437 = buffer.data(sgk + 437);
    const auto *sgk_441 = buffer.data(sgk + 441);
    const auto *sgk_446 = buffer.data(sgk + 446);
    const auto *sgk_452 = buffer.data(sgk + 452);
    const auto *sgk_460 = buffer.data(sgk + 460);
    const auto *sgk_462 = buffer.data(sgk + 462);
    const auto *sgk_463 = buffer.data(sgk + 463);
    const auto *sgk_464 = buffer.data(sgk + 464);
    const auto *sgk_465 = buffer.data(sgk + 465);
    const auto *sgk_466 = buffer.data(sgk + 466);
    const auto *sgk_467 = buffer.data(sgk + 467);
    const auto *sgk_468 = buffer.data(sgk + 468);
    const auto *sgk_470 = buffer.data(sgk + 470);

    const auto *sgl1_450 = buffer.data(sgl1 + 450);
    const auto *sgl1_453 = buffer.data(sgl1 + 453);
    const auto *sgl1_456 = buffer.data(sgl1 + 456);
    const auto *sgl1_460 = buffer.data(sgl1 + 460);
    const auto *sgl1_465 = buffer.data(sgl1 + 465);
    const auto *sgl1_471 = buffer.data(sgl1 + 471);
    const auto *sgl1_486 = buffer.data(sgl1 + 486);
    const auto *sgl1_488 = buffer.data(sgl1 + 488);
    const auto *sgl1_489 = buffer.data(sgl1 + 489);
    const auto *sgl1_490 = buffer.data(sgl1 + 490);
    const auto *sgl1_491 = buffer.data(sgl1 + 491);
    const auto *sgl1_492 = buffer.data(sgl1 + 492);

    const auto *shi0_438 = buffer.data(shi0 + 438);
    const auto *shi0_440 = buffer.data(shi0 + 440);
    const auto *shi0_441 = buffer.data(shi0 + 441);
    const auto *shi0_443 = buffer.data(shi0 + 443);
    const auto *shi0_444 = buffer.data(shi0 + 444);
    const auto *shi0_445 = buffer.data(shi0 + 445);
    const auto *shi0_446 = buffer.data(shi0 + 446);
    const auto *shi0_447 = buffer.data(shi0 + 447);
    const auto *shi0_453 = buffer.data(shi0 + 453);
    const auto *shi0_457 = buffer.data(shi0 + 457);
    const auto *shi0_460 = buffer.data(shi0 + 460);
    const auto *shi0_462 = buffer.data(shi0 + 462);
    const auto *shi0_465 = buffer.data(shi0 + 465);
    const auto *shi0_466 = buffer.data(shi0 + 466);
    const auto *shi0_468 = buffer.data(shi0 + 468);
    const auto *shi0_471 = buffer.data(shi0 + 471);
    const auto *shi0_472 = buffer.data(shi0 + 472);
    const auto *shi0_473 = buffer.data(shi0 + 473);
    const auto *shi0_475 = buffer.data(shi0 + 475);
    const auto *shi0_476 = buffer.data(shi0 + 476);
    const auto *shi0_479 = buffer.data(shi0 + 479);
    const auto *shi0_481 = buffer.data(shi0 + 481);
    const auto *shi0_482 = buffer.data(shi0 + 482);
    const auto *shi0_485 = buffer.data(shi0 + 485);
    const auto *shi0_486 = buffer.data(shi0 + 486);
    const auto *shi0_488 = buffer.data(shi0 + 488);
    const auto *shi0_490 = buffer.data(shi0 + 490);
    const auto *shi0_491 = buffer.data(shi0 + 491);
    const auto *shi0_493 = buffer.data(shi0 + 493);
    const auto *shi0_494 = buffer.data(shi0 + 494);
    const auto *shi0_496 = buffer.data(shi0 + 496);
    const auto *shi0_497 = buffer.data(shi0 + 497);
    const auto *shi0_499 = buffer.data(shi0 + 499);
    const auto *shi0_500 = buffer.data(shi0 + 500);
    const auto *shi0_501 = buffer.data(shi0 + 501);
    const auto *shi0_502 = buffer.data(shi0 + 502);
    const auto *shi0_503 = buffer.data(shi0 + 503);
    const auto *shi0_504 = buffer.data(shi0 + 504);
    const auto *shi0_507 = buffer.data(shi0 + 507);
    const auto *shi0_509 = buffer.data(shi0 + 509);

    const auto *shi1_438 = buffer.data(shi1 + 438);
    const auto *shi1_440 = buffer.data(shi1 + 440);
    const auto *shi1_441 = buffer.data(shi1 + 441);
    const auto *shi1_443 = buffer.data(shi1 + 443);
    const auto *shi1_444 = buffer.data(shi1 + 444);
    const auto *shi1_445 = buffer.data(shi1 + 445);
    const auto *shi1_446 = buffer.data(shi1 + 446);
    const auto *shi1_447 = buffer.data(shi1 + 447);
    const auto *shi1_453 = buffer.data(shi1 + 453);
    const auto *shi1_457 = buffer.data(shi1 + 457);
    const auto *shi1_460 = buffer.data(shi1 + 460);
    const auto *shi1_462 = buffer.data(shi1 + 462);
    const auto *shi1_465 = buffer.data(shi1 + 465);
    const auto *shi1_466 = buffer.data(shi1 + 466);
    const auto *shi1_468 = buffer.data(shi1 + 468);
    const auto *shi1_471 = buffer.data(shi1 + 471);
    const auto *shi1_472 = buffer.data(shi1 + 472);
    const auto *shi1_473 = buffer.data(shi1 + 473);
    const auto *shi1_475 = buffer.data(shi1 + 475);
    const auto *shi1_476 = buffer.data(shi1 + 476);
    const auto *shi1_479 = buffer.data(shi1 + 479);
    const auto *shi1_481 = buffer.data(shi1 + 481);
    const auto *shi1_482 = buffer.data(shi1 + 482);
    const auto *shi1_485 = buffer.data(shi1 + 485);
    const auto *shi1_486 = buffer.data(shi1 + 486);
    const auto *shi1_488 = buffer.data(shi1 + 488);
    const auto *shi1_490 = buffer.data(shi1 + 490);
    const auto *shi1_491 = buffer.data(shi1 + 491);
    const auto *shi1_493 = buffer.data(shi1 + 493);
    const auto *shi1_494 = buffer.data(shi1 + 494);
    const auto *shi1_496 = buffer.data(shi1 + 496);
    const auto *shi1_497 = buffer.data(shi1 + 497);
    const auto *shi1_499 = buffer.data(shi1 + 499);
    const auto *shi1_500 = buffer.data(shi1 + 500);
    const auto *shi1_501 = buffer.data(shi1 + 501);
    const auto *shi1_502 = buffer.data(shi1 + 502);
    const auto *shi1_503 = buffer.data(shi1 + 503);
    const auto *shi1_504 = buffer.data(shi1 + 504);
    const auto *shi1_507 = buffer.data(shi1 + 507);
    const auto *shi1_509 = buffer.data(shi1 + 509);

    const auto *shk_554 = buffer.data(shk + 554);
    const auto *shk_555 = buffer.data(shk + 555);
    const auto *shk_558 = buffer.data(shk + 558);
    const auto *shk_560 = buffer.data(shk + 560);
    const auto *shk_561 = buffer.data(shk + 561);
    const auto *shk_563 = buffer.data(shk + 563);
    const auto *shk_564 = buffer.data(shk + 564);
    const auto *shk_565 = buffer.data(shk + 565);
    const auto *shk_567 = buffer.data(shk + 567);
    const auto *shk_568 = buffer.data(shk + 568);
    const auto *shk_569 = buffer.data(shk + 569);
    const auto *shk_570 = buffer.data(shk + 570);
    const auto *shk_571 = buffer.data(shk + 571);
    const auto *shk_572 = buffer.data(shk + 572);
    const auto *shk_573 = buffer.data(shk + 573);
    const auto *shk_574 = buffer.data(shk + 574);
    const auto *shk_575 = buffer.data(shk + 575);
    const auto *shk_576 = buffer.data(shk + 576);
    const auto *shk_578 = buffer.data(shk + 578);
    const auto *shk_579 = buffer.data(shk + 579);
    const auto *shk_581 = buffer.data(shk + 581);
    const auto *shk_582 = buffer.data(shk + 582);
    const auto *shk_585 = buffer.data(shk + 585);
    const auto *shk_586 = buffer.data(shk + 586);
    const auto *shk_588 = buffer.data(shk + 588);
    const auto *shk_590 = buffer.data(shk + 590);
    const auto *shk_591 = buffer.data(shk + 591);
    const auto *shk_593 = buffer.data(shk + 593);
    const auto *shk_594 = buffer.data(shk + 594);
    const auto *shk_596 = buffer.data(shk + 596);
    const auto *shk_599 = buffer.data(shk + 599);
    const auto *shk_600 = buffer.data(shk + 600);
    const auto *shk_601 = buffer.data(shk + 601);
    const auto *shk_603 = buffer.data(shk + 603);
    const auto *shk_604 = buffer.data(shk + 604);
    const auto *shk_605 = buffer.data(shk + 605);
    const auto *shk_606 = buffer.data(shk + 606);
    const auto *shk_607 = buffer.data(shk + 607);
    const auto *shk_608 = buffer.data(shk + 608);
    const auto *shk_609 = buffer.data(shk + 609);
    const auto *shk_610 = buffer.data(shk + 610);
    const auto *shk_611 = buffer.data(shk + 611);
    const auto *shk_612 = buffer.data(shk + 612);
    const auto *shk_614 = buffer.data(shk + 614);
    const auto *shk_615 = buffer.data(shk + 615);
    const auto *shk_617 = buffer.data(shk + 617);
    const auto *shk_618 = buffer.data(shk + 618);
    const auto *shk_621 = buffer.data(shk + 621);
    const auto *shk_622 = buffer.data(shk + 622);
    const auto *shk_624 = buffer.data(shk + 624);
    const auto *shk_626 = buffer.data(shk + 626);
    const auto *shk_627 = buffer.data(shk + 627);
    const auto *shk_629 = buffer.data(shk + 629);
    const auto *shk_630 = buffer.data(shk + 630);
    const auto *shk_632 = buffer.data(shk + 632);
    const auto *shk_633 = buffer.data(shk + 633);
    const auto *shk_635 = buffer.data(shk + 635);
    const auto *shk_636 = buffer.data(shk + 636);
    const auto *shk_637 = buffer.data(shk + 637);
    const auto *shk_639 = buffer.data(shk + 639);
    const auto *shk_640 = buffer.data(shk + 640);
    const auto *shk_641 = buffer.data(shk + 641);
    const auto *shk_642 = buffer.data(shk + 642);
    const auto *shk_643 = buffer.data(shk + 643);
    const auto *shk_644 = buffer.data(shk + 644);
    const auto *shk_645 = buffer.data(shk + 645);
    const auto *shk_646 = buffer.data(shk + 646);
    const auto *shk_647 = buffer.data(shk + 647);
    const auto *shk_648 = buffer.data(shk + 648);
    const auto *shk_650 = buffer.data(shk + 650);
    const auto *shk_651 = buffer.data(shk + 651);
    const auto *shk_653 = buffer.data(shk + 653);

#pragma omp simd aligned(t_693, t_694, t_695, pc_x, pc_y, sgk_374, shi0_438, shi0_440, \
                         shi1_438, shi1_440, shk_554, shk_558, \
                         shk_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = f_10 * shi0_438[k]
                   - f_11 * shi1_438[k]
                   + f_3 * pc_x[k] * shk_558[k];

        t_694[k] = f_0 * sgk_374[k]
                   + f_3 * pc_y[k] * shk_554[k];

        t_695[k] = f_10 * shi0_440[k]
                   - f_11 * shi1_440[k]
                   + f_3 * pc_x[k] * shk_560[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, t_699, pc_x, pc_z, shi0_441, shi0_443, shi0_444, \
                         shi1_441, shi1_443, shi1_444, shk_555, shk_561, shk_563, \
                         shk_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = f_12 * shi0_441[k]
                   - f_13 * shi1_441[k]
                   + f_3 * pc_x[k] * shk_561[k];

        t_697[k] = f_3 * pc_z[k] * shk_555[k];

        t_698[k] = f_12 * shi0_443[k]
                   - f_13 * shi1_443[k]
                   + f_3 * pc_x[k] * shk_563[k];

        t_699[k] = f_12 * shi0_444[k]
                   - f_13 * shi1_444[k]
                   + f_3 * pc_x[k] * shk_564[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, pc_x, pc_y, sgk_380, shi0_445, shi0_447, \
                         shi1_445, shi1_447, shk_560, shk_565, shk_567, \
                         shk_568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_12 * shi0_445[k]
                   - f_13 * shi1_445[k]
                   + f_3 * pc_x[k] * shk_565[k];

        t_701[k] = f_0 * sgk_380[k]
                   + f_3 * pc_y[k] * shk_560[k];

        t_702[k] = f_12 * shi0_447[k]
                   - f_13 * shi1_447[k]
                   + f_3 * pc_x[k] * shk_567[k];

        t_703[k] = f_3 * pc_x[k] * shk_568[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, t_708, t_709, t_710, pc_x, shk_569, \
                         shk_570, shk_571, shk_572, shk_573, shk_574, \
                         shk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = f_3 * pc_x[k] * shk_569[k];

        t_705[k] = f_3 * pc_x[k] * shk_570[k];

        t_706[k] = f_3 * pc_x[k] * shk_571[k];

        t_707[k] = f_3 * pc_x[k] * shk_572[k];

        t_708[k] = f_3 * pc_x[k] * shk_573[k];

        t_709[k] = f_3 * pc_x[k] * shk_574[k];

        t_710[k] = f_3 * pc_x[k] * shk_575[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, pc_y, pc_z, sgk_388, sgk_390, shi0_441, \
                         shi0_443, shi1_441, shi1_443, shk_568, \
                         shk_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = f_0 * sgk_388[k]
                   + f_1 * shi0_441[k]
                   - f_2 * shi1_441[k]
                   + f_3 * pc_y[k] * shk_568[k];

        t_712[k] = f_3 * pc_z[k] * shk_568[k];

        t_713[k] = f_0 * sgk_390[k]
                   + f_4 * shi0_443[k]
                   - f_5 * shi1_443[k]
                   + f_3 * pc_y[k] * shk_570[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_y, sgk_391, sgk_392, sgk_393, shi0_444, \
                         shi0_445, shi0_446, shi1_444, shi1_445, shi1_446, shk_571, shk_572, \
                         shk_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_0 * sgk_391[k]
                   + f_6 * shi0_444[k]
                   - f_7 * shi1_444[k]
                   + f_3 * pc_y[k] * shk_571[k];

        t_715[k] = f_0 * sgk_392[k]
                   + f_8 * shi0_445[k]
                   - f_9 * shi1_445[k]
                   + f_3 * pc_y[k] * shk_572[k];

        t_716[k] = f_0 * sgk_393[k]
                   + f_10 * shi0_446[k]
                   - f_11 * shi1_446[k]
                   + f_3 * pc_y[k] * shk_573[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, t_720, pb_z, pc_y, pc_z, sgl0_450, sgk_394, \
                         sgk_395, sgl1_450, shi0_447, shi1_447, shk_574, \
                         shk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_0 * sgk_394[k]
                   + f_12 * shi0_447[k]
                   - f_13 * shi1_447[k]
                   + f_3 * pc_y[k] * shk_574[k];

        t_718[k] = f_0 * sgk_395[k]
                   + f_3 * pc_y[k] * shk_575[k];

        t_719[k] = f_1 * shi0_447[k]
                   - f_2 * shi1_447[k]
                   + f_3 * pc_z[k] * shk_575[k];

        t_720[k] = pb_z[k] * sgl0_450[k]
                   - f_14 * pc_z[k] * sgl1_450[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, t_724, pb_z, pc_y, pc_z, sgl0_453, sgk_360, \
                         sgk_396, sgk_398, sgl1_453, shk_576, shk_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = f_18 * sgk_396[k]
                   + f_3 * pc_y[k] * shk_576[k];

        t_722[k] = f_15 * sgk_360[k]
                   + f_3 * pc_z[k] * shk_576[k];

        t_723[k] = pb_z[k] * sgl0_453[k]
                   - f_14 * pc_z[k] * sgl1_453[k];

        t_724[k] = f_18 * sgk_398[k]
                   + f_3 * pc_y[k] * shk_578[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pb_z, pc_x, pc_y, pc_z, sgl0_456, \
                         sgk_363, sgk_401, sgl1_456, shi0_453, shi1_453, shk_579, \
                         shk_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_4 * shi0_453[k]
                   - f_5 * shi1_453[k]
                   + f_3 * pc_x[k] * shk_581[k];

        t_726[k] = pb_z[k] * sgl0_456[k]
                   - f_14 * pc_z[k] * sgl1_456[k];

        t_727[k] = f_15 * sgk_363[k]
                   + f_3 * pc_z[k] * shk_579[k];

        t_728[k] = f_18 * sgk_401[k]
                   + f_3 * pc_y[k] * shk_581[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pb_z, pc_x, pc_z, sgl0_460, sgk_366, sgl1_460, \
                         shi0_457, shi1_457, shk_582, shk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_6 * shi0_457[k]
                   - f_7 * shi1_457[k]
                   + f_3 * pc_x[k] * shk_585[k];

        t_730[k] = pb_z[k] * sgl0_460[k]
                   - f_14 * pc_z[k] * sgl1_460[k];

        t_731[k] = f_15 * sgk_366[k]
                   + f_3 * pc_z[k] * shk_582[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pc_x, pc_y, sgk_405, shi0_460, shi0_462, \
                         shi1_460, shi1_462, shk_585, shk_588, \
                         shk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_8 * shi0_460[k]
                   - f_9 * shi1_460[k]
                   + f_3 * pc_x[k] * shk_588[k];

        t_733[k] = f_18 * sgk_405[k]
                   + f_3 * pc_y[k] * shk_585[k];

        t_734[k] = f_8 * shi0_462[k]
                   - f_9 * shi1_462[k]
                   + f_3 * pc_x[k] * shk_590[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, pb_z, pc_x, pc_z, sgl0_465, sgk_370, sgl1_465, \
                         shi0_465, shi1_465, shk_586, shk_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = pb_z[k] * sgl0_465[k]
                   - f_14 * pc_z[k] * sgl1_465[k];

        t_736[k] = f_15 * sgk_370[k]
                   + f_3 * pc_z[k] * shk_586[k];

        t_737[k] = f_10 * shi0_465[k]
                   - f_11 * shi1_465[k]
                   + f_3 * pc_x[k] * shk_593[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pc_x, pc_y, sgk_410, shi0_466, shi0_468, \
                         shi1_466, shi1_468, shk_590, shk_594, \
                         shk_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = f_10 * shi0_466[k]
                   - f_11 * shi1_466[k]
                   + f_3 * pc_x[k] * shk_594[k];

        t_739[k] = f_18 * sgk_410[k]
                   + f_3 * pc_y[k] * shk_590[k];

        t_740[k] = f_10 * shi0_468[k]
                   - f_11 * shi1_468[k]
                   + f_3 * pc_x[k] * shk_596[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, pb_z, pc_x, pc_z, sgl0_471, sgk_375, sgl1_471, \
                         shi0_471, shi1_471, shk_591, shk_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = pb_z[k] * sgl0_471[k]
                   - f_14 * pc_z[k] * sgl1_471[k];

        t_742[k] = f_15 * sgk_375[k]
                   + f_3 * pc_z[k] * shk_591[k];

        t_743[k] = f_12 * shi0_471[k]
                   - f_13 * shi1_471[k]
                   + f_3 * pc_x[k] * shk_599[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, pc_x, pc_y, sgk_416, shi0_472, shi0_473, \
                         shi1_472, shi1_473, shk_596, shk_600, \
                         shk_601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_12 * shi0_472[k]
                   - f_13 * shi1_472[k]
                   + f_3 * pc_x[k] * shk_600[k];

        t_745[k] = f_12 * shi0_473[k]
                   - f_13 * shi1_473[k]
                   + f_3 * pc_x[k] * shk_601[k];

        t_746[k] = f_18 * sgk_416[k]
                   + f_3 * pc_y[k] * shk_596[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, t_750, t_751, t_752, pc_x, shi0_475, shi1_475, \
                         shk_603, shk_604, shk_605, shk_606, shk_607, \
                         shk_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_12 * shi0_475[k]
                   - f_13 * shi1_475[k]
                   + f_3 * pc_x[k] * shk_603[k];

        t_748[k] = f_3 * pc_x[k] * shk_604[k];

        t_749[k] = f_3 * pc_x[k] * shk_605[k];

        t_750[k] = f_3 * pc_x[k] * shk_606[k];

        t_751[k] = f_3 * pc_x[k] * shk_607[k];

        t_752[k] = f_3 * pc_x[k] * shk_608[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, t_756, t_757, pb_z, pc_x, pc_z, sgl0_486, \
                         sgk_388, sgl1_486, shk_604, shk_609, shk_610, \
                         shk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_3 * pc_x[k] * shk_609[k];

        t_754[k] = f_3 * pc_x[k] * shk_610[k];

        t_755[k] = f_3 * pc_x[k] * shk_611[k];

        t_756[k] = pb_z[k] * sgl0_486[k]
                   - f_14 * pc_z[k] * sgl1_486[k];

        t_757[k] = f_15 * sgk_388[k]
                   + f_3 * pc_z[k] * shk_604[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, pb_z, pc_z, sgl0_488, sgl0_489, sgl0_490, \
                         sgk_389, sgk_390, sgk_391, sgl1_488, sgl1_489, \
                         sgl1_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = pb_z[k] * sgl0_488[k]
                   + f_16 * sgk_389[k]
                   - f_14 * pc_z[k] * sgl1_488[k];

        t_759[k] = pb_z[k] * sgl0_489[k]
                   + f_17 * sgk_390[k]
                   - f_14 * pc_z[k] * sgl1_489[k];

        t_760[k] = pb_z[k] * sgl0_490[k]
                   + f_18 * sgk_391[k]
                   - f_14 * pc_z[k] * sgl1_490[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pb_z, pc_y, pc_z, sgl0_491, sgl0_492, sgk_392, \
                         sgk_393, sgk_431, sgl1_491, sgl1_492, \
                         shk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = pb_z[k] * sgl0_491[k]
                   + f_0 * sgk_392[k]
                   - f_14 * pc_z[k] * sgl1_491[k];

        t_762[k] = pb_z[k] * sgl0_492[k]
                   + f_19 * sgk_393[k]
                   - f_14 * pc_z[k] * sgl1_492[k];

        t_763[k] = f_18 * sgk_431[k]
                   + f_3 * pc_y[k] * shk_611[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pc_x, pc_y, pc_z, sgk_395, sgk_396, \
                         sgk_432, shi0_475, shi0_476, shi1_475, shi1_476, shk_611, \
                         shk_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_15 * sgk_395[k]
                   + f_1 * shi0_475[k]
                   - f_2 * shi1_475[k]
                   + f_3 * pc_z[k] * shk_611[k];

        t_765[k] = f_1 * shi0_476[k]
                   - f_2 * shi1_476[k]
                   + f_3 * pc_x[k] * shk_612[k];

        t_766[k] = f_17 * sgk_432[k]
                   + f_3 * pc_y[k] * shk_612[k];

        t_767[k] = f_16 * sgk_396[k]
                   + f_3 * pc_z[k] * shk_612[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, pc_x, pc_y, sgk_434, shi0_479, shi0_481, \
                         shi1_479, shi1_481, shk_614, shk_615, \
                         shk_617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_4 * shi0_479[k]
                   - f_5 * shi1_479[k]
                   + f_3 * pc_x[k] * shk_615[k];

        t_769[k] = f_17 * sgk_434[k]
                   + f_3 * pc_y[k] * shk_614[k];

        t_770[k] = f_4 * shi0_481[k]
                   - f_5 * shi1_481[k]
                   + f_3 * pc_x[k] * shk_617[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, pc_x, pc_y, pc_z, sgk_399, sgk_437, shi0_482, \
                         shi1_482, shk_615, shk_617, shk_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = f_6 * shi0_482[k]
                   - f_7 * shi1_482[k]
                   + f_3 * pc_x[k] * shk_618[k];

        t_772[k] = f_16 * sgk_399[k]
                   + f_3 * pc_z[k] * shk_615[k];

        t_773[k] = f_17 * sgk_437[k]
                   + f_3 * pc_y[k] * shk_617[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, pc_x, pc_z, sgk_402, shi0_485, shi0_486, \
                         shi1_485, shi1_486, shk_618, shk_621, \
                         shk_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_6 * shi0_485[k]
                   - f_7 * shi1_485[k]
                   + f_3 * pc_x[k] * shk_621[k];

        t_775[k] = f_8 * shi0_486[k]
                   - f_9 * shi1_486[k]
                   + f_3 * pc_x[k] * shk_622[k];

        t_776[k] = f_16 * sgk_402[k]
                   + f_3 * pc_z[k] * shk_618[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pc_x, pc_y, sgk_441, shi0_488, shi0_490, \
                         shi1_488, shi1_490, shk_621, shk_624, \
                         shk_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_8 * shi0_488[k]
                   - f_9 * shi1_488[k]
                   + f_3 * pc_x[k] * shk_624[k];

        t_778[k] = f_17 * sgk_441[k]
                   + f_3 * pc_y[k] * shk_621[k];

        t_779[k] = f_8 * shi0_490[k]
                   - f_9 * shi1_490[k]
                   + f_3 * pc_x[k] * shk_626[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, pc_x, pc_z, sgk_406, shi0_491, shi0_493, \
                         shi1_491, shi1_493, shk_622, shk_627, \
                         shk_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_10 * shi0_491[k]
                   - f_11 * shi1_491[k]
                   + f_3 * pc_x[k] * shk_627[k];

        t_781[k] = f_16 * sgk_406[k]
                   + f_3 * pc_z[k] * shk_622[k];

        t_782[k] = f_10 * shi0_493[k]
                   - f_11 * shi1_493[k]
                   + f_3 * pc_x[k] * shk_629[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, pc_x, pc_y, sgk_446, shi0_494, shi0_496, \
                         shi1_494, shi1_496, shk_626, shk_630, \
                         shk_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = f_10 * shi0_494[k]
                   - f_11 * shi1_494[k]
                   + f_3 * pc_x[k] * shk_630[k];

        t_784[k] = f_17 * sgk_446[k]
                   + f_3 * pc_y[k] * shk_626[k];

        t_785[k] = f_10 * shi0_496[k]
                   - f_11 * shi1_496[k]
                   + f_3 * pc_x[k] * shk_632[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, pc_x, pc_z, sgk_411, shi0_497, shi0_499, \
                         shi1_497, shi1_499, shk_627, shk_633, \
                         shk_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_12 * shi0_497[k]
                   - f_13 * shi1_497[k]
                   + f_3 * pc_x[k] * shk_633[k];

        t_787[k] = f_16 * sgk_411[k]
                   + f_3 * pc_z[k] * shk_627[k];

        t_788[k] = f_12 * shi0_499[k]
                   - f_13 * shi1_499[k]
                   + f_3 * pc_x[k] * shk_635[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, pc_x, pc_y, sgk_452, shi0_500, shi0_501, \
                         shi1_500, shi1_501, shk_632, shk_636, \
                         shk_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = f_12 * shi0_500[k]
                   - f_13 * shi1_500[k]
                   + f_3 * pc_x[k] * shk_636[k];

        t_790[k] = f_12 * shi0_501[k]
                   - f_13 * shi1_501[k]
                   + f_3 * pc_x[k] * shk_637[k];

        t_791[k] = f_17 * sgk_452[k]
                   + f_3 * pc_y[k] * shk_632[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, t_796, t_797, pc_x, shi0_503, shi1_503, \
                         shk_639, shk_640, shk_641, shk_642, shk_643, \
                         shk_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = f_12 * shi0_503[k]
                   - f_13 * shi1_503[k]
                   + f_3 * pc_x[k] * shk_639[k];

        t_793[k] = f_3 * pc_x[k] * shk_640[k];

        t_794[k] = f_3 * pc_x[k] * shk_641[k];

        t_795[k] = f_3 * pc_x[k] * shk_642[k];

        t_796[k] = f_3 * pc_x[k] * shk_643[k];

        t_797[k] = f_3 * pc_x[k] * shk_644[k];
    }

#pragma omp simd aligned(t_798, t_799, t_800, t_801, t_802, pc_x, pc_y, pc_z, sgk_424, \
                         sgk_460, shi0_497, shi1_497, shk_640, shk_645, shk_646, \
                         shk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_798[k] = f_3 * pc_x[k] * shk_645[k];

        t_799[k] = f_3 * pc_x[k] * shk_646[k];

        t_800[k] = f_3 * pc_x[k] * shk_647[k];

        t_801[k] = f_17 * sgk_460[k]
                   + f_1 * shi0_497[k]
                   - f_2 * shi1_497[k]
                   + f_3 * pc_y[k] * shk_640[k];

        t_802[k] = f_16 * sgk_424[k]
                   + f_3 * pc_z[k] * shk_640[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_y, sgk_462, sgk_463, sgk_464, shi0_499, \
                         shi0_500, shi0_501, shi1_499, shi1_500, shi1_501, shk_642, shk_643, \
                         shk_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_17 * sgk_462[k]
                   + f_4 * shi0_499[k]
                   - f_5 * shi1_499[k]
                   + f_3 * pc_y[k] * shk_642[k];

        t_804[k] = f_17 * sgk_463[k]
                   + f_6 * shi0_500[k]
                   - f_7 * shi1_500[k]
                   + f_3 * pc_y[k] * shk_643[k];

        t_805[k] = f_17 * sgk_464[k]
                   + f_8 * shi0_501[k]
                   - f_9 * shi1_501[k]
                   + f_3 * pc_y[k] * shk_644[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pc_y, sgk_465, sgk_466, sgk_467, shi0_502, \
                         shi0_503, shi1_502, shi1_503, shk_645, shk_646, \
                         shk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_17 * sgk_465[k]
                   + f_10 * shi0_502[k]
                   - f_11 * shi1_502[k]
                   + f_3 * pc_y[k] * shk_645[k];

        t_807[k] = f_17 * sgk_466[k]
                   + f_12 * shi0_503[k]
                   - f_13 * shi1_503[k]
                   + f_3 * pc_y[k] * shk_646[k];

        t_808[k] = f_17 * sgk_467[k]
                   + f_3 * pc_y[k] * shk_647[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, t_812, pc_x, pc_y, pc_z, sgk_431, sgk_432, \
                         sgk_468, shi0_503, shi0_504, shi1_503, shi1_504, shk_647, \
                         shk_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_16 * sgk_431[k]
                   + f_1 * shi0_503[k]
                   - f_2 * shi1_503[k]
                   + f_3 * pc_z[k] * shk_647[k];

        t_810[k] = f_1 * shi0_504[k]
                   - f_2 * shi1_504[k]
                   + f_3 * pc_x[k] * shk_648[k];

        t_811[k] = f_16 * sgk_468[k]
                   + f_3 * pc_y[k] * shk_648[k];

        t_812[k] = f_17 * sgk_432[k]
                   + f_3 * pc_z[k] * shk_648[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, pc_x, pc_y, sgk_470, shi0_507, shi0_509, \
                         shi1_507, shi1_509, shk_650, shk_651, \
                         shk_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = f_4 * shi0_507[k]
                   - f_5 * shi1_507[k]
                   + f_3 * pc_x[k] * shk_651[k];

        t_814[k] = f_16 * sgk_470[k]
                   + f_3 * pc_y[k] * shk_650[k];

        t_815[k] = f_4 * shi0_509[k]
                   - f_5 * shi1_509[k]
                   + f_3 * pc_x[k] * shk_653[k];
    }
}

static auto
compute_prim_shl_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgl0,
                                                          const size_t sgk, const size_t sgl1,
                                                          const size_t shi0, const size_t shi1,
                                                          const size_t shk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_18 = 2.0 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 4.0 / q;

    auto *t_816 = buffer.data(target + 816);
    auto *t_817 = buffer.data(target + 817);
    auto *t_818 = buffer.data(target + 818);
    auto *t_819 = buffer.data(target + 819);
    auto *t_820 = buffer.data(target + 820);
    auto *t_821 = buffer.data(target + 821);
    auto *t_822 = buffer.data(target + 822);
    auto *t_823 = buffer.data(target + 823);
    auto *t_824 = buffer.data(target + 824);
    auto *t_825 = buffer.data(target + 825);
    auto *t_826 = buffer.data(target + 826);
    auto *t_827 = buffer.data(target + 827);
    auto *t_828 = buffer.data(target + 828);
    auto *t_829 = buffer.data(target + 829);
    auto *t_830 = buffer.data(target + 830);
    auto *t_831 = buffer.data(target + 831);
    auto *t_832 = buffer.data(target + 832);
    auto *t_833 = buffer.data(target + 833);
    auto *t_834 = buffer.data(target + 834);
    auto *t_835 = buffer.data(target + 835);
    auto *t_836 = buffer.data(target + 836);
    auto *t_837 = buffer.data(target + 837);
    auto *t_838 = buffer.data(target + 838);
    auto *t_839 = buffer.data(target + 839);
    auto *t_840 = buffer.data(target + 840);
    auto *t_841 = buffer.data(target + 841);
    auto *t_842 = buffer.data(target + 842);
    auto *t_843 = buffer.data(target + 843);
    auto *t_844 = buffer.data(target + 844);
    auto *t_845 = buffer.data(target + 845);
    auto *t_846 = buffer.data(target + 846);
    auto *t_847 = buffer.data(target + 847);
    auto *t_848 = buffer.data(target + 848);
    auto *t_849 = buffer.data(target + 849);
    auto *t_850 = buffer.data(target + 850);
    auto *t_851 = buffer.data(target + 851);
    auto *t_852 = buffer.data(target + 852);
    auto *t_853 = buffer.data(target + 853);
    auto *t_854 = buffer.data(target + 854);
    auto *t_855 = buffer.data(target + 855);
    auto *t_856 = buffer.data(target + 856);
    auto *t_857 = buffer.data(target + 857);
    auto *t_858 = buffer.data(target + 858);
    auto *t_859 = buffer.data(target + 859);
    auto *t_860 = buffer.data(target + 860);
    auto *t_861 = buffer.data(target + 861);
    auto *t_862 = buffer.data(target + 862);
    auto *t_863 = buffer.data(target + 863);
    auto *t_864 = buffer.data(target + 864);
    auto *t_865 = buffer.data(target + 865);
    auto *t_866 = buffer.data(target + 866);
    auto *t_867 = buffer.data(target + 867);
    auto *t_868 = buffer.data(target + 868);
    auto *t_869 = buffer.data(target + 869);
    auto *t_870 = buffer.data(target + 870);
    auto *t_871 = buffer.data(target + 871);
    auto *t_872 = buffer.data(target + 872);
    auto *t_873 = buffer.data(target + 873);
    auto *t_874 = buffer.data(target + 874);
    auto *t_875 = buffer.data(target + 875);
    auto *t_876 = buffer.data(target + 876);
    auto *t_877 = buffer.data(target + 877);
    auto *t_878 = buffer.data(target + 878);
    auto *t_879 = buffer.data(target + 879);
    auto *t_880 = buffer.data(target + 880);
    auto *t_881 = buffer.data(target + 881);
    auto *t_882 = buffer.data(target + 882);
    auto *t_883 = buffer.data(target + 883);
    auto *t_884 = buffer.data(target + 884);
    auto *t_885 = buffer.data(target + 885);
    auto *t_886 = buffer.data(target + 886);
    auto *t_887 = buffer.data(target + 887);
    auto *t_888 = buffer.data(target + 888);
    auto *t_889 = buffer.data(target + 889);
    auto *t_890 = buffer.data(target + 890);
    auto *t_891 = buffer.data(target + 891);
    auto *t_892 = buffer.data(target + 892);
    auto *t_893 = buffer.data(target + 893);
    auto *t_894 = buffer.data(target + 894);
    auto *t_895 = buffer.data(target + 895);
    auto *t_896 = buffer.data(target + 896);
    auto *t_897 = buffer.data(target + 897);
    auto *t_898 = buffer.data(target + 898);
    auto *t_899 = buffer.data(target + 899);
    auto *t_900 = buffer.data(target + 900);
    auto *t_901 = buffer.data(target + 901);
    auto *t_902 = buffer.data(target + 902);
    auto *t_903 = buffer.data(target + 903);
    auto *t_904 = buffer.data(target + 904);
    auto *t_905 = buffer.data(target + 905);
    auto *t_906 = buffer.data(target + 906);
    auto *t_907 = buffer.data(target + 907);
    auto *t_908 = buffer.data(target + 908);
    auto *t_909 = buffer.data(target + 909);
    auto *t_910 = buffer.data(target + 910);
    auto *t_911 = buffer.data(target + 911);
    auto *t_912 = buffer.data(target + 912);
    auto *t_913 = buffer.data(target + 913);
    auto *t_914 = buffer.data(target + 914);
    auto *t_915 = buffer.data(target + 915);
    auto *t_916 = buffer.data(target + 916);
    auto *t_917 = buffer.data(target + 917);
    auto *t_918 = buffer.data(target + 918);
    auto *t_919 = buffer.data(target + 919);
    auto *t_920 = buffer.data(target + 920);
    auto *t_921 = buffer.data(target + 921);
    auto *t_922 = buffer.data(target + 922);
    auto *t_923 = buffer.data(target + 923);
    auto *t_924 = buffer.data(target + 924);
    auto *t_925 = buffer.data(target + 925);
    auto *t_926 = buffer.data(target + 926);
    auto *t_927 = buffer.data(target + 927);
    auto *t_928 = buffer.data(target + 928);
    auto *t_929 = buffer.data(target + 929);
    auto *t_930 = buffer.data(target + 930);
    auto *t_931 = buffer.data(target + 931);
    auto *t_932 = buffer.data(target + 932);
    auto *t_933 = buffer.data(target + 933);
    auto *t_934 = buffer.data(target + 934);
    auto *t_935 = buffer.data(target + 935);
    auto *t_936 = buffer.data(target + 936);
    auto *t_937 = buffer.data(target + 937);
    auto *t_938 = buffer.data(target + 938);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgl0_630 = buffer.data(sgl0 + 630);
    const auto *sgl0_635 = buffer.data(sgl0 + 635);
    const auto *sgl0_639 = buffer.data(sgl0 + 639);
    const auto *sgl0_644 = buffer.data(sgl0 + 644);
    const auto *sgl0_650 = buffer.data(sgl0 + 650);
    const auto *sgl0_657 = buffer.data(sgl0 + 657);
    const auto *sgl0_666 = buffer.data(sgl0 + 666);
    const auto *sgl0_668 = buffer.data(sgl0 + 668);
    const auto *sgl0_669 = buffer.data(sgl0 + 669);
    const auto *sgl0_670 = buffer.data(sgl0 + 670);
    const auto *sgl0_671 = buffer.data(sgl0 + 671);
    const auto *sgl0_672 = buffer.data(sgl0 + 672);
    const auto *sgl0_674 = buffer.data(sgl0 + 674);

    const auto *sgk_435 = buffer.data(sgk + 435);
    const auto *sgk_438 = buffer.data(sgk + 438);
    const auto *sgk_442 = buffer.data(sgk + 442);
    const auto *sgk_447 = buffer.data(sgk + 447);
    const auto *sgk_460 = buffer.data(sgk + 460);
    const auto *sgk_467 = buffer.data(sgk + 467);
    const auto *sgk_468 = buffer.data(sgk + 468);
    const auto *sgk_471 = buffer.data(sgk + 471);
    const auto *sgk_473 = buffer.data(sgk + 473);
    const auto *sgk_474 = buffer.data(sgk + 474);
    const auto *sgk_477 = buffer.data(sgk + 477);
    const auto *sgk_478 = buffer.data(sgk + 478);
    const auto *sgk_482 = buffer.data(sgk + 482);
    const auto *sgk_483 = buffer.data(sgk + 483);
    const auto *sgk_488 = buffer.data(sgk + 488);
    const auto *sgk_496 = buffer.data(sgk + 496);
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
    const auto *sgk_518 = buffer.data(sgk + 518);
    const auto *sgk_519 = buffer.data(sgk + 519);
    const auto *sgk_524 = buffer.data(sgk + 524);
    const auto *sgk_532 = buffer.data(sgk + 532);
    const auto *sgk_534 = buffer.data(sgk + 534);
    const auto *sgk_535 = buffer.data(sgk + 535);
    const auto *sgk_536 = buffer.data(sgk + 536);
    const auto *sgk_537 = buffer.data(sgk + 537);
    const auto *sgk_538 = buffer.data(sgk + 538);
    const auto *sgk_539 = buffer.data(sgk + 539);

    const auto *sgl1_630 = buffer.data(sgl1 + 630);
    const auto *sgl1_635 = buffer.data(sgl1 + 635);
    const auto *sgl1_639 = buffer.data(sgl1 + 639);
    const auto *sgl1_644 = buffer.data(sgl1 + 644);
    const auto *sgl1_650 = buffer.data(sgl1 + 650);
    const auto *sgl1_657 = buffer.data(sgl1 + 657);
    const auto *sgl1_666 = buffer.data(sgl1 + 666);
    const auto *sgl1_668 = buffer.data(sgl1 + 668);
    const auto *sgl1_669 = buffer.data(sgl1 + 669);
    const auto *sgl1_670 = buffer.data(sgl1 + 670);
    const auto *sgl1_671 = buffer.data(sgl1 + 671);
    const auto *sgl1_672 = buffer.data(sgl1 + 672);
    const auto *sgl1_674 = buffer.data(sgl1 + 674);

    const auto *shi0_510 = buffer.data(shi0 + 510);
    const auto *shi0_513 = buffer.data(shi0 + 513);
    const auto *shi0_514 = buffer.data(shi0 + 514);
    const auto *shi0_516 = buffer.data(shi0 + 516);
    const auto *shi0_518 = buffer.data(shi0 + 518);
    const auto *shi0_519 = buffer.data(shi0 + 519);
    const auto *shi0_521 = buffer.data(shi0 + 521);
    const auto *shi0_522 = buffer.data(shi0 + 522);
    const auto *shi0_524 = buffer.data(shi0 + 524);
    const auto *shi0_525 = buffer.data(shi0 + 525);
    const auto *shi0_527 = buffer.data(shi0 + 527);
    const auto *shi0_528 = buffer.data(shi0 + 528);
    const auto *shi0_529 = buffer.data(shi0 + 529);
    const auto *shi0_530 = buffer.data(shi0 + 530);
    const auto *shi0_531 = buffer.data(shi0 + 531);
    const auto *shi0_535 = buffer.data(shi0 + 535);
    const auto *shi0_538 = buffer.data(shi0 + 538);
    const auto *shi0_542 = buffer.data(shi0 + 542);
    const auto *shi0_544 = buffer.data(shi0 + 544);
    const auto *shi0_547 = buffer.data(shi0 + 547);
    const auto *shi0_549 = buffer.data(shi0 + 549);
    const auto *shi0_550 = buffer.data(shi0 + 550);
    const auto *shi0_553 = buffer.data(shi0 + 553);
    const auto *shi0_555 = buffer.data(shi0 + 555);
    const auto *shi0_556 = buffer.data(shi0 + 556);
    const auto *shi0_557 = buffer.data(shi0 + 557);
    const auto *shi0_560 = buffer.data(shi0 + 560);
    const auto *shi0_563 = buffer.data(shi0 + 563);
    const auto *shi0_565 = buffer.data(shi0 + 565);
    const auto *shi0_566 = buffer.data(shi0 + 566);
    const auto *shi0_569 = buffer.data(shi0 + 569);
    const auto *shi0_570 = buffer.data(shi0 + 570);
    const auto *shi0_572 = buffer.data(shi0 + 572);
    const auto *shi0_574 = buffer.data(shi0 + 574);
    const auto *shi0_575 = buffer.data(shi0 + 575);
    const auto *shi0_577 = buffer.data(shi0 + 577);
    const auto *shi0_578 = buffer.data(shi0 + 578);
    const auto *shi0_580 = buffer.data(shi0 + 580);
    const auto *shi0_581 = buffer.data(shi0 + 581);
    const auto *shi0_583 = buffer.data(shi0 + 583);
    const auto *shi0_584 = buffer.data(shi0 + 584);
    const auto *shi0_585 = buffer.data(shi0 + 585);
    const auto *shi0_587 = buffer.data(shi0 + 587);

    const auto *shi1_510 = buffer.data(shi1 + 510);
    const auto *shi1_513 = buffer.data(shi1 + 513);
    const auto *shi1_514 = buffer.data(shi1 + 514);
    const auto *shi1_516 = buffer.data(shi1 + 516);
    const auto *shi1_518 = buffer.data(shi1 + 518);
    const auto *shi1_519 = buffer.data(shi1 + 519);
    const auto *shi1_521 = buffer.data(shi1 + 521);
    const auto *shi1_522 = buffer.data(shi1 + 522);
    const auto *shi1_524 = buffer.data(shi1 + 524);
    const auto *shi1_525 = buffer.data(shi1 + 525);
    const auto *shi1_527 = buffer.data(shi1 + 527);
    const auto *shi1_528 = buffer.data(shi1 + 528);
    const auto *shi1_529 = buffer.data(shi1 + 529);
    const auto *shi1_530 = buffer.data(shi1 + 530);
    const auto *shi1_531 = buffer.data(shi1 + 531);
    const auto *shi1_535 = buffer.data(shi1 + 535);
    const auto *shi1_538 = buffer.data(shi1 + 538);
    const auto *shi1_542 = buffer.data(shi1 + 542);
    const auto *shi1_544 = buffer.data(shi1 + 544);
    const auto *shi1_547 = buffer.data(shi1 + 547);
    const auto *shi1_549 = buffer.data(shi1 + 549);
    const auto *shi1_550 = buffer.data(shi1 + 550);
    const auto *shi1_553 = buffer.data(shi1 + 553);
    const auto *shi1_555 = buffer.data(shi1 + 555);
    const auto *shi1_556 = buffer.data(shi1 + 556);
    const auto *shi1_557 = buffer.data(shi1 + 557);
    const auto *shi1_560 = buffer.data(shi1 + 560);
    const auto *shi1_563 = buffer.data(shi1 + 563);
    const auto *shi1_565 = buffer.data(shi1 + 565);
    const auto *shi1_566 = buffer.data(shi1 + 566);
    const auto *shi1_569 = buffer.data(shi1 + 569);
    const auto *shi1_570 = buffer.data(shi1 + 570);
    const auto *shi1_572 = buffer.data(shi1 + 572);
    const auto *shi1_574 = buffer.data(shi1 + 574);
    const auto *shi1_575 = buffer.data(shi1 + 575);
    const auto *shi1_577 = buffer.data(shi1 + 577);
    const auto *shi1_578 = buffer.data(shi1 + 578);
    const auto *shi1_580 = buffer.data(shi1 + 580);
    const auto *shi1_581 = buffer.data(shi1 + 581);
    const auto *shi1_583 = buffer.data(shi1 + 583);
    const auto *shi1_584 = buffer.data(shi1 + 584);
    const auto *shi1_585 = buffer.data(shi1 + 585);
    const auto *shi1_587 = buffer.data(shi1 + 587);

    const auto *shk_651 = buffer.data(shk + 651);
    const auto *shk_653 = buffer.data(shk + 653);
    const auto *shk_654 = buffer.data(shk + 654);
    const auto *shk_657 = buffer.data(shk + 657);
    const auto *shk_658 = buffer.data(shk + 658);
    const auto *shk_660 = buffer.data(shk + 660);
    const auto *shk_662 = buffer.data(shk + 662);
    const auto *shk_663 = buffer.data(shk + 663);
    const auto *shk_665 = buffer.data(shk + 665);
    const auto *shk_666 = buffer.data(shk + 666);
    const auto *shk_668 = buffer.data(shk + 668);
    const auto *shk_669 = buffer.data(shk + 669);
    const auto *shk_671 = buffer.data(shk + 671);
    const auto *shk_672 = buffer.data(shk + 672);
    const auto *shk_673 = buffer.data(shk + 673);
    const auto *shk_675 = buffer.data(shk + 675);
    const auto *shk_676 = buffer.data(shk + 676);
    const auto *shk_677 = buffer.data(shk + 677);
    const auto *shk_678 = buffer.data(shk + 678);
    const auto *shk_679 = buffer.data(shk + 679);
    const auto *shk_680 = buffer.data(shk + 680);
    const auto *shk_681 = buffer.data(shk + 681);
    const auto *shk_682 = buffer.data(shk + 682);
    const auto *shk_683 = buffer.data(shk + 683);
    const auto *shk_684 = buffer.data(shk + 684);
    const auto *shk_686 = buffer.data(shk + 686);
    const auto *shk_687 = buffer.data(shk + 687);
    const auto *shk_689 = buffer.data(shk + 689);
    const auto *shk_690 = buffer.data(shk + 690);
    const auto *shk_693 = buffer.data(shk + 693);
    const auto *shk_694 = buffer.data(shk + 694);
    const auto *shk_696 = buffer.data(shk + 696);
    const auto *shk_698 = buffer.data(shk + 698);
    const auto *shk_699 = buffer.data(shk + 699);
    const auto *shk_701 = buffer.data(shk + 701);
    const auto *shk_702 = buffer.data(shk + 702);
    const auto *shk_704 = buffer.data(shk + 704);
    const auto *shk_705 = buffer.data(shk + 705);
    const auto *shk_707 = buffer.data(shk + 707);
    const auto *shk_708 = buffer.data(shk + 708);
    const auto *shk_709 = buffer.data(shk + 709);
    const auto *shk_712 = buffer.data(shk + 712);
    const auto *shk_713 = buffer.data(shk + 713);
    const auto *shk_714 = buffer.data(shk + 714);
    const auto *shk_715 = buffer.data(shk + 715);
    const auto *shk_716 = buffer.data(shk + 716);
    const auto *shk_717 = buffer.data(shk + 717);
    const auto *shk_718 = buffer.data(shk + 718);
    const auto *shk_719 = buffer.data(shk + 719);
    const auto *shk_720 = buffer.data(shk + 720);
    const auto *shk_722 = buffer.data(shk + 722);
    const auto *shk_723 = buffer.data(shk + 723);
    const auto *shk_725 = buffer.data(shk + 725);
    const auto *shk_726 = buffer.data(shk + 726);
    const auto *shk_729 = buffer.data(shk + 729);
    const auto *shk_730 = buffer.data(shk + 730);
    const auto *shk_732 = buffer.data(shk + 732);
    const auto *shk_734 = buffer.data(shk + 734);
    const auto *shk_735 = buffer.data(shk + 735);
    const auto *shk_737 = buffer.data(shk + 737);
    const auto *shk_738 = buffer.data(shk + 738);
    const auto *shk_740 = buffer.data(shk + 740);
    const auto *shk_741 = buffer.data(shk + 741);
    const auto *shk_743 = buffer.data(shk + 743);
    const auto *shk_744 = buffer.data(shk + 744);
    const auto *shk_745 = buffer.data(shk + 745);
    const auto *shk_747 = buffer.data(shk + 747);
    const auto *shk_748 = buffer.data(shk + 748);
    const auto *shk_749 = buffer.data(shk + 749);
    const auto *shk_750 = buffer.data(shk + 750);
    const auto *shk_751 = buffer.data(shk + 751);
    const auto *shk_752 = buffer.data(shk + 752);
    const auto *shk_753 = buffer.data(shk + 753);
    const auto *shk_754 = buffer.data(shk + 754);
    const auto *shk_755 = buffer.data(shk + 755);

#pragma omp simd aligned(t_816, t_817, t_818, pc_x, pc_y, pc_z, sgk_435, sgk_473, shi0_510, \
                         shi1_510, shk_651, shk_653, shk_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_6 * shi0_510[k]
                   - f_7 * shi1_510[k]
                   + f_3 * pc_x[k] * shk_654[k];

        t_817[k] = f_17 * sgk_435[k]
                   + f_3 * pc_z[k] * shk_651[k];

        t_818[k] = f_16 * sgk_473[k]
                   + f_3 * pc_y[k] * shk_653[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pc_x, pc_z, sgk_438, shi0_513, shi0_514, \
                         shi1_513, shi1_514, shk_654, shk_657, \
                         shk_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_6 * shi0_513[k]
                   - f_7 * shi1_513[k]
                   + f_3 * pc_x[k] * shk_657[k];

        t_820[k] = f_8 * shi0_514[k]
                   - f_9 * shi1_514[k]
                   + f_3 * pc_x[k] * shk_658[k];

        t_821[k] = f_17 * sgk_438[k]
                   + f_3 * pc_z[k] * shk_654[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, pc_x, pc_y, sgk_477, shi0_516, shi0_518, \
                         shi1_516, shi1_518, shk_657, shk_660, \
                         shk_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = f_8 * shi0_516[k]
                   - f_9 * shi1_516[k]
                   + f_3 * pc_x[k] * shk_660[k];

        t_823[k] = f_16 * sgk_477[k]
                   + f_3 * pc_y[k] * shk_657[k];

        t_824[k] = f_8 * shi0_518[k]
                   - f_9 * shi1_518[k]
                   + f_3 * pc_x[k] * shk_662[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, pc_x, pc_z, sgk_442, shi0_519, shi0_521, \
                         shi1_519, shi1_521, shk_658, shk_663, \
                         shk_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = f_10 * shi0_519[k]
                   - f_11 * shi1_519[k]
                   + f_3 * pc_x[k] * shk_663[k];

        t_826[k] = f_17 * sgk_442[k]
                   + f_3 * pc_z[k] * shk_658[k];

        t_827[k] = f_10 * shi0_521[k]
                   - f_11 * shi1_521[k]
                   + f_3 * pc_x[k] * shk_665[k];
    }

#pragma omp simd aligned(t_828, t_829, t_830, pc_x, pc_y, sgk_482, shi0_522, shi0_524, \
                         shi1_522, shi1_524, shk_662, shk_666, \
                         shk_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_828[k] = f_10 * shi0_522[k]
                   - f_11 * shi1_522[k]
                   + f_3 * pc_x[k] * shk_666[k];

        t_829[k] = f_16 * sgk_482[k]
                   + f_3 * pc_y[k] * shk_662[k];

        t_830[k] = f_10 * shi0_524[k]
                   - f_11 * shi1_524[k]
                   + f_3 * pc_x[k] * shk_668[k];
    }

#pragma omp simd aligned(t_831, t_832, t_833, pc_x, pc_z, sgk_447, shi0_525, shi0_527, \
                         shi1_525, shi1_527, shk_663, shk_669, \
                         shk_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_831[k] = f_12 * shi0_525[k]
                   - f_13 * shi1_525[k]
                   + f_3 * pc_x[k] * shk_669[k];

        t_832[k] = f_17 * sgk_447[k]
                   + f_3 * pc_z[k] * shk_663[k];

        t_833[k] = f_12 * shi0_527[k]
                   - f_13 * shi1_527[k]
                   + f_3 * pc_x[k] * shk_671[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, pc_x, pc_y, sgk_488, shi0_528, shi0_529, \
                         shi1_528, shi1_529, shk_668, shk_672, \
                         shk_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_12 * shi0_528[k]
                   - f_13 * shi1_528[k]
                   + f_3 * pc_x[k] * shk_672[k];

        t_835[k] = f_12 * shi0_529[k]
                   - f_13 * shi1_529[k]
                   + f_3 * pc_x[k] * shk_673[k];

        t_836[k] = f_16 * sgk_488[k]
                   + f_3 * pc_y[k] * shk_668[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, t_840, t_841, t_842, pc_x, shi0_531, shi1_531, \
                         shk_675, shk_676, shk_677, shk_678, shk_679, \
                         shk_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_12 * shi0_531[k]
                   - f_13 * shi1_531[k]
                   + f_3 * pc_x[k] * shk_675[k];

        t_838[k] = f_3 * pc_x[k] * shk_676[k];

        t_839[k] = f_3 * pc_x[k] * shk_677[k];

        t_840[k] = f_3 * pc_x[k] * shk_678[k];

        t_841[k] = f_3 * pc_x[k] * shk_679[k];

        t_842[k] = f_3 * pc_x[k] * shk_680[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, t_846, t_847, pc_x, pc_y, pc_z, sgk_460, \
                         sgk_496, shi0_525, shi1_525, shk_676, shk_681, shk_682, \
                         shk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = f_3 * pc_x[k] * shk_681[k];

        t_844[k] = f_3 * pc_x[k] * shk_682[k];

        t_845[k] = f_3 * pc_x[k] * shk_683[k];

        t_846[k] = f_16 * sgk_496[k]
                   + f_1 * shi0_525[k]
                   - f_2 * shi1_525[k]
                   + f_3 * pc_y[k] * shk_676[k];

        t_847[k] = f_17 * sgk_460[k]
                   + f_3 * pc_z[k] * shk_676[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, pc_y, sgk_498, sgk_499, sgk_500, shi0_527, \
                         shi0_528, shi0_529, shi1_527, shi1_528, shi1_529, shk_678, shk_679, \
                         shk_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_16 * sgk_498[k]
                   + f_4 * shi0_527[k]
                   - f_5 * shi1_527[k]
                   + f_3 * pc_y[k] * shk_678[k];

        t_849[k] = f_16 * sgk_499[k]
                   + f_6 * shi0_528[k]
                   - f_7 * shi1_528[k]
                   + f_3 * pc_y[k] * shk_679[k];

        t_850[k] = f_16 * sgk_500[k]
                   + f_8 * shi0_529[k]
                   - f_9 * shi1_529[k]
                   + f_3 * pc_y[k] * shk_680[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, pc_y, sgk_501, sgk_502, sgk_503, shi0_530, \
                         shi0_531, shi1_530, shi1_531, shk_681, shk_682, \
                         shk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_16 * sgk_501[k]
                   + f_10 * shi0_530[k]
                   - f_11 * shi1_530[k]
                   + f_3 * pc_y[k] * shk_681[k];

        t_852[k] = f_16 * sgk_502[k]
                   + f_12 * shi0_531[k]
                   - f_13 * shi1_531[k]
                   + f_3 * pc_y[k] * shk_682[k];

        t_853[k] = f_16 * sgk_503[k]
                   + f_3 * pc_y[k] * shk_683[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pb_y, pc_y, pc_z, sgl0_630, sgk_467, \
                         sgk_468, sgk_504, sgl1_630, shi0_531, shi1_531, shk_683, \
                         shk_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = f_17 * sgk_467[k]
                   + f_1 * shi0_531[k]
                   - f_2 * shi1_531[k]
                   + f_3 * pc_z[k] * shk_683[k];

        t_855[k] = pb_y[k] * sgl0_630[k]
                   - f_14 * pc_y[k] * sgl1_630[k];

        t_856[k] = f_15 * sgk_504[k]
                   + f_3 * pc_y[k] * shk_684[k];

        t_857[k] = f_18 * sgk_468[k]
                   + f_3 * pc_z[k] * shk_684[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, pb_y, pc_x, pc_y, sgl0_635, sgk_506, sgl1_635, \
                         shi0_535, shi1_535, shk_686, shk_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_4 * shi0_535[k]
                   - f_5 * shi1_535[k]
                   + f_3 * pc_x[k] * shk_687[k];

        t_859[k] = f_15 * sgk_506[k]
                   + f_3 * pc_y[k] * shk_686[k];

        t_860[k] = pb_y[k] * sgl0_635[k]
                   - f_14 * pc_y[k] * sgl1_635[k];
    }

#pragma omp simd aligned(t_861, t_862, t_863, pc_x, pc_y, pc_z, sgk_471, sgk_509, shi0_538, \
                         shi1_538, shk_687, shk_689, shk_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_861[k] = f_6 * shi0_538[k]
                   - f_7 * shi1_538[k]
                   + f_3 * pc_x[k] * shk_690[k];

        t_862[k] = f_18 * sgk_471[k]
                   + f_3 * pc_z[k] * shk_687[k];

        t_863[k] = f_15 * sgk_509[k]
                   + f_3 * pc_y[k] * shk_689[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, pb_y, pc_x, pc_y, pc_z, sgl0_639, sgk_474, \
                         sgl1_639, shi0_542, shi1_542, shk_690, \
                         shk_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = pb_y[k] * sgl0_639[k]
                   - f_14 * pc_y[k] * sgl1_639[k];

        t_865[k] = f_8 * shi0_542[k]
                   - f_9 * shi1_542[k]
                   + f_3 * pc_x[k] * shk_694[k];

        t_866[k] = f_18 * sgk_474[k]
                   + f_3 * pc_z[k] * shk_690[k];
    }

#pragma omp simd aligned(t_867, t_868, t_869, pb_y, pc_x, pc_y, sgl0_644, sgk_513, sgl1_644, \
                         shi0_544, shi1_544, shk_693, shk_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_867[k] = f_8 * shi0_544[k]
                   - f_9 * shi1_544[k]
                   + f_3 * pc_x[k] * shk_696[k];

        t_868[k] = f_15 * sgk_513[k]
                   + f_3 * pc_y[k] * shk_693[k];

        t_869[k] = pb_y[k] * sgl0_644[k]
                   - f_14 * pc_y[k] * sgl1_644[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, pc_x, pc_z, sgk_478, shi0_547, shi0_549, \
                         shi1_547, shi1_549, shk_694, shk_699, \
                         shk_701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = f_10 * shi0_547[k]
                   - f_11 * shi1_547[k]
                   + f_3 * pc_x[k] * shk_699[k];

        t_871[k] = f_18 * sgk_478[k]
                   + f_3 * pc_z[k] * shk_694[k];

        t_872[k] = f_10 * shi0_549[k]
                   - f_11 * shi1_549[k]
                   + f_3 * pc_x[k] * shk_701[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pb_y, pc_x, pc_y, sgl0_650, sgk_518, sgl1_650, \
                         shi0_550, shi1_550, shk_698, shk_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_10 * shi0_550[k]
                   - f_11 * shi1_550[k]
                   + f_3 * pc_x[k] * shk_702[k];

        t_874[k] = f_15 * sgk_518[k]
                   + f_3 * pc_y[k] * shk_698[k];

        t_875[k] = pb_y[k] * sgl0_650[k]
                   - f_14 * pc_y[k] * sgl1_650[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pc_x, pc_z, sgk_483, shi0_553, shi0_555, \
                         shi1_553, shi1_555, shk_699, shk_705, \
                         shk_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_12 * shi0_553[k]
                   - f_13 * shi1_553[k]
                   + f_3 * pc_x[k] * shk_705[k];

        t_877[k] = f_18 * sgk_483[k]
                   + f_3 * pc_z[k] * shk_699[k];

        t_878[k] = f_12 * shi0_555[k]
                   - f_13 * shi1_555[k]
                   + f_3 * pc_x[k] * shk_707[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pc_x, pc_y, sgk_524, shi0_556, shi0_557, \
                         shi1_556, shi1_557, shk_704, shk_708, \
                         shk_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_12 * shi0_556[k]
                   - f_13 * shi1_556[k]
                   + f_3 * pc_x[k] * shk_708[k];

        t_880[k] = f_12 * shi0_557[k]
                   - f_13 * shi1_557[k]
                   + f_3 * pc_x[k] * shk_709[k];

        t_881[k] = f_15 * sgk_524[k]
                   + f_3 * pc_y[k] * shk_704[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, t_886, t_887, pb_y, pc_x, pc_y, sgl0_657, \
                         sgl1_657, shk_712, shk_713, shk_714, shk_715, \
                         shk_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = pb_y[k] * sgl0_657[k]
                   - f_14 * pc_y[k] * sgl1_657[k];

        t_883[k] = f_3 * pc_x[k] * shk_712[k];

        t_884[k] = f_3 * pc_x[k] * shk_713[k];

        t_885[k] = f_3 * pc_x[k] * shk_714[k];

        t_886[k] = f_3 * pc_x[k] * shk_715[k];

        t_887[k] = f_3 * pc_x[k] * shk_716[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, pb_y, pc_x, pc_y, sgl0_666, sgk_532, \
                         sgl1_666, shk_717, shk_718, shk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = f_3 * pc_x[k] * shk_717[k];

        t_889[k] = f_3 * pc_x[k] * shk_718[k];

        t_890[k] = f_3 * pc_x[k] * shk_719[k];

        t_891[k] = pb_y[k] * sgl0_666[k]
                   + f_20 * sgk_532[k]
                   - f_14 * pc_y[k] * sgl1_666[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, pb_y, pc_y, pc_z, sgl0_668, sgl0_669, sgk_496, \
                         sgk_534, sgk_535, sgl1_668, sgl1_669, \
                         shk_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = f_18 * sgk_496[k]
                   + f_3 * pc_z[k] * shk_712[k];

        t_893[k] = pb_y[k] * sgl0_668[k]
                   + f_19 * sgk_534[k]
                   - f_14 * pc_y[k] * sgl1_668[k];

        t_894[k] = pb_y[k] * sgl0_669[k]
                   + f_0 * sgk_535[k]
                   - f_14 * pc_y[k] * sgl1_669[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, pb_y, pc_y, sgl0_670, sgl0_671, sgl0_672, \
                         sgk_536, sgk_537, sgk_538, sgl1_670, sgl1_671, \
                         sgl1_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = pb_y[k] * sgl0_670[k]
                   + f_18 * sgk_536[k]
                   - f_14 * pc_y[k] * sgl1_670[k];

        t_896[k] = pb_y[k] * sgl0_671[k]
                   + f_17 * sgk_537[k]
                   - f_14 * pc_y[k] * sgl1_671[k];

        t_897[k] = pb_y[k] * sgl0_672[k]
                   + f_16 * sgk_538[k]
                   - f_14 * pc_y[k] * sgl1_672[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, t_901, pb_y, pc_x, pc_y, sgl0_674, sgk_539, \
                         sgl1_674, shi0_560, shi1_560, shk_719, \
                         shk_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_15 * sgk_539[k]
                   + f_3 * pc_y[k] * shk_719[k];

        t_899[k] = pb_y[k] * sgl0_674[k]
                   - f_14 * pc_y[k] * sgl1_674[k];

        t_900[k] = f_1 * shi0_560[k]
                   - f_2 * shi1_560[k]
                   + f_3 * pc_x[k] * shk_720[k];

        t_901[k] = f_3 * pc_y[k] * shk_720[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, t_905, pc_x, pc_y, pc_z, sgk_504, shi0_563, \
                         shi0_565, shi1_563, shi1_565, shk_720, shk_722, shk_723, \
                         shk_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_0 * sgk_504[k]
                   + f_3 * pc_z[k] * shk_720[k];

        t_903[k] = f_4 * shi0_563[k]
                   - f_5 * shi1_563[k]
                   + f_3 * pc_x[k] * shk_723[k];

        t_904[k] = f_3 * pc_y[k] * shk_722[k];

        t_905[k] = f_4 * shi0_565[k]
                   - f_5 * shi1_565[k]
                   + f_3 * pc_x[k] * shk_725[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, t_909, pc_x, pc_y, pc_z, sgk_507, shi0_566, \
                         shi0_569, shi1_566, shi1_569, shk_723, shk_725, shk_726, \
                         shk_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_6 * shi0_566[k]
                   - f_7 * shi1_566[k]
                   + f_3 * pc_x[k] * shk_726[k];

        t_907[k] = f_0 * sgk_507[k]
                   + f_3 * pc_z[k] * shk_723[k];

        t_908[k] = f_3 * pc_y[k] * shk_725[k];

        t_909[k] = f_6 * shi0_569[k]
                   - f_7 * shi1_569[k]
                   + f_3 * pc_x[k] * shk_729[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, pc_x, pc_y, pc_z, sgk_510, shi0_570, \
                         shi0_572, shi1_570, shi1_572, shk_726, shk_729, shk_730, \
                         shk_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = f_8 * shi0_570[k]
                   - f_9 * shi1_570[k]
                   + f_3 * pc_x[k] * shk_730[k];

        t_911[k] = f_0 * sgk_510[k]
                   + f_3 * pc_z[k] * shk_726[k];

        t_912[k] = f_8 * shi0_572[k]
                   - f_9 * shi1_572[k]
                   + f_3 * pc_x[k] * shk_732[k];

        t_913[k] = f_3 * pc_y[k] * shk_729[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, pc_x, pc_z, sgk_514, shi0_574, shi0_575, \
                         shi1_574, shi1_575, shk_730, shk_734, \
                         shk_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_8 * shi0_574[k]
                   - f_9 * shi1_574[k]
                   + f_3 * pc_x[k] * shk_734[k];

        t_915[k] = f_10 * shi0_575[k]
                   - f_11 * shi1_575[k]
                   + f_3 * pc_x[k] * shk_735[k];

        t_916[k] = f_0 * sgk_514[k]
                   + f_3 * pc_z[k] * shk_730[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, t_920, pc_x, pc_y, shi0_577, shi0_578, shi0_580, \
                         shi1_577, shi1_578, shi1_580, shk_734, shk_737, shk_738, \
                         shk_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = f_10 * shi0_577[k]
                   - f_11 * shi1_577[k]
                   + f_3 * pc_x[k] * shk_737[k];

        t_918[k] = f_10 * shi0_578[k]
                   - f_11 * shi1_578[k]
                   + f_3 * pc_x[k] * shk_738[k];

        t_919[k] = f_3 * pc_y[k] * shk_734[k];

        t_920[k] = f_10 * shi0_580[k]
                   - f_11 * shi1_580[k]
                   + f_3 * pc_x[k] * shk_740[k];
    }

#pragma omp simd aligned(t_921, t_922, t_923, pc_x, pc_z, sgk_519, shi0_581, shi0_583, \
                         shi1_581, shi1_583, shk_735, shk_741, \
                         shk_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_12 * shi0_581[k]
                   - f_13 * shi1_581[k]
                   + f_3 * pc_x[k] * shk_741[k];

        t_922[k] = f_0 * sgk_519[k]
                   + f_3 * pc_z[k] * shk_735[k];

        t_923[k] = f_12 * shi0_583[k]
                   - f_13 * shi1_583[k]
                   + f_3 * pc_x[k] * shk_743[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pc_x, pc_y, shi0_584, shi0_585, shi0_587, \
                         shi1_584, shi1_585, shi1_587, shk_740, shk_744, shk_745, \
                         shk_747 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_12 * shi0_584[k]
                   - f_13 * shi1_584[k]
                   + f_3 * pc_x[k] * shk_744[k];

        t_925[k] = f_12 * shi0_585[k]
                   - f_13 * shi1_585[k]
                   + f_3 * pc_x[k] * shk_745[k];

        t_926[k] = f_3 * pc_y[k] * shk_740[k];

        t_927[k] = f_12 * shi0_587[k]
                   - f_13 * shi1_587[k]
                   + f_3 * pc_x[k] * shk_747[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, t_931, t_932, t_933, t_934, pc_x, shk_748, \
                         shk_749, shk_750, shk_751, shk_752, shk_753, \
                         shk_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_3 * pc_x[k] * shk_748[k];

        t_929[k] = f_3 * pc_x[k] * shk_749[k];

        t_930[k] = f_3 * pc_x[k] * shk_750[k];

        t_931[k] = f_3 * pc_x[k] * shk_751[k];

        t_932[k] = f_3 * pc_x[k] * shk_752[k];

        t_933[k] = f_3 * pc_x[k] * shk_753[k];

        t_934[k] = f_3 * pc_x[k] * shk_754[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, pc_x, pc_y, pc_z, sgk_532, shi0_581, \
                         shi0_583, shi1_581, shi1_583, shk_748, shk_750, \
                         shk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_3 * pc_x[k] * shk_755[k];

        t_936[k] = f_1 * shi0_581[k]
                   - f_2 * shi1_581[k]
                   + f_3 * pc_y[k] * shk_748[k];

        t_937[k] = f_0 * sgk_532[k]
                   + f_3 * pc_z[k] * shk_748[k];

        t_938[k] = f_4 * shi0_583[k]
                   - f_5 * shi1_583[k]
                   + f_3 * pc_y[k] * shk_750[k];
    }
}

static auto
compute_prim_shl_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t sgk, const size_t shi0,
                                                          const size_t shi1, const size_t shk,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);

    auto *t_939 = buffer.data(target + 939);
    auto *t_940 = buffer.data(target + 940);
    auto *t_941 = buffer.data(target + 941);
    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);

    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgk_539 = buffer.data(sgk + 539);

    const auto *shi0_584 = buffer.data(shi0 + 584);
    const auto *shi0_585 = buffer.data(shi0 + 585);
    const auto *shi0_586 = buffer.data(shi0 + 586);
    const auto *shi0_587 = buffer.data(shi0 + 587);

    const auto *shi1_584 = buffer.data(shi1 + 584);
    const auto *shi1_585 = buffer.data(shi1 + 585);
    const auto *shi1_586 = buffer.data(shi1 + 586);
    const auto *shi1_587 = buffer.data(shi1 + 587);

    const auto *shk_751 = buffer.data(shk + 751);
    const auto *shk_752 = buffer.data(shk + 752);
    const auto *shk_753 = buffer.data(shk + 753);
    const auto *shk_754 = buffer.data(shk + 754);
    const auto *shk_755 = buffer.data(shk + 755);

#pragma omp simd aligned(t_939, t_940, t_941, pc_y, shi0_584, shi0_585, shi0_586, shi1_584, \
                         shi1_585, shi1_586, shk_751, shk_752, \
                         shk_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_6 * shi0_584[k]
                   - f_7 * shi1_584[k]
                   + f_3 * pc_y[k] * shk_751[k];

        t_940[k] = f_8 * shi0_585[k]
                   - f_9 * shi1_585[k]
                   + f_3 * pc_y[k] * shk_752[k];

        t_941[k] = f_10 * shi0_586[k]
                   - f_11 * shi1_586[k]
                   + f_3 * pc_y[k] * shk_753[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, pc_y, pc_z, sgk_539, shi0_587, shi1_587, \
                         shk_754, shk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_12 * shi0_587[k]
                   - f_13 * shi1_587[k]
                   + f_3 * pc_y[k] * shk_754[k];

        t_943[k] = f_3 * pc_y[k] * shk_755[k];

        t_944[k] = f_0 * sgk_539[k]
                   + f_1 * shi0_587[k]
                   - f_2 * shi1_587[k]
                   + f_3 * pc_z[k] * shk_755[k];
    }
}

auto
compute_prim_shl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sgl0, const size_t sgk,
                                                   const size_t sgl1, const size_t shi0,
                                                   const size_t shi1, const size_t shk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_shl_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sgl0, sgk,
                                                              sgl1, shi0, shi1, shk, ncols,
                                                              gamma, p, q);

    compute_prim_shl_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sgl0, sgk,
                                                              sgl1, shi0, shi1, shk, ncols,
                                                              gamma, p, q);

    compute_prim_shl_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sgl0, sgk,
                                                              sgl1, shi0, shi1, shk, ncols,
                                                              gamma, p, q);

    compute_prim_shl_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sgl0, sgk,
                                                              sgl1, shi0, shi1, shk, ncols,
                                                              gamma, p, q);

    compute_prim_shl_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sgl0, sgk,
                                                              sgl1, shk, ncols, gamma, p, q);

    compute_prim_shl_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, sgl0, sgk,
                                                              sgl1, shi0, shi1, shk, ncols,
                                                              gamma, p, q);

    compute_prim_shl_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, sgl0, sgk,
                                                              sgl1, shi0, shi1, shk, ncols,
                                                              gamma, p, q);

    compute_prim_shl_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, sgl0, sgk,
                                                              sgl1, shi0, shi1, shk, ncols,
                                                              gamma, p, q);

    compute_prim_shl_three_center_electron_repulsion_0_piece8(buffer, target, pc, sgk, shi0,
                                                              shi1, shk, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
