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


#include "SimdThreeCenterElectronRepulsionVrrRecSFL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sfl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sdl0,
                                                          const size_t sdk, const size_t sdl1,
                                                          const size_t sfi0, const size_t sfi1,
                                                          const size_t sfk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_17 = 2.0 / q;
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

    const auto *sdl0_0 = buffer.data(sdl0 + 0);
    const auto *sdl0_3 = buffer.data(sdl0 + 3);
    const auto *sdl0_5 = buffer.data(sdl0 + 5);
    const auto *sdl0_6 = buffer.data(sdl0 + 6);
    const auto *sdl0_9 = buffer.data(sdl0 + 9);
    const auto *sdl0_10 = buffer.data(sdl0 + 10);
    const auto *sdl0_12 = buffer.data(sdl0 + 12);
    const auto *sdl0_14 = buffer.data(sdl0 + 14);
    const auto *sdl0_15 = buffer.data(sdl0 + 15);
    const auto *sdl0_17 = buffer.data(sdl0 + 17);
    const auto *sdl0_18 = buffer.data(sdl0 + 18);
    const auto *sdl0_20 = buffer.data(sdl0 + 20);
    const auto *sdl0_21 = buffer.data(sdl0 + 21);
    const auto *sdl0_23 = buffer.data(sdl0 + 23);
    const auto *sdl0_24 = buffer.data(sdl0 + 24);
    const auto *sdl0_25 = buffer.data(sdl0 + 25);
    const auto *sdl0_27 = buffer.data(sdl0 + 27);
    const auto *sdl0_44 = buffer.data(sdl0 + 44);

    const auto *sdk_0 = buffer.data(sdk + 0);
    const auto *sdk_1 = buffer.data(sdk + 1);
    const auto *sdk_2 = buffer.data(sdk + 2);
    const auto *sdk_3 = buffer.data(sdk + 3);
    const auto *sdk_5 = buffer.data(sdk + 5);
    const auto *sdk_6 = buffer.data(sdk + 6);
    const auto *sdk_7 = buffer.data(sdk + 7);
    const auto *sdk_8 = buffer.data(sdk + 8);
    const auto *sdk_9 = buffer.data(sdk + 9);
    const auto *sdk_10 = buffer.data(sdk + 10);
    const auto *sdk_11 = buffer.data(sdk + 11);
    const auto *sdk_12 = buffer.data(sdk + 12);
    const auto *sdk_13 = buffer.data(sdk + 13);
    const auto *sdk_14 = buffer.data(sdk + 14);
    const auto *sdk_15 = buffer.data(sdk + 15);
    const auto *sdk_16 = buffer.data(sdk + 16);
    const auto *sdk_17 = buffer.data(sdk + 17);
    const auto *sdk_18 = buffer.data(sdk + 18);
    const auto *sdk_19 = buffer.data(sdk + 19);
    const auto *sdk_20 = buffer.data(sdk + 20);
    const auto *sdk_21 = buffer.data(sdk + 21);
    const auto *sdk_23 = buffer.data(sdk + 23);
    const auto *sdk_24 = buffer.data(sdk + 24);
    const auto *sdk_25 = buffer.data(sdk + 25);
    const auto *sdk_27 = buffer.data(sdk + 27);
    const auto *sdk_28 = buffer.data(sdk + 28);
    const auto *sdk_29 = buffer.data(sdk + 29);
    const auto *sdk_30 = buffer.data(sdk + 30);
    const auto *sdk_31 = buffer.data(sdk + 31);
    const auto *sdk_32 = buffer.data(sdk + 32);
    const auto *sdk_33 = buffer.data(sdk + 33);
    const auto *sdk_34 = buffer.data(sdk + 34);
    const auto *sdk_35 = buffer.data(sdk + 35);
    const auto *sdk_64 = buffer.data(sdk + 64);
    const auto *sdk_65 = buffer.data(sdk + 65);
    const auto *sdk_66 = buffer.data(sdk + 66);
    const auto *sdk_67 = buffer.data(sdk + 67);
    const auto *sdk_68 = buffer.data(sdk + 68);
    const auto *sdk_69 = buffer.data(sdk + 69);
    const auto *sdk_70 = buffer.data(sdk + 70);
    const auto *sdk_71 = buffer.data(sdk + 71);

    const auto *sdl1_0 = buffer.data(sdl1 + 0);
    const auto *sdl1_3 = buffer.data(sdl1 + 3);
    const auto *sdl1_5 = buffer.data(sdl1 + 5);
    const auto *sdl1_6 = buffer.data(sdl1 + 6);
    const auto *sdl1_9 = buffer.data(sdl1 + 9);
    const auto *sdl1_10 = buffer.data(sdl1 + 10);
    const auto *sdl1_12 = buffer.data(sdl1 + 12);
    const auto *sdl1_14 = buffer.data(sdl1 + 14);
    const auto *sdl1_15 = buffer.data(sdl1 + 15);
    const auto *sdl1_17 = buffer.data(sdl1 + 17);
    const auto *sdl1_18 = buffer.data(sdl1 + 18);
    const auto *sdl1_20 = buffer.data(sdl1 + 20);
    const auto *sdl1_21 = buffer.data(sdl1 + 21);
    const auto *sdl1_23 = buffer.data(sdl1 + 23);
    const auto *sdl1_24 = buffer.data(sdl1 + 24);
    const auto *sdl1_25 = buffer.data(sdl1 + 25);
    const auto *sdl1_27 = buffer.data(sdl1 + 27);
    const auto *sdl1_44 = buffer.data(sdl1 + 44);

    const auto *sfi0_0 = buffer.data(sfi0 + 0);
    const auto *sfi0_3 = buffer.data(sfi0 + 3);
    const auto *sfi0_5 = buffer.data(sfi0 + 5);
    const auto *sfi0_6 = buffer.data(sfi0 + 6);
    const auto *sfi0_9 = buffer.data(sfi0 + 9);
    const auto *sfi0_10 = buffer.data(sfi0 + 10);
    const auto *sfi0_12 = buffer.data(sfi0 + 12);
    const auto *sfi0_14 = buffer.data(sfi0 + 14);
    const auto *sfi0_15 = buffer.data(sfi0 + 15);
    const auto *sfi0_17 = buffer.data(sfi0 + 17);
    const auto *sfi0_18 = buffer.data(sfi0 + 18);
    const auto *sfi0_20 = buffer.data(sfi0 + 20);
    const auto *sfi0_21 = buffer.data(sfi0 + 21);
    const auto *sfi0_23 = buffer.data(sfi0 + 23);
    const auto *sfi0_24 = buffer.data(sfi0 + 24);
    const auto *sfi0_25 = buffer.data(sfi0 + 25);
    const auto *sfi0_26 = buffer.data(sfi0 + 26);
    const auto *sfi0_27 = buffer.data(sfi0 + 27);
    const auto *sfi0_49 = buffer.data(sfi0 + 49);
    const auto *sfi0_51 = buffer.data(sfi0 + 51);
    const auto *sfi0_52 = buffer.data(sfi0 + 52);
    const auto *sfi0_53 = buffer.data(sfi0 + 53);
    const auto *sfi0_54 = buffer.data(sfi0 + 54);
    const auto *sfi0_55 = buffer.data(sfi0 + 55);

    const auto *sfi1_0 = buffer.data(sfi1 + 0);
    const auto *sfi1_3 = buffer.data(sfi1 + 3);
    const auto *sfi1_5 = buffer.data(sfi1 + 5);
    const auto *sfi1_6 = buffer.data(sfi1 + 6);
    const auto *sfi1_9 = buffer.data(sfi1 + 9);
    const auto *sfi1_10 = buffer.data(sfi1 + 10);
    const auto *sfi1_12 = buffer.data(sfi1 + 12);
    const auto *sfi1_14 = buffer.data(sfi1 + 14);
    const auto *sfi1_15 = buffer.data(sfi1 + 15);
    const auto *sfi1_17 = buffer.data(sfi1 + 17);
    const auto *sfi1_18 = buffer.data(sfi1 + 18);
    const auto *sfi1_20 = buffer.data(sfi1 + 20);
    const auto *sfi1_21 = buffer.data(sfi1 + 21);
    const auto *sfi1_23 = buffer.data(sfi1 + 23);
    const auto *sfi1_24 = buffer.data(sfi1 + 24);
    const auto *sfi1_25 = buffer.data(sfi1 + 25);
    const auto *sfi1_26 = buffer.data(sfi1 + 26);
    const auto *sfi1_27 = buffer.data(sfi1 + 27);
    const auto *sfi1_49 = buffer.data(sfi1 + 49);
    const auto *sfi1_51 = buffer.data(sfi1 + 51);
    const auto *sfi1_52 = buffer.data(sfi1 + 52);
    const auto *sfi1_53 = buffer.data(sfi1 + 53);
    const auto *sfi1_54 = buffer.data(sfi1 + 54);
    const auto *sfi1_55 = buffer.data(sfi1 + 55);

    const auto *sfk_0 = buffer.data(sfk + 0);
    const auto *sfk_2 = buffer.data(sfk + 2);
    const auto *sfk_3 = buffer.data(sfk + 3);
    const auto *sfk_5 = buffer.data(sfk + 5);
    const auto *sfk_6 = buffer.data(sfk + 6);
    const auto *sfk_9 = buffer.data(sfk + 9);
    const auto *sfk_10 = buffer.data(sfk + 10);
    const auto *sfk_12 = buffer.data(sfk + 12);
    const auto *sfk_14 = buffer.data(sfk + 14);
    const auto *sfk_15 = buffer.data(sfk + 15);
    const auto *sfk_17 = buffer.data(sfk + 17);
    const auto *sfk_18 = buffer.data(sfk + 18);
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
    const auto *sfk_65 = buffer.data(sfk + 65);
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
    const auto *sfk_78 = buffer.data(sfk + 78);
    const auto *sfk_81 = buffer.data(sfk + 81);
    const auto *sfk_82 = buffer.data(sfk + 82);
    const auto *sfk_86 = buffer.data(sfk + 86);
    const auto *sfk_87 = buffer.data(sfk + 87);
    const auto *sfk_92 = buffer.data(sfk + 92);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sdk_0, sdk_3, sfi0_0, sfi0_3, \
                         sfi1_0, sfi1_3, sfk_0, sfk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sdk_0[k]
                 + f_1 * sfi0_0[k]
                 - f_2 * sfi1_0[k]
                 + f_3 * pc_x[k] * sfk_0[k];

        t_1[k] = f_3 * pc_y[k] * sfk_0[k];

        t_2[k] = f_3 * pc_z[k] * sfk_0[k];

        t_3[k] = f_0 * sdk_3[k]
                 + f_4 * sfi0_3[k]
                 - f_5 * sfi1_3[k]
                 + f_3 * pc_x[k] * sfk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sdk_5, sdk_6, sfi0_5, sfi0_6, sfi1_5, \
                         sfi1_6, sfk_2, sfk_5, sfk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sfk_2[k];

        t_5[k] = f_0 * sdk_5[k]
                 + f_4 * sfi0_5[k]
                 - f_5 * sfi1_5[k]
                 + f_3 * pc_x[k] * sfk_5[k];

        t_6[k] = f_0 * sdk_6[k]
                 + f_6 * sfi0_6[k]
                 - f_7 * sfi1_6[k]
                 + f_3 * pc_x[k] * sfk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sdk_9, sfi0_9, sfi1_9, sfk_3, sfk_5, \
                         sfk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sfk_3[k];

        t_8[k] = f_3 * pc_y[k] * sfk_5[k];

        t_9[k] = f_0 * sdk_9[k]
                 + f_6 * sfi0_9[k]
                 - f_7 * sfi1_9[k]
                 + f_3 * pc_x[k] * sfk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sdk_10, sdk_12, sfi0_10, sfi0_12, \
                         sfi1_10, sfi1_12, sfk_6, sfk_10, sfk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sdk_10[k]
                  + f_8 * sfi0_10[k]
                  - f_9 * sfi1_10[k]
                  + f_3 * pc_x[k] * sfk_10[k];

        t_11[k] = f_3 * pc_z[k] * sfk_6[k];

        t_12[k] = f_0 * sdk_12[k]
                  + f_8 * sfi0_12[k]
                  - f_9 * sfi1_12[k]
                  + f_3 * pc_x[k] * sfk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, sdk_14, sdk_15, sfi0_14, sfi0_15, \
                         sfi1_14, sfi1_15, sfk_9, sfk_14, sfk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sfk_9[k];

        t_14[k] = f_0 * sdk_14[k]
                  + f_8 * sfi0_14[k]
                  - f_9 * sfi1_14[k]
                  + f_3 * pc_x[k] * sfk_14[k];

        t_15[k] = f_0 * sdk_15[k]
                  + f_10 * sfi0_15[k]
                  - f_11 * sfi1_15[k]
                  + f_3 * pc_x[k] * sfk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, sdk_17, sdk_18, sfi0_17, sfi0_18, \
                         sfi1_17, sfi1_18, sfk_10, sfk_17, sfk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * sfk_10[k];

        t_17[k] = f_0 * sdk_17[k]
                  + f_10 * sfi0_17[k]
                  - f_11 * sfi1_17[k]
                  + f_3 * pc_x[k] * sfk_17[k];

        t_18[k] = f_0 * sdk_18[k]
                  + f_10 * sfi0_18[k]
                  - f_11 * sfi1_18[k]
                  + f_3 * pc_x[k] * sfk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pc_x, pc_y, sdk_20, sdk_21, sfi0_20, sfi0_21, \
                         sfi1_20, sfi1_21, sfk_14, sfk_20, sfk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sfk_14[k];

        t_20[k] = f_0 * sdk_20[k]
                  + f_10 * sfi0_20[k]
                  - f_11 * sfi1_20[k]
                  + f_3 * pc_x[k] * sfk_20[k];

        t_21[k] = f_0 * sdk_21[k]
                  + f_12 * sfi0_21[k]
                  - f_13 * sfi1_21[k]
                  + f_3 * pc_x[k] * sfk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pc_x, pc_z, sdk_23, sdk_24, sfi0_23, sfi0_24, \
                         sfi1_23, sfi1_24, sfk_15, sfk_23, sfk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * pc_z[k] * sfk_15[k];

        t_23[k] = f_0 * sdk_23[k]
                  + f_12 * sfi0_23[k]
                  - f_13 * sfi1_23[k]
                  + f_3 * pc_x[k] * sfk_23[k];

        t_24[k] = f_0 * sdk_24[k]
                  + f_12 * sfi0_24[k]
                  - f_13 * sfi1_24[k]
                  + f_3 * pc_x[k] * sfk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pc_x, pc_y, sdk_25, sdk_27, sfi0_25, sfi0_27, \
                         sfi1_25, sfi1_27, sfk_20, sfk_25, sfk_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * sdk_25[k]
                  + f_12 * sfi0_25[k]
                  - f_13 * sfi1_25[k]
                  + f_3 * pc_x[k] * sfk_25[k];

        t_26[k] = f_3 * pc_y[k] * sfk_20[k];

        t_27[k] = f_0 * sdk_27[k]
                  + f_12 * sfi0_27[k]
                  - f_13 * sfi1_27[k]
                  + f_3 * pc_x[k] * sfk_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pc_x, sdk_28, sdk_29, sdk_30, sdk_31, \
                         sdk_32, sfk_28, sfk_29, sfk_30, sfk_31, \
                         sfk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * sdk_28[k]
                  + f_3 * pc_x[k] * sfk_28[k];

        t_29[k] = f_0 * sdk_29[k]
                  + f_3 * pc_x[k] * sfk_29[k];

        t_30[k] = f_0 * sdk_30[k]
                  + f_3 * pc_x[k] * sfk_30[k];

        t_31[k] = f_0 * sdk_31[k]
                  + f_3 * pc_x[k] * sfk_31[k];

        t_32[k] = f_0 * sdk_32[k]
                  + f_3 * pc_x[k] * sfk_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pc_x, pc_y, sdk_33, sdk_34, sdk_35, sfi0_21, \
                         sfi1_21, sfk_28, sfk_33, sfk_34, sfk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * sdk_33[k]
                  + f_3 * pc_x[k] * sfk_33[k];

        t_34[k] = f_0 * sdk_34[k]
                  + f_3 * pc_x[k] * sfk_34[k];

        t_35[k] = f_0 * sdk_35[k]
                  + f_3 * pc_x[k] * sfk_35[k];

        t_36[k] = f_1 * sfi0_21[k]
                  - f_2 * sfi1_21[k]
                  + f_3 * pc_y[k] * sfk_28[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pc_y, pc_z, sfi0_23, sfi0_24, sfi0_25, \
                         sfi1_23, sfi1_24, sfi1_25, sfk_28, sfk_30, sfk_31, \
                         sfk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * pc_z[k] * sfk_28[k];

        t_38[k] = f_4 * sfi0_23[k]
                  - f_5 * sfi1_23[k]
                  + f_3 * pc_y[k] * sfk_30[k];

        t_39[k] = f_6 * sfi0_24[k]
                  - f_7 * sfi1_24[k]
                  + f_3 * pc_y[k] * sfk_31[k];

        t_40[k] = f_8 * sfi0_25[k]
                  - f_9 * sfi1_25[k]
                  + f_3 * pc_y[k] * sfk_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, sfi0_26, sfi0_27, sfi1_26, \
                         sfi1_27, sfk_33, sfk_34, sfk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * sfi0_26[k]
                  - f_11 * sfi1_26[k]
                  + f_3 * pc_y[k] * sfk_33[k];

        t_42[k] = f_12 * sfi0_27[k]
                  - f_13 * sfi1_27[k]
                  + f_3 * pc_y[k] * sfk_34[k];

        t_43[k] = f_3 * pc_y[k] * sfk_35[k];

        t_44[k] = f_1 * sfi0_27[k]
                  - f_2 * sfi1_27[k]
                  + f_3 * pc_z[k] * sfk_35[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_y, pc_y, pc_z, sdl0_0, sdl0_3, sdk_0, \
                         sdk_1, sdl1_0, sdl1_3, sfk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pb_y[k] * sdl0_0[k]
                  - f_14 * pc_y[k] * sdl1_0[k];

        t_46[k] = f_15 * sdk_0[k]
                  + f_3 * pc_y[k] * sfk_36[k];

        t_47[k] = f_3 * pc_z[k] * sfk_36[k];

        t_48[k] = pb_y[k] * sdl0_3[k]
                  + f_16 * sdk_1[k]
                  - f_14 * pc_y[k] * sdl1_3[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_y, pc_y, pc_z, sdl0_5, sdl0_6, sdk_2, \
                         sdk_3, sdl1_5, sdl1_6, sfk_38, sfk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_15 * sdk_2[k]
                  + f_3 * pc_y[k] * sfk_38[k];

        t_50[k] = pb_y[k] * sdl0_5[k]
                  - f_14 * pc_y[k] * sdl1_5[k];

        t_51[k] = pb_y[k] * sdl0_6[k]
                  + f_0 * sdk_3[k]
                  - f_14 * pc_y[k] * sdl1_6[k];

        t_52[k] = f_3 * pc_z[k] * sfk_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_y, pc_y, pc_z, sdl0_9, sdl0_10, sdk_5, \
                         sdk_6, sdl1_9, sdl1_10, sfk_41, sfk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_15 * sdk_5[k]
                  + f_3 * pc_y[k] * sfk_41[k];

        t_54[k] = pb_y[k] * sdl0_9[k]
                  - f_14 * pc_y[k] * sdl1_9[k];

        t_55[k] = pb_y[k] * sdl0_10[k]
                  + f_17 * sdk_6[k]
                  - f_14 * pc_y[k] * sdl1_10[k];

        t_56[k] = f_3 * pc_z[k] * sfk_42[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, sdl0_12, sdl0_14, sdl0_15, sdk_8, \
                         sdk_9, sdk_10, sdl1_12, sdl1_14, sdl1_15, \
                         sfk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_y[k] * sdl0_12[k]
                  + f_16 * sdk_8[k]
                  - f_14 * pc_y[k] * sdl1_12[k];

        t_58[k] = f_15 * sdk_9[k]
                  + f_3 * pc_y[k] * sfk_45[k];

        t_59[k] = pb_y[k] * sdl0_14[k]
                  - f_14 * pc_y[k] * sdl1_14[k];

        t_60[k] = pb_y[k] * sdl0_15[k]
                  + f_18 * sdk_10[k]
                  - f_14 * pc_y[k] * sdl1_15[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_y, pc_y, pc_z, sdl0_17, sdl0_18, sdk_12, \
                         sdk_13, sdk_14, sdl1_17, sdl1_18, sfk_46, \
                         sfk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * pc_z[k] * sfk_46[k];

        t_62[k] = pb_y[k] * sdl0_17[k]
                  + f_0 * sdk_12[k]
                  - f_14 * pc_y[k] * sdl1_17[k];

        t_63[k] = pb_y[k] * sdl0_18[k]
                  + f_16 * sdk_13[k]
                  - f_14 * pc_y[k] * sdl1_18[k];

        t_64[k] = f_15 * sdk_14[k]
                  + f_3 * pc_y[k] * sfk_50[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_y, pc_z, sdl0_20, sdl0_21, sdl0_23, \
                         sdk_15, sdk_17, sdl1_20, sdl1_21, sdl1_23, \
                         sfk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * sdl0_20[k]
                  - f_14 * pc_y[k] * sdl1_20[k];

        t_66[k] = pb_y[k] * sdl0_21[k]
                  + f_19 * sdk_15[k]
                  - f_14 * pc_y[k] * sdl1_21[k];

        t_67[k] = f_3 * pc_z[k] * sfk_51[k];

        t_68[k] = pb_y[k] * sdl0_23[k]
                  + f_17 * sdk_17[k]
                  - f_14 * pc_y[k] * sdl1_23[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_y, pc_y, sdl0_24, sdl0_25, sdl0_27, \
                         sdk_18, sdk_19, sdk_20, sdl1_24, sdl1_25, sdl1_27, \
                         sfk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_y[k] * sdl0_24[k]
                  + f_0 * sdk_18[k]
                  - f_14 * pc_y[k] * sdl1_24[k];

        t_70[k] = pb_y[k] * sdl0_25[k]
                  + f_16 * sdk_19[k]
                  - f_14 * pc_y[k] * sdl1_25[k];

        t_71[k] = f_15 * sdk_20[k]
                  + f_3 * pc_y[k] * sfk_56[k];

        t_72[k] = pb_y[k] * sdl0_27[k]
                  - f_14 * pc_y[k] * sdl1_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pc_x, sdk_64, sdk_65, sdk_66, sdk_67, \
                         sdk_68, sfk_64, sfk_65, sfk_66, sfk_67, \
                         sfk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_16 * sdk_64[k]
                  + f_3 * pc_x[k] * sfk_64[k];

        t_74[k] = f_16 * sdk_65[k]
                  + f_3 * pc_x[k] * sfk_65[k];

        t_75[k] = f_16 * sdk_66[k]
                  + f_3 * pc_x[k] * sfk_66[k];

        t_76[k] = f_16 * sdk_67[k]
                  + f_3 * pc_x[k] * sfk_67[k];

        t_77[k] = f_16 * sdk_68[k]
                  + f_3 * pc_x[k] * sfk_68[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pc_x, pc_y, sdk_28, sdk_69, sdk_70, sdk_71, \
                         sfi0_49, sfi1_49, sfk_64, sfk_69, sfk_70, \
                         sfk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_16 * sdk_69[k]
                  + f_3 * pc_x[k] * sfk_69[k];

        t_79[k] = f_16 * sdk_70[k]
                  + f_3 * pc_x[k] * sfk_70[k];

        t_80[k] = f_16 * sdk_71[k]
                  + f_3 * pc_x[k] * sfk_71[k];

        t_81[k] = f_15 * sdk_28[k]
                  + f_1 * sfi0_49[k]
                  - f_2 * sfi1_49[k]
                  + f_3 * pc_y[k] * sfk_64[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pc_y, pc_z, sdk_30, sdk_31, sfi0_51, sfi0_52, \
                         sfi1_51, sfi1_52, sfk_64, sfk_66, sfk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * pc_z[k] * sfk_64[k];

        t_83[k] = f_15 * sdk_30[k]
                  + f_4 * sfi0_51[k]
                  - f_5 * sfi1_51[k]
                  + f_3 * pc_y[k] * sfk_66[k];

        t_84[k] = f_15 * sdk_31[k]
                  + f_6 * sfi0_52[k]
                  - f_7 * sfi1_52[k]
                  + f_3 * pc_y[k] * sfk_67[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_y, sdk_32, sdk_33, sdk_34, sfi0_53, sfi0_54, \
                         sfi0_55, sfi1_53, sfi1_54, sfi1_55, sfk_68, sfk_69, \
                         sfk_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_15 * sdk_32[k]
                  + f_8 * sfi0_53[k]
                  - f_9 * sfi1_53[k]
                  + f_3 * pc_y[k] * sfk_68[k];

        t_86[k] = f_15 * sdk_33[k]
                  + f_10 * sfi0_54[k]
                  - f_11 * sfi1_54[k]
                  + f_3 * pc_y[k] * sfk_69[k];

        t_87[k] = f_15 * sdk_34[k]
                  + f_12 * sfi0_55[k]
                  - f_13 * sfi1_55[k]
                  + f_3 * pc_y[k] * sfk_70[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pb_y, pb_z, pc_y, pc_z, sdl0_0, sdl0_44, \
                         sdk_35, sdl1_0, sdl1_44, sfk_71, sfk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * sdk_35[k]
                  + f_3 * pc_y[k] * sfk_71[k];

        t_89[k] = pb_y[k] * sdl0_44[k]
                  - f_14 * pc_y[k] * sdl1_44[k];

        t_90[k] = pb_z[k] * sdl0_0[k]
                  - f_14 * pc_z[k] * sdl1_0[k];

        t_91[k] = f_3 * pc_y[k] * sfk_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_z, pc_y, pc_z, sdl0_3, sdl0_5, sdk_0, \
                         sdk_2, sdl1_3, sdl1_5, sfk_72, sfk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * sdk_0[k]
                  + f_3 * pc_z[k] * sfk_72[k];

        t_93[k] = pb_z[k] * sdl0_3[k]
                  - f_14 * pc_z[k] * sdl1_3[k];

        t_94[k] = f_3 * pc_y[k] * sfk_74[k];

        t_95[k] = pb_z[k] * sdl0_5[k]
                  + f_16 * sdk_2[k]
                  - f_14 * pc_z[k] * sdl1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_z, pc_y, pc_z, sdl0_6, sdl0_9, sdk_3, \
                         sdk_5, sdl1_6, sdl1_9, sfk_75, sfk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pb_z[k] * sdl0_6[k]
                  - f_14 * pc_z[k] * sdl1_6[k];

        t_97[k] = f_15 * sdk_3[k]
                  + f_3 * pc_z[k] * sfk_75[k];

        t_98[k] = f_3 * pc_y[k] * sfk_77[k];

        t_99[k] = pb_z[k] * sdl0_9[k]
                  + f_0 * sdk_5[k]
                  - f_14 * pc_z[k] * sdl1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pb_z, pc_y, pc_z, sdl0_10, sdl0_12, \
                         sdk_6, sdk_7, sdl1_10, sdl1_12, sfk_78, \
                         sfk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pb_z[k] * sdl0_10[k]
                   - f_14 * pc_z[k] * sdl1_10[k];

        t_101[k] = f_15 * sdk_6[k]
                   + f_3 * pc_z[k] * sfk_78[k];

        t_102[k] = pb_z[k] * sdl0_12[k]
                   + f_16 * sdk_7[k]
                   - f_14 * pc_z[k] * sdl1_12[k];

        t_103[k] = f_3 * pc_y[k] * sfk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_z, pc_z, sdl0_14, sdl0_15, sdl0_17, \
                         sdk_9, sdk_10, sdk_11, sdl1_14, sdl1_15, sdl1_17, \
                         sfk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_z[k] * sdl0_14[k]
                   + f_17 * sdk_9[k]
                   - f_14 * pc_z[k] * sdl1_14[k];

        t_105[k] = pb_z[k] * sdl0_15[k]
                   - f_14 * pc_z[k] * sdl1_15[k];

        t_106[k] = f_15 * sdk_10[k]
                   + f_3 * pc_z[k] * sfk_82[k];

        t_107[k] = pb_z[k] * sdl0_17[k]
                   + f_16 * sdk_11[k]
                   - f_14 * pc_z[k] * sdl1_17[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pb_z, pc_y, pc_z, sdl0_18, sdl0_20, \
                         sdl0_21, sdk_12, sdk_14, sdl1_18, sdl1_20, sdl1_21, \
                         sfk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_z[k] * sdl0_18[k]
                   + f_0 * sdk_12[k]
                   - f_14 * pc_z[k] * sdl1_18[k];

        t_109[k] = f_3 * pc_y[k] * sfk_86[k];

        t_110[k] = pb_z[k] * sdl0_20[k]
                   + f_18 * sdk_14[k]
                   - f_14 * pc_z[k] * sdl1_20[k];

        t_111[k] = pb_z[k] * sdl0_21[k]
                   - f_14 * pc_z[k] * sdl1_21[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_z, pc_z, sdl0_23, sdl0_24, sdk_15, sdk_16, \
                         sdk_17, sdl1_23, sdl1_24, sfk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_15 * sdk_15[k]
                   + f_3 * pc_z[k] * sfk_87[k];

        t_113[k] = pb_z[k] * sdl0_23[k]
                   + f_16 * sdk_16[k]
                   - f_14 * pc_z[k] * sdl1_23[k];

        t_114[k] = pb_z[k] * sdl0_24[k]
                   + f_0 * sdk_17[k]
                   - f_14 * pc_z[k] * sdl1_24[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pb_z, pc_y, pc_z, sdl0_25, sdl0_27, sdk_18, \
                         sdk_20, sdl1_25, sdl1_27, sfk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_z[k] * sdl0_25[k]
                   + f_17 * sdk_18[k]
                   - f_14 * pc_z[k] * sdl1_25[k];

        t_116[k] = f_3 * pc_y[k] * sfk_92[k];

        t_117[k] = pb_z[k] * sdl0_27[k]
                   + f_19 * sdk_20[k]
                   - f_14 * pc_z[k] * sdl1_27[k];
    }
}

static auto
compute_prim_sfl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sdl0,
                                                          const size_t sdk, const size_t sdl1,
                                                          const size_t sfi0, const size_t sfi1,
                                                          const size_t sfk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_17 = 2.0 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 4.0 / q;

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
    auto *t_233 = buffer.data(target + 233);
    auto *t_234 = buffer.data(target + 234);
    auto *t_235 = buffer.data(target + 235);
    auto *t_236 = buffer.data(target + 236);
    auto *t_237 = buffer.data(target + 237);
    auto *t_238 = buffer.data(target + 238);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdl0_36 = buffer.data(sdl0 + 36);
    const auto *sdl0_48 = buffer.data(sdl0 + 48);
    const auto *sdl0_51 = buffer.data(sdl0 + 51);
    const auto *sdl0_55 = buffer.data(sdl0 + 55);
    const auto *sdl0_60 = buffer.data(sdl0 + 60);
    const auto *sdl0_66 = buffer.data(sdl0 + 66);
    const auto *sdl0_90 = buffer.data(sdl0 + 90);
    const auto *sdl0_95 = buffer.data(sdl0 + 95);
    const auto *sdl0_99 = buffer.data(sdl0 + 99);
    const auto *sdl0_104 = buffer.data(sdl0 + 104);
    const auto *sdl0_110 = buffer.data(sdl0 + 110);
    const auto *sdl0_117 = buffer.data(sdl0 + 117);
    const auto *sdl0_135 = buffer.data(sdl0 + 135);
    const auto *sdl0_138 = buffer.data(sdl0 + 138);
    const auto *sdl0_140 = buffer.data(sdl0 + 140);
    const auto *sdl0_141 = buffer.data(sdl0 + 141);
    const auto *sdl0_144 = buffer.data(sdl0 + 144);
    const auto *sdl0_145 = buffer.data(sdl0 + 145);
    const auto *sdl0_147 = buffer.data(sdl0 + 147);
    const auto *sdl0_149 = buffer.data(sdl0 + 149);
    const auto *sdl0_150 = buffer.data(sdl0 + 150);
    const auto *sdl0_152 = buffer.data(sdl0 + 152);
    const auto *sdl0_153 = buffer.data(sdl0 + 153);
    const auto *sdl0_155 = buffer.data(sdl0 + 155);
    const auto *sdl0_156 = buffer.data(sdl0 + 156);
    const auto *sdl0_158 = buffer.data(sdl0 + 158);
    const auto *sdl0_159 = buffer.data(sdl0 + 159);
    const auto *sdl0_160 = buffer.data(sdl0 + 160);
    const auto *sdl0_162 = buffer.data(sdl0 + 162);
    const auto *sdl0_171 = buffer.data(sdl0 + 171);
    const auto *sdl0_173 = buffer.data(sdl0 + 173);
    const auto *sdl0_174 = buffer.data(sdl0 + 174);
    const auto *sdl0_175 = buffer.data(sdl0 + 175);
    const auto *sdl0_176 = buffer.data(sdl0 + 176);
    const auto *sdl0_177 = buffer.data(sdl0 + 177);
    const auto *sdl0_179 = buffer.data(sdl0 + 179);
    const auto *sdl0_192 = buffer.data(sdl0 + 192);
    const auto *sdl0_197 = buffer.data(sdl0 + 197);
    const auto *sdl0_198 = buffer.data(sdl0 + 198);
    const auto *sdl0_203 = buffer.data(sdl0 + 203);
    const auto *sdl0_204 = buffer.data(sdl0 + 204);
    const auto *sdl0_205 = buffer.data(sdl0 + 205);
    const auto *sdl0_216 = buffer.data(sdl0 + 216);
    const auto *sdl0_218 = buffer.data(sdl0 + 218);
    const auto *sdl0_219 = buffer.data(sdl0 + 219);
    const auto *sdl0_220 = buffer.data(sdl0 + 220);
    const auto *sdl0_221 = buffer.data(sdl0 + 221);
    const auto *sdl0_222 = buffer.data(sdl0 + 222);
    const auto *sdl0_224 = buffer.data(sdl0 + 224);
    const auto *sdl0_225 = buffer.data(sdl0 + 225);
    const auto *sdl0_228 = buffer.data(sdl0 + 228);
    const auto *sdl0_230 = buffer.data(sdl0 + 230);
    const auto *sdl0_231 = buffer.data(sdl0 + 231);
    const auto *sdl0_234 = buffer.data(sdl0 + 234);
    const auto *sdl0_235 = buffer.data(sdl0 + 235);
    const auto *sdl0_237 = buffer.data(sdl0 + 237);

    const auto *sdk_28 = buffer.data(sdk + 28);
    const auto *sdk_35 = buffer.data(sdk + 35);
    const auto *sdk_36 = buffer.data(sdk + 36);
    const auto *sdk_38 = buffer.data(sdk + 38);
    const auto *sdk_39 = buffer.data(sdk + 39);
    const auto *sdk_41 = buffer.data(sdk + 41);
    const auto *sdk_42 = buffer.data(sdk + 42);
    const auto *sdk_45 = buffer.data(sdk + 45);
    const auto *sdk_46 = buffer.data(sdk + 46);
    const auto *sdk_50 = buffer.data(sdk + 50);
    const auto *sdk_51 = buffer.data(sdk + 51);
    const auto *sdk_56 = buffer.data(sdk + 56);
    const auto *sdk_64 = buffer.data(sdk + 64);
    const auto *sdk_71 = buffer.data(sdk + 71);
    const auto *sdk_72 = buffer.data(sdk + 72);
    const auto *sdk_74 = buffer.data(sdk + 74);
    const auto *sdk_75 = buffer.data(sdk + 75);
    const auto *sdk_77 = buffer.data(sdk + 77);
    const auto *sdk_78 = buffer.data(sdk + 78);
    const auto *sdk_81 = buffer.data(sdk + 81);
    const auto *sdk_86 = buffer.data(sdk + 86);
    const auto *sdk_92 = buffer.data(sdk + 92);
    const auto *sdk_100 = buffer.data(sdk + 100);
    const auto *sdk_101 = buffer.data(sdk + 101);
    const auto *sdk_102 = buffer.data(sdk + 102);
    const auto *sdk_103 = buffer.data(sdk + 103);
    const auto *sdk_104 = buffer.data(sdk + 104);
    const auto *sdk_105 = buffer.data(sdk + 105);
    const auto *sdk_106 = buffer.data(sdk + 106);
    const auto *sdk_107 = buffer.data(sdk + 107);
    const auto *sdk_108 = buffer.data(sdk + 108);
    const auto *sdk_111 = buffer.data(sdk + 111);
    const auto *sdk_113 = buffer.data(sdk + 113);
    const auto *sdk_114 = buffer.data(sdk + 114);
    const auto *sdk_117 = buffer.data(sdk + 117);
    const auto *sdk_118 = buffer.data(sdk + 118);
    const auto *sdk_120 = buffer.data(sdk + 120);
    const auto *sdk_122 = buffer.data(sdk + 122);
    const auto *sdk_123 = buffer.data(sdk + 123);
    const auto *sdk_125 = buffer.data(sdk + 125);
    const auto *sdk_126 = buffer.data(sdk + 126);
    const auto *sdk_128 = buffer.data(sdk + 128);
    const auto *sdk_129 = buffer.data(sdk + 129);
    const auto *sdk_131 = buffer.data(sdk + 131);
    const auto *sdk_132 = buffer.data(sdk + 132);
    const auto *sdk_133 = buffer.data(sdk + 133);
    const auto *sdk_135 = buffer.data(sdk + 135);
    const auto *sdk_136 = buffer.data(sdk + 136);
    const auto *sdk_137 = buffer.data(sdk + 137);
    const auto *sdk_138 = buffer.data(sdk + 138);
    const auto *sdk_139 = buffer.data(sdk + 139);
    const auto *sdk_140 = buffer.data(sdk + 140);
    const auto *sdk_141 = buffer.data(sdk + 141);
    const auto *sdk_142 = buffer.data(sdk + 142);
    const auto *sdk_143 = buffer.data(sdk + 143);
    const auto *sdk_156 = buffer.data(sdk + 156);
    const auto *sdk_161 = buffer.data(sdk + 161);
    const auto *sdk_162 = buffer.data(sdk + 162);
    const auto *sdk_167 = buffer.data(sdk + 167);
    const auto *sdk_168 = buffer.data(sdk + 168);
    const auto *sdk_169 = buffer.data(sdk + 169);
    const auto *sdk_172 = buffer.data(sdk + 172);
    const auto *sdk_173 = buffer.data(sdk + 173);
    const auto *sdk_174 = buffer.data(sdk + 174);
    const auto *sdk_175 = buffer.data(sdk + 175);
    const auto *sdk_176 = buffer.data(sdk + 176);
    const auto *sdk_177 = buffer.data(sdk + 177);
    const auto *sdk_178 = buffer.data(sdk + 178);
    const auto *sdk_179 = buffer.data(sdk + 179);
    const auto *sdk_180 = buffer.data(sdk + 180);
    const auto *sdk_183 = buffer.data(sdk + 183);
    const auto *sdk_185 = buffer.data(sdk + 185);
    const auto *sdk_186 = buffer.data(sdk + 186);
    const auto *sdk_189 = buffer.data(sdk + 189);
    const auto *sdk_190 = buffer.data(sdk + 190);
    const auto *sdk_192 = buffer.data(sdk + 192);

    const auto *sdl1_36 = buffer.data(sdl1 + 36);
    const auto *sdl1_48 = buffer.data(sdl1 + 48);
    const auto *sdl1_51 = buffer.data(sdl1 + 51);
    const auto *sdl1_55 = buffer.data(sdl1 + 55);
    const auto *sdl1_60 = buffer.data(sdl1 + 60);
    const auto *sdl1_66 = buffer.data(sdl1 + 66);
    const auto *sdl1_90 = buffer.data(sdl1 + 90);
    const auto *sdl1_95 = buffer.data(sdl1 + 95);
    const auto *sdl1_99 = buffer.data(sdl1 + 99);
    const auto *sdl1_104 = buffer.data(sdl1 + 104);
    const auto *sdl1_110 = buffer.data(sdl1 + 110);
    const auto *sdl1_117 = buffer.data(sdl1 + 117);
    const auto *sdl1_135 = buffer.data(sdl1 + 135);
    const auto *sdl1_138 = buffer.data(sdl1 + 138);
    const auto *sdl1_140 = buffer.data(sdl1 + 140);
    const auto *sdl1_141 = buffer.data(sdl1 + 141);
    const auto *sdl1_144 = buffer.data(sdl1 + 144);
    const auto *sdl1_145 = buffer.data(sdl1 + 145);
    const auto *sdl1_147 = buffer.data(sdl1 + 147);
    const auto *sdl1_149 = buffer.data(sdl1 + 149);
    const auto *sdl1_150 = buffer.data(sdl1 + 150);
    const auto *sdl1_152 = buffer.data(sdl1 + 152);
    const auto *sdl1_153 = buffer.data(sdl1 + 153);
    const auto *sdl1_155 = buffer.data(sdl1 + 155);
    const auto *sdl1_156 = buffer.data(sdl1 + 156);
    const auto *sdl1_158 = buffer.data(sdl1 + 158);
    const auto *sdl1_159 = buffer.data(sdl1 + 159);
    const auto *sdl1_160 = buffer.data(sdl1 + 160);
    const auto *sdl1_162 = buffer.data(sdl1 + 162);
    const auto *sdl1_171 = buffer.data(sdl1 + 171);
    const auto *sdl1_173 = buffer.data(sdl1 + 173);
    const auto *sdl1_174 = buffer.data(sdl1 + 174);
    const auto *sdl1_175 = buffer.data(sdl1 + 175);
    const auto *sdl1_176 = buffer.data(sdl1 + 176);
    const auto *sdl1_177 = buffer.data(sdl1 + 177);
    const auto *sdl1_179 = buffer.data(sdl1 + 179);
    const auto *sdl1_192 = buffer.data(sdl1 + 192);
    const auto *sdl1_197 = buffer.data(sdl1 + 197);
    const auto *sdl1_198 = buffer.data(sdl1 + 198);
    const auto *sdl1_203 = buffer.data(sdl1 + 203);
    const auto *sdl1_204 = buffer.data(sdl1 + 204);
    const auto *sdl1_205 = buffer.data(sdl1 + 205);
    const auto *sdl1_216 = buffer.data(sdl1 + 216);
    const auto *sdl1_218 = buffer.data(sdl1 + 218);
    const auto *sdl1_219 = buffer.data(sdl1 + 219);
    const auto *sdl1_220 = buffer.data(sdl1 + 220);
    const auto *sdl1_221 = buffer.data(sdl1 + 221);
    const auto *sdl1_222 = buffer.data(sdl1 + 222);
    const auto *sdl1_224 = buffer.data(sdl1 + 224);
    const auto *sdl1_225 = buffer.data(sdl1 + 225);
    const auto *sdl1_228 = buffer.data(sdl1 + 228);
    const auto *sdl1_230 = buffer.data(sdl1 + 230);
    const auto *sdl1_231 = buffer.data(sdl1 + 231);
    const auto *sdl1_234 = buffer.data(sdl1 + 234);
    const auto *sdl1_235 = buffer.data(sdl1 + 235);
    const auto *sdl1_237 = buffer.data(sdl1 + 237);

    const auto *sfi0_79 = buffer.data(sfi0 + 79);
    const auto *sfi0_80 = buffer.data(sfi0 + 80);
    const auto *sfi0_81 = buffer.data(sfi0 + 81);
    const auto *sfi0_82 = buffer.data(sfi0 + 82);
    const auto *sfi0_83 = buffer.data(sfi0 + 83);

    const auto *sfi1_79 = buffer.data(sfi1 + 79);
    const auto *sfi1_80 = buffer.data(sfi1 + 80);
    const auto *sfi1_81 = buffer.data(sfi1 + 81);
    const auto *sfi1_82 = buffer.data(sfi1 + 82);
    const auto *sfi1_83 = buffer.data(sfi1 + 83);

    const auto *sfk_100 = buffer.data(sfk + 100);
    const auto *sfk_101 = buffer.data(sfk + 101);
    const auto *sfk_102 = buffer.data(sfk + 102);
    const auto *sfk_103 = buffer.data(sfk + 103);
    const auto *sfk_104 = buffer.data(sfk + 104);
    const auto *sfk_105 = buffer.data(sfk + 105);
    const auto *sfk_106 = buffer.data(sfk + 106);
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
    const auto *sfk_136 = buffer.data(sfk + 136);
    const auto *sfk_137 = buffer.data(sfk + 137);
    const auto *sfk_138 = buffer.data(sfk + 138);
    const auto *sfk_139 = buffer.data(sfk + 139);
    const auto *sfk_140 = buffer.data(sfk + 140);
    const auto *sfk_141 = buffer.data(sfk + 141);
    const auto *sfk_142 = buffer.data(sfk + 142);
    const auto *sfk_143 = buffer.data(sfk + 143);
    const auto *sfk_144 = buffer.data(sfk + 144);
    const auto *sfk_146 = buffer.data(sfk + 146);
    const auto *sfk_147 = buffer.data(sfk + 147);
    const auto *sfk_149 = buffer.data(sfk + 149);
    const auto *sfk_150 = buffer.data(sfk + 150);
    const auto *sfk_153 = buffer.data(sfk + 153);
    const auto *sfk_154 = buffer.data(sfk + 154);
    const auto *sfk_158 = buffer.data(sfk + 158);
    const auto *sfk_159 = buffer.data(sfk + 159);
    const auto *sfk_164 = buffer.data(sfk + 164);
    const auto *sfk_172 = buffer.data(sfk + 172);
    const auto *sfk_173 = buffer.data(sfk + 173);
    const auto *sfk_174 = buffer.data(sfk + 174);
    const auto *sfk_175 = buffer.data(sfk + 175);
    const auto *sfk_176 = buffer.data(sfk + 176);
    const auto *sfk_177 = buffer.data(sfk + 177);
    const auto *sfk_178 = buffer.data(sfk + 178);
    const auto *sfk_179 = buffer.data(sfk + 179);
    const auto *sfk_180 = buffer.data(sfk + 180);
    const auto *sfk_182 = buffer.data(sfk + 182);
    const auto *sfk_183 = buffer.data(sfk + 183);
    const auto *sfk_185 = buffer.data(sfk + 185);
    const auto *sfk_186 = buffer.data(sfk + 186);
    const auto *sfk_189 = buffer.data(sfk + 189);

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pc_x, sdk_100, sdk_101, sdk_102, \
                         sdk_103, sdk_104, sfk_100, sfk_101, sfk_102, sfk_103, \
                         sfk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_16 * sdk_100[k]
                   + f_3 * pc_x[k] * sfk_100[k];

        t_119[k] = f_16 * sdk_101[k]
                   + f_3 * pc_x[k] * sfk_101[k];

        t_120[k] = f_16 * sdk_102[k]
                   + f_3 * pc_x[k] * sfk_102[k];

        t_121[k] = f_16 * sdk_103[k]
                   + f_3 * pc_x[k] * sfk_103[k];

        t_122[k] = f_16 * sdk_104[k]
                   + f_3 * pc_x[k] * sfk_104[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pb_z, pc_x, pc_z, sdl0_36, sdk_105, \
                         sdk_106, sdk_107, sdl1_36, sfk_105, sfk_106, \
                         sfk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_16 * sdk_105[k]
                   + f_3 * pc_x[k] * sfk_105[k];

        t_124[k] = f_16 * sdk_106[k]
                   + f_3 * pc_x[k] * sfk_106[k];

        t_125[k] = f_16 * sdk_107[k]
                   + f_3 * pc_x[k] * sfk_107[k];

        t_126[k] = pb_z[k] * sdl0_36[k]
                   - f_14 * pc_z[k] * sdl1_36[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pc_y, pc_z, sdk_28, sfi0_79, sfi0_80, sfi1_79, \
                         sfi1_80, sfk_100, sfk_102, sfk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_15 * sdk_28[k]
                   + f_3 * pc_z[k] * sfk_100[k];

        t_128[k] = f_4 * sfi0_79[k]
                   - f_5 * sfi1_79[k]
                   + f_3 * pc_y[k] * sfk_102[k];

        t_129[k] = f_6 * sfi0_80[k]
                   - f_7 * sfi1_80[k]
                   + f_3 * pc_y[k] * sfk_103[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_y, sfi0_81, sfi0_82, sfi0_83, sfi1_81, \
                         sfi1_82, sfi1_83, sfk_104, sfk_105, sfk_106, \
                         sfk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_8 * sfi0_81[k]
                   - f_9 * sfi1_81[k]
                   + f_3 * pc_y[k] * sfk_104[k];

        t_131[k] = f_10 * sfi0_82[k]
                   - f_11 * sfi1_82[k]
                   + f_3 * pc_y[k] * sfk_105[k];

        t_132[k] = f_12 * sfi0_83[k]
                   - f_13 * sfi1_83[k]
                   + f_3 * pc_y[k] * sfk_106[k];

        t_133[k] = f_3 * pc_y[k] * sfk_107[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pb_x, pc_x, pc_y, pc_z, sdl0_135, sdk_35, \
                         sdk_36, sdk_108, sdl1_135, sfi0_83, sfi1_83, sfk_107, \
                         sfk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_15 * sdk_35[k]
                   + f_1 * sfi0_83[k]
                   - f_2 * sfi1_83[k]
                   + f_3 * pc_z[k] * sfk_107[k];

        t_135[k] = pb_x[k] * sdl0_135[k]
                   + f_20 * sdk_108[k]
                   - f_14 * pc_x[k] * sdl1_135[k];

        t_136[k] = f_16 * sdk_36[k]
                   + f_3 * pc_y[k] * sfk_108[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pb_x, pc_x, pc_y, pc_z, sdl0_138, sdk_38, \
                         sdk_111, sdl1_138, sfk_108, sfk_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_3 * pc_z[k] * sfk_108[k];

        t_138[k] = pb_x[k] * sdl0_138[k]
                   + f_19 * sdk_111[k]
                   - f_14 * pc_x[k] * sdl1_138[k];

        t_139[k] = f_16 * sdk_38[k]
                   + f_3 * pc_y[k] * sfk_110[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pb_x, pc_x, pc_z, sdl0_140, sdl0_141, sdk_113, \
                         sdk_114, sdl1_140, sdl1_141, sfk_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = pb_x[k] * sdl0_140[k]
                   + f_19 * sdk_113[k]
                   - f_14 * pc_x[k] * sdl1_140[k];

        t_141[k] = pb_x[k] * sdl0_141[k]
                   + f_18 * sdk_114[k]
                   - f_14 * pc_x[k] * sdl1_141[k];

        t_142[k] = f_3 * pc_z[k] * sfk_111[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pb_x, pc_x, pc_y, sdl0_144, sdl0_145, sdk_41, \
                         sdk_117, sdk_118, sdl1_144, sdl1_145, \
                         sfk_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_16 * sdk_41[k]
                   + f_3 * pc_y[k] * sfk_113[k];

        t_144[k] = pb_x[k] * sdl0_144[k]
                   + f_18 * sdk_117[k]
                   - f_14 * pc_x[k] * sdl1_144[k];

        t_145[k] = pb_x[k] * sdl0_145[k]
                   + f_17 * sdk_118[k]
                   - f_14 * pc_x[k] * sdl1_145[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pb_x, pc_x, pc_y, pc_z, sdl0_147, sdk_45, \
                         sdk_120, sdl1_147, sfk_114, sfk_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_3 * pc_z[k] * sfk_114[k];

        t_147[k] = pb_x[k] * sdl0_147[k]
                   + f_17 * sdk_120[k]
                   - f_14 * pc_x[k] * sdl1_147[k];

        t_148[k] = f_16 * sdk_45[k]
                   + f_3 * pc_y[k] * sfk_117[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pb_x, pc_x, pc_z, sdl0_149, sdl0_150, sdk_122, \
                         sdk_123, sdl1_149, sdl1_150, sfk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pb_x[k] * sdl0_149[k]
                   + f_17 * sdk_122[k]
                   - f_14 * pc_x[k] * sdl1_149[k];

        t_150[k] = pb_x[k] * sdl0_150[k]
                   + f_0 * sdk_123[k]
                   - f_14 * pc_x[k] * sdl1_150[k];

        t_151[k] = f_3 * pc_z[k] * sfk_118[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pb_x, pc_x, pc_y, sdl0_152, sdl0_153, sdk_50, \
                         sdk_125, sdk_126, sdl1_152, sdl1_153, \
                         sfk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = pb_x[k] * sdl0_152[k]
                   + f_0 * sdk_125[k]
                   - f_14 * pc_x[k] * sdl1_152[k];

        t_153[k] = pb_x[k] * sdl0_153[k]
                   + f_0 * sdk_126[k]
                   - f_14 * pc_x[k] * sdl1_153[k];

        t_154[k] = f_16 * sdk_50[k]
                   + f_3 * pc_y[k] * sfk_122[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, pb_x, pc_x, pc_z, sdl0_155, sdl0_156, sdk_128, \
                         sdk_129, sdl1_155, sdl1_156, sfk_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = pb_x[k] * sdl0_155[k]
                   + f_0 * sdk_128[k]
                   - f_14 * pc_x[k] * sdl1_155[k];

        t_156[k] = pb_x[k] * sdl0_156[k]
                   + f_16 * sdk_129[k]
                   - f_14 * pc_x[k] * sdl1_156[k];

        t_157[k] = f_3 * pc_z[k] * sfk_123[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pb_x, pc_x, sdl0_158, sdl0_159, sdl0_160, \
                         sdk_131, sdk_132, sdk_133, sdl1_158, sdl1_159, \
                         sdl1_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = pb_x[k] * sdl0_158[k]
                   + f_16 * sdk_131[k]
                   - f_14 * pc_x[k] * sdl1_158[k];

        t_159[k] = pb_x[k] * sdl0_159[k]
                   + f_16 * sdk_132[k]
                   - f_14 * pc_x[k] * sdl1_159[k];

        t_160[k] = pb_x[k] * sdl0_160[k]
                   + f_16 * sdk_133[k]
                   - f_14 * pc_x[k] * sdl1_160[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_x, pc_x, pc_y, sdl0_162, sdk_56, \
                         sdk_135, sdk_136, sdk_137, sdl1_162, sfk_128, sfk_136, \
                         sfk_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_16 * sdk_56[k]
                   + f_3 * pc_y[k] * sfk_128[k];

        t_162[k] = pb_x[k] * sdl0_162[k]
                   + f_16 * sdk_135[k]
                   - f_14 * pc_x[k] * sdl1_162[k];

        t_163[k] = f_15 * sdk_136[k]
                   + f_3 * pc_x[k] * sfk_136[k];

        t_164[k] = f_15 * sdk_137[k]
                   + f_3 * pc_x[k] * sfk_137[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pc_x, sdk_138, sdk_139, sdk_140, \
                         sdk_141, sdk_142, sfk_138, sfk_139, sfk_140, sfk_141, \
                         sfk_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_15 * sdk_138[k]
                   + f_3 * pc_x[k] * sfk_138[k];

        t_166[k] = f_15 * sdk_139[k]
                   + f_3 * pc_x[k] * sfk_139[k];

        t_167[k] = f_15 * sdk_140[k]
                   + f_3 * pc_x[k] * sfk_140[k];

        t_168[k] = f_15 * sdk_141[k]
                   + f_3 * pc_x[k] * sfk_141[k];

        t_169[k] = f_15 * sdk_142[k]
                   + f_3 * pc_x[k] * sfk_142[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_x, pc_x, pc_z, sdl0_171, sdl0_173, \
                         sdk_143, sdl1_171, sdl1_173, sfk_136, \
                         sfk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_15 * sdk_143[k]
                   + f_3 * pc_x[k] * sfk_143[k];

        t_171[k] = pb_x[k] * sdl0_171[k]
                   - f_14 * pc_x[k] * sdl1_171[k];

        t_172[k] = f_3 * pc_z[k] * sfk_136[k];

        t_173[k] = pb_x[k] * sdl0_173[k]
                   - f_14 * pc_x[k] * sdl1_173[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pb_x, pc_x, sdl0_174, sdl0_175, sdl0_176, \
                         sdl0_177, sdl1_174, sdl1_175, sdl1_176, \
                         sdl1_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pb_x[k] * sdl0_174[k]
                   - f_14 * pc_x[k] * sdl1_174[k];

        t_175[k] = pb_x[k] * sdl0_175[k]
                   - f_14 * pc_x[k] * sdl1_175[k];

        t_176[k] = pb_x[k] * sdl0_176[k]
                   - f_14 * pc_x[k] * sdl1_176[k];

        t_177[k] = pb_x[k] * sdl0_177[k]
                   - f_14 * pc_x[k] * sdl1_177[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, pb_x, pb_y, pc_x, pc_y, sdl0_90, \
                         sdl0_179, sdk_71, sdk_72, sdl1_90, sdl1_179, sfk_143, \
                         sfk_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_16 * sdk_71[k]
                   + f_3 * pc_y[k] * sfk_143[k];

        t_179[k] = pb_x[k] * sdl0_179[k]
                   - f_14 * pc_x[k] * sdl1_179[k];

        t_180[k] = pb_y[k] * sdl0_90[k]
                   - f_14 * pc_y[k] * sdl1_90[k];

        t_181[k] = f_15 * sdk_72[k]
                   + f_3 * pc_y[k] * sfk_144[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pb_y, pb_z, pc_y, pc_z, sdl0_48, sdl0_95, \
                         sdk_36, sdk_74, sdl1_48, sdl1_95, sfk_144, \
                         sfk_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_15 * sdk_36[k]
                   + f_3 * pc_z[k] * sfk_144[k];

        t_183[k] = pb_z[k] * sdl0_48[k]
                   - f_14 * pc_z[k] * sdl1_48[k];

        t_184[k] = f_15 * sdk_74[k]
                   + f_3 * pc_y[k] * sfk_146[k];

        t_185[k] = pb_y[k] * sdl0_95[k]
                   - f_14 * pc_y[k] * sdl1_95[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pb_y, pb_z, pc_y, pc_z, sdl0_51, sdl0_99, \
                         sdk_39, sdk_77, sdl1_51, sdl1_99, sfk_147, \
                         sfk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = pb_z[k] * sdl0_51[k]
                   - f_14 * pc_z[k] * sdl1_51[k];

        t_187[k] = f_15 * sdk_39[k]
                   + f_3 * pc_z[k] * sfk_147[k];

        t_188[k] = f_15 * sdk_77[k]
                   + f_3 * pc_y[k] * sfk_149[k];

        t_189[k] = pb_y[k] * sdl0_99[k]
                   - f_14 * pc_y[k] * sdl1_99[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pb_x, pb_z, pc_x, pc_z, sdl0_55, sdl0_192, \
                         sdk_42, sdk_156, sdl1_55, sdl1_192, sfk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = pb_z[k] * sdl0_55[k]
                   - f_14 * pc_z[k] * sdl1_55[k];

        t_191[k] = f_15 * sdk_42[k]
                   + f_3 * pc_z[k] * sfk_150[k];

        t_192[k] = pb_x[k] * sdl0_192[k]
                   + f_17 * sdk_156[k]
                   - f_14 * pc_x[k] * sdl1_192[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pb_y, pb_z, pc_y, pc_z, sdl0_60, \
                         sdl0_104, sdk_46, sdk_81, sdl1_60, sdl1_104, sfk_153, \
                         sfk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_15 * sdk_81[k]
                   + f_3 * pc_y[k] * sfk_153[k];

        t_194[k] = pb_y[k] * sdl0_104[k]
                   - f_14 * pc_y[k] * sdl1_104[k];

        t_195[k] = pb_z[k] * sdl0_60[k]
                   - f_14 * pc_z[k] * sdl1_60[k];

        t_196[k] = f_15 * sdk_46[k]
                   + f_3 * pc_z[k] * sfk_154[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pb_x, pc_x, pc_y, sdl0_197, sdl0_198, sdk_86, \
                         sdk_161, sdk_162, sdl1_197, sdl1_198, \
                         sfk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = pb_x[k] * sdl0_197[k]
                   + f_0 * sdk_161[k]
                   - f_14 * pc_x[k] * sdl1_197[k];

        t_198[k] = pb_x[k] * sdl0_198[k]
                   + f_0 * sdk_162[k]
                   - f_14 * pc_x[k] * sdl1_198[k];

        t_199[k] = f_15 * sdk_86[k]
                   + f_3 * pc_y[k] * sfk_158[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_y, pb_z, pc_y, pc_z, sdl0_66, sdl0_110, \
                         sdk_51, sdl1_66, sdl1_110, sfk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pb_y[k] * sdl0_110[k]
                   - f_14 * pc_y[k] * sdl1_110[k];

        t_201[k] = pb_z[k] * sdl0_66[k]
                   - f_14 * pc_z[k] * sdl1_66[k];

        t_202[k] = f_15 * sdk_51[k]
                   + f_3 * pc_z[k] * sfk_159[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pb_x, pc_x, sdl0_203, sdl0_204, sdl0_205, \
                         sdk_167, sdk_168, sdk_169, sdl1_203, sdl1_204, \
                         sdl1_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = pb_x[k] * sdl0_203[k]
                   + f_16 * sdk_167[k]
                   - f_14 * pc_x[k] * sdl1_203[k];

        t_204[k] = pb_x[k] * sdl0_204[k]
                   + f_16 * sdk_168[k]
                   - f_14 * pc_x[k] * sdl1_204[k];

        t_205[k] = pb_x[k] * sdl0_205[k]
                   + f_16 * sdk_169[k]
                   - f_14 * pc_x[k] * sdl1_205[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pb_y, pc_x, pc_y, sdl0_117, sdk_92, \
                         sdk_172, sdk_173, sdl1_117, sfk_164, sfk_172, \
                         sfk_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_15 * sdk_92[k]
                   + f_3 * pc_y[k] * sfk_164[k];

        t_207[k] = pb_y[k] * sdl0_117[k]
                   - f_14 * pc_y[k] * sdl1_117[k];

        t_208[k] = f_15 * sdk_172[k]
                   + f_3 * pc_x[k] * sfk_172[k];

        t_209[k] = f_15 * sdk_173[k]
                   + f_3 * pc_x[k] * sfk_173[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pc_x, sdk_174, sdk_175, sdk_176, \
                         sdk_177, sdk_178, sfk_174, sfk_175, sfk_176, sfk_177, \
                         sfk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_15 * sdk_174[k]
                   + f_3 * pc_x[k] * sfk_174[k];

        t_211[k] = f_15 * sdk_175[k]
                   + f_3 * pc_x[k] * sfk_175[k];

        t_212[k] = f_15 * sdk_176[k]
                   + f_3 * pc_x[k] * sfk_176[k];

        t_213[k] = f_15 * sdk_177[k]
                   + f_3 * pc_x[k] * sfk_177[k];

        t_214[k] = f_15 * sdk_178[k]
                   + f_3 * pc_x[k] * sfk_178[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pb_x, pc_x, pc_z, sdl0_216, sdl0_218, \
                         sdk_64, sdk_179, sdl1_216, sdl1_218, sfk_172, \
                         sfk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_15 * sdk_179[k]
                   + f_3 * pc_x[k] * sfk_179[k];

        t_216[k] = pb_x[k] * sdl0_216[k]
                   - f_14 * pc_x[k] * sdl1_216[k];

        t_217[k] = f_15 * sdk_64[k]
                   + f_3 * pc_z[k] * sfk_172[k];

        t_218[k] = pb_x[k] * sdl0_218[k]
                   - f_14 * pc_x[k] * sdl1_218[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pb_x, pc_x, sdl0_219, sdl0_220, sdl0_221, \
                         sdl0_222, sdl1_219, sdl1_220, sdl1_221, \
                         sdl1_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = pb_x[k] * sdl0_219[k]
                   - f_14 * pc_x[k] * sdl1_219[k];

        t_220[k] = pb_x[k] * sdl0_220[k]
                   - f_14 * pc_x[k] * sdl1_220[k];

        t_221[k] = pb_x[k] * sdl0_221[k]
                   - f_14 * pc_x[k] * sdl1_221[k];

        t_222[k] = pb_x[k] * sdl0_222[k]
                   - f_14 * pc_x[k] * sdl1_222[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pb_x, pc_x, pc_y, sdl0_224, sdl0_225, \
                         sdk_107, sdk_180, sdl1_224, sdl1_225, sfk_179, \
                         sfk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * sdk_107[k]
                   + f_3 * pc_y[k] * sfk_179[k];

        t_224[k] = pb_x[k] * sdl0_224[k]
                   - f_14 * pc_x[k] * sdl1_224[k];

        t_225[k] = pb_x[k] * sdl0_225[k]
                   + f_20 * sdk_180[k]
                   - f_14 * pc_x[k] * sdl1_225[k];

        t_226[k] = f_3 * pc_y[k] * sfk_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pb_x, pc_x, pc_y, pc_z, sdl0_228, sdk_72, \
                         sdk_183, sdl1_228, sfk_180, sfk_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_16 * sdk_72[k]
                   + f_3 * pc_z[k] * sfk_180[k];

        t_228[k] = pb_x[k] * sdl0_228[k]
                   + f_19 * sdk_183[k]
                   - f_14 * pc_x[k] * sdl1_228[k];

        t_229[k] = f_3 * pc_y[k] * sfk_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, pb_x, pc_x, pc_z, sdl0_230, sdl0_231, sdk_75, \
                         sdk_185, sdk_186, sdl1_230, sdl1_231, \
                         sfk_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = pb_x[k] * sdl0_230[k]
                   + f_19 * sdk_185[k]
                   - f_14 * pc_x[k] * sdl1_230[k];

        t_231[k] = pb_x[k] * sdl0_231[k]
                   + f_18 * sdk_186[k]
                   - f_14 * pc_x[k] * sdl1_231[k];

        t_232[k] = f_16 * sdk_75[k]
                   + f_3 * pc_z[k] * sfk_183[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_x, pc_x, pc_y, sdl0_234, sdl0_235, sdk_189, \
                         sdk_190, sdl1_234, sdl1_235, sfk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * pc_y[k] * sfk_185[k];

        t_234[k] = pb_x[k] * sdl0_234[k]
                   + f_18 * sdk_189[k]
                   - f_14 * pc_x[k] * sdl1_234[k];

        t_235[k] = pb_x[k] * sdl0_235[k]
                   + f_17 * sdk_190[k]
                   - f_14 * pc_x[k] * sdl1_235[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pb_x, pc_x, pc_y, pc_z, sdl0_237, sdk_78, \
                         sdk_192, sdl1_237, sfk_186, sfk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_16 * sdk_78[k]
                   + f_3 * pc_z[k] * sfk_186[k];

        t_237[k] = pb_x[k] * sdl0_237[k]
                   + f_17 * sdk_192[k]
                   - f_14 * pc_x[k] * sdl1_237[k];

        t_238[k] = f_3 * pc_y[k] * sfk_189[k];
    }
}

static auto
compute_prim_sfl_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sdl0,
                                                          const size_t sdk, const size_t sdl1,
                                                          const size_t sfi0, const size_t sfi1,
                                                          const size_t sfk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_17 = 2.0 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdl0_135 = buffer.data(sdl0 + 135);
    const auto *sdl0_138 = buffer.data(sdl0 + 138);
    const auto *sdl0_141 = buffer.data(sdl0 + 141);
    const auto *sdl0_145 = buffer.data(sdl0 + 145);
    const auto *sdl0_150 = buffer.data(sdl0 + 150);
    const auto *sdl0_156 = buffer.data(sdl0 + 156);
    const auto *sdl0_171 = buffer.data(sdl0 + 171);
    const auto *sdl0_173 = buffer.data(sdl0 + 173);
    const auto *sdl0_174 = buffer.data(sdl0 + 174);
    const auto *sdl0_175 = buffer.data(sdl0 + 175);
    const auto *sdl0_176 = buffer.data(sdl0 + 176);
    const auto *sdl0_177 = buffer.data(sdl0 + 177);
    const auto *sdl0_225 = buffer.data(sdl0 + 225);
    const auto *sdl0_239 = buffer.data(sdl0 + 239);
    const auto *sdl0_240 = buffer.data(sdl0 + 240);
    const auto *sdl0_242 = buffer.data(sdl0 + 242);
    const auto *sdl0_243 = buffer.data(sdl0 + 243);
    const auto *sdl0_245 = buffer.data(sdl0 + 245);
    const auto *sdl0_246 = buffer.data(sdl0 + 246);
    const auto *sdl0_248 = buffer.data(sdl0 + 248);
    const auto *sdl0_249 = buffer.data(sdl0 + 249);
    const auto *sdl0_250 = buffer.data(sdl0 + 250);
    const auto *sdl0_252 = buffer.data(sdl0 + 252);
    const auto *sdl0_261 = buffer.data(sdl0 + 261);
    const auto *sdl0_263 = buffer.data(sdl0 + 263);
    const auto *sdl0_264 = buffer.data(sdl0 + 264);
    const auto *sdl0_265 = buffer.data(sdl0 + 265);
    const auto *sdl0_266 = buffer.data(sdl0 + 266);
    const auto *sdl0_267 = buffer.data(sdl0 + 267);
    const auto *sdl0_269 = buffer.data(sdl0 + 269);

    const auto *sdk_82 = buffer.data(sdk + 82);
    const auto *sdk_87 = buffer.data(sdk + 87);
    const auto *sdk_100 = buffer.data(sdk + 100);
    const auto *sdk_108 = buffer.data(sdk + 108);
    const auto *sdk_110 = buffer.data(sdk + 110);
    const auto *sdk_111 = buffer.data(sdk + 111);
    const auto *sdk_113 = buffer.data(sdk + 113);
    const auto *sdk_114 = buffer.data(sdk + 114);
    const auto *sdk_117 = buffer.data(sdk + 117);
    const auto *sdk_118 = buffer.data(sdk + 118);
    const auto *sdk_122 = buffer.data(sdk + 122);
    const auto *sdk_123 = buffer.data(sdk + 123);
    const auto *sdk_128 = buffer.data(sdk + 128);
    const auto *sdk_136 = buffer.data(sdk + 136);
    const auto *sdk_137 = buffer.data(sdk + 137);
    const auto *sdk_138 = buffer.data(sdk + 138);
    const auto *sdk_139 = buffer.data(sdk + 139);
    const auto *sdk_140 = buffer.data(sdk + 140);
    const auto *sdk_141 = buffer.data(sdk + 141);
    const auto *sdk_142 = buffer.data(sdk + 142);
    const auto *sdk_143 = buffer.data(sdk + 143);
    const auto *sdk_144 = buffer.data(sdk + 144);
    const auto *sdk_146 = buffer.data(sdk + 146);
    const auto *sdk_149 = buffer.data(sdk + 149);
    const auto *sdk_153 = buffer.data(sdk + 153);
    const auto *sdk_158 = buffer.data(sdk + 158);
    const auto *sdk_164 = buffer.data(sdk + 164);
    const auto *sdk_179 = buffer.data(sdk + 179);
    const auto *sdk_180 = buffer.data(sdk + 180);
    const auto *sdk_194 = buffer.data(sdk + 194);
    const auto *sdk_195 = buffer.data(sdk + 195);
    const auto *sdk_197 = buffer.data(sdk + 197);
    const auto *sdk_198 = buffer.data(sdk + 198);
    const auto *sdk_200 = buffer.data(sdk + 200);
    const auto *sdk_201 = buffer.data(sdk + 201);
    const auto *sdk_203 = buffer.data(sdk + 203);
    const auto *sdk_204 = buffer.data(sdk + 204);
    const auto *sdk_205 = buffer.data(sdk + 205);
    const auto *sdk_207 = buffer.data(sdk + 207);
    const auto *sdk_208 = buffer.data(sdk + 208);
    const auto *sdk_209 = buffer.data(sdk + 209);
    const auto *sdk_210 = buffer.data(sdk + 210);
    const auto *sdk_211 = buffer.data(sdk + 211);
    const auto *sdk_212 = buffer.data(sdk + 212);
    const auto *sdk_213 = buffer.data(sdk + 213);
    const auto *sdk_214 = buffer.data(sdk + 214);
    const auto *sdk_215 = buffer.data(sdk + 215);

    const auto *sdl1_135 = buffer.data(sdl1 + 135);
    const auto *sdl1_138 = buffer.data(sdl1 + 138);
    const auto *sdl1_141 = buffer.data(sdl1 + 141);
    const auto *sdl1_145 = buffer.data(sdl1 + 145);
    const auto *sdl1_150 = buffer.data(sdl1 + 150);
    const auto *sdl1_156 = buffer.data(sdl1 + 156);
    const auto *sdl1_171 = buffer.data(sdl1 + 171);
    const auto *sdl1_173 = buffer.data(sdl1 + 173);
    const auto *sdl1_174 = buffer.data(sdl1 + 174);
    const auto *sdl1_175 = buffer.data(sdl1 + 175);
    const auto *sdl1_176 = buffer.data(sdl1 + 176);
    const auto *sdl1_177 = buffer.data(sdl1 + 177);
    const auto *sdl1_225 = buffer.data(sdl1 + 225);
    const auto *sdl1_239 = buffer.data(sdl1 + 239);
    const auto *sdl1_240 = buffer.data(sdl1 + 240);
    const auto *sdl1_242 = buffer.data(sdl1 + 242);
    const auto *sdl1_243 = buffer.data(sdl1 + 243);
    const auto *sdl1_245 = buffer.data(sdl1 + 245);
    const auto *sdl1_246 = buffer.data(sdl1 + 246);
    const auto *sdl1_248 = buffer.data(sdl1 + 248);
    const auto *sdl1_249 = buffer.data(sdl1 + 249);
    const auto *sdl1_250 = buffer.data(sdl1 + 250);
    const auto *sdl1_252 = buffer.data(sdl1 + 252);
    const auto *sdl1_261 = buffer.data(sdl1 + 261);
    const auto *sdl1_263 = buffer.data(sdl1 + 263);
    const auto *sdl1_264 = buffer.data(sdl1 + 264);
    const auto *sdl1_265 = buffer.data(sdl1 + 265);
    const auto *sdl1_266 = buffer.data(sdl1 + 266);
    const auto *sdl1_267 = buffer.data(sdl1 + 267);
    const auto *sdl1_269 = buffer.data(sdl1 + 269);

    const auto *sfi0_168 = buffer.data(sfi0 + 168);
    const auto *sfi0_171 = buffer.data(sfi0 + 171);
    const auto *sfi0_173 = buffer.data(sfi0 + 173);
    const auto *sfi0_174 = buffer.data(sfi0 + 174);
    const auto *sfi0_177 = buffer.data(sfi0 + 177);
    const auto *sfi0_178 = buffer.data(sfi0 + 178);
    const auto *sfi0_180 = buffer.data(sfi0 + 180);
    const auto *sfi0_182 = buffer.data(sfi0 + 182);
    const auto *sfi0_183 = buffer.data(sfi0 + 183);
    const auto *sfi0_185 = buffer.data(sfi0 + 185);
    const auto *sfi0_186 = buffer.data(sfi0 + 186);
    const auto *sfi0_188 = buffer.data(sfi0 + 188);
    const auto *sfi0_189 = buffer.data(sfi0 + 189);
    const auto *sfi0_191 = buffer.data(sfi0 + 191);
    const auto *sfi0_192 = buffer.data(sfi0 + 192);
    const auto *sfi0_193 = buffer.data(sfi0 + 193);
    const auto *sfi0_194 = buffer.data(sfi0 + 194);
    const auto *sfi0_195 = buffer.data(sfi0 + 195);
    const auto *sfi0_201 = buffer.data(sfi0 + 201);
    const auto *sfi0_205 = buffer.data(sfi0 + 205);
    const auto *sfi0_208 = buffer.data(sfi0 + 208);
    const auto *sfi0_210 = buffer.data(sfi0 + 210);
    const auto *sfi0_213 = buffer.data(sfi0 + 213);
    const auto *sfi0_214 = buffer.data(sfi0 + 214);
    const auto *sfi0_216 = buffer.data(sfi0 + 216);
    const auto *sfi0_219 = buffer.data(sfi0 + 219);
    const auto *sfi0_220 = buffer.data(sfi0 + 220);
    const auto *sfi0_221 = buffer.data(sfi0 + 221);
    const auto *sfi0_223 = buffer.data(sfi0 + 223);

    const auto *sfi1_168 = buffer.data(sfi1 + 168);
    const auto *sfi1_171 = buffer.data(sfi1 + 171);
    const auto *sfi1_173 = buffer.data(sfi1 + 173);
    const auto *sfi1_174 = buffer.data(sfi1 + 174);
    const auto *sfi1_177 = buffer.data(sfi1 + 177);
    const auto *sfi1_178 = buffer.data(sfi1 + 178);
    const auto *sfi1_180 = buffer.data(sfi1 + 180);
    const auto *sfi1_182 = buffer.data(sfi1 + 182);
    const auto *sfi1_183 = buffer.data(sfi1 + 183);
    const auto *sfi1_185 = buffer.data(sfi1 + 185);
    const auto *sfi1_186 = buffer.data(sfi1 + 186);
    const auto *sfi1_188 = buffer.data(sfi1 + 188);
    const auto *sfi1_189 = buffer.data(sfi1 + 189);
    const auto *sfi1_191 = buffer.data(sfi1 + 191);
    const auto *sfi1_192 = buffer.data(sfi1 + 192);
    const auto *sfi1_193 = buffer.data(sfi1 + 193);
    const auto *sfi1_194 = buffer.data(sfi1 + 194);
    const auto *sfi1_195 = buffer.data(sfi1 + 195);
    const auto *sfi1_201 = buffer.data(sfi1 + 201);
    const auto *sfi1_205 = buffer.data(sfi1 + 205);
    const auto *sfi1_208 = buffer.data(sfi1 + 208);
    const auto *sfi1_210 = buffer.data(sfi1 + 210);
    const auto *sfi1_213 = buffer.data(sfi1 + 213);
    const auto *sfi1_214 = buffer.data(sfi1 + 214);
    const auto *sfi1_216 = buffer.data(sfi1 + 216);
    const auto *sfi1_219 = buffer.data(sfi1 + 219);
    const auto *sfi1_220 = buffer.data(sfi1 + 220);
    const auto *sfi1_221 = buffer.data(sfi1 + 221);
    const auto *sfi1_223 = buffer.data(sfi1 + 223);

    const auto *sfk_190 = buffer.data(sfk + 190);
    const auto *sfk_194 = buffer.data(sfk + 194);
    const auto *sfk_195 = buffer.data(sfk + 195);
    const auto *sfk_200 = buffer.data(sfk + 200);
    const auto *sfk_208 = buffer.data(sfk + 208);
    const auto *sfk_209 = buffer.data(sfk + 209);
    const auto *sfk_210 = buffer.data(sfk + 210);
    const auto *sfk_211 = buffer.data(sfk + 211);
    const auto *sfk_212 = buffer.data(sfk + 212);
    const auto *sfk_213 = buffer.data(sfk + 213);
    const auto *sfk_214 = buffer.data(sfk + 214);
    const auto *sfk_215 = buffer.data(sfk + 215);
    const auto *sfk_216 = buffer.data(sfk + 216);
    const auto *sfk_218 = buffer.data(sfk + 218);
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
    const auto *sfk_252 = buffer.data(sfk + 252);
    const auto *sfk_254 = buffer.data(sfk + 254);
    const auto *sfk_255 = buffer.data(sfk + 255);
    const auto *sfk_257 = buffer.data(sfk + 257);
    const auto *sfk_258 = buffer.data(sfk + 258);
    const auto *sfk_261 = buffer.data(sfk + 261);
    const auto *sfk_262 = buffer.data(sfk + 262);
    const auto *sfk_264 = buffer.data(sfk + 264);
    const auto *sfk_266 = buffer.data(sfk + 266);
    const auto *sfk_267 = buffer.data(sfk + 267);
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
    const auto *sfk_287 = buffer.data(sfk + 287);
    const auto *sfk_288 = buffer.data(sfk + 288);

#pragma omp simd aligned(t_239, t_240, t_241, pb_x, pc_x, pc_z, sdl0_239, sdl0_240, sdk_82, \
                         sdk_194, sdk_195, sdl1_239, sdl1_240, \
                         sfk_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = pb_x[k] * sdl0_239[k]
                   + f_17 * sdk_194[k]
                   - f_14 * pc_x[k] * sdl1_239[k];

        t_240[k] = pb_x[k] * sdl0_240[k]
                   + f_0 * sdk_195[k]
                   - f_14 * pc_x[k] * sdl1_240[k];

        t_241[k] = f_16 * sdk_82[k]
                   + f_3 * pc_z[k] * sfk_190[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, pb_x, pc_x, pc_y, sdl0_242, sdl0_243, sdk_197, \
                         sdk_198, sdl1_242, sdl1_243, sfk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = pb_x[k] * sdl0_242[k]
                   + f_0 * sdk_197[k]
                   - f_14 * pc_x[k] * sdl1_242[k];

        t_243[k] = pb_x[k] * sdl0_243[k]
                   + f_0 * sdk_198[k]
                   - f_14 * pc_x[k] * sdl1_243[k];

        t_244[k] = f_3 * pc_y[k] * sfk_194[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pb_x, pc_x, pc_z, sdl0_245, sdl0_246, sdk_87, \
                         sdk_200, sdk_201, sdl1_245, sdl1_246, \
                         sfk_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = pb_x[k] * sdl0_245[k]
                   + f_0 * sdk_200[k]
                   - f_14 * pc_x[k] * sdl1_245[k];

        t_246[k] = pb_x[k] * sdl0_246[k]
                   + f_16 * sdk_201[k]
                   - f_14 * pc_x[k] * sdl1_246[k];

        t_247[k] = f_16 * sdk_87[k]
                   + f_3 * pc_z[k] * sfk_195[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pb_x, pc_x, sdl0_248, sdl0_249, sdl0_250, \
                         sdk_203, sdk_204, sdk_205, sdl1_248, sdl1_249, \
                         sdl1_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = pb_x[k] * sdl0_248[k]
                   + f_16 * sdk_203[k]
                   - f_14 * pc_x[k] * sdl1_248[k];

        t_249[k] = pb_x[k] * sdl0_249[k]
                   + f_16 * sdk_204[k]
                   - f_14 * pc_x[k] * sdl1_249[k];

        t_250[k] = pb_x[k] * sdl0_250[k]
                   + f_16 * sdk_205[k]
                   - f_14 * pc_x[k] * sdl1_250[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_x, pc_x, pc_y, sdl0_252, sdk_207, \
                         sdk_208, sdk_209, sdl1_252, sfk_200, sfk_208, \
                         sfk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_3 * pc_y[k] * sfk_200[k];

        t_252[k] = pb_x[k] * sdl0_252[k]
                   + f_16 * sdk_207[k]
                   - f_14 * pc_x[k] * sdl1_252[k];

        t_253[k] = f_15 * sdk_208[k]
                   + f_3 * pc_x[k] * sfk_208[k];

        t_254[k] = f_15 * sdk_209[k]
                   + f_3 * pc_x[k] * sfk_209[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pc_x, sdk_210, sdk_211, sdk_212, \
                         sdk_213, sdk_214, sfk_210, sfk_211, sfk_212, sfk_213, \
                         sfk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_15 * sdk_210[k]
                   + f_3 * pc_x[k] * sfk_210[k];

        t_256[k] = f_15 * sdk_211[k]
                   + f_3 * pc_x[k] * sfk_211[k];

        t_257[k] = f_15 * sdk_212[k]
                   + f_3 * pc_x[k] * sfk_212[k];

        t_258[k] = f_15 * sdk_213[k]
                   + f_3 * pc_x[k] * sfk_213[k];

        t_259[k] = f_15 * sdk_214[k]
                   + f_3 * pc_x[k] * sfk_214[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pb_x, pc_x, pc_z, sdl0_261, sdl0_263, \
                         sdk_100, sdk_215, sdl1_261, sdl1_263, sfk_208, \
                         sfk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_15 * sdk_215[k]
                   + f_3 * pc_x[k] * sfk_215[k];

        t_261[k] = pb_x[k] * sdl0_261[k]
                   - f_14 * pc_x[k] * sdl1_261[k];

        t_262[k] = f_16 * sdk_100[k]
                   + f_3 * pc_z[k] * sfk_208[k];

        t_263[k] = pb_x[k] * sdl0_263[k]
                   - f_14 * pc_x[k] * sdl1_263[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pb_x, pc_x, sdl0_264, sdl0_265, sdl0_266, \
                         sdl0_267, sdl1_264, sdl1_265, sdl1_266, \
                         sdl1_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = pb_x[k] * sdl0_264[k]
                   - f_14 * pc_x[k] * sdl1_264[k];

        t_265[k] = pb_x[k] * sdl0_265[k]
                   - f_14 * pc_x[k] * sdl1_265[k];

        t_266[k] = pb_x[k] * sdl0_266[k]
                   - f_14 * pc_x[k] * sdl1_266[k];

        t_267[k] = pb_x[k] * sdl0_267[k]
                   - f_14 * pc_x[k] * sdl1_267[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, pb_x, pc_x, pc_y, pc_z, sdl0_269, \
                         sdk_108, sdl1_269, sfi0_168, sfi1_168, sfk_215, \
                         sfk_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_3 * pc_y[k] * sfk_215[k];

        t_269[k] = pb_x[k] * sdl0_269[k]
                   - f_14 * pc_x[k] * sdl1_269[k];

        t_270[k] = f_1 * sfi0_168[k]
                   - f_2 * sfi1_168[k]
                   + f_3 * pc_x[k] * sfk_216[k];

        t_271[k] = f_0 * sdk_108[k]
                   + f_3 * pc_y[k] * sfk_216[k];

        t_272[k] = f_3 * pc_z[k] * sfk_216[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_x, pc_y, sdk_110, sfi0_171, sfi0_173, \
                         sfi1_171, sfi1_173, sfk_218, sfk_219, \
                         sfk_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_4 * sfi0_171[k]
                   - f_5 * sfi1_171[k]
                   + f_3 * pc_x[k] * sfk_219[k];

        t_274[k] = f_0 * sdk_110[k]
                   + f_3 * pc_y[k] * sfk_218[k];

        t_275[k] = f_4 * sfi0_173[k]
                   - f_5 * sfi1_173[k]
                   + f_3 * pc_x[k] * sfk_221[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, pc_y, pc_z, sdk_113, sfi0_174, \
                         sfi0_177, sfi1_174, sfi1_177, sfk_219, sfk_221, sfk_222, \
                         sfk_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * sfi0_174[k]
                   - f_7 * sfi1_174[k]
                   + f_3 * pc_x[k] * sfk_222[k];

        t_277[k] = f_3 * pc_z[k] * sfk_219[k];

        t_278[k] = f_0 * sdk_113[k]
                   + f_3 * pc_y[k] * sfk_221[k];

        t_279[k] = f_6 * sfi0_177[k]
                   - f_7 * sfi1_177[k]
                   + f_3 * pc_x[k] * sfk_225[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, sdk_117, sfi0_178, \
                         sfi0_180, sfi1_178, sfi1_180, sfk_222, sfk_225, sfk_226, \
                         sfk_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_8 * sfi0_178[k]
                   - f_9 * sfi1_178[k]
                   + f_3 * pc_x[k] * sfk_226[k];

        t_281[k] = f_3 * pc_z[k] * sfk_222[k];

        t_282[k] = f_8 * sfi0_180[k]
                   - f_9 * sfi1_180[k]
                   + f_3 * pc_x[k] * sfk_228[k];

        t_283[k] = f_0 * sdk_117[k]
                   + f_3 * pc_y[k] * sfk_225[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pc_x, pc_z, sfi0_182, sfi0_183, sfi0_185, \
                         sfi1_182, sfi1_183, sfi1_185, sfk_226, sfk_230, sfk_231, \
                         sfk_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_8 * sfi0_182[k]
                   - f_9 * sfi1_182[k]
                   + f_3 * pc_x[k] * sfk_230[k];

        t_285[k] = f_10 * sfi0_183[k]
                   - f_11 * sfi1_183[k]
                   + f_3 * pc_x[k] * sfk_231[k];

        t_286[k] = f_3 * pc_z[k] * sfk_226[k];

        t_287[k] = f_10 * sfi0_185[k]
                   - f_11 * sfi1_185[k]
                   + f_3 * pc_x[k] * sfk_233[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_x, pc_y, sdk_122, sfi0_186, sfi0_188, \
                         sfi1_186, sfi1_188, sfk_230, sfk_234, \
                         sfk_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_10 * sfi0_186[k]
                   - f_11 * sfi1_186[k]
                   + f_3 * pc_x[k] * sfk_234[k];

        t_289[k] = f_0 * sdk_122[k]
                   + f_3 * pc_y[k] * sfk_230[k];

        t_290[k] = f_10 * sfi0_188[k]
                   - f_11 * sfi1_188[k]
                   + f_3 * pc_x[k] * sfk_236[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pc_x, pc_z, sfi0_189, sfi0_191, sfi0_192, \
                         sfi1_189, sfi1_191, sfi1_192, sfk_231, sfk_237, sfk_239, \
                         sfk_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_12 * sfi0_189[k]
                   - f_13 * sfi1_189[k]
                   + f_3 * pc_x[k] * sfk_237[k];

        t_292[k] = f_3 * pc_z[k] * sfk_231[k];

        t_293[k] = f_12 * sfi0_191[k]
                   - f_13 * sfi1_191[k]
                   + f_3 * pc_x[k] * sfk_239[k];

        t_294[k] = f_12 * sfi0_192[k]
                   - f_13 * sfi1_192[k]
                   + f_3 * pc_x[k] * sfk_240[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pc_x, pc_y, sdk_128, sfi0_193, sfi0_195, \
                         sfi1_193, sfi1_195, sfk_236, sfk_241, sfk_243, \
                         sfk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_12 * sfi0_193[k]
                   - f_13 * sfi1_193[k]
                   + f_3 * pc_x[k] * sfk_241[k];

        t_296[k] = f_0 * sdk_128[k]
                   + f_3 * pc_y[k] * sfk_236[k];

        t_297[k] = f_12 * sfi0_195[k]
                   - f_13 * sfi1_195[k]
                   + f_3 * pc_x[k] * sfk_243[k];

        t_298[k] = f_3 * pc_x[k] * sfk_244[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, t_304, t_305, pc_x, sfk_245, \
                         sfk_246, sfk_247, sfk_248, sfk_249, sfk_250, \
                         sfk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_3 * pc_x[k] * sfk_245[k];

        t_300[k] = f_3 * pc_x[k] * sfk_246[k];

        t_301[k] = f_3 * pc_x[k] * sfk_247[k];

        t_302[k] = f_3 * pc_x[k] * sfk_248[k];

        t_303[k] = f_3 * pc_x[k] * sfk_249[k];

        t_304[k] = f_3 * pc_x[k] * sfk_250[k];

        t_305[k] = f_3 * pc_x[k] * sfk_251[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, pc_y, pc_z, sdk_136, sdk_138, sfi0_189, \
                         sfi0_191, sfi1_189, sfi1_191, sfk_244, \
                         sfk_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_0 * sdk_136[k]
                   + f_1 * sfi0_189[k]
                   - f_2 * sfi1_189[k]
                   + f_3 * pc_y[k] * sfk_244[k];

        t_307[k] = f_3 * pc_z[k] * sfk_244[k];

        t_308[k] = f_0 * sdk_138[k]
                   + f_4 * sfi0_191[k]
                   - f_5 * sfi1_191[k]
                   + f_3 * pc_y[k] * sfk_246[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pc_y, sdk_139, sdk_140, sdk_141, sfi0_192, \
                         sfi0_193, sfi0_194, sfi1_192, sfi1_193, sfi1_194, sfk_247, sfk_248, \
                         sfk_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_0 * sdk_139[k]
                   + f_6 * sfi0_192[k]
                   - f_7 * sfi1_192[k]
                   + f_3 * pc_y[k] * sfk_247[k];

        t_310[k] = f_0 * sdk_140[k]
                   + f_8 * sfi0_193[k]
                   - f_9 * sfi1_193[k]
                   + f_3 * pc_y[k] * sfk_248[k];

        t_311[k] = f_0 * sdk_141[k]
                   + f_10 * sfi0_194[k]
                   - f_11 * sfi1_194[k]
                   + f_3 * pc_y[k] * sfk_249[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pb_z, pc_y, pc_z, sdl0_135, sdk_142, \
                         sdk_143, sdl1_135, sfi0_195, sfi1_195, sfk_250, \
                         sfk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_0 * sdk_142[k]
                   + f_12 * sfi0_195[k]
                   - f_13 * sfi1_195[k]
                   + f_3 * pc_y[k] * sfk_250[k];

        t_313[k] = f_0 * sdk_143[k]
                   + f_3 * pc_y[k] * sfk_251[k];

        t_314[k] = f_1 * sfi0_195[k]
                   - f_2 * sfi1_195[k]
                   + f_3 * pc_z[k] * sfk_251[k];

        t_315[k] = pb_z[k] * sdl0_135[k]
                   - f_14 * pc_z[k] * sdl1_135[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pb_z, pc_y, pc_z, sdl0_138, sdk_108, \
                         sdk_144, sdk_146, sdl1_138, sfk_252, sfk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_16 * sdk_144[k]
                   + f_3 * pc_y[k] * sfk_252[k];

        t_317[k] = f_15 * sdk_108[k]
                   + f_3 * pc_z[k] * sfk_252[k];

        t_318[k] = pb_z[k] * sdl0_138[k]
                   - f_14 * pc_z[k] * sdl1_138[k];

        t_319[k] = f_16 * sdk_146[k]
                   + f_3 * pc_y[k] * sfk_254[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pb_z, pc_x, pc_y, pc_z, sdl0_141, \
                         sdk_111, sdk_149, sdl1_141, sfi0_201, sfi1_201, sfk_255, \
                         sfk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_4 * sfi0_201[k]
                   - f_5 * sfi1_201[k]
                   + f_3 * pc_x[k] * sfk_257[k];

        t_321[k] = pb_z[k] * sdl0_141[k]
                   - f_14 * pc_z[k] * sdl1_141[k];

        t_322[k] = f_15 * sdk_111[k]
                   + f_3 * pc_z[k] * sfk_255[k];

        t_323[k] = f_16 * sdk_149[k]
                   + f_3 * pc_y[k] * sfk_257[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, pb_z, pc_x, pc_z, sdl0_145, sdk_114, sdl1_145, \
                         sfi0_205, sfi1_205, sfk_258, sfk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_6 * sfi0_205[k]
                   - f_7 * sfi1_205[k]
                   + f_3 * pc_x[k] * sfk_261[k];

        t_325[k] = pb_z[k] * sdl0_145[k]
                   - f_14 * pc_z[k] * sdl1_145[k];

        t_326[k] = f_15 * sdk_114[k]
                   + f_3 * pc_z[k] * sfk_258[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, pc_x, pc_y, sdk_153, sfi0_208, sfi0_210, \
                         sfi1_208, sfi1_210, sfk_261, sfk_264, \
                         sfk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_8 * sfi0_208[k]
                   - f_9 * sfi1_208[k]
                   + f_3 * pc_x[k] * sfk_264[k];

        t_328[k] = f_16 * sdk_153[k]
                   + f_3 * pc_y[k] * sfk_261[k];

        t_329[k] = f_8 * sfi0_210[k]
                   - f_9 * sfi1_210[k]
                   + f_3 * pc_x[k] * sfk_266[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pb_z, pc_x, pc_z, sdl0_150, sdk_118, sdl1_150, \
                         sfi0_213, sfi1_213, sfk_262, sfk_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = pb_z[k] * sdl0_150[k]
                   - f_14 * pc_z[k] * sdl1_150[k];

        t_331[k] = f_15 * sdk_118[k]
                   + f_3 * pc_z[k] * sfk_262[k];

        t_332[k] = f_10 * sfi0_213[k]
                   - f_11 * sfi1_213[k]
                   + f_3 * pc_x[k] * sfk_269[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_x, pc_y, sdk_158, sfi0_214, sfi0_216, \
                         sfi1_214, sfi1_216, sfk_266, sfk_270, \
                         sfk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_10 * sfi0_214[k]
                   - f_11 * sfi1_214[k]
                   + f_3 * pc_x[k] * sfk_270[k];

        t_334[k] = f_16 * sdk_158[k]
                   + f_3 * pc_y[k] * sfk_266[k];

        t_335[k] = f_10 * sfi0_216[k]
                   - f_11 * sfi1_216[k]
                   + f_3 * pc_x[k] * sfk_272[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pb_z, pc_x, pc_z, sdl0_156, sdk_123, sdl1_156, \
                         sfi0_219, sfi1_219, sfk_267, sfk_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = pb_z[k] * sdl0_156[k]
                   - f_14 * pc_z[k] * sdl1_156[k];

        t_337[k] = f_15 * sdk_123[k]
                   + f_3 * pc_z[k] * sfk_267[k];

        t_338[k] = f_12 * sfi0_219[k]
                   - f_13 * sfi1_219[k]
                   + f_3 * pc_x[k] * sfk_275[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, sdk_164, sfi0_220, sfi0_221, \
                         sfi1_220, sfi1_221, sfk_272, sfk_276, \
                         sfk_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_12 * sfi0_220[k]
                   - f_13 * sfi1_220[k]
                   + f_3 * pc_x[k] * sfk_276[k];

        t_340[k] = f_12 * sfi0_221[k]
                   - f_13 * sfi1_221[k]
                   + f_3 * pc_x[k] * sfk_277[k];

        t_341[k] = f_16 * sdk_164[k]
                   + f_3 * pc_y[k] * sfk_272[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, t_347, pc_x, sfi0_223, sfi1_223, \
                         sfk_279, sfk_280, sfk_281, sfk_282, sfk_283, \
                         sfk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_12 * sfi0_223[k]
                   - f_13 * sfi1_223[k]
                   + f_3 * pc_x[k] * sfk_279[k];

        t_343[k] = f_3 * pc_x[k] * sfk_280[k];

        t_344[k] = f_3 * pc_x[k] * sfk_281[k];

        t_345[k] = f_3 * pc_x[k] * sfk_282[k];

        t_346[k] = f_3 * pc_x[k] * sfk_283[k];

        t_347[k] = f_3 * pc_x[k] * sfk_284[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, pb_z, pc_x, pc_z, sdl0_171, \
                         sdk_136, sdl1_171, sfk_280, sfk_285, sfk_286, \
                         sfk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_3 * pc_x[k] * sfk_285[k];

        t_349[k] = f_3 * pc_x[k] * sfk_286[k];

        t_350[k] = f_3 * pc_x[k] * sfk_287[k];

        t_351[k] = pb_z[k] * sdl0_171[k]
                   - f_14 * pc_z[k] * sdl1_171[k];

        t_352[k] = f_15 * sdk_136[k]
                   + f_3 * pc_z[k] * sfk_280[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pb_z, pc_z, sdl0_173, sdl0_174, sdl0_175, \
                         sdk_137, sdk_138, sdk_139, sdl1_173, sdl1_174, \
                         sdl1_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = pb_z[k] * sdl0_173[k]
                   + f_16 * sdk_137[k]
                   - f_14 * pc_z[k] * sdl1_173[k];

        t_354[k] = pb_z[k] * sdl0_174[k]
                   + f_0 * sdk_138[k]
                   - f_14 * pc_z[k] * sdl1_174[k];

        t_355[k] = pb_z[k] * sdl0_175[k]
                   + f_17 * sdk_139[k]
                   - f_14 * pc_z[k] * sdl1_175[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pb_z, pc_y, pc_z, sdl0_176, sdl0_177, sdk_140, \
                         sdk_141, sdk_179, sdl1_176, sdl1_177, \
                         sfk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = pb_z[k] * sdl0_176[k]
                   + f_18 * sdk_140[k]
                   - f_14 * pc_z[k] * sdl1_176[k];

        t_357[k] = pb_z[k] * sdl0_177[k]
                   + f_19 * sdk_141[k]
                   - f_14 * pc_z[k] * sdl1_177[k];

        t_358[k] = f_16 * sdk_179[k]
                   + f_3 * pc_y[k] * sfk_287[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pb_y, pc_y, pc_z, sdl0_225, sdk_143, \
                         sdk_144, sdk_180, sdl1_225, sfi0_223, sfi1_223, sfk_287, \
                         sfk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_15 * sdk_143[k]
                   + f_1 * sfi0_223[k]
                   - f_2 * sfi1_223[k]
                   + f_3 * pc_z[k] * sfk_287[k];

        t_360[k] = pb_y[k] * sdl0_225[k]
                   - f_14 * pc_y[k] * sdl1_225[k];

        t_361[k] = f_15 * sdk_180[k]
                   + f_3 * pc_y[k] * sfk_288[k];

        t_362[k] = f_16 * sdk_144[k]
                   + f_3 * pc_z[k] * sfk_288[k];
    }
}

static auto
compute_prim_sfl_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sdl0,
                                                          const size_t sdk, const size_t sdl1,
                                                          const size_t sfi0, const size_t sfi1,
                                                          const size_t sfk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_17 = 2.0 / q;
    const auto f_18 = 2.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 4.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdl0_230 = buffer.data(sdl0 + 230);
    const auto *sdl0_234 = buffer.data(sdl0 + 234);
    const auto *sdl0_239 = buffer.data(sdl0 + 239);
    const auto *sdl0_245 = buffer.data(sdl0 + 245);
    const auto *sdl0_252 = buffer.data(sdl0 + 252);
    const auto *sdl0_261 = buffer.data(sdl0 + 261);
    const auto *sdl0_263 = buffer.data(sdl0 + 263);
    const auto *sdl0_264 = buffer.data(sdl0 + 264);
    const auto *sdl0_265 = buffer.data(sdl0 + 265);
    const auto *sdl0_266 = buffer.data(sdl0 + 266);
    const auto *sdl0_267 = buffer.data(sdl0 + 267);
    const auto *sdl0_269 = buffer.data(sdl0 + 269);

    const auto *sdk_147 = buffer.data(sdk + 147);
    const auto *sdk_150 = buffer.data(sdk + 150);
    const auto *sdk_154 = buffer.data(sdk + 154);
    const auto *sdk_159 = buffer.data(sdk + 159);
    const auto *sdk_172 = buffer.data(sdk + 172);
    const auto *sdk_180 = buffer.data(sdk + 180);
    const auto *sdk_182 = buffer.data(sdk + 182);
    const auto *sdk_183 = buffer.data(sdk + 183);
    const auto *sdk_185 = buffer.data(sdk + 185);
    const auto *sdk_186 = buffer.data(sdk + 186);
    const auto *sdk_189 = buffer.data(sdk + 189);
    const auto *sdk_190 = buffer.data(sdk + 190);
    const auto *sdk_194 = buffer.data(sdk + 194);
    const auto *sdk_195 = buffer.data(sdk + 195);
    const auto *sdk_200 = buffer.data(sdk + 200);
    const auto *sdk_208 = buffer.data(sdk + 208);
    const auto *sdk_210 = buffer.data(sdk + 210);
    const auto *sdk_211 = buffer.data(sdk + 211);
    const auto *sdk_212 = buffer.data(sdk + 212);
    const auto *sdk_213 = buffer.data(sdk + 213);
    const auto *sdk_214 = buffer.data(sdk + 214);
    const auto *sdk_215 = buffer.data(sdk + 215);

    const auto *sdl1_230 = buffer.data(sdl1 + 230);
    const auto *sdl1_234 = buffer.data(sdl1 + 234);
    const auto *sdl1_239 = buffer.data(sdl1 + 239);
    const auto *sdl1_245 = buffer.data(sdl1 + 245);
    const auto *sdl1_252 = buffer.data(sdl1 + 252);
    const auto *sdl1_261 = buffer.data(sdl1 + 261);
    const auto *sdl1_263 = buffer.data(sdl1 + 263);
    const auto *sdl1_264 = buffer.data(sdl1 + 264);
    const auto *sdl1_265 = buffer.data(sdl1 + 265);
    const auto *sdl1_266 = buffer.data(sdl1 + 266);
    const auto *sdl1_267 = buffer.data(sdl1 + 267);
    const auto *sdl1_269 = buffer.data(sdl1 + 269);

    const auto *sfi0_227 = buffer.data(sfi0 + 227);
    const auto *sfi0_230 = buffer.data(sfi0 + 230);
    const auto *sfi0_234 = buffer.data(sfi0 + 234);
    const auto *sfi0_236 = buffer.data(sfi0 + 236);
    const auto *sfi0_239 = buffer.data(sfi0 + 239);
    const auto *sfi0_241 = buffer.data(sfi0 + 241);
    const auto *sfi0_242 = buffer.data(sfi0 + 242);
    const auto *sfi0_245 = buffer.data(sfi0 + 245);
    const auto *sfi0_247 = buffer.data(sfi0 + 247);
    const auto *sfi0_248 = buffer.data(sfi0 + 248);
    const auto *sfi0_249 = buffer.data(sfi0 + 249);
    const auto *sfi0_252 = buffer.data(sfi0 + 252);
    const auto *sfi0_255 = buffer.data(sfi0 + 255);
    const auto *sfi0_257 = buffer.data(sfi0 + 257);
    const auto *sfi0_258 = buffer.data(sfi0 + 258);
    const auto *sfi0_261 = buffer.data(sfi0 + 261);
    const auto *sfi0_262 = buffer.data(sfi0 + 262);
    const auto *sfi0_264 = buffer.data(sfi0 + 264);
    const auto *sfi0_266 = buffer.data(sfi0 + 266);
    const auto *sfi0_267 = buffer.data(sfi0 + 267);
    const auto *sfi0_269 = buffer.data(sfi0 + 269);
    const auto *sfi0_270 = buffer.data(sfi0 + 270);
    const auto *sfi0_272 = buffer.data(sfi0 + 272);
    const auto *sfi0_273 = buffer.data(sfi0 + 273);
    const auto *sfi0_275 = buffer.data(sfi0 + 275);
    const auto *sfi0_276 = buffer.data(sfi0 + 276);
    const auto *sfi0_277 = buffer.data(sfi0 + 277);
    const auto *sfi0_278 = buffer.data(sfi0 + 278);
    const auto *sfi0_279 = buffer.data(sfi0 + 279);

    const auto *sfi1_227 = buffer.data(sfi1 + 227);
    const auto *sfi1_230 = buffer.data(sfi1 + 230);
    const auto *sfi1_234 = buffer.data(sfi1 + 234);
    const auto *sfi1_236 = buffer.data(sfi1 + 236);
    const auto *sfi1_239 = buffer.data(sfi1 + 239);
    const auto *sfi1_241 = buffer.data(sfi1 + 241);
    const auto *sfi1_242 = buffer.data(sfi1 + 242);
    const auto *sfi1_245 = buffer.data(sfi1 + 245);
    const auto *sfi1_247 = buffer.data(sfi1 + 247);
    const auto *sfi1_248 = buffer.data(sfi1 + 248);
    const auto *sfi1_249 = buffer.data(sfi1 + 249);
    const auto *sfi1_252 = buffer.data(sfi1 + 252);
    const auto *sfi1_255 = buffer.data(sfi1 + 255);
    const auto *sfi1_257 = buffer.data(sfi1 + 257);
    const auto *sfi1_258 = buffer.data(sfi1 + 258);
    const auto *sfi1_261 = buffer.data(sfi1 + 261);
    const auto *sfi1_262 = buffer.data(sfi1 + 262);
    const auto *sfi1_264 = buffer.data(sfi1 + 264);
    const auto *sfi1_266 = buffer.data(sfi1 + 266);
    const auto *sfi1_267 = buffer.data(sfi1 + 267);
    const auto *sfi1_269 = buffer.data(sfi1 + 269);
    const auto *sfi1_270 = buffer.data(sfi1 + 270);
    const auto *sfi1_272 = buffer.data(sfi1 + 272);
    const auto *sfi1_273 = buffer.data(sfi1 + 273);
    const auto *sfi1_275 = buffer.data(sfi1 + 275);
    const auto *sfi1_276 = buffer.data(sfi1 + 276);
    const auto *sfi1_277 = buffer.data(sfi1 + 277);
    const auto *sfi1_278 = buffer.data(sfi1 + 278);
    const auto *sfi1_279 = buffer.data(sfi1 + 279);

    const auto *sfk_290 = buffer.data(sfk + 290);
    const auto *sfk_291 = buffer.data(sfk + 291);
    const auto *sfk_293 = buffer.data(sfk + 293);
    const auto *sfk_294 = buffer.data(sfk + 294);
    const auto *sfk_297 = buffer.data(sfk + 297);
    const auto *sfk_298 = buffer.data(sfk + 298);
    const auto *sfk_300 = buffer.data(sfk + 300);
    const auto *sfk_302 = buffer.data(sfk + 302);
    const auto *sfk_303 = buffer.data(sfk + 303);
    const auto *sfk_305 = buffer.data(sfk + 305);
    const auto *sfk_306 = buffer.data(sfk + 306);
    const auto *sfk_308 = buffer.data(sfk + 308);
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
    const auto *sfk_326 = buffer.data(sfk + 326);
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

#pragma omp simd aligned(t_363, t_364, t_365, pb_y, pc_x, pc_y, sdl0_230, sdk_182, sdl1_230, \
                         sfi0_227, sfi1_227, sfk_290, sfk_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_4 * sfi0_227[k]
                   - f_5 * sfi1_227[k]
                   + f_3 * pc_x[k] * sfk_291[k];

        t_364[k] = f_15 * sdk_182[k]
                   + f_3 * pc_y[k] * sfk_290[k];

        t_365[k] = pb_y[k] * sdl0_230[k]
                   - f_14 * pc_y[k] * sdl1_230[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, pc_x, pc_y, pc_z, sdk_147, sdk_185, sfi0_230, \
                         sfi1_230, sfk_291, sfk_293, sfk_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_6 * sfi0_230[k]
                   - f_7 * sfi1_230[k]
                   + f_3 * pc_x[k] * sfk_294[k];

        t_367[k] = f_16 * sdk_147[k]
                   + f_3 * pc_z[k] * sfk_291[k];

        t_368[k] = f_15 * sdk_185[k]
                   + f_3 * pc_y[k] * sfk_293[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pb_y, pc_x, pc_y, pc_z, sdl0_234, sdk_150, \
                         sdl1_234, sfi0_234, sfi1_234, sfk_294, \
                         sfk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = pb_y[k] * sdl0_234[k]
                   - f_14 * pc_y[k] * sdl1_234[k];

        t_370[k] = f_8 * sfi0_234[k]
                   - f_9 * sfi1_234[k]
                   + f_3 * pc_x[k] * sfk_298[k];

        t_371[k] = f_16 * sdk_150[k]
                   + f_3 * pc_z[k] * sfk_294[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pb_y, pc_x, pc_y, sdl0_239, sdk_189, sdl1_239, \
                         sfi0_236, sfi1_236, sfk_297, sfk_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_8 * sfi0_236[k]
                   - f_9 * sfi1_236[k]
                   + f_3 * pc_x[k] * sfk_300[k];

        t_373[k] = f_15 * sdk_189[k]
                   + f_3 * pc_y[k] * sfk_297[k];

        t_374[k] = pb_y[k] * sdl0_239[k]
                   - f_14 * pc_y[k] * sdl1_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, pc_z, sdk_154, sfi0_239, sfi0_241, \
                         sfi1_239, sfi1_241, sfk_298, sfk_303, \
                         sfk_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_10 * sfi0_239[k]
                   - f_11 * sfi1_239[k]
                   + f_3 * pc_x[k] * sfk_303[k];

        t_376[k] = f_16 * sdk_154[k]
                   + f_3 * pc_z[k] * sfk_298[k];

        t_377[k] = f_10 * sfi0_241[k]
                   - f_11 * sfi1_241[k]
                   + f_3 * pc_x[k] * sfk_305[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, pb_y, pc_x, pc_y, sdl0_245, sdk_194, sdl1_245, \
                         sfi0_242, sfi1_242, sfk_302, sfk_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_10 * sfi0_242[k]
                   - f_11 * sfi1_242[k]
                   + f_3 * pc_x[k] * sfk_306[k];

        t_379[k] = f_15 * sdk_194[k]
                   + f_3 * pc_y[k] * sfk_302[k];

        t_380[k] = pb_y[k] * sdl0_245[k]
                   - f_14 * pc_y[k] * sdl1_245[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, pc_x, pc_z, sdk_159, sfi0_245, sfi0_247, \
                         sfi1_245, sfi1_247, sfk_303, sfk_309, \
                         sfk_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_12 * sfi0_245[k]
                   - f_13 * sfi1_245[k]
                   + f_3 * pc_x[k] * sfk_309[k];

        t_382[k] = f_16 * sdk_159[k]
                   + f_3 * pc_z[k] * sfk_303[k];

        t_383[k] = f_12 * sfi0_247[k]
                   - f_13 * sfi1_247[k]
                   + f_3 * pc_x[k] * sfk_311[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pc_x, pc_y, sdk_200, sfi0_248, sfi0_249, \
                         sfi1_248, sfi1_249, sfk_308, sfk_312, \
                         sfk_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_12 * sfi0_248[k]
                   - f_13 * sfi1_248[k]
                   + f_3 * pc_x[k] * sfk_312[k];

        t_385[k] = f_12 * sfi0_249[k]
                   - f_13 * sfi1_249[k]
                   + f_3 * pc_x[k] * sfk_313[k];

        t_386[k] = f_15 * sdk_200[k]
                   + f_3 * pc_y[k] * sfk_308[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, t_392, pb_y, pc_x, pc_y, sdl0_252, \
                         sdl1_252, sfk_316, sfk_317, sfk_318, sfk_319, \
                         sfk_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = pb_y[k] * sdl0_252[k]
                   - f_14 * pc_y[k] * sdl1_252[k];

        t_388[k] = f_3 * pc_x[k] * sfk_316[k];

        t_389[k] = f_3 * pc_x[k] * sfk_317[k];

        t_390[k] = f_3 * pc_x[k] * sfk_318[k];

        t_391[k] = f_3 * pc_x[k] * sfk_319[k];

        t_392[k] = f_3 * pc_x[k] * sfk_320[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, pb_y, pc_x, pc_y, sdl0_261, sdk_208, \
                         sdl1_261, sfk_321, sfk_322, sfk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = f_3 * pc_x[k] * sfk_321[k];

        t_394[k] = f_3 * pc_x[k] * sfk_322[k];

        t_395[k] = f_3 * pc_x[k] * sfk_323[k];

        t_396[k] = pb_y[k] * sdl0_261[k]
                   + f_20 * sdk_208[k]
                   - f_14 * pc_y[k] * sdl1_261[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, pb_y, pc_y, pc_z, sdl0_263, sdl0_264, sdk_172, \
                         sdk_210, sdk_211, sdl1_263, sdl1_264, \
                         sfk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_16 * sdk_172[k]
                   + f_3 * pc_z[k] * sfk_316[k];

        t_398[k] = pb_y[k] * sdl0_263[k]
                   + f_19 * sdk_210[k]
                   - f_14 * pc_y[k] * sdl1_263[k];

        t_399[k] = pb_y[k] * sdl0_264[k]
                   + f_18 * sdk_211[k]
                   - f_14 * pc_y[k] * sdl1_264[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pb_y, pc_y, sdl0_265, sdl0_266, sdl0_267, \
                         sdk_212, sdk_213, sdk_214, sdl1_265, sdl1_266, \
                         sdl1_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = pb_y[k] * sdl0_265[k]
                   + f_17 * sdk_212[k]
                   - f_14 * pc_y[k] * sdl1_265[k];

        t_401[k] = pb_y[k] * sdl0_266[k]
                   + f_0 * sdk_213[k]
                   - f_14 * pc_y[k] * sdl1_266[k];

        t_402[k] = pb_y[k] * sdl0_267[k]
                   + f_16 * sdk_214[k]
                   - f_14 * pc_y[k] * sdl1_267[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, pb_y, pc_x, pc_y, sdl0_269, sdk_215, \
                         sdl1_269, sfi0_252, sfi1_252, sfk_323, \
                         sfk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_15 * sdk_215[k]
                   + f_3 * pc_y[k] * sfk_323[k];

        t_404[k] = pb_y[k] * sdl0_269[k]
                   - f_14 * pc_y[k] * sdl1_269[k];

        t_405[k] = f_1 * sfi0_252[k]
                   - f_2 * sfi1_252[k]
                   + f_3 * pc_x[k] * sfk_324[k];

        t_406[k] = f_3 * pc_y[k] * sfk_324[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, pc_x, pc_y, pc_z, sdk_180, sfi0_255, \
                         sfi0_257, sfi1_255, sfi1_257, sfk_324, sfk_326, sfk_327, \
                         sfk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_0 * sdk_180[k]
                   + f_3 * pc_z[k] * sfk_324[k];

        t_408[k] = f_4 * sfi0_255[k]
                   - f_5 * sfi1_255[k]
                   + f_3 * pc_x[k] * sfk_327[k];

        t_409[k] = f_3 * pc_y[k] * sfk_326[k];

        t_410[k] = f_4 * sfi0_257[k]
                   - f_5 * sfi1_257[k]
                   + f_3 * pc_x[k] * sfk_329[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, pc_x, pc_y, pc_z, sdk_183, sfi0_258, \
                         sfi0_261, sfi1_258, sfi1_261, sfk_327, sfk_329, sfk_330, \
                         sfk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_6 * sfi0_258[k]
                   - f_7 * sfi1_258[k]
                   + f_3 * pc_x[k] * sfk_330[k];

        t_412[k] = f_0 * sdk_183[k]
                   + f_3 * pc_z[k] * sfk_327[k];

        t_413[k] = f_3 * pc_y[k] * sfk_329[k];

        t_414[k] = f_6 * sfi0_261[k]
                   - f_7 * sfi1_261[k]
                   + f_3 * pc_x[k] * sfk_333[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, pc_x, pc_y, pc_z, sdk_186, sfi0_262, \
                         sfi0_264, sfi1_262, sfi1_264, sfk_330, sfk_333, sfk_334, \
                         sfk_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_8 * sfi0_262[k]
                   - f_9 * sfi1_262[k]
                   + f_3 * pc_x[k] * sfk_334[k];

        t_416[k] = f_0 * sdk_186[k]
                   + f_3 * pc_z[k] * sfk_330[k];

        t_417[k] = f_8 * sfi0_264[k]
                   - f_9 * sfi1_264[k]
                   + f_3 * pc_x[k] * sfk_336[k];

        t_418[k] = f_3 * pc_y[k] * sfk_333[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, pc_x, pc_z, sdk_190, sfi0_266, sfi0_267, \
                         sfi1_266, sfi1_267, sfk_334, sfk_338, \
                         sfk_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_8 * sfi0_266[k]
                   - f_9 * sfi1_266[k]
                   + f_3 * pc_x[k] * sfk_338[k];

        t_420[k] = f_10 * sfi0_267[k]
                   - f_11 * sfi1_267[k]
                   + f_3 * pc_x[k] * sfk_339[k];

        t_421[k] = f_0 * sdk_190[k]
                   + f_3 * pc_z[k] * sfk_334[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, pc_x, pc_y, sfi0_269, sfi0_270, sfi0_272, \
                         sfi1_269, sfi1_270, sfi1_272, sfk_338, sfk_341, sfk_342, \
                         sfk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_10 * sfi0_269[k]
                   - f_11 * sfi1_269[k]
                   + f_3 * pc_x[k] * sfk_341[k];

        t_423[k] = f_10 * sfi0_270[k]
                   - f_11 * sfi1_270[k]
                   + f_3 * pc_x[k] * sfk_342[k];

        t_424[k] = f_3 * pc_y[k] * sfk_338[k];

        t_425[k] = f_10 * sfi0_272[k]
                   - f_11 * sfi1_272[k]
                   + f_3 * pc_x[k] * sfk_344[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pc_x, pc_z, sdk_195, sfi0_273, sfi0_275, \
                         sfi1_273, sfi1_275, sfk_339, sfk_345, \
                         sfk_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_12 * sfi0_273[k]
                   - f_13 * sfi1_273[k]
                   + f_3 * pc_x[k] * sfk_345[k];

        t_427[k] = f_0 * sdk_195[k]
                   + f_3 * pc_z[k] * sfk_339[k];

        t_428[k] = f_12 * sfi0_275[k]
                   - f_13 * sfi1_275[k]
                   + f_3 * pc_x[k] * sfk_347[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pc_x, pc_y, sfi0_276, sfi0_277, sfi0_279, \
                         sfi1_276, sfi1_277, sfi1_279, sfk_344, sfk_348, sfk_349, \
                         sfk_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_12 * sfi0_276[k]
                   - f_13 * sfi1_276[k]
                   + f_3 * pc_x[k] * sfk_348[k];

        t_430[k] = f_12 * sfi0_277[k]
                   - f_13 * sfi1_277[k]
                   + f_3 * pc_x[k] * sfk_349[k];

        t_431[k] = f_3 * pc_y[k] * sfk_344[k];

        t_432[k] = f_12 * sfi0_279[k]
                   - f_13 * sfi1_279[k]
                   + f_3 * pc_x[k] * sfk_351[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, t_437, t_438, t_439, pc_x, sfk_352, \
                         sfk_353, sfk_354, sfk_355, sfk_356, sfk_357, \
                         sfk_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_3 * pc_x[k] * sfk_352[k];

        t_434[k] = f_3 * pc_x[k] * sfk_353[k];

        t_435[k] = f_3 * pc_x[k] * sfk_354[k];

        t_436[k] = f_3 * pc_x[k] * sfk_355[k];

        t_437[k] = f_3 * pc_x[k] * sfk_356[k];

        t_438[k] = f_3 * pc_x[k] * sfk_357[k];

        t_439[k] = f_3 * pc_x[k] * sfk_358[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, pc_x, pc_y, pc_z, sdk_208, sfi0_273, \
                         sfi0_275, sfi1_273, sfi1_275, sfk_352, sfk_354, \
                         sfk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_3 * pc_x[k] * sfk_359[k];

        t_441[k] = f_1 * sfi0_273[k]
                   - f_2 * sfi1_273[k]
                   + f_3 * pc_y[k] * sfk_352[k];

        t_442[k] = f_0 * sdk_208[k]
                   + f_3 * pc_z[k] * sfk_352[k];

        t_443[k] = f_4 * sfi0_275[k]
                   - f_5 * sfi1_275[k]
                   + f_3 * pc_y[k] * sfk_354[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pc_y, sfi0_276, sfi0_277, sfi0_278, sfi1_276, \
                         sfi1_277, sfi1_278, sfk_355, sfk_356, \
                         sfk_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_6 * sfi0_276[k]
                   - f_7 * sfi1_276[k]
                   + f_3 * pc_y[k] * sfk_355[k];

        t_445[k] = f_8 * sfi0_277[k]
                   - f_9 * sfi1_277[k]
                   + f_3 * pc_y[k] * sfk_356[k];

        t_446[k] = f_10 * sfi0_278[k]
                   - f_11 * sfi1_278[k]
                   + f_3 * pc_y[k] * sfk_357[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pc_y, pc_z, sdk_215, sfi0_279, sfi1_279, \
                         sfk_358, sfk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_12 * sfi0_279[k]
                   - f_13 * sfi1_279[k]
                   + f_3 * pc_y[k] * sfk_358[k];

        t_448[k] = f_3 * pc_y[k] * sfk_359[k];

        t_449[k] = f_0 * sdk_215[k]
                   + f_1 * sfi0_279[k]
                   - f_2 * sfi1_279[k]
                   + f_3 * pc_z[k] * sfk_359[k];
    }
}

auto
compute_prim_sfl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sdl0, const size_t sdk,
                                                   const size_t sdl1, const size_t sfi0,
                                                   const size_t sfi1, const size_t sfk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sfl_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sdl0, sdk,
                                                              sdl1, sfi0, sfi1, sfk, ncols,
                                                              gamma, p, q);

    compute_prim_sfl_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sdl0, sdk,
                                                              sdl1, sfi0, sfi1, sfk, ncols,
                                                              gamma, p, q);

    compute_prim_sfl_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sdl0, sdk,
                                                              sdl1, sfi0, sfi1, sfk, ncols,
                                                              gamma, p, q);

    compute_prim_sfl_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sdl0, sdk,
                                                              sdl1, sfi0, sfi1, sfk, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
