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


#include "SimdThreeCenterElectronRepulsionVrrRecSFK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sfk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sdk0,
                                                          const size_t sdi, const size_t sdk1,
                                                          const size_t sfh0, const size_t sfh1,
                                                          const size_t sfi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 2.0 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 3.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdk0_0 = buffer.data(sdk0 + 0);
    const auto *sdk0_3 = buffer.data(sdk0 + 3);
    const auto *sdk0_5 = buffer.data(sdk0 + 5);
    const auto *sdk0_6 = buffer.data(sdk0 + 6);
    const auto *sdk0_9 = buffer.data(sdk0 + 9);
    const auto *sdk0_10 = buffer.data(sdk0 + 10);
    const auto *sdk0_12 = buffer.data(sdk0 + 12);
    const auto *sdk0_14 = buffer.data(sdk0 + 14);
    const auto *sdk0_15 = buffer.data(sdk0 + 15);
    const auto *sdk0_17 = buffer.data(sdk0 + 17);
    const auto *sdk0_18 = buffer.data(sdk0 + 18);
    const auto *sdk0_20 = buffer.data(sdk0 + 20);
    const auto *sdk0_28 = buffer.data(sdk0 + 28);
    const auto *sdk0_35 = buffer.data(sdk0 + 35);
    const auto *sdk0_108 = buffer.data(sdk0 + 108);
    const auto *sdk0_111 = buffer.data(sdk0 + 111);
    const auto *sdk0_113 = buffer.data(sdk0 + 113);
    const auto *sdk0_114 = buffer.data(sdk0 + 114);
    const auto *sdk0_117 = buffer.data(sdk0 + 117);
    const auto *sdk0_118 = buffer.data(sdk0 + 118);
    const auto *sdk0_120 = buffer.data(sdk0 + 120);

    const auto *sdi_0 = buffer.data(sdi + 0);
    const auto *sdi_1 = buffer.data(sdi + 1);
    const auto *sdi_2 = buffer.data(sdi + 2);
    const auto *sdi_3 = buffer.data(sdi + 3);
    const auto *sdi_5 = buffer.data(sdi + 5);
    const auto *sdi_6 = buffer.data(sdi + 6);
    const auto *sdi_7 = buffer.data(sdi + 7);
    const auto *sdi_8 = buffer.data(sdi + 8);
    const auto *sdi_9 = buffer.data(sdi + 9);
    const auto *sdi_10 = buffer.data(sdi + 10);
    const auto *sdi_11 = buffer.data(sdi + 11);
    const auto *sdi_12 = buffer.data(sdi + 12);
    const auto *sdi_13 = buffer.data(sdi + 13);
    const auto *sdi_14 = buffer.data(sdi + 14);
    const auto *sdi_15 = buffer.data(sdi + 15);
    const auto *sdi_17 = buffer.data(sdi + 17);
    const auto *sdi_18 = buffer.data(sdi + 18);
    const auto *sdi_20 = buffer.data(sdi + 20);
    const auto *sdi_21 = buffer.data(sdi + 21);
    const auto *sdi_22 = buffer.data(sdi + 22);
    const auto *sdi_23 = buffer.data(sdi + 23);
    const auto *sdi_24 = buffer.data(sdi + 24);
    const auto *sdi_25 = buffer.data(sdi + 25);
    const auto *sdi_26 = buffer.data(sdi + 26);
    const auto *sdi_27 = buffer.data(sdi + 27);
    const auto *sdi_28 = buffer.data(sdi + 28);
    const auto *sdi_30 = buffer.data(sdi + 30);
    const auto *sdi_33 = buffer.data(sdi + 33);
    const auto *sdi_37 = buffer.data(sdi + 37);
    const auto *sdi_49 = buffer.data(sdi + 49);
    const auto *sdi_50 = buffer.data(sdi + 50);
    const auto *sdi_51 = buffer.data(sdi + 51);
    const auto *sdi_52 = buffer.data(sdi + 52);
    const auto *sdi_53 = buffer.data(sdi + 53);
    const auto *sdi_54 = buffer.data(sdi + 54);
    const auto *sdi_55 = buffer.data(sdi + 55);
    const auto *sdi_77 = buffer.data(sdi + 77);
    const auto *sdi_78 = buffer.data(sdi + 78);
    const auto *sdi_79 = buffer.data(sdi + 79);
    const auto *sdi_80 = buffer.data(sdi + 80);
    const auto *sdi_81 = buffer.data(sdi + 81);
    const auto *sdi_82 = buffer.data(sdi + 82);
    const auto *sdi_83 = buffer.data(sdi + 83);
    const auto *sdi_84 = buffer.data(sdi + 84);
    const auto *sdi_87 = buffer.data(sdi + 87);
    const auto *sdi_89 = buffer.data(sdi + 89);
    const auto *sdi_90 = buffer.data(sdi + 90);
    const auto *sdi_93 = buffer.data(sdi + 93);
    const auto *sdi_94 = buffer.data(sdi + 94);
    const auto *sdi_96 = buffer.data(sdi + 96);

    const auto *sdk1_0 = buffer.data(sdk1 + 0);
    const auto *sdk1_3 = buffer.data(sdk1 + 3);
    const auto *sdk1_5 = buffer.data(sdk1 + 5);
    const auto *sdk1_6 = buffer.data(sdk1 + 6);
    const auto *sdk1_9 = buffer.data(sdk1 + 9);
    const auto *sdk1_10 = buffer.data(sdk1 + 10);
    const auto *sdk1_12 = buffer.data(sdk1 + 12);
    const auto *sdk1_14 = buffer.data(sdk1 + 14);
    const auto *sdk1_15 = buffer.data(sdk1 + 15);
    const auto *sdk1_17 = buffer.data(sdk1 + 17);
    const auto *sdk1_18 = buffer.data(sdk1 + 18);
    const auto *sdk1_20 = buffer.data(sdk1 + 20);
    const auto *sdk1_28 = buffer.data(sdk1 + 28);
    const auto *sdk1_35 = buffer.data(sdk1 + 35);
    const auto *sdk1_108 = buffer.data(sdk1 + 108);
    const auto *sdk1_111 = buffer.data(sdk1 + 111);
    const auto *sdk1_113 = buffer.data(sdk1 + 113);
    const auto *sdk1_114 = buffer.data(sdk1 + 114);
    const auto *sdk1_117 = buffer.data(sdk1 + 117);
    const auto *sdk1_118 = buffer.data(sdk1 + 118);
    const auto *sdk1_120 = buffer.data(sdk1 + 120);

    const auto *sfh0_0 = buffer.data(sfh0 + 0);
    const auto *sfh0_3 = buffer.data(sfh0 + 3);
    const auto *sfh0_5 = buffer.data(sfh0 + 5);
    const auto *sfh0_6 = buffer.data(sfh0 + 6);
    const auto *sfh0_9 = buffer.data(sfh0 + 9);
    const auto *sfh0_10 = buffer.data(sfh0 + 10);
    const auto *sfh0_12 = buffer.data(sfh0 + 12);
    const auto *sfh0_14 = buffer.data(sfh0 + 14);
    const auto *sfh0_15 = buffer.data(sfh0 + 15);
    const auto *sfh0_17 = buffer.data(sfh0 + 17);
    const auto *sfh0_18 = buffer.data(sfh0 + 18);
    const auto *sfh0_19 = buffer.data(sfh0 + 19);
    const auto *sfh0_20 = buffer.data(sfh0 + 20);
    const auto *sfh0_36 = buffer.data(sfh0 + 36);
    const auto *sfh0_38 = buffer.data(sfh0 + 38);
    const auto *sfh0_39 = buffer.data(sfh0 + 39);
    const auto *sfh0_40 = buffer.data(sfh0 + 40);
    const auto *sfh0_41 = buffer.data(sfh0 + 41);
    const auto *sfh0_59 = buffer.data(sfh0 + 59);
    const auto *sfh0_60 = buffer.data(sfh0 + 60);
    const auto *sfh0_61 = buffer.data(sfh0 + 61);
    const auto *sfh0_62 = buffer.data(sfh0 + 62);

    const auto *sfh1_0 = buffer.data(sfh1 + 0);
    const auto *sfh1_3 = buffer.data(sfh1 + 3);
    const auto *sfh1_5 = buffer.data(sfh1 + 5);
    const auto *sfh1_6 = buffer.data(sfh1 + 6);
    const auto *sfh1_9 = buffer.data(sfh1 + 9);
    const auto *sfh1_10 = buffer.data(sfh1 + 10);
    const auto *sfh1_12 = buffer.data(sfh1 + 12);
    const auto *sfh1_14 = buffer.data(sfh1 + 14);
    const auto *sfh1_15 = buffer.data(sfh1 + 15);
    const auto *sfh1_17 = buffer.data(sfh1 + 17);
    const auto *sfh1_18 = buffer.data(sfh1 + 18);
    const auto *sfh1_19 = buffer.data(sfh1 + 19);
    const auto *sfh1_20 = buffer.data(sfh1 + 20);
    const auto *sfh1_36 = buffer.data(sfh1 + 36);
    const auto *sfh1_38 = buffer.data(sfh1 + 38);
    const auto *sfh1_39 = buffer.data(sfh1 + 39);
    const auto *sfh1_40 = buffer.data(sfh1 + 40);
    const auto *sfh1_41 = buffer.data(sfh1 + 41);
    const auto *sfh1_59 = buffer.data(sfh1 + 59);
    const auto *sfh1_60 = buffer.data(sfh1 + 60);
    const auto *sfh1_61 = buffer.data(sfh1 + 61);
    const auto *sfh1_62 = buffer.data(sfh1 + 62);

    const auto *sfi_0 = buffer.data(sfi + 0);
    const auto *sfi_2 = buffer.data(sfi + 2);
    const auto *sfi_3 = buffer.data(sfi + 3);
    const auto *sfi_5 = buffer.data(sfi + 5);
    const auto *sfi_6 = buffer.data(sfi + 6);
    const auto *sfi_9 = buffer.data(sfi + 9);
    const auto *sfi_10 = buffer.data(sfi + 10);
    const auto *sfi_12 = buffer.data(sfi + 12);
    const auto *sfi_14 = buffer.data(sfi + 14);
    const auto *sfi_15 = buffer.data(sfi + 15);
    const auto *sfi_17 = buffer.data(sfi + 17);
    const auto *sfi_18 = buffer.data(sfi + 18);
    const auto *sfi_20 = buffer.data(sfi + 20);
    const auto *sfi_21 = buffer.data(sfi + 21);
    const auto *sfi_22 = buffer.data(sfi + 22);
    const auto *sfi_23 = buffer.data(sfi + 23);
    const auto *sfi_24 = buffer.data(sfi + 24);
    const auto *sfi_25 = buffer.data(sfi + 25);
    const auto *sfi_26 = buffer.data(sfi + 26);
    const auto *sfi_27 = buffer.data(sfi + 27);
    const auto *sfi_28 = buffer.data(sfi + 28);
    const auto *sfi_30 = buffer.data(sfi + 30);
    const auto *sfi_31 = buffer.data(sfi + 31);
    const auto *sfi_33 = buffer.data(sfi + 33);
    const auto *sfi_34 = buffer.data(sfi + 34);
    const auto *sfi_37 = buffer.data(sfi + 37);
    const auto *sfi_38 = buffer.data(sfi + 38);
    const auto *sfi_42 = buffer.data(sfi + 42);
    const auto *sfi_49 = buffer.data(sfi + 49);
    const auto *sfi_50 = buffer.data(sfi + 50);
    const auto *sfi_51 = buffer.data(sfi + 51);
    const auto *sfi_52 = buffer.data(sfi + 52);
    const auto *sfi_53 = buffer.data(sfi + 53);
    const auto *sfi_54 = buffer.data(sfi + 54);
    const auto *sfi_55 = buffer.data(sfi + 55);
    const auto *sfi_56 = buffer.data(sfi + 56);
    const auto *sfi_58 = buffer.data(sfi + 58);
    const auto *sfi_59 = buffer.data(sfi + 59);
    const auto *sfi_61 = buffer.data(sfi + 61);
    const auto *sfi_62 = buffer.data(sfi + 62);
    const auto *sfi_65 = buffer.data(sfi + 65);
    const auto *sfi_66 = buffer.data(sfi + 66);
    const auto *sfi_70 = buffer.data(sfi + 70);
    const auto *sfi_77 = buffer.data(sfi + 77);
    const auto *sfi_78 = buffer.data(sfi + 78);
    const auto *sfi_79 = buffer.data(sfi + 79);
    const auto *sfi_80 = buffer.data(sfi + 80);
    const auto *sfi_81 = buffer.data(sfi + 81);
    const auto *sfi_82 = buffer.data(sfi + 82);
    const auto *sfi_83 = buffer.data(sfi + 83);
    const auto *sfi_84 = buffer.data(sfi + 84);
    const auto *sfi_86 = buffer.data(sfi + 86);
    const auto *sfi_87 = buffer.data(sfi + 87);
    const auto *sfi_89 = buffer.data(sfi + 89);
    const auto *sfi_90 = buffer.data(sfi + 90);
    const auto *sfi_93 = buffer.data(sfi + 93);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sdi_0, sdi_3, sfh0_0, sfh0_3, \
                         sfh1_0, sfh1_3, sfi_0, sfi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sdi_0[k]
                 + f_1 * sfh0_0[k]
                 - f_2 * sfh1_0[k]
                 + f_3 * pc_x[k] * sfi_0[k];

        t_1[k] = f_3 * pc_y[k] * sfi_0[k];

        t_2[k] = f_3 * pc_z[k] * sfi_0[k];

        t_3[k] = f_0 * sdi_3[k]
                 + f_4 * sfh0_3[k]
                 - f_5 * sfh1_3[k]
                 + f_3 * pc_x[k] * sfi_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sdi_5, sdi_6, sfh0_5, sfh0_6, sfh1_5, \
                         sfh1_6, sfi_2, sfi_5, sfi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sfi_2[k];

        t_5[k] = f_0 * sdi_5[k]
                 + f_4 * sfh0_5[k]
                 - f_5 * sfh1_5[k]
                 + f_3 * pc_x[k] * sfi_5[k];

        t_6[k] = f_0 * sdi_6[k]
                 + f_6 * sfh0_6[k]
                 - f_7 * sfh1_6[k]
                 + f_3 * pc_x[k] * sfi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sdi_9, sfh0_9, sfh1_9, sfi_3, sfi_5, \
                         sfi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sfi_3[k];

        t_8[k] = f_3 * pc_y[k] * sfi_5[k];

        t_9[k] = f_0 * sdi_9[k]
                 + f_6 * sfh0_9[k]
                 - f_7 * sfh1_9[k]
                 + f_3 * pc_x[k] * sfi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sdi_10, sdi_12, sfh0_10, sfh0_12, \
                         sfh1_10, sfh1_12, sfi_6, sfi_10, sfi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sdi_10[k]
                  + f_8 * sfh0_10[k]
                  - f_9 * sfh1_10[k]
                  + f_3 * pc_x[k] * sfi_10[k];

        t_11[k] = f_3 * pc_z[k] * sfi_6[k];

        t_12[k] = f_0 * sdi_12[k]
                  + f_8 * sfh0_12[k]
                  - f_9 * sfh1_12[k]
                  + f_3 * pc_x[k] * sfi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, sdi_14, sdi_15, sfh0_14, sfh0_15, \
                         sfh1_14, sfh1_15, sfi_9, sfi_14, sfi_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sfi_9[k];

        t_14[k] = f_0 * sdi_14[k]
                  + f_8 * sfh0_14[k]
                  - f_9 * sfh1_14[k]
                  + f_3 * pc_x[k] * sfi_14[k];

        t_15[k] = f_0 * sdi_15[k]
                  + f_10 * sfh0_15[k]
                  - f_11 * sfh1_15[k]
                  + f_3 * pc_x[k] * sfi_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, sdi_17, sdi_18, sfh0_17, sfh0_18, \
                         sfh1_17, sfh1_18, sfi_10, sfi_17, sfi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * sfi_10[k];

        t_17[k] = f_0 * sdi_17[k]
                  + f_10 * sfh0_17[k]
                  - f_11 * sfh1_17[k]
                  + f_3 * pc_x[k] * sfi_17[k];

        t_18[k] = f_0 * sdi_18[k]
                  + f_10 * sfh0_18[k]
                  - f_11 * sfh1_18[k]
                  + f_3 * pc_x[k] * sfi_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pc_x, pc_y, sdi_20, sdi_21, sdi_22, sfh0_20, \
                         sfh1_20, sfi_14, sfi_20, sfi_21, sfi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sfi_14[k];

        t_20[k] = f_0 * sdi_20[k]
                  + f_10 * sfh0_20[k]
                  - f_11 * sfh1_20[k]
                  + f_3 * pc_x[k] * sfi_20[k];

        t_21[k] = f_0 * sdi_21[k]
                  + f_3 * pc_x[k] * sfi_21[k];

        t_22[k] = f_0 * sdi_22[k]
                  + f_3 * pc_x[k] * sfi_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pc_x, sdi_23, sdi_24, sdi_25, sdi_26, \
                         sdi_27, sfi_23, sfi_24, sfi_25, sfi_26, \
                         sfi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * sdi_23[k]
                  + f_3 * pc_x[k] * sfi_23[k];

        t_24[k] = f_0 * sdi_24[k]
                  + f_3 * pc_x[k] * sfi_24[k];

        t_25[k] = f_0 * sdi_25[k]
                  + f_3 * pc_x[k] * sfi_25[k];

        t_26[k] = f_0 * sdi_26[k]
                  + f_3 * pc_x[k] * sfi_26[k];

        t_27[k] = f_0 * sdi_27[k]
                  + f_3 * pc_x[k] * sfi_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pc_y, pc_z, sfh0_15, sfh0_17, sfh0_18, \
                         sfh1_15, sfh1_17, sfh1_18, sfi_21, sfi_23, \
                         sfi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * sfh0_15[k]
                  - f_2 * sfh1_15[k]
                  + f_3 * pc_y[k] * sfi_21[k];

        t_29[k] = f_3 * pc_z[k] * sfi_21[k];

        t_30[k] = f_4 * sfh0_17[k]
                  - f_5 * sfh1_17[k]
                  + f_3 * pc_y[k] * sfi_23[k];

        t_31[k] = f_6 * sfh0_18[k]
                  - f_7 * sfh1_18[k]
                  + f_3 * pc_y[k] * sfi_24[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_y, pc_z, sfh0_19, sfh0_20, sfh1_19, \
                         sfh1_20, sfi_25, sfi_26, sfi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_8 * sfh0_19[k]
                  - f_9 * sfh1_19[k]
                  + f_3 * pc_y[k] * sfi_25[k];

        t_33[k] = f_10 * sfh0_20[k]
                  - f_11 * sfh1_20[k]
                  + f_3 * pc_y[k] * sfi_26[k];

        t_34[k] = f_3 * pc_y[k] * sfi_27[k];

        t_35[k] = f_1 * sfh0_20[k]
                  - f_2 * sfh1_20[k]
                  + f_3 * pc_z[k] * sfi_27[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, sdk0_0, sdk0_3, sdi_0, \
                         sdi_1, sdk1_0, sdk1_3, sfi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_y[k] * sdk0_0[k]
                  - f_12 * pc_y[k] * sdk1_0[k];

        t_37[k] = f_13 * sdi_0[k]
                  + f_3 * pc_y[k] * sfi_28[k];

        t_38[k] = f_3 * pc_z[k] * sfi_28[k];

        t_39[k] = pb_y[k] * sdk0_3[k]
                  + f_14 * sdi_1[k]
                  - f_12 * pc_y[k] * sdk1_3[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pc_y, pc_z, sdk0_5, sdk0_6, sdi_2, \
                         sdi_3, sdk1_5, sdk1_6, sfi_30, sfi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_13 * sdi_2[k]
                  + f_3 * pc_y[k] * sfi_30[k];

        t_41[k] = pb_y[k] * sdk0_5[k]
                  - f_12 * pc_y[k] * sdk1_5[k];

        t_42[k] = pb_y[k] * sdk0_6[k]
                  + f_0 * sdi_3[k]
                  - f_12 * pc_y[k] * sdk1_6[k];

        t_43[k] = f_3 * pc_z[k] * sfi_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pc_y, pc_z, sdk0_9, sdk0_10, sdi_5, \
                         sdi_6, sdk1_9, sdk1_10, sfi_33, sfi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_13 * sdi_5[k]
                  + f_3 * pc_y[k] * sfi_33[k];

        t_45[k] = pb_y[k] * sdk0_9[k]
                  - f_12 * pc_y[k] * sdk1_9[k];

        t_46[k] = pb_y[k] * sdk0_10[k]
                  + f_15 * sdi_6[k]
                  - f_12 * pc_y[k] * sdk1_10[k];

        t_47[k] = f_3 * pc_z[k] * sfi_34[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_y, pc_y, sdk0_12, sdk0_14, sdk0_15, sdi_8, \
                         sdi_9, sdi_10, sdk1_12, sdk1_14, sdk1_15, \
                         sfi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pb_y[k] * sdk0_12[k]
                  + f_14 * sdi_8[k]
                  - f_12 * pc_y[k] * sdk1_12[k];

        t_49[k] = f_13 * sdi_9[k]
                  + f_3 * pc_y[k] * sfi_37[k];

        t_50[k] = pb_y[k] * sdk0_14[k]
                  - f_12 * pc_y[k] * sdk1_14[k];

        t_51[k] = pb_y[k] * sdk0_15[k]
                  + f_16 * sdi_10[k]
                  - f_12 * pc_y[k] * sdk1_15[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_y, pc_y, pc_z, sdk0_17, sdk0_18, sdi_12, \
                         sdi_13, sdi_14, sdk1_17, sdk1_18, sfi_38, \
                         sfi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * pc_z[k] * sfi_38[k];

        t_53[k] = pb_y[k] * sdk0_17[k]
                  + f_0 * sdi_12[k]
                  - f_12 * pc_y[k] * sdk1_17[k];

        t_54[k] = pb_y[k] * sdk0_18[k]
                  + f_14 * sdi_13[k]
                  - f_12 * pc_y[k] * sdk1_18[k];

        t_55[k] = f_13 * sdi_14[k]
                  + f_3 * pc_y[k] * sfi_42[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pc_x, pc_y, sdk0_20, sdi_49, sdi_50, \
                         sdi_51, sdk1_20, sfi_49, sfi_50, sfi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_y[k] * sdk0_20[k]
                  - f_12 * pc_y[k] * sdk1_20[k];

        t_57[k] = f_14 * sdi_49[k]
                  + f_3 * pc_x[k] * sfi_49[k];

        t_58[k] = f_14 * sdi_50[k]
                  + f_3 * pc_x[k] * sfi_50[k];

        t_59[k] = f_14 * sdi_51[k]
                  + f_3 * pc_x[k] * sfi_51[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, sdi_52, sdi_53, sdi_54, sdi_55, sfi_52, \
                         sfi_53, sfi_54, sfi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_14 * sdi_52[k]
                  + f_3 * pc_x[k] * sfi_52[k];

        t_61[k] = f_14 * sdi_53[k]
                  + f_3 * pc_x[k] * sfi_53[k];

        t_62[k] = f_14 * sdi_54[k]
                  + f_3 * pc_x[k] * sfi_54[k];

        t_63[k] = f_14 * sdi_55[k]
                  + f_3 * pc_x[k] * sfi_55[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pc_y, pc_z, sdi_21, sdi_23, sfh0_36, sfh0_38, \
                         sfh1_36, sfh1_38, sfi_49, sfi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_13 * sdi_21[k]
                  + f_1 * sfh0_36[k]
                  - f_2 * sfh1_36[k]
                  + f_3 * pc_y[k] * sfi_49[k];

        t_65[k] = f_3 * pc_z[k] * sfi_49[k];

        t_66[k] = f_13 * sdi_23[k]
                  + f_4 * sfh0_38[k]
                  - f_5 * sfh1_38[k]
                  + f_3 * pc_y[k] * sfi_51[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_y, sdi_24, sdi_25, sdi_26, sfh0_39, sfh0_40, \
                         sfh0_41, sfh1_39, sfh1_40, sfh1_41, sfi_52, sfi_53, \
                         sfi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_13 * sdi_24[k]
                  + f_6 * sfh0_39[k]
                  - f_7 * sfh1_39[k]
                  + f_3 * pc_y[k] * sfi_52[k];

        t_68[k] = f_13 * sdi_25[k]
                  + f_8 * sfh0_40[k]
                  - f_9 * sfh1_40[k]
                  + f_3 * pc_y[k] * sfi_53[k];

        t_69[k] = f_13 * sdi_26[k]
                  + f_10 * sfh0_41[k]
                  - f_11 * sfh1_41[k]
                  + f_3 * pc_y[k] * sfi_54[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_y, pb_z, pc_y, pc_z, sdk0_0, sdk0_35, \
                         sdi_27, sdk1_0, sdk1_35, sfi_55, sfi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_13 * sdi_27[k]
                  + f_3 * pc_y[k] * sfi_55[k];

        t_71[k] = pb_y[k] * sdk0_35[k]
                  - f_12 * pc_y[k] * sdk1_35[k];

        t_72[k] = pb_z[k] * sdk0_0[k]
                  - f_12 * pc_z[k] * sdk1_0[k];

        t_73[k] = f_3 * pc_y[k] * sfi_56[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_z, pc_y, pc_z, sdk0_3, sdk0_5, sdi_0, \
                         sdi_2, sdk1_3, sdk1_5, sfi_56, sfi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * sdi_0[k]
                  + f_3 * pc_z[k] * sfi_56[k];

        t_75[k] = pb_z[k] * sdk0_3[k]
                  - f_12 * pc_z[k] * sdk1_3[k];

        t_76[k] = f_3 * pc_y[k] * sfi_58[k];

        t_77[k] = pb_z[k] * sdk0_5[k]
                  + f_14 * sdi_2[k]
                  - f_12 * pc_z[k] * sdk1_5[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_z, pc_y, pc_z, sdk0_6, sdk0_9, sdi_3, \
                         sdi_5, sdk1_6, sdk1_9, sfi_59, sfi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pb_z[k] * sdk0_6[k]
                  - f_12 * pc_z[k] * sdk1_6[k];

        t_79[k] = f_13 * sdi_3[k]
                  + f_3 * pc_z[k] * sfi_59[k];

        t_80[k] = f_3 * pc_y[k] * sfi_61[k];

        t_81[k] = pb_z[k] * sdk0_9[k]
                  + f_0 * sdi_5[k]
                  - f_12 * pc_z[k] * sdk1_9[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pb_z, pc_y, pc_z, sdk0_10, sdk0_12, sdi_6, \
                         sdi_7, sdk1_10, sdk1_12, sfi_62, sfi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pb_z[k] * sdk0_10[k]
                  - f_12 * pc_z[k] * sdk1_10[k];

        t_83[k] = f_13 * sdi_6[k]
                  + f_3 * pc_z[k] * sfi_62[k];

        t_84[k] = pb_z[k] * sdk0_12[k]
                  + f_14 * sdi_7[k]
                  - f_12 * pc_z[k] * sdk1_12[k];

        t_85[k] = f_3 * pc_y[k] * sfi_65[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_z, pc_z, sdk0_14, sdk0_15, sdk0_17, sdi_9, \
                         sdi_10, sdi_11, sdk1_14, sdk1_15, sdk1_17, \
                         sfi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_z[k] * sdk0_14[k]
                  + f_15 * sdi_9[k]
                  - f_12 * pc_z[k] * sdk1_14[k];

        t_87[k] = pb_z[k] * sdk0_15[k]
                  - f_12 * pc_z[k] * sdk1_15[k];

        t_88[k] = f_13 * sdi_10[k]
                  + f_3 * pc_z[k] * sfi_66[k];

        t_89[k] = pb_z[k] * sdk0_17[k]
                  + f_14 * sdi_11[k]
                  - f_12 * pc_z[k] * sdk1_17[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_z, pc_y, pc_z, sdk0_18, sdk0_20, sdi_12, sdi_14, \
                         sdk1_18, sdk1_20, sfi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_z[k] * sdk0_18[k]
                  + f_0 * sdi_12[k]
                  - f_12 * pc_z[k] * sdk1_18[k];

        t_91[k] = f_3 * pc_y[k] * sfi_70[k];

        t_92[k] = pb_z[k] * sdk0_20[k]
                  + f_16 * sdi_14[k]
                  - f_12 * pc_z[k] * sdk1_20[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pc_x, sdi_77, sdi_78, sdi_79, sdi_80, \
                         sdi_81, sfi_77, sfi_78, sfi_79, sfi_80, \
                         sfi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_14 * sdi_77[k]
                  + f_3 * pc_x[k] * sfi_77[k];

        t_94[k] = f_14 * sdi_78[k]
                  + f_3 * pc_x[k] * sfi_78[k];

        t_95[k] = f_14 * sdi_79[k]
                  + f_3 * pc_x[k] * sfi_79[k];

        t_96[k] = f_14 * sdi_80[k]
                  + f_3 * pc_x[k] * sfi_80[k];

        t_97[k] = f_14 * sdi_81[k]
                  + f_3 * pc_x[k] * sfi_81[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_z, pc_x, pc_z, sdk0_28, sdi_21, sdi_82, \
                         sdi_83, sdk1_28, sfi_77, sfi_82, sfi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_14 * sdi_82[k]
                  + f_3 * pc_x[k] * sfi_82[k];

        t_99[k] = f_14 * sdi_83[k]
                  + f_3 * pc_x[k] * sfi_83[k];

        t_100[k] = pb_z[k] * sdk0_28[k]
                   - f_12 * pc_z[k] * sdk1_28[k];

        t_101[k] = f_13 * sdi_21[k]
                   + f_3 * pc_z[k] * sfi_77[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pc_y, sfh0_59, sfh0_60, sfh0_61, sfh1_59, \
                         sfh1_60, sfh1_61, sfi_79, sfi_80, sfi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_4 * sfh0_59[k]
                   - f_5 * sfh1_59[k]
                   + f_3 * pc_y[k] * sfi_79[k];

        t_103[k] = f_6 * sfh0_60[k]
                   - f_7 * sfh1_60[k]
                   + f_3 * pc_y[k] * sfi_80[k];

        t_104[k] = f_8 * sfh0_61[k]
                   - f_9 * sfh1_61[k]
                   + f_3 * pc_y[k] * sfi_81[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pb_x, pc_x, pc_y, pc_z, sdk0_108, sdi_27, \
                         sdi_84, sdk1_108, sfh0_62, sfh1_62, sfi_82, \
                         sfi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_10 * sfh0_62[k]
                   - f_11 * sfh1_62[k]
                   + f_3 * pc_y[k] * sfi_82[k];

        t_106[k] = f_3 * pc_y[k] * sfi_83[k];

        t_107[k] = f_13 * sdi_27[k]
                   + f_1 * sfh0_62[k]
                   - f_2 * sfh1_62[k]
                   + f_3 * pc_z[k] * sfi_83[k];

        t_108[k] = pb_x[k] * sdk0_108[k]
                   + f_17 * sdi_84[k]
                   - f_12 * pc_x[k] * sdk1_108[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pb_x, pc_x, pc_y, pc_z, sdk0_111, sdi_28, \
                         sdi_30, sdi_87, sdk1_111, sfi_84, sfi_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_14 * sdi_28[k]
                   + f_3 * pc_y[k] * sfi_84[k];

        t_110[k] = f_3 * pc_z[k] * sfi_84[k];

        t_111[k] = pb_x[k] * sdk0_111[k]
                   + f_16 * sdi_87[k]
                   - f_12 * pc_x[k] * sdk1_111[k];

        t_112[k] = f_14 * sdi_30[k]
                   + f_3 * pc_y[k] * sfi_86[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pb_x, pc_x, pc_z, sdk0_113, sdk0_114, sdi_89, \
                         sdi_90, sdk1_113, sdk1_114, sfi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pb_x[k] * sdk0_113[k]
                   + f_16 * sdi_89[k]
                   - f_12 * pc_x[k] * sdk1_113[k];

        t_114[k] = pb_x[k] * sdk0_114[k]
                   + f_15 * sdi_90[k]
                   - f_12 * pc_x[k] * sdk1_114[k];

        t_115[k] = f_3 * pc_z[k] * sfi_87[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pb_x, pc_x, pc_y, sdk0_117, sdk0_118, sdi_33, \
                         sdi_93, sdi_94, sdk1_117, sdk1_118, sfi_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_14 * sdi_33[k]
                   + f_3 * pc_y[k] * sfi_89[k];

        t_117[k] = pb_x[k] * sdk0_117[k]
                   + f_15 * sdi_93[k]
                   - f_12 * pc_x[k] * sdk1_117[k];

        t_118[k] = pb_x[k] * sdk0_118[k]
                   + f_0 * sdi_94[k]
                   - f_12 * pc_x[k] * sdk1_118[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pb_x, pc_x, pc_y, pc_z, sdk0_120, sdi_37, \
                         sdi_96, sdk1_120, sfi_90, sfi_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_3 * pc_z[k] * sfi_90[k];

        t_120[k] = pb_x[k] * sdk0_120[k]
                   + f_0 * sdi_96[k]
                   - f_12 * pc_x[k] * sdk1_120[k];

        t_121[k] = f_14 * sdi_37[k]
                   + f_3 * pc_y[k] * sfi_93[k];
    }
}

static auto
compute_prim_sfk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sdk0,
                                                          const size_t sdi, const size_t sdk1,
                                                          const size_t sfh0, const size_t sfh1,
                                                          const size_t sfi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 2.0 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 3.5 / q;

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
    auto *t_239 = buffer.data(target + 239);
    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);
    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdk0_39 = buffer.data(sdk0 + 39);
    const auto *sdk0_42 = buffer.data(sdk0 + 42);
    const auto *sdk0_46 = buffer.data(sdk0 + 46);
    const auto *sdk0_51 = buffer.data(sdk0 + 51);
    const auto *sdk0_72 = buffer.data(sdk0 + 72);
    const auto *sdk0_77 = buffer.data(sdk0 + 77);
    const auto *sdk0_81 = buffer.data(sdk0 + 81);
    const auto *sdk0_86 = buffer.data(sdk0 + 86);
    const auto *sdk0_92 = buffer.data(sdk0 + 92);
    const auto *sdk0_122 = buffer.data(sdk0 + 122);
    const auto *sdk0_123 = buffer.data(sdk0 + 123);
    const auto *sdk0_125 = buffer.data(sdk0 + 125);
    const auto *sdk0_126 = buffer.data(sdk0 + 126);
    const auto *sdk0_128 = buffer.data(sdk0 + 128);
    const auto *sdk0_136 = buffer.data(sdk0 + 136);
    const auto *sdk0_138 = buffer.data(sdk0 + 138);
    const auto *sdk0_139 = buffer.data(sdk0 + 139);
    const auto *sdk0_140 = buffer.data(sdk0 + 140);
    const auto *sdk0_141 = buffer.data(sdk0 + 141);
    const auto *sdk0_143 = buffer.data(sdk0 + 143);
    const auto *sdk0_156 = buffer.data(sdk0 + 156);
    const auto *sdk0_161 = buffer.data(sdk0 + 161);
    const auto *sdk0_162 = buffer.data(sdk0 + 162);
    const auto *sdk0_172 = buffer.data(sdk0 + 172);
    const auto *sdk0_174 = buffer.data(sdk0 + 174);
    const auto *sdk0_175 = buffer.data(sdk0 + 175);
    const auto *sdk0_176 = buffer.data(sdk0 + 176);
    const auto *sdk0_177 = buffer.data(sdk0 + 177);
    const auto *sdk0_179 = buffer.data(sdk0 + 179);
    const auto *sdk0_180 = buffer.data(sdk0 + 180);
    const auto *sdk0_183 = buffer.data(sdk0 + 183);
    const auto *sdk0_185 = buffer.data(sdk0 + 185);
    const auto *sdk0_186 = buffer.data(sdk0 + 186);
    const auto *sdk0_189 = buffer.data(sdk0 + 189);
    const auto *sdk0_190 = buffer.data(sdk0 + 190);
    const auto *sdk0_192 = buffer.data(sdk0 + 192);
    const auto *sdk0_194 = buffer.data(sdk0 + 194);
    const auto *sdk0_195 = buffer.data(sdk0 + 195);
    const auto *sdk0_197 = buffer.data(sdk0 + 197);
    const auto *sdk0_198 = buffer.data(sdk0 + 198);
    const auto *sdk0_200 = buffer.data(sdk0 + 200);
    const auto *sdk0_208 = buffer.data(sdk0 + 208);
    const auto *sdk0_210 = buffer.data(sdk0 + 210);
    const auto *sdk0_211 = buffer.data(sdk0 + 211);
    const auto *sdk0_212 = buffer.data(sdk0 + 212);
    const auto *sdk0_213 = buffer.data(sdk0 + 213);
    const auto *sdk0_215 = buffer.data(sdk0 + 215);

    const auto *sdi_28 = buffer.data(sdi + 28);
    const auto *sdi_31 = buffer.data(sdi + 31);
    const auto *sdi_34 = buffer.data(sdi + 34);
    const auto *sdi_38 = buffer.data(sdi + 38);
    const auto *sdi_42 = buffer.data(sdi + 42);
    const auto *sdi_49 = buffer.data(sdi + 49);
    const auto *sdi_55 = buffer.data(sdi + 55);
    const auto *sdi_56 = buffer.data(sdi + 56);
    const auto *sdi_58 = buffer.data(sdi + 58);
    const auto *sdi_59 = buffer.data(sdi + 59);
    const auto *sdi_61 = buffer.data(sdi + 61);
    const auto *sdi_62 = buffer.data(sdi + 62);
    const auto *sdi_65 = buffer.data(sdi + 65);
    const auto *sdi_66 = buffer.data(sdi + 66);
    const auto *sdi_70 = buffer.data(sdi + 70);
    const auto *sdi_77 = buffer.data(sdi + 77);
    const auto *sdi_83 = buffer.data(sdi + 83);
    const auto *sdi_84 = buffer.data(sdi + 84);
    const auto *sdi_86 = buffer.data(sdi + 86);
    const auto *sdi_89 = buffer.data(sdi + 89);
    const auto *sdi_93 = buffer.data(sdi + 93);
    const auto *sdi_98 = buffer.data(sdi + 98);
    const auto *sdi_99 = buffer.data(sdi + 99);
    const auto *sdi_101 = buffer.data(sdi + 101);
    const auto *sdi_102 = buffer.data(sdi + 102);
    const auto *sdi_104 = buffer.data(sdi + 104);
    const auto *sdi_105 = buffer.data(sdi + 105);
    const auto *sdi_106 = buffer.data(sdi + 106);
    const auto *sdi_107 = buffer.data(sdi + 107);
    const auto *sdi_108 = buffer.data(sdi + 108);
    const auto *sdi_109 = buffer.data(sdi + 109);
    const auto *sdi_110 = buffer.data(sdi + 110);
    const auto *sdi_111 = buffer.data(sdi + 111);
    const auto *sdi_124 = buffer.data(sdi + 124);
    const auto *sdi_129 = buffer.data(sdi + 129);
    const auto *sdi_130 = buffer.data(sdi + 130);
    const auto *sdi_133 = buffer.data(sdi + 133);
    const auto *sdi_134 = buffer.data(sdi + 134);
    const auto *sdi_135 = buffer.data(sdi + 135);
    const auto *sdi_136 = buffer.data(sdi + 136);
    const auto *sdi_137 = buffer.data(sdi + 137);
    const auto *sdi_138 = buffer.data(sdi + 138);
    const auto *sdi_139 = buffer.data(sdi + 139);
    const auto *sdi_140 = buffer.data(sdi + 140);
    const auto *sdi_143 = buffer.data(sdi + 143);
    const auto *sdi_145 = buffer.data(sdi + 145);
    const auto *sdi_146 = buffer.data(sdi + 146);
    const auto *sdi_149 = buffer.data(sdi + 149);
    const auto *sdi_150 = buffer.data(sdi + 150);
    const auto *sdi_152 = buffer.data(sdi + 152);
    const auto *sdi_154 = buffer.data(sdi + 154);
    const auto *sdi_155 = buffer.data(sdi + 155);
    const auto *sdi_157 = buffer.data(sdi + 157);
    const auto *sdi_158 = buffer.data(sdi + 158);
    const auto *sdi_160 = buffer.data(sdi + 160);
    const auto *sdi_161 = buffer.data(sdi + 161);
    const auto *sdi_162 = buffer.data(sdi + 162);
    const auto *sdi_163 = buffer.data(sdi + 163);
    const auto *sdi_164 = buffer.data(sdi + 164);
    const auto *sdi_165 = buffer.data(sdi + 165);
    const auto *sdi_166 = buffer.data(sdi + 166);
    const auto *sdi_167 = buffer.data(sdi + 167);

    const auto *sdk1_39 = buffer.data(sdk1 + 39);
    const auto *sdk1_42 = buffer.data(sdk1 + 42);
    const auto *sdk1_46 = buffer.data(sdk1 + 46);
    const auto *sdk1_51 = buffer.data(sdk1 + 51);
    const auto *sdk1_72 = buffer.data(sdk1 + 72);
    const auto *sdk1_77 = buffer.data(sdk1 + 77);
    const auto *sdk1_81 = buffer.data(sdk1 + 81);
    const auto *sdk1_86 = buffer.data(sdk1 + 86);
    const auto *sdk1_92 = buffer.data(sdk1 + 92);
    const auto *sdk1_122 = buffer.data(sdk1 + 122);
    const auto *sdk1_123 = buffer.data(sdk1 + 123);
    const auto *sdk1_125 = buffer.data(sdk1 + 125);
    const auto *sdk1_126 = buffer.data(sdk1 + 126);
    const auto *sdk1_128 = buffer.data(sdk1 + 128);
    const auto *sdk1_136 = buffer.data(sdk1 + 136);
    const auto *sdk1_138 = buffer.data(sdk1 + 138);
    const auto *sdk1_139 = buffer.data(sdk1 + 139);
    const auto *sdk1_140 = buffer.data(sdk1 + 140);
    const auto *sdk1_141 = buffer.data(sdk1 + 141);
    const auto *sdk1_143 = buffer.data(sdk1 + 143);
    const auto *sdk1_156 = buffer.data(sdk1 + 156);
    const auto *sdk1_161 = buffer.data(sdk1 + 161);
    const auto *sdk1_162 = buffer.data(sdk1 + 162);
    const auto *sdk1_172 = buffer.data(sdk1 + 172);
    const auto *sdk1_174 = buffer.data(sdk1 + 174);
    const auto *sdk1_175 = buffer.data(sdk1 + 175);
    const auto *sdk1_176 = buffer.data(sdk1 + 176);
    const auto *sdk1_177 = buffer.data(sdk1 + 177);
    const auto *sdk1_179 = buffer.data(sdk1 + 179);
    const auto *sdk1_180 = buffer.data(sdk1 + 180);
    const auto *sdk1_183 = buffer.data(sdk1 + 183);
    const auto *sdk1_185 = buffer.data(sdk1 + 185);
    const auto *sdk1_186 = buffer.data(sdk1 + 186);
    const auto *sdk1_189 = buffer.data(sdk1 + 189);
    const auto *sdk1_190 = buffer.data(sdk1 + 190);
    const auto *sdk1_192 = buffer.data(sdk1 + 192);
    const auto *sdk1_194 = buffer.data(sdk1 + 194);
    const auto *sdk1_195 = buffer.data(sdk1 + 195);
    const auto *sdk1_197 = buffer.data(sdk1 + 197);
    const auto *sdk1_198 = buffer.data(sdk1 + 198);
    const auto *sdk1_200 = buffer.data(sdk1 + 200);
    const auto *sdk1_208 = buffer.data(sdk1 + 208);
    const auto *sdk1_210 = buffer.data(sdk1 + 210);
    const auto *sdk1_211 = buffer.data(sdk1 + 211);
    const auto *sdk1_212 = buffer.data(sdk1 + 212);
    const auto *sdk1_213 = buffer.data(sdk1 + 213);
    const auto *sdk1_215 = buffer.data(sdk1 + 215);

    const auto *sfh0_126 = buffer.data(sfh0 + 126);
    const auto *sfh0_129 = buffer.data(sfh0 + 129);
    const auto *sfh0_131 = buffer.data(sfh0 + 131);
    const auto *sfh0_132 = buffer.data(sfh0 + 132);
    const auto *sfh0_135 = buffer.data(sfh0 + 135);
    const auto *sfh0_136 = buffer.data(sfh0 + 136);
    const auto *sfh0_138 = buffer.data(sfh0 + 138);
    const auto *sfh0_140 = buffer.data(sfh0 + 140);
    const auto *sfh0_141 = buffer.data(sfh0 + 141);
    const auto *sfh0_143 = buffer.data(sfh0 + 143);
    const auto *sfh0_144 = buffer.data(sfh0 + 144);
    const auto *sfh0_146 = buffer.data(sfh0 + 146);

    const auto *sfh1_126 = buffer.data(sfh1 + 126);
    const auto *sfh1_129 = buffer.data(sfh1 + 129);
    const auto *sfh1_131 = buffer.data(sfh1 + 131);
    const auto *sfh1_132 = buffer.data(sfh1 + 132);
    const auto *sfh1_135 = buffer.data(sfh1 + 135);
    const auto *sfh1_136 = buffer.data(sfh1 + 136);
    const auto *sfh1_138 = buffer.data(sfh1 + 138);
    const auto *sfh1_140 = buffer.data(sfh1 + 140);
    const auto *sfh1_141 = buffer.data(sfh1 + 141);
    const auto *sfh1_143 = buffer.data(sfh1 + 143);
    const auto *sfh1_144 = buffer.data(sfh1 + 144);
    const auto *sfh1_146 = buffer.data(sfh1 + 146);

    const auto *sfi_94 = buffer.data(sfi + 94);
    const auto *sfi_98 = buffer.data(sfi + 98);
    const auto *sfi_105 = buffer.data(sfi + 105);
    const auto *sfi_106 = buffer.data(sfi + 106);
    const auto *sfi_107 = buffer.data(sfi + 107);
    const auto *sfi_108 = buffer.data(sfi + 108);
    const auto *sfi_109 = buffer.data(sfi + 109);
    const auto *sfi_110 = buffer.data(sfi + 110);
    const auto *sfi_111 = buffer.data(sfi + 111);
    const auto *sfi_112 = buffer.data(sfi + 112);
    const auto *sfi_114 = buffer.data(sfi + 114);
    const auto *sfi_115 = buffer.data(sfi + 115);
    const auto *sfi_117 = buffer.data(sfi + 117);
    const auto *sfi_118 = buffer.data(sfi + 118);
    const auto *sfi_121 = buffer.data(sfi + 121);
    const auto *sfi_122 = buffer.data(sfi + 122);
    const auto *sfi_126 = buffer.data(sfi + 126);
    const auto *sfi_133 = buffer.data(sfi + 133);
    const auto *sfi_134 = buffer.data(sfi + 134);
    const auto *sfi_135 = buffer.data(sfi + 135);
    const auto *sfi_136 = buffer.data(sfi + 136);
    const auto *sfi_137 = buffer.data(sfi + 137);
    const auto *sfi_138 = buffer.data(sfi + 138);
    const auto *sfi_139 = buffer.data(sfi + 139);
    const auto *sfi_140 = buffer.data(sfi + 140);
    const auto *sfi_142 = buffer.data(sfi + 142);
    const auto *sfi_143 = buffer.data(sfi + 143);
    const auto *sfi_145 = buffer.data(sfi + 145);
    const auto *sfi_146 = buffer.data(sfi + 146);
    const auto *sfi_149 = buffer.data(sfi + 149);
    const auto *sfi_150 = buffer.data(sfi + 150);
    const auto *sfi_154 = buffer.data(sfi + 154);
    const auto *sfi_161 = buffer.data(sfi + 161);
    const auto *sfi_162 = buffer.data(sfi + 162);
    const auto *sfi_163 = buffer.data(sfi + 163);
    const auto *sfi_164 = buffer.data(sfi + 164);
    const auto *sfi_165 = buffer.data(sfi + 165);
    const auto *sfi_166 = buffer.data(sfi + 166);
    const auto *sfi_167 = buffer.data(sfi + 167);
    const auto *sfi_168 = buffer.data(sfi + 168);
    const auto *sfi_170 = buffer.data(sfi + 170);
    const auto *sfi_171 = buffer.data(sfi + 171);
    const auto *sfi_173 = buffer.data(sfi + 173);
    const auto *sfi_174 = buffer.data(sfi + 174);
    const auto *sfi_177 = buffer.data(sfi + 177);
    const auto *sfi_178 = buffer.data(sfi + 178);
    const auto *sfi_180 = buffer.data(sfi + 180);
    const auto *sfi_182 = buffer.data(sfi + 182);
    const auto *sfi_183 = buffer.data(sfi + 183);
    const auto *sfi_185 = buffer.data(sfi + 185);
    const auto *sfi_186 = buffer.data(sfi + 186);
    const auto *sfi_188 = buffer.data(sfi + 188);
    const auto *sfi_189 = buffer.data(sfi + 189);
    const auto *sfi_190 = buffer.data(sfi + 190);
    const auto *sfi_191 = buffer.data(sfi + 191);
    const auto *sfi_192 = buffer.data(sfi + 192);
    const auto *sfi_193 = buffer.data(sfi + 193);
    const auto *sfi_194 = buffer.data(sfi + 194);
    const auto *sfi_195 = buffer.data(sfi + 195);

#pragma omp simd aligned(t_122, t_123, t_124, pb_x, pc_x, pc_z, sdk0_122, sdk0_123, sdi_98, \
                         sdi_99, sdk1_122, sdk1_123, sfi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = pb_x[k] * sdk0_122[k]
                   + f_0 * sdi_98[k]
                   - f_12 * pc_x[k] * sdk1_122[k];

        t_123[k] = pb_x[k] * sdk0_123[k]
                   + f_14 * sdi_99[k]
                   - f_12 * pc_x[k] * sdk1_123[k];

        t_124[k] = f_3 * pc_z[k] * sfi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pb_x, pc_x, pc_y, sdk0_125, sdk0_126, sdi_42, \
                         sdi_101, sdi_102, sdk1_125, sdk1_126, sfi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_x[k] * sdk0_125[k]
                   + f_14 * sdi_101[k]
                   - f_12 * pc_x[k] * sdk1_125[k];

        t_126[k] = pb_x[k] * sdk0_126[k]
                   + f_14 * sdi_102[k]
                   - f_12 * pc_x[k] * sdk1_126[k];

        t_127[k] = f_14 * sdi_42[k]
                   + f_3 * pc_y[k] * sfi_98[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pb_x, pc_x, sdk0_128, sdi_104, sdi_105, \
                         sdi_106, sdi_107, sdk1_128, sfi_105, sfi_106, \
                         sfi_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = pb_x[k] * sdk0_128[k]
                   + f_14 * sdi_104[k]
                   - f_12 * pc_x[k] * sdk1_128[k];

        t_129[k] = f_13 * sdi_105[k]
                   + f_3 * pc_x[k] * sfi_105[k];

        t_130[k] = f_13 * sdi_106[k]
                   + f_3 * pc_x[k] * sfi_106[k];

        t_131[k] = f_13 * sdi_107[k]
                   + f_3 * pc_x[k] * sfi_107[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, sdi_108, sdi_109, sdi_110, sdi_111, \
                         sfi_108, sfi_109, sfi_110, sfi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_13 * sdi_108[k]
                   + f_3 * pc_x[k] * sfi_108[k];

        t_133[k] = f_13 * sdi_109[k]
                   + f_3 * pc_x[k] * sfi_109[k];

        t_134[k] = f_13 * sdi_110[k]
                   + f_3 * pc_x[k] * sfi_110[k];

        t_135[k] = f_13 * sdi_111[k]
                   + f_3 * pc_x[k] * sfi_111[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pb_x, pc_x, pc_z, sdk0_136, sdk0_138, \
                         sdk0_139, sdk1_136, sdk1_138, sdk1_139, \
                         sfi_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_x[k] * sdk0_136[k]
                   - f_12 * pc_x[k] * sdk1_136[k];

        t_137[k] = f_3 * pc_z[k] * sfi_105[k];

        t_138[k] = pb_x[k] * sdk0_138[k]
                   - f_12 * pc_x[k] * sdk1_138[k];

        t_139[k] = pb_x[k] * sdk0_139[k]
                   - f_12 * pc_x[k] * sdk1_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pb_x, pc_x, pc_y, sdk0_140, sdk0_141, \
                         sdk0_143, sdi_55, sdk1_140, sdk1_141, sdk1_143, \
                         sfi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = pb_x[k] * sdk0_140[k]
                   - f_12 * pc_x[k] * sdk1_140[k];

        t_141[k] = pb_x[k] * sdk0_141[k]
                   - f_12 * pc_x[k] * sdk1_141[k];

        t_142[k] = f_14 * sdi_55[k]
                   + f_3 * pc_y[k] * sfi_111[k];

        t_143[k] = pb_x[k] * sdk0_143[k]
                   - f_12 * pc_x[k] * sdk1_143[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_y, pb_z, pc_y, pc_z, sdk0_39, sdk0_72, \
                         sdi_28, sdi_56, sdk1_39, sdk1_72, sfi_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pb_y[k] * sdk0_72[k]
                   - f_12 * pc_y[k] * sdk1_72[k];

        t_145[k] = f_13 * sdi_56[k]
                   + f_3 * pc_y[k] * sfi_112[k];

        t_146[k] = f_13 * sdi_28[k]
                   + f_3 * pc_z[k] * sfi_112[k];

        t_147[k] = pb_z[k] * sdk0_39[k]
                   - f_12 * pc_z[k] * sdk1_39[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_y, pb_z, pc_y, pc_z, sdk0_42, sdk0_77, \
                         sdi_31, sdi_58, sdk1_42, sdk1_77, sfi_114, \
                         sfi_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_13 * sdi_58[k]
                   + f_3 * pc_y[k] * sfi_114[k];

        t_149[k] = pb_y[k] * sdk0_77[k]
                   - f_12 * pc_y[k] * sdk1_77[k];

        t_150[k] = pb_z[k] * sdk0_42[k]
                   - f_12 * pc_z[k] * sdk1_42[k];

        t_151[k] = f_13 * sdi_31[k]
                   + f_3 * pc_z[k] * sfi_115[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_y, pb_z, pc_y, pc_z, sdk0_46, sdk0_81, \
                         sdi_34, sdi_61, sdk1_46, sdk1_81, sfi_117, \
                         sfi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_13 * sdi_61[k]
                   + f_3 * pc_y[k] * sfi_117[k];

        t_153[k] = pb_y[k] * sdk0_81[k]
                   - f_12 * pc_y[k] * sdk1_81[k];

        t_154[k] = pb_z[k] * sdk0_46[k]
                   - f_12 * pc_z[k] * sdk1_46[k];

        t_155[k] = f_13 * sdi_34[k]
                   + f_3 * pc_z[k] * sfi_118[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pb_x, pb_y, pc_x, pc_y, sdk0_86, sdk0_156, \
                         sdi_65, sdi_124, sdk1_86, sdk1_156, sfi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pb_x[k] * sdk0_156[k]
                   + f_0 * sdi_124[k]
                   - f_12 * pc_x[k] * sdk1_156[k];

        t_157[k] = f_13 * sdi_65[k]
                   + f_3 * pc_y[k] * sfi_121[k];

        t_158[k] = pb_y[k] * sdk0_86[k]
                   - f_12 * pc_y[k] * sdk1_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pb_x, pb_z, pc_x, pc_z, sdk0_51, sdk0_161, \
                         sdi_38, sdi_129, sdk1_51, sdk1_161, sfi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pb_z[k] * sdk0_51[k]
                   - f_12 * pc_z[k] * sdk1_51[k];

        t_160[k] = f_13 * sdi_38[k]
                   + f_3 * pc_z[k] * sfi_122[k];

        t_161[k] = pb_x[k] * sdk0_161[k]
                   + f_14 * sdi_129[k]
                   - f_12 * pc_x[k] * sdk1_161[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pb_x, pb_y, pc_x, pc_y, sdk0_92, sdk0_162, \
                         sdi_70, sdi_130, sdk1_92, sdk1_162, sfi_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pb_x[k] * sdk0_162[k]
                   + f_14 * sdi_130[k]
                   - f_12 * pc_x[k] * sdk1_162[k];

        t_163[k] = f_13 * sdi_70[k]
                   + f_3 * pc_y[k] * sfi_126[k];

        t_164[k] = pb_y[k] * sdk0_92[k]
                   - f_12 * pc_y[k] * sdk1_92[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pc_x, sdi_133, sdi_134, sdi_135, \
                         sdi_136, sdi_137, sfi_133, sfi_134, sfi_135, sfi_136, \
                         sfi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_13 * sdi_133[k]
                   + f_3 * pc_x[k] * sfi_133[k];

        t_166[k] = f_13 * sdi_134[k]
                   + f_3 * pc_x[k] * sfi_134[k];

        t_167[k] = f_13 * sdi_135[k]
                   + f_3 * pc_x[k] * sfi_135[k];

        t_168[k] = f_13 * sdi_136[k]
                   + f_3 * pc_x[k] * sfi_136[k];

        t_169[k] = f_13 * sdi_137[k]
                   + f_3 * pc_x[k] * sfi_137[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_x, pc_x, pc_z, sdk0_172, sdi_49, \
                         sdi_138, sdi_139, sdk1_172, sfi_133, sfi_138, \
                         sfi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_13 * sdi_138[k]
                   + f_3 * pc_x[k] * sfi_138[k];

        t_171[k] = f_13 * sdi_139[k]
                   + f_3 * pc_x[k] * sfi_139[k];

        t_172[k] = pb_x[k] * sdk0_172[k]
                   - f_12 * pc_x[k] * sdk1_172[k];

        t_173[k] = f_13 * sdi_49[k]
                   + f_3 * pc_z[k] * sfi_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pb_x, pc_x, sdk0_174, sdk0_175, sdk0_176, \
                         sdk0_177, sdk1_174, sdk1_175, sdk1_176, \
                         sdk1_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pb_x[k] * sdk0_174[k]
                   - f_12 * pc_x[k] * sdk1_174[k];

        t_175[k] = pb_x[k] * sdk0_175[k]
                   - f_12 * pc_x[k] * sdk1_175[k];

        t_176[k] = pb_x[k] * sdk0_176[k]
                   - f_12 * pc_x[k] * sdk1_176[k];

        t_177[k] = pb_x[k] * sdk0_177[k]
                   - f_12 * pc_x[k] * sdk1_177[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, pb_x, pc_x, pc_y, sdk0_179, sdk0_180, \
                         sdi_83, sdi_140, sdk1_179, sdk1_180, sfi_139, \
                         sfi_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_13 * sdi_83[k]
                   + f_3 * pc_y[k] * sfi_139[k];

        t_179[k] = pb_x[k] * sdk0_179[k]
                   - f_12 * pc_x[k] * sdk1_179[k];

        t_180[k] = pb_x[k] * sdk0_180[k]
                   + f_17 * sdi_140[k]
                   - f_12 * pc_x[k] * sdk1_180[k];

        t_181[k] = f_3 * pc_y[k] * sfi_140[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, pb_x, pc_x, pc_y, pc_z, sdk0_183, sdi_56, \
                         sdi_143, sdk1_183, sfi_140, sfi_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_14 * sdi_56[k]
                   + f_3 * pc_z[k] * sfi_140[k];

        t_183[k] = pb_x[k] * sdk0_183[k]
                   + f_16 * sdi_143[k]
                   - f_12 * pc_x[k] * sdk1_183[k];

        t_184[k] = f_3 * pc_y[k] * sfi_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, pb_x, pc_x, pc_z, sdk0_185, sdk0_186, sdi_59, \
                         sdi_145, sdi_146, sdk1_185, sdk1_186, \
                         sfi_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pb_x[k] * sdk0_185[k]
                   + f_16 * sdi_145[k]
                   - f_12 * pc_x[k] * sdk1_185[k];

        t_186[k] = pb_x[k] * sdk0_186[k]
                   + f_15 * sdi_146[k]
                   - f_12 * pc_x[k] * sdk1_186[k];

        t_187[k] = f_14 * sdi_59[k]
                   + f_3 * pc_z[k] * sfi_143[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, pb_x, pc_x, pc_y, sdk0_189, sdk0_190, sdi_149, \
                         sdi_150, sdk1_189, sdk1_190, sfi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_3 * pc_y[k] * sfi_145[k];

        t_189[k] = pb_x[k] * sdk0_189[k]
                   + f_15 * sdi_149[k]
                   - f_12 * pc_x[k] * sdk1_189[k];

        t_190[k] = pb_x[k] * sdk0_190[k]
                   + f_0 * sdi_150[k]
                   - f_12 * pc_x[k] * sdk1_190[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, pb_x, pc_x, pc_y, pc_z, sdk0_192, sdi_62, \
                         sdi_152, sdk1_192, sfi_146, sfi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_14 * sdi_62[k]
                   + f_3 * pc_z[k] * sfi_146[k];

        t_192[k] = pb_x[k] * sdk0_192[k]
                   + f_0 * sdi_152[k]
                   - f_12 * pc_x[k] * sdk1_192[k];

        t_193[k] = f_3 * pc_y[k] * sfi_149[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, pb_x, pc_x, pc_z, sdk0_194, sdk0_195, sdi_66, \
                         sdi_154, sdi_155, sdk1_194, sdk1_195, \
                         sfi_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = pb_x[k] * sdk0_194[k]
                   + f_0 * sdi_154[k]
                   - f_12 * pc_x[k] * sdk1_194[k];

        t_195[k] = pb_x[k] * sdk0_195[k]
                   + f_14 * sdi_155[k]
                   - f_12 * pc_x[k] * sdk1_195[k];

        t_196[k] = f_14 * sdi_66[k]
                   + f_3 * pc_z[k] * sfi_150[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pb_x, pc_x, pc_y, sdk0_197, sdk0_198, sdi_157, \
                         sdi_158, sdk1_197, sdk1_198, sfi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = pb_x[k] * sdk0_197[k]
                   + f_14 * sdi_157[k]
                   - f_12 * pc_x[k] * sdk1_197[k];

        t_198[k] = pb_x[k] * sdk0_198[k]
                   + f_14 * sdi_158[k]
                   - f_12 * pc_x[k] * sdk1_198[k];

        t_199[k] = f_3 * pc_y[k] * sfi_154[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pb_x, pc_x, sdk0_200, sdi_160, sdi_161, \
                         sdi_162, sdi_163, sdk1_200, sfi_161, sfi_162, \
                         sfi_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pb_x[k] * sdk0_200[k]
                   + f_14 * sdi_160[k]
                   - f_12 * pc_x[k] * sdk1_200[k];

        t_201[k] = f_13 * sdi_161[k]
                   + f_3 * pc_x[k] * sfi_161[k];

        t_202[k] = f_13 * sdi_162[k]
                   + f_3 * pc_x[k] * sfi_162[k];

        t_203[k] = f_13 * sdi_163[k]
                   + f_3 * pc_x[k] * sfi_163[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pc_x, sdi_164, sdi_165, sdi_166, sdi_167, \
                         sfi_164, sfi_165, sfi_166, sfi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_13 * sdi_164[k]
                   + f_3 * pc_x[k] * sfi_164[k];

        t_205[k] = f_13 * sdi_165[k]
                   + f_3 * pc_x[k] * sfi_165[k];

        t_206[k] = f_13 * sdi_166[k]
                   + f_3 * pc_x[k] * sfi_166[k];

        t_207[k] = f_13 * sdi_167[k]
                   + f_3 * pc_x[k] * sfi_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_x, pc_x, pc_z, sdk0_208, sdk0_210, \
                         sdk0_211, sdi_77, sdk1_208, sdk1_210, sdk1_211, \
                         sfi_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pb_x[k] * sdk0_208[k]
                   - f_12 * pc_x[k] * sdk1_208[k];

        t_209[k] = f_14 * sdi_77[k]
                   + f_3 * pc_z[k] * sfi_161[k];

        t_210[k] = pb_x[k] * sdk0_210[k]
                   - f_12 * pc_x[k] * sdk1_210[k];

        t_211[k] = pb_x[k] * sdk0_211[k]
                   - f_12 * pc_x[k] * sdk1_211[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pb_x, pc_x, pc_y, sdk0_212, sdk0_213, \
                         sdk0_215, sdk1_212, sdk1_213, sdk1_215, \
                         sfi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pb_x[k] * sdk0_212[k]
                   - f_12 * pc_x[k] * sdk1_212[k];

        t_213[k] = pb_x[k] * sdk0_213[k]
                   - f_12 * pc_x[k] * sdk1_213[k];

        t_214[k] = f_3 * pc_y[k] * sfi_167[k];

        t_215[k] = pb_x[k] * sdk0_215[k]
                   - f_12 * pc_x[k] * sdk1_215[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pc_x, pc_y, pc_z, sdi_84, sfh0_126, \
                         sfh0_129, sfh1_126, sfh1_129, sfi_168, \
                         sfi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_1 * sfh0_126[k]
                   - f_2 * sfh1_126[k]
                   + f_3 * pc_x[k] * sfi_168[k];

        t_217[k] = f_0 * sdi_84[k]
                   + f_3 * pc_y[k] * sfi_168[k];

        t_218[k] = f_3 * pc_z[k] * sfi_168[k];

        t_219[k] = f_4 * sfh0_129[k]
                   - f_5 * sfh1_129[k]
                   + f_3 * pc_x[k] * sfi_171[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pc_x, pc_y, pc_z, sdi_86, sfh0_131, \
                         sfh0_132, sfh1_131, sfh1_132, sfi_170, sfi_171, sfi_173, \
                         sfi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_0 * sdi_86[k]
                   + f_3 * pc_y[k] * sfi_170[k];

        t_221[k] = f_4 * sfh0_131[k]
                   - f_5 * sfh1_131[k]
                   + f_3 * pc_x[k] * sfi_173[k];

        t_222[k] = f_6 * sfh0_132[k]
                   - f_7 * sfh1_132[k]
                   + f_3 * pc_x[k] * sfi_174[k];

        t_223[k] = f_3 * pc_z[k] * sfi_171[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pc_x, pc_y, pc_z, sdi_89, sfh0_135, \
                         sfh0_136, sfh1_135, sfh1_136, sfi_173, sfi_174, sfi_177, \
                         sfi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_0 * sdi_89[k]
                   + f_3 * pc_y[k] * sfi_173[k];

        t_225[k] = f_6 * sfh0_135[k]
                   - f_7 * sfh1_135[k]
                   + f_3 * pc_x[k] * sfi_177[k];

        t_226[k] = f_8 * sfh0_136[k]
                   - f_9 * sfh1_136[k]
                   + f_3 * pc_x[k] * sfi_178[k];

        t_227[k] = f_3 * pc_z[k] * sfi_174[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pc_x, pc_y, sdi_93, sfh0_138, sfh0_140, \
                         sfh1_138, sfh1_140, sfi_177, sfi_180, \
                         sfi_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_8 * sfh0_138[k]
                   - f_9 * sfh1_138[k]
                   + f_3 * pc_x[k] * sfi_180[k];

        t_229[k] = f_0 * sdi_93[k]
                   + f_3 * pc_y[k] * sfi_177[k];

        t_230[k] = f_8 * sfh0_140[k]
                   - f_9 * sfh1_140[k]
                   + f_3 * pc_x[k] * sfi_182[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pc_x, pc_z, sfh0_141, sfh0_143, sfh0_144, \
                         sfh1_141, sfh1_143, sfh1_144, sfi_178, sfi_183, sfi_185, \
                         sfi_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_10 * sfh0_141[k]
                   - f_11 * sfh1_141[k]
                   + f_3 * pc_x[k] * sfi_183[k];

        t_232[k] = f_3 * pc_z[k] * sfi_178[k];

        t_233[k] = f_10 * sfh0_143[k]
                   - f_11 * sfh1_143[k]
                   + f_3 * pc_x[k] * sfi_185[k];

        t_234[k] = f_10 * sfh0_144[k]
                   - f_11 * sfh1_144[k]
                   + f_3 * pc_x[k] * sfi_186[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, pc_x, pc_y, sdi_98, sfh0_146, \
                         sfh1_146, sfi_182, sfi_188, sfi_189, sfi_190, \
                         sfi_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_0 * sdi_98[k]
                   + f_3 * pc_y[k] * sfi_182[k];

        t_236[k] = f_10 * sfh0_146[k]
                   - f_11 * sfh1_146[k]
                   + f_3 * pc_x[k] * sfi_188[k];

        t_237[k] = f_3 * pc_x[k] * sfi_189[k];

        t_238[k] = f_3 * pc_x[k] * sfi_190[k];

        t_239[k] = f_3 * pc_x[k] * sfi_191[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, pc_y, sdi_105, sfh0_141, \
                         sfh1_141, sfi_189, sfi_192, sfi_193, sfi_194, \
                         sfi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_3 * pc_x[k] * sfi_192[k];

        t_241[k] = f_3 * pc_x[k] * sfi_193[k];

        t_242[k] = f_3 * pc_x[k] * sfi_194[k];

        t_243[k] = f_3 * pc_x[k] * sfi_195[k];

        t_244[k] = f_0 * sdi_105[k]
                   + f_1 * sfh0_141[k]
                   - f_2 * sfh1_141[k]
                   + f_3 * pc_y[k] * sfi_189[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, sdi_107, sdi_108, sfh0_143, \
                         sfh0_144, sfh1_143, sfh1_144, sfi_189, sfi_191, \
                         sfi_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_3 * pc_z[k] * sfi_189[k];

        t_246[k] = f_0 * sdi_107[k]
                   + f_4 * sfh0_143[k]
                   - f_5 * sfh1_143[k]
                   + f_3 * pc_y[k] * sfi_191[k];

        t_247[k] = f_0 * sdi_108[k]
                   + f_6 * sfh0_144[k]
                   - f_7 * sfh1_144[k]
                   + f_3 * pc_y[k] * sfi_192[k];
    }
}

static auto
compute_prim_sfk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sdk0,
                                                          const size_t sdi, const size_t sdk1,
                                                          const size_t sfh0, const size_t sfh1,
                                                          const size_t sfi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 2.0 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 3.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdk0_108 = buffer.data(sdk0 + 108);
    const auto *sdk0_111 = buffer.data(sdk0 + 111);
    const auto *sdk0_114 = buffer.data(sdk0 + 114);
    const auto *sdk0_118 = buffer.data(sdk0 + 118);
    const auto *sdk0_123 = buffer.data(sdk0 + 123);
    const auto *sdk0_136 = buffer.data(sdk0 + 136);
    const auto *sdk0_138 = buffer.data(sdk0 + 138);
    const auto *sdk0_139 = buffer.data(sdk0 + 139);
    const auto *sdk0_140 = buffer.data(sdk0 + 140);
    const auto *sdk0_141 = buffer.data(sdk0 + 141);
    const auto *sdk0_180 = buffer.data(sdk0 + 180);
    const auto *sdk0_185 = buffer.data(sdk0 + 185);
    const auto *sdk0_189 = buffer.data(sdk0 + 189);
    const auto *sdk0_194 = buffer.data(sdk0 + 194);
    const auto *sdk0_200 = buffer.data(sdk0 + 200);
    const auto *sdk0_208 = buffer.data(sdk0 + 208);
    const auto *sdk0_210 = buffer.data(sdk0 + 210);
    const auto *sdk0_211 = buffer.data(sdk0 + 211);
    const auto *sdk0_212 = buffer.data(sdk0 + 212);
    const auto *sdk0_213 = buffer.data(sdk0 + 213);
    const auto *sdk0_215 = buffer.data(sdk0 + 215);

    const auto *sdi_84 = buffer.data(sdi + 84);
    const auto *sdi_87 = buffer.data(sdi + 87);
    const auto *sdi_90 = buffer.data(sdi + 90);
    const auto *sdi_94 = buffer.data(sdi + 94);
    const auto *sdi_105 = buffer.data(sdi + 105);
    const auto *sdi_106 = buffer.data(sdi + 106);
    const auto *sdi_107 = buffer.data(sdi + 107);
    const auto *sdi_108 = buffer.data(sdi + 108);
    const auto *sdi_109 = buffer.data(sdi + 109);
    const auto *sdi_110 = buffer.data(sdi + 110);
    const auto *sdi_111 = buffer.data(sdi + 111);
    const auto *sdi_112 = buffer.data(sdi + 112);
    const auto *sdi_114 = buffer.data(sdi + 114);
    const auto *sdi_115 = buffer.data(sdi + 115);
    const auto *sdi_117 = buffer.data(sdi + 117);
    const auto *sdi_118 = buffer.data(sdi + 118);
    const auto *sdi_121 = buffer.data(sdi + 121);
    const auto *sdi_122 = buffer.data(sdi + 122);
    const auto *sdi_126 = buffer.data(sdi + 126);
    const auto *sdi_133 = buffer.data(sdi + 133);
    const auto *sdi_139 = buffer.data(sdi + 139);
    const auto *sdi_140 = buffer.data(sdi + 140);
    const auto *sdi_142 = buffer.data(sdi + 142);
    const auto *sdi_143 = buffer.data(sdi + 143);
    const auto *sdi_145 = buffer.data(sdi + 145);
    const auto *sdi_146 = buffer.data(sdi + 146);
    const auto *sdi_149 = buffer.data(sdi + 149);
    const auto *sdi_150 = buffer.data(sdi + 150);
    const auto *sdi_154 = buffer.data(sdi + 154);
    const auto *sdi_161 = buffer.data(sdi + 161);
    const auto *sdi_163 = buffer.data(sdi + 163);
    const auto *sdi_164 = buffer.data(sdi + 164);
    const auto *sdi_165 = buffer.data(sdi + 165);
    const auto *sdi_166 = buffer.data(sdi + 166);
    const auto *sdi_167 = buffer.data(sdi + 167);

    const auto *sdk1_108 = buffer.data(sdk1 + 108);
    const auto *sdk1_111 = buffer.data(sdk1 + 111);
    const auto *sdk1_114 = buffer.data(sdk1 + 114);
    const auto *sdk1_118 = buffer.data(sdk1 + 118);
    const auto *sdk1_123 = buffer.data(sdk1 + 123);
    const auto *sdk1_136 = buffer.data(sdk1 + 136);
    const auto *sdk1_138 = buffer.data(sdk1 + 138);
    const auto *sdk1_139 = buffer.data(sdk1 + 139);
    const auto *sdk1_140 = buffer.data(sdk1 + 140);
    const auto *sdk1_141 = buffer.data(sdk1 + 141);
    const auto *sdk1_180 = buffer.data(sdk1 + 180);
    const auto *sdk1_185 = buffer.data(sdk1 + 185);
    const auto *sdk1_189 = buffer.data(sdk1 + 189);
    const auto *sdk1_194 = buffer.data(sdk1 + 194);
    const auto *sdk1_200 = buffer.data(sdk1 + 200);
    const auto *sdk1_208 = buffer.data(sdk1 + 208);
    const auto *sdk1_210 = buffer.data(sdk1 + 210);
    const auto *sdk1_211 = buffer.data(sdk1 + 211);
    const auto *sdk1_212 = buffer.data(sdk1 + 212);
    const auto *sdk1_213 = buffer.data(sdk1 + 213);
    const auto *sdk1_215 = buffer.data(sdk1 + 215);

    const auto *sfh0_145 = buffer.data(sfh0 + 145);
    const auto *sfh0_146 = buffer.data(sfh0 + 146);
    const auto *sfh0_152 = buffer.data(sfh0 + 152);
    const auto *sfh0_156 = buffer.data(sfh0 + 156);
    const auto *sfh0_159 = buffer.data(sfh0 + 159);
    const auto *sfh0_161 = buffer.data(sfh0 + 161);
    const auto *sfh0_164 = buffer.data(sfh0 + 164);
    const auto *sfh0_165 = buffer.data(sfh0 + 165);
    const auto *sfh0_167 = buffer.data(sfh0 + 167);
    const auto *sfh0_171 = buffer.data(sfh0 + 171);
    const auto *sfh0_174 = buffer.data(sfh0 + 174);
    const auto *sfh0_178 = buffer.data(sfh0 + 178);
    const auto *sfh0_180 = buffer.data(sfh0 + 180);
    const auto *sfh0_183 = buffer.data(sfh0 + 183);
    const auto *sfh0_185 = buffer.data(sfh0 + 185);
    const auto *sfh0_186 = buffer.data(sfh0 + 186);
    const auto *sfh0_189 = buffer.data(sfh0 + 189);
    const auto *sfh0_192 = buffer.data(sfh0 + 192);
    const auto *sfh0_194 = buffer.data(sfh0 + 194);
    const auto *sfh0_195 = buffer.data(sfh0 + 195);
    const auto *sfh0_198 = buffer.data(sfh0 + 198);
    const auto *sfh0_199 = buffer.data(sfh0 + 199);
    const auto *sfh0_201 = buffer.data(sfh0 + 201);
    const auto *sfh0_203 = buffer.data(sfh0 + 203);
    const auto *sfh0_204 = buffer.data(sfh0 + 204);
    const auto *sfh0_206 = buffer.data(sfh0 + 206);
    const auto *sfh0_207 = buffer.data(sfh0 + 207);
    const auto *sfh0_208 = buffer.data(sfh0 + 208);
    const auto *sfh0_209 = buffer.data(sfh0 + 209);

    const auto *sfh1_145 = buffer.data(sfh1 + 145);
    const auto *sfh1_146 = buffer.data(sfh1 + 146);
    const auto *sfh1_152 = buffer.data(sfh1 + 152);
    const auto *sfh1_156 = buffer.data(sfh1 + 156);
    const auto *sfh1_159 = buffer.data(sfh1 + 159);
    const auto *sfh1_161 = buffer.data(sfh1 + 161);
    const auto *sfh1_164 = buffer.data(sfh1 + 164);
    const auto *sfh1_165 = buffer.data(sfh1 + 165);
    const auto *sfh1_167 = buffer.data(sfh1 + 167);
    const auto *sfh1_171 = buffer.data(sfh1 + 171);
    const auto *sfh1_174 = buffer.data(sfh1 + 174);
    const auto *sfh1_178 = buffer.data(sfh1 + 178);
    const auto *sfh1_180 = buffer.data(sfh1 + 180);
    const auto *sfh1_183 = buffer.data(sfh1 + 183);
    const auto *sfh1_185 = buffer.data(sfh1 + 185);
    const auto *sfh1_186 = buffer.data(sfh1 + 186);
    const auto *sfh1_189 = buffer.data(sfh1 + 189);
    const auto *sfh1_192 = buffer.data(sfh1 + 192);
    const auto *sfh1_194 = buffer.data(sfh1 + 194);
    const auto *sfh1_195 = buffer.data(sfh1 + 195);
    const auto *sfh1_198 = buffer.data(sfh1 + 198);
    const auto *sfh1_199 = buffer.data(sfh1 + 199);
    const auto *sfh1_201 = buffer.data(sfh1 + 201);
    const auto *sfh1_203 = buffer.data(sfh1 + 203);
    const auto *sfh1_204 = buffer.data(sfh1 + 204);
    const auto *sfh1_206 = buffer.data(sfh1 + 206);
    const auto *sfh1_207 = buffer.data(sfh1 + 207);
    const auto *sfh1_208 = buffer.data(sfh1 + 208);
    const auto *sfh1_209 = buffer.data(sfh1 + 209);

    const auto *sfi_193 = buffer.data(sfi + 193);
    const auto *sfi_194 = buffer.data(sfi + 194);
    const auto *sfi_195 = buffer.data(sfi + 195);
    const auto *sfi_196 = buffer.data(sfi + 196);
    const auto *sfi_198 = buffer.data(sfi + 198);
    const auto *sfi_199 = buffer.data(sfi + 199);
    const auto *sfi_201 = buffer.data(sfi + 201);
    const auto *sfi_202 = buffer.data(sfi + 202);
    const auto *sfi_205 = buffer.data(sfi + 205);
    const auto *sfi_206 = buffer.data(sfi + 206);
    const auto *sfi_208 = buffer.data(sfi + 208);
    const auto *sfi_210 = buffer.data(sfi + 210);
    const auto *sfi_213 = buffer.data(sfi + 213);
    const auto *sfi_214 = buffer.data(sfi + 214);
    const auto *sfi_216 = buffer.data(sfi + 216);
    const auto *sfi_217 = buffer.data(sfi + 217);
    const auto *sfi_218 = buffer.data(sfi + 218);
    const auto *sfi_219 = buffer.data(sfi + 219);
    const auto *sfi_220 = buffer.data(sfi + 220);
    const auto *sfi_221 = buffer.data(sfi + 221);
    const auto *sfi_222 = buffer.data(sfi + 222);
    const auto *sfi_223 = buffer.data(sfi + 223);
    const auto *sfi_224 = buffer.data(sfi + 224);
    const auto *sfi_226 = buffer.data(sfi + 226);
    const auto *sfi_227 = buffer.data(sfi + 227);
    const auto *sfi_229 = buffer.data(sfi + 229);
    const auto *sfi_230 = buffer.data(sfi + 230);
    const auto *sfi_233 = buffer.data(sfi + 233);
    const auto *sfi_234 = buffer.data(sfi + 234);
    const auto *sfi_236 = buffer.data(sfi + 236);
    const auto *sfi_238 = buffer.data(sfi + 238);
    const auto *sfi_239 = buffer.data(sfi + 239);
    const auto *sfi_241 = buffer.data(sfi + 241);
    const auto *sfi_242 = buffer.data(sfi + 242);
    const auto *sfi_245 = buffer.data(sfi + 245);
    const auto *sfi_246 = buffer.data(sfi + 246);
    const auto *sfi_247 = buffer.data(sfi + 247);
    const auto *sfi_248 = buffer.data(sfi + 248);
    const auto *sfi_249 = buffer.data(sfi + 249);
    const auto *sfi_250 = buffer.data(sfi + 250);
    const auto *sfi_251 = buffer.data(sfi + 251);
    const auto *sfi_252 = buffer.data(sfi + 252);
    const auto *sfi_254 = buffer.data(sfi + 254);
    const auto *sfi_255 = buffer.data(sfi + 255);
    const auto *sfi_257 = buffer.data(sfi + 257);
    const auto *sfi_258 = buffer.data(sfi + 258);
    const auto *sfi_261 = buffer.data(sfi + 261);
    const auto *sfi_262 = buffer.data(sfi + 262);
    const auto *sfi_264 = buffer.data(sfi + 264);
    const auto *sfi_266 = buffer.data(sfi + 266);
    const auto *sfi_267 = buffer.data(sfi + 267);
    const auto *sfi_269 = buffer.data(sfi + 269);
    const auto *sfi_270 = buffer.data(sfi + 270);
    const auto *sfi_272 = buffer.data(sfi + 272);
    const auto *sfi_273 = buffer.data(sfi + 273);
    const auto *sfi_274 = buffer.data(sfi + 274);
    const auto *sfi_275 = buffer.data(sfi + 275);
    const auto *sfi_276 = buffer.data(sfi + 276);
    const auto *sfi_277 = buffer.data(sfi + 277);
    const auto *sfi_278 = buffer.data(sfi + 278);
    const auto *sfi_279 = buffer.data(sfi + 279);

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pc_y, pc_z, sdi_109, sdi_110, sdi_111, \
                         sfh0_145, sfh0_146, sfh1_145, sfh1_146, sfi_193, sfi_194, \
                         sfi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_0 * sdi_109[k]
                   + f_8 * sfh0_145[k]
                   - f_9 * sfh1_145[k]
                   + f_3 * pc_y[k] * sfi_193[k];

        t_249[k] = f_0 * sdi_110[k]
                   + f_10 * sfh0_146[k]
                   - f_11 * sfh1_146[k]
                   + f_3 * pc_y[k] * sfi_194[k];

        t_250[k] = f_0 * sdi_111[k]
                   + f_3 * pc_y[k] * sfi_195[k];

        t_251[k] = f_1 * sfh0_146[k]
                   - f_2 * sfh1_146[k]
                   + f_3 * pc_z[k] * sfi_195[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pb_z, pc_y, pc_z, sdk0_108, sdk0_111, \
                         sdi_84, sdi_112, sdk1_108, sdk1_111, sfi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = pb_z[k] * sdk0_108[k]
                   - f_12 * pc_z[k] * sdk1_108[k];

        t_253[k] = f_14 * sdi_112[k]
                   + f_3 * pc_y[k] * sfi_196[k];

        t_254[k] = f_13 * sdi_84[k]
                   + f_3 * pc_z[k] * sfi_196[k];

        t_255[k] = pb_z[k] * sdk0_111[k]
                   - f_12 * pc_z[k] * sdk1_111[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, pb_z, pc_x, pc_y, pc_z, sdk0_114, sdi_114, \
                         sdk1_114, sfh0_152, sfh1_152, sfi_198, \
                         sfi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_14 * sdi_114[k]
                   + f_3 * pc_y[k] * sfi_198[k];

        t_257[k] = f_4 * sfh0_152[k]
                   - f_5 * sfh1_152[k]
                   + f_3 * pc_x[k] * sfi_201[k];

        t_258[k] = pb_z[k] * sdk0_114[k]
                   - f_12 * pc_z[k] * sdk1_114[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, pc_x, pc_y, pc_z, sdi_87, sdi_117, sfh0_156, \
                         sfh1_156, sfi_199, sfi_201, sfi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_13 * sdi_87[k]
                   + f_3 * pc_z[k] * sfi_199[k];

        t_260[k] = f_14 * sdi_117[k]
                   + f_3 * pc_y[k] * sfi_201[k];

        t_261[k] = f_6 * sfh0_156[k]
                   - f_7 * sfh1_156[k]
                   + f_3 * pc_x[k] * sfi_205[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, pb_z, pc_x, pc_z, sdk0_118, sdi_90, sdk1_118, \
                         sfh0_159, sfh1_159, sfi_202, sfi_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = pb_z[k] * sdk0_118[k]
                   - f_12 * pc_z[k] * sdk1_118[k];

        t_263[k] = f_13 * sdi_90[k]
                   + f_3 * pc_z[k] * sfi_202[k];

        t_264[k] = f_8 * sfh0_159[k]
                   - f_9 * sfh1_159[k]
                   + f_3 * pc_x[k] * sfi_208[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, pb_z, pc_x, pc_y, pc_z, sdk0_123, sdi_121, \
                         sdk1_123, sfh0_161, sfh1_161, sfi_205, \
                         sfi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_14 * sdi_121[k]
                   + f_3 * pc_y[k] * sfi_205[k];

        t_266[k] = f_8 * sfh0_161[k]
                   - f_9 * sfh1_161[k]
                   + f_3 * pc_x[k] * sfi_210[k];

        t_267[k] = pb_z[k] * sdk0_123[k]
                   - f_12 * pc_z[k] * sdk1_123[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_x, pc_z, sdi_94, sfh0_164, sfh0_165, \
                         sfh1_164, sfh1_165, sfi_206, sfi_213, \
                         sfi_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_13 * sdi_94[k]
                   + f_3 * pc_z[k] * sfi_206[k];

        t_269[k] = f_10 * sfh0_164[k]
                   - f_11 * sfh1_164[k]
                   + f_3 * pc_x[k] * sfi_213[k];

        t_270[k] = f_10 * sfh0_165[k]
                   - f_11 * sfh1_165[k]
                   + f_3 * pc_x[k] * sfi_214[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, t_275, pc_x, pc_y, sdi_126, sfh0_167, \
                         sfh1_167, sfi_210, sfi_216, sfi_217, sfi_218, \
                         sfi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_14 * sdi_126[k]
                   + f_3 * pc_y[k] * sfi_210[k];

        t_272[k] = f_10 * sfh0_167[k]
                   - f_11 * sfh1_167[k]
                   + f_3 * pc_x[k] * sfi_216[k];

        t_273[k] = f_3 * pc_x[k] * sfi_217[k];

        t_274[k] = f_3 * pc_x[k] * sfi_218[k];

        t_275[k] = f_3 * pc_x[k] * sfi_219[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, pb_z, pc_x, pc_z, sdk0_136, \
                         sdk1_136, sfi_220, sfi_221, sfi_222, sfi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_3 * pc_x[k] * sfi_220[k];

        t_277[k] = f_3 * pc_x[k] * sfi_221[k];

        t_278[k] = f_3 * pc_x[k] * sfi_222[k];

        t_279[k] = f_3 * pc_x[k] * sfi_223[k];

        t_280[k] = pb_z[k] * sdk0_136[k]
                   - f_12 * pc_z[k] * sdk1_136[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, pb_z, pc_z, sdk0_138, sdk0_139, sdi_105, \
                         sdi_106, sdi_107, sdk1_138, sdk1_139, \
                         sfi_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_13 * sdi_105[k]
                   + f_3 * pc_z[k] * sfi_217[k];

        t_282[k] = pb_z[k] * sdk0_138[k]
                   + f_14 * sdi_106[k]
                   - f_12 * pc_z[k] * sdk1_138[k];

        t_283[k] = pb_z[k] * sdk0_139[k]
                   + f_0 * sdi_107[k]
                   - f_12 * pc_z[k] * sdk1_139[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pb_z, pc_y, pc_z, sdk0_140, sdk0_141, sdi_108, \
                         sdi_109, sdi_139, sdk1_140, sdk1_141, \
                         sfi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = pb_z[k] * sdk0_140[k]
                   + f_15 * sdi_108[k]
                   - f_12 * pc_z[k] * sdk1_140[k];

        t_285[k] = pb_z[k] * sdk0_141[k]
                   + f_16 * sdi_109[k]
                   - f_12 * pc_z[k] * sdk1_141[k];

        t_286[k] = f_14 * sdi_139[k]
                   + f_3 * pc_y[k] * sfi_223[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pb_y, pc_y, pc_z, sdk0_180, sdi_111, \
                         sdi_112, sdi_140, sdk1_180, sfh0_167, sfh1_167, sfi_223, \
                         sfi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_13 * sdi_111[k]
                   + f_1 * sfh0_167[k]
                   - f_2 * sfh1_167[k]
                   + f_3 * pc_z[k] * sfi_223[k];

        t_288[k] = pb_y[k] * sdk0_180[k]
                   - f_12 * pc_y[k] * sdk1_180[k];

        t_289[k] = f_13 * sdi_140[k]
                   + f_3 * pc_y[k] * sfi_224[k];

        t_290[k] = f_14 * sdi_112[k]
                   + f_3 * pc_z[k] * sfi_224[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pb_y, pc_x, pc_y, sdk0_185, sdi_142, sdk1_185, \
                         sfh0_171, sfh1_171, sfi_226, sfi_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_4 * sfh0_171[k]
                   - f_5 * sfh1_171[k]
                   + f_3 * pc_x[k] * sfi_227[k];

        t_292[k] = f_13 * sdi_142[k]
                   + f_3 * pc_y[k] * sfi_226[k];

        t_293[k] = pb_y[k] * sdk0_185[k]
                   - f_12 * pc_y[k] * sdk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pc_x, pc_y, pc_z, sdi_115, sdi_145, sfh0_174, \
                         sfh1_174, sfi_227, sfi_229, sfi_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_6 * sfh0_174[k]
                   - f_7 * sfh1_174[k]
                   + f_3 * pc_x[k] * sfi_230[k];

        t_295[k] = f_14 * sdi_115[k]
                   + f_3 * pc_z[k] * sfi_227[k];

        t_296[k] = f_13 * sdi_145[k]
                   + f_3 * pc_y[k] * sfi_229[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_y, pc_x, pc_y, pc_z, sdk0_189, sdi_118, \
                         sdk1_189, sfh0_178, sfh1_178, sfi_230, \
                         sfi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = pb_y[k] * sdk0_189[k]
                   - f_12 * pc_y[k] * sdk1_189[k];

        t_298[k] = f_8 * sfh0_178[k]
                   - f_9 * sfh1_178[k]
                   + f_3 * pc_x[k] * sfi_234[k];

        t_299[k] = f_14 * sdi_118[k]
                   + f_3 * pc_z[k] * sfi_230[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, pb_y, pc_x, pc_y, sdk0_194, sdi_149, sdk1_194, \
                         sfh0_180, sfh1_180, sfi_233, sfi_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_8 * sfh0_180[k]
                   - f_9 * sfh1_180[k]
                   + f_3 * pc_x[k] * sfi_236[k];

        t_301[k] = f_13 * sdi_149[k]
                   + f_3 * pc_y[k] * sfi_233[k];

        t_302[k] = pb_y[k] * sdk0_194[k]
                   - f_12 * pc_y[k] * sdk1_194[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, pc_x, pc_z, sdi_122, sfh0_183, sfh0_185, \
                         sfh1_183, sfh1_185, sfi_234, sfi_239, \
                         sfi_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_10 * sfh0_183[k]
                   - f_11 * sfh1_183[k]
                   + f_3 * pc_x[k] * sfi_239[k];

        t_304[k] = f_14 * sdi_122[k]
                   + f_3 * pc_z[k] * sfi_234[k];

        t_305[k] = f_10 * sfh0_185[k]
                   - f_11 * sfh1_185[k]
                   + f_3 * pc_x[k] * sfi_241[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pb_y, pc_x, pc_y, sdk0_200, sdi_154, \
                         sdk1_200, sfh0_186, sfh1_186, sfi_238, sfi_242, \
                         sfi_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_10 * sfh0_186[k]
                   - f_11 * sfh1_186[k]
                   + f_3 * pc_x[k] * sfi_242[k];

        t_307[k] = f_13 * sdi_154[k]
                   + f_3 * pc_y[k] * sfi_238[k];

        t_308[k] = pb_y[k] * sdk0_200[k]
                   - f_12 * pc_y[k] * sdk1_200[k];

        t_309[k] = f_3 * pc_x[k] * sfi_245[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, t_315, pc_x, sfi_246, sfi_247, \
                         sfi_248, sfi_249, sfi_250, sfi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_3 * pc_x[k] * sfi_246[k];

        t_311[k] = f_3 * pc_x[k] * sfi_247[k];

        t_312[k] = f_3 * pc_x[k] * sfi_248[k];

        t_313[k] = f_3 * pc_x[k] * sfi_249[k];

        t_314[k] = f_3 * pc_x[k] * sfi_250[k];

        t_315[k] = f_3 * pc_x[k] * sfi_251[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pb_y, pc_y, pc_z, sdk0_208, sdk0_210, sdi_133, \
                         sdi_161, sdi_163, sdk1_208, sdk1_210, \
                         sfi_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = pb_y[k] * sdk0_208[k]
                   + f_17 * sdi_161[k]
                   - f_12 * pc_y[k] * sdk1_208[k];

        t_317[k] = f_14 * sdi_133[k]
                   + f_3 * pc_z[k] * sfi_245[k];

        t_318[k] = pb_y[k] * sdk0_210[k]
                   + f_16 * sdi_163[k]
                   - f_12 * pc_y[k] * sdk1_210[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pb_y, pc_y, sdk0_211, sdk0_212, sdk0_213, \
                         sdi_164, sdi_165, sdi_166, sdk1_211, sdk1_212, \
                         sdk1_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = pb_y[k] * sdk0_211[k]
                   + f_15 * sdi_164[k]
                   - f_12 * pc_y[k] * sdk1_211[k];

        t_320[k] = pb_y[k] * sdk0_212[k]
                   + f_0 * sdi_165[k]
                   - f_12 * pc_y[k] * sdk1_212[k];

        t_321[k] = pb_y[k] * sdk0_213[k]
                   + f_14 * sdi_166[k]
                   - f_12 * pc_y[k] * sdk1_213[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pb_y, pc_x, pc_y, sdk0_215, sdi_167, \
                         sdk1_215, sfh0_189, sfh1_189, sfi_251, \
                         sfi_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_13 * sdi_167[k]
                   + f_3 * pc_y[k] * sfi_251[k];

        t_323[k] = pb_y[k] * sdk0_215[k]
                   - f_12 * pc_y[k] * sdk1_215[k];

        t_324[k] = f_1 * sfh0_189[k]
                   - f_2 * sfh1_189[k]
                   + f_3 * pc_x[k] * sfi_252[k];

        t_325[k] = f_3 * pc_y[k] * sfi_252[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pc_x, pc_y, pc_z, sdi_140, sfh0_192, \
                         sfh0_194, sfh1_192, sfh1_194, sfi_252, sfi_254, sfi_255, \
                         sfi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_0 * sdi_140[k]
                   + f_3 * pc_z[k] * sfi_252[k];

        t_327[k] = f_4 * sfh0_192[k]
                   - f_5 * sfh1_192[k]
                   + f_3 * pc_x[k] * sfi_255[k];

        t_328[k] = f_3 * pc_y[k] * sfi_254[k];

        t_329[k] = f_4 * sfh0_194[k]
                   - f_5 * sfh1_194[k]
                   + f_3 * pc_x[k] * sfi_257[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pc_x, pc_y, pc_z, sdi_143, sfh0_195, \
                         sfh0_198, sfh1_195, sfh1_198, sfi_255, sfi_257, sfi_258, \
                         sfi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_6 * sfh0_195[k]
                   - f_7 * sfh1_195[k]
                   + f_3 * pc_x[k] * sfi_258[k];

        t_331[k] = f_0 * sdi_143[k]
                   + f_3 * pc_z[k] * sfi_255[k];

        t_332[k] = f_3 * pc_y[k] * sfi_257[k];

        t_333[k] = f_6 * sfh0_198[k]
                   - f_7 * sfh1_198[k]
                   + f_3 * pc_x[k] * sfi_261[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pc_x, pc_y, pc_z, sdi_146, sfh0_199, \
                         sfh0_201, sfh1_199, sfh1_201, sfi_258, sfi_261, sfi_262, \
                         sfi_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_8 * sfh0_199[k]
                   - f_9 * sfh1_199[k]
                   + f_3 * pc_x[k] * sfi_262[k];

        t_335[k] = f_0 * sdi_146[k]
                   + f_3 * pc_z[k] * sfi_258[k];

        t_336[k] = f_8 * sfh0_201[k]
                   - f_9 * sfh1_201[k]
                   + f_3 * pc_x[k] * sfi_264[k];

        t_337[k] = f_3 * pc_y[k] * sfi_261[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pc_x, pc_z, sdi_150, sfh0_203, sfh0_204, \
                         sfh1_203, sfh1_204, sfi_262, sfi_266, \
                         sfi_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_8 * sfh0_203[k]
                   - f_9 * sfh1_203[k]
                   + f_3 * pc_x[k] * sfi_266[k];

        t_339[k] = f_10 * sfh0_204[k]
                   - f_11 * sfh1_204[k]
                   + f_3 * pc_x[k] * sfi_267[k];

        t_340[k] = f_0 * sdi_150[k]
                   + f_3 * pc_z[k] * sfi_262[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_y, sfh0_206, sfh0_207, sfh0_209, \
                         sfh1_206, sfh1_207, sfh1_209, sfi_266, sfi_269, sfi_270, \
                         sfi_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_10 * sfh0_206[k]
                   - f_11 * sfh1_206[k]
                   + f_3 * pc_x[k] * sfi_269[k];

        t_342[k] = f_10 * sfh0_207[k]
                   - f_11 * sfh1_207[k]
                   + f_3 * pc_x[k] * sfi_270[k];

        t_343[k] = f_3 * pc_y[k] * sfi_266[k];

        t_344[k] = f_10 * sfh0_209[k]
                   - f_11 * sfh1_209[k]
                   + f_3 * pc_x[k] * sfi_272[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, t_350, t_351, pc_x, sfi_273, \
                         sfi_274, sfi_275, sfi_276, sfi_277, sfi_278, \
                         sfi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_3 * pc_x[k] * sfi_273[k];

        t_346[k] = f_3 * pc_x[k] * sfi_274[k];

        t_347[k] = f_3 * pc_x[k] * sfi_275[k];

        t_348[k] = f_3 * pc_x[k] * sfi_276[k];

        t_349[k] = f_3 * pc_x[k] * sfi_277[k];

        t_350[k] = f_3 * pc_x[k] * sfi_278[k];

        t_351[k] = f_3 * pc_x[k] * sfi_279[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pc_y, pc_z, sdi_161, sfh0_204, sfh0_206, \
                         sfh0_207, sfh1_204, sfh1_206, sfh1_207, sfi_273, sfi_275, \
                         sfi_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_1 * sfh0_204[k]
                   - f_2 * sfh1_204[k]
                   + f_3 * pc_y[k] * sfi_273[k];

        t_353[k] = f_0 * sdi_161[k]
                   + f_3 * pc_z[k] * sfi_273[k];

        t_354[k] = f_4 * sfh0_206[k]
                   - f_5 * sfh1_206[k]
                   + f_3 * pc_y[k] * sfi_275[k];

        t_355[k] = f_6 * sfh0_207[k]
                   - f_7 * sfh1_207[k]
                   + f_3 * pc_y[k] * sfi_276[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pc_y, pc_z, sdi_167, sfh0_208, sfh0_209, \
                         sfh1_208, sfh1_209, sfi_277, sfi_278, \
                         sfi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_8 * sfh0_208[k]
                   - f_9 * sfh1_208[k]
                   + f_3 * pc_y[k] * sfi_277[k];

        t_357[k] = f_10 * sfh0_209[k]
                   - f_11 * sfh1_209[k]
                   + f_3 * pc_y[k] * sfi_278[k];

        t_358[k] = f_3 * pc_y[k] * sfi_279[k];

        t_359[k] = f_0 * sdi_167[k]
                   + f_1 * sfh0_209[k]
                   - f_2 * sfh1_209[k]
                   + f_3 * pc_z[k] * sfi_279[k];
    }
}

auto
compute_prim_sfk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sdk0, const size_t sdi,
                                                   const size_t sdk1, const size_t sfh0,
                                                   const size_t sfh1, const size_t sfi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sfk_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sdk0, sdi,
                                                              sdk1, sfh0, sfh1, sfi, ncols,
                                                              gamma, p, q);

    compute_prim_sfk_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sdk0, sdi,
                                                              sdk1, sfh0, sfh1, sfi, ncols,
                                                              gamma, p, q);

    compute_prim_sfk_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sdk0, sdi,
                                                              sdk1, sfh0, sfh1, sfi, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
