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


#include "SimdThreeCenterElectronRepulsionVrrRecSDI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sdi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t spi0,
                                                          const size_t sph, const size_t spi1,
                                                          const size_t sdg0, const size_t sdg1,
                                                          const size_t sdh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 1.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spi0_0 = buffer.data(spi0 + 0);
    const auto *spi0_3 = buffer.data(spi0 + 3);
    const auto *spi0_5 = buffer.data(spi0 + 5);
    const auto *spi0_6 = buffer.data(spi0 + 6);
    const auto *spi0_9 = buffer.data(spi0 + 9);
    const auto *spi0_10 = buffer.data(spi0 + 10);
    const auto *spi0_14 = buffer.data(spi0 + 14);
    const auto *spi0_31 = buffer.data(spi0 + 31);
    const auto *spi0_34 = buffer.data(spi0 + 34);
    const auto *spi0_38 = buffer.data(spi0 + 38);
    const auto *spi0_40 = buffer.data(spi0 + 40);
    const auto *spi0_49 = buffer.data(spi0 + 49);
    const auto *spi0_51 = buffer.data(spi0 + 51);
    const auto *spi0_52 = buffer.data(spi0 + 52);
    const auto *spi0_53 = buffer.data(spi0 + 53);
    const auto *spi0_55 = buffer.data(spi0 + 55);
    const auto *spi0_56 = buffer.data(spi0 + 56);
    const auto *spi0_61 = buffer.data(spi0 + 61);
    const auto *spi0_65 = buffer.data(spi0 + 65);
    const auto *spi0_68 = buffer.data(spi0 + 68);
    const auto *spi0_70 = buffer.data(spi0 + 70);
    const auto *spi0_77 = buffer.data(spi0 + 77);
    const auto *spi0_79 = buffer.data(spi0 + 79);
    const auto *spi0_80 = buffer.data(spi0 + 80);
    const auto *spi0_81 = buffer.data(spi0 + 81);
    const auto *spi0_83 = buffer.data(spi0 + 83);

    const auto *sph_0 = buffer.data(sph + 0);
    const auto *sph_2 = buffer.data(sph + 2);
    const auto *sph_3 = buffer.data(sph + 3);
    const auto *sph_5 = buffer.data(sph + 5);
    const auto *sph_6 = buffer.data(sph + 6);
    const auto *sph_9 = buffer.data(sph + 9);
    const auto *sph_10 = buffer.data(sph + 10);
    const auto *sph_12 = buffer.data(sph + 12);
    const auto *sph_14 = buffer.data(sph + 14);
    const auto *sph_15 = buffer.data(sph + 15);
    const auto *sph_16 = buffer.data(sph + 16);
    const auto *sph_17 = buffer.data(sph + 17);
    const auto *sph_18 = buffer.data(sph + 18);
    const auto *sph_19 = buffer.data(sph + 19);
    const auto *sph_20 = buffer.data(sph + 20);
    const auto *sph_21 = buffer.data(sph + 21);
    const auto *sph_23 = buffer.data(sph + 23);
    const auto *sph_24 = buffer.data(sph + 24);
    const auto *sph_26 = buffer.data(sph + 26);
    const auto *sph_27 = buffer.data(sph + 27);
    const auto *sph_30 = buffer.data(sph + 30);
    const auto *sph_31 = buffer.data(sph + 31);
    const auto *sph_33 = buffer.data(sph + 33);
    const auto *sph_36 = buffer.data(sph + 36);
    const auto *sph_37 = buffer.data(sph + 37);
    const auto *sph_38 = buffer.data(sph + 38);
    const auto *sph_39 = buffer.data(sph + 39);
    const auto *sph_40 = buffer.data(sph + 40);
    const auto *sph_41 = buffer.data(sph + 41);
    const auto *sph_42 = buffer.data(sph + 42);
    const auto *sph_44 = buffer.data(sph + 44);
    const auto *sph_47 = buffer.data(sph + 47);
    const auto *sph_51 = buffer.data(sph + 51);
    const auto *sph_54 = buffer.data(sph + 54);
    const auto *sph_56 = buffer.data(sph + 56);
    const auto *sph_57 = buffer.data(sph + 57);
    const auto *sph_58 = buffer.data(sph + 58);
    const auto *sph_59 = buffer.data(sph + 59);
    const auto *sph_60 = buffer.data(sph + 60);
    const auto *sph_61 = buffer.data(sph + 61);
    const auto *sph_62 = buffer.data(sph + 62);

    const auto *spi1_0 = buffer.data(spi1 + 0);
    const auto *spi1_3 = buffer.data(spi1 + 3);
    const auto *spi1_5 = buffer.data(spi1 + 5);
    const auto *spi1_6 = buffer.data(spi1 + 6);
    const auto *spi1_9 = buffer.data(spi1 + 9);
    const auto *spi1_10 = buffer.data(spi1 + 10);
    const auto *spi1_14 = buffer.data(spi1 + 14);
    const auto *spi1_31 = buffer.data(spi1 + 31);
    const auto *spi1_34 = buffer.data(spi1 + 34);
    const auto *spi1_38 = buffer.data(spi1 + 38);
    const auto *spi1_40 = buffer.data(spi1 + 40);
    const auto *spi1_49 = buffer.data(spi1 + 49);
    const auto *spi1_51 = buffer.data(spi1 + 51);
    const auto *spi1_52 = buffer.data(spi1 + 52);
    const auto *spi1_53 = buffer.data(spi1 + 53);
    const auto *spi1_55 = buffer.data(spi1 + 55);
    const auto *spi1_56 = buffer.data(spi1 + 56);
    const auto *spi1_61 = buffer.data(spi1 + 61);
    const auto *spi1_65 = buffer.data(spi1 + 65);
    const auto *spi1_68 = buffer.data(spi1 + 68);
    const auto *spi1_70 = buffer.data(spi1 + 70);
    const auto *spi1_77 = buffer.data(spi1 + 77);
    const auto *spi1_79 = buffer.data(spi1 + 79);
    const auto *spi1_80 = buffer.data(spi1 + 80);
    const auto *spi1_81 = buffer.data(spi1 + 81);
    const auto *spi1_83 = buffer.data(spi1 + 83);

    const auto *sdg0_0 = buffer.data(sdg0 + 0);
    const auto *sdg0_3 = buffer.data(sdg0 + 3);
    const auto *sdg0_5 = buffer.data(sdg0 + 5);
    const auto *sdg0_6 = buffer.data(sdg0 + 6);
    const auto *sdg0_9 = buffer.data(sdg0 + 9);
    const auto *sdg0_10 = buffer.data(sdg0 + 10);
    const auto *sdg0_12 = buffer.data(sdg0 + 12);
    const auto *sdg0_13 = buffer.data(sdg0 + 13);
    const auto *sdg0_14 = buffer.data(sdg0 + 14);
    const auto *sdg0_45 = buffer.data(sdg0 + 45);
    const auto *sdg0_48 = buffer.data(sdg0 + 48);
    const auto *sdg0_50 = buffer.data(sdg0 + 50);
    const auto *sdg0_51 = buffer.data(sdg0 + 51);
    const auto *sdg0_54 = buffer.data(sdg0 + 54);
    const auto *sdg0_55 = buffer.data(sdg0 + 55);
    const auto *sdg0_57 = buffer.data(sdg0 + 57);
    const auto *sdg0_58 = buffer.data(sdg0 + 58);
    const auto *sdg0_59 = buffer.data(sdg0 + 59);
    const auto *sdg0_72 = buffer.data(sdg0 + 72);

    const auto *sdg1_0 = buffer.data(sdg1 + 0);
    const auto *sdg1_3 = buffer.data(sdg1 + 3);
    const auto *sdg1_5 = buffer.data(sdg1 + 5);
    const auto *sdg1_6 = buffer.data(sdg1 + 6);
    const auto *sdg1_9 = buffer.data(sdg1 + 9);
    const auto *sdg1_10 = buffer.data(sdg1 + 10);
    const auto *sdg1_12 = buffer.data(sdg1 + 12);
    const auto *sdg1_13 = buffer.data(sdg1 + 13);
    const auto *sdg1_14 = buffer.data(sdg1 + 14);
    const auto *sdg1_45 = buffer.data(sdg1 + 45);
    const auto *sdg1_48 = buffer.data(sdg1 + 48);
    const auto *sdg1_50 = buffer.data(sdg1 + 50);
    const auto *sdg1_51 = buffer.data(sdg1 + 51);
    const auto *sdg1_54 = buffer.data(sdg1 + 54);
    const auto *sdg1_55 = buffer.data(sdg1 + 55);
    const auto *sdg1_57 = buffer.data(sdg1 + 57);
    const auto *sdg1_58 = buffer.data(sdg1 + 58);
    const auto *sdg1_59 = buffer.data(sdg1 + 59);
    const auto *sdg1_72 = buffer.data(sdg1 + 72);

    const auto *sdh_0 = buffer.data(sdh + 0);
    const auto *sdh_2 = buffer.data(sdh + 2);
    const auto *sdh_3 = buffer.data(sdh + 3);
    const auto *sdh_5 = buffer.data(sdh + 5);
    const auto *sdh_6 = buffer.data(sdh + 6);
    const auto *sdh_9 = buffer.data(sdh + 9);
    const auto *sdh_10 = buffer.data(sdh + 10);
    const auto *sdh_12 = buffer.data(sdh + 12);
    const auto *sdh_14 = buffer.data(sdh + 14);
    const auto *sdh_15 = buffer.data(sdh + 15);
    const auto *sdh_16 = buffer.data(sdh + 16);
    const auto *sdh_17 = buffer.data(sdh + 17);
    const auto *sdh_18 = buffer.data(sdh + 18);
    const auto *sdh_19 = buffer.data(sdh + 19);
    const auto *sdh_20 = buffer.data(sdh + 20);
    const auto *sdh_21 = buffer.data(sdh + 21);
    const auto *sdh_23 = buffer.data(sdh + 23);
    const auto *sdh_24 = buffer.data(sdh + 24);
    const auto *sdh_26 = buffer.data(sdh + 26);
    const auto *sdh_27 = buffer.data(sdh + 27);
    const auto *sdh_30 = buffer.data(sdh + 30);
    const auto *sdh_36 = buffer.data(sdh + 36);
    const auto *sdh_37 = buffer.data(sdh + 37);
    const auto *sdh_38 = buffer.data(sdh + 38);
    const auto *sdh_39 = buffer.data(sdh + 39);
    const auto *sdh_40 = buffer.data(sdh + 40);
    const auto *sdh_41 = buffer.data(sdh + 41);
    const auto *sdh_42 = buffer.data(sdh + 42);
    const auto *sdh_44 = buffer.data(sdh + 44);
    const auto *sdh_45 = buffer.data(sdh + 45);
    const auto *sdh_47 = buffer.data(sdh + 47);
    const auto *sdh_48 = buffer.data(sdh + 48);
    const auto *sdh_51 = buffer.data(sdh + 51);
    const auto *sdh_57 = buffer.data(sdh + 57);
    const auto *sdh_58 = buffer.data(sdh + 58);
    const auto *sdh_59 = buffer.data(sdh + 59);
    const auto *sdh_60 = buffer.data(sdh + 60);
    const auto *sdh_61 = buffer.data(sdh + 61);
    const auto *sdh_62 = buffer.data(sdh + 62);
    const auto *sdh_63 = buffer.data(sdh + 63);
    const auto *sdh_65 = buffer.data(sdh + 65);
    const auto *sdh_66 = buffer.data(sdh + 66);
    const auto *sdh_68 = buffer.data(sdh + 68);
    const auto *sdh_69 = buffer.data(sdh + 69);
    const auto *sdh_72 = buffer.data(sdh + 72);
    const auto *sdh_73 = buffer.data(sdh + 73);
    const auto *sdh_75 = buffer.data(sdh + 75);
    const auto *sdh_77 = buffer.data(sdh + 77);
    const auto *sdh_78 = buffer.data(sdh + 78);
    const auto *sdh_79 = buffer.data(sdh + 79);
    const auto *sdh_80 = buffer.data(sdh + 80);
    const auto *sdh_81 = buffer.data(sdh + 81);
    const auto *sdh_82 = buffer.data(sdh + 82);
    const auto *sdh_83 = buffer.data(sdh + 83);
    const auto *sdh_84 = buffer.data(sdh + 84);
    const auto *sdh_86 = buffer.data(sdh + 86);
    const auto *sdh_87 = buffer.data(sdh + 87);
    const auto *sdh_89 = buffer.data(sdh + 89);
    const auto *sdh_90 = buffer.data(sdh + 90);
    const auto *sdh_93 = buffer.data(sdh + 93);
    const auto *sdh_96 = buffer.data(sdh + 96);
    const auto *sdh_99 = buffer.data(sdh + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sph_0, sph_3, sdg0_0, sdg0_3, \
                         sdg1_0, sdg1_3, sdh_0, sdh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sph_0[k]
                 + f_1 * sdg0_0[k]
                 - f_2 * sdg1_0[k]
                 + f_3 * pc_x[k] * sdh_0[k];

        t_1[k] = f_3 * pc_y[k] * sdh_0[k];

        t_2[k] = f_3 * pc_z[k] * sdh_0[k];

        t_3[k] = f_0 * sph_3[k]
                 + f_4 * sdg0_3[k]
                 - f_5 * sdg1_3[k]
                 + f_3 * pc_x[k] * sdh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sph_5, sph_6, sdg0_5, sdg0_6, sdg1_5, \
                         sdg1_6, sdh_2, sdh_5, sdh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sdh_2[k];

        t_5[k] = f_0 * sph_5[k]
                 + f_4 * sdg0_5[k]
                 - f_5 * sdg1_5[k]
                 + f_3 * pc_x[k] * sdh_5[k];

        t_6[k] = f_0 * sph_6[k]
                 + f_6 * sdg0_6[k]
                 - f_7 * sdg1_6[k]
                 + f_3 * pc_x[k] * sdh_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sph_9, sdg0_9, sdg1_9, sdh_3, sdh_5, \
                         sdh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sdh_3[k];

        t_8[k] = f_3 * pc_y[k] * sdh_5[k];

        t_9[k] = f_0 * sph_9[k]
                 + f_6 * sdg0_9[k]
                 - f_7 * sdg1_9[k]
                 + f_3 * pc_x[k] * sdh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sph_10, sph_12, sdg0_10, sdg0_12, \
                         sdg1_10, sdg1_12, sdh_6, sdh_10, sdh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sph_10[k]
                  + f_8 * sdg0_10[k]
                  - f_9 * sdg1_10[k]
                  + f_3 * pc_x[k] * sdh_10[k];

        t_11[k] = f_3 * pc_z[k] * sdh_6[k];

        t_12[k] = f_0 * sph_12[k]
                  + f_8 * sdg0_12[k]
                  - f_9 * sdg1_12[k]
                  + f_3 * pc_x[k] * sdh_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pc_x, pc_y, sph_14, sph_15, sph_16, sdg0_14, \
                         sdg1_14, sdh_9, sdh_14, sdh_15, sdh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sdh_9[k];

        t_14[k] = f_0 * sph_14[k]
                  + f_8 * sdg0_14[k]
                  - f_9 * sdg1_14[k]
                  + f_3 * pc_x[k] * sdh_14[k];

        t_15[k] = f_0 * sph_15[k]
                  + f_3 * pc_x[k] * sdh_15[k];

        t_16[k] = f_0 * sph_16[k]
                  + f_3 * pc_x[k] * sdh_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_x, sph_17, sph_18, sph_19, sph_20, sdh_17, \
                         sdh_18, sdh_19, sdh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * sph_17[k]
                  + f_3 * pc_x[k] * sdh_17[k];

        t_18[k] = f_0 * sph_18[k]
                  + f_3 * pc_x[k] * sdh_18[k];

        t_19[k] = f_0 * sph_19[k]
                  + f_3 * pc_x[k] * sdh_19[k];

        t_20[k] = f_0 * sph_20[k]
                  + f_3 * pc_x[k] * sdh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, sdg0_10, sdg0_12, sdg0_13, \
                         sdg1_10, sdg1_12, sdg1_13, sdh_15, sdh_17, \
                         sdh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * sdg0_10[k]
                  - f_2 * sdg1_10[k]
                  + f_3 * pc_y[k] * sdh_15[k];

        t_22[k] = f_3 * pc_z[k] * sdh_15[k];

        t_23[k] = f_4 * sdg0_12[k]
                  - f_5 * sdg1_12[k]
                  + f_3 * pc_y[k] * sdh_17[k];

        t_24[k] = f_6 * sdg0_13[k]
                  - f_7 * sdg1_13[k]
                  + f_3 * pc_y[k] * sdh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_y, pc_y, pc_z, spi0_0, sph_0, \
                         spi1_0, sdg0_14, sdg1_14, sdh_19, sdh_20, \
                         sdh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * sdg0_14[k]
                  - f_9 * sdg1_14[k]
                  + f_3 * pc_y[k] * sdh_19[k];

        t_26[k] = f_3 * pc_y[k] * sdh_20[k];

        t_27[k] = f_1 * sdg0_14[k]
                  - f_2 * sdg1_14[k]
                  + f_3 * pc_z[k] * sdh_20[k];

        t_28[k] = pb_y[k] * spi0_0[k]
                  - f_10 * pc_y[k] * spi1_0[k];

        t_29[k] = f_11 * sph_0[k]
                  + f_3 * pc_y[k] * sdh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_x, pc_x, pc_y, pc_z, spi0_31, sph_2, sph_24, \
                         spi1_31, sdh_21, sdh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * sdh_21[k];

        t_31[k] = pb_x[k] * spi0_31[k]
                  + f_12 * sph_24[k]
                  - f_10 * pc_x[k] * spi1_31[k];

        t_32[k] = f_11 * sph_2[k]
                  + f_3 * pc_y[k] * sdh_23[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, pb_y, pc_x, pc_y, pc_z, spi0_5, spi0_34, \
                         sph_27, spi1_5, spi1_34, sdh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_y[k] * spi0_5[k]
                  - f_10 * pc_y[k] * spi1_5[k];

        t_34[k] = pb_x[k] * spi0_34[k]
                  + f_13 * sph_27[k]
                  - f_10 * pc_x[k] * spi1_34[k];

        t_35[k] = f_3 * pc_z[k] * sdh_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_x, pb_y, pc_x, pc_y, spi0_9, spi0_38, sph_5, \
                         sph_31, spi1_9, spi1_38, sdh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_11 * sph_5[k]
                  + f_3 * pc_y[k] * sdh_26[k];

        t_37[k] = pb_y[k] * spi0_9[k]
                  - f_10 * pc_y[k] * spi1_9[k];

        t_38[k] = pb_x[k] * spi0_38[k]
                  + f_0 * sph_31[k]
                  - f_10 * pc_x[k] * spi1_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, pc_x, pc_y, pc_z, spi0_40, sph_9, sph_33, \
                         spi1_40, sdh_27, sdh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_3 * pc_z[k] * sdh_27[k];

        t_40[k] = pb_x[k] * spi0_40[k]
                  + f_0 * sph_33[k]
                  - f_10 * pc_x[k] * spi1_40[k];

        t_41[k] = f_11 * sph_9[k]
                  + f_3 * pc_y[k] * sdh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_y, pc_x, pc_y, spi0_14, sph_36, sph_37, \
                         sph_38, spi1_14, sdh_36, sdh_37, sdh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * spi0_14[k]
                  - f_10 * pc_y[k] * spi1_14[k];

        t_43[k] = f_11 * sph_36[k]
                  + f_3 * pc_x[k] * sdh_36[k];

        t_44[k] = f_11 * sph_37[k]
                  + f_3 * pc_x[k] * sdh_37[k];

        t_45[k] = f_11 * sph_38[k]
                  + f_3 * pc_x[k] * sdh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_x, pc_x, spi0_49, sph_39, sph_40, sph_41, \
                         spi1_49, sdh_39, sdh_40, sdh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_11 * sph_39[k]
                  + f_3 * pc_x[k] * sdh_39[k];

        t_47[k] = f_11 * sph_40[k]
                  + f_3 * pc_x[k] * sdh_40[k];

        t_48[k] = f_11 * sph_41[k]
                  + f_3 * pc_x[k] * sdh_41[k];

        t_49[k] = pb_x[k] * spi0_49[k]
                  - f_10 * pc_x[k] * spi1_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pb_x, pc_x, pc_z, spi0_51, spi0_52, spi0_53, \
                         spi1_51, spi1_52, spi1_53, sdh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * sdh_36[k];

        t_51[k] = pb_x[k] * spi0_51[k]
                  - f_10 * pc_x[k] * spi1_51[k];

        t_52[k] = pb_x[k] * spi0_52[k]
                  - f_10 * pc_x[k] * spi1_52[k];

        t_53[k] = pb_x[k] * spi0_53[k]
                  - f_10 * pc_x[k] * spi1_53[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pb_x, pb_z, pc_x, pc_y, pc_z, spi0_0, \
                         spi0_55, sph_20, spi1_0, spi1_55, sdh_41, \
                         sdh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * sph_20[k]
                  + f_3 * pc_y[k] * sdh_41[k];

        t_55[k] = pb_x[k] * spi0_55[k]
                  - f_10 * pc_x[k] * spi1_55[k];

        t_56[k] = pb_z[k] * spi0_0[k]
                  - f_10 * pc_z[k] * spi1_0[k];

        t_57[k] = f_3 * pc_y[k] * sdh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_z, pc_y, pc_z, spi0_3, sph_0, spi1_3, sdh_42, \
                         sdh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_11 * sph_0[k]
                  + f_3 * pc_z[k] * sdh_42[k];

        t_59[k] = pb_z[k] * spi0_3[k]
                  - f_10 * pc_z[k] * spi1_3[k];

        t_60[k] = f_3 * pc_y[k] * sdh_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_x, pb_z, pc_x, pc_z, spi0_6, spi0_61, sph_3, \
                         sph_47, spi1_6, spi1_61, sdh_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_x[k] * spi0_61[k]
                  + f_12 * sph_47[k]
                  - f_10 * pc_x[k] * spi1_61[k];

        t_62[k] = pb_z[k] * spi0_6[k]
                  - f_10 * pc_z[k] * spi1_6[k];

        t_63[k] = f_11 * sph_3[k]
                  + f_3 * pc_z[k] * sdh_45[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, pb_z, pc_x, pc_y, pc_z, spi0_10, spi0_65, \
                         sph_51, spi1_10, spi1_65, sdh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * pc_y[k] * sdh_47[k];

        t_65[k] = pb_x[k] * spi0_65[k]
                  + f_13 * sph_51[k]
                  - f_10 * pc_x[k] * spi1_65[k];

        t_66[k] = pb_z[k] * spi0_10[k]
                  - f_10 * pc_z[k] * spi1_10[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_x, pc_x, pc_y, pc_z, spi0_68, sph_6, sph_54, \
                         spi1_68, sdh_48, sdh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_11 * sph_6[k]
                  + f_3 * pc_z[k] * sdh_48[k];

        t_68[k] = pb_x[k] * spi0_68[k]
                  + f_0 * sph_54[k]
                  - f_10 * pc_x[k] * spi1_68[k];

        t_69[k] = f_3 * pc_y[k] * sdh_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_x, pc_x, spi0_70, sph_56, sph_57, sph_58, \
                         sph_59, spi1_70, sdh_57, sdh_58, sdh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pb_x[k] * spi0_70[k]
                  + f_0 * sph_56[k]
                  - f_10 * pc_x[k] * spi1_70[k];

        t_71[k] = f_11 * sph_57[k]
                  + f_3 * pc_x[k] * sdh_57[k];

        t_72[k] = f_11 * sph_58[k]
                  + f_3 * pc_x[k] * sdh_58[k];

        t_73[k] = f_11 * sph_59[k]
                  + f_3 * pc_x[k] * sdh_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_x, pc_x, spi0_77, sph_60, sph_61, sph_62, \
                         spi1_77, sdh_60, sdh_61, sdh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_11 * sph_60[k]
                  + f_3 * pc_x[k] * sdh_60[k];

        t_75[k] = f_11 * sph_61[k]
                  + f_3 * pc_x[k] * sdh_61[k];

        t_76[k] = f_11 * sph_62[k]
                  + f_3 * pc_x[k] * sdh_62[k];

        t_77[k] = pb_x[k] * spi0_77[k]
                  - f_10 * pc_x[k] * spi1_77[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_x, pc_x, pc_z, spi0_79, spi0_80, spi0_81, \
                         sph_15, spi1_79, spi1_80, spi1_81, sdh_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_11 * sph_15[k]
                  + f_3 * pc_z[k] * sdh_57[k];

        t_79[k] = pb_x[k] * spi0_79[k]
                  - f_10 * pc_x[k] * spi1_79[k];

        t_80[k] = pb_x[k] * spi0_80[k]
                  - f_10 * pc_x[k] * spi1_80[k];

        t_81[k] = pb_x[k] * spi0_81[k]
                  - f_10 * pc_x[k] * spi1_81[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pb_x, pc_x, pc_y, pc_z, spi0_83, \
                         sph_21, spi1_83, sdg0_45, sdg1_45, sdh_62, \
                         sdh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * pc_y[k] * sdh_62[k];

        t_83[k] = pb_x[k] * spi0_83[k]
                  - f_10 * pc_x[k] * spi1_83[k];

        t_84[k] = f_1 * sdg0_45[k]
                  - f_2 * sdg1_45[k]
                  + f_3 * pc_x[k] * sdh_63[k];

        t_85[k] = f_0 * sph_21[k]
                  + f_3 * pc_y[k] * sdh_63[k];

        t_86[k] = f_3 * pc_z[k] * sdh_63[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pc_x, pc_y, sph_23, sdg0_48, sdg0_50, sdg1_48, \
                         sdg1_50, sdh_65, sdh_66, sdh_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * sdg0_48[k]
                  - f_5 * sdg1_48[k]
                  + f_3 * pc_x[k] * sdh_66[k];

        t_88[k] = f_0 * sph_23[k]
                  + f_3 * pc_y[k] * sdh_65[k];

        t_89[k] = f_4 * sdg0_50[k]
                  - f_5 * sdg1_50[k]
                  + f_3 * pc_x[k] * sdh_68[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, pc_z, sph_26, sdg0_51, sdg0_54, \
                         sdg1_51, sdg1_54, sdh_66, sdh_68, sdh_69, \
                         sdh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_6 * sdg0_51[k]
                  - f_7 * sdg1_51[k]
                  + f_3 * pc_x[k] * sdh_69[k];

        t_91[k] = f_3 * pc_z[k] * sdh_66[k];

        t_92[k] = f_0 * sph_26[k]
                  + f_3 * pc_y[k] * sdh_68[k];

        t_93[k] = f_6 * sdg0_54[k]
                  - f_7 * sdg1_54[k]
                  + f_3 * pc_x[k] * sdh_72[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pc_x, pc_y, pc_z, sph_30, sdg0_55, sdg0_57, \
                         sdg1_55, sdg1_57, sdh_69, sdh_72, sdh_73, \
                         sdh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_8 * sdg0_55[k]
                  - f_9 * sdg1_55[k]
                  + f_3 * pc_x[k] * sdh_73[k];

        t_95[k] = f_3 * pc_z[k] * sdh_69[k];

        t_96[k] = f_8 * sdg0_57[k]
                  - f_9 * sdg1_57[k]
                  + f_3 * pc_x[k] * sdh_75[k];

        t_97[k] = f_0 * sph_30[k]
                  + f_3 * pc_y[k] * sdh_72[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, pc_x, sdg0_59, sdg1_59, \
                         sdh_77, sdh_78, sdh_79, sdh_80, sdh_81, \
                         sdh_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_8 * sdg0_59[k]
                  - f_9 * sdg1_59[k]
                  + f_3 * pc_x[k] * sdh_77[k];

        t_99[k] = f_3 * pc_x[k] * sdh_78[k];

        t_100[k] = f_3 * pc_x[k] * sdh_79[k];

        t_101[k] = f_3 * pc_x[k] * sdh_80[k];

        t_102[k] = f_3 * pc_x[k] * sdh_81[k];

        t_103[k] = f_3 * pc_x[k] * sdh_82[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pc_x, pc_y, pc_z, sph_36, sph_38, \
                         sdg0_55, sdg0_57, sdg1_55, sdg1_57, sdh_78, sdh_80, \
                         sdh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_3 * pc_x[k] * sdh_83[k];

        t_105[k] = f_0 * sph_36[k]
                   + f_1 * sdg0_55[k]
                   - f_2 * sdg1_55[k]
                   + f_3 * pc_y[k] * sdh_78[k];

        t_106[k] = f_3 * pc_z[k] * sdh_78[k];

        t_107[k] = f_0 * sph_38[k]
                   + f_4 * sdg0_57[k]
                   - f_5 * sdg1_57[k]
                   + f_3 * pc_y[k] * sdh_80[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pc_y, pc_z, sph_39, sph_40, sph_41, \
                         sdg0_58, sdg0_59, sdg1_58, sdg1_59, sdh_81, sdh_82, \
                         sdh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_0 * sph_39[k]
                   + f_6 * sdg0_58[k]
                   - f_7 * sdg1_58[k]
                   + f_3 * pc_y[k] * sdh_81[k];

        t_109[k] = f_0 * sph_40[k]
                   + f_8 * sdg0_59[k]
                   - f_9 * sdg1_59[k]
                   + f_3 * pc_y[k] * sdh_82[k];

        t_110[k] = f_0 * sph_41[k]
                   + f_3 * pc_y[k] * sdh_83[k];

        t_111[k] = f_1 * sdg0_59[k]
                   - f_2 * sdg1_59[k]
                   + f_3 * pc_z[k] * sdh_83[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pb_y, pb_z, pc_y, pc_z, spi0_31, spi0_56, \
                         sph_21, sph_42, spi1_31, spi1_56, sdh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = pb_y[k] * spi0_56[k]
                   - f_10 * pc_y[k] * spi1_56[k];

        t_113[k] = f_11 * sph_42[k]
                   + f_3 * pc_y[k] * sdh_84[k];

        t_114[k] = f_11 * sph_21[k]
                   + f_3 * pc_z[k] * sdh_84[k];

        t_115[k] = pb_z[k] * spi0_31[k]
                   - f_10 * pc_z[k] * spi1_31[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pb_z, pc_y, pc_z, spi0_34, spi0_61, \
                         sph_24, sph_44, spi1_34, spi1_61, sdh_86, \
                         sdh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_11 * sph_44[k]
                   + f_3 * pc_y[k] * sdh_86[k];

        t_117[k] = pb_y[k] * spi0_61[k]
                   - f_10 * pc_y[k] * spi1_61[k];

        t_118[k] = pb_z[k] * spi0_34[k]
                   - f_10 * pc_z[k] * spi1_34[k];

        t_119[k] = f_11 * sph_24[k]
                   + f_3 * pc_z[k] * sdh_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pb_z, pc_y, pc_z, spi0_38, spi0_65, \
                         sph_27, sph_47, spi1_38, spi1_65, sdh_89, \
                         sdh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_11 * sph_47[k]
                   + f_3 * pc_y[k] * sdh_89[k];

        t_121[k] = pb_y[k] * spi0_65[k]
                   - f_10 * pc_y[k] * spi1_65[k];

        t_122[k] = pb_z[k] * spi0_38[k]
                   - f_10 * pc_z[k] * spi1_38[k];

        t_123[k] = f_11 * sph_27[k]
                   + f_3 * pc_z[k] * sdh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_y, pc_x, pc_y, spi0_70, sph_51, \
                         spi1_70, sdg0_72, sdg1_72, sdh_93, sdh_96, \
                         sdh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_8 * sdg0_72[k]
                   - f_9 * sdg1_72[k]
                   + f_3 * pc_x[k] * sdh_96[k];

        t_125[k] = f_11 * sph_51[k]
                   + f_3 * pc_y[k] * sdh_93[k];

        t_126[k] = pb_y[k] * spi0_70[k]
                   - f_10 * pc_y[k] * spi1_70[k];

        t_127[k] = f_3 * pc_x[k] * sdh_99[k];
    }
}

static auto
compute_prim_sdi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t spi0,
                                                          const size_t sph, const size_t spi1,
                                                          const size_t sdg0, const size_t sdg1,
                                                          const size_t sdh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spi0_49 = buffer.data(spi0 + 49);
    const auto *spi0_79 = buffer.data(spi0 + 79);
    const auto *spi0_80 = buffer.data(spi0 + 80);
    const auto *spi0_81 = buffer.data(spi0 + 81);
    const auto *spi0_83 = buffer.data(spi0 + 83);

    const auto *sph_36 = buffer.data(sph + 36);
    const auto *sph_42 = buffer.data(sph + 42);
    const auto *sph_45 = buffer.data(sph + 45);
    const auto *sph_48 = buffer.data(sph + 48);
    const auto *sph_57 = buffer.data(sph + 57);
    const auto *sph_59 = buffer.data(sph + 59);
    const auto *sph_60 = buffer.data(sph + 60);
    const auto *sph_61 = buffer.data(sph + 61);
    const auto *sph_62 = buffer.data(sph + 62);

    const auto *spi1_49 = buffer.data(spi1 + 49);
    const auto *spi1_79 = buffer.data(spi1 + 79);
    const auto *spi1_80 = buffer.data(spi1 + 80);
    const auto *spi1_81 = buffer.data(spi1 + 81);
    const auto *spi1_83 = buffer.data(spi1 + 83);

    const auto *sdg0_75 = buffer.data(sdg0 + 75);
    const auto *sdg0_78 = buffer.data(sdg0 + 78);
    const auto *sdg0_80 = buffer.data(sdg0 + 80);
    const auto *sdg0_81 = buffer.data(sdg0 + 81);
    const auto *sdg0_84 = buffer.data(sdg0 + 84);
    const auto *sdg0_85 = buffer.data(sdg0 + 85);
    const auto *sdg0_87 = buffer.data(sdg0 + 87);
    const auto *sdg0_88 = buffer.data(sdg0 + 88);
    const auto *sdg0_89 = buffer.data(sdg0 + 89);

    const auto *sdg1_75 = buffer.data(sdg1 + 75);
    const auto *sdg1_78 = buffer.data(sdg1 + 78);
    const auto *sdg1_80 = buffer.data(sdg1 + 80);
    const auto *sdg1_81 = buffer.data(sdg1 + 81);
    const auto *sdg1_84 = buffer.data(sdg1 + 84);
    const auto *sdg1_85 = buffer.data(sdg1 + 85);
    const auto *sdg1_87 = buffer.data(sdg1 + 87);
    const auto *sdg1_88 = buffer.data(sdg1 + 88);
    const auto *sdg1_89 = buffer.data(sdg1 + 89);

    const auto *sdh_99 = buffer.data(sdh + 99);
    const auto *sdh_100 = buffer.data(sdh + 100);
    const auto *sdh_101 = buffer.data(sdh + 101);
    const auto *sdh_102 = buffer.data(sdh + 102);
    const auto *sdh_103 = buffer.data(sdh + 103);
    const auto *sdh_104 = buffer.data(sdh + 104);
    const auto *sdh_105 = buffer.data(sdh + 105);
    const auto *sdh_107 = buffer.data(sdh + 107);
    const auto *sdh_108 = buffer.data(sdh + 108);
    const auto *sdh_110 = buffer.data(sdh + 110);
    const auto *sdh_111 = buffer.data(sdh + 111);
    const auto *sdh_114 = buffer.data(sdh + 114);
    const auto *sdh_115 = buffer.data(sdh + 115);
    const auto *sdh_117 = buffer.data(sdh + 117);
    const auto *sdh_119 = buffer.data(sdh + 119);
    const auto *sdh_120 = buffer.data(sdh + 120);
    const auto *sdh_121 = buffer.data(sdh + 121);
    const auto *sdh_122 = buffer.data(sdh + 122);
    const auto *sdh_123 = buffer.data(sdh + 123);
    const auto *sdh_124 = buffer.data(sdh + 124);
    const auto *sdh_125 = buffer.data(sdh + 125);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, pb_z, pc_x, pc_z, spi0_49, \
                         spi1_49, sdh_100, sdh_101, sdh_102, sdh_103, \
                         sdh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_3 * pc_x[k] * sdh_100[k];

        t_129[k] = f_3 * pc_x[k] * sdh_101[k];

        t_130[k] = f_3 * pc_x[k] * sdh_102[k];

        t_131[k] = f_3 * pc_x[k] * sdh_103[k];

        t_132[k] = f_3 * pc_x[k] * sdh_104[k];

        t_133[k] = pb_z[k] * spi0_49[k]
                   - f_10 * pc_z[k] * spi1_49[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pb_y, pc_y, pc_z, spi0_79, spi0_80, sph_36, \
                         sph_59, sph_60, spi1_79, spi1_80, sdh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * sph_36[k]
                   + f_3 * pc_z[k] * sdh_99[k];

        t_135[k] = pb_y[k] * spi0_79[k]
                   + f_12 * sph_59[k]
                   - f_10 * pc_y[k] * spi1_79[k];

        t_136[k] = pb_y[k] * spi0_80[k]
                   + f_13 * sph_60[k]
                   - f_10 * pc_y[k] * spi1_80[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pb_y, pc_y, spi0_81, spi0_83, sph_61, sph_62, \
                         spi1_81, spi1_83, sdh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = pb_y[k] * spi0_81[k]
                   + f_0 * sph_61[k]
                   - f_10 * pc_y[k] * spi1_81[k];

        t_138[k] = f_11 * sph_62[k]
                   + f_3 * pc_y[k] * sdh_104[k];

        t_139[k] = pb_y[k] * spi0_83[k]
                   - f_10 * pc_y[k] * spi1_83[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pc_x, pc_y, pc_z, sph_42, sdg0_75, \
                         sdg0_78, sdg1_75, sdg1_78, sdh_105, sdh_107, \
                         sdh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_1 * sdg0_75[k]
                   - f_2 * sdg1_75[k]
                   + f_3 * pc_x[k] * sdh_105[k];

        t_141[k] = f_3 * pc_y[k] * sdh_105[k];

        t_142[k] = f_0 * sph_42[k]
                   + f_3 * pc_z[k] * sdh_105[k];

        t_143[k] = f_4 * sdg0_78[k]
                   - f_5 * sdg1_78[k]
                   + f_3 * pc_x[k] * sdh_108[k];

        t_144[k] = f_3 * pc_y[k] * sdh_107[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pc_x, pc_y, pc_z, sph_45, sdg0_80, \
                         sdg0_81, sdg1_80, sdg1_81, sdh_108, sdh_110, \
                         sdh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_4 * sdg0_80[k]
                   - f_5 * sdg1_80[k]
                   + f_3 * pc_x[k] * sdh_110[k];

        t_146[k] = f_6 * sdg0_81[k]
                   - f_7 * sdg1_81[k]
                   + f_3 * pc_x[k] * sdh_111[k];

        t_147[k] = f_0 * sph_45[k]
                   + f_3 * pc_z[k] * sdh_108[k];

        t_148[k] = f_3 * pc_y[k] * sdh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_z, sph_48, sdg0_84, sdg0_85, sdg1_84, \
                         sdg1_85, sdh_111, sdh_114, sdh_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_6 * sdg0_84[k]
                   - f_7 * sdg1_84[k]
                   + f_3 * pc_x[k] * sdh_114[k];

        t_150[k] = f_8 * sdg0_85[k]
                   - f_9 * sdg1_85[k]
                   + f_3 * pc_x[k] * sdh_115[k];

        t_151[k] = f_0 * sph_48[k]
                   + f_3 * pc_z[k] * sdh_111[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, pc_x, pc_y, sdg0_87, sdg0_89, \
                         sdg1_87, sdg1_89, sdh_114, sdh_117, sdh_119, sdh_120, \
                         sdh_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_8 * sdg0_87[k]
                   - f_9 * sdg1_87[k]
                   + f_3 * pc_x[k] * sdh_117[k];

        t_153[k] = f_3 * pc_y[k] * sdh_114[k];

        t_154[k] = f_8 * sdg0_89[k]
                   - f_9 * sdg1_89[k]
                   + f_3 * pc_x[k] * sdh_119[k];

        t_155[k] = f_3 * pc_x[k] * sdh_120[k];

        t_156[k] = f_3 * pc_x[k] * sdh_121[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pc_x, pc_y, sdg0_85, sdg1_85, \
                         sdh_120, sdh_122, sdh_123, sdh_124, sdh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_3 * pc_x[k] * sdh_122[k];

        t_158[k] = f_3 * pc_x[k] * sdh_123[k];

        t_159[k] = f_3 * pc_x[k] * sdh_124[k];

        t_160[k] = f_3 * pc_x[k] * sdh_125[k];

        t_161[k] = f_1 * sdg0_85[k]
                   - f_2 * sdg1_85[k]
                   + f_3 * pc_y[k] * sdh_120[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pc_y, pc_z, sph_57, sdg0_87, sdg0_88, sdg1_87, \
                         sdg1_88, sdh_120, sdh_122, sdh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_0 * sph_57[k]
                   + f_3 * pc_z[k] * sdh_120[k];

        t_163[k] = f_4 * sdg0_87[k]
                   - f_5 * sdg1_87[k]
                   + f_3 * pc_y[k] * sdh_122[k];

        t_164[k] = f_6 * sdg0_88[k]
                   - f_7 * sdg1_88[k]
                   + f_3 * pc_y[k] * sdh_123[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pc_y, pc_z, sph_62, sdg0_89, sdg1_89, sdh_124, \
                         sdh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_8 * sdg0_89[k]
                   - f_9 * sdg1_89[k]
                   + f_3 * pc_y[k] * sdh_124[k];

        t_166[k] = f_3 * pc_y[k] * sdh_125[k];

        t_167[k] = f_0 * sph_62[k]
                   + f_1 * sdg0_89[k]
                   - f_2 * sdg1_89[k]
                   + f_3 * pc_z[k] * sdh_125[k];
    }
}

auto
compute_prim_sdi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t spi0, const size_t sph,
                                                   const size_t spi1, const size_t sdg0,
                                                   const size_t sdg1, const size_t sdh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sdi_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, spi0, sph,
                                                              spi1, sdg0, sdg1, sdh, ncols,
                                                              gamma, p, q);

    compute_prim_sdi_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, spi0, sph,
                                                              spi1, sdg0, sdg1, sdh, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
