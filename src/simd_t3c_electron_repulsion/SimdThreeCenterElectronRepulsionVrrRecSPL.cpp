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


#include "SimdThreeCenterElectronRepulsionVrrRecSPL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_spl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t ssl0,
                                                          const size_t ssk, const size_t ssl1,
                                                          const size_t spi0, const size_t spi1,
                                                          const size_t spk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = gamma / q;
    const auto f_2 = p / q;
    const auto f_3 = 3.0 / q;
    const auto f_4 = 2.5 / q;
    const auto f_5 = 2.0 / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 1.0 / q;
    const auto f_8 = 0.5 / q;
    const auto f_9 = 2.5 / gamma;
    const auto f_10 = 2.5 * p / (gamma * q);
    const auto f_11 = 2.0 / gamma;
    const auto f_12 = 2.0 * p / (gamma * q);
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 1.0 / gamma;
    const auto f_16 = p / (gamma * q);
    const auto f_17 = 0.5 / gamma;
    const auto f_18 = 0.5 * p / (gamma * q);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssl0_0 = buffer.data(ssl0 + 0);
    const auto *ssl0_3 = buffer.data(ssl0 + 3);
    const auto *ssl0_5 = buffer.data(ssl0 + 5);
    const auto *ssl0_6 = buffer.data(ssl0 + 6);
    const auto *ssl0_9 = buffer.data(ssl0 + 9);
    const auto *ssl0_10 = buffer.data(ssl0 + 10);
    const auto *ssl0_12 = buffer.data(ssl0 + 12);
    const auto *ssl0_14 = buffer.data(ssl0 + 14);
    const auto *ssl0_15 = buffer.data(ssl0 + 15);
    const auto *ssl0_17 = buffer.data(ssl0 + 17);
    const auto *ssl0_18 = buffer.data(ssl0 + 18);
    const auto *ssl0_20 = buffer.data(ssl0 + 20);
    const auto *ssl0_21 = buffer.data(ssl0 + 21);
    const auto *ssl0_23 = buffer.data(ssl0 + 23);
    const auto *ssl0_24 = buffer.data(ssl0 + 24);
    const auto *ssl0_25 = buffer.data(ssl0 + 25);
    const auto *ssl0_27 = buffer.data(ssl0 + 27);
    const auto *ssl0_36 = buffer.data(ssl0 + 36);
    const auto *ssl0_38 = buffer.data(ssl0 + 38);
    const auto *ssl0_39 = buffer.data(ssl0 + 39);
    const auto *ssl0_40 = buffer.data(ssl0 + 40);
    const auto *ssl0_41 = buffer.data(ssl0 + 41);
    const auto *ssl0_42 = buffer.data(ssl0 + 42);
    const auto *ssl0_44 = buffer.data(ssl0 + 44);

    const auto *ssk_0 = buffer.data(ssk + 0);
    const auto *ssk_2 = buffer.data(ssk + 2);
    const auto *ssk_3 = buffer.data(ssk + 3);
    const auto *ssk_5 = buffer.data(ssk + 5);
    const auto *ssk_6 = buffer.data(ssk + 6);
    const auto *ssk_9 = buffer.data(ssk + 9);
    const auto *ssk_10 = buffer.data(ssk + 10);
    const auto *ssk_12 = buffer.data(ssk + 12);
    const auto *ssk_14 = buffer.data(ssk + 14);
    const auto *ssk_15 = buffer.data(ssk + 15);
    const auto *ssk_17 = buffer.data(ssk + 17);
    const auto *ssk_18 = buffer.data(ssk + 18);
    const auto *ssk_20 = buffer.data(ssk + 20);
    const auto *ssk_21 = buffer.data(ssk + 21);
    const auto *ssk_23 = buffer.data(ssk + 23);
    const auto *ssk_24 = buffer.data(ssk + 24);
    const auto *ssk_25 = buffer.data(ssk + 25);
    const auto *ssk_27 = buffer.data(ssk + 27);
    const auto *ssk_28 = buffer.data(ssk + 28);
    const auto *ssk_29 = buffer.data(ssk + 29);
    const auto *ssk_30 = buffer.data(ssk + 30);
    const auto *ssk_31 = buffer.data(ssk + 31);
    const auto *ssk_32 = buffer.data(ssk + 32);
    const auto *ssk_33 = buffer.data(ssk + 33);
    const auto *ssk_34 = buffer.data(ssk + 34);
    const auto *ssk_35 = buffer.data(ssk + 35);

    const auto *ssl1_0 = buffer.data(ssl1 + 0);
    const auto *ssl1_3 = buffer.data(ssl1 + 3);
    const auto *ssl1_5 = buffer.data(ssl1 + 5);
    const auto *ssl1_6 = buffer.data(ssl1 + 6);
    const auto *ssl1_9 = buffer.data(ssl1 + 9);
    const auto *ssl1_10 = buffer.data(ssl1 + 10);
    const auto *ssl1_12 = buffer.data(ssl1 + 12);
    const auto *ssl1_14 = buffer.data(ssl1 + 14);
    const auto *ssl1_15 = buffer.data(ssl1 + 15);
    const auto *ssl1_17 = buffer.data(ssl1 + 17);
    const auto *ssl1_18 = buffer.data(ssl1 + 18);
    const auto *ssl1_20 = buffer.data(ssl1 + 20);
    const auto *ssl1_21 = buffer.data(ssl1 + 21);
    const auto *ssl1_23 = buffer.data(ssl1 + 23);
    const auto *ssl1_24 = buffer.data(ssl1 + 24);
    const auto *ssl1_25 = buffer.data(ssl1 + 25);
    const auto *ssl1_27 = buffer.data(ssl1 + 27);
    const auto *ssl1_36 = buffer.data(ssl1 + 36);
    const auto *ssl1_38 = buffer.data(ssl1 + 38);
    const auto *ssl1_39 = buffer.data(ssl1 + 39);
    const auto *ssl1_40 = buffer.data(ssl1 + 40);
    const auto *ssl1_41 = buffer.data(ssl1 + 41);
    const auto *ssl1_42 = buffer.data(ssl1 + 42);
    const auto *ssl1_44 = buffer.data(ssl1 + 44);

    const auto *spi0_31 = buffer.data(spi0 + 31);
    const auto *spi0_34 = buffer.data(spi0 + 34);
    const auto *spi0_38 = buffer.data(spi0 + 38);
    const auto *spi0_40 = buffer.data(spi0 + 40);
    const auto *spi0_43 = buffer.data(spi0 + 43);
    const auto *spi0_45 = buffer.data(spi0 + 45);
    const auto *spi0_46 = buffer.data(spi0 + 46);
    const auto *spi0_49 = buffer.data(spi0 + 49);
    const auto *spi0_51 = buffer.data(spi0 + 51);
    const auto *spi0_52 = buffer.data(spi0 + 52);
    const auto *spi0_53 = buffer.data(spi0 + 53);
    const auto *spi0_61 = buffer.data(spi0 + 61);
    const auto *spi0_65 = buffer.data(spi0 + 65);
    const auto *spi0_68 = buffer.data(spi0 + 68);
    const auto *spi0_70 = buffer.data(spi0 + 70);
    const auto *spi0_73 = buffer.data(spi0 + 73);
    const auto *spi0_74 = buffer.data(spi0 + 74);
    const auto *spi0_76 = buffer.data(spi0 + 76);
    const auto *spi0_79 = buffer.data(spi0 + 79);
    const auto *spi0_80 = buffer.data(spi0 + 80);
    const auto *spi0_81 = buffer.data(spi0 + 81);
    const auto *spi0_83 = buffer.data(spi0 + 83);

    const auto *spi1_31 = buffer.data(spi1 + 31);
    const auto *spi1_34 = buffer.data(spi1 + 34);
    const auto *spi1_38 = buffer.data(spi1 + 38);
    const auto *spi1_40 = buffer.data(spi1 + 40);
    const auto *spi1_43 = buffer.data(spi1 + 43);
    const auto *spi1_45 = buffer.data(spi1 + 45);
    const auto *spi1_46 = buffer.data(spi1 + 46);
    const auto *spi1_49 = buffer.data(spi1 + 49);
    const auto *spi1_51 = buffer.data(spi1 + 51);
    const auto *spi1_52 = buffer.data(spi1 + 52);
    const auto *spi1_53 = buffer.data(spi1 + 53);
    const auto *spi1_61 = buffer.data(spi1 + 61);
    const auto *spi1_65 = buffer.data(spi1 + 65);
    const auto *spi1_68 = buffer.data(spi1 + 68);
    const auto *spi1_70 = buffer.data(spi1 + 70);
    const auto *spi1_73 = buffer.data(spi1 + 73);
    const auto *spi1_74 = buffer.data(spi1 + 74);
    const auto *spi1_76 = buffer.data(spi1 + 76);
    const auto *spi1_79 = buffer.data(spi1 + 79);
    const auto *spi1_80 = buffer.data(spi1 + 80);
    const auto *spi1_81 = buffer.data(spi1 + 81);
    const auto *spi1_83 = buffer.data(spi1 + 83);

    const auto *spk_0 = buffer.data(spk + 0);
    const auto *spk_2 = buffer.data(spk + 2);
    const auto *spk_3 = buffer.data(spk + 3);
    const auto *spk_5 = buffer.data(spk + 5);
    const auto *spk_6 = buffer.data(spk + 6);
    const auto *spk_9 = buffer.data(spk + 9);
    const auto *spk_10 = buffer.data(spk + 10);
    const auto *spk_14 = buffer.data(spk + 14);
    const auto *spk_15 = buffer.data(spk + 15);
    const auto *spk_20 = buffer.data(spk + 20);
    const auto *spk_28 = buffer.data(spk + 28);
    const auto *spk_29 = buffer.data(spk + 29);
    const auto *spk_30 = buffer.data(spk + 30);
    const auto *spk_31 = buffer.data(spk + 31);
    const auto *spk_32 = buffer.data(spk + 32);
    const auto *spk_33 = buffer.data(spk + 33);
    const auto *spk_34 = buffer.data(spk + 34);
    const auto *spk_35 = buffer.data(spk + 35);
    const auto *spk_36 = buffer.data(spk + 36);
    const auto *spk_38 = buffer.data(spk + 38);
    const auto *spk_39 = buffer.data(spk + 39);
    const auto *spk_41 = buffer.data(spk + 41);
    const auto *spk_42 = buffer.data(spk + 42);
    const auto *spk_45 = buffer.data(spk + 45);
    const auto *spk_46 = buffer.data(spk + 46);
    const auto *spk_48 = buffer.data(spk + 48);
    const auto *spk_50 = buffer.data(spk + 50);
    const auto *spk_51 = buffer.data(spk + 51);
    const auto *spk_53 = buffer.data(spk + 53);
    const auto *spk_54 = buffer.data(spk + 54);
    const auto *spk_56 = buffer.data(spk + 56);
    const auto *spk_57 = buffer.data(spk + 57);
    const auto *spk_59 = buffer.data(spk + 59);
    const auto *spk_60 = buffer.data(spk + 60);
    const auto *spk_61 = buffer.data(spk + 61);
    const auto *spk_64 = buffer.data(spk + 64);
    const auto *spk_65 = buffer.data(spk + 65);
    const auto *spk_66 = buffer.data(spk + 66);
    const auto *spk_67 = buffer.data(spk + 67);
    const auto *spk_68 = buffer.data(spk + 68);
    const auto *spk_69 = buffer.data(spk + 69);
    const auto *spk_70 = buffer.data(spk + 70);
    const auto *spk_71 = buffer.data(spk + 71);
    const auto *spk_72 = buffer.data(spk + 72);
    const auto *spk_74 = buffer.data(spk + 74);
    const auto *spk_75 = buffer.data(spk + 75);
    const auto *spk_77 = buffer.data(spk + 77);
    const auto *spk_78 = buffer.data(spk + 78);
    const auto *spk_81 = buffer.data(spk + 81);
    const auto *spk_82 = buffer.data(spk + 82);
    const auto *spk_84 = buffer.data(spk + 84);
    const auto *spk_86 = buffer.data(spk + 86);
    const auto *spk_87 = buffer.data(spk + 87);
    const auto *spk_89 = buffer.data(spk + 89);
    const auto *spk_90 = buffer.data(spk + 90);
    const auto *spk_92 = buffer.data(spk + 92);
    const auto *spk_95 = buffer.data(spk + 95);
    const auto *spk_96 = buffer.data(spk + 96);
    const auto *spk_97 = buffer.data(spk + 97);
    const auto *spk_99 = buffer.data(spk + 99);
    const auto *spk_100 = buffer.data(spk + 100);
    const auto *spk_101 = buffer.data(spk + 101);
    const auto *spk_102 = buffer.data(spk + 102);
    const auto *spk_103 = buffer.data(spk + 103);
    const auto *spk_104 = buffer.data(spk + 104);
    const auto *spk_105 = buffer.data(spk + 105);
    const auto *spk_106 = buffer.data(spk + 106);
    const auto *spk_107 = buffer.data(spk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pc_x, pc_y, pc_z, ssl0_0, ssl0_3, ssk_0, \
                         ssk_3, ssl1_0, ssl1_3, spk_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = pb_x[k] * ssl0_0[k]
                 + f_0 * ssk_0[k]
                 - f_1 * pc_x[k] * ssl1_0[k];

        t_1[k] = f_2 * pc_y[k] * spk_0[k];

        t_2[k] = f_2 * pc_z[k] * spk_0[k];

        t_3[k] = pb_x[k] * ssl0_3[k]
                 + f_3 * ssk_3[k]
                 - f_1 * pc_x[k] * ssl1_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pc_x, pc_y, pc_z, ssl0_5, ssl0_6, ssk_5, \
                         ssk_6, ssl1_5, ssl1_6, spk_2, spk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * pc_y[k] * spk_2[k];

        t_5[k] = pb_x[k] * ssl0_5[k]
                 + f_3 * ssk_5[k]
                 - f_1 * pc_x[k] * ssl1_5[k];

        t_6[k] = pb_x[k] * ssl0_6[k]
                 + f_4 * ssk_6[k]
                 - f_1 * pc_x[k] * ssl1_6[k];

        t_7[k] = f_2 * pc_z[k] * spk_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pc_x, pc_y, pc_z, ssl0_9, ssl0_10, ssk_9, \
                         ssk_10, ssl1_9, ssl1_10, spk_5, spk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * pc_y[k] * spk_5[k];

        t_9[k] = pb_x[k] * ssl0_9[k]
                 + f_4 * ssk_9[k]
                 - f_1 * pc_x[k] * ssl1_9[k];

        t_10[k] = pb_x[k] * ssl0_10[k]
                  + f_5 * ssk_10[k]
                  - f_1 * pc_x[k] * ssl1_10[k];

        t_11[k] = f_2 * pc_z[k] * spk_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_x, pc_x, pc_y, ssl0_12, ssl0_14, ssk_12, ssk_14, \
                         ssl1_12, ssl1_14, spk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_x[k] * ssl0_12[k]
                  + f_5 * ssk_12[k]
                  - f_1 * pc_x[k] * ssl1_12[k];

        t_13[k] = f_2 * pc_y[k] * spk_9[k];

        t_14[k] = pb_x[k] * ssl0_14[k]
                  + f_5 * ssk_14[k]
                  - f_1 * pc_x[k] * ssl1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, pc_x, pc_z, ssl0_15, ssl0_17, ssk_15, ssk_17, \
                         ssl1_15, ssl1_17, spk_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_x[k] * ssl0_15[k]
                  + f_6 * ssk_15[k]
                  - f_1 * pc_x[k] * ssl1_15[k];

        t_16[k] = f_2 * pc_z[k] * spk_10[k];

        t_17[k] = pb_x[k] * ssl0_17[k]
                  + f_6 * ssk_17[k]
                  - f_1 * pc_x[k] * ssl1_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, pc_x, pc_y, ssl0_18, ssl0_20, ssk_18, ssk_20, \
                         ssl1_18, ssl1_20, spk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_x[k] * ssl0_18[k]
                  + f_6 * ssk_18[k]
                  - f_1 * pc_x[k] * ssl1_18[k];

        t_19[k] = f_2 * pc_y[k] * spk_14[k];

        t_20[k] = pb_x[k] * ssl0_20[k]
                  + f_6 * ssk_20[k]
                  - f_1 * pc_x[k] * ssl1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_x, pc_x, pc_z, ssl0_21, ssl0_23, ssk_21, ssk_23, \
                         ssl1_21, ssl1_23, spk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pb_x[k] * ssl0_21[k]
                  + f_7 * ssk_21[k]
                  - f_1 * pc_x[k] * ssl1_21[k];

        t_22[k] = f_2 * pc_z[k] * spk_15[k];

        t_23[k] = pb_x[k] * ssl0_23[k]
                  + f_7 * ssk_23[k]
                  - f_1 * pc_x[k] * ssl1_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, pc_x, pc_y, ssl0_24, ssl0_25, ssk_24, ssk_25, \
                         ssl1_24, ssl1_25, spk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_x[k] * ssl0_24[k]
                  + f_7 * ssk_24[k]
                  - f_1 * pc_x[k] * ssl1_24[k];

        t_25[k] = pb_x[k] * ssl0_25[k]
                  + f_7 * ssk_25[k]
                  - f_1 * pc_x[k] * ssl1_25[k];

        t_26[k] = f_2 * pc_y[k] * spk_20[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, pc_x, ssl0_27, ssk_27, ssk_28, ssk_29, \
                         ssk_30, ssl1_27, spk_28, spk_29, spk_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_x[k] * ssl0_27[k]
                  + f_7 * ssk_27[k]
                  - f_1 * pc_x[k] * ssl1_27[k];

        t_28[k] = f_8 * ssk_28[k]
                  + f_2 * pc_x[k] * spk_28[k];

        t_29[k] = f_8 * ssk_29[k]
                  + f_2 * pc_x[k] * spk_29[k];

        t_30[k] = f_8 * ssk_30[k]
                  + f_2 * pc_x[k] * spk_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, ssk_31, ssk_32, ssk_33, ssk_34, \
                         ssk_35, spk_31, spk_32, spk_33, spk_34, \
                         spk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_8 * ssk_31[k]
                  + f_2 * pc_x[k] * spk_31[k];

        t_32[k] = f_8 * ssk_32[k]
                  + f_2 * pc_x[k] * spk_32[k];

        t_33[k] = f_8 * ssk_33[k]
                  + f_2 * pc_x[k] * spk_33[k];

        t_34[k] = f_8 * ssk_34[k]
                  + f_2 * pc_x[k] * spk_34[k];

        t_35[k] = f_8 * ssk_35[k]
                  + f_2 * pc_x[k] * spk_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_x, pc_x, pc_z, ssl0_36, ssl0_38, ssl0_39, \
                         ssl1_36, ssl1_38, ssl1_39, spk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_x[k] * ssl0_36[k]
                  - f_1 * pc_x[k] * ssl1_36[k];

        t_37[k] = f_2 * pc_z[k] * spk_28[k];

        t_38[k] = pb_x[k] * ssl0_38[k]
                  - f_1 * pc_x[k] * ssl1_38[k];

        t_39[k] = pb_x[k] * ssl0_39[k]
                  - f_1 * pc_x[k] * ssl1_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_x, pc_x, pc_y, ssl0_40, ssl0_41, ssl0_42, \
                         ssl1_40, ssl1_41, ssl1_42, spk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_x[k] * ssl0_40[k]
                  - f_1 * pc_x[k] * ssl1_40[k];

        t_41[k] = pb_x[k] * ssl0_41[k]
                  - f_1 * pc_x[k] * ssl1_41[k];

        t_42[k] = pb_x[k] * ssl0_42[k]
                  - f_1 * pc_x[k] * ssl1_42[k];

        t_43[k] = f_2 * pc_y[k] * spk_35[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_y, pc_x, pc_y, pc_z, ssl0_0, \
                         ssl0_44, ssk_0, ssl1_0, ssl1_44, spk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * ssl0_44[k]
                  - f_1 * pc_x[k] * ssl1_44[k];

        t_45[k] = pb_y[k] * ssl0_0[k]
                  - f_1 * pc_y[k] * ssl1_0[k];

        t_46[k] = f_8 * ssk_0[k]
                  + f_2 * pc_y[k] * spk_36[k];

        t_47[k] = f_2 * pc_z[k] * spk_36[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_y, pc_x, pc_y, ssl0_5, ssk_2, ssl1_5, spi0_31, \
                         spi1_31, spk_38, spk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_9 * spi0_31[k]
                  - f_10 * spi1_31[k]
                  + f_2 * pc_x[k] * spk_39[k];

        t_49[k] = f_8 * ssk_2[k]
                  + f_2 * pc_y[k] * spk_38[k];

        t_50[k] = pb_y[k] * ssl0_5[k]
                  - f_1 * pc_y[k] * ssl1_5[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_y, pc_x, pc_y, pc_z, ssl0_9, ssk_5, \
                         ssl1_9, spi0_34, spi1_34, spk_39, spk_41, \
                         spk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_11 * spi0_34[k]
                  - f_12 * spi1_34[k]
                  + f_2 * pc_x[k] * spk_42[k];

        t_52[k] = f_2 * pc_z[k] * spk_39[k];

        t_53[k] = f_8 * ssk_5[k]
                  + f_2 * pc_y[k] * spk_41[k];

        t_54[k] = pb_y[k] * ssl0_9[k]
                  - f_1 * pc_y[k] * ssl1_9[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, ssk_9, spi0_38, spi0_40, \
                         spi1_38, spi1_40, spk_42, spk_45, spk_46, \
                         spk_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_13 * spi0_38[k]
                  - f_14 * spi1_38[k]
                  + f_2 * pc_x[k] * spk_46[k];

        t_56[k] = f_2 * pc_z[k] * spk_42[k];

        t_57[k] = f_13 * spi0_40[k]
                  - f_14 * spi1_40[k]
                  + f_2 * pc_x[k] * spk_48[k];

        t_58[k] = f_8 * ssk_9[k]
                  + f_2 * pc_y[k] * spk_45[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_y, pc_x, pc_y, pc_z, ssl0_14, ssl1_14, spi0_43, \
                         spi1_43, spk_46, spk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pb_y[k] * ssl0_14[k]
                  - f_1 * pc_y[k] * ssl1_14[k];

        t_60[k] = f_15 * spi0_43[k]
                  - f_16 * spi1_43[k]
                  + f_2 * pc_x[k] * spk_51[k];

        t_61[k] = f_2 * pc_z[k] * spk_46[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pc_x, pc_y, ssk_14, spi0_45, spi0_46, spi1_45, \
                         spi1_46, spk_50, spk_53, spk_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_15 * spi0_45[k]
                  - f_16 * spi1_45[k]
                  + f_2 * pc_x[k] * spk_53[k];

        t_63[k] = f_15 * spi0_46[k]
                  - f_16 * spi1_46[k]
                  + f_2 * pc_x[k] * spk_54[k];

        t_64[k] = f_8 * ssk_14[k]
                  + f_2 * pc_y[k] * spk_50[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pb_y, pc_x, pc_y, pc_z, ssl0_20, ssl1_20, spi0_49, \
                         spi1_49, spk_51, spk_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * ssl0_20[k]
                  - f_1 * pc_y[k] * ssl1_20[k];

        t_66[k] = f_17 * spi0_49[k]
                  - f_18 * spi1_49[k]
                  + f_2 * pc_x[k] * spk_57[k];

        t_67[k] = f_2 * pc_z[k] * spk_51[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pc_x, spi0_51, spi0_52, spi0_53, spi1_51, spi1_52, \
                         spi1_53, spk_59, spk_60, spk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_17 * spi0_51[k]
                  - f_18 * spi1_51[k]
                  + f_2 * pc_x[k] * spk_59[k];

        t_69[k] = f_17 * spi0_52[k]
                  - f_18 * spi1_52[k]
                  + f_2 * pc_x[k] * spk_60[k];

        t_70[k] = f_17 * spi0_53[k]
                  - f_18 * spi1_53[k]
                  + f_2 * pc_x[k] * spk_61[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, pb_y, pc_x, pc_y, ssl0_27, ssk_20, \
                         ssl1_27, spk_56, spk_64, spk_65, spk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_8 * ssk_20[k]
                  + f_2 * pc_y[k] * spk_56[k];

        t_72[k] = pb_y[k] * ssl0_27[k]
                  - f_1 * pc_y[k] * ssl1_27[k];

        t_73[k] = f_2 * pc_x[k] * spk_64[k];

        t_74[k] = f_2 * pc_x[k] * spk_65[k];

        t_75[k] = f_2 * pc_x[k] * spk_66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pc_x, spk_67, spk_68, spk_69, spk_70, \
                         spk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_2 * pc_x[k] * spk_67[k];

        t_77[k] = f_2 * pc_x[k] * spk_68[k];

        t_78[k] = f_2 * pc_x[k] * spk_69[k];

        t_79[k] = f_2 * pc_x[k] * spk_70[k];

        t_80[k] = f_2 * pc_x[k] * spk_71[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_y, pc_y, pc_z, ssl0_36, ssl0_38, ssk_28, ssk_30, \
                         ssl1_36, ssl1_38, spk_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pb_y[k] * ssl0_36[k]
                  + f_0 * ssk_28[k]
                  - f_1 * pc_y[k] * ssl1_36[k];

        t_82[k] = f_2 * pc_z[k] * spk_64[k];

        t_83[k] = pb_y[k] * ssl0_38[k]
                  + f_3 * ssk_30[k]
                  - f_1 * pc_y[k] * ssl1_38[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pb_y, pc_y, ssl0_39, ssl0_40, ssl0_41, ssk_31, \
                         ssk_32, ssk_33, ssl1_39, ssl1_40, ssl1_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pb_y[k] * ssl0_39[k]
                  + f_4 * ssk_31[k]
                  - f_1 * pc_y[k] * ssl1_39[k];

        t_85[k] = pb_y[k] * ssl0_40[k]
                  + f_5 * ssk_32[k]
                  - f_1 * pc_y[k] * ssl1_40[k];

        t_86[k] = pb_y[k] * ssl0_41[k]
                  + f_6 * ssk_33[k]
                  - f_1 * pc_y[k] * ssl1_41[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_y, pc_y, ssl0_42, ssl0_44, ssk_34, ssk_35, \
                         ssl1_42, ssl1_44, spk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pb_y[k] * ssl0_42[k]
                  + f_7 * ssk_34[k]
                  - f_1 * pc_y[k] * ssl1_42[k];

        t_88[k] = f_8 * ssk_35[k]
                  + f_2 * pc_y[k] * spk_71[k];

        t_89[k] = pb_y[k] * ssl0_44[k]
                  - f_1 * pc_y[k] * ssl1_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pb_z, pc_y, pc_z, ssl0_0, ssl0_3, \
                         ssk_0, ssl1_0, ssl1_3, spk_72, spk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_z[k] * ssl0_0[k]
                  - f_1 * pc_z[k] * ssl1_0[k];

        t_91[k] = f_2 * pc_y[k] * spk_72[k];

        t_92[k] = f_8 * ssk_0[k]
                  + f_2 * pc_z[k] * spk_72[k];

        t_93[k] = pb_z[k] * ssl0_3[k]
                  - f_1 * pc_z[k] * ssl1_3[k];

        t_94[k] = f_2 * pc_y[k] * spk_74[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pb_z, pc_x, pc_y, pc_z, ssl0_6, ssk_3, \
                         ssl1_6, spi0_61, spi1_61, spk_75, spk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_9 * spi0_61[k]
                  - f_10 * spi1_61[k]
                  + f_2 * pc_x[k] * spk_77[k];

        t_96[k] = pb_z[k] * ssl0_6[k]
                  - f_1 * pc_z[k] * ssl1_6[k];

        t_97[k] = f_8 * ssk_3[k]
                  + f_2 * pc_z[k] * spk_75[k];

        t_98[k] = f_2 * pc_y[k] * spk_77[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pb_z, pc_x, pc_z, ssl0_10, ssk_6, ssl1_10, \
                         spi0_65, spi1_65, spk_78, spk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_11 * spi0_65[k]
                  - f_12 * spi1_65[k]
                  + f_2 * pc_x[k] * spk_81[k];

        t_100[k] = pb_z[k] * ssl0_10[k]
                   - f_1 * pc_z[k] * ssl1_10[k];

        t_101[k] = f_8 * ssk_6[k]
                   + f_2 * pc_z[k] * spk_78[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pc_x, pc_y, spi0_68, spi0_70, spi1_68, spi1_70, \
                         spk_81, spk_84, spk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_13 * spi0_68[k]
                   - f_14 * spi1_68[k]
                   + f_2 * pc_x[k] * spk_84[k];

        t_103[k] = f_2 * pc_y[k] * spk_81[k];

        t_104[k] = f_13 * spi0_70[k]
                   - f_14 * spi1_70[k]
                   + f_2 * pc_x[k] * spk_86[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pb_z, pc_x, pc_z, ssl0_15, ssk_10, ssl1_15, \
                         spi0_73, spi1_73, spk_82, spk_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pb_z[k] * ssl0_15[k]
                   - f_1 * pc_z[k] * ssl1_15[k];

        t_106[k] = f_8 * ssk_10[k]
                   + f_2 * pc_z[k] * spk_82[k];

        t_107[k] = f_15 * spi0_73[k]
                   - f_16 * spi1_73[k]
                   + f_2 * pc_x[k] * spk_89[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pc_x, pc_y, spi0_74, spi0_76, spi1_74, spi1_76, \
                         spk_86, spk_90, spk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_15 * spi0_74[k]
                   - f_16 * spi1_74[k]
                   + f_2 * pc_x[k] * spk_90[k];

        t_109[k] = f_2 * pc_y[k] * spk_86[k];

        t_110[k] = f_15 * spi0_76[k]
                   - f_16 * spi1_76[k]
                   + f_2 * pc_x[k] * spk_92[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pb_z, pc_x, pc_z, ssl0_21, ssk_15, ssl1_21, \
                         spi0_79, spi1_79, spk_87, spk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pb_z[k] * ssl0_21[k]
                   - f_1 * pc_z[k] * ssl1_21[k];

        t_112[k] = f_8 * ssk_15[k]
                   + f_2 * pc_z[k] * spk_87[k];

        t_113[k] = f_17 * spi0_79[k]
                   - f_18 * spi1_79[k]
                   + f_2 * pc_x[k] * spk_95[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, pc_y, spi0_80, spi0_81, spi0_83, \
                         spi1_80, spi1_81, spi1_83, spk_92, spk_96, spk_97, \
                         spk_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_17 * spi0_80[k]
                   - f_18 * spi1_80[k]
                   + f_2 * pc_x[k] * spk_96[k];

        t_115[k] = f_17 * spi0_81[k]
                   - f_18 * spi1_81[k]
                   + f_2 * pc_x[k] * spk_97[k];

        t_116[k] = f_2 * pc_y[k] * spk_92[k];

        t_117[k] = f_17 * spi0_83[k]
                   - f_18 * spi1_83[k]
                   + f_2 * pc_x[k] * spk_99[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, t_123, t_124, pc_x, spk_100, \
                         spk_101, spk_102, spk_103, spk_104, spk_105, \
                         spk_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_2 * pc_x[k] * spk_100[k];

        t_119[k] = f_2 * pc_x[k] * spk_101[k];

        t_120[k] = f_2 * pc_x[k] * spk_102[k];

        t_121[k] = f_2 * pc_x[k] * spk_103[k];

        t_122[k] = f_2 * pc_x[k] * spk_104[k];

        t_123[k] = f_2 * pc_x[k] * spk_105[k];

        t_124[k] = f_2 * pc_x[k] * spk_106[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_z, pc_x, pc_y, pc_z, ssl0_36, ssk_28, \
                         ssl1_36, spi0_79, spi1_79, spk_100, spk_102, \
                         spk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_2 * pc_x[k] * spk_107[k];

        t_126[k] = pb_z[k] * ssl0_36[k]
                   - f_1 * pc_z[k] * ssl1_36[k];

        t_127[k] = f_8 * ssk_28[k]
                   + f_2 * pc_z[k] * spk_100[k];

        t_128[k] = f_9 * spi0_79[k]
                   - f_10 * spi1_79[k]
                   + f_2 * pc_y[k] * spk_102[k];
    }
}

static auto
compute_prim_spl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t ssl0,
                                                          const size_t ssk, const size_t ssl1,
                                                          const size_t spi0, const size_t spi1,
                                                          const size_t spk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = gamma / q;
    const auto f_2 = p / q;
    const auto f_11 = 2.0 / gamma;
    const auto f_12 = 2.0 * p / (gamma * q);
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 1.0 / gamma;
    const auto f_16 = p / (gamma * q);
    const auto f_17 = 0.5 / gamma;
    const auto f_18 = 0.5 * p / (gamma * q);

    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssl0_44 = buffer.data(ssl0 + 44);

    const auto *ssk_35 = buffer.data(ssk + 35);

    const auto *ssl1_44 = buffer.data(ssl1 + 44);

    const auto *spi0_80 = buffer.data(spi0 + 80);
    const auto *spi0_81 = buffer.data(spi0 + 81);
    const auto *spi0_82 = buffer.data(spi0 + 82);
    const auto *spi0_83 = buffer.data(spi0 + 83);

    const auto *spi1_80 = buffer.data(spi1 + 80);
    const auto *spi1_81 = buffer.data(spi1 + 81);
    const auto *spi1_82 = buffer.data(spi1 + 82);
    const auto *spi1_83 = buffer.data(spi1 + 83);

    const auto *spk_103 = buffer.data(spk + 103);
    const auto *spk_104 = buffer.data(spk + 104);
    const auto *spk_105 = buffer.data(spk + 105);
    const auto *spk_106 = buffer.data(spk + 106);
    const auto *spk_107 = buffer.data(spk + 107);

#pragma omp simd aligned(t_129, t_130, t_131, pc_y, spi0_80, spi0_81, spi0_82, spi1_80, \
                         spi1_81, spi1_82, spk_103, spk_104, spk_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_11 * spi0_80[k]
                   - f_12 * spi1_80[k]
                   + f_2 * pc_y[k] * spk_103[k];

        t_130[k] = f_13 * spi0_81[k]
                   - f_14 * spi1_81[k]
                   + f_2 * pc_y[k] * spk_104[k];

        t_131[k] = f_15 * spi0_82[k]
                   - f_16 * spi1_82[k]
                   + f_2 * pc_y[k] * spk_105[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pb_z, pc_y, pc_z, ssl0_44, ssk_35, ssl1_44, \
                         spi0_83, spi1_83, spk_106, spk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_17 * spi0_83[k]
                   - f_18 * spi1_83[k]
                   + f_2 * pc_y[k] * spk_106[k];

        t_133[k] = f_2 * pc_y[k] * spk_107[k];

        t_134[k] = pb_z[k] * ssl0_44[k]
                   + f_0 * ssk_35[k]
                   - f_1 * pc_z[k] * ssl1_44[k];
    }
}

auto
compute_prim_spl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t ssl0, const size_t ssk,
                                                   const size_t ssl1, const size_t spi0,
                                                   const size_t spi1, const size_t spk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_spl_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, ssl0, ssk,
                                                              ssl1, spi0, spi1, spk, ncols,
                                                              gamma, p, q);

    compute_prim_spl_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, ssl0, ssk,
                                                              ssl1, spi0, spi1, spk, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
