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


#include "SimdThreeCenterElectronRepulsionVrrRecSOH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snh0,
                                                          const size_t sng, const size_t snh1,
                                                          const size_t sof0, const size_t sof1,
                                                          const size_t sog, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 5.0 / q;
    const auto f_13 = 4.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_0 = buffer.data(snh0 + 0);
    const auto *snh0_3 = buffer.data(snh0 + 3);
    const auto *snh0_5 = buffer.data(snh0 + 5);
    const auto *snh0_6 = buffer.data(snh0 + 6);
    const auto *snh0_9 = buffer.data(snh0 + 9);
    const auto *snh0_15 = buffer.data(snh0 + 15);
    const auto *snh0_20 = buffer.data(snh0 + 20);
    const auto *snh0_24 = buffer.data(snh0 + 24);
    const auto *snh0_27 = buffer.data(snh0 + 27);
    const auto *snh0_36 = buffer.data(snh0 + 36);
    const auto *snh0_42 = buffer.data(snh0 + 42);
    const auto *snh0_47 = buffer.data(snh0 + 47);
    const auto *snh0_51 = buffer.data(snh0 + 51);
    const auto *snh0_62 = buffer.data(snh0 + 62);

    const auto *sng_0 = buffer.data(sng + 0);
    const auto *sng_1 = buffer.data(sng + 1);
    const auto *sng_2 = buffer.data(sng + 2);
    const auto *sng_3 = buffer.data(sng + 3);
    const auto *sng_5 = buffer.data(sng + 5);
    const auto *sng_6 = buffer.data(sng + 6);
    const auto *sng_9 = buffer.data(sng + 9);
    const auto *sng_10 = buffer.data(sng + 10);
    const auto *sng_11 = buffer.data(sng + 11);
    const auto *sng_12 = buffer.data(sng + 12);
    const auto *sng_13 = buffer.data(sng + 13);
    const auto *sng_14 = buffer.data(sng + 14);
    const auto *sng_15 = buffer.data(sng + 15);
    const auto *sng_17 = buffer.data(sng + 17);
    const auto *sng_18 = buffer.data(sng + 18);
    const auto *sng_20 = buffer.data(sng + 20);
    const auto *sng_25 = buffer.data(sng + 25);
    const auto *sng_26 = buffer.data(sng + 26);
    const auto *sng_27 = buffer.data(sng + 27);
    const auto *sng_28 = buffer.data(sng + 28);
    const auto *sng_29 = buffer.data(sng + 29);
    const auto *sng_30 = buffer.data(sng + 30);
    const auto *sng_32 = buffer.data(sng + 32);
    const auto *sng_33 = buffer.data(sng + 33);
    const auto *sng_35 = buffer.data(sng + 35);
    const auto *sng_40 = buffer.data(sng + 40);
    const auto *sng_41 = buffer.data(sng + 41);
    const auto *sng_42 = buffer.data(sng + 42);
    const auto *sng_43 = buffer.data(sng + 43);
    const auto *sng_44 = buffer.data(sng + 44);
    const auto *sng_45 = buffer.data(sng + 45);
    const auto *sng_48 = buffer.data(sng + 48);
    const auto *sng_50 = buffer.data(sng + 50);
    const auto *sng_51 = buffer.data(sng + 51);
    const auto *sng_54 = buffer.data(sng + 54);
    const auto *sng_55 = buffer.data(sng + 55);
    const auto *sng_56 = buffer.data(sng + 56);
    const auto *sng_57 = buffer.data(sng + 57);
    const auto *sng_58 = buffer.data(sng + 58);
    const auto *sng_59 = buffer.data(sng + 59);
    const auto *sng_70 = buffer.data(sng + 70);
    const auto *sng_71 = buffer.data(sng + 71);
    const auto *sng_72 = buffer.data(sng + 72);
    const auto *sng_73 = buffer.data(sng + 73);
    const auto *sng_74 = buffer.data(sng + 74);
    const auto *sng_75 = buffer.data(sng + 75);
    const auto *sng_78 = buffer.data(sng + 78);
    const auto *sng_80 = buffer.data(sng + 80);
    const auto *sng_81 = buffer.data(sng + 81);
    const auto *sng_84 = buffer.data(sng + 84);
    const auto *sng_85 = buffer.data(sng + 85);
    const auto *sng_86 = buffer.data(sng + 86);
    const auto *sng_87 = buffer.data(sng + 87);
    const auto *sng_88 = buffer.data(sng + 88);
    const auto *sng_89 = buffer.data(sng + 89);

    const auto *snh1_0 = buffer.data(snh1 + 0);
    const auto *snh1_3 = buffer.data(snh1 + 3);
    const auto *snh1_5 = buffer.data(snh1 + 5);
    const auto *snh1_6 = buffer.data(snh1 + 6);
    const auto *snh1_9 = buffer.data(snh1 + 9);
    const auto *snh1_15 = buffer.data(snh1 + 15);
    const auto *snh1_20 = buffer.data(snh1 + 20);
    const auto *snh1_24 = buffer.data(snh1 + 24);
    const auto *snh1_27 = buffer.data(snh1 + 27);
    const auto *snh1_36 = buffer.data(snh1 + 36);
    const auto *snh1_42 = buffer.data(snh1 + 42);
    const auto *snh1_47 = buffer.data(snh1 + 47);
    const auto *snh1_51 = buffer.data(snh1 + 51);
    const auto *snh1_62 = buffer.data(snh1 + 62);

    const auto *sof0_0 = buffer.data(sof0 + 0);
    const auto *sof0_3 = buffer.data(sof0 + 3);
    const auto *sof0_5 = buffer.data(sof0 + 5);
    const auto *sof0_6 = buffer.data(sof0 + 6);
    const auto *sof0_8 = buffer.data(sof0 + 8);
    const auto *sof0_9 = buffer.data(sof0 + 9);
    const auto *sof0_16 = buffer.data(sof0 + 16);
    const auto *sof0_18 = buffer.data(sof0 + 18);
    const auto *sof0_19 = buffer.data(sof0 + 19);
    const auto *sof0_28 = buffer.data(sof0 + 28);
    const auto *sof0_29 = buffer.data(sof0 + 29);
    const auto *sof0_30 = buffer.data(sof0 + 30);
    const auto *sof0_33 = buffer.data(sof0 + 33);
    const auto *sof0_35 = buffer.data(sof0 + 35);
    const auto *sof0_36 = buffer.data(sof0 + 36);
    const auto *sof0_38 = buffer.data(sof0 + 38);
    const auto *sof0_39 = buffer.data(sof0 + 39);
    const auto *sof0_48 = buffer.data(sof0 + 48);
    const auto *sof0_49 = buffer.data(sof0 + 49);
    const auto *sof0_50 = buffer.data(sof0 + 50);
    const auto *sof0_53 = buffer.data(sof0 + 53);
    const auto *sof0_55 = buffer.data(sof0 + 55);
    const auto *sof0_56 = buffer.data(sof0 + 56);
    const auto *sof0_58 = buffer.data(sof0 + 58);
    const auto *sof0_59 = buffer.data(sof0 + 59);

    const auto *sof1_0 = buffer.data(sof1 + 0);
    const auto *sof1_3 = buffer.data(sof1 + 3);
    const auto *sof1_5 = buffer.data(sof1 + 5);
    const auto *sof1_6 = buffer.data(sof1 + 6);
    const auto *sof1_8 = buffer.data(sof1 + 8);
    const auto *sof1_9 = buffer.data(sof1 + 9);
    const auto *sof1_16 = buffer.data(sof1 + 16);
    const auto *sof1_18 = buffer.data(sof1 + 18);
    const auto *sof1_19 = buffer.data(sof1 + 19);
    const auto *sof1_28 = buffer.data(sof1 + 28);
    const auto *sof1_29 = buffer.data(sof1 + 29);
    const auto *sof1_30 = buffer.data(sof1 + 30);
    const auto *sof1_33 = buffer.data(sof1 + 33);
    const auto *sof1_35 = buffer.data(sof1 + 35);
    const auto *sof1_36 = buffer.data(sof1 + 36);
    const auto *sof1_38 = buffer.data(sof1 + 38);
    const auto *sof1_39 = buffer.data(sof1 + 39);
    const auto *sof1_48 = buffer.data(sof1 + 48);
    const auto *sof1_49 = buffer.data(sof1 + 49);
    const auto *sof1_50 = buffer.data(sof1 + 50);
    const auto *sof1_53 = buffer.data(sof1 + 53);
    const auto *sof1_55 = buffer.data(sof1 + 55);
    const auto *sof1_56 = buffer.data(sof1 + 56);
    const auto *sof1_58 = buffer.data(sof1 + 58);
    const auto *sof1_59 = buffer.data(sof1 + 59);

    const auto *sog_0 = buffer.data(sog + 0);
    const auto *sog_2 = buffer.data(sog + 2);
    const auto *sog_3 = buffer.data(sog + 3);
    const auto *sog_5 = buffer.data(sog + 5);
    const auto *sog_6 = buffer.data(sog + 6);
    const auto *sog_9 = buffer.data(sog + 9);
    const auto *sog_10 = buffer.data(sog + 10);
    const auto *sog_11 = buffer.data(sog + 11);
    const auto *sog_12 = buffer.data(sog + 12);
    const auto *sog_13 = buffer.data(sog + 13);
    const auto *sog_14 = buffer.data(sog + 14);
    const auto *sog_15 = buffer.data(sog + 15);
    const auto *sog_17 = buffer.data(sog + 17);
    const auto *sog_18 = buffer.data(sog + 18);
    const auto *sog_20 = buffer.data(sog + 20);
    const auto *sog_25 = buffer.data(sog + 25);
    const auto *sog_26 = buffer.data(sog + 26);
    const auto *sog_27 = buffer.data(sog + 27);
    const auto *sog_28 = buffer.data(sog + 28);
    const auto *sog_29 = buffer.data(sog + 29);
    const auto *sog_30 = buffer.data(sog + 30);
    const auto *sog_32 = buffer.data(sog + 32);
    const auto *sog_33 = buffer.data(sog + 33);
    const auto *sog_35 = buffer.data(sog + 35);
    const auto *sog_40 = buffer.data(sog + 40);
    const auto *sog_41 = buffer.data(sog + 41);
    const auto *sog_42 = buffer.data(sog + 42);
    const auto *sog_43 = buffer.data(sog + 43);
    const auto *sog_44 = buffer.data(sog + 44);
    const auto *sog_45 = buffer.data(sog + 45);
    const auto *sog_47 = buffer.data(sog + 47);
    const auto *sog_48 = buffer.data(sog + 48);
    const auto *sog_50 = buffer.data(sog + 50);
    const auto *sog_51 = buffer.data(sog + 51);
    const auto *sog_54 = buffer.data(sog + 54);
    const auto *sog_55 = buffer.data(sog + 55);
    const auto *sog_56 = buffer.data(sog + 56);
    const auto *sog_57 = buffer.data(sog + 57);
    const auto *sog_58 = buffer.data(sog + 58);
    const auto *sog_59 = buffer.data(sog + 59);
    const auto *sog_60 = buffer.data(sog + 60);
    const auto *sog_62 = buffer.data(sog + 62);
    const auto *sog_63 = buffer.data(sog + 63);
    const auto *sog_65 = buffer.data(sog + 65);
    const auto *sog_70 = buffer.data(sog + 70);
    const auto *sog_71 = buffer.data(sog + 71);
    const auto *sog_72 = buffer.data(sog + 72);
    const auto *sog_73 = buffer.data(sog + 73);
    const auto *sog_74 = buffer.data(sog + 74);
    const auto *sog_75 = buffer.data(sog + 75);
    const auto *sog_77 = buffer.data(sog + 77);
    const auto *sog_78 = buffer.data(sog + 78);
    const auto *sog_80 = buffer.data(sog + 80);
    const auto *sog_81 = buffer.data(sog + 81);
    const auto *sog_84 = buffer.data(sog + 84);
    const auto *sog_85 = buffer.data(sog + 85);
    const auto *sog_86 = buffer.data(sog + 86);
    const auto *sog_87 = buffer.data(sog + 87);
    const auto *sog_88 = buffer.data(sog + 88);
    const auto *sog_89 = buffer.data(sog + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sng_0, sng_3, sof0_0, sof0_3, \
                         sof1_0, sof1_3, sog_0, sog_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sng_0[k]
                 + f_1 * sof0_0[k]
                 - f_2 * sof1_0[k]
                 + f_3 * pc_x[k] * sog_0[k];

        t_1[k] = f_3 * pc_y[k] * sog_0[k];

        t_2[k] = f_3 * pc_z[k] * sog_0[k];

        t_3[k] = f_0 * sng_3[k]
                 + f_4 * sof0_3[k]
                 - f_5 * sof1_3[k]
                 + f_3 * pc_x[k] * sog_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sng_5, sng_6, sof0_5, sof0_6, sof1_5, \
                         sof1_6, sog_2, sog_5, sog_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sog_2[k];

        t_5[k] = f_0 * sng_5[k]
                 + f_4 * sof0_5[k]
                 - f_5 * sof1_5[k]
                 + f_3 * pc_x[k] * sog_5[k];

        t_6[k] = f_0 * sng_6[k]
                 + f_6 * sof0_6[k]
                 - f_7 * sof1_6[k]
                 + f_3 * pc_x[k] * sog_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, sng_9, sng_10, sof0_9, sof1_9, \
                         sog_3, sog_5, sog_9, sog_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sog_3[k];

        t_8[k] = f_3 * pc_y[k] * sog_5[k];

        t_9[k] = f_0 * sng_9[k]
                 + f_6 * sof0_9[k]
                 - f_7 * sof1_9[k]
                 + f_3 * pc_x[k] * sog_9[k];

        t_10[k] = f_0 * sng_10[k]
                  + f_3 * pc_x[k] * sog_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, sng_11, sng_12, sng_13, sng_14, sog_11, \
                         sog_12, sog_13, sog_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * sng_11[k]
                  + f_3 * pc_x[k] * sog_11[k];

        t_12[k] = f_0 * sng_12[k]
                  + f_3 * pc_x[k] * sog_12[k];

        t_13[k] = f_0 * sng_13[k]
                  + f_3 * pc_x[k] * sog_13[k];

        t_14[k] = f_0 * sng_14[k]
                  + f_3 * pc_x[k] * sog_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, sof0_6, sof0_8, sof0_9, sof1_6, \
                         sof1_8, sof1_9, sog_10, sog_12, sog_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * sof0_6[k]
                  - f_2 * sof1_6[k]
                  + f_3 * pc_y[k] * sog_10[k];

        t_16[k] = f_3 * pc_z[k] * sog_10[k];

        t_17[k] = f_4 * sof0_8[k]
                  - f_5 * sof1_8[k]
                  + f_3 * pc_y[k] * sog_12[k];

        t_18[k] = f_6 * sof0_9[k]
                  - f_7 * sof1_9[k]
                  + f_3 * pc_y[k] * sog_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_y, pc_y, pc_z, snh0_0, sng_0, \
                         snh1_0, sof0_9, sof1_9, sog_14, sog_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sog_14[k];

        t_20[k] = f_1 * sof0_9[k]
                  - f_2 * sof1_9[k]
                  + f_3 * pc_z[k] * sog_14[k];

        t_21[k] = pb_y[k] * snh0_0[k]
                  - f_8 * pc_y[k] * snh1_0[k];

        t_22[k] = f_9 * sng_0[k]
                  + f_3 * pc_y[k] * sog_15[k];

        t_23[k] = f_3 * pc_z[k] * sog_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_y, pc_y, snh0_3, snh0_5, snh0_6, sng_1, \
                         sng_2, sng_3, snh1_3, snh1_5, snh1_6, sog_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * snh0_3[k]
                  + f_10 * sng_1[k]
                  - f_8 * pc_y[k] * snh1_3[k];

        t_25[k] = f_9 * sng_2[k]
                  + f_3 * pc_y[k] * sog_17[k];

        t_26[k] = pb_y[k] * snh0_5[k]
                  - f_8 * pc_y[k] * snh1_5[k];

        t_27[k] = pb_y[k] * snh0_6[k]
                  + f_11 * sng_3[k]
                  - f_8 * pc_y[k] * snh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_x, pc_y, pc_z, snh0_9, sng_5, \
                         sng_25, snh1_9, sog_18, sog_20, sog_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * sog_18[k];

        t_29[k] = f_9 * sng_5[k]
                  + f_3 * pc_y[k] * sog_20[k];

        t_30[k] = pb_y[k] * snh0_9[k]
                  - f_8 * pc_y[k] * snh1_9[k];

        t_31[k] = f_12 * sng_25[k]
                  + f_3 * pc_x[k] * sog_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, sng_26, sng_27, sng_28, sng_29, sog_26, \
                         sog_27, sog_28, sog_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_12 * sng_26[k]
                  + f_3 * pc_x[k] * sog_26[k];

        t_33[k] = f_12 * sng_27[k]
                  + f_3 * pc_x[k] * sog_27[k];

        t_34[k] = f_12 * sng_28[k]
                  + f_3 * pc_x[k] * sog_28[k];

        t_35[k] = f_12 * sng_29[k]
                  + f_3 * pc_x[k] * sog_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pc_y, pc_z, sng_10, sng_12, sof0_16, sof0_18, \
                         sof1_16, sof1_18, sog_25, sog_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * sng_10[k]
                  + f_1 * sof0_16[k]
                  - f_2 * sof1_16[k]
                  + f_3 * pc_y[k] * sog_25[k];

        t_37[k] = f_3 * pc_z[k] * sog_25[k];

        t_38[k] = f_9 * sng_12[k]
                  + f_4 * sof0_18[k]
                  - f_5 * sof1_18[k]
                  + f_3 * pc_y[k] * sog_27[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, pc_y, snh0_20, sng_13, sng_14, snh1_20, \
                         sof0_19, sof1_19, sog_28, sog_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * sng_13[k]
                  + f_6 * sof0_19[k]
                  - f_7 * sof1_19[k]
                  + f_3 * pc_y[k] * sog_28[k];

        t_40[k] = f_9 * sng_14[k]
                  + f_3 * pc_y[k] * sog_29[k];

        t_41[k] = pb_y[k] * snh0_20[k]
                  - f_8 * pc_y[k] * snh1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, snh0_0, snh0_3, \
                         sng_0, snh1_0, snh1_3, sog_30, sog_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * snh0_0[k]
                  - f_8 * pc_z[k] * snh1_0[k];

        t_43[k] = f_3 * pc_y[k] * sog_30[k];

        t_44[k] = f_9 * sng_0[k]
                  + f_3 * pc_z[k] * sog_30[k];

        t_45[k] = pb_z[k] * snh0_3[k]
                  - f_8 * pc_z[k] * snh1_3[k];

        t_46[k] = f_3 * pc_y[k] * sog_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_y, pc_z, snh0_5, snh0_6, sng_2, \
                         sng_3, snh1_5, snh1_6, sog_33, sog_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * snh0_5[k]
                  + f_10 * sng_2[k]
                  - f_8 * pc_z[k] * snh1_5[k];

        t_48[k] = pb_z[k] * snh0_6[k]
                  - f_8 * pc_z[k] * snh1_6[k];

        t_49[k] = f_9 * sng_3[k]
                  + f_3 * pc_z[k] * sog_33[k];

        t_50[k] = f_3 * pc_y[k] * sog_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_z, pc_x, pc_z, snh0_9, sng_5, sng_40, \
                         sng_41, sng_42, snh1_9, sog_40, sog_41, \
                         sog_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_z[k] * snh0_9[k]
                  + f_11 * sng_5[k]
                  - f_8 * pc_z[k] * snh1_9[k];

        t_52[k] = f_12 * sng_40[k]
                  + f_3 * pc_x[k] * sog_40[k];

        t_53[k] = f_12 * sng_41[k]
                  + f_3 * pc_x[k] * sog_41[k];

        t_54[k] = f_12 * sng_42[k]
                  + f_3 * pc_x[k] * sog_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_z, pc_x, pc_z, snh0_15, sng_10, sng_43, \
                         sng_44, snh1_15, sog_40, sog_43, sog_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_12 * sng_43[k]
                  + f_3 * pc_x[k] * sog_43[k];

        t_56[k] = f_12 * sng_44[k]
                  + f_3 * pc_x[k] * sog_44[k];

        t_57[k] = pb_z[k] * snh0_15[k]
                  - f_8 * pc_z[k] * snh1_15[k];

        t_58[k] = f_9 * sng_10[k]
                  + f_3 * pc_z[k] * sog_40[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_y, pc_z, sng_14, sof0_28, sof0_29, \
                         sof1_28, sof1_29, sog_42, sog_43, sog_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_4 * sof0_28[k]
                  - f_5 * sof1_28[k]
                  + f_3 * pc_y[k] * sog_42[k];

        t_60[k] = f_6 * sof0_29[k]
                  - f_7 * sof1_29[k]
                  + f_3 * pc_y[k] * sog_43[k];

        t_61[k] = f_3 * pc_y[k] * sog_44[k];

        t_62[k] = f_9 * sng_14[k]
                  + f_1 * sof0_29[k]
                  - f_2 * sof1_29[k]
                  + f_3 * pc_z[k] * sog_44[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, pc_z, sng_15, sng_45, sng_48, \
                         sof0_30, sof0_33, sof1_30, sof1_33, sog_45, \
                         sog_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_13 * sng_45[k]
                  + f_1 * sof0_30[k]
                  - f_2 * sof1_30[k]
                  + f_3 * pc_x[k] * sog_45[k];

        t_64[k] = f_10 * sng_15[k]
                  + f_3 * pc_y[k] * sog_45[k];

        t_65[k] = f_3 * pc_z[k] * sog_45[k];

        t_66[k] = f_13 * sng_48[k]
                  + f_4 * sof0_33[k]
                  - f_5 * sof1_33[k]
                  + f_3 * pc_x[k] * sog_48[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_x, pc_y, sng_17, sng_50, sng_51, sof0_35, \
                         sof0_36, sof1_35, sof1_36, sog_47, sog_50, \
                         sog_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_10 * sng_17[k]
                  + f_3 * pc_y[k] * sog_47[k];

        t_68[k] = f_13 * sng_50[k]
                  + f_4 * sof0_35[k]
                  - f_5 * sof1_35[k]
                  + f_3 * pc_x[k] * sog_50[k];

        t_69[k] = f_13 * sng_51[k]
                  + f_6 * sof0_36[k]
                  - f_7 * sof1_36[k]
                  + f_3 * pc_x[k] * sog_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pc_x, pc_y, pc_z, sng_20, sng_54, sng_55, \
                         sof0_39, sof1_39, sog_48, sog_50, sog_54, \
                         sog_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * pc_z[k] * sog_48[k];

        t_71[k] = f_10 * sng_20[k]
                  + f_3 * pc_y[k] * sog_50[k];

        t_72[k] = f_13 * sng_54[k]
                  + f_6 * sof0_39[k]
                  - f_7 * sof1_39[k]
                  + f_3 * pc_x[k] * sog_54[k];

        t_73[k] = f_13 * sng_55[k]
                  + f_3 * pc_x[k] * sog_55[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, sng_56, sng_57, sng_58, sng_59, sog_56, \
                         sog_57, sog_58, sog_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * sng_56[k]
                  + f_3 * pc_x[k] * sog_56[k];

        t_75[k] = f_13 * sng_57[k]
                  + f_3 * pc_x[k] * sog_57[k];

        t_76[k] = f_13 * sng_58[k]
                  + f_3 * pc_x[k] * sog_58[k];

        t_77[k] = f_13 * sng_59[k]
                  + f_3 * pc_x[k] * sog_59[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, sng_25, sng_27, sof0_36, sof0_38, \
                         sof1_36, sof1_38, sog_55, sog_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * sng_25[k]
                  + f_1 * sof0_36[k]
                  - f_2 * sof1_36[k]
                  + f_3 * pc_y[k] * sog_55[k];

        t_79[k] = f_3 * pc_z[k] * sog_55[k];

        t_80[k] = f_10 * sng_27[k]
                  + f_4 * sof0_38[k]
                  - f_5 * sof1_38[k]
                  + f_3 * pc_y[k] * sog_57[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pb_y, pc_y, pc_z, snh0_42, sng_28, sng_29, \
                         snh1_42, sof0_39, sof1_39, sog_58, sog_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_10 * sng_28[k]
                  + f_6 * sof0_39[k]
                  - f_7 * sof1_39[k]
                  + f_3 * pc_y[k] * sog_58[k];

        t_82[k] = f_10 * sng_29[k]
                  + f_3 * pc_y[k] * sog_59[k];

        t_83[k] = f_1 * sof0_39[k]
                  - f_2 * sof1_39[k]
                  + f_3 * pc_z[k] * sog_59[k];

        t_84[k] = pb_y[k] * snh0_42[k]
                  - f_8 * pc_y[k] * snh1_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_z, pc_y, pc_z, snh0_24, sng_15, sng_30, \
                         sng_32, snh1_24, sog_60, sog_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_9 * sng_30[k]
                  + f_3 * pc_y[k] * sog_60[k];

        t_86[k] = f_9 * sng_15[k]
                  + f_3 * pc_z[k] * sog_60[k];

        t_87[k] = pb_z[k] * snh0_24[k]
                  - f_8 * pc_z[k] * snh1_24[k];

        t_88[k] = f_9 * sng_32[k]
                  + f_3 * pc_y[k] * sog_62[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pb_y, pb_z, pc_y, pc_z, snh0_27, snh0_47, \
                         sng_18, sng_35, snh1_27, snh1_47, sog_63, \
                         sog_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pb_y[k] * snh0_47[k]
                  - f_8 * pc_y[k] * snh1_47[k];

        t_90[k] = pb_z[k] * snh0_27[k]
                  - f_8 * pc_z[k] * snh1_27[k];

        t_91[k] = f_9 * sng_18[k]
                  + f_3 * pc_z[k] * sog_63[k];

        t_92[k] = f_9 * sng_35[k]
                  + f_3 * pc_y[k] * sog_65[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_y, pc_x, pc_y, snh0_51, sng_70, sng_71, \
                         sng_72, snh1_51, sog_70, sog_71, sog_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pb_y[k] * snh0_51[k]
                  - f_8 * pc_y[k] * snh1_51[k];

        t_94[k] = f_13 * sng_70[k]
                  + f_3 * pc_x[k] * sog_70[k];

        t_95[k] = f_13 * sng_71[k]
                  + f_3 * pc_x[k] * sog_71[k];

        t_96[k] = f_13 * sng_72[k]
                  + f_3 * pc_x[k] * sog_72[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_z, pc_x, pc_z, snh0_36, sng_25, sng_73, \
                         sng_74, snh1_36, sog_70, sog_73, sog_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_13 * sng_73[k]
                  + f_3 * pc_x[k] * sog_73[k];

        t_98[k] = f_13 * sng_74[k]
                  + f_3 * pc_x[k] * sog_74[k];

        t_99[k] = pb_z[k] * snh0_36[k]
                  - f_8 * pc_z[k] * snh1_36[k];

        t_100[k] = f_9 * sng_25[k]
                   + f_3 * pc_z[k] * sog_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_y, sng_42, sng_43, sng_44, sof0_48, sof0_49, \
                         sof1_48, sof1_49, sog_72, sog_73, sog_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_9 * sng_42[k]
                   + f_4 * sof0_48[k]
                   - f_5 * sof1_48[k]
                   + f_3 * pc_y[k] * sog_72[k];

        t_102[k] = f_9 * sng_43[k]
                   + f_6 * sof0_49[k]
                   - f_7 * sof1_49[k]
                   + f_3 * pc_y[k] * sog_73[k];

        t_103[k] = f_9 * sng_44[k]
                   + f_3 * pc_y[k] * sog_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_y, pc_x, pc_y, pc_z, snh0_62, sng_30, \
                         sng_75, snh1_62, sof0_50, sof1_50, sog_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_y[k] * snh0_62[k]
                   - f_8 * pc_y[k] * snh1_62[k];

        t_105[k] = f_13 * sng_75[k]
                   + f_1 * sof0_50[k]
                   - f_2 * sof1_50[k]
                   + f_3 * pc_x[k] * sog_75[k];

        t_106[k] = f_3 * pc_y[k] * sog_75[k];

        t_107[k] = f_10 * sng_30[k]
                   + f_3 * pc_z[k] * sog_75[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pc_x, pc_y, sng_78, sng_80, sof0_53, sof0_55, \
                         sof1_53, sof1_55, sog_77, sog_78, sog_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_13 * sng_78[k]
                   + f_4 * sof0_53[k]
                   - f_5 * sof1_53[k]
                   + f_3 * pc_x[k] * sog_78[k];

        t_109[k] = f_3 * pc_y[k] * sog_77[k];

        t_110[k] = f_13 * sng_80[k]
                   + f_4 * sof0_55[k]
                   - f_5 * sof1_55[k]
                   + f_3 * pc_x[k] * sog_80[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_y, pc_z, sng_33, sng_81, sof0_56, \
                         sof1_56, sog_78, sog_80, sog_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * sng_81[k]
                   + f_6 * sof0_56[k]
                   - f_7 * sof1_56[k]
                   + f_3 * pc_x[k] * sog_81[k];

        t_112[k] = f_10 * sng_33[k]
                   + f_3 * pc_z[k] * sog_78[k];

        t_113[k] = f_3 * pc_y[k] * sog_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, sng_84, sng_85, sng_86, sng_87, \
                         sof0_59, sof1_59, sog_84, sog_85, sog_86, \
                         sog_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_13 * sng_84[k]
                   + f_6 * sof0_59[k]
                   - f_7 * sof1_59[k]
                   + f_3 * pc_x[k] * sog_84[k];

        t_115[k] = f_13 * sng_85[k]
                   + f_3 * pc_x[k] * sog_85[k];

        t_116[k] = f_13 * sng_86[k]
                   + f_3 * pc_x[k] * sog_86[k];

        t_117[k] = f_13 * sng_87[k]
                   + f_3 * pc_x[k] * sog_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, pc_z, sng_40, sng_88, sng_89, \
                         sof0_56, sof1_56, sog_85, sog_88, sog_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_13 * sng_88[k]
                   + f_3 * pc_x[k] * sog_88[k];

        t_119[k] = f_13 * sng_89[k]
                   + f_3 * pc_x[k] * sog_89[k];

        t_120[k] = f_1 * sof0_56[k]
                   - f_2 * sof1_56[k]
                   + f_3 * pc_y[k] * sog_85[k];

        t_121[k] = f_10 * sng_40[k]
                   + f_3 * pc_z[k] * sog_85[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, sng_44, sof0_58, sof0_59, \
                         sof1_58, sof1_59, sog_87, sog_88, sog_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_4 * sof0_58[k]
                   - f_5 * sof1_58[k]
                   + f_3 * pc_y[k] * sog_87[k];

        t_123[k] = f_6 * sof0_59[k]
                   - f_7 * sof1_59[k]
                   + f_3 * pc_y[k] * sog_88[k];

        t_124[k] = f_3 * pc_y[k] * sog_89[k];

        t_125[k] = f_10 * sng_44[k]
                   + f_1 * sof0_59[k]
                   - f_2 * sof1_59[k]
                   + f_3 * pc_z[k] * sog_89[k];
    }
}

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snh0,
                                                          const size_t sng, const size_t snh1,
                                                          const size_t sof0, const size_t sof1,
                                                          const size_t sog, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_14 = 4.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_63 = buffer.data(snh0 + 63);
    const auto *snh0_66 = buffer.data(snh0 + 66);
    const auto *snh0_69 = buffer.data(snh0 + 69);
    const auto *snh0_78 = buffer.data(snh0 + 78);
    const auto *snh0_105 = buffer.data(snh0 + 105);
    const auto *snh0_108 = buffer.data(snh0 + 108);
    const auto *snh0_110 = buffer.data(snh0 + 110);
    const auto *snh0_111 = buffer.data(snh0 + 111);
    const auto *snh0_114 = buffer.data(snh0 + 114);
    const auto *snh0_125 = buffer.data(snh0 + 125);
    const auto *snh0_126 = buffer.data(snh0 + 126);
    const auto *snh0_129 = buffer.data(snh0 + 129);
    const auto *snh0_132 = buffer.data(snh0 + 132);

    const auto *sng_45 = buffer.data(sng + 45);
    const auto *sng_47 = buffer.data(sng + 47);
    const auto *sng_48 = buffer.data(sng + 48);
    const auto *sng_50 = buffer.data(sng + 50);
    const auto *sng_55 = buffer.data(sng + 55);
    const auto *sng_57 = buffer.data(sng + 57);
    const auto *sng_58 = buffer.data(sng + 58);
    const auto *sng_59 = buffer.data(sng + 59);
    const auto *sng_60 = buffer.data(sng + 60);
    const auto *sng_62 = buffer.data(sng + 62);
    const auto *sng_63 = buffer.data(sng + 63);
    const auto *sng_65 = buffer.data(sng + 65);
    const auto *sng_70 = buffer.data(sng + 70);
    const auto *sng_72 = buffer.data(sng + 72);
    const auto *sng_73 = buffer.data(sng + 73);
    const auto *sng_74 = buffer.data(sng + 74);
    const auto *sng_75 = buffer.data(sng + 75);
    const auto *sng_76 = buffer.data(sng + 76);
    const auto *sng_77 = buffer.data(sng + 77);
    const auto *sng_78 = buffer.data(sng + 78);
    const auto *sng_80 = buffer.data(sng + 80);
    const auto *sng_85 = buffer.data(sng + 85);
    const auto *sng_87 = buffer.data(sng + 87);
    const auto *sng_88 = buffer.data(sng + 88);
    const auto *sng_89 = buffer.data(sng + 89);
    const auto *sng_90 = buffer.data(sng + 90);
    const auto *sng_92 = buffer.data(sng + 92);
    const auto *sng_93 = buffer.data(sng + 93);
    const auto *sng_95 = buffer.data(sng + 95);
    const auto *sng_96 = buffer.data(sng + 96);
    const auto *sng_99 = buffer.data(sng + 99);
    const auto *sng_100 = buffer.data(sng + 100);
    const auto *sng_101 = buffer.data(sng + 101);
    const auto *sng_102 = buffer.data(sng + 102);
    const auto *sng_103 = buffer.data(sng + 103);
    const auto *sng_104 = buffer.data(sng + 104);
    const auto *sng_105 = buffer.data(sng + 105);
    const auto *sng_107 = buffer.data(sng + 107);
    const auto *sng_110 = buffer.data(sng + 110);
    const auto *sng_114 = buffer.data(sng + 114);
    const auto *sng_115 = buffer.data(sng + 115);
    const auto *sng_116 = buffer.data(sng + 116);
    const auto *sng_117 = buffer.data(sng + 117);
    const auto *sng_118 = buffer.data(sng + 118);
    const auto *sng_119 = buffer.data(sng + 119);
    const auto *sng_130 = buffer.data(sng + 130);
    const auto *sng_131 = buffer.data(sng + 131);
    const auto *sng_132 = buffer.data(sng + 132);
    const auto *sng_133 = buffer.data(sng + 133);
    const auto *sng_134 = buffer.data(sng + 134);
    const auto *sng_135 = buffer.data(sng + 135);
    const auto *sng_138 = buffer.data(sng + 138);
    const auto *sng_140 = buffer.data(sng + 140);
    const auto *sng_141 = buffer.data(sng + 141);
    const auto *sng_144 = buffer.data(sng + 144);
    const auto *sng_145 = buffer.data(sng + 145);
    const auto *sng_146 = buffer.data(sng + 146);
    const auto *sng_147 = buffer.data(sng + 147);
    const auto *sng_148 = buffer.data(sng + 148);
    const auto *sng_149 = buffer.data(sng + 149);
    const auto *sng_150 = buffer.data(sng + 150);
    const auto *sng_153 = buffer.data(sng + 153);
    const auto *sng_155 = buffer.data(sng + 155);
    const auto *sng_156 = buffer.data(sng + 156);
    const auto *sng_159 = buffer.data(sng + 159);
    const auto *sng_160 = buffer.data(sng + 160);
    const auto *sng_161 = buffer.data(sng + 161);
    const auto *sng_162 = buffer.data(sng + 162);
    const auto *sng_163 = buffer.data(sng + 163);
    const auto *sng_164 = buffer.data(sng + 164);
    const auto *sng_170 = buffer.data(sng + 170);
    const auto *sng_174 = buffer.data(sng + 174);
    const auto *sng_175 = buffer.data(sng + 175);
    const auto *sng_176 = buffer.data(sng + 176);

    const auto *snh1_63 = buffer.data(snh1 + 63);
    const auto *snh1_66 = buffer.data(snh1 + 66);
    const auto *snh1_69 = buffer.data(snh1 + 69);
    const auto *snh1_78 = buffer.data(snh1 + 78);
    const auto *snh1_105 = buffer.data(snh1 + 105);
    const auto *snh1_108 = buffer.data(snh1 + 108);
    const auto *snh1_110 = buffer.data(snh1 + 110);
    const auto *snh1_111 = buffer.data(snh1 + 111);
    const auto *snh1_114 = buffer.data(snh1 + 114);
    const auto *snh1_125 = buffer.data(snh1 + 125);
    const auto *snh1_126 = buffer.data(snh1 + 126);
    const auto *snh1_129 = buffer.data(snh1 + 129);
    const auto *snh1_132 = buffer.data(snh1 + 132);

    const auto *sof0_60 = buffer.data(sof0 + 60);
    const auto *sof0_63 = buffer.data(sof0 + 63);
    const auto *sof0_65 = buffer.data(sof0 + 65);
    const auto *sof0_66 = buffer.data(sof0 + 66);
    const auto *sof0_68 = buffer.data(sof0 + 68);
    const auto *sof0_69 = buffer.data(sof0 + 69);
    const auto *sof0_75 = buffer.data(sof0 + 75);
    const auto *sof0_78 = buffer.data(sof0 + 78);
    const auto *sof0_79 = buffer.data(sof0 + 79);
    const auto *sof0_86 = buffer.data(sof0 + 86);
    const auto *sof0_88 = buffer.data(sof0 + 88);
    const auto *sof0_89 = buffer.data(sof0 + 89);
    const auto *sof0_90 = buffer.data(sof0 + 90);
    const auto *sof0_93 = buffer.data(sof0 + 93);
    const auto *sof0_95 = buffer.data(sof0 + 95);
    const auto *sof0_96 = buffer.data(sof0 + 96);
    const auto *sof0_98 = buffer.data(sof0 + 98);
    const auto *sof0_99 = buffer.data(sof0 + 99);
    const auto *sof0_100 = buffer.data(sof0 + 100);
    const auto *sof0_103 = buffer.data(sof0 + 103);
    const auto *sof0_105 = buffer.data(sof0 + 105);
    const auto *sof0_106 = buffer.data(sof0 + 106);
    const auto *sof0_108 = buffer.data(sof0 + 108);
    const auto *sof0_109 = buffer.data(sof0 + 109);
    const auto *sof0_115 = buffer.data(sof0 + 115);
    const auto *sof0_119 = buffer.data(sof0 + 119);

    const auto *sof1_60 = buffer.data(sof1 + 60);
    const auto *sof1_63 = buffer.data(sof1 + 63);
    const auto *sof1_65 = buffer.data(sof1 + 65);
    const auto *sof1_66 = buffer.data(sof1 + 66);
    const auto *sof1_68 = buffer.data(sof1 + 68);
    const auto *sof1_69 = buffer.data(sof1 + 69);
    const auto *sof1_75 = buffer.data(sof1 + 75);
    const auto *sof1_78 = buffer.data(sof1 + 78);
    const auto *sof1_79 = buffer.data(sof1 + 79);
    const auto *sof1_86 = buffer.data(sof1 + 86);
    const auto *sof1_88 = buffer.data(sof1 + 88);
    const auto *sof1_89 = buffer.data(sof1 + 89);
    const auto *sof1_90 = buffer.data(sof1 + 90);
    const auto *sof1_93 = buffer.data(sof1 + 93);
    const auto *sof1_95 = buffer.data(sof1 + 95);
    const auto *sof1_96 = buffer.data(sof1 + 96);
    const auto *sof1_98 = buffer.data(sof1 + 98);
    const auto *sof1_99 = buffer.data(sof1 + 99);
    const auto *sof1_100 = buffer.data(sof1 + 100);
    const auto *sof1_103 = buffer.data(sof1 + 103);
    const auto *sof1_105 = buffer.data(sof1 + 105);
    const auto *sof1_106 = buffer.data(sof1 + 106);
    const auto *sof1_108 = buffer.data(sof1 + 108);
    const auto *sof1_109 = buffer.data(sof1 + 109);
    const auto *sof1_115 = buffer.data(sof1 + 115);
    const auto *sof1_119 = buffer.data(sof1 + 119);

    const auto *sog_90 = buffer.data(sog + 90);
    const auto *sog_92 = buffer.data(sog + 92);
    const auto *sog_93 = buffer.data(sog + 93);
    const auto *sog_95 = buffer.data(sog + 95);
    const auto *sog_96 = buffer.data(sog + 96);
    const auto *sog_99 = buffer.data(sog + 99);
    const auto *sog_100 = buffer.data(sog + 100);
    const auto *sog_101 = buffer.data(sog + 101);
    const auto *sog_102 = buffer.data(sog + 102);
    const auto *sog_103 = buffer.data(sog + 103);
    const auto *sog_104 = buffer.data(sog + 104);
    const auto *sog_105 = buffer.data(sog + 105);
    const auto *sog_107 = buffer.data(sog + 107);
    const auto *sog_108 = buffer.data(sog + 108);
    const auto *sog_110 = buffer.data(sog + 110);
    const auto *sog_114 = buffer.data(sog + 114);
    const auto *sog_115 = buffer.data(sog + 115);
    const auto *sog_116 = buffer.data(sog + 116);
    const auto *sog_117 = buffer.data(sog + 117);
    const auto *sog_118 = buffer.data(sog + 118);
    const auto *sog_119 = buffer.data(sog + 119);
    const auto *sog_120 = buffer.data(sog + 120);
    const auto *sog_122 = buffer.data(sog + 122);
    const auto *sog_123 = buffer.data(sog + 123);
    const auto *sog_125 = buffer.data(sog + 125);
    const auto *sog_130 = buffer.data(sog + 130);
    const auto *sog_131 = buffer.data(sog + 131);
    const auto *sog_132 = buffer.data(sog + 132);
    const auto *sog_133 = buffer.data(sog + 133);
    const auto *sog_134 = buffer.data(sog + 134);
    const auto *sog_135 = buffer.data(sog + 135);
    const auto *sog_137 = buffer.data(sog + 137);
    const auto *sog_138 = buffer.data(sog + 138);
    const auto *sog_140 = buffer.data(sog + 140);
    const auto *sog_141 = buffer.data(sog + 141);
    const auto *sog_144 = buffer.data(sog + 144);
    const auto *sog_145 = buffer.data(sog + 145);
    const auto *sog_146 = buffer.data(sog + 146);
    const auto *sog_147 = buffer.data(sog + 147);
    const auto *sog_148 = buffer.data(sog + 148);
    const auto *sog_149 = buffer.data(sog + 149);
    const auto *sog_150 = buffer.data(sog + 150);
    const auto *sog_152 = buffer.data(sog + 152);
    const auto *sog_153 = buffer.data(sog + 153);
    const auto *sog_155 = buffer.data(sog + 155);
    const auto *sog_156 = buffer.data(sog + 156);
    const auto *sog_159 = buffer.data(sog + 159);
    const auto *sog_160 = buffer.data(sog + 160);
    const auto *sog_161 = buffer.data(sog + 161);
    const auto *sog_162 = buffer.data(sog + 162);
    const auto *sog_163 = buffer.data(sog + 163);
    const auto *sog_164 = buffer.data(sog + 164);
    const auto *sog_165 = buffer.data(sog + 165);
    const auto *sog_167 = buffer.data(sog + 167);
    const auto *sog_168 = buffer.data(sog + 168);
    const auto *sog_170 = buffer.data(sog + 170);
    const auto *sog_174 = buffer.data(sog + 174);
    const auto *sog_175 = buffer.data(sog + 175);
    const auto *sog_176 = buffer.data(sog + 176);

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, sng_45, sng_90, sng_93, \
                         sof0_60, sof0_63, sof1_60, sof1_63, sog_90, \
                         sog_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_14 * sng_90[k]
                   + f_1 * sof0_60[k]
                   - f_2 * sof1_60[k]
                   + f_3 * pc_x[k] * sog_90[k];

        t_127[k] = f_11 * sng_45[k]
                   + f_3 * pc_y[k] * sog_90[k];

        t_128[k] = f_3 * pc_z[k] * sog_90[k];

        t_129[k] = f_14 * sng_93[k]
                   + f_4 * sof0_63[k]
                   - f_5 * sof1_63[k]
                   + f_3 * pc_x[k] * sog_93[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pc_x, pc_y, sng_47, sng_95, sng_96, sof0_65, \
                         sof0_66, sof1_65, sof1_66, sog_92, sog_95, \
                         sog_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_11 * sng_47[k]
                   + f_3 * pc_y[k] * sog_92[k];

        t_131[k] = f_14 * sng_95[k]
                   + f_4 * sof0_65[k]
                   - f_5 * sof1_65[k]
                   + f_3 * pc_x[k] * sog_95[k];

        t_132[k] = f_14 * sng_96[k]
                   + f_6 * sof0_66[k]
                   - f_7 * sof1_66[k]
                   + f_3 * pc_x[k] * sog_96[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, pc_z, sng_50, sng_99, \
                         sng_100, sof0_69, sof1_69, sog_93, sog_95, sog_99, \
                         sog_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_3 * pc_z[k] * sog_93[k];

        t_134[k] = f_11 * sng_50[k]
                   + f_3 * pc_y[k] * sog_95[k];

        t_135[k] = f_14 * sng_99[k]
                   + f_6 * sof0_69[k]
                   - f_7 * sof1_69[k]
                   + f_3 * pc_x[k] * sog_99[k];

        t_136[k] = f_14 * sng_100[k]
                   + f_3 * pc_x[k] * sog_100[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pc_x, sng_101, sng_102, sng_103, sng_104, \
                         sog_101, sog_102, sog_103, sog_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_14 * sng_101[k]
                   + f_3 * pc_x[k] * sog_101[k];

        t_138[k] = f_14 * sng_102[k]
                   + f_3 * pc_x[k] * sog_102[k];

        t_139[k] = f_14 * sng_103[k]
                   + f_3 * pc_x[k] * sog_103[k];

        t_140[k] = f_14 * sng_104[k]
                   + f_3 * pc_x[k] * sog_104[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_y, pc_z, sng_55, sng_57, sof0_66, sof0_68, \
                         sof1_66, sof1_68, sog_100, sog_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_11 * sng_55[k]
                   + f_1 * sof0_66[k]
                   - f_2 * sof1_66[k]
                   + f_3 * pc_y[k] * sog_100[k];

        t_142[k] = f_3 * pc_z[k] * sog_100[k];

        t_143[k] = f_11 * sng_57[k]
                   + f_4 * sof0_68[k]
                   - f_5 * sof1_68[k]
                   + f_3 * pc_y[k] * sog_102[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_z, pc_y, pc_z, snh0_63, sng_58, \
                         sng_59, snh1_63, sof0_69, sof1_69, sog_103, \
                         sog_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_11 * sng_58[k]
                   + f_6 * sof0_69[k]
                   - f_7 * sof1_69[k]
                   + f_3 * pc_y[k] * sog_103[k];

        t_145[k] = f_11 * sng_59[k]
                   + f_3 * pc_y[k] * sog_104[k];

        t_146[k] = f_1 * sof0_69[k]
                   - f_2 * sof1_69[k]
                   + f_3 * pc_z[k] * sog_104[k];

        t_147[k] = pb_z[k] * snh0_63[k]
                   - f_8 * pc_z[k] * snh1_63[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_z, pc_y, pc_z, snh0_66, sng_45, \
                         sng_60, sng_62, snh1_66, sog_105, sog_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * sng_60[k]
                   + f_3 * pc_y[k] * sog_105[k];

        t_149[k] = f_9 * sng_45[k]
                   + f_3 * pc_z[k] * sog_105[k];

        t_150[k] = pb_z[k] * snh0_66[k]
                   - f_8 * pc_z[k] * snh1_66[k];

        t_151[k] = f_10 * sng_62[k]
                   + f_3 * pc_y[k] * sog_107[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pb_z, pc_x, pc_z, snh0_69, sng_48, sng_110, \
                         snh1_69, sof0_75, sof1_75, sog_108, sog_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_14 * sng_110[k]
                   + f_4 * sof0_75[k]
                   - f_5 * sof1_75[k]
                   + f_3 * pc_x[k] * sog_110[k];

        t_153[k] = pb_z[k] * snh0_69[k]
                   - f_8 * pc_z[k] * snh1_69[k];

        t_154[k] = f_9 * sng_48[k]
                   + f_3 * pc_z[k] * sog_108[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_y, sng_65, sng_114, sng_115, \
                         sng_116, sof0_79, sof1_79, sog_110, sog_114, sog_115, \
                         sog_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_10 * sng_65[k]
                   + f_3 * pc_y[k] * sog_110[k];

        t_156[k] = f_14 * sng_114[k]
                   + f_6 * sof0_79[k]
                   - f_7 * sof1_79[k]
                   + f_3 * pc_x[k] * sog_114[k];

        t_157[k] = f_14 * sng_115[k]
                   + f_3 * pc_x[k] * sog_115[k];

        t_158[k] = f_14 * sng_116[k]
                   + f_3 * pc_x[k] * sog_116[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_z, pc_x, pc_z, snh0_78, sng_117, \
                         sng_118, sng_119, snh1_78, sog_117, sog_118, \
                         sog_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_14 * sng_117[k]
                   + f_3 * pc_x[k] * sog_117[k];

        t_160[k] = f_14 * sng_118[k]
                   + f_3 * pc_x[k] * sog_118[k];

        t_161[k] = f_14 * sng_119[k]
                   + f_3 * pc_x[k] * sog_119[k];

        t_162[k] = pb_z[k] * snh0_78[k]
                   - f_8 * pc_z[k] * snh1_78[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pc_y, pc_z, sng_55, sng_72, sng_73, sof0_78, \
                         sof0_79, sof1_78, sof1_79, sog_115, sog_117, \
                         sog_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_9 * sng_55[k]
                   + f_3 * pc_z[k] * sog_115[k];

        t_164[k] = f_10 * sng_72[k]
                   + f_4 * sof0_78[k]
                   - f_5 * sof1_78[k]
                   + f_3 * pc_y[k] * sog_117[k];

        t_165[k] = f_10 * sng_73[k]
                   + f_6 * sof0_79[k]
                   - f_7 * sof1_79[k]
                   + f_3 * pc_y[k] * sog_118[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_y, pc_y, pc_z, snh0_105, sng_59, \
                         sng_74, sng_75, snh1_105, sof0_79, sof1_79, sog_119, \
                         sog_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_10 * sng_74[k]
                   + f_3 * pc_y[k] * sog_119[k];

        t_167[k] = f_9 * sng_59[k]
                   + f_1 * sof0_79[k]
                   - f_2 * sof1_79[k]
                   + f_3 * pc_z[k] * sog_119[k];

        t_168[k] = pb_y[k] * snh0_105[k]
                   - f_8 * pc_y[k] * snh1_105[k];

        t_169[k] = f_9 * sng_75[k]
                   + f_3 * pc_y[k] * sog_120[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_y, pc_y, pc_z, snh0_108, snh0_110, \
                         sng_60, sng_76, sng_77, snh1_108, snh1_110, sog_120, \
                         sog_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_10 * sng_60[k]
                   + f_3 * pc_z[k] * sog_120[k];

        t_171[k] = pb_y[k] * snh0_108[k]
                   + f_10 * sng_76[k]
                   - f_8 * pc_y[k] * snh1_108[k];

        t_172[k] = f_9 * sng_77[k]
                   + f_3 * pc_y[k] * sog_122[k];

        t_173[k] = pb_y[k] * snh0_110[k]
                   - f_8 * pc_y[k] * snh1_110[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pb_y, pc_y, pc_z, snh0_111, snh0_114, \
                         sng_63, sng_78, sng_80, snh1_111, snh1_114, sog_123, \
                         sog_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pb_y[k] * snh0_111[k]
                   + f_11 * sng_78[k]
                   - f_8 * pc_y[k] * snh1_111[k];

        t_175[k] = f_10 * sng_63[k]
                   + f_3 * pc_z[k] * sog_123[k];

        t_176[k] = f_9 * sng_80[k]
                   + f_3 * pc_y[k] * sog_125[k];

        t_177[k] = pb_y[k] * snh0_114[k]
                   - f_8 * pc_y[k] * snh1_114[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pc_x, sng_130, sng_131, sng_132, \
                         sng_133, sng_134, sog_130, sog_131, sog_132, sog_133, \
                         sog_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * sng_130[k]
                   + f_3 * pc_x[k] * sog_130[k];

        t_179[k] = f_14 * sng_131[k]
                   + f_3 * pc_x[k] * sog_131[k];

        t_180[k] = f_14 * sng_132[k]
                   + f_3 * pc_x[k] * sog_132[k];

        t_181[k] = f_14 * sng_133[k]
                   + f_3 * pc_x[k] * sog_133[k];

        t_182[k] = f_14 * sng_134[k]
                   + f_3 * pc_x[k] * sog_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, sng_70, sng_85, sng_87, sof0_86, \
                         sof0_88, sof1_86, sof1_88, sog_130, sog_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * sng_85[k]
                   + f_1 * sof0_86[k]
                   - f_2 * sof1_86[k]
                   + f_3 * pc_y[k] * sog_130[k];

        t_184[k] = f_10 * sng_70[k]
                   + f_3 * pc_z[k] * sog_130[k];

        t_185[k] = f_9 * sng_87[k]
                   + f_4 * sof0_88[k]
                   - f_5 * sof1_88[k]
                   + f_3 * pc_y[k] * sog_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pb_y, pc_y, snh0_125, sng_88, sng_89, snh1_125, \
                         sof0_89, sof1_89, sog_133, sog_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * sng_88[k]
                   + f_6 * sof0_89[k]
                   - f_7 * sof1_89[k]
                   + f_3 * pc_y[k] * sog_133[k];

        t_187[k] = f_9 * sng_89[k]
                   + f_3 * pc_y[k] * sog_134[k];

        t_188[k] = pb_y[k] * snh0_125[k]
                   - f_8 * pc_y[k] * snh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_x, pc_y, pc_z, sng_75, sng_135, \
                         sng_138, sof0_90, sof0_93, sof1_90, sof1_93, sog_135, \
                         sog_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_14 * sng_135[k]
                   + f_1 * sof0_90[k]
                   - f_2 * sof1_90[k]
                   + f_3 * pc_x[k] * sog_135[k];

        t_190[k] = f_3 * pc_y[k] * sog_135[k];

        t_191[k] = f_11 * sng_75[k]
                   + f_3 * pc_z[k] * sog_135[k];

        t_192[k] = f_14 * sng_138[k]
                   + f_4 * sof0_93[k]
                   - f_5 * sof1_93[k]
                   + f_3 * pc_x[k] * sog_138[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, sng_140, sng_141, sof0_95, sof0_96, \
                         sof1_95, sof1_96, sog_137, sog_140, sog_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * sog_137[k];

        t_194[k] = f_14 * sng_140[k]
                   + f_4 * sof0_95[k]
                   - f_5 * sof1_95[k]
                   + f_3 * pc_x[k] * sog_140[k];

        t_195[k] = f_14 * sng_141[k]
                   + f_6 * sof0_96[k]
                   - f_7 * sof1_96[k]
                   + f_3 * pc_x[k] * sog_141[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pc_x, pc_y, pc_z, sng_78, sng_144, \
                         sng_145, sof0_99, sof1_99, sog_138, sog_140, sog_144, \
                         sog_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_11 * sng_78[k]
                   + f_3 * pc_z[k] * sog_138[k];

        t_197[k] = f_3 * pc_y[k] * sog_140[k];

        t_198[k] = f_14 * sng_144[k]
                   + f_6 * sof0_99[k]
                   - f_7 * sof1_99[k]
                   + f_3 * pc_x[k] * sog_144[k];

        t_199[k] = f_14 * sng_145[k]
                   + f_3 * pc_x[k] * sog_145[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, sng_146, sng_147, sng_148, sng_149, \
                         sog_146, sog_147, sog_148, sog_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_14 * sng_146[k]
                   + f_3 * pc_x[k] * sog_146[k];

        t_201[k] = f_14 * sng_147[k]
                   + f_3 * pc_x[k] * sog_147[k];

        t_202[k] = f_14 * sng_148[k]
                   + f_3 * pc_x[k] * sog_148[k];

        t_203[k] = f_14 * sng_149[k]
                   + f_3 * pc_x[k] * sog_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pc_y, pc_z, sng_85, sof0_96, sof0_98, \
                         sof0_99, sof1_96, sof1_98, sof1_99, sog_145, sog_147, \
                         sog_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * sof0_96[k]
                   - f_2 * sof1_96[k]
                   + f_3 * pc_y[k] * sog_145[k];

        t_205[k] = f_11 * sng_85[k]
                   + f_3 * pc_z[k] * sog_145[k];

        t_206[k] = f_4 * sof0_98[k]
                   - f_5 * sof1_98[k]
                   + f_3 * pc_y[k] * sog_147[k];

        t_207[k] = f_6 * sof0_99[k]
                   - f_7 * sof1_99[k]
                   + f_3 * pc_y[k] * sog_148[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_x, pc_y, pc_z, sng_89, sng_90, \
                         sng_150, sof0_99, sof0_100, sof1_99, sof1_100, sog_149, \
                         sog_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_3 * pc_y[k] * sog_149[k];

        t_209[k] = f_11 * sng_89[k]
                   + f_1 * sof0_99[k]
                   - f_2 * sof1_99[k]
                   + f_3 * pc_z[k] * sog_149[k];

        t_210[k] = f_15 * sng_150[k]
                   + f_1 * sof0_100[k]
                   - f_2 * sof1_100[k]
                   + f_3 * pc_x[k] * sog_150[k];

        t_211[k] = f_16 * sng_90[k]
                   + f_3 * pc_y[k] * sog_150[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pc_x, pc_y, pc_z, sng_92, sng_153, sof0_103, \
                         sof1_103, sog_150, sog_152, sog_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_3 * pc_z[k] * sog_150[k];

        t_213[k] = f_15 * sng_153[k]
                   + f_4 * sof0_103[k]
                   - f_5 * sof1_103[k]
                   + f_3 * pc_x[k] * sog_153[k];

        t_214[k] = f_16 * sng_92[k]
                   + f_3 * pc_y[k] * sog_152[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pc_x, pc_z, sng_155, sng_156, sof0_105, \
                         sof0_106, sof1_105, sof1_106, sog_153, sog_155, \
                         sog_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_15 * sng_155[k]
                   + f_4 * sof0_105[k]
                   - f_5 * sof1_105[k]
                   + f_3 * pc_x[k] * sog_155[k];

        t_216[k] = f_15 * sng_156[k]
                   + f_6 * sof0_106[k]
                   - f_7 * sof1_106[k]
                   + f_3 * pc_x[k] * sog_156[k];

        t_217[k] = f_3 * pc_z[k] * sog_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, sng_95, sng_159, sng_160, \
                         sng_161, sof0_109, sof1_109, sog_155, sog_159, sog_160, \
                         sog_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_16 * sng_95[k]
                   + f_3 * pc_y[k] * sog_155[k];

        t_219[k] = f_15 * sng_159[k]
                   + f_6 * sof0_109[k]
                   - f_7 * sof1_109[k]
                   + f_3 * pc_x[k] * sog_159[k];

        t_220[k] = f_15 * sng_160[k]
                   + f_3 * pc_x[k] * sog_160[k];

        t_221[k] = f_15 * sng_161[k]
                   + f_3 * pc_x[k] * sog_161[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, sng_100, sng_162, sng_163, \
                         sng_164, sof0_106, sof1_106, sog_160, sog_162, sog_163, \
                         sog_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_15 * sng_162[k]
                   + f_3 * pc_x[k] * sog_162[k];

        t_223[k] = f_15 * sng_163[k]
                   + f_3 * pc_x[k] * sog_163[k];

        t_224[k] = f_15 * sng_164[k]
                   + f_3 * pc_x[k] * sog_164[k];

        t_225[k] = f_16 * sng_100[k]
                   + f_1 * sof0_106[k]
                   - f_2 * sof1_106[k]
                   + f_3 * pc_y[k] * sog_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_y, pc_z, sng_102, sng_103, sof0_108, \
                         sof0_109, sof1_108, sof1_109, sog_160, sog_162, \
                         sog_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * sog_160[k];

        t_227[k] = f_16 * sng_102[k]
                   + f_4 * sof0_108[k]
                   - f_5 * sof1_108[k]
                   + f_3 * pc_y[k] * sog_162[k];

        t_228[k] = f_16 * sng_103[k]
                   + f_6 * sof0_109[k]
                   - f_7 * sof1_109[k]
                   + f_3 * pc_y[k] * sog_163[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_z, pc_y, pc_z, snh0_126, sng_104, \
                         sng_105, snh1_126, sof0_109, sof1_109, sog_164, \
                         sog_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_16 * sng_104[k]
                   + f_3 * pc_y[k] * sog_164[k];

        t_230[k] = f_1 * sof0_109[k]
                   - f_2 * sof1_109[k]
                   + f_3 * pc_z[k] * sog_164[k];

        t_231[k] = pb_z[k] * snh0_126[k]
                   - f_8 * pc_z[k] * snh1_126[k];

        t_232[k] = f_11 * sng_105[k]
                   + f_3 * pc_y[k] * sog_165[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_z, pc_y, pc_z, snh0_129, sng_90, sng_107, \
                         snh1_129, sog_165, sog_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_9 * sng_90[k]
                   + f_3 * pc_z[k] * sog_165[k];

        t_234[k] = pb_z[k] * snh0_129[k]
                   - f_8 * pc_z[k] * snh1_129[k];

        t_235[k] = f_11 * sng_107[k]
                   + f_3 * pc_y[k] * sog_167[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pb_z, pc_x, pc_z, snh0_132, sng_93, sng_170, \
                         snh1_132, sof0_115, sof1_115, sog_168, \
                         sog_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_15 * sng_170[k]
                   + f_4 * sof0_115[k]
                   - f_5 * sof1_115[k]
                   + f_3 * pc_x[k] * sog_170[k];

        t_237[k] = pb_z[k] * snh0_132[k]
                   - f_8 * pc_z[k] * snh1_132[k];

        t_238[k] = f_9 * sng_93[k]
                   + f_3 * pc_z[k] * sog_168[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pc_x, pc_y, sng_110, sng_174, sng_175, \
                         sng_176, sof0_119, sof1_119, sog_170, sog_174, sog_175, \
                         sog_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_11 * sng_110[k]
                   + f_3 * pc_y[k] * sog_170[k];

        t_240[k] = f_15 * sng_174[k]
                   + f_6 * sof0_119[k]
                   - f_7 * sof1_119[k]
                   + f_3 * pc_x[k] * sog_174[k];

        t_241[k] = f_15 * sng_175[k]
                   + f_3 * pc_x[k] * sog_175[k];

        t_242[k] = f_15 * sng_176[k]
                   + f_3 * pc_x[k] * sog_176[k];
    }
}

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snh0,
                                                          const size_t sng, const size_t snh1,
                                                          const size_t sof0, const size_t sof1,
                                                          const size_t sog, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_141 = buffer.data(snh0 + 141);
    const auto *snh0_189 = buffer.data(snh0 + 189);
    const auto *snh0_192 = buffer.data(snh0 + 192);
    const auto *snh0_194 = buffer.data(snh0 + 194);
    const auto *snh0_195 = buffer.data(snh0 + 195);
    const auto *snh0_198 = buffer.data(snh0 + 198);
    const auto *snh0_209 = buffer.data(snh0 + 209);
    const auto *snh0_210 = buffer.data(snh0 + 210);
    const auto *snh0_213 = buffer.data(snh0 + 213);
    const auto *snh0_216 = buffer.data(snh0 + 216);
    const auto *snh0_225 = buffer.data(snh0 + 225);

    const auto *sng_100 = buffer.data(sng + 100);
    const auto *sng_104 = buffer.data(sng + 104);
    const auto *sng_105 = buffer.data(sng + 105);
    const auto *sng_108 = buffer.data(sng + 108);
    const auto *sng_115 = buffer.data(sng + 115);
    const auto *sng_117 = buffer.data(sng + 117);
    const auto *sng_118 = buffer.data(sng + 118);
    const auto *sng_119 = buffer.data(sng + 119);
    const auto *sng_120 = buffer.data(sng + 120);
    const auto *sng_122 = buffer.data(sng + 122);
    const auto *sng_123 = buffer.data(sng + 123);
    const auto *sng_125 = buffer.data(sng + 125);
    const auto *sng_130 = buffer.data(sng + 130);
    const auto *sng_132 = buffer.data(sng + 132);
    const auto *sng_133 = buffer.data(sng + 133);
    const auto *sng_134 = buffer.data(sng + 134);
    const auto *sng_135 = buffer.data(sng + 135);
    const auto *sng_136 = buffer.data(sng + 136);
    const auto *sng_137 = buffer.data(sng + 137);
    const auto *sng_138 = buffer.data(sng + 138);
    const auto *sng_140 = buffer.data(sng + 140);
    const auto *sng_145 = buffer.data(sng + 145);
    const auto *sng_147 = buffer.data(sng + 147);
    const auto *sng_148 = buffer.data(sng + 148);
    const auto *sng_149 = buffer.data(sng + 149);
    const auto *sng_150 = buffer.data(sng + 150);
    const auto *sng_152 = buffer.data(sng + 152);
    const auto *sng_153 = buffer.data(sng + 153);
    const auto *sng_155 = buffer.data(sng + 155);
    const auto *sng_160 = buffer.data(sng + 160);
    const auto *sng_162 = buffer.data(sng + 162);
    const auto *sng_163 = buffer.data(sng + 163);
    const auto *sng_164 = buffer.data(sng + 164);
    const auto *sng_165 = buffer.data(sng + 165);
    const auto *sng_167 = buffer.data(sng + 167);
    const auto *sng_170 = buffer.data(sng + 170);
    const auto *sng_177 = buffer.data(sng + 177);
    const auto *sng_178 = buffer.data(sng + 178);
    const auto *sng_179 = buffer.data(sng + 179);
    const auto *sng_180 = buffer.data(sng + 180);
    const auto *sng_183 = buffer.data(sng + 183);
    const auto *sng_185 = buffer.data(sng + 185);
    const auto *sng_186 = buffer.data(sng + 186);
    const auto *sng_189 = buffer.data(sng + 189);
    const auto *sng_190 = buffer.data(sng + 190);
    const auto *sng_191 = buffer.data(sng + 191);
    const auto *sng_192 = buffer.data(sng + 192);
    const auto *sng_193 = buffer.data(sng + 193);
    const auto *sng_194 = buffer.data(sng + 194);
    const auto *sng_205 = buffer.data(sng + 205);
    const auto *sng_206 = buffer.data(sng + 206);
    const auto *sng_207 = buffer.data(sng + 207);
    const auto *sng_208 = buffer.data(sng + 208);
    const auto *sng_209 = buffer.data(sng + 209);
    const auto *sng_210 = buffer.data(sng + 210);
    const auto *sng_213 = buffer.data(sng + 213);
    const auto *sng_215 = buffer.data(sng + 215);
    const auto *sng_216 = buffer.data(sng + 216);
    const auto *sng_219 = buffer.data(sng + 219);
    const auto *sng_220 = buffer.data(sng + 220);
    const auto *sng_221 = buffer.data(sng + 221);
    const auto *sng_222 = buffer.data(sng + 222);
    const auto *sng_223 = buffer.data(sng + 223);
    const auto *sng_224 = buffer.data(sng + 224);
    const auto *sng_225 = buffer.data(sng + 225);
    const auto *sng_228 = buffer.data(sng + 228);
    const auto *sng_230 = buffer.data(sng + 230);
    const auto *sng_231 = buffer.data(sng + 231);
    const auto *sng_234 = buffer.data(sng + 234);
    const auto *sng_235 = buffer.data(sng + 235);
    const auto *sng_236 = buffer.data(sng + 236);
    const auto *sng_237 = buffer.data(sng + 237);
    const auto *sng_238 = buffer.data(sng + 238);
    const auto *sng_239 = buffer.data(sng + 239);
    const auto *sng_245 = buffer.data(sng + 245);
    const auto *sng_249 = buffer.data(sng + 249);
    const auto *sng_250 = buffer.data(sng + 250);
    const auto *sng_251 = buffer.data(sng + 251);
    const auto *sng_252 = buffer.data(sng + 252);
    const auto *sng_253 = buffer.data(sng + 253);
    const auto *sng_254 = buffer.data(sng + 254);
    const auto *sng_255 = buffer.data(sng + 255);

    const auto *snh1_141 = buffer.data(snh1 + 141);
    const auto *snh1_189 = buffer.data(snh1 + 189);
    const auto *snh1_192 = buffer.data(snh1 + 192);
    const auto *snh1_194 = buffer.data(snh1 + 194);
    const auto *snh1_195 = buffer.data(snh1 + 195);
    const auto *snh1_198 = buffer.data(snh1 + 198);
    const auto *snh1_209 = buffer.data(snh1 + 209);
    const auto *snh1_210 = buffer.data(snh1 + 210);
    const auto *snh1_213 = buffer.data(snh1 + 213);
    const auto *snh1_216 = buffer.data(snh1 + 216);
    const auto *snh1_225 = buffer.data(snh1 + 225);

    const auto *sof0_118 = buffer.data(sof0 + 118);
    const auto *sof0_119 = buffer.data(sof0 + 119);
    const auto *sof0_120 = buffer.data(sof0 + 120);
    const auto *sof0_123 = buffer.data(sof0 + 123);
    const auto *sof0_125 = buffer.data(sof0 + 125);
    const auto *sof0_126 = buffer.data(sof0 + 126);
    const auto *sof0_128 = buffer.data(sof0 + 128);
    const auto *sof0_129 = buffer.data(sof0 + 129);
    const auto *sof0_136 = buffer.data(sof0 + 136);
    const auto *sof0_138 = buffer.data(sof0 + 138);
    const auto *sof0_139 = buffer.data(sof0 + 139);
    const auto *sof0_140 = buffer.data(sof0 + 140);
    const auto *sof0_143 = buffer.data(sof0 + 143);
    const auto *sof0_145 = buffer.data(sof0 + 145);
    const auto *sof0_146 = buffer.data(sof0 + 146);
    const auto *sof0_148 = buffer.data(sof0 + 148);
    const auto *sof0_149 = buffer.data(sof0 + 149);
    const auto *sof0_150 = buffer.data(sof0 + 150);
    const auto *sof0_153 = buffer.data(sof0 + 153);
    const auto *sof0_155 = buffer.data(sof0 + 155);
    const auto *sof0_156 = buffer.data(sof0 + 156);
    const auto *sof0_158 = buffer.data(sof0 + 158);
    const auto *sof0_159 = buffer.data(sof0 + 159);
    const auto *sof0_165 = buffer.data(sof0 + 165);
    const auto *sof0_168 = buffer.data(sof0 + 168);
    const auto *sof0_169 = buffer.data(sof0 + 169);
    const auto *sof0_170 = buffer.data(sof0 + 170);

    const auto *sof1_118 = buffer.data(sof1 + 118);
    const auto *sof1_119 = buffer.data(sof1 + 119);
    const auto *sof1_120 = buffer.data(sof1 + 120);
    const auto *sof1_123 = buffer.data(sof1 + 123);
    const auto *sof1_125 = buffer.data(sof1 + 125);
    const auto *sof1_126 = buffer.data(sof1 + 126);
    const auto *sof1_128 = buffer.data(sof1 + 128);
    const auto *sof1_129 = buffer.data(sof1 + 129);
    const auto *sof1_136 = buffer.data(sof1 + 136);
    const auto *sof1_138 = buffer.data(sof1 + 138);
    const auto *sof1_139 = buffer.data(sof1 + 139);
    const auto *sof1_140 = buffer.data(sof1 + 140);
    const auto *sof1_143 = buffer.data(sof1 + 143);
    const auto *sof1_145 = buffer.data(sof1 + 145);
    const auto *sof1_146 = buffer.data(sof1 + 146);
    const auto *sof1_148 = buffer.data(sof1 + 148);
    const auto *sof1_149 = buffer.data(sof1 + 149);
    const auto *sof1_150 = buffer.data(sof1 + 150);
    const auto *sof1_153 = buffer.data(sof1 + 153);
    const auto *sof1_155 = buffer.data(sof1 + 155);
    const auto *sof1_156 = buffer.data(sof1 + 156);
    const auto *sof1_158 = buffer.data(sof1 + 158);
    const auto *sof1_159 = buffer.data(sof1 + 159);
    const auto *sof1_165 = buffer.data(sof1 + 165);
    const auto *sof1_168 = buffer.data(sof1 + 168);
    const auto *sof1_169 = buffer.data(sof1 + 169);
    const auto *sof1_170 = buffer.data(sof1 + 170);

    const auto *sog_175 = buffer.data(sog + 175);
    const auto *sog_177 = buffer.data(sog + 177);
    const auto *sog_178 = buffer.data(sog + 178);
    const auto *sog_179 = buffer.data(sog + 179);
    const auto *sog_180 = buffer.data(sog + 180);
    const auto *sog_182 = buffer.data(sog + 182);
    const auto *sog_183 = buffer.data(sog + 183);
    const auto *sog_185 = buffer.data(sog + 185);
    const auto *sog_186 = buffer.data(sog + 186);
    const auto *sog_189 = buffer.data(sog + 189);
    const auto *sog_190 = buffer.data(sog + 190);
    const auto *sog_191 = buffer.data(sog + 191);
    const auto *sog_192 = buffer.data(sog + 192);
    const auto *sog_193 = buffer.data(sog + 193);
    const auto *sog_194 = buffer.data(sog + 194);
    const auto *sog_195 = buffer.data(sog + 195);
    const auto *sog_197 = buffer.data(sog + 197);
    const auto *sog_198 = buffer.data(sog + 198);
    const auto *sog_200 = buffer.data(sog + 200);
    const auto *sog_205 = buffer.data(sog + 205);
    const auto *sog_206 = buffer.data(sog + 206);
    const auto *sog_207 = buffer.data(sog + 207);
    const auto *sog_208 = buffer.data(sog + 208);
    const auto *sog_209 = buffer.data(sog + 209);
    const auto *sog_210 = buffer.data(sog + 210);
    const auto *sog_212 = buffer.data(sog + 212);
    const auto *sog_213 = buffer.data(sog + 213);
    const auto *sog_215 = buffer.data(sog + 215);
    const auto *sog_216 = buffer.data(sog + 216);
    const auto *sog_219 = buffer.data(sog + 219);
    const auto *sog_220 = buffer.data(sog + 220);
    const auto *sog_221 = buffer.data(sog + 221);
    const auto *sog_222 = buffer.data(sog + 222);
    const auto *sog_223 = buffer.data(sog + 223);
    const auto *sog_224 = buffer.data(sog + 224);
    const auto *sog_225 = buffer.data(sog + 225);
    const auto *sog_227 = buffer.data(sog + 227);
    const auto *sog_228 = buffer.data(sog + 228);
    const auto *sog_230 = buffer.data(sog + 230);
    const auto *sog_231 = buffer.data(sog + 231);
    const auto *sog_234 = buffer.data(sog + 234);
    const auto *sog_235 = buffer.data(sog + 235);
    const auto *sog_236 = buffer.data(sog + 236);
    const auto *sog_237 = buffer.data(sog + 237);
    const auto *sog_238 = buffer.data(sog + 238);
    const auto *sog_239 = buffer.data(sog + 239);
    const auto *sog_240 = buffer.data(sog + 240);
    const auto *sog_242 = buffer.data(sog + 242);
    const auto *sog_243 = buffer.data(sog + 243);
    const auto *sog_245 = buffer.data(sog + 245);
    const auto *sog_249 = buffer.data(sog + 249);
    const auto *sog_250 = buffer.data(sog + 250);
    const auto *sog_251 = buffer.data(sog + 251);
    const auto *sog_252 = buffer.data(sog + 252);
    const auto *sog_253 = buffer.data(sog + 253);
    const auto *sog_254 = buffer.data(sog + 254);
    const auto *sog_255 = buffer.data(sog + 255);

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pb_z, pc_x, pc_z, snh0_141, sng_177, \
                         sng_178, sng_179, snh1_141, sog_177, sog_178, \
                         sog_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_15 * sng_177[k]
                   + f_3 * pc_x[k] * sog_177[k];

        t_244[k] = f_15 * sng_178[k]
                   + f_3 * pc_x[k] * sog_178[k];

        t_245[k] = f_15 * sng_179[k]
                   + f_3 * pc_x[k] * sog_179[k];

        t_246[k] = pb_z[k] * snh0_141[k]
                   - f_8 * pc_z[k] * snh1_141[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_y, pc_z, sng_100, sng_117, sng_118, sof0_118, \
                         sof0_119, sof1_118, sof1_119, sog_175, sog_177, \
                         sog_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_9 * sng_100[k]
                   + f_3 * pc_z[k] * sog_175[k];

        t_248[k] = f_11 * sng_117[k]
                   + f_4 * sof0_118[k]
                   - f_5 * sof1_118[k]
                   + f_3 * pc_y[k] * sog_177[k];

        t_249[k] = f_11 * sng_118[k]
                   + f_6 * sof0_119[k]
                   - f_7 * sof1_119[k]
                   + f_3 * pc_y[k] * sog_178[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pc_x, pc_y, pc_z, sng_104, sng_119, sng_180, \
                         sof0_119, sof0_120, sof1_119, sof1_120, sog_179, \
                         sog_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_11 * sng_119[k]
                   + f_3 * pc_y[k] * sog_179[k];

        t_251[k] = f_9 * sng_104[k]
                   + f_1 * sof0_119[k]
                   - f_2 * sof1_119[k]
                   + f_3 * pc_z[k] * sog_179[k];

        t_252[k] = f_15 * sng_180[k]
                   + f_1 * sof0_120[k]
                   - f_2 * sof1_120[k]
                   + f_3 * pc_x[k] * sog_180[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pc_x, pc_y, pc_z, sng_105, sng_120, \
                         sng_122, sng_183, sof0_123, sof1_123, sog_180, sog_182, \
                         sog_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_10 * sng_120[k]
                   + f_3 * pc_y[k] * sog_180[k];

        t_254[k] = f_10 * sng_105[k]
                   + f_3 * pc_z[k] * sog_180[k];

        t_255[k] = f_15 * sng_183[k]
                   + f_4 * sof0_123[k]
                   - f_5 * sof1_123[k]
                   + f_3 * pc_x[k] * sog_183[k];

        t_256[k] = f_10 * sng_122[k]
                   + f_3 * pc_y[k] * sog_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, sng_108, sng_185, sng_186, sof0_125, \
                         sof0_126, sof1_125, sof1_126, sog_183, sog_185, \
                         sog_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_15 * sng_185[k]
                   + f_4 * sof0_125[k]
                   - f_5 * sof1_125[k]
                   + f_3 * pc_x[k] * sog_185[k];

        t_258[k] = f_15 * sng_186[k]
                   + f_6 * sof0_126[k]
                   - f_7 * sof1_126[k]
                   + f_3 * pc_x[k] * sog_186[k];

        t_259[k] = f_10 * sng_108[k]
                   + f_3 * pc_z[k] * sog_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, sng_125, sng_189, sng_190, \
                         sng_191, sof0_129, sof1_129, sog_185, sog_189, sog_190, \
                         sog_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * sng_125[k]
                   + f_3 * pc_y[k] * sog_185[k];

        t_261[k] = f_15 * sng_189[k]
                   + f_6 * sof0_129[k]
                   - f_7 * sof1_129[k]
                   + f_3 * pc_x[k] * sog_189[k];

        t_262[k] = f_15 * sng_190[k]
                   + f_3 * pc_x[k] * sog_190[k];

        t_263[k] = f_15 * sng_191[k]
                   + f_3 * pc_x[k] * sog_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, sng_130, sng_192, sng_193, \
                         sng_194, sof0_126, sof1_126, sog_190, sog_192, sog_193, \
                         sog_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_15 * sng_192[k]
                   + f_3 * pc_x[k] * sog_192[k];

        t_265[k] = f_15 * sng_193[k]
                   + f_3 * pc_x[k] * sog_193[k];

        t_266[k] = f_15 * sng_194[k]
                   + f_3 * pc_x[k] * sog_194[k];

        t_267[k] = f_10 * sng_130[k]
                   + f_1 * sof0_126[k]
                   - f_2 * sof1_126[k]
                   + f_3 * pc_y[k] * sog_190[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_y, pc_z, sng_115, sng_132, sng_133, sof0_128, \
                         sof0_129, sof1_128, sof1_129, sog_190, sog_192, \
                         sog_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * sng_115[k]
                   + f_3 * pc_z[k] * sog_190[k];

        t_269[k] = f_10 * sng_132[k]
                   + f_4 * sof0_128[k]
                   - f_5 * sof1_128[k]
                   + f_3 * pc_y[k] * sog_192[k];

        t_270[k] = f_10 * sng_133[k]
                   + f_6 * sof0_129[k]
                   - f_7 * sof1_129[k]
                   + f_3 * pc_y[k] * sog_193[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pb_y, pc_y, pc_z, snh0_189, sng_119, \
                         sng_134, sng_135, snh1_189, sof0_129, sof1_129, sog_194, \
                         sog_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * sng_134[k]
                   + f_3 * pc_y[k] * sog_194[k];

        t_272[k] = f_10 * sng_119[k]
                   + f_1 * sof0_129[k]
                   - f_2 * sof1_129[k]
                   + f_3 * pc_z[k] * sog_194[k];

        t_273[k] = pb_y[k] * snh0_189[k]
                   - f_8 * pc_y[k] * snh1_189[k];

        t_274[k] = f_9 * sng_135[k]
                   + f_3 * pc_y[k] * sog_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pb_y, pc_y, pc_z, snh0_192, snh0_194, \
                         sng_120, sng_136, sng_137, snh1_192, snh1_194, sog_195, \
                         sog_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * sng_120[k]
                   + f_3 * pc_z[k] * sog_195[k];

        t_276[k] = pb_y[k] * snh0_192[k]
                   + f_10 * sng_136[k]
                   - f_8 * pc_y[k] * snh1_192[k];

        t_277[k] = f_9 * sng_137[k]
                   + f_3 * pc_y[k] * sog_197[k];

        t_278[k] = pb_y[k] * snh0_194[k]
                   - f_8 * pc_y[k] * snh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pb_y, pc_y, pc_z, snh0_195, snh0_198, \
                         sng_123, sng_138, sng_140, snh1_195, snh1_198, sog_198, \
                         sog_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pb_y[k] * snh0_195[k]
                   + f_11 * sng_138[k]
                   - f_8 * pc_y[k] * snh1_195[k];

        t_280[k] = f_11 * sng_123[k]
                   + f_3 * pc_z[k] * sog_198[k];

        t_281[k] = f_9 * sng_140[k]
                   + f_3 * pc_y[k] * sog_200[k];

        t_282[k] = pb_y[k] * snh0_198[k]
                   - f_8 * pc_y[k] * snh1_198[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, sng_205, sng_206, sng_207, \
                         sng_208, sng_209, sog_205, sog_206, sog_207, sog_208, \
                         sog_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_15 * sng_205[k]
                   + f_3 * pc_x[k] * sog_205[k];

        t_284[k] = f_15 * sng_206[k]
                   + f_3 * pc_x[k] * sog_206[k];

        t_285[k] = f_15 * sng_207[k]
                   + f_3 * pc_x[k] * sog_207[k];

        t_286[k] = f_15 * sng_208[k]
                   + f_3 * pc_x[k] * sog_208[k];

        t_287[k] = f_15 * sng_209[k]
                   + f_3 * pc_x[k] * sog_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, sng_130, sng_145, sng_147, sof0_136, \
                         sof0_138, sof1_136, sof1_138, sog_205, \
                         sog_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * sng_145[k]
                   + f_1 * sof0_136[k]
                   - f_2 * sof1_136[k]
                   + f_3 * pc_y[k] * sog_205[k];

        t_289[k] = f_11 * sng_130[k]
                   + f_3 * pc_z[k] * sog_205[k];

        t_290[k] = f_9 * sng_147[k]
                   + f_4 * sof0_138[k]
                   - f_5 * sof1_138[k]
                   + f_3 * pc_y[k] * sog_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pb_y, pc_y, snh0_209, sng_148, sng_149, \
                         snh1_209, sof0_139, sof1_139, sog_208, \
                         sog_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * sng_148[k]
                   + f_6 * sof0_139[k]
                   - f_7 * sof1_139[k]
                   + f_3 * pc_y[k] * sog_208[k];

        t_292[k] = f_9 * sng_149[k]
                   + f_3 * pc_y[k] * sog_209[k];

        t_293[k] = pb_y[k] * snh0_209[k]
                   - f_8 * pc_y[k] * snh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pc_x, pc_y, pc_z, sng_135, sng_210, \
                         sng_213, sof0_140, sof0_143, sof1_140, sof1_143, sog_210, \
                         sog_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_15 * sng_210[k]
                   + f_1 * sof0_140[k]
                   - f_2 * sof1_140[k]
                   + f_3 * pc_x[k] * sog_210[k];

        t_295[k] = f_3 * pc_y[k] * sog_210[k];

        t_296[k] = f_16 * sng_135[k]
                   + f_3 * pc_z[k] * sog_210[k];

        t_297[k] = f_15 * sng_213[k]
                   + f_4 * sof0_143[k]
                   - f_5 * sof1_143[k]
                   + f_3 * pc_x[k] * sog_213[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pc_x, pc_y, sng_215, sng_216, sof0_145, \
                         sof0_146, sof1_145, sof1_146, sog_212, sog_215, \
                         sog_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_3 * pc_y[k] * sog_212[k];

        t_299[k] = f_15 * sng_215[k]
                   + f_4 * sof0_145[k]
                   - f_5 * sof1_145[k]
                   + f_3 * pc_x[k] * sog_215[k];

        t_300[k] = f_15 * sng_216[k]
                   + f_6 * sof0_146[k]
                   - f_7 * sof1_146[k]
                   + f_3 * pc_x[k] * sog_216[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_x, pc_y, pc_z, sng_138, sng_219, \
                         sng_220, sof0_149, sof1_149, sog_213, sog_215, sog_219, \
                         sog_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_16 * sng_138[k]
                   + f_3 * pc_z[k] * sog_213[k];

        t_302[k] = f_3 * pc_y[k] * sog_215[k];

        t_303[k] = f_15 * sng_219[k]
                   + f_6 * sof0_149[k]
                   - f_7 * sof1_149[k]
                   + f_3 * pc_x[k] * sog_219[k];

        t_304[k] = f_15 * sng_220[k]
                   + f_3 * pc_x[k] * sog_220[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, sng_221, sng_222, sng_223, sng_224, \
                         sog_221, sog_222, sog_223, sog_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_15 * sng_221[k]
                   + f_3 * pc_x[k] * sog_221[k];

        t_306[k] = f_15 * sng_222[k]
                   + f_3 * pc_x[k] * sog_222[k];

        t_307[k] = f_15 * sng_223[k]
                   + f_3 * pc_x[k] * sog_223[k];

        t_308[k] = f_15 * sng_224[k]
                   + f_3 * pc_x[k] * sog_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pc_y, pc_z, sng_145, sof0_146, sof0_148, \
                         sof0_149, sof1_146, sof1_148, sof1_149, sog_220, sog_222, \
                         sog_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_1 * sof0_146[k]
                   - f_2 * sof1_146[k]
                   + f_3 * pc_y[k] * sog_220[k];

        t_310[k] = f_16 * sng_145[k]
                   + f_3 * pc_z[k] * sog_220[k];

        t_311[k] = f_4 * sof0_148[k]
                   - f_5 * sof1_148[k]
                   + f_3 * pc_y[k] * sog_222[k];

        t_312[k] = f_6 * sof0_149[k]
                   - f_7 * sof1_149[k]
                   + f_3 * pc_y[k] * sog_223[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pc_x, pc_y, pc_z, sng_149, sng_150, \
                         sng_225, sof0_149, sof0_150, sof1_149, sof1_150, sog_224, \
                         sog_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_3 * pc_y[k] * sog_224[k];

        t_314[k] = f_16 * sng_149[k]
                   + f_1 * sof0_149[k]
                   - f_2 * sof1_149[k]
                   + f_3 * pc_z[k] * sog_224[k];

        t_315[k] = f_17 * sng_225[k]
                   + f_1 * sof0_150[k]
                   - f_2 * sof1_150[k]
                   + f_3 * pc_x[k] * sog_225[k];

        t_316[k] = f_18 * sng_150[k]
                   + f_3 * pc_y[k] * sog_225[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, pc_x, pc_y, pc_z, sng_152, sng_228, sof0_153, \
                         sof1_153, sog_225, sog_227, sog_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_3 * pc_z[k] * sog_225[k];

        t_318[k] = f_17 * sng_228[k]
                   + f_4 * sof0_153[k]
                   - f_5 * sof1_153[k]
                   + f_3 * pc_x[k] * sog_228[k];

        t_319[k] = f_18 * sng_152[k]
                   + f_3 * pc_y[k] * sog_227[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, pc_x, pc_z, sng_230, sng_231, sof0_155, \
                         sof0_156, sof1_155, sof1_156, sog_228, sog_230, \
                         sog_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_17 * sng_230[k]
                   + f_4 * sof0_155[k]
                   - f_5 * sof1_155[k]
                   + f_3 * pc_x[k] * sog_230[k];

        t_321[k] = f_17 * sng_231[k]
                   + f_6 * sof0_156[k]
                   - f_7 * sof1_156[k]
                   + f_3 * pc_x[k] * sog_231[k];

        t_322[k] = f_3 * pc_z[k] * sog_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pc_x, pc_y, sng_155, sng_234, sng_235, \
                         sng_236, sof0_159, sof1_159, sog_230, sog_234, sog_235, \
                         sog_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_18 * sng_155[k]
                   + f_3 * pc_y[k] * sog_230[k];

        t_324[k] = f_17 * sng_234[k]
                   + f_6 * sof0_159[k]
                   - f_7 * sof1_159[k]
                   + f_3 * pc_x[k] * sog_234[k];

        t_325[k] = f_17 * sng_235[k]
                   + f_3 * pc_x[k] * sog_235[k];

        t_326[k] = f_17 * sng_236[k]
                   + f_3 * pc_x[k] * sog_236[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pc_x, pc_y, sng_160, sng_237, sng_238, \
                         sng_239, sof0_156, sof1_156, sog_235, sog_237, sog_238, \
                         sog_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_17 * sng_237[k]
                   + f_3 * pc_x[k] * sog_237[k];

        t_328[k] = f_17 * sng_238[k]
                   + f_3 * pc_x[k] * sog_238[k];

        t_329[k] = f_17 * sng_239[k]
                   + f_3 * pc_x[k] * sog_239[k];

        t_330[k] = f_18 * sng_160[k]
                   + f_1 * sof0_156[k]
                   - f_2 * sof1_156[k]
                   + f_3 * pc_y[k] * sog_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pc_y, pc_z, sng_162, sng_163, sof0_158, \
                         sof0_159, sof1_158, sof1_159, sog_235, sog_237, \
                         sog_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * sog_235[k];

        t_332[k] = f_18 * sng_162[k]
                   + f_4 * sof0_158[k]
                   - f_5 * sof1_158[k]
                   + f_3 * pc_y[k] * sog_237[k];

        t_333[k] = f_18 * sng_163[k]
                   + f_6 * sof0_159[k]
                   - f_7 * sof1_159[k]
                   + f_3 * pc_y[k] * sog_238[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pb_z, pc_y, pc_z, snh0_210, sng_164, \
                         sng_165, snh1_210, sof0_159, sof1_159, sog_239, \
                         sog_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_18 * sng_164[k]
                   + f_3 * pc_y[k] * sog_239[k];

        t_335[k] = f_1 * sof0_159[k]
                   - f_2 * sof1_159[k]
                   + f_3 * pc_z[k] * sog_239[k];

        t_336[k] = pb_z[k] * snh0_210[k]
                   - f_8 * pc_z[k] * snh1_210[k];

        t_337[k] = f_16 * sng_165[k]
                   + f_3 * pc_y[k] * sog_240[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pb_z, pc_y, pc_z, snh0_213, sng_150, sng_167, \
                         snh1_213, sog_240, sog_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_9 * sng_150[k]
                   + f_3 * pc_z[k] * sog_240[k];

        t_339[k] = pb_z[k] * snh0_213[k]
                   - f_8 * pc_z[k] * snh1_213[k];

        t_340[k] = f_16 * sng_167[k]
                   + f_3 * pc_y[k] * sog_242[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pb_z, pc_x, pc_z, snh0_216, sng_153, sng_245, \
                         snh1_216, sof0_165, sof1_165, sog_243, \
                         sog_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_17 * sng_245[k]
                   + f_4 * sof0_165[k]
                   - f_5 * sof1_165[k]
                   + f_3 * pc_x[k] * sog_245[k];

        t_342[k] = pb_z[k] * snh0_216[k]
                   - f_8 * pc_z[k] * snh1_216[k];

        t_343[k] = f_9 * sng_153[k]
                   + f_3 * pc_z[k] * sog_243[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pc_x, pc_y, sng_170, sng_249, sng_250, \
                         sng_251, sof0_169, sof1_169, sog_245, sog_249, sog_250, \
                         sog_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_16 * sng_170[k]
                   + f_3 * pc_y[k] * sog_245[k];

        t_345[k] = f_17 * sng_249[k]
                   + f_6 * sof0_169[k]
                   - f_7 * sof1_169[k]
                   + f_3 * pc_x[k] * sog_249[k];

        t_346[k] = f_17 * sng_250[k]
                   + f_3 * pc_x[k] * sog_250[k];

        t_347[k] = f_17 * sng_251[k]
                   + f_3 * pc_x[k] * sog_251[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pb_z, pc_x, pc_z, snh0_225, sng_252, \
                         sng_253, sng_254, snh1_225, sog_252, sog_253, \
                         sog_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_17 * sng_252[k]
                   + f_3 * pc_x[k] * sog_252[k];

        t_349[k] = f_17 * sng_253[k]
                   + f_3 * pc_x[k] * sog_253[k];

        t_350[k] = f_17 * sng_254[k]
                   + f_3 * pc_x[k] * sog_254[k];

        t_351[k] = pb_z[k] * snh0_225[k]
                   - f_8 * pc_z[k] * snh1_225[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, pc_y, pc_z, sng_160, sng_177, sng_178, sof0_168, \
                         sof0_169, sof1_168, sof1_169, sog_250, sog_252, \
                         sog_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_9 * sng_160[k]
                   + f_3 * pc_z[k] * sog_250[k];

        t_353[k] = f_16 * sng_177[k]
                   + f_4 * sof0_168[k]
                   - f_5 * sof1_168[k]
                   + f_3 * pc_y[k] * sog_252[k];

        t_354[k] = f_16 * sng_178[k]
                   + f_6 * sof0_169[k]
                   - f_7 * sof1_169[k]
                   + f_3 * pc_y[k] * sog_253[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, pc_x, pc_y, pc_z, sng_164, sng_179, sng_255, \
                         sof0_169, sof0_170, sof1_169, sof1_170, sog_254, \
                         sog_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_16 * sng_179[k]
                   + f_3 * pc_y[k] * sog_254[k];

        t_356[k] = f_9 * sng_164[k]
                   + f_1 * sof0_169[k]
                   - f_2 * sof1_169[k]
                   + f_3 * pc_z[k] * sog_254[k];

        t_357[k] = f_17 * sng_255[k]
                   + f_1 * sof0_170[k]
                   - f_2 * sof1_170[k]
                   + f_3 * pc_x[k] * sog_255[k];
    }
}

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snh0,
                                                          const size_t sng, const size_t snh1,
                                                          const size_t sof0, const size_t sof1,
                                                          const size_t sog, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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
    auto *t_471 = buffer.data(target + 471);
    auto *t_472 = buffer.data(target + 472);
    auto *t_473 = buffer.data(target + 473);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_294 = buffer.data(snh0 + 294);
    const auto *snh0_297 = buffer.data(snh0 + 297);
    const auto *snh0_299 = buffer.data(snh0 + 299);
    const auto *snh0_300 = buffer.data(snh0 + 300);
    const auto *snh0_303 = buffer.data(snh0 + 303);
    const auto *snh0_314 = buffer.data(snh0 + 314);
    const auto *snh0_315 = buffer.data(snh0 + 315);
    const auto *snh0_318 = buffer.data(snh0 + 318);
    const auto *snh0_321 = buffer.data(snh0 + 321);

    const auto *sng_165 = buffer.data(sng + 165);
    const auto *sng_168 = buffer.data(sng + 168);
    const auto *sng_175 = buffer.data(sng + 175);
    const auto *sng_179 = buffer.data(sng + 179);
    const auto *sng_180 = buffer.data(sng + 180);
    const auto *sng_182 = buffer.data(sng + 182);
    const auto *sng_183 = buffer.data(sng + 183);
    const auto *sng_185 = buffer.data(sng + 185);
    const auto *sng_190 = buffer.data(sng + 190);
    const auto *sng_192 = buffer.data(sng + 192);
    const auto *sng_193 = buffer.data(sng + 193);
    const auto *sng_194 = buffer.data(sng + 194);
    const auto *sng_195 = buffer.data(sng + 195);
    const auto *sng_197 = buffer.data(sng + 197);
    const auto *sng_198 = buffer.data(sng + 198);
    const auto *sng_200 = buffer.data(sng + 200);
    const auto *sng_205 = buffer.data(sng + 205);
    const auto *sng_207 = buffer.data(sng + 207);
    const auto *sng_208 = buffer.data(sng + 208);
    const auto *sng_209 = buffer.data(sng + 209);
    const auto *sng_210 = buffer.data(sng + 210);
    const auto *sng_211 = buffer.data(sng + 211);
    const auto *sng_212 = buffer.data(sng + 212);
    const auto *sng_213 = buffer.data(sng + 213);
    const auto *sng_215 = buffer.data(sng + 215);
    const auto *sng_220 = buffer.data(sng + 220);
    const auto *sng_222 = buffer.data(sng + 222);
    const auto *sng_223 = buffer.data(sng + 223);
    const auto *sng_224 = buffer.data(sng + 224);
    const auto *sng_225 = buffer.data(sng + 225);
    const auto *sng_227 = buffer.data(sng + 227);
    const auto *sng_228 = buffer.data(sng + 228);
    const auto *sng_230 = buffer.data(sng + 230);
    const auto *sng_235 = buffer.data(sng + 235);
    const auto *sng_237 = buffer.data(sng + 237);
    const auto *sng_238 = buffer.data(sng + 238);
    const auto *sng_239 = buffer.data(sng + 239);
    const auto *sng_240 = buffer.data(sng + 240);
    const auto *sng_242 = buffer.data(sng + 242);
    const auto *sng_245 = buffer.data(sng + 245);
    const auto *sng_258 = buffer.data(sng + 258);
    const auto *sng_260 = buffer.data(sng + 260);
    const auto *sng_261 = buffer.data(sng + 261);
    const auto *sng_264 = buffer.data(sng + 264);
    const auto *sng_265 = buffer.data(sng + 265);
    const auto *sng_266 = buffer.data(sng + 266);
    const auto *sng_267 = buffer.data(sng + 267);
    const auto *sng_268 = buffer.data(sng + 268);
    const auto *sng_269 = buffer.data(sng + 269);
    const auto *sng_270 = buffer.data(sng + 270);
    const auto *sng_273 = buffer.data(sng + 273);
    const auto *sng_275 = buffer.data(sng + 275);
    const auto *sng_276 = buffer.data(sng + 276);
    const auto *sng_279 = buffer.data(sng + 279);
    const auto *sng_280 = buffer.data(sng + 280);
    const auto *sng_281 = buffer.data(sng + 281);
    const auto *sng_282 = buffer.data(sng + 282);
    const auto *sng_283 = buffer.data(sng + 283);
    const auto *sng_284 = buffer.data(sng + 284);
    const auto *sng_295 = buffer.data(sng + 295);
    const auto *sng_296 = buffer.data(sng + 296);
    const auto *sng_297 = buffer.data(sng + 297);
    const auto *sng_298 = buffer.data(sng + 298);
    const auto *sng_299 = buffer.data(sng + 299);
    const auto *sng_300 = buffer.data(sng + 300);
    const auto *sng_303 = buffer.data(sng + 303);
    const auto *sng_305 = buffer.data(sng + 305);
    const auto *sng_306 = buffer.data(sng + 306);
    const auto *sng_309 = buffer.data(sng + 309);
    const auto *sng_310 = buffer.data(sng + 310);
    const auto *sng_311 = buffer.data(sng + 311);
    const auto *sng_312 = buffer.data(sng + 312);
    const auto *sng_313 = buffer.data(sng + 313);
    const auto *sng_314 = buffer.data(sng + 314);
    const auto *sng_315 = buffer.data(sng + 315);
    const auto *sng_318 = buffer.data(sng + 318);
    const auto *sng_320 = buffer.data(sng + 320);
    const auto *sng_321 = buffer.data(sng + 321);
    const auto *sng_324 = buffer.data(sng + 324);
    const auto *sng_325 = buffer.data(sng + 325);
    const auto *sng_326 = buffer.data(sng + 326);
    const auto *sng_327 = buffer.data(sng + 327);
    const auto *sng_328 = buffer.data(sng + 328);
    const auto *sng_329 = buffer.data(sng + 329);
    const auto *sng_335 = buffer.data(sng + 335);
    const auto *sng_339 = buffer.data(sng + 339);
    const auto *sng_340 = buffer.data(sng + 340);
    const auto *sng_341 = buffer.data(sng + 341);

    const auto *snh1_294 = buffer.data(snh1 + 294);
    const auto *snh1_297 = buffer.data(snh1 + 297);
    const auto *snh1_299 = buffer.data(snh1 + 299);
    const auto *snh1_300 = buffer.data(snh1 + 300);
    const auto *snh1_303 = buffer.data(snh1 + 303);
    const auto *snh1_314 = buffer.data(snh1 + 314);
    const auto *snh1_315 = buffer.data(snh1 + 315);
    const auto *snh1_318 = buffer.data(snh1 + 318);
    const auto *snh1_321 = buffer.data(snh1 + 321);

    const auto *sof0_173 = buffer.data(sof0 + 173);
    const auto *sof0_175 = buffer.data(sof0 + 175);
    const auto *sof0_176 = buffer.data(sof0 + 176);
    const auto *sof0_178 = buffer.data(sof0 + 178);
    const auto *sof0_179 = buffer.data(sof0 + 179);
    const auto *sof0_180 = buffer.data(sof0 + 180);
    const auto *sof0_183 = buffer.data(sof0 + 183);
    const auto *sof0_185 = buffer.data(sof0 + 185);
    const auto *sof0_186 = buffer.data(sof0 + 186);
    const auto *sof0_188 = buffer.data(sof0 + 188);
    const auto *sof0_189 = buffer.data(sof0 + 189);
    const auto *sof0_196 = buffer.data(sof0 + 196);
    const auto *sof0_198 = buffer.data(sof0 + 198);
    const auto *sof0_199 = buffer.data(sof0 + 199);
    const auto *sof0_200 = buffer.data(sof0 + 200);
    const auto *sof0_203 = buffer.data(sof0 + 203);
    const auto *sof0_205 = buffer.data(sof0 + 205);
    const auto *sof0_206 = buffer.data(sof0 + 206);
    const auto *sof0_208 = buffer.data(sof0 + 208);
    const auto *sof0_209 = buffer.data(sof0 + 209);
    const auto *sof0_210 = buffer.data(sof0 + 210);
    const auto *sof0_213 = buffer.data(sof0 + 213);
    const auto *sof0_215 = buffer.data(sof0 + 215);
    const auto *sof0_216 = buffer.data(sof0 + 216);
    const auto *sof0_218 = buffer.data(sof0 + 218);
    const auto *sof0_219 = buffer.data(sof0 + 219);
    const auto *sof0_225 = buffer.data(sof0 + 225);
    const auto *sof0_229 = buffer.data(sof0 + 229);

    const auto *sof1_173 = buffer.data(sof1 + 173);
    const auto *sof1_175 = buffer.data(sof1 + 175);
    const auto *sof1_176 = buffer.data(sof1 + 176);
    const auto *sof1_178 = buffer.data(sof1 + 178);
    const auto *sof1_179 = buffer.data(sof1 + 179);
    const auto *sof1_180 = buffer.data(sof1 + 180);
    const auto *sof1_183 = buffer.data(sof1 + 183);
    const auto *sof1_185 = buffer.data(sof1 + 185);
    const auto *sof1_186 = buffer.data(sof1 + 186);
    const auto *sof1_188 = buffer.data(sof1 + 188);
    const auto *sof1_189 = buffer.data(sof1 + 189);
    const auto *sof1_196 = buffer.data(sof1 + 196);
    const auto *sof1_198 = buffer.data(sof1 + 198);
    const auto *sof1_199 = buffer.data(sof1 + 199);
    const auto *sof1_200 = buffer.data(sof1 + 200);
    const auto *sof1_203 = buffer.data(sof1 + 203);
    const auto *sof1_205 = buffer.data(sof1 + 205);
    const auto *sof1_206 = buffer.data(sof1 + 206);
    const auto *sof1_208 = buffer.data(sof1 + 208);
    const auto *sof1_209 = buffer.data(sof1 + 209);
    const auto *sof1_210 = buffer.data(sof1 + 210);
    const auto *sof1_213 = buffer.data(sof1 + 213);
    const auto *sof1_215 = buffer.data(sof1 + 215);
    const auto *sof1_216 = buffer.data(sof1 + 216);
    const auto *sof1_218 = buffer.data(sof1 + 218);
    const auto *sof1_219 = buffer.data(sof1 + 219);
    const auto *sof1_225 = buffer.data(sof1 + 225);
    const auto *sof1_229 = buffer.data(sof1 + 229);

    const auto *sog_255 = buffer.data(sog + 255);
    const auto *sog_257 = buffer.data(sog + 257);
    const auto *sog_258 = buffer.data(sog + 258);
    const auto *sog_260 = buffer.data(sog + 260);
    const auto *sog_261 = buffer.data(sog + 261);
    const auto *sog_264 = buffer.data(sog + 264);
    const auto *sog_265 = buffer.data(sog + 265);
    const auto *sog_266 = buffer.data(sog + 266);
    const auto *sog_267 = buffer.data(sog + 267);
    const auto *sog_268 = buffer.data(sog + 268);
    const auto *sog_269 = buffer.data(sog + 269);
    const auto *sog_270 = buffer.data(sog + 270);
    const auto *sog_272 = buffer.data(sog + 272);
    const auto *sog_273 = buffer.data(sog + 273);
    const auto *sog_275 = buffer.data(sog + 275);
    const auto *sog_276 = buffer.data(sog + 276);
    const auto *sog_279 = buffer.data(sog + 279);
    const auto *sog_280 = buffer.data(sog + 280);
    const auto *sog_281 = buffer.data(sog + 281);
    const auto *sog_282 = buffer.data(sog + 282);
    const auto *sog_283 = buffer.data(sog + 283);
    const auto *sog_284 = buffer.data(sog + 284);
    const auto *sog_285 = buffer.data(sog + 285);
    const auto *sog_287 = buffer.data(sog + 287);
    const auto *sog_288 = buffer.data(sog + 288);
    const auto *sog_290 = buffer.data(sog + 290);
    const auto *sog_295 = buffer.data(sog + 295);
    const auto *sog_296 = buffer.data(sog + 296);
    const auto *sog_297 = buffer.data(sog + 297);
    const auto *sog_298 = buffer.data(sog + 298);
    const auto *sog_299 = buffer.data(sog + 299);
    const auto *sog_300 = buffer.data(sog + 300);
    const auto *sog_302 = buffer.data(sog + 302);
    const auto *sog_303 = buffer.data(sog + 303);
    const auto *sog_305 = buffer.data(sog + 305);
    const auto *sog_306 = buffer.data(sog + 306);
    const auto *sog_309 = buffer.data(sog + 309);
    const auto *sog_310 = buffer.data(sog + 310);
    const auto *sog_311 = buffer.data(sog + 311);
    const auto *sog_312 = buffer.data(sog + 312);
    const auto *sog_313 = buffer.data(sog + 313);
    const auto *sog_314 = buffer.data(sog + 314);
    const auto *sog_315 = buffer.data(sog + 315);
    const auto *sog_317 = buffer.data(sog + 317);
    const auto *sog_318 = buffer.data(sog + 318);
    const auto *sog_320 = buffer.data(sog + 320);
    const auto *sog_321 = buffer.data(sog + 321);
    const auto *sog_324 = buffer.data(sog + 324);
    const auto *sog_325 = buffer.data(sog + 325);
    const auto *sog_326 = buffer.data(sog + 326);
    const auto *sog_327 = buffer.data(sog + 327);
    const auto *sog_328 = buffer.data(sog + 328);
    const auto *sog_329 = buffer.data(sog + 329);
    const auto *sog_330 = buffer.data(sog + 330);
    const auto *sog_332 = buffer.data(sog + 332);
    const auto *sog_333 = buffer.data(sog + 333);
    const auto *sog_335 = buffer.data(sog + 335);
    const auto *sog_339 = buffer.data(sog + 339);
    const auto *sog_340 = buffer.data(sog + 340);
    const auto *sog_341 = buffer.data(sog + 341);

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pc_x, pc_y, pc_z, sng_165, sng_180, \
                         sng_182, sng_258, sof0_173, sof1_173, sog_255, sog_257, \
                         sog_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_11 * sng_180[k]
                   + f_3 * pc_y[k] * sog_255[k];

        t_359[k] = f_10 * sng_165[k]
                   + f_3 * pc_z[k] * sog_255[k];

        t_360[k] = f_17 * sng_258[k]
                   + f_4 * sof0_173[k]
                   - f_5 * sof1_173[k]
                   + f_3 * pc_x[k] * sog_258[k];

        t_361[k] = f_11 * sng_182[k]
                   + f_3 * pc_y[k] * sog_257[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pc_x, pc_z, sng_168, sng_260, sng_261, sof0_175, \
                         sof0_176, sof1_175, sof1_176, sog_258, sog_260, \
                         sog_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_17 * sng_260[k]
                   + f_4 * sof0_175[k]
                   - f_5 * sof1_175[k]
                   + f_3 * pc_x[k] * sog_260[k];

        t_363[k] = f_17 * sng_261[k]
                   + f_6 * sof0_176[k]
                   - f_7 * sof1_176[k]
                   + f_3 * pc_x[k] * sog_261[k];

        t_364[k] = f_10 * sng_168[k]
                   + f_3 * pc_z[k] * sog_258[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, pc_y, sng_185, sng_264, sng_265, \
                         sng_266, sof0_179, sof1_179, sog_260, sog_264, sog_265, \
                         sog_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_11 * sng_185[k]
                   + f_3 * pc_y[k] * sog_260[k];

        t_366[k] = f_17 * sng_264[k]
                   + f_6 * sof0_179[k]
                   - f_7 * sof1_179[k]
                   + f_3 * pc_x[k] * sog_264[k];

        t_367[k] = f_17 * sng_265[k]
                   + f_3 * pc_x[k] * sog_265[k];

        t_368[k] = f_17 * sng_266[k]
                   + f_3 * pc_x[k] * sog_266[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pc_x, pc_y, sng_190, sng_267, sng_268, \
                         sng_269, sof0_176, sof1_176, sog_265, sog_267, sog_268, \
                         sog_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_17 * sng_267[k]
                   + f_3 * pc_x[k] * sog_267[k];

        t_370[k] = f_17 * sng_268[k]
                   + f_3 * pc_x[k] * sog_268[k];

        t_371[k] = f_17 * sng_269[k]
                   + f_3 * pc_x[k] * sog_269[k];

        t_372[k] = f_11 * sng_190[k]
                   + f_1 * sof0_176[k]
                   - f_2 * sof1_176[k]
                   + f_3 * pc_y[k] * sog_265[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_y, pc_z, sng_175, sng_192, sng_193, sof0_178, \
                         sof0_179, sof1_178, sof1_179, sog_265, sog_267, \
                         sog_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_10 * sng_175[k]
                   + f_3 * pc_z[k] * sog_265[k];

        t_374[k] = f_11 * sng_192[k]
                   + f_4 * sof0_178[k]
                   - f_5 * sof1_178[k]
                   + f_3 * pc_y[k] * sog_267[k];

        t_375[k] = f_11 * sng_193[k]
                   + f_6 * sof0_179[k]
                   - f_7 * sof1_179[k]
                   + f_3 * pc_y[k] * sog_268[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_y, pc_z, sng_179, sng_194, sng_270, \
                         sof0_179, sof0_180, sof1_179, sof1_180, sog_269, \
                         sog_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_11 * sng_194[k]
                   + f_3 * pc_y[k] * sog_269[k];

        t_377[k] = f_10 * sng_179[k]
                   + f_1 * sof0_179[k]
                   - f_2 * sof1_179[k]
                   + f_3 * pc_z[k] * sog_269[k];

        t_378[k] = f_17 * sng_270[k]
                   + f_1 * sof0_180[k]
                   - f_2 * sof1_180[k]
                   + f_3 * pc_x[k] * sog_270[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, pc_z, sng_180, sng_195, \
                         sng_197, sng_273, sof0_183, sof1_183, sog_270, sog_272, \
                         sog_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_10 * sng_195[k]
                   + f_3 * pc_y[k] * sog_270[k];

        t_380[k] = f_11 * sng_180[k]
                   + f_3 * pc_z[k] * sog_270[k];

        t_381[k] = f_17 * sng_273[k]
                   + f_4 * sof0_183[k]
                   - f_5 * sof1_183[k]
                   + f_3 * pc_x[k] * sog_273[k];

        t_382[k] = f_10 * sng_197[k]
                   + f_3 * pc_y[k] * sog_272[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pc_x, pc_z, sng_183, sng_275, sng_276, sof0_185, \
                         sof0_186, sof1_185, sof1_186, sog_273, sog_275, \
                         sog_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_17 * sng_275[k]
                   + f_4 * sof0_185[k]
                   - f_5 * sof1_185[k]
                   + f_3 * pc_x[k] * sog_275[k];

        t_384[k] = f_17 * sng_276[k]
                   + f_6 * sof0_186[k]
                   - f_7 * sof1_186[k]
                   + f_3 * pc_x[k] * sog_276[k];

        t_385[k] = f_11 * sng_183[k]
                   + f_3 * pc_z[k] * sog_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pc_x, pc_y, sng_200, sng_279, sng_280, \
                         sng_281, sof0_189, sof1_189, sog_275, sog_279, sog_280, \
                         sog_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_10 * sng_200[k]
                   + f_3 * pc_y[k] * sog_275[k];

        t_387[k] = f_17 * sng_279[k]
                   + f_6 * sof0_189[k]
                   - f_7 * sof1_189[k]
                   + f_3 * pc_x[k] * sog_279[k];

        t_388[k] = f_17 * sng_280[k]
                   + f_3 * pc_x[k] * sog_280[k];

        t_389[k] = f_17 * sng_281[k]
                   + f_3 * pc_x[k] * sog_281[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pc_x, pc_y, sng_205, sng_282, sng_283, \
                         sng_284, sof0_186, sof1_186, sog_280, sog_282, sog_283, \
                         sog_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_17 * sng_282[k]
                   + f_3 * pc_x[k] * sog_282[k];

        t_391[k] = f_17 * sng_283[k]
                   + f_3 * pc_x[k] * sog_283[k];

        t_392[k] = f_17 * sng_284[k]
                   + f_3 * pc_x[k] * sog_284[k];

        t_393[k] = f_10 * sng_205[k]
                   + f_1 * sof0_186[k]
                   - f_2 * sof1_186[k]
                   + f_3 * pc_y[k] * sog_280[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pc_y, pc_z, sng_190, sng_207, sng_208, sof0_188, \
                         sof0_189, sof1_188, sof1_189, sog_280, sog_282, \
                         sog_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_11 * sng_190[k]
                   + f_3 * pc_z[k] * sog_280[k];

        t_395[k] = f_10 * sng_207[k]
                   + f_4 * sof0_188[k]
                   - f_5 * sof1_188[k]
                   + f_3 * pc_y[k] * sog_282[k];

        t_396[k] = f_10 * sng_208[k]
                   + f_6 * sof0_189[k]
                   - f_7 * sof1_189[k]
                   + f_3 * pc_y[k] * sog_283[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pb_y, pc_y, pc_z, snh0_294, sng_194, \
                         sng_209, sng_210, snh1_294, sof0_189, sof1_189, sog_284, \
                         sog_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_10 * sng_209[k]
                   + f_3 * pc_y[k] * sog_284[k];

        t_398[k] = f_11 * sng_194[k]
                   + f_1 * sof0_189[k]
                   - f_2 * sof1_189[k]
                   + f_3 * pc_z[k] * sog_284[k];

        t_399[k] = pb_y[k] * snh0_294[k]
                   - f_8 * pc_y[k] * snh1_294[k];

        t_400[k] = f_9 * sng_210[k]
                   + f_3 * pc_y[k] * sog_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pb_y, pc_y, pc_z, snh0_297, snh0_299, \
                         sng_195, sng_211, sng_212, snh1_297, snh1_299, sog_285, \
                         sog_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_16 * sng_195[k]
                   + f_3 * pc_z[k] * sog_285[k];

        t_402[k] = pb_y[k] * snh0_297[k]
                   + f_10 * sng_211[k]
                   - f_8 * pc_y[k] * snh1_297[k];

        t_403[k] = f_9 * sng_212[k]
                   + f_3 * pc_y[k] * sog_287[k];

        t_404[k] = pb_y[k] * snh0_299[k]
                   - f_8 * pc_y[k] * snh1_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pb_y, pc_y, pc_z, snh0_300, snh0_303, \
                         sng_198, sng_213, sng_215, snh1_300, snh1_303, sog_288, \
                         sog_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pb_y[k] * snh0_300[k]
                   + f_11 * sng_213[k]
                   - f_8 * pc_y[k] * snh1_300[k];

        t_406[k] = f_16 * sng_198[k]
                   + f_3 * pc_z[k] * sog_288[k];

        t_407[k] = f_9 * sng_215[k]
                   + f_3 * pc_y[k] * sog_290[k];

        t_408[k] = pb_y[k] * snh0_303[k]
                   - f_8 * pc_y[k] * snh1_303[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, sng_295, sng_296, sng_297, \
                         sng_298, sng_299, sog_295, sog_296, sog_297, sog_298, \
                         sog_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_17 * sng_295[k]
                   + f_3 * pc_x[k] * sog_295[k];

        t_410[k] = f_17 * sng_296[k]
                   + f_3 * pc_x[k] * sog_296[k];

        t_411[k] = f_17 * sng_297[k]
                   + f_3 * pc_x[k] * sog_297[k];

        t_412[k] = f_17 * sng_298[k]
                   + f_3 * pc_x[k] * sog_298[k];

        t_413[k] = f_17 * sng_299[k]
                   + f_3 * pc_x[k] * sog_299[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, pc_z, sng_205, sng_220, sng_222, sof0_196, \
                         sof0_198, sof1_196, sof1_198, sog_295, \
                         sog_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_9 * sng_220[k]
                   + f_1 * sof0_196[k]
                   - f_2 * sof1_196[k]
                   + f_3 * pc_y[k] * sog_295[k];

        t_415[k] = f_16 * sng_205[k]
                   + f_3 * pc_z[k] * sog_295[k];

        t_416[k] = f_9 * sng_222[k]
                   + f_4 * sof0_198[k]
                   - f_5 * sof1_198[k]
                   + f_3 * pc_y[k] * sog_297[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pb_y, pc_y, snh0_314, sng_223, sng_224, \
                         snh1_314, sof0_199, sof1_199, sog_298, \
                         sog_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_9 * sng_223[k]
                   + f_6 * sof0_199[k]
                   - f_7 * sof1_199[k]
                   + f_3 * pc_y[k] * sog_298[k];

        t_418[k] = f_9 * sng_224[k]
                   + f_3 * pc_y[k] * sog_299[k];

        t_419[k] = pb_y[k] * snh0_314[k]
                   - f_8 * pc_y[k] * snh1_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, sng_210, sng_300, \
                         sng_303, sof0_200, sof0_203, sof1_200, sof1_203, sog_300, \
                         sog_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_17 * sng_300[k]
                   + f_1 * sof0_200[k]
                   - f_2 * sof1_200[k]
                   + f_3 * pc_x[k] * sog_300[k];

        t_421[k] = f_3 * pc_y[k] * sog_300[k];

        t_422[k] = f_18 * sng_210[k]
                   + f_3 * pc_z[k] * sog_300[k];

        t_423[k] = f_17 * sng_303[k]
                   + f_4 * sof0_203[k]
                   - f_5 * sof1_203[k]
                   + f_3 * pc_x[k] * sog_303[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pc_x, pc_y, sng_305, sng_306, sof0_205, \
                         sof0_206, sof1_205, sof1_206, sog_302, sog_305, \
                         sog_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_y[k] * sog_302[k];

        t_425[k] = f_17 * sng_305[k]
                   + f_4 * sof0_205[k]
                   - f_5 * sof1_205[k]
                   + f_3 * pc_x[k] * sog_305[k];

        t_426[k] = f_17 * sng_306[k]
                   + f_6 * sof0_206[k]
                   - f_7 * sof1_206[k]
                   + f_3 * pc_x[k] * sog_306[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, pc_x, pc_y, pc_z, sng_213, sng_309, \
                         sng_310, sof0_209, sof1_209, sog_303, sog_305, sog_309, \
                         sog_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_18 * sng_213[k]
                   + f_3 * pc_z[k] * sog_303[k];

        t_428[k] = f_3 * pc_y[k] * sog_305[k];

        t_429[k] = f_17 * sng_309[k]
                   + f_6 * sof0_209[k]
                   - f_7 * sof1_209[k]
                   + f_3 * pc_x[k] * sog_309[k];

        t_430[k] = f_17 * sng_310[k]
                   + f_3 * pc_x[k] * sog_310[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pc_x, sng_311, sng_312, sng_313, sng_314, \
                         sog_311, sog_312, sog_313, sog_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_17 * sng_311[k]
                   + f_3 * pc_x[k] * sog_311[k];

        t_432[k] = f_17 * sng_312[k]
                   + f_3 * pc_x[k] * sog_312[k];

        t_433[k] = f_17 * sng_313[k]
                   + f_3 * pc_x[k] * sog_313[k];

        t_434[k] = f_17 * sng_314[k]
                   + f_3 * pc_x[k] * sog_314[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, pc_y, pc_z, sng_220, sof0_206, sof0_208, \
                         sof0_209, sof1_206, sof1_208, sof1_209, sog_310, sog_312, \
                         sog_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_1 * sof0_206[k]
                   - f_2 * sof1_206[k]
                   + f_3 * pc_y[k] * sog_310[k];

        t_436[k] = f_18 * sng_220[k]
                   + f_3 * pc_z[k] * sog_310[k];

        t_437[k] = f_4 * sof0_208[k]
                   - f_5 * sof1_208[k]
                   + f_3 * pc_y[k] * sog_312[k];

        t_438[k] = f_6 * sof0_209[k]
                   - f_7 * sof1_209[k]
                   + f_3 * pc_y[k] * sog_313[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, pc_x, pc_y, pc_z, sng_224, sng_225, \
                         sng_315, sof0_209, sof0_210, sof1_209, sof1_210, sog_314, \
                         sog_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_3 * pc_y[k] * sog_314[k];

        t_440[k] = f_18 * sng_224[k]
                   + f_1 * sof0_209[k]
                   - f_2 * sof1_209[k]
                   + f_3 * pc_z[k] * sog_314[k];

        t_441[k] = f_18 * sng_315[k]
                   + f_1 * sof0_210[k]
                   - f_2 * sof1_210[k]
                   + f_3 * pc_x[k] * sog_315[k];

        t_442[k] = f_17 * sng_225[k]
                   + f_3 * pc_y[k] * sog_315[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, sng_227, sng_318, sof0_213, \
                         sof1_213, sog_315, sog_317, sog_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_3 * pc_z[k] * sog_315[k];

        t_444[k] = f_18 * sng_318[k]
                   + f_4 * sof0_213[k]
                   - f_5 * sof1_213[k]
                   + f_3 * pc_x[k] * sog_318[k];

        t_445[k] = f_17 * sng_227[k]
                   + f_3 * pc_y[k] * sog_317[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, sng_320, sng_321, sof0_215, \
                         sof0_216, sof1_215, sof1_216, sog_318, sog_320, \
                         sog_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_18 * sng_320[k]
                   + f_4 * sof0_215[k]
                   - f_5 * sof1_215[k]
                   + f_3 * pc_x[k] * sog_320[k];

        t_447[k] = f_18 * sng_321[k]
                   + f_6 * sof0_216[k]
                   - f_7 * sof1_216[k]
                   + f_3 * pc_x[k] * sog_321[k];

        t_448[k] = f_3 * pc_z[k] * sog_318[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, sng_230, sng_324, sng_325, \
                         sng_326, sof0_219, sof1_219, sog_320, sog_324, sog_325, \
                         sog_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_17 * sng_230[k]
                   + f_3 * pc_y[k] * sog_320[k];

        t_450[k] = f_18 * sng_324[k]
                   + f_6 * sof0_219[k]
                   - f_7 * sof1_219[k]
                   + f_3 * pc_x[k] * sog_324[k];

        t_451[k] = f_18 * sng_325[k]
                   + f_3 * pc_x[k] * sog_325[k];

        t_452[k] = f_18 * sng_326[k]
                   + f_3 * pc_x[k] * sog_326[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pc_x, pc_y, sng_235, sng_327, sng_328, \
                         sng_329, sof0_216, sof1_216, sog_325, sog_327, sog_328, \
                         sog_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_18 * sng_327[k]
                   + f_3 * pc_x[k] * sog_327[k];

        t_454[k] = f_18 * sng_328[k]
                   + f_3 * pc_x[k] * sog_328[k];

        t_455[k] = f_18 * sng_329[k]
                   + f_3 * pc_x[k] * sog_329[k];

        t_456[k] = f_17 * sng_235[k]
                   + f_1 * sof0_216[k]
                   - f_2 * sof1_216[k]
                   + f_3 * pc_y[k] * sog_325[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, pc_y, pc_z, sng_237, sng_238, sof0_218, \
                         sof0_219, sof1_218, sof1_219, sog_325, sog_327, \
                         sog_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * sog_325[k];

        t_458[k] = f_17 * sng_237[k]
                   + f_4 * sof0_218[k]
                   - f_5 * sof1_218[k]
                   + f_3 * pc_y[k] * sog_327[k];

        t_459[k] = f_17 * sng_238[k]
                   + f_6 * sof0_219[k]
                   - f_7 * sof1_219[k]
                   + f_3 * pc_y[k] * sog_328[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, pb_z, pc_y, pc_z, snh0_315, sng_239, \
                         sng_240, snh1_315, sof0_219, sof1_219, sog_329, \
                         sog_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_17 * sng_239[k]
                   + f_3 * pc_y[k] * sog_329[k];

        t_461[k] = f_1 * sof0_219[k]
                   - f_2 * sof1_219[k]
                   + f_3 * pc_z[k] * sog_329[k];

        t_462[k] = pb_z[k] * snh0_315[k]
                   - f_8 * pc_z[k] * snh1_315[k];

        t_463[k] = f_18 * sng_240[k]
                   + f_3 * pc_y[k] * sog_330[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, pb_z, pc_y, pc_z, snh0_318, sng_225, sng_242, \
                         snh1_318, sog_330, sog_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_9 * sng_225[k]
                   + f_3 * pc_z[k] * sog_330[k];

        t_465[k] = pb_z[k] * snh0_318[k]
                   - f_8 * pc_z[k] * snh1_318[k];

        t_466[k] = f_18 * sng_242[k]
                   + f_3 * pc_y[k] * sog_332[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, pb_z, pc_x, pc_z, snh0_321, sng_228, sng_335, \
                         snh1_321, sof0_225, sof1_225, sog_333, \
                         sog_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_18 * sng_335[k]
                   + f_4 * sof0_225[k]
                   - f_5 * sof1_225[k]
                   + f_3 * pc_x[k] * sog_335[k];

        t_468[k] = pb_z[k] * snh0_321[k]
                   - f_8 * pc_z[k] * snh1_321[k];

        t_469[k] = f_9 * sng_228[k]
                   + f_3 * pc_z[k] * sog_333[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pc_x, pc_y, sng_245, sng_339, sng_340, \
                         sng_341, sof0_229, sof1_229, sog_335, sog_339, sog_340, \
                         sog_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_18 * sng_245[k]
                   + f_3 * pc_y[k] * sog_335[k];

        t_471[k] = f_18 * sng_339[k]
                   + f_6 * sof0_229[k]
                   - f_7 * sof1_229[k]
                   + f_3 * pc_x[k] * sog_339[k];

        t_472[k] = f_18 * sng_340[k]
                   + f_3 * pc_x[k] * sog_340[k];

        t_473[k] = f_18 * sng_341[k]
                   + f_3 * pc_x[k] * sog_341[k];
    }
}

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snh0,
                                                          const size_t sng, const size_t snh1,
                                                          const size_t sof0, const size_t sof1,
                                                          const size_t sog, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_330 = buffer.data(snh0 + 330);
    const auto *snh0_420 = buffer.data(snh0 + 420);
    const auto *snh0_423 = buffer.data(snh0 + 423);
    const auto *snh0_425 = buffer.data(snh0 + 425);
    const auto *snh0_426 = buffer.data(snh0 + 426);
    const auto *snh0_429 = buffer.data(snh0 + 429);
    const auto *snh0_440 = buffer.data(snh0 + 440);

    const auto *sng_235 = buffer.data(sng + 235);
    const auto *sng_239 = buffer.data(sng + 239);
    const auto *sng_240 = buffer.data(sng + 240);
    const auto *sng_243 = buffer.data(sng + 243);
    const auto *sng_250 = buffer.data(sng + 250);
    const auto *sng_252 = buffer.data(sng + 252);
    const auto *sng_253 = buffer.data(sng + 253);
    const auto *sng_254 = buffer.data(sng + 254);
    const auto *sng_255 = buffer.data(sng + 255);
    const auto *sng_257 = buffer.data(sng + 257);
    const auto *sng_258 = buffer.data(sng + 258);
    const auto *sng_260 = buffer.data(sng + 260);
    const auto *sng_265 = buffer.data(sng + 265);
    const auto *sng_267 = buffer.data(sng + 267);
    const auto *sng_268 = buffer.data(sng + 268);
    const auto *sng_269 = buffer.data(sng + 269);
    const auto *sng_270 = buffer.data(sng + 270);
    const auto *sng_272 = buffer.data(sng + 272);
    const auto *sng_273 = buffer.data(sng + 273);
    const auto *sng_275 = buffer.data(sng + 275);
    const auto *sng_280 = buffer.data(sng + 280);
    const auto *sng_282 = buffer.data(sng + 282);
    const auto *sng_283 = buffer.data(sng + 283);
    const auto *sng_284 = buffer.data(sng + 284);
    const auto *sng_285 = buffer.data(sng + 285);
    const auto *sng_287 = buffer.data(sng + 287);
    const auto *sng_288 = buffer.data(sng + 288);
    const auto *sng_290 = buffer.data(sng + 290);
    const auto *sng_295 = buffer.data(sng + 295);
    const auto *sng_297 = buffer.data(sng + 297);
    const auto *sng_298 = buffer.data(sng + 298);
    const auto *sng_299 = buffer.data(sng + 299);
    const auto *sng_300 = buffer.data(sng + 300);
    const auto *sng_301 = buffer.data(sng + 301);
    const auto *sng_302 = buffer.data(sng + 302);
    const auto *sng_303 = buffer.data(sng + 303);
    const auto *sng_305 = buffer.data(sng + 305);
    const auto *sng_310 = buffer.data(sng + 310);
    const auto *sng_312 = buffer.data(sng + 312);
    const auto *sng_313 = buffer.data(sng + 313);
    const auto *sng_314 = buffer.data(sng + 314);
    const auto *sng_342 = buffer.data(sng + 342);
    const auto *sng_343 = buffer.data(sng + 343);
    const auto *sng_344 = buffer.data(sng + 344);
    const auto *sng_345 = buffer.data(sng + 345);
    const auto *sng_348 = buffer.data(sng + 348);
    const auto *sng_350 = buffer.data(sng + 350);
    const auto *sng_351 = buffer.data(sng + 351);
    const auto *sng_354 = buffer.data(sng + 354);
    const auto *sng_355 = buffer.data(sng + 355);
    const auto *sng_356 = buffer.data(sng + 356);
    const auto *sng_357 = buffer.data(sng + 357);
    const auto *sng_358 = buffer.data(sng + 358);
    const auto *sng_359 = buffer.data(sng + 359);
    const auto *sng_360 = buffer.data(sng + 360);
    const auto *sng_363 = buffer.data(sng + 363);
    const auto *sng_365 = buffer.data(sng + 365);
    const auto *sng_366 = buffer.data(sng + 366);
    const auto *sng_369 = buffer.data(sng + 369);
    const auto *sng_370 = buffer.data(sng + 370);
    const auto *sng_371 = buffer.data(sng + 371);
    const auto *sng_372 = buffer.data(sng + 372);
    const auto *sng_373 = buffer.data(sng + 373);
    const auto *sng_374 = buffer.data(sng + 374);
    const auto *sng_375 = buffer.data(sng + 375);
    const auto *sng_378 = buffer.data(sng + 378);
    const auto *sng_380 = buffer.data(sng + 380);
    const auto *sng_381 = buffer.data(sng + 381);
    const auto *sng_384 = buffer.data(sng + 384);
    const auto *sng_385 = buffer.data(sng + 385);
    const auto *sng_386 = buffer.data(sng + 386);
    const auto *sng_387 = buffer.data(sng + 387);
    const auto *sng_388 = buffer.data(sng + 388);
    const auto *sng_389 = buffer.data(sng + 389);
    const auto *sng_400 = buffer.data(sng + 400);
    const auto *sng_401 = buffer.data(sng + 401);
    const auto *sng_402 = buffer.data(sng + 402);
    const auto *sng_403 = buffer.data(sng + 403);
    const auto *sng_404 = buffer.data(sng + 404);
    const auto *sng_405 = buffer.data(sng + 405);
    const auto *sng_408 = buffer.data(sng + 408);
    const auto *sng_410 = buffer.data(sng + 410);
    const auto *sng_411 = buffer.data(sng + 411);
    const auto *sng_414 = buffer.data(sng + 414);
    const auto *sng_415 = buffer.data(sng + 415);
    const auto *sng_416 = buffer.data(sng + 416);
    const auto *sng_417 = buffer.data(sng + 417);
    const auto *sng_418 = buffer.data(sng + 418);
    const auto *sng_419 = buffer.data(sng + 419);

    const auto *snh1_330 = buffer.data(snh1 + 330);
    const auto *snh1_420 = buffer.data(snh1 + 420);
    const auto *snh1_423 = buffer.data(snh1 + 423);
    const auto *snh1_425 = buffer.data(snh1 + 425);
    const auto *snh1_426 = buffer.data(snh1 + 426);
    const auto *snh1_429 = buffer.data(snh1 + 429);
    const auto *snh1_440 = buffer.data(snh1 + 440);

    const auto *sof0_228 = buffer.data(sof0 + 228);
    const auto *sof0_229 = buffer.data(sof0 + 229);
    const auto *sof0_230 = buffer.data(sof0 + 230);
    const auto *sof0_233 = buffer.data(sof0 + 233);
    const auto *sof0_235 = buffer.data(sof0 + 235);
    const auto *sof0_236 = buffer.data(sof0 + 236);
    const auto *sof0_238 = buffer.data(sof0 + 238);
    const auto *sof0_239 = buffer.data(sof0 + 239);
    const auto *sof0_240 = buffer.data(sof0 + 240);
    const auto *sof0_243 = buffer.data(sof0 + 243);
    const auto *sof0_245 = buffer.data(sof0 + 245);
    const auto *sof0_246 = buffer.data(sof0 + 246);
    const auto *sof0_248 = buffer.data(sof0 + 248);
    const auto *sof0_249 = buffer.data(sof0 + 249);
    const auto *sof0_250 = buffer.data(sof0 + 250);
    const auto *sof0_253 = buffer.data(sof0 + 253);
    const auto *sof0_255 = buffer.data(sof0 + 255);
    const auto *sof0_256 = buffer.data(sof0 + 256);
    const auto *sof0_258 = buffer.data(sof0 + 258);
    const auto *sof0_259 = buffer.data(sof0 + 259);
    const auto *sof0_266 = buffer.data(sof0 + 266);
    const auto *sof0_268 = buffer.data(sof0 + 268);
    const auto *sof0_269 = buffer.data(sof0 + 269);
    const auto *sof0_270 = buffer.data(sof0 + 270);
    const auto *sof0_273 = buffer.data(sof0 + 273);
    const auto *sof0_275 = buffer.data(sof0 + 275);
    const auto *sof0_276 = buffer.data(sof0 + 276);
    const auto *sof0_278 = buffer.data(sof0 + 278);
    const auto *sof0_279 = buffer.data(sof0 + 279);

    const auto *sof1_228 = buffer.data(sof1 + 228);
    const auto *sof1_229 = buffer.data(sof1 + 229);
    const auto *sof1_230 = buffer.data(sof1 + 230);
    const auto *sof1_233 = buffer.data(sof1 + 233);
    const auto *sof1_235 = buffer.data(sof1 + 235);
    const auto *sof1_236 = buffer.data(sof1 + 236);
    const auto *sof1_238 = buffer.data(sof1 + 238);
    const auto *sof1_239 = buffer.data(sof1 + 239);
    const auto *sof1_240 = buffer.data(sof1 + 240);
    const auto *sof1_243 = buffer.data(sof1 + 243);
    const auto *sof1_245 = buffer.data(sof1 + 245);
    const auto *sof1_246 = buffer.data(sof1 + 246);
    const auto *sof1_248 = buffer.data(sof1 + 248);
    const auto *sof1_249 = buffer.data(sof1 + 249);
    const auto *sof1_250 = buffer.data(sof1 + 250);
    const auto *sof1_253 = buffer.data(sof1 + 253);
    const auto *sof1_255 = buffer.data(sof1 + 255);
    const auto *sof1_256 = buffer.data(sof1 + 256);
    const auto *sof1_258 = buffer.data(sof1 + 258);
    const auto *sof1_259 = buffer.data(sof1 + 259);
    const auto *sof1_266 = buffer.data(sof1 + 266);
    const auto *sof1_268 = buffer.data(sof1 + 268);
    const auto *sof1_269 = buffer.data(sof1 + 269);
    const auto *sof1_270 = buffer.data(sof1 + 270);
    const auto *sof1_273 = buffer.data(sof1 + 273);
    const auto *sof1_275 = buffer.data(sof1 + 275);
    const auto *sof1_276 = buffer.data(sof1 + 276);
    const auto *sof1_278 = buffer.data(sof1 + 278);
    const auto *sof1_279 = buffer.data(sof1 + 279);

    const auto *sog_340 = buffer.data(sog + 340);
    const auto *sog_342 = buffer.data(sog + 342);
    const auto *sog_343 = buffer.data(sog + 343);
    const auto *sog_344 = buffer.data(sog + 344);
    const auto *sog_345 = buffer.data(sog + 345);
    const auto *sog_347 = buffer.data(sog + 347);
    const auto *sog_348 = buffer.data(sog + 348);
    const auto *sog_350 = buffer.data(sog + 350);
    const auto *sog_351 = buffer.data(sog + 351);
    const auto *sog_354 = buffer.data(sog + 354);
    const auto *sog_355 = buffer.data(sog + 355);
    const auto *sog_356 = buffer.data(sog + 356);
    const auto *sog_357 = buffer.data(sog + 357);
    const auto *sog_358 = buffer.data(sog + 358);
    const auto *sog_359 = buffer.data(sog + 359);
    const auto *sog_360 = buffer.data(sog + 360);
    const auto *sog_362 = buffer.data(sog + 362);
    const auto *sog_363 = buffer.data(sog + 363);
    const auto *sog_365 = buffer.data(sog + 365);
    const auto *sog_366 = buffer.data(sog + 366);
    const auto *sog_369 = buffer.data(sog + 369);
    const auto *sog_370 = buffer.data(sog + 370);
    const auto *sog_371 = buffer.data(sog + 371);
    const auto *sog_372 = buffer.data(sog + 372);
    const auto *sog_373 = buffer.data(sog + 373);
    const auto *sog_374 = buffer.data(sog + 374);
    const auto *sog_375 = buffer.data(sog + 375);
    const auto *sog_377 = buffer.data(sog + 377);
    const auto *sog_378 = buffer.data(sog + 378);
    const auto *sog_380 = buffer.data(sog + 380);
    const auto *sog_381 = buffer.data(sog + 381);
    const auto *sog_384 = buffer.data(sog + 384);
    const auto *sog_385 = buffer.data(sog + 385);
    const auto *sog_386 = buffer.data(sog + 386);
    const auto *sog_387 = buffer.data(sog + 387);
    const auto *sog_388 = buffer.data(sog + 388);
    const auto *sog_389 = buffer.data(sog + 389);
    const auto *sog_390 = buffer.data(sog + 390);
    const auto *sog_392 = buffer.data(sog + 392);
    const auto *sog_393 = buffer.data(sog + 393);
    const auto *sog_395 = buffer.data(sog + 395);
    const auto *sog_400 = buffer.data(sog + 400);
    const auto *sog_401 = buffer.data(sog + 401);
    const auto *sog_402 = buffer.data(sog + 402);
    const auto *sog_403 = buffer.data(sog + 403);
    const auto *sog_404 = buffer.data(sog + 404);
    const auto *sog_405 = buffer.data(sog + 405);
    const auto *sog_407 = buffer.data(sog + 407);
    const auto *sog_408 = buffer.data(sog + 408);
    const auto *sog_410 = buffer.data(sog + 410);
    const auto *sog_411 = buffer.data(sog + 411);
    const auto *sog_414 = buffer.data(sog + 414);
    const auto *sog_415 = buffer.data(sog + 415);
    const auto *sog_416 = buffer.data(sog + 416);
    const auto *sog_417 = buffer.data(sog + 417);
    const auto *sog_418 = buffer.data(sog + 418);
    const auto *sog_419 = buffer.data(sog + 419);

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_z, pc_x, pc_z, snh0_330, sng_342, \
                         sng_343, sng_344, snh1_330, sog_342, sog_343, \
                         sog_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_18 * sng_342[k]
                   + f_3 * pc_x[k] * sog_342[k];

        t_475[k] = f_18 * sng_343[k]
                   + f_3 * pc_x[k] * sog_343[k];

        t_476[k] = f_18 * sng_344[k]
                   + f_3 * pc_x[k] * sog_344[k];

        t_477[k] = pb_z[k] * snh0_330[k]
                   - f_8 * pc_z[k] * snh1_330[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pc_y, pc_z, sng_235, sng_252, sng_253, sof0_228, \
                         sof0_229, sof1_228, sof1_229, sog_340, sog_342, \
                         sog_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_9 * sng_235[k]
                   + f_3 * pc_z[k] * sog_340[k];

        t_479[k] = f_18 * sng_252[k]
                   + f_4 * sof0_228[k]
                   - f_5 * sof1_228[k]
                   + f_3 * pc_y[k] * sog_342[k];

        t_480[k] = f_18 * sng_253[k]
                   + f_6 * sof0_229[k]
                   - f_7 * sof1_229[k]
                   + f_3 * pc_y[k] * sog_343[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, pc_x, pc_y, pc_z, sng_239, sng_254, sng_345, \
                         sof0_229, sof0_230, sof1_229, sof1_230, sog_344, \
                         sog_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_18 * sng_254[k]
                   + f_3 * pc_y[k] * sog_344[k];

        t_482[k] = f_9 * sng_239[k]
                   + f_1 * sof0_229[k]
                   - f_2 * sof1_229[k]
                   + f_3 * pc_z[k] * sog_344[k];

        t_483[k] = f_18 * sng_345[k]
                   + f_1 * sof0_230[k]
                   - f_2 * sof1_230[k]
                   + f_3 * pc_x[k] * sog_345[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, pc_x, pc_y, pc_z, sng_240, sng_255, \
                         sng_257, sng_348, sof0_233, sof1_233, sog_345, sog_347, \
                         sog_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_16 * sng_255[k]
                   + f_3 * pc_y[k] * sog_345[k];

        t_485[k] = f_10 * sng_240[k]
                   + f_3 * pc_z[k] * sog_345[k];

        t_486[k] = f_18 * sng_348[k]
                   + f_4 * sof0_233[k]
                   - f_5 * sof1_233[k]
                   + f_3 * pc_x[k] * sog_348[k];

        t_487[k] = f_16 * sng_257[k]
                   + f_3 * pc_y[k] * sog_347[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_z, sng_243, sng_350, sng_351, sof0_235, \
                         sof0_236, sof1_235, sof1_236, sog_348, sog_350, \
                         sog_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_18 * sng_350[k]
                   + f_4 * sof0_235[k]
                   - f_5 * sof1_235[k]
                   + f_3 * pc_x[k] * sog_350[k];

        t_489[k] = f_18 * sng_351[k]
                   + f_6 * sof0_236[k]
                   - f_7 * sof1_236[k]
                   + f_3 * pc_x[k] * sog_351[k];

        t_490[k] = f_10 * sng_243[k]
                   + f_3 * pc_z[k] * sog_348[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_x, pc_y, sng_260, sng_354, sng_355, \
                         sng_356, sof0_239, sof1_239, sog_350, sog_354, sog_355, \
                         sog_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_16 * sng_260[k]
                   + f_3 * pc_y[k] * sog_350[k];

        t_492[k] = f_18 * sng_354[k]
                   + f_6 * sof0_239[k]
                   - f_7 * sof1_239[k]
                   + f_3 * pc_x[k] * sog_354[k];

        t_493[k] = f_18 * sng_355[k]
                   + f_3 * pc_x[k] * sog_355[k];

        t_494[k] = f_18 * sng_356[k]
                   + f_3 * pc_x[k] * sog_356[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pc_x, pc_y, sng_265, sng_357, sng_358, \
                         sng_359, sof0_236, sof1_236, sog_355, sog_357, sog_358, \
                         sog_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_18 * sng_357[k]
                   + f_3 * pc_x[k] * sog_357[k];

        t_496[k] = f_18 * sng_358[k]
                   + f_3 * pc_x[k] * sog_358[k];

        t_497[k] = f_18 * sng_359[k]
                   + f_3 * pc_x[k] * sog_359[k];

        t_498[k] = f_16 * sng_265[k]
                   + f_1 * sof0_236[k]
                   - f_2 * sof1_236[k]
                   + f_3 * pc_y[k] * sog_355[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, pc_z, sng_250, sng_267, sng_268, sof0_238, \
                         sof0_239, sof1_238, sof1_239, sog_355, sog_357, \
                         sog_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_10 * sng_250[k]
                   + f_3 * pc_z[k] * sog_355[k];

        t_500[k] = f_16 * sng_267[k]
                   + f_4 * sof0_238[k]
                   - f_5 * sof1_238[k]
                   + f_3 * pc_y[k] * sog_357[k];

        t_501[k] = f_16 * sng_268[k]
                   + f_6 * sof0_239[k]
                   - f_7 * sof1_239[k]
                   + f_3 * pc_y[k] * sog_358[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, sng_254, sng_269, sng_360, \
                         sof0_239, sof0_240, sof1_239, sof1_240, sog_359, \
                         sog_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_16 * sng_269[k]
                   + f_3 * pc_y[k] * sog_359[k];

        t_503[k] = f_10 * sng_254[k]
                   + f_1 * sof0_239[k]
                   - f_2 * sof1_239[k]
                   + f_3 * pc_z[k] * sog_359[k];

        t_504[k] = f_18 * sng_360[k]
                   + f_1 * sof0_240[k]
                   - f_2 * sof1_240[k]
                   + f_3 * pc_x[k] * sog_360[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, sng_255, sng_270, \
                         sng_272, sng_363, sof0_243, sof1_243, sog_360, sog_362, \
                         sog_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_11 * sng_270[k]
                   + f_3 * pc_y[k] * sog_360[k];

        t_506[k] = f_11 * sng_255[k]
                   + f_3 * pc_z[k] * sog_360[k];

        t_507[k] = f_18 * sng_363[k]
                   + f_4 * sof0_243[k]
                   - f_5 * sof1_243[k]
                   + f_3 * pc_x[k] * sog_363[k];

        t_508[k] = f_11 * sng_272[k]
                   + f_3 * pc_y[k] * sog_362[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, sng_258, sng_365, sng_366, sof0_245, \
                         sof0_246, sof1_245, sof1_246, sog_363, sog_365, \
                         sog_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_18 * sng_365[k]
                   + f_4 * sof0_245[k]
                   - f_5 * sof1_245[k]
                   + f_3 * pc_x[k] * sog_365[k];

        t_510[k] = f_18 * sng_366[k]
                   + f_6 * sof0_246[k]
                   - f_7 * sof1_246[k]
                   + f_3 * pc_x[k] * sog_366[k];

        t_511[k] = f_11 * sng_258[k]
                   + f_3 * pc_z[k] * sog_363[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pc_x, pc_y, sng_275, sng_369, sng_370, \
                         sng_371, sof0_249, sof1_249, sog_365, sog_369, sog_370, \
                         sog_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_11 * sng_275[k]
                   + f_3 * pc_y[k] * sog_365[k];

        t_513[k] = f_18 * sng_369[k]
                   + f_6 * sof0_249[k]
                   - f_7 * sof1_249[k]
                   + f_3 * pc_x[k] * sog_369[k];

        t_514[k] = f_18 * sng_370[k]
                   + f_3 * pc_x[k] * sog_370[k];

        t_515[k] = f_18 * sng_371[k]
                   + f_3 * pc_x[k] * sog_371[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, pc_y, sng_280, sng_372, sng_373, \
                         sng_374, sof0_246, sof1_246, sog_370, sog_372, sog_373, \
                         sog_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_18 * sng_372[k]
                   + f_3 * pc_x[k] * sog_372[k];

        t_517[k] = f_18 * sng_373[k]
                   + f_3 * pc_x[k] * sog_373[k];

        t_518[k] = f_18 * sng_374[k]
                   + f_3 * pc_x[k] * sog_374[k];

        t_519[k] = f_11 * sng_280[k]
                   + f_1 * sof0_246[k]
                   - f_2 * sof1_246[k]
                   + f_3 * pc_y[k] * sog_370[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, sng_265, sng_282, sng_283, sof0_248, \
                         sof0_249, sof1_248, sof1_249, sog_370, sog_372, \
                         sog_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_11 * sng_265[k]
                   + f_3 * pc_z[k] * sog_370[k];

        t_521[k] = f_11 * sng_282[k]
                   + f_4 * sof0_248[k]
                   - f_5 * sof1_248[k]
                   + f_3 * pc_y[k] * sog_372[k];

        t_522[k] = f_11 * sng_283[k]
                   + f_6 * sof0_249[k]
                   - f_7 * sof1_249[k]
                   + f_3 * pc_y[k] * sog_373[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pc_x, pc_y, pc_z, sng_269, sng_284, sng_375, \
                         sof0_249, sof0_250, sof1_249, sof1_250, sog_374, \
                         sog_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_11 * sng_284[k]
                   + f_3 * pc_y[k] * sog_374[k];

        t_524[k] = f_11 * sng_269[k]
                   + f_1 * sof0_249[k]
                   - f_2 * sof1_249[k]
                   + f_3 * pc_z[k] * sog_374[k];

        t_525[k] = f_18 * sng_375[k]
                   + f_1 * sof0_250[k]
                   - f_2 * sof1_250[k]
                   + f_3 * pc_x[k] * sog_375[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pc_x, pc_y, pc_z, sng_270, sng_285, \
                         sng_287, sng_378, sof0_253, sof1_253, sog_375, sog_377, \
                         sog_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_10 * sng_285[k]
                   + f_3 * pc_y[k] * sog_375[k];

        t_527[k] = f_16 * sng_270[k]
                   + f_3 * pc_z[k] * sog_375[k];

        t_528[k] = f_18 * sng_378[k]
                   + f_4 * sof0_253[k]
                   - f_5 * sof1_253[k]
                   + f_3 * pc_x[k] * sog_378[k];

        t_529[k] = f_10 * sng_287[k]
                   + f_3 * pc_y[k] * sog_377[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pc_x, pc_z, sng_273, sng_380, sng_381, sof0_255, \
                         sof0_256, sof1_255, sof1_256, sog_378, sog_380, \
                         sog_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_18 * sng_380[k]
                   + f_4 * sof0_255[k]
                   - f_5 * sof1_255[k]
                   + f_3 * pc_x[k] * sog_380[k];

        t_531[k] = f_18 * sng_381[k]
                   + f_6 * sof0_256[k]
                   - f_7 * sof1_256[k]
                   + f_3 * pc_x[k] * sog_381[k];

        t_532[k] = f_16 * sng_273[k]
                   + f_3 * pc_z[k] * sog_378[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_y, sng_290, sng_384, sng_385, \
                         sng_386, sof0_259, sof1_259, sog_380, sog_384, sog_385, \
                         sog_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_10 * sng_290[k]
                   + f_3 * pc_y[k] * sog_380[k];

        t_534[k] = f_18 * sng_384[k]
                   + f_6 * sof0_259[k]
                   - f_7 * sof1_259[k]
                   + f_3 * pc_x[k] * sog_384[k];

        t_535[k] = f_18 * sng_385[k]
                   + f_3 * pc_x[k] * sog_385[k];

        t_536[k] = f_18 * sng_386[k]
                   + f_3 * pc_x[k] * sog_386[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, sng_295, sng_387, sng_388, \
                         sng_389, sof0_256, sof1_256, sog_385, sog_387, sog_388, \
                         sog_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_18 * sng_387[k]
                   + f_3 * pc_x[k] * sog_387[k];

        t_538[k] = f_18 * sng_388[k]
                   + f_3 * pc_x[k] * sog_388[k];

        t_539[k] = f_18 * sng_389[k]
                   + f_3 * pc_x[k] * sog_389[k];

        t_540[k] = f_10 * sng_295[k]
                   + f_1 * sof0_256[k]
                   - f_2 * sof1_256[k]
                   + f_3 * pc_y[k] * sog_385[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pc_y, pc_z, sng_280, sng_297, sng_298, sof0_258, \
                         sof0_259, sof1_258, sof1_259, sog_385, sog_387, \
                         sog_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_16 * sng_280[k]
                   + f_3 * pc_z[k] * sog_385[k];

        t_542[k] = f_10 * sng_297[k]
                   + f_4 * sof0_258[k]
                   - f_5 * sof1_258[k]
                   + f_3 * pc_y[k] * sog_387[k];

        t_543[k] = f_10 * sng_298[k]
                   + f_6 * sof0_259[k]
                   - f_7 * sof1_259[k]
                   + f_3 * pc_y[k] * sog_388[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pb_y, pc_y, pc_z, snh0_420, sng_284, \
                         sng_299, sng_300, snh1_420, sof0_259, sof1_259, sog_389, \
                         sog_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_10 * sng_299[k]
                   + f_3 * pc_y[k] * sog_389[k];

        t_545[k] = f_16 * sng_284[k]
                   + f_1 * sof0_259[k]
                   - f_2 * sof1_259[k]
                   + f_3 * pc_z[k] * sog_389[k];

        t_546[k] = pb_y[k] * snh0_420[k]
                   - f_8 * pc_y[k] * snh1_420[k];

        t_547[k] = f_9 * sng_300[k]
                   + f_3 * pc_y[k] * sog_390[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pb_y, pc_y, pc_z, snh0_423, snh0_425, \
                         sng_285, sng_301, sng_302, snh1_423, snh1_425, sog_390, \
                         sog_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_18 * sng_285[k]
                   + f_3 * pc_z[k] * sog_390[k];

        t_549[k] = pb_y[k] * snh0_423[k]
                   + f_10 * sng_301[k]
                   - f_8 * pc_y[k] * snh1_423[k];

        t_550[k] = f_9 * sng_302[k]
                   + f_3 * pc_y[k] * sog_392[k];

        t_551[k] = pb_y[k] * snh0_425[k]
                   - f_8 * pc_y[k] * snh1_425[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pb_y, pc_y, pc_z, snh0_426, snh0_429, \
                         sng_288, sng_303, sng_305, snh1_426, snh1_429, sog_393, \
                         sog_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pb_y[k] * snh0_426[k]
                   + f_11 * sng_303[k]
                   - f_8 * pc_y[k] * snh1_426[k];

        t_553[k] = f_18 * sng_288[k]
                   + f_3 * pc_z[k] * sog_393[k];

        t_554[k] = f_9 * sng_305[k]
                   + f_3 * pc_y[k] * sog_395[k];

        t_555[k] = pb_y[k] * snh0_429[k]
                   - f_8 * pc_y[k] * snh1_429[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, pc_x, sng_400, sng_401, sng_402, \
                         sng_403, sng_404, sog_400, sog_401, sog_402, sog_403, \
                         sog_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_18 * sng_400[k]
                   + f_3 * pc_x[k] * sog_400[k];

        t_557[k] = f_18 * sng_401[k]
                   + f_3 * pc_x[k] * sog_401[k];

        t_558[k] = f_18 * sng_402[k]
                   + f_3 * pc_x[k] * sog_402[k];

        t_559[k] = f_18 * sng_403[k]
                   + f_3 * pc_x[k] * sog_403[k];

        t_560[k] = f_18 * sng_404[k]
                   + f_3 * pc_x[k] * sog_404[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pc_y, pc_z, sng_295, sng_310, sng_312, sof0_266, \
                         sof0_268, sof1_266, sof1_268, sog_400, \
                         sog_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_9 * sng_310[k]
                   + f_1 * sof0_266[k]
                   - f_2 * sof1_266[k]
                   + f_3 * pc_y[k] * sog_400[k];

        t_562[k] = f_18 * sng_295[k]
                   + f_3 * pc_z[k] * sog_400[k];

        t_563[k] = f_9 * sng_312[k]
                   + f_4 * sof0_268[k]
                   - f_5 * sof1_268[k]
                   + f_3 * pc_y[k] * sog_402[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pb_y, pc_y, snh0_440, sng_313, sng_314, \
                         snh1_440, sof0_269, sof1_269, sog_403, \
                         sog_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_9 * sng_313[k]
                   + f_6 * sof0_269[k]
                   - f_7 * sof1_269[k]
                   + f_3 * pc_y[k] * sog_403[k];

        t_565[k] = f_9 * sng_314[k]
                   + f_3 * pc_y[k] * sog_404[k];

        t_566[k] = pb_y[k] * snh0_440[k]
                   - f_8 * pc_y[k] * snh1_440[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pc_x, pc_y, pc_z, sng_300, sng_405, \
                         sng_408, sof0_270, sof0_273, sof1_270, sof1_273, sog_405, \
                         sog_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_18 * sng_405[k]
                   + f_1 * sof0_270[k]
                   - f_2 * sof1_270[k]
                   + f_3 * pc_x[k] * sog_405[k];

        t_568[k] = f_3 * pc_y[k] * sog_405[k];

        t_569[k] = f_17 * sng_300[k]
                   + f_3 * pc_z[k] * sog_405[k];

        t_570[k] = f_18 * sng_408[k]
                   + f_4 * sof0_273[k]
                   - f_5 * sof1_273[k]
                   + f_3 * pc_x[k] * sog_408[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pc_x, pc_y, sng_410, sng_411, sof0_275, \
                         sof0_276, sof1_275, sof1_276, sog_407, sog_410, \
                         sog_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_3 * pc_y[k] * sog_407[k];

        t_572[k] = f_18 * sng_410[k]
                   + f_4 * sof0_275[k]
                   - f_5 * sof1_275[k]
                   + f_3 * pc_x[k] * sog_410[k];

        t_573[k] = f_18 * sng_411[k]
                   + f_6 * sof0_276[k]
                   - f_7 * sof1_276[k]
                   + f_3 * pc_x[k] * sog_411[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pc_x, pc_y, pc_z, sng_303, sng_414, \
                         sng_415, sof0_279, sof1_279, sog_408, sog_410, sog_414, \
                         sog_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_17 * sng_303[k]
                   + f_3 * pc_z[k] * sog_408[k];

        t_575[k] = f_3 * pc_y[k] * sog_410[k];

        t_576[k] = f_18 * sng_414[k]
                   + f_6 * sof0_279[k]
                   - f_7 * sof1_279[k]
                   + f_3 * pc_x[k] * sog_414[k];

        t_577[k] = f_18 * sng_415[k]
                   + f_3 * pc_x[k] * sog_415[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, pc_x, sng_416, sng_417, sng_418, sng_419, \
                         sog_416, sog_417, sog_418, sog_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_18 * sng_416[k]
                   + f_3 * pc_x[k] * sog_416[k];

        t_579[k] = f_18 * sng_417[k]
                   + f_3 * pc_x[k] * sog_417[k];

        t_580[k] = f_18 * sng_418[k]
                   + f_3 * pc_x[k] * sog_418[k];

        t_581[k] = f_18 * sng_419[k]
                   + f_3 * pc_x[k] * sog_419[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pc_y, pc_z, sng_310, sof0_276, sof0_278, \
                         sof0_279, sof1_276, sof1_278, sof1_279, sog_415, sog_417, \
                         sog_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_1 * sof0_276[k]
                   - f_2 * sof1_276[k]
                   + f_3 * pc_y[k] * sog_415[k];

        t_583[k] = f_17 * sng_310[k]
                   + f_3 * pc_z[k] * sog_415[k];

        t_584[k] = f_4 * sof0_278[k]
                   - f_5 * sof1_278[k]
                   + f_3 * pc_y[k] * sog_417[k];

        t_585[k] = f_6 * sof0_279[k]
                   - f_7 * sof1_279[k]
                   + f_3 * pc_y[k] * sog_418[k];
    }
}

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snh0,
                                                          const size_t sng, const size_t snh1,
                                                          const size_t sof0, const size_t sof1,
                                                          const size_t sog, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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
    auto *t_693 = buffer.data(target + 693);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_441 = buffer.data(snh0 + 441);
    const auto *snh0_444 = buffer.data(snh0 + 444);
    const auto *snh0_447 = buffer.data(snh0 + 447);
    const auto *snh0_456 = buffer.data(snh0 + 456);

    const auto *sng_314 = buffer.data(sng + 314);
    const auto *sng_315 = buffer.data(sng + 315);
    const auto *sng_317 = buffer.data(sng + 317);
    const auto *sng_318 = buffer.data(sng + 318);
    const auto *sng_320 = buffer.data(sng + 320);
    const auto *sng_325 = buffer.data(sng + 325);
    const auto *sng_327 = buffer.data(sng + 327);
    const auto *sng_328 = buffer.data(sng + 328);
    const auto *sng_329 = buffer.data(sng + 329);
    const auto *sng_330 = buffer.data(sng + 330);
    const auto *sng_332 = buffer.data(sng + 332);
    const auto *sng_333 = buffer.data(sng + 333);
    const auto *sng_335 = buffer.data(sng + 335);
    const auto *sng_340 = buffer.data(sng + 340);
    const auto *sng_342 = buffer.data(sng + 342);
    const auto *sng_343 = buffer.data(sng + 343);
    const auto *sng_344 = buffer.data(sng + 344);
    const auto *sng_345 = buffer.data(sng + 345);
    const auto *sng_347 = buffer.data(sng + 347);
    const auto *sng_348 = buffer.data(sng + 348);
    const auto *sng_350 = buffer.data(sng + 350);
    const auto *sng_355 = buffer.data(sng + 355);
    const auto *sng_357 = buffer.data(sng + 357);
    const auto *sng_358 = buffer.data(sng + 358);
    const auto *sng_359 = buffer.data(sng + 359);
    const auto *sng_360 = buffer.data(sng + 360);
    const auto *sng_362 = buffer.data(sng + 362);
    const auto *sng_363 = buffer.data(sng + 363);
    const auto *sng_365 = buffer.data(sng + 365);
    const auto *sng_370 = buffer.data(sng + 370);
    const auto *sng_372 = buffer.data(sng + 372);
    const auto *sng_373 = buffer.data(sng + 373);
    const auto *sng_374 = buffer.data(sng + 374);
    const auto *sng_375 = buffer.data(sng + 375);
    const auto *sng_377 = buffer.data(sng + 377);
    const auto *sng_380 = buffer.data(sng + 380);
    const auto *sng_385 = buffer.data(sng + 385);
    const auto *sng_387 = buffer.data(sng + 387);
    const auto *sng_388 = buffer.data(sng + 388);
    const auto *sng_389 = buffer.data(sng + 389);
    const auto *sng_420 = buffer.data(sng + 420);
    const auto *sng_423 = buffer.data(sng + 423);
    const auto *sng_425 = buffer.data(sng + 425);
    const auto *sng_426 = buffer.data(sng + 426);
    const auto *sng_429 = buffer.data(sng + 429);
    const auto *sng_430 = buffer.data(sng + 430);
    const auto *sng_431 = buffer.data(sng + 431);
    const auto *sng_432 = buffer.data(sng + 432);
    const auto *sng_433 = buffer.data(sng + 433);
    const auto *sng_434 = buffer.data(sng + 434);
    const auto *sng_440 = buffer.data(sng + 440);
    const auto *sng_444 = buffer.data(sng + 444);
    const auto *sng_445 = buffer.data(sng + 445);
    const auto *sng_446 = buffer.data(sng + 446);
    const auto *sng_447 = buffer.data(sng + 447);
    const auto *sng_448 = buffer.data(sng + 448);
    const auto *sng_449 = buffer.data(sng + 449);
    const auto *sng_450 = buffer.data(sng + 450);
    const auto *sng_453 = buffer.data(sng + 453);
    const auto *sng_455 = buffer.data(sng + 455);
    const auto *sng_456 = buffer.data(sng + 456);
    const auto *sng_459 = buffer.data(sng + 459);
    const auto *sng_460 = buffer.data(sng + 460);
    const auto *sng_461 = buffer.data(sng + 461);
    const auto *sng_462 = buffer.data(sng + 462);
    const auto *sng_463 = buffer.data(sng + 463);
    const auto *sng_464 = buffer.data(sng + 464);
    const auto *sng_465 = buffer.data(sng + 465);
    const auto *sng_468 = buffer.data(sng + 468);
    const auto *sng_470 = buffer.data(sng + 470);
    const auto *sng_471 = buffer.data(sng + 471);
    const auto *sng_474 = buffer.data(sng + 474);
    const auto *sng_475 = buffer.data(sng + 475);
    const auto *sng_476 = buffer.data(sng + 476);
    const auto *sng_477 = buffer.data(sng + 477);
    const auto *sng_478 = buffer.data(sng + 478);
    const auto *sng_479 = buffer.data(sng + 479);
    const auto *sng_480 = buffer.data(sng + 480);
    const auto *sng_483 = buffer.data(sng + 483);
    const auto *sng_485 = buffer.data(sng + 485);
    const auto *sng_486 = buffer.data(sng + 486);
    const auto *sng_489 = buffer.data(sng + 489);
    const auto *sng_490 = buffer.data(sng + 490);
    const auto *sng_491 = buffer.data(sng + 491);
    const auto *sng_492 = buffer.data(sng + 492);
    const auto *sng_493 = buffer.data(sng + 493);
    const auto *sng_494 = buffer.data(sng + 494);
    const auto *sng_495 = buffer.data(sng + 495);

    const auto *snh1_441 = buffer.data(snh1 + 441);
    const auto *snh1_444 = buffer.data(snh1 + 444);
    const auto *snh1_447 = buffer.data(snh1 + 447);
    const auto *snh1_456 = buffer.data(snh1 + 456);

    const auto *sof0_279 = buffer.data(sof0 + 279);
    const auto *sof0_280 = buffer.data(sof0 + 280);
    const auto *sof0_283 = buffer.data(sof0 + 283);
    const auto *sof0_285 = buffer.data(sof0 + 285);
    const auto *sof0_286 = buffer.data(sof0 + 286);
    const auto *sof0_288 = buffer.data(sof0 + 288);
    const auto *sof0_289 = buffer.data(sof0 + 289);
    const auto *sof0_295 = buffer.data(sof0 + 295);
    const auto *sof0_298 = buffer.data(sof0 + 298);
    const auto *sof0_299 = buffer.data(sof0 + 299);
    const auto *sof0_300 = buffer.data(sof0 + 300);
    const auto *sof0_303 = buffer.data(sof0 + 303);
    const auto *sof0_305 = buffer.data(sof0 + 305);
    const auto *sof0_306 = buffer.data(sof0 + 306);
    const auto *sof0_308 = buffer.data(sof0 + 308);
    const auto *sof0_309 = buffer.data(sof0 + 309);
    const auto *sof0_310 = buffer.data(sof0 + 310);
    const auto *sof0_313 = buffer.data(sof0 + 313);
    const auto *sof0_315 = buffer.data(sof0 + 315);
    const auto *sof0_316 = buffer.data(sof0 + 316);
    const auto *sof0_318 = buffer.data(sof0 + 318);
    const auto *sof0_319 = buffer.data(sof0 + 319);
    const auto *sof0_320 = buffer.data(sof0 + 320);
    const auto *sof0_323 = buffer.data(sof0 + 323);
    const auto *sof0_325 = buffer.data(sof0 + 325);
    const auto *sof0_326 = buffer.data(sof0 + 326);
    const auto *sof0_328 = buffer.data(sof0 + 328);
    const auto *sof0_329 = buffer.data(sof0 + 329);
    const auto *sof0_330 = buffer.data(sof0 + 330);

    const auto *sof1_279 = buffer.data(sof1 + 279);
    const auto *sof1_280 = buffer.data(sof1 + 280);
    const auto *sof1_283 = buffer.data(sof1 + 283);
    const auto *sof1_285 = buffer.data(sof1 + 285);
    const auto *sof1_286 = buffer.data(sof1 + 286);
    const auto *sof1_288 = buffer.data(sof1 + 288);
    const auto *sof1_289 = buffer.data(sof1 + 289);
    const auto *sof1_295 = buffer.data(sof1 + 295);
    const auto *sof1_298 = buffer.data(sof1 + 298);
    const auto *sof1_299 = buffer.data(sof1 + 299);
    const auto *sof1_300 = buffer.data(sof1 + 300);
    const auto *sof1_303 = buffer.data(sof1 + 303);
    const auto *sof1_305 = buffer.data(sof1 + 305);
    const auto *sof1_306 = buffer.data(sof1 + 306);
    const auto *sof1_308 = buffer.data(sof1 + 308);
    const auto *sof1_309 = buffer.data(sof1 + 309);
    const auto *sof1_310 = buffer.data(sof1 + 310);
    const auto *sof1_313 = buffer.data(sof1 + 313);
    const auto *sof1_315 = buffer.data(sof1 + 315);
    const auto *sof1_316 = buffer.data(sof1 + 316);
    const auto *sof1_318 = buffer.data(sof1 + 318);
    const auto *sof1_319 = buffer.data(sof1 + 319);
    const auto *sof1_320 = buffer.data(sof1 + 320);
    const auto *sof1_323 = buffer.data(sof1 + 323);
    const auto *sof1_325 = buffer.data(sof1 + 325);
    const auto *sof1_326 = buffer.data(sof1 + 326);
    const auto *sof1_328 = buffer.data(sof1 + 328);
    const auto *sof1_329 = buffer.data(sof1 + 329);
    const auto *sof1_330 = buffer.data(sof1 + 330);

    const auto *sog_419 = buffer.data(sog + 419);
    const auto *sog_420 = buffer.data(sog + 420);
    const auto *sog_422 = buffer.data(sog + 422);
    const auto *sog_423 = buffer.data(sog + 423);
    const auto *sog_425 = buffer.data(sog + 425);
    const auto *sog_426 = buffer.data(sog + 426);
    const auto *sog_429 = buffer.data(sog + 429);
    const auto *sog_430 = buffer.data(sog + 430);
    const auto *sog_431 = buffer.data(sog + 431);
    const auto *sog_432 = buffer.data(sog + 432);
    const auto *sog_433 = buffer.data(sog + 433);
    const auto *sog_434 = buffer.data(sog + 434);
    const auto *sog_435 = buffer.data(sog + 435);
    const auto *sog_437 = buffer.data(sog + 437);
    const auto *sog_438 = buffer.data(sog + 438);
    const auto *sog_440 = buffer.data(sog + 440);
    const auto *sog_444 = buffer.data(sog + 444);
    const auto *sog_445 = buffer.data(sog + 445);
    const auto *sog_446 = buffer.data(sog + 446);
    const auto *sog_447 = buffer.data(sog + 447);
    const auto *sog_448 = buffer.data(sog + 448);
    const auto *sog_449 = buffer.data(sog + 449);
    const auto *sog_450 = buffer.data(sog + 450);
    const auto *sog_452 = buffer.data(sog + 452);
    const auto *sog_453 = buffer.data(sog + 453);
    const auto *sog_455 = buffer.data(sog + 455);
    const auto *sog_456 = buffer.data(sog + 456);
    const auto *sog_459 = buffer.data(sog + 459);
    const auto *sog_460 = buffer.data(sog + 460);
    const auto *sog_461 = buffer.data(sog + 461);
    const auto *sog_462 = buffer.data(sog + 462);
    const auto *sog_463 = buffer.data(sog + 463);
    const auto *sog_464 = buffer.data(sog + 464);
    const auto *sog_465 = buffer.data(sog + 465);
    const auto *sog_467 = buffer.data(sog + 467);
    const auto *sog_468 = buffer.data(sog + 468);
    const auto *sog_470 = buffer.data(sog + 470);
    const auto *sog_471 = buffer.data(sog + 471);
    const auto *sog_474 = buffer.data(sog + 474);
    const auto *sog_475 = buffer.data(sog + 475);
    const auto *sog_476 = buffer.data(sog + 476);
    const auto *sog_477 = buffer.data(sog + 477);
    const auto *sog_478 = buffer.data(sog + 478);
    const auto *sog_479 = buffer.data(sog + 479);
    const auto *sog_480 = buffer.data(sog + 480);
    const auto *sog_482 = buffer.data(sog + 482);
    const auto *sog_483 = buffer.data(sog + 483);
    const auto *sog_485 = buffer.data(sog + 485);
    const auto *sog_486 = buffer.data(sog + 486);
    const auto *sog_489 = buffer.data(sog + 489);
    const auto *sog_490 = buffer.data(sog + 490);
    const auto *sog_491 = buffer.data(sog + 491);
    const auto *sog_492 = buffer.data(sog + 492);
    const auto *sog_493 = buffer.data(sog + 493);
    const auto *sog_494 = buffer.data(sog + 494);
    const auto *sog_495 = buffer.data(sog + 495);

#pragma omp simd aligned(t_586, t_587, t_588, t_589, pc_x, pc_y, pc_z, sng_314, sng_315, \
                         sng_420, sof0_279, sof0_280, sof1_279, sof1_280, sog_419, \
                         sog_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_3 * pc_y[k] * sog_419[k];

        t_587[k] = f_17 * sng_314[k]
                   + f_1 * sof0_279[k]
                   - f_2 * sof1_279[k]
                   + f_3 * pc_z[k] * sog_419[k];

        t_588[k] = f_16 * sng_420[k]
                   + f_1 * sof0_280[k]
                   - f_2 * sof1_280[k]
                   + f_3 * pc_x[k] * sog_420[k];

        t_589[k] = f_15 * sng_315[k]
                   + f_3 * pc_y[k] * sog_420[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, pc_x, pc_y, pc_z, sng_317, sng_423, sof0_283, \
                         sof1_283, sog_420, sog_422, sog_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_3 * pc_z[k] * sog_420[k];

        t_591[k] = f_16 * sng_423[k]
                   + f_4 * sof0_283[k]
                   - f_5 * sof1_283[k]
                   + f_3 * pc_x[k] * sog_423[k];

        t_592[k] = f_15 * sng_317[k]
                   + f_3 * pc_y[k] * sog_422[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pc_x, pc_z, sng_425, sng_426, sof0_285, \
                         sof0_286, sof1_285, sof1_286, sog_423, sog_425, \
                         sog_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = f_16 * sng_425[k]
                   + f_4 * sof0_285[k]
                   - f_5 * sof1_285[k]
                   + f_3 * pc_x[k] * sog_425[k];

        t_594[k] = f_16 * sng_426[k]
                   + f_6 * sof0_286[k]
                   - f_7 * sof1_286[k]
                   + f_3 * pc_x[k] * sog_426[k];

        t_595[k] = f_3 * pc_z[k] * sog_423[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, pc_y, sng_320, sng_429, sng_430, \
                         sng_431, sof0_289, sof1_289, sog_425, sog_429, sog_430, \
                         sog_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_15 * sng_320[k]
                   + f_3 * pc_y[k] * sog_425[k];

        t_597[k] = f_16 * sng_429[k]
                   + f_6 * sof0_289[k]
                   - f_7 * sof1_289[k]
                   + f_3 * pc_x[k] * sog_429[k];

        t_598[k] = f_16 * sng_430[k]
                   + f_3 * pc_x[k] * sog_430[k];

        t_599[k] = f_16 * sng_431[k]
                   + f_3 * pc_x[k] * sog_431[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, sng_325, sng_432, sng_433, \
                         sng_434, sof0_286, sof1_286, sog_430, sog_432, sog_433, \
                         sog_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_16 * sng_432[k]
                   + f_3 * pc_x[k] * sog_432[k];

        t_601[k] = f_16 * sng_433[k]
                   + f_3 * pc_x[k] * sog_433[k];

        t_602[k] = f_16 * sng_434[k]
                   + f_3 * pc_x[k] * sog_434[k];

        t_603[k] = f_15 * sng_325[k]
                   + f_1 * sof0_286[k]
                   - f_2 * sof1_286[k]
                   + f_3 * pc_y[k] * sog_430[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pc_y, pc_z, sng_327, sng_328, sof0_288, \
                         sof0_289, sof1_288, sof1_289, sog_430, sog_432, \
                         sog_433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * sog_430[k];

        t_605[k] = f_15 * sng_327[k]
                   + f_4 * sof0_288[k]
                   - f_5 * sof1_288[k]
                   + f_3 * pc_y[k] * sog_432[k];

        t_606[k] = f_15 * sng_328[k]
                   + f_6 * sof0_289[k]
                   - f_7 * sof1_289[k]
                   + f_3 * pc_y[k] * sog_433[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, pb_z, pc_y, pc_z, snh0_441, sng_329, \
                         sng_330, snh1_441, sof0_289, sof1_289, sog_434, \
                         sog_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_15 * sng_329[k]
                   + f_3 * pc_y[k] * sog_434[k];

        t_608[k] = f_1 * sof0_289[k]
                   - f_2 * sof1_289[k]
                   + f_3 * pc_z[k] * sog_434[k];

        t_609[k] = pb_z[k] * snh0_441[k]
                   - f_8 * pc_z[k] * snh1_441[k];

        t_610[k] = f_17 * sng_330[k]
                   + f_3 * pc_y[k] * sog_435[k];
    }

#pragma omp simd aligned(t_611, t_612, t_613, pb_z, pc_y, pc_z, snh0_444, sng_315, sng_332, \
                         snh1_444, sog_435, sog_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = f_9 * sng_315[k]
                   + f_3 * pc_z[k] * sog_435[k];

        t_612[k] = pb_z[k] * snh0_444[k]
                   - f_8 * pc_z[k] * snh1_444[k];

        t_613[k] = f_17 * sng_332[k]
                   + f_3 * pc_y[k] * sog_437[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, pb_z, pc_x, pc_z, snh0_447, sng_318, sng_440, \
                         snh1_447, sof0_295, sof1_295, sog_438, \
                         sog_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_16 * sng_440[k]
                   + f_4 * sof0_295[k]
                   - f_5 * sof1_295[k]
                   + f_3 * pc_x[k] * sog_440[k];

        t_615[k] = pb_z[k] * snh0_447[k]
                   - f_8 * pc_z[k] * snh1_447[k];

        t_616[k] = f_9 * sng_318[k]
                   + f_3 * pc_z[k] * sog_438[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, pc_x, pc_y, sng_335, sng_444, sng_445, \
                         sng_446, sof0_299, sof1_299, sog_440, sog_444, sog_445, \
                         sog_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_17 * sng_335[k]
                   + f_3 * pc_y[k] * sog_440[k];

        t_618[k] = f_16 * sng_444[k]
                   + f_6 * sof0_299[k]
                   - f_7 * sof1_299[k]
                   + f_3 * pc_x[k] * sog_444[k];

        t_619[k] = f_16 * sng_445[k]
                   + f_3 * pc_x[k] * sog_445[k];

        t_620[k] = f_16 * sng_446[k]
                   + f_3 * pc_x[k] * sog_446[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, t_624, pb_z, pc_x, pc_z, snh0_456, sng_447, \
                         sng_448, sng_449, snh1_456, sog_447, sog_448, \
                         sog_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_16 * sng_447[k]
                   + f_3 * pc_x[k] * sog_447[k];

        t_622[k] = f_16 * sng_448[k]
                   + f_3 * pc_x[k] * sog_448[k];

        t_623[k] = f_16 * sng_449[k]
                   + f_3 * pc_x[k] * sog_449[k];

        t_624[k] = pb_z[k] * snh0_456[k]
                   - f_8 * pc_z[k] * snh1_456[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, pc_y, pc_z, sng_325, sng_342, sng_343, sof0_298, \
                         sof0_299, sof1_298, sof1_299, sog_445, sog_447, \
                         sog_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = f_9 * sng_325[k]
                   + f_3 * pc_z[k] * sog_445[k];

        t_626[k] = f_17 * sng_342[k]
                   + f_4 * sof0_298[k]
                   - f_5 * sof1_298[k]
                   + f_3 * pc_y[k] * sog_447[k];

        t_627[k] = f_17 * sng_343[k]
                   + f_6 * sof0_299[k]
                   - f_7 * sof1_299[k]
                   + f_3 * pc_y[k] * sog_448[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, pc_x, pc_y, pc_z, sng_329, sng_344, sng_450, \
                         sof0_299, sof0_300, sof1_299, sof1_300, sog_449, \
                         sog_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_17 * sng_344[k]
                   + f_3 * pc_y[k] * sog_449[k];

        t_629[k] = f_9 * sng_329[k]
                   + f_1 * sof0_299[k]
                   - f_2 * sof1_299[k]
                   + f_3 * pc_z[k] * sog_449[k];

        t_630[k] = f_16 * sng_450[k]
                   + f_1 * sof0_300[k]
                   - f_2 * sof1_300[k]
                   + f_3 * pc_x[k] * sog_450[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, pc_x, pc_y, pc_z, sng_330, sng_345, \
                         sng_347, sng_453, sof0_303, sof1_303, sog_450, sog_452, \
                         sog_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = f_18 * sng_345[k]
                   + f_3 * pc_y[k] * sog_450[k];

        t_632[k] = f_10 * sng_330[k]
                   + f_3 * pc_z[k] * sog_450[k];

        t_633[k] = f_16 * sng_453[k]
                   + f_4 * sof0_303[k]
                   - f_5 * sof1_303[k]
                   + f_3 * pc_x[k] * sog_453[k];

        t_634[k] = f_18 * sng_347[k]
                   + f_3 * pc_y[k] * sog_452[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, pc_x, pc_z, sng_333, sng_455, sng_456, sof0_305, \
                         sof0_306, sof1_305, sof1_306, sog_453, sog_455, \
                         sog_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_16 * sng_455[k]
                   + f_4 * sof0_305[k]
                   - f_5 * sof1_305[k]
                   + f_3 * pc_x[k] * sog_455[k];

        t_636[k] = f_16 * sng_456[k]
                   + f_6 * sof0_306[k]
                   - f_7 * sof1_306[k]
                   + f_3 * pc_x[k] * sog_456[k];

        t_637[k] = f_10 * sng_333[k]
                   + f_3 * pc_z[k] * sog_453[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pc_x, pc_y, sng_350, sng_459, sng_460, \
                         sng_461, sof0_309, sof1_309, sog_455, sog_459, sog_460, \
                         sog_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_18 * sng_350[k]
                   + f_3 * pc_y[k] * sog_455[k];

        t_639[k] = f_16 * sng_459[k]
                   + f_6 * sof0_309[k]
                   - f_7 * sof1_309[k]
                   + f_3 * pc_x[k] * sog_459[k];

        t_640[k] = f_16 * sng_460[k]
                   + f_3 * pc_x[k] * sog_460[k];

        t_641[k] = f_16 * sng_461[k]
                   + f_3 * pc_x[k] * sog_461[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pc_x, pc_y, sng_355, sng_462, sng_463, \
                         sng_464, sof0_306, sof1_306, sog_460, sog_462, sog_463, \
                         sog_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_16 * sng_462[k]
                   + f_3 * pc_x[k] * sog_462[k];

        t_643[k] = f_16 * sng_463[k]
                   + f_3 * pc_x[k] * sog_463[k];

        t_644[k] = f_16 * sng_464[k]
                   + f_3 * pc_x[k] * sog_464[k];

        t_645[k] = f_18 * sng_355[k]
                   + f_1 * sof0_306[k]
                   - f_2 * sof1_306[k]
                   + f_3 * pc_y[k] * sog_460[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_y, pc_z, sng_340, sng_357, sng_358, sof0_308, \
                         sof0_309, sof1_308, sof1_309, sog_460, sog_462, \
                         sog_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_10 * sng_340[k]
                   + f_3 * pc_z[k] * sog_460[k];

        t_647[k] = f_18 * sng_357[k]
                   + f_4 * sof0_308[k]
                   - f_5 * sof1_308[k]
                   + f_3 * pc_y[k] * sog_462[k];

        t_648[k] = f_18 * sng_358[k]
                   + f_6 * sof0_309[k]
                   - f_7 * sof1_309[k]
                   + f_3 * pc_y[k] * sog_463[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, pc_x, pc_y, pc_z, sng_344, sng_359, sng_465, \
                         sof0_309, sof0_310, sof1_309, sof1_310, sog_464, \
                         sog_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_18 * sng_359[k]
                   + f_3 * pc_y[k] * sog_464[k];

        t_650[k] = f_10 * sng_344[k]
                   + f_1 * sof0_309[k]
                   - f_2 * sof1_309[k]
                   + f_3 * pc_z[k] * sog_464[k];

        t_651[k] = f_16 * sng_465[k]
                   + f_1 * sof0_310[k]
                   - f_2 * sof1_310[k]
                   + f_3 * pc_x[k] * sog_465[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, pc_x, pc_y, pc_z, sng_345, sng_360, \
                         sng_362, sng_468, sof0_313, sof1_313, sog_465, sog_467, \
                         sog_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_16 * sng_360[k]
                   + f_3 * pc_y[k] * sog_465[k];

        t_653[k] = f_11 * sng_345[k]
                   + f_3 * pc_z[k] * sog_465[k];

        t_654[k] = f_16 * sng_468[k]
                   + f_4 * sof0_313[k]
                   - f_5 * sof1_313[k]
                   + f_3 * pc_x[k] * sog_468[k];

        t_655[k] = f_16 * sng_362[k]
                   + f_3 * pc_y[k] * sog_467[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_z, sng_348, sng_470, sng_471, sof0_315, \
                         sof0_316, sof1_315, sof1_316, sog_468, sog_470, \
                         sog_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_16 * sng_470[k]
                   + f_4 * sof0_315[k]
                   - f_5 * sof1_315[k]
                   + f_3 * pc_x[k] * sog_470[k];

        t_657[k] = f_16 * sng_471[k]
                   + f_6 * sof0_316[k]
                   - f_7 * sof1_316[k]
                   + f_3 * pc_x[k] * sog_471[k];

        t_658[k] = f_11 * sng_348[k]
                   + f_3 * pc_z[k] * sog_468[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, pc_x, pc_y, sng_365, sng_474, sng_475, \
                         sng_476, sof0_319, sof1_319, sog_470, sog_474, sog_475, \
                         sog_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_16 * sng_365[k]
                   + f_3 * pc_y[k] * sog_470[k];

        t_660[k] = f_16 * sng_474[k]
                   + f_6 * sof0_319[k]
                   - f_7 * sof1_319[k]
                   + f_3 * pc_x[k] * sog_474[k];

        t_661[k] = f_16 * sng_475[k]
                   + f_3 * pc_x[k] * sog_475[k];

        t_662[k] = f_16 * sng_476[k]
                   + f_3 * pc_x[k] * sog_476[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, t_666, pc_x, pc_y, sng_370, sng_477, sng_478, \
                         sng_479, sof0_316, sof1_316, sog_475, sog_477, sog_478, \
                         sog_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_16 * sng_477[k]
                   + f_3 * pc_x[k] * sog_477[k];

        t_664[k] = f_16 * sng_478[k]
                   + f_3 * pc_x[k] * sog_478[k];

        t_665[k] = f_16 * sng_479[k]
                   + f_3 * pc_x[k] * sog_479[k];

        t_666[k] = f_16 * sng_370[k]
                   + f_1 * sof0_316[k]
                   - f_2 * sof1_316[k]
                   + f_3 * pc_y[k] * sog_475[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, pc_z, sng_355, sng_372, sng_373, sof0_318, \
                         sof0_319, sof1_318, sof1_319, sog_475, sog_477, \
                         sog_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_11 * sng_355[k]
                   + f_3 * pc_z[k] * sog_475[k];

        t_668[k] = f_16 * sng_372[k]
                   + f_4 * sof0_318[k]
                   - f_5 * sof1_318[k]
                   + f_3 * pc_y[k] * sog_477[k];

        t_669[k] = f_16 * sng_373[k]
                   + f_6 * sof0_319[k]
                   - f_7 * sof1_319[k]
                   + f_3 * pc_y[k] * sog_478[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, sng_359, sng_374, sng_480, \
                         sof0_319, sof0_320, sof1_319, sof1_320, sog_479, \
                         sog_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_16 * sng_374[k]
                   + f_3 * pc_y[k] * sog_479[k];

        t_671[k] = f_11 * sng_359[k]
                   + f_1 * sof0_319[k]
                   - f_2 * sof1_319[k]
                   + f_3 * pc_z[k] * sog_479[k];

        t_672[k] = f_16 * sng_480[k]
                   + f_1 * sof0_320[k]
                   - f_2 * sof1_320[k]
                   + f_3 * pc_x[k] * sog_480[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, sng_360, sng_375, \
                         sng_377, sng_483, sof0_323, sof1_323, sog_480, sog_482, \
                         sog_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_11 * sng_375[k]
                   + f_3 * pc_y[k] * sog_480[k];

        t_674[k] = f_16 * sng_360[k]
                   + f_3 * pc_z[k] * sog_480[k];

        t_675[k] = f_16 * sng_483[k]
                   + f_4 * sof0_323[k]
                   - f_5 * sof1_323[k]
                   + f_3 * pc_x[k] * sog_483[k];

        t_676[k] = f_11 * sng_377[k]
                   + f_3 * pc_y[k] * sog_482[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, sng_363, sng_485, sng_486, sof0_325, \
                         sof0_326, sof1_325, sof1_326, sog_483, sog_485, \
                         sog_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_16 * sng_485[k]
                   + f_4 * sof0_325[k]
                   - f_5 * sof1_325[k]
                   + f_3 * pc_x[k] * sog_485[k];

        t_678[k] = f_16 * sng_486[k]
                   + f_6 * sof0_326[k]
                   - f_7 * sof1_326[k]
                   + f_3 * pc_x[k] * sog_486[k];

        t_679[k] = f_16 * sng_363[k]
                   + f_3 * pc_z[k] * sog_483[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pc_x, pc_y, sng_380, sng_489, sng_490, \
                         sng_491, sof0_329, sof1_329, sog_485, sog_489, sog_490, \
                         sog_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_11 * sng_380[k]
                   + f_3 * pc_y[k] * sog_485[k];

        t_681[k] = f_16 * sng_489[k]
                   + f_6 * sof0_329[k]
                   - f_7 * sof1_329[k]
                   + f_3 * pc_x[k] * sog_489[k];

        t_682[k] = f_16 * sng_490[k]
                   + f_3 * pc_x[k] * sog_490[k];

        t_683[k] = f_16 * sng_491[k]
                   + f_3 * pc_x[k] * sog_491[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, pc_x, pc_y, sng_385, sng_492, sng_493, \
                         sng_494, sof0_326, sof1_326, sog_490, sog_492, sog_493, \
                         sog_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_16 * sng_492[k]
                   + f_3 * pc_x[k] * sog_492[k];

        t_685[k] = f_16 * sng_493[k]
                   + f_3 * pc_x[k] * sog_493[k];

        t_686[k] = f_16 * sng_494[k]
                   + f_3 * pc_x[k] * sog_494[k];

        t_687[k] = f_11 * sng_385[k]
                   + f_1 * sof0_326[k]
                   - f_2 * sof1_326[k]
                   + f_3 * pc_y[k] * sog_490[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, pc_y, pc_z, sng_370, sng_387, sng_388, sof0_328, \
                         sof0_329, sof1_328, sof1_329, sog_490, sog_492, \
                         sog_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_16 * sng_370[k]
                   + f_3 * pc_z[k] * sog_490[k];

        t_689[k] = f_11 * sng_387[k]
                   + f_4 * sof0_328[k]
                   - f_5 * sof1_328[k]
                   + f_3 * pc_y[k] * sog_492[k];

        t_690[k] = f_11 * sng_388[k]
                   + f_6 * sof0_329[k]
                   - f_7 * sof1_329[k]
                   + f_3 * pc_y[k] * sog_493[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, pc_x, pc_y, pc_z, sng_374, sng_389, sng_495, \
                         sof0_329, sof0_330, sof1_329, sof1_330, sog_494, \
                         sog_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_11 * sng_389[k]
                   + f_3 * pc_y[k] * sog_494[k];

        t_692[k] = f_16 * sng_374[k]
                   + f_1 * sof0_329[k]
                   - f_2 * sof1_329[k]
                   + f_3 * pc_z[k] * sog_494[k];

        t_693[k] = f_16 * sng_495[k]
                   + f_1 * sof0_330[k]
                   - f_2 * sof1_330[k]
                   + f_3 * pc_x[k] * sog_495[k];
    }
}

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snh0,
                                                          const size_t sng, const size_t snh1,
                                                          const size_t sof0, const size_t sof1,
                                                          const size_t sog, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_14 = 4.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_567 = buffer.data(snh0 + 567);
    const auto *snh0_570 = buffer.data(snh0 + 570);
    const auto *snh0_572 = buffer.data(snh0 + 572);
    const auto *snh0_573 = buffer.data(snh0 + 573);
    const auto *snh0_576 = buffer.data(snh0 + 576);
    const auto *snh0_587 = buffer.data(snh0 + 587);
    const auto *snh0_588 = buffer.data(snh0 + 588);
    const auto *snh0_591 = buffer.data(snh0 + 591);
    const auto *snh0_594 = buffer.data(snh0 + 594);
    const auto *snh0_603 = buffer.data(snh0 + 603);

    const auto *sng_375 = buffer.data(sng + 375);
    const auto *sng_378 = buffer.data(sng + 378);
    const auto *sng_385 = buffer.data(sng + 385);
    const auto *sng_389 = buffer.data(sng + 389);
    const auto *sng_390 = buffer.data(sng + 390);
    const auto *sng_392 = buffer.data(sng + 392);
    const auto *sng_393 = buffer.data(sng + 393);
    const auto *sng_395 = buffer.data(sng + 395);
    const auto *sng_400 = buffer.data(sng + 400);
    const auto *sng_402 = buffer.data(sng + 402);
    const auto *sng_403 = buffer.data(sng + 403);
    const auto *sng_404 = buffer.data(sng + 404);
    const auto *sng_405 = buffer.data(sng + 405);
    const auto *sng_406 = buffer.data(sng + 406);
    const auto *sng_407 = buffer.data(sng + 407);
    const auto *sng_408 = buffer.data(sng + 408);
    const auto *sng_410 = buffer.data(sng + 410);
    const auto *sng_415 = buffer.data(sng + 415);
    const auto *sng_417 = buffer.data(sng + 417);
    const auto *sng_418 = buffer.data(sng + 418);
    const auto *sng_419 = buffer.data(sng + 419);
    const auto *sng_420 = buffer.data(sng + 420);
    const auto *sng_422 = buffer.data(sng + 422);
    const auto *sng_423 = buffer.data(sng + 423);
    const auto *sng_425 = buffer.data(sng + 425);
    const auto *sng_430 = buffer.data(sng + 430);
    const auto *sng_432 = buffer.data(sng + 432);
    const auto *sng_433 = buffer.data(sng + 433);
    const auto *sng_434 = buffer.data(sng + 434);
    const auto *sng_435 = buffer.data(sng + 435);
    const auto *sng_437 = buffer.data(sng + 437);
    const auto *sng_438 = buffer.data(sng + 438);
    const auto *sng_440 = buffer.data(sng + 440);
    const auto *sng_447 = buffer.data(sng + 447);
    const auto *sng_448 = buffer.data(sng + 448);
    const auto *sng_449 = buffer.data(sng + 449);
    const auto *sng_450 = buffer.data(sng + 450);
    const auto *sng_452 = buffer.data(sng + 452);
    const auto *sng_455 = buffer.data(sng + 455);
    const auto *sng_498 = buffer.data(sng + 498);
    const auto *sng_500 = buffer.data(sng + 500);
    const auto *sng_501 = buffer.data(sng + 501);
    const auto *sng_504 = buffer.data(sng + 504);
    const auto *sng_505 = buffer.data(sng + 505);
    const auto *sng_506 = buffer.data(sng + 506);
    const auto *sng_507 = buffer.data(sng + 507);
    const auto *sng_508 = buffer.data(sng + 508);
    const auto *sng_509 = buffer.data(sng + 509);
    const auto *sng_520 = buffer.data(sng + 520);
    const auto *sng_521 = buffer.data(sng + 521);
    const auto *sng_522 = buffer.data(sng + 522);
    const auto *sng_523 = buffer.data(sng + 523);
    const auto *sng_524 = buffer.data(sng + 524);
    const auto *sng_525 = buffer.data(sng + 525);
    const auto *sng_528 = buffer.data(sng + 528);
    const auto *sng_530 = buffer.data(sng + 530);
    const auto *sng_531 = buffer.data(sng + 531);
    const auto *sng_534 = buffer.data(sng + 534);
    const auto *sng_535 = buffer.data(sng + 535);
    const auto *sng_536 = buffer.data(sng + 536);
    const auto *sng_537 = buffer.data(sng + 537);
    const auto *sng_538 = buffer.data(sng + 538);
    const auto *sng_539 = buffer.data(sng + 539);
    const auto *sng_540 = buffer.data(sng + 540);
    const auto *sng_543 = buffer.data(sng + 543);
    const auto *sng_545 = buffer.data(sng + 545);
    const auto *sng_546 = buffer.data(sng + 546);
    const auto *sng_549 = buffer.data(sng + 549);
    const auto *sng_550 = buffer.data(sng + 550);
    const auto *sng_551 = buffer.data(sng + 551);
    const auto *sng_552 = buffer.data(sng + 552);
    const auto *sng_553 = buffer.data(sng + 553);
    const auto *sng_554 = buffer.data(sng + 554);
    const auto *sng_560 = buffer.data(sng + 560);
    const auto *sng_564 = buffer.data(sng + 564);
    const auto *sng_565 = buffer.data(sng + 565);
    const auto *sng_566 = buffer.data(sng + 566);
    const auto *sng_567 = buffer.data(sng + 567);
    const auto *sng_568 = buffer.data(sng + 568);
    const auto *sng_569 = buffer.data(sng + 569);
    const auto *sng_570 = buffer.data(sng + 570);
    const auto *sng_573 = buffer.data(sng + 573);
    const auto *sng_575 = buffer.data(sng + 575);
    const auto *sng_576 = buffer.data(sng + 576);
    const auto *sng_579 = buffer.data(sng + 579);
    const auto *sng_580 = buffer.data(sng + 580);
    const auto *sng_581 = buffer.data(sng + 581);

    const auto *snh1_567 = buffer.data(snh1 + 567);
    const auto *snh1_570 = buffer.data(snh1 + 570);
    const auto *snh1_572 = buffer.data(snh1 + 572);
    const auto *snh1_573 = buffer.data(snh1 + 573);
    const auto *snh1_576 = buffer.data(snh1 + 576);
    const auto *snh1_587 = buffer.data(snh1 + 587);
    const auto *snh1_588 = buffer.data(snh1 + 588);
    const auto *snh1_591 = buffer.data(snh1 + 591);
    const auto *snh1_594 = buffer.data(snh1 + 594);
    const auto *snh1_603 = buffer.data(snh1 + 603);

    const auto *sof0_333 = buffer.data(sof0 + 333);
    const auto *sof0_335 = buffer.data(sof0 + 335);
    const auto *sof0_336 = buffer.data(sof0 + 336);
    const auto *sof0_338 = buffer.data(sof0 + 338);
    const auto *sof0_339 = buffer.data(sof0 + 339);
    const auto *sof0_346 = buffer.data(sof0 + 346);
    const auto *sof0_348 = buffer.data(sof0 + 348);
    const auto *sof0_349 = buffer.data(sof0 + 349);
    const auto *sof0_350 = buffer.data(sof0 + 350);
    const auto *sof0_353 = buffer.data(sof0 + 353);
    const auto *sof0_355 = buffer.data(sof0 + 355);
    const auto *sof0_356 = buffer.data(sof0 + 356);
    const auto *sof0_358 = buffer.data(sof0 + 358);
    const auto *sof0_359 = buffer.data(sof0 + 359);
    const auto *sof0_360 = buffer.data(sof0 + 360);
    const auto *sof0_363 = buffer.data(sof0 + 363);
    const auto *sof0_365 = buffer.data(sof0 + 365);
    const auto *sof0_366 = buffer.data(sof0 + 366);
    const auto *sof0_368 = buffer.data(sof0 + 368);
    const auto *sof0_369 = buffer.data(sof0 + 369);
    const auto *sof0_375 = buffer.data(sof0 + 375);
    const auto *sof0_378 = buffer.data(sof0 + 378);
    const auto *sof0_379 = buffer.data(sof0 + 379);
    const auto *sof0_380 = buffer.data(sof0 + 380);
    const auto *sof0_383 = buffer.data(sof0 + 383);
    const auto *sof0_385 = buffer.data(sof0 + 385);
    const auto *sof0_386 = buffer.data(sof0 + 386);
    const auto *sof0_389 = buffer.data(sof0 + 389);

    const auto *sof1_333 = buffer.data(sof1 + 333);
    const auto *sof1_335 = buffer.data(sof1 + 335);
    const auto *sof1_336 = buffer.data(sof1 + 336);
    const auto *sof1_338 = buffer.data(sof1 + 338);
    const auto *sof1_339 = buffer.data(sof1 + 339);
    const auto *sof1_346 = buffer.data(sof1 + 346);
    const auto *sof1_348 = buffer.data(sof1 + 348);
    const auto *sof1_349 = buffer.data(sof1 + 349);
    const auto *sof1_350 = buffer.data(sof1 + 350);
    const auto *sof1_353 = buffer.data(sof1 + 353);
    const auto *sof1_355 = buffer.data(sof1 + 355);
    const auto *sof1_356 = buffer.data(sof1 + 356);
    const auto *sof1_358 = buffer.data(sof1 + 358);
    const auto *sof1_359 = buffer.data(sof1 + 359);
    const auto *sof1_360 = buffer.data(sof1 + 360);
    const auto *sof1_363 = buffer.data(sof1 + 363);
    const auto *sof1_365 = buffer.data(sof1 + 365);
    const auto *sof1_366 = buffer.data(sof1 + 366);
    const auto *sof1_368 = buffer.data(sof1 + 368);
    const auto *sof1_369 = buffer.data(sof1 + 369);
    const auto *sof1_375 = buffer.data(sof1 + 375);
    const auto *sof1_378 = buffer.data(sof1 + 378);
    const auto *sof1_379 = buffer.data(sof1 + 379);
    const auto *sof1_380 = buffer.data(sof1 + 380);
    const auto *sof1_383 = buffer.data(sof1 + 383);
    const auto *sof1_385 = buffer.data(sof1 + 385);
    const auto *sof1_386 = buffer.data(sof1 + 386);
    const auto *sof1_389 = buffer.data(sof1 + 389);

    const auto *sog_495 = buffer.data(sog + 495);
    const auto *sog_497 = buffer.data(sog + 497);
    const auto *sog_498 = buffer.data(sog + 498);
    const auto *sog_500 = buffer.data(sog + 500);
    const auto *sog_501 = buffer.data(sog + 501);
    const auto *sog_504 = buffer.data(sog + 504);
    const auto *sog_505 = buffer.data(sog + 505);
    const auto *sog_506 = buffer.data(sog + 506);
    const auto *sog_507 = buffer.data(sog + 507);
    const auto *sog_508 = buffer.data(sog + 508);
    const auto *sog_509 = buffer.data(sog + 509);
    const auto *sog_510 = buffer.data(sog + 510);
    const auto *sog_512 = buffer.data(sog + 512);
    const auto *sog_513 = buffer.data(sog + 513);
    const auto *sog_515 = buffer.data(sog + 515);
    const auto *sog_520 = buffer.data(sog + 520);
    const auto *sog_521 = buffer.data(sog + 521);
    const auto *sog_522 = buffer.data(sog + 522);
    const auto *sog_523 = buffer.data(sog + 523);
    const auto *sog_524 = buffer.data(sog + 524);
    const auto *sog_525 = buffer.data(sog + 525);
    const auto *sog_527 = buffer.data(sog + 527);
    const auto *sog_528 = buffer.data(sog + 528);
    const auto *sog_530 = buffer.data(sog + 530);
    const auto *sog_531 = buffer.data(sog + 531);
    const auto *sog_534 = buffer.data(sog + 534);
    const auto *sog_535 = buffer.data(sog + 535);
    const auto *sog_536 = buffer.data(sog + 536);
    const auto *sog_537 = buffer.data(sog + 537);
    const auto *sog_538 = buffer.data(sog + 538);
    const auto *sog_539 = buffer.data(sog + 539);
    const auto *sog_540 = buffer.data(sog + 540);
    const auto *sog_542 = buffer.data(sog + 542);
    const auto *sog_543 = buffer.data(sog + 543);
    const auto *sog_545 = buffer.data(sog + 545);
    const auto *sog_546 = buffer.data(sog + 546);
    const auto *sog_549 = buffer.data(sog + 549);
    const auto *sog_550 = buffer.data(sog + 550);
    const auto *sog_551 = buffer.data(sog + 551);
    const auto *sog_552 = buffer.data(sog + 552);
    const auto *sog_553 = buffer.data(sog + 553);
    const auto *sog_554 = buffer.data(sog + 554);
    const auto *sog_555 = buffer.data(sog + 555);
    const auto *sog_557 = buffer.data(sog + 557);
    const auto *sog_558 = buffer.data(sog + 558);
    const auto *sog_560 = buffer.data(sog + 560);
    const auto *sog_564 = buffer.data(sog + 564);
    const auto *sog_565 = buffer.data(sog + 565);
    const auto *sog_566 = buffer.data(sog + 566);
    const auto *sog_567 = buffer.data(sog + 567);
    const auto *sog_568 = buffer.data(sog + 568);
    const auto *sog_569 = buffer.data(sog + 569);
    const auto *sog_570 = buffer.data(sog + 570);
    const auto *sog_572 = buffer.data(sog + 572);
    const auto *sog_573 = buffer.data(sog + 573);
    const auto *sog_575 = buffer.data(sog + 575);
    const auto *sog_576 = buffer.data(sog + 576);
    const auto *sog_579 = buffer.data(sog + 579);
    const auto *sog_580 = buffer.data(sog + 580);
    const auto *sog_581 = buffer.data(sog + 581);

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pc_x, pc_y, pc_z, sng_375, sng_390, \
                         sng_392, sng_498, sof0_333, sof1_333, sog_495, sog_497, \
                         sog_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_10 * sng_390[k]
                   + f_3 * pc_y[k] * sog_495[k];

        t_695[k] = f_18 * sng_375[k]
                   + f_3 * pc_z[k] * sog_495[k];

        t_696[k] = f_16 * sng_498[k]
                   + f_4 * sof0_333[k]
                   - f_5 * sof1_333[k]
                   + f_3 * pc_x[k] * sog_498[k];

        t_697[k] = f_10 * sng_392[k]
                   + f_3 * pc_y[k] * sog_497[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pc_x, pc_z, sng_378, sng_500, sng_501, sof0_335, \
                         sof0_336, sof1_335, sof1_336, sog_498, sog_500, \
                         sog_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_16 * sng_500[k]
                   + f_4 * sof0_335[k]
                   - f_5 * sof1_335[k]
                   + f_3 * pc_x[k] * sog_500[k];

        t_699[k] = f_16 * sng_501[k]
                   + f_6 * sof0_336[k]
                   - f_7 * sof1_336[k]
                   + f_3 * pc_x[k] * sog_501[k];

        t_700[k] = f_18 * sng_378[k]
                   + f_3 * pc_z[k] * sog_498[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pc_x, pc_y, sng_395, sng_504, sng_505, \
                         sng_506, sof0_339, sof1_339, sog_500, sog_504, sog_505, \
                         sog_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_10 * sng_395[k]
                   + f_3 * pc_y[k] * sog_500[k];

        t_702[k] = f_16 * sng_504[k]
                   + f_6 * sof0_339[k]
                   - f_7 * sof1_339[k]
                   + f_3 * pc_x[k] * sog_504[k];

        t_703[k] = f_16 * sng_505[k]
                   + f_3 * pc_x[k] * sog_505[k];

        t_704[k] = f_16 * sng_506[k]
                   + f_3 * pc_x[k] * sog_506[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, pc_x, pc_y, sng_400, sng_507, sng_508, \
                         sng_509, sof0_336, sof1_336, sog_505, sog_507, sog_508, \
                         sog_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_16 * sng_507[k]
                   + f_3 * pc_x[k] * sog_507[k];

        t_706[k] = f_16 * sng_508[k]
                   + f_3 * pc_x[k] * sog_508[k];

        t_707[k] = f_16 * sng_509[k]
                   + f_3 * pc_x[k] * sog_509[k];

        t_708[k] = f_10 * sng_400[k]
                   + f_1 * sof0_336[k]
                   - f_2 * sof1_336[k]
                   + f_3 * pc_y[k] * sog_505[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_y, pc_z, sng_385, sng_402, sng_403, sof0_338, \
                         sof0_339, sof1_338, sof1_339, sog_505, sog_507, \
                         sog_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_18 * sng_385[k]
                   + f_3 * pc_z[k] * sog_505[k];

        t_710[k] = f_10 * sng_402[k]
                   + f_4 * sof0_338[k]
                   - f_5 * sof1_338[k]
                   + f_3 * pc_y[k] * sog_507[k];

        t_711[k] = f_10 * sng_403[k]
                   + f_6 * sof0_339[k]
                   - f_7 * sof1_339[k]
                   + f_3 * pc_y[k] * sog_508[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pb_y, pc_y, pc_z, snh0_567, sng_389, \
                         sng_404, sng_405, snh1_567, sof0_339, sof1_339, sog_509, \
                         sog_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_10 * sng_404[k]
                   + f_3 * pc_y[k] * sog_509[k];

        t_713[k] = f_18 * sng_389[k]
                   + f_1 * sof0_339[k]
                   - f_2 * sof1_339[k]
                   + f_3 * pc_z[k] * sog_509[k];

        t_714[k] = pb_y[k] * snh0_567[k]
                   - f_8 * pc_y[k] * snh1_567[k];

        t_715[k] = f_9 * sng_405[k]
                   + f_3 * pc_y[k] * sog_510[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pb_y, pc_y, pc_z, snh0_570, snh0_572, \
                         sng_390, sng_406, sng_407, snh1_570, snh1_572, sog_510, \
                         sog_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_17 * sng_390[k]
                   + f_3 * pc_z[k] * sog_510[k];

        t_717[k] = pb_y[k] * snh0_570[k]
                   + f_10 * sng_406[k]
                   - f_8 * pc_y[k] * snh1_570[k];

        t_718[k] = f_9 * sng_407[k]
                   + f_3 * pc_y[k] * sog_512[k];

        t_719[k] = pb_y[k] * snh0_572[k]
                   - f_8 * pc_y[k] * snh1_572[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pb_y, pc_y, pc_z, snh0_573, snh0_576, \
                         sng_393, sng_408, sng_410, snh1_573, snh1_576, sog_513, \
                         sog_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pb_y[k] * snh0_573[k]
                   + f_11 * sng_408[k]
                   - f_8 * pc_y[k] * snh1_573[k];

        t_721[k] = f_17 * sng_393[k]
                   + f_3 * pc_z[k] * sog_513[k];

        t_722[k] = f_9 * sng_410[k]
                   + f_3 * pc_y[k] * sog_515[k];

        t_723[k] = pb_y[k] * snh0_576[k]
                   - f_8 * pc_y[k] * snh1_576[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, pc_x, sng_520, sng_521, sng_522, \
                         sng_523, sng_524, sog_520, sog_521, sog_522, sog_523, \
                         sog_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_16 * sng_520[k]
                   + f_3 * pc_x[k] * sog_520[k];

        t_725[k] = f_16 * sng_521[k]
                   + f_3 * pc_x[k] * sog_521[k];

        t_726[k] = f_16 * sng_522[k]
                   + f_3 * pc_x[k] * sog_522[k];

        t_727[k] = f_16 * sng_523[k]
                   + f_3 * pc_x[k] * sog_523[k];

        t_728[k] = f_16 * sng_524[k]
                   + f_3 * pc_x[k] * sog_524[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_y, pc_z, sng_400, sng_415, sng_417, sof0_346, \
                         sof0_348, sof1_346, sof1_348, sog_520, \
                         sog_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_9 * sng_415[k]
                   + f_1 * sof0_346[k]
                   - f_2 * sof1_346[k]
                   + f_3 * pc_y[k] * sog_520[k];

        t_730[k] = f_17 * sng_400[k]
                   + f_3 * pc_z[k] * sog_520[k];

        t_731[k] = f_9 * sng_417[k]
                   + f_4 * sof0_348[k]
                   - f_5 * sof1_348[k]
                   + f_3 * pc_y[k] * sog_522[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pb_y, pc_y, snh0_587, sng_418, sng_419, \
                         snh1_587, sof0_349, sof1_349, sog_523, \
                         sog_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_9 * sng_418[k]
                   + f_6 * sof0_349[k]
                   - f_7 * sof1_349[k]
                   + f_3 * pc_y[k] * sog_523[k];

        t_733[k] = f_9 * sng_419[k]
                   + f_3 * pc_y[k] * sog_524[k];

        t_734[k] = pb_y[k] * snh0_587[k]
                   - f_8 * pc_y[k] * snh1_587[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, pc_x, pc_y, pc_z, sng_405, sng_525, \
                         sng_528, sof0_350, sof0_353, sof1_350, sof1_353, sog_525, \
                         sog_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_16 * sng_525[k]
                   + f_1 * sof0_350[k]
                   - f_2 * sof1_350[k]
                   + f_3 * pc_x[k] * sog_525[k];

        t_736[k] = f_3 * pc_y[k] * sog_525[k];

        t_737[k] = f_15 * sng_405[k]
                   + f_3 * pc_z[k] * sog_525[k];

        t_738[k] = f_16 * sng_528[k]
                   + f_4 * sof0_353[k]
                   - f_5 * sof1_353[k]
                   + f_3 * pc_x[k] * sog_528[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, pc_x, pc_y, sng_530, sng_531, sof0_355, \
                         sof0_356, sof1_355, sof1_356, sog_527, sog_530, \
                         sog_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_3 * pc_y[k] * sog_527[k];

        t_740[k] = f_16 * sng_530[k]
                   + f_4 * sof0_355[k]
                   - f_5 * sof1_355[k]
                   + f_3 * pc_x[k] * sog_530[k];

        t_741[k] = f_16 * sng_531[k]
                   + f_6 * sof0_356[k]
                   - f_7 * sof1_356[k]
                   + f_3 * pc_x[k] * sog_531[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, pc_x, pc_y, pc_z, sng_408, sng_534, \
                         sng_535, sof0_359, sof1_359, sog_528, sog_530, sog_534, \
                         sog_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_15 * sng_408[k]
                   + f_3 * pc_z[k] * sog_528[k];

        t_743[k] = f_3 * pc_y[k] * sog_530[k];

        t_744[k] = f_16 * sng_534[k]
                   + f_6 * sof0_359[k]
                   - f_7 * sof1_359[k]
                   + f_3 * pc_x[k] * sog_534[k];

        t_745[k] = f_16 * sng_535[k]
                   + f_3 * pc_x[k] * sog_535[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pc_x, sng_536, sng_537, sng_538, sng_539, \
                         sog_536, sog_537, sog_538, sog_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_16 * sng_536[k]
                   + f_3 * pc_x[k] * sog_536[k];

        t_747[k] = f_16 * sng_537[k]
                   + f_3 * pc_x[k] * sog_537[k];

        t_748[k] = f_16 * sng_538[k]
                   + f_3 * pc_x[k] * sog_538[k];

        t_749[k] = f_16 * sng_539[k]
                   + f_3 * pc_x[k] * sog_539[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, pc_y, pc_z, sng_415, sof0_356, sof0_358, \
                         sof0_359, sof1_356, sof1_358, sof1_359, sog_535, sog_537, \
                         sog_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_1 * sof0_356[k]
                   - f_2 * sof1_356[k]
                   + f_3 * pc_y[k] * sog_535[k];

        t_751[k] = f_15 * sng_415[k]
                   + f_3 * pc_z[k] * sog_535[k];

        t_752[k] = f_4 * sof0_358[k]
                   - f_5 * sof1_358[k]
                   + f_3 * pc_y[k] * sog_537[k];

        t_753[k] = f_6 * sof0_359[k]
                   - f_7 * sof1_359[k]
                   + f_3 * pc_y[k] * sog_538[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, pc_x, pc_y, pc_z, sng_419, sng_420, \
                         sng_540, sof0_359, sof0_360, sof1_359, sof1_360, sog_539, \
                         sog_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = f_3 * pc_y[k] * sog_539[k];

        t_755[k] = f_15 * sng_419[k]
                   + f_1 * sof0_359[k]
                   - f_2 * sof1_359[k]
                   + f_3 * pc_z[k] * sog_539[k];

        t_756[k] = f_11 * sng_540[k]
                   + f_1 * sof0_360[k]
                   - f_2 * sof1_360[k]
                   + f_3 * pc_x[k] * sog_540[k];

        t_757[k] = f_14 * sng_420[k]
                   + f_3 * pc_y[k] * sog_540[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, pc_x, pc_y, pc_z, sng_422, sng_543, sof0_363, \
                         sof1_363, sog_540, sog_542, sog_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = f_3 * pc_z[k] * sog_540[k];

        t_759[k] = f_11 * sng_543[k]
                   + f_4 * sof0_363[k]
                   - f_5 * sof1_363[k]
                   + f_3 * pc_x[k] * sog_543[k];

        t_760[k] = f_14 * sng_422[k]
                   + f_3 * pc_y[k] * sog_542[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pc_x, pc_z, sng_545, sng_546, sof0_365, \
                         sof0_366, sof1_365, sof1_366, sog_543, sog_545, \
                         sog_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_11 * sng_545[k]
                   + f_4 * sof0_365[k]
                   - f_5 * sof1_365[k]
                   + f_3 * pc_x[k] * sog_545[k];

        t_762[k] = f_11 * sng_546[k]
                   + f_6 * sof0_366[k]
                   - f_7 * sof1_366[k]
                   + f_3 * pc_x[k] * sog_546[k];

        t_763[k] = f_3 * pc_z[k] * sog_543[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pc_x, pc_y, sng_425, sng_549, sng_550, \
                         sng_551, sof0_369, sof1_369, sog_545, sog_549, sog_550, \
                         sog_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_14 * sng_425[k]
                   + f_3 * pc_y[k] * sog_545[k];

        t_765[k] = f_11 * sng_549[k]
                   + f_6 * sof0_369[k]
                   - f_7 * sof1_369[k]
                   + f_3 * pc_x[k] * sog_549[k];

        t_766[k] = f_11 * sng_550[k]
                   + f_3 * pc_x[k] * sog_550[k];

        t_767[k] = f_11 * sng_551[k]
                   + f_3 * pc_x[k] * sog_551[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pc_x, pc_y, sng_430, sng_552, sng_553, \
                         sng_554, sof0_366, sof1_366, sog_550, sog_552, sog_553, \
                         sog_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_11 * sng_552[k]
                   + f_3 * pc_x[k] * sog_552[k];

        t_769[k] = f_11 * sng_553[k]
                   + f_3 * pc_x[k] * sog_553[k];

        t_770[k] = f_11 * sng_554[k]
                   + f_3 * pc_x[k] * sog_554[k];

        t_771[k] = f_14 * sng_430[k]
                   + f_1 * sof0_366[k]
                   - f_2 * sof1_366[k]
                   + f_3 * pc_y[k] * sog_550[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, pc_y, pc_z, sng_432, sng_433, sof0_368, \
                         sof0_369, sof1_368, sof1_369, sog_550, sog_552, \
                         sog_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_z[k] * sog_550[k];

        t_773[k] = f_14 * sng_432[k]
                   + f_4 * sof0_368[k]
                   - f_5 * sof1_368[k]
                   + f_3 * pc_y[k] * sog_552[k];

        t_774[k] = f_14 * sng_433[k]
                   + f_6 * sof0_369[k]
                   - f_7 * sof1_369[k]
                   + f_3 * pc_y[k] * sog_553[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, pb_z, pc_y, pc_z, snh0_588, sng_434, \
                         sng_435, snh1_588, sof0_369, sof1_369, sog_554, \
                         sog_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_14 * sng_434[k]
                   + f_3 * pc_y[k] * sog_554[k];

        t_776[k] = f_1 * sof0_369[k]
                   - f_2 * sof1_369[k]
                   + f_3 * pc_z[k] * sog_554[k];

        t_777[k] = pb_z[k] * snh0_588[k]
                   - f_8 * pc_z[k] * snh1_588[k];

        t_778[k] = f_15 * sng_435[k]
                   + f_3 * pc_y[k] * sog_555[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pb_z, pc_y, pc_z, snh0_591, sng_420, sng_437, \
                         snh1_591, sog_555, sog_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_9 * sng_420[k]
                   + f_3 * pc_z[k] * sog_555[k];

        t_780[k] = pb_z[k] * snh0_591[k]
                   - f_8 * pc_z[k] * snh1_591[k];

        t_781[k] = f_15 * sng_437[k]
                   + f_3 * pc_y[k] * sog_557[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, pb_z, pc_x, pc_z, snh0_594, sng_423, sng_560, \
                         snh1_594, sof0_375, sof1_375, sog_558, \
                         sog_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_11 * sng_560[k]
                   + f_4 * sof0_375[k]
                   - f_5 * sof1_375[k]
                   + f_3 * pc_x[k] * sog_560[k];

        t_783[k] = pb_z[k] * snh0_594[k]
                   - f_8 * pc_z[k] * snh1_594[k];

        t_784[k] = f_9 * sng_423[k]
                   + f_3 * pc_z[k] * sog_558[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, pc_x, pc_y, sng_440, sng_564, sng_565, \
                         sng_566, sof0_379, sof1_379, sog_560, sog_564, sog_565, \
                         sog_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_15 * sng_440[k]
                   + f_3 * pc_y[k] * sog_560[k];

        t_786[k] = f_11 * sng_564[k]
                   + f_6 * sof0_379[k]
                   - f_7 * sof1_379[k]
                   + f_3 * pc_x[k] * sog_564[k];

        t_787[k] = f_11 * sng_565[k]
                   + f_3 * pc_x[k] * sog_565[k];

        t_788[k] = f_11 * sng_566[k]
                   + f_3 * pc_x[k] * sog_566[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, t_792, pb_z, pc_x, pc_z, snh0_603, sng_567, \
                         sng_568, sng_569, snh1_603, sog_567, sog_568, \
                         sog_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = f_11 * sng_567[k]
                   + f_3 * pc_x[k] * sog_567[k];

        t_790[k] = f_11 * sng_568[k]
                   + f_3 * pc_x[k] * sog_568[k];

        t_791[k] = f_11 * sng_569[k]
                   + f_3 * pc_x[k] * sog_569[k];

        t_792[k] = pb_z[k] * snh0_603[k]
                   - f_8 * pc_z[k] * snh1_603[k];
    }

#pragma omp simd aligned(t_793, t_794, t_795, pc_y, pc_z, sng_430, sng_447, sng_448, sof0_378, \
                         sof0_379, sof1_378, sof1_379, sog_565, sog_567, \
                         sog_568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_793[k] = f_9 * sng_430[k]
                   + f_3 * pc_z[k] * sog_565[k];

        t_794[k] = f_15 * sng_447[k]
                   + f_4 * sof0_378[k]
                   - f_5 * sof1_378[k]
                   + f_3 * pc_y[k] * sog_567[k];

        t_795[k] = f_15 * sng_448[k]
                   + f_6 * sof0_379[k]
                   - f_7 * sof1_379[k]
                   + f_3 * pc_y[k] * sog_568[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, pc_x, pc_y, pc_z, sng_434, sng_449, sng_570, \
                         sof0_379, sof0_380, sof1_379, sof1_380, sog_569, \
                         sog_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_15 * sng_449[k]
                   + f_3 * pc_y[k] * sog_569[k];

        t_797[k] = f_9 * sng_434[k]
                   + f_1 * sof0_379[k]
                   - f_2 * sof1_379[k]
                   + f_3 * pc_z[k] * sog_569[k];

        t_798[k] = f_11 * sng_570[k]
                   + f_1 * sof0_380[k]
                   - f_2 * sof1_380[k]
                   + f_3 * pc_x[k] * sog_570[k];
    }

#pragma omp simd aligned(t_799, t_800, t_801, t_802, pc_x, pc_y, pc_z, sng_435, sng_450, \
                         sng_452, sng_573, sof0_383, sof1_383, sog_570, sog_572, \
                         sog_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_799[k] = f_17 * sng_450[k]
                   + f_3 * pc_y[k] * sog_570[k];

        t_800[k] = f_10 * sng_435[k]
                   + f_3 * pc_z[k] * sog_570[k];

        t_801[k] = f_11 * sng_573[k]
                   + f_4 * sof0_383[k]
                   - f_5 * sof1_383[k]
                   + f_3 * pc_x[k] * sog_573[k];

        t_802[k] = f_17 * sng_452[k]
                   + f_3 * pc_y[k] * sog_572[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_x, pc_z, sng_438, sng_575, sng_576, sof0_385, \
                         sof0_386, sof1_385, sof1_386, sog_573, sog_575, \
                         sog_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_11 * sng_575[k]
                   + f_4 * sof0_385[k]
                   - f_5 * sof1_385[k]
                   + f_3 * pc_x[k] * sog_575[k];

        t_804[k] = f_11 * sng_576[k]
                   + f_6 * sof0_386[k]
                   - f_7 * sof1_386[k]
                   + f_3 * pc_x[k] * sog_576[k];

        t_805[k] = f_10 * sng_438[k]
                   + f_3 * pc_z[k] * sog_573[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, t_809, pc_x, pc_y, sng_455, sng_579, sng_580, \
                         sng_581, sof0_389, sof1_389, sog_575, sog_579, sog_580, \
                         sog_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_17 * sng_455[k]
                   + f_3 * pc_y[k] * sog_575[k];

        t_807[k] = f_11 * sng_579[k]
                   + f_6 * sof0_389[k]
                   - f_7 * sof1_389[k]
                   + f_3 * pc_x[k] * sog_579[k];

        t_808[k] = f_11 * sng_580[k]
                   + f_3 * pc_x[k] * sog_580[k];

        t_809[k] = f_11 * sng_581[k]
                   + f_3 * pc_x[k] * sog_581[k];
    }
}

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snh0,
                                                          const size_t sng, const size_t snh1,
                                                          const size_t sof0, const size_t sof1,
                                                          const size_t sog, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_810 = buffer.data(target + 810);
    auto *t_811 = buffer.data(target + 811);
    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);
    auto *t_815 = buffer.data(target + 815);
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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_735 = buffer.data(snh0 + 735);
    const auto *snh0_738 = buffer.data(snh0 + 738);
    const auto *snh0_740 = buffer.data(snh0 + 740);
    const auto *snh0_741 = buffer.data(snh0 + 741);
    const auto *snh0_744 = buffer.data(snh0 + 744);

    const auto *sng_445 = buffer.data(sng + 445);
    const auto *sng_449 = buffer.data(sng + 449);
    const auto *sng_450 = buffer.data(sng + 450);
    const auto *sng_453 = buffer.data(sng + 453);
    const auto *sng_460 = buffer.data(sng + 460);
    const auto *sng_462 = buffer.data(sng + 462);
    const auto *sng_463 = buffer.data(sng + 463);
    const auto *sng_464 = buffer.data(sng + 464);
    const auto *sng_465 = buffer.data(sng + 465);
    const auto *sng_467 = buffer.data(sng + 467);
    const auto *sng_468 = buffer.data(sng + 468);
    const auto *sng_470 = buffer.data(sng + 470);
    const auto *sng_475 = buffer.data(sng + 475);
    const auto *sng_477 = buffer.data(sng + 477);
    const auto *sng_478 = buffer.data(sng + 478);
    const auto *sng_479 = buffer.data(sng + 479);
    const auto *sng_480 = buffer.data(sng + 480);
    const auto *sng_482 = buffer.data(sng + 482);
    const auto *sng_483 = buffer.data(sng + 483);
    const auto *sng_485 = buffer.data(sng + 485);
    const auto *sng_490 = buffer.data(sng + 490);
    const auto *sng_492 = buffer.data(sng + 492);
    const auto *sng_493 = buffer.data(sng + 493);
    const auto *sng_494 = buffer.data(sng + 494);
    const auto *sng_495 = buffer.data(sng + 495);
    const auto *sng_497 = buffer.data(sng + 497);
    const auto *sng_498 = buffer.data(sng + 498);
    const auto *sng_500 = buffer.data(sng + 500);
    const auto *sng_505 = buffer.data(sng + 505);
    const auto *sng_507 = buffer.data(sng + 507);
    const auto *sng_508 = buffer.data(sng + 508);
    const auto *sng_509 = buffer.data(sng + 509);
    const auto *sng_510 = buffer.data(sng + 510);
    const auto *sng_512 = buffer.data(sng + 512);
    const auto *sng_513 = buffer.data(sng + 513);
    const auto *sng_515 = buffer.data(sng + 515);
    const auto *sng_520 = buffer.data(sng + 520);
    const auto *sng_522 = buffer.data(sng + 522);
    const auto *sng_523 = buffer.data(sng + 523);
    const auto *sng_524 = buffer.data(sng + 524);
    const auto *sng_525 = buffer.data(sng + 525);
    const auto *sng_526 = buffer.data(sng + 526);
    const auto *sng_527 = buffer.data(sng + 527);
    const auto *sng_528 = buffer.data(sng + 528);
    const auto *sng_530 = buffer.data(sng + 530);
    const auto *sng_535 = buffer.data(sng + 535);
    const auto *sng_537 = buffer.data(sng + 537);
    const auto *sng_582 = buffer.data(sng + 582);
    const auto *sng_583 = buffer.data(sng + 583);
    const auto *sng_584 = buffer.data(sng + 584);
    const auto *sng_585 = buffer.data(sng + 585);
    const auto *sng_588 = buffer.data(sng + 588);
    const auto *sng_590 = buffer.data(sng + 590);
    const auto *sng_591 = buffer.data(sng + 591);
    const auto *sng_594 = buffer.data(sng + 594);
    const auto *sng_595 = buffer.data(sng + 595);
    const auto *sng_596 = buffer.data(sng + 596);
    const auto *sng_597 = buffer.data(sng + 597);
    const auto *sng_598 = buffer.data(sng + 598);
    const auto *sng_599 = buffer.data(sng + 599);
    const auto *sng_600 = buffer.data(sng + 600);
    const auto *sng_603 = buffer.data(sng + 603);
    const auto *sng_605 = buffer.data(sng + 605);
    const auto *sng_606 = buffer.data(sng + 606);
    const auto *sng_609 = buffer.data(sng + 609);
    const auto *sng_610 = buffer.data(sng + 610);
    const auto *sng_611 = buffer.data(sng + 611);
    const auto *sng_612 = buffer.data(sng + 612);
    const auto *sng_613 = buffer.data(sng + 613);
    const auto *sng_614 = buffer.data(sng + 614);
    const auto *sng_615 = buffer.data(sng + 615);
    const auto *sng_618 = buffer.data(sng + 618);
    const auto *sng_620 = buffer.data(sng + 620);
    const auto *sng_621 = buffer.data(sng + 621);
    const auto *sng_624 = buffer.data(sng + 624);
    const auto *sng_625 = buffer.data(sng + 625);
    const auto *sng_626 = buffer.data(sng + 626);
    const auto *sng_627 = buffer.data(sng + 627);
    const auto *sng_628 = buffer.data(sng + 628);
    const auto *sng_629 = buffer.data(sng + 629);
    const auto *sng_630 = buffer.data(sng + 630);
    const auto *sng_633 = buffer.data(sng + 633);
    const auto *sng_635 = buffer.data(sng + 635);
    const auto *sng_636 = buffer.data(sng + 636);
    const auto *sng_639 = buffer.data(sng + 639);
    const auto *sng_640 = buffer.data(sng + 640);
    const auto *sng_641 = buffer.data(sng + 641);
    const auto *sng_642 = buffer.data(sng + 642);
    const auto *sng_643 = buffer.data(sng + 643);
    const auto *sng_644 = buffer.data(sng + 644);
    const auto *sng_655 = buffer.data(sng + 655);
    const auto *sng_656 = buffer.data(sng + 656);
    const auto *sng_657 = buffer.data(sng + 657);
    const auto *sng_658 = buffer.data(sng + 658);
    const auto *sng_659 = buffer.data(sng + 659);

    const auto *snh1_735 = buffer.data(snh1 + 735);
    const auto *snh1_738 = buffer.data(snh1 + 738);
    const auto *snh1_740 = buffer.data(snh1 + 740);
    const auto *snh1_741 = buffer.data(snh1 + 741);
    const auto *snh1_744 = buffer.data(snh1 + 744);

    const auto *sof0_386 = buffer.data(sof0 + 386);
    const auto *sof0_388 = buffer.data(sof0 + 388);
    const auto *sof0_389 = buffer.data(sof0 + 389);
    const auto *sof0_390 = buffer.data(sof0 + 390);
    const auto *sof0_393 = buffer.data(sof0 + 393);
    const auto *sof0_395 = buffer.data(sof0 + 395);
    const auto *sof0_396 = buffer.data(sof0 + 396);
    const auto *sof0_398 = buffer.data(sof0 + 398);
    const auto *sof0_399 = buffer.data(sof0 + 399);
    const auto *sof0_400 = buffer.data(sof0 + 400);
    const auto *sof0_403 = buffer.data(sof0 + 403);
    const auto *sof0_405 = buffer.data(sof0 + 405);
    const auto *sof0_406 = buffer.data(sof0 + 406);
    const auto *sof0_408 = buffer.data(sof0 + 408);
    const auto *sof0_409 = buffer.data(sof0 + 409);
    const auto *sof0_410 = buffer.data(sof0 + 410);
    const auto *sof0_413 = buffer.data(sof0 + 413);
    const auto *sof0_415 = buffer.data(sof0 + 415);
    const auto *sof0_416 = buffer.data(sof0 + 416);
    const auto *sof0_418 = buffer.data(sof0 + 418);
    const auto *sof0_419 = buffer.data(sof0 + 419);
    const auto *sof0_420 = buffer.data(sof0 + 420);
    const auto *sof0_423 = buffer.data(sof0 + 423);
    const auto *sof0_425 = buffer.data(sof0 + 425);
    const auto *sof0_426 = buffer.data(sof0 + 426);
    const auto *sof0_428 = buffer.data(sof0 + 428);
    const auto *sof0_429 = buffer.data(sof0 + 429);
    const auto *sof0_436 = buffer.data(sof0 + 436);
    const auto *sof0_438 = buffer.data(sof0 + 438);

    const auto *sof1_386 = buffer.data(sof1 + 386);
    const auto *sof1_388 = buffer.data(sof1 + 388);
    const auto *sof1_389 = buffer.data(sof1 + 389);
    const auto *sof1_390 = buffer.data(sof1 + 390);
    const auto *sof1_393 = buffer.data(sof1 + 393);
    const auto *sof1_395 = buffer.data(sof1 + 395);
    const auto *sof1_396 = buffer.data(sof1 + 396);
    const auto *sof1_398 = buffer.data(sof1 + 398);
    const auto *sof1_399 = buffer.data(sof1 + 399);
    const auto *sof1_400 = buffer.data(sof1 + 400);
    const auto *sof1_403 = buffer.data(sof1 + 403);
    const auto *sof1_405 = buffer.data(sof1 + 405);
    const auto *sof1_406 = buffer.data(sof1 + 406);
    const auto *sof1_408 = buffer.data(sof1 + 408);
    const auto *sof1_409 = buffer.data(sof1 + 409);
    const auto *sof1_410 = buffer.data(sof1 + 410);
    const auto *sof1_413 = buffer.data(sof1 + 413);
    const auto *sof1_415 = buffer.data(sof1 + 415);
    const auto *sof1_416 = buffer.data(sof1 + 416);
    const auto *sof1_418 = buffer.data(sof1 + 418);
    const auto *sof1_419 = buffer.data(sof1 + 419);
    const auto *sof1_420 = buffer.data(sof1 + 420);
    const auto *sof1_423 = buffer.data(sof1 + 423);
    const auto *sof1_425 = buffer.data(sof1 + 425);
    const auto *sof1_426 = buffer.data(sof1 + 426);
    const auto *sof1_428 = buffer.data(sof1 + 428);
    const auto *sof1_429 = buffer.data(sof1 + 429);
    const auto *sof1_436 = buffer.data(sof1 + 436);
    const auto *sof1_438 = buffer.data(sof1 + 438);

    const auto *sog_580 = buffer.data(sog + 580);
    const auto *sog_582 = buffer.data(sog + 582);
    const auto *sog_583 = buffer.data(sog + 583);
    const auto *sog_584 = buffer.data(sog + 584);
    const auto *sog_585 = buffer.data(sog + 585);
    const auto *sog_587 = buffer.data(sog + 587);
    const auto *sog_588 = buffer.data(sog + 588);
    const auto *sog_590 = buffer.data(sog + 590);
    const auto *sog_591 = buffer.data(sog + 591);
    const auto *sog_594 = buffer.data(sog + 594);
    const auto *sog_595 = buffer.data(sog + 595);
    const auto *sog_596 = buffer.data(sog + 596);
    const auto *sog_597 = buffer.data(sog + 597);
    const auto *sog_598 = buffer.data(sog + 598);
    const auto *sog_599 = buffer.data(sog + 599);
    const auto *sog_600 = buffer.data(sog + 600);
    const auto *sog_602 = buffer.data(sog + 602);
    const auto *sog_603 = buffer.data(sog + 603);
    const auto *sog_605 = buffer.data(sog + 605);
    const auto *sog_606 = buffer.data(sog + 606);
    const auto *sog_609 = buffer.data(sog + 609);
    const auto *sog_610 = buffer.data(sog + 610);
    const auto *sog_611 = buffer.data(sog + 611);
    const auto *sog_612 = buffer.data(sog + 612);
    const auto *sog_613 = buffer.data(sog + 613);
    const auto *sog_614 = buffer.data(sog + 614);
    const auto *sog_615 = buffer.data(sog + 615);
    const auto *sog_617 = buffer.data(sog + 617);
    const auto *sog_618 = buffer.data(sog + 618);
    const auto *sog_620 = buffer.data(sog + 620);
    const auto *sog_621 = buffer.data(sog + 621);
    const auto *sog_624 = buffer.data(sog + 624);
    const auto *sog_625 = buffer.data(sog + 625);
    const auto *sog_626 = buffer.data(sog + 626);
    const auto *sog_627 = buffer.data(sog + 627);
    const auto *sog_628 = buffer.data(sog + 628);
    const auto *sog_629 = buffer.data(sog + 629);
    const auto *sog_630 = buffer.data(sog + 630);
    const auto *sog_632 = buffer.data(sog + 632);
    const auto *sog_633 = buffer.data(sog + 633);
    const auto *sog_635 = buffer.data(sog + 635);
    const auto *sog_636 = buffer.data(sog + 636);
    const auto *sog_639 = buffer.data(sog + 639);
    const auto *sog_640 = buffer.data(sog + 640);
    const auto *sog_641 = buffer.data(sog + 641);
    const auto *sog_642 = buffer.data(sog + 642);
    const auto *sog_643 = buffer.data(sog + 643);
    const auto *sog_644 = buffer.data(sog + 644);
    const auto *sog_645 = buffer.data(sog + 645);
    const auto *sog_647 = buffer.data(sog + 647);
    const auto *sog_648 = buffer.data(sog + 648);
    const auto *sog_650 = buffer.data(sog + 650);
    const auto *sog_655 = buffer.data(sog + 655);
    const auto *sog_656 = buffer.data(sog + 656);
    const auto *sog_657 = buffer.data(sog + 657);
    const auto *sog_658 = buffer.data(sog + 658);
    const auto *sog_659 = buffer.data(sog + 659);

#pragma omp simd aligned(t_810, t_811, t_812, t_813, pc_x, pc_y, sng_460, sng_582, sng_583, \
                         sng_584, sof0_386, sof1_386, sog_580, sog_582, sog_583, \
                         sog_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_11 * sng_582[k]
                   + f_3 * pc_x[k] * sog_582[k];

        t_811[k] = f_11 * sng_583[k]
                   + f_3 * pc_x[k] * sog_583[k];

        t_812[k] = f_11 * sng_584[k]
                   + f_3 * pc_x[k] * sog_584[k];

        t_813[k] = f_17 * sng_460[k]
                   + f_1 * sof0_386[k]
                   - f_2 * sof1_386[k]
                   + f_3 * pc_y[k] * sog_580[k];
    }

#pragma omp simd aligned(t_814, t_815, t_816, pc_y, pc_z, sng_445, sng_462, sng_463, sof0_388, \
                         sof0_389, sof1_388, sof1_389, sog_580, sog_582, \
                         sog_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_814[k] = f_10 * sng_445[k]
                   + f_3 * pc_z[k] * sog_580[k];

        t_815[k] = f_17 * sng_462[k]
                   + f_4 * sof0_388[k]
                   - f_5 * sof1_388[k]
                   + f_3 * pc_y[k] * sog_582[k];

        t_816[k] = f_17 * sng_463[k]
                   + f_6 * sof0_389[k]
                   - f_7 * sof1_389[k]
                   + f_3 * pc_y[k] * sog_583[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, pc_x, pc_y, pc_z, sng_449, sng_464, sng_585, \
                         sof0_389, sof0_390, sof1_389, sof1_390, sog_584, \
                         sog_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = f_17 * sng_464[k]
                   + f_3 * pc_y[k] * sog_584[k];

        t_818[k] = f_10 * sng_449[k]
                   + f_1 * sof0_389[k]
                   - f_2 * sof1_389[k]
                   + f_3 * pc_z[k] * sog_584[k];

        t_819[k] = f_11 * sng_585[k]
                   + f_1 * sof0_390[k]
                   - f_2 * sof1_390[k]
                   + f_3 * pc_x[k] * sog_585[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, pc_x, pc_y, pc_z, sng_450, sng_465, \
                         sng_467, sng_588, sof0_393, sof1_393, sog_585, sog_587, \
                         sog_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = f_18 * sng_465[k]
                   + f_3 * pc_y[k] * sog_585[k];

        t_821[k] = f_11 * sng_450[k]
                   + f_3 * pc_z[k] * sog_585[k];

        t_822[k] = f_11 * sng_588[k]
                   + f_4 * sof0_393[k]
                   - f_5 * sof1_393[k]
                   + f_3 * pc_x[k] * sog_588[k];

        t_823[k] = f_18 * sng_467[k]
                   + f_3 * pc_y[k] * sog_587[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pc_x, pc_z, sng_453, sng_590, sng_591, sof0_395, \
                         sof0_396, sof1_395, sof1_396, sog_588, sog_590, \
                         sog_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_11 * sng_590[k]
                   + f_4 * sof0_395[k]
                   - f_5 * sof1_395[k]
                   + f_3 * pc_x[k] * sog_590[k];

        t_825[k] = f_11 * sng_591[k]
                   + f_6 * sof0_396[k]
                   - f_7 * sof1_396[k]
                   + f_3 * pc_x[k] * sog_591[k];

        t_826[k] = f_11 * sng_453[k]
                   + f_3 * pc_z[k] * sog_588[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, t_830, pc_x, pc_y, sng_470, sng_594, sng_595, \
                         sng_596, sof0_399, sof1_399, sog_590, sog_594, sog_595, \
                         sog_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_18 * sng_470[k]
                   + f_3 * pc_y[k] * sog_590[k];

        t_828[k] = f_11 * sng_594[k]
                   + f_6 * sof0_399[k]
                   - f_7 * sof1_399[k]
                   + f_3 * pc_x[k] * sog_594[k];

        t_829[k] = f_11 * sng_595[k]
                   + f_3 * pc_x[k] * sog_595[k];

        t_830[k] = f_11 * sng_596[k]
                   + f_3 * pc_x[k] * sog_596[k];
    }

#pragma omp simd aligned(t_831, t_832, t_833, t_834, pc_x, pc_y, sng_475, sng_597, sng_598, \
                         sng_599, sof0_396, sof1_396, sog_595, sog_597, sog_598, \
                         sog_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_831[k] = f_11 * sng_597[k]
                   + f_3 * pc_x[k] * sog_597[k];

        t_832[k] = f_11 * sng_598[k]
                   + f_3 * pc_x[k] * sog_598[k];

        t_833[k] = f_11 * sng_599[k]
                   + f_3 * pc_x[k] * sog_599[k];

        t_834[k] = f_18 * sng_475[k]
                   + f_1 * sof0_396[k]
                   - f_2 * sof1_396[k]
                   + f_3 * pc_y[k] * sog_595[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, pc_y, pc_z, sng_460, sng_477, sng_478, sof0_398, \
                         sof0_399, sof1_398, sof1_399, sog_595, sog_597, \
                         sog_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = f_11 * sng_460[k]
                   + f_3 * pc_z[k] * sog_595[k];

        t_836[k] = f_18 * sng_477[k]
                   + f_4 * sof0_398[k]
                   - f_5 * sof1_398[k]
                   + f_3 * pc_y[k] * sog_597[k];

        t_837[k] = f_18 * sng_478[k]
                   + f_6 * sof0_399[k]
                   - f_7 * sof1_399[k]
                   + f_3 * pc_y[k] * sog_598[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, pc_x, pc_y, pc_z, sng_464, sng_479, sng_600, \
                         sof0_399, sof0_400, sof1_399, sof1_400, sog_599, \
                         sog_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_18 * sng_479[k]
                   + f_3 * pc_y[k] * sog_599[k];

        t_839[k] = f_11 * sng_464[k]
                   + f_1 * sof0_399[k]
                   - f_2 * sof1_399[k]
                   + f_3 * pc_z[k] * sog_599[k];

        t_840[k] = f_11 * sng_600[k]
                   + f_1 * sof0_400[k]
                   - f_2 * sof1_400[k]
                   + f_3 * pc_x[k] * sog_600[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, t_844, pc_x, pc_y, pc_z, sng_465, sng_480, \
                         sng_482, sng_603, sof0_403, sof1_403, sog_600, sog_602, \
                         sog_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_16 * sng_480[k]
                   + f_3 * pc_y[k] * sog_600[k];

        t_842[k] = f_16 * sng_465[k]
                   + f_3 * pc_z[k] * sog_600[k];

        t_843[k] = f_11 * sng_603[k]
                   + f_4 * sof0_403[k]
                   - f_5 * sof1_403[k]
                   + f_3 * pc_x[k] * sog_603[k];

        t_844[k] = f_16 * sng_482[k]
                   + f_3 * pc_y[k] * sog_602[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_z, sng_468, sng_605, sng_606, sof0_405, \
                         sof0_406, sof1_405, sof1_406, sog_603, sog_605, \
                         sog_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_11 * sng_605[k]
                   + f_4 * sof0_405[k]
                   - f_5 * sof1_405[k]
                   + f_3 * pc_x[k] * sog_605[k];

        t_846[k] = f_11 * sng_606[k]
                   + f_6 * sof0_406[k]
                   - f_7 * sof1_406[k]
                   + f_3 * pc_x[k] * sog_606[k];

        t_847[k] = f_16 * sng_468[k]
                   + f_3 * pc_z[k] * sog_603[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, pc_x, pc_y, sng_485, sng_609, sng_610, \
                         sng_611, sof0_409, sof1_409, sog_605, sog_609, sog_610, \
                         sog_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_16 * sng_485[k]
                   + f_3 * pc_y[k] * sog_605[k];

        t_849[k] = f_11 * sng_609[k]
                   + f_6 * sof0_409[k]
                   - f_7 * sof1_409[k]
                   + f_3 * pc_x[k] * sog_609[k];

        t_850[k] = f_11 * sng_610[k]
                   + f_3 * pc_x[k] * sog_610[k];

        t_851[k] = f_11 * sng_611[k]
                   + f_3 * pc_x[k] * sog_611[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pc_x, pc_y, sng_490, sng_612, sng_613, \
                         sng_614, sof0_406, sof1_406, sog_610, sog_612, sog_613, \
                         sog_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_11 * sng_612[k]
                   + f_3 * pc_x[k] * sog_612[k];

        t_853[k] = f_11 * sng_613[k]
                   + f_3 * pc_x[k] * sog_613[k];

        t_854[k] = f_11 * sng_614[k]
                   + f_3 * pc_x[k] * sog_614[k];

        t_855[k] = f_16 * sng_490[k]
                   + f_1 * sof0_406[k]
                   - f_2 * sof1_406[k]
                   + f_3 * pc_y[k] * sog_610[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pc_y, pc_z, sng_475, sng_492, sng_493, sof0_408, \
                         sof0_409, sof1_408, sof1_409, sog_610, sog_612, \
                         sog_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_16 * sng_475[k]
                   + f_3 * pc_z[k] * sog_610[k];

        t_857[k] = f_16 * sng_492[k]
                   + f_4 * sof0_408[k]
                   - f_5 * sof1_408[k]
                   + f_3 * pc_y[k] * sog_612[k];

        t_858[k] = f_16 * sng_493[k]
                   + f_6 * sof0_409[k]
                   - f_7 * sof1_409[k]
                   + f_3 * pc_y[k] * sog_613[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pc_x, pc_y, pc_z, sng_479, sng_494, sng_615, \
                         sof0_409, sof0_410, sof1_409, sof1_410, sog_614, \
                         sog_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_16 * sng_494[k]
                   + f_3 * pc_y[k] * sog_614[k];

        t_860[k] = f_16 * sng_479[k]
                   + f_1 * sof0_409[k]
                   - f_2 * sof1_409[k]
                   + f_3 * pc_z[k] * sog_614[k];

        t_861[k] = f_11 * sng_615[k]
                   + f_1 * sof0_410[k]
                   - f_2 * sof1_410[k]
                   + f_3 * pc_x[k] * sog_615[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pc_x, pc_y, pc_z, sng_480, sng_495, \
                         sng_497, sng_618, sof0_413, sof1_413, sog_615, sog_617, \
                         sog_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_11 * sng_495[k]
                   + f_3 * pc_y[k] * sog_615[k];

        t_863[k] = f_18 * sng_480[k]
                   + f_3 * pc_z[k] * sog_615[k];

        t_864[k] = f_11 * sng_618[k]
                   + f_4 * sof0_413[k]
                   - f_5 * sof1_413[k]
                   + f_3 * pc_x[k] * sog_618[k];

        t_865[k] = f_11 * sng_497[k]
                   + f_3 * pc_y[k] * sog_617[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_z, sng_483, sng_620, sng_621, sof0_415, \
                         sof0_416, sof1_415, sof1_416, sog_618, sog_620, \
                         sog_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_11 * sng_620[k]
                   + f_4 * sof0_415[k]
                   - f_5 * sof1_415[k]
                   + f_3 * pc_x[k] * sog_620[k];

        t_867[k] = f_11 * sng_621[k]
                   + f_6 * sof0_416[k]
                   - f_7 * sof1_416[k]
                   + f_3 * pc_x[k] * sog_621[k];

        t_868[k] = f_18 * sng_483[k]
                   + f_3 * pc_z[k] * sog_618[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_y, sng_500, sng_624, sng_625, \
                         sng_626, sof0_419, sof1_419, sog_620, sog_624, sog_625, \
                         sog_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_11 * sng_500[k]
                   + f_3 * pc_y[k] * sog_620[k];

        t_870[k] = f_11 * sng_624[k]
                   + f_6 * sof0_419[k]
                   - f_7 * sof1_419[k]
                   + f_3 * pc_x[k] * sog_624[k];

        t_871[k] = f_11 * sng_625[k]
                   + f_3 * pc_x[k] * sog_625[k];

        t_872[k] = f_11 * sng_626[k]
                   + f_3 * pc_x[k] * sog_626[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, pc_x, pc_y, sng_505, sng_627, sng_628, \
                         sng_629, sof0_416, sof1_416, sog_625, sog_627, sog_628, \
                         sog_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_11 * sng_627[k]
                   + f_3 * pc_x[k] * sog_627[k];

        t_874[k] = f_11 * sng_628[k]
                   + f_3 * pc_x[k] * sog_628[k];

        t_875[k] = f_11 * sng_629[k]
                   + f_3 * pc_x[k] * sog_629[k];

        t_876[k] = f_11 * sng_505[k]
                   + f_1 * sof0_416[k]
                   - f_2 * sof1_416[k]
                   + f_3 * pc_y[k] * sog_625[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pc_y, pc_z, sng_490, sng_507, sng_508, sof0_418, \
                         sof0_419, sof1_418, sof1_419, sog_625, sog_627, \
                         sog_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_18 * sng_490[k]
                   + f_3 * pc_z[k] * sog_625[k];

        t_878[k] = f_11 * sng_507[k]
                   + f_4 * sof0_418[k]
                   - f_5 * sof1_418[k]
                   + f_3 * pc_y[k] * sog_627[k];

        t_879[k] = f_11 * sng_508[k]
                   + f_6 * sof0_419[k]
                   - f_7 * sof1_419[k]
                   + f_3 * pc_y[k] * sog_628[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, pc_x, pc_y, pc_z, sng_494, sng_509, sng_630, \
                         sof0_419, sof0_420, sof1_419, sof1_420, sog_629, \
                         sog_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = f_11 * sng_509[k]
                   + f_3 * pc_y[k] * sog_629[k];

        t_881[k] = f_18 * sng_494[k]
                   + f_1 * sof0_419[k]
                   - f_2 * sof1_419[k]
                   + f_3 * pc_z[k] * sog_629[k];

        t_882[k] = f_11 * sng_630[k]
                   + f_1 * sof0_420[k]
                   - f_2 * sof1_420[k]
                   + f_3 * pc_x[k] * sog_630[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, pc_x, pc_y, pc_z, sng_495, sng_510, \
                         sng_512, sng_633, sof0_423, sof1_423, sog_630, sog_632, \
                         sog_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = f_10 * sng_510[k]
                   + f_3 * pc_y[k] * sog_630[k];

        t_884[k] = f_17 * sng_495[k]
                   + f_3 * pc_z[k] * sog_630[k];

        t_885[k] = f_11 * sng_633[k]
                   + f_4 * sof0_423[k]
                   - f_5 * sof1_423[k]
                   + f_3 * pc_x[k] * sog_633[k];

        t_886[k] = f_10 * sng_512[k]
                   + f_3 * pc_y[k] * sog_632[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, pc_x, pc_z, sng_498, sng_635, sng_636, sof0_425, \
                         sof0_426, sof1_425, sof1_426, sog_633, sog_635, \
                         sog_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = f_11 * sng_635[k]
                   + f_4 * sof0_425[k]
                   - f_5 * sof1_425[k]
                   + f_3 * pc_x[k] * sog_635[k];

        t_888[k] = f_11 * sng_636[k]
                   + f_6 * sof0_426[k]
                   - f_7 * sof1_426[k]
                   + f_3 * pc_x[k] * sog_636[k];

        t_889[k] = f_17 * sng_498[k]
                   + f_3 * pc_z[k] * sog_633[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, pc_x, pc_y, sng_515, sng_639, sng_640, \
                         sng_641, sof0_429, sof1_429, sog_635, sog_639, sog_640, \
                         sog_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_10 * sng_515[k]
                   + f_3 * pc_y[k] * sog_635[k];

        t_891[k] = f_11 * sng_639[k]
                   + f_6 * sof0_429[k]
                   - f_7 * sof1_429[k]
                   + f_3 * pc_x[k] * sog_639[k];

        t_892[k] = f_11 * sng_640[k]
                   + f_3 * pc_x[k] * sog_640[k];

        t_893[k] = f_11 * sng_641[k]
                   + f_3 * pc_x[k] * sog_641[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, pc_x, pc_y, sng_520, sng_642, sng_643, \
                         sng_644, sof0_426, sof1_426, sog_640, sog_642, sog_643, \
                         sog_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_11 * sng_642[k]
                   + f_3 * pc_x[k] * sog_642[k];

        t_895[k] = f_11 * sng_643[k]
                   + f_3 * pc_x[k] * sog_643[k];

        t_896[k] = f_11 * sng_644[k]
                   + f_3 * pc_x[k] * sog_644[k];

        t_897[k] = f_10 * sng_520[k]
                   + f_1 * sof0_426[k]
                   - f_2 * sof1_426[k]
                   + f_3 * pc_y[k] * sog_640[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, pc_y, pc_z, sng_505, sng_522, sng_523, sof0_428, \
                         sof0_429, sof1_428, sof1_429, sog_640, sog_642, \
                         sog_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_17 * sng_505[k]
                   + f_3 * pc_z[k] * sog_640[k];

        t_899[k] = f_10 * sng_522[k]
                   + f_4 * sof0_428[k]
                   - f_5 * sof1_428[k]
                   + f_3 * pc_y[k] * sog_642[k];

        t_900[k] = f_10 * sng_523[k]
                   + f_6 * sof0_429[k]
                   - f_7 * sof1_429[k]
                   + f_3 * pc_y[k] * sog_643[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, t_904, pb_y, pc_y, pc_z, snh0_735, sng_509, \
                         sng_524, sng_525, snh1_735, sof0_429, sof1_429, sog_644, \
                         sog_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = f_10 * sng_524[k]
                   + f_3 * pc_y[k] * sog_644[k];

        t_902[k] = f_17 * sng_509[k]
                   + f_1 * sof0_429[k]
                   - f_2 * sof1_429[k]
                   + f_3 * pc_z[k] * sog_644[k];

        t_903[k] = pb_y[k] * snh0_735[k]
                   - f_8 * pc_y[k] * snh1_735[k];

        t_904[k] = f_9 * sng_525[k]
                   + f_3 * pc_y[k] * sog_645[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, pb_y, pc_y, pc_z, snh0_738, snh0_740, \
                         sng_510, sng_526, sng_527, snh1_738, snh1_740, sog_645, \
                         sog_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_15 * sng_510[k]
                   + f_3 * pc_z[k] * sog_645[k];

        t_906[k] = pb_y[k] * snh0_738[k]
                   + f_10 * sng_526[k]
                   - f_8 * pc_y[k] * snh1_738[k];

        t_907[k] = f_9 * sng_527[k]
                   + f_3 * pc_y[k] * sog_647[k];

        t_908[k] = pb_y[k] * snh0_740[k]
                   - f_8 * pc_y[k] * snh1_740[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, t_912, pb_y, pc_y, pc_z, snh0_741, snh0_744, \
                         sng_513, sng_528, sng_530, snh1_741, snh1_744, sog_648, \
                         sog_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = pb_y[k] * snh0_741[k]
                   + f_11 * sng_528[k]
                   - f_8 * pc_y[k] * snh1_741[k];

        t_910[k] = f_15 * sng_513[k]
                   + f_3 * pc_z[k] * sog_648[k];

        t_911[k] = f_9 * sng_530[k]
                   + f_3 * pc_y[k] * sog_650[k];

        t_912[k] = pb_y[k] * snh0_744[k]
                   - f_8 * pc_y[k] * snh1_744[k];
    }

#pragma omp simd aligned(t_913, t_914, t_915, t_916, t_917, pc_x, sng_655, sng_656, sng_657, \
                         sng_658, sng_659, sog_655, sog_656, sog_657, sog_658, \
                         sog_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_913[k] = f_11 * sng_655[k]
                   + f_3 * pc_x[k] * sog_655[k];

        t_914[k] = f_11 * sng_656[k]
                   + f_3 * pc_x[k] * sog_656[k];

        t_915[k] = f_11 * sng_657[k]
                   + f_3 * pc_x[k] * sog_657[k];

        t_916[k] = f_11 * sng_658[k]
                   + f_3 * pc_x[k] * sog_658[k];

        t_917[k] = f_11 * sng_659[k]
                   + f_3 * pc_x[k] * sog_659[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pc_y, pc_z, sng_520, sng_535, sng_537, sof0_436, \
                         sof0_438, sof1_436, sof1_438, sog_655, \
                         sog_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_9 * sng_535[k]
                   + f_1 * sof0_436[k]
                   - f_2 * sof1_436[k]
                   + f_3 * pc_y[k] * sog_655[k];

        t_919[k] = f_15 * sng_520[k]
                   + f_3 * pc_z[k] * sog_655[k];

        t_920[k] = f_9 * sng_537[k]
                   + f_4 * sof0_438[k]
                   - f_5 * sof1_438[k]
                   + f_3 * pc_y[k] * sog_657[k];
    }
}

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snh0,
                                                          const size_t sng, const size_t snh1,
                                                          const size_t sof0, const size_t sof1,
                                                          const size_t sog, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 4.5 / q;
    const auto f_14 = 4.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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
    auto *t_939 = buffer.data(target + 939);
    auto *t_940 = buffer.data(target + 940);
    auto *t_941 = buffer.data(target + 941);
    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);
    auto *t_945 = buffer.data(target + 945);
    auto *t_946 = buffer.data(target + 946);
    auto *t_947 = buffer.data(target + 947);
    auto *t_948 = buffer.data(target + 948);
    auto *t_949 = buffer.data(target + 949);
    auto *t_950 = buffer.data(target + 950);
    auto *t_951 = buffer.data(target + 951);
    auto *t_952 = buffer.data(target + 952);
    auto *t_953 = buffer.data(target + 953);
    auto *t_954 = buffer.data(target + 954);
    auto *t_955 = buffer.data(target + 955);
    auto *t_956 = buffer.data(target + 956);
    auto *t_957 = buffer.data(target + 957);
    auto *t_958 = buffer.data(target + 958);
    auto *t_959 = buffer.data(target + 959);
    auto *t_960 = buffer.data(target + 960);
    auto *t_961 = buffer.data(target + 961);
    auto *t_962 = buffer.data(target + 962);
    auto *t_963 = buffer.data(target + 963);
    auto *t_964 = buffer.data(target + 964);
    auto *t_965 = buffer.data(target + 965);
    auto *t_966 = buffer.data(target + 966);
    auto *t_967 = buffer.data(target + 967);
    auto *t_968 = buffer.data(target + 968);
    auto *t_969 = buffer.data(target + 969);
    auto *t_970 = buffer.data(target + 970);
    auto *t_971 = buffer.data(target + 971);
    auto *t_972 = buffer.data(target + 972);
    auto *t_973 = buffer.data(target + 973);
    auto *t_974 = buffer.data(target + 974);
    auto *t_975 = buffer.data(target + 975);
    auto *t_976 = buffer.data(target + 976);
    auto *t_977 = buffer.data(target + 977);
    auto *t_978 = buffer.data(target + 978);
    auto *t_979 = buffer.data(target + 979);
    auto *t_980 = buffer.data(target + 980);
    auto *t_981 = buffer.data(target + 981);
    auto *t_982 = buffer.data(target + 982);
    auto *t_983 = buffer.data(target + 983);
    auto *t_984 = buffer.data(target + 984);
    auto *t_985 = buffer.data(target + 985);
    auto *t_986 = buffer.data(target + 986);
    auto *t_987 = buffer.data(target + 987);
    auto *t_988 = buffer.data(target + 988);
    auto *t_989 = buffer.data(target + 989);
    auto *t_990 = buffer.data(target + 990);
    auto *t_991 = buffer.data(target + 991);
    auto *t_992 = buffer.data(target + 992);
    auto *t_993 = buffer.data(target + 993);
    auto *t_994 = buffer.data(target + 994);
    auto *t_995 = buffer.data(target + 995);
    auto *t_996 = buffer.data(target + 996);
    auto *t_997 = buffer.data(target + 997);
    auto *t_998 = buffer.data(target + 998);
    auto *t_999 = buffer.data(target + 999);
    auto *t_1000 = buffer.data(target + 1000);
    auto *t_1001 = buffer.data(target + 1001);
    auto *t_1002 = buffer.data(target + 1002);
    auto *t_1003 = buffer.data(target + 1003);
    auto *t_1004 = buffer.data(target + 1004);
    auto *t_1005 = buffer.data(target + 1005);
    auto *t_1006 = buffer.data(target + 1006);
    auto *t_1007 = buffer.data(target + 1007);
    auto *t_1008 = buffer.data(target + 1008);
    auto *t_1009 = buffer.data(target + 1009);
    auto *t_1010 = buffer.data(target + 1010);
    auto *t_1011 = buffer.data(target + 1011);
    auto *t_1012 = buffer.data(target + 1012);
    auto *t_1013 = buffer.data(target + 1013);
    auto *t_1014 = buffer.data(target + 1014);
    auto *t_1015 = buffer.data(target + 1015);
    auto *t_1016 = buffer.data(target + 1016);
    auto *t_1017 = buffer.data(target + 1017);
    auto *t_1018 = buffer.data(target + 1018);
    auto *t_1019 = buffer.data(target + 1019);
    auto *t_1020 = buffer.data(target + 1020);
    auto *t_1021 = buffer.data(target + 1021);
    auto *t_1022 = buffer.data(target + 1022);
    auto *t_1023 = buffer.data(target + 1023);
    auto *t_1024 = buffer.data(target + 1024);
    auto *t_1025 = buffer.data(target + 1025);
    auto *t_1026 = buffer.data(target + 1026);
    auto *t_1027 = buffer.data(target + 1027);
    auto *t_1028 = buffer.data(target + 1028);
    auto *t_1029 = buffer.data(target + 1029);
    auto *t_1030 = buffer.data(target + 1030);
    auto *t_1031 = buffer.data(target + 1031);
    auto *t_1032 = buffer.data(target + 1032);
    auto *t_1033 = buffer.data(target + 1033);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_755 = buffer.data(snh0 + 755);
    const auto *snh0_756 = buffer.data(snh0 + 756);
    const auto *snh0_759 = buffer.data(snh0 + 759);
    const auto *snh0_762 = buffer.data(snh0 + 762);
    const auto *snh0_771 = buffer.data(snh0 + 771);

    const auto *sng_525 = buffer.data(sng + 525);
    const auto *sng_528 = buffer.data(sng + 528);
    const auto *sng_535 = buffer.data(sng + 535);
    const auto *sng_538 = buffer.data(sng + 538);
    const auto *sng_539 = buffer.data(sng + 539);
    const auto *sng_540 = buffer.data(sng + 540);
    const auto *sng_542 = buffer.data(sng + 542);
    const auto *sng_543 = buffer.data(sng + 543);
    const auto *sng_545 = buffer.data(sng + 545);
    const auto *sng_550 = buffer.data(sng + 550);
    const auto *sng_552 = buffer.data(sng + 552);
    const auto *sng_553 = buffer.data(sng + 553);
    const auto *sng_554 = buffer.data(sng + 554);
    const auto *sng_555 = buffer.data(sng + 555);
    const auto *sng_557 = buffer.data(sng + 557);
    const auto *sng_558 = buffer.data(sng + 558);
    const auto *sng_560 = buffer.data(sng + 560);
    const auto *sng_565 = buffer.data(sng + 565);
    const auto *sng_567 = buffer.data(sng + 567);
    const auto *sng_568 = buffer.data(sng + 568);
    const auto *sng_569 = buffer.data(sng + 569);
    const auto *sng_570 = buffer.data(sng + 570);
    const auto *sng_572 = buffer.data(sng + 572);
    const auto *sng_573 = buffer.data(sng + 573);
    const auto *sng_575 = buffer.data(sng + 575);
    const auto *sng_580 = buffer.data(sng + 580);
    const auto *sng_582 = buffer.data(sng + 582);
    const auto *sng_583 = buffer.data(sng + 583);
    const auto *sng_584 = buffer.data(sng + 584);
    const auto *sng_585 = buffer.data(sng + 585);
    const auto *sng_587 = buffer.data(sng + 587);
    const auto *sng_590 = buffer.data(sng + 590);
    const auto *sng_595 = buffer.data(sng + 595);
    const auto *sng_597 = buffer.data(sng + 597);
    const auto *sng_598 = buffer.data(sng + 598);
    const auto *sng_599 = buffer.data(sng + 599);
    const auto *sng_600 = buffer.data(sng + 600);
    const auto *sng_602 = buffer.data(sng + 602);
    const auto *sng_660 = buffer.data(sng + 660);
    const auto *sng_663 = buffer.data(sng + 663);
    const auto *sng_665 = buffer.data(sng + 665);
    const auto *sng_666 = buffer.data(sng + 666);
    const auto *sng_669 = buffer.data(sng + 669);
    const auto *sng_670 = buffer.data(sng + 670);
    const auto *sng_671 = buffer.data(sng + 671);
    const auto *sng_672 = buffer.data(sng + 672);
    const auto *sng_673 = buffer.data(sng + 673);
    const auto *sng_674 = buffer.data(sng + 674);
    const auto *sng_675 = buffer.data(sng + 675);
    const auto *sng_678 = buffer.data(sng + 678);
    const auto *sng_680 = buffer.data(sng + 680);
    const auto *sng_681 = buffer.data(sng + 681);
    const auto *sng_684 = buffer.data(sng + 684);
    const auto *sng_685 = buffer.data(sng + 685);
    const auto *sng_686 = buffer.data(sng + 686);
    const auto *sng_687 = buffer.data(sng + 687);
    const auto *sng_688 = buffer.data(sng + 688);
    const auto *sng_689 = buffer.data(sng + 689);
    const auto *sng_695 = buffer.data(sng + 695);
    const auto *sng_699 = buffer.data(sng + 699);
    const auto *sng_700 = buffer.data(sng + 700);
    const auto *sng_701 = buffer.data(sng + 701);
    const auto *sng_702 = buffer.data(sng + 702);
    const auto *sng_703 = buffer.data(sng + 703);
    const auto *sng_704 = buffer.data(sng + 704);
    const auto *sng_705 = buffer.data(sng + 705);
    const auto *sng_708 = buffer.data(sng + 708);
    const auto *sng_710 = buffer.data(sng + 710);
    const auto *sng_711 = buffer.data(sng + 711);
    const auto *sng_714 = buffer.data(sng + 714);
    const auto *sng_715 = buffer.data(sng + 715);
    const auto *sng_716 = buffer.data(sng + 716);
    const auto *sng_717 = buffer.data(sng + 717);
    const auto *sng_718 = buffer.data(sng + 718);
    const auto *sng_719 = buffer.data(sng + 719);
    const auto *sng_720 = buffer.data(sng + 720);
    const auto *sng_723 = buffer.data(sng + 723);
    const auto *sng_725 = buffer.data(sng + 725);
    const auto *sng_726 = buffer.data(sng + 726);
    const auto *sng_729 = buffer.data(sng + 729);
    const auto *sng_730 = buffer.data(sng + 730);
    const auto *sng_731 = buffer.data(sng + 731);
    const auto *sng_732 = buffer.data(sng + 732);
    const auto *sng_733 = buffer.data(sng + 733);
    const auto *sng_734 = buffer.data(sng + 734);
    const auto *sng_735 = buffer.data(sng + 735);
    const auto *sng_738 = buffer.data(sng + 738);

    const auto *snh1_755 = buffer.data(snh1 + 755);
    const auto *snh1_756 = buffer.data(snh1 + 756);
    const auto *snh1_759 = buffer.data(snh1 + 759);
    const auto *snh1_762 = buffer.data(snh1 + 762);
    const auto *snh1_771 = buffer.data(snh1 + 771);

    const auto *sof0_439 = buffer.data(sof0 + 439);
    const auto *sof0_440 = buffer.data(sof0 + 440);
    const auto *sof0_443 = buffer.data(sof0 + 443);
    const auto *sof0_445 = buffer.data(sof0 + 445);
    const auto *sof0_446 = buffer.data(sof0 + 446);
    const auto *sof0_448 = buffer.data(sof0 + 448);
    const auto *sof0_449 = buffer.data(sof0 + 449);
    const auto *sof0_450 = buffer.data(sof0 + 450);
    const auto *sof0_453 = buffer.data(sof0 + 453);
    const auto *sof0_455 = buffer.data(sof0 + 455);
    const auto *sof0_456 = buffer.data(sof0 + 456);
    const auto *sof0_458 = buffer.data(sof0 + 458);
    const auto *sof0_459 = buffer.data(sof0 + 459);
    const auto *sof0_465 = buffer.data(sof0 + 465);
    const auto *sof0_468 = buffer.data(sof0 + 468);
    const auto *sof0_469 = buffer.data(sof0 + 469);
    const auto *sof0_470 = buffer.data(sof0 + 470);
    const auto *sof0_473 = buffer.data(sof0 + 473);
    const auto *sof0_475 = buffer.data(sof0 + 475);
    const auto *sof0_476 = buffer.data(sof0 + 476);
    const auto *sof0_478 = buffer.data(sof0 + 478);
    const auto *sof0_479 = buffer.data(sof0 + 479);
    const auto *sof0_480 = buffer.data(sof0 + 480);
    const auto *sof0_483 = buffer.data(sof0 + 483);
    const auto *sof0_485 = buffer.data(sof0 + 485);
    const auto *sof0_486 = buffer.data(sof0 + 486);
    const auto *sof0_488 = buffer.data(sof0 + 488);
    const auto *sof0_489 = buffer.data(sof0 + 489);
    const auto *sof0_490 = buffer.data(sof0 + 490);
    const auto *sof0_493 = buffer.data(sof0 + 493);

    const auto *sof1_439 = buffer.data(sof1 + 439);
    const auto *sof1_440 = buffer.data(sof1 + 440);
    const auto *sof1_443 = buffer.data(sof1 + 443);
    const auto *sof1_445 = buffer.data(sof1 + 445);
    const auto *sof1_446 = buffer.data(sof1 + 446);
    const auto *sof1_448 = buffer.data(sof1 + 448);
    const auto *sof1_449 = buffer.data(sof1 + 449);
    const auto *sof1_450 = buffer.data(sof1 + 450);
    const auto *sof1_453 = buffer.data(sof1 + 453);
    const auto *sof1_455 = buffer.data(sof1 + 455);
    const auto *sof1_456 = buffer.data(sof1 + 456);
    const auto *sof1_458 = buffer.data(sof1 + 458);
    const auto *sof1_459 = buffer.data(sof1 + 459);
    const auto *sof1_465 = buffer.data(sof1 + 465);
    const auto *sof1_468 = buffer.data(sof1 + 468);
    const auto *sof1_469 = buffer.data(sof1 + 469);
    const auto *sof1_470 = buffer.data(sof1 + 470);
    const auto *sof1_473 = buffer.data(sof1 + 473);
    const auto *sof1_475 = buffer.data(sof1 + 475);
    const auto *sof1_476 = buffer.data(sof1 + 476);
    const auto *sof1_478 = buffer.data(sof1 + 478);
    const auto *sof1_479 = buffer.data(sof1 + 479);
    const auto *sof1_480 = buffer.data(sof1 + 480);
    const auto *sof1_483 = buffer.data(sof1 + 483);
    const auto *sof1_485 = buffer.data(sof1 + 485);
    const auto *sof1_486 = buffer.data(sof1 + 486);
    const auto *sof1_488 = buffer.data(sof1 + 488);
    const auto *sof1_489 = buffer.data(sof1 + 489);
    const auto *sof1_490 = buffer.data(sof1 + 490);
    const auto *sof1_493 = buffer.data(sof1 + 493);

    const auto *sog_658 = buffer.data(sog + 658);
    const auto *sog_659 = buffer.data(sog + 659);
    const auto *sog_660 = buffer.data(sog + 660);
    const auto *sog_662 = buffer.data(sog + 662);
    const auto *sog_663 = buffer.data(sog + 663);
    const auto *sog_665 = buffer.data(sog + 665);
    const auto *sog_666 = buffer.data(sog + 666);
    const auto *sog_669 = buffer.data(sog + 669);
    const auto *sog_670 = buffer.data(sog + 670);
    const auto *sog_671 = buffer.data(sog + 671);
    const auto *sog_672 = buffer.data(sog + 672);
    const auto *sog_673 = buffer.data(sog + 673);
    const auto *sog_674 = buffer.data(sog + 674);
    const auto *sog_675 = buffer.data(sog + 675);
    const auto *sog_677 = buffer.data(sog + 677);
    const auto *sog_678 = buffer.data(sog + 678);
    const auto *sog_680 = buffer.data(sog + 680);
    const auto *sog_681 = buffer.data(sog + 681);
    const auto *sog_684 = buffer.data(sog + 684);
    const auto *sog_685 = buffer.data(sog + 685);
    const auto *sog_686 = buffer.data(sog + 686);
    const auto *sog_687 = buffer.data(sog + 687);
    const auto *sog_688 = buffer.data(sog + 688);
    const auto *sog_689 = buffer.data(sog + 689);
    const auto *sog_690 = buffer.data(sog + 690);
    const auto *sog_692 = buffer.data(sog + 692);
    const auto *sog_693 = buffer.data(sog + 693);
    const auto *sog_695 = buffer.data(sog + 695);
    const auto *sog_699 = buffer.data(sog + 699);
    const auto *sog_700 = buffer.data(sog + 700);
    const auto *sog_701 = buffer.data(sog + 701);
    const auto *sog_702 = buffer.data(sog + 702);
    const auto *sog_703 = buffer.data(sog + 703);
    const auto *sog_704 = buffer.data(sog + 704);
    const auto *sog_705 = buffer.data(sog + 705);
    const auto *sog_707 = buffer.data(sog + 707);
    const auto *sog_708 = buffer.data(sog + 708);
    const auto *sog_710 = buffer.data(sog + 710);
    const auto *sog_711 = buffer.data(sog + 711);
    const auto *sog_714 = buffer.data(sog + 714);
    const auto *sog_715 = buffer.data(sog + 715);
    const auto *sog_716 = buffer.data(sog + 716);
    const auto *sog_717 = buffer.data(sog + 717);
    const auto *sog_718 = buffer.data(sog + 718);
    const auto *sog_719 = buffer.data(sog + 719);
    const auto *sog_720 = buffer.data(sog + 720);
    const auto *sog_722 = buffer.data(sog + 722);
    const auto *sog_723 = buffer.data(sog + 723);
    const auto *sog_725 = buffer.data(sog + 725);
    const auto *sog_726 = buffer.data(sog + 726);
    const auto *sog_729 = buffer.data(sog + 729);
    const auto *sog_730 = buffer.data(sog + 730);
    const auto *sog_731 = buffer.data(sog + 731);
    const auto *sog_732 = buffer.data(sog + 732);
    const auto *sog_733 = buffer.data(sog + 733);
    const auto *sog_734 = buffer.data(sog + 734);
    const auto *sog_735 = buffer.data(sog + 735);
    const auto *sog_737 = buffer.data(sog + 737);
    const auto *sog_738 = buffer.data(sog + 738);

#pragma omp simd aligned(t_921, t_922, t_923, pb_y, pc_y, snh0_755, sng_538, sng_539, \
                         snh1_755, sof0_439, sof1_439, sog_658, \
                         sog_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_9 * sng_538[k]
                   + f_6 * sof0_439[k]
                   - f_7 * sof1_439[k]
                   + f_3 * pc_y[k] * sog_658[k];

        t_922[k] = f_9 * sng_539[k]
                   + f_3 * pc_y[k] * sog_659[k];

        t_923[k] = pb_y[k] * snh0_755[k]
                   - f_8 * pc_y[k] * snh1_755[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pc_x, pc_y, pc_z, sng_525, sng_660, \
                         sng_663, sof0_440, sof0_443, sof1_440, sof1_443, sog_660, \
                         sog_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_11 * sng_660[k]
                   + f_1 * sof0_440[k]
                   - f_2 * sof1_440[k]
                   + f_3 * pc_x[k] * sog_660[k];

        t_925[k] = f_3 * pc_y[k] * sog_660[k];

        t_926[k] = f_14 * sng_525[k]
                   + f_3 * pc_z[k] * sog_660[k];

        t_927[k] = f_11 * sng_663[k]
                   + f_4 * sof0_443[k]
                   - f_5 * sof1_443[k]
                   + f_3 * pc_x[k] * sog_663[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, pc_x, pc_y, sng_665, sng_666, sof0_445, \
                         sof0_446, sof1_445, sof1_446, sog_662, sog_665, \
                         sog_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_3 * pc_y[k] * sog_662[k];

        t_929[k] = f_11 * sng_665[k]
                   + f_4 * sof0_445[k]
                   - f_5 * sof1_445[k]
                   + f_3 * pc_x[k] * sog_665[k];

        t_930[k] = f_11 * sng_666[k]
                   + f_6 * sof0_446[k]
                   - f_7 * sof1_446[k]
                   + f_3 * pc_x[k] * sog_666[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, pc_x, pc_y, pc_z, sng_528, sng_669, \
                         sng_670, sof0_449, sof1_449, sog_663, sog_665, sog_669, \
                         sog_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_14 * sng_528[k]
                   + f_3 * pc_z[k] * sog_663[k];

        t_932[k] = f_3 * pc_y[k] * sog_665[k];

        t_933[k] = f_11 * sng_669[k]
                   + f_6 * sof0_449[k]
                   - f_7 * sof1_449[k]
                   + f_3 * pc_x[k] * sog_669[k];

        t_934[k] = f_11 * sng_670[k]
                   + f_3 * pc_x[k] * sog_670[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, pc_x, sng_671, sng_672, sng_673, sng_674, \
                         sog_671, sog_672, sog_673, sog_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_11 * sng_671[k]
                   + f_3 * pc_x[k] * sog_671[k];

        t_936[k] = f_11 * sng_672[k]
                   + f_3 * pc_x[k] * sog_672[k];

        t_937[k] = f_11 * sng_673[k]
                   + f_3 * pc_x[k] * sog_673[k];

        t_938[k] = f_11 * sng_674[k]
                   + f_3 * pc_x[k] * sog_674[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, t_942, pc_y, pc_z, sng_535, sof0_446, sof0_448, \
                         sof0_449, sof1_446, sof1_448, sof1_449, sog_670, sog_672, \
                         sog_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_1 * sof0_446[k]
                   - f_2 * sof1_446[k]
                   + f_3 * pc_y[k] * sog_670[k];

        t_940[k] = f_14 * sng_535[k]
                   + f_3 * pc_z[k] * sog_670[k];

        t_941[k] = f_4 * sof0_448[k]
                   - f_5 * sof1_448[k]
                   + f_3 * pc_y[k] * sog_672[k];

        t_942[k] = f_6 * sof0_449[k]
                   - f_7 * sof1_449[k]
                   + f_3 * pc_y[k] * sog_673[k];
    }

#pragma omp simd aligned(t_943, t_944, t_945, t_946, pc_x, pc_y, pc_z, sng_539, sng_540, \
                         sng_675, sof0_449, sof0_450, sof1_449, sof1_450, sog_674, \
                         sog_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = f_3 * pc_y[k] * sog_674[k];

        t_944[k] = f_14 * sng_539[k]
                   + f_1 * sof0_449[k]
                   - f_2 * sof1_449[k]
                   + f_3 * pc_z[k] * sog_674[k];

        t_945[k] = f_10 * sng_675[k]
                   + f_1 * sof0_450[k]
                   - f_2 * sof1_450[k]
                   + f_3 * pc_x[k] * sog_675[k];

        t_946[k] = f_13 * sng_540[k]
                   + f_3 * pc_y[k] * sog_675[k];
    }

#pragma omp simd aligned(t_947, t_948, t_949, pc_x, pc_y, pc_z, sng_542, sng_678, sof0_453, \
                         sof1_453, sog_675, sog_677, sog_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_947[k] = f_3 * pc_z[k] * sog_675[k];

        t_948[k] = f_10 * sng_678[k]
                   + f_4 * sof0_453[k]
                   - f_5 * sof1_453[k]
                   + f_3 * pc_x[k] * sog_678[k];

        t_949[k] = f_13 * sng_542[k]
                   + f_3 * pc_y[k] * sog_677[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, pc_x, pc_z, sng_680, sng_681, sof0_455, \
                         sof0_456, sof1_455, sof1_456, sog_678, sog_680, \
                         sog_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = f_10 * sng_680[k]
                   + f_4 * sof0_455[k]
                   - f_5 * sof1_455[k]
                   + f_3 * pc_x[k] * sog_680[k];

        t_951[k] = f_10 * sng_681[k]
                   + f_6 * sof0_456[k]
                   - f_7 * sof1_456[k]
                   + f_3 * pc_x[k] * sog_681[k];

        t_952[k] = f_3 * pc_z[k] * sog_678[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pc_x, pc_y, sng_545, sng_684, sng_685, \
                         sng_686, sof0_459, sof1_459, sog_680, sog_684, sog_685, \
                         sog_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_13 * sng_545[k]
                   + f_3 * pc_y[k] * sog_680[k];

        t_954[k] = f_10 * sng_684[k]
                   + f_6 * sof0_459[k]
                   - f_7 * sof1_459[k]
                   + f_3 * pc_x[k] * sog_684[k];

        t_955[k] = f_10 * sng_685[k]
                   + f_3 * pc_x[k] * sog_685[k];

        t_956[k] = f_10 * sng_686[k]
                   + f_3 * pc_x[k] * sog_686[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pc_x, pc_y, sng_550, sng_687, sng_688, \
                         sng_689, sof0_456, sof1_456, sog_685, sog_687, sog_688, \
                         sog_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_10 * sng_687[k]
                   + f_3 * pc_x[k] * sog_687[k];

        t_958[k] = f_10 * sng_688[k]
                   + f_3 * pc_x[k] * sog_688[k];

        t_959[k] = f_10 * sng_689[k]
                   + f_3 * pc_x[k] * sog_689[k];

        t_960[k] = f_13 * sng_550[k]
                   + f_1 * sof0_456[k]
                   - f_2 * sof1_456[k]
                   + f_3 * pc_y[k] * sog_685[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, pc_y, pc_z, sng_552, sng_553, sof0_458, \
                         sof0_459, sof1_458, sof1_459, sog_685, sog_687, \
                         sog_688 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = f_3 * pc_z[k] * sog_685[k];

        t_962[k] = f_13 * sng_552[k]
                   + f_4 * sof0_458[k]
                   - f_5 * sof1_458[k]
                   + f_3 * pc_y[k] * sog_687[k];

        t_963[k] = f_13 * sng_553[k]
                   + f_6 * sof0_459[k]
                   - f_7 * sof1_459[k]
                   + f_3 * pc_y[k] * sog_688[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pb_z, pc_y, pc_z, snh0_756, sng_554, \
                         sng_555, snh1_756, sof0_459, sof1_459, sog_689, \
                         sog_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = f_13 * sng_554[k]
                   + f_3 * pc_y[k] * sog_689[k];

        t_965[k] = f_1 * sof0_459[k]
                   - f_2 * sof1_459[k]
                   + f_3 * pc_z[k] * sog_689[k];

        t_966[k] = pb_z[k] * snh0_756[k]
                   - f_8 * pc_z[k] * snh1_756[k];

        t_967[k] = f_14 * sng_555[k]
                   + f_3 * pc_y[k] * sog_690[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, pb_z, pc_y, pc_z, snh0_759, sng_540, sng_557, \
                         snh1_759, sog_690, sog_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_9 * sng_540[k]
                   + f_3 * pc_z[k] * sog_690[k];

        t_969[k] = pb_z[k] * snh0_759[k]
                   - f_8 * pc_z[k] * snh1_759[k];

        t_970[k] = f_14 * sng_557[k]
                   + f_3 * pc_y[k] * sog_692[k];
    }

#pragma omp simd aligned(t_971, t_972, t_973, pb_z, pc_x, pc_z, snh0_762, sng_543, sng_695, \
                         snh1_762, sof0_465, sof1_465, sog_693, \
                         sog_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_971[k] = f_10 * sng_695[k]
                   + f_4 * sof0_465[k]
                   - f_5 * sof1_465[k]
                   + f_3 * pc_x[k] * sog_695[k];

        t_972[k] = pb_z[k] * snh0_762[k]
                   - f_8 * pc_z[k] * snh1_762[k];

        t_973[k] = f_9 * sng_543[k]
                   + f_3 * pc_z[k] * sog_693[k];
    }

#pragma omp simd aligned(t_974, t_975, t_976, t_977, pc_x, pc_y, sng_560, sng_699, sng_700, \
                         sng_701, sof0_469, sof1_469, sog_695, sog_699, sog_700, \
                         sog_701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_974[k] = f_14 * sng_560[k]
                   + f_3 * pc_y[k] * sog_695[k];

        t_975[k] = f_10 * sng_699[k]
                   + f_6 * sof0_469[k]
                   - f_7 * sof1_469[k]
                   + f_3 * pc_x[k] * sog_699[k];

        t_976[k] = f_10 * sng_700[k]
                   + f_3 * pc_x[k] * sog_700[k];

        t_977[k] = f_10 * sng_701[k]
                   + f_3 * pc_x[k] * sog_701[k];
    }

#pragma omp simd aligned(t_978, t_979, t_980, t_981, pb_z, pc_x, pc_z, snh0_771, sng_702, \
                         sng_703, sng_704, snh1_771, sog_702, sog_703, \
                         sog_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_978[k] = f_10 * sng_702[k]
                   + f_3 * pc_x[k] * sog_702[k];

        t_979[k] = f_10 * sng_703[k]
                   + f_3 * pc_x[k] * sog_703[k];

        t_980[k] = f_10 * sng_704[k]
                   + f_3 * pc_x[k] * sog_704[k];

        t_981[k] = pb_z[k] * snh0_771[k]
                   - f_8 * pc_z[k] * snh1_771[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, pc_y, pc_z, sng_550, sng_567, sng_568, sof0_468, \
                         sof0_469, sof1_468, sof1_469, sog_700, sog_702, \
                         sog_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = f_9 * sng_550[k]
                   + f_3 * pc_z[k] * sog_700[k];

        t_983[k] = f_14 * sng_567[k]
                   + f_4 * sof0_468[k]
                   - f_5 * sof1_468[k]
                   + f_3 * pc_y[k] * sog_702[k];

        t_984[k] = f_14 * sng_568[k]
                   + f_6 * sof0_469[k]
                   - f_7 * sof1_469[k]
                   + f_3 * pc_y[k] * sog_703[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, pc_x, pc_y, pc_z, sng_554, sng_569, sng_705, \
                         sof0_469, sof0_470, sof1_469, sof1_470, sog_704, \
                         sog_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_14 * sng_569[k]
                   + f_3 * pc_y[k] * sog_704[k];

        t_986[k] = f_9 * sng_554[k]
                   + f_1 * sof0_469[k]
                   - f_2 * sof1_469[k]
                   + f_3 * pc_z[k] * sog_704[k];

        t_987[k] = f_10 * sng_705[k]
                   + f_1 * sof0_470[k]
                   - f_2 * sof1_470[k]
                   + f_3 * pc_x[k] * sog_705[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, t_991, pc_x, pc_y, pc_z, sng_555, sng_570, \
                         sng_572, sng_708, sof0_473, sof1_473, sog_705, sog_707, \
                         sog_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_15 * sng_570[k]
                   + f_3 * pc_y[k] * sog_705[k];

        t_989[k] = f_10 * sng_555[k]
                   + f_3 * pc_z[k] * sog_705[k];

        t_990[k] = f_10 * sng_708[k]
                   + f_4 * sof0_473[k]
                   - f_5 * sof1_473[k]
                   + f_3 * pc_x[k] * sog_708[k];

        t_991[k] = f_15 * sng_572[k]
                   + f_3 * pc_y[k] * sog_707[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, pc_x, pc_z, sng_558, sng_710, sng_711, sof0_475, \
                         sof0_476, sof1_475, sof1_476, sog_708, sog_710, \
                         sog_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_10 * sng_710[k]
                   + f_4 * sof0_475[k]
                   - f_5 * sof1_475[k]
                   + f_3 * pc_x[k] * sog_710[k];

        t_993[k] = f_10 * sng_711[k]
                   + f_6 * sof0_476[k]
                   - f_7 * sof1_476[k]
                   + f_3 * pc_x[k] * sog_711[k];

        t_994[k] = f_10 * sng_558[k]
                   + f_3 * pc_z[k] * sog_708[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, pc_x, pc_y, sng_575, sng_714, sng_715, \
                         sng_716, sof0_479, sof1_479, sog_710, sog_714, sog_715, \
                         sog_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_15 * sng_575[k]
                   + f_3 * pc_y[k] * sog_710[k];

        t_996[k] = f_10 * sng_714[k]
                   + f_6 * sof0_479[k]
                   - f_7 * sof1_479[k]
                   + f_3 * pc_x[k] * sog_714[k];

        t_997[k] = f_10 * sng_715[k]
                   + f_3 * pc_x[k] * sog_715[k];

        t_998[k] = f_10 * sng_716[k]
                   + f_3 * pc_x[k] * sog_716[k];
    }

#pragma omp simd aligned(t_999, t_1000, t_1001, t_1002, pc_x, pc_y, sng_580, sng_717, sng_718, \
                         sng_719, sof0_476, sof1_476, sog_715, sog_717, sog_718, \
                         sog_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_999[k] = f_10 * sng_717[k]
                   + f_3 * pc_x[k] * sog_717[k];

        t_1000[k] = f_10 * sng_718[k]
                    + f_3 * pc_x[k] * sog_718[k];

        t_1001[k] = f_10 * sng_719[k]
                    + f_3 * pc_x[k] * sog_719[k];

        t_1002[k] = f_15 * sng_580[k]
                    + f_1 * sof0_476[k]
                    - f_2 * sof1_476[k]
                    + f_3 * pc_y[k] * sog_715[k];
    }

#pragma omp simd aligned(t_1003, t_1004, t_1005, pc_y, pc_z, sng_565, sng_582, sng_583, \
                         sof0_478, sof0_479, sof1_478, sof1_479, sog_715, sog_717, \
                         sog_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1003[k] = f_10 * sng_565[k]
                    + f_3 * pc_z[k] * sog_715[k];

        t_1004[k] = f_15 * sng_582[k]
                    + f_4 * sof0_478[k]
                    - f_5 * sof1_478[k]
                    + f_3 * pc_y[k] * sog_717[k];

        t_1005[k] = f_15 * sng_583[k]
                    + f_6 * sof0_479[k]
                    - f_7 * sof1_479[k]
                    + f_3 * pc_y[k] * sog_718[k];
    }

#pragma omp simd aligned(t_1006, t_1007, t_1008, pc_x, pc_y, pc_z, sng_569, sng_584, sng_720, \
                         sof0_479, sof0_480, sof1_479, sof1_480, sog_719, \
                         sog_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1006[k] = f_15 * sng_584[k]
                    + f_3 * pc_y[k] * sog_719[k];

        t_1007[k] = f_10 * sng_569[k]
                    + f_1 * sof0_479[k]
                    - f_2 * sof1_479[k]
                    + f_3 * pc_z[k] * sog_719[k];

        t_1008[k] = f_10 * sng_720[k]
                    + f_1 * sof0_480[k]
                    - f_2 * sof1_480[k]
                    + f_3 * pc_x[k] * sog_720[k];
    }

#pragma omp simd aligned(t_1009, t_1010, t_1011, t_1012, pc_x, pc_y, pc_z, sng_570, sng_585, \
                         sng_587, sng_723, sof0_483, sof1_483, sog_720, sog_722, \
                         sog_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1009[k] = f_17 * sng_585[k]
                    + f_3 * pc_y[k] * sog_720[k];

        t_1010[k] = f_11 * sng_570[k]
                    + f_3 * pc_z[k] * sog_720[k];

        t_1011[k] = f_10 * sng_723[k]
                    + f_4 * sof0_483[k]
                    - f_5 * sof1_483[k]
                    + f_3 * pc_x[k] * sog_723[k];

        t_1012[k] = f_17 * sng_587[k]
                    + f_3 * pc_y[k] * sog_722[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, pc_x, pc_z, sng_573, sng_725, sng_726, \
                         sof0_485, sof0_486, sof1_485, sof1_486, sog_723, sog_725, \
                         sog_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = f_10 * sng_725[k]
                    + f_4 * sof0_485[k]
                    - f_5 * sof1_485[k]
                    + f_3 * pc_x[k] * sog_725[k];

        t_1014[k] = f_10 * sng_726[k]
                    + f_6 * sof0_486[k]
                    - f_7 * sof1_486[k]
                    + f_3 * pc_x[k] * sog_726[k];

        t_1015[k] = f_11 * sng_573[k]
                    + f_3 * pc_z[k] * sog_723[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pc_x, pc_y, sng_590, sng_729, \
                         sng_730, sng_731, sof0_489, sof1_489, sog_725, sog_729, sog_730, \
                         sog_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_17 * sng_590[k]
                    + f_3 * pc_y[k] * sog_725[k];

        t_1017[k] = f_10 * sng_729[k]
                    + f_6 * sof0_489[k]
                    - f_7 * sof1_489[k]
                    + f_3 * pc_x[k] * sog_729[k];

        t_1018[k] = f_10 * sng_730[k]
                    + f_3 * pc_x[k] * sog_730[k];

        t_1019[k] = f_10 * sng_731[k]
                    + f_3 * pc_x[k] * sog_731[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, pc_x, pc_y, sng_595, sng_732, \
                         sng_733, sng_734, sof0_486, sof1_486, sog_730, sog_732, sog_733, \
                         sog_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_10 * sng_732[k]
                    + f_3 * pc_x[k] * sog_732[k];

        t_1021[k] = f_10 * sng_733[k]
                    + f_3 * pc_x[k] * sog_733[k];

        t_1022[k] = f_10 * sng_734[k]
                    + f_3 * pc_x[k] * sog_734[k];

        t_1023[k] = f_17 * sng_595[k]
                    + f_1 * sof0_486[k]
                    - f_2 * sof1_486[k]
                    + f_3 * pc_y[k] * sog_730[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, pc_y, pc_z, sng_580, sng_597, sng_598, \
                         sof0_488, sof0_489, sof1_488, sof1_489, sog_730, sog_732, \
                         sog_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_11 * sng_580[k]
                    + f_3 * pc_z[k] * sog_730[k];

        t_1025[k] = f_17 * sng_597[k]
                    + f_4 * sof0_488[k]
                    - f_5 * sof1_488[k]
                    + f_3 * pc_y[k] * sog_732[k];

        t_1026[k] = f_17 * sng_598[k]
                    + f_6 * sof0_489[k]
                    - f_7 * sof1_489[k]
                    + f_3 * pc_y[k] * sog_733[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, pc_x, pc_y, pc_z, sng_584, sng_599, sng_735, \
                         sof0_489, sof0_490, sof1_489, sof1_490, sog_734, \
                         sog_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_17 * sng_599[k]
                    + f_3 * pc_y[k] * sog_734[k];

        t_1028[k] = f_11 * sng_584[k]
                    + f_1 * sof0_489[k]
                    - f_2 * sof1_489[k]
                    + f_3 * pc_z[k] * sog_734[k];

        t_1029[k] = f_10 * sng_735[k]
                    + f_1 * sof0_490[k]
                    - f_2 * sof1_490[k]
                    + f_3 * pc_x[k] * sog_735[k];
    }

#pragma omp simd aligned(t_1030, t_1031, t_1032, t_1033, pc_x, pc_y, pc_z, sng_585, sng_600, \
                         sng_602, sng_738, sof0_493, sof1_493, sog_735, sog_737, \
                         sog_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1030[k] = f_18 * sng_600[k]
                    + f_3 * pc_y[k] * sog_735[k];

        t_1031[k] = f_16 * sng_585[k]
                    + f_3 * pc_z[k] * sog_735[k];

        t_1032[k] = f_10 * sng_738[k]
                    + f_4 * sof0_493[k]
                    - f_5 * sof1_493[k]
                    + f_3 * pc_x[k] * sog_738[k];

        t_1033[k] = f_18 * sng_602[k]
                    + f_3 * pc_y[k] * sog_737[k];
    }
}

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snh0,
                                                          const size_t sng, const size_t snh1,
                                                          const size_t sof0, const size_t sof1,
                                                          const size_t sog, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 4.5 / q;
    const auto f_14 = 4.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_1034 = buffer.data(target + 1034);
    auto *t_1035 = buffer.data(target + 1035);
    auto *t_1036 = buffer.data(target + 1036);
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);
    auto *t_1041 = buffer.data(target + 1041);
    auto *t_1042 = buffer.data(target + 1042);
    auto *t_1043 = buffer.data(target + 1043);
    auto *t_1044 = buffer.data(target + 1044);
    auto *t_1045 = buffer.data(target + 1045);
    auto *t_1046 = buffer.data(target + 1046);
    auto *t_1047 = buffer.data(target + 1047);
    auto *t_1048 = buffer.data(target + 1048);
    auto *t_1049 = buffer.data(target + 1049);
    auto *t_1050 = buffer.data(target + 1050);
    auto *t_1051 = buffer.data(target + 1051);
    auto *t_1052 = buffer.data(target + 1052);
    auto *t_1053 = buffer.data(target + 1053);
    auto *t_1054 = buffer.data(target + 1054);
    auto *t_1055 = buffer.data(target + 1055);
    auto *t_1056 = buffer.data(target + 1056);
    auto *t_1057 = buffer.data(target + 1057);
    auto *t_1058 = buffer.data(target + 1058);
    auto *t_1059 = buffer.data(target + 1059);
    auto *t_1060 = buffer.data(target + 1060);
    auto *t_1061 = buffer.data(target + 1061);
    auto *t_1062 = buffer.data(target + 1062);
    auto *t_1063 = buffer.data(target + 1063);
    auto *t_1064 = buffer.data(target + 1064);
    auto *t_1065 = buffer.data(target + 1065);
    auto *t_1066 = buffer.data(target + 1066);
    auto *t_1067 = buffer.data(target + 1067);
    auto *t_1068 = buffer.data(target + 1068);
    auto *t_1069 = buffer.data(target + 1069);
    auto *t_1070 = buffer.data(target + 1070);
    auto *t_1071 = buffer.data(target + 1071);
    auto *t_1072 = buffer.data(target + 1072);
    auto *t_1073 = buffer.data(target + 1073);
    auto *t_1074 = buffer.data(target + 1074);
    auto *t_1075 = buffer.data(target + 1075);
    auto *t_1076 = buffer.data(target + 1076);
    auto *t_1077 = buffer.data(target + 1077);
    auto *t_1078 = buffer.data(target + 1078);
    auto *t_1079 = buffer.data(target + 1079);
    auto *t_1080 = buffer.data(target + 1080);
    auto *t_1081 = buffer.data(target + 1081);
    auto *t_1082 = buffer.data(target + 1082);
    auto *t_1083 = buffer.data(target + 1083);
    auto *t_1084 = buffer.data(target + 1084);
    auto *t_1085 = buffer.data(target + 1085);
    auto *t_1086 = buffer.data(target + 1086);
    auto *t_1087 = buffer.data(target + 1087);
    auto *t_1088 = buffer.data(target + 1088);
    auto *t_1089 = buffer.data(target + 1089);
    auto *t_1090 = buffer.data(target + 1090);
    auto *t_1091 = buffer.data(target + 1091);
    auto *t_1092 = buffer.data(target + 1092);
    auto *t_1093 = buffer.data(target + 1093);
    auto *t_1094 = buffer.data(target + 1094);
    auto *t_1095 = buffer.data(target + 1095);
    auto *t_1096 = buffer.data(target + 1096);
    auto *t_1097 = buffer.data(target + 1097);
    auto *t_1098 = buffer.data(target + 1098);
    auto *t_1099 = buffer.data(target + 1099);
    auto *t_1100 = buffer.data(target + 1100);
    auto *t_1101 = buffer.data(target + 1101);
    auto *t_1102 = buffer.data(target + 1102);
    auto *t_1103 = buffer.data(target + 1103);
    auto *t_1104 = buffer.data(target + 1104);
    auto *t_1105 = buffer.data(target + 1105);
    auto *t_1106 = buffer.data(target + 1106);
    auto *t_1107 = buffer.data(target + 1107);
    auto *t_1108 = buffer.data(target + 1108);
    auto *t_1109 = buffer.data(target + 1109);
    auto *t_1110 = buffer.data(target + 1110);
    auto *t_1111 = buffer.data(target + 1111);
    auto *t_1112 = buffer.data(target + 1112);
    auto *t_1113 = buffer.data(target + 1113);
    auto *t_1114 = buffer.data(target + 1114);
    auto *t_1115 = buffer.data(target + 1115);
    auto *t_1116 = buffer.data(target + 1116);
    auto *t_1117 = buffer.data(target + 1117);
    auto *t_1118 = buffer.data(target + 1118);
    auto *t_1119 = buffer.data(target + 1119);
    auto *t_1120 = buffer.data(target + 1120);
    auto *t_1121 = buffer.data(target + 1121);
    auto *t_1122 = buffer.data(target + 1122);
    auto *t_1123 = buffer.data(target + 1123);
    auto *t_1124 = buffer.data(target + 1124);
    auto *t_1125 = buffer.data(target + 1125);
    auto *t_1126 = buffer.data(target + 1126);
    auto *t_1127 = buffer.data(target + 1127);
    auto *t_1128 = buffer.data(target + 1128);
    auto *t_1129 = buffer.data(target + 1129);
    auto *t_1130 = buffer.data(target + 1130);
    auto *t_1131 = buffer.data(target + 1131);
    auto *t_1132 = buffer.data(target + 1132);
    auto *t_1133 = buffer.data(target + 1133);
    auto *t_1134 = buffer.data(target + 1134);
    auto *t_1135 = buffer.data(target + 1135);
    auto *t_1136 = buffer.data(target + 1136);
    auto *t_1137 = buffer.data(target + 1137);
    auto *t_1138 = buffer.data(target + 1138);
    auto *t_1139 = buffer.data(target + 1139);
    auto *t_1140 = buffer.data(target + 1140);
    auto *t_1141 = buffer.data(target + 1141);
    auto *t_1142 = buffer.data(target + 1142);
    auto *t_1143 = buffer.data(target + 1143);
    auto *t_1144 = buffer.data(target + 1144);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_924 = buffer.data(snh0 + 924);
    const auto *snh0_927 = buffer.data(snh0 + 927);
    const auto *snh0_929 = buffer.data(snh0 + 929);
    const auto *snh0_930 = buffer.data(snh0 + 930);
    const auto *snh0_933 = buffer.data(snh0 + 933);
    const auto *snh0_944 = buffer.data(snh0 + 944);

    const auto *sng_588 = buffer.data(sng + 588);
    const auto *sng_595 = buffer.data(sng + 595);
    const auto *sng_599 = buffer.data(sng + 599);
    const auto *sng_600 = buffer.data(sng + 600);
    const auto *sng_603 = buffer.data(sng + 603);
    const auto *sng_605 = buffer.data(sng + 605);
    const auto *sng_610 = buffer.data(sng + 610);
    const auto *sng_612 = buffer.data(sng + 612);
    const auto *sng_613 = buffer.data(sng + 613);
    const auto *sng_614 = buffer.data(sng + 614);
    const auto *sng_615 = buffer.data(sng + 615);
    const auto *sng_617 = buffer.data(sng + 617);
    const auto *sng_618 = buffer.data(sng + 618);
    const auto *sng_620 = buffer.data(sng + 620);
    const auto *sng_625 = buffer.data(sng + 625);
    const auto *sng_627 = buffer.data(sng + 627);
    const auto *sng_628 = buffer.data(sng + 628);
    const auto *sng_629 = buffer.data(sng + 629);
    const auto *sng_630 = buffer.data(sng + 630);
    const auto *sng_632 = buffer.data(sng + 632);
    const auto *sng_633 = buffer.data(sng + 633);
    const auto *sng_635 = buffer.data(sng + 635);
    const auto *sng_640 = buffer.data(sng + 640);
    const auto *sng_642 = buffer.data(sng + 642);
    const auto *sng_643 = buffer.data(sng + 643);
    const auto *sng_644 = buffer.data(sng + 644);
    const auto *sng_645 = buffer.data(sng + 645);
    const auto *sng_647 = buffer.data(sng + 647);
    const auto *sng_648 = buffer.data(sng + 648);
    const auto *sng_650 = buffer.data(sng + 650);
    const auto *sng_655 = buffer.data(sng + 655);
    const auto *sng_657 = buffer.data(sng + 657);
    const auto *sng_658 = buffer.data(sng + 658);
    const auto *sng_659 = buffer.data(sng + 659);
    const auto *sng_660 = buffer.data(sng + 660);
    const auto *sng_661 = buffer.data(sng + 661);
    const auto *sng_662 = buffer.data(sng + 662);
    const auto *sng_663 = buffer.data(sng + 663);
    const auto *sng_665 = buffer.data(sng + 665);
    const auto *sng_670 = buffer.data(sng + 670);
    const auto *sng_672 = buffer.data(sng + 672);
    const auto *sng_673 = buffer.data(sng + 673);
    const auto *sng_674 = buffer.data(sng + 674);
    const auto *sng_740 = buffer.data(sng + 740);
    const auto *sng_741 = buffer.data(sng + 741);
    const auto *sng_744 = buffer.data(sng + 744);
    const auto *sng_745 = buffer.data(sng + 745);
    const auto *sng_746 = buffer.data(sng + 746);
    const auto *sng_747 = buffer.data(sng + 747);
    const auto *sng_748 = buffer.data(sng + 748);
    const auto *sng_749 = buffer.data(sng + 749);
    const auto *sng_750 = buffer.data(sng + 750);
    const auto *sng_753 = buffer.data(sng + 753);
    const auto *sng_755 = buffer.data(sng + 755);
    const auto *sng_756 = buffer.data(sng + 756);
    const auto *sng_759 = buffer.data(sng + 759);
    const auto *sng_760 = buffer.data(sng + 760);
    const auto *sng_761 = buffer.data(sng + 761);
    const auto *sng_762 = buffer.data(sng + 762);
    const auto *sng_763 = buffer.data(sng + 763);
    const auto *sng_764 = buffer.data(sng + 764);
    const auto *sng_765 = buffer.data(sng + 765);
    const auto *sng_768 = buffer.data(sng + 768);
    const auto *sng_770 = buffer.data(sng + 770);
    const auto *sng_771 = buffer.data(sng + 771);
    const auto *sng_774 = buffer.data(sng + 774);
    const auto *sng_775 = buffer.data(sng + 775);
    const auto *sng_776 = buffer.data(sng + 776);
    const auto *sng_777 = buffer.data(sng + 777);
    const auto *sng_778 = buffer.data(sng + 778);
    const auto *sng_779 = buffer.data(sng + 779);
    const auto *sng_780 = buffer.data(sng + 780);
    const auto *sng_783 = buffer.data(sng + 783);
    const auto *sng_785 = buffer.data(sng + 785);
    const auto *sng_786 = buffer.data(sng + 786);
    const auto *sng_789 = buffer.data(sng + 789);
    const auto *sng_790 = buffer.data(sng + 790);
    const auto *sng_791 = buffer.data(sng + 791);
    const auto *sng_792 = buffer.data(sng + 792);
    const auto *sng_793 = buffer.data(sng + 793);
    const auto *sng_794 = buffer.data(sng + 794);
    const auto *sng_805 = buffer.data(sng + 805);
    const auto *sng_806 = buffer.data(sng + 806);
    const auto *sng_807 = buffer.data(sng + 807);
    const auto *sng_808 = buffer.data(sng + 808);
    const auto *sng_809 = buffer.data(sng + 809);
    const auto *sng_810 = buffer.data(sng + 810);
    const auto *sng_813 = buffer.data(sng + 813);
    const auto *sng_815 = buffer.data(sng + 815);
    const auto *sng_816 = buffer.data(sng + 816);
    const auto *sng_819 = buffer.data(sng + 819);
    const auto *sng_820 = buffer.data(sng + 820);

    const auto *snh1_924 = buffer.data(snh1 + 924);
    const auto *snh1_927 = buffer.data(snh1 + 927);
    const auto *snh1_929 = buffer.data(snh1 + 929);
    const auto *snh1_930 = buffer.data(snh1 + 930);
    const auto *snh1_933 = buffer.data(snh1 + 933);
    const auto *snh1_944 = buffer.data(snh1 + 944);

    const auto *sof0_495 = buffer.data(sof0 + 495);
    const auto *sof0_496 = buffer.data(sof0 + 496);
    const auto *sof0_498 = buffer.data(sof0 + 498);
    const auto *sof0_499 = buffer.data(sof0 + 499);
    const auto *sof0_500 = buffer.data(sof0 + 500);
    const auto *sof0_503 = buffer.data(sof0 + 503);
    const auto *sof0_505 = buffer.data(sof0 + 505);
    const auto *sof0_506 = buffer.data(sof0 + 506);
    const auto *sof0_508 = buffer.data(sof0 + 508);
    const auto *sof0_509 = buffer.data(sof0 + 509);
    const auto *sof0_510 = buffer.data(sof0 + 510);
    const auto *sof0_513 = buffer.data(sof0 + 513);
    const auto *sof0_515 = buffer.data(sof0 + 515);
    const auto *sof0_516 = buffer.data(sof0 + 516);
    const auto *sof0_518 = buffer.data(sof0 + 518);
    const auto *sof0_519 = buffer.data(sof0 + 519);
    const auto *sof0_520 = buffer.data(sof0 + 520);
    const auto *sof0_523 = buffer.data(sof0 + 523);
    const auto *sof0_525 = buffer.data(sof0 + 525);
    const auto *sof0_526 = buffer.data(sof0 + 526);
    const auto *sof0_528 = buffer.data(sof0 + 528);
    const auto *sof0_529 = buffer.data(sof0 + 529);
    const auto *sof0_536 = buffer.data(sof0 + 536);
    const auto *sof0_538 = buffer.data(sof0 + 538);
    const auto *sof0_539 = buffer.data(sof0 + 539);
    const auto *sof0_540 = buffer.data(sof0 + 540);
    const auto *sof0_543 = buffer.data(sof0 + 543);
    const auto *sof0_545 = buffer.data(sof0 + 545);
    const auto *sof0_546 = buffer.data(sof0 + 546);
    const auto *sof0_549 = buffer.data(sof0 + 549);

    const auto *sof1_495 = buffer.data(sof1 + 495);
    const auto *sof1_496 = buffer.data(sof1 + 496);
    const auto *sof1_498 = buffer.data(sof1 + 498);
    const auto *sof1_499 = buffer.data(sof1 + 499);
    const auto *sof1_500 = buffer.data(sof1 + 500);
    const auto *sof1_503 = buffer.data(sof1 + 503);
    const auto *sof1_505 = buffer.data(sof1 + 505);
    const auto *sof1_506 = buffer.data(sof1 + 506);
    const auto *sof1_508 = buffer.data(sof1 + 508);
    const auto *sof1_509 = buffer.data(sof1 + 509);
    const auto *sof1_510 = buffer.data(sof1 + 510);
    const auto *sof1_513 = buffer.data(sof1 + 513);
    const auto *sof1_515 = buffer.data(sof1 + 515);
    const auto *sof1_516 = buffer.data(sof1 + 516);
    const auto *sof1_518 = buffer.data(sof1 + 518);
    const auto *sof1_519 = buffer.data(sof1 + 519);
    const auto *sof1_520 = buffer.data(sof1 + 520);
    const auto *sof1_523 = buffer.data(sof1 + 523);
    const auto *sof1_525 = buffer.data(sof1 + 525);
    const auto *sof1_526 = buffer.data(sof1 + 526);
    const auto *sof1_528 = buffer.data(sof1 + 528);
    const auto *sof1_529 = buffer.data(sof1 + 529);
    const auto *sof1_536 = buffer.data(sof1 + 536);
    const auto *sof1_538 = buffer.data(sof1 + 538);
    const auto *sof1_539 = buffer.data(sof1 + 539);
    const auto *sof1_540 = buffer.data(sof1 + 540);
    const auto *sof1_543 = buffer.data(sof1 + 543);
    const auto *sof1_545 = buffer.data(sof1 + 545);
    const auto *sof1_546 = buffer.data(sof1 + 546);
    const auto *sof1_549 = buffer.data(sof1 + 549);

    const auto *sog_738 = buffer.data(sog + 738);
    const auto *sog_740 = buffer.data(sog + 740);
    const auto *sog_741 = buffer.data(sog + 741);
    const auto *sog_744 = buffer.data(sog + 744);
    const auto *sog_745 = buffer.data(sog + 745);
    const auto *sog_746 = buffer.data(sog + 746);
    const auto *sog_747 = buffer.data(sog + 747);
    const auto *sog_748 = buffer.data(sog + 748);
    const auto *sog_749 = buffer.data(sog + 749);
    const auto *sog_750 = buffer.data(sog + 750);
    const auto *sog_752 = buffer.data(sog + 752);
    const auto *sog_753 = buffer.data(sog + 753);
    const auto *sog_755 = buffer.data(sog + 755);
    const auto *sog_756 = buffer.data(sog + 756);
    const auto *sog_759 = buffer.data(sog + 759);
    const auto *sog_760 = buffer.data(sog + 760);
    const auto *sog_761 = buffer.data(sog + 761);
    const auto *sog_762 = buffer.data(sog + 762);
    const auto *sog_763 = buffer.data(sog + 763);
    const auto *sog_764 = buffer.data(sog + 764);
    const auto *sog_765 = buffer.data(sog + 765);
    const auto *sog_767 = buffer.data(sog + 767);
    const auto *sog_768 = buffer.data(sog + 768);
    const auto *sog_770 = buffer.data(sog + 770);
    const auto *sog_771 = buffer.data(sog + 771);
    const auto *sog_774 = buffer.data(sog + 774);
    const auto *sog_775 = buffer.data(sog + 775);
    const auto *sog_776 = buffer.data(sog + 776);
    const auto *sog_777 = buffer.data(sog + 777);
    const auto *sog_778 = buffer.data(sog + 778);
    const auto *sog_779 = buffer.data(sog + 779);
    const auto *sog_780 = buffer.data(sog + 780);
    const auto *sog_782 = buffer.data(sog + 782);
    const auto *sog_783 = buffer.data(sog + 783);
    const auto *sog_785 = buffer.data(sog + 785);
    const auto *sog_786 = buffer.data(sog + 786);
    const auto *sog_789 = buffer.data(sog + 789);
    const auto *sog_790 = buffer.data(sog + 790);
    const auto *sog_791 = buffer.data(sog + 791);
    const auto *sog_792 = buffer.data(sog + 792);
    const auto *sog_793 = buffer.data(sog + 793);
    const auto *sog_794 = buffer.data(sog + 794);
    const auto *sog_795 = buffer.data(sog + 795);
    const auto *sog_797 = buffer.data(sog + 797);
    const auto *sog_798 = buffer.data(sog + 798);
    const auto *sog_800 = buffer.data(sog + 800);
    const auto *sog_805 = buffer.data(sog + 805);
    const auto *sog_806 = buffer.data(sog + 806);
    const auto *sog_807 = buffer.data(sog + 807);
    const auto *sog_808 = buffer.data(sog + 808);
    const auto *sog_809 = buffer.data(sog + 809);
    const auto *sog_810 = buffer.data(sog + 810);
    const auto *sog_812 = buffer.data(sog + 812);
    const auto *sog_813 = buffer.data(sog + 813);
    const auto *sog_815 = buffer.data(sog + 815);
    const auto *sog_816 = buffer.data(sog + 816);
    const auto *sog_819 = buffer.data(sog + 819);
    const auto *sog_820 = buffer.data(sog + 820);

#pragma omp simd aligned(t_1034, t_1035, t_1036, pc_x, pc_z, sng_588, sng_740, sng_741, \
                         sof0_495, sof0_496, sof1_495, sof1_496, sog_738, sog_740, \
                         sog_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = f_10 * sng_740[k]
                    + f_4 * sof0_495[k]
                    - f_5 * sof1_495[k]
                    + f_3 * pc_x[k] * sog_740[k];

        t_1035[k] = f_10 * sng_741[k]
                    + f_6 * sof0_496[k]
                    - f_7 * sof1_496[k]
                    + f_3 * pc_x[k] * sog_741[k];

        t_1036[k] = f_16 * sng_588[k]
                    + f_3 * pc_z[k] * sog_738[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, t_1040, pc_x, pc_y, sng_605, sng_744, \
                         sng_745, sng_746, sof0_499, sof1_499, sog_740, sog_744, sog_745, \
                         sog_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_18 * sng_605[k]
                    + f_3 * pc_y[k] * sog_740[k];

        t_1038[k] = f_10 * sng_744[k]
                    + f_6 * sof0_499[k]
                    - f_7 * sof1_499[k]
                    + f_3 * pc_x[k] * sog_744[k];

        t_1039[k] = f_10 * sng_745[k]
                    + f_3 * pc_x[k] * sog_745[k];

        t_1040[k] = f_10 * sng_746[k]
                    + f_3 * pc_x[k] * sog_746[k];
    }

#pragma omp simd aligned(t_1041, t_1042, t_1043, t_1044, pc_x, pc_y, sng_610, sng_747, \
                         sng_748, sng_749, sof0_496, sof1_496, sog_745, sog_747, sog_748, \
                         sog_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1041[k] = f_10 * sng_747[k]
                    + f_3 * pc_x[k] * sog_747[k];

        t_1042[k] = f_10 * sng_748[k]
                    + f_3 * pc_x[k] * sog_748[k];

        t_1043[k] = f_10 * sng_749[k]
                    + f_3 * pc_x[k] * sog_749[k];

        t_1044[k] = f_18 * sng_610[k]
                    + f_1 * sof0_496[k]
                    - f_2 * sof1_496[k]
                    + f_3 * pc_y[k] * sog_745[k];
    }

#pragma omp simd aligned(t_1045, t_1046, t_1047, pc_y, pc_z, sng_595, sng_612, sng_613, \
                         sof0_498, sof0_499, sof1_498, sof1_499, sog_745, sog_747, \
                         sog_748 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1045[k] = f_16 * sng_595[k]
                    + f_3 * pc_z[k] * sog_745[k];

        t_1046[k] = f_18 * sng_612[k]
                    + f_4 * sof0_498[k]
                    - f_5 * sof1_498[k]
                    + f_3 * pc_y[k] * sog_747[k];

        t_1047[k] = f_18 * sng_613[k]
                    + f_6 * sof0_499[k]
                    - f_7 * sof1_499[k]
                    + f_3 * pc_y[k] * sog_748[k];
    }

#pragma omp simd aligned(t_1048, t_1049, t_1050, pc_x, pc_y, pc_z, sng_599, sng_614, sng_750, \
                         sof0_499, sof0_500, sof1_499, sof1_500, sog_749, \
                         sog_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1048[k] = f_18 * sng_614[k]
                    + f_3 * pc_y[k] * sog_749[k];

        t_1049[k] = f_16 * sng_599[k]
                    + f_1 * sof0_499[k]
                    - f_2 * sof1_499[k]
                    + f_3 * pc_z[k] * sog_749[k];

        t_1050[k] = f_10 * sng_750[k]
                    + f_1 * sof0_500[k]
                    - f_2 * sof1_500[k]
                    + f_3 * pc_x[k] * sog_750[k];
    }

#pragma omp simd aligned(t_1051, t_1052, t_1053, t_1054, pc_x, pc_y, pc_z, sng_600, sng_615, \
                         sng_617, sng_753, sof0_503, sof1_503, sog_750, sog_752, \
                         sog_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1051[k] = f_16 * sng_615[k]
                    + f_3 * pc_y[k] * sog_750[k];

        t_1052[k] = f_18 * sng_600[k]
                    + f_3 * pc_z[k] * sog_750[k];

        t_1053[k] = f_10 * sng_753[k]
                    + f_4 * sof0_503[k]
                    - f_5 * sof1_503[k]
                    + f_3 * pc_x[k] * sog_753[k];

        t_1054[k] = f_16 * sng_617[k]
                    + f_3 * pc_y[k] * sog_752[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pc_x, pc_z, sng_603, sng_755, sng_756, \
                         sof0_505, sof0_506, sof1_505, sof1_506, sog_753, sog_755, \
                         sog_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_10 * sng_755[k]
                    + f_4 * sof0_505[k]
                    - f_5 * sof1_505[k]
                    + f_3 * pc_x[k] * sog_755[k];

        t_1056[k] = f_10 * sng_756[k]
                    + f_6 * sof0_506[k]
                    - f_7 * sof1_506[k]
                    + f_3 * pc_x[k] * sog_756[k];

        t_1057[k] = f_18 * sng_603[k]
                    + f_3 * pc_z[k] * sog_753[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, t_1061, pc_x, pc_y, sng_620, sng_759, \
                         sng_760, sng_761, sof0_509, sof1_509, sog_755, sog_759, sog_760, \
                         sog_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_16 * sng_620[k]
                    + f_3 * pc_y[k] * sog_755[k];

        t_1059[k] = f_10 * sng_759[k]
                    + f_6 * sof0_509[k]
                    - f_7 * sof1_509[k]
                    + f_3 * pc_x[k] * sog_759[k];

        t_1060[k] = f_10 * sng_760[k]
                    + f_3 * pc_x[k] * sog_760[k];

        t_1061[k] = f_10 * sng_761[k]
                    + f_3 * pc_x[k] * sog_761[k];
    }

#pragma omp simd aligned(t_1062, t_1063, t_1064, t_1065, pc_x, pc_y, sng_625, sng_762, \
                         sng_763, sng_764, sof0_506, sof1_506, sog_760, sog_762, sog_763, \
                         sog_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1062[k] = f_10 * sng_762[k]
                    + f_3 * pc_x[k] * sog_762[k];

        t_1063[k] = f_10 * sng_763[k]
                    + f_3 * pc_x[k] * sog_763[k];

        t_1064[k] = f_10 * sng_764[k]
                    + f_3 * pc_x[k] * sog_764[k];

        t_1065[k] = f_16 * sng_625[k]
                    + f_1 * sof0_506[k]
                    - f_2 * sof1_506[k]
                    + f_3 * pc_y[k] * sog_760[k];
    }

#pragma omp simd aligned(t_1066, t_1067, t_1068, pc_y, pc_z, sng_610, sng_627, sng_628, \
                         sof0_508, sof0_509, sof1_508, sof1_509, sog_760, sog_762, \
                         sog_763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1066[k] = f_18 * sng_610[k]
                    + f_3 * pc_z[k] * sog_760[k];

        t_1067[k] = f_16 * sng_627[k]
                    + f_4 * sof0_508[k]
                    - f_5 * sof1_508[k]
                    + f_3 * pc_y[k] * sog_762[k];

        t_1068[k] = f_16 * sng_628[k]
                    + f_6 * sof0_509[k]
                    - f_7 * sof1_509[k]
                    + f_3 * pc_y[k] * sog_763[k];
    }

#pragma omp simd aligned(t_1069, t_1070, t_1071, pc_x, pc_y, pc_z, sng_614, sng_629, sng_765, \
                         sof0_509, sof0_510, sof1_509, sof1_510, sog_764, \
                         sog_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1069[k] = f_16 * sng_629[k]
                    + f_3 * pc_y[k] * sog_764[k];

        t_1070[k] = f_18 * sng_614[k]
                    + f_1 * sof0_509[k]
                    - f_2 * sof1_509[k]
                    + f_3 * pc_z[k] * sog_764[k];

        t_1071[k] = f_10 * sng_765[k]
                    + f_1 * sof0_510[k]
                    - f_2 * sof1_510[k]
                    + f_3 * pc_x[k] * sog_765[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, t_1075, pc_x, pc_y, pc_z, sng_615, sng_630, \
                         sng_632, sng_768, sof0_513, sof1_513, sog_765, sog_767, \
                         sog_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = f_11 * sng_630[k]
                    + f_3 * pc_y[k] * sog_765[k];

        t_1073[k] = f_17 * sng_615[k]
                    + f_3 * pc_z[k] * sog_765[k];

        t_1074[k] = f_10 * sng_768[k]
                    + f_4 * sof0_513[k]
                    - f_5 * sof1_513[k]
                    + f_3 * pc_x[k] * sog_768[k];

        t_1075[k] = f_11 * sng_632[k]
                    + f_3 * pc_y[k] * sog_767[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_x, pc_z, sng_618, sng_770, sng_771, \
                         sof0_515, sof0_516, sof1_515, sof1_516, sog_768, sog_770, \
                         sog_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_10 * sng_770[k]
                    + f_4 * sof0_515[k]
                    - f_5 * sof1_515[k]
                    + f_3 * pc_x[k] * sog_770[k];

        t_1077[k] = f_10 * sng_771[k]
                    + f_6 * sof0_516[k]
                    - f_7 * sof1_516[k]
                    + f_3 * pc_x[k] * sog_771[k];

        t_1078[k] = f_17 * sng_618[k]
                    + f_3 * pc_z[k] * sog_768[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, pc_x, pc_y, sng_635, sng_774, \
                         sng_775, sng_776, sof0_519, sof1_519, sog_770, sog_774, sog_775, \
                         sog_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_11 * sng_635[k]
                    + f_3 * pc_y[k] * sog_770[k];

        t_1080[k] = f_10 * sng_774[k]
                    + f_6 * sof0_519[k]
                    - f_7 * sof1_519[k]
                    + f_3 * pc_x[k] * sog_774[k];

        t_1081[k] = f_10 * sng_775[k]
                    + f_3 * pc_x[k] * sog_775[k];

        t_1082[k] = f_10 * sng_776[k]
                    + f_3 * pc_x[k] * sog_776[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, t_1086, pc_x, pc_y, sng_640, sng_777, \
                         sng_778, sng_779, sof0_516, sof1_516, sog_775, sog_777, sog_778, \
                         sog_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = f_10 * sng_777[k]
                    + f_3 * pc_x[k] * sog_777[k];

        t_1084[k] = f_10 * sng_778[k]
                    + f_3 * pc_x[k] * sog_778[k];

        t_1085[k] = f_10 * sng_779[k]
                    + f_3 * pc_x[k] * sog_779[k];

        t_1086[k] = f_11 * sng_640[k]
                    + f_1 * sof0_516[k]
                    - f_2 * sof1_516[k]
                    + f_3 * pc_y[k] * sog_775[k];
    }

#pragma omp simd aligned(t_1087, t_1088, t_1089, pc_y, pc_z, sng_625, sng_642, sng_643, \
                         sof0_518, sof0_519, sof1_518, sof1_519, sog_775, sog_777, \
                         sog_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1087[k] = f_17 * sng_625[k]
                    + f_3 * pc_z[k] * sog_775[k];

        t_1088[k] = f_11 * sng_642[k]
                    + f_4 * sof0_518[k]
                    - f_5 * sof1_518[k]
                    + f_3 * pc_y[k] * sog_777[k];

        t_1089[k] = f_11 * sng_643[k]
                    + f_6 * sof0_519[k]
                    - f_7 * sof1_519[k]
                    + f_3 * pc_y[k] * sog_778[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, pc_x, pc_y, pc_z, sng_629, sng_644, sng_780, \
                         sof0_519, sof0_520, sof1_519, sof1_520, sog_779, \
                         sog_780 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = f_11 * sng_644[k]
                    + f_3 * pc_y[k] * sog_779[k];

        t_1091[k] = f_17 * sng_629[k]
                    + f_1 * sof0_519[k]
                    - f_2 * sof1_519[k]
                    + f_3 * pc_z[k] * sog_779[k];

        t_1092[k] = f_10 * sng_780[k]
                    + f_1 * sof0_520[k]
                    - f_2 * sof1_520[k]
                    + f_3 * pc_x[k] * sog_780[k];
    }

#pragma omp simd aligned(t_1093, t_1094, t_1095, t_1096, pc_x, pc_y, pc_z, sng_630, sng_645, \
                         sng_647, sng_783, sof0_523, sof1_523, sog_780, sog_782, \
                         sog_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1093[k] = f_10 * sng_645[k]
                    + f_3 * pc_y[k] * sog_780[k];

        t_1094[k] = f_15 * sng_630[k]
                    + f_3 * pc_z[k] * sog_780[k];

        t_1095[k] = f_10 * sng_783[k]
                    + f_4 * sof0_523[k]
                    - f_5 * sof1_523[k]
                    + f_3 * pc_x[k] * sog_783[k];

        t_1096[k] = f_10 * sng_647[k]
                    + f_3 * pc_y[k] * sog_782[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_z, sng_633, sng_785, sng_786, \
                         sof0_525, sof0_526, sof1_525, sof1_526, sog_783, sog_785, \
                         sog_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_10 * sng_785[k]
                    + f_4 * sof0_525[k]
                    - f_5 * sof1_525[k]
                    + f_3 * pc_x[k] * sog_785[k];

        t_1098[k] = f_10 * sng_786[k]
                    + f_6 * sof0_526[k]
                    - f_7 * sof1_526[k]
                    + f_3 * pc_x[k] * sog_786[k];

        t_1099[k] = f_15 * sng_633[k]
                    + f_3 * pc_z[k] * sog_783[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, t_1103, pc_x, pc_y, sng_650, sng_789, \
                         sng_790, sng_791, sof0_529, sof1_529, sog_785, sog_789, sog_790, \
                         sog_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_10 * sng_650[k]
                    + f_3 * pc_y[k] * sog_785[k];

        t_1101[k] = f_10 * sng_789[k]
                    + f_6 * sof0_529[k]
                    - f_7 * sof1_529[k]
                    + f_3 * pc_x[k] * sog_789[k];

        t_1102[k] = f_10 * sng_790[k]
                    + f_3 * pc_x[k] * sog_790[k];

        t_1103[k] = f_10 * sng_791[k]
                    + f_3 * pc_x[k] * sog_791[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, t_1107, pc_x, pc_y, sng_655, sng_792, \
                         sng_793, sng_794, sof0_526, sof1_526, sog_790, sog_792, sog_793, \
                         sog_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_10 * sng_792[k]
                    + f_3 * pc_x[k] * sog_792[k];

        t_1105[k] = f_10 * sng_793[k]
                    + f_3 * pc_x[k] * sog_793[k];

        t_1106[k] = f_10 * sng_794[k]
                    + f_3 * pc_x[k] * sog_794[k];

        t_1107[k] = f_10 * sng_655[k]
                    + f_1 * sof0_526[k]
                    - f_2 * sof1_526[k]
                    + f_3 * pc_y[k] * sog_790[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, pc_y, pc_z, sng_640, sng_657, sng_658, \
                         sof0_528, sof0_529, sof1_528, sof1_529, sog_790, sog_792, \
                         sog_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_15 * sng_640[k]
                    + f_3 * pc_z[k] * sog_790[k];

        t_1109[k] = f_10 * sng_657[k]
                    + f_4 * sof0_528[k]
                    - f_5 * sof1_528[k]
                    + f_3 * pc_y[k] * sog_792[k];

        t_1110[k] = f_10 * sng_658[k]
                    + f_6 * sof0_529[k]
                    - f_7 * sof1_529[k]
                    + f_3 * pc_y[k] * sog_793[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, t_1114, pb_y, pc_y, pc_z, snh0_924, sng_644, \
                         sng_659, sng_660, snh1_924, sof0_529, sof1_529, sog_794, \
                         sog_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_10 * sng_659[k]
                    + f_3 * pc_y[k] * sog_794[k];

        t_1112[k] = f_15 * sng_644[k]
                    + f_1 * sof0_529[k]
                    - f_2 * sof1_529[k]
                    + f_3 * pc_z[k] * sog_794[k];

        t_1113[k] = pb_y[k] * snh0_924[k]
                    - f_8 * pc_y[k] * snh1_924[k];

        t_1114[k] = f_9 * sng_660[k]
                    + f_3 * pc_y[k] * sog_795[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, t_1118, pb_y, pc_y, pc_z, snh0_927, snh0_929, \
                         sng_645, sng_661, sng_662, snh1_927, snh1_929, sog_795, \
                         sog_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = f_14 * sng_645[k]
                    + f_3 * pc_z[k] * sog_795[k];

        t_1116[k] = pb_y[k] * snh0_927[k]
                    + f_10 * sng_661[k]
                    - f_8 * pc_y[k] * snh1_927[k];

        t_1117[k] = f_9 * sng_662[k]
                    + f_3 * pc_y[k] * sog_797[k];

        t_1118[k] = pb_y[k] * snh0_929[k]
                    - f_8 * pc_y[k] * snh1_929[k];
    }

#pragma omp simd aligned(t_1119, t_1120, t_1121, t_1122, pb_y, pc_y, pc_z, snh0_930, snh0_933, \
                         sng_648, sng_663, sng_665, snh1_930, snh1_933, sog_798, \
                         sog_800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1119[k] = pb_y[k] * snh0_930[k]
                    + f_11 * sng_663[k]
                    - f_8 * pc_y[k] * snh1_930[k];

        t_1120[k] = f_14 * sng_648[k]
                    + f_3 * pc_z[k] * sog_798[k];

        t_1121[k] = f_9 * sng_665[k]
                    + f_3 * pc_y[k] * sog_800[k];

        t_1122[k] = pb_y[k] * snh0_933[k]
                    - f_8 * pc_y[k] * snh1_933[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, t_1126, t_1127, pc_x, sng_805, sng_806, \
                         sng_807, sng_808, sng_809, sog_805, sog_806, sog_807, sog_808, \
                         sog_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_10 * sng_805[k]
                    + f_3 * pc_x[k] * sog_805[k];

        t_1124[k] = f_10 * sng_806[k]
                    + f_3 * pc_x[k] * sog_806[k];

        t_1125[k] = f_10 * sng_807[k]
                    + f_3 * pc_x[k] * sog_807[k];

        t_1126[k] = f_10 * sng_808[k]
                    + f_3 * pc_x[k] * sog_808[k];

        t_1127[k] = f_10 * sng_809[k]
                    + f_3 * pc_x[k] * sog_809[k];
    }

#pragma omp simd aligned(t_1128, t_1129, t_1130, pc_y, pc_z, sng_655, sng_670, sng_672, \
                         sof0_536, sof0_538, sof1_536, sof1_538, sog_805, \
                         sog_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1128[k] = f_9 * sng_670[k]
                    + f_1 * sof0_536[k]
                    - f_2 * sof1_536[k]
                    + f_3 * pc_y[k] * sog_805[k];

        t_1129[k] = f_14 * sng_655[k]
                    + f_3 * pc_z[k] * sog_805[k];

        t_1130[k] = f_9 * sng_672[k]
                    + f_4 * sof0_538[k]
                    - f_5 * sof1_538[k]
                    + f_3 * pc_y[k] * sog_807[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, pb_y, pc_y, snh0_944, sng_673, sng_674, \
                         snh1_944, sof0_539, sof1_539, sog_808, \
                         sog_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = f_9 * sng_673[k]
                    + f_6 * sof0_539[k]
                    - f_7 * sof1_539[k]
                    + f_3 * pc_y[k] * sog_808[k];

        t_1132[k] = f_9 * sng_674[k]
                    + f_3 * pc_y[k] * sog_809[k];

        t_1133[k] = pb_y[k] * snh0_944[k]
                    - f_8 * pc_y[k] * snh1_944[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, pc_x, pc_y, pc_z, sng_660, sng_810, \
                         sng_813, sof0_540, sof0_543, sof1_540, sof1_543, sog_810, \
                         sog_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = f_10 * sng_810[k]
                    + f_1 * sof0_540[k]
                    - f_2 * sof1_540[k]
                    + f_3 * pc_x[k] * sog_810[k];

        t_1135[k] = f_3 * pc_y[k] * sog_810[k];

        t_1136[k] = f_13 * sng_660[k]
                    + f_3 * pc_z[k] * sog_810[k];

        t_1137[k] = f_10 * sng_813[k]
                    + f_4 * sof0_543[k]
                    - f_5 * sof1_543[k]
                    + f_3 * pc_x[k] * sog_813[k];
    }

#pragma omp simd aligned(t_1138, t_1139, t_1140, pc_x, pc_y, sng_815, sng_816, sof0_545, \
                         sof0_546, sof1_545, sof1_546, sog_812, sog_815, \
                         sog_816 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1138[k] = f_3 * pc_y[k] * sog_812[k];

        t_1139[k] = f_10 * sng_815[k]
                    + f_4 * sof0_545[k]
                    - f_5 * sof1_545[k]
                    + f_3 * pc_x[k] * sog_815[k];

        t_1140[k] = f_10 * sng_816[k]
                    + f_6 * sof0_546[k]
                    - f_7 * sof1_546[k]
                    + f_3 * pc_x[k] * sog_816[k];
    }

#pragma omp simd aligned(t_1141, t_1142, t_1143, t_1144, pc_x, pc_y, pc_z, sng_663, sng_819, \
                         sng_820, sof0_549, sof1_549, sog_813, sog_815, sog_819, \
                         sog_820 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1141[k] = f_13 * sng_663[k]
                    + f_3 * pc_z[k] * sog_813[k];

        t_1142[k] = f_3 * pc_y[k] * sog_815[k];

        t_1143[k] = f_10 * sng_819[k]
                    + f_6 * sof0_549[k]
                    - f_7 * sof1_549[k]
                    + f_3 * pc_x[k] * sog_819[k];

        t_1144[k] = f_10 * sng_820[k]
                    + f_3 * pc_x[k] * sog_820[k];
    }
}

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t snh0,
                                                           const size_t sng, const size_t snh1,
                                                           const size_t sof0, const size_t sof1,
                                                           const size_t sog, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 5.0 / q;
    const auto f_13 = 4.5 / q;
    const auto f_14 = 4.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_1145 = buffer.data(target + 1145);
    auto *t_1146 = buffer.data(target + 1146);
    auto *t_1147 = buffer.data(target + 1147);
    auto *t_1148 = buffer.data(target + 1148);
    auto *t_1149 = buffer.data(target + 1149);
    auto *t_1150 = buffer.data(target + 1150);
    auto *t_1151 = buffer.data(target + 1151);
    auto *t_1152 = buffer.data(target + 1152);
    auto *t_1153 = buffer.data(target + 1153);
    auto *t_1154 = buffer.data(target + 1154);
    auto *t_1155 = buffer.data(target + 1155);
    auto *t_1156 = buffer.data(target + 1156);
    auto *t_1157 = buffer.data(target + 1157);
    auto *t_1158 = buffer.data(target + 1158);
    auto *t_1159 = buffer.data(target + 1159);
    auto *t_1160 = buffer.data(target + 1160);
    auto *t_1161 = buffer.data(target + 1161);
    auto *t_1162 = buffer.data(target + 1162);
    auto *t_1163 = buffer.data(target + 1163);
    auto *t_1164 = buffer.data(target + 1164);
    auto *t_1165 = buffer.data(target + 1165);
    auto *t_1166 = buffer.data(target + 1166);
    auto *t_1167 = buffer.data(target + 1167);
    auto *t_1168 = buffer.data(target + 1168);
    auto *t_1169 = buffer.data(target + 1169);
    auto *t_1170 = buffer.data(target + 1170);
    auto *t_1171 = buffer.data(target + 1171);
    auto *t_1172 = buffer.data(target + 1172);
    auto *t_1173 = buffer.data(target + 1173);
    auto *t_1174 = buffer.data(target + 1174);
    auto *t_1175 = buffer.data(target + 1175);
    auto *t_1176 = buffer.data(target + 1176);
    auto *t_1177 = buffer.data(target + 1177);
    auto *t_1178 = buffer.data(target + 1178);
    auto *t_1179 = buffer.data(target + 1179);
    auto *t_1180 = buffer.data(target + 1180);
    auto *t_1181 = buffer.data(target + 1181);
    auto *t_1182 = buffer.data(target + 1182);
    auto *t_1183 = buffer.data(target + 1183);
    auto *t_1184 = buffer.data(target + 1184);
    auto *t_1185 = buffer.data(target + 1185);
    auto *t_1186 = buffer.data(target + 1186);
    auto *t_1187 = buffer.data(target + 1187);
    auto *t_1188 = buffer.data(target + 1188);
    auto *t_1189 = buffer.data(target + 1189);
    auto *t_1190 = buffer.data(target + 1190);
    auto *t_1191 = buffer.data(target + 1191);
    auto *t_1192 = buffer.data(target + 1192);
    auto *t_1193 = buffer.data(target + 1193);
    auto *t_1194 = buffer.data(target + 1194);
    auto *t_1195 = buffer.data(target + 1195);
    auto *t_1196 = buffer.data(target + 1196);
    auto *t_1197 = buffer.data(target + 1197);
    auto *t_1198 = buffer.data(target + 1198);
    auto *t_1199 = buffer.data(target + 1199);
    auto *t_1200 = buffer.data(target + 1200);
    auto *t_1201 = buffer.data(target + 1201);
    auto *t_1202 = buffer.data(target + 1202);
    auto *t_1203 = buffer.data(target + 1203);
    auto *t_1204 = buffer.data(target + 1204);
    auto *t_1205 = buffer.data(target + 1205);
    auto *t_1206 = buffer.data(target + 1206);
    auto *t_1207 = buffer.data(target + 1207);
    auto *t_1208 = buffer.data(target + 1208);
    auto *t_1209 = buffer.data(target + 1209);
    auto *t_1210 = buffer.data(target + 1210);
    auto *t_1211 = buffer.data(target + 1211);
    auto *t_1212 = buffer.data(target + 1212);
    auto *t_1213 = buffer.data(target + 1213);
    auto *t_1214 = buffer.data(target + 1214);
    auto *t_1215 = buffer.data(target + 1215);
    auto *t_1216 = buffer.data(target + 1216);
    auto *t_1217 = buffer.data(target + 1217);
    auto *t_1218 = buffer.data(target + 1218);
    auto *t_1219 = buffer.data(target + 1219);
    auto *t_1220 = buffer.data(target + 1220);
    auto *t_1221 = buffer.data(target + 1221);
    auto *t_1222 = buffer.data(target + 1222);
    auto *t_1223 = buffer.data(target + 1223);
    auto *t_1224 = buffer.data(target + 1224);
    auto *t_1225 = buffer.data(target + 1225);
    auto *t_1226 = buffer.data(target + 1226);
    auto *t_1227 = buffer.data(target + 1227);
    auto *t_1228 = buffer.data(target + 1228);
    auto *t_1229 = buffer.data(target + 1229);
    auto *t_1230 = buffer.data(target + 1230);
    auto *t_1231 = buffer.data(target + 1231);
    auto *t_1232 = buffer.data(target + 1232);
    auto *t_1233 = buffer.data(target + 1233);
    auto *t_1234 = buffer.data(target + 1234);
    auto *t_1235 = buffer.data(target + 1235);
    auto *t_1236 = buffer.data(target + 1236);
    auto *t_1237 = buffer.data(target + 1237);
    auto *t_1238 = buffer.data(target + 1238);
    auto *t_1239 = buffer.data(target + 1239);
    auto *t_1240 = buffer.data(target + 1240);
    auto *t_1241 = buffer.data(target + 1241);
    auto *t_1242 = buffer.data(target + 1242);
    auto *t_1243 = buffer.data(target + 1243);
    auto *t_1244 = buffer.data(target + 1244);
    auto *t_1245 = buffer.data(target + 1245);
    auto *t_1246 = buffer.data(target + 1246);
    auto *t_1247 = buffer.data(target + 1247);
    auto *t_1248 = buffer.data(target + 1248);
    auto *t_1249 = buffer.data(target + 1249);
    auto *t_1250 = buffer.data(target + 1250);
    auto *t_1251 = buffer.data(target + 1251);
    auto *t_1252 = buffer.data(target + 1252);
    auto *t_1253 = buffer.data(target + 1253);
    auto *t_1254 = buffer.data(target + 1254);
    auto *t_1255 = buffer.data(target + 1255);
    auto *t_1256 = buffer.data(target + 1256);
    auto *t_1257 = buffer.data(target + 1257);
    auto *t_1258 = buffer.data(target + 1258);
    auto *t_1259 = buffer.data(target + 1259);
    auto *t_1260 = buffer.data(target + 1260);
    auto *t_1261 = buffer.data(target + 1261);
    auto *t_1262 = buffer.data(target + 1262);
    auto *t_1263 = buffer.data(target + 1263);
    auto *t_1264 = buffer.data(target + 1264);
    auto *t_1265 = buffer.data(target + 1265);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_945 = buffer.data(snh0 + 945);
    const auto *snh0_948 = buffer.data(snh0 + 948);
    const auto *snh0_951 = buffer.data(snh0 + 951);
    const auto *snh0_1155 = buffer.data(snh0 + 1155);
    const auto *snh0_1158 = buffer.data(snh0 + 1158);
    const auto *snh0_1160 = buffer.data(snh0 + 1160);
    const auto *snh0_1161 = buffer.data(snh0 + 1161);
    const auto *snh0_1164 = buffer.data(snh0 + 1164);
    const auto *snh0_1170 = buffer.data(snh0 + 1170);
    const auto *snh0_1172 = buffer.data(snh0 + 1172);
    const auto *snh0_1173 = buffer.data(snh0 + 1173);
    const auto *snh0_1175 = buffer.data(snh0 + 1175);
    const auto *snh0_1181 = buffer.data(snh0 + 1181);
    const auto *snh0_1185 = buffer.data(snh0 + 1185);
    const auto *snh0_1191 = buffer.data(snh0 + 1191);
    const auto *snh0_1193 = buffer.data(snh0 + 1193);
    const auto *snh0_1194 = buffer.data(snh0 + 1194);
    const auto *snh0_1196 = buffer.data(snh0 + 1196);
    const auto *snh0_1197 = buffer.data(snh0 + 1197);
    const auto *snh0_1200 = buffer.data(snh0 + 1200);
    const auto *snh0_1202 = buffer.data(snh0 + 1202);
    const auto *snh0_1203 = buffer.data(snh0 + 1203);
    const auto *snh0_1206 = buffer.data(snh0 + 1206);
    const auto *snh0_1212 = buffer.data(snh0 + 1212);
    const auto *snh0_1214 = buffer.data(snh0 + 1214);
    const auto *snh0_1215 = buffer.data(snh0 + 1215);
    const auto *snh0_1217 = buffer.data(snh0 + 1217);
    const auto *snh0_1218 = buffer.data(snh0 + 1218);
    const auto *snh0_1221 = buffer.data(snh0 + 1221);
    const auto *snh0_1223 = buffer.data(snh0 + 1223);
    const auto *snh0_1224 = buffer.data(snh0 + 1224);
    const auto *snh0_1227 = buffer.data(snh0 + 1227);
    const auto *snh0_1233 = buffer.data(snh0 + 1233);
    const auto *snh0_1235 = buffer.data(snh0 + 1235);
    const auto *snh0_1236 = buffer.data(snh0 + 1236);
    const auto *snh0_1238 = buffer.data(snh0 + 1238);
    const auto *snh0_1239 = buffer.data(snh0 + 1239);
    const auto *snh0_1242 = buffer.data(snh0 + 1242);
    const auto *snh0_1244 = buffer.data(snh0 + 1244);
    const auto *snh0_1245 = buffer.data(snh0 + 1245);
    const auto *snh0_1248 = buffer.data(snh0 + 1248);
    const auto *snh0_1254 = buffer.data(snh0 + 1254);
    const auto *snh0_1256 = buffer.data(snh0 + 1256);
    const auto *snh0_1257 = buffer.data(snh0 + 1257);
    const auto *snh0_1259 = buffer.data(snh0 + 1259);
    const auto *snh0_1260 = buffer.data(snh0 + 1260);
    const auto *snh0_1263 = buffer.data(snh0 + 1263);
    const auto *snh0_1265 = buffer.data(snh0 + 1265);

    const auto *sng_670 = buffer.data(sng + 670);
    const auto *sng_674 = buffer.data(sng + 674);
    const auto *sng_675 = buffer.data(sng + 675);
    const auto *sng_677 = buffer.data(sng + 677);
    const auto *sng_678 = buffer.data(sng + 678);
    const auto *sng_680 = buffer.data(sng + 680);
    const auto *sng_685 = buffer.data(sng + 685);
    const auto *sng_689 = buffer.data(sng + 689);
    const auto *sng_690 = buffer.data(sng + 690);
    const auto *sng_692 = buffer.data(sng + 692);
    const auto *sng_693 = buffer.data(sng + 693);
    const auto *sng_695 = buffer.data(sng + 695);
    const auto *sng_700 = buffer.data(sng + 700);
    const auto *sng_704 = buffer.data(sng + 704);
    const auto *sng_705 = buffer.data(sng + 705);
    const auto *sng_707 = buffer.data(sng + 707);
    const auto *sng_708 = buffer.data(sng + 708);
    const auto *sng_710 = buffer.data(sng + 710);
    const auto *sng_715 = buffer.data(sng + 715);
    const auto *sng_719 = buffer.data(sng + 719);
    const auto *sng_720 = buffer.data(sng + 720);
    const auto *sng_722 = buffer.data(sng + 722);
    const auto *sng_723 = buffer.data(sng + 723);
    const auto *sng_725 = buffer.data(sng + 725);
    const auto *sng_730 = buffer.data(sng + 730);
    const auto *sng_734 = buffer.data(sng + 734);
    const auto *sng_735 = buffer.data(sng + 735);
    const auto *sng_737 = buffer.data(sng + 737);
    const auto *sng_740 = buffer.data(sng + 740);
    const auto *sng_749 = buffer.data(sng + 749);
    const auto *sng_750 = buffer.data(sng + 750);
    const auto *sng_752 = buffer.data(sng + 752);
    const auto *sng_821 = buffer.data(sng + 821);
    const auto *sng_822 = buffer.data(sng + 822);
    const auto *sng_823 = buffer.data(sng + 823);
    const auto *sng_824 = buffer.data(sng + 824);
    const auto *sng_825 = buffer.data(sng + 825);
    const auto *sng_828 = buffer.data(sng + 828);
    const auto *sng_830 = buffer.data(sng + 830);
    const auto *sng_831 = buffer.data(sng + 831);
    const auto *sng_834 = buffer.data(sng + 834);
    const auto *sng_835 = buffer.data(sng + 835);
    const auto *sng_836 = buffer.data(sng + 836);
    const auto *sng_837 = buffer.data(sng + 837);
    const auto *sng_838 = buffer.data(sng + 838);
    const auto *sng_839 = buffer.data(sng + 839);
    const auto *sng_845 = buffer.data(sng + 845);
    const auto *sng_849 = buffer.data(sng + 849);
    const auto *sng_850 = buffer.data(sng + 850);
    const auto *sng_851 = buffer.data(sng + 851);
    const auto *sng_852 = buffer.data(sng + 852);
    const auto *sng_853 = buffer.data(sng + 853);
    const auto *sng_854 = buffer.data(sng + 854);
    const auto *sng_855 = buffer.data(sng + 855);
    const auto *sng_858 = buffer.data(sng + 858);
    const auto *sng_860 = buffer.data(sng + 860);
    const auto *sng_861 = buffer.data(sng + 861);
    const auto *sng_864 = buffer.data(sng + 864);
    const auto *sng_865 = buffer.data(sng + 865);
    const auto *sng_866 = buffer.data(sng + 866);
    const auto *sng_867 = buffer.data(sng + 867);
    const auto *sng_868 = buffer.data(sng + 868);
    const auto *sng_869 = buffer.data(sng + 869);
    const auto *sng_870 = buffer.data(sng + 870);
    const auto *sng_873 = buffer.data(sng + 873);
    const auto *sng_875 = buffer.data(sng + 875);
    const auto *sng_876 = buffer.data(sng + 876);
    const auto *sng_879 = buffer.data(sng + 879);
    const auto *sng_880 = buffer.data(sng + 880);
    const auto *sng_881 = buffer.data(sng + 881);
    const auto *sng_882 = buffer.data(sng + 882);
    const auto *sng_883 = buffer.data(sng + 883);
    const auto *sng_884 = buffer.data(sng + 884);
    const auto *sng_885 = buffer.data(sng + 885);
    const auto *sng_888 = buffer.data(sng + 888);
    const auto *sng_890 = buffer.data(sng + 890);
    const auto *sng_891 = buffer.data(sng + 891);
    const auto *sng_894 = buffer.data(sng + 894);
    const auto *sng_895 = buffer.data(sng + 895);
    const auto *sng_896 = buffer.data(sng + 896);
    const auto *sng_897 = buffer.data(sng + 897);
    const auto *sng_898 = buffer.data(sng + 898);
    const auto *sng_899 = buffer.data(sng + 899);
    const auto *sng_900 = buffer.data(sng + 900);
    const auto *sng_903 = buffer.data(sng + 903);
    const auto *sng_905 = buffer.data(sng + 905);

    const auto *snh1_945 = buffer.data(snh1 + 945);
    const auto *snh1_948 = buffer.data(snh1 + 948);
    const auto *snh1_951 = buffer.data(snh1 + 951);
    const auto *snh1_1155 = buffer.data(snh1 + 1155);
    const auto *snh1_1158 = buffer.data(snh1 + 1158);
    const auto *snh1_1160 = buffer.data(snh1 + 1160);
    const auto *snh1_1161 = buffer.data(snh1 + 1161);
    const auto *snh1_1164 = buffer.data(snh1 + 1164);
    const auto *snh1_1170 = buffer.data(snh1 + 1170);
    const auto *snh1_1172 = buffer.data(snh1 + 1172);
    const auto *snh1_1173 = buffer.data(snh1 + 1173);
    const auto *snh1_1175 = buffer.data(snh1 + 1175);
    const auto *snh1_1181 = buffer.data(snh1 + 1181);
    const auto *snh1_1185 = buffer.data(snh1 + 1185);
    const auto *snh1_1191 = buffer.data(snh1 + 1191);
    const auto *snh1_1193 = buffer.data(snh1 + 1193);
    const auto *snh1_1194 = buffer.data(snh1 + 1194);
    const auto *snh1_1196 = buffer.data(snh1 + 1196);
    const auto *snh1_1197 = buffer.data(snh1 + 1197);
    const auto *snh1_1200 = buffer.data(snh1 + 1200);
    const auto *snh1_1202 = buffer.data(snh1 + 1202);
    const auto *snh1_1203 = buffer.data(snh1 + 1203);
    const auto *snh1_1206 = buffer.data(snh1 + 1206);
    const auto *snh1_1212 = buffer.data(snh1 + 1212);
    const auto *snh1_1214 = buffer.data(snh1 + 1214);
    const auto *snh1_1215 = buffer.data(snh1 + 1215);
    const auto *snh1_1217 = buffer.data(snh1 + 1217);
    const auto *snh1_1218 = buffer.data(snh1 + 1218);
    const auto *snh1_1221 = buffer.data(snh1 + 1221);
    const auto *snh1_1223 = buffer.data(snh1 + 1223);
    const auto *snh1_1224 = buffer.data(snh1 + 1224);
    const auto *snh1_1227 = buffer.data(snh1 + 1227);
    const auto *snh1_1233 = buffer.data(snh1 + 1233);
    const auto *snh1_1235 = buffer.data(snh1 + 1235);
    const auto *snh1_1236 = buffer.data(snh1 + 1236);
    const auto *snh1_1238 = buffer.data(snh1 + 1238);
    const auto *snh1_1239 = buffer.data(snh1 + 1239);
    const auto *snh1_1242 = buffer.data(snh1 + 1242);
    const auto *snh1_1244 = buffer.data(snh1 + 1244);
    const auto *snh1_1245 = buffer.data(snh1 + 1245);
    const auto *snh1_1248 = buffer.data(snh1 + 1248);
    const auto *snh1_1254 = buffer.data(snh1 + 1254);
    const auto *snh1_1256 = buffer.data(snh1 + 1256);
    const auto *snh1_1257 = buffer.data(snh1 + 1257);
    const auto *snh1_1259 = buffer.data(snh1 + 1259);
    const auto *snh1_1260 = buffer.data(snh1 + 1260);
    const auto *snh1_1263 = buffer.data(snh1 + 1263);
    const auto *snh1_1265 = buffer.data(snh1 + 1265);

    const auto *sof0_546 = buffer.data(sof0 + 546);
    const auto *sof0_548 = buffer.data(sof0 + 548);
    const auto *sof0_549 = buffer.data(sof0 + 549);

    const auto *sof1_546 = buffer.data(sof1 + 546);
    const auto *sof1_548 = buffer.data(sof1 + 548);
    const auto *sof1_549 = buffer.data(sof1 + 549);

    const auto *sog_820 = buffer.data(sog + 820);
    const auto *sog_821 = buffer.data(sog + 821);
    const auto *sog_822 = buffer.data(sog + 822);
    const auto *sog_823 = buffer.data(sog + 823);
    const auto *sog_824 = buffer.data(sog + 824);
    const auto *sog_825 = buffer.data(sog + 825);
    const auto *sog_827 = buffer.data(sog + 827);
    const auto *sog_828 = buffer.data(sog + 828);
    const auto *sog_830 = buffer.data(sog + 830);
    const auto *sog_835 = buffer.data(sog + 835);
    const auto *sog_836 = buffer.data(sog + 836);
    const auto *sog_837 = buffer.data(sog + 837);
    const auto *sog_838 = buffer.data(sog + 838);
    const auto *sog_839 = buffer.data(sog + 839);
    const auto *sog_840 = buffer.data(sog + 840);
    const auto *sog_842 = buffer.data(sog + 842);
    const auto *sog_843 = buffer.data(sog + 843);
    const auto *sog_845 = buffer.data(sog + 845);
    const auto *sog_850 = buffer.data(sog + 850);
    const auto *sog_851 = buffer.data(sog + 851);
    const auto *sog_852 = buffer.data(sog + 852);
    const auto *sog_853 = buffer.data(sog + 853);
    const auto *sog_854 = buffer.data(sog + 854);
    const auto *sog_855 = buffer.data(sog + 855);
    const auto *sog_857 = buffer.data(sog + 857);
    const auto *sog_858 = buffer.data(sog + 858);
    const auto *sog_860 = buffer.data(sog + 860);
    const auto *sog_865 = buffer.data(sog + 865);
    const auto *sog_866 = buffer.data(sog + 866);
    const auto *sog_867 = buffer.data(sog + 867);
    const auto *sog_868 = buffer.data(sog + 868);
    const auto *sog_869 = buffer.data(sog + 869);
    const auto *sog_870 = buffer.data(sog + 870);
    const auto *sog_872 = buffer.data(sog + 872);
    const auto *sog_873 = buffer.data(sog + 873);
    const auto *sog_875 = buffer.data(sog + 875);
    const auto *sog_880 = buffer.data(sog + 880);
    const auto *sog_881 = buffer.data(sog + 881);
    const auto *sog_882 = buffer.data(sog + 882);
    const auto *sog_883 = buffer.data(sog + 883);
    const auto *sog_884 = buffer.data(sog + 884);
    const auto *sog_885 = buffer.data(sog + 885);
    const auto *sog_887 = buffer.data(sog + 887);
    const auto *sog_888 = buffer.data(sog + 888);
    const auto *sog_890 = buffer.data(sog + 890);
    const auto *sog_895 = buffer.data(sog + 895);
    const auto *sog_896 = buffer.data(sog + 896);
    const auto *sog_897 = buffer.data(sog + 897);
    const auto *sog_898 = buffer.data(sog + 898);
    const auto *sog_899 = buffer.data(sog + 899);
    const auto *sog_900 = buffer.data(sog + 900);
    const auto *sog_902 = buffer.data(sog + 902);

#pragma omp simd aligned(t_1145, t_1146, t_1147, t_1148, pc_x, sng_821, sng_822, sng_823, \
                         sng_824, sog_821, sog_822, sog_823, sog_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = f_10 * sng_821[k]
                    + f_3 * pc_x[k] * sog_821[k];

        t_1146[k] = f_10 * sng_822[k]
                    + f_3 * pc_x[k] * sog_822[k];

        t_1147[k] = f_10 * sng_823[k]
                    + f_3 * pc_x[k] * sog_823[k];

        t_1148[k] = f_10 * sng_824[k]
                    + f_3 * pc_x[k] * sog_824[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, t_1152, pc_y, pc_z, sng_670, sof0_546, \
                         sof0_548, sof0_549, sof1_546, sof1_548, sof1_549, sog_820, sog_822, \
                         sog_823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_1 * sof0_546[k]
                    - f_2 * sof1_546[k]
                    + f_3 * pc_y[k] * sog_820[k];

        t_1150[k] = f_13 * sng_670[k]
                    + f_3 * pc_z[k] * sog_820[k];

        t_1151[k] = f_4 * sof0_548[k]
                    - f_5 * sof1_548[k]
                    + f_3 * pc_y[k] * sog_822[k];

        t_1152[k] = f_6 * sof0_549[k]
                    - f_7 * sof1_549[k]
                    + f_3 * pc_y[k] * sog_823[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, pb_x, pc_x, pc_y, pc_z, snh0_1155, sng_674, \
                         sng_825, snh1_1155, sof0_549, sof1_549, \
                         sog_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_3 * pc_y[k] * sog_824[k];

        t_1154[k] = f_13 * sng_674[k]
                    + f_1 * sof0_549[k]
                    - f_2 * sof1_549[k]
                    + f_3 * pc_z[k] * sog_824[k];

        t_1155[k] = pb_x[k] * snh0_1155[k]
                    + f_18 * sng_825[k]
                    - f_8 * pc_x[k] * snh1_1155[k];
    }

#pragma omp simd aligned(t_1156, t_1157, t_1158, t_1159, pb_x, pc_x, pc_y, pc_z, snh0_1158, \
                         sng_675, sng_677, sng_828, snh1_1158, sog_825, \
                         sog_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1156[k] = f_12 * sng_675[k]
                    + f_3 * pc_y[k] * sog_825[k];

        t_1157[k] = f_3 * pc_z[k] * sog_825[k];

        t_1158[k] = pb_x[k] * snh0_1158[k]
                    + f_11 * sng_828[k]
                    - f_8 * pc_x[k] * snh1_1158[k];

        t_1159[k] = f_12 * sng_677[k]
                    + f_3 * pc_y[k] * sog_827[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, pb_x, pc_x, pc_z, snh0_1160, snh0_1161, \
                         sng_830, sng_831, snh1_1160, snh1_1161, \
                         sog_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = pb_x[k] * snh0_1160[k]
                    + f_11 * sng_830[k]
                    - f_8 * pc_x[k] * snh1_1160[k];

        t_1161[k] = pb_x[k] * snh0_1161[k]
                    + f_10 * sng_831[k]
                    - f_8 * pc_x[k] * snh1_1161[k];

        t_1162[k] = f_3 * pc_z[k] * sog_828[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, t_1166, pb_x, pc_x, pc_y, snh0_1164, sng_680, \
                         sng_834, sng_835, sng_836, snh1_1164, sog_830, sog_835, \
                         sog_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_12 * sng_680[k]
                    + f_3 * pc_y[k] * sog_830[k];

        t_1164[k] = pb_x[k] * snh0_1164[k]
                    + f_10 * sng_834[k]
                    - f_8 * pc_x[k] * snh1_1164[k];

        t_1165[k] = f_9 * sng_835[k]
                    + f_3 * pc_x[k] * sog_835[k];

        t_1166[k] = f_9 * sng_836[k]
                    + f_3 * pc_x[k] * sog_836[k];
    }

#pragma omp simd aligned(t_1167, t_1168, t_1169, t_1170, pb_x, pc_x, snh0_1170, sng_837, \
                         sng_838, sng_839, snh1_1170, sog_837, sog_838, \
                         sog_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1167[k] = f_9 * sng_837[k]
                    + f_3 * pc_x[k] * sog_837[k];

        t_1168[k] = f_9 * sng_838[k]
                    + f_3 * pc_x[k] * sog_838[k];

        t_1169[k] = f_9 * sng_839[k]
                    + f_3 * pc_x[k] * sog_839[k];

        t_1170[k] = pb_x[k] * snh0_1170[k]
                    - f_8 * pc_x[k] * snh1_1170[k];
    }

#pragma omp simd aligned(t_1171, t_1172, t_1173, t_1174, pb_x, pc_x, pc_y, pc_z, snh0_1172, \
                         snh0_1173, sng_689, snh1_1172, snh1_1173, sog_835, \
                         sog_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1171[k] = f_3 * pc_z[k] * sog_835[k];

        t_1172[k] = pb_x[k] * snh0_1172[k]
                    - f_8 * pc_x[k] * snh1_1172[k];

        t_1173[k] = pb_x[k] * snh0_1173[k]
                    - f_8 * pc_x[k] * snh1_1173[k];

        t_1174[k] = f_12 * sng_689[k]
                    + f_3 * pc_y[k] * sog_839[k];
    }

#pragma omp simd aligned(t_1175, t_1176, t_1177, t_1178, pb_x, pb_z, pc_x, pc_y, pc_z, \
                         snh0_945, snh0_1175, sng_675, sng_690, snh1_945, snh1_1175, \
                         sog_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1175[k] = pb_x[k] * snh0_1175[k]
                    - f_8 * pc_x[k] * snh1_1175[k];

        t_1176[k] = pb_z[k] * snh0_945[k]
                    - f_8 * pc_z[k] * snh1_945[k];

        t_1177[k] = f_13 * sng_690[k]
                    + f_3 * pc_y[k] * sog_840[k];

        t_1178[k] = f_9 * sng_675[k]
                    + f_3 * pc_z[k] * sog_840[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, pb_x, pb_z, pc_x, pc_y, pc_z, snh0_948, \
                         snh0_1181, sng_692, sng_845, snh1_948, snh1_1181, \
                         sog_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = pb_z[k] * snh0_948[k]
                    - f_8 * pc_z[k] * snh1_948[k];

        t_1180[k] = f_13 * sng_692[k]
                    + f_3 * pc_y[k] * sog_842[k];

        t_1181[k] = pb_x[k] * snh0_1181[k]
                    + f_11 * sng_845[k]
                    - f_8 * pc_x[k] * snh1_1181[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, pb_z, pc_y, pc_z, snh0_951, sng_678, sng_695, \
                         snh1_951, sog_843, sog_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = pb_z[k] * snh0_951[k]
                    - f_8 * pc_z[k] * snh1_951[k];

        t_1183[k] = f_9 * sng_678[k]
                    + f_3 * pc_z[k] * sog_843[k];

        t_1184[k] = f_13 * sng_695[k]
                    + f_3 * pc_y[k] * sog_845[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, t_1188, pb_x, pc_x, snh0_1185, sng_849, \
                         sng_850, sng_851, sng_852, snh1_1185, sog_850, sog_851, \
                         sog_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = pb_x[k] * snh0_1185[k]
                    + f_10 * sng_849[k]
                    - f_8 * pc_x[k] * snh1_1185[k];

        t_1186[k] = f_9 * sng_850[k]
                    + f_3 * pc_x[k] * sog_850[k];

        t_1187[k] = f_9 * sng_851[k]
                    + f_3 * pc_x[k] * sog_851[k];

        t_1188[k] = f_9 * sng_852[k]
                    + f_3 * pc_x[k] * sog_852[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, t_1192, pb_x, pc_x, pc_z, snh0_1191, sng_685, \
                         sng_853, sng_854, snh1_1191, sog_850, sog_853, \
                         sog_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_9 * sng_853[k]
                    + f_3 * pc_x[k] * sog_853[k];

        t_1190[k] = f_9 * sng_854[k]
                    + f_3 * pc_x[k] * sog_854[k];

        t_1191[k] = pb_x[k] * snh0_1191[k]
                    - f_8 * pc_x[k] * snh1_1191[k];

        t_1192[k] = f_9 * sng_685[k]
                    + f_3 * pc_z[k] * sog_850[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, t_1196, pb_x, pc_x, pc_y, snh0_1193, \
                         snh0_1194, snh0_1196, sng_704, snh1_1193, snh1_1194, snh1_1196, \
                         sog_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = pb_x[k] * snh0_1193[k]
                    - f_8 * pc_x[k] * snh1_1193[k];

        t_1194[k] = pb_x[k] * snh0_1194[k]
                    - f_8 * pc_x[k] * snh1_1194[k];

        t_1195[k] = f_13 * sng_704[k]
                    + f_3 * pc_y[k] * sog_854[k];

        t_1196[k] = pb_x[k] * snh0_1196[k]
                    - f_8 * pc_x[k] * snh1_1196[k];
    }

#pragma omp simd aligned(t_1197, t_1198, t_1199, pb_x, pc_x, pc_y, pc_z, snh0_1197, sng_690, \
                         sng_705, sng_855, snh1_1197, sog_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1197[k] = pb_x[k] * snh0_1197[k]
                    + f_18 * sng_855[k]
                    - f_8 * pc_x[k] * snh1_1197[k];

        t_1198[k] = f_14 * sng_705[k]
                    + f_3 * pc_y[k] * sog_855[k];

        t_1199[k] = f_10 * sng_690[k]
                    + f_3 * pc_z[k] * sog_855[k];
    }

#pragma omp simd aligned(t_1200, t_1201, t_1202, pb_x, pc_x, pc_y, snh0_1200, snh0_1202, \
                         sng_707, sng_858, sng_860, snh1_1200, snh1_1202, \
                         sog_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1200[k] = pb_x[k] * snh0_1200[k]
                    + f_11 * sng_858[k]
                    - f_8 * pc_x[k] * snh1_1200[k];

        t_1201[k] = f_14 * sng_707[k]
                    + f_3 * pc_y[k] * sog_857[k];

        t_1202[k] = pb_x[k] * snh0_1202[k]
                    + f_11 * sng_860[k]
                    - f_8 * pc_x[k] * snh1_1202[k];
    }

#pragma omp simd aligned(t_1203, t_1204, t_1205, pb_x, pc_x, pc_y, pc_z, snh0_1203, sng_693, \
                         sng_710, sng_861, snh1_1203, sog_858, \
                         sog_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1203[k] = pb_x[k] * snh0_1203[k]
                    + f_10 * sng_861[k]
                    - f_8 * pc_x[k] * snh1_1203[k];

        t_1204[k] = f_10 * sng_693[k]
                    + f_3 * pc_z[k] * sog_858[k];

        t_1205[k] = f_14 * sng_710[k]
                    + f_3 * pc_y[k] * sog_860[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, t_1209, pb_x, pc_x, snh0_1206, sng_864, \
                         sng_865, sng_866, sng_867, snh1_1206, sog_865, sog_866, \
                         sog_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = pb_x[k] * snh0_1206[k]
                    + f_10 * sng_864[k]
                    - f_8 * pc_x[k] * snh1_1206[k];

        t_1207[k] = f_9 * sng_865[k]
                    + f_3 * pc_x[k] * sog_865[k];

        t_1208[k] = f_9 * sng_866[k]
                    + f_3 * pc_x[k] * sog_866[k];

        t_1209[k] = f_9 * sng_867[k]
                    + f_3 * pc_x[k] * sog_867[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, pb_x, pc_x, pc_z, snh0_1212, sng_700, \
                         sng_868, sng_869, snh1_1212, sog_865, sog_868, \
                         sog_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = f_9 * sng_868[k]
                    + f_3 * pc_x[k] * sog_868[k];

        t_1211[k] = f_9 * sng_869[k]
                    + f_3 * pc_x[k] * sog_869[k];

        t_1212[k] = pb_x[k] * snh0_1212[k]
                    - f_8 * pc_x[k] * snh1_1212[k];

        t_1213[k] = f_10 * sng_700[k]
                    + f_3 * pc_z[k] * sog_865[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, t_1217, pb_x, pc_x, pc_y, snh0_1214, \
                         snh0_1215, snh0_1217, sng_719, snh1_1214, snh1_1215, snh1_1217, \
                         sog_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = pb_x[k] * snh0_1214[k]
                    - f_8 * pc_x[k] * snh1_1214[k];

        t_1215[k] = pb_x[k] * snh0_1215[k]
                    - f_8 * pc_x[k] * snh1_1215[k];

        t_1216[k] = f_14 * sng_719[k]
                    + f_3 * pc_y[k] * sog_869[k];

        t_1217[k] = pb_x[k] * snh0_1217[k]
                    - f_8 * pc_x[k] * snh1_1217[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, pb_x, pc_x, pc_y, pc_z, snh0_1218, sng_705, \
                         sng_720, sng_870, snh1_1218, sog_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = pb_x[k] * snh0_1218[k]
                    + f_18 * sng_870[k]
                    - f_8 * pc_x[k] * snh1_1218[k];

        t_1219[k] = f_15 * sng_720[k]
                    + f_3 * pc_y[k] * sog_870[k];

        t_1220[k] = f_11 * sng_705[k]
                    + f_3 * pc_z[k] * sog_870[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, pb_x, pc_x, pc_y, snh0_1221, snh0_1223, \
                         sng_722, sng_873, sng_875, snh1_1221, snh1_1223, \
                         sog_872 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = pb_x[k] * snh0_1221[k]
                    + f_11 * sng_873[k]
                    - f_8 * pc_x[k] * snh1_1221[k];

        t_1222[k] = f_15 * sng_722[k]
                    + f_3 * pc_y[k] * sog_872[k];

        t_1223[k] = pb_x[k] * snh0_1223[k]
                    + f_11 * sng_875[k]
                    - f_8 * pc_x[k] * snh1_1223[k];
    }

#pragma omp simd aligned(t_1224, t_1225, t_1226, pb_x, pc_x, pc_y, pc_z, snh0_1224, sng_708, \
                         sng_725, sng_876, snh1_1224, sog_873, \
                         sog_875 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1224[k] = pb_x[k] * snh0_1224[k]
                    + f_10 * sng_876[k]
                    - f_8 * pc_x[k] * snh1_1224[k];

        t_1225[k] = f_11 * sng_708[k]
                    + f_3 * pc_z[k] * sog_873[k];

        t_1226[k] = f_15 * sng_725[k]
                    + f_3 * pc_y[k] * sog_875[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, t_1230, pb_x, pc_x, snh0_1227, sng_879, \
                         sng_880, sng_881, sng_882, snh1_1227, sog_880, sog_881, \
                         sog_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = pb_x[k] * snh0_1227[k]
                    + f_10 * sng_879[k]
                    - f_8 * pc_x[k] * snh1_1227[k];

        t_1228[k] = f_9 * sng_880[k]
                    + f_3 * pc_x[k] * sog_880[k];

        t_1229[k] = f_9 * sng_881[k]
                    + f_3 * pc_x[k] * sog_881[k];

        t_1230[k] = f_9 * sng_882[k]
                    + f_3 * pc_x[k] * sog_882[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, t_1234, pb_x, pc_x, pc_z, snh0_1233, sng_715, \
                         sng_883, sng_884, snh1_1233, sog_880, sog_883, \
                         sog_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = f_9 * sng_883[k]
                    + f_3 * pc_x[k] * sog_883[k];

        t_1232[k] = f_9 * sng_884[k]
                    + f_3 * pc_x[k] * sog_884[k];

        t_1233[k] = pb_x[k] * snh0_1233[k]
                    - f_8 * pc_x[k] * snh1_1233[k];

        t_1234[k] = f_11 * sng_715[k]
                    + f_3 * pc_z[k] * sog_880[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, t_1238, pb_x, pc_x, pc_y, snh0_1235, \
                         snh0_1236, snh0_1238, sng_734, snh1_1235, snh1_1236, snh1_1238, \
                         sog_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = pb_x[k] * snh0_1235[k]
                    - f_8 * pc_x[k] * snh1_1235[k];

        t_1236[k] = pb_x[k] * snh0_1236[k]
                    - f_8 * pc_x[k] * snh1_1236[k];

        t_1237[k] = f_15 * sng_734[k]
                    + f_3 * pc_y[k] * sog_884[k];

        t_1238[k] = pb_x[k] * snh0_1238[k]
                    - f_8 * pc_x[k] * snh1_1238[k];
    }

#pragma omp simd aligned(t_1239, t_1240, t_1241, pb_x, pc_x, pc_y, pc_z, snh0_1239, sng_720, \
                         sng_735, sng_885, snh1_1239, sog_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1239[k] = pb_x[k] * snh0_1239[k]
                    + f_18 * sng_885[k]
                    - f_8 * pc_x[k] * snh1_1239[k];

        t_1240[k] = f_17 * sng_735[k]
                    + f_3 * pc_y[k] * sog_885[k];

        t_1241[k] = f_16 * sng_720[k]
                    + f_3 * pc_z[k] * sog_885[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, pb_x, pc_x, pc_y, snh0_1242, snh0_1244, \
                         sng_737, sng_888, sng_890, snh1_1242, snh1_1244, \
                         sog_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = pb_x[k] * snh0_1242[k]
                    + f_11 * sng_888[k]
                    - f_8 * pc_x[k] * snh1_1242[k];

        t_1243[k] = f_17 * sng_737[k]
                    + f_3 * pc_y[k] * sog_887[k];

        t_1244[k] = pb_x[k] * snh0_1244[k]
                    + f_11 * sng_890[k]
                    - f_8 * pc_x[k] * snh1_1244[k];
    }

#pragma omp simd aligned(t_1245, t_1246, t_1247, pb_x, pc_x, pc_y, pc_z, snh0_1245, sng_723, \
                         sng_740, sng_891, snh1_1245, sog_888, \
                         sog_890 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1245[k] = pb_x[k] * snh0_1245[k]
                    + f_10 * sng_891[k]
                    - f_8 * pc_x[k] * snh1_1245[k];

        t_1246[k] = f_16 * sng_723[k]
                    + f_3 * pc_z[k] * sog_888[k];

        t_1247[k] = f_17 * sng_740[k]
                    + f_3 * pc_y[k] * sog_890[k];
    }

#pragma omp simd aligned(t_1248, t_1249, t_1250, t_1251, pb_x, pc_x, snh0_1248, sng_894, \
                         sng_895, sng_896, sng_897, snh1_1248, sog_895, sog_896, \
                         sog_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1248[k] = pb_x[k] * snh0_1248[k]
                    + f_10 * sng_894[k]
                    - f_8 * pc_x[k] * snh1_1248[k];

        t_1249[k] = f_9 * sng_895[k]
                    + f_3 * pc_x[k] * sog_895[k];

        t_1250[k] = f_9 * sng_896[k]
                    + f_3 * pc_x[k] * sog_896[k];

        t_1251[k] = f_9 * sng_897[k]
                    + f_3 * pc_x[k] * sog_897[k];
    }

#pragma omp simd aligned(t_1252, t_1253, t_1254, t_1255, pb_x, pc_x, pc_z, snh0_1254, sng_730, \
                         sng_898, sng_899, snh1_1254, sog_895, sog_898, \
                         sog_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1252[k] = f_9 * sng_898[k]
                    + f_3 * pc_x[k] * sog_898[k];

        t_1253[k] = f_9 * sng_899[k]
                    + f_3 * pc_x[k] * sog_899[k];

        t_1254[k] = pb_x[k] * snh0_1254[k]
                    - f_8 * pc_x[k] * snh1_1254[k];

        t_1255[k] = f_16 * sng_730[k]
                    + f_3 * pc_z[k] * sog_895[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, t_1259, pb_x, pc_x, pc_y, snh0_1256, \
                         snh0_1257, snh0_1259, sng_749, snh1_1256, snh1_1257, snh1_1259, \
                         sog_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = pb_x[k] * snh0_1256[k]
                    - f_8 * pc_x[k] * snh1_1256[k];

        t_1257[k] = pb_x[k] * snh0_1257[k]
                    - f_8 * pc_x[k] * snh1_1257[k];

        t_1258[k] = f_17 * sng_749[k]
                    + f_3 * pc_y[k] * sog_899[k];

        t_1259[k] = pb_x[k] * snh0_1259[k]
                    - f_8 * pc_x[k] * snh1_1259[k];
    }

#pragma omp simd aligned(t_1260, t_1261, t_1262, pb_x, pc_x, pc_y, pc_z, snh0_1260, sng_735, \
                         sng_750, sng_900, snh1_1260, sog_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = pb_x[k] * snh0_1260[k]
                    + f_18 * sng_900[k]
                    - f_8 * pc_x[k] * snh1_1260[k];

        t_1261[k] = f_18 * sng_750[k]
                    + f_3 * pc_y[k] * sog_900[k];

        t_1262[k] = f_18 * sng_735[k]
                    + f_3 * pc_z[k] * sog_900[k];
    }

#pragma omp simd aligned(t_1263, t_1264, t_1265, pb_x, pc_x, pc_y, snh0_1263, snh0_1265, \
                         sng_752, sng_903, sng_905, snh1_1263, snh1_1265, \
                         sog_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1263[k] = pb_x[k] * snh0_1263[k]
                    + f_11 * sng_903[k]
                    - f_8 * pc_x[k] * snh1_1263[k];

        t_1264[k] = f_18 * sng_752[k]
                    + f_3 * pc_y[k] * sog_902[k];

        t_1265[k] = pb_x[k] * snh0_1265[k]
                    + f_11 * sng_905[k]
                    - f_8 * pc_x[k] * snh1_1265[k];
    }
}

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t snh0,
                                                           const size_t sng, const size_t snh1,
                                                           const size_t sof0, const size_t sof1,
                                                           const size_t sog, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 5.0 / q;
    const auto f_13 = 4.5 / q;
    const auto f_14 = 4.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_1266 = buffer.data(target + 1266);
    auto *t_1267 = buffer.data(target + 1267);
    auto *t_1268 = buffer.data(target + 1268);
    auto *t_1269 = buffer.data(target + 1269);
    auto *t_1270 = buffer.data(target + 1270);
    auto *t_1271 = buffer.data(target + 1271);
    auto *t_1272 = buffer.data(target + 1272);
    auto *t_1273 = buffer.data(target + 1273);
    auto *t_1274 = buffer.data(target + 1274);
    auto *t_1275 = buffer.data(target + 1275);
    auto *t_1276 = buffer.data(target + 1276);
    auto *t_1277 = buffer.data(target + 1277);
    auto *t_1278 = buffer.data(target + 1278);
    auto *t_1279 = buffer.data(target + 1279);
    auto *t_1280 = buffer.data(target + 1280);
    auto *t_1281 = buffer.data(target + 1281);
    auto *t_1282 = buffer.data(target + 1282);
    auto *t_1283 = buffer.data(target + 1283);
    auto *t_1284 = buffer.data(target + 1284);
    auto *t_1285 = buffer.data(target + 1285);
    auto *t_1286 = buffer.data(target + 1286);
    auto *t_1287 = buffer.data(target + 1287);
    auto *t_1288 = buffer.data(target + 1288);
    auto *t_1289 = buffer.data(target + 1289);
    auto *t_1290 = buffer.data(target + 1290);
    auto *t_1291 = buffer.data(target + 1291);
    auto *t_1292 = buffer.data(target + 1292);
    auto *t_1293 = buffer.data(target + 1293);
    auto *t_1294 = buffer.data(target + 1294);
    auto *t_1295 = buffer.data(target + 1295);
    auto *t_1296 = buffer.data(target + 1296);
    auto *t_1297 = buffer.data(target + 1297);
    auto *t_1298 = buffer.data(target + 1298);
    auto *t_1299 = buffer.data(target + 1299);
    auto *t_1300 = buffer.data(target + 1300);
    auto *t_1301 = buffer.data(target + 1301);
    auto *t_1302 = buffer.data(target + 1302);
    auto *t_1303 = buffer.data(target + 1303);
    auto *t_1304 = buffer.data(target + 1304);
    auto *t_1305 = buffer.data(target + 1305);
    auto *t_1306 = buffer.data(target + 1306);
    auto *t_1307 = buffer.data(target + 1307);
    auto *t_1308 = buffer.data(target + 1308);
    auto *t_1309 = buffer.data(target + 1309);
    auto *t_1310 = buffer.data(target + 1310);
    auto *t_1311 = buffer.data(target + 1311);
    auto *t_1312 = buffer.data(target + 1312);
    auto *t_1313 = buffer.data(target + 1313);
    auto *t_1314 = buffer.data(target + 1314);
    auto *t_1315 = buffer.data(target + 1315);
    auto *t_1316 = buffer.data(target + 1316);
    auto *t_1317 = buffer.data(target + 1317);
    auto *t_1318 = buffer.data(target + 1318);
    auto *t_1319 = buffer.data(target + 1319);
    auto *t_1320 = buffer.data(target + 1320);
    auto *t_1321 = buffer.data(target + 1321);
    auto *t_1322 = buffer.data(target + 1322);
    auto *t_1323 = buffer.data(target + 1323);
    auto *t_1324 = buffer.data(target + 1324);
    auto *t_1325 = buffer.data(target + 1325);
    auto *t_1326 = buffer.data(target + 1326);
    auto *t_1327 = buffer.data(target + 1327);
    auto *t_1328 = buffer.data(target + 1328);
    auto *t_1329 = buffer.data(target + 1329);
    auto *t_1330 = buffer.data(target + 1330);
    auto *t_1331 = buffer.data(target + 1331);
    auto *t_1332 = buffer.data(target + 1332);
    auto *t_1333 = buffer.data(target + 1333);
    auto *t_1334 = buffer.data(target + 1334);
    auto *t_1335 = buffer.data(target + 1335);
    auto *t_1336 = buffer.data(target + 1336);
    auto *t_1337 = buffer.data(target + 1337);
    auto *t_1338 = buffer.data(target + 1338);
    auto *t_1339 = buffer.data(target + 1339);
    auto *t_1340 = buffer.data(target + 1340);
    auto *t_1341 = buffer.data(target + 1341);
    auto *t_1342 = buffer.data(target + 1342);
    auto *t_1343 = buffer.data(target + 1343);
    auto *t_1344 = buffer.data(target + 1344);
    auto *t_1345 = buffer.data(target + 1345);
    auto *t_1346 = buffer.data(target + 1346);
    auto *t_1347 = buffer.data(target + 1347);
    auto *t_1348 = buffer.data(target + 1348);
    auto *t_1349 = buffer.data(target + 1349);
    auto *t_1350 = buffer.data(target + 1350);
    auto *t_1351 = buffer.data(target + 1351);
    auto *t_1352 = buffer.data(target + 1352);
    auto *t_1353 = buffer.data(target + 1353);
    auto *t_1354 = buffer.data(target + 1354);
    auto *t_1355 = buffer.data(target + 1355);
    auto *t_1356 = buffer.data(target + 1356);
    auto *t_1357 = buffer.data(target + 1357);
    auto *t_1358 = buffer.data(target + 1358);
    auto *t_1359 = buffer.data(target + 1359);
    auto *t_1360 = buffer.data(target + 1360);
    auto *t_1361 = buffer.data(target + 1361);
    auto *t_1362 = buffer.data(target + 1362);
    auto *t_1363 = buffer.data(target + 1363);
    auto *t_1364 = buffer.data(target + 1364);
    auto *t_1365 = buffer.data(target + 1365);
    auto *t_1366 = buffer.data(target + 1366);
    auto *t_1367 = buffer.data(target + 1367);
    auto *t_1368 = buffer.data(target + 1368);
    auto *t_1369 = buffer.data(target + 1369);
    auto *t_1370 = buffer.data(target + 1370);
    auto *t_1371 = buffer.data(target + 1371);
    auto *t_1372 = buffer.data(target + 1372);
    auto *t_1373 = buffer.data(target + 1373);
    auto *t_1374 = buffer.data(target + 1374);
    auto *t_1375 = buffer.data(target + 1375);
    auto *t_1376 = buffer.data(target + 1376);
    auto *t_1377 = buffer.data(target + 1377);
    auto *t_1378 = buffer.data(target + 1378);
    auto *t_1379 = buffer.data(target + 1379);
    auto *t_1380 = buffer.data(target + 1380);
    auto *t_1381 = buffer.data(target + 1381);
    auto *t_1382 = buffer.data(target + 1382);
    auto *t_1383 = buffer.data(target + 1383);
    auto *t_1384 = buffer.data(target + 1384);
    auto *t_1385 = buffer.data(target + 1385);
    auto *t_1386 = buffer.data(target + 1386);
    auto *t_1387 = buffer.data(target + 1387);
    auto *t_1388 = buffer.data(target + 1388);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_1134 = buffer.data(snh0 + 1134);
    const auto *snh0_1139 = buffer.data(snh0 + 1139);
    const auto *snh0_1143 = buffer.data(snh0 + 1143);
    const auto *snh0_1266 = buffer.data(snh0 + 1266);
    const auto *snh0_1269 = buffer.data(snh0 + 1269);
    const auto *snh0_1275 = buffer.data(snh0 + 1275);
    const auto *snh0_1277 = buffer.data(snh0 + 1277);
    const auto *snh0_1278 = buffer.data(snh0 + 1278);
    const auto *snh0_1280 = buffer.data(snh0 + 1280);
    const auto *snh0_1281 = buffer.data(snh0 + 1281);
    const auto *snh0_1284 = buffer.data(snh0 + 1284);
    const auto *snh0_1286 = buffer.data(snh0 + 1286);
    const auto *snh0_1287 = buffer.data(snh0 + 1287);
    const auto *snh0_1290 = buffer.data(snh0 + 1290);
    const auto *snh0_1296 = buffer.data(snh0 + 1296);
    const auto *snh0_1298 = buffer.data(snh0 + 1298);
    const auto *snh0_1299 = buffer.data(snh0 + 1299);
    const auto *snh0_1301 = buffer.data(snh0 + 1301);
    const auto *snh0_1302 = buffer.data(snh0 + 1302);
    const auto *snh0_1305 = buffer.data(snh0 + 1305);
    const auto *snh0_1307 = buffer.data(snh0 + 1307);
    const auto *snh0_1308 = buffer.data(snh0 + 1308);
    const auto *snh0_1311 = buffer.data(snh0 + 1311);
    const auto *snh0_1317 = buffer.data(snh0 + 1317);
    const auto *snh0_1319 = buffer.data(snh0 + 1319);
    const auto *snh0_1320 = buffer.data(snh0 + 1320);
    const auto *snh0_1322 = buffer.data(snh0 + 1322);
    const auto *snh0_1323 = buffer.data(snh0 + 1323);
    const auto *snh0_1326 = buffer.data(snh0 + 1326);
    const auto *snh0_1328 = buffer.data(snh0 + 1328);
    const auto *snh0_1329 = buffer.data(snh0 + 1329);
    const auto *snh0_1332 = buffer.data(snh0 + 1332);
    const auto *snh0_1338 = buffer.data(snh0 + 1338);
    const auto *snh0_1340 = buffer.data(snh0 + 1340);
    const auto *snh0_1341 = buffer.data(snh0 + 1341);
    const auto *snh0_1343 = buffer.data(snh0 + 1343);
    const auto *snh0_1347 = buffer.data(snh0 + 1347);
    const auto *snh0_1350 = buffer.data(snh0 + 1350);
    const auto *snh0_1359 = buffer.data(snh0 + 1359);
    const auto *snh0_1361 = buffer.data(snh0 + 1361);
    const auto *snh0_1362 = buffer.data(snh0 + 1362);
    const auto *snh0_1364 = buffer.data(snh0 + 1364);
    const auto *snh0_1365 = buffer.data(snh0 + 1365);
    const auto *snh0_1368 = buffer.data(snh0 + 1368);
    const auto *snh0_1370 = buffer.data(snh0 + 1370);
    const auto *snh0_1371 = buffer.data(snh0 + 1371);
    const auto *snh0_1374 = buffer.data(snh0 + 1374);
    const auto *snh0_1380 = buffer.data(snh0 + 1380);
    const auto *snh0_1382 = buffer.data(snh0 + 1382);
    const auto *snh0_1383 = buffer.data(snh0 + 1383);
    const auto *snh0_1385 = buffer.data(snh0 + 1385);

    const auto *sng_738 = buffer.data(sng + 738);
    const auto *sng_745 = buffer.data(sng + 745);
    const auto *sng_750 = buffer.data(sng + 750);
    const auto *sng_753 = buffer.data(sng + 753);
    const auto *sng_755 = buffer.data(sng + 755);
    const auto *sng_760 = buffer.data(sng + 760);
    const auto *sng_764 = buffer.data(sng + 764);
    const auto *sng_765 = buffer.data(sng + 765);
    const auto *sng_767 = buffer.data(sng + 767);
    const auto *sng_768 = buffer.data(sng + 768);
    const auto *sng_770 = buffer.data(sng + 770);
    const auto *sng_775 = buffer.data(sng + 775);
    const auto *sng_779 = buffer.data(sng + 779);
    const auto *sng_780 = buffer.data(sng + 780);
    const auto *sng_782 = buffer.data(sng + 782);
    const auto *sng_783 = buffer.data(sng + 783);
    const auto *sng_785 = buffer.data(sng + 785);
    const auto *sng_790 = buffer.data(sng + 790);
    const auto *sng_794 = buffer.data(sng + 794);
    const auto *sng_795 = buffer.data(sng + 795);
    const auto *sng_797 = buffer.data(sng + 797);
    const auto *sng_798 = buffer.data(sng + 798);
    const auto *sng_800 = buffer.data(sng + 800);
    const auto *sng_805 = buffer.data(sng + 805);
    const auto *sng_809 = buffer.data(sng + 809);
    const auto *sng_810 = buffer.data(sng + 810);
    const auto *sng_812 = buffer.data(sng + 812);
    const auto *sng_813 = buffer.data(sng + 813);
    const auto *sng_815 = buffer.data(sng + 815);
    const auto *sng_820 = buffer.data(sng + 820);
    const auto *sng_824 = buffer.data(sng + 824);
    const auto *sng_825 = buffer.data(sng + 825);
    const auto *sng_906 = buffer.data(sng + 906);
    const auto *sng_909 = buffer.data(sng + 909);
    const auto *sng_910 = buffer.data(sng + 910);
    const auto *sng_911 = buffer.data(sng + 911);
    const auto *sng_912 = buffer.data(sng + 912);
    const auto *sng_913 = buffer.data(sng + 913);
    const auto *sng_914 = buffer.data(sng + 914);
    const auto *sng_915 = buffer.data(sng + 915);
    const auto *sng_918 = buffer.data(sng + 918);
    const auto *sng_920 = buffer.data(sng + 920);
    const auto *sng_921 = buffer.data(sng + 921);
    const auto *sng_924 = buffer.data(sng + 924);
    const auto *sng_925 = buffer.data(sng + 925);
    const auto *sng_926 = buffer.data(sng + 926);
    const auto *sng_927 = buffer.data(sng + 927);
    const auto *sng_928 = buffer.data(sng + 928);
    const auto *sng_929 = buffer.data(sng + 929);
    const auto *sng_930 = buffer.data(sng + 930);
    const auto *sng_933 = buffer.data(sng + 933);
    const auto *sng_935 = buffer.data(sng + 935);
    const auto *sng_936 = buffer.data(sng + 936);
    const auto *sng_939 = buffer.data(sng + 939);
    const auto *sng_940 = buffer.data(sng + 940);
    const auto *sng_941 = buffer.data(sng + 941);
    const auto *sng_942 = buffer.data(sng + 942);
    const auto *sng_943 = buffer.data(sng + 943);
    const auto *sng_944 = buffer.data(sng + 944);
    const auto *sng_945 = buffer.data(sng + 945);
    const auto *sng_948 = buffer.data(sng + 948);
    const auto *sng_950 = buffer.data(sng + 950);
    const auto *sng_951 = buffer.data(sng + 951);
    const auto *sng_954 = buffer.data(sng + 954);
    const auto *sng_955 = buffer.data(sng + 955);
    const auto *sng_956 = buffer.data(sng + 956);
    const auto *sng_957 = buffer.data(sng + 957);
    const auto *sng_958 = buffer.data(sng + 958);
    const auto *sng_959 = buffer.data(sng + 959);
    const auto *sng_963 = buffer.data(sng + 963);
    const auto *sng_966 = buffer.data(sng + 966);
    const auto *sng_970 = buffer.data(sng + 970);
    const auto *sng_971 = buffer.data(sng + 971);
    const auto *sng_972 = buffer.data(sng + 972);
    const auto *sng_973 = buffer.data(sng + 973);
    const auto *sng_974 = buffer.data(sng + 974);
    const auto *sng_975 = buffer.data(sng + 975);
    const auto *sng_978 = buffer.data(sng + 978);
    const auto *sng_980 = buffer.data(sng + 980);
    const auto *sng_981 = buffer.data(sng + 981);
    const auto *sng_984 = buffer.data(sng + 984);
    const auto *sng_985 = buffer.data(sng + 985);
    const auto *sng_986 = buffer.data(sng + 986);
    const auto *sng_987 = buffer.data(sng + 987);
    const auto *sng_988 = buffer.data(sng + 988);
    const auto *sng_989 = buffer.data(sng + 989);

    const auto *snh1_1134 = buffer.data(snh1 + 1134);
    const auto *snh1_1139 = buffer.data(snh1 + 1139);
    const auto *snh1_1143 = buffer.data(snh1 + 1143);
    const auto *snh1_1266 = buffer.data(snh1 + 1266);
    const auto *snh1_1269 = buffer.data(snh1 + 1269);
    const auto *snh1_1275 = buffer.data(snh1 + 1275);
    const auto *snh1_1277 = buffer.data(snh1 + 1277);
    const auto *snh1_1278 = buffer.data(snh1 + 1278);
    const auto *snh1_1280 = buffer.data(snh1 + 1280);
    const auto *snh1_1281 = buffer.data(snh1 + 1281);
    const auto *snh1_1284 = buffer.data(snh1 + 1284);
    const auto *snh1_1286 = buffer.data(snh1 + 1286);
    const auto *snh1_1287 = buffer.data(snh1 + 1287);
    const auto *snh1_1290 = buffer.data(snh1 + 1290);
    const auto *snh1_1296 = buffer.data(snh1 + 1296);
    const auto *snh1_1298 = buffer.data(snh1 + 1298);
    const auto *snh1_1299 = buffer.data(snh1 + 1299);
    const auto *snh1_1301 = buffer.data(snh1 + 1301);
    const auto *snh1_1302 = buffer.data(snh1 + 1302);
    const auto *snh1_1305 = buffer.data(snh1 + 1305);
    const auto *snh1_1307 = buffer.data(snh1 + 1307);
    const auto *snh1_1308 = buffer.data(snh1 + 1308);
    const auto *snh1_1311 = buffer.data(snh1 + 1311);
    const auto *snh1_1317 = buffer.data(snh1 + 1317);
    const auto *snh1_1319 = buffer.data(snh1 + 1319);
    const auto *snh1_1320 = buffer.data(snh1 + 1320);
    const auto *snh1_1322 = buffer.data(snh1 + 1322);
    const auto *snh1_1323 = buffer.data(snh1 + 1323);
    const auto *snh1_1326 = buffer.data(snh1 + 1326);
    const auto *snh1_1328 = buffer.data(snh1 + 1328);
    const auto *snh1_1329 = buffer.data(snh1 + 1329);
    const auto *snh1_1332 = buffer.data(snh1 + 1332);
    const auto *snh1_1338 = buffer.data(snh1 + 1338);
    const auto *snh1_1340 = buffer.data(snh1 + 1340);
    const auto *snh1_1341 = buffer.data(snh1 + 1341);
    const auto *snh1_1343 = buffer.data(snh1 + 1343);
    const auto *snh1_1347 = buffer.data(snh1 + 1347);
    const auto *snh1_1350 = buffer.data(snh1 + 1350);
    const auto *snh1_1359 = buffer.data(snh1 + 1359);
    const auto *snh1_1361 = buffer.data(snh1 + 1361);
    const auto *snh1_1362 = buffer.data(snh1 + 1362);
    const auto *snh1_1364 = buffer.data(snh1 + 1364);
    const auto *snh1_1365 = buffer.data(snh1 + 1365);
    const auto *snh1_1368 = buffer.data(snh1 + 1368);
    const auto *snh1_1370 = buffer.data(snh1 + 1370);
    const auto *snh1_1371 = buffer.data(snh1 + 1371);
    const auto *snh1_1374 = buffer.data(snh1 + 1374);
    const auto *snh1_1380 = buffer.data(snh1 + 1380);
    const auto *snh1_1382 = buffer.data(snh1 + 1382);
    const auto *snh1_1383 = buffer.data(snh1 + 1383);
    const auto *snh1_1385 = buffer.data(snh1 + 1385);

    const auto *sof0_660 = buffer.data(sof0 + 660);

    const auto *sof1_660 = buffer.data(sof1 + 660);

    const auto *sog_903 = buffer.data(sog + 903);
    const auto *sog_905 = buffer.data(sog + 905);
    const auto *sog_910 = buffer.data(sog + 910);
    const auto *sog_911 = buffer.data(sog + 911);
    const auto *sog_912 = buffer.data(sog + 912);
    const auto *sog_913 = buffer.data(sog + 913);
    const auto *sog_914 = buffer.data(sog + 914);
    const auto *sog_915 = buffer.data(sog + 915);
    const auto *sog_917 = buffer.data(sog + 917);
    const auto *sog_918 = buffer.data(sog + 918);
    const auto *sog_920 = buffer.data(sog + 920);
    const auto *sog_925 = buffer.data(sog + 925);
    const auto *sog_926 = buffer.data(sog + 926);
    const auto *sog_927 = buffer.data(sog + 927);
    const auto *sog_928 = buffer.data(sog + 928);
    const auto *sog_929 = buffer.data(sog + 929);
    const auto *sog_930 = buffer.data(sog + 930);
    const auto *sog_932 = buffer.data(sog + 932);
    const auto *sog_933 = buffer.data(sog + 933);
    const auto *sog_935 = buffer.data(sog + 935);
    const auto *sog_940 = buffer.data(sog + 940);
    const auto *sog_941 = buffer.data(sog + 941);
    const auto *sog_942 = buffer.data(sog + 942);
    const auto *sog_943 = buffer.data(sog + 943);
    const auto *sog_944 = buffer.data(sog + 944);
    const auto *sog_945 = buffer.data(sog + 945);
    const auto *sog_947 = buffer.data(sog + 947);
    const auto *sog_948 = buffer.data(sog + 948);
    const auto *sog_950 = buffer.data(sog + 950);
    const auto *sog_955 = buffer.data(sog + 955);
    const auto *sog_956 = buffer.data(sog + 956);
    const auto *sog_957 = buffer.data(sog + 957);
    const auto *sog_958 = buffer.data(sog + 958);
    const auto *sog_959 = buffer.data(sog + 959);
    const auto *sog_960 = buffer.data(sog + 960);
    const auto *sog_962 = buffer.data(sog + 962);
    const auto *sog_963 = buffer.data(sog + 963);
    const auto *sog_965 = buffer.data(sog + 965);
    const auto *sog_970 = buffer.data(sog + 970);
    const auto *sog_971 = buffer.data(sog + 971);
    const auto *sog_972 = buffer.data(sog + 972);
    const auto *sog_973 = buffer.data(sog + 973);
    const auto *sog_974 = buffer.data(sog + 974);
    const auto *sog_975 = buffer.data(sog + 975);
    const auto *sog_977 = buffer.data(sog + 977);
    const auto *sog_978 = buffer.data(sog + 978);
    const auto *sog_980 = buffer.data(sog + 980);
    const auto *sog_985 = buffer.data(sog + 985);
    const auto *sog_986 = buffer.data(sog + 986);
    const auto *sog_987 = buffer.data(sog + 987);
    const auto *sog_988 = buffer.data(sog + 988);
    const auto *sog_989 = buffer.data(sog + 989);
    const auto *sog_990 = buffer.data(sog + 990);

#pragma omp simd aligned(t_1266, t_1267, t_1268, pb_x, pc_x, pc_y, pc_z, snh0_1266, sng_738, \
                         sng_755, sng_906, snh1_1266, sog_903, \
                         sog_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1266[k] = pb_x[k] * snh0_1266[k]
                    + f_10 * sng_906[k]
                    - f_8 * pc_x[k] * snh1_1266[k];

        t_1267[k] = f_18 * sng_738[k]
                    + f_3 * pc_z[k] * sog_903[k];

        t_1268[k] = f_18 * sng_755[k]
                    + f_3 * pc_y[k] * sog_905[k];
    }

#pragma omp simd aligned(t_1269, t_1270, t_1271, t_1272, pb_x, pc_x, snh0_1269, sng_909, \
                         sng_910, sng_911, sng_912, snh1_1269, sog_910, sog_911, \
                         sog_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1269[k] = pb_x[k] * snh0_1269[k]
                    + f_10 * sng_909[k]
                    - f_8 * pc_x[k] * snh1_1269[k];

        t_1270[k] = f_9 * sng_910[k]
                    + f_3 * pc_x[k] * sog_910[k];

        t_1271[k] = f_9 * sng_911[k]
                    + f_3 * pc_x[k] * sog_911[k];

        t_1272[k] = f_9 * sng_912[k]
                    + f_3 * pc_x[k] * sog_912[k];
    }

#pragma omp simd aligned(t_1273, t_1274, t_1275, t_1276, pb_x, pc_x, pc_z, snh0_1275, sng_745, \
                         sng_913, sng_914, snh1_1275, sog_910, sog_913, \
                         sog_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1273[k] = f_9 * sng_913[k]
                    + f_3 * pc_x[k] * sog_913[k];

        t_1274[k] = f_9 * sng_914[k]
                    + f_3 * pc_x[k] * sog_914[k];

        t_1275[k] = pb_x[k] * snh0_1275[k]
                    - f_8 * pc_x[k] * snh1_1275[k];

        t_1276[k] = f_18 * sng_745[k]
                    + f_3 * pc_z[k] * sog_910[k];
    }

#pragma omp simd aligned(t_1277, t_1278, t_1279, t_1280, pb_x, pc_x, pc_y, snh0_1277, \
                         snh0_1278, snh0_1280, sng_764, snh1_1277, snh1_1278, snh1_1280, \
                         sog_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1277[k] = pb_x[k] * snh0_1277[k]
                    - f_8 * pc_x[k] * snh1_1277[k];

        t_1278[k] = pb_x[k] * snh0_1278[k]
                    - f_8 * pc_x[k] * snh1_1278[k];

        t_1279[k] = f_18 * sng_764[k]
                    + f_3 * pc_y[k] * sog_914[k];

        t_1280[k] = pb_x[k] * snh0_1280[k]
                    - f_8 * pc_x[k] * snh1_1280[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, pb_x, pc_x, pc_y, pc_z, snh0_1281, sng_750, \
                         sng_765, sng_915, snh1_1281, sog_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = pb_x[k] * snh0_1281[k]
                    + f_18 * sng_915[k]
                    - f_8 * pc_x[k] * snh1_1281[k];

        t_1282[k] = f_16 * sng_765[k]
                    + f_3 * pc_y[k] * sog_915[k];

        t_1283[k] = f_17 * sng_750[k]
                    + f_3 * pc_z[k] * sog_915[k];
    }

#pragma omp simd aligned(t_1284, t_1285, t_1286, pb_x, pc_x, pc_y, snh0_1284, snh0_1286, \
                         sng_767, sng_918, sng_920, snh1_1284, snh1_1286, \
                         sog_917 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1284[k] = pb_x[k] * snh0_1284[k]
                    + f_11 * sng_918[k]
                    - f_8 * pc_x[k] * snh1_1284[k];

        t_1285[k] = f_16 * sng_767[k]
                    + f_3 * pc_y[k] * sog_917[k];

        t_1286[k] = pb_x[k] * snh0_1286[k]
                    + f_11 * sng_920[k]
                    - f_8 * pc_x[k] * snh1_1286[k];
    }

#pragma omp simd aligned(t_1287, t_1288, t_1289, pb_x, pc_x, pc_y, pc_z, snh0_1287, sng_753, \
                         sng_770, sng_921, snh1_1287, sog_918, \
                         sog_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1287[k] = pb_x[k] * snh0_1287[k]
                    + f_10 * sng_921[k]
                    - f_8 * pc_x[k] * snh1_1287[k];

        t_1288[k] = f_17 * sng_753[k]
                    + f_3 * pc_z[k] * sog_918[k];

        t_1289[k] = f_16 * sng_770[k]
                    + f_3 * pc_y[k] * sog_920[k];
    }

#pragma omp simd aligned(t_1290, t_1291, t_1292, t_1293, pb_x, pc_x, snh0_1290, sng_924, \
                         sng_925, sng_926, sng_927, snh1_1290, sog_925, sog_926, \
                         sog_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1290[k] = pb_x[k] * snh0_1290[k]
                    + f_10 * sng_924[k]
                    - f_8 * pc_x[k] * snh1_1290[k];

        t_1291[k] = f_9 * sng_925[k]
                    + f_3 * pc_x[k] * sog_925[k];

        t_1292[k] = f_9 * sng_926[k]
                    + f_3 * pc_x[k] * sog_926[k];

        t_1293[k] = f_9 * sng_927[k]
                    + f_3 * pc_x[k] * sog_927[k];
    }

#pragma omp simd aligned(t_1294, t_1295, t_1296, t_1297, pb_x, pc_x, pc_z, snh0_1296, sng_760, \
                         sng_928, sng_929, snh1_1296, sog_925, sog_928, \
                         sog_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1294[k] = f_9 * sng_928[k]
                    + f_3 * pc_x[k] * sog_928[k];

        t_1295[k] = f_9 * sng_929[k]
                    + f_3 * pc_x[k] * sog_929[k];

        t_1296[k] = pb_x[k] * snh0_1296[k]
                    - f_8 * pc_x[k] * snh1_1296[k];

        t_1297[k] = f_17 * sng_760[k]
                    + f_3 * pc_z[k] * sog_925[k];
    }

#pragma omp simd aligned(t_1298, t_1299, t_1300, t_1301, pb_x, pc_x, pc_y, snh0_1298, \
                         snh0_1299, snh0_1301, sng_779, snh1_1298, snh1_1299, snh1_1301, \
                         sog_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1298[k] = pb_x[k] * snh0_1298[k]
                    - f_8 * pc_x[k] * snh1_1298[k];

        t_1299[k] = pb_x[k] * snh0_1299[k]
                    - f_8 * pc_x[k] * snh1_1299[k];

        t_1300[k] = f_16 * sng_779[k]
                    + f_3 * pc_y[k] * sog_929[k];

        t_1301[k] = pb_x[k] * snh0_1301[k]
                    - f_8 * pc_x[k] * snh1_1301[k];
    }

#pragma omp simd aligned(t_1302, t_1303, t_1304, pb_x, pc_x, pc_y, pc_z, snh0_1302, sng_765, \
                         sng_780, sng_930, snh1_1302, sog_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1302[k] = pb_x[k] * snh0_1302[k]
                    + f_18 * sng_930[k]
                    - f_8 * pc_x[k] * snh1_1302[k];

        t_1303[k] = f_11 * sng_780[k]
                    + f_3 * pc_y[k] * sog_930[k];

        t_1304[k] = f_15 * sng_765[k]
                    + f_3 * pc_z[k] * sog_930[k];
    }

#pragma omp simd aligned(t_1305, t_1306, t_1307, pb_x, pc_x, pc_y, snh0_1305, snh0_1307, \
                         sng_782, sng_933, sng_935, snh1_1305, snh1_1307, \
                         sog_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1305[k] = pb_x[k] * snh0_1305[k]
                    + f_11 * sng_933[k]
                    - f_8 * pc_x[k] * snh1_1305[k];

        t_1306[k] = f_11 * sng_782[k]
                    + f_3 * pc_y[k] * sog_932[k];

        t_1307[k] = pb_x[k] * snh0_1307[k]
                    + f_11 * sng_935[k]
                    - f_8 * pc_x[k] * snh1_1307[k];
    }

#pragma omp simd aligned(t_1308, t_1309, t_1310, pb_x, pc_x, pc_y, pc_z, snh0_1308, sng_768, \
                         sng_785, sng_936, snh1_1308, sog_933, \
                         sog_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1308[k] = pb_x[k] * snh0_1308[k]
                    + f_10 * sng_936[k]
                    - f_8 * pc_x[k] * snh1_1308[k];

        t_1309[k] = f_15 * sng_768[k]
                    + f_3 * pc_z[k] * sog_933[k];

        t_1310[k] = f_11 * sng_785[k]
                    + f_3 * pc_y[k] * sog_935[k];
    }

#pragma omp simd aligned(t_1311, t_1312, t_1313, t_1314, pb_x, pc_x, snh0_1311, sng_939, \
                         sng_940, sng_941, sng_942, snh1_1311, sog_940, sog_941, \
                         sog_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1311[k] = pb_x[k] * snh0_1311[k]
                    + f_10 * sng_939[k]
                    - f_8 * pc_x[k] * snh1_1311[k];

        t_1312[k] = f_9 * sng_940[k]
                    + f_3 * pc_x[k] * sog_940[k];

        t_1313[k] = f_9 * sng_941[k]
                    + f_3 * pc_x[k] * sog_941[k];

        t_1314[k] = f_9 * sng_942[k]
                    + f_3 * pc_x[k] * sog_942[k];
    }

#pragma omp simd aligned(t_1315, t_1316, t_1317, t_1318, pb_x, pc_x, pc_z, snh0_1317, sng_775, \
                         sng_943, sng_944, snh1_1317, sog_940, sog_943, \
                         sog_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1315[k] = f_9 * sng_943[k]
                    + f_3 * pc_x[k] * sog_943[k];

        t_1316[k] = f_9 * sng_944[k]
                    + f_3 * pc_x[k] * sog_944[k];

        t_1317[k] = pb_x[k] * snh0_1317[k]
                    - f_8 * pc_x[k] * snh1_1317[k];

        t_1318[k] = f_15 * sng_775[k]
                    + f_3 * pc_z[k] * sog_940[k];
    }

#pragma omp simd aligned(t_1319, t_1320, t_1321, t_1322, pb_x, pc_x, pc_y, snh0_1319, \
                         snh0_1320, snh0_1322, sng_794, snh1_1319, snh1_1320, snh1_1322, \
                         sog_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1319[k] = pb_x[k] * snh0_1319[k]
                    - f_8 * pc_x[k] * snh1_1319[k];

        t_1320[k] = pb_x[k] * snh0_1320[k]
                    - f_8 * pc_x[k] * snh1_1320[k];

        t_1321[k] = f_11 * sng_794[k]
                    + f_3 * pc_y[k] * sog_944[k];

        t_1322[k] = pb_x[k] * snh0_1322[k]
                    - f_8 * pc_x[k] * snh1_1322[k];
    }

#pragma omp simd aligned(t_1323, t_1324, t_1325, pb_x, pc_x, pc_y, pc_z, snh0_1323, sng_780, \
                         sng_795, sng_945, snh1_1323, sog_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1323[k] = pb_x[k] * snh0_1323[k]
                    + f_18 * sng_945[k]
                    - f_8 * pc_x[k] * snh1_1323[k];

        t_1324[k] = f_10 * sng_795[k]
                    + f_3 * pc_y[k] * sog_945[k];

        t_1325[k] = f_14 * sng_780[k]
                    + f_3 * pc_z[k] * sog_945[k];
    }

#pragma omp simd aligned(t_1326, t_1327, t_1328, pb_x, pc_x, pc_y, snh0_1326, snh0_1328, \
                         sng_797, sng_948, sng_950, snh1_1326, snh1_1328, \
                         sog_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1326[k] = pb_x[k] * snh0_1326[k]
                    + f_11 * sng_948[k]
                    - f_8 * pc_x[k] * snh1_1326[k];

        t_1327[k] = f_10 * sng_797[k]
                    + f_3 * pc_y[k] * sog_947[k];

        t_1328[k] = pb_x[k] * snh0_1328[k]
                    + f_11 * sng_950[k]
                    - f_8 * pc_x[k] * snh1_1328[k];
    }

#pragma omp simd aligned(t_1329, t_1330, t_1331, pb_x, pc_x, pc_y, pc_z, snh0_1329, sng_783, \
                         sng_800, sng_951, snh1_1329, sog_948, \
                         sog_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1329[k] = pb_x[k] * snh0_1329[k]
                    + f_10 * sng_951[k]
                    - f_8 * pc_x[k] * snh1_1329[k];

        t_1330[k] = f_14 * sng_783[k]
                    + f_3 * pc_z[k] * sog_948[k];

        t_1331[k] = f_10 * sng_800[k]
                    + f_3 * pc_y[k] * sog_950[k];
    }

#pragma omp simd aligned(t_1332, t_1333, t_1334, t_1335, pb_x, pc_x, snh0_1332, sng_954, \
                         sng_955, sng_956, sng_957, snh1_1332, sog_955, sog_956, \
                         sog_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1332[k] = pb_x[k] * snh0_1332[k]
                    + f_10 * sng_954[k]
                    - f_8 * pc_x[k] * snh1_1332[k];

        t_1333[k] = f_9 * sng_955[k]
                    + f_3 * pc_x[k] * sog_955[k];

        t_1334[k] = f_9 * sng_956[k]
                    + f_3 * pc_x[k] * sog_956[k];

        t_1335[k] = f_9 * sng_957[k]
                    + f_3 * pc_x[k] * sog_957[k];
    }

#pragma omp simd aligned(t_1336, t_1337, t_1338, t_1339, pb_x, pc_x, pc_z, snh0_1338, sng_790, \
                         sng_958, sng_959, snh1_1338, sog_955, sog_958, \
                         sog_959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1336[k] = f_9 * sng_958[k]
                    + f_3 * pc_x[k] * sog_958[k];

        t_1337[k] = f_9 * sng_959[k]
                    + f_3 * pc_x[k] * sog_959[k];

        t_1338[k] = pb_x[k] * snh0_1338[k]
                    - f_8 * pc_x[k] * snh1_1338[k];

        t_1339[k] = f_14 * sng_790[k]
                    + f_3 * pc_z[k] * sog_955[k];
    }

#pragma omp simd aligned(t_1340, t_1341, t_1342, t_1343, pb_x, pc_x, pc_y, snh0_1340, \
                         snh0_1341, snh0_1343, sng_809, snh1_1340, snh1_1341, snh1_1343, \
                         sog_959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1340[k] = pb_x[k] * snh0_1340[k]
                    - f_8 * pc_x[k] * snh1_1340[k];

        t_1341[k] = pb_x[k] * snh0_1341[k]
                    - f_8 * pc_x[k] * snh1_1341[k];

        t_1342[k] = f_10 * sng_809[k]
                    + f_3 * pc_y[k] * sog_959[k];

        t_1343[k] = pb_x[k] * snh0_1343[k]
                    - f_8 * pc_x[k] * snh1_1343[k];
    }

#pragma omp simd aligned(t_1344, t_1345, t_1346, pb_y, pc_y, pc_z, snh0_1134, sng_795, \
                         sng_810, snh1_1134, sog_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1344[k] = pb_y[k] * snh0_1134[k]
                    - f_8 * pc_y[k] * snh1_1134[k];

        t_1345[k] = f_9 * sng_810[k]
                    + f_3 * pc_y[k] * sog_960[k];

        t_1346[k] = f_13 * sng_795[k]
                    + f_3 * pc_z[k] * sog_960[k];
    }

#pragma omp simd aligned(t_1347, t_1348, t_1349, pb_x, pb_y, pc_x, pc_y, snh0_1139, snh0_1347, \
                         sng_812, sng_963, snh1_1139, snh1_1347, \
                         sog_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1347[k] = pb_x[k] * snh0_1347[k]
                    + f_11 * sng_963[k]
                    - f_8 * pc_x[k] * snh1_1347[k];

        t_1348[k] = f_9 * sng_812[k]
                    + f_3 * pc_y[k] * sog_962[k];

        t_1349[k] = pb_y[k] * snh0_1139[k]
                    - f_8 * pc_y[k] * snh1_1139[k];
    }

#pragma omp simd aligned(t_1350, t_1351, t_1352, pb_x, pc_x, pc_y, pc_z, snh0_1350, sng_798, \
                         sng_815, sng_966, snh1_1350, sog_963, \
                         sog_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1350[k] = pb_x[k] * snh0_1350[k]
                    + f_10 * sng_966[k]
                    - f_8 * pc_x[k] * snh1_1350[k];

        t_1351[k] = f_13 * sng_798[k]
                    + f_3 * pc_z[k] * sog_963[k];

        t_1352[k] = f_9 * sng_815[k]
                    + f_3 * pc_y[k] * sog_965[k];
    }

#pragma omp simd aligned(t_1353, t_1354, t_1355, t_1356, pb_y, pc_x, pc_y, snh0_1143, sng_970, \
                         sng_971, sng_972, snh1_1143, sog_970, sog_971, \
                         sog_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1353[k] = pb_y[k] * snh0_1143[k]
                    - f_8 * pc_y[k] * snh1_1143[k];

        t_1354[k] = f_9 * sng_970[k]
                    + f_3 * pc_x[k] * sog_970[k];

        t_1355[k] = f_9 * sng_971[k]
                    + f_3 * pc_x[k] * sog_971[k];

        t_1356[k] = f_9 * sng_972[k]
                    + f_3 * pc_x[k] * sog_972[k];
    }

#pragma omp simd aligned(t_1357, t_1358, t_1359, t_1360, pb_x, pc_x, pc_z, snh0_1359, sng_805, \
                         sng_973, sng_974, snh1_1359, sog_970, sog_973, \
                         sog_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1357[k] = f_9 * sng_973[k]
                    + f_3 * pc_x[k] * sog_973[k];

        t_1358[k] = f_9 * sng_974[k]
                    + f_3 * pc_x[k] * sog_974[k];

        t_1359[k] = pb_x[k] * snh0_1359[k]
                    - f_8 * pc_x[k] * snh1_1359[k];

        t_1360[k] = f_13 * sng_805[k]
                    + f_3 * pc_z[k] * sog_970[k];
    }

#pragma omp simd aligned(t_1361, t_1362, t_1363, t_1364, pb_x, pc_x, pc_y, snh0_1361, \
                         snh0_1362, snh0_1364, sng_824, snh1_1361, snh1_1362, snh1_1364, \
                         sog_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1361[k] = pb_x[k] * snh0_1361[k]
                    - f_8 * pc_x[k] * snh1_1361[k];

        t_1362[k] = pb_x[k] * snh0_1362[k]
                    - f_8 * pc_x[k] * snh1_1362[k];

        t_1363[k] = f_9 * sng_824[k]
                    + f_3 * pc_y[k] * sog_974[k];

        t_1364[k] = pb_x[k] * snh0_1364[k]
                    - f_8 * pc_x[k] * snh1_1364[k];
    }

#pragma omp simd aligned(t_1365, t_1366, t_1367, t_1368, pb_x, pc_x, pc_y, pc_z, snh0_1365, \
                         snh0_1368, sng_810, sng_975, sng_978, snh1_1365, snh1_1368, \
                         sog_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1365[k] = pb_x[k] * snh0_1365[k]
                    + f_18 * sng_975[k]
                    - f_8 * pc_x[k] * snh1_1365[k];

        t_1366[k] = f_3 * pc_y[k] * sog_975[k];

        t_1367[k] = f_12 * sng_810[k]
                    + f_3 * pc_z[k] * sog_975[k];

        t_1368[k] = pb_x[k] * snh0_1368[k]
                    + f_11 * sng_978[k]
                    - f_8 * pc_x[k] * snh1_1368[k];
    }

#pragma omp simd aligned(t_1369, t_1370, t_1371, pb_x, pc_x, pc_y, snh0_1370, snh0_1371, \
                         sng_980, sng_981, snh1_1370, snh1_1371, \
                         sog_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1369[k] = f_3 * pc_y[k] * sog_977[k];

        t_1370[k] = pb_x[k] * snh0_1370[k]
                    + f_11 * sng_980[k]
                    - f_8 * pc_x[k] * snh1_1370[k];

        t_1371[k] = pb_x[k] * snh0_1371[k]
                    + f_10 * sng_981[k]
                    - f_8 * pc_x[k] * snh1_1371[k];
    }

#pragma omp simd aligned(t_1372, t_1373, t_1374, t_1375, pb_x, pc_x, pc_y, pc_z, snh0_1374, \
                         sng_813, sng_984, sng_985, snh1_1374, sog_978, sog_980, \
                         sog_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1372[k] = f_12 * sng_813[k]
                    + f_3 * pc_z[k] * sog_978[k];

        t_1373[k] = f_3 * pc_y[k] * sog_980[k];

        t_1374[k] = pb_x[k] * snh0_1374[k]
                    + f_10 * sng_984[k]
                    - f_8 * pc_x[k] * snh1_1374[k];

        t_1375[k] = f_9 * sng_985[k]
                    + f_3 * pc_x[k] * sog_985[k];
    }

#pragma omp simd aligned(t_1376, t_1377, t_1378, t_1379, pc_x, sng_986, sng_987, sng_988, \
                         sng_989, sog_986, sog_987, sog_988, sog_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1376[k] = f_9 * sng_986[k]
                    + f_3 * pc_x[k] * sog_986[k];

        t_1377[k] = f_9 * sng_987[k]
                    + f_3 * pc_x[k] * sog_987[k];

        t_1378[k] = f_9 * sng_988[k]
                    + f_3 * pc_x[k] * sog_988[k];

        t_1379[k] = f_9 * sng_989[k]
                    + f_3 * pc_x[k] * sog_989[k];
    }

#pragma omp simd aligned(t_1380, t_1381, t_1382, t_1383, pb_x, pc_x, pc_z, snh0_1380, \
                         snh0_1382, snh0_1383, sng_820, snh1_1380, snh1_1382, snh1_1383, \
                         sog_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1380[k] = pb_x[k] * snh0_1380[k]
                    - f_8 * pc_x[k] * snh1_1380[k];

        t_1381[k] = f_12 * sng_820[k]
                    + f_3 * pc_z[k] * sog_985[k];

        t_1382[k] = pb_x[k] * snh0_1382[k]
                    - f_8 * pc_x[k] * snh1_1382[k];

        t_1383[k] = pb_x[k] * snh0_1383[k]
                    - f_8 * pc_x[k] * snh1_1383[k];
    }

#pragma omp simd aligned(t_1384, t_1385, t_1386, t_1387, t_1388, pb_x, pc_x, pc_y, pc_z, \
                         snh0_1385, sng_825, snh1_1385, sof0_660, sof1_660, sog_989, \
                         sog_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1384[k] = f_3 * pc_y[k] * sog_989[k];

        t_1385[k] = pb_x[k] * snh0_1385[k]
                    - f_8 * pc_x[k] * snh1_1385[k];

        t_1386[k] = f_1 * sof0_660[k]
                    - f_2 * sof1_660[k]
                    + f_3 * pc_x[k] * sog_990[k];

        t_1387[k] = f_0 * sng_825[k]
                    + f_3 * pc_y[k] * sog_990[k];

        t_1388[k] = f_3 * pc_z[k] * sog_990[k];
    }
}

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t snh0,
                                                           const size_t sng, const size_t snh1,
                                                           const size_t sof0, const size_t sof1,
                                                           const size_t sog, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 5.0 / q;
    const auto f_13 = 4.5 / q;
    const auto f_14 = 4.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_1389 = buffer.data(target + 1389);
    auto *t_1390 = buffer.data(target + 1390);
    auto *t_1391 = buffer.data(target + 1391);
    auto *t_1392 = buffer.data(target + 1392);
    auto *t_1393 = buffer.data(target + 1393);
    auto *t_1394 = buffer.data(target + 1394);
    auto *t_1395 = buffer.data(target + 1395);
    auto *t_1396 = buffer.data(target + 1396);
    auto *t_1397 = buffer.data(target + 1397);
    auto *t_1398 = buffer.data(target + 1398);
    auto *t_1399 = buffer.data(target + 1399);
    auto *t_1400 = buffer.data(target + 1400);
    auto *t_1401 = buffer.data(target + 1401);
    auto *t_1402 = buffer.data(target + 1402);
    auto *t_1403 = buffer.data(target + 1403);
    auto *t_1404 = buffer.data(target + 1404);
    auto *t_1405 = buffer.data(target + 1405);
    auto *t_1406 = buffer.data(target + 1406);
    auto *t_1407 = buffer.data(target + 1407);
    auto *t_1408 = buffer.data(target + 1408);
    auto *t_1409 = buffer.data(target + 1409);
    auto *t_1410 = buffer.data(target + 1410);
    auto *t_1411 = buffer.data(target + 1411);
    auto *t_1412 = buffer.data(target + 1412);
    auto *t_1413 = buffer.data(target + 1413);
    auto *t_1414 = buffer.data(target + 1414);
    auto *t_1415 = buffer.data(target + 1415);
    auto *t_1416 = buffer.data(target + 1416);
    auto *t_1417 = buffer.data(target + 1417);
    auto *t_1418 = buffer.data(target + 1418);
    auto *t_1419 = buffer.data(target + 1419);
    auto *t_1420 = buffer.data(target + 1420);
    auto *t_1421 = buffer.data(target + 1421);
    auto *t_1422 = buffer.data(target + 1422);
    auto *t_1423 = buffer.data(target + 1423);
    auto *t_1424 = buffer.data(target + 1424);
    auto *t_1425 = buffer.data(target + 1425);
    auto *t_1426 = buffer.data(target + 1426);
    auto *t_1427 = buffer.data(target + 1427);
    auto *t_1428 = buffer.data(target + 1428);
    auto *t_1429 = buffer.data(target + 1429);
    auto *t_1430 = buffer.data(target + 1430);
    auto *t_1431 = buffer.data(target + 1431);
    auto *t_1432 = buffer.data(target + 1432);
    auto *t_1433 = buffer.data(target + 1433);
    auto *t_1434 = buffer.data(target + 1434);
    auto *t_1435 = buffer.data(target + 1435);
    auto *t_1436 = buffer.data(target + 1436);
    auto *t_1437 = buffer.data(target + 1437);
    auto *t_1438 = buffer.data(target + 1438);
    auto *t_1439 = buffer.data(target + 1439);
    auto *t_1440 = buffer.data(target + 1440);
    auto *t_1441 = buffer.data(target + 1441);
    auto *t_1442 = buffer.data(target + 1442);
    auto *t_1443 = buffer.data(target + 1443);
    auto *t_1444 = buffer.data(target + 1444);
    auto *t_1445 = buffer.data(target + 1445);
    auto *t_1446 = buffer.data(target + 1446);
    auto *t_1447 = buffer.data(target + 1447);
    auto *t_1448 = buffer.data(target + 1448);
    auto *t_1449 = buffer.data(target + 1449);
    auto *t_1450 = buffer.data(target + 1450);
    auto *t_1451 = buffer.data(target + 1451);
    auto *t_1452 = buffer.data(target + 1452);
    auto *t_1453 = buffer.data(target + 1453);
    auto *t_1454 = buffer.data(target + 1454);
    auto *t_1455 = buffer.data(target + 1455);
    auto *t_1456 = buffer.data(target + 1456);
    auto *t_1457 = buffer.data(target + 1457);
    auto *t_1458 = buffer.data(target + 1458);
    auto *t_1459 = buffer.data(target + 1459);
    auto *t_1460 = buffer.data(target + 1460);
    auto *t_1461 = buffer.data(target + 1461);
    auto *t_1462 = buffer.data(target + 1462);
    auto *t_1463 = buffer.data(target + 1463);
    auto *t_1464 = buffer.data(target + 1464);
    auto *t_1465 = buffer.data(target + 1465);
    auto *t_1466 = buffer.data(target + 1466);
    auto *t_1467 = buffer.data(target + 1467);
    auto *t_1468 = buffer.data(target + 1468);
    auto *t_1469 = buffer.data(target + 1469);
    auto *t_1470 = buffer.data(target + 1470);
    auto *t_1471 = buffer.data(target + 1471);
    auto *t_1472 = buffer.data(target + 1472);
    auto *t_1473 = buffer.data(target + 1473);
    auto *t_1474 = buffer.data(target + 1474);
    auto *t_1475 = buffer.data(target + 1475);
    auto *t_1476 = buffer.data(target + 1476);
    auto *t_1477 = buffer.data(target + 1477);
    auto *t_1478 = buffer.data(target + 1478);
    auto *t_1479 = buffer.data(target + 1479);
    auto *t_1480 = buffer.data(target + 1480);
    auto *t_1481 = buffer.data(target + 1481);
    auto *t_1482 = buffer.data(target + 1482);
    auto *t_1483 = buffer.data(target + 1483);
    auto *t_1484 = buffer.data(target + 1484);
    auto *t_1485 = buffer.data(target + 1485);
    auto *t_1486 = buffer.data(target + 1486);
    auto *t_1487 = buffer.data(target + 1487);
    auto *t_1488 = buffer.data(target + 1488);
    auto *t_1489 = buffer.data(target + 1489);
    auto *t_1490 = buffer.data(target + 1490);
    auto *t_1491 = buffer.data(target + 1491);
    auto *t_1492 = buffer.data(target + 1492);
    auto *t_1493 = buffer.data(target + 1493);
    auto *t_1494 = buffer.data(target + 1494);
    auto *t_1495 = buffer.data(target + 1495);
    auto *t_1496 = buffer.data(target + 1496);
    auto *t_1497 = buffer.data(target + 1497);
    auto *t_1498 = buffer.data(target + 1498);
    auto *t_1499 = buffer.data(target + 1499);
    auto *t_1500 = buffer.data(target + 1500);
    auto *t_1501 = buffer.data(target + 1501);
    auto *t_1502 = buffer.data(target + 1502);
    auto *t_1503 = buffer.data(target + 1503);
    auto *t_1504 = buffer.data(target + 1504);
    auto *t_1505 = buffer.data(target + 1505);
    auto *t_1506 = buffer.data(target + 1506);
    auto *t_1507 = buffer.data(target + 1507);
    auto *t_1508 = buffer.data(target + 1508);
    auto *t_1509 = buffer.data(target + 1509);
    auto *t_1510 = buffer.data(target + 1510);
    auto *t_1511 = buffer.data(target + 1511);
    auto *t_1512 = buffer.data(target + 1512);
    auto *t_1513 = buffer.data(target + 1513);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_1155 = buffer.data(snh0 + 1155);
    const auto *snh0_1158 = buffer.data(snh0 + 1158);
    const auto *snh0_1161 = buffer.data(snh0 + 1161);
    const auto *snh0_1170 = buffer.data(snh0 + 1170);
    const auto *snh0_1172 = buffer.data(snh0 + 1172);
    const auto *snh0_1173 = buffer.data(snh0 + 1173);

    const auto *sng_825 = buffer.data(sng + 825);
    const auto *sng_827 = buffer.data(sng + 827);
    const auto *sng_828 = buffer.data(sng + 828);
    const auto *sng_830 = buffer.data(sng + 830);
    const auto *sng_835 = buffer.data(sng + 835);
    const auto *sng_836 = buffer.data(sng + 836);
    const auto *sng_837 = buffer.data(sng + 837);
    const auto *sng_838 = buffer.data(sng + 838);
    const auto *sng_839 = buffer.data(sng + 839);
    const auto *sng_840 = buffer.data(sng + 840);
    const auto *sng_842 = buffer.data(sng + 842);
    const auto *sng_843 = buffer.data(sng + 843);
    const auto *sng_845 = buffer.data(sng + 845);
    const auto *sng_850 = buffer.data(sng + 850);
    const auto *sng_854 = buffer.data(sng + 854);
    const auto *sng_855 = buffer.data(sng + 855);
    const auto *sng_857 = buffer.data(sng + 857);
    const auto *sng_858 = buffer.data(sng + 858);
    const auto *sng_860 = buffer.data(sng + 860);
    const auto *sng_865 = buffer.data(sng + 865);
    const auto *sng_867 = buffer.data(sng + 867);
    const auto *sng_868 = buffer.data(sng + 868);
    const auto *sng_869 = buffer.data(sng + 869);
    const auto *sng_870 = buffer.data(sng + 870);
    const auto *sng_872 = buffer.data(sng + 872);
    const auto *sng_873 = buffer.data(sng + 873);
    const auto *sng_875 = buffer.data(sng + 875);
    const auto *sng_880 = buffer.data(sng + 880);
    const auto *sng_882 = buffer.data(sng + 882);
    const auto *sng_883 = buffer.data(sng + 883);
    const auto *sng_884 = buffer.data(sng + 884);
    const auto *sng_885 = buffer.data(sng + 885);
    const auto *sng_887 = buffer.data(sng + 887);
    const auto *sng_888 = buffer.data(sng + 888);
    const auto *sng_890 = buffer.data(sng + 890);
    const auto *sng_895 = buffer.data(sng + 895);
    const auto *sng_897 = buffer.data(sng + 897);
    const auto *sng_898 = buffer.data(sng + 898);
    const auto *sng_899 = buffer.data(sng + 899);
    const auto *sng_900 = buffer.data(sng + 900);
    const auto *sng_902 = buffer.data(sng + 902);
    const auto *sng_905 = buffer.data(sng + 905);
    const auto *sng_910 = buffer.data(sng + 910);
    const auto *sng_912 = buffer.data(sng + 912);
    const auto *sng_913 = buffer.data(sng + 913);
    const auto *sng_914 = buffer.data(sng + 914);
    const auto *sng_915 = buffer.data(sng + 915);

    const auto *snh1_1155 = buffer.data(snh1 + 1155);
    const auto *snh1_1158 = buffer.data(snh1 + 1158);
    const auto *snh1_1161 = buffer.data(snh1 + 1161);
    const auto *snh1_1170 = buffer.data(snh1 + 1170);
    const auto *snh1_1172 = buffer.data(snh1 + 1172);
    const auto *snh1_1173 = buffer.data(snh1 + 1173);

    const auto *sof0_663 = buffer.data(sof0 + 663);
    const auto *sof0_665 = buffer.data(sof0 + 665);
    const auto *sof0_666 = buffer.data(sof0 + 666);
    const auto *sof0_668 = buffer.data(sof0 + 668);
    const auto *sof0_669 = buffer.data(sof0 + 669);
    const auto *sof0_675 = buffer.data(sof0 + 675);
    const auto *sof0_679 = buffer.data(sof0 + 679);
    const auto *sof0_680 = buffer.data(sof0 + 680);
    const auto *sof0_683 = buffer.data(sof0 + 683);
    const auto *sof0_685 = buffer.data(sof0 + 685);
    const auto *sof0_686 = buffer.data(sof0 + 686);
    const auto *sof0_688 = buffer.data(sof0 + 688);
    const auto *sof0_689 = buffer.data(sof0 + 689);
    const auto *sof0_690 = buffer.data(sof0 + 690);
    const auto *sof0_693 = buffer.data(sof0 + 693);
    const auto *sof0_695 = buffer.data(sof0 + 695);
    const auto *sof0_696 = buffer.data(sof0 + 696);
    const auto *sof0_698 = buffer.data(sof0 + 698);
    const auto *sof0_699 = buffer.data(sof0 + 699);
    const auto *sof0_700 = buffer.data(sof0 + 700);
    const auto *sof0_703 = buffer.data(sof0 + 703);
    const auto *sof0_705 = buffer.data(sof0 + 705);
    const auto *sof0_706 = buffer.data(sof0 + 706);
    const auto *sof0_708 = buffer.data(sof0 + 708);
    const auto *sof0_709 = buffer.data(sof0 + 709);
    const auto *sof0_710 = buffer.data(sof0 + 710);
    const auto *sof0_713 = buffer.data(sof0 + 713);
    const auto *sof0_715 = buffer.data(sof0 + 715);
    const auto *sof0_716 = buffer.data(sof0 + 716);
    const auto *sof0_718 = buffer.data(sof0 + 718);
    const auto *sof0_719 = buffer.data(sof0 + 719);
    const auto *sof0_720 = buffer.data(sof0 + 720);

    const auto *sof1_663 = buffer.data(sof1 + 663);
    const auto *sof1_665 = buffer.data(sof1 + 665);
    const auto *sof1_666 = buffer.data(sof1 + 666);
    const auto *sof1_668 = buffer.data(sof1 + 668);
    const auto *sof1_669 = buffer.data(sof1 + 669);
    const auto *sof1_675 = buffer.data(sof1 + 675);
    const auto *sof1_679 = buffer.data(sof1 + 679);
    const auto *sof1_680 = buffer.data(sof1 + 680);
    const auto *sof1_683 = buffer.data(sof1 + 683);
    const auto *sof1_685 = buffer.data(sof1 + 685);
    const auto *sof1_686 = buffer.data(sof1 + 686);
    const auto *sof1_688 = buffer.data(sof1 + 688);
    const auto *sof1_689 = buffer.data(sof1 + 689);
    const auto *sof1_690 = buffer.data(sof1 + 690);
    const auto *sof1_693 = buffer.data(sof1 + 693);
    const auto *sof1_695 = buffer.data(sof1 + 695);
    const auto *sof1_696 = buffer.data(sof1 + 696);
    const auto *sof1_698 = buffer.data(sof1 + 698);
    const auto *sof1_699 = buffer.data(sof1 + 699);
    const auto *sof1_700 = buffer.data(sof1 + 700);
    const auto *sof1_703 = buffer.data(sof1 + 703);
    const auto *sof1_705 = buffer.data(sof1 + 705);
    const auto *sof1_706 = buffer.data(sof1 + 706);
    const auto *sof1_708 = buffer.data(sof1 + 708);
    const auto *sof1_709 = buffer.data(sof1 + 709);
    const auto *sof1_710 = buffer.data(sof1 + 710);
    const auto *sof1_713 = buffer.data(sof1 + 713);
    const auto *sof1_715 = buffer.data(sof1 + 715);
    const auto *sof1_716 = buffer.data(sof1 + 716);
    const auto *sof1_718 = buffer.data(sof1 + 718);
    const auto *sof1_719 = buffer.data(sof1 + 719);
    const auto *sof1_720 = buffer.data(sof1 + 720);

    const auto *sog_992 = buffer.data(sog + 992);
    const auto *sog_993 = buffer.data(sog + 993);
    const auto *sog_995 = buffer.data(sog + 995);
    const auto *sog_996 = buffer.data(sog + 996);
    const auto *sog_999 = buffer.data(sog + 999);
    const auto *sog_1000 = buffer.data(sog + 1000);
    const auto *sog_1001 = buffer.data(sog + 1001);
    const auto *sog_1002 = buffer.data(sog + 1002);
    const auto *sog_1003 = buffer.data(sog + 1003);
    const auto *sog_1004 = buffer.data(sog + 1004);
    const auto *sog_1005 = buffer.data(sog + 1005);
    const auto *sog_1007 = buffer.data(sog + 1007);
    const auto *sog_1008 = buffer.data(sog + 1008);
    const auto *sog_1010 = buffer.data(sog + 1010);
    const auto *sog_1014 = buffer.data(sog + 1014);
    const auto *sog_1015 = buffer.data(sog + 1015);
    const auto *sog_1016 = buffer.data(sog + 1016);
    const auto *sog_1017 = buffer.data(sog + 1017);
    const auto *sog_1018 = buffer.data(sog + 1018);
    const auto *sog_1019 = buffer.data(sog + 1019);
    const auto *sog_1020 = buffer.data(sog + 1020);
    const auto *sog_1022 = buffer.data(sog + 1022);
    const auto *sog_1023 = buffer.data(sog + 1023);
    const auto *sog_1025 = buffer.data(sog + 1025);
    const auto *sog_1026 = buffer.data(sog + 1026);
    const auto *sog_1029 = buffer.data(sog + 1029);
    const auto *sog_1030 = buffer.data(sog + 1030);
    const auto *sog_1031 = buffer.data(sog + 1031);
    const auto *sog_1032 = buffer.data(sog + 1032);
    const auto *sog_1033 = buffer.data(sog + 1033);
    const auto *sog_1034 = buffer.data(sog + 1034);
    const auto *sog_1035 = buffer.data(sog + 1035);
    const auto *sog_1037 = buffer.data(sog + 1037);
    const auto *sog_1038 = buffer.data(sog + 1038);
    const auto *sog_1040 = buffer.data(sog + 1040);
    const auto *sog_1041 = buffer.data(sog + 1041);
    const auto *sog_1044 = buffer.data(sog + 1044);
    const auto *sog_1045 = buffer.data(sog + 1045);
    const auto *sog_1046 = buffer.data(sog + 1046);
    const auto *sog_1047 = buffer.data(sog + 1047);
    const auto *sog_1048 = buffer.data(sog + 1048);
    const auto *sog_1049 = buffer.data(sog + 1049);
    const auto *sog_1050 = buffer.data(sog + 1050);
    const auto *sog_1052 = buffer.data(sog + 1052);
    const auto *sog_1053 = buffer.data(sog + 1053);
    const auto *sog_1055 = buffer.data(sog + 1055);
    const auto *sog_1056 = buffer.data(sog + 1056);
    const auto *sog_1059 = buffer.data(sog + 1059);
    const auto *sog_1060 = buffer.data(sog + 1060);
    const auto *sog_1061 = buffer.data(sog + 1061);
    const auto *sog_1062 = buffer.data(sog + 1062);
    const auto *sog_1063 = buffer.data(sog + 1063);
    const auto *sog_1064 = buffer.data(sog + 1064);
    const auto *sog_1065 = buffer.data(sog + 1065);
    const auto *sog_1067 = buffer.data(sog + 1067);
    const auto *sog_1068 = buffer.data(sog + 1068);
    const auto *sog_1070 = buffer.data(sog + 1070);
    const auto *sog_1071 = buffer.data(sog + 1071);
    const auto *sog_1074 = buffer.data(sog + 1074);
    const auto *sog_1075 = buffer.data(sog + 1075);
    const auto *sog_1076 = buffer.data(sog + 1076);
    const auto *sog_1077 = buffer.data(sog + 1077);
    const auto *sog_1078 = buffer.data(sog + 1078);
    const auto *sog_1079 = buffer.data(sog + 1079);
    const auto *sog_1080 = buffer.data(sog + 1080);

#pragma omp simd aligned(t_1389, t_1390, t_1391, pc_x, pc_y, sng_827, sof0_663, sof0_665, \
                         sof1_663, sof1_665, sog_992, sog_993, \
                         sog_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1389[k] = f_4 * sof0_663[k]
                    - f_5 * sof1_663[k]
                    + f_3 * pc_x[k] * sog_993[k];

        t_1390[k] = f_0 * sng_827[k]
                    + f_3 * pc_y[k] * sog_992[k];

        t_1391[k] = f_4 * sof0_665[k]
                    - f_5 * sof1_665[k]
                    + f_3 * pc_x[k] * sog_995[k];
    }

#pragma omp simd aligned(t_1392, t_1393, t_1394, t_1395, pc_x, pc_y, pc_z, sng_830, sof0_666, \
                         sof0_669, sof1_666, sof1_669, sog_993, sog_995, sog_996, \
                         sog_999 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1392[k] = f_6 * sof0_666[k]
                    - f_7 * sof1_666[k]
                    + f_3 * pc_x[k] * sog_996[k];

        t_1393[k] = f_3 * pc_z[k] * sog_993[k];

        t_1394[k] = f_0 * sng_830[k]
                    + f_3 * pc_y[k] * sog_995[k];

        t_1395[k] = f_6 * sof0_669[k]
                    - f_7 * sof1_669[k]
                    + f_3 * pc_x[k] * sog_999[k];
    }

#pragma omp simd aligned(t_1396, t_1397, t_1398, t_1399, t_1400, t_1401, pc_x, pc_y, sng_835, \
                         sof0_666, sof1_666, sog_1000, sog_1001, sog_1002, sog_1003, \
                         sog_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1396[k] = f_3 * pc_x[k] * sog_1000[k];

        t_1397[k] = f_3 * pc_x[k] * sog_1001[k];

        t_1398[k] = f_3 * pc_x[k] * sog_1002[k];

        t_1399[k] = f_3 * pc_x[k] * sog_1003[k];

        t_1400[k] = f_3 * pc_x[k] * sog_1004[k];

        t_1401[k] = f_0 * sng_835[k]
                    + f_1 * sof0_666[k]
                    - f_2 * sof1_666[k]
                    + f_3 * pc_y[k] * sog_1000[k];
    }

#pragma omp simd aligned(t_1402, t_1403, t_1404, pc_y, pc_z, sng_837, sng_838, sof0_668, \
                         sof0_669, sof1_668, sof1_669, sog_1000, sog_1002, \
                         sog_1003 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1402[k] = f_3 * pc_z[k] * sog_1000[k];

        t_1403[k] = f_0 * sng_837[k]
                    + f_4 * sof0_668[k]
                    - f_5 * sof1_668[k]
                    + f_3 * pc_y[k] * sog_1002[k];

        t_1404[k] = f_0 * sng_838[k]
                    + f_6 * sof0_669[k]
                    - f_7 * sof1_669[k]
                    + f_3 * pc_y[k] * sog_1003[k];
    }

#pragma omp simd aligned(t_1405, t_1406, t_1407, t_1408, pb_z, pc_y, pc_z, snh0_1155, sng_839, \
                         sng_840, snh1_1155, sof0_669, sof1_669, sog_1004, \
                         sog_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1405[k] = f_0 * sng_839[k]
                    + f_3 * pc_y[k] * sog_1004[k];

        t_1406[k] = f_1 * sof0_669[k]
                    - f_2 * sof1_669[k]
                    + f_3 * pc_z[k] * sog_1004[k];

        t_1407[k] = pb_z[k] * snh0_1155[k]
                    - f_8 * pc_z[k] * snh1_1155[k];

        t_1408[k] = f_12 * sng_840[k]
                    + f_3 * pc_y[k] * sog_1005[k];
    }

#pragma omp simd aligned(t_1409, t_1410, t_1411, pb_z, pc_y, pc_z, snh0_1158, sng_825, \
                         sng_842, snh1_1158, sog_1005, sog_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1409[k] = f_9 * sng_825[k]
                    + f_3 * pc_z[k] * sog_1005[k];

        t_1410[k] = pb_z[k] * snh0_1158[k]
                    - f_8 * pc_z[k] * snh1_1158[k];

        t_1411[k] = f_12 * sng_842[k]
                    + f_3 * pc_y[k] * sog_1007[k];
    }

#pragma omp simd aligned(t_1412, t_1413, t_1414, t_1415, pb_z, pc_x, pc_y, pc_z, snh0_1161, \
                         sng_828, sng_845, snh1_1161, sof0_675, sof1_675, sog_1008, \
                         sog_1010 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1412[k] = f_4 * sof0_675[k]
                    - f_5 * sof1_675[k]
                    + f_3 * pc_x[k] * sog_1010[k];

        t_1413[k] = pb_z[k] * snh0_1161[k]
                    - f_8 * pc_z[k] * snh1_1161[k];

        t_1414[k] = f_9 * sng_828[k]
                    + f_3 * pc_z[k] * sog_1008[k];

        t_1415[k] = f_12 * sng_845[k]
                    + f_3 * pc_y[k] * sog_1010[k];
    }

#pragma omp simd aligned(t_1416, t_1417, t_1418, t_1419, t_1420, t_1421, pc_x, sof0_679, \
                         sof1_679, sog_1014, sog_1015, sog_1016, sog_1017, sog_1018, \
                         sog_1019 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1416[k] = f_6 * sof0_679[k]
                    - f_7 * sof1_679[k]
                    + f_3 * pc_x[k] * sog_1014[k];

        t_1417[k] = f_3 * pc_x[k] * sog_1015[k];

        t_1418[k] = f_3 * pc_x[k] * sog_1016[k];

        t_1419[k] = f_3 * pc_x[k] * sog_1017[k];

        t_1420[k] = f_3 * pc_x[k] * sog_1018[k];

        t_1421[k] = f_3 * pc_x[k] * sog_1019[k];
    }

#pragma omp simd aligned(t_1422, t_1423, t_1424, t_1425, pb_z, pc_z, snh0_1170, snh0_1172, \
                         snh0_1173, sng_835, sng_836, sng_837, snh1_1170, snh1_1172, \
                         snh1_1173, sog_1015 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1422[k] = pb_z[k] * snh0_1170[k]
                    - f_8 * pc_z[k] * snh1_1170[k];

        t_1423[k] = f_9 * sng_835[k]
                    + f_3 * pc_z[k] * sog_1015[k];

        t_1424[k] = pb_z[k] * snh0_1172[k]
                    + f_10 * sng_836[k]
                    - f_8 * pc_z[k] * snh1_1172[k];

        t_1425[k] = pb_z[k] * snh0_1173[k]
                    + f_11 * sng_837[k]
                    - f_8 * pc_z[k] * snh1_1173[k];
    }

#pragma omp simd aligned(t_1426, t_1427, t_1428, t_1429, pc_x, pc_y, pc_z, sng_839, sng_854, \
                         sng_855, sof0_679, sof0_680, sof1_679, sof1_680, sog_1019, \
                         sog_1020 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1426[k] = f_12 * sng_854[k]
                    + f_3 * pc_y[k] * sog_1019[k];

        t_1427[k] = f_9 * sng_839[k]
                    + f_1 * sof0_679[k]
                    - f_2 * sof1_679[k]
                    + f_3 * pc_z[k] * sog_1019[k];

        t_1428[k] = f_1 * sof0_680[k]
                    - f_2 * sof1_680[k]
                    + f_3 * pc_x[k] * sog_1020[k];

        t_1429[k] = f_13 * sng_855[k]
                    + f_3 * pc_y[k] * sog_1020[k];
    }

#pragma omp simd aligned(t_1430, t_1431, t_1432, pc_x, pc_y, pc_z, sng_840, sng_857, sof0_683, \
                         sof1_683, sog_1020, sog_1022, sog_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1430[k] = f_10 * sng_840[k]
                    + f_3 * pc_z[k] * sog_1020[k];

        t_1431[k] = f_4 * sof0_683[k]
                    - f_5 * sof1_683[k]
                    + f_3 * pc_x[k] * sog_1023[k];

        t_1432[k] = f_13 * sng_857[k]
                    + f_3 * pc_y[k] * sog_1022[k];
    }

#pragma omp simd aligned(t_1433, t_1434, t_1435, t_1436, pc_x, pc_y, pc_z, sng_843, sng_860, \
                         sof0_685, sof0_686, sof1_685, sof1_686, sog_1023, sog_1025, \
                         sog_1026 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1433[k] = f_4 * sof0_685[k]
                    - f_5 * sof1_685[k]
                    + f_3 * pc_x[k] * sog_1025[k];

        t_1434[k] = f_6 * sof0_686[k]
                    - f_7 * sof1_686[k]
                    + f_3 * pc_x[k] * sog_1026[k];

        t_1435[k] = f_10 * sng_843[k]
                    + f_3 * pc_z[k] * sog_1023[k];

        t_1436[k] = f_13 * sng_860[k]
                    + f_3 * pc_y[k] * sog_1025[k];
    }

#pragma omp simd aligned(t_1437, t_1438, t_1439, t_1440, t_1441, t_1442, pc_x, sof0_689, \
                         sof1_689, sog_1029, sog_1030, sog_1031, sog_1032, sog_1033, \
                         sog_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1437[k] = f_6 * sof0_689[k]
                    - f_7 * sof1_689[k]
                    + f_3 * pc_x[k] * sog_1029[k];

        t_1438[k] = f_3 * pc_x[k] * sog_1030[k];

        t_1439[k] = f_3 * pc_x[k] * sog_1031[k];

        t_1440[k] = f_3 * pc_x[k] * sog_1032[k];

        t_1441[k] = f_3 * pc_x[k] * sog_1033[k];

        t_1442[k] = f_3 * pc_x[k] * sog_1034[k];
    }

#pragma omp simd aligned(t_1443, t_1444, t_1445, pc_y, pc_z, sng_850, sng_865, sng_867, \
                         sof0_686, sof0_688, sof1_686, sof1_688, sog_1030, \
                         sog_1032 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1443[k] = f_13 * sng_865[k]
                    + f_1 * sof0_686[k]
                    - f_2 * sof1_686[k]
                    + f_3 * pc_y[k] * sog_1030[k];

        t_1444[k] = f_10 * sng_850[k]
                    + f_3 * pc_z[k] * sog_1030[k];

        t_1445[k] = f_13 * sng_867[k]
                    + f_4 * sof0_688[k]
                    - f_5 * sof1_688[k]
                    + f_3 * pc_y[k] * sog_1032[k];
    }

#pragma omp simd aligned(t_1446, t_1447, t_1448, pc_y, pc_z, sng_854, sng_868, sng_869, \
                         sof0_689, sof1_689, sog_1033, sog_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1446[k] = f_13 * sng_868[k]
                    + f_6 * sof0_689[k]
                    - f_7 * sof1_689[k]
                    + f_3 * pc_y[k] * sog_1033[k];

        t_1447[k] = f_13 * sng_869[k]
                    + f_3 * pc_y[k] * sog_1034[k];

        t_1448[k] = f_10 * sng_854[k]
                    + f_1 * sof0_689[k]
                    - f_2 * sof1_689[k]
                    + f_3 * pc_z[k] * sog_1034[k];
    }

#pragma omp simd aligned(t_1449, t_1450, t_1451, t_1452, pc_x, pc_y, pc_z, sng_855, sng_870, \
                         sof0_690, sof0_693, sof1_690, sof1_693, sog_1035, \
                         sog_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1449[k] = f_1 * sof0_690[k]
                    - f_2 * sof1_690[k]
                    + f_3 * pc_x[k] * sog_1035[k];

        t_1450[k] = f_14 * sng_870[k]
                    + f_3 * pc_y[k] * sog_1035[k];

        t_1451[k] = f_11 * sng_855[k]
                    + f_3 * pc_z[k] * sog_1035[k];

        t_1452[k] = f_4 * sof0_693[k]
                    - f_5 * sof1_693[k]
                    + f_3 * pc_x[k] * sog_1038[k];
    }

#pragma omp simd aligned(t_1453, t_1454, t_1455, pc_x, pc_y, sng_872, sof0_695, sof0_696, \
                         sof1_695, sof1_696, sog_1037, sog_1040, \
                         sog_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1453[k] = f_14 * sng_872[k]
                    + f_3 * pc_y[k] * sog_1037[k];

        t_1454[k] = f_4 * sof0_695[k]
                    - f_5 * sof1_695[k]
                    + f_3 * pc_x[k] * sog_1040[k];

        t_1455[k] = f_6 * sof0_696[k]
                    - f_7 * sof1_696[k]
                    + f_3 * pc_x[k] * sog_1041[k];
    }

#pragma omp simd aligned(t_1456, t_1457, t_1458, t_1459, pc_x, pc_y, pc_z, sng_858, sng_875, \
                         sof0_699, sof1_699, sog_1038, sog_1040, sog_1044, \
                         sog_1045 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1456[k] = f_11 * sng_858[k]
                    + f_3 * pc_z[k] * sog_1038[k];

        t_1457[k] = f_14 * sng_875[k]
                    + f_3 * pc_y[k] * sog_1040[k];

        t_1458[k] = f_6 * sof0_699[k]
                    - f_7 * sof1_699[k]
                    + f_3 * pc_x[k] * sog_1044[k];

        t_1459[k] = f_3 * pc_x[k] * sog_1045[k];
    }

#pragma omp simd aligned(t_1460, t_1461, t_1462, t_1463, t_1464, pc_x, pc_y, sng_880, \
                         sof0_696, sof1_696, sog_1045, sog_1046, sog_1047, sog_1048, \
                         sog_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1460[k] = f_3 * pc_x[k] * sog_1046[k];

        t_1461[k] = f_3 * pc_x[k] * sog_1047[k];

        t_1462[k] = f_3 * pc_x[k] * sog_1048[k];

        t_1463[k] = f_3 * pc_x[k] * sog_1049[k];

        t_1464[k] = f_14 * sng_880[k]
                    + f_1 * sof0_696[k]
                    - f_2 * sof1_696[k]
                    + f_3 * pc_y[k] * sog_1045[k];
    }

#pragma omp simd aligned(t_1465, t_1466, t_1467, pc_y, pc_z, sng_865, sng_882, sng_883, \
                         sof0_698, sof0_699, sof1_698, sof1_699, sog_1045, sog_1047, \
                         sog_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1465[k] = f_11 * sng_865[k]
                    + f_3 * pc_z[k] * sog_1045[k];

        t_1466[k] = f_14 * sng_882[k]
                    + f_4 * sof0_698[k]
                    - f_5 * sof1_698[k]
                    + f_3 * pc_y[k] * sog_1047[k];

        t_1467[k] = f_14 * sng_883[k]
                    + f_6 * sof0_699[k]
                    - f_7 * sof1_699[k]
                    + f_3 * pc_y[k] * sog_1048[k];
    }

#pragma omp simd aligned(t_1468, t_1469, t_1470, t_1471, pc_x, pc_y, pc_z, sng_869, sng_884, \
                         sng_885, sof0_699, sof0_700, sof1_699, sof1_700, sog_1049, \
                         sog_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1468[k] = f_14 * sng_884[k]
                    + f_3 * pc_y[k] * sog_1049[k];

        t_1469[k] = f_11 * sng_869[k]
                    + f_1 * sof0_699[k]
                    - f_2 * sof1_699[k]
                    + f_3 * pc_z[k] * sog_1049[k];

        t_1470[k] = f_1 * sof0_700[k]
                    - f_2 * sof1_700[k]
                    + f_3 * pc_x[k] * sog_1050[k];

        t_1471[k] = f_15 * sng_885[k]
                    + f_3 * pc_y[k] * sog_1050[k];
    }

#pragma omp simd aligned(t_1472, t_1473, t_1474, pc_x, pc_y, pc_z, sng_870, sng_887, sof0_703, \
                         sof1_703, sog_1050, sog_1052, sog_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1472[k] = f_16 * sng_870[k]
                    + f_3 * pc_z[k] * sog_1050[k];

        t_1473[k] = f_4 * sof0_703[k]
                    - f_5 * sof1_703[k]
                    + f_3 * pc_x[k] * sog_1053[k];

        t_1474[k] = f_15 * sng_887[k]
                    + f_3 * pc_y[k] * sog_1052[k];
    }

#pragma omp simd aligned(t_1475, t_1476, t_1477, t_1478, pc_x, pc_y, pc_z, sng_873, sng_890, \
                         sof0_705, sof0_706, sof1_705, sof1_706, sog_1053, sog_1055, \
                         sog_1056 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1475[k] = f_4 * sof0_705[k]
                    - f_5 * sof1_705[k]
                    + f_3 * pc_x[k] * sog_1055[k];

        t_1476[k] = f_6 * sof0_706[k]
                    - f_7 * sof1_706[k]
                    + f_3 * pc_x[k] * sog_1056[k];

        t_1477[k] = f_16 * sng_873[k]
                    + f_3 * pc_z[k] * sog_1053[k];

        t_1478[k] = f_15 * sng_890[k]
                    + f_3 * pc_y[k] * sog_1055[k];
    }

#pragma omp simd aligned(t_1479, t_1480, t_1481, t_1482, t_1483, t_1484, pc_x, sof0_709, \
                         sof1_709, sog_1059, sog_1060, sog_1061, sog_1062, sog_1063, \
                         sog_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1479[k] = f_6 * sof0_709[k]
                    - f_7 * sof1_709[k]
                    + f_3 * pc_x[k] * sog_1059[k];

        t_1480[k] = f_3 * pc_x[k] * sog_1060[k];

        t_1481[k] = f_3 * pc_x[k] * sog_1061[k];

        t_1482[k] = f_3 * pc_x[k] * sog_1062[k];

        t_1483[k] = f_3 * pc_x[k] * sog_1063[k];

        t_1484[k] = f_3 * pc_x[k] * sog_1064[k];
    }

#pragma omp simd aligned(t_1485, t_1486, t_1487, pc_y, pc_z, sng_880, sng_895, sng_897, \
                         sof0_706, sof0_708, sof1_706, sof1_708, sog_1060, \
                         sog_1062 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1485[k] = f_15 * sng_895[k]
                    + f_1 * sof0_706[k]
                    - f_2 * sof1_706[k]
                    + f_3 * pc_y[k] * sog_1060[k];

        t_1486[k] = f_16 * sng_880[k]
                    + f_3 * pc_z[k] * sog_1060[k];

        t_1487[k] = f_15 * sng_897[k]
                    + f_4 * sof0_708[k]
                    - f_5 * sof1_708[k]
                    + f_3 * pc_y[k] * sog_1062[k];
    }

#pragma omp simd aligned(t_1488, t_1489, t_1490, pc_y, pc_z, sng_884, sng_898, sng_899, \
                         sof0_709, sof1_709, sog_1063, sog_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1488[k] = f_15 * sng_898[k]
                    + f_6 * sof0_709[k]
                    - f_7 * sof1_709[k]
                    + f_3 * pc_y[k] * sog_1063[k];

        t_1489[k] = f_15 * sng_899[k]
                    + f_3 * pc_y[k] * sog_1064[k];

        t_1490[k] = f_16 * sng_884[k]
                    + f_1 * sof0_709[k]
                    - f_2 * sof1_709[k]
                    + f_3 * pc_z[k] * sog_1064[k];
    }

#pragma omp simd aligned(t_1491, t_1492, t_1493, t_1494, pc_x, pc_y, pc_z, sng_885, sng_900, \
                         sof0_710, sof0_713, sof1_710, sof1_713, sog_1065, \
                         sog_1068 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1491[k] = f_1 * sof0_710[k]
                    - f_2 * sof1_710[k]
                    + f_3 * pc_x[k] * sog_1065[k];

        t_1492[k] = f_17 * sng_900[k]
                    + f_3 * pc_y[k] * sog_1065[k];

        t_1493[k] = f_18 * sng_885[k]
                    + f_3 * pc_z[k] * sog_1065[k];

        t_1494[k] = f_4 * sof0_713[k]
                    - f_5 * sof1_713[k]
                    + f_3 * pc_x[k] * sog_1068[k];
    }

#pragma omp simd aligned(t_1495, t_1496, t_1497, pc_x, pc_y, sng_902, sof0_715, sof0_716, \
                         sof1_715, sof1_716, sog_1067, sog_1070, \
                         sog_1071 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1495[k] = f_17 * sng_902[k]
                    + f_3 * pc_y[k] * sog_1067[k];

        t_1496[k] = f_4 * sof0_715[k]
                    - f_5 * sof1_715[k]
                    + f_3 * pc_x[k] * sog_1070[k];

        t_1497[k] = f_6 * sof0_716[k]
                    - f_7 * sof1_716[k]
                    + f_3 * pc_x[k] * sog_1071[k];
    }

#pragma omp simd aligned(t_1498, t_1499, t_1500, t_1501, pc_x, pc_y, pc_z, sng_888, sng_905, \
                         sof0_719, sof1_719, sog_1068, sog_1070, sog_1074, \
                         sog_1075 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1498[k] = f_18 * sng_888[k]
                    + f_3 * pc_z[k] * sog_1068[k];

        t_1499[k] = f_17 * sng_905[k]
                    + f_3 * pc_y[k] * sog_1070[k];

        t_1500[k] = f_6 * sof0_719[k]
                    - f_7 * sof1_719[k]
                    + f_3 * pc_x[k] * sog_1074[k];

        t_1501[k] = f_3 * pc_x[k] * sog_1075[k];
    }

#pragma omp simd aligned(t_1502, t_1503, t_1504, t_1505, t_1506, pc_x, pc_y, sng_910, \
                         sof0_716, sof1_716, sog_1075, sog_1076, sog_1077, sog_1078, \
                         sog_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1502[k] = f_3 * pc_x[k] * sog_1076[k];

        t_1503[k] = f_3 * pc_x[k] * sog_1077[k];

        t_1504[k] = f_3 * pc_x[k] * sog_1078[k];

        t_1505[k] = f_3 * pc_x[k] * sog_1079[k];

        t_1506[k] = f_17 * sng_910[k]
                    + f_1 * sof0_716[k]
                    - f_2 * sof1_716[k]
                    + f_3 * pc_y[k] * sog_1075[k];
    }

#pragma omp simd aligned(t_1507, t_1508, t_1509, pc_y, pc_z, sng_895, sng_912, sng_913, \
                         sof0_718, sof0_719, sof1_718, sof1_719, sog_1075, sog_1077, \
                         sog_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1507[k] = f_18 * sng_895[k]
                    + f_3 * pc_z[k] * sog_1075[k];

        t_1508[k] = f_17 * sng_912[k]
                    + f_4 * sof0_718[k]
                    - f_5 * sof1_718[k]
                    + f_3 * pc_y[k] * sog_1077[k];

        t_1509[k] = f_17 * sng_913[k]
                    + f_6 * sof0_719[k]
                    - f_7 * sof1_719[k]
                    + f_3 * pc_y[k] * sog_1078[k];
    }

#pragma omp simd aligned(t_1510, t_1511, t_1512, t_1513, pc_x, pc_y, pc_z, sng_899, sng_914, \
                         sng_915, sof0_719, sof0_720, sof1_719, sof1_720, sog_1079, \
                         sog_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1510[k] = f_17 * sng_914[k]
                    + f_3 * pc_y[k] * sog_1079[k];

        t_1511[k] = f_18 * sng_899[k]
                    + f_1 * sof0_719[k]
                    - f_2 * sof1_719[k]
                    + f_3 * pc_z[k] * sog_1079[k];

        t_1512[k] = f_1 * sof0_720[k]
                    - f_2 * sof1_720[k]
                    + f_3 * pc_x[k] * sog_1080[k];

        t_1513[k] = f_18 * sng_915[k]
                    + f_3 * pc_y[k] * sog_1080[k];
    }
}

static auto
compute_prim_soh_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t snh0,
                                                           const size_t sng, const size_t snh1,
                                                           const size_t sof0, const size_t sof1,
                                                           const size_t sog, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 5.0 / q;
    const auto f_13 = 4.5 / q;
    const auto f_14 = 4.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_1514 = buffer.data(target + 1514);
    auto *t_1515 = buffer.data(target + 1515);
    auto *t_1516 = buffer.data(target + 1516);
    auto *t_1517 = buffer.data(target + 1517);
    auto *t_1518 = buffer.data(target + 1518);
    auto *t_1519 = buffer.data(target + 1519);
    auto *t_1520 = buffer.data(target + 1520);
    auto *t_1521 = buffer.data(target + 1521);
    auto *t_1522 = buffer.data(target + 1522);
    auto *t_1523 = buffer.data(target + 1523);
    auto *t_1524 = buffer.data(target + 1524);
    auto *t_1525 = buffer.data(target + 1525);
    auto *t_1526 = buffer.data(target + 1526);
    auto *t_1527 = buffer.data(target + 1527);
    auto *t_1528 = buffer.data(target + 1528);
    auto *t_1529 = buffer.data(target + 1529);
    auto *t_1530 = buffer.data(target + 1530);
    auto *t_1531 = buffer.data(target + 1531);
    auto *t_1532 = buffer.data(target + 1532);
    auto *t_1533 = buffer.data(target + 1533);
    auto *t_1534 = buffer.data(target + 1534);
    auto *t_1535 = buffer.data(target + 1535);
    auto *t_1536 = buffer.data(target + 1536);
    auto *t_1537 = buffer.data(target + 1537);
    auto *t_1538 = buffer.data(target + 1538);
    auto *t_1539 = buffer.data(target + 1539);
    auto *t_1540 = buffer.data(target + 1540);
    auto *t_1541 = buffer.data(target + 1541);
    auto *t_1542 = buffer.data(target + 1542);
    auto *t_1543 = buffer.data(target + 1543);
    auto *t_1544 = buffer.data(target + 1544);
    auto *t_1545 = buffer.data(target + 1545);
    auto *t_1546 = buffer.data(target + 1546);
    auto *t_1547 = buffer.data(target + 1547);
    auto *t_1548 = buffer.data(target + 1548);
    auto *t_1549 = buffer.data(target + 1549);
    auto *t_1550 = buffer.data(target + 1550);
    auto *t_1551 = buffer.data(target + 1551);
    auto *t_1552 = buffer.data(target + 1552);
    auto *t_1553 = buffer.data(target + 1553);
    auto *t_1554 = buffer.data(target + 1554);
    auto *t_1555 = buffer.data(target + 1555);
    auto *t_1556 = buffer.data(target + 1556);
    auto *t_1557 = buffer.data(target + 1557);
    auto *t_1558 = buffer.data(target + 1558);
    auto *t_1559 = buffer.data(target + 1559);
    auto *t_1560 = buffer.data(target + 1560);
    auto *t_1561 = buffer.data(target + 1561);
    auto *t_1562 = buffer.data(target + 1562);
    auto *t_1563 = buffer.data(target + 1563);
    auto *t_1564 = buffer.data(target + 1564);
    auto *t_1565 = buffer.data(target + 1565);
    auto *t_1566 = buffer.data(target + 1566);
    auto *t_1567 = buffer.data(target + 1567);
    auto *t_1568 = buffer.data(target + 1568);
    auto *t_1569 = buffer.data(target + 1569);
    auto *t_1570 = buffer.data(target + 1570);
    auto *t_1571 = buffer.data(target + 1571);
    auto *t_1572 = buffer.data(target + 1572);
    auto *t_1573 = buffer.data(target + 1573);
    auto *t_1574 = buffer.data(target + 1574);
    auto *t_1575 = buffer.data(target + 1575);
    auto *t_1576 = buffer.data(target + 1576);
    auto *t_1577 = buffer.data(target + 1577);
    auto *t_1578 = buffer.data(target + 1578);
    auto *t_1579 = buffer.data(target + 1579);
    auto *t_1580 = buffer.data(target + 1580);
    auto *t_1581 = buffer.data(target + 1581);
    auto *t_1582 = buffer.data(target + 1582);
    auto *t_1583 = buffer.data(target + 1583);
    auto *t_1584 = buffer.data(target + 1584);
    auto *t_1585 = buffer.data(target + 1585);
    auto *t_1586 = buffer.data(target + 1586);
    auto *t_1587 = buffer.data(target + 1587);
    auto *t_1588 = buffer.data(target + 1588);
    auto *t_1589 = buffer.data(target + 1589);
    auto *t_1590 = buffer.data(target + 1590);
    auto *t_1591 = buffer.data(target + 1591);
    auto *t_1592 = buffer.data(target + 1592);
    auto *t_1593 = buffer.data(target + 1593);
    auto *t_1594 = buffer.data(target + 1594);
    auto *t_1595 = buffer.data(target + 1595);
    auto *t_1596 = buffer.data(target + 1596);
    auto *t_1597 = buffer.data(target + 1597);
    auto *t_1598 = buffer.data(target + 1598);
    auto *t_1599 = buffer.data(target + 1599);
    auto *t_1600 = buffer.data(target + 1600);
    auto *t_1601 = buffer.data(target + 1601);
    auto *t_1602 = buffer.data(target + 1602);
    auto *t_1603 = buffer.data(target + 1603);
    auto *t_1604 = buffer.data(target + 1604);
    auto *t_1605 = buffer.data(target + 1605);
    auto *t_1606 = buffer.data(target + 1606);
    auto *t_1607 = buffer.data(target + 1607);
    auto *t_1608 = buffer.data(target + 1608);
    auto *t_1609 = buffer.data(target + 1609);
    auto *t_1610 = buffer.data(target + 1610);
    auto *t_1611 = buffer.data(target + 1611);
    auto *t_1612 = buffer.data(target + 1612);
    auto *t_1613 = buffer.data(target + 1613);
    auto *t_1614 = buffer.data(target + 1614);
    auto *t_1615 = buffer.data(target + 1615);
    auto *t_1616 = buffer.data(target + 1616);
    auto *t_1617 = buffer.data(target + 1617);
    auto *t_1618 = buffer.data(target + 1618);
    auto *t_1619 = buffer.data(target + 1619);
    auto *t_1620 = buffer.data(target + 1620);
    auto *t_1621 = buffer.data(target + 1621);
    auto *t_1622 = buffer.data(target + 1622);
    auto *t_1623 = buffer.data(target + 1623);
    auto *t_1624 = buffer.data(target + 1624);
    auto *t_1625 = buffer.data(target + 1625);
    auto *t_1626 = buffer.data(target + 1626);
    auto *t_1627 = buffer.data(target + 1627);
    auto *t_1628 = buffer.data(target + 1628);
    auto *t_1629 = buffer.data(target + 1629);
    auto *t_1630 = buffer.data(target + 1630);
    auto *t_1631 = buffer.data(target + 1631);
    auto *t_1632 = buffer.data(target + 1632);
    auto *t_1633 = buffer.data(target + 1633);
    auto *t_1634 = buffer.data(target + 1634);
    auto *t_1635 = buffer.data(target + 1635);
    auto *t_1636 = buffer.data(target + 1636);
    auto *t_1637 = buffer.data(target + 1637);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh0_1365 = buffer.data(snh0 + 1365);
    const auto *snh0_1370 = buffer.data(snh0 + 1370);
    const auto *snh0_1374 = buffer.data(snh0 + 1374);
    const auto *snh0_1380 = buffer.data(snh0 + 1380);
    const auto *snh0_1382 = buffer.data(snh0 + 1382);
    const auto *snh0_1383 = buffer.data(snh0 + 1383);
    const auto *snh0_1385 = buffer.data(snh0 + 1385);

    const auto *sng_900 = buffer.data(sng + 900);
    const auto *sng_903 = buffer.data(sng + 903);
    const auto *sng_910 = buffer.data(sng + 910);
    const auto *sng_914 = buffer.data(sng + 914);
    const auto *sng_915 = buffer.data(sng + 915);
    const auto *sng_917 = buffer.data(sng + 917);
    const auto *sng_918 = buffer.data(sng + 918);
    const auto *sng_920 = buffer.data(sng + 920);
    const auto *sng_925 = buffer.data(sng + 925);
    const auto *sng_927 = buffer.data(sng + 927);
    const auto *sng_928 = buffer.data(sng + 928);
    const auto *sng_929 = buffer.data(sng + 929);
    const auto *sng_930 = buffer.data(sng + 930);
    const auto *sng_932 = buffer.data(sng + 932);
    const auto *sng_933 = buffer.data(sng + 933);
    const auto *sng_935 = buffer.data(sng + 935);
    const auto *sng_940 = buffer.data(sng + 940);
    const auto *sng_942 = buffer.data(sng + 942);
    const auto *sng_943 = buffer.data(sng + 943);
    const auto *sng_944 = buffer.data(sng + 944);
    const auto *sng_945 = buffer.data(sng + 945);
    const auto *sng_947 = buffer.data(sng + 947);
    const auto *sng_948 = buffer.data(sng + 948);
    const auto *sng_950 = buffer.data(sng + 950);
    const auto *sng_955 = buffer.data(sng + 955);
    const auto *sng_957 = buffer.data(sng + 957);
    const auto *sng_958 = buffer.data(sng + 958);
    const auto *sng_959 = buffer.data(sng + 959);
    const auto *sng_960 = buffer.data(sng + 960);
    const auto *sng_962 = buffer.data(sng + 962);
    const auto *sng_963 = buffer.data(sng + 963);
    const auto *sng_965 = buffer.data(sng + 965);
    const auto *sng_970 = buffer.data(sng + 970);
    const auto *sng_972 = buffer.data(sng + 972);
    const auto *sng_973 = buffer.data(sng + 973);
    const auto *sng_974 = buffer.data(sng + 974);
    const auto *sng_975 = buffer.data(sng + 975);
    const auto *sng_977 = buffer.data(sng + 977);
    const auto *sng_978 = buffer.data(sng + 978);
    const auto *sng_980 = buffer.data(sng + 980);
    const auto *sng_985 = buffer.data(sng + 985);
    const auto *sng_987 = buffer.data(sng + 987);
    const auto *sng_988 = buffer.data(sng + 988);
    const auto *sng_989 = buffer.data(sng + 989);

    const auto *snh1_1365 = buffer.data(snh1 + 1365);
    const auto *snh1_1370 = buffer.data(snh1 + 1370);
    const auto *snh1_1374 = buffer.data(snh1 + 1374);
    const auto *snh1_1380 = buffer.data(snh1 + 1380);
    const auto *snh1_1382 = buffer.data(snh1 + 1382);
    const auto *snh1_1383 = buffer.data(snh1 + 1383);
    const auto *snh1_1385 = buffer.data(snh1 + 1385);

    const auto *sof0_723 = buffer.data(sof0 + 723);
    const auto *sof0_725 = buffer.data(sof0 + 725);
    const auto *sof0_726 = buffer.data(sof0 + 726);
    const auto *sof0_728 = buffer.data(sof0 + 728);
    const auto *sof0_729 = buffer.data(sof0 + 729);
    const auto *sof0_730 = buffer.data(sof0 + 730);
    const auto *sof0_733 = buffer.data(sof0 + 733);
    const auto *sof0_735 = buffer.data(sof0 + 735);
    const auto *sof0_736 = buffer.data(sof0 + 736);
    const auto *sof0_738 = buffer.data(sof0 + 738);
    const auto *sof0_739 = buffer.data(sof0 + 739);
    const auto *sof0_740 = buffer.data(sof0 + 740);
    const auto *sof0_743 = buffer.data(sof0 + 743);
    const auto *sof0_745 = buffer.data(sof0 + 745);
    const auto *sof0_746 = buffer.data(sof0 + 746);
    const auto *sof0_748 = buffer.data(sof0 + 748);
    const auto *sof0_749 = buffer.data(sof0 + 749);
    const auto *sof0_750 = buffer.data(sof0 + 750);
    const auto *sof0_753 = buffer.data(sof0 + 753);
    const auto *sof0_755 = buffer.data(sof0 + 755);
    const auto *sof0_756 = buffer.data(sof0 + 756);
    const auto *sof0_758 = buffer.data(sof0 + 758);
    const auto *sof0_759 = buffer.data(sof0 + 759);
    const auto *sof0_763 = buffer.data(sof0 + 763);
    const auto *sof0_766 = buffer.data(sof0 + 766);
    const auto *sof0_770 = buffer.data(sof0 + 770);
    const auto *sof0_773 = buffer.data(sof0 + 773);
    const auto *sof0_775 = buffer.data(sof0 + 775);
    const auto *sof0_776 = buffer.data(sof0 + 776);
    const auto *sof0_778 = buffer.data(sof0 + 778);
    const auto *sof0_779 = buffer.data(sof0 + 779);

    const auto *sof1_723 = buffer.data(sof1 + 723);
    const auto *sof1_725 = buffer.data(sof1 + 725);
    const auto *sof1_726 = buffer.data(sof1 + 726);
    const auto *sof1_728 = buffer.data(sof1 + 728);
    const auto *sof1_729 = buffer.data(sof1 + 729);
    const auto *sof1_730 = buffer.data(sof1 + 730);
    const auto *sof1_733 = buffer.data(sof1 + 733);
    const auto *sof1_735 = buffer.data(sof1 + 735);
    const auto *sof1_736 = buffer.data(sof1 + 736);
    const auto *sof1_738 = buffer.data(sof1 + 738);
    const auto *sof1_739 = buffer.data(sof1 + 739);
    const auto *sof1_740 = buffer.data(sof1 + 740);
    const auto *sof1_743 = buffer.data(sof1 + 743);
    const auto *sof1_745 = buffer.data(sof1 + 745);
    const auto *sof1_746 = buffer.data(sof1 + 746);
    const auto *sof1_748 = buffer.data(sof1 + 748);
    const auto *sof1_749 = buffer.data(sof1 + 749);
    const auto *sof1_750 = buffer.data(sof1 + 750);
    const auto *sof1_753 = buffer.data(sof1 + 753);
    const auto *sof1_755 = buffer.data(sof1 + 755);
    const auto *sof1_756 = buffer.data(sof1 + 756);
    const auto *sof1_758 = buffer.data(sof1 + 758);
    const auto *sof1_759 = buffer.data(sof1 + 759);
    const auto *sof1_763 = buffer.data(sof1 + 763);
    const auto *sof1_766 = buffer.data(sof1 + 766);
    const auto *sof1_770 = buffer.data(sof1 + 770);
    const auto *sof1_773 = buffer.data(sof1 + 773);
    const auto *sof1_775 = buffer.data(sof1 + 775);
    const auto *sof1_776 = buffer.data(sof1 + 776);
    const auto *sof1_778 = buffer.data(sof1 + 778);
    const auto *sof1_779 = buffer.data(sof1 + 779);

    const auto *sog_1080 = buffer.data(sog + 1080);
    const auto *sog_1082 = buffer.data(sog + 1082);
    const auto *sog_1083 = buffer.data(sog + 1083);
    const auto *sog_1085 = buffer.data(sog + 1085);
    const auto *sog_1086 = buffer.data(sog + 1086);
    const auto *sog_1089 = buffer.data(sog + 1089);
    const auto *sog_1090 = buffer.data(sog + 1090);
    const auto *sog_1091 = buffer.data(sog + 1091);
    const auto *sog_1092 = buffer.data(sog + 1092);
    const auto *sog_1093 = buffer.data(sog + 1093);
    const auto *sog_1094 = buffer.data(sog + 1094);
    const auto *sog_1095 = buffer.data(sog + 1095);
    const auto *sog_1097 = buffer.data(sog + 1097);
    const auto *sog_1098 = buffer.data(sog + 1098);
    const auto *sog_1100 = buffer.data(sog + 1100);
    const auto *sog_1101 = buffer.data(sog + 1101);
    const auto *sog_1104 = buffer.data(sog + 1104);
    const auto *sog_1105 = buffer.data(sog + 1105);
    const auto *sog_1106 = buffer.data(sog + 1106);
    const auto *sog_1107 = buffer.data(sog + 1107);
    const auto *sog_1108 = buffer.data(sog + 1108);
    const auto *sog_1109 = buffer.data(sog + 1109);
    const auto *sog_1110 = buffer.data(sog + 1110);
    const auto *sog_1112 = buffer.data(sog + 1112);
    const auto *sog_1113 = buffer.data(sog + 1113);
    const auto *sog_1115 = buffer.data(sog + 1115);
    const auto *sog_1116 = buffer.data(sog + 1116);
    const auto *sog_1119 = buffer.data(sog + 1119);
    const auto *sog_1120 = buffer.data(sog + 1120);
    const auto *sog_1121 = buffer.data(sog + 1121);
    const auto *sog_1122 = buffer.data(sog + 1122);
    const auto *sog_1123 = buffer.data(sog + 1123);
    const auto *sog_1124 = buffer.data(sog + 1124);
    const auto *sog_1125 = buffer.data(sog + 1125);
    const auto *sog_1127 = buffer.data(sog + 1127);
    const auto *sog_1128 = buffer.data(sog + 1128);
    const auto *sog_1130 = buffer.data(sog + 1130);
    const auto *sog_1131 = buffer.data(sog + 1131);
    const auto *sog_1134 = buffer.data(sog + 1134);
    const auto *sog_1135 = buffer.data(sog + 1135);
    const auto *sog_1136 = buffer.data(sog + 1136);
    const auto *sog_1137 = buffer.data(sog + 1137);
    const auto *sog_1138 = buffer.data(sog + 1138);
    const auto *sog_1139 = buffer.data(sog + 1139);
    const auto *sog_1140 = buffer.data(sog + 1140);
    const auto *sog_1142 = buffer.data(sog + 1142);
    const auto *sog_1143 = buffer.data(sog + 1143);
    const auto *sog_1145 = buffer.data(sog + 1145);
    const auto *sog_1146 = buffer.data(sog + 1146);
    const auto *sog_1150 = buffer.data(sog + 1150);
    const auto *sog_1151 = buffer.data(sog + 1151);
    const auto *sog_1152 = buffer.data(sog + 1152);
    const auto *sog_1153 = buffer.data(sog + 1153);
    const auto *sog_1154 = buffer.data(sog + 1154);
    const auto *sog_1155 = buffer.data(sog + 1155);
    const auto *sog_1157 = buffer.data(sog + 1157);
    const auto *sog_1158 = buffer.data(sog + 1158);
    const auto *sog_1160 = buffer.data(sog + 1160);
    const auto *sog_1161 = buffer.data(sog + 1161);
    const auto *sog_1164 = buffer.data(sog + 1164);
    const auto *sog_1165 = buffer.data(sog + 1165);
    const auto *sog_1166 = buffer.data(sog + 1166);
    const auto *sog_1167 = buffer.data(sog + 1167);
    const auto *sog_1168 = buffer.data(sog + 1168);
    const auto *sog_1169 = buffer.data(sog + 1169);

#pragma omp simd aligned(t_1514, t_1515, t_1516, pc_x, pc_y, pc_z, sng_900, sng_917, sof0_723, \
                         sof1_723, sog_1080, sog_1082, sog_1083 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1514[k] = f_17 * sng_900[k]
                    + f_3 * pc_z[k] * sog_1080[k];

        t_1515[k] = f_4 * sof0_723[k]
                    - f_5 * sof1_723[k]
                    + f_3 * pc_x[k] * sog_1083[k];

        t_1516[k] = f_18 * sng_917[k]
                    + f_3 * pc_y[k] * sog_1082[k];
    }

#pragma omp simd aligned(t_1517, t_1518, t_1519, t_1520, pc_x, pc_y, pc_z, sng_903, sng_920, \
                         sof0_725, sof0_726, sof1_725, sof1_726, sog_1083, sog_1085, \
                         sog_1086 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1517[k] = f_4 * sof0_725[k]
                    - f_5 * sof1_725[k]
                    + f_3 * pc_x[k] * sog_1085[k];

        t_1518[k] = f_6 * sof0_726[k]
                    - f_7 * sof1_726[k]
                    + f_3 * pc_x[k] * sog_1086[k];

        t_1519[k] = f_17 * sng_903[k]
                    + f_3 * pc_z[k] * sog_1083[k];

        t_1520[k] = f_18 * sng_920[k]
                    + f_3 * pc_y[k] * sog_1085[k];
    }

#pragma omp simd aligned(t_1521, t_1522, t_1523, t_1524, t_1525, t_1526, pc_x, sof0_729, \
                         sof1_729, sog_1089, sog_1090, sog_1091, sog_1092, sog_1093, \
                         sog_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1521[k] = f_6 * sof0_729[k]
                    - f_7 * sof1_729[k]
                    + f_3 * pc_x[k] * sog_1089[k];

        t_1522[k] = f_3 * pc_x[k] * sog_1090[k];

        t_1523[k] = f_3 * pc_x[k] * sog_1091[k];

        t_1524[k] = f_3 * pc_x[k] * sog_1092[k];

        t_1525[k] = f_3 * pc_x[k] * sog_1093[k];

        t_1526[k] = f_3 * pc_x[k] * sog_1094[k];
    }

#pragma omp simd aligned(t_1527, t_1528, t_1529, pc_y, pc_z, sng_910, sng_925, sng_927, \
                         sof0_726, sof0_728, sof1_726, sof1_728, sog_1090, \
                         sog_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1527[k] = f_18 * sng_925[k]
                    + f_1 * sof0_726[k]
                    - f_2 * sof1_726[k]
                    + f_3 * pc_y[k] * sog_1090[k];

        t_1528[k] = f_17 * sng_910[k]
                    + f_3 * pc_z[k] * sog_1090[k];

        t_1529[k] = f_18 * sng_927[k]
                    + f_4 * sof0_728[k]
                    - f_5 * sof1_728[k]
                    + f_3 * pc_y[k] * sog_1092[k];
    }

#pragma omp simd aligned(t_1530, t_1531, t_1532, pc_y, pc_z, sng_914, sng_928, sng_929, \
                         sof0_729, sof1_729, sog_1093, sog_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1530[k] = f_18 * sng_928[k]
                    + f_6 * sof0_729[k]
                    - f_7 * sof1_729[k]
                    + f_3 * pc_y[k] * sog_1093[k];

        t_1531[k] = f_18 * sng_929[k]
                    + f_3 * pc_y[k] * sog_1094[k];

        t_1532[k] = f_17 * sng_914[k]
                    + f_1 * sof0_729[k]
                    - f_2 * sof1_729[k]
                    + f_3 * pc_z[k] * sog_1094[k];
    }

#pragma omp simd aligned(t_1533, t_1534, t_1535, t_1536, pc_x, pc_y, pc_z, sng_915, sng_930, \
                         sof0_730, sof0_733, sof1_730, sof1_733, sog_1095, \
                         sog_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1533[k] = f_1 * sof0_730[k]
                    - f_2 * sof1_730[k]
                    + f_3 * pc_x[k] * sog_1095[k];

        t_1534[k] = f_16 * sng_930[k]
                    + f_3 * pc_y[k] * sog_1095[k];

        t_1535[k] = f_15 * sng_915[k]
                    + f_3 * pc_z[k] * sog_1095[k];

        t_1536[k] = f_4 * sof0_733[k]
                    - f_5 * sof1_733[k]
                    + f_3 * pc_x[k] * sog_1098[k];
    }

#pragma omp simd aligned(t_1537, t_1538, t_1539, pc_x, pc_y, sng_932, sof0_735, sof0_736, \
                         sof1_735, sof1_736, sog_1097, sog_1100, \
                         sog_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1537[k] = f_16 * sng_932[k]
                    + f_3 * pc_y[k] * sog_1097[k];

        t_1538[k] = f_4 * sof0_735[k]
                    - f_5 * sof1_735[k]
                    + f_3 * pc_x[k] * sog_1100[k];

        t_1539[k] = f_6 * sof0_736[k]
                    - f_7 * sof1_736[k]
                    + f_3 * pc_x[k] * sog_1101[k];
    }

#pragma omp simd aligned(t_1540, t_1541, t_1542, t_1543, pc_x, pc_y, pc_z, sng_918, sng_935, \
                         sof0_739, sof1_739, sog_1098, sog_1100, sog_1104, \
                         sog_1105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1540[k] = f_15 * sng_918[k]
                    + f_3 * pc_z[k] * sog_1098[k];

        t_1541[k] = f_16 * sng_935[k]
                    + f_3 * pc_y[k] * sog_1100[k];

        t_1542[k] = f_6 * sof0_739[k]
                    - f_7 * sof1_739[k]
                    + f_3 * pc_x[k] * sog_1104[k];

        t_1543[k] = f_3 * pc_x[k] * sog_1105[k];
    }

#pragma omp simd aligned(t_1544, t_1545, t_1546, t_1547, t_1548, pc_x, pc_y, sng_940, \
                         sof0_736, sof1_736, sog_1105, sog_1106, sog_1107, sog_1108, \
                         sog_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1544[k] = f_3 * pc_x[k] * sog_1106[k];

        t_1545[k] = f_3 * pc_x[k] * sog_1107[k];

        t_1546[k] = f_3 * pc_x[k] * sog_1108[k];

        t_1547[k] = f_3 * pc_x[k] * sog_1109[k];

        t_1548[k] = f_16 * sng_940[k]
                    + f_1 * sof0_736[k]
                    - f_2 * sof1_736[k]
                    + f_3 * pc_y[k] * sog_1105[k];
    }

#pragma omp simd aligned(t_1549, t_1550, t_1551, pc_y, pc_z, sng_925, sng_942, sng_943, \
                         sof0_738, sof0_739, sof1_738, sof1_739, sog_1105, sog_1107, \
                         sog_1108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1549[k] = f_15 * sng_925[k]
                    + f_3 * pc_z[k] * sog_1105[k];

        t_1550[k] = f_16 * sng_942[k]
                    + f_4 * sof0_738[k]
                    - f_5 * sof1_738[k]
                    + f_3 * pc_y[k] * sog_1107[k];

        t_1551[k] = f_16 * sng_943[k]
                    + f_6 * sof0_739[k]
                    - f_7 * sof1_739[k]
                    + f_3 * pc_y[k] * sog_1108[k];
    }

#pragma omp simd aligned(t_1552, t_1553, t_1554, t_1555, pc_x, pc_y, pc_z, sng_929, sng_944, \
                         sng_945, sof0_739, sof0_740, sof1_739, sof1_740, sog_1109, \
                         sog_1110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1552[k] = f_16 * sng_944[k]
                    + f_3 * pc_y[k] * sog_1109[k];

        t_1553[k] = f_15 * sng_929[k]
                    + f_1 * sof0_739[k]
                    - f_2 * sof1_739[k]
                    + f_3 * pc_z[k] * sog_1109[k];

        t_1554[k] = f_1 * sof0_740[k]
                    - f_2 * sof1_740[k]
                    + f_3 * pc_x[k] * sog_1110[k];

        t_1555[k] = f_11 * sng_945[k]
                    + f_3 * pc_y[k] * sog_1110[k];
    }

#pragma omp simd aligned(t_1556, t_1557, t_1558, pc_x, pc_y, pc_z, sng_930, sng_947, sof0_743, \
                         sof1_743, sog_1110, sog_1112, sog_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1556[k] = f_14 * sng_930[k]
                    + f_3 * pc_z[k] * sog_1110[k];

        t_1557[k] = f_4 * sof0_743[k]
                    - f_5 * sof1_743[k]
                    + f_3 * pc_x[k] * sog_1113[k];

        t_1558[k] = f_11 * sng_947[k]
                    + f_3 * pc_y[k] * sog_1112[k];
    }

#pragma omp simd aligned(t_1559, t_1560, t_1561, t_1562, pc_x, pc_y, pc_z, sng_933, sng_950, \
                         sof0_745, sof0_746, sof1_745, sof1_746, sog_1113, sog_1115, \
                         sog_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1559[k] = f_4 * sof0_745[k]
                    - f_5 * sof1_745[k]
                    + f_3 * pc_x[k] * sog_1115[k];

        t_1560[k] = f_6 * sof0_746[k]
                    - f_7 * sof1_746[k]
                    + f_3 * pc_x[k] * sog_1116[k];

        t_1561[k] = f_14 * sng_933[k]
                    + f_3 * pc_z[k] * sog_1113[k];

        t_1562[k] = f_11 * sng_950[k]
                    + f_3 * pc_y[k] * sog_1115[k];
    }

#pragma omp simd aligned(t_1563, t_1564, t_1565, t_1566, t_1567, t_1568, pc_x, sof0_749, \
                         sof1_749, sog_1119, sog_1120, sog_1121, sog_1122, sog_1123, \
                         sog_1124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1563[k] = f_6 * sof0_749[k]
                    - f_7 * sof1_749[k]
                    + f_3 * pc_x[k] * sog_1119[k];

        t_1564[k] = f_3 * pc_x[k] * sog_1120[k];

        t_1565[k] = f_3 * pc_x[k] * sog_1121[k];

        t_1566[k] = f_3 * pc_x[k] * sog_1122[k];

        t_1567[k] = f_3 * pc_x[k] * sog_1123[k];

        t_1568[k] = f_3 * pc_x[k] * sog_1124[k];
    }

#pragma omp simd aligned(t_1569, t_1570, t_1571, pc_y, pc_z, sng_940, sng_955, sng_957, \
                         sof0_746, sof0_748, sof1_746, sof1_748, sog_1120, \
                         sog_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1569[k] = f_11 * sng_955[k]
                    + f_1 * sof0_746[k]
                    - f_2 * sof1_746[k]
                    + f_3 * pc_y[k] * sog_1120[k];

        t_1570[k] = f_14 * sng_940[k]
                    + f_3 * pc_z[k] * sog_1120[k];

        t_1571[k] = f_11 * sng_957[k]
                    + f_4 * sof0_748[k]
                    - f_5 * sof1_748[k]
                    + f_3 * pc_y[k] * sog_1122[k];
    }

#pragma omp simd aligned(t_1572, t_1573, t_1574, pc_y, pc_z, sng_944, sng_958, sng_959, \
                         sof0_749, sof1_749, sog_1123, sog_1124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1572[k] = f_11 * sng_958[k]
                    + f_6 * sof0_749[k]
                    - f_7 * sof1_749[k]
                    + f_3 * pc_y[k] * sog_1123[k];

        t_1573[k] = f_11 * sng_959[k]
                    + f_3 * pc_y[k] * sog_1124[k];

        t_1574[k] = f_14 * sng_944[k]
                    + f_1 * sof0_749[k]
                    - f_2 * sof1_749[k]
                    + f_3 * pc_z[k] * sog_1124[k];
    }

#pragma omp simd aligned(t_1575, t_1576, t_1577, t_1578, pc_x, pc_y, pc_z, sng_945, sng_960, \
                         sof0_750, sof0_753, sof1_750, sof1_753, sog_1125, \
                         sog_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1575[k] = f_1 * sof0_750[k]
                    - f_2 * sof1_750[k]
                    + f_3 * pc_x[k] * sog_1125[k];

        t_1576[k] = f_10 * sng_960[k]
                    + f_3 * pc_y[k] * sog_1125[k];

        t_1577[k] = f_13 * sng_945[k]
                    + f_3 * pc_z[k] * sog_1125[k];

        t_1578[k] = f_4 * sof0_753[k]
                    - f_5 * sof1_753[k]
                    + f_3 * pc_x[k] * sog_1128[k];
    }

#pragma omp simd aligned(t_1579, t_1580, t_1581, pc_x, pc_y, sng_962, sof0_755, sof0_756, \
                         sof1_755, sof1_756, sog_1127, sog_1130, \
                         sog_1131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1579[k] = f_10 * sng_962[k]
                    + f_3 * pc_y[k] * sog_1127[k];

        t_1580[k] = f_4 * sof0_755[k]
                    - f_5 * sof1_755[k]
                    + f_3 * pc_x[k] * sog_1130[k];

        t_1581[k] = f_6 * sof0_756[k]
                    - f_7 * sof1_756[k]
                    + f_3 * pc_x[k] * sog_1131[k];
    }

#pragma omp simd aligned(t_1582, t_1583, t_1584, t_1585, pc_x, pc_y, pc_z, sng_948, sng_965, \
                         sof0_759, sof1_759, sog_1128, sog_1130, sog_1134, \
                         sog_1135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1582[k] = f_13 * sng_948[k]
                    + f_3 * pc_z[k] * sog_1128[k];

        t_1583[k] = f_10 * sng_965[k]
                    + f_3 * pc_y[k] * sog_1130[k];

        t_1584[k] = f_6 * sof0_759[k]
                    - f_7 * sof1_759[k]
                    + f_3 * pc_x[k] * sog_1134[k];

        t_1585[k] = f_3 * pc_x[k] * sog_1135[k];
    }

#pragma omp simd aligned(t_1586, t_1587, t_1588, t_1589, t_1590, pc_x, pc_y, sng_970, \
                         sof0_756, sof1_756, sog_1135, sog_1136, sog_1137, sog_1138, \
                         sog_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1586[k] = f_3 * pc_x[k] * sog_1136[k];

        t_1587[k] = f_3 * pc_x[k] * sog_1137[k];

        t_1588[k] = f_3 * pc_x[k] * sog_1138[k];

        t_1589[k] = f_3 * pc_x[k] * sog_1139[k];

        t_1590[k] = f_10 * sng_970[k]
                    + f_1 * sof0_756[k]
                    - f_2 * sof1_756[k]
                    + f_3 * pc_y[k] * sog_1135[k];
    }

#pragma omp simd aligned(t_1591, t_1592, t_1593, pc_y, pc_z, sng_955, sng_972, sng_973, \
                         sof0_758, sof0_759, sof1_758, sof1_759, sog_1135, sog_1137, \
                         sog_1138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1591[k] = f_13 * sng_955[k]
                    + f_3 * pc_z[k] * sog_1135[k];

        t_1592[k] = f_10 * sng_972[k]
                    + f_4 * sof0_758[k]
                    - f_5 * sof1_758[k]
                    + f_3 * pc_y[k] * sog_1137[k];

        t_1593[k] = f_10 * sng_973[k]
                    + f_6 * sof0_759[k]
                    - f_7 * sof1_759[k]
                    + f_3 * pc_y[k] * sog_1138[k];
    }

#pragma omp simd aligned(t_1594, t_1595, t_1596, t_1597, pb_y, pc_y, pc_z, snh0_1365, sng_959, \
                         sng_974, sng_975, snh1_1365, sof0_759, sof1_759, sog_1139, \
                         sog_1140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1594[k] = f_10 * sng_974[k]
                    + f_3 * pc_y[k] * sog_1139[k];

        t_1595[k] = f_13 * sng_959[k]
                    + f_1 * sof0_759[k]
                    - f_2 * sof1_759[k]
                    + f_3 * pc_z[k] * sog_1139[k];

        t_1596[k] = pb_y[k] * snh0_1365[k]
                    - f_8 * pc_y[k] * snh1_1365[k];

        t_1597[k] = f_9 * sng_975[k]
                    + f_3 * pc_y[k] * sog_1140[k];
    }

#pragma omp simd aligned(t_1598, t_1599, t_1600, pc_x, pc_y, pc_z, sng_960, sng_977, sof0_763, \
                         sof1_763, sog_1140, sog_1142, sog_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1598[k] = f_12 * sng_960[k]
                    + f_3 * pc_z[k] * sog_1140[k];

        t_1599[k] = f_4 * sof0_763[k]
                    - f_5 * sof1_763[k]
                    + f_3 * pc_x[k] * sog_1143[k];

        t_1600[k] = f_9 * sng_977[k]
                    + f_3 * pc_y[k] * sog_1142[k];
    }

#pragma omp simd aligned(t_1601, t_1602, t_1603, pb_y, pc_x, pc_y, pc_z, snh0_1370, sng_963, \
                         snh1_1370, sof0_766, sof1_766, sog_1143, \
                         sog_1146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1601[k] = pb_y[k] * snh0_1370[k]
                    - f_8 * pc_y[k] * snh1_1370[k];

        t_1602[k] = f_6 * sof0_766[k]
                    - f_7 * sof1_766[k]
                    + f_3 * pc_x[k] * sog_1146[k];

        t_1603[k] = f_12 * sng_963[k]
                    + f_3 * pc_z[k] * sog_1143[k];
    }

#pragma omp simd aligned(t_1604, t_1605, t_1606, t_1607, t_1608, pb_y, pc_x, pc_y, snh0_1374, \
                         sng_980, snh1_1374, sog_1145, sog_1150, sog_1151, \
                         sog_1152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1604[k] = f_9 * sng_980[k]
                    + f_3 * pc_y[k] * sog_1145[k];

        t_1605[k] = pb_y[k] * snh0_1374[k]
                    - f_8 * pc_y[k] * snh1_1374[k];

        t_1606[k] = f_3 * pc_x[k] * sog_1150[k];

        t_1607[k] = f_3 * pc_x[k] * sog_1151[k];

        t_1608[k] = f_3 * pc_x[k] * sog_1152[k];
    }

#pragma omp simd aligned(t_1609, t_1610, t_1611, t_1612, pb_y, pc_x, pc_y, pc_z, snh0_1380, \
                         sng_970, sng_985, snh1_1380, sog_1150, sog_1153, \
                         sog_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1609[k] = f_3 * pc_x[k] * sog_1153[k];

        t_1610[k] = f_3 * pc_x[k] * sog_1154[k];

        t_1611[k] = pb_y[k] * snh0_1380[k]
                    + f_18 * sng_985[k]
                    - f_8 * pc_y[k] * snh1_1380[k];

        t_1612[k] = f_12 * sng_970[k]
                    + f_3 * pc_z[k] * sog_1150[k];
    }

#pragma omp simd aligned(t_1613, t_1614, t_1615, t_1616, pb_y, pc_y, snh0_1382, snh0_1383, \
                         snh0_1385, sng_987, sng_988, sng_989, snh1_1382, snh1_1383, \
                         snh1_1385, sog_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1613[k] = pb_y[k] * snh0_1382[k]
                    + f_11 * sng_987[k]
                    - f_8 * pc_y[k] * snh1_1382[k];

        t_1614[k] = pb_y[k] * snh0_1383[k]
                    + f_10 * sng_988[k]
                    - f_8 * pc_y[k] * snh1_1383[k];

        t_1615[k] = f_9 * sng_989[k]
                    + f_3 * pc_y[k] * sog_1154[k];

        t_1616[k] = pb_y[k] * snh0_1385[k]
                    - f_8 * pc_y[k] * snh1_1385[k];
    }

#pragma omp simd aligned(t_1617, t_1618, t_1619, t_1620, t_1621, pc_x, pc_y, pc_z, sng_975, \
                         sof0_770, sof0_773, sof1_770, sof1_773, sog_1155, sog_1157, \
                         sog_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1617[k] = f_1 * sof0_770[k]
                    - f_2 * sof1_770[k]
                    + f_3 * pc_x[k] * sog_1155[k];

        t_1618[k] = f_3 * pc_y[k] * sog_1155[k];

        t_1619[k] = f_0 * sng_975[k]
                    + f_3 * pc_z[k] * sog_1155[k];

        t_1620[k] = f_4 * sof0_773[k]
                    - f_5 * sof1_773[k]
                    + f_3 * pc_x[k] * sog_1158[k];

        t_1621[k] = f_3 * pc_y[k] * sog_1157[k];
    }

#pragma omp simd aligned(t_1622, t_1623, t_1624, t_1625, pc_x, pc_y, pc_z, sng_978, sof0_775, \
                         sof0_776, sof1_775, sof1_776, sog_1158, sog_1160, \
                         sog_1161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1622[k] = f_4 * sof0_775[k]
                    - f_5 * sof1_775[k]
                    + f_3 * pc_x[k] * sog_1160[k];

        t_1623[k] = f_6 * sof0_776[k]
                    - f_7 * sof1_776[k]
                    + f_3 * pc_x[k] * sog_1161[k];

        t_1624[k] = f_0 * sng_978[k]
                    + f_3 * pc_z[k] * sog_1158[k];

        t_1625[k] = f_3 * pc_y[k] * sog_1160[k];
    }

#pragma omp simd aligned(t_1626, t_1627, t_1628, t_1629, t_1630, t_1631, pc_x, sof0_779, \
                         sof1_779, sog_1164, sog_1165, sog_1166, sog_1167, sog_1168, \
                         sog_1169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1626[k] = f_6 * sof0_779[k]
                    - f_7 * sof1_779[k]
                    + f_3 * pc_x[k] * sog_1164[k];

        t_1627[k] = f_3 * pc_x[k] * sog_1165[k];

        t_1628[k] = f_3 * pc_x[k] * sog_1166[k];

        t_1629[k] = f_3 * pc_x[k] * sog_1167[k];

        t_1630[k] = f_3 * pc_x[k] * sog_1168[k];

        t_1631[k] = f_3 * pc_x[k] * sog_1169[k];
    }

#pragma omp simd aligned(t_1632, t_1633, t_1634, t_1635, pc_y, pc_z, sng_985, sof0_776, \
                         sof0_778, sof0_779, sof1_776, sof1_778, sof1_779, sog_1165, sog_1167, \
                         sog_1168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1632[k] = f_1 * sof0_776[k]
                    - f_2 * sof1_776[k]
                    + f_3 * pc_y[k] * sog_1165[k];

        t_1633[k] = f_0 * sng_985[k]
                    + f_3 * pc_z[k] * sog_1165[k];

        t_1634[k] = f_4 * sof0_778[k]
                    - f_5 * sof1_778[k]
                    + f_3 * pc_y[k] * sog_1167[k];

        t_1635[k] = f_6 * sof0_779[k]
                    - f_7 * sof1_779[k]
                    + f_3 * pc_y[k] * sog_1168[k];
    }

#pragma omp simd aligned(t_1636, t_1637, pc_y, pc_z, sng_989, sof0_779, sof1_779, \
                         sog_1169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1636[k] = f_3 * pc_y[k] * sog_1169[k];

        t_1637[k] = f_0 * sng_989[k]
                    + f_1 * sof0_779[k]
                    - f_2 * sof1_779[k]
                    + f_3 * pc_z[k] * sog_1169[k];
    }
}

auto
compute_prim_soh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t snh0, const size_t sng,
                                                   const size_t snh1, const size_t sof0,
                                                   const size_t sof1, const size_t sog,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_soh_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, snh0, sng,
                                                              snh1, sof0, sof1, sog, ncols,
                                                              gamma, p, q);

    compute_prim_soh_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, snh0, sng,
                                                              snh1, sof0, sof1, sog, ncols,
                                                              gamma, p, q);

    compute_prim_soh_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, snh0, sng,
                                                              snh1, sof0, sof1, sog, ncols,
                                                              gamma, p, q);

    compute_prim_soh_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, snh0, sng,
                                                              snh1, sof0, sof1, sog, ncols,
                                                              gamma, p, q);

    compute_prim_soh_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, snh0, sng,
                                                              snh1, sof0, sof1, sog, ncols,
                                                              gamma, p, q);

    compute_prim_soh_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, snh0, sng,
                                                              snh1, sof0, sof1, sog, ncols,
                                                              gamma, p, q);

    compute_prim_soh_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, snh0, sng,
                                                              snh1, sof0, sof1, sog, ncols,
                                                              gamma, p, q);

    compute_prim_soh_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, snh0, sng,
                                                              snh1, sof0, sof1, sog, ncols,
                                                              gamma, p, q);

    compute_prim_soh_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, snh0, sng,
                                                              snh1, sof0, sof1, sog, ncols,
                                                              gamma, p, q);

    compute_prim_soh_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, snh0, sng,
                                                              snh1, sof0, sof1, sog, ncols,
                                                              gamma, p, q);

    compute_prim_soh_three_center_electron_repulsion_0_piece10(buffer, target, pb, pc, snh0,
                                                               sng, snh1, sof0, sof1, sog,
                                                               ncols, gamma, p, q);

    compute_prim_soh_three_center_electron_repulsion_0_piece11(buffer, target, pb, pc, snh0,
                                                               sng, snh1, sof0, sof1, sog,
                                                               ncols, gamma, p, q);

    compute_prim_soh_three_center_electron_repulsion_0_piece12(buffer, target, pb, pc, snh0,
                                                               sng, snh1, sof0, sof1, sog,
                                                               ncols, gamma, p, q);

    compute_prim_soh_three_center_electron_repulsion_0_piece13(buffer, target, pb, pc, snh0,
                                                               sng, snh1, sof0, sof1, sog,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
