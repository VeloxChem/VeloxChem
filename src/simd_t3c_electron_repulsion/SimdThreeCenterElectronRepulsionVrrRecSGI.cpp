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


#include "SimdThreeCenterElectronRepulsionVrrRecSGI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sgi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfi0,
                                                          const size_t sfh, const size_t sfi1,
                                                          const size_t sgg0, const size_t sgg1,
                                                          const size_t sgh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_12 = 1.0 / q;
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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfi0_0 = buffer.data(sfi0 + 0);
    const auto *sfi0_3 = buffer.data(sfi0 + 3);
    const auto *sfi0_5 = buffer.data(sfi0 + 5);
    const auto *sfi0_6 = buffer.data(sfi0 + 6);
    const auto *sfi0_9 = buffer.data(sfi0 + 9);
    const auto *sfi0_10 = buffer.data(sfi0 + 10);
    const auto *sfi0_12 = buffer.data(sfi0 + 12);
    const auto *sfi0_14 = buffer.data(sfi0 + 14);
    const auto *sfi0_21 = buffer.data(sfi0 + 21);
    const auto *sfi0_27 = buffer.data(sfi0 + 27);
    const auto *sfi0_31 = buffer.data(sfi0 + 31);
    const auto *sfi0_34 = buffer.data(sfi0 + 34);
    const auto *sfi0_38 = buffer.data(sfi0 + 38);
    const auto *sfi0_56 = buffer.data(sfi0 + 56);
    const auto *sfi0_61 = buffer.data(sfi0 + 61);
    const auto *sfi0_65 = buffer.data(sfi0 + 65);

    const auto *sfh_0 = buffer.data(sfh + 0);
    const auto *sfh_1 = buffer.data(sfh + 1);
    const auto *sfh_2 = buffer.data(sfh + 2);
    const auto *sfh_3 = buffer.data(sfh + 3);
    const auto *sfh_5 = buffer.data(sfh + 5);
    const auto *sfh_6 = buffer.data(sfh + 6);
    const auto *sfh_7 = buffer.data(sfh + 7);
    const auto *sfh_8 = buffer.data(sfh + 8);
    const auto *sfh_9 = buffer.data(sfh + 9);
    const auto *sfh_10 = buffer.data(sfh + 10);
    const auto *sfh_12 = buffer.data(sfh + 12);
    const auto *sfh_14 = buffer.data(sfh + 14);
    const auto *sfh_15 = buffer.data(sfh + 15);
    const auto *sfh_16 = buffer.data(sfh + 16);
    const auto *sfh_17 = buffer.data(sfh + 17);
    const auto *sfh_18 = buffer.data(sfh + 18);
    const auto *sfh_19 = buffer.data(sfh + 19);
    const auto *sfh_20 = buffer.data(sfh + 20);
    const auto *sfh_21 = buffer.data(sfh + 21);
    const auto *sfh_23 = buffer.data(sfh + 23);
    const auto *sfh_24 = buffer.data(sfh + 24);
    const auto *sfh_26 = buffer.data(sfh + 26);
    const auto *sfh_27 = buffer.data(sfh + 27);
    const auto *sfh_30 = buffer.data(sfh + 30);
    const auto *sfh_36 = buffer.data(sfh + 36);
    const auto *sfh_37 = buffer.data(sfh + 37);
    const auto *sfh_38 = buffer.data(sfh + 38);
    const auto *sfh_39 = buffer.data(sfh + 39);
    const auto *sfh_40 = buffer.data(sfh + 40);
    const auto *sfh_41 = buffer.data(sfh + 41);
    const auto *sfh_42 = buffer.data(sfh + 42);
    const auto *sfh_44 = buffer.data(sfh + 44);
    const auto *sfh_47 = buffer.data(sfh + 47);
    const auto *sfh_57 = buffer.data(sfh + 57);
    const auto *sfh_58 = buffer.data(sfh + 58);
    const auto *sfh_59 = buffer.data(sfh + 59);
    const auto *sfh_60 = buffer.data(sfh + 60);
    const auto *sfh_61 = buffer.data(sfh + 61);
    const auto *sfh_62 = buffer.data(sfh + 62);
    const auto *sfh_63 = buffer.data(sfh + 63);
    const auto *sfh_66 = buffer.data(sfh + 66);
    const auto *sfh_68 = buffer.data(sfh + 68);
    const auto *sfh_69 = buffer.data(sfh + 69);
    const auto *sfh_72 = buffer.data(sfh + 72);
    const auto *sfh_73 = buffer.data(sfh + 73);
    const auto *sfh_75 = buffer.data(sfh + 75);
    const auto *sfh_77 = buffer.data(sfh + 77);
    const auto *sfh_78 = buffer.data(sfh + 78);
    const auto *sfh_79 = buffer.data(sfh + 79);
    const auto *sfh_80 = buffer.data(sfh + 80);
    const auto *sfh_81 = buffer.data(sfh + 81);
    const auto *sfh_82 = buffer.data(sfh + 82);
    const auto *sfh_83 = buffer.data(sfh + 83);

    const auto *sfi1_0 = buffer.data(sfi1 + 0);
    const auto *sfi1_3 = buffer.data(sfi1 + 3);
    const auto *sfi1_5 = buffer.data(sfi1 + 5);
    const auto *sfi1_6 = buffer.data(sfi1 + 6);
    const auto *sfi1_9 = buffer.data(sfi1 + 9);
    const auto *sfi1_10 = buffer.data(sfi1 + 10);
    const auto *sfi1_12 = buffer.data(sfi1 + 12);
    const auto *sfi1_14 = buffer.data(sfi1 + 14);
    const auto *sfi1_21 = buffer.data(sfi1 + 21);
    const auto *sfi1_27 = buffer.data(sfi1 + 27);
    const auto *sfi1_31 = buffer.data(sfi1 + 31);
    const auto *sfi1_34 = buffer.data(sfi1 + 34);
    const auto *sfi1_38 = buffer.data(sfi1 + 38);
    const auto *sfi1_56 = buffer.data(sfi1 + 56);
    const auto *sfi1_61 = buffer.data(sfi1 + 61);
    const auto *sfi1_65 = buffer.data(sfi1 + 65);

    const auto *sgg0_0 = buffer.data(sgg0 + 0);
    const auto *sgg0_3 = buffer.data(sgg0 + 3);
    const auto *sgg0_5 = buffer.data(sgg0 + 5);
    const auto *sgg0_6 = buffer.data(sgg0 + 6);
    const auto *sgg0_9 = buffer.data(sgg0 + 9);
    const auto *sgg0_10 = buffer.data(sgg0 + 10);
    const auto *sgg0_12 = buffer.data(sgg0 + 12);
    const auto *sgg0_13 = buffer.data(sgg0 + 13);
    const auto *sgg0_14 = buffer.data(sgg0 + 14);
    const auto *sgg0_25 = buffer.data(sgg0 + 25);
    const auto *sgg0_27 = buffer.data(sgg0 + 27);
    const auto *sgg0_28 = buffer.data(sgg0 + 28);
    const auto *sgg0_29 = buffer.data(sgg0 + 29);
    const auto *sgg0_42 = buffer.data(sgg0 + 42);
    const auto *sgg0_43 = buffer.data(sgg0 + 43);
    const auto *sgg0_44 = buffer.data(sgg0 + 44);
    const auto *sgg0_45 = buffer.data(sgg0 + 45);
    const auto *sgg0_48 = buffer.data(sgg0 + 48);
    const auto *sgg0_50 = buffer.data(sgg0 + 50);
    const auto *sgg0_51 = buffer.data(sgg0 + 51);
    const auto *sgg0_54 = buffer.data(sgg0 + 54);
    const auto *sgg0_55 = buffer.data(sgg0 + 55);
    const auto *sgg0_57 = buffer.data(sgg0 + 57);
    const auto *sgg0_58 = buffer.data(sgg0 + 58);
    const auto *sgg0_59 = buffer.data(sgg0 + 59);

    const auto *sgg1_0 = buffer.data(sgg1 + 0);
    const auto *sgg1_3 = buffer.data(sgg1 + 3);
    const auto *sgg1_5 = buffer.data(sgg1 + 5);
    const auto *sgg1_6 = buffer.data(sgg1 + 6);
    const auto *sgg1_9 = buffer.data(sgg1 + 9);
    const auto *sgg1_10 = buffer.data(sgg1 + 10);
    const auto *sgg1_12 = buffer.data(sgg1 + 12);
    const auto *sgg1_13 = buffer.data(sgg1 + 13);
    const auto *sgg1_14 = buffer.data(sgg1 + 14);
    const auto *sgg1_25 = buffer.data(sgg1 + 25);
    const auto *sgg1_27 = buffer.data(sgg1 + 27);
    const auto *sgg1_28 = buffer.data(sgg1 + 28);
    const auto *sgg1_29 = buffer.data(sgg1 + 29);
    const auto *sgg1_42 = buffer.data(sgg1 + 42);
    const auto *sgg1_43 = buffer.data(sgg1 + 43);
    const auto *sgg1_44 = buffer.data(sgg1 + 44);
    const auto *sgg1_45 = buffer.data(sgg1 + 45);
    const auto *sgg1_48 = buffer.data(sgg1 + 48);
    const auto *sgg1_50 = buffer.data(sgg1 + 50);
    const auto *sgg1_51 = buffer.data(sgg1 + 51);
    const auto *sgg1_54 = buffer.data(sgg1 + 54);
    const auto *sgg1_55 = buffer.data(sgg1 + 55);
    const auto *sgg1_57 = buffer.data(sgg1 + 57);
    const auto *sgg1_58 = buffer.data(sgg1 + 58);
    const auto *sgg1_59 = buffer.data(sgg1 + 59);

    const auto *sgh_0 = buffer.data(sgh + 0);
    const auto *sgh_2 = buffer.data(sgh + 2);
    const auto *sgh_3 = buffer.data(sgh + 3);
    const auto *sgh_5 = buffer.data(sgh + 5);
    const auto *sgh_6 = buffer.data(sgh + 6);
    const auto *sgh_9 = buffer.data(sgh + 9);
    const auto *sgh_10 = buffer.data(sgh + 10);
    const auto *sgh_12 = buffer.data(sgh + 12);
    const auto *sgh_14 = buffer.data(sgh + 14);
    const auto *sgh_15 = buffer.data(sgh + 15);
    const auto *sgh_16 = buffer.data(sgh + 16);
    const auto *sgh_17 = buffer.data(sgh + 17);
    const auto *sgh_18 = buffer.data(sgh + 18);
    const auto *sgh_19 = buffer.data(sgh + 19);
    const auto *sgh_20 = buffer.data(sgh + 20);
    const auto *sgh_21 = buffer.data(sgh + 21);
    const auto *sgh_23 = buffer.data(sgh + 23);
    const auto *sgh_24 = buffer.data(sgh + 24);
    const auto *sgh_26 = buffer.data(sgh + 26);
    const auto *sgh_27 = buffer.data(sgh + 27);
    const auto *sgh_30 = buffer.data(sgh + 30);
    const auto *sgh_36 = buffer.data(sgh + 36);
    const auto *sgh_37 = buffer.data(sgh + 37);
    const auto *sgh_38 = buffer.data(sgh + 38);
    const auto *sgh_39 = buffer.data(sgh + 39);
    const auto *sgh_40 = buffer.data(sgh + 40);
    const auto *sgh_41 = buffer.data(sgh + 41);
    const auto *sgh_42 = buffer.data(sgh + 42);
    const auto *sgh_44 = buffer.data(sgh + 44);
    const auto *sgh_45 = buffer.data(sgh + 45);
    const auto *sgh_47 = buffer.data(sgh + 47);
    const auto *sgh_48 = buffer.data(sgh + 48);
    const auto *sgh_51 = buffer.data(sgh + 51);
    const auto *sgh_57 = buffer.data(sgh + 57);
    const auto *sgh_58 = buffer.data(sgh + 58);
    const auto *sgh_59 = buffer.data(sgh + 59);
    const auto *sgh_60 = buffer.data(sgh + 60);
    const auto *sgh_61 = buffer.data(sgh + 61);
    const auto *sgh_62 = buffer.data(sgh + 62);
    const auto *sgh_63 = buffer.data(sgh + 63);
    const auto *sgh_65 = buffer.data(sgh + 65);
    const auto *sgh_66 = buffer.data(sgh + 66);
    const auto *sgh_68 = buffer.data(sgh + 68);
    const auto *sgh_69 = buffer.data(sgh + 69);
    const auto *sgh_72 = buffer.data(sgh + 72);
    const auto *sgh_73 = buffer.data(sgh + 73);
    const auto *sgh_75 = buffer.data(sgh + 75);
    const auto *sgh_77 = buffer.data(sgh + 77);
    const auto *sgh_78 = buffer.data(sgh + 78);
    const auto *sgh_79 = buffer.data(sgh + 79);
    const auto *sgh_80 = buffer.data(sgh + 80);
    const auto *sgh_81 = buffer.data(sgh + 81);
    const auto *sgh_82 = buffer.data(sgh + 82);
    const auto *sgh_83 = buffer.data(sgh + 83);
    const auto *sgh_84 = buffer.data(sgh + 84);
    const auto *sgh_86 = buffer.data(sgh + 86);
    const auto *sgh_87 = buffer.data(sgh + 87);
    const auto *sgh_89 = buffer.data(sgh + 89);
    const auto *sgh_90 = buffer.data(sgh + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sfh_0, sfh_3, sgg0_0, sgg0_3, \
                         sgg1_0, sgg1_3, sgh_0, sgh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfh_0[k]
                 + f_1 * sgg0_0[k]
                 - f_2 * sgg1_0[k]
                 + f_3 * pc_x[k] * sgh_0[k];

        t_1[k] = f_3 * pc_y[k] * sgh_0[k];

        t_2[k] = f_3 * pc_z[k] * sgh_0[k];

        t_3[k] = f_0 * sfh_3[k]
                 + f_4 * sgg0_3[k]
                 - f_5 * sgg1_3[k]
                 + f_3 * pc_x[k] * sgh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sfh_5, sfh_6, sgg0_5, sgg0_6, sgg1_5, \
                         sgg1_6, sgh_2, sgh_5, sgh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sgh_2[k];

        t_5[k] = f_0 * sfh_5[k]
                 + f_4 * sgg0_5[k]
                 - f_5 * sgg1_5[k]
                 + f_3 * pc_x[k] * sgh_5[k];

        t_6[k] = f_0 * sfh_6[k]
                 + f_6 * sgg0_6[k]
                 - f_7 * sgg1_6[k]
                 + f_3 * pc_x[k] * sgh_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sfh_9, sgg0_9, sgg1_9, sgh_3, sgh_5, \
                         sgh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sgh_3[k];

        t_8[k] = f_3 * pc_y[k] * sgh_5[k];

        t_9[k] = f_0 * sfh_9[k]
                 + f_6 * sgg0_9[k]
                 - f_7 * sgg1_9[k]
                 + f_3 * pc_x[k] * sgh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sfh_10, sfh_12, sgg0_10, sgg0_12, \
                         sgg1_10, sgg1_12, sgh_6, sgh_10, sgh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sfh_10[k]
                  + f_8 * sgg0_10[k]
                  - f_9 * sgg1_10[k]
                  + f_3 * pc_x[k] * sgh_10[k];

        t_11[k] = f_3 * pc_z[k] * sgh_6[k];

        t_12[k] = f_0 * sfh_12[k]
                  + f_8 * sgg0_12[k]
                  - f_9 * sgg1_12[k]
                  + f_3 * pc_x[k] * sgh_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pc_x, pc_y, sfh_14, sfh_15, sfh_16, sgg0_14, \
                         sgg1_14, sgh_9, sgh_14, sgh_15, sgh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sgh_9[k];

        t_14[k] = f_0 * sfh_14[k]
                  + f_8 * sgg0_14[k]
                  - f_9 * sgg1_14[k]
                  + f_3 * pc_x[k] * sgh_14[k];

        t_15[k] = f_0 * sfh_15[k]
                  + f_3 * pc_x[k] * sgh_15[k];

        t_16[k] = f_0 * sfh_16[k]
                  + f_3 * pc_x[k] * sgh_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_x, sfh_17, sfh_18, sfh_19, sfh_20, sgh_17, \
                         sgh_18, sgh_19, sgh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * sfh_17[k]
                  + f_3 * pc_x[k] * sgh_17[k];

        t_18[k] = f_0 * sfh_18[k]
                  + f_3 * pc_x[k] * sgh_18[k];

        t_19[k] = f_0 * sfh_19[k]
                  + f_3 * pc_x[k] * sgh_19[k];

        t_20[k] = f_0 * sfh_20[k]
                  + f_3 * pc_x[k] * sgh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, sgg0_10, sgg0_12, sgg0_13, \
                         sgg1_10, sgg1_12, sgg1_13, sgh_15, sgh_17, \
                         sgh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * sgg0_10[k]
                  - f_2 * sgg1_10[k]
                  + f_3 * pc_y[k] * sgh_15[k];

        t_22[k] = f_3 * pc_z[k] * sgh_15[k];

        t_23[k] = f_4 * sgg0_12[k]
                  - f_5 * sgg1_12[k]
                  + f_3 * pc_y[k] * sgh_17[k];

        t_24[k] = f_6 * sgg0_13[k]
                  - f_7 * sgg1_13[k]
                  + f_3 * pc_y[k] * sgh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_y, pc_y, pc_z, sfi0_0, sfh_0, \
                         sfi1_0, sgg0_14, sgg1_14, sgh_19, sgh_20, \
                         sgh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * sgg0_14[k]
                  - f_9 * sgg1_14[k]
                  + f_3 * pc_y[k] * sgh_19[k];

        t_26[k] = f_3 * pc_y[k] * sgh_20[k];

        t_27[k] = f_1 * sgg0_14[k]
                  - f_2 * sgg1_14[k]
                  + f_3 * pc_z[k] * sgh_20[k];

        t_28[k] = pb_y[k] * sfi0_0[k]
                  - f_10 * pc_y[k] * sfi1_0[k];

        t_29[k] = f_11 * sfh_0[k]
                  + f_3 * pc_y[k] * sgh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_y, pc_y, pc_z, sfi0_3, sfi0_5, sfh_1, \
                         sfh_2, sfi1_3, sfi1_5, sgh_21, sgh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * sgh_21[k];

        t_31[k] = pb_y[k] * sfi0_3[k]
                  + f_12 * sfh_1[k]
                  - f_10 * pc_y[k] * sfi1_3[k];

        t_32[k] = f_11 * sfh_2[k]
                  + f_3 * pc_y[k] * sgh_23[k];

        t_33[k] = pb_y[k] * sfi0_5[k]
                  - f_10 * pc_y[k] * sfi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_y, pc_y, pc_z, sfi0_6, sfi0_9, sfh_3, \
                         sfh_5, sfi1_6, sfi1_9, sgh_24, sgh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * sfi0_6[k]
                  + f_13 * sfh_3[k]
                  - f_10 * pc_y[k] * sfi1_6[k];

        t_35[k] = f_3 * pc_z[k] * sgh_24[k];

        t_36[k] = f_11 * sfh_5[k]
                  + f_3 * pc_y[k] * sgh_26[k];

        t_37[k] = pb_y[k] * sfi0_9[k]
                  - f_10 * pc_y[k] * sfi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, sfi0_10, sfi0_12, sfh_6, \
                         sfh_8, sfh_9, sfi1_10, sfi1_12, sgh_27, \
                         sgh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * sfi0_10[k]
                  + f_0 * sfh_6[k]
                  - f_10 * pc_y[k] * sfi1_10[k];

        t_39[k] = f_3 * pc_z[k] * sgh_27[k];

        t_40[k] = pb_y[k] * sfi0_12[k]
                  + f_12 * sfh_8[k]
                  - f_10 * pc_y[k] * sfi1_12[k];

        t_41[k] = f_11 * sfh_9[k]
                  + f_3 * pc_y[k] * sgh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_y, pc_x, pc_y, sfi0_14, sfh_36, sfh_37, \
                         sfh_38, sfi1_14, sgh_36, sgh_37, sgh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * sfi0_14[k]
                  - f_10 * pc_y[k] * sfi1_14[k];

        t_43[k] = f_13 * sfh_36[k]
                  + f_3 * pc_x[k] * sgh_36[k];

        t_44[k] = f_13 * sfh_37[k]
                  + f_3 * pc_x[k] * sgh_37[k];

        t_45[k] = f_13 * sfh_38[k]
                  + f_3 * pc_x[k] * sgh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, sfh_15, sfh_39, sfh_40, sfh_41, \
                         sgg0_25, sgg1_25, sgh_36, sgh_39, sgh_40, \
                         sgh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_13 * sfh_39[k]
                  + f_3 * pc_x[k] * sgh_39[k];

        t_47[k] = f_13 * sfh_40[k]
                  + f_3 * pc_x[k] * sgh_40[k];

        t_48[k] = f_13 * sfh_41[k]
                  + f_3 * pc_x[k] * sgh_41[k];

        t_49[k] = f_11 * sfh_15[k]
                  + f_1 * sgg0_25[k]
                  - f_2 * sgg1_25[k]
                  + f_3 * pc_y[k] * sgh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pc_y, pc_z, sfh_17, sfh_18, sgg0_27, sgg0_28, \
                         sgg1_27, sgg1_28, sgh_36, sgh_38, sgh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * sgh_36[k];

        t_51[k] = f_11 * sfh_17[k]
                  + f_4 * sgg0_27[k]
                  - f_5 * sgg1_27[k]
                  + f_3 * pc_y[k] * sgh_38[k];

        t_52[k] = f_11 * sfh_18[k]
                  + f_6 * sgg0_28[k]
                  - f_7 * sgg1_28[k]
                  + f_3 * pc_y[k] * sgh_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_y, sfi0_27, sfh_19, sfh_20, sfi1_27, \
                         sgg0_29, sgg1_29, sgh_40, sgh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_11 * sfh_19[k]
                  + f_8 * sgg0_29[k]
                  - f_9 * sgg1_29[k]
                  + f_3 * pc_y[k] * sgh_40[k];

        t_54[k] = f_11 * sfh_20[k]
                  + f_3 * pc_y[k] * sgh_41[k];

        t_55[k] = pb_y[k] * sfi0_27[k]
                  - f_10 * pc_y[k] * sfi1_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pb_z, pc_y, pc_z, sfi0_0, sfi0_3, \
                         sfh_0, sfi1_0, sfi1_3, sgh_42, sgh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * sfi0_0[k]
                  - f_10 * pc_z[k] * sfi1_0[k];

        t_57[k] = f_3 * pc_y[k] * sgh_42[k];

        t_58[k] = f_11 * sfh_0[k]
                  + f_3 * pc_z[k] * sgh_42[k];

        t_59[k] = pb_z[k] * sfi0_3[k]
                  - f_10 * pc_z[k] * sfi1_3[k];

        t_60[k] = f_3 * pc_y[k] * sgh_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, sfi0_5, sfi0_6, sfh_2, \
                         sfh_3, sfi1_5, sfi1_6, sgh_45, sgh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_z[k] * sfi0_5[k]
                  + f_12 * sfh_2[k]
                  - f_10 * pc_z[k] * sfi1_5[k];

        t_62[k] = pb_z[k] * sfi0_6[k]
                  - f_10 * pc_z[k] * sfi1_6[k];

        t_63[k] = f_11 * sfh_3[k]
                  + f_3 * pc_z[k] * sgh_45[k];

        t_64[k] = f_3 * pc_y[k] * sgh_47[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_z, pc_z, sfi0_9, sfi0_10, sfi0_12, sfh_5, \
                         sfh_6, sfh_7, sfi1_9, sfi1_10, sfi1_12, \
                         sgh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_z[k] * sfi0_9[k]
                  + f_13 * sfh_5[k]
                  - f_10 * pc_z[k] * sfi1_9[k];

        t_66[k] = pb_z[k] * sfi0_10[k]
                  - f_10 * pc_z[k] * sfi1_10[k];

        t_67[k] = f_11 * sfh_6[k]
                  + f_3 * pc_z[k] * sgh_48[k];

        t_68[k] = pb_z[k] * sfi0_12[k]
                  + f_12 * sfh_7[k]
                  - f_10 * pc_z[k] * sfi1_12[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_z, pc_x, pc_y, pc_z, sfi0_14, sfh_9, \
                         sfh_57, sfh_58, sfi1_14, sgh_51, sgh_57, \
                         sgh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * pc_y[k] * sgh_51[k];

        t_70[k] = pb_z[k] * sfi0_14[k]
                  + f_0 * sfh_9[k]
                  - f_10 * pc_z[k] * sfi1_14[k];

        t_71[k] = f_13 * sfh_57[k]
                  + f_3 * pc_x[k] * sgh_57[k];

        t_72[k] = f_13 * sfh_58[k]
                  + f_3 * pc_x[k] * sgh_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pc_x, sfh_59, sfh_60, sfh_61, sfh_62, sgh_59, \
                         sgh_60, sgh_61, sgh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_13 * sfh_59[k]
                  + f_3 * pc_x[k] * sgh_59[k];

        t_74[k] = f_13 * sfh_60[k]
                  + f_3 * pc_x[k] * sgh_60[k];

        t_75[k] = f_13 * sfh_61[k]
                  + f_3 * pc_x[k] * sgh_61[k];

        t_76[k] = f_13 * sfh_62[k]
                  + f_3 * pc_x[k] * sgh_62[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_z, pc_y, pc_z, sfi0_21, sfh_15, sfi1_21, \
                         sgg0_42, sgg1_42, sgh_57, sgh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_z[k] * sfi0_21[k]
                  - f_10 * pc_z[k] * sfi1_21[k];

        t_78[k] = f_11 * sfh_15[k]
                  + f_3 * pc_z[k] * sgh_57[k];

        t_79[k] = f_4 * sgg0_42[k]
                  - f_5 * sgg1_42[k]
                  + f_3 * pc_y[k] * sgh_59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_y, pc_z, sfh_20, sgg0_43, sgg0_44, \
                         sgg1_43, sgg1_44, sgh_60, sgh_61, sgh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_6 * sgg0_43[k]
                  - f_7 * sgg1_43[k]
                  + f_3 * pc_y[k] * sgh_60[k];

        t_81[k] = f_8 * sgg0_44[k]
                  - f_9 * sgg1_44[k]
                  + f_3 * pc_y[k] * sgh_61[k];

        t_82[k] = f_3 * pc_y[k] * sgh_62[k];

        t_83[k] = f_11 * sfh_20[k]
                  + f_1 * sgg0_44[k]
                  - f_2 * sgg1_44[k]
                  + f_3 * pc_z[k] * sgh_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pc_x, pc_y, pc_z, sfh_21, sfh_63, sfh_66, \
                         sgg0_45, sgg0_48, sgg1_45, sgg1_48, sgh_63, \
                         sgh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * sfh_63[k]
                  + f_1 * sgg0_45[k]
                  - f_2 * sgg1_45[k]
                  + f_3 * pc_x[k] * sgh_63[k];

        t_85[k] = f_12 * sfh_21[k]
                  + f_3 * pc_y[k] * sgh_63[k];

        t_86[k] = f_3 * pc_z[k] * sgh_63[k];

        t_87[k] = f_12 * sfh_66[k]
                  + f_4 * sgg0_48[k]
                  - f_5 * sgg1_48[k]
                  + f_3 * pc_x[k] * sgh_66[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pc_x, pc_y, sfh_23, sfh_68, sfh_69, sgg0_50, \
                         sgg0_51, sgg1_50, sgg1_51, sgh_65, sgh_68, \
                         sgh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_12 * sfh_23[k]
                  + f_3 * pc_y[k] * sgh_65[k];

        t_89[k] = f_12 * sfh_68[k]
                  + f_4 * sgg0_50[k]
                  - f_5 * sgg1_50[k]
                  + f_3 * pc_x[k] * sgh_68[k];

        t_90[k] = f_12 * sfh_69[k]
                  + f_6 * sgg0_51[k]
                  - f_7 * sgg1_51[k]
                  + f_3 * pc_x[k] * sgh_69[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pc_x, pc_y, pc_z, sfh_26, sfh_72, sgg0_54, sgg1_54, \
                         sgh_66, sgh_68, sgh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * pc_z[k] * sgh_66[k];

        t_92[k] = f_12 * sfh_26[k]
                  + f_3 * pc_y[k] * sgh_68[k];

        t_93[k] = f_12 * sfh_72[k]
                  + f_6 * sgg0_54[k]
                  - f_7 * sgg1_54[k]
                  + f_3 * pc_x[k] * sgh_72[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pc_x, pc_z, sfh_73, sfh_75, sgg0_55, sgg0_57, \
                         sgg1_55, sgg1_57, sgh_69, sgh_73, sgh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_12 * sfh_73[k]
                  + f_8 * sgg0_55[k]
                  - f_9 * sgg1_55[k]
                  + f_3 * pc_x[k] * sgh_73[k];

        t_95[k] = f_3 * pc_z[k] * sgh_69[k];

        t_96[k] = f_12 * sfh_75[k]
                  + f_8 * sgg0_57[k]
                  - f_9 * sgg1_57[k]
                  + f_3 * pc_x[k] * sgh_75[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, sfh_30, sfh_77, sfh_78, sfh_79, \
                         sgg0_59, sgg1_59, sgh_72, sgh_77, sgh_78, \
                         sgh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_12 * sfh_30[k]
                  + f_3 * pc_y[k] * sgh_72[k];

        t_98[k] = f_12 * sfh_77[k]
                  + f_8 * sgg0_59[k]
                  - f_9 * sgg1_59[k]
                  + f_3 * pc_x[k] * sgh_77[k];

        t_99[k] = f_12 * sfh_78[k]
                  + f_3 * pc_x[k] * sgh_78[k];

        t_100[k] = f_12 * sfh_79[k]
                   + f_3 * pc_x[k] * sgh_79[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_x, sfh_80, sfh_81, sfh_82, sfh_83, \
                         sgh_80, sgh_81, sgh_82, sgh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_12 * sfh_80[k]
                   + f_3 * pc_x[k] * sgh_80[k];

        t_102[k] = f_12 * sfh_81[k]
                   + f_3 * pc_x[k] * sgh_81[k];

        t_103[k] = f_12 * sfh_82[k]
                   + f_3 * pc_x[k] * sgh_82[k];

        t_104[k] = f_12 * sfh_83[k]
                   + f_3 * pc_x[k] * sgh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pc_y, pc_z, sfh_36, sfh_38, sgg0_55, sgg0_57, \
                         sgg1_55, sgg1_57, sgh_78, sgh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * sfh_36[k]
                   + f_1 * sgg0_55[k]
                   - f_2 * sgg1_55[k]
                   + f_3 * pc_y[k] * sgh_78[k];

        t_106[k] = f_3 * pc_z[k] * sgh_78[k];

        t_107[k] = f_12 * sfh_38[k]
                   + f_4 * sgg0_57[k]
                   - f_5 * sgg1_57[k]
                   + f_3 * pc_y[k] * sgh_80[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pc_y, pc_z, sfh_39, sfh_40, sfh_41, \
                         sgg0_58, sgg0_59, sgg1_58, sgg1_59, sgh_81, sgh_82, \
                         sgh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_12 * sfh_39[k]
                   + f_6 * sgg0_58[k]
                   - f_7 * sgg1_58[k]
                   + f_3 * pc_y[k] * sgh_81[k];

        t_109[k] = f_12 * sfh_40[k]
                   + f_8 * sgg0_59[k]
                   - f_9 * sgg1_59[k]
                   + f_3 * pc_y[k] * sgh_82[k];

        t_110[k] = f_12 * sfh_41[k]
                   + f_3 * pc_y[k] * sgh_83[k];

        t_111[k] = f_1 * sgg0_59[k]
                   - f_2 * sgg1_59[k]
                   + f_3 * pc_z[k] * sgh_83[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pb_y, pb_z, pc_y, pc_z, sfi0_31, sfi0_56, \
                         sfh_21, sfh_42, sfi1_31, sfi1_56, sgh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = pb_y[k] * sfi0_56[k]
                   - f_10 * pc_y[k] * sfi1_56[k];

        t_113[k] = f_11 * sfh_42[k]
                   + f_3 * pc_y[k] * sgh_84[k];

        t_114[k] = f_11 * sfh_21[k]
                   + f_3 * pc_z[k] * sgh_84[k];

        t_115[k] = pb_z[k] * sfi0_31[k]
                   - f_10 * pc_z[k] * sfi1_31[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pb_z, pc_y, pc_z, sfi0_34, sfi0_61, \
                         sfh_24, sfh_44, sfi1_34, sfi1_61, sgh_86, \
                         sgh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_11 * sfh_44[k]
                   + f_3 * pc_y[k] * sgh_86[k];

        t_117[k] = pb_y[k] * sfi0_61[k]
                   - f_10 * pc_y[k] * sfi1_61[k];

        t_118[k] = pb_z[k] * sfi0_34[k]
                   - f_10 * pc_z[k] * sfi1_34[k];

        t_119[k] = f_11 * sfh_24[k]
                   + f_3 * pc_z[k] * sgh_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pb_z, pc_y, pc_z, sfi0_38, sfi0_65, \
                         sfh_27, sfh_47, sfi1_38, sfi1_65, sgh_89, \
                         sgh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_11 * sfh_47[k]
                   + f_3 * pc_y[k] * sgh_89[k];

        t_121[k] = pb_y[k] * sfi0_65[k]
                   - f_10 * pc_y[k] * sfi1_65[k];

        t_122[k] = pb_z[k] * sfi0_38[k]
                   - f_10 * pc_z[k] * sfi1_38[k];

        t_123[k] = f_11 * sfh_27[k]
                   + f_3 * pc_z[k] * sgh_90[k];
    }
}

static auto
compute_prim_sgi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfi0,
                                                          const size_t sfh, const size_t sfi1,
                                                          const size_t sgg0, const size_t sgg1,
                                                          const size_t sgh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfi0_49 = buffer.data(sfi0 + 49);
    const auto *sfi0_68 = buffer.data(sfi0 + 68);
    const auto *sfi0_70 = buffer.data(sfi0 + 70);
    const auto *sfi0_83 = buffer.data(sfi0 + 83);
    const auto *sfi0_84 = buffer.data(sfi0 + 84);
    const auto *sfi0_87 = buffer.data(sfi0 + 87);
    const auto *sfi0_90 = buffer.data(sfi0 + 90);
    const auto *sfi0_94 = buffer.data(sfi0 + 94);
    const auto *sfi0_140 = buffer.data(sfi0 + 140);
    const auto *sfi0_145 = buffer.data(sfi0 + 145);
    const auto *sfi0_149 = buffer.data(sfi0 + 149);
    const auto *sfi0_154 = buffer.data(sfi0 + 154);
    const auto *sfi0_168 = buffer.data(sfi0 + 168);
    const auto *sfi0_171 = buffer.data(sfi0 + 171);
    const auto *sfi0_173 = buffer.data(sfi0 + 173);
    const auto *sfi0_174 = buffer.data(sfi0 + 174);
    const auto *sfi0_177 = buffer.data(sfi0 + 177);
    const auto *sfi0_178 = buffer.data(sfi0 + 178);
    const auto *sfi0_180 = buffer.data(sfi0 + 180);
    const auto *sfi0_182 = buffer.data(sfi0 + 182);
    const auto *sfi0_189 = buffer.data(sfi0 + 189);
    const auto *sfi0_191 = buffer.data(sfi0 + 191);
    const auto *sfi0_192 = buffer.data(sfi0 + 192);
    const auto *sfi0_193 = buffer.data(sfi0 + 193);
    const auto *sfi0_195 = buffer.data(sfi0 + 195);
    const auto *sfi0_201 = buffer.data(sfi0 + 201);
    const auto *sfi0_205 = buffer.data(sfi0 + 205);
    const auto *sfi0_208 = buffer.data(sfi0 + 208);
    const auto *sfi0_210 = buffer.data(sfi0 + 210);
    const auto *sfi0_217 = buffer.data(sfi0 + 217);
    const auto *sfi0_219 = buffer.data(sfi0 + 219);
    const auto *sfi0_220 = buffer.data(sfi0 + 220);
    const auto *sfi0_221 = buffer.data(sfi0 + 221);
    const auto *sfi0_223 = buffer.data(sfi0 + 223);
    const auto *sfi0_227 = buffer.data(sfi0 + 227);
    const auto *sfi0_230 = buffer.data(sfi0 + 230);
    const auto *sfi0_234 = buffer.data(sfi0 + 234);
    const auto *sfi0_236 = buffer.data(sfi0 + 236);
    const auto *sfi0_245 = buffer.data(sfi0 + 245);

    const auto *sfh_36 = buffer.data(sfh + 36);
    const auto *sfh_42 = buffer.data(sfh + 42);
    const auto *sfh_45 = buffer.data(sfh + 45);
    const auto *sfh_48 = buffer.data(sfh + 48);
    const auto *sfh_50 = buffer.data(sfh + 50);
    const auto *sfh_51 = buffer.data(sfh + 51);
    const auto *sfh_57 = buffer.data(sfh + 57);
    const auto *sfh_59 = buffer.data(sfh + 59);
    const auto *sfh_60 = buffer.data(sfh + 60);
    const auto *sfh_61 = buffer.data(sfh + 61);
    const auto *sfh_62 = buffer.data(sfh + 62);
    const auto *sfh_63 = buffer.data(sfh + 63);
    const auto *sfh_65 = buffer.data(sfh + 65);
    const auto *sfh_66 = buffer.data(sfh + 66);
    const auto *sfh_68 = buffer.data(sfh + 68);
    const auto *sfh_69 = buffer.data(sfh + 69);
    const auto *sfh_72 = buffer.data(sfh + 72);
    const auto *sfh_78 = buffer.data(sfh + 78);
    const auto *sfh_83 = buffer.data(sfh + 83);
    const auto *sfh_84 = buffer.data(sfh + 84);
    const auto *sfh_86 = buffer.data(sfh + 86);
    const auto *sfh_87 = buffer.data(sfh + 87);
    const auto *sfh_89 = buffer.data(sfh + 89);
    const auto *sfh_90 = buffer.data(sfh + 90);
    const auto *sfh_93 = buffer.data(sfh + 93);
    const auto *sfh_99 = buffer.data(sfh + 99);
    const auto *sfh_100 = buffer.data(sfh + 100);
    const auto *sfh_101 = buffer.data(sfh + 101);
    const auto *sfh_102 = buffer.data(sfh + 102);
    const auto *sfh_103 = buffer.data(sfh + 103);
    const auto *sfh_104 = buffer.data(sfh + 104);
    const auto *sfh_105 = buffer.data(sfh + 105);
    const auto *sfh_107 = buffer.data(sfh + 107);
    const auto *sfh_108 = buffer.data(sfh + 108);
    const auto *sfh_110 = buffer.data(sfh + 110);
    const auto *sfh_111 = buffer.data(sfh + 111);
    const auto *sfh_114 = buffer.data(sfh + 114);
    const auto *sfh_115 = buffer.data(sfh + 115);
    const auto *sfh_117 = buffer.data(sfh + 117);
    const auto *sfh_119 = buffer.data(sfh + 119);
    const auto *sfh_120 = buffer.data(sfh + 120);
    const auto *sfh_121 = buffer.data(sfh + 121);
    const auto *sfh_122 = buffer.data(sfh + 122);
    const auto *sfh_123 = buffer.data(sfh + 123);
    const auto *sfh_124 = buffer.data(sfh + 124);
    const auto *sfh_125 = buffer.data(sfh + 125);
    const auto *sfh_126 = buffer.data(sfh + 126);
    const auto *sfh_129 = buffer.data(sfh + 129);
    const auto *sfh_131 = buffer.data(sfh + 131);
    const auto *sfh_132 = buffer.data(sfh + 132);
    const auto *sfh_135 = buffer.data(sfh + 135);
    const auto *sfh_136 = buffer.data(sfh + 136);
    const auto *sfh_138 = buffer.data(sfh + 138);
    const auto *sfh_140 = buffer.data(sfh + 140);
    const auto *sfh_141 = buffer.data(sfh + 141);
    const auto *sfh_142 = buffer.data(sfh + 142);
    const auto *sfh_143 = buffer.data(sfh + 143);
    const auto *sfh_144 = buffer.data(sfh + 144);
    const auto *sfh_145 = buffer.data(sfh + 145);
    const auto *sfh_146 = buffer.data(sfh + 146);
    const auto *sfh_152 = buffer.data(sfh + 152);
    const auto *sfh_156 = buffer.data(sfh + 156);
    const auto *sfh_159 = buffer.data(sfh + 159);
    const auto *sfh_161 = buffer.data(sfh + 161);
    const auto *sfh_162 = buffer.data(sfh + 162);
    const auto *sfh_163 = buffer.data(sfh + 163);
    const auto *sfh_164 = buffer.data(sfh + 164);
    const auto *sfh_165 = buffer.data(sfh + 165);
    const auto *sfh_166 = buffer.data(sfh + 166);
    const auto *sfh_167 = buffer.data(sfh + 167);
    const auto *sfh_171 = buffer.data(sfh + 171);
    const auto *sfh_174 = buffer.data(sfh + 174);
    const auto *sfh_178 = buffer.data(sfh + 178);
    const auto *sfh_180 = buffer.data(sfh + 180);
    const auto *sfh_183 = buffer.data(sfh + 183);
    const auto *sfh_184 = buffer.data(sfh + 184);
    const auto *sfh_185 = buffer.data(sfh + 185);
    const auto *sfh_186 = buffer.data(sfh + 186);
    const auto *sfh_187 = buffer.data(sfh + 187);
    const auto *sfh_188 = buffer.data(sfh + 188);

    const auto *sfi1_49 = buffer.data(sfi1 + 49);
    const auto *sfi1_68 = buffer.data(sfi1 + 68);
    const auto *sfi1_70 = buffer.data(sfi1 + 70);
    const auto *sfi1_83 = buffer.data(sfi1 + 83);
    const auto *sfi1_84 = buffer.data(sfi1 + 84);
    const auto *sfi1_87 = buffer.data(sfi1 + 87);
    const auto *sfi1_90 = buffer.data(sfi1 + 90);
    const auto *sfi1_94 = buffer.data(sfi1 + 94);
    const auto *sfi1_140 = buffer.data(sfi1 + 140);
    const auto *sfi1_145 = buffer.data(sfi1 + 145);
    const auto *sfi1_149 = buffer.data(sfi1 + 149);
    const auto *sfi1_154 = buffer.data(sfi1 + 154);
    const auto *sfi1_168 = buffer.data(sfi1 + 168);
    const auto *sfi1_171 = buffer.data(sfi1 + 171);
    const auto *sfi1_173 = buffer.data(sfi1 + 173);
    const auto *sfi1_174 = buffer.data(sfi1 + 174);
    const auto *sfi1_177 = buffer.data(sfi1 + 177);
    const auto *sfi1_178 = buffer.data(sfi1 + 178);
    const auto *sfi1_180 = buffer.data(sfi1 + 180);
    const auto *sfi1_182 = buffer.data(sfi1 + 182);
    const auto *sfi1_189 = buffer.data(sfi1 + 189);
    const auto *sfi1_191 = buffer.data(sfi1 + 191);
    const auto *sfi1_192 = buffer.data(sfi1 + 192);
    const auto *sfi1_193 = buffer.data(sfi1 + 193);
    const auto *sfi1_195 = buffer.data(sfi1 + 195);
    const auto *sfi1_201 = buffer.data(sfi1 + 201);
    const auto *sfi1_205 = buffer.data(sfi1 + 205);
    const auto *sfi1_208 = buffer.data(sfi1 + 208);
    const auto *sfi1_210 = buffer.data(sfi1 + 210);
    const auto *sfi1_217 = buffer.data(sfi1 + 217);
    const auto *sfi1_219 = buffer.data(sfi1 + 219);
    const auto *sfi1_220 = buffer.data(sfi1 + 220);
    const auto *sfi1_221 = buffer.data(sfi1 + 221);
    const auto *sfi1_223 = buffer.data(sfi1 + 223);
    const auto *sfi1_227 = buffer.data(sfi1 + 227);
    const auto *sfi1_230 = buffer.data(sfi1 + 230);
    const auto *sfi1_234 = buffer.data(sfi1 + 234);
    const auto *sfi1_236 = buffer.data(sfi1 + 236);
    const auto *sfi1_245 = buffer.data(sfi1 + 245);

    const auto *sgg0_72 = buffer.data(sgg0 + 72);
    const auto *sgg0_73 = buffer.data(sgg0 + 73);
    const auto *sgg0_74 = buffer.data(sgg0 + 74);
    const auto *sgg0_75 = buffer.data(sgg0 + 75);
    const auto *sgg0_78 = buffer.data(sgg0 + 78);
    const auto *sgg0_80 = buffer.data(sgg0 + 80);
    const auto *sgg0_81 = buffer.data(sgg0 + 81);
    const auto *sgg0_84 = buffer.data(sgg0 + 84);
    const auto *sgg0_85 = buffer.data(sgg0 + 85);
    const auto *sgg0_87 = buffer.data(sgg0 + 87);
    const auto *sgg0_88 = buffer.data(sgg0 + 88);
    const auto *sgg0_89 = buffer.data(sgg0 + 89);

    const auto *sgg1_72 = buffer.data(sgg1 + 72);
    const auto *sgg1_73 = buffer.data(sgg1 + 73);
    const auto *sgg1_74 = buffer.data(sgg1 + 74);
    const auto *sgg1_75 = buffer.data(sgg1 + 75);
    const auto *sgg1_78 = buffer.data(sgg1 + 78);
    const auto *sgg1_80 = buffer.data(sgg1 + 80);
    const auto *sgg1_81 = buffer.data(sgg1 + 81);
    const auto *sgg1_84 = buffer.data(sgg1 + 84);
    const auto *sgg1_85 = buffer.data(sgg1 + 85);
    const auto *sgg1_87 = buffer.data(sgg1 + 87);
    const auto *sgg1_88 = buffer.data(sgg1 + 88);
    const auto *sgg1_89 = buffer.data(sgg1 + 89);

    const auto *sgh_93 = buffer.data(sgh + 93);
    const auto *sgh_99 = buffer.data(sgh + 99);
    const auto *sgh_100 = buffer.data(sgh + 100);
    const auto *sgh_101 = buffer.data(sgh + 101);
    const auto *sgh_102 = buffer.data(sgh + 102);
    const auto *sgh_103 = buffer.data(sgh + 103);
    const auto *sgh_104 = buffer.data(sgh + 104);
    const auto *sgh_105 = buffer.data(sgh + 105);
    const auto *sgh_107 = buffer.data(sgh + 107);
    const auto *sgh_108 = buffer.data(sgh + 108);
    const auto *sgh_110 = buffer.data(sgh + 110);
    const auto *sgh_111 = buffer.data(sgh + 111);
    const auto *sgh_114 = buffer.data(sgh + 114);
    const auto *sgh_115 = buffer.data(sgh + 115);
    const auto *sgh_117 = buffer.data(sgh + 117);
    const auto *sgh_119 = buffer.data(sgh + 119);
    const auto *sgh_120 = buffer.data(sgh + 120);
    const auto *sgh_121 = buffer.data(sgh + 121);
    const auto *sgh_122 = buffer.data(sgh + 122);
    const auto *sgh_123 = buffer.data(sgh + 123);
    const auto *sgh_124 = buffer.data(sgh + 124);
    const auto *sgh_125 = buffer.data(sgh + 125);
    const auto *sgh_126 = buffer.data(sgh + 126);
    const auto *sgh_128 = buffer.data(sgh + 128);
    const auto *sgh_129 = buffer.data(sgh + 129);
    const auto *sgh_131 = buffer.data(sgh + 131);
    const auto *sgh_132 = buffer.data(sgh + 132);
    const auto *sgh_135 = buffer.data(sgh + 135);
    const auto *sgh_141 = buffer.data(sgh + 141);
    const auto *sgh_142 = buffer.data(sgh + 142);
    const auto *sgh_143 = buffer.data(sgh + 143);
    const auto *sgh_144 = buffer.data(sgh + 144);
    const auto *sgh_145 = buffer.data(sgh + 145);
    const auto *sgh_146 = buffer.data(sgh + 146);
    const auto *sgh_147 = buffer.data(sgh + 147);
    const auto *sgh_149 = buffer.data(sgh + 149);
    const auto *sgh_150 = buffer.data(sgh + 150);
    const auto *sgh_152 = buffer.data(sgh + 152);
    const auto *sgh_153 = buffer.data(sgh + 153);
    const auto *sgh_156 = buffer.data(sgh + 156);
    const auto *sgh_162 = buffer.data(sgh + 162);
    const auto *sgh_163 = buffer.data(sgh + 163);
    const auto *sgh_164 = buffer.data(sgh + 164);
    const auto *sgh_165 = buffer.data(sgh + 165);
    const auto *sgh_166 = buffer.data(sgh + 166);
    const auto *sgh_167 = buffer.data(sgh + 167);
    const auto *sgh_168 = buffer.data(sgh + 168);
    const auto *sgh_170 = buffer.data(sgh + 170);
    const auto *sgh_171 = buffer.data(sgh + 171);
    const auto *sgh_173 = buffer.data(sgh + 173);
    const auto *sgh_174 = buffer.data(sgh + 174);
    const auto *sgh_177 = buffer.data(sgh + 177);
    const auto *sgh_183 = buffer.data(sgh + 183);
    const auto *sgh_184 = buffer.data(sgh + 184);
    const auto *sgh_185 = buffer.data(sgh + 185);
    const auto *sgh_186 = buffer.data(sgh + 186);
    const auto *sgh_187 = buffer.data(sgh + 187);
    const auto *sgh_188 = buffer.data(sgh + 188);

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_y, pc_x, pc_y, sfi0_68, sfi0_70, \
                         sfh_50, sfh_51, sfh_99, sfi1_68, sfi1_70, sgh_93, \
                         sgh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * sfi0_68[k]
                   + f_12 * sfh_50[k]
                   - f_10 * pc_y[k] * sfi1_68[k];

        t_125[k] = f_11 * sfh_51[k]
                   + f_3 * pc_y[k] * sgh_93[k];

        t_126[k] = pb_y[k] * sfi0_70[k]
                   - f_10 * pc_y[k] * sfi1_70[k];

        t_127[k] = f_12 * sfh_99[k]
                   + f_3 * pc_x[k] * sgh_99[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, sfh_100, sfh_101, sfh_102, \
                         sfh_103, sfh_104, sgh_100, sgh_101, sgh_102, sgh_103, \
                         sgh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_12 * sfh_100[k]
                   + f_3 * pc_x[k] * sgh_100[k];

        t_129[k] = f_12 * sfh_101[k]
                   + f_3 * pc_x[k] * sgh_101[k];

        t_130[k] = f_12 * sfh_102[k]
                   + f_3 * pc_x[k] * sgh_102[k];

        t_131[k] = f_12 * sfh_103[k]
                   + f_3 * pc_x[k] * sgh_103[k];

        t_132[k] = f_12 * sfh_104[k]
                   + f_3 * pc_x[k] * sgh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_z, pc_y, pc_z, sfi0_49, sfh_36, sfh_59, \
                         sfi1_49, sgg0_72, sgg1_72, sgh_99, sgh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pb_z[k] * sfi0_49[k]
                   - f_10 * pc_z[k] * sfi1_49[k];

        t_134[k] = f_11 * sfh_36[k]
                   + f_3 * pc_z[k] * sgh_99[k];

        t_135[k] = f_11 * sfh_59[k]
                   + f_4 * sgg0_72[k]
                   - f_5 * sgg1_72[k]
                   + f_3 * pc_y[k] * sgh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, sfh_60, sfh_61, sfh_62, sgg0_73, sgg0_74, \
                         sgg1_73, sgg1_74, sgh_102, sgh_103, sgh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * sfh_60[k]
                   + f_6 * sgg0_73[k]
                   - f_7 * sgg1_73[k]
                   + f_3 * pc_y[k] * sgh_102[k];

        t_137[k] = f_11 * sfh_61[k]
                   + f_8 * sgg0_74[k]
                   - f_9 * sgg1_74[k]
                   + f_3 * pc_y[k] * sgh_103[k];

        t_138[k] = f_11 * sfh_62[k]
                   + f_3 * pc_y[k] * sgh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pb_y, pc_x, pc_y, pc_z, sfi0_83, sfh_42, \
                         sfh_105, sfi1_83, sgg0_75, sgg1_75, sgh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * sfi0_83[k]
                   - f_10 * pc_y[k] * sfi1_83[k];

        t_140[k] = f_12 * sfh_105[k]
                   + f_1 * sgg0_75[k]
                   - f_2 * sgg1_75[k]
                   + f_3 * pc_x[k] * sgh_105[k];

        t_141[k] = f_3 * pc_y[k] * sgh_105[k];

        t_142[k] = f_12 * sfh_42[k]
                   + f_3 * pc_z[k] * sgh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, sfh_108, sfh_110, sgg0_78, sgg0_80, \
                         sgg1_78, sgg1_80, sgh_107, sgh_108, sgh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_12 * sfh_108[k]
                   + f_4 * sgg0_78[k]
                   - f_5 * sgg1_78[k]
                   + f_3 * pc_x[k] * sgh_108[k];

        t_144[k] = f_3 * pc_y[k] * sgh_107[k];

        t_145[k] = f_12 * sfh_110[k]
                   + f_4 * sgg0_80[k]
                   - f_5 * sgg1_80[k]
                   + f_3 * pc_x[k] * sgh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_x, pc_y, pc_z, sfh_45, sfh_111, sgg0_81, \
                         sgg1_81, sgh_108, sgh_110, sgh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_12 * sfh_111[k]
                   + f_6 * sgg0_81[k]
                   - f_7 * sgg1_81[k]
                   + f_3 * pc_x[k] * sgh_111[k];

        t_147[k] = f_12 * sfh_45[k]
                   + f_3 * pc_z[k] * sgh_108[k];

        t_148[k] = f_3 * pc_y[k] * sgh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_z, sfh_48, sfh_114, sfh_115, sgg0_84, \
                         sgg0_85, sgg1_84, sgg1_85, sgh_111, sgh_114, \
                         sgh_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_12 * sfh_114[k]
                   + f_6 * sgg0_84[k]
                   - f_7 * sgg1_84[k]
                   + f_3 * pc_x[k] * sgh_114[k];

        t_150[k] = f_12 * sfh_115[k]
                   + f_8 * sgg0_85[k]
                   - f_9 * sgg1_85[k]
                   + f_3 * pc_x[k] * sgh_115[k];

        t_151[k] = f_12 * sfh_48[k]
                   + f_3 * pc_z[k] * sgh_111[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pc_x, pc_y, sfh_117, sfh_119, sgg0_87, sgg0_89, \
                         sgg1_87, sgg1_89, sgh_114, sgh_117, sgh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_12 * sfh_117[k]
                   + f_8 * sgg0_87[k]
                   - f_9 * sgg1_87[k]
                   + f_3 * pc_x[k] * sgh_117[k];

        t_153[k] = f_3 * pc_y[k] * sgh_114[k];

        t_154[k] = f_12 * sfh_119[k]
                   + f_8 * sgg0_89[k]
                   - f_9 * sgg1_89[k]
                   + f_3 * pc_x[k] * sgh_119[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pc_x, sfh_120, sfh_121, sfh_122, \
                         sfh_123, sfh_124, sgh_120, sgh_121, sgh_122, sgh_123, \
                         sgh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_12 * sfh_120[k]
                   + f_3 * pc_x[k] * sgh_120[k];

        t_156[k] = f_12 * sfh_121[k]
                   + f_3 * pc_x[k] * sgh_121[k];

        t_157[k] = f_12 * sfh_122[k]
                   + f_3 * pc_x[k] * sgh_122[k];

        t_158[k] = f_12 * sfh_123[k]
                   + f_3 * pc_x[k] * sgh_123[k];

        t_159[k] = f_12 * sfh_124[k]
                   + f_3 * pc_x[k] * sgh_124[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, sfh_57, sfh_125, \
                         sgg0_85, sgg0_87, sgg1_85, sgg1_87, sgh_120, sgh_122, \
                         sgh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_12 * sfh_125[k]
                   + f_3 * pc_x[k] * sgh_125[k];

        t_161[k] = f_1 * sgg0_85[k]
                   - f_2 * sgg1_85[k]
                   + f_3 * pc_y[k] * sgh_120[k];

        t_162[k] = f_12 * sfh_57[k]
                   + f_3 * pc_z[k] * sgh_120[k];

        t_163[k] = f_4 * sgg0_87[k]
                   - f_5 * sgg1_87[k]
                   + f_3 * pc_y[k] * sgh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, sfh_62, sgg0_88, sgg0_89, \
                         sgg1_88, sgg1_89, sgh_123, sgh_124, sgh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * sgg0_88[k]
                   - f_7 * sgg1_88[k]
                   + f_3 * pc_y[k] * sgh_123[k];

        t_165[k] = f_8 * sgg0_89[k]
                   - f_9 * sgg1_89[k]
                   + f_3 * pc_y[k] * sgh_124[k];

        t_166[k] = f_3 * pc_y[k] * sgh_125[k];

        t_167[k] = f_12 * sfh_62[k]
                   + f_1 * sgg0_89[k]
                   - f_2 * sgg1_89[k]
                   + f_3 * pc_z[k] * sgh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pb_x, pc_x, pc_y, pc_z, sfi0_168, \
                         sfi0_171, sfh_63, sfh_126, sfh_129, sfi1_168, sfi1_171, \
                         sgh_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pb_x[k] * sfi0_168[k]
                   + f_14 * sfh_126[k]
                   - f_10 * pc_x[k] * sfi1_168[k];

        t_169[k] = f_13 * sfh_63[k]
                   + f_3 * pc_y[k] * sgh_126[k];

        t_170[k] = f_3 * pc_z[k] * sgh_126[k];

        t_171[k] = pb_x[k] * sfi0_171[k]
                   + f_0 * sfh_129[k]
                   - f_10 * pc_x[k] * sfi1_171[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pb_x, pc_x, pc_y, sfi0_173, sfi0_174, sfh_65, \
                         sfh_131, sfh_132, sfi1_173, sfi1_174, \
                         sgh_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_13 * sfh_65[k]
                   + f_3 * pc_y[k] * sgh_128[k];

        t_173[k] = pb_x[k] * sfi0_173[k]
                   + f_0 * sfh_131[k]
                   - f_10 * pc_x[k] * sfi1_173[k];

        t_174[k] = pb_x[k] * sfi0_174[k]
                   + f_13 * sfh_132[k]
                   - f_10 * pc_x[k] * sfi1_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pb_x, pc_x, pc_y, pc_z, sfi0_177, sfh_68, \
                         sfh_135, sfi1_177, sgh_129, sgh_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_3 * pc_z[k] * sgh_129[k];

        t_176[k] = f_13 * sfh_68[k]
                   + f_3 * pc_y[k] * sgh_131[k];

        t_177[k] = pb_x[k] * sfi0_177[k]
                   + f_13 * sfh_135[k]
                   - f_10 * pc_x[k] * sfi1_177[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pb_x, pc_x, pc_z, sfi0_178, sfi0_180, sfh_136, \
                         sfh_138, sfi1_178, sfi1_180, sgh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = pb_x[k] * sfi0_178[k]
                   + f_12 * sfh_136[k]
                   - f_10 * pc_x[k] * sfi1_178[k];

        t_179[k] = f_3 * pc_z[k] * sgh_132[k];

        t_180[k] = pb_x[k] * sfi0_180[k]
                   + f_12 * sfh_138[k]
                   - f_10 * pc_x[k] * sfi1_180[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_x, pc_x, pc_y, sfi0_182, sfh_72, \
                         sfh_140, sfh_141, sfh_142, sfi1_182, sgh_135, sgh_141, \
                         sgh_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_13 * sfh_72[k]
                   + f_3 * pc_y[k] * sgh_135[k];

        t_182[k] = pb_x[k] * sfi0_182[k]
                   + f_12 * sfh_140[k]
                   - f_10 * pc_x[k] * sfi1_182[k];

        t_183[k] = f_11 * sfh_141[k]
                   + f_3 * pc_x[k] * sgh_141[k];

        t_184[k] = f_11 * sfh_142[k]
                   + f_3 * pc_x[k] * sgh_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, sfh_143, sfh_144, sfh_145, sfh_146, \
                         sgh_143, sgh_144, sgh_145, sgh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_11 * sfh_143[k]
                   + f_3 * pc_x[k] * sgh_143[k];

        t_186[k] = f_11 * sfh_144[k]
                   + f_3 * pc_x[k] * sgh_144[k];

        t_187[k] = f_11 * sfh_145[k]
                   + f_3 * pc_x[k] * sgh_145[k];

        t_188[k] = f_11 * sfh_146[k]
                   + f_3 * pc_x[k] * sgh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_x, pc_x, pc_z, sfi0_189, sfi0_191, \
                         sfi0_192, sfi1_189, sfi1_191, sfi1_192, \
                         sgh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pb_x[k] * sfi0_189[k]
                   - f_10 * pc_x[k] * sfi1_189[k];

        t_190[k] = f_3 * pc_z[k] * sgh_141[k];

        t_191[k] = pb_x[k] * sfi0_191[k]
                   - f_10 * pc_x[k] * sfi1_191[k];

        t_192[k] = pb_x[k] * sfi0_192[k]
                   - f_10 * pc_x[k] * sfi1_192[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pb_x, pc_x, pc_y, sfi0_193, sfi0_195, sfh_83, \
                         sfi1_193, sfi1_195, sgh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_x[k] * sfi0_193[k]
                   - f_10 * pc_x[k] * sfi1_193[k];

        t_194[k] = f_13 * sfh_83[k]
                   + f_3 * pc_y[k] * sgh_146[k];

        t_195[k] = pb_x[k] * sfi0_195[k]
                   - f_10 * pc_x[k] * sfi1_195[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pb_z, pc_y, pc_z, sfi0_84, sfi0_87, \
                         sfh_63, sfh_84, sfi1_84, sfi1_87, sgh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pb_z[k] * sfi0_84[k]
                   - f_10 * pc_z[k] * sfi1_84[k];

        t_197[k] = f_12 * sfh_84[k]
                   + f_3 * pc_y[k] * sgh_147[k];

        t_198[k] = f_11 * sfh_63[k]
                   + f_3 * pc_z[k] * sgh_147[k];

        t_199[k] = pb_z[k] * sfi0_87[k]
                   - f_10 * pc_z[k] * sfi1_87[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_x, pb_z, pc_x, pc_y, pc_z, sfi0_90, sfi0_201, \
                         sfh_86, sfh_152, sfi1_90, sfi1_201, sgh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_12 * sfh_86[k]
                   + f_3 * pc_y[k] * sgh_149[k];

        t_201[k] = pb_x[k] * sfi0_201[k]
                   + f_0 * sfh_152[k]
                   - f_10 * pc_x[k] * sfi1_201[k];

        t_202[k] = pb_z[k] * sfi0_90[k]
                   - f_10 * pc_z[k] * sfi1_90[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pb_x, pc_x, pc_y, pc_z, sfi0_205, sfh_66, \
                         sfh_89, sfh_156, sfi1_205, sgh_150, sgh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_11 * sfh_66[k]
                   + f_3 * pc_z[k] * sgh_150[k];

        t_204[k] = f_12 * sfh_89[k]
                   + f_3 * pc_y[k] * sgh_152[k];

        t_205[k] = pb_x[k] * sfi0_205[k]
                   + f_13 * sfh_156[k]
                   - f_10 * pc_x[k] * sfi1_205[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, pb_x, pb_z, pc_x, pc_z, sfi0_94, sfi0_208, \
                         sfh_69, sfh_159, sfi1_94, sfi1_208, sgh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pb_z[k] * sfi0_94[k]
                   - f_10 * pc_z[k] * sfi1_94[k];

        t_207[k] = f_11 * sfh_69[k]
                   + f_3 * pc_z[k] * sgh_153[k];

        t_208[k] = pb_x[k] * sfi0_208[k]
                   + f_12 * sfh_159[k]
                   - f_10 * pc_x[k] * sfi1_208[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pb_x, pc_x, pc_y, sfi0_210, sfh_93, \
                         sfh_161, sfh_162, sfh_163, sfi1_210, sgh_156, sgh_162, \
                         sgh_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_12 * sfh_93[k]
                   + f_3 * pc_y[k] * sgh_156[k];

        t_210[k] = pb_x[k] * sfi0_210[k]
                   + f_12 * sfh_161[k]
                   - f_10 * pc_x[k] * sfi1_210[k];

        t_211[k] = f_11 * sfh_162[k]
                   + f_3 * pc_x[k] * sgh_162[k];

        t_212[k] = f_11 * sfh_163[k]
                   + f_3 * pc_x[k] * sgh_163[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, sfh_164, sfh_165, sfh_166, sfh_167, \
                         sgh_164, sgh_165, sgh_166, sgh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_11 * sfh_164[k]
                   + f_3 * pc_x[k] * sgh_164[k];

        t_214[k] = f_11 * sfh_165[k]
                   + f_3 * pc_x[k] * sgh_165[k];

        t_215[k] = f_11 * sfh_166[k]
                   + f_3 * pc_x[k] * sgh_166[k];

        t_216[k] = f_11 * sfh_167[k]
                   + f_3 * pc_x[k] * sgh_167[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pb_x, pc_x, pc_z, sfi0_217, sfi0_219, \
                         sfi0_220, sfh_78, sfi1_217, sfi1_219, sfi1_220, \
                         sgh_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pb_x[k] * sfi0_217[k]
                   - f_10 * pc_x[k] * sfi1_217[k];

        t_218[k] = f_11 * sfh_78[k]
                   + f_3 * pc_z[k] * sgh_162[k];

        t_219[k] = pb_x[k] * sfi0_219[k]
                   - f_10 * pc_x[k] * sfi1_219[k];

        t_220[k] = pb_x[k] * sfi0_220[k]
                   - f_10 * pc_x[k] * sfi1_220[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pb_x, pb_y, pc_x, pc_y, sfi0_140, \
                         sfi0_221, sfi0_223, sfh_104, sfi1_140, sfi1_221, sfi1_223, \
                         sgh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = pb_x[k] * sfi0_221[k]
                   - f_10 * pc_x[k] * sfi1_221[k];

        t_222[k] = f_12 * sfh_104[k]
                   + f_3 * pc_y[k] * sgh_167[k];

        t_223[k] = pb_x[k] * sfi0_223[k]
                   - f_10 * pc_x[k] * sfi1_223[k];

        t_224[k] = pb_y[k] * sfi0_140[k]
                   - f_10 * pc_y[k] * sfi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_x, pc_x, pc_y, pc_z, sfi0_227, sfh_84, \
                         sfh_105, sfh_107, sfh_171, sfi1_227, sgh_168, \
                         sgh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * sfh_105[k]
                   + f_3 * pc_y[k] * sgh_168[k];

        t_226[k] = f_12 * sfh_84[k]
                   + f_3 * pc_z[k] * sgh_168[k];

        t_227[k] = pb_x[k] * sfi0_227[k]
                   + f_0 * sfh_171[k]
                   - f_10 * pc_x[k] * sfi1_227[k];

        t_228[k] = f_11 * sfh_107[k]
                   + f_3 * pc_y[k] * sgh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pb_x, pb_y, pc_x, pc_y, pc_z, sfi0_145, \
                         sfi0_230, sfh_87, sfh_174, sfi1_145, sfi1_230, \
                         sgh_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * sfi0_145[k]
                   - f_10 * pc_y[k] * sfi1_145[k];

        t_230[k] = pb_x[k] * sfi0_230[k]
                   + f_13 * sfh_174[k]
                   - f_10 * pc_x[k] * sfi1_230[k];

        t_231[k] = f_12 * sfh_87[k]
                   + f_3 * pc_z[k] * sgh_171[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pb_x, pb_y, pc_x, pc_y, sfi0_149, sfi0_234, \
                         sfh_110, sfh_178, sfi1_149, sfi1_234, \
                         sgh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_11 * sfh_110[k]
                   + f_3 * pc_y[k] * sgh_173[k];

        t_233[k] = pb_y[k] * sfi0_149[k]
                   - f_10 * pc_y[k] * sfi1_149[k];

        t_234[k] = pb_x[k] * sfi0_234[k]
                   + f_12 * sfh_178[k]
                   - f_10 * pc_x[k] * sfi1_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, pb_x, pc_x, pc_y, pc_z, sfi0_236, sfh_90, \
                         sfh_114, sfh_180, sfi1_236, sgh_174, sgh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_12 * sfh_90[k]
                   + f_3 * pc_z[k] * sgh_174[k];

        t_236[k] = pb_x[k] * sfi0_236[k]
                   + f_12 * sfh_180[k]
                   - f_10 * pc_x[k] * sfi1_236[k];

        t_237[k] = f_11 * sfh_114[k]
                   + f_3 * pc_y[k] * sgh_177[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, pb_y, pc_x, pc_y, sfi0_154, sfh_183, \
                         sfh_184, sfh_185, sfi1_154, sgh_183, sgh_184, \
                         sgh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = pb_y[k] * sfi0_154[k]
                   - f_10 * pc_y[k] * sfi1_154[k];

        t_239[k] = f_11 * sfh_183[k]
                   + f_3 * pc_x[k] * sgh_183[k];

        t_240[k] = f_11 * sfh_184[k]
                   + f_3 * pc_x[k] * sgh_184[k];

        t_241[k] = f_11 * sfh_185[k]
                   + f_3 * pc_x[k] * sgh_185[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pb_x, pc_x, sfi0_245, sfh_186, sfh_187, \
                         sfh_188, sfi1_245, sgh_186, sgh_187, sgh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_11 * sfh_186[k]
                   + f_3 * pc_x[k] * sgh_186[k];

        t_243[k] = f_11 * sfh_187[k]
                   + f_3 * pc_x[k] * sgh_187[k];

        t_244[k] = f_11 * sfh_188[k]
                   + f_3 * pc_x[k] * sgh_188[k];

        t_245[k] = pb_x[k] * sfi0_245[k]
                   - f_10 * pc_x[k] * sfi1_245[k];
    }
}

static auto
compute_prim_sgi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfi0,
                                                          const size_t sfh, const size_t sfi1,
                                                          const size_t sgg0, const size_t sgg1,
                                                          const size_t sgh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfi0_168 = buffer.data(sfi0 + 168);
    const auto *sfi0_171 = buffer.data(sfi0 + 171);
    const auto *sfi0_174 = buffer.data(sfi0 + 174);
    const auto *sfi0_178 = buffer.data(sfi0 + 178);
    const auto *sfi0_189 = buffer.data(sfi0 + 189);
    const auto *sfi0_191 = buffer.data(sfi0 + 191);
    const auto *sfi0_192 = buffer.data(sfi0 + 192);
    const auto *sfi0_193 = buffer.data(sfi0 + 193);
    const auto *sfi0_247 = buffer.data(sfi0 + 247);
    const auto *sfi0_248 = buffer.data(sfi0 + 248);
    const auto *sfi0_249 = buffer.data(sfi0 + 249);
    const auto *sfi0_251 = buffer.data(sfi0 + 251);
    const auto *sfi0_252 = buffer.data(sfi0 + 252);
    const auto *sfi0_255 = buffer.data(sfi0 + 255);
    const auto *sfi0_257 = buffer.data(sfi0 + 257);
    const auto *sfi0_258 = buffer.data(sfi0 + 258);
    const auto *sfi0_261 = buffer.data(sfi0 + 261);
    const auto *sfi0_262 = buffer.data(sfi0 + 262);
    const auto *sfi0_264 = buffer.data(sfi0 + 264);
    const auto *sfi0_266 = buffer.data(sfi0 + 266);
    const auto *sfi0_273 = buffer.data(sfi0 + 273);
    const auto *sfi0_275 = buffer.data(sfi0 + 275);
    const auto *sfi0_276 = buffer.data(sfi0 + 276);
    const auto *sfi0_277 = buffer.data(sfi0 + 277);
    const auto *sfi0_279 = buffer.data(sfi0 + 279);

    const auto *sfh_99 = buffer.data(sfh + 99);
    const auto *sfh_105 = buffer.data(sfh + 105);
    const auto *sfh_108 = buffer.data(sfh + 108);
    const auto *sfh_111 = buffer.data(sfh + 111);
    const auto *sfh_120 = buffer.data(sfh + 120);
    const auto *sfh_125 = buffer.data(sfh + 125);
    const auto *sfh_126 = buffer.data(sfh + 126);
    const auto *sfh_128 = buffer.data(sfh + 128);
    const auto *sfh_129 = buffer.data(sfh + 129);
    const auto *sfh_131 = buffer.data(sfh + 131);
    const auto *sfh_132 = buffer.data(sfh + 132);
    const auto *sfh_135 = buffer.data(sfh + 135);
    const auto *sfh_141 = buffer.data(sfh + 141);
    const auto *sfh_142 = buffer.data(sfh + 142);
    const auto *sfh_143 = buffer.data(sfh + 143);
    const auto *sfh_144 = buffer.data(sfh + 144);
    const auto *sfh_145 = buffer.data(sfh + 145);
    const auto *sfh_146 = buffer.data(sfh + 146);
    const auto *sfh_147 = buffer.data(sfh + 147);
    const auto *sfh_149 = buffer.data(sfh + 149);
    const auto *sfh_150 = buffer.data(sfh + 150);
    const auto *sfh_152 = buffer.data(sfh + 152);
    const auto *sfh_153 = buffer.data(sfh + 153);
    const auto *sfh_156 = buffer.data(sfh + 156);
    const auto *sfh_162 = buffer.data(sfh + 162);
    const auto *sfh_167 = buffer.data(sfh + 167);
    const auto *sfh_168 = buffer.data(sfh + 168);
    const auto *sfh_170 = buffer.data(sfh + 170);
    const auto *sfh_171 = buffer.data(sfh + 171);
    const auto *sfh_173 = buffer.data(sfh + 173);
    const auto *sfh_177 = buffer.data(sfh + 177);
    const auto *sfh_183 = buffer.data(sfh + 183);
    const auto *sfh_185 = buffer.data(sfh + 185);
    const auto *sfh_186 = buffer.data(sfh + 186);
    const auto *sfh_187 = buffer.data(sfh + 187);
    const auto *sfh_188 = buffer.data(sfh + 188);
    const auto *sfh_189 = buffer.data(sfh + 189);
    const auto *sfh_191 = buffer.data(sfh + 191);
    const auto *sfh_192 = buffer.data(sfh + 192);
    const auto *sfh_194 = buffer.data(sfh + 194);
    const auto *sfh_195 = buffer.data(sfh + 195);
    const auto *sfh_198 = buffer.data(sfh + 198);
    const auto *sfh_199 = buffer.data(sfh + 199);
    const auto *sfh_201 = buffer.data(sfh + 201);
    const auto *sfh_203 = buffer.data(sfh + 203);
    const auto *sfh_204 = buffer.data(sfh + 204);
    const auto *sfh_205 = buffer.data(sfh + 205);
    const auto *sfh_206 = buffer.data(sfh + 206);
    const auto *sfh_207 = buffer.data(sfh + 207);
    const auto *sfh_208 = buffer.data(sfh + 208);
    const auto *sfh_209 = buffer.data(sfh + 209);

    const auto *sfi1_168 = buffer.data(sfi1 + 168);
    const auto *sfi1_171 = buffer.data(sfi1 + 171);
    const auto *sfi1_174 = buffer.data(sfi1 + 174);
    const auto *sfi1_178 = buffer.data(sfi1 + 178);
    const auto *sfi1_189 = buffer.data(sfi1 + 189);
    const auto *sfi1_191 = buffer.data(sfi1 + 191);
    const auto *sfi1_192 = buffer.data(sfi1 + 192);
    const auto *sfi1_193 = buffer.data(sfi1 + 193);
    const auto *sfi1_247 = buffer.data(sfi1 + 247);
    const auto *sfi1_248 = buffer.data(sfi1 + 248);
    const auto *sfi1_249 = buffer.data(sfi1 + 249);
    const auto *sfi1_251 = buffer.data(sfi1 + 251);
    const auto *sfi1_252 = buffer.data(sfi1 + 252);
    const auto *sfi1_255 = buffer.data(sfi1 + 255);
    const auto *sfi1_257 = buffer.data(sfi1 + 257);
    const auto *sfi1_258 = buffer.data(sfi1 + 258);
    const auto *sfi1_261 = buffer.data(sfi1 + 261);
    const auto *sfi1_262 = buffer.data(sfi1 + 262);
    const auto *sfi1_264 = buffer.data(sfi1 + 264);
    const auto *sfi1_266 = buffer.data(sfi1 + 266);
    const auto *sfi1_273 = buffer.data(sfi1 + 273);
    const auto *sfi1_275 = buffer.data(sfi1 + 275);
    const auto *sfi1_276 = buffer.data(sfi1 + 276);
    const auto *sfi1_277 = buffer.data(sfi1 + 277);
    const auto *sfi1_279 = buffer.data(sfi1 + 279);

    const auto *sgg0_150 = buffer.data(sgg0 + 150);
    const auto *sgg0_153 = buffer.data(sgg0 + 153);
    const auto *sgg0_155 = buffer.data(sgg0 + 155);
    const auto *sgg0_156 = buffer.data(sgg0 + 156);
    const auto *sgg0_159 = buffer.data(sgg0 + 159);
    const auto *sgg0_160 = buffer.data(sgg0 + 160);
    const auto *sgg0_162 = buffer.data(sgg0 + 162);
    const auto *sgg0_163 = buffer.data(sgg0 + 163);
    const auto *sgg0_164 = buffer.data(sgg0 + 164);
    const auto *sgg0_170 = buffer.data(sgg0 + 170);
    const auto *sgg0_174 = buffer.data(sgg0 + 174);
    const auto *sgg0_177 = buffer.data(sgg0 + 177);
    const auto *sgg0_179 = buffer.data(sgg0 + 179);
    const auto *sgg0_180 = buffer.data(sgg0 + 180);
    const auto *sgg0_183 = buffer.data(sgg0 + 183);
    const auto *sgg0_185 = buffer.data(sgg0 + 185);
    const auto *sgg0_186 = buffer.data(sgg0 + 186);
    const auto *sgg0_189 = buffer.data(sgg0 + 189);
    const auto *sgg0_190 = buffer.data(sgg0 + 190);
    const auto *sgg0_192 = buffer.data(sgg0 + 192);
    const auto *sgg0_193 = buffer.data(sgg0 + 193);
    const auto *sgg0_194 = buffer.data(sgg0 + 194);
    const auto *sgg0_198 = buffer.data(sgg0 + 198);
    const auto *sgg0_201 = buffer.data(sgg0 + 201);

    const auto *sgg1_150 = buffer.data(sgg1 + 150);
    const auto *sgg1_153 = buffer.data(sgg1 + 153);
    const auto *sgg1_155 = buffer.data(sgg1 + 155);
    const auto *sgg1_156 = buffer.data(sgg1 + 156);
    const auto *sgg1_159 = buffer.data(sgg1 + 159);
    const auto *sgg1_160 = buffer.data(sgg1 + 160);
    const auto *sgg1_162 = buffer.data(sgg1 + 162);
    const auto *sgg1_163 = buffer.data(sgg1 + 163);
    const auto *sgg1_164 = buffer.data(sgg1 + 164);
    const auto *sgg1_170 = buffer.data(sgg1 + 170);
    const auto *sgg1_174 = buffer.data(sgg1 + 174);
    const auto *sgg1_177 = buffer.data(sgg1 + 177);
    const auto *sgg1_179 = buffer.data(sgg1 + 179);
    const auto *sgg1_180 = buffer.data(sgg1 + 180);
    const auto *sgg1_183 = buffer.data(sgg1 + 183);
    const auto *sgg1_185 = buffer.data(sgg1 + 185);
    const auto *sgg1_186 = buffer.data(sgg1 + 186);
    const auto *sgg1_189 = buffer.data(sgg1 + 189);
    const auto *sgg1_190 = buffer.data(sgg1 + 190);
    const auto *sgg1_192 = buffer.data(sgg1 + 192);
    const auto *sgg1_193 = buffer.data(sgg1 + 193);
    const auto *sgg1_194 = buffer.data(sgg1 + 194);
    const auto *sgg1_198 = buffer.data(sgg1 + 198);
    const auto *sgg1_201 = buffer.data(sgg1 + 201);

    const auto *sgh_183 = buffer.data(sgh + 183);
    const auto *sgh_188 = buffer.data(sgh + 188);
    const auto *sgh_189 = buffer.data(sgh + 189);
    const auto *sgh_191 = buffer.data(sgh + 191);
    const auto *sgh_192 = buffer.data(sgh + 192);
    const auto *sgh_194 = buffer.data(sgh + 194);
    const auto *sgh_195 = buffer.data(sgh + 195);
    const auto *sgh_198 = buffer.data(sgh + 198);
    const auto *sgh_204 = buffer.data(sgh + 204);
    const auto *sgh_205 = buffer.data(sgh + 205);
    const auto *sgh_206 = buffer.data(sgh + 206);
    const auto *sgh_207 = buffer.data(sgh + 207);
    const auto *sgh_208 = buffer.data(sgh + 208);
    const auto *sgh_209 = buffer.data(sgh + 209);
    const auto *sgh_210 = buffer.data(sgh + 210);
    const auto *sgh_212 = buffer.data(sgh + 212);
    const auto *sgh_213 = buffer.data(sgh + 213);
    const auto *sgh_215 = buffer.data(sgh + 215);
    const auto *sgh_216 = buffer.data(sgh + 216);
    const auto *sgh_219 = buffer.data(sgh + 219);
    const auto *sgh_220 = buffer.data(sgh + 220);
    const auto *sgh_222 = buffer.data(sgh + 222);
    const auto *sgh_224 = buffer.data(sgh + 224);
    const auto *sgh_225 = buffer.data(sgh + 225);
    const auto *sgh_226 = buffer.data(sgh + 226);
    const auto *sgh_227 = buffer.data(sgh + 227);
    const auto *sgh_228 = buffer.data(sgh + 228);
    const auto *sgh_229 = buffer.data(sgh + 229);
    const auto *sgh_230 = buffer.data(sgh + 230);
    const auto *sgh_231 = buffer.data(sgh + 231);
    const auto *sgh_233 = buffer.data(sgh + 233);
    const auto *sgh_234 = buffer.data(sgh + 234);
    const auto *sgh_236 = buffer.data(sgh + 236);
    const auto *sgh_237 = buffer.data(sgh + 237);
    const auto *sgh_240 = buffer.data(sgh + 240);
    const auto *sgh_243 = buffer.data(sgh + 243);
    const auto *sgh_245 = buffer.data(sgh + 245);
    const auto *sgh_246 = buffer.data(sgh + 246);
    const auto *sgh_247 = buffer.data(sgh + 247);
    const auto *sgh_248 = buffer.data(sgh + 248);
    const auto *sgh_249 = buffer.data(sgh + 249);
    const auto *sgh_250 = buffer.data(sgh + 250);
    const auto *sgh_251 = buffer.data(sgh + 251);
    const auto *sgh_252 = buffer.data(sgh + 252);
    const auto *sgh_254 = buffer.data(sgh + 254);
    const auto *sgh_255 = buffer.data(sgh + 255);
    const auto *sgh_257 = buffer.data(sgh + 257);
    const auto *sgh_258 = buffer.data(sgh + 258);
    const auto *sgh_261 = buffer.data(sgh + 261);
    const auto *sgh_262 = buffer.data(sgh + 262);
    const auto *sgh_264 = buffer.data(sgh + 264);
    const auto *sgh_266 = buffer.data(sgh + 266);
    const auto *sgh_267 = buffer.data(sgh + 267);
    const auto *sgh_268 = buffer.data(sgh + 268);
    const auto *sgh_269 = buffer.data(sgh + 269);
    const auto *sgh_270 = buffer.data(sgh + 270);
    const auto *sgh_271 = buffer.data(sgh + 271);
    const auto *sgh_272 = buffer.data(sgh + 272);
    const auto *sgh_273 = buffer.data(sgh + 273);
    const auto *sgh_275 = buffer.data(sgh + 275);
    const auto *sgh_276 = buffer.data(sgh + 276);
    const auto *sgh_278 = buffer.data(sgh + 278);
    const auto *sgh_279 = buffer.data(sgh + 279);

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pb_x, pc_x, pc_z, sfi0_247, sfi0_248, \
                         sfi0_249, sfh_99, sfi1_247, sfi1_248, sfi1_249, \
                         sgh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * sfh_99[k]
                   + f_3 * pc_z[k] * sgh_183[k];

        t_247[k] = pb_x[k] * sfi0_247[k]
                   - f_10 * pc_x[k] * sfi1_247[k];

        t_248[k] = pb_x[k] * sfi0_248[k]
                   - f_10 * pc_x[k] * sfi1_248[k];

        t_249[k] = pb_x[k] * sfi0_249[k]
                   - f_10 * pc_x[k] * sfi1_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pb_x, pc_x, pc_y, sfi0_251, sfi0_252, \
                         sfh_125, sfh_189, sfi1_251, sfi1_252, sgh_188, \
                         sgh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_11 * sfh_125[k]
                   + f_3 * pc_y[k] * sgh_188[k];

        t_251[k] = pb_x[k] * sfi0_251[k]
                   - f_10 * pc_x[k] * sfi1_251[k];

        t_252[k] = pb_x[k] * sfi0_252[k]
                   + f_14 * sfh_189[k]
                   - f_10 * pc_x[k] * sfi1_252[k];

        t_253[k] = f_3 * pc_y[k] * sgh_189[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pb_x, pc_x, pc_y, pc_z, sfi0_255, sfh_105, \
                         sfh_192, sfi1_255, sgh_189, sgh_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * sfh_105[k]
                   + f_3 * pc_z[k] * sgh_189[k];

        t_255[k] = pb_x[k] * sfi0_255[k]
                   + f_0 * sfh_192[k]
                   - f_10 * pc_x[k] * sfi1_255[k];

        t_256[k] = f_3 * pc_y[k] * sgh_191[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pb_x, pc_x, pc_z, sfi0_257, sfi0_258, sfh_108, \
                         sfh_194, sfh_195, sfi1_257, sfi1_258, \
                         sgh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = pb_x[k] * sfi0_257[k]
                   + f_0 * sfh_194[k]
                   - f_10 * pc_x[k] * sfi1_257[k];

        t_258[k] = pb_x[k] * sfi0_258[k]
                   + f_13 * sfh_195[k]
                   - f_10 * pc_x[k] * sfi1_258[k];

        t_259[k] = f_13 * sfh_108[k]
                   + f_3 * pc_z[k] * sgh_192[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pb_x, pc_x, pc_y, sfi0_261, sfi0_262, sfh_198, \
                         sfh_199, sfi1_261, sfi1_262, sgh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_3 * pc_y[k] * sgh_194[k];

        t_261[k] = pb_x[k] * sfi0_261[k]
                   + f_13 * sfh_198[k]
                   - f_10 * pc_x[k] * sfi1_261[k];

        t_262[k] = pb_x[k] * sfi0_262[k]
                   + f_12 * sfh_199[k]
                   - f_10 * pc_x[k] * sfi1_262[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pb_x, pc_x, pc_y, pc_z, sfi0_264, sfh_111, \
                         sfh_201, sfi1_264, sgh_195, sgh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_13 * sfh_111[k]
                   + f_3 * pc_z[k] * sgh_195[k];

        t_264[k] = pb_x[k] * sfi0_264[k]
                   + f_12 * sfh_201[k]
                   - f_10 * pc_x[k] * sfi1_264[k];

        t_265[k] = f_3 * pc_y[k] * sgh_198[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, pb_x, pc_x, sfi0_266, sfh_203, sfh_204, \
                         sfh_205, sfh_206, sfi1_266, sgh_204, sgh_205, \
                         sgh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = pb_x[k] * sfi0_266[k]
                   + f_12 * sfh_203[k]
                   - f_10 * pc_x[k] * sfi1_266[k];

        t_267[k] = f_11 * sfh_204[k]
                   + f_3 * pc_x[k] * sgh_204[k];

        t_268[k] = f_11 * sfh_205[k]
                   + f_3 * pc_x[k] * sgh_205[k];

        t_269[k] = f_11 * sfh_206[k]
                   + f_3 * pc_x[k] * sgh_206[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pb_x, pc_x, sfi0_273, sfh_207, sfh_208, \
                         sfh_209, sfi1_273, sgh_207, sgh_208, sgh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_11 * sfh_207[k]
                   + f_3 * pc_x[k] * sgh_207[k];

        t_271[k] = f_11 * sfh_208[k]
                   + f_3 * pc_x[k] * sgh_208[k];

        t_272[k] = f_11 * sfh_209[k]
                   + f_3 * pc_x[k] * sgh_209[k];

        t_273[k] = pb_x[k] * sfi0_273[k]
                   - f_10 * pc_x[k] * sfi1_273[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, pb_x, pc_x, pc_z, sfi0_275, sfi0_276, \
                         sfi0_277, sfh_120, sfi1_275, sfi1_276, sfi1_277, \
                         sgh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_13 * sfh_120[k]
                   + f_3 * pc_z[k] * sgh_204[k];

        t_275[k] = pb_x[k] * sfi0_275[k]
                   - f_10 * pc_x[k] * sfi1_275[k];

        t_276[k] = pb_x[k] * sfi0_276[k]
                   - f_10 * pc_x[k] * sfi1_276[k];

        t_277[k] = pb_x[k] * sfi0_277[k]
                   - f_10 * pc_x[k] * sfi1_277[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, pb_x, pc_x, pc_y, pc_z, sfi0_279, \
                         sfh_126, sfi1_279, sgg0_150, sgg1_150, sgh_209, \
                         sgh_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_3 * pc_y[k] * sgh_209[k];

        t_279[k] = pb_x[k] * sfi0_279[k]
                   - f_10 * pc_x[k] * sfi1_279[k];

        t_280[k] = f_1 * sgg0_150[k]
                   - f_2 * sgg1_150[k]
                   + f_3 * pc_x[k] * sgh_210[k];

        t_281[k] = f_0 * sfh_126[k]
                   + f_3 * pc_y[k] * sgh_210[k];

        t_282[k] = f_3 * pc_z[k] * sgh_210[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_x, pc_y, sfh_128, sgg0_153, sgg0_155, \
                         sgg1_153, sgg1_155, sgh_212, sgh_213, \
                         sgh_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_4 * sgg0_153[k]
                   - f_5 * sgg1_153[k]
                   + f_3 * pc_x[k] * sgh_213[k];

        t_284[k] = f_0 * sfh_128[k]
                   + f_3 * pc_y[k] * sgh_212[k];

        t_285[k] = f_4 * sgg0_155[k]
                   - f_5 * sgg1_155[k]
                   + f_3 * pc_x[k] * sgh_215[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pc_x, pc_y, pc_z, sfh_131, sgg0_156, \
                         sgg0_159, sgg1_156, sgg1_159, sgh_213, sgh_215, sgh_216, \
                         sgh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_6 * sgg0_156[k]
                   - f_7 * sgg1_156[k]
                   + f_3 * pc_x[k] * sgh_216[k];

        t_287[k] = f_3 * pc_z[k] * sgh_213[k];

        t_288[k] = f_0 * sfh_131[k]
                   + f_3 * pc_y[k] * sgh_215[k];

        t_289[k] = f_6 * sgg0_159[k]
                   - f_7 * sgg1_159[k]
                   + f_3 * pc_x[k] * sgh_219[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pc_x, pc_y, pc_z, sfh_135, sgg0_160, \
                         sgg0_162, sgg1_160, sgg1_162, sgh_216, sgh_219, sgh_220, \
                         sgh_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_8 * sgg0_160[k]
                   - f_9 * sgg1_160[k]
                   + f_3 * pc_x[k] * sgh_220[k];

        t_291[k] = f_3 * pc_z[k] * sgh_216[k];

        t_292[k] = f_8 * sgg0_162[k]
                   - f_9 * sgg1_162[k]
                   + f_3 * pc_x[k] * sgh_222[k];

        t_293[k] = f_0 * sfh_135[k]
                   + f_3 * pc_y[k] * sgh_219[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, t_299, pc_x, sgg0_164, sgg1_164, \
                         sgh_224, sgh_225, sgh_226, sgh_227, sgh_228, \
                         sgh_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_8 * sgg0_164[k]
                   - f_9 * sgg1_164[k]
                   + f_3 * pc_x[k] * sgh_224[k];

        t_295[k] = f_3 * pc_x[k] * sgh_225[k];

        t_296[k] = f_3 * pc_x[k] * sgh_226[k];

        t_297[k] = f_3 * pc_x[k] * sgh_227[k];

        t_298[k] = f_3 * pc_x[k] * sgh_228[k];

        t_299[k] = f_3 * pc_x[k] * sgh_229[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pc_x, pc_y, pc_z, sfh_141, sfh_143, \
                         sgg0_160, sgg0_162, sgg1_160, sgg1_162, sgh_225, sgh_227, \
                         sgh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_3 * pc_x[k] * sgh_230[k];

        t_301[k] = f_0 * sfh_141[k]
                   + f_1 * sgg0_160[k]
                   - f_2 * sgg1_160[k]
                   + f_3 * pc_y[k] * sgh_225[k];

        t_302[k] = f_3 * pc_z[k] * sgh_225[k];

        t_303[k] = f_0 * sfh_143[k]
                   + f_4 * sgg0_162[k]
                   - f_5 * sgg1_162[k]
                   + f_3 * pc_y[k] * sgh_227[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_y, pc_z, sfh_144, sfh_145, sfh_146, \
                         sgg0_163, sgg0_164, sgg1_163, sgg1_164, sgh_228, sgh_229, \
                         sgh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_0 * sfh_144[k]
                   + f_6 * sgg0_163[k]
                   - f_7 * sgg1_163[k]
                   + f_3 * pc_y[k] * sgh_228[k];

        t_305[k] = f_0 * sfh_145[k]
                   + f_8 * sgg0_164[k]
                   - f_9 * sgg1_164[k]
                   + f_3 * pc_y[k] * sgh_229[k];

        t_306[k] = f_0 * sfh_146[k]
                   + f_3 * pc_y[k] * sgh_230[k];

        t_307[k] = f_1 * sgg0_164[k]
                   - f_2 * sgg1_164[k]
                   + f_3 * pc_z[k] * sgh_230[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pb_z, pc_y, pc_z, sfi0_168, sfi0_171, \
                         sfh_126, sfh_147, sfi1_168, sfi1_171, \
                         sgh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pb_z[k] * sfi0_168[k]
                   - f_10 * pc_z[k] * sfi1_168[k];

        t_309[k] = f_13 * sfh_147[k]
                   + f_3 * pc_y[k] * sgh_231[k];

        t_310[k] = f_11 * sfh_126[k]
                   + f_3 * pc_z[k] * sgh_231[k];

        t_311[k] = pb_z[k] * sfi0_171[k]
                   - f_10 * pc_z[k] * sfi1_171[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pb_z, pc_x, pc_y, pc_z, sfi0_174, sfh_149, \
                         sfi1_174, sgg0_170, sgg1_170, sgh_233, \
                         sgh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_13 * sfh_149[k]
                   + f_3 * pc_y[k] * sgh_233[k];

        t_313[k] = f_4 * sgg0_170[k]
                   - f_5 * sgg1_170[k]
                   + f_3 * pc_x[k] * sgh_236[k];

        t_314[k] = pb_z[k] * sfi0_174[k]
                   - f_10 * pc_z[k] * sfi1_174[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pc_x, pc_y, pc_z, sfh_129, sfh_152, sgg0_174, \
                         sgg1_174, sgh_234, sgh_236, sgh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_11 * sfh_129[k]
                   + f_3 * pc_z[k] * sgh_234[k];

        t_316[k] = f_13 * sfh_152[k]
                   + f_3 * pc_y[k] * sgh_236[k];

        t_317[k] = f_6 * sgg0_174[k]
                   - f_7 * sgg1_174[k]
                   + f_3 * pc_x[k] * sgh_240[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pb_z, pc_x, pc_z, sfi0_178, sfh_132, sfi1_178, \
                         sgg0_177, sgg1_177, sgh_237, sgh_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pb_z[k] * sfi0_178[k]
                   - f_10 * pc_z[k] * sfi1_178[k];

        t_319[k] = f_11 * sfh_132[k]
                   + f_3 * pc_z[k] * sgh_237[k];

        t_320[k] = f_8 * sgg0_177[k]
                   - f_9 * sgg1_177[k]
                   + f_3 * pc_x[k] * sgh_243[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, t_325, pc_x, pc_y, sfh_156, sgg0_179, \
                         sgg1_179, sgh_240, sgh_245, sgh_246, sgh_247, \
                         sgh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * sfh_156[k]
                   + f_3 * pc_y[k] * sgh_240[k];

        t_322[k] = f_8 * sgg0_179[k]
                   - f_9 * sgg1_179[k]
                   + f_3 * pc_x[k] * sgh_245[k];

        t_323[k] = f_3 * pc_x[k] * sgh_246[k];

        t_324[k] = f_3 * pc_x[k] * sgh_247[k];

        t_325[k] = f_3 * pc_x[k] * sgh_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, pb_z, pc_x, pc_z, sfi0_189, \
                         sfh_141, sfi1_189, sgh_246, sgh_249, sgh_250, \
                         sgh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_3 * pc_x[k] * sgh_249[k];

        t_327[k] = f_3 * pc_x[k] * sgh_250[k];

        t_328[k] = f_3 * pc_x[k] * sgh_251[k];

        t_329[k] = pb_z[k] * sfi0_189[k]
                   - f_10 * pc_z[k] * sfi1_189[k];

        t_330[k] = f_11 * sfh_141[k]
                   + f_3 * pc_z[k] * sgh_246[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pb_z, pc_z, sfi0_191, sfi0_192, sfi0_193, \
                         sfh_142, sfh_143, sfh_144, sfi1_191, sfi1_192, \
                         sfi1_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = pb_z[k] * sfi0_191[k]
                   + f_12 * sfh_142[k]
                   - f_10 * pc_z[k] * sfi1_191[k];

        t_332[k] = pb_z[k] * sfi0_192[k]
                   + f_13 * sfh_143[k]
                   - f_10 * pc_z[k] * sfi1_192[k];

        t_333[k] = pb_z[k] * sfi0_193[k]
                   + f_0 * sfh_144[k]
                   - f_10 * pc_z[k] * sfi1_193[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pc_x, pc_y, pc_z, sfh_146, sfh_167, \
                         sfh_168, sgg0_179, sgg0_180, sgg1_179, sgg1_180, sgh_251, \
                         sgh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_13 * sfh_167[k]
                   + f_3 * pc_y[k] * sgh_251[k];

        t_335[k] = f_11 * sfh_146[k]
                   + f_1 * sgg0_179[k]
                   - f_2 * sgg1_179[k]
                   + f_3 * pc_z[k] * sgh_251[k];

        t_336[k] = f_1 * sgg0_180[k]
                   - f_2 * sgg1_180[k]
                   + f_3 * pc_x[k] * sgh_252[k];

        t_337[k] = f_12 * sfh_168[k]
                   + f_3 * pc_y[k] * sgh_252[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pc_x, pc_y, pc_z, sfh_147, sfh_170, sgg0_183, \
                         sgg1_183, sgh_252, sgh_254, sgh_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_12 * sfh_147[k]
                   + f_3 * pc_z[k] * sgh_252[k];

        t_339[k] = f_4 * sgg0_183[k]
                   - f_5 * sgg1_183[k]
                   + f_3 * pc_x[k] * sgh_255[k];

        t_340[k] = f_12 * sfh_170[k]
                   + f_3 * pc_y[k] * sgh_254[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_y, pc_z, sfh_150, sfh_173, \
                         sgg0_185, sgg0_186, sgg1_185, sgg1_186, sgh_255, sgh_257, \
                         sgh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_4 * sgg0_185[k]
                   - f_5 * sgg1_185[k]
                   + f_3 * pc_x[k] * sgh_257[k];

        t_342[k] = f_6 * sgg0_186[k]
                   - f_7 * sgg1_186[k]
                   + f_3 * pc_x[k] * sgh_258[k];

        t_343[k] = f_12 * sfh_150[k]
                   + f_3 * pc_z[k] * sgh_255[k];

        t_344[k] = f_12 * sfh_173[k]
                   + f_3 * pc_y[k] * sgh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, sfh_153, sgg0_189, sgg0_190, \
                         sgg1_189, sgg1_190, sgh_258, sgh_261, \
                         sgh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_6 * sgg0_189[k]
                   - f_7 * sgg1_189[k]
                   + f_3 * pc_x[k] * sgh_261[k];

        t_346[k] = f_8 * sgg0_190[k]
                   - f_9 * sgg1_190[k]
                   + f_3 * pc_x[k] * sgh_262[k];

        t_347[k] = f_12 * sfh_153[k]
                   + f_3 * pc_z[k] * sgh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pc_x, pc_y, sfh_177, sgg0_192, sgg0_194, \
                         sgg1_192, sgg1_194, sgh_261, sgh_264, sgh_266, \
                         sgh_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_8 * sgg0_192[k]
                   - f_9 * sgg1_192[k]
                   + f_3 * pc_x[k] * sgh_264[k];

        t_349[k] = f_12 * sfh_177[k]
                   + f_3 * pc_y[k] * sgh_261[k];

        t_350[k] = f_8 * sgg0_194[k]
                   - f_9 * sgg1_194[k]
                   + f_3 * pc_x[k] * sgh_266[k];

        t_351[k] = f_3 * pc_x[k] * sgh_267[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, pc_x, sgh_268, sgh_269, sgh_270, \
                         sgh_271, sgh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_3 * pc_x[k] * sgh_268[k];

        t_353[k] = f_3 * pc_x[k] * sgh_269[k];

        t_354[k] = f_3 * pc_x[k] * sgh_270[k];

        t_355[k] = f_3 * pc_x[k] * sgh_271[k];

        t_356[k] = f_3 * pc_x[k] * sgh_272[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pc_y, pc_z, sfh_162, sfh_183, sfh_185, sgg0_190, \
                         sgg0_192, sgg1_190, sgg1_192, sgh_267, \
                         sgh_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_12 * sfh_183[k]
                   + f_1 * sgg0_190[k]
                   - f_2 * sgg1_190[k]
                   + f_3 * pc_y[k] * sgh_267[k];

        t_358[k] = f_12 * sfh_162[k]
                   + f_3 * pc_z[k] * sgh_267[k];

        t_359[k] = f_12 * sfh_185[k]
                   + f_4 * sgg0_192[k]
                   - f_5 * sgg1_192[k]
                   + f_3 * pc_y[k] * sgh_269[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pc_y, sfh_186, sfh_187, sfh_188, sgg0_193, \
                         sgg0_194, sgg1_193, sgg1_194, sgh_270, sgh_271, \
                         sgh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_12 * sfh_186[k]
                   + f_6 * sgg0_193[k]
                   - f_7 * sgg1_193[k]
                   + f_3 * pc_y[k] * sgh_270[k];

        t_361[k] = f_12 * sfh_187[k]
                   + f_8 * sgg0_194[k]
                   - f_9 * sgg1_194[k]
                   + f_3 * pc_y[k] * sgh_271[k];

        t_362[k] = f_12 * sfh_188[k]
                   + f_3 * pc_y[k] * sgh_272[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pb_y, pc_y, pc_z, sfi0_252, sfh_167, \
                         sfh_168, sfh_189, sfi1_252, sgg0_194, sgg1_194, sgh_272, \
                         sgh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_12 * sfh_167[k]
                   + f_1 * sgg0_194[k]
                   - f_2 * sgg1_194[k]
                   + f_3 * pc_z[k] * sgh_272[k];

        t_364[k] = pb_y[k] * sfi0_252[k]
                   - f_10 * pc_y[k] * sfi1_252[k];

        t_365[k] = f_11 * sfh_189[k]
                   + f_3 * pc_y[k] * sgh_273[k];

        t_366[k] = f_13 * sfh_168[k]
                   + f_3 * pc_z[k] * sgh_273[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, pb_y, pc_x, pc_y, sfi0_257, sfh_191, sfi1_257, \
                         sgg0_198, sgg1_198, sgh_275, sgh_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_4 * sgg0_198[k]
                   - f_5 * sgg1_198[k]
                   + f_3 * pc_x[k] * sgh_276[k];

        t_368[k] = f_11 * sfh_191[k]
                   + f_3 * pc_y[k] * sgh_275[k];

        t_369[k] = pb_y[k] * sfi0_257[k]
                   - f_10 * pc_y[k] * sfi1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pc_x, pc_y, pc_z, sfh_171, sfh_194, sgg0_201, \
                         sgg1_201, sgh_276, sgh_278, sgh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_6 * sgg0_201[k]
                   - f_7 * sgg1_201[k]
                   + f_3 * pc_x[k] * sgh_279[k];

        t_371[k] = f_13 * sfh_171[k]
                   + f_3 * pc_z[k] * sgh_276[k];

        t_372[k] = f_11 * sfh_194[k]
                   + f_3 * pc_y[k] * sgh_278[k];
    }
}

static auto
compute_prim_sgi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfi0,
                                                          const size_t sfh, const size_t sfi1,
                                                          const size_t sgg0, const size_t sgg1,
                                                          const size_t sgh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfi0_261 = buffer.data(sfi0 + 261);
    const auto *sfi0_266 = buffer.data(sfi0 + 266);
    const auto *sfi0_273 = buffer.data(sfi0 + 273);
    const auto *sfi0_275 = buffer.data(sfi0 + 275);
    const auto *sfi0_276 = buffer.data(sfi0 + 276);
    const auto *sfi0_277 = buffer.data(sfi0 + 277);
    const auto *sfi0_279 = buffer.data(sfi0 + 279);

    const auto *sfh_174 = buffer.data(sfh + 174);
    const auto *sfh_183 = buffer.data(sfh + 183);
    const auto *sfh_189 = buffer.data(sfh + 189);
    const auto *sfh_192 = buffer.data(sfh + 192);
    const auto *sfh_195 = buffer.data(sfh + 195);
    const auto *sfh_198 = buffer.data(sfh + 198);
    const auto *sfh_204 = buffer.data(sfh + 204);
    const auto *sfh_206 = buffer.data(sfh + 206);
    const auto *sfh_207 = buffer.data(sfh + 207);
    const auto *sfh_208 = buffer.data(sfh + 208);
    const auto *sfh_209 = buffer.data(sfh + 209);

    const auto *sfi1_261 = buffer.data(sfi1 + 261);
    const auto *sfi1_266 = buffer.data(sfi1 + 266);
    const auto *sfi1_273 = buffer.data(sfi1 + 273);
    const auto *sfi1_275 = buffer.data(sfi1 + 275);
    const auto *sfi1_276 = buffer.data(sfi1 + 276);
    const auto *sfi1_277 = buffer.data(sfi1 + 277);
    const auto *sfi1_279 = buffer.data(sfi1 + 279);

    const auto *sgg0_205 = buffer.data(sgg0 + 205);
    const auto *sgg0_207 = buffer.data(sgg0 + 207);
    const auto *sgg0_210 = buffer.data(sgg0 + 210);
    const auto *sgg0_213 = buffer.data(sgg0 + 213);
    const auto *sgg0_215 = buffer.data(sgg0 + 215);
    const auto *sgg0_216 = buffer.data(sgg0 + 216);
    const auto *sgg0_219 = buffer.data(sgg0 + 219);
    const auto *sgg0_220 = buffer.data(sgg0 + 220);
    const auto *sgg0_222 = buffer.data(sgg0 + 222);
    const auto *sgg0_223 = buffer.data(sgg0 + 223);
    const auto *sgg0_224 = buffer.data(sgg0 + 224);

    const auto *sgg1_205 = buffer.data(sgg1 + 205);
    const auto *sgg1_207 = buffer.data(sgg1 + 207);
    const auto *sgg1_210 = buffer.data(sgg1 + 210);
    const auto *sgg1_213 = buffer.data(sgg1 + 213);
    const auto *sgg1_215 = buffer.data(sgg1 + 215);
    const auto *sgg1_216 = buffer.data(sgg1 + 216);
    const auto *sgg1_219 = buffer.data(sgg1 + 219);
    const auto *sgg1_220 = buffer.data(sgg1 + 220);
    const auto *sgg1_222 = buffer.data(sgg1 + 222);
    const auto *sgg1_223 = buffer.data(sgg1 + 223);
    const auto *sgg1_224 = buffer.data(sgg1 + 224);

    const auto *sgh_279 = buffer.data(sgh + 279);
    const auto *sgh_282 = buffer.data(sgh + 282);
    const auto *sgh_283 = buffer.data(sgh + 283);
    const auto *sgh_285 = buffer.data(sgh + 285);
    const auto *sgh_288 = buffer.data(sgh + 288);
    const auto *sgh_289 = buffer.data(sgh + 289);
    const auto *sgh_290 = buffer.data(sgh + 290);
    const auto *sgh_291 = buffer.data(sgh + 291);
    const auto *sgh_292 = buffer.data(sgh + 292);
    const auto *sgh_293 = buffer.data(sgh + 293);
    const auto *sgh_294 = buffer.data(sgh + 294);
    const auto *sgh_296 = buffer.data(sgh + 296);
    const auto *sgh_297 = buffer.data(sgh + 297);
    const auto *sgh_299 = buffer.data(sgh + 299);
    const auto *sgh_300 = buffer.data(sgh + 300);
    const auto *sgh_303 = buffer.data(sgh + 303);
    const auto *sgh_304 = buffer.data(sgh + 304);
    const auto *sgh_306 = buffer.data(sgh + 306);
    const auto *sgh_308 = buffer.data(sgh + 308);
    const auto *sgh_309 = buffer.data(sgh + 309);
    const auto *sgh_310 = buffer.data(sgh + 310);
    const auto *sgh_311 = buffer.data(sgh + 311);
    const auto *sgh_312 = buffer.data(sgh + 312);
    const auto *sgh_313 = buffer.data(sgh + 313);
    const auto *sgh_314 = buffer.data(sgh + 314);

#pragma omp simd aligned(t_373, t_374, t_375, pb_y, pc_x, pc_y, pc_z, sfi0_261, sfh_174, \
                         sfi1_261, sgg0_205, sgg1_205, sgh_279, \
                         sgh_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = pb_y[k] * sfi0_261[k]
                   - f_10 * pc_y[k] * sfi1_261[k];

        t_374[k] = f_8 * sgg0_205[k]
                   - f_9 * sgg1_205[k]
                   + f_3 * pc_x[k] * sgh_283[k];

        t_375[k] = f_13 * sfh_174[k]
                   + f_3 * pc_z[k] * sgh_279[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pb_y, pc_x, pc_y, sfi0_266, sfh_198, \
                         sfi1_266, sgg0_207, sgg1_207, sgh_282, sgh_285, \
                         sgh_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_8 * sgg0_207[k]
                   - f_9 * sgg1_207[k]
                   + f_3 * pc_x[k] * sgh_285[k];

        t_377[k] = f_11 * sfh_198[k]
                   + f_3 * pc_y[k] * sgh_282[k];

        t_378[k] = pb_y[k] * sfi0_266[k]
                   - f_10 * pc_y[k] * sfi1_266[k];

        t_379[k] = f_3 * pc_x[k] * sgh_288[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, pc_x, sgh_289, sgh_290, sgh_291, \
                         sgh_292, sgh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_3 * pc_x[k] * sgh_289[k];

        t_381[k] = f_3 * pc_x[k] * sgh_290[k];

        t_382[k] = f_3 * pc_x[k] * sgh_291[k];

        t_383[k] = f_3 * pc_x[k] * sgh_292[k];

        t_384[k] = f_3 * pc_x[k] * sgh_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pb_y, pc_y, pc_z, sfi0_273, sfi0_275, sfh_183, \
                         sfh_204, sfh_206, sfi1_273, sfi1_275, \
                         sgh_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pb_y[k] * sfi0_273[k]
                   + f_14 * sfh_204[k]
                   - f_10 * pc_y[k] * sfi1_273[k];

        t_386[k] = f_13 * sfh_183[k]
                   + f_3 * pc_z[k] * sgh_288[k];

        t_387[k] = pb_y[k] * sfi0_275[k]
                   + f_0 * sfh_206[k]
                   - f_10 * pc_y[k] * sfi1_275[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pb_y, pc_y, sfi0_276, sfi0_277, sfi0_279, \
                         sfh_207, sfh_208, sfh_209, sfi1_276, sfi1_277, sfi1_279, \
                         sgh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = pb_y[k] * sfi0_276[k]
                   + f_13 * sfh_207[k]
                   - f_10 * pc_y[k] * sfi1_276[k];

        t_389[k] = pb_y[k] * sfi0_277[k]
                   + f_12 * sfh_208[k]
                   - f_10 * pc_y[k] * sfi1_277[k];

        t_390[k] = f_11 * sfh_209[k]
                   + f_3 * pc_y[k] * sgh_293[k];

        t_391[k] = pb_y[k] * sfi0_279[k]
                   - f_10 * pc_y[k] * sfi1_279[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, pc_x, pc_y, pc_z, sfh_189, \
                         sgg0_210, sgg0_213, sgg1_210, sgg1_213, sgh_294, sgh_296, \
                         sgh_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_1 * sgg0_210[k]
                   - f_2 * sgg1_210[k]
                   + f_3 * pc_x[k] * sgh_294[k];

        t_393[k] = f_3 * pc_y[k] * sgh_294[k];

        t_394[k] = f_0 * sfh_189[k]
                   + f_3 * pc_z[k] * sgh_294[k];

        t_395[k] = f_4 * sgg0_213[k]
                   - f_5 * sgg1_213[k]
                   + f_3 * pc_x[k] * sgh_297[k];

        t_396[k] = f_3 * pc_y[k] * sgh_296[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pc_x, pc_y, pc_z, sfh_192, sgg0_215, \
                         sgg0_216, sgg1_215, sgg1_216, sgh_297, sgh_299, \
                         sgh_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_4 * sgg0_215[k]
                   - f_5 * sgg1_215[k]
                   + f_3 * pc_x[k] * sgh_299[k];

        t_398[k] = f_6 * sgg0_216[k]
                   - f_7 * sgg1_216[k]
                   + f_3 * pc_x[k] * sgh_300[k];

        t_399[k] = f_0 * sfh_192[k]
                   + f_3 * pc_z[k] * sgh_297[k];

        t_400[k] = f_3 * pc_y[k] * sgh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_z, sfh_195, sgg0_219, sgg0_220, \
                         sgg1_219, sgg1_220, sgh_300, sgh_303, \
                         sgh_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_6 * sgg0_219[k]
                   - f_7 * sgg1_219[k]
                   + f_3 * pc_x[k] * sgh_303[k];

        t_402[k] = f_8 * sgg0_220[k]
                   - f_9 * sgg1_220[k]
                   + f_3 * pc_x[k] * sgh_304[k];

        t_403[k] = f_0 * sfh_195[k]
                   + f_3 * pc_z[k] * sgh_300[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, t_408, pc_x, pc_y, sgg0_222, sgg0_224, \
                         sgg1_222, sgg1_224, sgh_303, sgh_306, sgh_308, sgh_309, \
                         sgh_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_8 * sgg0_222[k]
                   - f_9 * sgg1_222[k]
                   + f_3 * pc_x[k] * sgh_306[k];

        t_405[k] = f_3 * pc_y[k] * sgh_303[k];

        t_406[k] = f_8 * sgg0_224[k]
                   - f_9 * sgg1_224[k]
                   + f_3 * pc_x[k] * sgh_308[k];

        t_407[k] = f_3 * pc_x[k] * sgh_309[k];

        t_408[k] = f_3 * pc_x[k] * sgh_310[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, pc_y, sgg0_220, sgg1_220, \
                         sgh_309, sgh_311, sgh_312, sgh_313, sgh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_3 * pc_x[k] * sgh_311[k];

        t_410[k] = f_3 * pc_x[k] * sgh_312[k];

        t_411[k] = f_3 * pc_x[k] * sgh_313[k];

        t_412[k] = f_3 * pc_x[k] * sgh_314[k];

        t_413[k] = f_1 * sgg0_220[k]
                   - f_2 * sgg1_220[k]
                   + f_3 * pc_y[k] * sgh_309[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, pc_z, sfh_204, sgg0_222, sgg0_223, \
                         sgg1_222, sgg1_223, sgh_309, sgh_311, \
                         sgh_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_0 * sfh_204[k]
                   + f_3 * pc_z[k] * sgh_309[k];

        t_415[k] = f_4 * sgg0_222[k]
                   - f_5 * sgg1_222[k]
                   + f_3 * pc_y[k] * sgh_311[k];

        t_416[k] = f_6 * sgg0_223[k]
                   - f_7 * sgg1_223[k]
                   + f_3 * pc_y[k] * sgh_312[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_y, pc_z, sfh_209, sgg0_224, sgg1_224, \
                         sgh_313, sgh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_8 * sgg0_224[k]
                   - f_9 * sgg1_224[k]
                   + f_3 * pc_y[k] * sgh_313[k];

        t_418[k] = f_3 * pc_y[k] * sgh_314[k];

        t_419[k] = f_0 * sfh_209[k]
                   + f_1 * sgg0_224[k]
                   - f_2 * sgg1_224[k]
                   + f_3 * pc_z[k] * sgh_314[k];
    }
}

auto
compute_prim_sgi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sfi0, const size_t sfh,
                                                   const size_t sfi1, const size_t sgg0,
                                                   const size_t sgg1, const size_t sgh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sgi_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sfi0, sfh,
                                                              sfi1, sgg0, sgg1, sgh, ncols,
                                                              gamma, p, q);

    compute_prim_sgi_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sfi0, sfh,
                                                              sfi1, sgg0, sgg1, sgh, ncols,
                                                              gamma, p, q);

    compute_prim_sgi_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sfi0, sfh,
                                                              sfi1, sgg0, sgg1, sgh, ncols,
                                                              gamma, p, q);

    compute_prim_sgi_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sfi0, sfh,
                                                              sfi1, sgg0, sgg1, sgh, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
