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


#include "SimdThreeCenterElectronRepulsionVrrRecSHI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_shi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgi0,
                                                          const size_t sgh, const size_t sgi1,
                                                          const size_t shg0, const size_t shg1,
                                                          const size_t shh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_14 = 2.0 / q;

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

    const auto *sgi0_0 = buffer.data(sgi0 + 0);
    const auto *sgi0_3 = buffer.data(sgi0 + 3);
    const auto *sgi0_5 = buffer.data(sgi0 + 5);
    const auto *sgi0_6 = buffer.data(sgi0 + 6);
    const auto *sgi0_9 = buffer.data(sgi0 + 9);
    const auto *sgi0_10 = buffer.data(sgi0 + 10);
    const auto *sgi0_12 = buffer.data(sgi0 + 12);
    const auto *sgi0_14 = buffer.data(sgi0 + 14);
    const auto *sgi0_21 = buffer.data(sgi0 + 21);
    const auto *sgi0_27 = buffer.data(sgi0 + 27);
    const auto *sgi0_31 = buffer.data(sgi0 + 31);
    const auto *sgi0_34 = buffer.data(sgi0 + 34);
    const auto *sgi0_38 = buffer.data(sgi0 + 38);
    const auto *sgi0_56 = buffer.data(sgi0 + 56);
    const auto *sgi0_61 = buffer.data(sgi0 + 61);
    const auto *sgi0_65 = buffer.data(sgi0 + 65);

    const auto *sgh_0 = buffer.data(sgh + 0);
    const auto *sgh_1 = buffer.data(sgh + 1);
    const auto *sgh_2 = buffer.data(sgh + 2);
    const auto *sgh_3 = buffer.data(sgh + 3);
    const auto *sgh_5 = buffer.data(sgh + 5);
    const auto *sgh_6 = buffer.data(sgh + 6);
    const auto *sgh_7 = buffer.data(sgh + 7);
    const auto *sgh_8 = buffer.data(sgh + 8);
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
    const auto *sgh_47 = buffer.data(sgh + 47);
    const auto *sgh_57 = buffer.data(sgh + 57);
    const auto *sgh_58 = buffer.data(sgh + 58);
    const auto *sgh_59 = buffer.data(sgh + 59);
    const auto *sgh_60 = buffer.data(sgh + 60);
    const auto *sgh_61 = buffer.data(sgh + 61);
    const auto *sgh_62 = buffer.data(sgh + 62);
    const auto *sgh_63 = buffer.data(sgh + 63);
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

    const auto *sgi1_0 = buffer.data(sgi1 + 0);
    const auto *sgi1_3 = buffer.data(sgi1 + 3);
    const auto *sgi1_5 = buffer.data(sgi1 + 5);
    const auto *sgi1_6 = buffer.data(sgi1 + 6);
    const auto *sgi1_9 = buffer.data(sgi1 + 9);
    const auto *sgi1_10 = buffer.data(sgi1 + 10);
    const auto *sgi1_12 = buffer.data(sgi1 + 12);
    const auto *sgi1_14 = buffer.data(sgi1 + 14);
    const auto *sgi1_21 = buffer.data(sgi1 + 21);
    const auto *sgi1_27 = buffer.data(sgi1 + 27);
    const auto *sgi1_31 = buffer.data(sgi1 + 31);
    const auto *sgi1_34 = buffer.data(sgi1 + 34);
    const auto *sgi1_38 = buffer.data(sgi1 + 38);
    const auto *sgi1_56 = buffer.data(sgi1 + 56);
    const auto *sgi1_61 = buffer.data(sgi1 + 61);
    const auto *sgi1_65 = buffer.data(sgi1 + 65);

    const auto *shg0_0 = buffer.data(shg0 + 0);
    const auto *shg0_3 = buffer.data(shg0 + 3);
    const auto *shg0_5 = buffer.data(shg0 + 5);
    const auto *shg0_6 = buffer.data(shg0 + 6);
    const auto *shg0_9 = buffer.data(shg0 + 9);
    const auto *shg0_10 = buffer.data(shg0 + 10);
    const auto *shg0_12 = buffer.data(shg0 + 12);
    const auto *shg0_13 = buffer.data(shg0 + 13);
    const auto *shg0_14 = buffer.data(shg0 + 14);
    const auto *shg0_25 = buffer.data(shg0 + 25);
    const auto *shg0_27 = buffer.data(shg0 + 27);
    const auto *shg0_28 = buffer.data(shg0 + 28);
    const auto *shg0_29 = buffer.data(shg0 + 29);
    const auto *shg0_42 = buffer.data(shg0 + 42);
    const auto *shg0_43 = buffer.data(shg0 + 43);
    const auto *shg0_44 = buffer.data(shg0 + 44);
    const auto *shg0_45 = buffer.data(shg0 + 45);
    const auto *shg0_48 = buffer.data(shg0 + 48);
    const auto *shg0_50 = buffer.data(shg0 + 50);
    const auto *shg0_51 = buffer.data(shg0 + 51);
    const auto *shg0_54 = buffer.data(shg0 + 54);
    const auto *shg0_55 = buffer.data(shg0 + 55);
    const auto *shg0_57 = buffer.data(shg0 + 57);
    const auto *shg0_58 = buffer.data(shg0 + 58);
    const auto *shg0_59 = buffer.data(shg0 + 59);

    const auto *shg1_0 = buffer.data(shg1 + 0);
    const auto *shg1_3 = buffer.data(shg1 + 3);
    const auto *shg1_5 = buffer.data(shg1 + 5);
    const auto *shg1_6 = buffer.data(shg1 + 6);
    const auto *shg1_9 = buffer.data(shg1 + 9);
    const auto *shg1_10 = buffer.data(shg1 + 10);
    const auto *shg1_12 = buffer.data(shg1 + 12);
    const auto *shg1_13 = buffer.data(shg1 + 13);
    const auto *shg1_14 = buffer.data(shg1 + 14);
    const auto *shg1_25 = buffer.data(shg1 + 25);
    const auto *shg1_27 = buffer.data(shg1 + 27);
    const auto *shg1_28 = buffer.data(shg1 + 28);
    const auto *shg1_29 = buffer.data(shg1 + 29);
    const auto *shg1_42 = buffer.data(shg1 + 42);
    const auto *shg1_43 = buffer.data(shg1 + 43);
    const auto *shg1_44 = buffer.data(shg1 + 44);
    const auto *shg1_45 = buffer.data(shg1 + 45);
    const auto *shg1_48 = buffer.data(shg1 + 48);
    const auto *shg1_50 = buffer.data(shg1 + 50);
    const auto *shg1_51 = buffer.data(shg1 + 51);
    const auto *shg1_54 = buffer.data(shg1 + 54);
    const auto *shg1_55 = buffer.data(shg1 + 55);
    const auto *shg1_57 = buffer.data(shg1 + 57);
    const auto *shg1_58 = buffer.data(shg1 + 58);
    const auto *shg1_59 = buffer.data(shg1 + 59);

    const auto *shh_0 = buffer.data(shh + 0);
    const auto *shh_2 = buffer.data(shh + 2);
    const auto *shh_3 = buffer.data(shh + 3);
    const auto *shh_5 = buffer.data(shh + 5);
    const auto *shh_6 = buffer.data(shh + 6);
    const auto *shh_9 = buffer.data(shh + 9);
    const auto *shh_10 = buffer.data(shh + 10);
    const auto *shh_12 = buffer.data(shh + 12);
    const auto *shh_14 = buffer.data(shh + 14);
    const auto *shh_15 = buffer.data(shh + 15);
    const auto *shh_16 = buffer.data(shh + 16);
    const auto *shh_17 = buffer.data(shh + 17);
    const auto *shh_18 = buffer.data(shh + 18);
    const auto *shh_19 = buffer.data(shh + 19);
    const auto *shh_20 = buffer.data(shh + 20);
    const auto *shh_21 = buffer.data(shh + 21);
    const auto *shh_23 = buffer.data(shh + 23);
    const auto *shh_24 = buffer.data(shh + 24);
    const auto *shh_26 = buffer.data(shh + 26);
    const auto *shh_27 = buffer.data(shh + 27);
    const auto *shh_30 = buffer.data(shh + 30);
    const auto *shh_36 = buffer.data(shh + 36);
    const auto *shh_37 = buffer.data(shh + 37);
    const auto *shh_38 = buffer.data(shh + 38);
    const auto *shh_39 = buffer.data(shh + 39);
    const auto *shh_40 = buffer.data(shh + 40);
    const auto *shh_41 = buffer.data(shh + 41);
    const auto *shh_42 = buffer.data(shh + 42);
    const auto *shh_44 = buffer.data(shh + 44);
    const auto *shh_45 = buffer.data(shh + 45);
    const auto *shh_47 = buffer.data(shh + 47);
    const auto *shh_48 = buffer.data(shh + 48);
    const auto *shh_51 = buffer.data(shh + 51);
    const auto *shh_57 = buffer.data(shh + 57);
    const auto *shh_58 = buffer.data(shh + 58);
    const auto *shh_59 = buffer.data(shh + 59);
    const auto *shh_60 = buffer.data(shh + 60);
    const auto *shh_61 = buffer.data(shh + 61);
    const auto *shh_62 = buffer.data(shh + 62);
    const auto *shh_63 = buffer.data(shh + 63);
    const auto *shh_65 = buffer.data(shh + 65);
    const auto *shh_66 = buffer.data(shh + 66);
    const auto *shh_68 = buffer.data(shh + 68);
    const auto *shh_69 = buffer.data(shh + 69);
    const auto *shh_72 = buffer.data(shh + 72);
    const auto *shh_73 = buffer.data(shh + 73);
    const auto *shh_75 = buffer.data(shh + 75);
    const auto *shh_77 = buffer.data(shh + 77);
    const auto *shh_78 = buffer.data(shh + 78);
    const auto *shh_79 = buffer.data(shh + 79);
    const auto *shh_80 = buffer.data(shh + 80);
    const auto *shh_81 = buffer.data(shh + 81);
    const auto *shh_82 = buffer.data(shh + 82);
    const auto *shh_83 = buffer.data(shh + 83);
    const auto *shh_84 = buffer.data(shh + 84);
    const auto *shh_86 = buffer.data(shh + 86);
    const auto *shh_87 = buffer.data(shh + 87);
    const auto *shh_89 = buffer.data(shh + 89);
    const auto *shh_90 = buffer.data(shh + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sgh_0, sgh_3, shg0_0, shg0_3, \
                         shg1_0, shg1_3, shh_0, shh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sgh_0[k]
                 + f_1 * shg0_0[k]
                 - f_2 * shg1_0[k]
                 + f_3 * pc_x[k] * shh_0[k];

        t_1[k] = f_3 * pc_y[k] * shh_0[k];

        t_2[k] = f_3 * pc_z[k] * shh_0[k];

        t_3[k] = f_0 * sgh_3[k]
                 + f_4 * shg0_3[k]
                 - f_5 * shg1_3[k]
                 + f_3 * pc_x[k] * shh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sgh_5, sgh_6, shg0_5, shg0_6, shg1_5, \
                         shg1_6, shh_2, shh_5, shh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * shh_2[k];

        t_5[k] = f_0 * sgh_5[k]
                 + f_4 * shg0_5[k]
                 - f_5 * shg1_5[k]
                 + f_3 * pc_x[k] * shh_5[k];

        t_6[k] = f_0 * sgh_6[k]
                 + f_6 * shg0_6[k]
                 - f_7 * shg1_6[k]
                 + f_3 * pc_x[k] * shh_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sgh_9, shg0_9, shg1_9, shh_3, shh_5, \
                         shh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * shh_3[k];

        t_8[k] = f_3 * pc_y[k] * shh_5[k];

        t_9[k] = f_0 * sgh_9[k]
                 + f_6 * shg0_9[k]
                 - f_7 * shg1_9[k]
                 + f_3 * pc_x[k] * shh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sgh_10, sgh_12, shg0_10, shg0_12, \
                         shg1_10, shg1_12, shh_6, shh_10, shh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sgh_10[k]
                  + f_8 * shg0_10[k]
                  - f_9 * shg1_10[k]
                  + f_3 * pc_x[k] * shh_10[k];

        t_11[k] = f_3 * pc_z[k] * shh_6[k];

        t_12[k] = f_0 * sgh_12[k]
                  + f_8 * shg0_12[k]
                  - f_9 * shg1_12[k]
                  + f_3 * pc_x[k] * shh_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pc_x, pc_y, sgh_14, sgh_15, sgh_16, shg0_14, \
                         shg1_14, shh_9, shh_14, shh_15, shh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * shh_9[k];

        t_14[k] = f_0 * sgh_14[k]
                  + f_8 * shg0_14[k]
                  - f_9 * shg1_14[k]
                  + f_3 * pc_x[k] * shh_14[k];

        t_15[k] = f_0 * sgh_15[k]
                  + f_3 * pc_x[k] * shh_15[k];

        t_16[k] = f_0 * sgh_16[k]
                  + f_3 * pc_x[k] * shh_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_x, sgh_17, sgh_18, sgh_19, sgh_20, shh_17, \
                         shh_18, shh_19, shh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * sgh_17[k]
                  + f_3 * pc_x[k] * shh_17[k];

        t_18[k] = f_0 * sgh_18[k]
                  + f_3 * pc_x[k] * shh_18[k];

        t_19[k] = f_0 * sgh_19[k]
                  + f_3 * pc_x[k] * shh_19[k];

        t_20[k] = f_0 * sgh_20[k]
                  + f_3 * pc_x[k] * shh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, shg0_10, shg0_12, shg0_13, \
                         shg1_10, shg1_12, shg1_13, shh_15, shh_17, \
                         shh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * shg0_10[k]
                  - f_2 * shg1_10[k]
                  + f_3 * pc_y[k] * shh_15[k];

        t_22[k] = f_3 * pc_z[k] * shh_15[k];

        t_23[k] = f_4 * shg0_12[k]
                  - f_5 * shg1_12[k]
                  + f_3 * pc_y[k] * shh_17[k];

        t_24[k] = f_6 * shg0_13[k]
                  - f_7 * shg1_13[k]
                  + f_3 * pc_y[k] * shh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_y, pc_y, pc_z, sgi0_0, sgh_0, \
                         sgi1_0, shg0_14, shg1_14, shh_19, shh_20, \
                         shh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * shg0_14[k]
                  - f_9 * shg1_14[k]
                  + f_3 * pc_y[k] * shh_19[k];

        t_26[k] = f_3 * pc_y[k] * shh_20[k];

        t_27[k] = f_1 * shg0_14[k]
                  - f_2 * shg1_14[k]
                  + f_3 * pc_z[k] * shh_20[k];

        t_28[k] = pb_y[k] * sgi0_0[k]
                  - f_10 * pc_y[k] * sgi1_0[k];

        t_29[k] = f_11 * sgh_0[k]
                  + f_3 * pc_y[k] * shh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_y, pc_y, pc_z, sgi0_3, sgi0_5, sgh_1, \
                         sgh_2, sgi1_3, sgi1_5, shh_21, shh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * shh_21[k];

        t_31[k] = pb_y[k] * sgi0_3[k]
                  + f_12 * sgh_1[k]
                  - f_10 * pc_y[k] * sgi1_3[k];

        t_32[k] = f_11 * sgh_2[k]
                  + f_3 * pc_y[k] * shh_23[k];

        t_33[k] = pb_y[k] * sgi0_5[k]
                  - f_10 * pc_y[k] * sgi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_y, pc_y, pc_z, sgi0_6, sgi0_9, sgh_3, \
                         sgh_5, sgi1_6, sgi1_9, shh_24, shh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * sgi0_6[k]
                  + f_13 * sgh_3[k]
                  - f_10 * pc_y[k] * sgi1_6[k];

        t_35[k] = f_3 * pc_z[k] * shh_24[k];

        t_36[k] = f_11 * sgh_5[k]
                  + f_3 * pc_y[k] * shh_26[k];

        t_37[k] = pb_y[k] * sgi0_9[k]
                  - f_10 * pc_y[k] * sgi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, sgi0_10, sgi0_12, sgh_6, \
                         sgh_8, sgh_9, sgi1_10, sgi1_12, shh_27, \
                         shh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * sgi0_10[k]
                  + f_14 * sgh_6[k]
                  - f_10 * pc_y[k] * sgi1_10[k];

        t_39[k] = f_3 * pc_z[k] * shh_27[k];

        t_40[k] = pb_y[k] * sgi0_12[k]
                  + f_12 * sgh_8[k]
                  - f_10 * pc_y[k] * sgi1_12[k];

        t_41[k] = f_11 * sgh_9[k]
                  + f_3 * pc_y[k] * shh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_y, pc_x, pc_y, sgi0_14, sgh_36, sgh_37, \
                         sgh_38, sgi1_14, shh_36, shh_37, shh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * sgi0_14[k]
                  - f_10 * pc_y[k] * sgi1_14[k];

        t_43[k] = f_14 * sgh_36[k]
                  + f_3 * pc_x[k] * shh_36[k];

        t_44[k] = f_14 * sgh_37[k]
                  + f_3 * pc_x[k] * shh_37[k];

        t_45[k] = f_14 * sgh_38[k]
                  + f_3 * pc_x[k] * shh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, sgh_15, sgh_39, sgh_40, sgh_41, \
                         shg0_25, shg1_25, shh_36, shh_39, shh_40, \
                         shh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_14 * sgh_39[k]
                  + f_3 * pc_x[k] * shh_39[k];

        t_47[k] = f_14 * sgh_40[k]
                  + f_3 * pc_x[k] * shh_40[k];

        t_48[k] = f_14 * sgh_41[k]
                  + f_3 * pc_x[k] * shh_41[k];

        t_49[k] = f_11 * sgh_15[k]
                  + f_1 * shg0_25[k]
                  - f_2 * shg1_25[k]
                  + f_3 * pc_y[k] * shh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pc_y, pc_z, sgh_17, sgh_18, shg0_27, shg0_28, \
                         shg1_27, shg1_28, shh_36, shh_38, shh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * shh_36[k];

        t_51[k] = f_11 * sgh_17[k]
                  + f_4 * shg0_27[k]
                  - f_5 * shg1_27[k]
                  + f_3 * pc_y[k] * shh_38[k];

        t_52[k] = f_11 * sgh_18[k]
                  + f_6 * shg0_28[k]
                  - f_7 * shg1_28[k]
                  + f_3 * pc_y[k] * shh_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_y, sgi0_27, sgh_19, sgh_20, sgi1_27, \
                         shg0_29, shg1_29, shh_40, shh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_11 * sgh_19[k]
                  + f_8 * shg0_29[k]
                  - f_9 * shg1_29[k]
                  + f_3 * pc_y[k] * shh_40[k];

        t_54[k] = f_11 * sgh_20[k]
                  + f_3 * pc_y[k] * shh_41[k];

        t_55[k] = pb_y[k] * sgi0_27[k]
                  - f_10 * pc_y[k] * sgi1_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pb_z, pc_y, pc_z, sgi0_0, sgi0_3, \
                         sgh_0, sgi1_0, sgi1_3, shh_42, shh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * sgi0_0[k]
                  - f_10 * pc_z[k] * sgi1_0[k];

        t_57[k] = f_3 * pc_y[k] * shh_42[k];

        t_58[k] = f_11 * sgh_0[k]
                  + f_3 * pc_z[k] * shh_42[k];

        t_59[k] = pb_z[k] * sgi0_3[k]
                  - f_10 * pc_z[k] * sgi1_3[k];

        t_60[k] = f_3 * pc_y[k] * shh_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, sgi0_5, sgi0_6, sgh_2, \
                         sgh_3, sgi1_5, sgi1_6, shh_45, shh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_z[k] * sgi0_5[k]
                  + f_12 * sgh_2[k]
                  - f_10 * pc_z[k] * sgi1_5[k];

        t_62[k] = pb_z[k] * sgi0_6[k]
                  - f_10 * pc_z[k] * sgi1_6[k];

        t_63[k] = f_11 * sgh_3[k]
                  + f_3 * pc_z[k] * shh_45[k];

        t_64[k] = f_3 * pc_y[k] * shh_47[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_z, pc_z, sgi0_9, sgi0_10, sgi0_12, sgh_5, \
                         sgh_6, sgh_7, sgi1_9, sgi1_10, sgi1_12, \
                         shh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_z[k] * sgi0_9[k]
                  + f_13 * sgh_5[k]
                  - f_10 * pc_z[k] * sgi1_9[k];

        t_66[k] = pb_z[k] * sgi0_10[k]
                  - f_10 * pc_z[k] * sgi1_10[k];

        t_67[k] = f_11 * sgh_6[k]
                  + f_3 * pc_z[k] * shh_48[k];

        t_68[k] = pb_z[k] * sgi0_12[k]
                  + f_12 * sgh_7[k]
                  - f_10 * pc_z[k] * sgi1_12[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_z, pc_x, pc_y, pc_z, sgi0_14, sgh_9, \
                         sgh_57, sgh_58, sgi1_14, shh_51, shh_57, \
                         shh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * pc_y[k] * shh_51[k];

        t_70[k] = pb_z[k] * sgi0_14[k]
                  + f_14 * sgh_9[k]
                  - f_10 * pc_z[k] * sgi1_14[k];

        t_71[k] = f_14 * sgh_57[k]
                  + f_3 * pc_x[k] * shh_57[k];

        t_72[k] = f_14 * sgh_58[k]
                  + f_3 * pc_x[k] * shh_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pc_x, sgh_59, sgh_60, sgh_61, sgh_62, shh_59, \
                         shh_60, shh_61, shh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_14 * sgh_59[k]
                  + f_3 * pc_x[k] * shh_59[k];

        t_74[k] = f_14 * sgh_60[k]
                  + f_3 * pc_x[k] * shh_60[k];

        t_75[k] = f_14 * sgh_61[k]
                  + f_3 * pc_x[k] * shh_61[k];

        t_76[k] = f_14 * sgh_62[k]
                  + f_3 * pc_x[k] * shh_62[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_z, pc_y, pc_z, sgi0_21, sgh_15, sgi1_21, \
                         shg0_42, shg1_42, shh_57, shh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_z[k] * sgi0_21[k]
                  - f_10 * pc_z[k] * sgi1_21[k];

        t_78[k] = f_11 * sgh_15[k]
                  + f_3 * pc_z[k] * shh_57[k];

        t_79[k] = f_4 * shg0_42[k]
                  - f_5 * shg1_42[k]
                  + f_3 * pc_y[k] * shh_59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_y, pc_z, sgh_20, shg0_43, shg0_44, \
                         shg1_43, shg1_44, shh_60, shh_61, shh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_6 * shg0_43[k]
                  - f_7 * shg1_43[k]
                  + f_3 * pc_y[k] * shh_60[k];

        t_81[k] = f_8 * shg0_44[k]
                  - f_9 * shg1_44[k]
                  + f_3 * pc_y[k] * shh_61[k];

        t_82[k] = f_3 * pc_y[k] * shh_62[k];

        t_83[k] = f_11 * sgh_20[k]
                  + f_1 * shg0_44[k]
                  - f_2 * shg1_44[k]
                  + f_3 * pc_z[k] * shh_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pc_x, pc_y, pc_z, sgh_21, sgh_63, sgh_66, \
                         shg0_45, shg0_48, shg1_45, shg1_48, shh_63, \
                         shh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_13 * sgh_63[k]
                  + f_1 * shg0_45[k]
                  - f_2 * shg1_45[k]
                  + f_3 * pc_x[k] * shh_63[k];

        t_85[k] = f_12 * sgh_21[k]
                  + f_3 * pc_y[k] * shh_63[k];

        t_86[k] = f_3 * pc_z[k] * shh_63[k];

        t_87[k] = f_13 * sgh_66[k]
                  + f_4 * shg0_48[k]
                  - f_5 * shg1_48[k]
                  + f_3 * pc_x[k] * shh_66[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pc_x, pc_y, sgh_23, sgh_68, sgh_69, shg0_50, \
                         shg0_51, shg1_50, shg1_51, shh_65, shh_68, \
                         shh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_12 * sgh_23[k]
                  + f_3 * pc_y[k] * shh_65[k];

        t_89[k] = f_13 * sgh_68[k]
                  + f_4 * shg0_50[k]
                  - f_5 * shg1_50[k]
                  + f_3 * pc_x[k] * shh_68[k];

        t_90[k] = f_13 * sgh_69[k]
                  + f_6 * shg0_51[k]
                  - f_7 * shg1_51[k]
                  + f_3 * pc_x[k] * shh_69[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pc_x, pc_y, pc_z, sgh_26, sgh_72, shg0_54, shg1_54, \
                         shh_66, shh_68, shh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * pc_z[k] * shh_66[k];

        t_92[k] = f_12 * sgh_26[k]
                  + f_3 * pc_y[k] * shh_68[k];

        t_93[k] = f_13 * sgh_72[k]
                  + f_6 * shg0_54[k]
                  - f_7 * shg1_54[k]
                  + f_3 * pc_x[k] * shh_72[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pc_x, pc_z, sgh_73, sgh_75, shg0_55, shg0_57, \
                         shg1_55, shg1_57, shh_69, shh_73, shh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_13 * sgh_73[k]
                  + f_8 * shg0_55[k]
                  - f_9 * shg1_55[k]
                  + f_3 * pc_x[k] * shh_73[k];

        t_95[k] = f_3 * pc_z[k] * shh_69[k];

        t_96[k] = f_13 * sgh_75[k]
                  + f_8 * shg0_57[k]
                  - f_9 * shg1_57[k]
                  + f_3 * pc_x[k] * shh_75[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, sgh_30, sgh_77, sgh_78, sgh_79, \
                         shg0_59, shg1_59, shh_72, shh_77, shh_78, \
                         shh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_12 * sgh_30[k]
                  + f_3 * pc_y[k] * shh_72[k];

        t_98[k] = f_13 * sgh_77[k]
                  + f_8 * shg0_59[k]
                  - f_9 * shg1_59[k]
                  + f_3 * pc_x[k] * shh_77[k];

        t_99[k] = f_13 * sgh_78[k]
                  + f_3 * pc_x[k] * shh_78[k];

        t_100[k] = f_13 * sgh_79[k]
                   + f_3 * pc_x[k] * shh_79[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_x, sgh_80, sgh_81, sgh_82, sgh_83, \
                         shh_80, shh_81, shh_82, shh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_13 * sgh_80[k]
                   + f_3 * pc_x[k] * shh_80[k];

        t_102[k] = f_13 * sgh_81[k]
                   + f_3 * pc_x[k] * shh_81[k];

        t_103[k] = f_13 * sgh_82[k]
                   + f_3 * pc_x[k] * shh_82[k];

        t_104[k] = f_13 * sgh_83[k]
                   + f_3 * pc_x[k] * shh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pc_y, pc_z, sgh_36, sgh_38, shg0_55, shg0_57, \
                         shg1_55, shg1_57, shh_78, shh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * sgh_36[k]
                   + f_1 * shg0_55[k]
                   - f_2 * shg1_55[k]
                   + f_3 * pc_y[k] * shh_78[k];

        t_106[k] = f_3 * pc_z[k] * shh_78[k];

        t_107[k] = f_12 * sgh_38[k]
                   + f_4 * shg0_57[k]
                   - f_5 * shg1_57[k]
                   + f_3 * pc_y[k] * shh_80[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pc_y, pc_z, sgh_39, sgh_40, sgh_41, \
                         shg0_58, shg0_59, shg1_58, shg1_59, shh_81, shh_82, \
                         shh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_12 * sgh_39[k]
                   + f_6 * shg0_58[k]
                   - f_7 * shg1_58[k]
                   + f_3 * pc_y[k] * shh_81[k];

        t_109[k] = f_12 * sgh_40[k]
                   + f_8 * shg0_59[k]
                   - f_9 * shg1_59[k]
                   + f_3 * pc_y[k] * shh_82[k];

        t_110[k] = f_12 * sgh_41[k]
                   + f_3 * pc_y[k] * shh_83[k];

        t_111[k] = f_1 * shg0_59[k]
                   - f_2 * shg1_59[k]
                   + f_3 * pc_z[k] * shh_83[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pb_y, pb_z, pc_y, pc_z, sgi0_31, sgi0_56, \
                         sgh_21, sgh_42, sgi1_31, sgi1_56, shh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = pb_y[k] * sgi0_56[k]
                   - f_10 * pc_y[k] * sgi1_56[k];

        t_113[k] = f_11 * sgh_42[k]
                   + f_3 * pc_y[k] * shh_84[k];

        t_114[k] = f_11 * sgh_21[k]
                   + f_3 * pc_z[k] * shh_84[k];

        t_115[k] = pb_z[k] * sgi0_31[k]
                   - f_10 * pc_z[k] * sgi1_31[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pb_z, pc_y, pc_z, sgi0_34, sgi0_61, \
                         sgh_24, sgh_44, sgi1_34, sgi1_61, shh_86, \
                         shh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_11 * sgh_44[k]
                   + f_3 * pc_y[k] * shh_86[k];

        t_117[k] = pb_y[k] * sgi0_61[k]
                   - f_10 * pc_y[k] * sgi1_61[k];

        t_118[k] = pb_z[k] * sgi0_34[k]
                   - f_10 * pc_z[k] * sgi1_34[k];

        t_119[k] = f_11 * sgh_24[k]
                   + f_3 * pc_z[k] * shh_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pb_z, pc_y, pc_z, sgi0_38, sgi0_65, \
                         sgh_27, sgh_47, sgi1_38, sgi1_65, shh_89, \
                         shh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_11 * sgh_47[k]
                   + f_3 * pc_y[k] * shh_89[k];

        t_121[k] = pb_y[k] * sgi0_65[k]
                   - f_10 * pc_y[k] * sgi1_65[k];

        t_122[k] = pb_z[k] * sgi0_38[k]
                   - f_10 * pc_z[k] * sgi1_38[k];

        t_123[k] = f_11 * sgh_27[k]
                   + f_3 * pc_z[k] * shh_90[k];
    }
}

static auto
compute_prim_shi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgi0,
                                                          const size_t sgh, const size_t sgi1,
                                                          const size_t shg0, const size_t shg1,
                                                          const size_t shh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_14 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgi0_49 = buffer.data(sgi0 + 49);
    const auto *sgi0_68 = buffer.data(sgi0 + 68);
    const auto *sgi0_70 = buffer.data(sgi0 + 70);
    const auto *sgi0_83 = buffer.data(sgi0 + 83);
    const auto *sgi0_84 = buffer.data(sgi0 + 84);
    const auto *sgi0_87 = buffer.data(sgi0 + 87);
    const auto *sgi0_90 = buffer.data(sgi0 + 90);
    const auto *sgi0_94 = buffer.data(sgi0 + 94);
    const auto *sgi0_96 = buffer.data(sgi0 + 96);
    const auto *sgi0_105 = buffer.data(sgi0 + 105);
    const auto *sgi0_140 = buffer.data(sgi0 + 140);
    const auto *sgi0_143 = buffer.data(sgi0 + 143);
    const auto *sgi0_145 = buffer.data(sgi0 + 145);
    const auto *sgi0_146 = buffer.data(sgi0 + 146);
    const auto *sgi0_149 = buffer.data(sgi0 + 149);
    const auto *sgi0_150 = buffer.data(sgi0 + 150);
    const auto *sgi0_152 = buffer.data(sgi0 + 152);
    const auto *sgi0_154 = buffer.data(sgi0 + 154);

    const auto *sgh_36 = buffer.data(sgh + 36);
    const auto *sgh_42 = buffer.data(sgh + 42);
    const auto *sgh_45 = buffer.data(sgh + 45);
    const auto *sgh_48 = buffer.data(sgh + 48);
    const auto *sgh_50 = buffer.data(sgh + 50);
    const auto *sgh_51 = buffer.data(sgh + 51);
    const auto *sgh_57 = buffer.data(sgh + 57);
    const auto *sgh_59 = buffer.data(sgh + 59);
    const auto *sgh_60 = buffer.data(sgh + 60);
    const auto *sgh_61 = buffer.data(sgh + 61);
    const auto *sgh_62 = buffer.data(sgh + 62);
    const auto *sgh_63 = buffer.data(sgh + 63);
    const auto *sgh_65 = buffer.data(sgh + 65);
    const auto *sgh_66 = buffer.data(sgh + 66);
    const auto *sgh_68 = buffer.data(sgh + 68);
    const auto *sgh_69 = buffer.data(sgh + 69);
    const auto *sgh_70 = buffer.data(sgh + 70);
    const auto *sgh_72 = buffer.data(sgh + 72);
    const auto *sgh_78 = buffer.data(sgh + 78);
    const auto *sgh_80 = buffer.data(sgh + 80);
    const auto *sgh_81 = buffer.data(sgh + 81);
    const auto *sgh_82 = buffer.data(sgh + 82);
    const auto *sgh_83 = buffer.data(sgh + 83);
    const auto *sgh_84 = buffer.data(sgh + 84);
    const auto *sgh_86 = buffer.data(sgh + 86);
    const auto *sgh_87 = buffer.data(sgh + 87);
    const auto *sgh_89 = buffer.data(sgh + 89);
    const auto *sgh_90 = buffer.data(sgh + 90);
    const auto *sgh_93 = buffer.data(sgh + 93);
    const auto *sgh_99 = buffer.data(sgh + 99);
    const auto *sgh_100 = buffer.data(sgh + 100);
    const auto *sgh_101 = buffer.data(sgh + 101);
    const auto *sgh_102 = buffer.data(sgh + 102);
    const auto *sgh_103 = buffer.data(sgh + 103);
    const auto *sgh_104 = buffer.data(sgh + 104);
    const auto *sgh_105 = buffer.data(sgh + 105);
    const auto *sgh_106 = buffer.data(sgh + 106);
    const auto *sgh_107 = buffer.data(sgh + 107);
    const auto *sgh_108 = buffer.data(sgh + 108);
    const auto *sgh_110 = buffer.data(sgh + 110);
    const auto *sgh_111 = buffer.data(sgh + 111);
    const auto *sgh_113 = buffer.data(sgh + 113);
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
    const auto *sgh_129 = buffer.data(sgh + 129);
    const auto *sgh_131 = buffer.data(sgh + 131);
    const auto *sgh_132 = buffer.data(sgh + 132);
    const auto *sgh_135 = buffer.data(sgh + 135);
    const auto *sgh_136 = buffer.data(sgh + 136);
    const auto *sgh_138 = buffer.data(sgh + 138);
    const auto *sgh_140 = buffer.data(sgh + 140);
    const auto *sgh_141 = buffer.data(sgh + 141);
    const auto *sgh_142 = buffer.data(sgh + 142);
    const auto *sgh_143 = buffer.data(sgh + 143);
    const auto *sgh_144 = buffer.data(sgh + 144);
    const auto *sgh_145 = buffer.data(sgh + 145);
    const auto *sgh_146 = buffer.data(sgh + 146);
    const auto *sgh_152 = buffer.data(sgh + 152);
    const auto *sgh_156 = buffer.data(sgh + 156);
    const auto *sgh_161 = buffer.data(sgh + 161);
    const auto *sgh_162 = buffer.data(sgh + 162);
    const auto *sgh_163 = buffer.data(sgh + 163);
    const auto *sgh_164 = buffer.data(sgh + 164);
    const auto *sgh_165 = buffer.data(sgh + 165);
    const auto *sgh_166 = buffer.data(sgh + 166);
    const auto *sgh_167 = buffer.data(sgh + 167);
    const auto *sgh_183 = buffer.data(sgh + 183);

    const auto *sgi1_49 = buffer.data(sgi1 + 49);
    const auto *sgi1_68 = buffer.data(sgi1 + 68);
    const auto *sgi1_70 = buffer.data(sgi1 + 70);
    const auto *sgi1_83 = buffer.data(sgi1 + 83);
    const auto *sgi1_84 = buffer.data(sgi1 + 84);
    const auto *sgi1_87 = buffer.data(sgi1 + 87);
    const auto *sgi1_90 = buffer.data(sgi1 + 90);
    const auto *sgi1_94 = buffer.data(sgi1 + 94);
    const auto *sgi1_96 = buffer.data(sgi1 + 96);
    const auto *sgi1_105 = buffer.data(sgi1 + 105);
    const auto *sgi1_140 = buffer.data(sgi1 + 140);
    const auto *sgi1_143 = buffer.data(sgi1 + 143);
    const auto *sgi1_145 = buffer.data(sgi1 + 145);
    const auto *sgi1_146 = buffer.data(sgi1 + 146);
    const auto *sgi1_149 = buffer.data(sgi1 + 149);
    const auto *sgi1_150 = buffer.data(sgi1 + 150);
    const auto *sgi1_152 = buffer.data(sgi1 + 152);
    const auto *sgi1_154 = buffer.data(sgi1 + 154);

    const auto *shg0_72 = buffer.data(shg0 + 72);
    const auto *shg0_73 = buffer.data(shg0 + 73);
    const auto *shg0_74 = buffer.data(shg0 + 74);
    const auto *shg0_75 = buffer.data(shg0 + 75);
    const auto *shg0_78 = buffer.data(shg0 + 78);
    const auto *shg0_80 = buffer.data(shg0 + 80);
    const auto *shg0_81 = buffer.data(shg0 + 81);
    const auto *shg0_84 = buffer.data(shg0 + 84);
    const auto *shg0_85 = buffer.data(shg0 + 85);
    const auto *shg0_87 = buffer.data(shg0 + 87);
    const auto *shg0_88 = buffer.data(shg0 + 88);
    const auto *shg0_89 = buffer.data(shg0 + 89);
    const auto *shg0_90 = buffer.data(shg0 + 90);
    const auto *shg0_93 = buffer.data(shg0 + 93);
    const auto *shg0_95 = buffer.data(shg0 + 95);
    const auto *shg0_96 = buffer.data(shg0 + 96);
    const auto *shg0_99 = buffer.data(shg0 + 99);
    const auto *shg0_100 = buffer.data(shg0 + 100);
    const auto *shg0_102 = buffer.data(shg0 + 102);
    const auto *shg0_103 = buffer.data(shg0 + 103);
    const auto *shg0_104 = buffer.data(shg0 + 104);
    const auto *shg0_110 = buffer.data(shg0 + 110);
    const auto *shg0_114 = buffer.data(shg0 + 114);
    const auto *shg0_117 = buffer.data(shg0 + 117);
    const auto *shg0_118 = buffer.data(shg0 + 118);
    const auto *shg0_119 = buffer.data(shg0 + 119);

    const auto *shg1_72 = buffer.data(shg1 + 72);
    const auto *shg1_73 = buffer.data(shg1 + 73);
    const auto *shg1_74 = buffer.data(shg1 + 74);
    const auto *shg1_75 = buffer.data(shg1 + 75);
    const auto *shg1_78 = buffer.data(shg1 + 78);
    const auto *shg1_80 = buffer.data(shg1 + 80);
    const auto *shg1_81 = buffer.data(shg1 + 81);
    const auto *shg1_84 = buffer.data(shg1 + 84);
    const auto *shg1_85 = buffer.data(shg1 + 85);
    const auto *shg1_87 = buffer.data(shg1 + 87);
    const auto *shg1_88 = buffer.data(shg1 + 88);
    const auto *shg1_89 = buffer.data(shg1 + 89);
    const auto *shg1_90 = buffer.data(shg1 + 90);
    const auto *shg1_93 = buffer.data(shg1 + 93);
    const auto *shg1_95 = buffer.data(shg1 + 95);
    const auto *shg1_96 = buffer.data(shg1 + 96);
    const auto *shg1_99 = buffer.data(shg1 + 99);
    const auto *shg1_100 = buffer.data(shg1 + 100);
    const auto *shg1_102 = buffer.data(shg1 + 102);
    const auto *shg1_103 = buffer.data(shg1 + 103);
    const auto *shg1_104 = buffer.data(shg1 + 104);
    const auto *shg1_110 = buffer.data(shg1 + 110);
    const auto *shg1_114 = buffer.data(shg1 + 114);
    const auto *shg1_117 = buffer.data(shg1 + 117);
    const auto *shg1_118 = buffer.data(shg1 + 118);
    const auto *shg1_119 = buffer.data(shg1 + 119);

    const auto *shh_93 = buffer.data(shh + 93);
    const auto *shh_99 = buffer.data(shh + 99);
    const auto *shh_100 = buffer.data(shh + 100);
    const auto *shh_101 = buffer.data(shh + 101);
    const auto *shh_102 = buffer.data(shh + 102);
    const auto *shh_103 = buffer.data(shh + 103);
    const auto *shh_104 = buffer.data(shh + 104);
    const auto *shh_105 = buffer.data(shh + 105);
    const auto *shh_107 = buffer.data(shh + 107);
    const auto *shh_108 = buffer.data(shh + 108);
    const auto *shh_110 = buffer.data(shh + 110);
    const auto *shh_111 = buffer.data(shh + 111);
    const auto *shh_114 = buffer.data(shh + 114);
    const auto *shh_115 = buffer.data(shh + 115);
    const auto *shh_117 = buffer.data(shh + 117);
    const auto *shh_119 = buffer.data(shh + 119);
    const auto *shh_120 = buffer.data(shh + 120);
    const auto *shh_121 = buffer.data(shh + 121);
    const auto *shh_122 = buffer.data(shh + 122);
    const auto *shh_123 = buffer.data(shh + 123);
    const auto *shh_124 = buffer.data(shh + 124);
    const auto *shh_125 = buffer.data(shh + 125);
    const auto *shh_126 = buffer.data(shh + 126);
    const auto *shh_128 = buffer.data(shh + 128);
    const auto *shh_129 = buffer.data(shh + 129);
    const auto *shh_131 = buffer.data(shh + 131);
    const auto *shh_132 = buffer.data(shh + 132);
    const auto *shh_135 = buffer.data(shh + 135);
    const auto *shh_136 = buffer.data(shh + 136);
    const auto *shh_138 = buffer.data(shh + 138);
    const auto *shh_140 = buffer.data(shh + 140);
    const auto *shh_141 = buffer.data(shh + 141);
    const auto *shh_142 = buffer.data(shh + 142);
    const auto *shh_143 = buffer.data(shh + 143);
    const auto *shh_144 = buffer.data(shh + 144);
    const auto *shh_145 = buffer.data(shh + 145);
    const auto *shh_146 = buffer.data(shh + 146);
    const auto *shh_147 = buffer.data(shh + 147);
    const auto *shh_149 = buffer.data(shh + 149);
    const auto *shh_150 = buffer.data(shh + 150);
    const auto *shh_152 = buffer.data(shh + 152);
    const auto *shh_153 = buffer.data(shh + 153);
    const auto *shh_156 = buffer.data(shh + 156);
    const auto *shh_161 = buffer.data(shh + 161);
    const auto *shh_162 = buffer.data(shh + 162);
    const auto *shh_163 = buffer.data(shh + 163);
    const auto *shh_164 = buffer.data(shh + 164);
    const auto *shh_165 = buffer.data(shh + 165);
    const auto *shh_166 = buffer.data(shh + 166);
    const auto *shh_167 = buffer.data(shh + 167);
    const auto *shh_168 = buffer.data(shh + 168);
    const auto *shh_170 = buffer.data(shh + 170);
    const auto *shh_171 = buffer.data(shh + 171);
    const auto *shh_173 = buffer.data(shh + 173);
    const auto *shh_174 = buffer.data(shh + 174);
    const auto *shh_177 = buffer.data(shh + 177);
    const auto *shh_183 = buffer.data(shh + 183);

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_y, pc_x, pc_y, sgi0_68, sgi0_70, \
                         sgh_50, sgh_51, sgh_99, sgi1_68, sgi1_70, shh_93, \
                         shh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * sgi0_68[k]
                   + f_12 * sgh_50[k]
                   - f_10 * pc_y[k] * sgi1_68[k];

        t_125[k] = f_11 * sgh_51[k]
                   + f_3 * pc_y[k] * shh_93[k];

        t_126[k] = pb_y[k] * sgi0_70[k]
                   - f_10 * pc_y[k] * sgi1_70[k];

        t_127[k] = f_13 * sgh_99[k]
                   + f_3 * pc_x[k] * shh_99[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, sgh_100, sgh_101, sgh_102, \
                         sgh_103, sgh_104, shh_100, shh_101, shh_102, shh_103, \
                         shh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_13 * sgh_100[k]
                   + f_3 * pc_x[k] * shh_100[k];

        t_129[k] = f_13 * sgh_101[k]
                   + f_3 * pc_x[k] * shh_101[k];

        t_130[k] = f_13 * sgh_102[k]
                   + f_3 * pc_x[k] * shh_102[k];

        t_131[k] = f_13 * sgh_103[k]
                   + f_3 * pc_x[k] * shh_103[k];

        t_132[k] = f_13 * sgh_104[k]
                   + f_3 * pc_x[k] * shh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_z, pc_y, pc_z, sgi0_49, sgh_36, sgh_59, \
                         sgi1_49, shg0_72, shg1_72, shh_99, shh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pb_z[k] * sgi0_49[k]
                   - f_10 * pc_z[k] * sgi1_49[k];

        t_134[k] = f_11 * sgh_36[k]
                   + f_3 * pc_z[k] * shh_99[k];

        t_135[k] = f_11 * sgh_59[k]
                   + f_4 * shg0_72[k]
                   - f_5 * shg1_72[k]
                   + f_3 * pc_y[k] * shh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, sgh_60, sgh_61, sgh_62, shg0_73, shg0_74, \
                         shg1_73, shg1_74, shh_102, shh_103, shh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * sgh_60[k]
                   + f_6 * shg0_73[k]
                   - f_7 * shg1_73[k]
                   + f_3 * pc_y[k] * shh_102[k];

        t_137[k] = f_11 * sgh_61[k]
                   + f_8 * shg0_74[k]
                   - f_9 * shg1_74[k]
                   + f_3 * pc_y[k] * shh_103[k];

        t_138[k] = f_11 * sgh_62[k]
                   + f_3 * pc_y[k] * shh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pb_y, pc_x, pc_y, pc_z, sgi0_83, sgh_42, \
                         sgh_105, sgi1_83, shg0_75, shg1_75, shh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * sgi0_83[k]
                   - f_10 * pc_y[k] * sgi1_83[k];

        t_140[k] = f_13 * sgh_105[k]
                   + f_1 * shg0_75[k]
                   - f_2 * shg1_75[k]
                   + f_3 * pc_x[k] * shh_105[k];

        t_141[k] = f_3 * pc_y[k] * shh_105[k];

        t_142[k] = f_12 * sgh_42[k]
                   + f_3 * pc_z[k] * shh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, sgh_108, sgh_110, shg0_78, shg0_80, \
                         shg1_78, shg1_80, shh_107, shh_108, shh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_13 * sgh_108[k]
                   + f_4 * shg0_78[k]
                   - f_5 * shg1_78[k]
                   + f_3 * pc_x[k] * shh_108[k];

        t_144[k] = f_3 * pc_y[k] * shh_107[k];

        t_145[k] = f_13 * sgh_110[k]
                   + f_4 * shg0_80[k]
                   - f_5 * shg1_80[k]
                   + f_3 * pc_x[k] * shh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_x, pc_y, pc_z, sgh_45, sgh_111, shg0_81, \
                         shg1_81, shh_108, shh_110, shh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_13 * sgh_111[k]
                   + f_6 * shg0_81[k]
                   - f_7 * shg1_81[k]
                   + f_3 * pc_x[k] * shh_111[k];

        t_147[k] = f_12 * sgh_45[k]
                   + f_3 * pc_z[k] * shh_108[k];

        t_148[k] = f_3 * pc_y[k] * shh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_z, sgh_48, sgh_114, sgh_115, shg0_84, \
                         shg0_85, shg1_84, shg1_85, shh_111, shh_114, \
                         shh_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_13 * sgh_114[k]
                   + f_6 * shg0_84[k]
                   - f_7 * shg1_84[k]
                   + f_3 * pc_x[k] * shh_114[k];

        t_150[k] = f_13 * sgh_115[k]
                   + f_8 * shg0_85[k]
                   - f_9 * shg1_85[k]
                   + f_3 * pc_x[k] * shh_115[k];

        t_151[k] = f_12 * sgh_48[k]
                   + f_3 * pc_z[k] * shh_111[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pc_x, pc_y, sgh_117, sgh_119, shg0_87, shg0_89, \
                         shg1_87, shg1_89, shh_114, shh_117, shh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_13 * sgh_117[k]
                   + f_8 * shg0_87[k]
                   - f_9 * shg1_87[k]
                   + f_3 * pc_x[k] * shh_117[k];

        t_153[k] = f_3 * pc_y[k] * shh_114[k];

        t_154[k] = f_13 * sgh_119[k]
                   + f_8 * shg0_89[k]
                   - f_9 * shg1_89[k]
                   + f_3 * pc_x[k] * shh_119[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pc_x, sgh_120, sgh_121, sgh_122, \
                         sgh_123, sgh_124, shh_120, shh_121, shh_122, shh_123, \
                         shh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_13 * sgh_120[k]
                   + f_3 * pc_x[k] * shh_120[k];

        t_156[k] = f_13 * sgh_121[k]
                   + f_3 * pc_x[k] * shh_121[k];

        t_157[k] = f_13 * sgh_122[k]
                   + f_3 * pc_x[k] * shh_122[k];

        t_158[k] = f_13 * sgh_123[k]
                   + f_3 * pc_x[k] * shh_123[k];

        t_159[k] = f_13 * sgh_124[k]
                   + f_3 * pc_x[k] * shh_124[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, sgh_57, sgh_125, \
                         shg0_85, shg0_87, shg1_85, shg1_87, shh_120, shh_122, \
                         shh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_13 * sgh_125[k]
                   + f_3 * pc_x[k] * shh_125[k];

        t_161[k] = f_1 * shg0_85[k]
                   - f_2 * shg1_85[k]
                   + f_3 * pc_y[k] * shh_120[k];

        t_162[k] = f_12 * sgh_57[k]
                   + f_3 * pc_z[k] * shh_120[k];

        t_163[k] = f_4 * shg0_87[k]
                   - f_5 * shg1_87[k]
                   + f_3 * pc_y[k] * shh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, sgh_62, shg0_88, shg0_89, \
                         shg1_88, shg1_89, shh_123, shh_124, shh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * shg0_88[k]
                   - f_7 * shg1_88[k]
                   + f_3 * pc_y[k] * shh_123[k];

        t_165[k] = f_8 * shg0_89[k]
                   - f_9 * shg1_89[k]
                   + f_3 * pc_y[k] * shh_124[k];

        t_166[k] = f_3 * pc_y[k] * shh_125[k];

        t_167[k] = f_12 * sgh_62[k]
                   + f_1 * shg0_89[k]
                   - f_2 * shg1_89[k]
                   + f_3 * pc_z[k] * shh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, sgh_63, sgh_126, \
                         sgh_129, shg0_90, shg0_93, shg1_90, shg1_93, shh_126, \
                         shh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_12 * sgh_126[k]
                   + f_1 * shg0_90[k]
                   - f_2 * shg1_90[k]
                   + f_3 * pc_x[k] * shh_126[k];

        t_169[k] = f_13 * sgh_63[k]
                   + f_3 * pc_y[k] * shh_126[k];

        t_170[k] = f_3 * pc_z[k] * shh_126[k];

        t_171[k] = f_12 * sgh_129[k]
                   + f_4 * shg0_93[k]
                   - f_5 * shg1_93[k]
                   + f_3 * pc_x[k] * shh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pc_x, pc_y, sgh_65, sgh_131, sgh_132, shg0_95, \
                         shg0_96, shg1_95, shg1_96, shh_128, shh_131, \
                         shh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_13 * sgh_65[k]
                   + f_3 * pc_y[k] * shh_128[k];

        t_173[k] = f_12 * sgh_131[k]
                   + f_4 * shg0_95[k]
                   - f_5 * shg1_95[k]
                   + f_3 * pc_x[k] * shh_131[k];

        t_174[k] = f_12 * sgh_132[k]
                   + f_6 * shg0_96[k]
                   - f_7 * shg1_96[k]
                   + f_3 * pc_x[k] * shh_132[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pc_x, pc_y, pc_z, sgh_68, sgh_135, shg0_99, \
                         shg1_99, shh_129, shh_131, shh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_3 * pc_z[k] * shh_129[k];

        t_176[k] = f_13 * sgh_68[k]
                   + f_3 * pc_y[k] * shh_131[k];

        t_177[k] = f_12 * sgh_135[k]
                   + f_6 * shg0_99[k]
                   - f_7 * shg1_99[k]
                   + f_3 * pc_x[k] * shh_135[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_z, sgh_136, sgh_138, shg0_100, \
                         shg0_102, shg1_100, shg1_102, shh_132, shh_136, \
                         shh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_12 * sgh_136[k]
                   + f_8 * shg0_100[k]
                   - f_9 * shg1_100[k]
                   + f_3 * pc_x[k] * shh_136[k];

        t_179[k] = f_3 * pc_z[k] * shh_132[k];

        t_180[k] = f_12 * sgh_138[k]
                   + f_8 * shg0_102[k]
                   - f_9 * shg1_102[k]
                   + f_3 * pc_x[k] * shh_138[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, sgh_72, sgh_140, sgh_141, \
                         sgh_142, shg0_104, shg1_104, shh_135, shh_140, shh_141, \
                         shh_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_13 * sgh_72[k]
                   + f_3 * pc_y[k] * shh_135[k];

        t_182[k] = f_12 * sgh_140[k]
                   + f_8 * shg0_104[k]
                   - f_9 * shg1_104[k]
                   + f_3 * pc_x[k] * shh_140[k];

        t_183[k] = f_12 * sgh_141[k]
                   + f_3 * pc_x[k] * shh_141[k];

        t_184[k] = f_12 * sgh_142[k]
                   + f_3 * pc_x[k] * shh_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, sgh_143, sgh_144, sgh_145, sgh_146, \
                         shh_143, shh_144, shh_145, shh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_12 * sgh_143[k]
                   + f_3 * pc_x[k] * shh_143[k];

        t_186[k] = f_12 * sgh_144[k]
                   + f_3 * pc_x[k] * shh_144[k];

        t_187[k] = f_12 * sgh_145[k]
                   + f_3 * pc_x[k] * shh_145[k];

        t_188[k] = f_12 * sgh_146[k]
                   + f_3 * pc_x[k] * shh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_y, pc_z, sgh_78, sgh_80, shg0_100, shg0_102, \
                         shg1_100, shg1_102, shh_141, shh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * sgh_78[k]
                   + f_1 * shg0_100[k]
                   - f_2 * shg1_100[k]
                   + f_3 * pc_y[k] * shh_141[k];

        t_190[k] = f_3 * pc_z[k] * shh_141[k];

        t_191[k] = f_13 * sgh_80[k]
                   + f_4 * shg0_102[k]
                   - f_5 * shg1_102[k]
                   + f_3 * pc_y[k] * shh_143[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_y, pc_z, sgh_81, sgh_82, sgh_83, \
                         shg0_103, shg0_104, shg1_103, shg1_104, shh_144, shh_145, \
                         shh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_13 * sgh_81[k]
                   + f_6 * shg0_103[k]
                   - f_7 * shg1_103[k]
                   + f_3 * pc_y[k] * shh_144[k];

        t_193[k] = f_13 * sgh_82[k]
                   + f_8 * shg0_104[k]
                   - f_9 * shg1_104[k]
                   + f_3 * pc_y[k] * shh_145[k];

        t_194[k] = f_13 * sgh_83[k]
                   + f_3 * pc_y[k] * shh_146[k];

        t_195[k] = f_1 * shg0_104[k]
                   - f_2 * shg1_104[k]
                   + f_3 * pc_z[k] * shh_146[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pb_z, pc_y, pc_z, sgi0_84, sgi0_87, \
                         sgh_63, sgh_84, sgi1_84, sgi1_87, shh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pb_z[k] * sgi0_84[k]
                   - f_10 * pc_z[k] * sgi1_84[k];

        t_197[k] = f_12 * sgh_84[k]
                   + f_3 * pc_y[k] * shh_147[k];

        t_198[k] = f_11 * sgh_63[k]
                   + f_3 * pc_z[k] * shh_147[k];

        t_199[k] = pb_z[k] * sgi0_87[k]
                   - f_10 * pc_z[k] * sgi1_87[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_z, pc_x, pc_y, pc_z, sgi0_90, sgh_86, \
                         sgh_152, sgi1_90, shg0_110, shg1_110, shh_149, \
                         shh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_12 * sgh_86[k]
                   + f_3 * pc_y[k] * shh_149[k];

        t_201[k] = f_12 * sgh_152[k]
                   + f_4 * shg0_110[k]
                   - f_5 * shg1_110[k]
                   + f_3 * pc_x[k] * shh_152[k];

        t_202[k] = pb_z[k] * sgi0_90[k]
                   - f_10 * pc_z[k] * sgi1_90[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pc_x, pc_y, pc_z, sgh_66, sgh_89, sgh_156, \
                         shg0_114, shg1_114, shh_150, shh_152, \
                         shh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_11 * sgh_66[k]
                   + f_3 * pc_z[k] * shh_150[k];

        t_204[k] = f_12 * sgh_89[k]
                   + f_3 * pc_y[k] * shh_152[k];

        t_205[k] = f_12 * sgh_156[k]
                   + f_6 * shg0_114[k]
                   - f_7 * shg1_114[k]
                   + f_3 * pc_x[k] * shh_156[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pb_z, pc_y, pc_z, sgi0_94, sgi0_96, \
                         sgh_69, sgh_70, sgh_93, sgi1_94, sgi1_96, shh_153, \
                         shh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pb_z[k] * sgi0_94[k]
                   - f_10 * pc_z[k] * sgi1_94[k];

        t_207[k] = f_11 * sgh_69[k]
                   + f_3 * pc_z[k] * shh_153[k];

        t_208[k] = pb_z[k] * sgi0_96[k]
                   + f_12 * sgh_70[k]
                   - f_10 * pc_z[k] * sgi1_96[k];

        t_209[k] = f_12 * sgh_93[k]
                   + f_3 * pc_y[k] * shh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, sgh_161, sgh_162, sgh_163, sgh_164, \
                         shg0_119, shg1_119, shh_161, shh_162, shh_163, \
                         shh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_12 * sgh_161[k]
                   + f_8 * shg0_119[k]
                   - f_9 * shg1_119[k]
                   + f_3 * pc_x[k] * shh_161[k];

        t_211[k] = f_12 * sgh_162[k]
                   + f_3 * pc_x[k] * shh_162[k];

        t_212[k] = f_12 * sgh_163[k]
                   + f_3 * pc_x[k] * shh_163[k];

        t_213[k] = f_12 * sgh_164[k]
                   + f_3 * pc_x[k] * shh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pb_z, pc_x, pc_z, sgi0_105, sgh_165, \
                         sgh_166, sgh_167, sgi1_105, shh_165, shh_166, \
                         shh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_12 * sgh_165[k]
                   + f_3 * pc_x[k] * shh_165[k];

        t_215[k] = f_12 * sgh_166[k]
                   + f_3 * pc_x[k] * shh_166[k];

        t_216[k] = f_12 * sgh_167[k]
                   + f_3 * pc_x[k] * shh_167[k];

        t_217[k] = pb_z[k] * sgi0_105[k]
                   - f_10 * pc_z[k] * sgi1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, sgh_78, sgh_101, sgh_102, shg0_117, \
                         shg0_118, shg1_117, shg1_118, shh_162, shh_164, \
                         shh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * sgh_78[k]
                   + f_3 * pc_z[k] * shh_162[k];

        t_219[k] = f_12 * sgh_101[k]
                   + f_4 * shg0_117[k]
                   - f_5 * shg1_117[k]
                   + f_3 * pc_y[k] * shh_164[k];

        t_220[k] = f_12 * sgh_102[k]
                   + f_6 * shg0_118[k]
                   - f_7 * shg1_118[k]
                   + f_3 * pc_y[k] * shh_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pb_y, pc_y, pc_z, sgi0_140, sgh_83, \
                         sgh_103, sgh_104, sgi1_140, shg0_119, shg1_119, shh_166, \
                         shh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * sgh_103[k]
                   + f_8 * shg0_119[k]
                   - f_9 * shg1_119[k]
                   + f_3 * pc_y[k] * shh_166[k];

        t_222[k] = f_12 * sgh_104[k]
                   + f_3 * pc_y[k] * shh_167[k];

        t_223[k] = f_11 * sgh_83[k]
                   + f_1 * shg0_119[k]
                   - f_2 * shg1_119[k]
                   + f_3 * pc_z[k] * shh_167[k];

        t_224[k] = pb_y[k] * sgi0_140[k]
                   - f_10 * pc_y[k] * sgi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_y, pc_y, pc_z, sgi0_143, sgh_84, \
                         sgh_105, sgh_106, sgh_107, sgi1_143, shh_168, \
                         shh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * sgh_105[k]
                   + f_3 * pc_y[k] * shh_168[k];

        t_226[k] = f_12 * sgh_84[k]
                   + f_3 * pc_z[k] * shh_168[k];

        t_227[k] = pb_y[k] * sgi0_143[k]
                   + f_12 * sgh_106[k]
                   - f_10 * pc_y[k] * sgi1_143[k];

        t_228[k] = f_11 * sgh_107[k]
                   + f_3 * pc_y[k] * shh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_y, pc_y, pc_z, sgi0_145, sgi0_146, \
                         sgh_87, sgh_108, sgh_110, sgi1_145, sgi1_146, shh_171, \
                         shh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * sgi0_145[k]
                   - f_10 * pc_y[k] * sgi1_145[k];

        t_230[k] = pb_y[k] * sgi0_146[k]
                   + f_13 * sgh_108[k]
                   - f_10 * pc_y[k] * sgi1_146[k];

        t_231[k] = f_12 * sgh_87[k]
                   + f_3 * pc_z[k] * shh_171[k];

        t_232[k] = f_11 * sgh_110[k]
                   + f_3 * pc_y[k] * shh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_y, pc_y, pc_z, sgi0_149, sgi0_150, sgh_90, \
                         sgh_111, sgi1_149, sgi1_150, shh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * sgi0_149[k]
                   - f_10 * pc_y[k] * sgi1_149[k];

        t_234[k] = pb_y[k] * sgi0_150[k]
                   + f_14 * sgh_111[k]
                   - f_10 * pc_y[k] * sgi1_150[k];

        t_235[k] = f_12 * sgh_90[k]
                   + f_3 * pc_z[k] * shh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pb_y, pc_x, pc_y, sgi0_152, sgi0_154, \
                         sgh_113, sgh_114, sgh_183, sgi1_152, sgi1_154, shh_177, \
                         shh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pb_y[k] * sgi0_152[k]
                   + f_12 * sgh_113[k]
                   - f_10 * pc_y[k] * sgi1_152[k];

        t_237[k] = f_11 * sgh_114[k]
                   + f_3 * pc_y[k] * shh_177[k];

        t_238[k] = pb_y[k] * sgi0_154[k]
                   - f_10 * pc_y[k] * sgi1_154[k];

        t_239[k] = f_12 * sgh_183[k]
                   + f_3 * pc_x[k] * shh_183[k];
    }
}

static auto
compute_prim_shi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgi0,
                                                          const size_t sgh, const size_t sgi1,
                                                          const size_t shg0, const size_t shg1,
                                                          const size_t shh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgi0_167 = buffer.data(sgi0 + 167);
    const auto *sgi0_168 = buffer.data(sgi0 + 168);
    const auto *sgi0_171 = buffer.data(sgi0 + 171);
    const auto *sgi0_174 = buffer.data(sgi0 + 174);
    const auto *sgi0_178 = buffer.data(sgi0 + 178);
    const auto *sgi0_280 = buffer.data(sgi0 + 280);
    const auto *sgi0_283 = buffer.data(sgi0 + 283);
    const auto *sgi0_285 = buffer.data(sgi0 + 285);
    const auto *sgi0_286 = buffer.data(sgi0 + 286);
    const auto *sgi0_289 = buffer.data(sgi0 + 289);
    const auto *sgi0_290 = buffer.data(sgi0 + 290);
    const auto *sgi0_292 = buffer.data(sgi0 + 292);
    const auto *sgi0_294 = buffer.data(sgi0 + 294);
    const auto *sgi0_301 = buffer.data(sgi0 + 301);
    const auto *sgi0_303 = buffer.data(sgi0 + 303);
    const auto *sgi0_304 = buffer.data(sgi0 + 304);
    const auto *sgi0_305 = buffer.data(sgi0 + 305);
    const auto *sgi0_307 = buffer.data(sgi0 + 307);
    const auto *sgi0_313 = buffer.data(sgi0 + 313);
    const auto *sgi0_317 = buffer.data(sgi0 + 317);
    const auto *sgi0_320 = buffer.data(sgi0 + 320);
    const auto *sgi0_322 = buffer.data(sgi0 + 322);
    const auto *sgi0_329 = buffer.data(sgi0 + 329);
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
    const auto *sgi0_357 = buffer.data(sgi0 + 357);
    const auto *sgi0_359 = buffer.data(sgi0 + 359);
    const auto *sgi0_360 = buffer.data(sgi0 + 360);
    const auto *sgi0_361 = buffer.data(sgi0 + 361);

    const auto *sgh_99 = buffer.data(sgh + 99);
    const auto *sgh_105 = buffer.data(sgh + 105);
    const auto *sgh_108 = buffer.data(sgh + 108);
    const auto *sgh_111 = buffer.data(sgh + 111);
    const auto *sgh_120 = buffer.data(sgh + 120);
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
    const auto *sgh_146 = buffer.data(sgh + 146);
    const auto *sgh_147 = buffer.data(sgh + 147);
    const auto *sgh_149 = buffer.data(sgh + 149);
    const auto *sgh_150 = buffer.data(sgh + 150);
    const auto *sgh_152 = buffer.data(sgh + 152);
    const auto *sgh_153 = buffer.data(sgh + 153);
    const auto *sgh_156 = buffer.data(sgh + 156);
    const auto *sgh_162 = buffer.data(sgh + 162);
    const auto *sgh_167 = buffer.data(sgh + 167);
    const auto *sgh_168 = buffer.data(sgh + 168);
    const auto *sgh_170 = buffer.data(sgh + 170);
    const auto *sgh_173 = buffer.data(sgh + 173);
    const auto *sgh_177 = buffer.data(sgh + 177);
    const auto *sgh_184 = buffer.data(sgh + 184);
    const auto *sgh_185 = buffer.data(sgh + 185);
    const auto *sgh_186 = buffer.data(sgh + 186);
    const auto *sgh_187 = buffer.data(sgh + 187);
    const auto *sgh_188 = buffer.data(sgh + 188);
    const auto *sgh_189 = buffer.data(sgh + 189);
    const auto *sgh_192 = buffer.data(sgh + 192);
    const auto *sgh_194 = buffer.data(sgh + 194);
    const auto *sgh_195 = buffer.data(sgh + 195);
    const auto *sgh_198 = buffer.data(sgh + 198);
    const auto *sgh_199 = buffer.data(sgh + 199);
    const auto *sgh_201 = buffer.data(sgh + 201);
    const auto *sgh_203 = buffer.data(sgh + 203);
    const auto *sgh_204 = buffer.data(sgh + 204);
    const auto *sgh_205 = buffer.data(sgh + 205);
    const auto *sgh_206 = buffer.data(sgh + 206);
    const auto *sgh_207 = buffer.data(sgh + 207);
    const auto *sgh_208 = buffer.data(sgh + 208);
    const auto *sgh_209 = buffer.data(sgh + 209);
    const auto *sgh_210 = buffer.data(sgh + 210);
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
    const auto *sgh_236 = buffer.data(sgh + 236);
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

    const auto *sgi1_167 = buffer.data(sgi1 + 167);
    const auto *sgi1_168 = buffer.data(sgi1 + 168);
    const auto *sgi1_171 = buffer.data(sgi1 + 171);
    const auto *sgi1_174 = buffer.data(sgi1 + 174);
    const auto *sgi1_178 = buffer.data(sgi1 + 178);
    const auto *sgi1_280 = buffer.data(sgi1 + 280);
    const auto *sgi1_283 = buffer.data(sgi1 + 283);
    const auto *sgi1_285 = buffer.data(sgi1 + 285);
    const auto *sgi1_286 = buffer.data(sgi1 + 286);
    const auto *sgi1_289 = buffer.data(sgi1 + 289);
    const auto *sgi1_290 = buffer.data(sgi1 + 290);
    const auto *sgi1_292 = buffer.data(sgi1 + 292);
    const auto *sgi1_294 = buffer.data(sgi1 + 294);
    const auto *sgi1_301 = buffer.data(sgi1 + 301);
    const auto *sgi1_303 = buffer.data(sgi1 + 303);
    const auto *sgi1_304 = buffer.data(sgi1 + 304);
    const auto *sgi1_305 = buffer.data(sgi1 + 305);
    const auto *sgi1_307 = buffer.data(sgi1 + 307);
    const auto *sgi1_313 = buffer.data(sgi1 + 313);
    const auto *sgi1_317 = buffer.data(sgi1 + 317);
    const auto *sgi1_320 = buffer.data(sgi1 + 320);
    const auto *sgi1_322 = buffer.data(sgi1 + 322);
    const auto *sgi1_329 = buffer.data(sgi1 + 329);
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
    const auto *sgi1_357 = buffer.data(sgi1 + 357);
    const auto *sgi1_359 = buffer.data(sgi1 + 359);
    const auto *sgi1_360 = buffer.data(sgi1 + 360);
    const auto *sgi1_361 = buffer.data(sgi1 + 361);

    const auto *shg0_130 = buffer.data(shg0 + 130);
    const auto *shg0_132 = buffer.data(shg0 + 132);
    const auto *shg0_133 = buffer.data(shg0 + 133);
    const auto *shg0_134 = buffer.data(shg0 + 134);
    const auto *shg0_135 = buffer.data(shg0 + 135);
    const auto *shg0_138 = buffer.data(shg0 + 138);
    const auto *shg0_140 = buffer.data(shg0 + 140);
    const auto *shg0_141 = buffer.data(shg0 + 141);
    const auto *shg0_144 = buffer.data(shg0 + 144);
    const auto *shg0_145 = buffer.data(shg0 + 145);
    const auto *shg0_147 = buffer.data(shg0 + 147);
    const auto *shg0_148 = buffer.data(shg0 + 148);
    const auto *shg0_149 = buffer.data(shg0 + 149);

    const auto *shg1_130 = buffer.data(shg1 + 130);
    const auto *shg1_132 = buffer.data(shg1 + 132);
    const auto *shg1_133 = buffer.data(shg1 + 133);
    const auto *shg1_134 = buffer.data(shg1 + 134);
    const auto *shg1_135 = buffer.data(shg1 + 135);
    const auto *shg1_138 = buffer.data(shg1 + 138);
    const auto *shg1_140 = buffer.data(shg1 + 140);
    const auto *shg1_141 = buffer.data(shg1 + 141);
    const auto *shg1_144 = buffer.data(shg1 + 144);
    const auto *shg1_145 = buffer.data(shg1 + 145);
    const auto *shg1_147 = buffer.data(shg1 + 147);
    const auto *shg1_148 = buffer.data(shg1 + 148);
    const auto *shg1_149 = buffer.data(shg1 + 149);

    const auto *shh_183 = buffer.data(shh + 183);
    const auto *shh_184 = buffer.data(shh + 184);
    const auto *shh_185 = buffer.data(shh + 185);
    const auto *shh_186 = buffer.data(shh + 186);
    const auto *shh_187 = buffer.data(shh + 187);
    const auto *shh_188 = buffer.data(shh + 188);
    const auto *shh_189 = buffer.data(shh + 189);
    const auto *shh_191 = buffer.data(shh + 191);
    const auto *shh_192 = buffer.data(shh + 192);
    const auto *shh_194 = buffer.data(shh + 194);
    const auto *shh_195 = buffer.data(shh + 195);
    const auto *shh_198 = buffer.data(shh + 198);
    const auto *shh_199 = buffer.data(shh + 199);
    const auto *shh_201 = buffer.data(shh + 201);
    const auto *shh_203 = buffer.data(shh + 203);
    const auto *shh_204 = buffer.data(shh + 204);
    const auto *shh_205 = buffer.data(shh + 205);
    const auto *shh_206 = buffer.data(shh + 206);
    const auto *shh_207 = buffer.data(shh + 207);
    const auto *shh_208 = buffer.data(shh + 208);
    const auto *shh_209 = buffer.data(shh + 209);
    const auto *shh_210 = buffer.data(shh + 210);
    const auto *shh_212 = buffer.data(shh + 212);
    const auto *shh_213 = buffer.data(shh + 213);
    const auto *shh_215 = buffer.data(shh + 215);
    const auto *shh_216 = buffer.data(shh + 216);
    const auto *shh_219 = buffer.data(shh + 219);
    const auto *shh_225 = buffer.data(shh + 225);
    const auto *shh_226 = buffer.data(shh + 226);
    const auto *shh_227 = buffer.data(shh + 227);
    const auto *shh_228 = buffer.data(shh + 228);
    const auto *shh_229 = buffer.data(shh + 229);
    const auto *shh_230 = buffer.data(shh + 230);
    const auto *shh_231 = buffer.data(shh + 231);
    const auto *shh_233 = buffer.data(shh + 233);
    const auto *shh_234 = buffer.data(shh + 234);
    const auto *shh_236 = buffer.data(shh + 236);
    const auto *shh_237 = buffer.data(shh + 237);
    const auto *shh_240 = buffer.data(shh + 240);
    const auto *shh_246 = buffer.data(shh + 246);
    const auto *shh_247 = buffer.data(shh + 247);
    const auto *shh_248 = buffer.data(shh + 248);
    const auto *shh_249 = buffer.data(shh + 249);
    const auto *shh_250 = buffer.data(shh + 250);
    const auto *shh_251 = buffer.data(shh + 251);
    const auto *shh_252 = buffer.data(shh + 252);
    const auto *shh_254 = buffer.data(shh + 254);
    const auto *shh_255 = buffer.data(shh + 255);
    const auto *shh_257 = buffer.data(shh + 257);
    const auto *shh_258 = buffer.data(shh + 258);
    const auto *shh_261 = buffer.data(shh + 261);
    const auto *shh_267 = buffer.data(shh + 267);
    const auto *shh_268 = buffer.data(shh + 268);
    const auto *shh_269 = buffer.data(shh + 269);
    const auto *shh_270 = buffer.data(shh + 270);
    const auto *shh_271 = buffer.data(shh + 271);
    const auto *shh_272 = buffer.data(shh + 272);

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, sgh_184, sgh_185, sgh_186, \
                         sgh_187, sgh_188, shh_184, shh_185, shh_186, shh_187, \
                         shh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_12 * sgh_184[k]
                   + f_3 * pc_x[k] * shh_184[k];

        t_241[k] = f_12 * sgh_185[k]
                   + f_3 * pc_x[k] * shh_185[k];

        t_242[k] = f_12 * sgh_186[k]
                   + f_3 * pc_x[k] * shh_186[k];

        t_243[k] = f_12 * sgh_187[k]
                   + f_3 * pc_x[k] * shh_187[k];

        t_244[k] = f_12 * sgh_188[k]
                   + f_3 * pc_x[k] * shh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, sgh_99, sgh_120, sgh_122, shg0_130, \
                         shg0_132, shg1_130, shg1_132, shh_183, \
                         shh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * sgh_120[k]
                   + f_1 * shg0_130[k]
                   - f_2 * shg1_130[k]
                   + f_3 * pc_y[k] * shh_183[k];

        t_246[k] = f_12 * sgh_99[k]
                   + f_3 * pc_z[k] * shh_183[k];

        t_247[k] = f_11 * sgh_122[k]
                   + f_4 * shg0_132[k]
                   - f_5 * shg1_132[k]
                   + f_3 * pc_y[k] * shh_185[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, sgh_123, sgh_124, sgh_125, shg0_133, \
                         shg0_134, shg1_133, shg1_134, shh_186, shh_187, \
                         shh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * sgh_123[k]
                   + f_6 * shg0_133[k]
                   - f_7 * shg1_133[k]
                   + f_3 * pc_y[k] * shh_186[k];

        t_249[k] = f_11 * sgh_124[k]
                   + f_8 * shg0_134[k]
                   - f_9 * shg1_134[k]
                   + f_3 * pc_y[k] * shh_187[k];

        t_250[k] = f_11 * sgh_125[k]
                   + f_3 * pc_y[k] * shh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_y, pc_x, pc_y, pc_z, sgi0_167, \
                         sgh_105, sgh_189, sgi1_167, shg0_135, shg1_135, \
                         shh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pb_y[k] * sgi0_167[k]
                   - f_10 * pc_y[k] * sgi1_167[k];

        t_252[k] = f_12 * sgh_189[k]
                   + f_1 * shg0_135[k]
                   - f_2 * shg1_135[k]
                   + f_3 * pc_x[k] * shh_189[k];

        t_253[k] = f_3 * pc_y[k] * shh_189[k];

        t_254[k] = f_13 * sgh_105[k]
                   + f_3 * pc_z[k] * shh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, sgh_192, sgh_194, shg0_138, \
                         shg0_140, shg1_138, shg1_140, shh_191, shh_192, \
                         shh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_12 * sgh_192[k]
                   + f_4 * shg0_138[k]
                   - f_5 * shg1_138[k]
                   + f_3 * pc_x[k] * shh_192[k];

        t_256[k] = f_3 * pc_y[k] * shh_191[k];

        t_257[k] = f_12 * sgh_194[k]
                   + f_4 * shg0_140[k]
                   - f_5 * shg1_140[k]
                   + f_3 * pc_x[k] * shh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_x, pc_y, pc_z, sgh_108, sgh_195, shg0_141, \
                         shg1_141, shh_192, shh_194, shh_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_12 * sgh_195[k]
                   + f_6 * shg0_141[k]
                   - f_7 * shg1_141[k]
                   + f_3 * pc_x[k] * shh_195[k];

        t_259[k] = f_13 * sgh_108[k]
                   + f_3 * pc_z[k] * shh_192[k];

        t_260[k] = f_3 * pc_y[k] * shh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_z, sgh_111, sgh_198, sgh_199, shg0_144, \
                         shg0_145, shg1_144, shg1_145, shh_195, shh_198, \
                         shh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_12 * sgh_198[k]
                   + f_6 * shg0_144[k]
                   - f_7 * shg1_144[k]
                   + f_3 * pc_x[k] * shh_198[k];

        t_262[k] = f_12 * sgh_199[k]
                   + f_8 * shg0_145[k]
                   - f_9 * shg1_145[k]
                   + f_3 * pc_x[k] * shh_199[k];

        t_263[k] = f_13 * sgh_111[k]
                   + f_3 * pc_z[k] * shh_195[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, sgh_201, sgh_203, shg0_147, \
                         shg0_149, shg1_147, shg1_149, shh_198, shh_201, \
                         shh_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_12 * sgh_201[k]
                   + f_8 * shg0_147[k]
                   - f_9 * shg1_147[k]
                   + f_3 * pc_x[k] * shh_201[k];

        t_265[k] = f_3 * pc_y[k] * shh_198[k];

        t_266[k] = f_12 * sgh_203[k]
                   + f_8 * shg0_149[k]
                   - f_9 * shg1_149[k]
                   + f_3 * pc_x[k] * shh_203[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, pc_x, sgh_204, sgh_205, sgh_206, \
                         sgh_207, sgh_208, shh_204, shh_205, shh_206, shh_207, \
                         shh_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_12 * sgh_204[k]
                   + f_3 * pc_x[k] * shh_204[k];

        t_268[k] = f_12 * sgh_205[k]
                   + f_3 * pc_x[k] * shh_205[k];

        t_269[k] = f_12 * sgh_206[k]
                   + f_3 * pc_x[k] * shh_206[k];

        t_270[k] = f_12 * sgh_207[k]
                   + f_3 * pc_x[k] * shh_207[k];

        t_271[k] = f_12 * sgh_208[k]
                   + f_3 * pc_x[k] * shh_208[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, sgh_120, sgh_209, \
                         shg0_145, shg0_147, shg1_145, shg1_147, shh_204, shh_206, \
                         shh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_12 * sgh_209[k]
                   + f_3 * pc_x[k] * shh_209[k];

        t_273[k] = f_1 * shg0_145[k]
                   - f_2 * shg1_145[k]
                   + f_3 * pc_y[k] * shh_204[k];

        t_274[k] = f_13 * sgh_120[k]
                   + f_3 * pc_z[k] * shh_204[k];

        t_275[k] = f_4 * shg0_147[k]
                   - f_5 * shg1_147[k]
                   + f_3 * pc_y[k] * shh_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, sgh_125, shg0_148, shg0_149, \
                         shg1_148, shg1_149, shh_207, shh_208, \
                         shh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * shg0_148[k]
                   - f_7 * shg1_148[k]
                   + f_3 * pc_y[k] * shh_207[k];

        t_277[k] = f_8 * shg0_149[k]
                   - f_9 * shg1_149[k]
                   + f_3 * pc_y[k] * shh_208[k];

        t_278[k] = f_3 * pc_y[k] * shh_209[k];

        t_279[k] = f_13 * sgh_125[k]
                   + f_1 * shg0_149[k]
                   - f_2 * shg1_149[k]
                   + f_3 * pc_z[k] * shh_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pb_x, pc_x, pc_y, pc_z, sgi0_280, \
                         sgi0_283, sgh_126, sgh_210, sgh_213, sgi1_280, sgi1_283, \
                         shh_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pb_x[k] * sgi0_280[k]
                   + f_15 * sgh_210[k]
                   - f_10 * pc_x[k] * sgi1_280[k];

        t_281[k] = f_14 * sgh_126[k]
                   + f_3 * pc_y[k] * shh_210[k];

        t_282[k] = f_3 * pc_z[k] * shh_210[k];

        t_283[k] = pb_x[k] * sgi0_283[k]
                   + f_14 * sgh_213[k]
                   - f_10 * pc_x[k] * sgi1_283[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pb_x, pc_x, pc_y, sgi0_285, sgi0_286, sgh_128, \
                         sgh_215, sgh_216, sgi1_285, sgi1_286, \
                         shh_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_14 * sgh_128[k]
                   + f_3 * pc_y[k] * shh_212[k];

        t_285[k] = pb_x[k] * sgi0_285[k]
                   + f_14 * sgh_215[k]
                   - f_10 * pc_x[k] * sgi1_285[k];

        t_286[k] = pb_x[k] * sgi0_286[k]
                   + f_13 * sgh_216[k]
                   - f_10 * pc_x[k] * sgi1_286[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pb_x, pc_x, pc_y, pc_z, sgi0_289, sgh_131, \
                         sgh_219, sgi1_289, shh_213, shh_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_3 * pc_z[k] * shh_213[k];

        t_288[k] = f_14 * sgh_131[k]
                   + f_3 * pc_y[k] * shh_215[k];

        t_289[k] = pb_x[k] * sgi0_289[k]
                   + f_13 * sgh_219[k]
                   - f_10 * pc_x[k] * sgi1_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pb_x, pc_x, pc_z, sgi0_290, sgi0_292, sgh_220, \
                         sgh_222, sgi1_290, sgi1_292, shh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_x[k] * sgi0_290[k]
                   + f_12 * sgh_220[k]
                   - f_10 * pc_x[k] * sgi1_290[k];

        t_291[k] = f_3 * pc_z[k] * shh_216[k];

        t_292[k] = pb_x[k] * sgi0_292[k]
                   + f_12 * sgh_222[k]
                   - f_10 * pc_x[k] * sgi1_292[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pb_x, pc_x, pc_y, sgi0_294, sgh_135, \
                         sgh_224, sgh_225, sgh_226, sgi1_294, shh_219, shh_225, \
                         shh_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_14 * sgh_135[k]
                   + f_3 * pc_y[k] * shh_219[k];

        t_294[k] = pb_x[k] * sgi0_294[k]
                   + f_12 * sgh_224[k]
                   - f_10 * pc_x[k] * sgi1_294[k];

        t_295[k] = f_11 * sgh_225[k]
                   + f_3 * pc_x[k] * shh_225[k];

        t_296[k] = f_11 * sgh_226[k]
                   + f_3 * pc_x[k] * shh_226[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, pc_x, sgh_227, sgh_228, sgh_229, sgh_230, \
                         shh_227, shh_228, shh_229, shh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_11 * sgh_227[k]
                   + f_3 * pc_x[k] * shh_227[k];

        t_298[k] = f_11 * sgh_228[k]
                   + f_3 * pc_x[k] * shh_228[k];

        t_299[k] = f_11 * sgh_229[k]
                   + f_3 * pc_x[k] * shh_229[k];

        t_300[k] = f_11 * sgh_230[k]
                   + f_3 * pc_x[k] * shh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pb_x, pc_x, pc_z, sgi0_301, sgi0_303, \
                         sgi0_304, sgi1_301, sgi1_303, sgi1_304, \
                         shh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = pb_x[k] * sgi0_301[k]
                   - f_10 * pc_x[k] * sgi1_301[k];

        t_302[k] = f_3 * pc_z[k] * shh_225[k];

        t_303[k] = pb_x[k] * sgi0_303[k]
                   - f_10 * pc_x[k] * sgi1_303[k];

        t_304[k] = pb_x[k] * sgi0_304[k]
                   - f_10 * pc_x[k] * sgi1_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, pb_x, pc_x, pc_y, sgi0_305, sgi0_307, sgh_146, \
                         sgi1_305, sgi1_307, shh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pb_x[k] * sgi0_305[k]
                   - f_10 * pc_x[k] * sgi1_305[k];

        t_306[k] = f_14 * sgh_146[k]
                   + f_3 * pc_y[k] * shh_230[k];

        t_307[k] = pb_x[k] * sgi0_307[k]
                   - f_10 * pc_x[k] * sgi1_307[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pb_z, pc_y, pc_z, sgi0_168, sgi0_171, \
                         sgh_126, sgh_147, sgi1_168, sgi1_171, \
                         shh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pb_z[k] * sgi0_168[k]
                   - f_10 * pc_z[k] * sgi1_168[k];

        t_309[k] = f_13 * sgh_147[k]
                   + f_3 * pc_y[k] * shh_231[k];

        t_310[k] = f_11 * sgh_126[k]
                   + f_3 * pc_z[k] * shh_231[k];

        t_311[k] = pb_z[k] * sgi0_171[k]
                   - f_10 * pc_z[k] * sgi1_171[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pb_x, pb_z, pc_x, pc_y, pc_z, sgi0_174, \
                         sgi0_313, sgh_149, sgh_236, sgi1_174, sgi1_313, \
                         shh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_13 * sgh_149[k]
                   + f_3 * pc_y[k] * shh_233[k];

        t_313[k] = pb_x[k] * sgi0_313[k]
                   + f_14 * sgh_236[k]
                   - f_10 * pc_x[k] * sgi1_313[k];

        t_314[k] = pb_z[k] * sgi0_174[k]
                   - f_10 * pc_z[k] * sgi1_174[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pb_x, pc_x, pc_y, pc_z, sgi0_317, sgh_129, \
                         sgh_152, sgh_240, sgi1_317, shh_234, shh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_11 * sgh_129[k]
                   + f_3 * pc_z[k] * shh_234[k];

        t_316[k] = f_13 * sgh_152[k]
                   + f_3 * pc_y[k] * shh_236[k];

        t_317[k] = pb_x[k] * sgi0_317[k]
                   + f_13 * sgh_240[k]
                   - f_10 * pc_x[k] * sgi1_317[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pb_x, pb_z, pc_x, pc_z, sgi0_178, sgi0_320, \
                         sgh_132, sgh_243, sgi1_178, sgi1_320, \
                         shh_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pb_z[k] * sgi0_178[k]
                   - f_10 * pc_z[k] * sgi1_178[k];

        t_319[k] = f_11 * sgh_132[k]
                   + f_3 * pc_z[k] * shh_237[k];

        t_320[k] = pb_x[k] * sgi0_320[k]
                   + f_12 * sgh_243[k]
                   - f_10 * pc_x[k] * sgi1_320[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pb_x, pc_x, pc_y, sgi0_322, sgh_156, \
                         sgh_245, sgh_246, sgh_247, sgi1_322, shh_240, shh_246, \
                         shh_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * sgh_156[k]
                   + f_3 * pc_y[k] * shh_240[k];

        t_322[k] = pb_x[k] * sgi0_322[k]
                   + f_12 * sgh_245[k]
                   - f_10 * pc_x[k] * sgi1_322[k];

        t_323[k] = f_11 * sgh_246[k]
                   + f_3 * pc_x[k] * shh_246[k];

        t_324[k] = f_11 * sgh_247[k]
                   + f_3 * pc_x[k] * shh_247[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pc_x, sgh_248, sgh_249, sgh_250, sgh_251, \
                         shh_248, shh_249, shh_250, shh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_11 * sgh_248[k]
                   + f_3 * pc_x[k] * shh_248[k];

        t_326[k] = f_11 * sgh_249[k]
                   + f_3 * pc_x[k] * shh_249[k];

        t_327[k] = f_11 * sgh_250[k]
                   + f_3 * pc_x[k] * shh_250[k];

        t_328[k] = f_11 * sgh_251[k]
                   + f_3 * pc_x[k] * shh_251[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pb_x, pc_x, pc_z, sgi0_329, sgi0_331, \
                         sgi0_332, sgh_141, sgi1_329, sgi1_331, sgi1_332, \
                         shh_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = pb_x[k] * sgi0_329[k]
                   - f_10 * pc_x[k] * sgi1_329[k];

        t_330[k] = f_11 * sgh_141[k]
                   + f_3 * pc_z[k] * shh_246[k];

        t_331[k] = pb_x[k] * sgi0_331[k]
                   - f_10 * pc_x[k] * sgi1_331[k];

        t_332[k] = pb_x[k] * sgi0_332[k]
                   - f_10 * pc_x[k] * sgi1_332[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, pb_x, pc_x, pc_y, sgi0_333, sgi0_335, \
                         sgi0_336, sgh_167, sgh_252, sgi1_333, sgi1_335, sgi1_336, \
                         shh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = pb_x[k] * sgi0_333[k]
                   - f_10 * pc_x[k] * sgi1_333[k];

        t_334[k] = f_13 * sgh_167[k]
                   + f_3 * pc_y[k] * shh_251[k];

        t_335[k] = pb_x[k] * sgi0_335[k]
                   - f_10 * pc_x[k] * sgi1_335[k];

        t_336[k] = pb_x[k] * sgi0_336[k]
                   + f_15 * sgh_252[k]
                   - f_10 * pc_x[k] * sgi1_336[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pb_x, pc_x, pc_y, pc_z, sgi0_339, \
                         sgh_147, sgh_168, sgh_170, sgh_255, sgi1_339, shh_252, \
                         shh_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_12 * sgh_168[k]
                   + f_3 * pc_y[k] * shh_252[k];

        t_338[k] = f_12 * sgh_147[k]
                   + f_3 * pc_z[k] * shh_252[k];

        t_339[k] = pb_x[k] * sgi0_339[k]
                   + f_14 * sgh_255[k]
                   - f_10 * pc_x[k] * sgi1_339[k];

        t_340[k] = f_12 * sgh_170[k]
                   + f_3 * pc_y[k] * shh_254[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pb_x, pc_x, pc_z, sgi0_341, sgi0_342, sgh_150, \
                         sgh_257, sgh_258, sgi1_341, sgi1_342, \
                         shh_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = pb_x[k] * sgi0_341[k]
                   + f_14 * sgh_257[k]
                   - f_10 * pc_x[k] * sgi1_341[k];

        t_342[k] = pb_x[k] * sgi0_342[k]
                   + f_13 * sgh_258[k]
                   - f_10 * pc_x[k] * sgi1_342[k];

        t_343[k] = f_12 * sgh_150[k]
                   + f_3 * pc_z[k] * shh_255[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, pb_x, pc_x, pc_y, sgi0_345, sgi0_346, sgh_173, \
                         sgh_261, sgh_262, sgi1_345, sgi1_346, \
                         shh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_12 * sgh_173[k]
                   + f_3 * pc_y[k] * shh_257[k];

        t_345[k] = pb_x[k] * sgi0_345[k]
                   + f_13 * sgh_261[k]
                   - f_10 * pc_x[k] * sgi1_345[k];

        t_346[k] = pb_x[k] * sgi0_346[k]
                   + f_12 * sgh_262[k]
                   - f_10 * pc_x[k] * sgi1_346[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, pb_x, pc_x, pc_y, pc_z, sgi0_348, sgh_153, \
                         sgh_177, sgh_264, sgi1_348, shh_258, shh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_12 * sgh_153[k]
                   + f_3 * pc_z[k] * shh_258[k];

        t_348[k] = pb_x[k] * sgi0_348[k]
                   + f_12 * sgh_264[k]
                   - f_10 * pc_x[k] * sgi1_348[k];

        t_349[k] = f_12 * sgh_177[k]
                   + f_3 * pc_y[k] * shh_261[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pb_x, pc_x, sgi0_350, sgh_266, sgh_267, \
                         sgh_268, sgh_269, sgi1_350, shh_267, shh_268, \
                         shh_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = pb_x[k] * sgi0_350[k]
                   + f_12 * sgh_266[k]
                   - f_10 * pc_x[k] * sgi1_350[k];

        t_351[k] = f_11 * sgh_267[k]
                   + f_3 * pc_x[k] * shh_267[k];

        t_352[k] = f_11 * sgh_268[k]
                   + f_3 * pc_x[k] * shh_268[k];

        t_353[k] = f_11 * sgh_269[k]
                   + f_3 * pc_x[k] * shh_269[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pb_x, pc_x, sgi0_357, sgh_270, sgh_271, \
                         sgh_272, sgi1_357, shh_270, shh_271, shh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_11 * sgh_270[k]
                   + f_3 * pc_x[k] * shh_270[k];

        t_355[k] = f_11 * sgh_271[k]
                   + f_3 * pc_x[k] * shh_271[k];

        t_356[k] = f_11 * sgh_272[k]
                   + f_3 * pc_x[k] * shh_272[k];

        t_357[k] = pb_x[k] * sgi0_357[k]
                   - f_10 * pc_x[k] * sgi1_357[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pb_x, pc_x, pc_z, sgi0_359, sgi0_360, \
                         sgi0_361, sgh_162, sgi1_359, sgi1_360, sgi1_361, \
                         shh_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_12 * sgh_162[k]
                   + f_3 * pc_z[k] * shh_267[k];

        t_359[k] = pb_x[k] * sgi0_359[k]
                   - f_10 * pc_x[k] * sgi1_359[k];

        t_360[k] = pb_x[k] * sgi0_360[k]
                   - f_10 * pc_x[k] * sgi1_360[k];

        t_361[k] = pb_x[k] * sgi0_361[k]
                   - f_10 * pc_x[k] * sgi1_361[k];
    }
}

static auto
compute_prim_shi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgi0,
                                                          const size_t sgh, const size_t sgi1,
                                                          const size_t shg0, const size_t shg1,
                                                          const size_t shh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;

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
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgi0_252 = buffer.data(sgi0 + 252);
    const auto *sgi0_257 = buffer.data(sgi0 + 257);
    const auto *sgi0_261 = buffer.data(sgi0 + 261);
    const auto *sgi0_266 = buffer.data(sgi0 + 266);
    const auto *sgi0_280 = buffer.data(sgi0 + 280);
    const auto *sgi0_283 = buffer.data(sgi0 + 283);
    const auto *sgi0_286 = buffer.data(sgi0 + 286);
    const auto *sgi0_290 = buffer.data(sgi0 + 290);
    const auto *sgi0_301 = buffer.data(sgi0 + 301);
    const auto *sgi0_303 = buffer.data(sgi0 + 303);
    const auto *sgi0_304 = buffer.data(sgi0 + 304);
    const auto *sgi0_305 = buffer.data(sgi0 + 305);
    const auto *sgi0_363 = buffer.data(sgi0 + 363);
    const auto *sgi0_367 = buffer.data(sgi0 + 367);
    const auto *sgi0_370 = buffer.data(sgi0 + 370);
    const auto *sgi0_374 = buffer.data(sgi0 + 374);
    const auto *sgi0_376 = buffer.data(sgi0 + 376);
    const auto *sgi0_385 = buffer.data(sgi0 + 385);
    const auto *sgi0_387 = buffer.data(sgi0 + 387);
    const auto *sgi0_388 = buffer.data(sgi0 + 388);
    const auto *sgi0_389 = buffer.data(sgi0 + 389);
    const auto *sgi0_391 = buffer.data(sgi0 + 391);
    const auto *sgi0_392 = buffer.data(sgi0 + 392);
    const auto *sgi0_395 = buffer.data(sgi0 + 395);
    const auto *sgi0_397 = buffer.data(sgi0 + 397);
    const auto *sgi0_398 = buffer.data(sgi0 + 398);
    const auto *sgi0_401 = buffer.data(sgi0 + 401);
    const auto *sgi0_402 = buffer.data(sgi0 + 402);
    const auto *sgi0_404 = buffer.data(sgi0 + 404);
    const auto *sgi0_406 = buffer.data(sgi0 + 406);
    const auto *sgi0_413 = buffer.data(sgi0 + 413);
    const auto *sgi0_415 = buffer.data(sgi0 + 415);
    const auto *sgi0_416 = buffer.data(sgi0 + 416);
    const auto *sgi0_417 = buffer.data(sgi0 + 417);
    const auto *sgi0_419 = buffer.data(sgi0 + 419);

    const auto *sgh_168 = buffer.data(sgh + 168);
    const auto *sgh_171 = buffer.data(sgh + 171);
    const auto *sgh_174 = buffer.data(sgh + 174);
    const auto *sgh_183 = buffer.data(sgh + 183);
    const auto *sgh_188 = buffer.data(sgh + 188);
    const auto *sgh_189 = buffer.data(sgh + 189);
    const auto *sgh_191 = buffer.data(sgh + 191);
    const auto *sgh_192 = buffer.data(sgh + 192);
    const auto *sgh_194 = buffer.data(sgh + 194);
    const auto *sgh_195 = buffer.data(sgh + 195);
    const auto *sgh_198 = buffer.data(sgh + 198);
    const auto *sgh_204 = buffer.data(sgh + 204);
    const auto *sgh_209 = buffer.data(sgh + 209);
    const auto *sgh_210 = buffer.data(sgh + 210);
    const auto *sgh_212 = buffer.data(sgh + 212);
    const auto *sgh_213 = buffer.data(sgh + 213);
    const auto *sgh_215 = buffer.data(sgh + 215);
    const auto *sgh_216 = buffer.data(sgh + 216);
    const auto *sgh_219 = buffer.data(sgh + 219);
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
    const auto *sgh_251 = buffer.data(sgh + 251);
    const auto *sgh_252 = buffer.data(sgh + 252);
    const auto *sgh_254 = buffer.data(sgh + 254);
    const auto *sgh_257 = buffer.data(sgh + 257);
    const auto *sgh_276 = buffer.data(sgh + 276);
    const auto *sgh_279 = buffer.data(sgh + 279);
    const auto *sgh_283 = buffer.data(sgh + 283);
    const auto *sgh_285 = buffer.data(sgh + 285);
    const auto *sgh_288 = buffer.data(sgh + 288);
    const auto *sgh_289 = buffer.data(sgh + 289);
    const auto *sgh_290 = buffer.data(sgh + 290);
    const auto *sgh_291 = buffer.data(sgh + 291);
    const auto *sgh_292 = buffer.data(sgh + 292);
    const auto *sgh_293 = buffer.data(sgh + 293);
    const auto *sgh_294 = buffer.data(sgh + 294);
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

    const auto *sgi1_252 = buffer.data(sgi1 + 252);
    const auto *sgi1_257 = buffer.data(sgi1 + 257);
    const auto *sgi1_261 = buffer.data(sgi1 + 261);
    const auto *sgi1_266 = buffer.data(sgi1 + 266);
    const auto *sgi1_280 = buffer.data(sgi1 + 280);
    const auto *sgi1_283 = buffer.data(sgi1 + 283);
    const auto *sgi1_286 = buffer.data(sgi1 + 286);
    const auto *sgi1_290 = buffer.data(sgi1 + 290);
    const auto *sgi1_301 = buffer.data(sgi1 + 301);
    const auto *sgi1_303 = buffer.data(sgi1 + 303);
    const auto *sgi1_304 = buffer.data(sgi1 + 304);
    const auto *sgi1_305 = buffer.data(sgi1 + 305);
    const auto *sgi1_363 = buffer.data(sgi1 + 363);
    const auto *sgi1_367 = buffer.data(sgi1 + 367);
    const auto *sgi1_370 = buffer.data(sgi1 + 370);
    const auto *sgi1_374 = buffer.data(sgi1 + 374);
    const auto *sgi1_376 = buffer.data(sgi1 + 376);
    const auto *sgi1_385 = buffer.data(sgi1 + 385);
    const auto *sgi1_387 = buffer.data(sgi1 + 387);
    const auto *sgi1_388 = buffer.data(sgi1 + 388);
    const auto *sgi1_389 = buffer.data(sgi1 + 389);
    const auto *sgi1_391 = buffer.data(sgi1 + 391);
    const auto *sgi1_392 = buffer.data(sgi1 + 392);
    const auto *sgi1_395 = buffer.data(sgi1 + 395);
    const auto *sgi1_397 = buffer.data(sgi1 + 397);
    const auto *sgi1_398 = buffer.data(sgi1 + 398);
    const auto *sgi1_401 = buffer.data(sgi1 + 401);
    const auto *sgi1_402 = buffer.data(sgi1 + 402);
    const auto *sgi1_404 = buffer.data(sgi1 + 404);
    const auto *sgi1_406 = buffer.data(sgi1 + 406);
    const auto *sgi1_413 = buffer.data(sgi1 + 413);
    const auto *sgi1_415 = buffer.data(sgi1 + 415);
    const auto *sgi1_416 = buffer.data(sgi1 + 416);
    const auto *sgi1_417 = buffer.data(sgi1 + 417);
    const auto *sgi1_419 = buffer.data(sgi1 + 419);

    const auto *shg0_225 = buffer.data(shg0 + 225);
    const auto *shg0_228 = buffer.data(shg0 + 228);
    const auto *shg0_230 = buffer.data(shg0 + 230);
    const auto *shg0_231 = buffer.data(shg0 + 231);
    const auto *shg0_234 = buffer.data(shg0 + 234);
    const auto *shg0_235 = buffer.data(shg0 + 235);
    const auto *shg0_237 = buffer.data(shg0 + 237);
    const auto *shg0_238 = buffer.data(shg0 + 238);
    const auto *shg0_239 = buffer.data(shg0 + 239);
    const auto *shg0_245 = buffer.data(shg0 + 245);
    const auto *shg0_249 = buffer.data(shg0 + 249);
    const auto *shg0_252 = buffer.data(shg0 + 252);
    const auto *shg0_254 = buffer.data(shg0 + 254);
    const auto *shg0_255 = buffer.data(shg0 + 255);
    const auto *shg0_258 = buffer.data(shg0 + 258);
    const auto *shg0_260 = buffer.data(shg0 + 260);
    const auto *shg0_261 = buffer.data(shg0 + 261);
    const auto *shg0_264 = buffer.data(shg0 + 264);
    const auto *shg0_265 = buffer.data(shg0 + 265);

    const auto *shg1_225 = buffer.data(shg1 + 225);
    const auto *shg1_228 = buffer.data(shg1 + 228);
    const auto *shg1_230 = buffer.data(shg1 + 230);
    const auto *shg1_231 = buffer.data(shg1 + 231);
    const auto *shg1_234 = buffer.data(shg1 + 234);
    const auto *shg1_235 = buffer.data(shg1 + 235);
    const auto *shg1_237 = buffer.data(shg1 + 237);
    const auto *shg1_238 = buffer.data(shg1 + 238);
    const auto *shg1_239 = buffer.data(shg1 + 239);
    const auto *shg1_245 = buffer.data(shg1 + 245);
    const auto *shg1_249 = buffer.data(shg1 + 249);
    const auto *shg1_252 = buffer.data(shg1 + 252);
    const auto *shg1_254 = buffer.data(shg1 + 254);
    const auto *shg1_255 = buffer.data(shg1 + 255);
    const auto *shg1_258 = buffer.data(shg1 + 258);
    const auto *shg1_260 = buffer.data(shg1 + 260);
    const auto *shg1_261 = buffer.data(shg1 + 261);
    const auto *shg1_264 = buffer.data(shg1 + 264);
    const auto *shg1_265 = buffer.data(shg1 + 265);

    const auto *shh_272 = buffer.data(shh + 272);
    const auto *shh_273 = buffer.data(shh + 273);
    const auto *shh_275 = buffer.data(shh + 275);
    const auto *shh_276 = buffer.data(shh + 276);
    const auto *shh_278 = buffer.data(shh + 278);
    const auto *shh_279 = buffer.data(shh + 279);
    const auto *shh_282 = buffer.data(shh + 282);
    const auto *shh_288 = buffer.data(shh + 288);
    const auto *shh_289 = buffer.data(shh + 289);
    const auto *shh_290 = buffer.data(shh + 290);
    const auto *shh_291 = buffer.data(shh + 291);
    const auto *shh_292 = buffer.data(shh + 292);
    const auto *shh_293 = buffer.data(shh + 293);
    const auto *shh_294 = buffer.data(shh + 294);
    const auto *shh_296 = buffer.data(shh + 296);
    const auto *shh_297 = buffer.data(shh + 297);
    const auto *shh_299 = buffer.data(shh + 299);
    const auto *shh_300 = buffer.data(shh + 300);
    const auto *shh_303 = buffer.data(shh + 303);
    const auto *shh_309 = buffer.data(shh + 309);
    const auto *shh_310 = buffer.data(shh + 310);
    const auto *shh_311 = buffer.data(shh + 311);
    const auto *shh_312 = buffer.data(shh + 312);
    const auto *shh_313 = buffer.data(shh + 313);
    const auto *shh_314 = buffer.data(shh + 314);
    const auto *shh_315 = buffer.data(shh + 315);
    const auto *shh_317 = buffer.data(shh + 317);
    const auto *shh_318 = buffer.data(shh + 318);
    const auto *shh_320 = buffer.data(shh + 320);
    const auto *shh_321 = buffer.data(shh + 321);
    const auto *shh_324 = buffer.data(shh + 324);
    const auto *shh_325 = buffer.data(shh + 325);
    const auto *shh_327 = buffer.data(shh + 327);
    const auto *shh_329 = buffer.data(shh + 329);
    const auto *shh_330 = buffer.data(shh + 330);
    const auto *shh_331 = buffer.data(shh + 331);
    const auto *shh_332 = buffer.data(shh + 332);
    const auto *shh_333 = buffer.data(shh + 333);
    const auto *shh_334 = buffer.data(shh + 334);
    const auto *shh_335 = buffer.data(shh + 335);
    const auto *shh_336 = buffer.data(shh + 336);
    const auto *shh_338 = buffer.data(shh + 338);
    const auto *shh_339 = buffer.data(shh + 339);
    const auto *shh_341 = buffer.data(shh + 341);
    const auto *shh_342 = buffer.data(shh + 342);
    const auto *shh_345 = buffer.data(shh + 345);
    const auto *shh_348 = buffer.data(shh + 348);
    const auto *shh_350 = buffer.data(shh + 350);
    const auto *shh_351 = buffer.data(shh + 351);
    const auto *shh_352 = buffer.data(shh + 352);
    const auto *shh_353 = buffer.data(shh + 353);
    const auto *shh_354 = buffer.data(shh + 354);
    const auto *shh_355 = buffer.data(shh + 355);
    const auto *shh_356 = buffer.data(shh + 356);
    const auto *shh_357 = buffer.data(shh + 357);
    const auto *shh_359 = buffer.data(shh + 359);
    const auto *shh_360 = buffer.data(shh + 360);
    const auto *shh_362 = buffer.data(shh + 362);
    const auto *shh_363 = buffer.data(shh + 363);
    const auto *shh_366 = buffer.data(shh + 366);
    const auto *shh_367 = buffer.data(shh + 367);

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pb_x, pb_y, pc_x, pc_y, sgi0_252, \
                         sgi0_363, sgh_188, sgh_189, sgi1_252, sgi1_363, shh_272, \
                         shh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * sgh_188[k]
                   + f_3 * pc_y[k] * shh_272[k];

        t_363[k] = pb_x[k] * sgi0_363[k]
                   - f_10 * pc_x[k] * sgi1_363[k];

        t_364[k] = pb_y[k] * sgi0_252[k]
                   - f_10 * pc_y[k] * sgi1_252[k];

        t_365[k] = f_11 * sgh_189[k]
                   + f_3 * pc_y[k] * shh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, pb_x, pc_x, pc_y, pc_z, sgi0_367, sgh_168, \
                         sgh_191, sgh_276, sgi1_367, shh_273, shh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * sgh_168[k]
                   + f_3 * pc_z[k] * shh_273[k];

        t_367[k] = pb_x[k] * sgi0_367[k]
                   + f_14 * sgh_276[k]
                   - f_10 * pc_x[k] * sgi1_367[k];

        t_368[k] = f_11 * sgh_191[k]
                   + f_3 * pc_y[k] * shh_275[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pb_x, pb_y, pc_x, pc_y, pc_z, sgi0_257, \
                         sgi0_370, sgh_171, sgh_279, sgi1_257, sgi1_370, \
                         shh_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = pb_y[k] * sgi0_257[k]
                   - f_10 * pc_y[k] * sgi1_257[k];

        t_370[k] = pb_x[k] * sgi0_370[k]
                   + f_13 * sgh_279[k]
                   - f_10 * pc_x[k] * sgi1_370[k];

        t_371[k] = f_13 * sgh_171[k]
                   + f_3 * pc_z[k] * shh_276[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pb_x, pb_y, pc_x, pc_y, sgi0_261, sgi0_374, \
                         sgh_194, sgh_283, sgi1_261, sgi1_374, \
                         shh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_11 * sgh_194[k]
                   + f_3 * pc_y[k] * shh_278[k];

        t_373[k] = pb_y[k] * sgi0_261[k]
                   - f_10 * pc_y[k] * sgi1_261[k];

        t_374[k] = pb_x[k] * sgi0_374[k]
                   + f_12 * sgh_283[k]
                   - f_10 * pc_x[k] * sgi1_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pb_x, pc_x, pc_y, pc_z, sgi0_376, sgh_174, \
                         sgh_198, sgh_285, sgi1_376, shh_279, shh_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_13 * sgh_174[k]
                   + f_3 * pc_z[k] * shh_279[k];

        t_376[k] = pb_x[k] * sgi0_376[k]
                   + f_12 * sgh_285[k]
                   - f_10 * pc_x[k] * sgi1_376[k];

        t_377[k] = f_11 * sgh_198[k]
                   + f_3 * pc_y[k] * shh_282[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pb_y, pc_x, pc_y, sgi0_266, sgh_288, \
                         sgh_289, sgh_290, sgi1_266, shh_288, shh_289, \
                         shh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pb_y[k] * sgi0_266[k]
                   - f_10 * pc_y[k] * sgi1_266[k];

        t_379[k] = f_11 * sgh_288[k]
                   + f_3 * pc_x[k] * shh_288[k];

        t_380[k] = f_11 * sgh_289[k]
                   + f_3 * pc_x[k] * shh_289[k];

        t_381[k] = f_11 * sgh_290[k]
                   + f_3 * pc_x[k] * shh_290[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, pb_x, pc_x, sgi0_385, sgh_291, sgh_292, \
                         sgh_293, sgi1_385, shh_291, shh_292, shh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_11 * sgh_291[k]
                   + f_3 * pc_x[k] * shh_291[k];

        t_383[k] = f_11 * sgh_292[k]
                   + f_3 * pc_x[k] * shh_292[k];

        t_384[k] = f_11 * sgh_293[k]
                   + f_3 * pc_x[k] * shh_293[k];

        t_385[k] = pb_x[k] * sgi0_385[k]
                   - f_10 * pc_x[k] * sgi1_385[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pb_x, pc_x, pc_z, sgi0_387, sgi0_388, \
                         sgi0_389, sgh_183, sgi1_387, sgi1_388, sgi1_389, \
                         shh_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_13 * sgh_183[k]
                   + f_3 * pc_z[k] * shh_288[k];

        t_387[k] = pb_x[k] * sgi0_387[k]
                   - f_10 * pc_x[k] * sgi1_387[k];

        t_388[k] = pb_x[k] * sgi0_388[k]
                   - f_10 * pc_x[k] * sgi1_388[k];

        t_389[k] = pb_x[k] * sgi0_389[k]
                   - f_10 * pc_x[k] * sgi1_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pb_x, pc_x, pc_y, sgi0_391, sgi0_392, \
                         sgh_209, sgh_294, sgi1_391, sgi1_392, shh_293, \
                         shh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_11 * sgh_209[k]
                   + f_3 * pc_y[k] * shh_293[k];

        t_391[k] = pb_x[k] * sgi0_391[k]
                   - f_10 * pc_x[k] * sgi1_391[k];

        t_392[k] = pb_x[k] * sgi0_392[k]
                   + f_15 * sgh_294[k]
                   - f_10 * pc_x[k] * sgi1_392[k];

        t_393[k] = f_3 * pc_y[k] * shh_294[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pb_x, pc_x, pc_y, pc_z, sgi0_395, sgh_189, \
                         sgh_297, sgi1_395, shh_294, shh_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_14 * sgh_189[k]
                   + f_3 * pc_z[k] * shh_294[k];

        t_395[k] = pb_x[k] * sgi0_395[k]
                   + f_14 * sgh_297[k]
                   - f_10 * pc_x[k] * sgi1_395[k];

        t_396[k] = f_3 * pc_y[k] * shh_296[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, pb_x, pc_x, pc_z, sgi0_397, sgi0_398, sgh_192, \
                         sgh_299, sgh_300, sgi1_397, sgi1_398, \
                         shh_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = pb_x[k] * sgi0_397[k]
                   + f_14 * sgh_299[k]
                   - f_10 * pc_x[k] * sgi1_397[k];

        t_398[k] = pb_x[k] * sgi0_398[k]
                   + f_13 * sgh_300[k]
                   - f_10 * pc_x[k] * sgi1_398[k];

        t_399[k] = f_14 * sgh_192[k]
                   + f_3 * pc_z[k] * shh_297[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pb_x, pc_x, pc_y, sgi0_401, sgi0_402, sgh_303, \
                         sgh_304, sgi1_401, sgi1_402, shh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_3 * pc_y[k] * shh_299[k];

        t_401[k] = pb_x[k] * sgi0_401[k]
                   + f_13 * sgh_303[k]
                   - f_10 * pc_x[k] * sgi1_401[k];

        t_402[k] = pb_x[k] * sgi0_402[k]
                   + f_12 * sgh_304[k]
                   - f_10 * pc_x[k] * sgi1_402[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, pb_x, pc_x, pc_y, pc_z, sgi0_404, sgh_195, \
                         sgh_306, sgi1_404, shh_300, shh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_14 * sgh_195[k]
                   + f_3 * pc_z[k] * shh_300[k];

        t_404[k] = pb_x[k] * sgi0_404[k]
                   + f_12 * sgh_306[k]
                   - f_10 * pc_x[k] * sgi1_404[k];

        t_405[k] = f_3 * pc_y[k] * shh_303[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pb_x, pc_x, sgi0_406, sgh_308, sgh_309, \
                         sgh_310, sgh_311, sgi1_406, shh_309, shh_310, \
                         shh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = pb_x[k] * sgi0_406[k]
                   + f_12 * sgh_308[k]
                   - f_10 * pc_x[k] * sgi1_406[k];

        t_407[k] = f_11 * sgh_309[k]
                   + f_3 * pc_x[k] * shh_309[k];

        t_408[k] = f_11 * sgh_310[k]
                   + f_3 * pc_x[k] * shh_310[k];

        t_409[k] = f_11 * sgh_311[k]
                   + f_3 * pc_x[k] * shh_311[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pb_x, pc_x, sgi0_413, sgh_312, sgh_313, \
                         sgh_314, sgi1_413, shh_312, shh_313, shh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_11 * sgh_312[k]
                   + f_3 * pc_x[k] * shh_312[k];

        t_411[k] = f_11 * sgh_313[k]
                   + f_3 * pc_x[k] * shh_313[k];

        t_412[k] = f_11 * sgh_314[k]
                   + f_3 * pc_x[k] * shh_314[k];

        t_413[k] = pb_x[k] * sgi0_413[k]
                   - f_10 * pc_x[k] * sgi1_413[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pb_x, pc_x, pc_z, sgi0_415, sgi0_416, \
                         sgi0_417, sgh_204, sgi1_415, sgi1_416, sgi1_417, \
                         shh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_14 * sgh_204[k]
                   + f_3 * pc_z[k] * shh_309[k];

        t_415[k] = pb_x[k] * sgi0_415[k]
                   - f_10 * pc_x[k] * sgi1_415[k];

        t_416[k] = pb_x[k] * sgi0_416[k]
                   - f_10 * pc_x[k] * sgi1_416[k];

        t_417[k] = pb_x[k] * sgi0_417[k]
                   - f_10 * pc_x[k] * sgi1_417[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, t_422, pb_x, pc_x, pc_y, pc_z, sgi0_419, \
                         sgh_210, sgi1_419, shg0_225, shg1_225, shh_314, \
                         shh_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_3 * pc_y[k] * shh_314[k];

        t_419[k] = pb_x[k] * sgi0_419[k]
                   - f_10 * pc_x[k] * sgi1_419[k];

        t_420[k] = f_1 * shg0_225[k]
                   - f_2 * shg1_225[k]
                   + f_3 * pc_x[k] * shh_315[k];

        t_421[k] = f_0 * sgh_210[k]
                   + f_3 * pc_y[k] * shh_315[k];

        t_422[k] = f_3 * pc_z[k] * shh_315[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pc_x, pc_y, sgh_212, shg0_228, shg0_230, \
                         shg1_228, shg1_230, shh_317, shh_318, \
                         shh_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_4 * shg0_228[k]
                   - f_5 * shg1_228[k]
                   + f_3 * pc_x[k] * shh_318[k];

        t_424[k] = f_0 * sgh_212[k]
                   + f_3 * pc_y[k] * shh_317[k];

        t_425[k] = f_4 * shg0_230[k]
                   - f_5 * shg1_230[k]
                   + f_3 * pc_x[k] * shh_320[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, pc_x, pc_y, pc_z, sgh_215, shg0_231, \
                         shg0_234, shg1_231, shg1_234, shh_318, shh_320, shh_321, \
                         shh_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_6 * shg0_231[k]
                   - f_7 * shg1_231[k]
                   + f_3 * pc_x[k] * shh_321[k];

        t_427[k] = f_3 * pc_z[k] * shh_318[k];

        t_428[k] = f_0 * sgh_215[k]
                   + f_3 * pc_y[k] * shh_320[k];

        t_429[k] = f_6 * shg0_234[k]
                   - f_7 * shg1_234[k]
                   + f_3 * pc_x[k] * shh_324[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, pc_x, pc_y, pc_z, sgh_219, shg0_235, \
                         shg0_237, shg1_235, shg1_237, shh_321, shh_324, shh_325, \
                         shh_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_8 * shg0_235[k]
                   - f_9 * shg1_235[k]
                   + f_3 * pc_x[k] * shh_325[k];

        t_431[k] = f_3 * pc_z[k] * shh_321[k];

        t_432[k] = f_8 * shg0_237[k]
                   - f_9 * shg1_237[k]
                   + f_3 * pc_x[k] * shh_327[k];

        t_433[k] = f_0 * sgh_219[k]
                   + f_3 * pc_y[k] * shh_324[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, t_438, t_439, pc_x, shg0_239, shg1_239, \
                         shh_329, shh_330, shh_331, shh_332, shh_333, \
                         shh_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_8 * shg0_239[k]
                   - f_9 * shg1_239[k]
                   + f_3 * pc_x[k] * shh_329[k];

        t_435[k] = f_3 * pc_x[k] * shh_330[k];

        t_436[k] = f_3 * pc_x[k] * shh_331[k];

        t_437[k] = f_3 * pc_x[k] * shh_332[k];

        t_438[k] = f_3 * pc_x[k] * shh_333[k];

        t_439[k] = f_3 * pc_x[k] * shh_334[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, pc_x, pc_y, pc_z, sgh_225, sgh_227, \
                         shg0_235, shg0_237, shg1_235, shg1_237, shh_330, shh_332, \
                         shh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_3 * pc_x[k] * shh_335[k];

        t_441[k] = f_0 * sgh_225[k]
                   + f_1 * shg0_235[k]
                   - f_2 * shg1_235[k]
                   + f_3 * pc_y[k] * shh_330[k];

        t_442[k] = f_3 * pc_z[k] * shh_330[k];

        t_443[k] = f_0 * sgh_227[k]
                   + f_4 * shg0_237[k]
                   - f_5 * shg1_237[k]
                   + f_3 * pc_y[k] * shh_332[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, pc_y, pc_z, sgh_228, sgh_229, sgh_230, \
                         shg0_238, shg0_239, shg1_238, shg1_239, shh_333, shh_334, \
                         shh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_0 * sgh_228[k]
                   + f_6 * shg0_238[k]
                   - f_7 * shg1_238[k]
                   + f_3 * pc_y[k] * shh_333[k];

        t_445[k] = f_0 * sgh_229[k]
                   + f_8 * shg0_239[k]
                   - f_9 * shg1_239[k]
                   + f_3 * pc_y[k] * shh_334[k];

        t_446[k] = f_0 * sgh_230[k]
                   + f_3 * pc_y[k] * shh_335[k];

        t_447[k] = f_1 * shg0_239[k]
                   - f_2 * shg1_239[k]
                   + f_3 * pc_z[k] * shh_335[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, pb_z, pc_y, pc_z, sgi0_280, sgi0_283, \
                         sgh_210, sgh_231, sgi1_280, sgi1_283, \
                         shh_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = pb_z[k] * sgi0_280[k]
                   - f_10 * pc_z[k] * sgi1_280[k];

        t_449[k] = f_14 * sgh_231[k]
                   + f_3 * pc_y[k] * shh_336[k];

        t_450[k] = f_11 * sgh_210[k]
                   + f_3 * pc_z[k] * shh_336[k];

        t_451[k] = pb_z[k] * sgi0_283[k]
                   - f_10 * pc_z[k] * sgi1_283[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, pb_z, pc_x, pc_y, pc_z, sgi0_286, sgh_233, \
                         sgi1_286, shg0_245, shg1_245, shh_338, \
                         shh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_14 * sgh_233[k]
                   + f_3 * pc_y[k] * shh_338[k];

        t_453[k] = f_4 * shg0_245[k]
                   - f_5 * shg1_245[k]
                   + f_3 * pc_x[k] * shh_341[k];

        t_454[k] = pb_z[k] * sgi0_286[k]
                   - f_10 * pc_z[k] * sgi1_286[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pc_x, pc_y, pc_z, sgh_213, sgh_236, shg0_249, \
                         shg1_249, shh_339, shh_341, shh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_11 * sgh_213[k]
                   + f_3 * pc_z[k] * shh_339[k];

        t_456[k] = f_14 * sgh_236[k]
                   + f_3 * pc_y[k] * shh_341[k];

        t_457[k] = f_6 * shg0_249[k]
                   - f_7 * shg1_249[k]
                   + f_3 * pc_x[k] * shh_345[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, pb_z, pc_x, pc_z, sgi0_290, sgh_216, sgi1_290, \
                         shg0_252, shg1_252, shh_342, shh_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = pb_z[k] * sgi0_290[k]
                   - f_10 * pc_z[k] * sgi1_290[k];

        t_459[k] = f_11 * sgh_216[k]
                   + f_3 * pc_z[k] * shh_342[k];

        t_460[k] = f_8 * shg0_252[k]
                   - f_9 * shg1_252[k]
                   + f_3 * pc_x[k] * shh_348[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, t_465, pc_x, pc_y, sgh_240, shg0_254, \
                         shg1_254, shh_345, shh_350, shh_351, shh_352, \
                         shh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_14 * sgh_240[k]
                   + f_3 * pc_y[k] * shh_345[k];

        t_462[k] = f_8 * shg0_254[k]
                   - f_9 * shg1_254[k]
                   + f_3 * pc_x[k] * shh_350[k];

        t_463[k] = f_3 * pc_x[k] * shh_351[k];

        t_464[k] = f_3 * pc_x[k] * shh_352[k];

        t_465[k] = f_3 * pc_x[k] * shh_353[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, t_470, pb_z, pc_x, pc_z, sgi0_301, \
                         sgh_225, sgi1_301, shh_351, shh_354, shh_355, \
                         shh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_3 * pc_x[k] * shh_354[k];

        t_467[k] = f_3 * pc_x[k] * shh_355[k];

        t_468[k] = f_3 * pc_x[k] * shh_356[k];

        t_469[k] = pb_z[k] * sgi0_301[k]
                   - f_10 * pc_z[k] * sgi1_301[k];

        t_470[k] = f_11 * sgh_225[k]
                   + f_3 * pc_z[k] * shh_351[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pb_z, pc_z, sgi0_303, sgi0_304, sgi0_305, \
                         sgh_226, sgh_227, sgh_228, sgi1_303, sgi1_304, \
                         sgi1_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = pb_z[k] * sgi0_303[k]
                   + f_12 * sgh_226[k]
                   - f_10 * pc_z[k] * sgi1_303[k];

        t_472[k] = pb_z[k] * sgi0_304[k]
                   + f_13 * sgh_227[k]
                   - f_10 * pc_z[k] * sgi1_304[k];

        t_473[k] = pb_z[k] * sgi0_305[k]
                   + f_14 * sgh_228[k]
                   - f_10 * pc_z[k] * sgi1_305[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pc_x, pc_y, pc_z, sgh_230, sgh_251, \
                         sgh_252, shg0_254, shg0_255, shg1_254, shg1_255, shh_356, \
                         shh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_14 * sgh_251[k]
                   + f_3 * pc_y[k] * shh_356[k];

        t_475[k] = f_11 * sgh_230[k]
                   + f_1 * shg0_254[k]
                   - f_2 * shg1_254[k]
                   + f_3 * pc_z[k] * shh_356[k];

        t_476[k] = f_1 * shg0_255[k]
                   - f_2 * shg1_255[k]
                   + f_3 * pc_x[k] * shh_357[k];

        t_477[k] = f_13 * sgh_252[k]
                   + f_3 * pc_y[k] * shh_357[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pc_x, pc_y, pc_z, sgh_231, sgh_254, shg0_258, \
                         shg1_258, shh_357, shh_359, shh_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_12 * sgh_231[k]
                   + f_3 * pc_z[k] * shh_357[k];

        t_479[k] = f_4 * shg0_258[k]
                   - f_5 * shg1_258[k]
                   + f_3 * pc_x[k] * shh_360[k];

        t_480[k] = f_13 * sgh_254[k]
                   + f_3 * pc_y[k] * shh_359[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pc_x, pc_y, pc_z, sgh_234, sgh_257, \
                         shg0_260, shg0_261, shg1_260, shg1_261, shh_360, shh_362, \
                         shh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_4 * shg0_260[k]
                   - f_5 * shg1_260[k]
                   + f_3 * pc_x[k] * shh_362[k];

        t_482[k] = f_6 * shg0_261[k]
                   - f_7 * shg1_261[k]
                   + f_3 * pc_x[k] * shh_363[k];

        t_483[k] = f_12 * sgh_234[k]
                   + f_3 * pc_z[k] * shh_360[k];

        t_484[k] = f_13 * sgh_257[k]
                   + f_3 * pc_y[k] * shh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_z, sgh_237, shg0_264, shg0_265, \
                         shg1_264, shg1_265, shh_363, shh_366, \
                         shh_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_6 * shg0_264[k]
                   - f_7 * shg1_264[k]
                   + f_3 * pc_x[k] * shh_366[k];

        t_486[k] = f_8 * shg0_265[k]
                   - f_9 * shg1_265[k]
                   + f_3 * pc_x[k] * shh_367[k];

        t_487[k] = f_12 * sgh_237[k]
                   + f_3 * pc_z[k] * shh_363[k];
    }
}

static auto
compute_prim_shi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgi0,
                                                          const size_t sgh, const size_t sgi1,
                                                          const size_t shg0, const size_t shg1,
                                                          const size_t shh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgi0_392 = buffer.data(sgi0 + 392);
    const auto *sgi0_397 = buffer.data(sgi0 + 397);
    const auto *sgi0_401 = buffer.data(sgi0 + 401);
    const auto *sgi0_406 = buffer.data(sgi0 + 406);
    const auto *sgi0_413 = buffer.data(sgi0 + 413);
    const auto *sgi0_415 = buffer.data(sgi0 + 415);
    const auto *sgi0_416 = buffer.data(sgi0 + 416);
    const auto *sgi0_417 = buffer.data(sgi0 + 417);
    const auto *sgi0_419 = buffer.data(sgi0 + 419);

    const auto *sgh_246 = buffer.data(sgh + 246);
    const auto *sgh_251 = buffer.data(sgh + 251);
    const auto *sgh_252 = buffer.data(sgh + 252);
    const auto *sgh_255 = buffer.data(sgh + 255);
    const auto *sgh_258 = buffer.data(sgh + 258);
    const auto *sgh_261 = buffer.data(sgh + 261);
    const auto *sgh_267 = buffer.data(sgh + 267);
    const auto *sgh_269 = buffer.data(sgh + 269);
    const auto *sgh_270 = buffer.data(sgh + 270);
    const auto *sgh_271 = buffer.data(sgh + 271);
    const auto *sgh_272 = buffer.data(sgh + 272);
    const auto *sgh_273 = buffer.data(sgh + 273);
    const auto *sgh_275 = buffer.data(sgh + 275);
    const auto *sgh_276 = buffer.data(sgh + 276);
    const auto *sgh_278 = buffer.data(sgh + 278);
    const auto *sgh_279 = buffer.data(sgh + 279);
    const auto *sgh_282 = buffer.data(sgh + 282);
    const auto *sgh_288 = buffer.data(sgh + 288);
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
    const auto *sgh_309 = buffer.data(sgh + 309);
    const auto *sgh_311 = buffer.data(sgh + 311);
    const auto *sgh_312 = buffer.data(sgh + 312);
    const auto *sgh_313 = buffer.data(sgh + 313);
    const auto *sgh_314 = buffer.data(sgh + 314);

    const auto *sgi1_392 = buffer.data(sgi1 + 392);
    const auto *sgi1_397 = buffer.data(sgi1 + 397);
    const auto *sgi1_401 = buffer.data(sgi1 + 401);
    const auto *sgi1_406 = buffer.data(sgi1 + 406);
    const auto *sgi1_413 = buffer.data(sgi1 + 413);
    const auto *sgi1_415 = buffer.data(sgi1 + 415);
    const auto *sgi1_416 = buffer.data(sgi1 + 416);
    const auto *sgi1_417 = buffer.data(sgi1 + 417);
    const auto *sgi1_419 = buffer.data(sgi1 + 419);

    const auto *shg0_265 = buffer.data(shg0 + 265);
    const auto *shg0_267 = buffer.data(shg0 + 267);
    const auto *shg0_268 = buffer.data(shg0 + 268);
    const auto *shg0_269 = buffer.data(shg0 + 269);
    const auto *shg0_270 = buffer.data(shg0 + 270);
    const auto *shg0_273 = buffer.data(shg0 + 273);
    const auto *shg0_275 = buffer.data(shg0 + 275);
    const auto *shg0_276 = buffer.data(shg0 + 276);
    const auto *shg0_279 = buffer.data(shg0 + 279);
    const auto *shg0_280 = buffer.data(shg0 + 280);
    const auto *shg0_282 = buffer.data(shg0 + 282);
    const auto *shg0_283 = buffer.data(shg0 + 283);
    const auto *shg0_284 = buffer.data(shg0 + 284);
    const auto *shg0_288 = buffer.data(shg0 + 288);
    const auto *shg0_291 = buffer.data(shg0 + 291);
    const auto *shg0_295 = buffer.data(shg0 + 295);
    const auto *shg0_297 = buffer.data(shg0 + 297);
    const auto *shg0_300 = buffer.data(shg0 + 300);
    const auto *shg0_303 = buffer.data(shg0 + 303);
    const auto *shg0_305 = buffer.data(shg0 + 305);
    const auto *shg0_306 = buffer.data(shg0 + 306);
    const auto *shg0_309 = buffer.data(shg0 + 309);
    const auto *shg0_310 = buffer.data(shg0 + 310);
    const auto *shg0_312 = buffer.data(shg0 + 312);
    const auto *shg0_313 = buffer.data(shg0 + 313);
    const auto *shg0_314 = buffer.data(shg0 + 314);

    const auto *shg1_265 = buffer.data(shg1 + 265);
    const auto *shg1_267 = buffer.data(shg1 + 267);
    const auto *shg1_268 = buffer.data(shg1 + 268);
    const auto *shg1_269 = buffer.data(shg1 + 269);
    const auto *shg1_270 = buffer.data(shg1 + 270);
    const auto *shg1_273 = buffer.data(shg1 + 273);
    const auto *shg1_275 = buffer.data(shg1 + 275);
    const auto *shg1_276 = buffer.data(shg1 + 276);
    const auto *shg1_279 = buffer.data(shg1 + 279);
    const auto *shg1_280 = buffer.data(shg1 + 280);
    const auto *shg1_282 = buffer.data(shg1 + 282);
    const auto *shg1_283 = buffer.data(shg1 + 283);
    const auto *shg1_284 = buffer.data(shg1 + 284);
    const auto *shg1_288 = buffer.data(shg1 + 288);
    const auto *shg1_291 = buffer.data(shg1 + 291);
    const auto *shg1_295 = buffer.data(shg1 + 295);
    const auto *shg1_297 = buffer.data(shg1 + 297);
    const auto *shg1_300 = buffer.data(shg1 + 300);
    const auto *shg1_303 = buffer.data(shg1 + 303);
    const auto *shg1_305 = buffer.data(shg1 + 305);
    const auto *shg1_306 = buffer.data(shg1 + 306);
    const auto *shg1_309 = buffer.data(shg1 + 309);
    const auto *shg1_310 = buffer.data(shg1 + 310);
    const auto *shg1_312 = buffer.data(shg1 + 312);
    const auto *shg1_313 = buffer.data(shg1 + 313);
    const auto *shg1_314 = buffer.data(shg1 + 314);

    const auto *shh_366 = buffer.data(shh + 366);
    const auto *shh_369 = buffer.data(shh + 369);
    const auto *shh_371 = buffer.data(shh + 371);
    const auto *shh_372 = buffer.data(shh + 372);
    const auto *shh_373 = buffer.data(shh + 373);
    const auto *shh_374 = buffer.data(shh + 374);
    const auto *shh_375 = buffer.data(shh + 375);
    const auto *shh_376 = buffer.data(shh + 376);
    const auto *shh_377 = buffer.data(shh + 377);
    const auto *shh_378 = buffer.data(shh + 378);
    const auto *shh_380 = buffer.data(shh + 380);
    const auto *shh_381 = buffer.data(shh + 381);
    const auto *shh_383 = buffer.data(shh + 383);
    const auto *shh_384 = buffer.data(shh + 384);
    const auto *shh_387 = buffer.data(shh + 387);
    const auto *shh_388 = buffer.data(shh + 388);
    const auto *shh_390 = buffer.data(shh + 390);
    const auto *shh_392 = buffer.data(shh + 392);
    const auto *shh_393 = buffer.data(shh + 393);
    const auto *shh_394 = buffer.data(shh + 394);
    const auto *shh_395 = buffer.data(shh + 395);
    const auto *shh_396 = buffer.data(shh + 396);
    const auto *shh_397 = buffer.data(shh + 397);
    const auto *shh_398 = buffer.data(shh + 398);
    const auto *shh_399 = buffer.data(shh + 399);
    const auto *shh_401 = buffer.data(shh + 401);
    const auto *shh_402 = buffer.data(shh + 402);
    const auto *shh_404 = buffer.data(shh + 404);
    const auto *shh_405 = buffer.data(shh + 405);
    const auto *shh_408 = buffer.data(shh + 408);
    const auto *shh_409 = buffer.data(shh + 409);
    const auto *shh_411 = buffer.data(shh + 411);
    const auto *shh_414 = buffer.data(shh + 414);
    const auto *shh_415 = buffer.data(shh + 415);
    const auto *shh_416 = buffer.data(shh + 416);
    const auto *shh_417 = buffer.data(shh + 417);
    const auto *shh_418 = buffer.data(shh + 418);
    const auto *shh_419 = buffer.data(shh + 419);
    const auto *shh_420 = buffer.data(shh + 420);
    const auto *shh_422 = buffer.data(shh + 422);
    const auto *shh_423 = buffer.data(shh + 423);
    const auto *shh_425 = buffer.data(shh + 425);
    const auto *shh_426 = buffer.data(shh + 426);
    const auto *shh_429 = buffer.data(shh + 429);
    const auto *shh_430 = buffer.data(shh + 430);
    const auto *shh_432 = buffer.data(shh + 432);
    const auto *shh_434 = buffer.data(shh + 434);
    const auto *shh_435 = buffer.data(shh + 435);
    const auto *shh_436 = buffer.data(shh + 436);
    const auto *shh_437 = buffer.data(shh + 437);
    const auto *shh_438 = buffer.data(shh + 438);
    const auto *shh_439 = buffer.data(shh + 439);
    const auto *shh_440 = buffer.data(shh + 440);

#pragma omp simd aligned(t_488, t_489, t_490, t_491, pc_x, pc_y, sgh_261, shg0_267, shg0_269, \
                         shg1_267, shg1_269, shh_366, shh_369, shh_371, \
                         shh_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_8 * shg0_267[k]
                   - f_9 * shg1_267[k]
                   + f_3 * pc_x[k] * shh_369[k];

        t_489[k] = f_13 * sgh_261[k]
                   + f_3 * pc_y[k] * shh_366[k];

        t_490[k] = f_8 * shg0_269[k]
                   - f_9 * shg1_269[k]
                   + f_3 * pc_x[k] * shh_371[k];

        t_491[k] = f_3 * pc_x[k] * shh_372[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, pc_x, shh_373, shh_374, shh_375, \
                         shh_376, shh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_3 * pc_x[k] * shh_373[k];

        t_493[k] = f_3 * pc_x[k] * shh_374[k];

        t_494[k] = f_3 * pc_x[k] * shh_375[k];

        t_495[k] = f_3 * pc_x[k] * shh_376[k];

        t_496[k] = f_3 * pc_x[k] * shh_377[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, pc_y, pc_z, sgh_246, sgh_267, sgh_269, shg0_265, \
                         shg0_267, shg1_265, shg1_267, shh_372, \
                         shh_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = f_13 * sgh_267[k]
                   + f_1 * shg0_265[k]
                   - f_2 * shg1_265[k]
                   + f_3 * pc_y[k] * shh_372[k];

        t_498[k] = f_12 * sgh_246[k]
                   + f_3 * pc_z[k] * shh_372[k];

        t_499[k] = f_13 * sgh_269[k]
                   + f_4 * shg0_267[k]
                   - f_5 * shg1_267[k]
                   + f_3 * pc_y[k] * shh_374[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, pc_y, sgh_270, sgh_271, sgh_272, shg0_268, \
                         shg0_269, shg1_268, shg1_269, shh_375, shh_376, \
                         shh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_13 * sgh_270[k]
                   + f_6 * shg0_268[k]
                   - f_7 * shg1_268[k]
                   + f_3 * pc_y[k] * shh_375[k];

        t_501[k] = f_13 * sgh_271[k]
                   + f_8 * shg0_269[k]
                   - f_9 * shg1_269[k]
                   + f_3 * pc_y[k] * shh_376[k];

        t_502[k] = f_13 * sgh_272[k]
                   + f_3 * pc_y[k] * shh_377[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, t_506, pc_x, pc_y, pc_z, sgh_251, sgh_252, \
                         sgh_273, shg0_269, shg0_270, shg1_269, shg1_270, shh_377, \
                         shh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_12 * sgh_251[k]
                   + f_1 * shg0_269[k]
                   - f_2 * shg1_269[k]
                   + f_3 * pc_z[k] * shh_377[k];

        t_504[k] = f_1 * shg0_270[k]
                   - f_2 * shg1_270[k]
                   + f_3 * pc_x[k] * shh_378[k];

        t_505[k] = f_12 * sgh_273[k]
                   + f_3 * pc_y[k] * shh_378[k];

        t_506[k] = f_13 * sgh_252[k]
                   + f_3 * pc_z[k] * shh_378[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, pc_x, pc_y, sgh_275, shg0_273, shg0_275, \
                         shg1_273, shg1_275, shh_380, shh_381, \
                         shh_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_4 * shg0_273[k]
                   - f_5 * shg1_273[k]
                   + f_3 * pc_x[k] * shh_381[k];

        t_508[k] = f_12 * sgh_275[k]
                   + f_3 * pc_y[k] * shh_380[k];

        t_509[k] = f_4 * shg0_275[k]
                   - f_5 * shg1_275[k]
                   + f_3 * pc_x[k] * shh_383[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, pc_x, pc_y, pc_z, sgh_255, sgh_278, shg0_276, \
                         shg1_276, shh_381, shh_383, shh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = f_6 * shg0_276[k]
                   - f_7 * shg1_276[k]
                   + f_3 * pc_x[k] * shh_384[k];

        t_511[k] = f_13 * sgh_255[k]
                   + f_3 * pc_z[k] * shh_381[k];

        t_512[k] = f_12 * sgh_278[k]
                   + f_3 * pc_y[k] * shh_383[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, pc_x, pc_z, sgh_258, shg0_279, shg0_280, \
                         shg1_279, shg1_280, shh_384, shh_387, \
                         shh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = f_6 * shg0_279[k]
                   - f_7 * shg1_279[k]
                   + f_3 * pc_x[k] * shh_387[k];

        t_514[k] = f_8 * shg0_280[k]
                   - f_9 * shg1_280[k]
                   + f_3 * pc_x[k] * shh_388[k];

        t_515[k] = f_13 * sgh_258[k]
                   + f_3 * pc_z[k] * shh_384[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, pc_y, sgh_282, shg0_282, shg0_284, \
                         shg1_282, shg1_284, shh_387, shh_390, shh_392, \
                         shh_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_8 * shg0_282[k]
                   - f_9 * shg1_282[k]
                   + f_3 * pc_x[k] * shh_390[k];

        t_517[k] = f_12 * sgh_282[k]
                   + f_3 * pc_y[k] * shh_387[k];

        t_518[k] = f_8 * shg0_284[k]
                   - f_9 * shg1_284[k]
                   + f_3 * pc_x[k] * shh_392[k];

        t_519[k] = f_3 * pc_x[k] * shh_393[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, pc_x, shh_394, shh_395, shh_396, \
                         shh_397, shh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_3 * pc_x[k] * shh_394[k];

        t_521[k] = f_3 * pc_x[k] * shh_395[k];

        t_522[k] = f_3 * pc_x[k] * shh_396[k];

        t_523[k] = f_3 * pc_x[k] * shh_397[k];

        t_524[k] = f_3 * pc_x[k] * shh_398[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, pc_y, pc_z, sgh_267, sgh_288, sgh_290, shg0_280, \
                         shg0_282, shg1_280, shg1_282, shh_393, \
                         shh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_12 * sgh_288[k]
                   + f_1 * shg0_280[k]
                   - f_2 * shg1_280[k]
                   + f_3 * pc_y[k] * shh_393[k];

        t_526[k] = f_13 * sgh_267[k]
                   + f_3 * pc_z[k] * shh_393[k];

        t_527[k] = f_12 * sgh_290[k]
                   + f_4 * shg0_282[k]
                   - f_5 * shg1_282[k]
                   + f_3 * pc_y[k] * shh_395[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, pc_y, sgh_291, sgh_292, sgh_293, shg0_283, \
                         shg0_284, shg1_283, shg1_284, shh_396, shh_397, \
                         shh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = f_12 * sgh_291[k]
                   + f_6 * shg0_283[k]
                   - f_7 * shg1_283[k]
                   + f_3 * pc_y[k] * shh_396[k];

        t_529[k] = f_12 * sgh_292[k]
                   + f_8 * shg0_284[k]
                   - f_9 * shg1_284[k]
                   + f_3 * pc_y[k] * shh_397[k];

        t_530[k] = f_12 * sgh_293[k]
                   + f_3 * pc_y[k] * shh_398[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, pb_y, pc_y, pc_z, sgi0_392, sgh_272, \
                         sgh_273, sgh_294, sgi1_392, shg0_284, shg1_284, shh_398, \
                         shh_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = f_13 * sgh_272[k]
                   + f_1 * shg0_284[k]
                   - f_2 * shg1_284[k]
                   + f_3 * pc_z[k] * shh_398[k];

        t_532[k] = pb_y[k] * sgi0_392[k]
                   - f_10 * pc_y[k] * sgi1_392[k];

        t_533[k] = f_11 * sgh_294[k]
                   + f_3 * pc_y[k] * shh_399[k];

        t_534[k] = f_14 * sgh_273[k]
                   + f_3 * pc_z[k] * shh_399[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, pb_y, pc_x, pc_y, sgi0_397, sgh_296, sgi1_397, \
                         shg0_288, shg1_288, shh_401, shh_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = f_4 * shg0_288[k]
                   - f_5 * shg1_288[k]
                   + f_3 * pc_x[k] * shh_402[k];

        t_536[k] = f_11 * sgh_296[k]
                   + f_3 * pc_y[k] * shh_401[k];

        t_537[k] = pb_y[k] * sgi0_397[k]
                   - f_10 * pc_y[k] * sgi1_397[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, pc_x, pc_y, pc_z, sgh_276, sgh_299, shg0_291, \
                         shg1_291, shh_402, shh_404, shh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_6 * shg0_291[k]
                   - f_7 * shg1_291[k]
                   + f_3 * pc_x[k] * shh_405[k];

        t_539[k] = f_14 * sgh_276[k]
                   + f_3 * pc_z[k] * shh_402[k];

        t_540[k] = f_11 * sgh_299[k]
                   + f_3 * pc_y[k] * shh_404[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pb_y, pc_x, pc_y, pc_z, sgi0_401, sgh_279, \
                         sgi1_401, shg0_295, shg1_295, shh_405, \
                         shh_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = pb_y[k] * sgi0_401[k]
                   - f_10 * pc_y[k] * sgi1_401[k];

        t_542[k] = f_8 * shg0_295[k]
                   - f_9 * shg1_295[k]
                   + f_3 * pc_x[k] * shh_409[k];

        t_543[k] = f_14 * sgh_279[k]
                   + f_3 * pc_z[k] * shh_405[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pb_y, pc_x, pc_y, sgi0_406, sgh_303, \
                         sgi1_406, shg0_297, shg1_297, shh_408, shh_411, \
                         shh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_8 * shg0_297[k]
                   - f_9 * shg1_297[k]
                   + f_3 * pc_x[k] * shh_411[k];

        t_545[k] = f_11 * sgh_303[k]
                   + f_3 * pc_y[k] * shh_408[k];

        t_546[k] = pb_y[k] * sgi0_406[k]
                   - f_10 * pc_y[k] * sgi1_406[k];

        t_547[k] = f_3 * pc_x[k] * shh_414[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pc_x, shh_415, shh_416, shh_417, \
                         shh_418, shh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_3 * pc_x[k] * shh_415[k];

        t_549[k] = f_3 * pc_x[k] * shh_416[k];

        t_550[k] = f_3 * pc_x[k] * shh_417[k];

        t_551[k] = f_3 * pc_x[k] * shh_418[k];

        t_552[k] = f_3 * pc_x[k] * shh_419[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pb_y, pc_y, pc_z, sgi0_413, sgi0_415, sgh_288, \
                         sgh_309, sgh_311, sgi1_413, sgi1_415, \
                         shh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = pb_y[k] * sgi0_413[k]
                   + f_15 * sgh_309[k]
                   - f_10 * pc_y[k] * sgi1_413[k];

        t_554[k] = f_14 * sgh_288[k]
                   + f_3 * pc_z[k] * shh_414[k];

        t_555[k] = pb_y[k] * sgi0_415[k]
                   + f_14 * sgh_311[k]
                   - f_10 * pc_y[k] * sgi1_415[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, pb_y, pc_y, sgi0_416, sgi0_417, sgi0_419, \
                         sgh_312, sgh_313, sgh_314, sgi1_416, sgi1_417, sgi1_419, \
                         shh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = pb_y[k] * sgi0_416[k]
                   + f_13 * sgh_312[k]
                   - f_10 * pc_y[k] * sgi1_416[k];

        t_557[k] = pb_y[k] * sgi0_417[k]
                   + f_12 * sgh_313[k]
                   - f_10 * pc_y[k] * sgi1_417[k];

        t_558[k] = f_11 * sgh_314[k]
                   + f_3 * pc_y[k] * shh_419[k];

        t_559[k] = pb_y[k] * sgi0_419[k]
                   - f_10 * pc_y[k] * sgi1_419[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, pc_x, pc_y, pc_z, sgh_294, \
                         shg0_300, shg0_303, shg1_300, shg1_303, shh_420, shh_422, \
                         shh_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_1 * shg0_300[k]
                   - f_2 * shg1_300[k]
                   + f_3 * pc_x[k] * shh_420[k];

        t_561[k] = f_3 * pc_y[k] * shh_420[k];

        t_562[k] = f_0 * sgh_294[k]
                   + f_3 * pc_z[k] * shh_420[k];

        t_563[k] = f_4 * shg0_303[k]
                   - f_5 * shg1_303[k]
                   + f_3 * pc_x[k] * shh_423[k];

        t_564[k] = f_3 * pc_y[k] * shh_422[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, pc_x, pc_y, pc_z, sgh_297, shg0_305, \
                         shg0_306, shg1_305, shg1_306, shh_423, shh_425, \
                         shh_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = f_4 * shg0_305[k]
                   - f_5 * shg1_305[k]
                   + f_3 * pc_x[k] * shh_425[k];

        t_566[k] = f_6 * shg0_306[k]
                   - f_7 * shg1_306[k]
                   + f_3 * pc_x[k] * shh_426[k];

        t_567[k] = f_0 * sgh_297[k]
                   + f_3 * pc_z[k] * shh_423[k];

        t_568[k] = f_3 * pc_y[k] * shh_425[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, pc_z, sgh_300, shg0_309, shg0_310, \
                         shg1_309, shg1_310, shh_426, shh_429, \
                         shh_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_6 * shg0_309[k]
                   - f_7 * shg1_309[k]
                   + f_3 * pc_x[k] * shh_429[k];

        t_570[k] = f_8 * shg0_310[k]
                   - f_9 * shg1_310[k]
                   + f_3 * pc_x[k] * shh_430[k];

        t_571[k] = f_0 * sgh_300[k]
                   + f_3 * pc_z[k] * shh_426[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, pc_x, pc_y, shg0_312, shg0_314, \
                         shg1_312, shg1_314, shh_429, shh_432, shh_434, shh_435, \
                         shh_436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_8 * shg0_312[k]
                   - f_9 * shg1_312[k]
                   + f_3 * pc_x[k] * shh_432[k];

        t_573[k] = f_3 * pc_y[k] * shh_429[k];

        t_574[k] = f_8 * shg0_314[k]
                   - f_9 * shg1_314[k]
                   + f_3 * pc_x[k] * shh_434[k];

        t_575[k] = f_3 * pc_x[k] * shh_435[k];

        t_576[k] = f_3 * pc_x[k] * shh_436[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, pc_x, pc_y, shg0_310, shg1_310, \
                         shh_435, shh_437, shh_438, shh_439, shh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = f_3 * pc_x[k] * shh_437[k];

        t_578[k] = f_3 * pc_x[k] * shh_438[k];

        t_579[k] = f_3 * pc_x[k] * shh_439[k];

        t_580[k] = f_3 * pc_x[k] * shh_440[k];

        t_581[k] = f_1 * shg0_310[k]
                   - f_2 * shg1_310[k]
                   + f_3 * pc_y[k] * shh_435[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, pc_y, pc_z, sgh_309, shg0_312, shg0_313, \
                         shg1_312, shg1_313, shh_435, shh_437, \
                         shh_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_0 * sgh_309[k]
                   + f_3 * pc_z[k] * shh_435[k];

        t_583[k] = f_4 * shg0_312[k]
                   - f_5 * shg1_312[k]
                   + f_3 * pc_y[k] * shh_437[k];

        t_584[k] = f_6 * shg0_313[k]
                   - f_7 * shg1_313[k]
                   + f_3 * pc_y[k] * shh_438[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, pc_y, pc_z, sgh_314, shg0_314, shg1_314, \
                         shh_439, shh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_8 * shg0_314[k]
                   - f_9 * shg1_314[k]
                   + f_3 * pc_y[k] * shh_439[k];

        t_586[k] = f_3 * pc_y[k] * shh_440[k];

        t_587[k] = f_0 * sgh_314[k]
                   + f_1 * shg0_314[k]
                   - f_2 * shg1_314[k]
                   + f_3 * pc_z[k] * shh_440[k];
    }
}

auto
compute_prim_shi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sgi0, const size_t sgh,
                                                   const size_t sgi1, const size_t shg0,
                                                   const size_t shg1, const size_t shh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_shi_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sgi0, sgh,
                                                              sgi1, shg0, shg1, shh, ncols,
                                                              gamma, p, q);

    compute_prim_shi_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sgi0, sgh,
                                                              sgi1, shg0, shg1, shh, ncols,
                                                              gamma, p, q);

    compute_prim_shi_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sgi0, sgh,
                                                              sgi1, shg0, shg1, shh, ncols,
                                                              gamma, p, q);

    compute_prim_shi_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sgi0, sgh,
                                                              sgi1, shg0, shg1, shh, ncols,
                                                              gamma, p, q);

    compute_prim_shi_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sgi0, sgh,
                                                              sgi1, shg0, shg1, shh, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
