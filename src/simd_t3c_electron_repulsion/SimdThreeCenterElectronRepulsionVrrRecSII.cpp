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


#include "SimdThreeCenterElectronRepulsionVrrRecSII.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sii_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shi0,
                                                          const size_t shh, const size_t shi1,
                                                          const size_t sig0, const size_t sig1,
                                                          const size_t sih, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_15 = 2.5 / q;

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

    const auto *shi0_0 = buffer.data(shi0 + 0);
    const auto *shi0_3 = buffer.data(shi0 + 3);
    const auto *shi0_5 = buffer.data(shi0 + 5);
    const auto *shi0_6 = buffer.data(shi0 + 6);
    const auto *shi0_9 = buffer.data(shi0 + 9);
    const auto *shi0_10 = buffer.data(shi0 + 10);
    const auto *shi0_12 = buffer.data(shi0 + 12);
    const auto *shi0_14 = buffer.data(shi0 + 14);
    const auto *shi0_21 = buffer.data(shi0 + 21);
    const auto *shi0_27 = buffer.data(shi0 + 27);
    const auto *shi0_31 = buffer.data(shi0 + 31);
    const auto *shi0_34 = buffer.data(shi0 + 34);
    const auto *shi0_38 = buffer.data(shi0 + 38);
    const auto *shi0_56 = buffer.data(shi0 + 56);
    const auto *shi0_61 = buffer.data(shi0 + 61);
    const auto *shi0_65 = buffer.data(shi0 + 65);

    const auto *shh_0 = buffer.data(shh + 0);
    const auto *shh_1 = buffer.data(shh + 1);
    const auto *shh_2 = buffer.data(shh + 2);
    const auto *shh_3 = buffer.data(shh + 3);
    const auto *shh_5 = buffer.data(shh + 5);
    const auto *shh_6 = buffer.data(shh + 6);
    const auto *shh_7 = buffer.data(shh + 7);
    const auto *shh_8 = buffer.data(shh + 8);
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
    const auto *shh_47 = buffer.data(shh + 47);
    const auto *shh_57 = buffer.data(shh + 57);
    const auto *shh_58 = buffer.data(shh + 58);
    const auto *shh_59 = buffer.data(shh + 59);
    const auto *shh_60 = buffer.data(shh + 60);
    const auto *shh_61 = buffer.data(shh + 61);
    const auto *shh_62 = buffer.data(shh + 62);
    const auto *shh_63 = buffer.data(shh + 63);
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

    const auto *shi1_0 = buffer.data(shi1 + 0);
    const auto *shi1_3 = buffer.data(shi1 + 3);
    const auto *shi1_5 = buffer.data(shi1 + 5);
    const auto *shi1_6 = buffer.data(shi1 + 6);
    const auto *shi1_9 = buffer.data(shi1 + 9);
    const auto *shi1_10 = buffer.data(shi1 + 10);
    const auto *shi1_12 = buffer.data(shi1 + 12);
    const auto *shi1_14 = buffer.data(shi1 + 14);
    const auto *shi1_21 = buffer.data(shi1 + 21);
    const auto *shi1_27 = buffer.data(shi1 + 27);
    const auto *shi1_31 = buffer.data(shi1 + 31);
    const auto *shi1_34 = buffer.data(shi1 + 34);
    const auto *shi1_38 = buffer.data(shi1 + 38);
    const auto *shi1_56 = buffer.data(shi1 + 56);
    const auto *shi1_61 = buffer.data(shi1 + 61);
    const auto *shi1_65 = buffer.data(shi1 + 65);

    const auto *sig0_0 = buffer.data(sig0 + 0);
    const auto *sig0_3 = buffer.data(sig0 + 3);
    const auto *sig0_5 = buffer.data(sig0 + 5);
    const auto *sig0_6 = buffer.data(sig0 + 6);
    const auto *sig0_9 = buffer.data(sig0 + 9);
    const auto *sig0_10 = buffer.data(sig0 + 10);
    const auto *sig0_12 = buffer.data(sig0 + 12);
    const auto *sig0_13 = buffer.data(sig0 + 13);
    const auto *sig0_14 = buffer.data(sig0 + 14);
    const auto *sig0_25 = buffer.data(sig0 + 25);
    const auto *sig0_27 = buffer.data(sig0 + 27);
    const auto *sig0_28 = buffer.data(sig0 + 28);
    const auto *sig0_29 = buffer.data(sig0 + 29);
    const auto *sig0_42 = buffer.data(sig0 + 42);
    const auto *sig0_43 = buffer.data(sig0 + 43);
    const auto *sig0_44 = buffer.data(sig0 + 44);
    const auto *sig0_45 = buffer.data(sig0 + 45);
    const auto *sig0_48 = buffer.data(sig0 + 48);
    const auto *sig0_50 = buffer.data(sig0 + 50);
    const auto *sig0_51 = buffer.data(sig0 + 51);
    const auto *sig0_54 = buffer.data(sig0 + 54);
    const auto *sig0_55 = buffer.data(sig0 + 55);
    const auto *sig0_57 = buffer.data(sig0 + 57);
    const auto *sig0_58 = buffer.data(sig0 + 58);
    const auto *sig0_59 = buffer.data(sig0 + 59);

    const auto *sig1_0 = buffer.data(sig1 + 0);
    const auto *sig1_3 = buffer.data(sig1 + 3);
    const auto *sig1_5 = buffer.data(sig1 + 5);
    const auto *sig1_6 = buffer.data(sig1 + 6);
    const auto *sig1_9 = buffer.data(sig1 + 9);
    const auto *sig1_10 = buffer.data(sig1 + 10);
    const auto *sig1_12 = buffer.data(sig1 + 12);
    const auto *sig1_13 = buffer.data(sig1 + 13);
    const auto *sig1_14 = buffer.data(sig1 + 14);
    const auto *sig1_25 = buffer.data(sig1 + 25);
    const auto *sig1_27 = buffer.data(sig1 + 27);
    const auto *sig1_28 = buffer.data(sig1 + 28);
    const auto *sig1_29 = buffer.data(sig1 + 29);
    const auto *sig1_42 = buffer.data(sig1 + 42);
    const auto *sig1_43 = buffer.data(sig1 + 43);
    const auto *sig1_44 = buffer.data(sig1 + 44);
    const auto *sig1_45 = buffer.data(sig1 + 45);
    const auto *sig1_48 = buffer.data(sig1 + 48);
    const auto *sig1_50 = buffer.data(sig1 + 50);
    const auto *sig1_51 = buffer.data(sig1 + 51);
    const auto *sig1_54 = buffer.data(sig1 + 54);
    const auto *sig1_55 = buffer.data(sig1 + 55);
    const auto *sig1_57 = buffer.data(sig1 + 57);
    const auto *sig1_58 = buffer.data(sig1 + 58);
    const auto *sig1_59 = buffer.data(sig1 + 59);

    const auto *sih_0 = buffer.data(sih + 0);
    const auto *sih_2 = buffer.data(sih + 2);
    const auto *sih_3 = buffer.data(sih + 3);
    const auto *sih_5 = buffer.data(sih + 5);
    const auto *sih_6 = buffer.data(sih + 6);
    const auto *sih_9 = buffer.data(sih + 9);
    const auto *sih_10 = buffer.data(sih + 10);
    const auto *sih_12 = buffer.data(sih + 12);
    const auto *sih_14 = buffer.data(sih + 14);
    const auto *sih_15 = buffer.data(sih + 15);
    const auto *sih_16 = buffer.data(sih + 16);
    const auto *sih_17 = buffer.data(sih + 17);
    const auto *sih_18 = buffer.data(sih + 18);
    const auto *sih_19 = buffer.data(sih + 19);
    const auto *sih_20 = buffer.data(sih + 20);
    const auto *sih_21 = buffer.data(sih + 21);
    const auto *sih_23 = buffer.data(sih + 23);
    const auto *sih_24 = buffer.data(sih + 24);
    const auto *sih_26 = buffer.data(sih + 26);
    const auto *sih_27 = buffer.data(sih + 27);
    const auto *sih_30 = buffer.data(sih + 30);
    const auto *sih_36 = buffer.data(sih + 36);
    const auto *sih_37 = buffer.data(sih + 37);
    const auto *sih_38 = buffer.data(sih + 38);
    const auto *sih_39 = buffer.data(sih + 39);
    const auto *sih_40 = buffer.data(sih + 40);
    const auto *sih_41 = buffer.data(sih + 41);
    const auto *sih_42 = buffer.data(sih + 42);
    const auto *sih_44 = buffer.data(sih + 44);
    const auto *sih_45 = buffer.data(sih + 45);
    const auto *sih_47 = buffer.data(sih + 47);
    const auto *sih_48 = buffer.data(sih + 48);
    const auto *sih_51 = buffer.data(sih + 51);
    const auto *sih_57 = buffer.data(sih + 57);
    const auto *sih_58 = buffer.data(sih + 58);
    const auto *sih_59 = buffer.data(sih + 59);
    const auto *sih_60 = buffer.data(sih + 60);
    const auto *sih_61 = buffer.data(sih + 61);
    const auto *sih_62 = buffer.data(sih + 62);
    const auto *sih_63 = buffer.data(sih + 63);
    const auto *sih_65 = buffer.data(sih + 65);
    const auto *sih_66 = buffer.data(sih + 66);
    const auto *sih_68 = buffer.data(sih + 68);
    const auto *sih_69 = buffer.data(sih + 69);
    const auto *sih_72 = buffer.data(sih + 72);
    const auto *sih_73 = buffer.data(sih + 73);
    const auto *sih_75 = buffer.data(sih + 75);
    const auto *sih_77 = buffer.data(sih + 77);
    const auto *sih_78 = buffer.data(sih + 78);
    const auto *sih_79 = buffer.data(sih + 79);
    const auto *sih_80 = buffer.data(sih + 80);
    const auto *sih_81 = buffer.data(sih + 81);
    const auto *sih_82 = buffer.data(sih + 82);
    const auto *sih_83 = buffer.data(sih + 83);
    const auto *sih_84 = buffer.data(sih + 84);
    const auto *sih_86 = buffer.data(sih + 86);
    const auto *sih_87 = buffer.data(sih + 87);
    const auto *sih_89 = buffer.data(sih + 89);
    const auto *sih_90 = buffer.data(sih + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, shh_0, shh_3, sig0_0, sig0_3, \
                         sig1_0, sig1_3, sih_0, sih_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * shh_0[k]
                 + f_1 * sig0_0[k]
                 - f_2 * sig1_0[k]
                 + f_3 * pc_x[k] * sih_0[k];

        t_1[k] = f_3 * pc_y[k] * sih_0[k];

        t_2[k] = f_3 * pc_z[k] * sih_0[k];

        t_3[k] = f_0 * shh_3[k]
                 + f_4 * sig0_3[k]
                 - f_5 * sig1_3[k]
                 + f_3 * pc_x[k] * sih_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, shh_5, shh_6, sig0_5, sig0_6, sig1_5, \
                         sig1_6, sih_2, sih_5, sih_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sih_2[k];

        t_5[k] = f_0 * shh_5[k]
                 + f_4 * sig0_5[k]
                 - f_5 * sig1_5[k]
                 + f_3 * pc_x[k] * sih_5[k];

        t_6[k] = f_0 * shh_6[k]
                 + f_6 * sig0_6[k]
                 - f_7 * sig1_6[k]
                 + f_3 * pc_x[k] * sih_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, shh_9, sig0_9, sig1_9, sih_3, sih_5, \
                         sih_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sih_3[k];

        t_8[k] = f_3 * pc_y[k] * sih_5[k];

        t_9[k] = f_0 * shh_9[k]
                 + f_6 * sig0_9[k]
                 - f_7 * sig1_9[k]
                 + f_3 * pc_x[k] * sih_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, shh_10, shh_12, sig0_10, sig0_12, \
                         sig1_10, sig1_12, sih_6, sih_10, sih_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * shh_10[k]
                  + f_8 * sig0_10[k]
                  - f_9 * sig1_10[k]
                  + f_3 * pc_x[k] * sih_10[k];

        t_11[k] = f_3 * pc_z[k] * sih_6[k];

        t_12[k] = f_0 * shh_12[k]
                  + f_8 * sig0_12[k]
                  - f_9 * sig1_12[k]
                  + f_3 * pc_x[k] * sih_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pc_x, pc_y, shh_14, shh_15, shh_16, sig0_14, \
                         sig1_14, sih_9, sih_14, sih_15, sih_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sih_9[k];

        t_14[k] = f_0 * shh_14[k]
                  + f_8 * sig0_14[k]
                  - f_9 * sig1_14[k]
                  + f_3 * pc_x[k] * sih_14[k];

        t_15[k] = f_0 * shh_15[k]
                  + f_3 * pc_x[k] * sih_15[k];

        t_16[k] = f_0 * shh_16[k]
                  + f_3 * pc_x[k] * sih_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_x, shh_17, shh_18, shh_19, shh_20, sih_17, \
                         sih_18, sih_19, sih_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * shh_17[k]
                  + f_3 * pc_x[k] * sih_17[k];

        t_18[k] = f_0 * shh_18[k]
                  + f_3 * pc_x[k] * sih_18[k];

        t_19[k] = f_0 * shh_19[k]
                  + f_3 * pc_x[k] * sih_19[k];

        t_20[k] = f_0 * shh_20[k]
                  + f_3 * pc_x[k] * sih_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, sig0_10, sig0_12, sig0_13, \
                         sig1_10, sig1_12, sig1_13, sih_15, sih_17, \
                         sih_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * sig0_10[k]
                  - f_2 * sig1_10[k]
                  + f_3 * pc_y[k] * sih_15[k];

        t_22[k] = f_3 * pc_z[k] * sih_15[k];

        t_23[k] = f_4 * sig0_12[k]
                  - f_5 * sig1_12[k]
                  + f_3 * pc_y[k] * sih_17[k];

        t_24[k] = f_6 * sig0_13[k]
                  - f_7 * sig1_13[k]
                  + f_3 * pc_y[k] * sih_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_y, pc_y, pc_z, shi0_0, shh_0, \
                         shi1_0, sig0_14, sig1_14, sih_19, sih_20, \
                         sih_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * sig0_14[k]
                  - f_9 * sig1_14[k]
                  + f_3 * pc_y[k] * sih_19[k];

        t_26[k] = f_3 * pc_y[k] * sih_20[k];

        t_27[k] = f_1 * sig0_14[k]
                  - f_2 * sig1_14[k]
                  + f_3 * pc_z[k] * sih_20[k];

        t_28[k] = pb_y[k] * shi0_0[k]
                  - f_10 * pc_y[k] * shi1_0[k];

        t_29[k] = f_11 * shh_0[k]
                  + f_3 * pc_y[k] * sih_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_y, pc_y, pc_z, shi0_3, shi0_5, shh_1, \
                         shh_2, shi1_3, shi1_5, sih_21, sih_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * sih_21[k];

        t_31[k] = pb_y[k] * shi0_3[k]
                  + f_12 * shh_1[k]
                  - f_10 * pc_y[k] * shi1_3[k];

        t_32[k] = f_11 * shh_2[k]
                  + f_3 * pc_y[k] * sih_23[k];

        t_33[k] = pb_y[k] * shi0_5[k]
                  - f_10 * pc_y[k] * shi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_y, pc_y, pc_z, shi0_6, shi0_9, shh_3, \
                         shh_5, shi1_6, shi1_9, sih_24, sih_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * shi0_6[k]
                  + f_13 * shh_3[k]
                  - f_10 * pc_y[k] * shi1_6[k];

        t_35[k] = f_3 * pc_z[k] * sih_24[k];

        t_36[k] = f_11 * shh_5[k]
                  + f_3 * pc_y[k] * sih_26[k];

        t_37[k] = pb_y[k] * shi0_9[k]
                  - f_10 * pc_y[k] * shi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, shi0_10, shi0_12, shh_6, \
                         shh_8, shh_9, shi1_10, shi1_12, sih_27, \
                         sih_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * shi0_10[k]
                  + f_14 * shh_6[k]
                  - f_10 * pc_y[k] * shi1_10[k];

        t_39[k] = f_3 * pc_z[k] * sih_27[k];

        t_40[k] = pb_y[k] * shi0_12[k]
                  + f_12 * shh_8[k]
                  - f_10 * pc_y[k] * shi1_12[k];

        t_41[k] = f_11 * shh_9[k]
                  + f_3 * pc_y[k] * sih_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_y, pc_x, pc_y, shi0_14, shh_36, shh_37, \
                         shh_38, shi1_14, sih_36, sih_37, sih_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * shi0_14[k]
                  - f_10 * pc_y[k] * shi1_14[k];

        t_43[k] = f_15 * shh_36[k]
                  + f_3 * pc_x[k] * sih_36[k];

        t_44[k] = f_15 * shh_37[k]
                  + f_3 * pc_x[k] * sih_37[k];

        t_45[k] = f_15 * shh_38[k]
                  + f_3 * pc_x[k] * sih_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, shh_15, shh_39, shh_40, shh_41, \
                         sig0_25, sig1_25, sih_36, sih_39, sih_40, \
                         sih_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_15 * shh_39[k]
                  + f_3 * pc_x[k] * sih_39[k];

        t_47[k] = f_15 * shh_40[k]
                  + f_3 * pc_x[k] * sih_40[k];

        t_48[k] = f_15 * shh_41[k]
                  + f_3 * pc_x[k] * sih_41[k];

        t_49[k] = f_11 * shh_15[k]
                  + f_1 * sig0_25[k]
                  - f_2 * sig1_25[k]
                  + f_3 * pc_y[k] * sih_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pc_y, pc_z, shh_17, shh_18, sig0_27, sig0_28, \
                         sig1_27, sig1_28, sih_36, sih_38, sih_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * sih_36[k];

        t_51[k] = f_11 * shh_17[k]
                  + f_4 * sig0_27[k]
                  - f_5 * sig1_27[k]
                  + f_3 * pc_y[k] * sih_38[k];

        t_52[k] = f_11 * shh_18[k]
                  + f_6 * sig0_28[k]
                  - f_7 * sig1_28[k]
                  + f_3 * pc_y[k] * sih_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_y, shi0_27, shh_19, shh_20, shi1_27, \
                         sig0_29, sig1_29, sih_40, sih_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_11 * shh_19[k]
                  + f_8 * sig0_29[k]
                  - f_9 * sig1_29[k]
                  + f_3 * pc_y[k] * sih_40[k];

        t_54[k] = f_11 * shh_20[k]
                  + f_3 * pc_y[k] * sih_41[k];

        t_55[k] = pb_y[k] * shi0_27[k]
                  - f_10 * pc_y[k] * shi1_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pb_z, pc_y, pc_z, shi0_0, shi0_3, \
                         shh_0, shi1_0, shi1_3, sih_42, sih_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * shi0_0[k]
                  - f_10 * pc_z[k] * shi1_0[k];

        t_57[k] = f_3 * pc_y[k] * sih_42[k];

        t_58[k] = f_11 * shh_0[k]
                  + f_3 * pc_z[k] * sih_42[k];

        t_59[k] = pb_z[k] * shi0_3[k]
                  - f_10 * pc_z[k] * shi1_3[k];

        t_60[k] = f_3 * pc_y[k] * sih_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, shi0_5, shi0_6, shh_2, \
                         shh_3, shi1_5, shi1_6, sih_45, sih_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_z[k] * shi0_5[k]
                  + f_12 * shh_2[k]
                  - f_10 * pc_z[k] * shi1_5[k];

        t_62[k] = pb_z[k] * shi0_6[k]
                  - f_10 * pc_z[k] * shi1_6[k];

        t_63[k] = f_11 * shh_3[k]
                  + f_3 * pc_z[k] * sih_45[k];

        t_64[k] = f_3 * pc_y[k] * sih_47[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_z, pc_z, shi0_9, shi0_10, shi0_12, shh_5, \
                         shh_6, shh_7, shi1_9, shi1_10, shi1_12, \
                         sih_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_z[k] * shi0_9[k]
                  + f_13 * shh_5[k]
                  - f_10 * pc_z[k] * shi1_9[k];

        t_66[k] = pb_z[k] * shi0_10[k]
                  - f_10 * pc_z[k] * shi1_10[k];

        t_67[k] = f_11 * shh_6[k]
                  + f_3 * pc_z[k] * sih_48[k];

        t_68[k] = pb_z[k] * shi0_12[k]
                  + f_12 * shh_7[k]
                  - f_10 * pc_z[k] * shi1_12[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_z, pc_x, pc_y, pc_z, shi0_14, shh_9, \
                         shh_57, shh_58, shi1_14, sih_51, sih_57, \
                         sih_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * pc_y[k] * sih_51[k];

        t_70[k] = pb_z[k] * shi0_14[k]
                  + f_14 * shh_9[k]
                  - f_10 * pc_z[k] * shi1_14[k];

        t_71[k] = f_15 * shh_57[k]
                  + f_3 * pc_x[k] * sih_57[k];

        t_72[k] = f_15 * shh_58[k]
                  + f_3 * pc_x[k] * sih_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pc_x, shh_59, shh_60, shh_61, shh_62, sih_59, \
                         sih_60, sih_61, sih_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_15 * shh_59[k]
                  + f_3 * pc_x[k] * sih_59[k];

        t_74[k] = f_15 * shh_60[k]
                  + f_3 * pc_x[k] * sih_60[k];

        t_75[k] = f_15 * shh_61[k]
                  + f_3 * pc_x[k] * sih_61[k];

        t_76[k] = f_15 * shh_62[k]
                  + f_3 * pc_x[k] * sih_62[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_z, pc_y, pc_z, shi0_21, shh_15, shi1_21, \
                         sig0_42, sig1_42, sih_57, sih_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_z[k] * shi0_21[k]
                  - f_10 * pc_z[k] * shi1_21[k];

        t_78[k] = f_11 * shh_15[k]
                  + f_3 * pc_z[k] * sih_57[k];

        t_79[k] = f_4 * sig0_42[k]
                  - f_5 * sig1_42[k]
                  + f_3 * pc_y[k] * sih_59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_y, pc_z, shh_20, sig0_43, sig0_44, \
                         sig1_43, sig1_44, sih_60, sih_61, sih_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_6 * sig0_43[k]
                  - f_7 * sig1_43[k]
                  + f_3 * pc_y[k] * sih_60[k];

        t_81[k] = f_8 * sig0_44[k]
                  - f_9 * sig1_44[k]
                  + f_3 * pc_y[k] * sih_61[k];

        t_82[k] = f_3 * pc_y[k] * sih_62[k];

        t_83[k] = f_11 * shh_20[k]
                  + f_1 * sig0_44[k]
                  - f_2 * sig1_44[k]
                  + f_3 * pc_z[k] * sih_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pc_x, pc_y, pc_z, shh_21, shh_63, shh_66, \
                         sig0_45, sig0_48, sig1_45, sig1_48, sih_63, \
                         sih_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_14 * shh_63[k]
                  + f_1 * sig0_45[k]
                  - f_2 * sig1_45[k]
                  + f_3 * pc_x[k] * sih_63[k];

        t_85[k] = f_12 * shh_21[k]
                  + f_3 * pc_y[k] * sih_63[k];

        t_86[k] = f_3 * pc_z[k] * sih_63[k];

        t_87[k] = f_14 * shh_66[k]
                  + f_4 * sig0_48[k]
                  - f_5 * sig1_48[k]
                  + f_3 * pc_x[k] * sih_66[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pc_x, pc_y, shh_23, shh_68, shh_69, sig0_50, \
                         sig0_51, sig1_50, sig1_51, sih_65, sih_68, \
                         sih_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_12 * shh_23[k]
                  + f_3 * pc_y[k] * sih_65[k];

        t_89[k] = f_14 * shh_68[k]
                  + f_4 * sig0_50[k]
                  - f_5 * sig1_50[k]
                  + f_3 * pc_x[k] * sih_68[k];

        t_90[k] = f_14 * shh_69[k]
                  + f_6 * sig0_51[k]
                  - f_7 * sig1_51[k]
                  + f_3 * pc_x[k] * sih_69[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pc_x, pc_y, pc_z, shh_26, shh_72, sig0_54, sig1_54, \
                         sih_66, sih_68, sih_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * pc_z[k] * sih_66[k];

        t_92[k] = f_12 * shh_26[k]
                  + f_3 * pc_y[k] * sih_68[k];

        t_93[k] = f_14 * shh_72[k]
                  + f_6 * sig0_54[k]
                  - f_7 * sig1_54[k]
                  + f_3 * pc_x[k] * sih_72[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pc_x, pc_z, shh_73, shh_75, sig0_55, sig0_57, \
                         sig1_55, sig1_57, sih_69, sih_73, sih_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_14 * shh_73[k]
                  + f_8 * sig0_55[k]
                  - f_9 * sig1_55[k]
                  + f_3 * pc_x[k] * sih_73[k];

        t_95[k] = f_3 * pc_z[k] * sih_69[k];

        t_96[k] = f_14 * shh_75[k]
                  + f_8 * sig0_57[k]
                  - f_9 * sig1_57[k]
                  + f_3 * pc_x[k] * sih_75[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, shh_30, shh_77, shh_78, shh_79, \
                         sig0_59, sig1_59, sih_72, sih_77, sih_78, \
                         sih_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_12 * shh_30[k]
                  + f_3 * pc_y[k] * sih_72[k];

        t_98[k] = f_14 * shh_77[k]
                  + f_8 * sig0_59[k]
                  - f_9 * sig1_59[k]
                  + f_3 * pc_x[k] * sih_77[k];

        t_99[k] = f_14 * shh_78[k]
                  + f_3 * pc_x[k] * sih_78[k];

        t_100[k] = f_14 * shh_79[k]
                   + f_3 * pc_x[k] * sih_79[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_x, shh_80, shh_81, shh_82, shh_83, \
                         sih_80, sih_81, sih_82, sih_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_14 * shh_80[k]
                   + f_3 * pc_x[k] * sih_80[k];

        t_102[k] = f_14 * shh_81[k]
                   + f_3 * pc_x[k] * sih_81[k];

        t_103[k] = f_14 * shh_82[k]
                   + f_3 * pc_x[k] * sih_82[k];

        t_104[k] = f_14 * shh_83[k]
                   + f_3 * pc_x[k] * sih_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pc_y, pc_z, shh_36, shh_38, sig0_55, sig0_57, \
                         sig1_55, sig1_57, sih_78, sih_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * shh_36[k]
                   + f_1 * sig0_55[k]
                   - f_2 * sig1_55[k]
                   + f_3 * pc_y[k] * sih_78[k];

        t_106[k] = f_3 * pc_z[k] * sih_78[k];

        t_107[k] = f_12 * shh_38[k]
                   + f_4 * sig0_57[k]
                   - f_5 * sig1_57[k]
                   + f_3 * pc_y[k] * sih_80[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pc_y, pc_z, shh_39, shh_40, shh_41, \
                         sig0_58, sig0_59, sig1_58, sig1_59, sih_81, sih_82, \
                         sih_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_12 * shh_39[k]
                   + f_6 * sig0_58[k]
                   - f_7 * sig1_58[k]
                   + f_3 * pc_y[k] * sih_81[k];

        t_109[k] = f_12 * shh_40[k]
                   + f_8 * sig0_59[k]
                   - f_9 * sig1_59[k]
                   + f_3 * pc_y[k] * sih_82[k];

        t_110[k] = f_12 * shh_41[k]
                   + f_3 * pc_y[k] * sih_83[k];

        t_111[k] = f_1 * sig0_59[k]
                   - f_2 * sig1_59[k]
                   + f_3 * pc_z[k] * sih_83[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pb_y, pb_z, pc_y, pc_z, shi0_31, shi0_56, \
                         shh_21, shh_42, shi1_31, shi1_56, sih_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = pb_y[k] * shi0_56[k]
                   - f_10 * pc_y[k] * shi1_56[k];

        t_113[k] = f_11 * shh_42[k]
                   + f_3 * pc_y[k] * sih_84[k];

        t_114[k] = f_11 * shh_21[k]
                   + f_3 * pc_z[k] * sih_84[k];

        t_115[k] = pb_z[k] * shi0_31[k]
                   - f_10 * pc_z[k] * shi1_31[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pb_z, pc_y, pc_z, shi0_34, shi0_61, \
                         shh_24, shh_44, shi1_34, shi1_61, sih_86, \
                         sih_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_11 * shh_44[k]
                   + f_3 * pc_y[k] * sih_86[k];

        t_117[k] = pb_y[k] * shi0_61[k]
                   - f_10 * pc_y[k] * shi1_61[k];

        t_118[k] = pb_z[k] * shi0_34[k]
                   - f_10 * pc_z[k] * shi1_34[k];

        t_119[k] = f_11 * shh_24[k]
                   + f_3 * pc_z[k] * sih_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pb_z, pc_y, pc_z, shi0_38, shi0_65, \
                         shh_27, shh_47, shi1_38, shi1_65, sih_89, \
                         sih_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_11 * shh_47[k]
                   + f_3 * pc_y[k] * sih_89[k];

        t_121[k] = pb_y[k] * shi0_65[k]
                   - f_10 * pc_y[k] * shi1_65[k];

        t_122[k] = pb_z[k] * shi0_38[k]
                   - f_10 * pc_z[k] * shi1_38[k];

        t_123[k] = f_11 * shh_27[k]
                   + f_3 * pc_z[k] * sih_90[k];
    }
}

static auto
compute_prim_sii_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shi0,
                                                          const size_t shh, const size_t shi1,
                                                          const size_t sig0, const size_t sig1,
                                                          const size_t sih, const size_t ncols,
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

    const auto *shi0_49 = buffer.data(shi0 + 49);
    const auto *shi0_68 = buffer.data(shi0 + 68);
    const auto *shi0_70 = buffer.data(shi0 + 70);
    const auto *shi0_83 = buffer.data(shi0 + 83);
    const auto *shi0_84 = buffer.data(shi0 + 84);
    const auto *shi0_87 = buffer.data(shi0 + 87);
    const auto *shi0_90 = buffer.data(shi0 + 90);
    const auto *shi0_94 = buffer.data(shi0 + 94);
    const auto *shi0_96 = buffer.data(shi0 + 96);
    const auto *shi0_105 = buffer.data(shi0 + 105);
    const auto *shi0_140 = buffer.data(shi0 + 140);
    const auto *shi0_143 = buffer.data(shi0 + 143);
    const auto *shi0_145 = buffer.data(shi0 + 145);
    const auto *shi0_146 = buffer.data(shi0 + 146);
    const auto *shi0_149 = buffer.data(shi0 + 149);
    const auto *shi0_150 = buffer.data(shi0 + 150);
    const auto *shi0_152 = buffer.data(shi0 + 152);
    const auto *shi0_154 = buffer.data(shi0 + 154);

    const auto *shh_36 = buffer.data(shh + 36);
    const auto *shh_42 = buffer.data(shh + 42);
    const auto *shh_45 = buffer.data(shh + 45);
    const auto *shh_48 = buffer.data(shh + 48);
    const auto *shh_50 = buffer.data(shh + 50);
    const auto *shh_51 = buffer.data(shh + 51);
    const auto *shh_57 = buffer.data(shh + 57);
    const auto *shh_59 = buffer.data(shh + 59);
    const auto *shh_60 = buffer.data(shh + 60);
    const auto *shh_61 = buffer.data(shh + 61);
    const auto *shh_62 = buffer.data(shh + 62);
    const auto *shh_63 = buffer.data(shh + 63);
    const auto *shh_65 = buffer.data(shh + 65);
    const auto *shh_66 = buffer.data(shh + 66);
    const auto *shh_68 = buffer.data(shh + 68);
    const auto *shh_69 = buffer.data(shh + 69);
    const auto *shh_70 = buffer.data(shh + 70);
    const auto *shh_72 = buffer.data(shh + 72);
    const auto *shh_78 = buffer.data(shh + 78);
    const auto *shh_80 = buffer.data(shh + 80);
    const auto *shh_81 = buffer.data(shh + 81);
    const auto *shh_82 = buffer.data(shh + 82);
    const auto *shh_83 = buffer.data(shh + 83);
    const auto *shh_84 = buffer.data(shh + 84);
    const auto *shh_86 = buffer.data(shh + 86);
    const auto *shh_87 = buffer.data(shh + 87);
    const auto *shh_89 = buffer.data(shh + 89);
    const auto *shh_90 = buffer.data(shh + 90);
    const auto *shh_93 = buffer.data(shh + 93);
    const auto *shh_99 = buffer.data(shh + 99);
    const auto *shh_100 = buffer.data(shh + 100);
    const auto *shh_101 = buffer.data(shh + 101);
    const auto *shh_102 = buffer.data(shh + 102);
    const auto *shh_103 = buffer.data(shh + 103);
    const auto *shh_104 = buffer.data(shh + 104);
    const auto *shh_105 = buffer.data(shh + 105);
    const auto *shh_106 = buffer.data(shh + 106);
    const auto *shh_107 = buffer.data(shh + 107);
    const auto *shh_108 = buffer.data(shh + 108);
    const auto *shh_110 = buffer.data(shh + 110);
    const auto *shh_111 = buffer.data(shh + 111);
    const auto *shh_113 = buffer.data(shh + 113);
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
    const auto *shh_152 = buffer.data(shh + 152);
    const auto *shh_156 = buffer.data(shh + 156);
    const auto *shh_161 = buffer.data(shh + 161);
    const auto *shh_162 = buffer.data(shh + 162);
    const auto *shh_163 = buffer.data(shh + 163);
    const auto *shh_164 = buffer.data(shh + 164);
    const auto *shh_165 = buffer.data(shh + 165);
    const auto *shh_166 = buffer.data(shh + 166);
    const auto *shh_167 = buffer.data(shh + 167);
    const auto *shh_183 = buffer.data(shh + 183);

    const auto *shi1_49 = buffer.data(shi1 + 49);
    const auto *shi1_68 = buffer.data(shi1 + 68);
    const auto *shi1_70 = buffer.data(shi1 + 70);
    const auto *shi1_83 = buffer.data(shi1 + 83);
    const auto *shi1_84 = buffer.data(shi1 + 84);
    const auto *shi1_87 = buffer.data(shi1 + 87);
    const auto *shi1_90 = buffer.data(shi1 + 90);
    const auto *shi1_94 = buffer.data(shi1 + 94);
    const auto *shi1_96 = buffer.data(shi1 + 96);
    const auto *shi1_105 = buffer.data(shi1 + 105);
    const auto *shi1_140 = buffer.data(shi1 + 140);
    const auto *shi1_143 = buffer.data(shi1 + 143);
    const auto *shi1_145 = buffer.data(shi1 + 145);
    const auto *shi1_146 = buffer.data(shi1 + 146);
    const auto *shi1_149 = buffer.data(shi1 + 149);
    const auto *shi1_150 = buffer.data(shi1 + 150);
    const auto *shi1_152 = buffer.data(shi1 + 152);
    const auto *shi1_154 = buffer.data(shi1 + 154);

    const auto *sig0_72 = buffer.data(sig0 + 72);
    const auto *sig0_73 = buffer.data(sig0 + 73);
    const auto *sig0_74 = buffer.data(sig0 + 74);
    const auto *sig0_75 = buffer.data(sig0 + 75);
    const auto *sig0_78 = buffer.data(sig0 + 78);
    const auto *sig0_80 = buffer.data(sig0 + 80);
    const auto *sig0_81 = buffer.data(sig0 + 81);
    const auto *sig0_84 = buffer.data(sig0 + 84);
    const auto *sig0_85 = buffer.data(sig0 + 85);
    const auto *sig0_87 = buffer.data(sig0 + 87);
    const auto *sig0_88 = buffer.data(sig0 + 88);
    const auto *sig0_89 = buffer.data(sig0 + 89);
    const auto *sig0_90 = buffer.data(sig0 + 90);
    const auto *sig0_93 = buffer.data(sig0 + 93);
    const auto *sig0_95 = buffer.data(sig0 + 95);
    const auto *sig0_96 = buffer.data(sig0 + 96);
    const auto *sig0_99 = buffer.data(sig0 + 99);
    const auto *sig0_100 = buffer.data(sig0 + 100);
    const auto *sig0_102 = buffer.data(sig0 + 102);
    const auto *sig0_103 = buffer.data(sig0 + 103);
    const auto *sig0_104 = buffer.data(sig0 + 104);
    const auto *sig0_110 = buffer.data(sig0 + 110);
    const auto *sig0_114 = buffer.data(sig0 + 114);
    const auto *sig0_117 = buffer.data(sig0 + 117);
    const auto *sig0_118 = buffer.data(sig0 + 118);
    const auto *sig0_119 = buffer.data(sig0 + 119);

    const auto *sig1_72 = buffer.data(sig1 + 72);
    const auto *sig1_73 = buffer.data(sig1 + 73);
    const auto *sig1_74 = buffer.data(sig1 + 74);
    const auto *sig1_75 = buffer.data(sig1 + 75);
    const auto *sig1_78 = buffer.data(sig1 + 78);
    const auto *sig1_80 = buffer.data(sig1 + 80);
    const auto *sig1_81 = buffer.data(sig1 + 81);
    const auto *sig1_84 = buffer.data(sig1 + 84);
    const auto *sig1_85 = buffer.data(sig1 + 85);
    const auto *sig1_87 = buffer.data(sig1 + 87);
    const auto *sig1_88 = buffer.data(sig1 + 88);
    const auto *sig1_89 = buffer.data(sig1 + 89);
    const auto *sig1_90 = buffer.data(sig1 + 90);
    const auto *sig1_93 = buffer.data(sig1 + 93);
    const auto *sig1_95 = buffer.data(sig1 + 95);
    const auto *sig1_96 = buffer.data(sig1 + 96);
    const auto *sig1_99 = buffer.data(sig1 + 99);
    const auto *sig1_100 = buffer.data(sig1 + 100);
    const auto *sig1_102 = buffer.data(sig1 + 102);
    const auto *sig1_103 = buffer.data(sig1 + 103);
    const auto *sig1_104 = buffer.data(sig1 + 104);
    const auto *sig1_110 = buffer.data(sig1 + 110);
    const auto *sig1_114 = buffer.data(sig1 + 114);
    const auto *sig1_117 = buffer.data(sig1 + 117);
    const auto *sig1_118 = buffer.data(sig1 + 118);
    const auto *sig1_119 = buffer.data(sig1 + 119);

    const auto *sih_93 = buffer.data(sih + 93);
    const auto *sih_99 = buffer.data(sih + 99);
    const auto *sih_100 = buffer.data(sih + 100);
    const auto *sih_101 = buffer.data(sih + 101);
    const auto *sih_102 = buffer.data(sih + 102);
    const auto *sih_103 = buffer.data(sih + 103);
    const auto *sih_104 = buffer.data(sih + 104);
    const auto *sih_105 = buffer.data(sih + 105);
    const auto *sih_107 = buffer.data(sih + 107);
    const auto *sih_108 = buffer.data(sih + 108);
    const auto *sih_110 = buffer.data(sih + 110);
    const auto *sih_111 = buffer.data(sih + 111);
    const auto *sih_114 = buffer.data(sih + 114);
    const auto *sih_115 = buffer.data(sih + 115);
    const auto *sih_117 = buffer.data(sih + 117);
    const auto *sih_119 = buffer.data(sih + 119);
    const auto *sih_120 = buffer.data(sih + 120);
    const auto *sih_121 = buffer.data(sih + 121);
    const auto *sih_122 = buffer.data(sih + 122);
    const auto *sih_123 = buffer.data(sih + 123);
    const auto *sih_124 = buffer.data(sih + 124);
    const auto *sih_125 = buffer.data(sih + 125);
    const auto *sih_126 = buffer.data(sih + 126);
    const auto *sih_128 = buffer.data(sih + 128);
    const auto *sih_129 = buffer.data(sih + 129);
    const auto *sih_131 = buffer.data(sih + 131);
    const auto *sih_132 = buffer.data(sih + 132);
    const auto *sih_135 = buffer.data(sih + 135);
    const auto *sih_136 = buffer.data(sih + 136);
    const auto *sih_138 = buffer.data(sih + 138);
    const auto *sih_140 = buffer.data(sih + 140);
    const auto *sih_141 = buffer.data(sih + 141);
    const auto *sih_142 = buffer.data(sih + 142);
    const auto *sih_143 = buffer.data(sih + 143);
    const auto *sih_144 = buffer.data(sih + 144);
    const auto *sih_145 = buffer.data(sih + 145);
    const auto *sih_146 = buffer.data(sih + 146);
    const auto *sih_147 = buffer.data(sih + 147);
    const auto *sih_149 = buffer.data(sih + 149);
    const auto *sih_150 = buffer.data(sih + 150);
    const auto *sih_152 = buffer.data(sih + 152);
    const auto *sih_153 = buffer.data(sih + 153);
    const auto *sih_156 = buffer.data(sih + 156);
    const auto *sih_161 = buffer.data(sih + 161);
    const auto *sih_162 = buffer.data(sih + 162);
    const auto *sih_163 = buffer.data(sih + 163);
    const auto *sih_164 = buffer.data(sih + 164);
    const auto *sih_165 = buffer.data(sih + 165);
    const auto *sih_166 = buffer.data(sih + 166);
    const auto *sih_167 = buffer.data(sih + 167);
    const auto *sih_168 = buffer.data(sih + 168);
    const auto *sih_170 = buffer.data(sih + 170);
    const auto *sih_171 = buffer.data(sih + 171);
    const auto *sih_173 = buffer.data(sih + 173);
    const auto *sih_174 = buffer.data(sih + 174);
    const auto *sih_177 = buffer.data(sih + 177);
    const auto *sih_183 = buffer.data(sih + 183);

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_y, pc_x, pc_y, shi0_68, shi0_70, \
                         shh_50, shh_51, shh_99, shi1_68, shi1_70, sih_93, \
                         sih_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * shi0_68[k]
                   + f_12 * shh_50[k]
                   - f_10 * pc_y[k] * shi1_68[k];

        t_125[k] = f_11 * shh_51[k]
                   + f_3 * pc_y[k] * sih_93[k];

        t_126[k] = pb_y[k] * shi0_70[k]
                   - f_10 * pc_y[k] * shi1_70[k];

        t_127[k] = f_14 * shh_99[k]
                   + f_3 * pc_x[k] * sih_99[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, shh_100, shh_101, shh_102, \
                         shh_103, shh_104, sih_100, sih_101, sih_102, sih_103, \
                         sih_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_14 * shh_100[k]
                   + f_3 * pc_x[k] * sih_100[k];

        t_129[k] = f_14 * shh_101[k]
                   + f_3 * pc_x[k] * sih_101[k];

        t_130[k] = f_14 * shh_102[k]
                   + f_3 * pc_x[k] * sih_102[k];

        t_131[k] = f_14 * shh_103[k]
                   + f_3 * pc_x[k] * sih_103[k];

        t_132[k] = f_14 * shh_104[k]
                   + f_3 * pc_x[k] * sih_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_z, pc_y, pc_z, shi0_49, shh_36, shh_59, \
                         shi1_49, sig0_72, sig1_72, sih_99, sih_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pb_z[k] * shi0_49[k]
                   - f_10 * pc_z[k] * shi1_49[k];

        t_134[k] = f_11 * shh_36[k]
                   + f_3 * pc_z[k] * sih_99[k];

        t_135[k] = f_11 * shh_59[k]
                   + f_4 * sig0_72[k]
                   - f_5 * sig1_72[k]
                   + f_3 * pc_y[k] * sih_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, shh_60, shh_61, shh_62, sig0_73, sig0_74, \
                         sig1_73, sig1_74, sih_102, sih_103, sih_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * shh_60[k]
                   + f_6 * sig0_73[k]
                   - f_7 * sig1_73[k]
                   + f_3 * pc_y[k] * sih_102[k];

        t_137[k] = f_11 * shh_61[k]
                   + f_8 * sig0_74[k]
                   - f_9 * sig1_74[k]
                   + f_3 * pc_y[k] * sih_103[k];

        t_138[k] = f_11 * shh_62[k]
                   + f_3 * pc_y[k] * sih_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pb_y, pc_x, pc_y, pc_z, shi0_83, shh_42, \
                         shh_105, shi1_83, sig0_75, sig1_75, sih_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * shi0_83[k]
                   - f_10 * pc_y[k] * shi1_83[k];

        t_140[k] = f_14 * shh_105[k]
                   + f_1 * sig0_75[k]
                   - f_2 * sig1_75[k]
                   + f_3 * pc_x[k] * sih_105[k];

        t_141[k] = f_3 * pc_y[k] * sih_105[k];

        t_142[k] = f_12 * shh_42[k]
                   + f_3 * pc_z[k] * sih_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, shh_108, shh_110, sig0_78, sig0_80, \
                         sig1_78, sig1_80, sih_107, sih_108, sih_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_14 * shh_108[k]
                   + f_4 * sig0_78[k]
                   - f_5 * sig1_78[k]
                   + f_3 * pc_x[k] * sih_108[k];

        t_144[k] = f_3 * pc_y[k] * sih_107[k];

        t_145[k] = f_14 * shh_110[k]
                   + f_4 * sig0_80[k]
                   - f_5 * sig1_80[k]
                   + f_3 * pc_x[k] * sih_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_x, pc_y, pc_z, shh_45, shh_111, sig0_81, \
                         sig1_81, sih_108, sih_110, sih_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_14 * shh_111[k]
                   + f_6 * sig0_81[k]
                   - f_7 * sig1_81[k]
                   + f_3 * pc_x[k] * sih_111[k];

        t_147[k] = f_12 * shh_45[k]
                   + f_3 * pc_z[k] * sih_108[k];

        t_148[k] = f_3 * pc_y[k] * sih_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_z, shh_48, shh_114, shh_115, sig0_84, \
                         sig0_85, sig1_84, sig1_85, sih_111, sih_114, \
                         sih_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_14 * shh_114[k]
                   + f_6 * sig0_84[k]
                   - f_7 * sig1_84[k]
                   + f_3 * pc_x[k] * sih_114[k];

        t_150[k] = f_14 * shh_115[k]
                   + f_8 * sig0_85[k]
                   - f_9 * sig1_85[k]
                   + f_3 * pc_x[k] * sih_115[k];

        t_151[k] = f_12 * shh_48[k]
                   + f_3 * pc_z[k] * sih_111[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pc_x, pc_y, shh_117, shh_119, sig0_87, sig0_89, \
                         sig1_87, sig1_89, sih_114, sih_117, sih_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_14 * shh_117[k]
                   + f_8 * sig0_87[k]
                   - f_9 * sig1_87[k]
                   + f_3 * pc_x[k] * sih_117[k];

        t_153[k] = f_3 * pc_y[k] * sih_114[k];

        t_154[k] = f_14 * shh_119[k]
                   + f_8 * sig0_89[k]
                   - f_9 * sig1_89[k]
                   + f_3 * pc_x[k] * sih_119[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pc_x, shh_120, shh_121, shh_122, \
                         shh_123, shh_124, sih_120, sih_121, sih_122, sih_123, \
                         sih_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_14 * shh_120[k]
                   + f_3 * pc_x[k] * sih_120[k];

        t_156[k] = f_14 * shh_121[k]
                   + f_3 * pc_x[k] * sih_121[k];

        t_157[k] = f_14 * shh_122[k]
                   + f_3 * pc_x[k] * sih_122[k];

        t_158[k] = f_14 * shh_123[k]
                   + f_3 * pc_x[k] * sih_123[k];

        t_159[k] = f_14 * shh_124[k]
                   + f_3 * pc_x[k] * sih_124[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, shh_57, shh_125, \
                         sig0_85, sig0_87, sig1_85, sig1_87, sih_120, sih_122, \
                         sih_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_14 * shh_125[k]
                   + f_3 * pc_x[k] * sih_125[k];

        t_161[k] = f_1 * sig0_85[k]
                   - f_2 * sig1_85[k]
                   + f_3 * pc_y[k] * sih_120[k];

        t_162[k] = f_12 * shh_57[k]
                   + f_3 * pc_z[k] * sih_120[k];

        t_163[k] = f_4 * sig0_87[k]
                   - f_5 * sig1_87[k]
                   + f_3 * pc_y[k] * sih_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, shh_62, sig0_88, sig0_89, \
                         sig1_88, sig1_89, sih_123, sih_124, sih_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * sig0_88[k]
                   - f_7 * sig1_88[k]
                   + f_3 * pc_y[k] * sih_123[k];

        t_165[k] = f_8 * sig0_89[k]
                   - f_9 * sig1_89[k]
                   + f_3 * pc_y[k] * sih_124[k];

        t_166[k] = f_3 * pc_y[k] * sih_125[k];

        t_167[k] = f_12 * shh_62[k]
                   + f_1 * sig0_89[k]
                   - f_2 * sig1_89[k]
                   + f_3 * pc_z[k] * sih_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, shh_63, shh_126, \
                         shh_129, sig0_90, sig0_93, sig1_90, sig1_93, sih_126, \
                         sih_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_13 * shh_126[k]
                   + f_1 * sig0_90[k]
                   - f_2 * sig1_90[k]
                   + f_3 * pc_x[k] * sih_126[k];

        t_169[k] = f_13 * shh_63[k]
                   + f_3 * pc_y[k] * sih_126[k];

        t_170[k] = f_3 * pc_z[k] * sih_126[k];

        t_171[k] = f_13 * shh_129[k]
                   + f_4 * sig0_93[k]
                   - f_5 * sig1_93[k]
                   + f_3 * pc_x[k] * sih_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pc_x, pc_y, shh_65, shh_131, shh_132, sig0_95, \
                         sig0_96, sig1_95, sig1_96, sih_128, sih_131, \
                         sih_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_13 * shh_65[k]
                   + f_3 * pc_y[k] * sih_128[k];

        t_173[k] = f_13 * shh_131[k]
                   + f_4 * sig0_95[k]
                   - f_5 * sig1_95[k]
                   + f_3 * pc_x[k] * sih_131[k];

        t_174[k] = f_13 * shh_132[k]
                   + f_6 * sig0_96[k]
                   - f_7 * sig1_96[k]
                   + f_3 * pc_x[k] * sih_132[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pc_x, pc_y, pc_z, shh_68, shh_135, sig0_99, \
                         sig1_99, sih_129, sih_131, sih_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_3 * pc_z[k] * sih_129[k];

        t_176[k] = f_13 * shh_68[k]
                   + f_3 * pc_y[k] * sih_131[k];

        t_177[k] = f_13 * shh_135[k]
                   + f_6 * sig0_99[k]
                   - f_7 * sig1_99[k]
                   + f_3 * pc_x[k] * sih_135[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_z, shh_136, shh_138, sig0_100, \
                         sig0_102, sig1_100, sig1_102, sih_132, sih_136, \
                         sih_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_13 * shh_136[k]
                   + f_8 * sig0_100[k]
                   - f_9 * sig1_100[k]
                   + f_3 * pc_x[k] * sih_136[k];

        t_179[k] = f_3 * pc_z[k] * sih_132[k];

        t_180[k] = f_13 * shh_138[k]
                   + f_8 * sig0_102[k]
                   - f_9 * sig1_102[k]
                   + f_3 * pc_x[k] * sih_138[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, shh_72, shh_140, shh_141, \
                         shh_142, sig0_104, sig1_104, sih_135, sih_140, sih_141, \
                         sih_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_13 * shh_72[k]
                   + f_3 * pc_y[k] * sih_135[k];

        t_182[k] = f_13 * shh_140[k]
                   + f_8 * sig0_104[k]
                   - f_9 * sig1_104[k]
                   + f_3 * pc_x[k] * sih_140[k];

        t_183[k] = f_13 * shh_141[k]
                   + f_3 * pc_x[k] * sih_141[k];

        t_184[k] = f_13 * shh_142[k]
                   + f_3 * pc_x[k] * sih_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, shh_143, shh_144, shh_145, shh_146, \
                         sih_143, sih_144, sih_145, sih_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_13 * shh_143[k]
                   + f_3 * pc_x[k] * sih_143[k];

        t_186[k] = f_13 * shh_144[k]
                   + f_3 * pc_x[k] * sih_144[k];

        t_187[k] = f_13 * shh_145[k]
                   + f_3 * pc_x[k] * sih_145[k];

        t_188[k] = f_13 * shh_146[k]
                   + f_3 * pc_x[k] * sih_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_y, pc_z, shh_78, shh_80, sig0_100, sig0_102, \
                         sig1_100, sig1_102, sih_141, sih_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * shh_78[k]
                   + f_1 * sig0_100[k]
                   - f_2 * sig1_100[k]
                   + f_3 * pc_y[k] * sih_141[k];

        t_190[k] = f_3 * pc_z[k] * sih_141[k];

        t_191[k] = f_13 * shh_80[k]
                   + f_4 * sig0_102[k]
                   - f_5 * sig1_102[k]
                   + f_3 * pc_y[k] * sih_143[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_y, pc_z, shh_81, shh_82, shh_83, \
                         sig0_103, sig0_104, sig1_103, sig1_104, sih_144, sih_145, \
                         sih_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_13 * shh_81[k]
                   + f_6 * sig0_103[k]
                   - f_7 * sig1_103[k]
                   + f_3 * pc_y[k] * sih_144[k];

        t_193[k] = f_13 * shh_82[k]
                   + f_8 * sig0_104[k]
                   - f_9 * sig1_104[k]
                   + f_3 * pc_y[k] * sih_145[k];

        t_194[k] = f_13 * shh_83[k]
                   + f_3 * pc_y[k] * sih_146[k];

        t_195[k] = f_1 * sig0_104[k]
                   - f_2 * sig1_104[k]
                   + f_3 * pc_z[k] * sih_146[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pb_z, pc_y, pc_z, shi0_84, shi0_87, \
                         shh_63, shh_84, shi1_84, shi1_87, sih_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pb_z[k] * shi0_84[k]
                   - f_10 * pc_z[k] * shi1_84[k];

        t_197[k] = f_12 * shh_84[k]
                   + f_3 * pc_y[k] * sih_147[k];

        t_198[k] = f_11 * shh_63[k]
                   + f_3 * pc_z[k] * sih_147[k];

        t_199[k] = pb_z[k] * shi0_87[k]
                   - f_10 * pc_z[k] * shi1_87[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_z, pc_x, pc_y, pc_z, shi0_90, shh_86, \
                         shh_152, shi1_90, sig0_110, sig1_110, sih_149, \
                         sih_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_12 * shh_86[k]
                   + f_3 * pc_y[k] * sih_149[k];

        t_201[k] = f_13 * shh_152[k]
                   + f_4 * sig0_110[k]
                   - f_5 * sig1_110[k]
                   + f_3 * pc_x[k] * sih_152[k];

        t_202[k] = pb_z[k] * shi0_90[k]
                   - f_10 * pc_z[k] * shi1_90[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pc_x, pc_y, pc_z, shh_66, shh_89, shh_156, \
                         sig0_114, sig1_114, sih_150, sih_152, \
                         sih_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_11 * shh_66[k]
                   + f_3 * pc_z[k] * sih_150[k];

        t_204[k] = f_12 * shh_89[k]
                   + f_3 * pc_y[k] * sih_152[k];

        t_205[k] = f_13 * shh_156[k]
                   + f_6 * sig0_114[k]
                   - f_7 * sig1_114[k]
                   + f_3 * pc_x[k] * sih_156[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pb_z, pc_y, pc_z, shi0_94, shi0_96, \
                         shh_69, shh_70, shh_93, shi1_94, shi1_96, sih_153, \
                         sih_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pb_z[k] * shi0_94[k]
                   - f_10 * pc_z[k] * shi1_94[k];

        t_207[k] = f_11 * shh_69[k]
                   + f_3 * pc_z[k] * sih_153[k];

        t_208[k] = pb_z[k] * shi0_96[k]
                   + f_12 * shh_70[k]
                   - f_10 * pc_z[k] * shi1_96[k];

        t_209[k] = f_12 * shh_93[k]
                   + f_3 * pc_y[k] * sih_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, shh_161, shh_162, shh_163, shh_164, \
                         sig0_119, sig1_119, sih_161, sih_162, sih_163, \
                         sih_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_13 * shh_161[k]
                   + f_8 * sig0_119[k]
                   - f_9 * sig1_119[k]
                   + f_3 * pc_x[k] * sih_161[k];

        t_211[k] = f_13 * shh_162[k]
                   + f_3 * pc_x[k] * sih_162[k];

        t_212[k] = f_13 * shh_163[k]
                   + f_3 * pc_x[k] * sih_163[k];

        t_213[k] = f_13 * shh_164[k]
                   + f_3 * pc_x[k] * sih_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pb_z, pc_x, pc_z, shi0_105, shh_165, \
                         shh_166, shh_167, shi1_105, sih_165, sih_166, \
                         sih_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_13 * shh_165[k]
                   + f_3 * pc_x[k] * sih_165[k];

        t_215[k] = f_13 * shh_166[k]
                   + f_3 * pc_x[k] * sih_166[k];

        t_216[k] = f_13 * shh_167[k]
                   + f_3 * pc_x[k] * sih_167[k];

        t_217[k] = pb_z[k] * shi0_105[k]
                   - f_10 * pc_z[k] * shi1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, shh_78, shh_101, shh_102, sig0_117, \
                         sig0_118, sig1_117, sig1_118, sih_162, sih_164, \
                         sih_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * shh_78[k]
                   + f_3 * pc_z[k] * sih_162[k];

        t_219[k] = f_12 * shh_101[k]
                   + f_4 * sig0_117[k]
                   - f_5 * sig1_117[k]
                   + f_3 * pc_y[k] * sih_164[k];

        t_220[k] = f_12 * shh_102[k]
                   + f_6 * sig0_118[k]
                   - f_7 * sig1_118[k]
                   + f_3 * pc_y[k] * sih_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pb_y, pc_y, pc_z, shi0_140, shh_83, \
                         shh_103, shh_104, shi1_140, sig0_119, sig1_119, sih_166, \
                         sih_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * shh_103[k]
                   + f_8 * sig0_119[k]
                   - f_9 * sig1_119[k]
                   + f_3 * pc_y[k] * sih_166[k];

        t_222[k] = f_12 * shh_104[k]
                   + f_3 * pc_y[k] * sih_167[k];

        t_223[k] = f_11 * shh_83[k]
                   + f_1 * sig0_119[k]
                   - f_2 * sig1_119[k]
                   + f_3 * pc_z[k] * sih_167[k];

        t_224[k] = pb_y[k] * shi0_140[k]
                   - f_10 * pc_y[k] * shi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_y, pc_y, pc_z, shi0_143, shh_84, \
                         shh_105, shh_106, shh_107, shi1_143, sih_168, \
                         sih_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * shh_105[k]
                   + f_3 * pc_y[k] * sih_168[k];

        t_226[k] = f_12 * shh_84[k]
                   + f_3 * pc_z[k] * sih_168[k];

        t_227[k] = pb_y[k] * shi0_143[k]
                   + f_12 * shh_106[k]
                   - f_10 * pc_y[k] * shi1_143[k];

        t_228[k] = f_11 * shh_107[k]
                   + f_3 * pc_y[k] * sih_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_y, pc_y, pc_z, shi0_145, shi0_146, \
                         shh_87, shh_108, shh_110, shi1_145, shi1_146, sih_171, \
                         sih_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * shi0_145[k]
                   - f_10 * pc_y[k] * shi1_145[k];

        t_230[k] = pb_y[k] * shi0_146[k]
                   + f_13 * shh_108[k]
                   - f_10 * pc_y[k] * shi1_146[k];

        t_231[k] = f_12 * shh_87[k]
                   + f_3 * pc_z[k] * sih_171[k];

        t_232[k] = f_11 * shh_110[k]
                   + f_3 * pc_y[k] * sih_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_y, pc_y, pc_z, shi0_149, shi0_150, shh_90, \
                         shh_111, shi1_149, shi1_150, sih_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * shi0_149[k]
                   - f_10 * pc_y[k] * shi1_149[k];

        t_234[k] = pb_y[k] * shi0_150[k]
                   + f_14 * shh_111[k]
                   - f_10 * pc_y[k] * shi1_150[k];

        t_235[k] = f_12 * shh_90[k]
                   + f_3 * pc_z[k] * sih_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pb_y, pc_x, pc_y, shi0_152, shi0_154, \
                         shh_113, shh_114, shh_183, shi1_152, shi1_154, sih_177, \
                         sih_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pb_y[k] * shi0_152[k]
                   + f_12 * shh_113[k]
                   - f_10 * pc_y[k] * shi1_152[k];

        t_237[k] = f_11 * shh_114[k]
                   + f_3 * pc_y[k] * sih_177[k];

        t_238[k] = pb_y[k] * shi0_154[k]
                   - f_10 * pc_y[k] * shi1_154[k];

        t_239[k] = f_13 * shh_183[k]
                   + f_3 * pc_x[k] * sih_183[k];
    }
}

static auto
compute_prim_sii_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shi0,
                                                          const size_t shh, const size_t shi1,
                                                          const size_t sig0, const size_t sig1,
                                                          const size_t sih, const size_t ncols,
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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shi0_167 = buffer.data(shi0 + 167);
    const auto *shi0_168 = buffer.data(shi0 + 168);
    const auto *shi0_171 = buffer.data(shi0 + 171);
    const auto *shi0_174 = buffer.data(shi0 + 174);
    const auto *shi0_178 = buffer.data(shi0 + 178);
    const auto *shi0_180 = buffer.data(shi0 + 180);
    const auto *shi0_189 = buffer.data(shi0 + 189);

    const auto *shh_99 = buffer.data(shh + 99);
    const auto *shh_105 = buffer.data(shh + 105);
    const auto *shh_108 = buffer.data(shh + 108);
    const auto *shh_111 = buffer.data(shh + 111);
    const auto *shh_120 = buffer.data(shh + 120);
    const auto *shh_122 = buffer.data(shh + 122);
    const auto *shh_123 = buffer.data(shh + 123);
    const auto *shh_124 = buffer.data(shh + 124);
    const auto *shh_125 = buffer.data(shh + 125);
    const auto *shh_126 = buffer.data(shh + 126);
    const auto *shh_128 = buffer.data(shh + 128);
    const auto *shh_129 = buffer.data(shh + 129);
    const auto *shh_131 = buffer.data(shh + 131);
    const auto *shh_132 = buffer.data(shh + 132);
    const auto *shh_133 = buffer.data(shh + 133);
    const auto *shh_135 = buffer.data(shh + 135);
    const auto *shh_141 = buffer.data(shh + 141);
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
    const auto *shh_164 = buffer.data(shh + 164);
    const auto *shh_165 = buffer.data(shh + 165);
    const auto *shh_166 = buffer.data(shh + 166);
    const auto *shh_167 = buffer.data(shh + 167);
    const auto *shh_168 = buffer.data(shh + 168);
    const auto *shh_170 = buffer.data(shh + 170);
    const auto *shh_173 = buffer.data(shh + 173);
    const auto *shh_177 = buffer.data(shh + 177);
    const auto *shh_184 = buffer.data(shh + 184);
    const auto *shh_185 = buffer.data(shh + 185);
    const auto *shh_186 = buffer.data(shh + 186);
    const auto *shh_187 = buffer.data(shh + 187);
    const auto *shh_188 = buffer.data(shh + 188);
    const auto *shh_189 = buffer.data(shh + 189);
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
    const auto *shh_213 = buffer.data(shh + 213);
    const auto *shh_215 = buffer.data(shh + 215);
    const auto *shh_216 = buffer.data(shh + 216);
    const auto *shh_219 = buffer.data(shh + 219);
    const auto *shh_220 = buffer.data(shh + 220);
    const auto *shh_222 = buffer.data(shh + 222);
    const auto *shh_224 = buffer.data(shh + 224);
    const auto *shh_225 = buffer.data(shh + 225);
    const auto *shh_226 = buffer.data(shh + 226);
    const auto *shh_227 = buffer.data(shh + 227);
    const auto *shh_228 = buffer.data(shh + 228);
    const auto *shh_229 = buffer.data(shh + 229);
    const auto *shh_230 = buffer.data(shh + 230);
    const auto *shh_236 = buffer.data(shh + 236);
    const auto *shh_240 = buffer.data(shh + 240);
    const auto *shh_245 = buffer.data(shh + 245);
    const auto *shh_246 = buffer.data(shh + 246);
    const auto *shh_247 = buffer.data(shh + 247);
    const auto *shh_248 = buffer.data(shh + 248);
    const auto *shh_249 = buffer.data(shh + 249);
    const auto *shh_250 = buffer.data(shh + 250);
    const auto *shh_251 = buffer.data(shh + 251);
    const auto *shh_252 = buffer.data(shh + 252);
    const auto *shh_255 = buffer.data(shh + 255);
    const auto *shh_257 = buffer.data(shh + 257);
    const auto *shh_258 = buffer.data(shh + 258);
    const auto *shh_261 = buffer.data(shh + 261);
    const auto *shh_262 = buffer.data(shh + 262);
    const auto *shh_264 = buffer.data(shh + 264);
    const auto *shh_266 = buffer.data(shh + 266);

    const auto *shi1_167 = buffer.data(shi1 + 167);
    const auto *shi1_168 = buffer.data(shi1 + 168);
    const auto *shi1_171 = buffer.data(shi1 + 171);
    const auto *shi1_174 = buffer.data(shi1 + 174);
    const auto *shi1_178 = buffer.data(shi1 + 178);
    const auto *shi1_180 = buffer.data(shi1 + 180);
    const auto *shi1_189 = buffer.data(shi1 + 189);

    const auto *sig0_130 = buffer.data(sig0 + 130);
    const auto *sig0_132 = buffer.data(sig0 + 132);
    const auto *sig0_133 = buffer.data(sig0 + 133);
    const auto *sig0_134 = buffer.data(sig0 + 134);
    const auto *sig0_135 = buffer.data(sig0 + 135);
    const auto *sig0_138 = buffer.data(sig0 + 138);
    const auto *sig0_140 = buffer.data(sig0 + 140);
    const auto *sig0_141 = buffer.data(sig0 + 141);
    const auto *sig0_144 = buffer.data(sig0 + 144);
    const auto *sig0_145 = buffer.data(sig0 + 145);
    const auto *sig0_147 = buffer.data(sig0 + 147);
    const auto *sig0_148 = buffer.data(sig0 + 148);
    const auto *sig0_149 = buffer.data(sig0 + 149);
    const auto *sig0_150 = buffer.data(sig0 + 150);
    const auto *sig0_153 = buffer.data(sig0 + 153);
    const auto *sig0_155 = buffer.data(sig0 + 155);
    const auto *sig0_156 = buffer.data(sig0 + 156);
    const auto *sig0_159 = buffer.data(sig0 + 159);
    const auto *sig0_160 = buffer.data(sig0 + 160);
    const auto *sig0_162 = buffer.data(sig0 + 162);
    const auto *sig0_163 = buffer.data(sig0 + 163);
    const auto *sig0_164 = buffer.data(sig0 + 164);
    const auto *sig0_170 = buffer.data(sig0 + 170);
    const auto *sig0_174 = buffer.data(sig0 + 174);
    const auto *sig0_177 = buffer.data(sig0 + 177);
    const auto *sig0_178 = buffer.data(sig0 + 178);
    const auto *sig0_179 = buffer.data(sig0 + 179);
    const auto *sig0_180 = buffer.data(sig0 + 180);
    const auto *sig0_183 = buffer.data(sig0 + 183);
    const auto *sig0_185 = buffer.data(sig0 + 185);
    const auto *sig0_186 = buffer.data(sig0 + 186);
    const auto *sig0_189 = buffer.data(sig0 + 189);
    const auto *sig0_190 = buffer.data(sig0 + 190);
    const auto *sig0_192 = buffer.data(sig0 + 192);
    const auto *sig0_194 = buffer.data(sig0 + 194);

    const auto *sig1_130 = buffer.data(sig1 + 130);
    const auto *sig1_132 = buffer.data(sig1 + 132);
    const auto *sig1_133 = buffer.data(sig1 + 133);
    const auto *sig1_134 = buffer.data(sig1 + 134);
    const auto *sig1_135 = buffer.data(sig1 + 135);
    const auto *sig1_138 = buffer.data(sig1 + 138);
    const auto *sig1_140 = buffer.data(sig1 + 140);
    const auto *sig1_141 = buffer.data(sig1 + 141);
    const auto *sig1_144 = buffer.data(sig1 + 144);
    const auto *sig1_145 = buffer.data(sig1 + 145);
    const auto *sig1_147 = buffer.data(sig1 + 147);
    const auto *sig1_148 = buffer.data(sig1 + 148);
    const auto *sig1_149 = buffer.data(sig1 + 149);
    const auto *sig1_150 = buffer.data(sig1 + 150);
    const auto *sig1_153 = buffer.data(sig1 + 153);
    const auto *sig1_155 = buffer.data(sig1 + 155);
    const auto *sig1_156 = buffer.data(sig1 + 156);
    const auto *sig1_159 = buffer.data(sig1 + 159);
    const auto *sig1_160 = buffer.data(sig1 + 160);
    const auto *sig1_162 = buffer.data(sig1 + 162);
    const auto *sig1_163 = buffer.data(sig1 + 163);
    const auto *sig1_164 = buffer.data(sig1 + 164);
    const auto *sig1_170 = buffer.data(sig1 + 170);
    const auto *sig1_174 = buffer.data(sig1 + 174);
    const auto *sig1_177 = buffer.data(sig1 + 177);
    const auto *sig1_178 = buffer.data(sig1 + 178);
    const auto *sig1_179 = buffer.data(sig1 + 179);
    const auto *sig1_180 = buffer.data(sig1 + 180);
    const auto *sig1_183 = buffer.data(sig1 + 183);
    const auto *sig1_185 = buffer.data(sig1 + 185);
    const auto *sig1_186 = buffer.data(sig1 + 186);
    const auto *sig1_189 = buffer.data(sig1 + 189);
    const auto *sig1_190 = buffer.data(sig1 + 190);
    const auto *sig1_192 = buffer.data(sig1 + 192);
    const auto *sig1_194 = buffer.data(sig1 + 194);

    const auto *sih_183 = buffer.data(sih + 183);
    const auto *sih_184 = buffer.data(sih + 184);
    const auto *sih_185 = buffer.data(sih + 185);
    const auto *sih_186 = buffer.data(sih + 186);
    const auto *sih_187 = buffer.data(sih + 187);
    const auto *sih_188 = buffer.data(sih + 188);
    const auto *sih_189 = buffer.data(sih + 189);
    const auto *sih_191 = buffer.data(sih + 191);
    const auto *sih_192 = buffer.data(sih + 192);
    const auto *sih_194 = buffer.data(sih + 194);
    const auto *sih_195 = buffer.data(sih + 195);
    const auto *sih_198 = buffer.data(sih + 198);
    const auto *sih_199 = buffer.data(sih + 199);
    const auto *sih_201 = buffer.data(sih + 201);
    const auto *sih_203 = buffer.data(sih + 203);
    const auto *sih_204 = buffer.data(sih + 204);
    const auto *sih_205 = buffer.data(sih + 205);
    const auto *sih_206 = buffer.data(sih + 206);
    const auto *sih_207 = buffer.data(sih + 207);
    const auto *sih_208 = buffer.data(sih + 208);
    const auto *sih_209 = buffer.data(sih + 209);
    const auto *sih_210 = buffer.data(sih + 210);
    const auto *sih_212 = buffer.data(sih + 212);
    const auto *sih_213 = buffer.data(sih + 213);
    const auto *sih_215 = buffer.data(sih + 215);
    const auto *sih_216 = buffer.data(sih + 216);
    const auto *sih_219 = buffer.data(sih + 219);
    const auto *sih_220 = buffer.data(sih + 220);
    const auto *sih_222 = buffer.data(sih + 222);
    const auto *sih_224 = buffer.data(sih + 224);
    const auto *sih_225 = buffer.data(sih + 225);
    const auto *sih_226 = buffer.data(sih + 226);
    const auto *sih_227 = buffer.data(sih + 227);
    const auto *sih_228 = buffer.data(sih + 228);
    const auto *sih_229 = buffer.data(sih + 229);
    const auto *sih_230 = buffer.data(sih + 230);
    const auto *sih_231 = buffer.data(sih + 231);
    const auto *sih_233 = buffer.data(sih + 233);
    const auto *sih_234 = buffer.data(sih + 234);
    const auto *sih_236 = buffer.data(sih + 236);
    const auto *sih_237 = buffer.data(sih + 237);
    const auto *sih_240 = buffer.data(sih + 240);
    const auto *sih_245 = buffer.data(sih + 245);
    const auto *sih_246 = buffer.data(sih + 246);
    const auto *sih_247 = buffer.data(sih + 247);
    const auto *sih_248 = buffer.data(sih + 248);
    const auto *sih_249 = buffer.data(sih + 249);
    const auto *sih_250 = buffer.data(sih + 250);
    const auto *sih_251 = buffer.data(sih + 251);
    const auto *sih_252 = buffer.data(sih + 252);
    const auto *sih_254 = buffer.data(sih + 254);
    const auto *sih_255 = buffer.data(sih + 255);
    const auto *sih_257 = buffer.data(sih + 257);
    const auto *sih_258 = buffer.data(sih + 258);
    const auto *sih_261 = buffer.data(sih + 261);
    const auto *sih_262 = buffer.data(sih + 262);
    const auto *sih_264 = buffer.data(sih + 264);
    const auto *sih_266 = buffer.data(sih + 266);

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, shh_184, shh_185, shh_186, \
                         shh_187, shh_188, sih_184, sih_185, sih_186, sih_187, \
                         sih_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_13 * shh_184[k]
                   + f_3 * pc_x[k] * sih_184[k];

        t_241[k] = f_13 * shh_185[k]
                   + f_3 * pc_x[k] * sih_185[k];

        t_242[k] = f_13 * shh_186[k]
                   + f_3 * pc_x[k] * sih_186[k];

        t_243[k] = f_13 * shh_187[k]
                   + f_3 * pc_x[k] * sih_187[k];

        t_244[k] = f_13 * shh_188[k]
                   + f_3 * pc_x[k] * sih_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, shh_99, shh_120, shh_122, sig0_130, \
                         sig0_132, sig1_130, sig1_132, sih_183, \
                         sih_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * shh_120[k]
                   + f_1 * sig0_130[k]
                   - f_2 * sig1_130[k]
                   + f_3 * pc_y[k] * sih_183[k];

        t_246[k] = f_12 * shh_99[k]
                   + f_3 * pc_z[k] * sih_183[k];

        t_247[k] = f_11 * shh_122[k]
                   + f_4 * sig0_132[k]
                   - f_5 * sig1_132[k]
                   + f_3 * pc_y[k] * sih_185[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, shh_123, shh_124, shh_125, sig0_133, \
                         sig0_134, sig1_133, sig1_134, sih_186, sih_187, \
                         sih_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * shh_123[k]
                   + f_6 * sig0_133[k]
                   - f_7 * sig1_133[k]
                   + f_3 * pc_y[k] * sih_186[k];

        t_249[k] = f_11 * shh_124[k]
                   + f_8 * sig0_134[k]
                   - f_9 * sig1_134[k]
                   + f_3 * pc_y[k] * sih_187[k];

        t_250[k] = f_11 * shh_125[k]
                   + f_3 * pc_y[k] * sih_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_y, pc_x, pc_y, pc_z, shi0_167, \
                         shh_105, shh_189, shi1_167, sig0_135, sig1_135, \
                         sih_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pb_y[k] * shi0_167[k]
                   - f_10 * pc_y[k] * shi1_167[k];

        t_252[k] = f_13 * shh_189[k]
                   + f_1 * sig0_135[k]
                   - f_2 * sig1_135[k]
                   + f_3 * pc_x[k] * sih_189[k];

        t_253[k] = f_3 * pc_y[k] * sih_189[k];

        t_254[k] = f_13 * shh_105[k]
                   + f_3 * pc_z[k] * sih_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, shh_192, shh_194, sig0_138, \
                         sig0_140, sig1_138, sig1_140, sih_191, sih_192, \
                         sih_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_13 * shh_192[k]
                   + f_4 * sig0_138[k]
                   - f_5 * sig1_138[k]
                   + f_3 * pc_x[k] * sih_192[k];

        t_256[k] = f_3 * pc_y[k] * sih_191[k];

        t_257[k] = f_13 * shh_194[k]
                   + f_4 * sig0_140[k]
                   - f_5 * sig1_140[k]
                   + f_3 * pc_x[k] * sih_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_x, pc_y, pc_z, shh_108, shh_195, sig0_141, \
                         sig1_141, sih_192, sih_194, sih_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_13 * shh_195[k]
                   + f_6 * sig0_141[k]
                   - f_7 * sig1_141[k]
                   + f_3 * pc_x[k] * sih_195[k];

        t_259[k] = f_13 * shh_108[k]
                   + f_3 * pc_z[k] * sih_192[k];

        t_260[k] = f_3 * pc_y[k] * sih_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_z, shh_111, shh_198, shh_199, sig0_144, \
                         sig0_145, sig1_144, sig1_145, sih_195, sih_198, \
                         sih_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_13 * shh_198[k]
                   + f_6 * sig0_144[k]
                   - f_7 * sig1_144[k]
                   + f_3 * pc_x[k] * sih_198[k];

        t_262[k] = f_13 * shh_199[k]
                   + f_8 * sig0_145[k]
                   - f_9 * sig1_145[k]
                   + f_3 * pc_x[k] * sih_199[k];

        t_263[k] = f_13 * shh_111[k]
                   + f_3 * pc_z[k] * sih_195[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, shh_201, shh_203, sig0_147, \
                         sig0_149, sig1_147, sig1_149, sih_198, sih_201, \
                         sih_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_13 * shh_201[k]
                   + f_8 * sig0_147[k]
                   - f_9 * sig1_147[k]
                   + f_3 * pc_x[k] * sih_201[k];

        t_265[k] = f_3 * pc_y[k] * sih_198[k];

        t_266[k] = f_13 * shh_203[k]
                   + f_8 * sig0_149[k]
                   - f_9 * sig1_149[k]
                   + f_3 * pc_x[k] * sih_203[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, pc_x, shh_204, shh_205, shh_206, \
                         shh_207, shh_208, sih_204, sih_205, sih_206, sih_207, \
                         sih_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_13 * shh_204[k]
                   + f_3 * pc_x[k] * sih_204[k];

        t_268[k] = f_13 * shh_205[k]
                   + f_3 * pc_x[k] * sih_205[k];

        t_269[k] = f_13 * shh_206[k]
                   + f_3 * pc_x[k] * sih_206[k];

        t_270[k] = f_13 * shh_207[k]
                   + f_3 * pc_x[k] * sih_207[k];

        t_271[k] = f_13 * shh_208[k]
                   + f_3 * pc_x[k] * sih_208[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, shh_120, shh_209, \
                         sig0_145, sig0_147, sig1_145, sig1_147, sih_204, sih_206, \
                         sih_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_13 * shh_209[k]
                   + f_3 * pc_x[k] * sih_209[k];

        t_273[k] = f_1 * sig0_145[k]
                   - f_2 * sig1_145[k]
                   + f_3 * pc_y[k] * sih_204[k];

        t_274[k] = f_13 * shh_120[k]
                   + f_3 * pc_z[k] * sih_204[k];

        t_275[k] = f_4 * sig0_147[k]
                   - f_5 * sig1_147[k]
                   + f_3 * pc_y[k] * sih_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, shh_125, sig0_148, sig0_149, \
                         sig1_148, sig1_149, sih_207, sih_208, \
                         sih_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * sig0_148[k]
                   - f_7 * sig1_148[k]
                   + f_3 * pc_y[k] * sih_207[k];

        t_277[k] = f_8 * sig0_149[k]
                   - f_9 * sig1_149[k]
                   + f_3 * pc_y[k] * sih_208[k];

        t_278[k] = f_3 * pc_y[k] * sih_209[k];

        t_279[k] = f_13 * shh_125[k]
                   + f_1 * sig0_149[k]
                   - f_2 * sig1_149[k]
                   + f_3 * pc_z[k] * sih_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, shh_126, shh_210, \
                         shh_213, sig0_150, sig0_153, sig1_150, sig1_153, sih_210, \
                         sih_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_12 * shh_210[k]
                   + f_1 * sig0_150[k]
                   - f_2 * sig1_150[k]
                   + f_3 * pc_x[k] * sih_210[k];

        t_281[k] = f_14 * shh_126[k]
                   + f_3 * pc_y[k] * sih_210[k];

        t_282[k] = f_3 * pc_z[k] * sih_210[k];

        t_283[k] = f_12 * shh_213[k]
                   + f_4 * sig0_153[k]
                   - f_5 * sig1_153[k]
                   + f_3 * pc_x[k] * sih_213[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pc_x, pc_y, shh_128, shh_215, shh_216, sig0_155, \
                         sig0_156, sig1_155, sig1_156, sih_212, sih_215, \
                         sih_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_14 * shh_128[k]
                   + f_3 * pc_y[k] * sih_212[k];

        t_285[k] = f_12 * shh_215[k]
                   + f_4 * sig0_155[k]
                   - f_5 * sig1_155[k]
                   + f_3 * pc_x[k] * sih_215[k];

        t_286[k] = f_12 * shh_216[k]
                   + f_6 * sig0_156[k]
                   - f_7 * sig1_156[k]
                   + f_3 * pc_x[k] * sih_216[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pc_x, pc_y, pc_z, shh_131, shh_219, sig0_159, \
                         sig1_159, sih_213, sih_215, sih_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_3 * pc_z[k] * sih_213[k];

        t_288[k] = f_14 * shh_131[k]
                   + f_3 * pc_y[k] * sih_215[k];

        t_289[k] = f_12 * shh_219[k]
                   + f_6 * sig0_159[k]
                   - f_7 * sig1_159[k]
                   + f_3 * pc_x[k] * sih_219[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pc_x, pc_z, shh_220, shh_222, sig0_160, \
                         sig0_162, sig1_160, sig1_162, sih_216, sih_220, \
                         sih_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_12 * shh_220[k]
                   + f_8 * sig0_160[k]
                   - f_9 * sig1_160[k]
                   + f_3 * pc_x[k] * sih_220[k];

        t_291[k] = f_3 * pc_z[k] * sih_216[k];

        t_292[k] = f_12 * shh_222[k]
                   + f_8 * sig0_162[k]
                   - f_9 * sig1_162[k]
                   + f_3 * pc_x[k] * sih_222[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pc_x, pc_y, shh_135, shh_224, shh_225, \
                         shh_226, sig0_164, sig1_164, sih_219, sih_224, sih_225, \
                         sih_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_14 * shh_135[k]
                   + f_3 * pc_y[k] * sih_219[k];

        t_294[k] = f_12 * shh_224[k]
                   + f_8 * sig0_164[k]
                   - f_9 * sig1_164[k]
                   + f_3 * pc_x[k] * sih_224[k];

        t_295[k] = f_12 * shh_225[k]
                   + f_3 * pc_x[k] * sih_225[k];

        t_296[k] = f_12 * shh_226[k]
                   + f_3 * pc_x[k] * sih_226[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, pc_x, shh_227, shh_228, shh_229, shh_230, \
                         sih_227, sih_228, sih_229, sih_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_12 * shh_227[k]
                   + f_3 * pc_x[k] * sih_227[k];

        t_298[k] = f_12 * shh_228[k]
                   + f_3 * pc_x[k] * sih_228[k];

        t_299[k] = f_12 * shh_229[k]
                   + f_3 * pc_x[k] * sih_229[k];

        t_300[k] = f_12 * shh_230[k]
                   + f_3 * pc_x[k] * sih_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, pc_y, pc_z, shh_141, shh_143, sig0_160, \
                         sig0_162, sig1_160, sig1_162, sih_225, \
                         sih_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * shh_141[k]
                   + f_1 * sig0_160[k]
                   - f_2 * sig1_160[k]
                   + f_3 * pc_y[k] * sih_225[k];

        t_302[k] = f_3 * pc_z[k] * sih_225[k];

        t_303[k] = f_14 * shh_143[k]
                   + f_4 * sig0_162[k]
                   - f_5 * sig1_162[k]
                   + f_3 * pc_y[k] * sih_227[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_y, pc_z, shh_144, shh_145, shh_146, \
                         sig0_163, sig0_164, sig1_163, sig1_164, sih_228, sih_229, \
                         sih_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_14 * shh_144[k]
                   + f_6 * sig0_163[k]
                   - f_7 * sig1_163[k]
                   + f_3 * pc_y[k] * sih_228[k];

        t_305[k] = f_14 * shh_145[k]
                   + f_8 * sig0_164[k]
                   - f_9 * sig1_164[k]
                   + f_3 * pc_y[k] * sih_229[k];

        t_306[k] = f_14 * shh_146[k]
                   + f_3 * pc_y[k] * sih_230[k];

        t_307[k] = f_1 * sig0_164[k]
                   - f_2 * sig1_164[k]
                   + f_3 * pc_z[k] * sih_230[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pb_z, pc_y, pc_z, shi0_168, shi0_171, \
                         shh_126, shh_147, shi1_168, shi1_171, \
                         sih_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pb_z[k] * shi0_168[k]
                   - f_10 * pc_z[k] * shi1_168[k];

        t_309[k] = f_13 * shh_147[k]
                   + f_3 * pc_y[k] * sih_231[k];

        t_310[k] = f_11 * shh_126[k]
                   + f_3 * pc_z[k] * sih_231[k];

        t_311[k] = pb_z[k] * shi0_171[k]
                   - f_10 * pc_z[k] * shi1_171[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pb_z, pc_x, pc_y, pc_z, shi0_174, shh_149, \
                         shh_236, shi1_174, sig0_170, sig1_170, sih_233, \
                         sih_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_13 * shh_149[k]
                   + f_3 * pc_y[k] * sih_233[k];

        t_313[k] = f_12 * shh_236[k]
                   + f_4 * sig0_170[k]
                   - f_5 * sig1_170[k]
                   + f_3 * pc_x[k] * sih_236[k];

        t_314[k] = pb_z[k] * shi0_174[k]
                   - f_10 * pc_z[k] * shi1_174[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pc_x, pc_y, pc_z, shh_129, shh_152, shh_240, \
                         sig0_174, sig1_174, sih_234, sih_236, \
                         sih_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_11 * shh_129[k]
                   + f_3 * pc_z[k] * sih_234[k];

        t_316[k] = f_13 * shh_152[k]
                   + f_3 * pc_y[k] * sih_236[k];

        t_317[k] = f_12 * shh_240[k]
                   + f_6 * sig0_174[k]
                   - f_7 * sig1_174[k]
                   + f_3 * pc_x[k] * sih_240[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, pb_z, pc_y, pc_z, shi0_178, shi0_180, \
                         shh_132, shh_133, shh_156, shi1_178, shi1_180, sih_237, \
                         sih_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pb_z[k] * shi0_178[k]
                   - f_10 * pc_z[k] * shi1_178[k];

        t_319[k] = f_11 * shh_132[k]
                   + f_3 * pc_z[k] * sih_237[k];

        t_320[k] = pb_z[k] * shi0_180[k]
                   + f_12 * shh_133[k]
                   - f_10 * pc_z[k] * shi1_180[k];

        t_321[k] = f_13 * shh_156[k]
                   + f_3 * pc_y[k] * sih_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, shh_245, shh_246, shh_247, shh_248, \
                         sig0_179, sig1_179, sih_245, sih_246, sih_247, \
                         sih_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_12 * shh_245[k]
                   + f_8 * sig0_179[k]
                   - f_9 * sig1_179[k]
                   + f_3 * pc_x[k] * sih_245[k];

        t_323[k] = f_12 * shh_246[k]
                   + f_3 * pc_x[k] * sih_246[k];

        t_324[k] = f_12 * shh_247[k]
                   + f_3 * pc_x[k] * sih_247[k];

        t_325[k] = f_12 * shh_248[k]
                   + f_3 * pc_x[k] * sih_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pb_z, pc_x, pc_z, shi0_189, shh_249, \
                         shh_250, shh_251, shi1_189, sih_249, sih_250, \
                         sih_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_12 * shh_249[k]
                   + f_3 * pc_x[k] * sih_249[k];

        t_327[k] = f_12 * shh_250[k]
                   + f_3 * pc_x[k] * sih_250[k];

        t_328[k] = f_12 * shh_251[k]
                   + f_3 * pc_x[k] * sih_251[k];

        t_329[k] = pb_z[k] * shi0_189[k]
                   - f_10 * pc_z[k] * shi1_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, pc_z, shh_141, shh_164, shh_165, sig0_177, \
                         sig0_178, sig1_177, sig1_178, sih_246, sih_248, \
                         sih_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_11 * shh_141[k]
                   + f_3 * pc_z[k] * sih_246[k];

        t_331[k] = f_13 * shh_164[k]
                   + f_4 * sig0_177[k]
                   - f_5 * sig1_177[k]
                   + f_3 * pc_y[k] * sih_248[k];

        t_332[k] = f_13 * shh_165[k]
                   + f_6 * sig0_178[k]
                   - f_7 * sig1_178[k]
                   + f_3 * pc_y[k] * sih_249[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_y, pc_z, shh_146, shh_166, shh_167, sig0_179, \
                         sig1_179, sih_250, sih_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_13 * shh_166[k]
                   + f_8 * sig0_179[k]
                   - f_9 * sig1_179[k]
                   + f_3 * pc_y[k] * sih_250[k];

        t_334[k] = f_13 * shh_167[k]
                   + f_3 * pc_y[k] * sih_251[k];

        t_335[k] = f_11 * shh_146[k]
                   + f_1 * sig0_179[k]
                   - f_2 * sig1_179[k]
                   + f_3 * pc_z[k] * sih_251[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, pc_z, shh_147, shh_168, shh_252, \
                         sig0_180, sig1_180, sih_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_12 * shh_252[k]
                   + f_1 * sig0_180[k]
                   - f_2 * sig1_180[k]
                   + f_3 * pc_x[k] * sih_252[k];

        t_337[k] = f_12 * shh_168[k]
                   + f_3 * pc_y[k] * sih_252[k];

        t_338[k] = f_12 * shh_147[k]
                   + f_3 * pc_z[k] * sih_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, shh_170, shh_255, shh_257, sig0_183, \
                         sig0_185, sig1_183, sig1_185, sih_254, sih_255, \
                         sih_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_12 * shh_255[k]
                   + f_4 * sig0_183[k]
                   - f_5 * sig1_183[k]
                   + f_3 * pc_x[k] * sih_255[k];

        t_340[k] = f_12 * shh_170[k]
                   + f_3 * pc_y[k] * sih_254[k];

        t_341[k] = f_12 * shh_257[k]
                   + f_4 * sig0_185[k]
                   - f_5 * sig1_185[k]
                   + f_3 * pc_x[k] * sih_257[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_y, pc_z, shh_150, shh_173, shh_258, \
                         sig0_186, sig1_186, sih_255, sih_257, \
                         sih_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_12 * shh_258[k]
                   + f_6 * sig0_186[k]
                   - f_7 * sig1_186[k]
                   + f_3 * pc_x[k] * sih_258[k];

        t_343[k] = f_12 * shh_150[k]
                   + f_3 * pc_z[k] * sih_255[k];

        t_344[k] = f_12 * shh_173[k]
                   + f_3 * pc_y[k] * sih_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, shh_153, shh_261, shh_262, sig0_189, \
                         sig0_190, sig1_189, sig1_190, sih_258, sih_261, \
                         sih_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_12 * shh_261[k]
                   + f_6 * sig0_189[k]
                   - f_7 * sig1_189[k]
                   + f_3 * pc_x[k] * sih_261[k];

        t_346[k] = f_12 * shh_262[k]
                   + f_8 * sig0_190[k]
                   - f_9 * sig1_190[k]
                   + f_3 * pc_x[k] * sih_262[k];

        t_347[k] = f_12 * shh_153[k]
                   + f_3 * pc_z[k] * sih_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, shh_177, shh_264, shh_266, sig0_192, \
                         sig0_194, sig1_192, sig1_194, sih_261, sih_264, \
                         sih_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_12 * shh_264[k]
                   + f_8 * sig0_192[k]
                   - f_9 * sig1_192[k]
                   + f_3 * pc_x[k] * sih_264[k];

        t_349[k] = f_12 * shh_177[k]
                   + f_3 * pc_y[k] * sih_261[k];

        t_350[k] = f_12 * shh_266[k]
                   + f_8 * sig0_194[k]
                   - f_9 * sig1_194[k]
                   + f_3 * pc_x[k] * sih_266[k];
    }
}

static auto
compute_prim_sii_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shi0,
                                                          const size_t shh, const size_t shi1,
                                                          const size_t sig0, const size_t sig1,
                                                          const size_t sih, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_15 = 2.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shi0_252 = buffer.data(shi0 + 252);
    const auto *shi0_255 = buffer.data(shi0 + 255);
    const auto *shi0_257 = buffer.data(shi0 + 257);
    const auto *shi0_258 = buffer.data(shi0 + 258);
    const auto *shi0_261 = buffer.data(shi0 + 261);
    const auto *shi0_262 = buffer.data(shi0 + 262);
    const auto *shi0_264 = buffer.data(shi0 + 264);
    const auto *shi0_266 = buffer.data(shi0 + 266);
    const auto *shi0_279 = buffer.data(shi0 + 279);
    const auto *shi0_280 = buffer.data(shi0 + 280);
    const auto *shi0_283 = buffer.data(shi0 + 283);
    const auto *shi0_286 = buffer.data(shi0 + 286);
    const auto *shi0_290 = buffer.data(shi0 + 290);
    const auto *shi0_420 = buffer.data(shi0 + 420);
    const auto *shi0_423 = buffer.data(shi0 + 423);
    const auto *shi0_425 = buffer.data(shi0 + 425);
    const auto *shi0_426 = buffer.data(shi0 + 426);
    const auto *shi0_429 = buffer.data(shi0 + 429);
    const auto *shi0_430 = buffer.data(shi0 + 430);
    const auto *shi0_432 = buffer.data(shi0 + 432);
    const auto *shi0_434 = buffer.data(shi0 + 434);
    const auto *shi0_441 = buffer.data(shi0 + 441);
    const auto *shi0_443 = buffer.data(shi0 + 443);
    const auto *shi0_444 = buffer.data(shi0 + 444);
    const auto *shi0_445 = buffer.data(shi0 + 445);
    const auto *shi0_447 = buffer.data(shi0 + 447);
    const auto *shi0_453 = buffer.data(shi0 + 453);
    const auto *shi0_457 = buffer.data(shi0 + 457);
    const auto *shi0_460 = buffer.data(shi0 + 460);
    const auto *shi0_462 = buffer.data(shi0 + 462);

    const auto *shh_162 = buffer.data(shh + 162);
    const auto *shh_167 = buffer.data(shh + 167);
    const auto *shh_168 = buffer.data(shh + 168);
    const auto *shh_171 = buffer.data(shh + 171);
    const auto *shh_174 = buffer.data(shh + 174);
    const auto *shh_183 = buffer.data(shh + 183);
    const auto *shh_185 = buffer.data(shh + 185);
    const auto *shh_186 = buffer.data(shh + 186);
    const auto *shh_187 = buffer.data(shh + 187);
    const auto *shh_188 = buffer.data(shh + 188);
    const auto *shh_189 = buffer.data(shh + 189);
    const auto *shh_190 = buffer.data(shh + 190);
    const auto *shh_191 = buffer.data(shh + 191);
    const auto *shh_192 = buffer.data(shh + 192);
    const auto *shh_194 = buffer.data(shh + 194);
    const auto *shh_195 = buffer.data(shh + 195);
    const auto *shh_197 = buffer.data(shh + 197);
    const auto *shh_198 = buffer.data(shh + 198);
    const auto *shh_204 = buffer.data(shh + 204);
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
    const auto *shh_230 = buffer.data(shh + 230);
    const auto *shh_231 = buffer.data(shh + 231);
    const auto *shh_233 = buffer.data(shh + 233);
    const auto *shh_236 = buffer.data(shh + 236);
    const auto *shh_240 = buffer.data(shh + 240);
    const auto *shh_267 = buffer.data(shh + 267);
    const auto *shh_268 = buffer.data(shh + 268);
    const auto *shh_269 = buffer.data(shh + 269);
    const auto *shh_270 = buffer.data(shh + 270);
    const auto *shh_271 = buffer.data(shh + 271);
    const auto *shh_272 = buffer.data(shh + 272);
    const auto *shh_288 = buffer.data(shh + 288);
    const auto *shh_289 = buffer.data(shh + 289);
    const auto *shh_290 = buffer.data(shh + 290);
    const auto *shh_291 = buffer.data(shh + 291);
    const auto *shh_292 = buffer.data(shh + 292);
    const auto *shh_293 = buffer.data(shh + 293);
    const auto *shh_294 = buffer.data(shh + 294);
    const auto *shh_297 = buffer.data(shh + 297);
    const auto *shh_299 = buffer.data(shh + 299);
    const auto *shh_300 = buffer.data(shh + 300);
    const auto *shh_303 = buffer.data(shh + 303);
    const auto *shh_304 = buffer.data(shh + 304);
    const auto *shh_306 = buffer.data(shh + 306);
    const auto *shh_308 = buffer.data(shh + 308);
    const auto *shh_309 = buffer.data(shh + 309);
    const auto *shh_310 = buffer.data(shh + 310);
    const auto *shh_311 = buffer.data(shh + 311);
    const auto *shh_312 = buffer.data(shh + 312);
    const auto *shh_313 = buffer.data(shh + 313);
    const auto *shh_314 = buffer.data(shh + 314);
    const auto *shh_315 = buffer.data(shh + 315);
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
    const auto *shh_341 = buffer.data(shh + 341);
    const auto *shh_345 = buffer.data(shh + 345);
    const auto *shh_348 = buffer.data(shh + 348);
    const auto *shh_350 = buffer.data(shh + 350);
    const auto *shh_351 = buffer.data(shh + 351);
    const auto *shh_352 = buffer.data(shh + 352);
    const auto *shh_353 = buffer.data(shh + 353);
    const auto *shh_354 = buffer.data(shh + 354);
    const auto *shh_355 = buffer.data(shh + 355);
    const auto *shh_356 = buffer.data(shh + 356);

    const auto *shi1_252 = buffer.data(shi1 + 252);
    const auto *shi1_255 = buffer.data(shi1 + 255);
    const auto *shi1_257 = buffer.data(shi1 + 257);
    const auto *shi1_258 = buffer.data(shi1 + 258);
    const auto *shi1_261 = buffer.data(shi1 + 261);
    const auto *shi1_262 = buffer.data(shi1 + 262);
    const auto *shi1_264 = buffer.data(shi1 + 264);
    const auto *shi1_266 = buffer.data(shi1 + 266);
    const auto *shi1_279 = buffer.data(shi1 + 279);
    const auto *shi1_280 = buffer.data(shi1 + 280);
    const auto *shi1_283 = buffer.data(shi1 + 283);
    const auto *shi1_286 = buffer.data(shi1 + 286);
    const auto *shi1_290 = buffer.data(shi1 + 290);
    const auto *shi1_420 = buffer.data(shi1 + 420);
    const auto *shi1_423 = buffer.data(shi1 + 423);
    const auto *shi1_425 = buffer.data(shi1 + 425);
    const auto *shi1_426 = buffer.data(shi1 + 426);
    const auto *shi1_429 = buffer.data(shi1 + 429);
    const auto *shi1_430 = buffer.data(shi1 + 430);
    const auto *shi1_432 = buffer.data(shi1 + 432);
    const auto *shi1_434 = buffer.data(shi1 + 434);
    const auto *shi1_441 = buffer.data(shi1 + 441);
    const auto *shi1_443 = buffer.data(shi1 + 443);
    const auto *shi1_444 = buffer.data(shi1 + 444);
    const auto *shi1_445 = buffer.data(shi1 + 445);
    const auto *shi1_447 = buffer.data(shi1 + 447);
    const auto *shi1_453 = buffer.data(shi1 + 453);
    const auto *shi1_457 = buffer.data(shi1 + 457);
    const auto *shi1_460 = buffer.data(shi1 + 460);
    const auto *shi1_462 = buffer.data(shi1 + 462);

    const auto *sig0_190 = buffer.data(sig0 + 190);
    const auto *sig0_192 = buffer.data(sig0 + 192);
    const auto *sig0_193 = buffer.data(sig0 + 193);
    const auto *sig0_194 = buffer.data(sig0 + 194);
    const auto *sig0_205 = buffer.data(sig0 + 205);
    const auto *sig0_207 = buffer.data(sig0 + 207);
    const auto *sig0_208 = buffer.data(sig0 + 208);
    const auto *sig0_209 = buffer.data(sig0 + 209);
    const auto *sig0_210 = buffer.data(sig0 + 210);
    const auto *sig0_213 = buffer.data(sig0 + 213);
    const auto *sig0_215 = buffer.data(sig0 + 215);
    const auto *sig0_216 = buffer.data(sig0 + 216);
    const auto *sig0_219 = buffer.data(sig0 + 219);
    const auto *sig0_220 = buffer.data(sig0 + 220);
    const auto *sig0_222 = buffer.data(sig0 + 222);
    const auto *sig0_223 = buffer.data(sig0 + 223);
    const auto *sig0_224 = buffer.data(sig0 + 224);

    const auto *sig1_190 = buffer.data(sig1 + 190);
    const auto *sig1_192 = buffer.data(sig1 + 192);
    const auto *sig1_193 = buffer.data(sig1 + 193);
    const auto *sig1_194 = buffer.data(sig1 + 194);
    const auto *sig1_205 = buffer.data(sig1 + 205);
    const auto *sig1_207 = buffer.data(sig1 + 207);
    const auto *sig1_208 = buffer.data(sig1 + 208);
    const auto *sig1_209 = buffer.data(sig1 + 209);
    const auto *sig1_210 = buffer.data(sig1 + 210);
    const auto *sig1_213 = buffer.data(sig1 + 213);
    const auto *sig1_215 = buffer.data(sig1 + 215);
    const auto *sig1_216 = buffer.data(sig1 + 216);
    const auto *sig1_219 = buffer.data(sig1 + 219);
    const auto *sig1_220 = buffer.data(sig1 + 220);
    const auto *sig1_222 = buffer.data(sig1 + 222);
    const auto *sig1_223 = buffer.data(sig1 + 223);
    const auto *sig1_224 = buffer.data(sig1 + 224);

    const auto *sih_267 = buffer.data(sih + 267);
    const auto *sih_268 = buffer.data(sih + 268);
    const auto *sih_269 = buffer.data(sih + 269);
    const auto *sih_270 = buffer.data(sih + 270);
    const auto *sih_271 = buffer.data(sih + 271);
    const auto *sih_272 = buffer.data(sih + 272);
    const auto *sih_273 = buffer.data(sih + 273);
    const auto *sih_275 = buffer.data(sih + 275);
    const auto *sih_276 = buffer.data(sih + 276);
    const auto *sih_278 = buffer.data(sih + 278);
    const auto *sih_279 = buffer.data(sih + 279);
    const auto *sih_282 = buffer.data(sih + 282);
    const auto *sih_288 = buffer.data(sih + 288);
    const auto *sih_289 = buffer.data(sih + 289);
    const auto *sih_290 = buffer.data(sih + 290);
    const auto *sih_291 = buffer.data(sih + 291);
    const auto *sih_292 = buffer.data(sih + 292);
    const auto *sih_293 = buffer.data(sih + 293);
    const auto *sih_294 = buffer.data(sih + 294);
    const auto *sih_296 = buffer.data(sih + 296);
    const auto *sih_297 = buffer.data(sih + 297);
    const auto *sih_299 = buffer.data(sih + 299);
    const auto *sih_300 = buffer.data(sih + 300);
    const auto *sih_303 = buffer.data(sih + 303);
    const auto *sih_304 = buffer.data(sih + 304);
    const auto *sih_306 = buffer.data(sih + 306);
    const auto *sih_308 = buffer.data(sih + 308);
    const auto *sih_309 = buffer.data(sih + 309);
    const auto *sih_310 = buffer.data(sih + 310);
    const auto *sih_311 = buffer.data(sih + 311);
    const auto *sih_312 = buffer.data(sih + 312);
    const auto *sih_313 = buffer.data(sih + 313);
    const auto *sih_314 = buffer.data(sih + 314);
    const auto *sih_315 = buffer.data(sih + 315);
    const auto *sih_317 = buffer.data(sih + 317);
    const auto *sih_318 = buffer.data(sih + 318);
    const auto *sih_320 = buffer.data(sih + 320);
    const auto *sih_321 = buffer.data(sih + 321);
    const auto *sih_324 = buffer.data(sih + 324);
    const auto *sih_330 = buffer.data(sih + 330);
    const auto *sih_331 = buffer.data(sih + 331);
    const auto *sih_332 = buffer.data(sih + 332);
    const auto *sih_333 = buffer.data(sih + 333);
    const auto *sih_334 = buffer.data(sih + 334);
    const auto *sih_335 = buffer.data(sih + 335);
    const auto *sih_336 = buffer.data(sih + 336);
    const auto *sih_338 = buffer.data(sih + 338);
    const auto *sih_339 = buffer.data(sih + 339);
    const auto *sih_341 = buffer.data(sih + 341);
    const auto *sih_342 = buffer.data(sih + 342);
    const auto *sih_345 = buffer.data(sih + 345);
    const auto *sih_351 = buffer.data(sih + 351);
    const auto *sih_352 = buffer.data(sih + 352);
    const auto *sih_353 = buffer.data(sih + 353);
    const auto *sih_354 = buffer.data(sih + 354);
    const auto *sih_355 = buffer.data(sih + 355);
    const auto *sih_356 = buffer.data(sih + 356);

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pc_x, shh_267, shh_268, shh_269, \
                         shh_270, shh_271, sih_267, sih_268, sih_269, sih_270, \
                         sih_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_12 * shh_267[k]
                   + f_3 * pc_x[k] * sih_267[k];

        t_352[k] = f_12 * shh_268[k]
                   + f_3 * pc_x[k] * sih_268[k];

        t_353[k] = f_12 * shh_269[k]
                   + f_3 * pc_x[k] * sih_269[k];

        t_354[k] = f_12 * shh_270[k]
                   + f_3 * pc_x[k] * sih_270[k];

        t_355[k] = f_12 * shh_271[k]
                   + f_3 * pc_x[k] * sih_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, shh_162, shh_183, shh_272, \
                         sig0_190, sig1_190, sih_267, sih_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_12 * shh_272[k]
                   + f_3 * pc_x[k] * sih_272[k];

        t_357[k] = f_12 * shh_183[k]
                   + f_1 * sig0_190[k]
                   - f_2 * sig1_190[k]
                   + f_3 * pc_y[k] * sih_267[k];

        t_358[k] = f_12 * shh_162[k]
                   + f_3 * pc_z[k] * sih_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, shh_185, shh_186, shh_187, sig0_192, \
                         sig0_193, sig0_194, sig1_192, sig1_193, sig1_194, sih_269, sih_270, \
                         sih_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * shh_185[k]
                   + f_4 * sig0_192[k]
                   - f_5 * sig1_192[k]
                   + f_3 * pc_y[k] * sih_269[k];

        t_360[k] = f_12 * shh_186[k]
                   + f_6 * sig0_193[k]
                   - f_7 * sig1_193[k]
                   + f_3 * pc_y[k] * sih_270[k];

        t_361[k] = f_12 * shh_187[k]
                   + f_8 * sig0_194[k]
                   - f_9 * sig1_194[k]
                   + f_3 * pc_y[k] * sih_271[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pb_y, pc_y, pc_z, shi0_252, shh_167, \
                         shh_188, shh_189, shi1_252, sig0_194, sig1_194, sih_272, \
                         sih_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * shh_188[k]
                   + f_3 * pc_y[k] * sih_272[k];

        t_363[k] = f_12 * shh_167[k]
                   + f_1 * sig0_194[k]
                   - f_2 * sig1_194[k]
                   + f_3 * pc_z[k] * sih_272[k];

        t_364[k] = pb_y[k] * shi0_252[k]
                   - f_10 * pc_y[k] * shi1_252[k];

        t_365[k] = f_11 * shh_189[k]
                   + f_3 * pc_y[k] * sih_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_y, pc_y, pc_z, shi0_255, shi0_257, \
                         shh_168, shh_190, shh_191, shi1_255, shi1_257, sih_273, \
                         sih_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * shh_168[k]
                   + f_3 * pc_z[k] * sih_273[k];

        t_367[k] = pb_y[k] * shi0_255[k]
                   + f_12 * shh_190[k]
                   - f_10 * pc_y[k] * shi1_255[k];

        t_368[k] = f_11 * shh_191[k]
                   + f_3 * pc_y[k] * sih_275[k];

        t_369[k] = pb_y[k] * shi0_257[k]
                   - f_10 * pc_y[k] * shi1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pb_y, pc_y, pc_z, shi0_258, shi0_261, \
                         shh_171, shh_192, shh_194, shi1_258, shi1_261, sih_276, \
                         sih_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pb_y[k] * shi0_258[k]
                   + f_13 * shh_192[k]
                   - f_10 * pc_y[k] * shi1_258[k];

        t_371[k] = f_13 * shh_171[k]
                   + f_3 * pc_z[k] * sih_276[k];

        t_372[k] = f_11 * shh_194[k]
                   + f_3 * pc_y[k] * sih_278[k];

        t_373[k] = pb_y[k] * shi0_261[k]
                   - f_10 * pc_y[k] * shi1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pb_y, pc_y, pc_z, shi0_262, shi0_264, shh_174, \
                         shh_195, shh_197, shi1_262, shi1_264, \
                         sih_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pb_y[k] * shi0_262[k]
                   + f_14 * shh_195[k]
                   - f_10 * pc_y[k] * shi1_262[k];

        t_375[k] = f_13 * shh_174[k]
                   + f_3 * pc_z[k] * sih_279[k];

        t_376[k] = pb_y[k] * shi0_264[k]
                   + f_12 * shh_197[k]
                   - f_10 * pc_y[k] * shi1_264[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pb_y, pc_x, pc_y, shi0_266, shh_198, \
                         shh_288, shh_289, shi1_266, sih_282, sih_288, \
                         sih_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * shh_198[k]
                   + f_3 * pc_y[k] * sih_282[k];

        t_378[k] = pb_y[k] * shi0_266[k]
                   - f_10 * pc_y[k] * shi1_266[k];

        t_379[k] = f_12 * shh_288[k]
                   + f_3 * pc_x[k] * sih_288[k];

        t_380[k] = f_12 * shh_289[k]
                   + f_3 * pc_x[k] * sih_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, shh_290, shh_291, shh_292, shh_293, \
                         sih_290, sih_291, sih_292, sih_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_12 * shh_290[k]
                   + f_3 * pc_x[k] * sih_290[k];

        t_382[k] = f_12 * shh_291[k]
                   + f_3 * pc_x[k] * sih_291[k];

        t_383[k] = f_12 * shh_292[k]
                   + f_3 * pc_x[k] * sih_292[k];

        t_384[k] = f_12 * shh_293[k]
                   + f_3 * pc_x[k] * sih_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, shh_183, shh_204, shh_206, sig0_205, \
                         sig0_207, sig1_205, sig1_207, sih_288, \
                         sih_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * shh_204[k]
                   + f_1 * sig0_205[k]
                   - f_2 * sig1_205[k]
                   + f_3 * pc_y[k] * sih_288[k];

        t_386[k] = f_13 * shh_183[k]
                   + f_3 * pc_z[k] * sih_288[k];

        t_387[k] = f_11 * shh_206[k]
                   + f_4 * sig0_207[k]
                   - f_5 * sig1_207[k]
                   + f_3 * pc_y[k] * sih_290[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, shh_207, shh_208, shh_209, sig0_208, \
                         sig0_209, sig1_208, sig1_209, sih_291, sih_292, \
                         sih_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * shh_207[k]
                   + f_6 * sig0_208[k]
                   - f_7 * sig1_208[k]
                   + f_3 * pc_y[k] * sih_291[k];

        t_389[k] = f_11 * shh_208[k]
                   + f_8 * sig0_209[k]
                   - f_9 * sig1_209[k]
                   + f_3 * pc_y[k] * sih_292[k];

        t_390[k] = f_11 * shh_209[k]
                   + f_3 * pc_y[k] * sih_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pb_y, pc_x, pc_y, pc_z, shi0_279, \
                         shh_189, shh_294, shi1_279, sig0_210, sig1_210, \
                         sih_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pb_y[k] * shi0_279[k]
                   - f_10 * pc_y[k] * shi1_279[k];

        t_392[k] = f_12 * shh_294[k]
                   + f_1 * sig0_210[k]
                   - f_2 * sig1_210[k]
                   + f_3 * pc_x[k] * sih_294[k];

        t_393[k] = f_3 * pc_y[k] * sih_294[k];

        t_394[k] = f_14 * shh_189[k]
                   + f_3 * pc_z[k] * sih_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_x, pc_y, shh_297, shh_299, sig0_213, \
                         sig0_215, sig1_213, sig1_215, sih_296, sih_297, \
                         sih_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_12 * shh_297[k]
                   + f_4 * sig0_213[k]
                   - f_5 * sig1_213[k]
                   + f_3 * pc_x[k] * sih_297[k];

        t_396[k] = f_3 * pc_y[k] * sih_296[k];

        t_397[k] = f_12 * shh_299[k]
                   + f_4 * sig0_215[k]
                   - f_5 * sig1_215[k]
                   + f_3 * pc_x[k] * sih_299[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_x, pc_y, pc_z, shh_192, shh_300, sig0_216, \
                         sig1_216, sih_297, sih_299, sih_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_12 * shh_300[k]
                   + f_6 * sig0_216[k]
                   - f_7 * sig1_216[k]
                   + f_3 * pc_x[k] * sih_300[k];

        t_399[k] = f_14 * shh_192[k]
                   + f_3 * pc_z[k] * sih_297[k];

        t_400[k] = f_3 * pc_y[k] * sih_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_z, shh_195, shh_303, shh_304, sig0_219, \
                         sig0_220, sig1_219, sig1_220, sih_300, sih_303, \
                         sih_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_12 * shh_303[k]
                   + f_6 * sig0_219[k]
                   - f_7 * sig1_219[k]
                   + f_3 * pc_x[k] * sih_303[k];

        t_402[k] = f_12 * shh_304[k]
                   + f_8 * sig0_220[k]
                   - f_9 * sig1_220[k]
                   + f_3 * pc_x[k] * sih_304[k];

        t_403[k] = f_14 * shh_195[k]
                   + f_3 * pc_z[k] * sih_300[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pc_x, pc_y, shh_306, shh_308, sig0_222, \
                         sig0_224, sig1_222, sig1_224, sih_303, sih_306, \
                         sih_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_12 * shh_306[k]
                   + f_8 * sig0_222[k]
                   - f_9 * sig1_222[k]
                   + f_3 * pc_x[k] * sih_306[k];

        t_405[k] = f_3 * pc_y[k] * sih_303[k];

        t_406[k] = f_12 * shh_308[k]
                   + f_8 * sig0_224[k]
                   - f_9 * sig1_224[k]
                   + f_3 * pc_x[k] * sih_308[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, pc_x, shh_309, shh_310, shh_311, \
                         shh_312, shh_313, sih_309, sih_310, sih_311, sih_312, \
                         sih_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_12 * shh_309[k]
                   + f_3 * pc_x[k] * sih_309[k];

        t_408[k] = f_12 * shh_310[k]
                   + f_3 * pc_x[k] * sih_310[k];

        t_409[k] = f_12 * shh_311[k]
                   + f_3 * pc_x[k] * sih_311[k];

        t_410[k] = f_12 * shh_312[k]
                   + f_3 * pc_x[k] * sih_312[k];

        t_411[k] = f_12 * shh_313[k]
                   + f_3 * pc_x[k] * sih_313[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pc_x, pc_y, pc_z, shh_204, shh_314, \
                         sig0_220, sig0_222, sig1_220, sig1_222, sih_309, sih_311, \
                         sih_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_12 * shh_314[k]
                   + f_3 * pc_x[k] * sih_314[k];

        t_413[k] = f_1 * sig0_220[k]
                   - f_2 * sig1_220[k]
                   + f_3 * pc_y[k] * sih_309[k];

        t_414[k] = f_14 * shh_204[k]
                   + f_3 * pc_z[k] * sih_309[k];

        t_415[k] = f_4 * sig0_222[k]
                   - f_5 * sig1_222[k]
                   + f_3 * pc_y[k] * sih_311[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_y, pc_z, shh_209, sig0_223, sig0_224, \
                         sig1_223, sig1_224, sih_312, sih_313, \
                         sih_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * sig0_223[k]
                   - f_7 * sig1_223[k]
                   + f_3 * pc_y[k] * sih_312[k];

        t_417[k] = f_8 * sig0_224[k]
                   - f_9 * sig1_224[k]
                   + f_3 * pc_y[k] * sih_313[k];

        t_418[k] = f_3 * pc_y[k] * sih_314[k];

        t_419[k] = f_14 * shh_209[k]
                   + f_1 * sig0_224[k]
                   - f_2 * sig1_224[k]
                   + f_3 * pc_z[k] * sih_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pb_x, pc_x, pc_y, pc_z, shi0_420, \
                         shi0_423, shh_210, shh_315, shh_318, shi1_420, shi1_423, \
                         sih_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = pb_x[k] * shi0_420[k]
                   + f_0 * shh_315[k]
                   - f_10 * pc_x[k] * shi1_420[k];

        t_421[k] = f_15 * shh_210[k]
                   + f_3 * pc_y[k] * sih_315[k];

        t_422[k] = f_3 * pc_z[k] * sih_315[k];

        t_423[k] = pb_x[k] * shi0_423[k]
                   + f_14 * shh_318[k]
                   - f_10 * pc_x[k] * shi1_423[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pb_x, pc_x, pc_y, shi0_425, shi0_426, shh_212, \
                         shh_320, shh_321, shi1_425, shi1_426, \
                         sih_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_15 * shh_212[k]
                   + f_3 * pc_y[k] * sih_317[k];

        t_425[k] = pb_x[k] * shi0_425[k]
                   + f_14 * shh_320[k]
                   - f_10 * pc_x[k] * shi1_425[k];

        t_426[k] = pb_x[k] * shi0_426[k]
                   + f_13 * shh_321[k]
                   - f_10 * pc_x[k] * shi1_426[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pb_x, pc_x, pc_y, pc_z, shi0_429, shh_215, \
                         shh_324, shi1_429, sih_318, sih_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_3 * pc_z[k] * sih_318[k];

        t_428[k] = f_15 * shh_215[k]
                   + f_3 * pc_y[k] * sih_320[k];

        t_429[k] = pb_x[k] * shi0_429[k]
                   + f_13 * shh_324[k]
                   - f_10 * pc_x[k] * shi1_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pb_x, pc_x, pc_z, shi0_430, shi0_432, shh_325, \
                         shh_327, shi1_430, shi1_432, sih_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = pb_x[k] * shi0_430[k]
                   + f_12 * shh_325[k]
                   - f_10 * pc_x[k] * shi1_430[k];

        t_431[k] = f_3 * pc_z[k] * sih_321[k];

        t_432[k] = pb_x[k] * shi0_432[k]
                   + f_12 * shh_327[k]
                   - f_10 * pc_x[k] * shi1_432[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pb_x, pc_x, pc_y, shi0_434, shh_219, \
                         shh_329, shh_330, shh_331, shi1_434, sih_324, sih_330, \
                         sih_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_15 * shh_219[k]
                   + f_3 * pc_y[k] * sih_324[k];

        t_434[k] = pb_x[k] * shi0_434[k]
                   + f_12 * shh_329[k]
                   - f_10 * pc_x[k] * shi1_434[k];

        t_435[k] = f_11 * shh_330[k]
                   + f_3 * pc_x[k] * sih_330[k];

        t_436[k] = f_11 * shh_331[k]
                   + f_3 * pc_x[k] * sih_331[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_x, shh_332, shh_333, shh_334, shh_335, \
                         sih_332, sih_333, sih_334, sih_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_11 * shh_332[k]
                   + f_3 * pc_x[k] * sih_332[k];

        t_438[k] = f_11 * shh_333[k]
                   + f_3 * pc_x[k] * sih_333[k];

        t_439[k] = f_11 * shh_334[k]
                   + f_3 * pc_x[k] * sih_334[k];

        t_440[k] = f_11 * shh_335[k]
                   + f_3 * pc_x[k] * sih_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pb_x, pc_x, pc_z, shi0_441, shi0_443, \
                         shi0_444, shi1_441, shi1_443, shi1_444, \
                         sih_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = pb_x[k] * shi0_441[k]
                   - f_10 * pc_x[k] * shi1_441[k];

        t_442[k] = f_3 * pc_z[k] * sih_330[k];

        t_443[k] = pb_x[k] * shi0_443[k]
                   - f_10 * pc_x[k] * shi1_443[k];

        t_444[k] = pb_x[k] * shi0_444[k]
                   - f_10 * pc_x[k] * shi1_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, pb_x, pc_x, pc_y, shi0_445, shi0_447, shh_230, \
                         shi1_445, shi1_447, sih_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = pb_x[k] * shi0_445[k]
                   - f_10 * pc_x[k] * shi1_445[k];

        t_446[k] = f_15 * shh_230[k]
                   + f_3 * pc_y[k] * sih_335[k];

        t_447[k] = pb_x[k] * shi0_447[k]
                   - f_10 * pc_x[k] * shi1_447[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, pb_z, pc_y, pc_z, shi0_280, shi0_283, \
                         shh_210, shh_231, shi1_280, shi1_283, \
                         sih_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = pb_z[k] * shi0_280[k]
                   - f_10 * pc_z[k] * shi1_280[k];

        t_449[k] = f_14 * shh_231[k]
                   + f_3 * pc_y[k] * sih_336[k];

        t_450[k] = f_11 * shh_210[k]
                   + f_3 * pc_z[k] * sih_336[k];

        t_451[k] = pb_z[k] * shi0_283[k]
                   - f_10 * pc_z[k] * shi1_283[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, pb_x, pb_z, pc_x, pc_y, pc_z, shi0_286, \
                         shi0_453, shh_233, shh_341, shi1_286, shi1_453, \
                         sih_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_14 * shh_233[k]
                   + f_3 * pc_y[k] * sih_338[k];

        t_453[k] = pb_x[k] * shi0_453[k]
                   + f_14 * shh_341[k]
                   - f_10 * pc_x[k] * shi1_453[k];

        t_454[k] = pb_z[k] * shi0_286[k]
                   - f_10 * pc_z[k] * shi1_286[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pb_x, pc_x, pc_y, pc_z, shi0_457, shh_213, \
                         shh_236, shh_345, shi1_457, sih_339, sih_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_11 * shh_213[k]
                   + f_3 * pc_z[k] * sih_339[k];

        t_456[k] = f_14 * shh_236[k]
                   + f_3 * pc_y[k] * sih_341[k];

        t_457[k] = pb_x[k] * shi0_457[k]
                   + f_13 * shh_345[k]
                   - f_10 * pc_x[k] * shi1_457[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, pb_x, pb_z, pc_x, pc_z, shi0_290, shi0_460, \
                         shh_216, shh_348, shi1_290, shi1_460, \
                         sih_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = pb_z[k] * shi0_290[k]
                   - f_10 * pc_z[k] * shi1_290[k];

        t_459[k] = f_11 * shh_216[k]
                   + f_3 * pc_z[k] * sih_342[k];

        t_460[k] = pb_x[k] * shi0_460[k]
                   + f_12 * shh_348[k]
                   - f_10 * pc_x[k] * shi1_460[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pb_x, pc_x, pc_y, shi0_462, shh_240, \
                         shh_350, shh_351, shh_352, shi1_462, sih_345, sih_351, \
                         sih_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_14 * shh_240[k]
                   + f_3 * pc_y[k] * sih_345[k];

        t_462[k] = pb_x[k] * shi0_462[k]
                   + f_12 * shh_350[k]
                   - f_10 * pc_x[k] * shi1_462[k];

        t_463[k] = f_11 * shh_351[k]
                   + f_3 * pc_x[k] * sih_351[k];

        t_464[k] = f_11 * shh_352[k]
                   + f_3 * pc_x[k] * sih_352[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, pc_x, shh_353, shh_354, shh_355, shh_356, \
                         sih_353, sih_354, sih_355, sih_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_11 * shh_353[k]
                   + f_3 * pc_x[k] * sih_353[k];

        t_466[k] = f_11 * shh_354[k]
                   + f_3 * pc_x[k] * sih_354[k];

        t_467[k] = f_11 * shh_355[k]
                   + f_3 * pc_x[k] * sih_355[k];

        t_468[k] = f_11 * shh_356[k]
                   + f_3 * pc_x[k] * sih_356[k];
    }
}

static auto
compute_prim_sii_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shi0,
                                                          const size_t shh, const size_t shi1,
                                                          const size_t sig0, const size_t sig1,
                                                          const size_t sih, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 2.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shi0_392 = buffer.data(shi0 + 392);
    const auto *shi0_397 = buffer.data(shi0 + 397);
    const auto *shi0_401 = buffer.data(shi0 + 401);
    const auto *shi0_406 = buffer.data(shi0 + 406);
    const auto *shi0_469 = buffer.data(shi0 + 469);
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
    const auto *shi0_497 = buffer.data(shi0 + 497);
    const auto *shi0_499 = buffer.data(shi0 + 499);
    const auto *shi0_500 = buffer.data(shi0 + 500);
    const auto *shi0_501 = buffer.data(shi0 + 501);
    const auto *shi0_503 = buffer.data(shi0 + 503);
    const auto *shi0_504 = buffer.data(shi0 + 504);
    const auto *shi0_507 = buffer.data(shi0 + 507);
    const auto *shi0_509 = buffer.data(shi0 + 509);
    const auto *shi0_510 = buffer.data(shi0 + 510);
    const auto *shi0_513 = buffer.data(shi0 + 513);
    const auto *shi0_514 = buffer.data(shi0 + 514);
    const auto *shi0_516 = buffer.data(shi0 + 516);
    const auto *shi0_518 = buffer.data(shi0 + 518);
    const auto *shi0_525 = buffer.data(shi0 + 525);
    const auto *shi0_527 = buffer.data(shi0 + 527);
    const auto *shi0_528 = buffer.data(shi0 + 528);
    const auto *shi0_529 = buffer.data(shi0 + 529);
    const auto *shi0_531 = buffer.data(shi0 + 531);
    const auto *shi0_535 = buffer.data(shi0 + 535);
    const auto *shi0_538 = buffer.data(shi0 + 538);
    const auto *shi0_542 = buffer.data(shi0 + 542);
    const auto *shi0_544 = buffer.data(shi0 + 544);
    const auto *shi0_553 = buffer.data(shi0 + 553);
    const auto *shi0_555 = buffer.data(shi0 + 555);
    const auto *shi0_556 = buffer.data(shi0 + 556);
    const auto *shi0_557 = buffer.data(shi0 + 557);
    const auto *shi0_559 = buffer.data(shi0 + 559);
    const auto *shi0_560 = buffer.data(shi0 + 560);
    const auto *shi0_563 = buffer.data(shi0 + 563);
    const auto *shi0_565 = buffer.data(shi0 + 565);
    const auto *shi0_566 = buffer.data(shi0 + 566);
    const auto *shi0_569 = buffer.data(shi0 + 569);
    const auto *shi0_570 = buffer.data(shi0 + 570);
    const auto *shi0_572 = buffer.data(shi0 + 572);
    const auto *shi0_574 = buffer.data(shi0 + 574);
    const auto *shi0_581 = buffer.data(shi0 + 581);
    const auto *shi0_583 = buffer.data(shi0 + 583);
    const auto *shi0_584 = buffer.data(shi0 + 584);
    const auto *shi0_585 = buffer.data(shi0 + 585);
    const auto *shi0_587 = buffer.data(shi0 + 587);

    const auto *shh_225 = buffer.data(shh + 225);
    const auto *shh_231 = buffer.data(shh + 231);
    const auto *shh_234 = buffer.data(shh + 234);
    const auto *shh_237 = buffer.data(shh + 237);
    const auto *shh_246 = buffer.data(shh + 246);
    const auto *shh_251 = buffer.data(shh + 251);
    const auto *shh_252 = buffer.data(shh + 252);
    const auto *shh_254 = buffer.data(shh + 254);
    const auto *shh_255 = buffer.data(shh + 255);
    const auto *shh_257 = buffer.data(shh + 257);
    const auto *shh_258 = buffer.data(shh + 258);
    const auto *shh_261 = buffer.data(shh + 261);
    const auto *shh_267 = buffer.data(shh + 267);
    const auto *shh_272 = buffer.data(shh + 272);
    const auto *shh_273 = buffer.data(shh + 273);
    const auto *shh_275 = buffer.data(shh + 275);
    const auto *shh_276 = buffer.data(shh + 276);
    const auto *shh_278 = buffer.data(shh + 278);
    const auto *shh_279 = buffer.data(shh + 279);
    const auto *shh_282 = buffer.data(shh + 282);
    const auto *shh_288 = buffer.data(shh + 288);
    const auto *shh_293 = buffer.data(shh + 293);
    const auto *shh_294 = buffer.data(shh + 294);
    const auto *shh_296 = buffer.data(shh + 296);
    const auto *shh_297 = buffer.data(shh + 297);
    const auto *shh_299 = buffer.data(shh + 299);
    const auto *shh_300 = buffer.data(shh + 300);
    const auto *shh_303 = buffer.data(shh + 303);
    const auto *shh_309 = buffer.data(shh + 309);
    const auto *shh_314 = buffer.data(shh + 314);
    const auto *shh_315 = buffer.data(shh + 315);
    const auto *shh_357 = buffer.data(shh + 357);
    const auto *shh_360 = buffer.data(shh + 360);
    const auto *shh_362 = buffer.data(shh + 362);
    const auto *shh_363 = buffer.data(shh + 363);
    const auto *shh_366 = buffer.data(shh + 366);
    const auto *shh_367 = buffer.data(shh + 367);
    const auto *shh_369 = buffer.data(shh + 369);
    const auto *shh_371 = buffer.data(shh + 371);
    const auto *shh_372 = buffer.data(shh + 372);
    const auto *shh_373 = buffer.data(shh + 373);
    const auto *shh_374 = buffer.data(shh + 374);
    const auto *shh_375 = buffer.data(shh + 375);
    const auto *shh_376 = buffer.data(shh + 376);
    const auto *shh_377 = buffer.data(shh + 377);
    const auto *shh_378 = buffer.data(shh + 378);
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
    const auto *shh_402 = buffer.data(shh + 402);
    const auto *shh_405 = buffer.data(shh + 405);
    const auto *shh_409 = buffer.data(shh + 409);
    const auto *shh_411 = buffer.data(shh + 411);
    const auto *shh_414 = buffer.data(shh + 414);
    const auto *shh_415 = buffer.data(shh + 415);
    const auto *shh_416 = buffer.data(shh + 416);
    const auto *shh_417 = buffer.data(shh + 417);
    const auto *shh_418 = buffer.data(shh + 418);
    const auto *shh_419 = buffer.data(shh + 419);
    const auto *shh_420 = buffer.data(shh + 420);
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

    const auto *shi1_392 = buffer.data(shi1 + 392);
    const auto *shi1_397 = buffer.data(shi1 + 397);
    const auto *shi1_401 = buffer.data(shi1 + 401);
    const auto *shi1_406 = buffer.data(shi1 + 406);
    const auto *shi1_469 = buffer.data(shi1 + 469);
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
    const auto *shi1_497 = buffer.data(shi1 + 497);
    const auto *shi1_499 = buffer.data(shi1 + 499);
    const auto *shi1_500 = buffer.data(shi1 + 500);
    const auto *shi1_501 = buffer.data(shi1 + 501);
    const auto *shi1_503 = buffer.data(shi1 + 503);
    const auto *shi1_504 = buffer.data(shi1 + 504);
    const auto *shi1_507 = buffer.data(shi1 + 507);
    const auto *shi1_509 = buffer.data(shi1 + 509);
    const auto *shi1_510 = buffer.data(shi1 + 510);
    const auto *shi1_513 = buffer.data(shi1 + 513);
    const auto *shi1_514 = buffer.data(shi1 + 514);
    const auto *shi1_516 = buffer.data(shi1 + 516);
    const auto *shi1_518 = buffer.data(shi1 + 518);
    const auto *shi1_525 = buffer.data(shi1 + 525);
    const auto *shi1_527 = buffer.data(shi1 + 527);
    const auto *shi1_528 = buffer.data(shi1 + 528);
    const auto *shi1_529 = buffer.data(shi1 + 529);
    const auto *shi1_531 = buffer.data(shi1 + 531);
    const auto *shi1_535 = buffer.data(shi1 + 535);
    const auto *shi1_538 = buffer.data(shi1 + 538);
    const auto *shi1_542 = buffer.data(shi1 + 542);
    const auto *shi1_544 = buffer.data(shi1 + 544);
    const auto *shi1_553 = buffer.data(shi1 + 553);
    const auto *shi1_555 = buffer.data(shi1 + 555);
    const auto *shi1_556 = buffer.data(shi1 + 556);
    const auto *shi1_557 = buffer.data(shi1 + 557);
    const auto *shi1_559 = buffer.data(shi1 + 559);
    const auto *shi1_560 = buffer.data(shi1 + 560);
    const auto *shi1_563 = buffer.data(shi1 + 563);
    const auto *shi1_565 = buffer.data(shi1 + 565);
    const auto *shi1_566 = buffer.data(shi1 + 566);
    const auto *shi1_569 = buffer.data(shi1 + 569);
    const auto *shi1_570 = buffer.data(shi1 + 570);
    const auto *shi1_572 = buffer.data(shi1 + 572);
    const auto *shi1_574 = buffer.data(shi1 + 574);
    const auto *shi1_581 = buffer.data(shi1 + 581);
    const auto *shi1_583 = buffer.data(shi1 + 583);
    const auto *shi1_584 = buffer.data(shi1 + 584);
    const auto *shi1_585 = buffer.data(shi1 + 585);
    const auto *shi1_587 = buffer.data(shi1 + 587);

    const auto *sig0_315 = buffer.data(sig0 + 315);

    const auto *sig1_315 = buffer.data(sig1 + 315);

    const auto *sih_351 = buffer.data(sih + 351);
    const auto *sih_356 = buffer.data(sih + 356);
    const auto *sih_357 = buffer.data(sih + 357);
    const auto *sih_359 = buffer.data(sih + 359);
    const auto *sih_360 = buffer.data(sih + 360);
    const auto *sih_362 = buffer.data(sih + 362);
    const auto *sih_363 = buffer.data(sih + 363);
    const auto *sih_366 = buffer.data(sih + 366);
    const auto *sih_372 = buffer.data(sih + 372);
    const auto *sih_373 = buffer.data(sih + 373);
    const auto *sih_374 = buffer.data(sih + 374);
    const auto *sih_375 = buffer.data(sih + 375);
    const auto *sih_376 = buffer.data(sih + 376);
    const auto *sih_377 = buffer.data(sih + 377);
    const auto *sih_378 = buffer.data(sih + 378);
    const auto *sih_380 = buffer.data(sih + 380);
    const auto *sih_381 = buffer.data(sih + 381);
    const auto *sih_383 = buffer.data(sih + 383);
    const auto *sih_384 = buffer.data(sih + 384);
    const auto *sih_387 = buffer.data(sih + 387);
    const auto *sih_393 = buffer.data(sih + 393);
    const auto *sih_394 = buffer.data(sih + 394);
    const auto *sih_395 = buffer.data(sih + 395);
    const auto *sih_396 = buffer.data(sih + 396);
    const auto *sih_397 = buffer.data(sih + 397);
    const auto *sih_398 = buffer.data(sih + 398);
    const auto *sih_399 = buffer.data(sih + 399);
    const auto *sih_401 = buffer.data(sih + 401);
    const auto *sih_402 = buffer.data(sih + 402);
    const auto *sih_404 = buffer.data(sih + 404);
    const auto *sih_405 = buffer.data(sih + 405);
    const auto *sih_408 = buffer.data(sih + 408);
    const auto *sih_414 = buffer.data(sih + 414);
    const auto *sih_415 = buffer.data(sih + 415);
    const auto *sih_416 = buffer.data(sih + 416);
    const auto *sih_417 = buffer.data(sih + 417);
    const auto *sih_418 = buffer.data(sih + 418);
    const auto *sih_419 = buffer.data(sih + 419);
    const auto *sih_420 = buffer.data(sih + 420);
    const auto *sih_422 = buffer.data(sih + 422);
    const auto *sih_423 = buffer.data(sih + 423);
    const auto *sih_425 = buffer.data(sih + 425);
    const auto *sih_426 = buffer.data(sih + 426);
    const auto *sih_429 = buffer.data(sih + 429);
    const auto *sih_435 = buffer.data(sih + 435);
    const auto *sih_436 = buffer.data(sih + 436);
    const auto *sih_437 = buffer.data(sih + 437);
    const auto *sih_438 = buffer.data(sih + 438);
    const auto *sih_439 = buffer.data(sih + 439);
    const auto *sih_440 = buffer.data(sih + 440);
    const auto *sih_441 = buffer.data(sih + 441);

#pragma omp simd aligned(t_469, t_470, t_471, t_472, pb_x, pc_x, pc_z, shi0_469, shi0_471, \
                         shi0_472, shh_225, shi1_469, shi1_471, shi1_472, \
                         sih_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = pb_x[k] * shi0_469[k]
                   - f_10 * pc_x[k] * shi1_469[k];

        t_470[k] = f_11 * shh_225[k]
                   + f_3 * pc_z[k] * sih_351[k];

        t_471[k] = pb_x[k] * shi0_471[k]
                   - f_10 * pc_x[k] * shi1_471[k];

        t_472[k] = pb_x[k] * shi0_472[k]
                   - f_10 * pc_x[k] * shi1_472[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, t_476, pb_x, pc_x, pc_y, shi0_473, shi0_475, \
                         shi0_476, shh_251, shh_357, shi1_473, shi1_475, shi1_476, \
                         sih_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = pb_x[k] * shi0_473[k]
                   - f_10 * pc_x[k] * shi1_473[k];

        t_474[k] = f_14 * shh_251[k]
                   + f_3 * pc_y[k] * sih_356[k];

        t_475[k] = pb_x[k] * shi0_475[k]
                   - f_10 * pc_x[k] * shi1_475[k];

        t_476[k] = pb_x[k] * shi0_476[k]
                   + f_0 * shh_357[k]
                   - f_10 * pc_x[k] * shi1_476[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, pb_x, pc_x, pc_y, pc_z, shi0_479, \
                         shh_231, shh_252, shh_254, shh_360, shi1_479, sih_357, \
                         sih_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_13 * shh_252[k]
                   + f_3 * pc_y[k] * sih_357[k];

        t_478[k] = f_12 * shh_231[k]
                   + f_3 * pc_z[k] * sih_357[k];

        t_479[k] = pb_x[k] * shi0_479[k]
                   + f_14 * shh_360[k]
                   - f_10 * pc_x[k] * shi1_479[k];

        t_480[k] = f_13 * shh_254[k]
                   + f_3 * pc_y[k] * sih_359[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, pb_x, pc_x, pc_z, shi0_481, shi0_482, shh_234, \
                         shh_362, shh_363, shi1_481, shi1_482, \
                         sih_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = pb_x[k] * shi0_481[k]
                   + f_14 * shh_362[k]
                   - f_10 * pc_x[k] * shi1_481[k];

        t_482[k] = pb_x[k] * shi0_482[k]
                   + f_13 * shh_363[k]
                   - f_10 * pc_x[k] * shi1_482[k];

        t_483[k] = f_12 * shh_234[k]
                   + f_3 * pc_z[k] * sih_360[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, pb_x, pc_x, pc_y, shi0_485, shi0_486, shh_257, \
                         shh_366, shh_367, shi1_485, shi1_486, \
                         sih_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_13 * shh_257[k]
                   + f_3 * pc_y[k] * sih_362[k];

        t_485[k] = pb_x[k] * shi0_485[k]
                   + f_13 * shh_366[k]
                   - f_10 * pc_x[k] * shi1_485[k];

        t_486[k] = pb_x[k] * shi0_486[k]
                   + f_12 * shh_367[k]
                   - f_10 * pc_x[k] * shi1_486[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, pb_x, pc_x, pc_y, pc_z, shi0_488, shh_237, \
                         shh_261, shh_369, shi1_488, sih_363, sih_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_12 * shh_237[k]
                   + f_3 * pc_z[k] * sih_363[k];

        t_488[k] = pb_x[k] * shi0_488[k]
                   + f_12 * shh_369[k]
                   - f_10 * pc_x[k] * shi1_488[k];

        t_489[k] = f_13 * shh_261[k]
                   + f_3 * pc_y[k] * sih_366[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, pb_x, pc_x, shi0_490, shh_371, shh_372, \
                         shh_373, shh_374, shi1_490, sih_372, sih_373, \
                         sih_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = pb_x[k] * shi0_490[k]
                   + f_12 * shh_371[k]
                   - f_10 * pc_x[k] * shi1_490[k];

        t_491[k] = f_11 * shh_372[k]
                   + f_3 * pc_x[k] * sih_372[k];

        t_492[k] = f_11 * shh_373[k]
                   + f_3 * pc_x[k] * sih_373[k];

        t_493[k] = f_11 * shh_374[k]
                   + f_3 * pc_x[k] * sih_374[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pb_x, pc_x, shi0_497, shh_375, shh_376, \
                         shh_377, shi1_497, sih_375, sih_376, sih_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_11 * shh_375[k]
                   + f_3 * pc_x[k] * sih_375[k];

        t_495[k] = f_11 * shh_376[k]
                   + f_3 * pc_x[k] * sih_376[k];

        t_496[k] = f_11 * shh_377[k]
                   + f_3 * pc_x[k] * sih_377[k];

        t_497[k] = pb_x[k] * shi0_497[k]
                   - f_10 * pc_x[k] * shi1_497[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, t_501, pb_x, pc_x, pc_z, shi0_499, shi0_500, \
                         shi0_501, shh_246, shi1_499, shi1_500, shi1_501, \
                         sih_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_12 * shh_246[k]
                   + f_3 * pc_z[k] * sih_372[k];

        t_499[k] = pb_x[k] * shi0_499[k]
                   - f_10 * pc_x[k] * shi1_499[k];

        t_500[k] = pb_x[k] * shi0_500[k]
                   - f_10 * pc_x[k] * shi1_500[k];

        t_501[k] = pb_x[k] * shi0_501[k]
                   - f_10 * pc_x[k] * shi1_501[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, pb_x, pc_x, pc_y, shi0_503, shi0_504, \
                         shh_272, shh_273, shh_378, shi1_503, shi1_504, sih_377, \
                         sih_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * shh_272[k]
                   + f_3 * pc_y[k] * sih_377[k];

        t_503[k] = pb_x[k] * shi0_503[k]
                   - f_10 * pc_x[k] * shi1_503[k];

        t_504[k] = pb_x[k] * shi0_504[k]
                   + f_0 * shh_378[k]
                   - f_10 * pc_x[k] * shi1_504[k];

        t_505[k] = f_12 * shh_273[k]
                   + f_3 * pc_y[k] * sih_378[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pb_x, pc_x, pc_y, pc_z, shi0_507, shh_252, \
                         shh_275, shh_381, shi1_507, sih_378, sih_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_13 * shh_252[k]
                   + f_3 * pc_z[k] * sih_378[k];

        t_507[k] = pb_x[k] * shi0_507[k]
                   + f_14 * shh_381[k]
                   - f_10 * pc_x[k] * shi1_507[k];

        t_508[k] = f_12 * shh_275[k]
                   + f_3 * pc_y[k] * sih_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pb_x, pc_x, pc_z, shi0_509, shi0_510, shh_255, \
                         shh_383, shh_384, shi1_509, shi1_510, \
                         sih_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = pb_x[k] * shi0_509[k]
                   + f_14 * shh_383[k]
                   - f_10 * pc_x[k] * shi1_509[k];

        t_510[k] = pb_x[k] * shi0_510[k]
                   + f_13 * shh_384[k]
                   - f_10 * pc_x[k] * shi1_510[k];

        t_511[k] = f_13 * shh_255[k]
                   + f_3 * pc_z[k] * sih_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pb_x, pc_x, pc_y, shi0_513, shi0_514, shh_278, \
                         shh_387, shh_388, shi1_513, shi1_514, \
                         sih_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * shh_278[k]
                   + f_3 * pc_y[k] * sih_383[k];

        t_513[k] = pb_x[k] * shi0_513[k]
                   + f_13 * shh_387[k]
                   - f_10 * pc_x[k] * shi1_513[k];

        t_514[k] = pb_x[k] * shi0_514[k]
                   + f_12 * shh_388[k]
                   - f_10 * pc_x[k] * shi1_514[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pb_x, pc_x, pc_y, pc_z, shi0_516, shh_258, \
                         shh_282, shh_390, shi1_516, sih_384, sih_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * shh_258[k]
                   + f_3 * pc_z[k] * sih_384[k];

        t_516[k] = pb_x[k] * shi0_516[k]
                   + f_12 * shh_390[k]
                   - f_10 * pc_x[k] * shi1_516[k];

        t_517[k] = f_12 * shh_282[k]
                   + f_3 * pc_y[k] * sih_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pb_x, pc_x, shi0_518, shh_392, shh_393, \
                         shh_394, shh_395, shi1_518, sih_393, sih_394, \
                         sih_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = pb_x[k] * shi0_518[k]
                   + f_12 * shh_392[k]
                   - f_10 * pc_x[k] * shi1_518[k];

        t_519[k] = f_11 * shh_393[k]
                   + f_3 * pc_x[k] * sih_393[k];

        t_520[k] = f_11 * shh_394[k]
                   + f_3 * pc_x[k] * sih_394[k];

        t_521[k] = f_11 * shh_395[k]
                   + f_3 * pc_x[k] * sih_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pb_x, pc_x, shi0_525, shh_396, shh_397, \
                         shh_398, shi1_525, sih_396, sih_397, sih_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_11 * shh_396[k]
                   + f_3 * pc_x[k] * sih_396[k];

        t_523[k] = f_11 * shh_397[k]
                   + f_3 * pc_x[k] * sih_397[k];

        t_524[k] = f_11 * shh_398[k]
                   + f_3 * pc_x[k] * sih_398[k];

        t_525[k] = pb_x[k] * shi0_525[k]
                   - f_10 * pc_x[k] * shi1_525[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pb_x, pc_x, pc_z, shi0_527, shi0_528, \
                         shi0_529, shh_267, shi1_527, shi1_528, shi1_529, \
                         sih_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_13 * shh_267[k]
                   + f_3 * pc_z[k] * sih_393[k];

        t_527[k] = pb_x[k] * shi0_527[k]
                   - f_10 * pc_x[k] * shi1_527[k];

        t_528[k] = pb_x[k] * shi0_528[k]
                   - f_10 * pc_x[k] * shi1_528[k];

        t_529[k] = pb_x[k] * shi0_529[k]
                   - f_10 * pc_x[k] * shi1_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pb_x, pb_y, pc_x, pc_y, shi0_392, \
                         shi0_531, shh_293, shh_294, shi1_392, shi1_531, sih_398, \
                         sih_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_12 * shh_293[k]
                   + f_3 * pc_y[k] * sih_398[k];

        t_531[k] = pb_x[k] * shi0_531[k]
                   - f_10 * pc_x[k] * shi1_531[k];

        t_532[k] = pb_y[k] * shi0_392[k]
                   - f_10 * pc_y[k] * shi1_392[k];

        t_533[k] = f_11 * shh_294[k]
                   + f_3 * pc_y[k] * sih_399[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pb_x, pc_x, pc_y, pc_z, shi0_535, shh_273, \
                         shh_296, shh_402, shi1_535, sih_399, sih_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_14 * shh_273[k]
                   + f_3 * pc_z[k] * sih_399[k];

        t_535[k] = pb_x[k] * shi0_535[k]
                   + f_14 * shh_402[k]
                   - f_10 * pc_x[k] * shi1_535[k];

        t_536[k] = f_11 * shh_296[k]
                   + f_3 * pc_y[k] * sih_401[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, pb_x, pb_y, pc_x, pc_y, pc_z, shi0_397, \
                         shi0_538, shh_276, shh_405, shi1_397, shi1_538, \
                         sih_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pb_y[k] * shi0_397[k]
                   - f_10 * pc_y[k] * shi1_397[k];

        t_538[k] = pb_x[k] * shi0_538[k]
                   + f_13 * shh_405[k]
                   - f_10 * pc_x[k] * shi1_538[k];

        t_539[k] = f_14 * shh_276[k]
                   + f_3 * pc_z[k] * sih_402[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, pb_x, pb_y, pc_x, pc_y, shi0_401, shi0_542, \
                         shh_299, shh_409, shi1_401, shi1_542, \
                         sih_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_11 * shh_299[k]
                   + f_3 * pc_y[k] * sih_404[k];

        t_541[k] = pb_y[k] * shi0_401[k]
                   - f_10 * pc_y[k] * shi1_401[k];

        t_542[k] = pb_x[k] * shi0_542[k]
                   + f_12 * shh_409[k]
                   - f_10 * pc_x[k] * shi1_542[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, pb_x, pc_x, pc_y, pc_z, shi0_544, shh_279, \
                         shh_303, shh_411, shi1_544, sih_405, sih_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_14 * shh_279[k]
                   + f_3 * pc_z[k] * sih_405[k];

        t_544[k] = pb_x[k] * shi0_544[k]
                   + f_12 * shh_411[k]
                   - f_10 * pc_x[k] * shi1_544[k];

        t_545[k] = f_11 * shh_303[k]
                   + f_3 * pc_y[k] * sih_408[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, pb_y, pc_x, pc_y, shi0_406, shh_414, \
                         shh_415, shh_416, shi1_406, sih_414, sih_415, \
                         sih_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = pb_y[k] * shi0_406[k]
                   - f_10 * pc_y[k] * shi1_406[k];

        t_547[k] = f_11 * shh_414[k]
                   + f_3 * pc_x[k] * sih_414[k];

        t_548[k] = f_11 * shh_415[k]
                   + f_3 * pc_x[k] * sih_415[k];

        t_549[k] = f_11 * shh_416[k]
                   + f_3 * pc_x[k] * sih_416[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, pb_x, pc_x, shi0_553, shh_417, shh_418, \
                         shh_419, shi1_553, sih_417, sih_418, sih_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_11 * shh_417[k]
                   + f_3 * pc_x[k] * sih_417[k];

        t_551[k] = f_11 * shh_418[k]
                   + f_3 * pc_x[k] * sih_418[k];

        t_552[k] = f_11 * shh_419[k]
                   + f_3 * pc_x[k] * sih_419[k];

        t_553[k] = pb_x[k] * shi0_553[k]
                   - f_10 * pc_x[k] * shi1_553[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, t_557, pb_x, pc_x, pc_z, shi0_555, shi0_556, \
                         shi0_557, shh_288, shi1_555, shi1_556, shi1_557, \
                         sih_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_14 * shh_288[k]
                   + f_3 * pc_z[k] * sih_414[k];

        t_555[k] = pb_x[k] * shi0_555[k]
                   - f_10 * pc_x[k] * shi1_555[k];

        t_556[k] = pb_x[k] * shi0_556[k]
                   - f_10 * pc_x[k] * shi1_556[k];

        t_557[k] = pb_x[k] * shi0_557[k]
                   - f_10 * pc_x[k] * shi1_557[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, pb_x, pc_x, pc_y, shi0_559, shi0_560, \
                         shh_314, shh_420, shi1_559, shi1_560, sih_419, \
                         sih_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_11 * shh_314[k]
                   + f_3 * pc_y[k] * sih_419[k];

        t_559[k] = pb_x[k] * shi0_559[k]
                   - f_10 * pc_x[k] * shi1_559[k];

        t_560[k] = pb_x[k] * shi0_560[k]
                   + f_0 * shh_420[k]
                   - f_10 * pc_x[k] * shi1_560[k];

        t_561[k] = f_3 * pc_y[k] * sih_420[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, pb_x, pc_x, pc_y, pc_z, shi0_563, shh_294, \
                         shh_423, shi1_563, sih_420, sih_422 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_15 * shh_294[k]
                   + f_3 * pc_z[k] * sih_420[k];

        t_563[k] = pb_x[k] * shi0_563[k]
                   + f_14 * shh_423[k]
                   - f_10 * pc_x[k] * shi1_563[k];

        t_564[k] = f_3 * pc_y[k] * sih_422[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, pb_x, pc_x, pc_z, shi0_565, shi0_566, shh_297, \
                         shh_425, shh_426, shi1_565, shi1_566, \
                         sih_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = pb_x[k] * shi0_565[k]
                   + f_14 * shh_425[k]
                   - f_10 * pc_x[k] * shi1_565[k];

        t_566[k] = pb_x[k] * shi0_566[k]
                   + f_13 * shh_426[k]
                   - f_10 * pc_x[k] * shi1_566[k];

        t_567[k] = f_15 * shh_297[k]
                   + f_3 * pc_z[k] * sih_423[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, pb_x, pc_x, pc_y, shi0_569, shi0_570, shh_429, \
                         shh_430, shi1_569, shi1_570, sih_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_3 * pc_y[k] * sih_425[k];

        t_569[k] = pb_x[k] * shi0_569[k]
                   + f_13 * shh_429[k]
                   - f_10 * pc_x[k] * shi1_569[k];

        t_570[k] = pb_x[k] * shi0_570[k]
                   + f_12 * shh_430[k]
                   - f_10 * pc_x[k] * shi1_570[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pb_x, pc_x, pc_y, pc_z, shi0_572, shh_300, \
                         shh_432, shi1_572, sih_426, sih_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_15 * shh_300[k]
                   + f_3 * pc_z[k] * sih_426[k];

        t_572[k] = pb_x[k] * shi0_572[k]
                   + f_12 * shh_432[k]
                   - f_10 * pc_x[k] * shi1_572[k];

        t_573[k] = f_3 * pc_y[k] * sih_429[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pb_x, pc_x, shi0_574, shh_434, shh_435, \
                         shh_436, shh_437, shi1_574, sih_435, sih_436, \
                         sih_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = pb_x[k] * shi0_574[k]
                   + f_12 * shh_434[k]
                   - f_10 * pc_x[k] * shi1_574[k];

        t_575[k] = f_11 * shh_435[k]
                   + f_3 * pc_x[k] * sih_435[k];

        t_576[k] = f_11 * shh_436[k]
                   + f_3 * pc_x[k] * sih_436[k];

        t_577[k] = f_11 * shh_437[k]
                   + f_3 * pc_x[k] * sih_437[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, pb_x, pc_x, shi0_581, shh_438, shh_439, \
                         shh_440, shi1_581, sih_438, sih_439, sih_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_11 * shh_438[k]
                   + f_3 * pc_x[k] * sih_438[k];

        t_579[k] = f_11 * shh_439[k]
                   + f_3 * pc_x[k] * sih_439[k];

        t_580[k] = f_11 * shh_440[k]
                   + f_3 * pc_x[k] * sih_440[k];

        t_581[k] = pb_x[k] * shi0_581[k]
                   - f_10 * pc_x[k] * shi1_581[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pb_x, pc_x, pc_z, shi0_583, shi0_584, \
                         shi0_585, shh_309, shi1_583, shi1_584, shi1_585, \
                         sih_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_15 * shh_309[k]
                   + f_3 * pc_z[k] * sih_435[k];

        t_583[k] = pb_x[k] * shi0_583[k]
                   - f_10 * pc_x[k] * shi1_583[k];

        t_584[k] = pb_x[k] * shi0_584[k]
                   - f_10 * pc_x[k] * shi1_584[k];

        t_585[k] = pb_x[k] * shi0_585[k]
                   - f_10 * pc_x[k] * shi1_585[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, t_590, pb_x, pc_x, pc_y, pc_z, shi0_587, \
                         shh_315, shi1_587, sig0_315, sig1_315, sih_440, \
                         sih_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_3 * pc_y[k] * sih_440[k];

        t_587[k] = pb_x[k] * shi0_587[k]
                   - f_10 * pc_x[k] * shi1_587[k];

        t_588[k] = f_1 * sig0_315[k]
                   - f_2 * sig1_315[k]
                   + f_3 * pc_x[k] * sih_441[k];

        t_589[k] = f_0 * shh_315[k]
                   + f_3 * pc_y[k] * sih_441[k];

        t_590[k] = f_3 * pc_z[k] * sih_441[k];
    }
}

static auto
compute_prim_sii_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shi0,
                                                          const size_t shh, const size_t shi1,
                                                          const size_t sig0, const size_t sig1,
                                                          const size_t sih, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_15 = 2.5 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shi0_420 = buffer.data(shi0 + 420);
    const auto *shi0_423 = buffer.data(shi0 + 423);
    const auto *shi0_426 = buffer.data(shi0 + 426);
    const auto *shi0_430 = buffer.data(shi0 + 430);
    const auto *shi0_441 = buffer.data(shi0 + 441);
    const auto *shi0_443 = buffer.data(shi0 + 443);
    const auto *shi0_444 = buffer.data(shi0 + 444);
    const auto *shi0_445 = buffer.data(shi0 + 445);

    const auto *shh_315 = buffer.data(shh + 315);
    const auto *shh_317 = buffer.data(shh + 317);
    const auto *shh_318 = buffer.data(shh + 318);
    const auto *shh_320 = buffer.data(shh + 320);
    const auto *shh_321 = buffer.data(shh + 321);
    const auto *shh_324 = buffer.data(shh + 324);
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
    const auto *shh_351 = buffer.data(shh + 351);
    const auto *shh_356 = buffer.data(shh + 356);
    const auto *shh_357 = buffer.data(shh + 357);
    const auto *shh_359 = buffer.data(shh + 359);
    const auto *shh_360 = buffer.data(shh + 360);
    const auto *shh_362 = buffer.data(shh + 362);
    const auto *shh_363 = buffer.data(shh + 363);
    const auto *shh_366 = buffer.data(shh + 366);
    const auto *shh_372 = buffer.data(shh + 372);
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
    const auto *shh_393 = buffer.data(shh + 393);
    const auto *shh_395 = buffer.data(shh + 395);
    const auto *shh_396 = buffer.data(shh + 396);
    const auto *shh_397 = buffer.data(shh + 397);
    const auto *shh_398 = buffer.data(shh + 398);
    const auto *shh_399 = buffer.data(shh + 399);
    const auto *shh_401 = buffer.data(shh + 401);
    const auto *shh_404 = buffer.data(shh + 404);
    const auto *shh_408 = buffer.data(shh + 408);

    const auto *shi1_420 = buffer.data(shi1 + 420);
    const auto *shi1_423 = buffer.data(shi1 + 423);
    const auto *shi1_426 = buffer.data(shi1 + 426);
    const auto *shi1_430 = buffer.data(shi1 + 430);
    const auto *shi1_441 = buffer.data(shi1 + 441);
    const auto *shi1_443 = buffer.data(shi1 + 443);
    const auto *shi1_444 = buffer.data(shi1 + 444);
    const auto *shi1_445 = buffer.data(shi1 + 445);

    const auto *sig0_318 = buffer.data(sig0 + 318);
    const auto *sig0_320 = buffer.data(sig0 + 320);
    const auto *sig0_321 = buffer.data(sig0 + 321);
    const auto *sig0_324 = buffer.data(sig0 + 324);
    const auto *sig0_325 = buffer.data(sig0 + 325);
    const auto *sig0_327 = buffer.data(sig0 + 327);
    const auto *sig0_328 = buffer.data(sig0 + 328);
    const auto *sig0_329 = buffer.data(sig0 + 329);
    const auto *sig0_335 = buffer.data(sig0 + 335);
    const auto *sig0_339 = buffer.data(sig0 + 339);
    const auto *sig0_342 = buffer.data(sig0 + 342);
    const auto *sig0_344 = buffer.data(sig0 + 344);
    const auto *sig0_345 = buffer.data(sig0 + 345);
    const auto *sig0_348 = buffer.data(sig0 + 348);
    const auto *sig0_350 = buffer.data(sig0 + 350);
    const auto *sig0_351 = buffer.data(sig0 + 351);
    const auto *sig0_354 = buffer.data(sig0 + 354);
    const auto *sig0_355 = buffer.data(sig0 + 355);
    const auto *sig0_357 = buffer.data(sig0 + 357);
    const auto *sig0_358 = buffer.data(sig0 + 358);
    const auto *sig0_359 = buffer.data(sig0 + 359);
    const auto *sig0_360 = buffer.data(sig0 + 360);
    const auto *sig0_363 = buffer.data(sig0 + 363);
    const auto *sig0_365 = buffer.data(sig0 + 365);
    const auto *sig0_366 = buffer.data(sig0 + 366);
    const auto *sig0_369 = buffer.data(sig0 + 369);
    const auto *sig0_370 = buffer.data(sig0 + 370);
    const auto *sig0_372 = buffer.data(sig0 + 372);
    const auto *sig0_373 = buffer.data(sig0 + 373);
    const auto *sig0_374 = buffer.data(sig0 + 374);
    const auto *sig0_375 = buffer.data(sig0 + 375);
    const auto *sig0_378 = buffer.data(sig0 + 378);
    const auto *sig0_380 = buffer.data(sig0 + 380);
    const auto *sig0_381 = buffer.data(sig0 + 381);
    const auto *sig0_384 = buffer.data(sig0 + 384);
    const auto *sig0_385 = buffer.data(sig0 + 385);
    const auto *sig0_387 = buffer.data(sig0 + 387);
    const auto *sig0_389 = buffer.data(sig0 + 389);

    const auto *sig1_318 = buffer.data(sig1 + 318);
    const auto *sig1_320 = buffer.data(sig1 + 320);
    const auto *sig1_321 = buffer.data(sig1 + 321);
    const auto *sig1_324 = buffer.data(sig1 + 324);
    const auto *sig1_325 = buffer.data(sig1 + 325);
    const auto *sig1_327 = buffer.data(sig1 + 327);
    const auto *sig1_328 = buffer.data(sig1 + 328);
    const auto *sig1_329 = buffer.data(sig1 + 329);
    const auto *sig1_335 = buffer.data(sig1 + 335);
    const auto *sig1_339 = buffer.data(sig1 + 339);
    const auto *sig1_342 = buffer.data(sig1 + 342);
    const auto *sig1_344 = buffer.data(sig1 + 344);
    const auto *sig1_345 = buffer.data(sig1 + 345);
    const auto *sig1_348 = buffer.data(sig1 + 348);
    const auto *sig1_350 = buffer.data(sig1 + 350);
    const auto *sig1_351 = buffer.data(sig1 + 351);
    const auto *sig1_354 = buffer.data(sig1 + 354);
    const auto *sig1_355 = buffer.data(sig1 + 355);
    const auto *sig1_357 = buffer.data(sig1 + 357);
    const auto *sig1_358 = buffer.data(sig1 + 358);
    const auto *sig1_359 = buffer.data(sig1 + 359);
    const auto *sig1_360 = buffer.data(sig1 + 360);
    const auto *sig1_363 = buffer.data(sig1 + 363);
    const auto *sig1_365 = buffer.data(sig1 + 365);
    const auto *sig1_366 = buffer.data(sig1 + 366);
    const auto *sig1_369 = buffer.data(sig1 + 369);
    const auto *sig1_370 = buffer.data(sig1 + 370);
    const auto *sig1_372 = buffer.data(sig1 + 372);
    const auto *sig1_373 = buffer.data(sig1 + 373);
    const auto *sig1_374 = buffer.data(sig1 + 374);
    const auto *sig1_375 = buffer.data(sig1 + 375);
    const auto *sig1_378 = buffer.data(sig1 + 378);
    const auto *sig1_380 = buffer.data(sig1 + 380);
    const auto *sig1_381 = buffer.data(sig1 + 381);
    const auto *sig1_384 = buffer.data(sig1 + 384);
    const auto *sig1_385 = buffer.data(sig1 + 385);
    const auto *sig1_387 = buffer.data(sig1 + 387);
    const auto *sig1_389 = buffer.data(sig1 + 389);

    const auto *sih_443 = buffer.data(sih + 443);
    const auto *sih_444 = buffer.data(sih + 444);
    const auto *sih_446 = buffer.data(sih + 446);
    const auto *sih_447 = buffer.data(sih + 447);
    const auto *sih_450 = buffer.data(sih + 450);
    const auto *sih_451 = buffer.data(sih + 451);
    const auto *sih_453 = buffer.data(sih + 453);
    const auto *sih_455 = buffer.data(sih + 455);
    const auto *sih_456 = buffer.data(sih + 456);
    const auto *sih_457 = buffer.data(sih + 457);
    const auto *sih_458 = buffer.data(sih + 458);
    const auto *sih_459 = buffer.data(sih + 459);
    const auto *sih_460 = buffer.data(sih + 460);
    const auto *sih_461 = buffer.data(sih + 461);
    const auto *sih_462 = buffer.data(sih + 462);
    const auto *sih_464 = buffer.data(sih + 464);
    const auto *sih_465 = buffer.data(sih + 465);
    const auto *sih_467 = buffer.data(sih + 467);
    const auto *sih_468 = buffer.data(sih + 468);
    const auto *sih_471 = buffer.data(sih + 471);
    const auto *sih_474 = buffer.data(sih + 474);
    const auto *sih_476 = buffer.data(sih + 476);
    const auto *sih_477 = buffer.data(sih + 477);
    const auto *sih_478 = buffer.data(sih + 478);
    const auto *sih_479 = buffer.data(sih + 479);
    const auto *sih_480 = buffer.data(sih + 480);
    const auto *sih_481 = buffer.data(sih + 481);
    const auto *sih_482 = buffer.data(sih + 482);
    const auto *sih_483 = buffer.data(sih + 483);
    const auto *sih_485 = buffer.data(sih + 485);
    const auto *sih_486 = buffer.data(sih + 486);
    const auto *sih_488 = buffer.data(sih + 488);
    const auto *sih_489 = buffer.data(sih + 489);
    const auto *sih_492 = buffer.data(sih + 492);
    const auto *sih_493 = buffer.data(sih + 493);
    const auto *sih_495 = buffer.data(sih + 495);
    const auto *sih_497 = buffer.data(sih + 497);
    const auto *sih_498 = buffer.data(sih + 498);
    const auto *sih_499 = buffer.data(sih + 499);
    const auto *sih_500 = buffer.data(sih + 500);
    const auto *sih_501 = buffer.data(sih + 501);
    const auto *sih_502 = buffer.data(sih + 502);
    const auto *sih_503 = buffer.data(sih + 503);
    const auto *sih_504 = buffer.data(sih + 504);
    const auto *sih_506 = buffer.data(sih + 506);
    const auto *sih_507 = buffer.data(sih + 507);
    const auto *sih_509 = buffer.data(sih + 509);
    const auto *sih_510 = buffer.data(sih + 510);
    const auto *sih_513 = buffer.data(sih + 513);
    const auto *sih_514 = buffer.data(sih + 514);
    const auto *sih_516 = buffer.data(sih + 516);
    const auto *sih_518 = buffer.data(sih + 518);
    const auto *sih_519 = buffer.data(sih + 519);
    const auto *sih_520 = buffer.data(sih + 520);
    const auto *sih_521 = buffer.data(sih + 521);
    const auto *sih_522 = buffer.data(sih + 522);
    const auto *sih_523 = buffer.data(sih + 523);
    const auto *sih_524 = buffer.data(sih + 524);
    const auto *sih_525 = buffer.data(sih + 525);
    const auto *sih_527 = buffer.data(sih + 527);
    const auto *sih_528 = buffer.data(sih + 528);
    const auto *sih_530 = buffer.data(sih + 530);
    const auto *sih_531 = buffer.data(sih + 531);
    const auto *sih_534 = buffer.data(sih + 534);
    const auto *sih_535 = buffer.data(sih + 535);
    const auto *sih_537 = buffer.data(sih + 537);
    const auto *sih_539 = buffer.data(sih + 539);
    const auto *sih_540 = buffer.data(sih + 540);

#pragma omp simd aligned(t_591, t_592, t_593, pc_x, pc_y, shh_317, sig0_318, sig0_320, \
                         sig1_318, sig1_320, sih_443, sih_444, \
                         sih_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_4 * sig0_318[k]
                   - f_5 * sig1_318[k]
                   + f_3 * pc_x[k] * sih_444[k];

        t_592[k] = f_0 * shh_317[k]
                   + f_3 * pc_y[k] * sih_443[k];

        t_593[k] = f_4 * sig0_320[k]
                   - f_5 * sig1_320[k]
                   + f_3 * pc_x[k] * sih_446[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, pc_x, pc_y, pc_z, shh_320, sig0_321, \
                         sig0_324, sig1_321, sig1_324, sih_444, sih_446, sih_447, \
                         sih_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_6 * sig0_321[k]
                   - f_7 * sig1_321[k]
                   + f_3 * pc_x[k] * sih_447[k];

        t_595[k] = f_3 * pc_z[k] * sih_444[k];

        t_596[k] = f_0 * shh_320[k]
                   + f_3 * pc_y[k] * sih_446[k];

        t_597[k] = f_6 * sig0_324[k]
                   - f_7 * sig1_324[k]
                   + f_3 * pc_x[k] * sih_450[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, t_601, pc_x, pc_y, pc_z, shh_324, sig0_325, \
                         sig0_327, sig1_325, sig1_327, sih_447, sih_450, sih_451, \
                         sih_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_8 * sig0_325[k]
                   - f_9 * sig1_325[k]
                   + f_3 * pc_x[k] * sih_451[k];

        t_599[k] = f_3 * pc_z[k] * sih_447[k];

        t_600[k] = f_8 * sig0_327[k]
                   - f_9 * sig1_327[k]
                   + f_3 * pc_x[k] * sih_453[k];

        t_601[k] = f_0 * shh_324[k]
                   + f_3 * pc_y[k] * sih_450[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, t_607, pc_x, sig0_329, sig1_329, \
                         sih_455, sih_456, sih_457, sih_458, sih_459, \
                         sih_460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_8 * sig0_329[k]
                   - f_9 * sig1_329[k]
                   + f_3 * pc_x[k] * sih_455[k];

        t_603[k] = f_3 * pc_x[k] * sih_456[k];

        t_604[k] = f_3 * pc_x[k] * sih_457[k];

        t_605[k] = f_3 * pc_x[k] * sih_458[k];

        t_606[k] = f_3 * pc_x[k] * sih_459[k];

        t_607[k] = f_3 * pc_x[k] * sih_460[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, pc_x, pc_y, pc_z, shh_330, shh_332, \
                         sig0_325, sig0_327, sig1_325, sig1_327, sih_456, sih_458, \
                         sih_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_3 * pc_x[k] * sih_461[k];

        t_609[k] = f_0 * shh_330[k]
                   + f_1 * sig0_325[k]
                   - f_2 * sig1_325[k]
                   + f_3 * pc_y[k] * sih_456[k];

        t_610[k] = f_3 * pc_z[k] * sih_456[k];

        t_611[k] = f_0 * shh_332[k]
                   + f_4 * sig0_327[k]
                   - f_5 * sig1_327[k]
                   + f_3 * pc_y[k] * sih_458[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, pc_y, pc_z, shh_333, shh_334, shh_335, \
                         sig0_328, sig0_329, sig1_328, sig1_329, sih_459, sih_460, \
                         sih_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_0 * shh_333[k]
                   + f_6 * sig0_328[k]
                   - f_7 * sig1_328[k]
                   + f_3 * pc_y[k] * sih_459[k];

        t_613[k] = f_0 * shh_334[k]
                   + f_8 * sig0_329[k]
                   - f_9 * sig1_329[k]
                   + f_3 * pc_y[k] * sih_460[k];

        t_614[k] = f_0 * shh_335[k]
                   + f_3 * pc_y[k] * sih_461[k];

        t_615[k] = f_1 * sig0_329[k]
                   - f_2 * sig1_329[k]
                   + f_3 * pc_z[k] * sih_461[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, pb_z, pc_y, pc_z, shi0_420, shi0_423, \
                         shh_315, shh_336, shi1_420, shi1_423, \
                         sih_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = pb_z[k] * shi0_420[k]
                   - f_10 * pc_z[k] * shi1_420[k];

        t_617[k] = f_15 * shh_336[k]
                   + f_3 * pc_y[k] * sih_462[k];

        t_618[k] = f_11 * shh_315[k]
                   + f_3 * pc_z[k] * sih_462[k];

        t_619[k] = pb_z[k] * shi0_423[k]
                   - f_10 * pc_z[k] * shi1_423[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pb_z, pc_x, pc_y, pc_z, shi0_426, shh_338, \
                         shi1_426, sig0_335, sig1_335, sih_464, \
                         sih_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_15 * shh_338[k]
                   + f_3 * pc_y[k] * sih_464[k];

        t_621[k] = f_4 * sig0_335[k]
                   - f_5 * sig1_335[k]
                   + f_3 * pc_x[k] * sih_467[k];

        t_622[k] = pb_z[k] * shi0_426[k]
                   - f_10 * pc_z[k] * shi1_426[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, pc_y, pc_z, shh_318, shh_341, sig0_339, \
                         sig1_339, sih_465, sih_467, sih_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_11 * shh_318[k]
                   + f_3 * pc_z[k] * sih_465[k];

        t_624[k] = f_15 * shh_341[k]
                   + f_3 * pc_y[k] * sih_467[k];

        t_625[k] = f_6 * sig0_339[k]
                   - f_7 * sig1_339[k]
                   + f_3 * pc_x[k] * sih_471[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pb_z, pc_x, pc_z, shi0_430, shh_321, shi1_430, \
                         sig0_342, sig1_342, sih_468, sih_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pb_z[k] * shi0_430[k]
                   - f_10 * pc_z[k] * shi1_430[k];

        t_627[k] = f_11 * shh_321[k]
                   + f_3 * pc_z[k] * sih_468[k];

        t_628[k] = f_8 * sig0_342[k]
                   - f_9 * sig1_342[k]
                   + f_3 * pc_x[k] * sih_474[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, t_633, pc_x, pc_y, shh_345, sig0_344, \
                         sig1_344, sih_471, sih_476, sih_477, sih_478, \
                         sih_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_15 * shh_345[k]
                   + f_3 * pc_y[k] * sih_471[k];

        t_630[k] = f_8 * sig0_344[k]
                   - f_9 * sig1_344[k]
                   + f_3 * pc_x[k] * sih_476[k];

        t_631[k] = f_3 * pc_x[k] * sih_477[k];

        t_632[k] = f_3 * pc_x[k] * sih_478[k];

        t_633[k] = f_3 * pc_x[k] * sih_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, pb_z, pc_x, pc_z, shi0_441, \
                         shh_330, shi1_441, sih_477, sih_480, sih_481, \
                         sih_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_3 * pc_x[k] * sih_480[k];

        t_635[k] = f_3 * pc_x[k] * sih_481[k];

        t_636[k] = f_3 * pc_x[k] * sih_482[k];

        t_637[k] = pb_z[k] * shi0_441[k]
                   - f_10 * pc_z[k] * shi1_441[k];

        t_638[k] = f_11 * shh_330[k]
                   + f_3 * pc_z[k] * sih_477[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pb_z, pc_z, shi0_443, shi0_444, shi0_445, \
                         shh_331, shh_332, shh_333, shi1_443, shi1_444, \
                         shi1_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = pb_z[k] * shi0_443[k]
                   + f_12 * shh_331[k]
                   - f_10 * pc_z[k] * shi1_443[k];

        t_640[k] = pb_z[k] * shi0_444[k]
                   + f_13 * shh_332[k]
                   - f_10 * pc_z[k] * shi1_444[k];

        t_641[k] = pb_z[k] * shi0_445[k]
                   + f_14 * shh_333[k]
                   - f_10 * pc_z[k] * shi1_445[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pc_x, pc_y, pc_z, shh_335, shh_356, \
                         shh_357, sig0_344, sig0_345, sig1_344, sig1_345, sih_482, \
                         sih_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_15 * shh_356[k]
                   + f_3 * pc_y[k] * sih_482[k];

        t_643[k] = f_11 * shh_335[k]
                   + f_1 * sig0_344[k]
                   - f_2 * sig1_344[k]
                   + f_3 * pc_z[k] * sih_482[k];

        t_644[k] = f_1 * sig0_345[k]
                   - f_2 * sig1_345[k]
                   + f_3 * pc_x[k] * sih_483[k];

        t_645[k] = f_14 * shh_357[k]
                   + f_3 * pc_y[k] * sih_483[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_x, pc_y, pc_z, shh_336, shh_359, sig0_348, \
                         sig1_348, sih_483, sih_485, sih_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_12 * shh_336[k]
                   + f_3 * pc_z[k] * sih_483[k];

        t_647[k] = f_4 * sig0_348[k]
                   - f_5 * sig1_348[k]
                   + f_3 * pc_x[k] * sih_486[k];

        t_648[k] = f_14 * shh_359[k]
                   + f_3 * pc_y[k] * sih_485[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pc_x, pc_y, pc_z, shh_339, shh_362, \
                         sig0_350, sig0_351, sig1_350, sig1_351, sih_486, sih_488, \
                         sih_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_4 * sig0_350[k]
                   - f_5 * sig1_350[k]
                   + f_3 * pc_x[k] * sih_488[k];

        t_650[k] = f_6 * sig0_351[k]
                   - f_7 * sig1_351[k]
                   + f_3 * pc_x[k] * sih_489[k];

        t_651[k] = f_12 * shh_339[k]
                   + f_3 * pc_z[k] * sih_486[k];

        t_652[k] = f_14 * shh_362[k]
                   + f_3 * pc_y[k] * sih_488[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, shh_342, sig0_354, sig0_355, \
                         sig1_354, sig1_355, sih_489, sih_492, \
                         sih_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_6 * sig0_354[k]
                   - f_7 * sig1_354[k]
                   + f_3 * pc_x[k] * sih_492[k];

        t_654[k] = f_8 * sig0_355[k]
                   - f_9 * sig1_355[k]
                   + f_3 * pc_x[k] * sih_493[k];

        t_655[k] = f_12 * shh_342[k]
                   + f_3 * pc_z[k] * sih_489[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, t_659, pc_x, pc_y, shh_366, sig0_357, sig0_359, \
                         sig1_357, sig1_359, sih_492, sih_495, sih_497, \
                         sih_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_8 * sig0_357[k]
                   - f_9 * sig1_357[k]
                   + f_3 * pc_x[k] * sih_495[k];

        t_657[k] = f_14 * shh_366[k]
                   + f_3 * pc_y[k] * sih_492[k];

        t_658[k] = f_8 * sig0_359[k]
                   - f_9 * sig1_359[k]
                   + f_3 * pc_x[k] * sih_497[k];

        t_659[k] = f_3 * pc_x[k] * sih_498[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, pc_x, sih_499, sih_500, sih_501, \
                         sih_502, sih_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_3 * pc_x[k] * sih_499[k];

        t_661[k] = f_3 * pc_x[k] * sih_500[k];

        t_662[k] = f_3 * pc_x[k] * sih_501[k];

        t_663[k] = f_3 * pc_x[k] * sih_502[k];

        t_664[k] = f_3 * pc_x[k] * sih_503[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, pc_y, pc_z, shh_351, shh_372, shh_374, sig0_355, \
                         sig0_357, sig1_355, sig1_357, sih_498, \
                         sih_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_14 * shh_372[k]
                   + f_1 * sig0_355[k]
                   - f_2 * sig1_355[k]
                   + f_3 * pc_y[k] * sih_498[k];

        t_666[k] = f_12 * shh_351[k]
                   + f_3 * pc_z[k] * sih_498[k];

        t_667[k] = f_14 * shh_374[k]
                   + f_4 * sig0_357[k]
                   - f_5 * sig1_357[k]
                   + f_3 * pc_y[k] * sih_500[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, pc_y, shh_375, shh_376, shh_377, sig0_358, \
                         sig0_359, sig1_358, sig1_359, sih_501, sih_502, \
                         sih_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_14 * shh_375[k]
                   + f_6 * sig0_358[k]
                   - f_7 * sig1_358[k]
                   + f_3 * pc_y[k] * sih_501[k];

        t_669[k] = f_14 * shh_376[k]
                   + f_8 * sig0_359[k]
                   - f_9 * sig1_359[k]
                   + f_3 * pc_y[k] * sih_502[k];

        t_670[k] = f_14 * shh_377[k]
                   + f_3 * pc_y[k] * sih_503[k];
    }

#pragma omp simd aligned(t_671, t_672, t_673, t_674, pc_x, pc_y, pc_z, shh_356, shh_357, \
                         shh_378, sig0_359, sig0_360, sig1_359, sig1_360, sih_503, \
                         sih_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_671[k] = f_12 * shh_356[k]
                   + f_1 * sig0_359[k]
                   - f_2 * sig1_359[k]
                   + f_3 * pc_z[k] * sih_503[k];

        t_672[k] = f_1 * sig0_360[k]
                   - f_2 * sig1_360[k]
                   + f_3 * pc_x[k] * sih_504[k];

        t_673[k] = f_13 * shh_378[k]
                   + f_3 * pc_y[k] * sih_504[k];

        t_674[k] = f_13 * shh_357[k]
                   + f_3 * pc_z[k] * sih_504[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, pc_x, pc_y, shh_380, sig0_363, sig0_365, \
                         sig1_363, sig1_365, sih_506, sih_507, \
                         sih_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = f_4 * sig0_363[k]
                   - f_5 * sig1_363[k]
                   + f_3 * pc_x[k] * sih_507[k];

        t_676[k] = f_13 * shh_380[k]
                   + f_3 * pc_y[k] * sih_506[k];

        t_677[k] = f_4 * sig0_365[k]
                   - f_5 * sig1_365[k]
                   + f_3 * pc_x[k] * sih_509[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, pc_x, pc_y, pc_z, shh_360, shh_383, sig0_366, \
                         sig1_366, sih_507, sih_509, sih_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = f_6 * sig0_366[k]
                   - f_7 * sig1_366[k]
                   + f_3 * pc_x[k] * sih_510[k];

        t_679[k] = f_13 * shh_360[k]
                   + f_3 * pc_z[k] * sih_507[k];

        t_680[k] = f_13 * shh_383[k]
                   + f_3 * pc_y[k] * sih_509[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, pc_x, pc_z, shh_363, sig0_369, sig0_370, \
                         sig1_369, sig1_370, sih_510, sih_513, \
                         sih_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = f_6 * sig0_369[k]
                   - f_7 * sig1_369[k]
                   + f_3 * pc_x[k] * sih_513[k];

        t_682[k] = f_8 * sig0_370[k]
                   - f_9 * sig1_370[k]
                   + f_3 * pc_x[k] * sih_514[k];

        t_683[k] = f_13 * shh_363[k]
                   + f_3 * pc_z[k] * sih_510[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, pc_x, pc_y, shh_387, sig0_372, sig0_374, \
                         sig1_372, sig1_374, sih_513, sih_516, sih_518, \
                         sih_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_8 * sig0_372[k]
                   - f_9 * sig1_372[k]
                   + f_3 * pc_x[k] * sih_516[k];

        t_685[k] = f_13 * shh_387[k]
                   + f_3 * pc_y[k] * sih_513[k];

        t_686[k] = f_8 * sig0_374[k]
                   - f_9 * sig1_374[k]
                   + f_3 * pc_x[k] * sih_518[k];

        t_687[k] = f_3 * pc_x[k] * sih_519[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, t_691, t_692, pc_x, sih_520, sih_521, sih_522, \
                         sih_523, sih_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_3 * pc_x[k] * sih_520[k];

        t_689[k] = f_3 * pc_x[k] * sih_521[k];

        t_690[k] = f_3 * pc_x[k] * sih_522[k];

        t_691[k] = f_3 * pc_x[k] * sih_523[k];

        t_692[k] = f_3 * pc_x[k] * sih_524[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, pc_y, pc_z, shh_372, shh_393, shh_395, sig0_370, \
                         sig0_372, sig1_370, sig1_372, sih_519, \
                         sih_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = f_13 * shh_393[k]
                   + f_1 * sig0_370[k]
                   - f_2 * sig1_370[k]
                   + f_3 * pc_y[k] * sih_519[k];

        t_694[k] = f_13 * shh_372[k]
                   + f_3 * pc_z[k] * sih_519[k];

        t_695[k] = f_13 * shh_395[k]
                   + f_4 * sig0_372[k]
                   - f_5 * sig1_372[k]
                   + f_3 * pc_y[k] * sih_521[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, pc_y, shh_396, shh_397, shh_398, sig0_373, \
                         sig0_374, sig1_373, sig1_374, sih_522, sih_523, \
                         sih_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = f_13 * shh_396[k]
                   + f_6 * sig0_373[k]
                   - f_7 * sig1_373[k]
                   + f_3 * pc_y[k] * sih_522[k];

        t_697[k] = f_13 * shh_397[k]
                   + f_8 * sig0_374[k]
                   - f_9 * sig1_374[k]
                   + f_3 * pc_y[k] * sih_523[k];

        t_698[k] = f_13 * shh_398[k]
                   + f_3 * pc_y[k] * sih_524[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, pc_x, pc_y, pc_z, shh_377, shh_378, \
                         shh_399, sig0_374, sig0_375, sig1_374, sig1_375, sih_524, \
                         sih_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_13 * shh_377[k]
                   + f_1 * sig0_374[k]
                   - f_2 * sig1_374[k]
                   + f_3 * pc_z[k] * sih_524[k];

        t_700[k] = f_1 * sig0_375[k]
                   - f_2 * sig1_375[k]
                   + f_3 * pc_x[k] * sih_525[k];

        t_701[k] = f_12 * shh_399[k]
                   + f_3 * pc_y[k] * sih_525[k];

        t_702[k] = f_14 * shh_378[k]
                   + f_3 * pc_z[k] * sih_525[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, pc_x, pc_y, shh_401, sig0_378, sig0_380, \
                         sig1_378, sig1_380, sih_527, sih_528, \
                         sih_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_4 * sig0_378[k]
                   - f_5 * sig1_378[k]
                   + f_3 * pc_x[k] * sih_528[k];

        t_704[k] = f_12 * shh_401[k]
                   + f_3 * pc_y[k] * sih_527[k];

        t_705[k] = f_4 * sig0_380[k]
                   - f_5 * sig1_380[k]
                   + f_3 * pc_x[k] * sih_530[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, pc_x, pc_y, pc_z, shh_381, shh_404, sig0_381, \
                         sig1_381, sih_528, sih_530, sih_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_6 * sig0_381[k]
                   - f_7 * sig1_381[k]
                   + f_3 * pc_x[k] * sih_531[k];

        t_707[k] = f_14 * shh_381[k]
                   + f_3 * pc_z[k] * sih_528[k];

        t_708[k] = f_12 * shh_404[k]
                   + f_3 * pc_y[k] * sih_530[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_x, pc_z, shh_384, sig0_384, sig0_385, \
                         sig1_384, sig1_385, sih_531, sih_534, \
                         sih_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_6 * sig0_384[k]
                   - f_7 * sig1_384[k]
                   + f_3 * pc_x[k] * sih_534[k];

        t_710[k] = f_8 * sig0_385[k]
                   - f_9 * sig1_385[k]
                   + f_3 * pc_x[k] * sih_535[k];

        t_711[k] = f_14 * shh_384[k]
                   + f_3 * pc_z[k] * sih_531[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pc_x, pc_y, shh_408, sig0_387, sig0_389, \
                         sig1_387, sig1_389, sih_534, sih_537, sih_539, \
                         sih_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_8 * sig0_387[k]
                   - f_9 * sig1_387[k]
                   + f_3 * pc_x[k] * sih_537[k];

        t_713[k] = f_12 * shh_408[k]
                   + f_3 * pc_y[k] * sih_534[k];

        t_714[k] = f_8 * sig0_389[k]
                   - f_9 * sig1_389[k]
                   + f_3 * pc_x[k] * sih_539[k];

        t_715[k] = f_3 * pc_x[k] * sih_540[k];
    }
}

static auto
compute_prim_sii_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shi0,
                                                          const size_t shh, const size_t shi1,
                                                          const size_t sig0, const size_t sig1,
                                                          const size_t sih, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_15 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shi0_560 = buffer.data(shi0 + 560);
    const auto *shi0_565 = buffer.data(shi0 + 565);
    const auto *shi0_569 = buffer.data(shi0 + 569);
    const auto *shi0_574 = buffer.data(shi0 + 574);
    const auto *shi0_581 = buffer.data(shi0 + 581);
    const auto *shi0_583 = buffer.data(shi0 + 583);
    const auto *shi0_584 = buffer.data(shi0 + 584);
    const auto *shi0_585 = buffer.data(shi0 + 585);
    const auto *shi0_587 = buffer.data(shi0 + 587);

    const auto *shh_393 = buffer.data(shh + 393);
    const auto *shh_398 = buffer.data(shh + 398);
    const auto *shh_399 = buffer.data(shh + 399);
    const auto *shh_402 = buffer.data(shh + 402);
    const auto *shh_405 = buffer.data(shh + 405);
    const auto *shh_414 = buffer.data(shh + 414);
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
    const auto *shh_435 = buffer.data(shh + 435);
    const auto *shh_437 = buffer.data(shh + 437);
    const auto *shh_438 = buffer.data(shh + 438);
    const auto *shh_439 = buffer.data(shh + 439);
    const auto *shh_440 = buffer.data(shh + 440);

    const auto *shi1_560 = buffer.data(shi1 + 560);
    const auto *shi1_565 = buffer.data(shi1 + 565);
    const auto *shi1_569 = buffer.data(shi1 + 569);
    const auto *shi1_574 = buffer.data(shi1 + 574);
    const auto *shi1_581 = buffer.data(shi1 + 581);
    const auto *shi1_583 = buffer.data(shi1 + 583);
    const auto *shi1_584 = buffer.data(shi1 + 584);
    const auto *shi1_585 = buffer.data(shi1 + 585);
    const auto *shi1_587 = buffer.data(shi1 + 587);

    const auto *sig0_385 = buffer.data(sig0 + 385);
    const auto *sig0_387 = buffer.data(sig0 + 387);
    const auto *sig0_388 = buffer.data(sig0 + 388);
    const auto *sig0_389 = buffer.data(sig0 + 389);
    const auto *sig0_393 = buffer.data(sig0 + 393);
    const auto *sig0_396 = buffer.data(sig0 + 396);
    const auto *sig0_400 = buffer.data(sig0 + 400);
    const auto *sig0_402 = buffer.data(sig0 + 402);
    const auto *sig0_405 = buffer.data(sig0 + 405);
    const auto *sig0_408 = buffer.data(sig0 + 408);
    const auto *sig0_410 = buffer.data(sig0 + 410);
    const auto *sig0_411 = buffer.data(sig0 + 411);
    const auto *sig0_414 = buffer.data(sig0 + 414);
    const auto *sig0_415 = buffer.data(sig0 + 415);
    const auto *sig0_417 = buffer.data(sig0 + 417);
    const auto *sig0_418 = buffer.data(sig0 + 418);
    const auto *sig0_419 = buffer.data(sig0 + 419);

    const auto *sig1_385 = buffer.data(sig1 + 385);
    const auto *sig1_387 = buffer.data(sig1 + 387);
    const auto *sig1_388 = buffer.data(sig1 + 388);
    const auto *sig1_389 = buffer.data(sig1 + 389);
    const auto *sig1_393 = buffer.data(sig1 + 393);
    const auto *sig1_396 = buffer.data(sig1 + 396);
    const auto *sig1_400 = buffer.data(sig1 + 400);
    const auto *sig1_402 = buffer.data(sig1 + 402);
    const auto *sig1_405 = buffer.data(sig1 + 405);
    const auto *sig1_408 = buffer.data(sig1 + 408);
    const auto *sig1_410 = buffer.data(sig1 + 410);
    const auto *sig1_411 = buffer.data(sig1 + 411);
    const auto *sig1_414 = buffer.data(sig1 + 414);
    const auto *sig1_415 = buffer.data(sig1 + 415);
    const auto *sig1_417 = buffer.data(sig1 + 417);
    const auto *sig1_418 = buffer.data(sig1 + 418);
    const auto *sig1_419 = buffer.data(sig1 + 419);

    const auto *sih_540 = buffer.data(sih + 540);
    const auto *sih_541 = buffer.data(sih + 541);
    const auto *sih_542 = buffer.data(sih + 542);
    const auto *sih_543 = buffer.data(sih + 543);
    const auto *sih_544 = buffer.data(sih + 544);
    const auto *sih_545 = buffer.data(sih + 545);
    const auto *sih_546 = buffer.data(sih + 546);
    const auto *sih_548 = buffer.data(sih + 548);
    const auto *sih_549 = buffer.data(sih + 549);
    const auto *sih_551 = buffer.data(sih + 551);
    const auto *sih_552 = buffer.data(sih + 552);
    const auto *sih_555 = buffer.data(sih + 555);
    const auto *sih_556 = buffer.data(sih + 556);
    const auto *sih_558 = buffer.data(sih + 558);
    const auto *sih_561 = buffer.data(sih + 561);
    const auto *sih_562 = buffer.data(sih + 562);
    const auto *sih_563 = buffer.data(sih + 563);
    const auto *sih_564 = buffer.data(sih + 564);
    const auto *sih_565 = buffer.data(sih + 565);
    const auto *sih_566 = buffer.data(sih + 566);
    const auto *sih_567 = buffer.data(sih + 567);
    const auto *sih_569 = buffer.data(sih + 569);
    const auto *sih_570 = buffer.data(sih + 570);
    const auto *sih_572 = buffer.data(sih + 572);
    const auto *sih_573 = buffer.data(sih + 573);
    const auto *sih_576 = buffer.data(sih + 576);
    const auto *sih_577 = buffer.data(sih + 577);
    const auto *sih_579 = buffer.data(sih + 579);
    const auto *sih_581 = buffer.data(sih + 581);
    const auto *sih_582 = buffer.data(sih + 582);
    const auto *sih_583 = buffer.data(sih + 583);
    const auto *sih_584 = buffer.data(sih + 584);
    const auto *sih_585 = buffer.data(sih + 585);
    const auto *sih_586 = buffer.data(sih + 586);
    const auto *sih_587 = buffer.data(sih + 587);

#pragma omp simd aligned(t_716, t_717, t_718, t_719, t_720, pc_x, sih_541, sih_542, sih_543, \
                         sih_544, sih_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_3 * pc_x[k] * sih_541[k];

        t_717[k] = f_3 * pc_x[k] * sih_542[k];

        t_718[k] = f_3 * pc_x[k] * sih_543[k];

        t_719[k] = f_3 * pc_x[k] * sih_544[k];

        t_720[k] = f_3 * pc_x[k] * sih_545[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, pc_y, pc_z, shh_393, shh_414, shh_416, sig0_385, \
                         sig0_387, sig1_385, sig1_387, sih_540, \
                         sih_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = f_12 * shh_414[k]
                   + f_1 * sig0_385[k]
                   - f_2 * sig1_385[k]
                   + f_3 * pc_y[k] * sih_540[k];

        t_722[k] = f_14 * shh_393[k]
                   + f_3 * pc_z[k] * sih_540[k];

        t_723[k] = f_12 * shh_416[k]
                   + f_4 * sig0_387[k]
                   - f_5 * sig1_387[k]
                   + f_3 * pc_y[k] * sih_542[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, pc_y, shh_417, shh_418, shh_419, sig0_388, \
                         sig0_389, sig1_388, sig1_389, sih_543, sih_544, \
                         sih_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_12 * shh_417[k]
                   + f_6 * sig0_388[k]
                   - f_7 * sig1_388[k]
                   + f_3 * pc_y[k] * sih_543[k];

        t_725[k] = f_12 * shh_418[k]
                   + f_8 * sig0_389[k]
                   - f_9 * sig1_389[k]
                   + f_3 * pc_y[k] * sih_544[k];

        t_726[k] = f_12 * shh_419[k]
                   + f_3 * pc_y[k] * sih_545[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, t_730, pb_y, pc_y, pc_z, shi0_560, shh_398, \
                         shh_399, shh_420, shi1_560, sig0_389, sig1_389, sih_545, \
                         sih_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_14 * shh_398[k]
                   + f_1 * sig0_389[k]
                   - f_2 * sig1_389[k]
                   + f_3 * pc_z[k] * sih_545[k];

        t_728[k] = pb_y[k] * shi0_560[k]
                   - f_10 * pc_y[k] * shi1_560[k];

        t_729[k] = f_11 * shh_420[k]
                   + f_3 * pc_y[k] * sih_546[k];

        t_730[k] = f_15 * shh_399[k]
                   + f_3 * pc_z[k] * sih_546[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, pb_y, pc_x, pc_y, shi0_565, shh_422, shi1_565, \
                         sig0_393, sig1_393, sih_548, sih_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = f_4 * sig0_393[k]
                   - f_5 * sig1_393[k]
                   + f_3 * pc_x[k] * sih_549[k];

        t_732[k] = f_11 * shh_422[k]
                   + f_3 * pc_y[k] * sih_548[k];

        t_733[k] = pb_y[k] * shi0_565[k]
                   - f_10 * pc_y[k] * shi1_565[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, pc_x, pc_y, pc_z, shh_402, shh_425, sig0_396, \
                         sig1_396, sih_549, sih_551, sih_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = f_6 * sig0_396[k]
                   - f_7 * sig1_396[k]
                   + f_3 * pc_x[k] * sih_552[k];

        t_735[k] = f_15 * shh_402[k]
                   + f_3 * pc_z[k] * sih_549[k];

        t_736[k] = f_11 * shh_425[k]
                   + f_3 * pc_y[k] * sih_551[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pb_y, pc_x, pc_y, pc_z, shi0_569, shh_405, \
                         shi1_569, sig0_400, sig1_400, sih_552, \
                         sih_556 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = pb_y[k] * shi0_569[k]
                   - f_10 * pc_y[k] * shi1_569[k];

        t_738[k] = f_8 * sig0_400[k]
                   - f_9 * sig1_400[k]
                   + f_3 * pc_x[k] * sih_556[k];

        t_739[k] = f_15 * shh_405[k]
                   + f_3 * pc_z[k] * sih_552[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pb_y, pc_x, pc_y, shi0_574, shh_429, \
                         shi1_574, sig0_402, sig1_402, sih_555, sih_558, \
                         sih_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_8 * sig0_402[k]
                   - f_9 * sig1_402[k]
                   + f_3 * pc_x[k] * sih_558[k];

        t_741[k] = f_11 * shh_429[k]
                   + f_3 * pc_y[k] * sih_555[k];

        t_742[k] = pb_y[k] * shi0_574[k]
                   - f_10 * pc_y[k] * shi1_574[k];

        t_743[k] = f_3 * pc_x[k] * sih_561[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, t_748, pc_x, sih_562, sih_563, sih_564, \
                         sih_565, sih_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_3 * pc_x[k] * sih_562[k];

        t_745[k] = f_3 * pc_x[k] * sih_563[k];

        t_746[k] = f_3 * pc_x[k] * sih_564[k];

        t_747[k] = f_3 * pc_x[k] * sih_565[k];

        t_748[k] = f_3 * pc_x[k] * sih_566[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, pb_y, pc_y, pc_z, shi0_581, shi0_583, shh_414, \
                         shh_435, shh_437, shi1_581, shi1_583, \
                         sih_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = pb_y[k] * shi0_581[k]
                   + f_0 * shh_435[k]
                   - f_10 * pc_y[k] * shi1_581[k];

        t_750[k] = f_15 * shh_414[k]
                   + f_3 * pc_z[k] * sih_561[k];

        t_751[k] = pb_y[k] * shi0_583[k]
                   + f_14 * shh_437[k]
                   - f_10 * pc_y[k] * shi1_583[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pb_y, pc_y, shi0_584, shi0_585, shi0_587, \
                         shh_438, shh_439, shh_440, shi1_584, shi1_585, shi1_587, \
                         sih_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = pb_y[k] * shi0_584[k]
                   + f_13 * shh_438[k]
                   - f_10 * pc_y[k] * shi1_584[k];

        t_753[k] = pb_y[k] * shi0_585[k]
                   + f_12 * shh_439[k]
                   - f_10 * pc_y[k] * shi1_585[k];

        t_754[k] = f_11 * shh_440[k]
                   + f_3 * pc_y[k] * sih_566[k];

        t_755[k] = pb_y[k] * shi0_587[k]
                   - f_10 * pc_y[k] * shi1_587[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, t_760, pc_x, pc_y, pc_z, shh_420, \
                         sig0_405, sig0_408, sig1_405, sig1_408, sih_567, sih_569, \
                         sih_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_1 * sig0_405[k]
                   - f_2 * sig1_405[k]
                   + f_3 * pc_x[k] * sih_567[k];

        t_757[k] = f_3 * pc_y[k] * sih_567[k];

        t_758[k] = f_0 * shh_420[k]
                   + f_3 * pc_z[k] * sih_567[k];

        t_759[k] = f_4 * sig0_408[k]
                   - f_5 * sig1_408[k]
                   + f_3 * pc_x[k] * sih_570[k];

        t_760[k] = f_3 * pc_y[k] * sih_569[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, t_764, pc_x, pc_y, pc_z, shh_423, sig0_410, \
                         sig0_411, sig1_410, sig1_411, sih_570, sih_572, \
                         sih_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_4 * sig0_410[k]
                   - f_5 * sig1_410[k]
                   + f_3 * pc_x[k] * sih_572[k];

        t_762[k] = f_6 * sig0_411[k]
                   - f_7 * sig1_411[k]
                   + f_3 * pc_x[k] * sih_573[k];

        t_763[k] = f_0 * shh_423[k]
                   + f_3 * pc_z[k] * sih_570[k];

        t_764[k] = f_3 * pc_y[k] * sih_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, pc_z, shh_426, sig0_414, sig0_415, \
                         sig1_414, sig1_415, sih_573, sih_576, \
                         sih_577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_6 * sig0_414[k]
                   - f_7 * sig1_414[k]
                   + f_3 * pc_x[k] * sih_576[k];

        t_766[k] = f_8 * sig0_415[k]
                   - f_9 * sig1_415[k]
                   + f_3 * pc_x[k] * sih_577[k];

        t_767[k] = f_0 * shh_426[k]
                   + f_3 * pc_z[k] * sih_573[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, t_772, pc_x, pc_y, sig0_417, sig0_419, \
                         sig1_417, sig1_419, sih_576, sih_579, sih_581, sih_582, \
                         sih_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_8 * sig0_417[k]
                   - f_9 * sig1_417[k]
                   + f_3 * pc_x[k] * sih_579[k];

        t_769[k] = f_3 * pc_y[k] * sih_576[k];

        t_770[k] = f_8 * sig0_419[k]
                   - f_9 * sig1_419[k]
                   + f_3 * pc_x[k] * sih_581[k];

        t_771[k] = f_3 * pc_x[k] * sih_582[k];

        t_772[k] = f_3 * pc_x[k] * sih_583[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, t_776, t_777, pc_x, pc_y, sig0_415, sig1_415, \
                         sih_582, sih_584, sih_585, sih_586, sih_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_3 * pc_x[k] * sih_584[k];

        t_774[k] = f_3 * pc_x[k] * sih_585[k];

        t_775[k] = f_3 * pc_x[k] * sih_586[k];

        t_776[k] = f_3 * pc_x[k] * sih_587[k];

        t_777[k] = f_1 * sig0_415[k]
                   - f_2 * sig1_415[k]
                   + f_3 * pc_y[k] * sih_582[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, pc_y, pc_z, shh_435, sig0_417, sig0_418, \
                         sig1_417, sig1_418, sih_582, sih_584, \
                         sih_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_0 * shh_435[k]
                   + f_3 * pc_z[k] * sih_582[k];

        t_779[k] = f_4 * sig0_417[k]
                   - f_5 * sig1_417[k]
                   + f_3 * pc_y[k] * sih_584[k];

        t_780[k] = f_6 * sig0_418[k]
                   - f_7 * sig1_418[k]
                   + f_3 * pc_y[k] * sih_585[k];
    }

#pragma omp simd aligned(t_781, t_782, t_783, pc_y, pc_z, shh_440, sig0_419, sig1_419, \
                         sih_586, sih_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_781[k] = f_8 * sig0_419[k]
                   - f_9 * sig1_419[k]
                   + f_3 * pc_y[k] * sih_586[k];

        t_782[k] = f_3 * pc_y[k] * sih_587[k];

        t_783[k] = f_0 * shh_440[k]
                   + f_1 * sig0_419[k]
                   - f_2 * sig1_419[k]
                   + f_3 * pc_z[k] * sih_587[k];
    }
}

auto
compute_prim_sii_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t shi0, const size_t shh,
                                                   const size_t shi1, const size_t sig0,
                                                   const size_t sig1, const size_t sih,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sii_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, shi0, shh,
                                                              shi1, sig0, sig1, sih, ncols,
                                                              gamma, p, q);

    compute_prim_sii_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, shi0, shh,
                                                              shi1, sig0, sig1, sih, ncols,
                                                              gamma, p, q);

    compute_prim_sii_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, shi0, shh,
                                                              shi1, sig0, sig1, sih, ncols,
                                                              gamma, p, q);

    compute_prim_sii_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, shi0, shh,
                                                              shi1, sig0, sig1, sih, ncols,
                                                              gamma, p, q);

    compute_prim_sii_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, shi0, shh,
                                                              shi1, sig0, sig1, sih, ncols,
                                                              gamma, p, q);

    compute_prim_sii_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, shi0, shh,
                                                              shi1, sig0, sig1, sih, ncols,
                                                              gamma, p, q);

    compute_prim_sii_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, shi0, shh,
                                                              shi1, sig0, sig1, sih, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
