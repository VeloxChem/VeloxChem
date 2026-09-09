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


#include "SimdThreeCenterElectronRepulsionVrrRecSHK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_shk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgk0,
                                                          const size_t sgi, const size_t sgk1,
                                                          const size_t shh0, const size_t shh1,
                                                          const size_t shi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgk0_0 = buffer.data(sgk0 + 0);
    const auto *sgk0_3 = buffer.data(sgk0 + 3);
    const auto *sgk0_5 = buffer.data(sgk0 + 5);
    const auto *sgk0_6 = buffer.data(sgk0 + 6);
    const auto *sgk0_9 = buffer.data(sgk0 + 9);
    const auto *sgk0_10 = buffer.data(sgk0 + 10);
    const auto *sgk0_12 = buffer.data(sgk0 + 12);
    const auto *sgk0_14 = buffer.data(sgk0 + 14);
    const auto *sgk0_15 = buffer.data(sgk0 + 15);
    const auto *sgk0_17 = buffer.data(sgk0 + 17);
    const auto *sgk0_18 = buffer.data(sgk0 + 18);
    const auto *sgk0_20 = buffer.data(sgk0 + 20);
    const auto *sgk0_28 = buffer.data(sgk0 + 28);
    const auto *sgk0_35 = buffer.data(sgk0 + 35);

    const auto *sgi_0 = buffer.data(sgi + 0);
    const auto *sgi_1 = buffer.data(sgi + 1);
    const auto *sgi_2 = buffer.data(sgi + 2);
    const auto *sgi_3 = buffer.data(sgi + 3);
    const auto *sgi_5 = buffer.data(sgi + 5);
    const auto *sgi_6 = buffer.data(sgi + 6);
    const auto *sgi_7 = buffer.data(sgi + 7);
    const auto *sgi_8 = buffer.data(sgi + 8);
    const auto *sgi_9 = buffer.data(sgi + 9);
    const auto *sgi_10 = buffer.data(sgi + 10);
    const auto *sgi_11 = buffer.data(sgi + 11);
    const auto *sgi_12 = buffer.data(sgi + 12);
    const auto *sgi_13 = buffer.data(sgi + 13);
    const auto *sgi_14 = buffer.data(sgi + 14);
    const auto *sgi_15 = buffer.data(sgi + 15);
    const auto *sgi_17 = buffer.data(sgi + 17);
    const auto *sgi_18 = buffer.data(sgi + 18);
    const auto *sgi_20 = buffer.data(sgi + 20);
    const auto *sgi_21 = buffer.data(sgi + 21);
    const auto *sgi_22 = buffer.data(sgi + 22);
    const auto *sgi_23 = buffer.data(sgi + 23);
    const auto *sgi_24 = buffer.data(sgi + 24);
    const auto *sgi_25 = buffer.data(sgi + 25);
    const auto *sgi_26 = buffer.data(sgi + 26);
    const auto *sgi_27 = buffer.data(sgi + 27);
    const auto *sgi_28 = buffer.data(sgi + 28);
    const auto *sgi_30 = buffer.data(sgi + 30);
    const auto *sgi_33 = buffer.data(sgi + 33);
    const auto *sgi_37 = buffer.data(sgi + 37);
    const auto *sgi_49 = buffer.data(sgi + 49);
    const auto *sgi_50 = buffer.data(sgi + 50);
    const auto *sgi_51 = buffer.data(sgi + 51);
    const auto *sgi_52 = buffer.data(sgi + 52);
    const auto *sgi_53 = buffer.data(sgi + 53);
    const auto *sgi_54 = buffer.data(sgi + 54);
    const auto *sgi_55 = buffer.data(sgi + 55);
    const auto *sgi_77 = buffer.data(sgi + 77);
    const auto *sgi_78 = buffer.data(sgi + 78);
    const auto *sgi_79 = buffer.data(sgi + 79);
    const auto *sgi_80 = buffer.data(sgi + 80);
    const auto *sgi_81 = buffer.data(sgi + 81);
    const auto *sgi_82 = buffer.data(sgi + 82);
    const auto *sgi_83 = buffer.data(sgi + 83);
    const auto *sgi_84 = buffer.data(sgi + 84);
    const auto *sgi_87 = buffer.data(sgi + 87);
    const auto *sgi_89 = buffer.data(sgi + 89);
    const auto *sgi_90 = buffer.data(sgi + 90);
    const auto *sgi_93 = buffer.data(sgi + 93);
    const auto *sgi_94 = buffer.data(sgi + 94);
    const auto *sgi_96 = buffer.data(sgi + 96);

    const auto *sgk1_0 = buffer.data(sgk1 + 0);
    const auto *sgk1_3 = buffer.data(sgk1 + 3);
    const auto *sgk1_5 = buffer.data(sgk1 + 5);
    const auto *sgk1_6 = buffer.data(sgk1 + 6);
    const auto *sgk1_9 = buffer.data(sgk1 + 9);
    const auto *sgk1_10 = buffer.data(sgk1 + 10);
    const auto *sgk1_12 = buffer.data(sgk1 + 12);
    const auto *sgk1_14 = buffer.data(sgk1 + 14);
    const auto *sgk1_15 = buffer.data(sgk1 + 15);
    const auto *sgk1_17 = buffer.data(sgk1 + 17);
    const auto *sgk1_18 = buffer.data(sgk1 + 18);
    const auto *sgk1_20 = buffer.data(sgk1 + 20);
    const auto *sgk1_28 = buffer.data(sgk1 + 28);
    const auto *sgk1_35 = buffer.data(sgk1 + 35);

    const auto *shh0_0 = buffer.data(shh0 + 0);
    const auto *shh0_3 = buffer.data(shh0 + 3);
    const auto *shh0_5 = buffer.data(shh0 + 5);
    const auto *shh0_6 = buffer.data(shh0 + 6);
    const auto *shh0_9 = buffer.data(shh0 + 9);
    const auto *shh0_10 = buffer.data(shh0 + 10);
    const auto *shh0_12 = buffer.data(shh0 + 12);
    const auto *shh0_14 = buffer.data(shh0 + 14);
    const auto *shh0_15 = buffer.data(shh0 + 15);
    const auto *shh0_17 = buffer.data(shh0 + 17);
    const auto *shh0_18 = buffer.data(shh0 + 18);
    const auto *shh0_19 = buffer.data(shh0 + 19);
    const auto *shh0_20 = buffer.data(shh0 + 20);
    const auto *shh0_36 = buffer.data(shh0 + 36);
    const auto *shh0_38 = buffer.data(shh0 + 38);
    const auto *shh0_39 = buffer.data(shh0 + 39);
    const auto *shh0_40 = buffer.data(shh0 + 40);
    const auto *shh0_41 = buffer.data(shh0 + 41);
    const auto *shh0_59 = buffer.data(shh0 + 59);
    const auto *shh0_60 = buffer.data(shh0 + 60);
    const auto *shh0_61 = buffer.data(shh0 + 61);
    const auto *shh0_62 = buffer.data(shh0 + 62);
    const auto *shh0_63 = buffer.data(shh0 + 63);
    const auto *shh0_66 = buffer.data(shh0 + 66);
    const auto *shh0_68 = buffer.data(shh0 + 68);
    const auto *shh0_69 = buffer.data(shh0 + 69);
    const auto *shh0_72 = buffer.data(shh0 + 72);
    const auto *shh0_73 = buffer.data(shh0 + 73);
    const auto *shh0_75 = buffer.data(shh0 + 75);

    const auto *shh1_0 = buffer.data(shh1 + 0);
    const auto *shh1_3 = buffer.data(shh1 + 3);
    const auto *shh1_5 = buffer.data(shh1 + 5);
    const auto *shh1_6 = buffer.data(shh1 + 6);
    const auto *shh1_9 = buffer.data(shh1 + 9);
    const auto *shh1_10 = buffer.data(shh1 + 10);
    const auto *shh1_12 = buffer.data(shh1 + 12);
    const auto *shh1_14 = buffer.data(shh1 + 14);
    const auto *shh1_15 = buffer.data(shh1 + 15);
    const auto *shh1_17 = buffer.data(shh1 + 17);
    const auto *shh1_18 = buffer.data(shh1 + 18);
    const auto *shh1_19 = buffer.data(shh1 + 19);
    const auto *shh1_20 = buffer.data(shh1 + 20);
    const auto *shh1_36 = buffer.data(shh1 + 36);
    const auto *shh1_38 = buffer.data(shh1 + 38);
    const auto *shh1_39 = buffer.data(shh1 + 39);
    const auto *shh1_40 = buffer.data(shh1 + 40);
    const auto *shh1_41 = buffer.data(shh1 + 41);
    const auto *shh1_59 = buffer.data(shh1 + 59);
    const auto *shh1_60 = buffer.data(shh1 + 60);
    const auto *shh1_61 = buffer.data(shh1 + 61);
    const auto *shh1_62 = buffer.data(shh1 + 62);
    const auto *shh1_63 = buffer.data(shh1 + 63);
    const auto *shh1_66 = buffer.data(shh1 + 66);
    const auto *shh1_68 = buffer.data(shh1 + 68);
    const auto *shh1_69 = buffer.data(shh1 + 69);
    const auto *shh1_72 = buffer.data(shh1 + 72);
    const auto *shh1_73 = buffer.data(shh1 + 73);
    const auto *shh1_75 = buffer.data(shh1 + 75);

    const auto *shi_0 = buffer.data(shi + 0);
    const auto *shi_2 = buffer.data(shi + 2);
    const auto *shi_3 = buffer.data(shi + 3);
    const auto *shi_5 = buffer.data(shi + 5);
    const auto *shi_6 = buffer.data(shi + 6);
    const auto *shi_9 = buffer.data(shi + 9);
    const auto *shi_10 = buffer.data(shi + 10);
    const auto *shi_12 = buffer.data(shi + 12);
    const auto *shi_14 = buffer.data(shi + 14);
    const auto *shi_15 = buffer.data(shi + 15);
    const auto *shi_17 = buffer.data(shi + 17);
    const auto *shi_18 = buffer.data(shi + 18);
    const auto *shi_20 = buffer.data(shi + 20);
    const auto *shi_21 = buffer.data(shi + 21);
    const auto *shi_22 = buffer.data(shi + 22);
    const auto *shi_23 = buffer.data(shi + 23);
    const auto *shi_24 = buffer.data(shi + 24);
    const auto *shi_25 = buffer.data(shi + 25);
    const auto *shi_26 = buffer.data(shi + 26);
    const auto *shi_27 = buffer.data(shi + 27);
    const auto *shi_28 = buffer.data(shi + 28);
    const auto *shi_30 = buffer.data(shi + 30);
    const auto *shi_31 = buffer.data(shi + 31);
    const auto *shi_33 = buffer.data(shi + 33);
    const auto *shi_34 = buffer.data(shi + 34);
    const auto *shi_37 = buffer.data(shi + 37);
    const auto *shi_38 = buffer.data(shi + 38);
    const auto *shi_42 = buffer.data(shi + 42);
    const auto *shi_49 = buffer.data(shi + 49);
    const auto *shi_50 = buffer.data(shi + 50);
    const auto *shi_51 = buffer.data(shi + 51);
    const auto *shi_52 = buffer.data(shi + 52);
    const auto *shi_53 = buffer.data(shi + 53);
    const auto *shi_54 = buffer.data(shi + 54);
    const auto *shi_55 = buffer.data(shi + 55);
    const auto *shi_56 = buffer.data(shi + 56);
    const auto *shi_58 = buffer.data(shi + 58);
    const auto *shi_59 = buffer.data(shi + 59);
    const auto *shi_61 = buffer.data(shi + 61);
    const auto *shi_62 = buffer.data(shi + 62);
    const auto *shi_65 = buffer.data(shi + 65);
    const auto *shi_66 = buffer.data(shi + 66);
    const auto *shi_70 = buffer.data(shi + 70);
    const auto *shi_77 = buffer.data(shi + 77);
    const auto *shi_78 = buffer.data(shi + 78);
    const auto *shi_79 = buffer.data(shi + 79);
    const auto *shi_80 = buffer.data(shi + 80);
    const auto *shi_81 = buffer.data(shi + 81);
    const auto *shi_82 = buffer.data(shi + 82);
    const auto *shi_83 = buffer.data(shi + 83);
    const auto *shi_84 = buffer.data(shi + 84);
    const auto *shi_86 = buffer.data(shi + 86);
    const auto *shi_87 = buffer.data(shi + 87);
    const auto *shi_89 = buffer.data(shi + 89);
    const auto *shi_90 = buffer.data(shi + 90);
    const auto *shi_93 = buffer.data(shi + 93);
    const auto *shi_94 = buffer.data(shi + 94);
    const auto *shi_96 = buffer.data(shi + 96);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sgi_0, sgi_3, shh0_0, shh0_3, \
                         shh1_0, shh1_3, shi_0, shi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sgi_0[k]
                 + f_1 * shh0_0[k]
                 - f_2 * shh1_0[k]
                 + f_3 * pc_x[k] * shi_0[k];

        t_1[k] = f_3 * pc_y[k] * shi_0[k];

        t_2[k] = f_3 * pc_z[k] * shi_0[k];

        t_3[k] = f_0 * sgi_3[k]
                 + f_4 * shh0_3[k]
                 - f_5 * shh1_3[k]
                 + f_3 * pc_x[k] * shi_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sgi_5, sgi_6, shh0_5, shh0_6, shh1_5, \
                         shh1_6, shi_2, shi_5, shi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * shi_2[k];

        t_5[k] = f_0 * sgi_5[k]
                 + f_4 * shh0_5[k]
                 - f_5 * shh1_5[k]
                 + f_3 * pc_x[k] * shi_5[k];

        t_6[k] = f_0 * sgi_6[k]
                 + f_6 * shh0_6[k]
                 - f_7 * shh1_6[k]
                 + f_3 * pc_x[k] * shi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sgi_9, shh0_9, shh1_9, shi_3, shi_5, \
                         shi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * shi_3[k];

        t_8[k] = f_3 * pc_y[k] * shi_5[k];

        t_9[k] = f_0 * sgi_9[k]
                 + f_6 * shh0_9[k]
                 - f_7 * shh1_9[k]
                 + f_3 * pc_x[k] * shi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sgi_10, sgi_12, shh0_10, shh0_12, \
                         shh1_10, shh1_12, shi_6, shi_10, shi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sgi_10[k]
                  + f_8 * shh0_10[k]
                  - f_9 * shh1_10[k]
                  + f_3 * pc_x[k] * shi_10[k];

        t_11[k] = f_3 * pc_z[k] * shi_6[k];

        t_12[k] = f_0 * sgi_12[k]
                  + f_8 * shh0_12[k]
                  - f_9 * shh1_12[k]
                  + f_3 * pc_x[k] * shi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, sgi_14, sgi_15, shh0_14, shh0_15, \
                         shh1_14, shh1_15, shi_9, shi_14, shi_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * shi_9[k];

        t_14[k] = f_0 * sgi_14[k]
                  + f_8 * shh0_14[k]
                  - f_9 * shh1_14[k]
                  + f_3 * pc_x[k] * shi_14[k];

        t_15[k] = f_0 * sgi_15[k]
                  + f_10 * shh0_15[k]
                  - f_11 * shh1_15[k]
                  + f_3 * pc_x[k] * shi_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, sgi_17, sgi_18, shh0_17, shh0_18, \
                         shh1_17, shh1_18, shi_10, shi_17, shi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * shi_10[k];

        t_17[k] = f_0 * sgi_17[k]
                  + f_10 * shh0_17[k]
                  - f_11 * shh1_17[k]
                  + f_3 * pc_x[k] * shi_17[k];

        t_18[k] = f_0 * sgi_18[k]
                  + f_10 * shh0_18[k]
                  - f_11 * shh1_18[k]
                  + f_3 * pc_x[k] * shi_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pc_x, pc_y, sgi_20, sgi_21, sgi_22, shh0_20, \
                         shh1_20, shi_14, shi_20, shi_21, shi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * shi_14[k];

        t_20[k] = f_0 * sgi_20[k]
                  + f_10 * shh0_20[k]
                  - f_11 * shh1_20[k]
                  + f_3 * pc_x[k] * shi_20[k];

        t_21[k] = f_0 * sgi_21[k]
                  + f_3 * pc_x[k] * shi_21[k];

        t_22[k] = f_0 * sgi_22[k]
                  + f_3 * pc_x[k] * shi_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pc_x, sgi_23, sgi_24, sgi_25, sgi_26, \
                         sgi_27, shi_23, shi_24, shi_25, shi_26, \
                         shi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * sgi_23[k]
                  + f_3 * pc_x[k] * shi_23[k];

        t_24[k] = f_0 * sgi_24[k]
                  + f_3 * pc_x[k] * shi_24[k];

        t_25[k] = f_0 * sgi_25[k]
                  + f_3 * pc_x[k] * shi_25[k];

        t_26[k] = f_0 * sgi_26[k]
                  + f_3 * pc_x[k] * shi_26[k];

        t_27[k] = f_0 * sgi_27[k]
                  + f_3 * pc_x[k] * shi_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pc_y, pc_z, shh0_15, shh0_17, shh0_18, \
                         shh1_15, shh1_17, shh1_18, shi_21, shi_23, \
                         shi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * shh0_15[k]
                  - f_2 * shh1_15[k]
                  + f_3 * pc_y[k] * shi_21[k];

        t_29[k] = f_3 * pc_z[k] * shi_21[k];

        t_30[k] = f_4 * shh0_17[k]
                  - f_5 * shh1_17[k]
                  + f_3 * pc_y[k] * shi_23[k];

        t_31[k] = f_6 * shh0_18[k]
                  - f_7 * shh1_18[k]
                  + f_3 * pc_y[k] * shi_24[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_y, pc_z, shh0_19, shh0_20, shh1_19, \
                         shh1_20, shi_25, shi_26, shi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_8 * shh0_19[k]
                  - f_9 * shh1_19[k]
                  + f_3 * pc_y[k] * shi_25[k];

        t_33[k] = f_10 * shh0_20[k]
                  - f_11 * shh1_20[k]
                  + f_3 * pc_y[k] * shi_26[k];

        t_34[k] = f_3 * pc_y[k] * shi_27[k];

        t_35[k] = f_1 * shh0_20[k]
                  - f_2 * shh1_20[k]
                  + f_3 * pc_z[k] * shi_27[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, sgk0_0, sgk0_3, sgi_0, \
                         sgi_1, sgk1_0, sgk1_3, shi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_y[k] * sgk0_0[k]
                  - f_12 * pc_y[k] * sgk1_0[k];

        t_37[k] = f_13 * sgi_0[k]
                  + f_3 * pc_y[k] * shi_28[k];

        t_38[k] = f_3 * pc_z[k] * shi_28[k];

        t_39[k] = pb_y[k] * sgk0_3[k]
                  + f_14 * sgi_1[k]
                  - f_12 * pc_y[k] * sgk1_3[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pc_y, pc_z, sgk0_5, sgk0_6, sgi_2, \
                         sgi_3, sgk1_5, sgk1_6, shi_30, shi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_13 * sgi_2[k]
                  + f_3 * pc_y[k] * shi_30[k];

        t_41[k] = pb_y[k] * sgk0_5[k]
                  - f_12 * pc_y[k] * sgk1_5[k];

        t_42[k] = pb_y[k] * sgk0_6[k]
                  + f_15 * sgi_3[k]
                  - f_12 * pc_y[k] * sgk1_6[k];

        t_43[k] = f_3 * pc_z[k] * shi_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pc_y, pc_z, sgk0_9, sgk0_10, sgi_5, \
                         sgi_6, sgk1_9, sgk1_10, shi_33, shi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_13 * sgi_5[k]
                  + f_3 * pc_y[k] * shi_33[k];

        t_45[k] = pb_y[k] * sgk0_9[k]
                  - f_12 * pc_y[k] * sgk1_9[k];

        t_46[k] = pb_y[k] * sgk0_10[k]
                  + f_16 * sgi_6[k]
                  - f_12 * pc_y[k] * sgk1_10[k];

        t_47[k] = f_3 * pc_z[k] * shi_34[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_y, pc_y, sgk0_12, sgk0_14, sgk0_15, sgi_8, \
                         sgi_9, sgi_10, sgk1_12, sgk1_14, sgk1_15, \
                         shi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pb_y[k] * sgk0_12[k]
                  + f_14 * sgi_8[k]
                  - f_12 * pc_y[k] * sgk1_12[k];

        t_49[k] = f_13 * sgi_9[k]
                  + f_3 * pc_y[k] * shi_37[k];

        t_50[k] = pb_y[k] * sgk0_14[k]
                  - f_12 * pc_y[k] * sgk1_14[k];

        t_51[k] = pb_y[k] * sgk0_15[k]
                  + f_0 * sgi_10[k]
                  - f_12 * pc_y[k] * sgk1_15[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_y, pc_y, pc_z, sgk0_17, sgk0_18, sgi_12, \
                         sgi_13, sgi_14, sgk1_17, sgk1_18, shi_38, \
                         shi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * pc_z[k] * shi_38[k];

        t_53[k] = pb_y[k] * sgk0_17[k]
                  + f_15 * sgi_12[k]
                  - f_12 * pc_y[k] * sgk1_17[k];

        t_54[k] = pb_y[k] * sgk0_18[k]
                  + f_14 * sgi_13[k]
                  - f_12 * pc_y[k] * sgk1_18[k];

        t_55[k] = f_13 * sgi_14[k]
                  + f_3 * pc_y[k] * shi_42[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pc_x, pc_y, sgk0_20, sgi_49, sgi_50, \
                         sgi_51, sgk1_20, shi_49, shi_50, shi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_y[k] * sgk0_20[k]
                  - f_12 * pc_y[k] * sgk1_20[k];

        t_57[k] = f_16 * sgi_49[k]
                  + f_3 * pc_x[k] * shi_49[k];

        t_58[k] = f_16 * sgi_50[k]
                  + f_3 * pc_x[k] * shi_50[k];

        t_59[k] = f_16 * sgi_51[k]
                  + f_3 * pc_x[k] * shi_51[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, sgi_52, sgi_53, sgi_54, sgi_55, shi_52, \
                         shi_53, shi_54, shi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_16 * sgi_52[k]
                  + f_3 * pc_x[k] * shi_52[k];

        t_61[k] = f_16 * sgi_53[k]
                  + f_3 * pc_x[k] * shi_53[k];

        t_62[k] = f_16 * sgi_54[k]
                  + f_3 * pc_x[k] * shi_54[k];

        t_63[k] = f_16 * sgi_55[k]
                  + f_3 * pc_x[k] * shi_55[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pc_y, pc_z, sgi_21, sgi_23, shh0_36, shh0_38, \
                         shh1_36, shh1_38, shi_49, shi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_13 * sgi_21[k]
                  + f_1 * shh0_36[k]
                  - f_2 * shh1_36[k]
                  + f_3 * pc_y[k] * shi_49[k];

        t_65[k] = f_3 * pc_z[k] * shi_49[k];

        t_66[k] = f_13 * sgi_23[k]
                  + f_4 * shh0_38[k]
                  - f_5 * shh1_38[k]
                  + f_3 * pc_y[k] * shi_51[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_y, sgi_24, sgi_25, sgi_26, shh0_39, shh0_40, \
                         shh0_41, shh1_39, shh1_40, shh1_41, shi_52, shi_53, \
                         shi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_13 * sgi_24[k]
                  + f_6 * shh0_39[k]
                  - f_7 * shh1_39[k]
                  + f_3 * pc_y[k] * shi_52[k];

        t_68[k] = f_13 * sgi_25[k]
                  + f_8 * shh0_40[k]
                  - f_9 * shh1_40[k]
                  + f_3 * pc_y[k] * shi_53[k];

        t_69[k] = f_13 * sgi_26[k]
                  + f_10 * shh0_41[k]
                  - f_11 * shh1_41[k]
                  + f_3 * pc_y[k] * shi_54[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_y, pb_z, pc_y, pc_z, sgk0_0, sgk0_35, \
                         sgi_27, sgk1_0, sgk1_35, shi_55, shi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_13 * sgi_27[k]
                  + f_3 * pc_y[k] * shi_55[k];

        t_71[k] = pb_y[k] * sgk0_35[k]
                  - f_12 * pc_y[k] * sgk1_35[k];

        t_72[k] = pb_z[k] * sgk0_0[k]
                  - f_12 * pc_z[k] * sgk1_0[k];

        t_73[k] = f_3 * pc_y[k] * shi_56[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_z, pc_y, pc_z, sgk0_3, sgk0_5, sgi_0, \
                         sgi_2, sgk1_3, sgk1_5, shi_56, shi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * sgi_0[k]
                  + f_3 * pc_z[k] * shi_56[k];

        t_75[k] = pb_z[k] * sgk0_3[k]
                  - f_12 * pc_z[k] * sgk1_3[k];

        t_76[k] = f_3 * pc_y[k] * shi_58[k];

        t_77[k] = pb_z[k] * sgk0_5[k]
                  + f_14 * sgi_2[k]
                  - f_12 * pc_z[k] * sgk1_5[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_z, pc_y, pc_z, sgk0_6, sgk0_9, sgi_3, \
                         sgi_5, sgk1_6, sgk1_9, shi_59, shi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pb_z[k] * sgk0_6[k]
                  - f_12 * pc_z[k] * sgk1_6[k];

        t_79[k] = f_13 * sgi_3[k]
                  + f_3 * pc_z[k] * shi_59[k];

        t_80[k] = f_3 * pc_y[k] * shi_61[k];

        t_81[k] = pb_z[k] * sgk0_9[k]
                  + f_15 * sgi_5[k]
                  - f_12 * pc_z[k] * sgk1_9[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pb_z, pc_y, pc_z, sgk0_10, sgk0_12, sgi_6, \
                         sgi_7, sgk1_10, sgk1_12, shi_62, shi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pb_z[k] * sgk0_10[k]
                  - f_12 * pc_z[k] * sgk1_10[k];

        t_83[k] = f_13 * sgi_6[k]
                  + f_3 * pc_z[k] * shi_62[k];

        t_84[k] = pb_z[k] * sgk0_12[k]
                  + f_14 * sgi_7[k]
                  - f_12 * pc_z[k] * sgk1_12[k];

        t_85[k] = f_3 * pc_y[k] * shi_65[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_z, pc_z, sgk0_14, sgk0_15, sgk0_17, sgi_9, \
                         sgi_10, sgi_11, sgk1_14, sgk1_15, sgk1_17, \
                         shi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_z[k] * sgk0_14[k]
                  + f_16 * sgi_9[k]
                  - f_12 * pc_z[k] * sgk1_14[k];

        t_87[k] = pb_z[k] * sgk0_15[k]
                  - f_12 * pc_z[k] * sgk1_15[k];

        t_88[k] = f_13 * sgi_10[k]
                  + f_3 * pc_z[k] * shi_66[k];

        t_89[k] = pb_z[k] * sgk0_17[k]
                  + f_14 * sgi_11[k]
                  - f_12 * pc_z[k] * sgk1_17[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_z, pc_y, pc_z, sgk0_18, sgk0_20, sgi_12, sgi_14, \
                         sgk1_18, sgk1_20, shi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_z[k] * sgk0_18[k]
                  + f_15 * sgi_12[k]
                  - f_12 * pc_z[k] * sgk1_18[k];

        t_91[k] = f_3 * pc_y[k] * shi_70[k];

        t_92[k] = pb_z[k] * sgk0_20[k]
                  + f_0 * sgi_14[k]
                  - f_12 * pc_z[k] * sgk1_20[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pc_x, sgi_77, sgi_78, sgi_79, sgi_80, \
                         sgi_81, shi_77, shi_78, shi_79, shi_80, \
                         shi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_16 * sgi_77[k]
                  + f_3 * pc_x[k] * shi_77[k];

        t_94[k] = f_16 * sgi_78[k]
                  + f_3 * pc_x[k] * shi_78[k];

        t_95[k] = f_16 * sgi_79[k]
                  + f_3 * pc_x[k] * shi_79[k];

        t_96[k] = f_16 * sgi_80[k]
                  + f_3 * pc_x[k] * shi_80[k];

        t_97[k] = f_16 * sgi_81[k]
                  + f_3 * pc_x[k] * shi_81[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_z, pc_x, pc_z, sgk0_28, sgi_21, sgi_82, \
                         sgi_83, sgk1_28, shi_77, shi_82, shi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_16 * sgi_82[k]
                  + f_3 * pc_x[k] * shi_82[k];

        t_99[k] = f_16 * sgi_83[k]
                  + f_3 * pc_x[k] * shi_83[k];

        t_100[k] = pb_z[k] * sgk0_28[k]
                   - f_12 * pc_z[k] * sgk1_28[k];

        t_101[k] = f_13 * sgi_21[k]
                   + f_3 * pc_z[k] * shi_77[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pc_y, shh0_59, shh0_60, shh0_61, shh1_59, \
                         shh1_60, shh1_61, shi_79, shi_80, shi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_4 * shh0_59[k]
                   - f_5 * shh1_59[k]
                   + f_3 * pc_y[k] * shi_79[k];

        t_103[k] = f_6 * shh0_60[k]
                   - f_7 * shh1_60[k]
                   + f_3 * pc_y[k] * shi_80[k];

        t_104[k] = f_8 * shh0_61[k]
                   - f_9 * shh1_61[k]
                   + f_3 * pc_y[k] * shi_81[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_x, pc_y, pc_z, sgi_27, sgi_84, \
                         shh0_62, shh0_63, shh1_62, shh1_63, shi_82, shi_83, \
                         shi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_10 * shh0_62[k]
                   - f_11 * shh1_62[k]
                   + f_3 * pc_y[k] * shi_82[k];

        t_106[k] = f_3 * pc_y[k] * shi_83[k];

        t_107[k] = f_13 * sgi_27[k]
                   + f_1 * shh0_62[k]
                   - f_2 * shh1_62[k]
                   + f_3 * pc_z[k] * shi_83[k];

        t_108[k] = f_15 * sgi_84[k]
                   + f_1 * shh0_63[k]
                   - f_2 * shh1_63[k]
                   + f_3 * pc_x[k] * shi_84[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, pc_z, sgi_28, sgi_30, sgi_87, \
                         shh0_66, shh1_66, shi_84, shi_86, shi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_14 * sgi_28[k]
                   + f_3 * pc_y[k] * shi_84[k];

        t_110[k] = f_3 * pc_z[k] * shi_84[k];

        t_111[k] = f_15 * sgi_87[k]
                   + f_4 * shh0_66[k]
                   - f_5 * shh1_66[k]
                   + f_3 * pc_x[k] * shi_87[k];

        t_112[k] = f_14 * sgi_30[k]
                   + f_3 * pc_y[k] * shi_86[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pc_x, pc_z, sgi_89, sgi_90, shh0_68, shh0_69, \
                         shh1_68, shh1_69, shi_87, shi_89, shi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_15 * sgi_89[k]
                   + f_4 * shh0_68[k]
                   - f_5 * shh1_68[k]
                   + f_3 * pc_x[k] * shi_89[k];

        t_114[k] = f_15 * sgi_90[k]
                   + f_6 * shh0_69[k]
                   - f_7 * shh1_69[k]
                   + f_3 * pc_x[k] * shi_90[k];

        t_115[k] = f_3 * pc_z[k] * shi_87[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pc_x, pc_y, sgi_33, sgi_93, sgi_94, shh0_72, \
                         shh0_73, shh1_72, shh1_73, shi_89, shi_93, \
                         shi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_14 * sgi_33[k]
                   + f_3 * pc_y[k] * shi_89[k];

        t_117[k] = f_15 * sgi_93[k]
                   + f_6 * shh0_72[k]
                   - f_7 * shh1_72[k]
                   + f_3 * pc_x[k] * shi_93[k];

        t_118[k] = f_15 * sgi_94[k]
                   + f_8 * shh0_73[k]
                   - f_9 * shh1_73[k]
                   + f_3 * pc_x[k] * shi_94[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_y, pc_z, sgi_37, sgi_96, shh0_75, \
                         shh1_75, shi_90, shi_93, shi_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_3 * pc_z[k] * shi_90[k];

        t_120[k] = f_15 * sgi_96[k]
                   + f_8 * shh0_75[k]
                   - f_9 * shh1_75[k]
                   + f_3 * pc_x[k] * shi_96[k];

        t_121[k] = f_14 * sgi_37[k]
                   + f_3 * pc_y[k] * shi_93[k];
    }
}

static auto
compute_prim_shk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgk0,
                                                          const size_t sgi, const size_t sgk1,
                                                          const size_t shh0, const size_t shh1,
                                                          const size_t shi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_15 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgk0_39 = buffer.data(sgk0 + 39);
    const auto *sgk0_42 = buffer.data(sgk0 + 42);
    const auto *sgk0_46 = buffer.data(sgk0 + 46);
    const auto *sgk0_51 = buffer.data(sgk0 + 51);
    const auto *sgk0_64 = buffer.data(sgk0 + 64);
    const auto *sgk0_72 = buffer.data(sgk0 + 72);
    const auto *sgk0_77 = buffer.data(sgk0 + 77);
    const auto *sgk0_81 = buffer.data(sgk0 + 81);
    const auto *sgk0_84 = buffer.data(sgk0 + 84);
    const auto *sgk0_86 = buffer.data(sgk0 + 86);
    const auto *sgk0_89 = buffer.data(sgk0 + 89);
    const auto *sgk0_90 = buffer.data(sgk0 + 90);
    const auto *sgk0_92 = buffer.data(sgk0 + 92);
    const auto *sgk0_107 = buffer.data(sgk0 + 107);

    const auto *sgi_28 = buffer.data(sgi + 28);
    const auto *sgi_31 = buffer.data(sgi + 31);
    const auto *sgi_34 = buffer.data(sgi + 34);
    const auto *sgi_38 = buffer.data(sgi + 38);
    const auto *sgi_42 = buffer.data(sgi + 42);
    const auto *sgi_49 = buffer.data(sgi + 49);
    const auto *sgi_51 = buffer.data(sgi + 51);
    const auto *sgi_52 = buffer.data(sgi + 52);
    const auto *sgi_53 = buffer.data(sgi + 53);
    const auto *sgi_54 = buffer.data(sgi + 54);
    const auto *sgi_55 = buffer.data(sgi + 55);
    const auto *sgi_56 = buffer.data(sgi + 56);
    const auto *sgi_58 = buffer.data(sgi + 58);
    const auto *sgi_59 = buffer.data(sgi + 59);
    const auto *sgi_61 = buffer.data(sgi + 61);
    const auto *sgi_62 = buffer.data(sgi + 62);
    const auto *sgi_64 = buffer.data(sgi + 64);
    const auto *sgi_65 = buffer.data(sgi + 65);
    const auto *sgi_66 = buffer.data(sgi + 66);
    const auto *sgi_68 = buffer.data(sgi + 68);
    const auto *sgi_69 = buffer.data(sgi + 69);
    const auto *sgi_70 = buffer.data(sgi + 70);
    const auto *sgi_77 = buffer.data(sgi + 77);
    const auto *sgi_79 = buffer.data(sgi + 79);
    const auto *sgi_80 = buffer.data(sgi + 80);
    const auto *sgi_81 = buffer.data(sgi + 81);
    const auto *sgi_82 = buffer.data(sgi + 82);
    const auto *sgi_83 = buffer.data(sgi + 83);
    const auto *sgi_84 = buffer.data(sgi + 84);
    const auto *sgi_86 = buffer.data(sgi + 86);
    const auto *sgi_89 = buffer.data(sgi + 89);
    const auto *sgi_93 = buffer.data(sgi + 93);
    const auto *sgi_98 = buffer.data(sgi + 98);
    const auto *sgi_99 = buffer.data(sgi + 99);
    const auto *sgi_101 = buffer.data(sgi + 101);
    const auto *sgi_102 = buffer.data(sgi + 102);
    const auto *sgi_104 = buffer.data(sgi + 104);
    const auto *sgi_105 = buffer.data(sgi + 105);
    const auto *sgi_106 = buffer.data(sgi + 106);
    const auto *sgi_107 = buffer.data(sgi + 107);
    const auto *sgi_108 = buffer.data(sgi + 108);
    const auto *sgi_109 = buffer.data(sgi + 109);
    const auto *sgi_110 = buffer.data(sgi + 110);
    const auto *sgi_111 = buffer.data(sgi + 111);
    const auto *sgi_133 = buffer.data(sgi + 133);
    const auto *sgi_134 = buffer.data(sgi + 134);
    const auto *sgi_135 = buffer.data(sgi + 135);
    const auto *sgi_136 = buffer.data(sgi + 136);
    const auto *sgi_137 = buffer.data(sgi + 137);
    const auto *sgi_138 = buffer.data(sgi + 138);
    const auto *sgi_139 = buffer.data(sgi + 139);
    const auto *sgi_140 = buffer.data(sgi + 140);
    const auto *sgi_143 = buffer.data(sgi + 143);
    const auto *sgi_145 = buffer.data(sgi + 145);
    const auto *sgi_146 = buffer.data(sgi + 146);
    const auto *sgi_149 = buffer.data(sgi + 149);
    const auto *sgi_150 = buffer.data(sgi + 150);
    const auto *sgi_152 = buffer.data(sgi + 152);
    const auto *sgi_154 = buffer.data(sgi + 154);
    const auto *sgi_155 = buffer.data(sgi + 155);
    const auto *sgi_157 = buffer.data(sgi + 157);
    const auto *sgi_158 = buffer.data(sgi + 158);
    const auto *sgi_160 = buffer.data(sgi + 160);
    const auto *sgi_161 = buffer.data(sgi + 161);
    const auto *sgi_162 = buffer.data(sgi + 162);
    const auto *sgi_163 = buffer.data(sgi + 163);
    const auto *sgi_164 = buffer.data(sgi + 164);
    const auto *sgi_165 = buffer.data(sgi + 165);
    const auto *sgi_166 = buffer.data(sgi + 166);
    const auto *sgi_167 = buffer.data(sgi + 167);
    const auto *sgi_168 = buffer.data(sgi + 168);
    const auto *sgi_171 = buffer.data(sgi + 171);
    const auto *sgi_173 = buffer.data(sgi + 173);
    const auto *sgi_174 = buffer.data(sgi + 174);
    const auto *sgi_177 = buffer.data(sgi + 177);
    const auto *sgi_178 = buffer.data(sgi + 178);
    const auto *sgi_180 = buffer.data(sgi + 180);
    const auto *sgi_182 = buffer.data(sgi + 182);
    const auto *sgi_183 = buffer.data(sgi + 183);
    const auto *sgi_185 = buffer.data(sgi + 185);
    const auto *sgi_186 = buffer.data(sgi + 186);

    const auto *sgk1_39 = buffer.data(sgk1 + 39);
    const auto *sgk1_42 = buffer.data(sgk1 + 42);
    const auto *sgk1_46 = buffer.data(sgk1 + 46);
    const auto *sgk1_51 = buffer.data(sgk1 + 51);
    const auto *sgk1_64 = buffer.data(sgk1 + 64);
    const auto *sgk1_72 = buffer.data(sgk1 + 72);
    const auto *sgk1_77 = buffer.data(sgk1 + 77);
    const auto *sgk1_81 = buffer.data(sgk1 + 81);
    const auto *sgk1_84 = buffer.data(sgk1 + 84);
    const auto *sgk1_86 = buffer.data(sgk1 + 86);
    const auto *sgk1_89 = buffer.data(sgk1 + 89);
    const auto *sgk1_90 = buffer.data(sgk1 + 90);
    const auto *sgk1_92 = buffer.data(sgk1 + 92);
    const auto *sgk1_107 = buffer.data(sgk1 + 107);

    const auto *shh0_77 = buffer.data(shh0 + 77);
    const auto *shh0_78 = buffer.data(shh0 + 78);
    const auto *shh0_80 = buffer.data(shh0 + 80);
    const auto *shh0_81 = buffer.data(shh0 + 81);
    const auto *shh0_82 = buffer.data(shh0 + 82);
    const auto *shh0_83 = buffer.data(shh0 + 83);
    const auto *shh0_101 = buffer.data(shh0 + 101);
    const auto *shh0_102 = buffer.data(shh0 + 102);
    const auto *shh0_103 = buffer.data(shh0 + 103);
    const auto *shh0_104 = buffer.data(shh0 + 104);
    const auto *shh0_105 = buffer.data(shh0 + 105);
    const auto *shh0_108 = buffer.data(shh0 + 108);
    const auto *shh0_110 = buffer.data(shh0 + 110);
    const auto *shh0_111 = buffer.data(shh0 + 111);
    const auto *shh0_114 = buffer.data(shh0 + 114);
    const auto *shh0_115 = buffer.data(shh0 + 115);
    const auto *shh0_117 = buffer.data(shh0 + 117);
    const auto *shh0_119 = buffer.data(shh0 + 119);
    const auto *shh0_120 = buffer.data(shh0 + 120);
    const auto *shh0_122 = buffer.data(shh0 + 122);
    const auto *shh0_123 = buffer.data(shh0 + 123);
    const auto *shh0_124 = buffer.data(shh0 + 124);
    const auto *shh0_125 = buffer.data(shh0 + 125);
    const auto *shh0_126 = buffer.data(shh0 + 126);
    const auto *shh0_129 = buffer.data(shh0 + 129);
    const auto *shh0_131 = buffer.data(shh0 + 131);
    const auto *shh0_132 = buffer.data(shh0 + 132);
    const auto *shh0_135 = buffer.data(shh0 + 135);
    const auto *shh0_136 = buffer.data(shh0 + 136);
    const auto *shh0_138 = buffer.data(shh0 + 138);
    const auto *shh0_140 = buffer.data(shh0 + 140);
    const auto *shh0_141 = buffer.data(shh0 + 141);
    const auto *shh0_143 = buffer.data(shh0 + 143);
    const auto *shh0_144 = buffer.data(shh0 + 144);

    const auto *shh1_77 = buffer.data(shh1 + 77);
    const auto *shh1_78 = buffer.data(shh1 + 78);
    const auto *shh1_80 = buffer.data(shh1 + 80);
    const auto *shh1_81 = buffer.data(shh1 + 81);
    const auto *shh1_82 = buffer.data(shh1 + 82);
    const auto *shh1_83 = buffer.data(shh1 + 83);
    const auto *shh1_101 = buffer.data(shh1 + 101);
    const auto *shh1_102 = buffer.data(shh1 + 102);
    const auto *shh1_103 = buffer.data(shh1 + 103);
    const auto *shh1_104 = buffer.data(shh1 + 104);
    const auto *shh1_105 = buffer.data(shh1 + 105);
    const auto *shh1_108 = buffer.data(shh1 + 108);
    const auto *shh1_110 = buffer.data(shh1 + 110);
    const auto *shh1_111 = buffer.data(shh1 + 111);
    const auto *shh1_114 = buffer.data(shh1 + 114);
    const auto *shh1_115 = buffer.data(shh1 + 115);
    const auto *shh1_117 = buffer.data(shh1 + 117);
    const auto *shh1_119 = buffer.data(shh1 + 119);
    const auto *shh1_120 = buffer.data(shh1 + 120);
    const auto *shh1_122 = buffer.data(shh1 + 122);
    const auto *shh1_123 = buffer.data(shh1 + 123);
    const auto *shh1_124 = buffer.data(shh1 + 124);
    const auto *shh1_125 = buffer.data(shh1 + 125);
    const auto *shh1_126 = buffer.data(shh1 + 126);
    const auto *shh1_129 = buffer.data(shh1 + 129);
    const auto *shh1_131 = buffer.data(shh1 + 131);
    const auto *shh1_132 = buffer.data(shh1 + 132);
    const auto *shh1_135 = buffer.data(shh1 + 135);
    const auto *shh1_136 = buffer.data(shh1 + 136);
    const auto *shh1_138 = buffer.data(shh1 + 138);
    const auto *shh1_140 = buffer.data(shh1 + 140);
    const auto *shh1_141 = buffer.data(shh1 + 141);
    const auto *shh1_143 = buffer.data(shh1 + 143);
    const auto *shh1_144 = buffer.data(shh1 + 144);

    const auto *shi_94 = buffer.data(shi + 94);
    const auto *shi_98 = buffer.data(shi + 98);
    const auto *shi_99 = buffer.data(shi + 99);
    const auto *shi_101 = buffer.data(shi + 101);
    const auto *shi_102 = buffer.data(shi + 102);
    const auto *shi_104 = buffer.data(shi + 104);
    const auto *shi_105 = buffer.data(shi + 105);
    const auto *shi_106 = buffer.data(shi + 106);
    const auto *shi_107 = buffer.data(shi + 107);
    const auto *shi_108 = buffer.data(shi + 108);
    const auto *shi_109 = buffer.data(shi + 109);
    const auto *shi_110 = buffer.data(shi + 110);
    const auto *shi_111 = buffer.data(shi + 111);
    const auto *shi_112 = buffer.data(shi + 112);
    const auto *shi_114 = buffer.data(shi + 114);
    const auto *shi_115 = buffer.data(shi + 115);
    const auto *shi_117 = buffer.data(shi + 117);
    const auto *shi_118 = buffer.data(shi + 118);
    const auto *shi_121 = buffer.data(shi + 121);
    const auto *shi_122 = buffer.data(shi + 122);
    const auto *shi_126 = buffer.data(shi + 126);
    const auto *shi_133 = buffer.data(shi + 133);
    const auto *shi_134 = buffer.data(shi + 134);
    const auto *shi_135 = buffer.data(shi + 135);
    const auto *shi_136 = buffer.data(shi + 136);
    const auto *shi_137 = buffer.data(shi + 137);
    const auto *shi_138 = buffer.data(shi + 138);
    const auto *shi_139 = buffer.data(shi + 139);
    const auto *shi_140 = buffer.data(shi + 140);
    const auto *shi_142 = buffer.data(shi + 142);
    const auto *shi_143 = buffer.data(shi + 143);
    const auto *shi_145 = buffer.data(shi + 145);
    const auto *shi_146 = buffer.data(shi + 146);
    const auto *shi_149 = buffer.data(shi + 149);
    const auto *shi_150 = buffer.data(shi + 150);
    const auto *shi_152 = buffer.data(shi + 152);
    const auto *shi_154 = buffer.data(shi + 154);
    const auto *shi_155 = buffer.data(shi + 155);
    const auto *shi_157 = buffer.data(shi + 157);
    const auto *shi_158 = buffer.data(shi + 158);
    const auto *shi_160 = buffer.data(shi + 160);
    const auto *shi_161 = buffer.data(shi + 161);
    const auto *shi_162 = buffer.data(shi + 162);
    const auto *shi_163 = buffer.data(shi + 163);
    const auto *shi_164 = buffer.data(shi + 164);
    const auto *shi_165 = buffer.data(shi + 165);
    const auto *shi_166 = buffer.data(shi + 166);
    const auto *shi_167 = buffer.data(shi + 167);
    const auto *shi_168 = buffer.data(shi + 168);
    const auto *shi_170 = buffer.data(shi + 170);
    const auto *shi_171 = buffer.data(shi + 171);
    const auto *shi_173 = buffer.data(shi + 173);
    const auto *shi_174 = buffer.data(shi + 174);
    const auto *shi_177 = buffer.data(shi + 177);
    const auto *shi_178 = buffer.data(shi + 178);
    const auto *shi_180 = buffer.data(shi + 180);
    const auto *shi_182 = buffer.data(shi + 182);
    const auto *shi_183 = buffer.data(shi + 183);
    const auto *shi_185 = buffer.data(shi + 185);
    const auto *shi_186 = buffer.data(shi + 186);

#pragma omp simd aligned(t_122, t_123, t_124, pc_x, pc_z, sgi_98, sgi_99, shh0_77, shh0_78, \
                         shh1_77, shh1_78, shi_94, shi_98, shi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_15 * sgi_98[k]
                   + f_8 * shh0_77[k]
                   - f_9 * shh1_77[k]
                   + f_3 * pc_x[k] * shi_98[k];

        t_123[k] = f_15 * sgi_99[k]
                   + f_10 * shh0_78[k]
                   - f_11 * shh1_78[k]
                   + f_3 * pc_x[k] * shi_99[k];

        t_124[k] = f_3 * pc_z[k] * shi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pc_x, pc_y, sgi_42, sgi_101, sgi_102, shh0_80, \
                         shh0_81, shh1_80, shh1_81, shi_98, shi_101, \
                         shi_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_15 * sgi_101[k]
                   + f_10 * shh0_80[k]
                   - f_11 * shh1_80[k]
                   + f_3 * pc_x[k] * shi_101[k];

        t_126[k] = f_15 * sgi_102[k]
                   + f_10 * shh0_81[k]
                   - f_11 * shh1_81[k]
                   + f_3 * pc_x[k] * shi_102[k];

        t_127[k] = f_14 * sgi_42[k]
                   + f_3 * pc_y[k] * shi_98[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_x, sgi_104, sgi_105, sgi_106, sgi_107, \
                         shh0_83, shh1_83, shi_104, shi_105, shi_106, \
                         shi_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_15 * sgi_104[k]
                   + f_10 * shh0_83[k]
                   - f_11 * shh1_83[k]
                   + f_3 * pc_x[k] * shi_104[k];

        t_129[k] = f_15 * sgi_105[k]
                   + f_3 * pc_x[k] * shi_105[k];

        t_130[k] = f_15 * sgi_106[k]
                   + f_3 * pc_x[k] * shi_106[k];

        t_131[k] = f_15 * sgi_107[k]
                   + f_3 * pc_x[k] * shi_107[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, sgi_108, sgi_109, sgi_110, sgi_111, \
                         shi_108, shi_109, shi_110, shi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_15 * sgi_108[k]
                   + f_3 * pc_x[k] * shi_108[k];

        t_133[k] = f_15 * sgi_109[k]
                   + f_3 * pc_x[k] * shi_109[k];

        t_134[k] = f_15 * sgi_110[k]
                   + f_3 * pc_x[k] * shi_110[k];

        t_135[k] = f_15 * sgi_111[k]
                   + f_3 * pc_x[k] * shi_111[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, pc_z, sgi_49, sgi_51, shh0_78, shh0_80, \
                         shh1_78, shh1_80, shi_105, shi_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_14 * sgi_49[k]
                   + f_1 * shh0_78[k]
                   - f_2 * shh1_78[k]
                   + f_3 * pc_y[k] * shi_105[k];

        t_137[k] = f_3 * pc_z[k] * shi_105[k];

        t_138[k] = f_14 * sgi_51[k]
                   + f_4 * shh0_80[k]
                   - f_5 * shh1_80[k]
                   + f_3 * pc_y[k] * shi_107[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pc_y, sgi_52, sgi_53, sgi_54, shh0_81, shh0_82, \
                         shh0_83, shh1_81, shh1_82, shh1_83, shi_108, shi_109, \
                         shi_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_14 * sgi_52[k]
                   + f_6 * shh0_81[k]
                   - f_7 * shh1_81[k]
                   + f_3 * pc_y[k] * shi_108[k];

        t_140[k] = f_14 * sgi_53[k]
                   + f_8 * shh0_82[k]
                   - f_9 * shh1_82[k]
                   + f_3 * pc_y[k] * shi_109[k];

        t_141[k] = f_14 * sgi_54[k]
                   + f_10 * shh0_83[k]
                   - f_11 * shh1_83[k]
                   + f_3 * pc_y[k] * shi_110[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_y, pc_y, pc_z, sgk0_72, sgi_55, \
                         sgi_56, sgk1_72, shh0_83, shh1_83, shi_111, \
                         shi_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_14 * sgi_55[k]
                   + f_3 * pc_y[k] * shi_111[k];

        t_143[k] = f_1 * shh0_83[k]
                   - f_2 * shh1_83[k]
                   + f_3 * pc_z[k] * shi_111[k];

        t_144[k] = pb_y[k] * sgk0_72[k]
                   - f_12 * pc_y[k] * sgk1_72[k];

        t_145[k] = f_13 * sgi_56[k]
                   + f_3 * pc_y[k] * shi_112[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pb_y, pb_z, pc_y, pc_z, sgk0_39, sgk0_77, \
                         sgi_28, sgi_58, sgk1_39, sgk1_77, shi_112, \
                         shi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_13 * sgi_28[k]
                   + f_3 * pc_z[k] * shi_112[k];

        t_147[k] = pb_z[k] * sgk0_39[k]
                   - f_12 * pc_z[k] * sgk1_39[k];

        t_148[k] = f_13 * sgi_58[k]
                   + f_3 * pc_y[k] * shi_114[k];

        t_149[k] = pb_y[k] * sgk0_77[k]
                   - f_12 * pc_y[k] * sgk1_77[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_y, pb_z, pc_y, pc_z, sgk0_42, sgk0_81, \
                         sgi_31, sgi_61, sgk1_42, sgk1_81, shi_115, \
                         shi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pb_z[k] * sgk0_42[k]
                   - f_12 * pc_z[k] * sgk1_42[k];

        t_151[k] = f_13 * sgi_31[k]
                   + f_3 * pc_z[k] * shi_115[k];

        t_152[k] = f_13 * sgi_61[k]
                   + f_3 * pc_y[k] * shi_117[k];

        t_153[k] = pb_y[k] * sgk0_81[k]
                   - f_12 * pc_y[k] * sgk1_81[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_y, pb_z, pc_y, pc_z, sgk0_46, sgk0_84, \
                         sgi_34, sgi_64, sgk1_46, sgk1_84, shi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_z[k] * sgk0_46[k]
                   - f_12 * pc_z[k] * sgk1_46[k];

        t_155[k] = f_13 * sgi_34[k]
                   + f_3 * pc_z[k] * shi_118[k];

        t_156[k] = pb_y[k] * sgk0_84[k]
                   + f_14 * sgi_64[k]
                   - f_12 * pc_y[k] * sgk1_84[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pb_y, pb_z, pc_y, pc_z, sgk0_51, sgk0_86, \
                         sgi_38, sgi_65, sgk1_51, sgk1_86, shi_121, \
                         shi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_13 * sgi_65[k]
                   + f_3 * pc_y[k] * shi_121[k];

        t_158[k] = pb_y[k] * sgk0_86[k]
                   - f_12 * pc_y[k] * sgk1_86[k];

        t_159[k] = pb_z[k] * sgk0_51[k]
                   - f_12 * pc_z[k] * sgk1_51[k];

        t_160[k] = f_13 * sgi_38[k]
                   + f_3 * pc_z[k] * shi_122[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_y, pc_y, sgk0_89, sgk0_90, sgk0_92, \
                         sgi_68, sgi_69, sgi_70, sgk1_89, sgk1_90, sgk1_92, \
                         shi_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_y[k] * sgk0_89[k]
                   + f_15 * sgi_68[k]
                   - f_12 * pc_y[k] * sgk1_89[k];

        t_162[k] = pb_y[k] * sgk0_90[k]
                   + f_14 * sgi_69[k]
                   - f_12 * pc_y[k] * sgk1_90[k];

        t_163[k] = f_13 * sgi_70[k]
                   + f_3 * pc_y[k] * shi_126[k];

        t_164[k] = pb_y[k] * sgk0_92[k]
                   - f_12 * pc_y[k] * sgk1_92[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pc_x, sgi_133, sgi_134, sgi_135, \
                         sgi_136, sgi_137, shi_133, shi_134, shi_135, shi_136, \
                         shi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_15 * sgi_133[k]
                   + f_3 * pc_x[k] * shi_133[k];

        t_166[k] = f_15 * sgi_134[k]
                   + f_3 * pc_x[k] * shi_134[k];

        t_167[k] = f_15 * sgi_135[k]
                   + f_3 * pc_x[k] * shi_135[k];

        t_168[k] = f_15 * sgi_136[k]
                   + f_3 * pc_x[k] * shi_136[k];

        t_169[k] = f_15 * sgi_137[k]
                   + f_3 * pc_x[k] * shi_137[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_z, pc_x, pc_z, sgk0_64, sgi_49, \
                         sgi_138, sgi_139, sgk1_64, shi_133, shi_138, \
                         shi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_15 * sgi_138[k]
                   + f_3 * pc_x[k] * shi_138[k];

        t_171[k] = f_15 * sgi_139[k]
                   + f_3 * pc_x[k] * shi_139[k];

        t_172[k] = pb_z[k] * sgk0_64[k]
                   - f_12 * pc_z[k] * sgk1_64[k];

        t_173[k] = f_13 * sgi_49[k]
                   + f_3 * pc_z[k] * shi_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, sgi_79, sgi_80, sgi_81, shh0_101, \
                         shh0_102, shh0_103, shh1_101, shh1_102, shh1_103, shi_135, shi_136, \
                         shi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * sgi_79[k]
                   + f_4 * shh0_101[k]
                   - f_5 * shh1_101[k]
                   + f_3 * pc_y[k] * shi_135[k];

        t_175[k] = f_13 * sgi_80[k]
                   + f_6 * shh0_102[k]
                   - f_7 * shh1_102[k]
                   + f_3 * pc_y[k] * shi_136[k];

        t_176[k] = f_13 * sgi_81[k]
                   + f_8 * shh0_103[k]
                   - f_9 * shh1_103[k]
                   + f_3 * pc_y[k] * shi_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_y, pc_y, sgk0_107, sgi_82, sgi_83, sgk1_107, \
                         shh0_104, shh1_104, shi_138, shi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * sgi_82[k]
                   + f_10 * shh0_104[k]
                   - f_11 * shh1_104[k]
                   + f_3 * pc_y[k] * shi_138[k];

        t_178[k] = f_13 * sgi_83[k]
                   + f_3 * pc_y[k] * shi_139[k];

        t_179[k] = pb_y[k] * sgk0_107[k]
                   - f_12 * pc_y[k] * sgk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, sgi_56, sgi_140, \
                         sgi_143, shh0_105, shh0_108, shh1_105, shh1_108, shi_140, \
                         shi_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_15 * sgi_140[k]
                   + f_1 * shh0_105[k]
                   - f_2 * shh1_105[k]
                   + f_3 * pc_x[k] * shi_140[k];

        t_181[k] = f_3 * pc_y[k] * shi_140[k];

        t_182[k] = f_14 * sgi_56[k]
                   + f_3 * pc_z[k] * shi_140[k];

        t_183[k] = f_15 * sgi_143[k]
                   + f_4 * shh0_108[k]
                   - f_5 * shh1_108[k]
                   + f_3 * pc_x[k] * shi_143[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pc_x, pc_y, sgi_145, sgi_146, shh0_110, \
                         shh0_111, shh1_110, shh1_111, shi_142, shi_145, \
                         shi_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_y[k] * shi_142[k];

        t_185[k] = f_15 * sgi_145[k]
                   + f_4 * shh0_110[k]
                   - f_5 * shh1_110[k]
                   + f_3 * pc_x[k] * shi_145[k];

        t_186[k] = f_15 * sgi_146[k]
                   + f_6 * shh0_111[k]
                   - f_7 * shh1_111[k]
                   + f_3 * pc_x[k] * shi_146[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pc_x, pc_y, pc_z, sgi_59, sgi_149, shh0_114, \
                         shh1_114, shi_143, shi_145, shi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_14 * sgi_59[k]
                   + f_3 * pc_z[k] * shi_143[k];

        t_188[k] = f_3 * pc_y[k] * shi_145[k];

        t_189[k] = f_15 * sgi_149[k]
                   + f_6 * shh0_114[k]
                   - f_7 * shh1_114[k]
                   + f_3 * pc_x[k] * shi_149[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_x, pc_z, sgi_62, sgi_150, sgi_152, shh0_115, \
                         shh0_117, shh1_115, shh1_117, shi_146, shi_150, \
                         shi_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_15 * sgi_150[k]
                   + f_8 * shh0_115[k]
                   - f_9 * shh1_115[k]
                   + f_3 * pc_x[k] * shi_150[k];

        t_191[k] = f_14 * sgi_62[k]
                   + f_3 * pc_z[k] * shi_146[k];

        t_192[k] = f_15 * sgi_152[k]
                   + f_8 * shh0_117[k]
                   - f_9 * shh1_117[k]
                   + f_3 * pc_x[k] * shi_152[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, sgi_154, sgi_155, shh0_119, \
                         shh0_120, shh1_119, shh1_120, shi_149, shi_154, \
                         shi_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * shi_149[k];

        t_194[k] = f_15 * sgi_154[k]
                   + f_8 * shh0_119[k]
                   - f_9 * shh1_119[k]
                   + f_3 * pc_x[k] * shi_154[k];

        t_195[k] = f_15 * sgi_155[k]
                   + f_10 * shh0_120[k]
                   - f_11 * shh1_120[k]
                   + f_3 * pc_x[k] * shi_155[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, pc_x, pc_z, sgi_66, sgi_157, sgi_158, shh0_122, \
                         shh0_123, shh1_122, shh1_123, shi_150, shi_157, \
                         shi_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_14 * sgi_66[k]
                   + f_3 * pc_z[k] * shi_150[k];

        t_197[k] = f_15 * sgi_157[k]
                   + f_10 * shh0_122[k]
                   - f_11 * shh1_122[k]
                   + f_3 * pc_x[k] * shi_157[k];

        t_198[k] = f_15 * sgi_158[k]
                   + f_10 * shh0_123[k]
                   - f_11 * shh1_123[k]
                   + f_3 * pc_x[k] * shi_158[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pc_x, pc_y, sgi_160, sgi_161, sgi_162, \
                         shh0_125, shh1_125, shi_154, shi_160, shi_161, \
                         shi_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_3 * pc_y[k] * shi_154[k];

        t_200[k] = f_15 * sgi_160[k]
                   + f_10 * shh0_125[k]
                   - f_11 * shh1_125[k]
                   + f_3 * pc_x[k] * shi_160[k];

        t_201[k] = f_15 * sgi_161[k]
                   + f_3 * pc_x[k] * shi_161[k];

        t_202[k] = f_15 * sgi_162[k]
                   + f_3 * pc_x[k] * shi_162[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pc_x, sgi_163, sgi_164, sgi_165, \
                         sgi_166, sgi_167, shi_163, shi_164, shi_165, shi_166, \
                         shi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_15 * sgi_163[k]
                   + f_3 * pc_x[k] * shi_163[k];

        t_204[k] = f_15 * sgi_164[k]
                   + f_3 * pc_x[k] * shi_164[k];

        t_205[k] = f_15 * sgi_165[k]
                   + f_3 * pc_x[k] * shi_165[k];

        t_206[k] = f_15 * sgi_166[k]
                   + f_3 * pc_x[k] * shi_166[k];

        t_207[k] = f_15 * sgi_167[k]
                   + f_3 * pc_x[k] * shi_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_y, pc_z, sgi_77, shh0_120, shh0_122, \
                         shh0_123, shh1_120, shh1_122, shh1_123, shi_161, shi_163, \
                         shi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * shh0_120[k]
                   - f_2 * shh1_120[k]
                   + f_3 * pc_y[k] * shi_161[k];

        t_209[k] = f_14 * sgi_77[k]
                   + f_3 * pc_z[k] * shi_161[k];

        t_210[k] = f_4 * shh0_122[k]
                   - f_5 * shh1_122[k]
                   + f_3 * pc_y[k] * shi_163[k];

        t_211[k] = f_6 * shh0_123[k]
                   - f_7 * shh1_123[k]
                   + f_3 * pc_y[k] * shi_164[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_y, pc_z, sgi_83, shh0_124, shh0_125, \
                         shh1_124, shh1_125, shi_165, shi_166, \
                         shi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_8 * shh0_124[k]
                   - f_9 * shh1_124[k]
                   + f_3 * pc_y[k] * shi_165[k];

        t_213[k] = f_10 * shh0_125[k]
                   - f_11 * shh1_125[k]
                   + f_3 * pc_y[k] * shi_166[k];

        t_214[k] = f_3 * pc_y[k] * shi_167[k];

        t_215[k] = f_14 * sgi_83[k]
                   + f_1 * shh0_125[k]
                   - f_2 * shh1_125[k]
                   + f_3 * pc_z[k] * shi_167[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pc_x, pc_y, pc_z, sgi_84, sgi_168, \
                         sgi_171, shh0_126, shh0_129, shh1_126, shh1_129, shi_168, \
                         shi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_14 * sgi_168[k]
                   + f_1 * shh0_126[k]
                   - f_2 * shh1_126[k]
                   + f_3 * pc_x[k] * shi_168[k];

        t_217[k] = f_15 * sgi_84[k]
                   + f_3 * pc_y[k] * shi_168[k];

        t_218[k] = f_3 * pc_z[k] * shi_168[k];

        t_219[k] = f_14 * sgi_171[k]
                   + f_4 * shh0_129[k]
                   - f_5 * shh1_129[k]
                   + f_3 * pc_x[k] * shi_171[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_x, pc_y, sgi_86, sgi_173, sgi_174, shh0_131, \
                         shh0_132, shh1_131, shh1_132, shi_170, shi_173, \
                         shi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * sgi_86[k]
                   + f_3 * pc_y[k] * shi_170[k];

        t_221[k] = f_14 * sgi_173[k]
                   + f_4 * shh0_131[k]
                   - f_5 * shh1_131[k]
                   + f_3 * pc_x[k] * shi_173[k];

        t_222[k] = f_14 * sgi_174[k]
                   + f_6 * shh0_132[k]
                   - f_7 * shh1_132[k]
                   + f_3 * pc_x[k] * shi_174[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, pc_x, pc_y, pc_z, sgi_89, sgi_177, shh0_135, \
                         shh1_135, shi_171, shi_173, shi_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_z[k] * shi_171[k];

        t_224[k] = f_15 * sgi_89[k]
                   + f_3 * pc_y[k] * shi_173[k];

        t_225[k] = f_14 * sgi_177[k]
                   + f_6 * shh0_135[k]
                   - f_7 * shh1_135[k]
                   + f_3 * pc_x[k] * shi_177[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_x, pc_z, sgi_178, sgi_180, shh0_136, \
                         shh0_138, shh1_136, shh1_138, shi_174, shi_178, \
                         shi_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_14 * sgi_178[k]
                   + f_8 * shh0_136[k]
                   - f_9 * shh1_136[k]
                   + f_3 * pc_x[k] * shi_178[k];

        t_227[k] = f_3 * pc_z[k] * shi_174[k];

        t_228[k] = f_14 * sgi_180[k]
                   + f_8 * shh0_138[k]
                   - f_9 * shh1_138[k]
                   + f_3 * pc_x[k] * shi_180[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pc_x, pc_y, sgi_93, sgi_182, sgi_183, shh0_140, \
                         shh0_141, shh1_140, shh1_141, shi_177, shi_182, \
                         shi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_15 * sgi_93[k]
                   + f_3 * pc_y[k] * shi_177[k];

        t_230[k] = f_14 * sgi_182[k]
                   + f_8 * shh0_140[k]
                   - f_9 * shh1_140[k]
                   + f_3 * pc_x[k] * shi_182[k];

        t_231[k] = f_14 * sgi_183[k]
                   + f_10 * shh0_141[k]
                   - f_11 * shh1_141[k]
                   + f_3 * pc_x[k] * shi_183[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pc_x, pc_z, sgi_185, sgi_186, shh0_143, \
                         shh0_144, shh1_143, shh1_144, shi_178, shi_185, \
                         shi_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_3 * pc_z[k] * shi_178[k];

        t_233[k] = f_14 * sgi_185[k]
                   + f_10 * shh0_143[k]
                   - f_11 * shh1_143[k]
                   + f_3 * pc_x[k] * shi_185[k];

        t_234[k] = f_14 * sgi_186[k]
                   + f_10 * shh0_144[k]
                   - f_11 * shh1_144[k]
                   + f_3 * pc_x[k] * shi_186[k];
    }
}

static auto
compute_prim_shk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgk0,
                                                          const size_t sgi, const size_t sgk1,
                                                          const size_t shh0, const size_t shh1,
                                                          const size_t shi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgk0_108 = buffer.data(sgk0 + 108);
    const auto *sgk0_111 = buffer.data(sgk0 + 111);
    const auto *sgk0_114 = buffer.data(sgk0 + 114);
    const auto *sgk0_118 = buffer.data(sgk0 + 118);
    const auto *sgk0_120 = buffer.data(sgk0 + 120);
    const auto *sgk0_123 = buffer.data(sgk0 + 123);
    const auto *sgk0_125 = buffer.data(sgk0 + 125);
    const auto *sgk0_126 = buffer.data(sgk0 + 126);
    const auto *sgk0_136 = buffer.data(sgk0 + 136);
    const auto *sgk0_180 = buffer.data(sgk0 + 180);
    const auto *sgk0_183 = buffer.data(sgk0 + 183);
    const auto *sgk0_185 = buffer.data(sgk0 + 185);
    const auto *sgk0_186 = buffer.data(sgk0 + 186);
    const auto *sgk0_189 = buffer.data(sgk0 + 189);
    const auto *sgk0_190 = buffer.data(sgk0 + 190);
    const auto *sgk0_192 = buffer.data(sgk0 + 192);
    const auto *sgk0_194 = buffer.data(sgk0 + 194);
    const auto *sgk0_195 = buffer.data(sgk0 + 195);
    const auto *sgk0_197 = buffer.data(sgk0 + 197);
    const auto *sgk0_198 = buffer.data(sgk0 + 198);
    const auto *sgk0_200 = buffer.data(sgk0 + 200);
    const auto *sgk0_215 = buffer.data(sgk0 + 215);

    const auto *sgi_84 = buffer.data(sgi + 84);
    const auto *sgi_87 = buffer.data(sgi + 87);
    const auto *sgi_90 = buffer.data(sgi + 90);
    const auto *sgi_91 = buffer.data(sgi + 91);
    const auto *sgi_94 = buffer.data(sgi + 94);
    const auto *sgi_95 = buffer.data(sgi + 95);
    const auto *sgi_96 = buffer.data(sgi + 96);
    const auto *sgi_98 = buffer.data(sgi + 98);
    const auto *sgi_105 = buffer.data(sgi + 105);
    const auto *sgi_107 = buffer.data(sgi + 107);
    const auto *sgi_108 = buffer.data(sgi + 108);
    const auto *sgi_109 = buffer.data(sgi + 109);
    const auto *sgi_110 = buffer.data(sgi + 110);
    const auto *sgi_111 = buffer.data(sgi + 111);
    const auto *sgi_112 = buffer.data(sgi + 112);
    const auto *sgi_114 = buffer.data(sgi + 114);
    const auto *sgi_115 = buffer.data(sgi + 115);
    const auto *sgi_117 = buffer.data(sgi + 117);
    const auto *sgi_118 = buffer.data(sgi + 118);
    const auto *sgi_121 = buffer.data(sgi + 121);
    const auto *sgi_122 = buffer.data(sgi + 122);
    const auto *sgi_126 = buffer.data(sgi + 126);
    const auto *sgi_133 = buffer.data(sgi + 133);
    const auto *sgi_135 = buffer.data(sgi + 135);
    const auto *sgi_136 = buffer.data(sgi + 136);
    const auto *sgi_137 = buffer.data(sgi + 137);
    const auto *sgi_138 = buffer.data(sgi + 138);
    const auto *sgi_139 = buffer.data(sgi + 139);
    const auto *sgi_140 = buffer.data(sgi + 140);
    const auto *sgi_141 = buffer.data(sgi + 141);
    const auto *sgi_142 = buffer.data(sgi + 142);
    const auto *sgi_143 = buffer.data(sgi + 143);
    const auto *sgi_145 = buffer.data(sgi + 145);
    const auto *sgi_146 = buffer.data(sgi + 146);
    const auto *sgi_148 = buffer.data(sgi + 148);
    const auto *sgi_149 = buffer.data(sgi + 149);
    const auto *sgi_150 = buffer.data(sgi + 150);
    const auto *sgi_152 = buffer.data(sgi + 152);
    const auto *sgi_153 = buffer.data(sgi + 153);
    const auto *sgi_154 = buffer.data(sgi + 154);
    const auto *sgi_161 = buffer.data(sgi + 161);
    const auto *sgi_163 = buffer.data(sgi + 163);
    const auto *sgi_164 = buffer.data(sgi + 164);
    const auto *sgi_165 = buffer.data(sgi + 165);
    const auto *sgi_166 = buffer.data(sgi + 166);
    const auto *sgi_167 = buffer.data(sgi + 167);
    const auto *sgi_188 = buffer.data(sgi + 188);
    const auto *sgi_189 = buffer.data(sgi + 189);
    const auto *sgi_190 = buffer.data(sgi + 190);
    const auto *sgi_191 = buffer.data(sgi + 191);
    const auto *sgi_192 = buffer.data(sgi + 192);
    const auto *sgi_193 = buffer.data(sgi + 193);
    const auto *sgi_194 = buffer.data(sgi + 194);
    const auto *sgi_195 = buffer.data(sgi + 195);
    const auto *sgi_201 = buffer.data(sgi + 201);
    const auto *sgi_205 = buffer.data(sgi + 205);
    const auto *sgi_210 = buffer.data(sgi + 210);
    const auto *sgi_216 = buffer.data(sgi + 216);
    const auto *sgi_217 = buffer.data(sgi + 217);
    const auto *sgi_218 = buffer.data(sgi + 218);
    const auto *sgi_219 = buffer.data(sgi + 219);
    const auto *sgi_220 = buffer.data(sgi + 220);
    const auto *sgi_221 = buffer.data(sgi + 221);
    const auto *sgi_222 = buffer.data(sgi + 222);
    const auto *sgi_223 = buffer.data(sgi + 223);
    const auto *sgi_245 = buffer.data(sgi + 245);
    const auto *sgi_246 = buffer.data(sgi + 246);
    const auto *sgi_247 = buffer.data(sgi + 247);
    const auto *sgi_248 = buffer.data(sgi + 248);
    const auto *sgi_249 = buffer.data(sgi + 249);
    const auto *sgi_250 = buffer.data(sgi + 250);
    const auto *sgi_251 = buffer.data(sgi + 251);
    const auto *sgi_252 = buffer.data(sgi + 252);
    const auto *sgi_255 = buffer.data(sgi + 255);
    const auto *sgi_257 = buffer.data(sgi + 257);
    const auto *sgi_258 = buffer.data(sgi + 258);
    const auto *sgi_261 = buffer.data(sgi + 261);
    const auto *sgi_262 = buffer.data(sgi + 262);
    const auto *sgi_264 = buffer.data(sgi + 264);
    const auto *sgi_266 = buffer.data(sgi + 266);
    const auto *sgi_267 = buffer.data(sgi + 267);
    const auto *sgi_269 = buffer.data(sgi + 269);
    const auto *sgi_270 = buffer.data(sgi + 270);
    const auto *sgi_272 = buffer.data(sgi + 272);
    const auto *sgi_273 = buffer.data(sgi + 273);
    const auto *sgi_274 = buffer.data(sgi + 274);

    const auto *sgk1_108 = buffer.data(sgk1 + 108);
    const auto *sgk1_111 = buffer.data(sgk1 + 111);
    const auto *sgk1_114 = buffer.data(sgk1 + 114);
    const auto *sgk1_118 = buffer.data(sgk1 + 118);
    const auto *sgk1_120 = buffer.data(sgk1 + 120);
    const auto *sgk1_123 = buffer.data(sgk1 + 123);
    const auto *sgk1_125 = buffer.data(sgk1 + 125);
    const auto *sgk1_126 = buffer.data(sgk1 + 126);
    const auto *sgk1_136 = buffer.data(sgk1 + 136);
    const auto *sgk1_180 = buffer.data(sgk1 + 180);
    const auto *sgk1_183 = buffer.data(sgk1 + 183);
    const auto *sgk1_185 = buffer.data(sgk1 + 185);
    const auto *sgk1_186 = buffer.data(sgk1 + 186);
    const auto *sgk1_189 = buffer.data(sgk1 + 189);
    const auto *sgk1_190 = buffer.data(sgk1 + 190);
    const auto *sgk1_192 = buffer.data(sgk1 + 192);
    const auto *sgk1_194 = buffer.data(sgk1 + 194);
    const auto *sgk1_195 = buffer.data(sgk1 + 195);
    const auto *sgk1_197 = buffer.data(sgk1 + 197);
    const auto *sgk1_198 = buffer.data(sgk1 + 198);
    const auto *sgk1_200 = buffer.data(sgk1 + 200);
    const auto *sgk1_215 = buffer.data(sgk1 + 215);

    const auto *shh0_141 = buffer.data(shh0 + 141);
    const auto *shh0_143 = buffer.data(shh0 + 143);
    const auto *shh0_144 = buffer.data(shh0 + 144);
    const auto *shh0_145 = buffer.data(shh0 + 145);
    const auto *shh0_146 = buffer.data(shh0 + 146);
    const auto *shh0_152 = buffer.data(shh0 + 152);
    const auto *shh0_156 = buffer.data(shh0 + 156);
    const auto *shh0_161 = buffer.data(shh0 + 161);
    const auto *shh0_164 = buffer.data(shh0 + 164);
    const auto *shh0_165 = buffer.data(shh0 + 165);
    const auto *shh0_166 = buffer.data(shh0 + 166);
    const auto *shh0_167 = buffer.data(shh0 + 167);
    const auto *shh0_183 = buffer.data(shh0 + 183);
    const auto *shh0_185 = buffer.data(shh0 + 185);
    const auto *shh0_186 = buffer.data(shh0 + 186);
    const auto *shh0_187 = buffer.data(shh0 + 187);
    const auto *shh0_188 = buffer.data(shh0 + 188);
    const auto *shh0_189 = buffer.data(shh0 + 189);
    const auto *shh0_192 = buffer.data(shh0 + 192);
    const auto *shh0_194 = buffer.data(shh0 + 194);
    const auto *shh0_195 = buffer.data(shh0 + 195);
    const auto *shh0_198 = buffer.data(shh0 + 198);
    const auto *shh0_199 = buffer.data(shh0 + 199);
    const auto *shh0_201 = buffer.data(shh0 + 201);
    const auto *shh0_203 = buffer.data(shh0 + 203);
    const auto *shh0_204 = buffer.data(shh0 + 204);
    const auto *shh0_206 = buffer.data(shh0 + 206);
    const auto *shh0_207 = buffer.data(shh0 + 207);
    const auto *shh0_209 = buffer.data(shh0 + 209);

    const auto *shh1_141 = buffer.data(shh1 + 141);
    const auto *shh1_143 = buffer.data(shh1 + 143);
    const auto *shh1_144 = buffer.data(shh1 + 144);
    const auto *shh1_145 = buffer.data(shh1 + 145);
    const auto *shh1_146 = buffer.data(shh1 + 146);
    const auto *shh1_152 = buffer.data(shh1 + 152);
    const auto *shh1_156 = buffer.data(shh1 + 156);
    const auto *shh1_161 = buffer.data(shh1 + 161);
    const auto *shh1_164 = buffer.data(shh1 + 164);
    const auto *shh1_165 = buffer.data(shh1 + 165);
    const auto *shh1_166 = buffer.data(shh1 + 166);
    const auto *shh1_167 = buffer.data(shh1 + 167);
    const auto *shh1_183 = buffer.data(shh1 + 183);
    const auto *shh1_185 = buffer.data(shh1 + 185);
    const auto *shh1_186 = buffer.data(shh1 + 186);
    const auto *shh1_187 = buffer.data(shh1 + 187);
    const auto *shh1_188 = buffer.data(shh1 + 188);
    const auto *shh1_189 = buffer.data(shh1 + 189);
    const auto *shh1_192 = buffer.data(shh1 + 192);
    const auto *shh1_194 = buffer.data(shh1 + 194);
    const auto *shh1_195 = buffer.data(shh1 + 195);
    const auto *shh1_198 = buffer.data(shh1 + 198);
    const auto *shh1_199 = buffer.data(shh1 + 199);
    const auto *shh1_201 = buffer.data(shh1 + 201);
    const auto *shh1_203 = buffer.data(shh1 + 203);
    const auto *shh1_204 = buffer.data(shh1 + 204);
    const auto *shh1_206 = buffer.data(shh1 + 206);
    const auto *shh1_207 = buffer.data(shh1 + 207);
    const auto *shh1_209 = buffer.data(shh1 + 209);

    const auto *shi_182 = buffer.data(shi + 182);
    const auto *shi_188 = buffer.data(shi + 188);
    const auto *shi_189 = buffer.data(shi + 189);
    const auto *shi_190 = buffer.data(shi + 190);
    const auto *shi_191 = buffer.data(shi + 191);
    const auto *shi_192 = buffer.data(shi + 192);
    const auto *shi_193 = buffer.data(shi + 193);
    const auto *shi_194 = buffer.data(shi + 194);
    const auto *shi_195 = buffer.data(shi + 195);
    const auto *shi_196 = buffer.data(shi + 196);
    const auto *shi_198 = buffer.data(shi + 198);
    const auto *shi_199 = buffer.data(shi + 199);
    const auto *shi_201 = buffer.data(shi + 201);
    const auto *shi_202 = buffer.data(shi + 202);
    const auto *shi_205 = buffer.data(shi + 205);
    const auto *shi_206 = buffer.data(shi + 206);
    const auto *shi_210 = buffer.data(shi + 210);
    const auto *shi_216 = buffer.data(shi + 216);
    const auto *shi_217 = buffer.data(shi + 217);
    const auto *shi_218 = buffer.data(shi + 218);
    const auto *shi_219 = buffer.data(shi + 219);
    const auto *shi_220 = buffer.data(shi + 220);
    const auto *shi_221 = buffer.data(shi + 221);
    const auto *shi_222 = buffer.data(shi + 222);
    const auto *shi_223 = buffer.data(shi + 223);
    const auto *shi_224 = buffer.data(shi + 224);
    const auto *shi_226 = buffer.data(shi + 226);
    const auto *shi_227 = buffer.data(shi + 227);
    const auto *shi_229 = buffer.data(shi + 229);
    const auto *shi_230 = buffer.data(shi + 230);
    const auto *shi_233 = buffer.data(shi + 233);
    const auto *shi_234 = buffer.data(shi + 234);
    const auto *shi_238 = buffer.data(shi + 238);
    const auto *shi_245 = buffer.data(shi + 245);
    const auto *shi_246 = buffer.data(shi + 246);
    const auto *shi_247 = buffer.data(shi + 247);
    const auto *shi_248 = buffer.data(shi + 248);
    const auto *shi_249 = buffer.data(shi + 249);
    const auto *shi_250 = buffer.data(shi + 250);
    const auto *shi_251 = buffer.data(shi + 251);
    const auto *shi_252 = buffer.data(shi + 252);
    const auto *shi_254 = buffer.data(shi + 254);
    const auto *shi_255 = buffer.data(shi + 255);
    const auto *shi_257 = buffer.data(shi + 257);
    const auto *shi_258 = buffer.data(shi + 258);
    const auto *shi_261 = buffer.data(shi + 261);
    const auto *shi_262 = buffer.data(shi + 262);
    const auto *shi_264 = buffer.data(shi + 264);
    const auto *shi_266 = buffer.data(shi + 266);
    const auto *shi_267 = buffer.data(shi + 267);
    const auto *shi_269 = buffer.data(shi + 269);
    const auto *shi_270 = buffer.data(shi + 270);
    const auto *shi_272 = buffer.data(shi + 272);
    const auto *shi_273 = buffer.data(shi + 273);
    const auto *shi_274 = buffer.data(shi + 274);

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pc_x, pc_y, sgi_98, sgi_188, sgi_189, \
                         sgi_190, shh0_146, shh1_146, shi_182, shi_188, shi_189, \
                         shi_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_15 * sgi_98[k]
                   + f_3 * pc_y[k] * shi_182[k];

        t_236[k] = f_14 * sgi_188[k]
                   + f_10 * shh0_146[k]
                   - f_11 * shh1_146[k]
                   + f_3 * pc_x[k] * shi_188[k];

        t_237[k] = f_14 * sgi_189[k]
                   + f_3 * pc_x[k] * shi_189[k];

        t_238[k] = f_14 * sgi_190[k]
                   + f_3 * pc_x[k] * shi_190[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pc_x, sgi_191, sgi_192, sgi_193, \
                         sgi_194, sgi_195, shi_191, shi_192, shi_193, shi_194, \
                         shi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_14 * sgi_191[k]
                   + f_3 * pc_x[k] * shi_191[k];

        t_240[k] = f_14 * sgi_192[k]
                   + f_3 * pc_x[k] * shi_192[k];

        t_241[k] = f_14 * sgi_193[k]
                   + f_3 * pc_x[k] * shi_193[k];

        t_242[k] = f_14 * sgi_194[k]
                   + f_3 * pc_x[k] * shi_194[k];

        t_243[k] = f_14 * sgi_195[k]
                   + f_3 * pc_x[k] * shi_195[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, pc_y, pc_z, sgi_105, sgi_107, shh0_141, \
                         shh0_143, shh1_141, shh1_143, shi_189, \
                         shi_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_15 * sgi_105[k]
                   + f_1 * shh0_141[k]
                   - f_2 * shh1_141[k]
                   + f_3 * pc_y[k] * shi_189[k];

        t_245[k] = f_3 * pc_z[k] * shi_189[k];

        t_246[k] = f_15 * sgi_107[k]
                   + f_4 * shh0_143[k]
                   - f_5 * shh1_143[k]
                   + f_3 * pc_y[k] * shi_191[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_y, sgi_108, sgi_109, sgi_110, shh0_144, \
                         shh0_145, shh0_146, shh1_144, shh1_145, shh1_146, shi_192, shi_193, \
                         shi_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_15 * sgi_108[k]
                   + f_6 * shh0_144[k]
                   - f_7 * shh1_144[k]
                   + f_3 * pc_y[k] * shi_192[k];

        t_248[k] = f_15 * sgi_109[k]
                   + f_8 * shh0_145[k]
                   - f_9 * shh1_145[k]
                   + f_3 * pc_y[k] * shi_193[k];

        t_249[k] = f_15 * sgi_110[k]
                   + f_10 * shh0_146[k]
                   - f_11 * shh1_146[k]
                   + f_3 * pc_y[k] * shi_194[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pb_z, pc_y, pc_z, sgk0_108, sgi_111, \
                         sgi_112, sgk1_108, shh0_146, shh1_146, shi_195, \
                         shi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_15 * sgi_111[k]
                   + f_3 * pc_y[k] * shi_195[k];

        t_251[k] = f_1 * shh0_146[k]
                   - f_2 * shh1_146[k]
                   + f_3 * pc_z[k] * shi_195[k];

        t_252[k] = pb_z[k] * sgk0_108[k]
                   - f_12 * pc_z[k] * sgk1_108[k];

        t_253[k] = f_14 * sgi_112[k]
                   + f_3 * pc_y[k] * shi_196[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pb_z, pc_y, pc_z, sgk0_111, sgi_84, sgi_114, \
                         sgk1_111, shi_196, shi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * sgi_84[k]
                   + f_3 * pc_z[k] * shi_196[k];

        t_255[k] = pb_z[k] * sgk0_111[k]
                   - f_12 * pc_z[k] * sgk1_111[k];

        t_256[k] = f_14 * sgi_114[k]
                   + f_3 * pc_y[k] * shi_198[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pb_z, pc_x, pc_z, sgk0_114, sgi_87, sgi_201, \
                         sgk1_114, shh0_152, shh1_152, shi_199, \
                         shi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_14 * sgi_201[k]
                   + f_4 * shh0_152[k]
                   - f_5 * shh1_152[k]
                   + f_3 * pc_x[k] * shi_201[k];

        t_258[k] = pb_z[k] * sgk0_114[k]
                   - f_12 * pc_z[k] * sgk1_114[k];

        t_259[k] = f_13 * sgi_87[k]
                   + f_3 * pc_z[k] * shi_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pb_z, pc_x, pc_y, pc_z, sgk0_118, sgi_117, \
                         sgi_205, sgk1_118, shh0_156, shh1_156, shi_201, \
                         shi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * sgi_117[k]
                   + f_3 * pc_y[k] * shi_201[k];

        t_261[k] = f_14 * sgi_205[k]
                   + f_6 * shh0_156[k]
                   - f_7 * shh1_156[k]
                   + f_3 * pc_x[k] * shi_205[k];

        t_262[k] = pb_z[k] * sgk0_118[k]
                   - f_12 * pc_z[k] * sgk1_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pb_z, pc_y, pc_z, sgk0_120, sgi_90, sgi_91, \
                         sgi_121, sgk1_120, shi_202, shi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_13 * sgi_90[k]
                   + f_3 * pc_z[k] * shi_202[k];

        t_264[k] = pb_z[k] * sgk0_120[k]
                   + f_14 * sgi_91[k]
                   - f_12 * pc_z[k] * sgk1_120[k];

        t_265[k] = f_14 * sgi_121[k]
                   + f_3 * pc_y[k] * shi_205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pb_z, pc_x, pc_z, sgk0_123, sgi_94, sgi_210, \
                         sgk1_123, shh0_161, shh1_161, shi_206, \
                         shi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_14 * sgi_210[k]
                   + f_8 * shh0_161[k]
                   - f_9 * shh1_161[k]
                   + f_3 * pc_x[k] * shi_210[k];

        t_267[k] = pb_z[k] * sgk0_123[k]
                   - f_12 * pc_z[k] * sgk1_123[k];

        t_268[k] = f_13 * sgi_94[k]
                   + f_3 * pc_z[k] * shi_206[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pb_z, pc_y, pc_z, sgk0_125, sgk0_126, sgi_95, \
                         sgi_96, sgi_126, sgk1_125, sgk1_126, shi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pb_z[k] * sgk0_125[k]
                   + f_14 * sgi_95[k]
                   - f_12 * pc_z[k] * sgk1_125[k];

        t_270[k] = pb_z[k] * sgk0_126[k]
                   + f_15 * sgi_96[k]
                   - f_12 * pc_z[k] * sgk1_126[k];

        t_271[k] = f_14 * sgi_126[k]
                   + f_3 * pc_y[k] * shi_210[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, sgi_216, sgi_217, sgi_218, sgi_219, \
                         shh0_167, shh1_167, shi_216, shi_217, shi_218, \
                         shi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_14 * sgi_216[k]
                   + f_10 * shh0_167[k]
                   - f_11 * shh1_167[k]
                   + f_3 * pc_x[k] * shi_216[k];

        t_273[k] = f_14 * sgi_217[k]
                   + f_3 * pc_x[k] * shi_217[k];

        t_274[k] = f_14 * sgi_218[k]
                   + f_3 * pc_x[k] * shi_218[k];

        t_275[k] = f_14 * sgi_219[k]
                   + f_3 * pc_x[k] * shi_219[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, sgi_220, sgi_221, sgi_222, sgi_223, \
                         shi_220, shi_221, shi_222, shi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_14 * sgi_220[k]
                   + f_3 * pc_x[k] * shi_220[k];

        t_277[k] = f_14 * sgi_221[k]
                   + f_3 * pc_x[k] * shi_221[k];

        t_278[k] = f_14 * sgi_222[k]
                   + f_3 * pc_x[k] * shi_222[k];

        t_279[k] = f_14 * sgi_223[k]
                   + f_3 * pc_x[k] * shi_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pb_z, pc_y, pc_z, sgk0_136, sgi_105, sgi_135, \
                         sgk1_136, shh0_164, shh1_164, shi_217, \
                         shi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pb_z[k] * sgk0_136[k]
                   - f_12 * pc_z[k] * sgk1_136[k];

        t_281[k] = f_13 * sgi_105[k]
                   + f_3 * pc_z[k] * shi_217[k];

        t_282[k] = f_14 * sgi_135[k]
                   + f_4 * shh0_164[k]
                   - f_5 * shh1_164[k]
                   + f_3 * pc_y[k] * shi_219[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_y, sgi_136, sgi_137, sgi_138, shh0_165, \
                         shh0_166, shh0_167, shh1_165, shh1_166, shh1_167, shi_220, shi_221, \
                         shi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_14 * sgi_136[k]
                   + f_6 * shh0_165[k]
                   - f_7 * shh1_165[k]
                   + f_3 * pc_y[k] * shi_220[k];

        t_284[k] = f_14 * sgi_137[k]
                   + f_8 * shh0_166[k]
                   - f_9 * shh1_166[k]
                   + f_3 * pc_y[k] * shi_221[k];

        t_285[k] = f_14 * sgi_138[k]
                   + f_10 * shh0_167[k]
                   - f_11 * shh1_167[k]
                   + f_3 * pc_y[k] * shi_222[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pc_y, pc_z, sgk0_180, sgi_111, \
                         sgi_139, sgi_140, sgk1_180, shh0_167, shh1_167, shi_223, \
                         shi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * sgi_139[k]
                   + f_3 * pc_y[k] * shi_223[k];

        t_287[k] = f_13 * sgi_111[k]
                   + f_1 * shh0_167[k]
                   - f_2 * shh1_167[k]
                   + f_3 * pc_z[k] * shi_223[k];

        t_288[k] = pb_y[k] * sgk0_180[k]
                   - f_12 * pc_y[k] * sgk1_180[k];

        t_289[k] = f_13 * sgi_140[k]
                   + f_3 * pc_y[k] * shi_224[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_y, pc_y, pc_z, sgk0_183, sgk0_185, \
                         sgi_112, sgi_141, sgi_142, sgk1_183, sgk1_185, shi_224, \
                         shi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * sgi_112[k]
                   + f_3 * pc_z[k] * shi_224[k];

        t_291[k] = pb_y[k] * sgk0_183[k]
                   + f_14 * sgi_141[k]
                   - f_12 * pc_y[k] * sgk1_183[k];

        t_292[k] = f_13 * sgi_142[k]
                   + f_3 * pc_y[k] * shi_226[k];

        t_293[k] = pb_y[k] * sgk0_185[k]
                   - f_12 * pc_y[k] * sgk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pb_y, pc_y, pc_z, sgk0_186, sgk0_189, \
                         sgi_115, sgi_143, sgi_145, sgk1_186, sgk1_189, shi_227, \
                         shi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pb_y[k] * sgk0_186[k]
                   + f_15 * sgi_143[k]
                   - f_12 * pc_y[k] * sgk1_186[k];

        t_295[k] = f_14 * sgi_115[k]
                   + f_3 * pc_z[k] * shi_227[k];

        t_296[k] = f_13 * sgi_145[k]
                   + f_3 * pc_y[k] * shi_229[k];

        t_297[k] = pb_y[k] * sgk0_189[k]
                   - f_12 * pc_y[k] * sgk1_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pb_y, pc_y, pc_z, sgk0_190, sgk0_192, sgi_118, \
                         sgi_146, sgi_148, sgk1_190, sgk1_192, \
                         shi_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = pb_y[k] * sgk0_190[k]
                   + f_16 * sgi_146[k]
                   - f_12 * pc_y[k] * sgk1_190[k];

        t_299[k] = f_14 * sgi_118[k]
                   + f_3 * pc_z[k] * shi_230[k];

        t_300[k] = pb_y[k] * sgk0_192[k]
                   + f_14 * sgi_148[k]
                   - f_12 * pc_y[k] * sgk1_192[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pb_y, pc_y, pc_z, sgk0_194, sgk0_195, \
                         sgi_122, sgi_149, sgi_150, sgk1_194, sgk1_195, shi_233, \
                         shi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_13 * sgi_149[k]
                   + f_3 * pc_y[k] * shi_233[k];

        t_302[k] = pb_y[k] * sgk0_194[k]
                   - f_12 * pc_y[k] * sgk1_194[k];

        t_303[k] = pb_y[k] * sgk0_195[k]
                   + f_0 * sgi_150[k]
                   - f_12 * pc_y[k] * sgk1_195[k];

        t_304[k] = f_14 * sgi_122[k]
                   + f_3 * pc_z[k] * shi_234[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pb_y, pc_y, sgk0_197, sgk0_198, sgk0_200, \
                         sgi_152, sgi_153, sgi_154, sgk1_197, sgk1_198, sgk1_200, \
                         shi_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pb_y[k] * sgk0_197[k]
                   + f_15 * sgi_152[k]
                   - f_12 * pc_y[k] * sgk1_197[k];

        t_306[k] = pb_y[k] * sgk0_198[k]
                   + f_14 * sgi_153[k]
                   - f_12 * pc_y[k] * sgk1_198[k];

        t_307[k] = f_13 * sgi_154[k]
                   + f_3 * pc_y[k] * shi_238[k];

        t_308[k] = pb_y[k] * sgk0_200[k]
                   - f_12 * pc_y[k] * sgk1_200[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pc_x, sgi_245, sgi_246, sgi_247, \
                         sgi_248, sgi_249, shi_245, shi_246, shi_247, shi_248, \
                         shi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_14 * sgi_245[k]
                   + f_3 * pc_x[k] * shi_245[k];

        t_310[k] = f_14 * sgi_246[k]
                   + f_3 * pc_x[k] * shi_246[k];

        t_311[k] = f_14 * sgi_247[k]
                   + f_3 * pc_x[k] * shi_247[k];

        t_312[k] = f_14 * sgi_248[k]
                   + f_3 * pc_x[k] * shi_248[k];

        t_313[k] = f_14 * sgi_249[k]
                   + f_3 * pc_x[k] * shi_249[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, sgi_133, sgi_161, \
                         sgi_250, sgi_251, shh0_183, shh1_183, shi_245, shi_250, \
                         shi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_14 * sgi_250[k]
                   + f_3 * pc_x[k] * shi_250[k];

        t_315[k] = f_14 * sgi_251[k]
                   + f_3 * pc_x[k] * shi_251[k];

        t_316[k] = f_13 * sgi_161[k]
                   + f_1 * shh0_183[k]
                   - f_2 * shh1_183[k]
                   + f_3 * pc_y[k] * shi_245[k];

        t_317[k] = f_14 * sgi_133[k]
                   + f_3 * pc_z[k] * shi_245[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_y, sgi_163, sgi_164, sgi_165, shh0_185, \
                         shh0_186, shh0_187, shh1_185, shh1_186, shh1_187, shi_247, shi_248, \
                         shi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_13 * sgi_163[k]
                   + f_4 * shh0_185[k]
                   - f_5 * shh1_185[k]
                   + f_3 * pc_y[k] * shi_247[k];

        t_319[k] = f_13 * sgi_164[k]
                   + f_6 * shh0_186[k]
                   - f_7 * shh1_186[k]
                   + f_3 * pc_y[k] * shi_248[k];

        t_320[k] = f_13 * sgi_165[k]
                   + f_8 * shh0_187[k]
                   - f_9 * shh1_187[k]
                   + f_3 * pc_y[k] * shi_249[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pb_y, pc_y, sgk0_215, sgi_166, sgi_167, \
                         sgk1_215, shh0_188, shh1_188, shi_250, \
                         shi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * sgi_166[k]
                   + f_10 * shh0_188[k]
                   - f_11 * shh1_188[k]
                   + f_3 * pc_y[k] * shi_250[k];

        t_322[k] = f_13 * sgi_167[k]
                   + f_3 * pc_y[k] * shi_251[k];

        t_323[k] = pb_y[k] * sgk0_215[k]
                   - f_12 * pc_y[k] * sgk1_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pc_x, pc_y, pc_z, sgi_140, sgi_252, \
                         sgi_255, shh0_189, shh0_192, shh1_189, shh1_192, shi_252, \
                         shi_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_14 * sgi_252[k]
                   + f_1 * shh0_189[k]
                   - f_2 * shh1_189[k]
                   + f_3 * pc_x[k] * shi_252[k];

        t_325[k] = f_3 * pc_y[k] * shi_252[k];

        t_326[k] = f_15 * sgi_140[k]
                   + f_3 * pc_z[k] * shi_252[k];

        t_327[k] = f_14 * sgi_255[k]
                   + f_4 * shh0_192[k]
                   - f_5 * shh1_192[k]
                   + f_3 * pc_x[k] * shi_255[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pc_x, pc_y, sgi_257, sgi_258, shh0_194, \
                         shh0_195, shh1_194, shh1_195, shi_254, shi_257, \
                         shi_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_3 * pc_y[k] * shi_254[k];

        t_329[k] = f_14 * sgi_257[k]
                   + f_4 * shh0_194[k]
                   - f_5 * shh1_194[k]
                   + f_3 * pc_x[k] * shi_257[k];

        t_330[k] = f_14 * sgi_258[k]
                   + f_6 * shh0_195[k]
                   - f_7 * shh1_195[k]
                   + f_3 * pc_x[k] * shi_258[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pc_x, pc_y, pc_z, sgi_143, sgi_261, shh0_198, \
                         shh1_198, shi_255, shi_257, shi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_15 * sgi_143[k]
                   + f_3 * pc_z[k] * shi_255[k];

        t_332[k] = f_3 * pc_y[k] * shi_257[k];

        t_333[k] = f_14 * sgi_261[k]
                   + f_6 * shh0_198[k]
                   - f_7 * shh1_198[k]
                   + f_3 * pc_x[k] * shi_261[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pc_x, pc_z, sgi_146, sgi_262, sgi_264, shh0_199, \
                         shh0_201, shh1_199, shh1_201, shi_258, shi_262, \
                         shi_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_14 * sgi_262[k]
                   + f_8 * shh0_199[k]
                   - f_9 * shh1_199[k]
                   + f_3 * pc_x[k] * shi_262[k];

        t_335[k] = f_15 * sgi_146[k]
                   + f_3 * pc_z[k] * shi_258[k];

        t_336[k] = f_14 * sgi_264[k]
                   + f_8 * shh0_201[k]
                   - f_9 * shh1_201[k]
                   + f_3 * pc_x[k] * shi_264[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, pc_x, pc_y, sgi_266, sgi_267, shh0_203, \
                         shh0_204, shh1_203, shh1_204, shi_261, shi_266, \
                         shi_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_3 * pc_y[k] * shi_261[k];

        t_338[k] = f_14 * sgi_266[k]
                   + f_8 * shh0_203[k]
                   - f_9 * shh1_203[k]
                   + f_3 * pc_x[k] * shi_266[k];

        t_339[k] = f_14 * sgi_267[k]
                   + f_10 * shh0_204[k]
                   - f_11 * shh1_204[k]
                   + f_3 * pc_x[k] * shi_267[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pc_x, pc_z, sgi_150, sgi_269, sgi_270, shh0_206, \
                         shh0_207, shh1_206, shh1_207, shi_262, shi_269, \
                         shi_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_15 * sgi_150[k]
                   + f_3 * pc_z[k] * shi_262[k];

        t_341[k] = f_14 * sgi_269[k]
                   + f_10 * shh0_206[k]
                   - f_11 * shh1_206[k]
                   + f_3 * pc_x[k] * shi_269[k];

        t_342[k] = f_14 * sgi_270[k]
                   + f_10 * shh0_207[k]
                   - f_11 * shh1_207[k]
                   + f_3 * pc_x[k] * shi_270[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, pc_x, pc_y, sgi_272, sgi_273, sgi_274, \
                         shh0_209, shh1_209, shi_266, shi_272, shi_273, \
                         shi_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_3 * pc_y[k] * shi_266[k];

        t_344[k] = f_14 * sgi_272[k]
                   + f_10 * shh0_209[k]
                   - f_11 * shh1_209[k]
                   + f_3 * pc_x[k] * shi_272[k];

        t_345[k] = f_14 * sgi_273[k]
                   + f_3 * pc_x[k] * shi_273[k];

        t_346[k] = f_14 * sgi_274[k]
                   + f_3 * pc_x[k] * shi_274[k];
    }
}

static auto
compute_prim_shk_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgk0,
                                                          const size_t sgi, const size_t sgk1,
                                                          const size_t shh0, const size_t shh1,
                                                          const size_t shi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgk0_216 = buffer.data(sgk0 + 216);
    const auto *sgk0_219 = buffer.data(sgk0 + 219);
    const auto *sgk0_222 = buffer.data(sgk0 + 222);
    const auto *sgk0_226 = buffer.data(sgk0 + 226);
    const auto *sgk0_231 = buffer.data(sgk0 + 231);
    const auto *sgk0_324 = buffer.data(sgk0 + 324);
    const auto *sgk0_360 = buffer.data(sgk0 + 360);
    const auto *sgk0_363 = buffer.data(sgk0 + 363);
    const auto *sgk0_365 = buffer.data(sgk0 + 365);
    const auto *sgk0_366 = buffer.data(sgk0 + 366);
    const auto *sgk0_369 = buffer.data(sgk0 + 369);
    const auto *sgk0_370 = buffer.data(sgk0 + 370);
    const auto *sgk0_372 = buffer.data(sgk0 + 372);
    const auto *sgk0_374 = buffer.data(sgk0 + 374);
    const auto *sgk0_375 = buffer.data(sgk0 + 375);
    const auto *sgk0_377 = buffer.data(sgk0 + 377);
    const auto *sgk0_378 = buffer.data(sgk0 + 378);
    const auto *sgk0_380 = buffer.data(sgk0 + 380);
    const auto *sgk0_388 = buffer.data(sgk0 + 388);
    const auto *sgk0_390 = buffer.data(sgk0 + 390);
    const auto *sgk0_391 = buffer.data(sgk0 + 391);
    const auto *sgk0_392 = buffer.data(sgk0 + 392);
    const auto *sgk0_393 = buffer.data(sgk0 + 393);
    const auto *sgk0_395 = buffer.data(sgk0 + 395);
    const auto *sgk0_401 = buffer.data(sgk0 + 401);
    const auto *sgk0_405 = buffer.data(sgk0 + 405);
    const auto *sgk0_408 = buffer.data(sgk0 + 408);
    const auto *sgk0_410 = buffer.data(sgk0 + 410);
    const auto *sgk0_413 = buffer.data(sgk0 + 413);
    const auto *sgk0_414 = buffer.data(sgk0 + 414);
    const auto *sgk0_416 = buffer.data(sgk0 + 416);
    const auto *sgk0_424 = buffer.data(sgk0 + 424);
    const auto *sgk0_426 = buffer.data(sgk0 + 426);
    const auto *sgk0_427 = buffer.data(sgk0 + 427);
    const auto *sgk0_428 = buffer.data(sgk0 + 428);
    const auto *sgk0_429 = buffer.data(sgk0 + 429);
    const auto *sgk0_431 = buffer.data(sgk0 + 431);
    const auto *sgk0_432 = buffer.data(sgk0 + 432);
    const auto *sgk0_435 = buffer.data(sgk0 + 435);
    const auto *sgk0_437 = buffer.data(sgk0 + 437);
    const auto *sgk0_438 = buffer.data(sgk0 + 438);
    const auto *sgk0_441 = buffer.data(sgk0 + 441);
    const auto *sgk0_442 = buffer.data(sgk0 + 442);
    const auto *sgk0_444 = buffer.data(sgk0 + 444);
    const auto *sgk0_446 = buffer.data(sgk0 + 446);
    const auto *sgk0_447 = buffer.data(sgk0 + 447);
    const auto *sgk0_449 = buffer.data(sgk0 + 449);
    const auto *sgk0_450 = buffer.data(sgk0 + 450);
    const auto *sgk0_452 = buffer.data(sgk0 + 452);
    const auto *sgk0_460 = buffer.data(sgk0 + 460);
    const auto *sgk0_462 = buffer.data(sgk0 + 462);
    const auto *sgk0_463 = buffer.data(sgk0 + 463);
    const auto *sgk0_464 = buffer.data(sgk0 + 464);
    const auto *sgk0_465 = buffer.data(sgk0 + 465);
    const auto *sgk0_467 = buffer.data(sgk0 + 467);

    const auto *sgi_161 = buffer.data(sgi + 161);
    const auto *sgi_167 = buffer.data(sgi + 167);
    const auto *sgi_168 = buffer.data(sgi + 168);
    const auto *sgi_170 = buffer.data(sgi + 170);
    const auto *sgi_171 = buffer.data(sgi + 171);
    const auto *sgi_173 = buffer.data(sgi + 173);
    const auto *sgi_174 = buffer.data(sgi + 174);
    const auto *sgi_177 = buffer.data(sgi + 177);
    const auto *sgi_178 = buffer.data(sgi + 178);
    const auto *sgi_182 = buffer.data(sgi + 182);
    const auto *sgi_189 = buffer.data(sgi + 189);
    const auto *sgi_195 = buffer.data(sgi + 195);
    const auto *sgi_196 = buffer.data(sgi + 196);
    const auto *sgi_198 = buffer.data(sgi + 198);
    const auto *sgi_199 = buffer.data(sgi + 199);
    const auto *sgi_201 = buffer.data(sgi + 201);
    const auto *sgi_202 = buffer.data(sgi + 202);
    const auto *sgi_205 = buffer.data(sgi + 205);
    const auto *sgi_206 = buffer.data(sgi + 206);
    const auto *sgi_210 = buffer.data(sgi + 210);
    const auto *sgi_217 = buffer.data(sgi + 217);
    const auto *sgi_223 = buffer.data(sgi + 223);
    const auto *sgi_224 = buffer.data(sgi + 224);
    const auto *sgi_226 = buffer.data(sgi + 226);
    const auto *sgi_229 = buffer.data(sgi + 229);
    const auto *sgi_233 = buffer.data(sgi + 233);
    const auto *sgi_238 = buffer.data(sgi + 238);
    const auto *sgi_251 = buffer.data(sgi + 251);
    const auto *sgi_252 = buffer.data(sgi + 252);
    const auto *sgi_275 = buffer.data(sgi + 275);
    const auto *sgi_276 = buffer.data(sgi + 276);
    const auto *sgi_277 = buffer.data(sgi + 277);
    const auto *sgi_278 = buffer.data(sgi + 278);
    const auto *sgi_279 = buffer.data(sgi + 279);
    const auto *sgi_280 = buffer.data(sgi + 280);
    const auto *sgi_283 = buffer.data(sgi + 283);
    const auto *sgi_285 = buffer.data(sgi + 285);
    const auto *sgi_286 = buffer.data(sgi + 286);
    const auto *sgi_289 = buffer.data(sgi + 289);
    const auto *sgi_290 = buffer.data(sgi + 290);
    const auto *sgi_292 = buffer.data(sgi + 292);
    const auto *sgi_294 = buffer.data(sgi + 294);
    const auto *sgi_295 = buffer.data(sgi + 295);
    const auto *sgi_297 = buffer.data(sgi + 297);
    const auto *sgi_298 = buffer.data(sgi + 298);
    const auto *sgi_300 = buffer.data(sgi + 300);
    const auto *sgi_301 = buffer.data(sgi + 301);
    const auto *sgi_302 = buffer.data(sgi + 302);
    const auto *sgi_303 = buffer.data(sgi + 303);
    const auto *sgi_304 = buffer.data(sgi + 304);
    const auto *sgi_305 = buffer.data(sgi + 305);
    const auto *sgi_306 = buffer.data(sgi + 306);
    const auto *sgi_307 = buffer.data(sgi + 307);
    const auto *sgi_313 = buffer.data(sgi + 313);
    const auto *sgi_317 = buffer.data(sgi + 317);
    const auto *sgi_320 = buffer.data(sgi + 320);
    const auto *sgi_322 = buffer.data(sgi + 322);
    const auto *sgi_325 = buffer.data(sgi + 325);
    const auto *sgi_326 = buffer.data(sgi + 326);
    const auto *sgi_328 = buffer.data(sgi + 328);
    const auto *sgi_329 = buffer.data(sgi + 329);
    const auto *sgi_330 = buffer.data(sgi + 330);
    const auto *sgi_331 = buffer.data(sgi + 331);
    const auto *sgi_332 = buffer.data(sgi + 332);
    const auto *sgi_333 = buffer.data(sgi + 333);
    const auto *sgi_334 = buffer.data(sgi + 334);
    const auto *sgi_335 = buffer.data(sgi + 335);
    const auto *sgi_336 = buffer.data(sgi + 336);
    const auto *sgi_339 = buffer.data(sgi + 339);
    const auto *sgi_341 = buffer.data(sgi + 341);
    const auto *sgi_342 = buffer.data(sgi + 342);
    const auto *sgi_345 = buffer.data(sgi + 345);
    const auto *sgi_346 = buffer.data(sgi + 346);
    const auto *sgi_348 = buffer.data(sgi + 348);
    const auto *sgi_350 = buffer.data(sgi + 350);
    const auto *sgi_351 = buffer.data(sgi + 351);
    const auto *sgi_353 = buffer.data(sgi + 353);
    const auto *sgi_354 = buffer.data(sgi + 354);
    const auto *sgi_356 = buffer.data(sgi + 356);
    const auto *sgi_357 = buffer.data(sgi + 357);
    const auto *sgi_358 = buffer.data(sgi + 358);
    const auto *sgi_359 = buffer.data(sgi + 359);
    const auto *sgi_360 = buffer.data(sgi + 360);
    const auto *sgi_361 = buffer.data(sgi + 361);
    const auto *sgi_362 = buffer.data(sgi + 362);
    const auto *sgi_363 = buffer.data(sgi + 363);

    const auto *sgk1_216 = buffer.data(sgk1 + 216);
    const auto *sgk1_219 = buffer.data(sgk1 + 219);
    const auto *sgk1_222 = buffer.data(sgk1 + 222);
    const auto *sgk1_226 = buffer.data(sgk1 + 226);
    const auto *sgk1_231 = buffer.data(sgk1 + 231);
    const auto *sgk1_324 = buffer.data(sgk1 + 324);
    const auto *sgk1_360 = buffer.data(sgk1 + 360);
    const auto *sgk1_363 = buffer.data(sgk1 + 363);
    const auto *sgk1_365 = buffer.data(sgk1 + 365);
    const auto *sgk1_366 = buffer.data(sgk1 + 366);
    const auto *sgk1_369 = buffer.data(sgk1 + 369);
    const auto *sgk1_370 = buffer.data(sgk1 + 370);
    const auto *sgk1_372 = buffer.data(sgk1 + 372);
    const auto *sgk1_374 = buffer.data(sgk1 + 374);
    const auto *sgk1_375 = buffer.data(sgk1 + 375);
    const auto *sgk1_377 = buffer.data(sgk1 + 377);
    const auto *sgk1_378 = buffer.data(sgk1 + 378);
    const auto *sgk1_380 = buffer.data(sgk1 + 380);
    const auto *sgk1_388 = buffer.data(sgk1 + 388);
    const auto *sgk1_390 = buffer.data(sgk1 + 390);
    const auto *sgk1_391 = buffer.data(sgk1 + 391);
    const auto *sgk1_392 = buffer.data(sgk1 + 392);
    const auto *sgk1_393 = buffer.data(sgk1 + 393);
    const auto *sgk1_395 = buffer.data(sgk1 + 395);
    const auto *sgk1_401 = buffer.data(sgk1 + 401);
    const auto *sgk1_405 = buffer.data(sgk1 + 405);
    const auto *sgk1_408 = buffer.data(sgk1 + 408);
    const auto *sgk1_410 = buffer.data(sgk1 + 410);
    const auto *sgk1_413 = buffer.data(sgk1 + 413);
    const auto *sgk1_414 = buffer.data(sgk1 + 414);
    const auto *sgk1_416 = buffer.data(sgk1 + 416);
    const auto *sgk1_424 = buffer.data(sgk1 + 424);
    const auto *sgk1_426 = buffer.data(sgk1 + 426);
    const auto *sgk1_427 = buffer.data(sgk1 + 427);
    const auto *sgk1_428 = buffer.data(sgk1 + 428);
    const auto *sgk1_429 = buffer.data(sgk1 + 429);
    const auto *sgk1_431 = buffer.data(sgk1 + 431);
    const auto *sgk1_432 = buffer.data(sgk1 + 432);
    const auto *sgk1_435 = buffer.data(sgk1 + 435);
    const auto *sgk1_437 = buffer.data(sgk1 + 437);
    const auto *sgk1_438 = buffer.data(sgk1 + 438);
    const auto *sgk1_441 = buffer.data(sgk1 + 441);
    const auto *sgk1_442 = buffer.data(sgk1 + 442);
    const auto *sgk1_444 = buffer.data(sgk1 + 444);
    const auto *sgk1_446 = buffer.data(sgk1 + 446);
    const auto *sgk1_447 = buffer.data(sgk1 + 447);
    const auto *sgk1_449 = buffer.data(sgk1 + 449);
    const auto *sgk1_450 = buffer.data(sgk1 + 450);
    const auto *sgk1_452 = buffer.data(sgk1 + 452);
    const auto *sgk1_460 = buffer.data(sgk1 + 460);
    const auto *sgk1_462 = buffer.data(sgk1 + 462);
    const auto *sgk1_463 = buffer.data(sgk1 + 463);
    const auto *sgk1_464 = buffer.data(sgk1 + 464);
    const auto *sgk1_465 = buffer.data(sgk1 + 465);
    const auto *sgk1_467 = buffer.data(sgk1 + 467);

    const auto *shh0_204 = buffer.data(shh0 + 204);
    const auto *shh0_206 = buffer.data(shh0 + 206);
    const auto *shh0_207 = buffer.data(shh0 + 207);
    const auto *shh0_208 = buffer.data(shh0 + 208);
    const auto *shh0_209 = buffer.data(shh0 + 209);

    const auto *shh1_204 = buffer.data(shh1 + 204);
    const auto *shh1_206 = buffer.data(shh1 + 206);
    const auto *shh1_207 = buffer.data(shh1 + 207);
    const auto *shh1_208 = buffer.data(shh1 + 208);
    const auto *shh1_209 = buffer.data(shh1 + 209);

    const auto *shi_273 = buffer.data(shi + 273);
    const auto *shi_275 = buffer.data(shi + 275);
    const auto *shi_276 = buffer.data(shi + 276);
    const auto *shi_277 = buffer.data(shi + 277);
    const auto *shi_278 = buffer.data(shi + 278);
    const auto *shi_279 = buffer.data(shi + 279);
    const auto *shi_280 = buffer.data(shi + 280);
    const auto *shi_282 = buffer.data(shi + 282);
    const auto *shi_283 = buffer.data(shi + 283);
    const auto *shi_285 = buffer.data(shi + 285);
    const auto *shi_286 = buffer.data(shi + 286);
    const auto *shi_289 = buffer.data(shi + 289);
    const auto *shi_290 = buffer.data(shi + 290);
    const auto *shi_294 = buffer.data(shi + 294);
    const auto *shi_301 = buffer.data(shi + 301);
    const auto *shi_302 = buffer.data(shi + 302);
    const auto *shi_303 = buffer.data(shi + 303);
    const auto *shi_304 = buffer.data(shi + 304);
    const auto *shi_305 = buffer.data(shi + 305);
    const auto *shi_306 = buffer.data(shi + 306);
    const auto *shi_307 = buffer.data(shi + 307);
    const auto *shi_308 = buffer.data(shi + 308);
    const auto *shi_310 = buffer.data(shi + 310);
    const auto *shi_311 = buffer.data(shi + 311);
    const auto *shi_313 = buffer.data(shi + 313);
    const auto *shi_314 = buffer.data(shi + 314);
    const auto *shi_317 = buffer.data(shi + 317);
    const auto *shi_318 = buffer.data(shi + 318);
    const auto *shi_322 = buffer.data(shi + 322);
    const auto *shi_329 = buffer.data(shi + 329);
    const auto *shi_330 = buffer.data(shi + 330);
    const auto *shi_331 = buffer.data(shi + 331);
    const auto *shi_332 = buffer.data(shi + 332);
    const auto *shi_333 = buffer.data(shi + 333);
    const auto *shi_334 = buffer.data(shi + 334);
    const auto *shi_335 = buffer.data(shi + 335);
    const auto *shi_336 = buffer.data(shi + 336);
    const auto *shi_338 = buffer.data(shi + 338);
    const auto *shi_339 = buffer.data(shi + 339);
    const auto *shi_341 = buffer.data(shi + 341);
    const auto *shi_342 = buffer.data(shi + 342);
    const auto *shi_345 = buffer.data(shi + 345);
    const auto *shi_346 = buffer.data(shi + 346);
    const auto *shi_350 = buffer.data(shi + 350);
    const auto *shi_357 = buffer.data(shi + 357);
    const auto *shi_358 = buffer.data(shi + 358);
    const auto *shi_359 = buffer.data(shi + 359);
    const auto *shi_360 = buffer.data(shi + 360);
    const auto *shi_361 = buffer.data(shi + 361);
    const auto *shi_362 = buffer.data(shi + 362);
    const auto *shi_363 = buffer.data(shi + 363);
    const auto *shi_364 = buffer.data(shi + 364);

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, pc_x, sgi_275, sgi_276, sgi_277, \
                         sgi_278, sgi_279, shi_275, shi_276, shi_277, shi_278, \
                         shi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_14 * sgi_275[k]
                   + f_3 * pc_x[k] * shi_275[k];

        t_348[k] = f_14 * sgi_276[k]
                   + f_3 * pc_x[k] * shi_276[k];

        t_349[k] = f_14 * sgi_277[k]
                   + f_3 * pc_x[k] * shi_277[k];

        t_350[k] = f_14 * sgi_278[k]
                   + f_3 * pc_x[k] * shi_278[k];

        t_351[k] = f_14 * sgi_279[k]
                   + f_3 * pc_x[k] * shi_279[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pc_y, pc_z, sgi_161, shh0_204, shh0_206, \
                         shh0_207, shh1_204, shh1_206, shh1_207, shi_273, shi_275, \
                         shi_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_1 * shh0_204[k]
                   - f_2 * shh1_204[k]
                   + f_3 * pc_y[k] * shi_273[k];

        t_353[k] = f_15 * sgi_161[k]
                   + f_3 * pc_z[k] * shi_273[k];

        t_354[k] = f_4 * shh0_206[k]
                   - f_5 * shh1_206[k]
                   + f_3 * pc_y[k] * shi_275[k];

        t_355[k] = f_6 * shh0_207[k]
                   - f_7 * shh1_207[k]
                   + f_3 * pc_y[k] * shi_276[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pc_y, pc_z, sgi_167, shh0_208, shh0_209, \
                         shh1_208, shh1_209, shi_277, shi_278, \
                         shi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_8 * shh0_208[k]
                   - f_9 * shh1_208[k]
                   + f_3 * pc_y[k] * shi_277[k];

        t_357[k] = f_10 * shh0_209[k]
                   - f_11 * shh1_209[k]
                   + f_3 * pc_y[k] * shi_278[k];

        t_358[k] = f_3 * pc_y[k] * shi_279[k];

        t_359[k] = f_15 * sgi_167[k]
                   + f_1 * shh0_209[k]
                   - f_2 * shh1_209[k]
                   + f_3 * pc_z[k] * shi_279[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pb_x, pc_x, pc_y, pc_z, sgk0_360, \
                         sgk0_363, sgi_168, sgi_280, sgi_283, sgk1_360, sgk1_363, \
                         shi_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = pb_x[k] * sgk0_360[k]
                   + f_17 * sgi_280[k]
                   - f_12 * pc_x[k] * sgk1_360[k];

        t_361[k] = f_16 * sgi_168[k]
                   + f_3 * pc_y[k] * shi_280[k];

        t_362[k] = f_3 * pc_z[k] * shi_280[k];

        t_363[k] = pb_x[k] * sgk0_363[k]
                   + f_0 * sgi_283[k]
                   - f_12 * pc_x[k] * sgk1_363[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, pb_x, pc_x, pc_y, sgk0_365, sgk0_366, sgi_170, \
                         sgi_285, sgi_286, sgk1_365, sgk1_366, \
                         shi_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_16 * sgi_170[k]
                   + f_3 * pc_y[k] * shi_282[k];

        t_365[k] = pb_x[k] * sgk0_365[k]
                   + f_0 * sgi_285[k]
                   - f_12 * pc_x[k] * sgk1_365[k];

        t_366[k] = pb_x[k] * sgk0_366[k]
                   + f_16 * sgi_286[k]
                   - f_12 * pc_x[k] * sgk1_366[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, pb_x, pc_x, pc_y, pc_z, sgk0_369, sgi_173, \
                         sgi_289, sgk1_369, shi_283, shi_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_3 * pc_z[k] * shi_283[k];

        t_368[k] = f_16 * sgi_173[k]
                   + f_3 * pc_y[k] * shi_285[k];

        t_369[k] = pb_x[k] * sgk0_369[k]
                   + f_16 * sgi_289[k]
                   - f_12 * pc_x[k] * sgk1_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pb_x, pc_x, pc_z, sgk0_370, sgk0_372, sgi_290, \
                         sgi_292, sgk1_370, sgk1_372, shi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pb_x[k] * sgk0_370[k]
                   + f_15 * sgi_290[k]
                   - f_12 * pc_x[k] * sgk1_370[k];

        t_371[k] = f_3 * pc_z[k] * shi_286[k];

        t_372[k] = pb_x[k] * sgk0_372[k]
                   + f_15 * sgi_292[k]
                   - f_12 * pc_x[k] * sgk1_372[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pb_x, pc_x, pc_y, sgk0_374, sgk0_375, sgi_177, \
                         sgi_294, sgi_295, sgk1_374, sgk1_375, \
                         shi_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_16 * sgi_177[k]
                   + f_3 * pc_y[k] * shi_289[k];

        t_374[k] = pb_x[k] * sgk0_374[k]
                   + f_15 * sgi_294[k]
                   - f_12 * pc_x[k] * sgk1_374[k];

        t_375[k] = pb_x[k] * sgk0_375[k]
                   + f_14 * sgi_295[k]
                   - f_12 * pc_x[k] * sgk1_375[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pb_x, pc_x, pc_z, sgk0_377, sgk0_378, sgi_297, \
                         sgi_298, sgk1_377, sgk1_378, shi_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_3 * pc_z[k] * shi_290[k];

        t_377[k] = pb_x[k] * sgk0_377[k]
                   + f_14 * sgi_297[k]
                   - f_12 * pc_x[k] * sgk1_377[k];

        t_378[k] = pb_x[k] * sgk0_378[k]
                   + f_14 * sgi_298[k]
                   - f_12 * pc_x[k] * sgk1_378[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pb_x, pc_x, pc_y, sgk0_380, sgi_182, \
                         sgi_300, sgi_301, sgi_302, sgk1_380, shi_294, shi_301, \
                         shi_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_16 * sgi_182[k]
                   + f_3 * pc_y[k] * shi_294[k];

        t_380[k] = pb_x[k] * sgk0_380[k]
                   + f_14 * sgi_300[k]
                   - f_12 * pc_x[k] * sgk1_380[k];

        t_381[k] = f_13 * sgi_301[k]
                   + f_3 * pc_x[k] * shi_301[k];

        t_382[k] = f_13 * sgi_302[k]
                   + f_3 * pc_x[k] * shi_302[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, pc_x, sgi_303, sgi_304, sgi_305, \
                         sgi_306, sgi_307, shi_303, shi_304, shi_305, shi_306, \
                         shi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_13 * sgi_303[k]
                   + f_3 * pc_x[k] * shi_303[k];

        t_384[k] = f_13 * sgi_304[k]
                   + f_3 * pc_x[k] * shi_304[k];

        t_385[k] = f_13 * sgi_305[k]
                   + f_3 * pc_x[k] * shi_305[k];

        t_386[k] = f_13 * sgi_306[k]
                   + f_3 * pc_x[k] * shi_306[k];

        t_387[k] = f_13 * sgi_307[k]
                   + f_3 * pc_x[k] * shi_307[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pb_x, pc_x, pc_z, sgk0_388, sgk0_390, \
                         sgk0_391, sgk1_388, sgk1_390, sgk1_391, \
                         shi_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = pb_x[k] * sgk0_388[k]
                   - f_12 * pc_x[k] * sgk1_388[k];

        t_389[k] = f_3 * pc_z[k] * shi_301[k];

        t_390[k] = pb_x[k] * sgk0_390[k]
                   - f_12 * pc_x[k] * sgk1_390[k];

        t_391[k] = pb_x[k] * sgk0_391[k]
                   - f_12 * pc_x[k] * sgk1_391[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, pb_x, pc_x, pc_y, sgk0_392, sgk0_393, \
                         sgk0_395, sgi_195, sgk1_392, sgk1_393, sgk1_395, \
                         shi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = pb_x[k] * sgk0_392[k]
                   - f_12 * pc_x[k] * sgk1_392[k];

        t_393[k] = pb_x[k] * sgk0_393[k]
                   - f_12 * pc_x[k] * sgk1_393[k];

        t_394[k] = f_16 * sgi_195[k]
                   + f_3 * pc_y[k] * shi_307[k];

        t_395[k] = pb_x[k] * sgk0_395[k]
                   - f_12 * pc_x[k] * sgk1_395[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, pb_z, pc_y, pc_z, sgk0_216, sgk0_219, \
                         sgi_168, sgi_196, sgk1_216, sgk1_219, \
                         shi_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = pb_z[k] * sgk0_216[k]
                   - f_12 * pc_z[k] * sgk1_216[k];

        t_397[k] = f_15 * sgi_196[k]
                   + f_3 * pc_y[k] * shi_308[k];

        t_398[k] = f_13 * sgi_168[k]
                   + f_3 * pc_z[k] * shi_308[k];

        t_399[k] = pb_z[k] * sgk0_219[k]
                   - f_12 * pc_z[k] * sgk1_219[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pb_x, pb_z, pc_x, pc_y, pc_z, sgk0_222, \
                         sgk0_401, sgi_198, sgi_313, sgk1_222, sgk1_401, \
                         shi_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_15 * sgi_198[k]
                   + f_3 * pc_y[k] * shi_310[k];

        t_401[k] = pb_x[k] * sgk0_401[k]
                   + f_0 * sgi_313[k]
                   - f_12 * pc_x[k] * sgk1_401[k];

        t_402[k] = pb_z[k] * sgk0_222[k]
                   - f_12 * pc_z[k] * sgk1_222[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, pb_x, pc_x, pc_y, pc_z, sgk0_405, sgi_171, \
                         sgi_201, sgi_317, sgk1_405, shi_311, shi_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_13 * sgi_171[k]
                   + f_3 * pc_z[k] * shi_311[k];

        t_404[k] = f_15 * sgi_201[k]
                   + f_3 * pc_y[k] * shi_313[k];

        t_405[k] = pb_x[k] * sgk0_405[k]
                   + f_16 * sgi_317[k]
                   - f_12 * pc_x[k] * sgk1_405[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, pb_x, pb_z, pc_x, pc_z, sgk0_226, sgk0_408, \
                         sgi_174, sgi_320, sgk1_226, sgk1_408, \
                         shi_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = pb_z[k] * sgk0_226[k]
                   - f_12 * pc_z[k] * sgk1_226[k];

        t_407[k] = f_13 * sgi_174[k]
                   + f_3 * pc_z[k] * shi_314[k];

        t_408[k] = pb_x[k] * sgk0_408[k]
                   + f_15 * sgi_320[k]
                   - f_12 * pc_x[k] * sgk1_408[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, pb_x, pb_z, pc_x, pc_y, pc_z, sgk0_231, \
                         sgk0_410, sgi_205, sgi_322, sgk1_231, sgk1_410, \
                         shi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_15 * sgi_205[k]
                   + f_3 * pc_y[k] * shi_317[k];

        t_410[k] = pb_x[k] * sgk0_410[k]
                   + f_15 * sgi_322[k]
                   - f_12 * pc_x[k] * sgk1_410[k];

        t_411[k] = pb_z[k] * sgk0_231[k]
                   - f_12 * pc_z[k] * sgk1_231[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, pb_x, pc_x, pc_z, sgk0_413, sgk0_414, sgi_178, \
                         sgi_325, sgi_326, sgk1_413, sgk1_414, \
                         shi_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_13 * sgi_178[k]
                   + f_3 * pc_z[k] * shi_318[k];

        t_413[k] = pb_x[k] * sgk0_413[k]
                   + f_14 * sgi_325[k]
                   - f_12 * pc_x[k] * sgk1_413[k];

        t_414[k] = pb_x[k] * sgk0_414[k]
                   + f_14 * sgi_326[k]
                   - f_12 * pc_x[k] * sgk1_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, pb_x, pc_x, pc_y, sgk0_416, sgi_210, \
                         sgi_328, sgi_329, sgi_330, sgk1_416, shi_322, shi_329, \
                         shi_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_15 * sgi_210[k]
                   + f_3 * pc_y[k] * shi_322[k];

        t_416[k] = pb_x[k] * sgk0_416[k]
                   + f_14 * sgi_328[k]
                   - f_12 * pc_x[k] * sgk1_416[k];

        t_417[k] = f_13 * sgi_329[k]
                   + f_3 * pc_x[k] * shi_329[k];

        t_418[k] = f_13 * sgi_330[k]
                   + f_3 * pc_x[k] * shi_330[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, t_423, pc_x, sgi_331, sgi_332, sgi_333, \
                         sgi_334, sgi_335, shi_331, shi_332, shi_333, shi_334, \
                         shi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_13 * sgi_331[k]
                   + f_3 * pc_x[k] * shi_331[k];

        t_420[k] = f_13 * sgi_332[k]
                   + f_3 * pc_x[k] * shi_332[k];

        t_421[k] = f_13 * sgi_333[k]
                   + f_3 * pc_x[k] * shi_333[k];

        t_422[k] = f_13 * sgi_334[k]
                   + f_3 * pc_x[k] * shi_334[k];

        t_423[k] = f_13 * sgi_335[k]
                   + f_3 * pc_x[k] * shi_335[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pb_x, pc_x, pc_z, sgk0_424, sgk0_426, \
                         sgk0_427, sgi_189, sgk1_424, sgk1_426, sgk1_427, \
                         shi_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pb_x[k] * sgk0_424[k]
                   - f_12 * pc_x[k] * sgk1_424[k];

        t_425[k] = f_13 * sgi_189[k]
                   + f_3 * pc_z[k] * shi_329[k];

        t_426[k] = pb_x[k] * sgk0_426[k]
                   - f_12 * pc_x[k] * sgk1_426[k];

        t_427[k] = pb_x[k] * sgk0_427[k]
                   - f_12 * pc_x[k] * sgk1_427[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pb_x, pc_x, pc_y, sgk0_428, sgk0_429, \
                         sgk0_431, sgi_223, sgk1_428, sgk1_429, sgk1_431, \
                         shi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = pb_x[k] * sgk0_428[k]
                   - f_12 * pc_x[k] * sgk1_428[k];

        t_429[k] = pb_x[k] * sgk0_429[k]
                   - f_12 * pc_x[k] * sgk1_429[k];

        t_430[k] = f_15 * sgi_223[k]
                   + f_3 * pc_y[k] * shi_335[k];

        t_431[k] = pb_x[k] * sgk0_431[k]
                   - f_12 * pc_x[k] * sgk1_431[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, pb_x, pc_x, pc_y, pc_z, sgk0_432, sgi_196, \
                         sgi_224, sgi_336, sgk1_432, shi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = pb_x[k] * sgk0_432[k]
                   + f_17 * sgi_336[k]
                   - f_12 * pc_x[k] * sgk1_432[k];

        t_433[k] = f_14 * sgi_224[k]
                   + f_3 * pc_y[k] * shi_336[k];

        t_434[k] = f_14 * sgi_196[k]
                   + f_3 * pc_z[k] * shi_336[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, pb_x, pc_x, pc_y, sgk0_435, sgk0_437, sgi_226, \
                         sgi_339, sgi_341, sgk1_435, sgk1_437, \
                         shi_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = pb_x[k] * sgk0_435[k]
                   + f_0 * sgi_339[k]
                   - f_12 * pc_x[k] * sgk1_435[k];

        t_436[k] = f_14 * sgi_226[k]
                   + f_3 * pc_y[k] * shi_338[k];

        t_437[k] = pb_x[k] * sgk0_437[k]
                   + f_0 * sgi_341[k]
                   - f_12 * pc_x[k] * sgk1_437[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, pb_x, pc_x, pc_y, pc_z, sgk0_438, sgi_199, \
                         sgi_229, sgi_342, sgk1_438, shi_339, shi_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = pb_x[k] * sgk0_438[k]
                   + f_16 * sgi_342[k]
                   - f_12 * pc_x[k] * sgk1_438[k];

        t_439[k] = f_14 * sgi_199[k]
                   + f_3 * pc_z[k] * shi_339[k];

        t_440[k] = f_14 * sgi_229[k]
                   + f_3 * pc_y[k] * shi_341[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pb_x, pc_x, pc_z, sgk0_441, sgk0_442, sgi_202, \
                         sgi_345, sgi_346, sgk1_441, sgk1_442, \
                         shi_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = pb_x[k] * sgk0_441[k]
                   + f_16 * sgi_345[k]
                   - f_12 * pc_x[k] * sgk1_441[k];

        t_442[k] = pb_x[k] * sgk0_442[k]
                   + f_15 * sgi_346[k]
                   - f_12 * pc_x[k] * sgk1_442[k];

        t_443[k] = f_14 * sgi_202[k]
                   + f_3 * pc_z[k] * shi_342[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pb_x, pc_x, pc_y, sgk0_444, sgk0_446, sgi_233, \
                         sgi_348, sgi_350, sgk1_444, sgk1_446, \
                         shi_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = pb_x[k] * sgk0_444[k]
                   + f_15 * sgi_348[k]
                   - f_12 * pc_x[k] * sgk1_444[k];

        t_445[k] = f_14 * sgi_233[k]
                   + f_3 * pc_y[k] * shi_345[k];

        t_446[k] = pb_x[k] * sgk0_446[k]
                   + f_15 * sgi_350[k]
                   - f_12 * pc_x[k] * sgk1_446[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pb_x, pc_x, pc_z, sgk0_447, sgk0_449, sgi_206, \
                         sgi_351, sgi_353, sgk1_447, sgk1_449, \
                         shi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = pb_x[k] * sgk0_447[k]
                   + f_14 * sgi_351[k]
                   - f_12 * pc_x[k] * sgk1_447[k];

        t_448[k] = f_14 * sgi_206[k]
                   + f_3 * pc_z[k] * shi_346[k];

        t_449[k] = pb_x[k] * sgk0_449[k]
                   + f_14 * sgi_353[k]
                   - f_12 * pc_x[k] * sgk1_449[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, pb_x, pc_x, pc_y, sgk0_450, sgk0_452, sgi_238, \
                         sgi_354, sgi_356, sgk1_450, sgk1_452, \
                         shi_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = pb_x[k] * sgk0_450[k]
                   + f_14 * sgi_354[k]
                   - f_12 * pc_x[k] * sgk1_450[k];

        t_451[k] = f_14 * sgi_238[k]
                   + f_3 * pc_y[k] * shi_350[k];

        t_452[k] = pb_x[k] * sgk0_452[k]
                   + f_14 * sgi_356[k]
                   - f_12 * pc_x[k] * sgk1_452[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, t_457, pc_x, sgi_357, sgi_358, sgi_359, \
                         sgi_360, sgi_361, shi_357, shi_358, shi_359, shi_360, \
                         shi_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_13 * sgi_357[k]
                   + f_3 * pc_x[k] * shi_357[k];

        t_454[k] = f_13 * sgi_358[k]
                   + f_3 * pc_x[k] * shi_358[k];

        t_455[k] = f_13 * sgi_359[k]
                   + f_3 * pc_x[k] * shi_359[k];

        t_456[k] = f_13 * sgi_360[k]
                   + f_3 * pc_x[k] * shi_360[k];

        t_457[k] = f_13 * sgi_361[k]
                   + f_3 * pc_x[k] * shi_361[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pb_x, pc_x, pc_z, sgk0_460, sgi_217, \
                         sgi_362, sgi_363, sgk1_460, shi_357, shi_362, \
                         shi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = f_13 * sgi_362[k]
                   + f_3 * pc_x[k] * shi_362[k];

        t_459[k] = f_13 * sgi_363[k]
                   + f_3 * pc_x[k] * shi_363[k];

        t_460[k] = pb_x[k] * sgk0_460[k]
                   - f_12 * pc_x[k] * sgk1_460[k];

        t_461[k] = f_14 * sgi_217[k]
                   + f_3 * pc_z[k] * shi_357[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pb_x, pc_x, sgk0_462, sgk0_463, sgk0_464, \
                         sgk0_465, sgk1_462, sgk1_463, sgk1_464, \
                         sgk1_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = pb_x[k] * sgk0_462[k]
                   - f_12 * pc_x[k] * sgk1_462[k];

        t_463[k] = pb_x[k] * sgk0_463[k]
                   - f_12 * pc_x[k] * sgk1_463[k];

        t_464[k] = pb_x[k] * sgk0_464[k]
                   - f_12 * pc_x[k] * sgk1_464[k];

        t_465[k] = pb_x[k] * sgk0_465[k]
                   - f_12 * pc_x[k] * sgk1_465[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pb_x, pb_y, pc_x, pc_y, sgk0_324, \
                         sgk0_467, sgi_251, sgi_252, sgk1_324, sgk1_467, shi_363, \
                         shi_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * sgi_251[k]
                   + f_3 * pc_y[k] * shi_363[k];

        t_467[k] = pb_x[k] * sgk0_467[k]
                   - f_12 * pc_x[k] * sgk1_467[k];

        t_468[k] = pb_y[k] * sgk0_324[k]
                   - f_12 * pc_y[k] * sgk1_324[k];

        t_469[k] = f_13 * sgi_252[k]
                   + f_3 * pc_y[k] * shi_364[k];
    }
}

static auto
compute_prim_shk_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgk0,
                                                          const size_t sgi, const size_t sgk1,
                                                          const size_t shh0, const size_t shh1,
                                                          const size_t shi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.5 / q;

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
    auto *t_591 = buffer.data(target + 591);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgk0_329 = buffer.data(sgk0 + 329);
    const auto *sgk0_333 = buffer.data(sgk0 + 333);
    const auto *sgk0_338 = buffer.data(sgk0 + 338);
    const auto *sgk0_344 = buffer.data(sgk0 + 344);
    const auto *sgk0_360 = buffer.data(sgk0 + 360);
    const auto *sgk0_363 = buffer.data(sgk0 + 363);
    const auto *sgk0_366 = buffer.data(sgk0 + 366);
    const auto *sgk0_370 = buffer.data(sgk0 + 370);
    const auto *sgk0_375 = buffer.data(sgk0 + 375);
    const auto *sgk0_471 = buffer.data(sgk0 + 471);
    const auto *sgk0_474 = buffer.data(sgk0 + 474);
    const auto *sgk0_478 = buffer.data(sgk0 + 478);
    const auto *sgk0_480 = buffer.data(sgk0 + 480);
    const auto *sgk0_483 = buffer.data(sgk0 + 483);
    const auto *sgk0_485 = buffer.data(sgk0 + 485);
    const auto *sgk0_486 = buffer.data(sgk0 + 486);
    const auto *sgk0_496 = buffer.data(sgk0 + 496);
    const auto *sgk0_498 = buffer.data(sgk0 + 498);
    const auto *sgk0_499 = buffer.data(sgk0 + 499);
    const auto *sgk0_500 = buffer.data(sgk0 + 500);
    const auto *sgk0_501 = buffer.data(sgk0 + 501);
    const auto *sgk0_503 = buffer.data(sgk0 + 503);
    const auto *sgk0_504 = buffer.data(sgk0 + 504);
    const auto *sgk0_507 = buffer.data(sgk0 + 507);
    const auto *sgk0_509 = buffer.data(sgk0 + 509);
    const auto *sgk0_510 = buffer.data(sgk0 + 510);
    const auto *sgk0_513 = buffer.data(sgk0 + 513);
    const auto *sgk0_514 = buffer.data(sgk0 + 514);
    const auto *sgk0_516 = buffer.data(sgk0 + 516);
    const auto *sgk0_518 = buffer.data(sgk0 + 518);
    const auto *sgk0_519 = buffer.data(sgk0 + 519);
    const auto *sgk0_521 = buffer.data(sgk0 + 521);
    const auto *sgk0_522 = buffer.data(sgk0 + 522);
    const auto *sgk0_524 = buffer.data(sgk0 + 524);
    const auto *sgk0_532 = buffer.data(sgk0 + 532);
    const auto *sgk0_534 = buffer.data(sgk0 + 534);
    const auto *sgk0_535 = buffer.data(sgk0 + 535);
    const auto *sgk0_536 = buffer.data(sgk0 + 536);
    const auto *sgk0_537 = buffer.data(sgk0 + 537);
    const auto *sgk0_539 = buffer.data(sgk0 + 539);

    const auto *sgi_224 = buffer.data(sgi + 224);
    const auto *sgi_227 = buffer.data(sgi + 227);
    const auto *sgi_230 = buffer.data(sgi + 230);
    const auto *sgi_234 = buffer.data(sgi + 234);
    const auto *sgi_245 = buffer.data(sgi + 245);
    const auto *sgi_252 = buffer.data(sgi + 252);
    const auto *sgi_254 = buffer.data(sgi + 254);
    const auto *sgi_255 = buffer.data(sgi + 255);
    const auto *sgi_257 = buffer.data(sgi + 257);
    const auto *sgi_258 = buffer.data(sgi + 258);
    const auto *sgi_261 = buffer.data(sgi + 261);
    const auto *sgi_262 = buffer.data(sgi + 262);
    const auto *sgi_266 = buffer.data(sgi + 266);
    const auto *sgi_273 = buffer.data(sgi + 273);
    const auto *sgi_279 = buffer.data(sgi + 279);
    const auto *sgi_280 = buffer.data(sgi + 280);
    const auto *sgi_282 = buffer.data(sgi + 282);
    const auto *sgi_283 = buffer.data(sgi + 283);
    const auto *sgi_285 = buffer.data(sgi + 285);
    const auto *sgi_286 = buffer.data(sgi + 286);
    const auto *sgi_289 = buffer.data(sgi + 289);
    const auto *sgi_294 = buffer.data(sgi + 294);
    const auto *sgi_301 = buffer.data(sgi + 301);
    const auto *sgi_303 = buffer.data(sgi + 303);
    const auto *sgi_304 = buffer.data(sgi + 304);
    const auto *sgi_305 = buffer.data(sgi + 305);
    const auto *sgi_306 = buffer.data(sgi + 306);
    const auto *sgi_307 = buffer.data(sgi + 307);
    const auto *sgi_308 = buffer.data(sgi + 308);
    const auto *sgi_310 = buffer.data(sgi + 310);
    const auto *sgi_313 = buffer.data(sgi + 313);
    const auto *sgi_317 = buffer.data(sgi + 317);
    const auto *sgi_367 = buffer.data(sgi + 367);
    const auto *sgi_370 = buffer.data(sgi + 370);
    const auto *sgi_374 = buffer.data(sgi + 374);
    const auto *sgi_376 = buffer.data(sgi + 376);
    const auto *sgi_379 = buffer.data(sgi + 379);
    const auto *sgi_381 = buffer.data(sgi + 381);
    const auto *sgi_382 = buffer.data(sgi + 382);
    const auto *sgi_385 = buffer.data(sgi + 385);
    const auto *sgi_386 = buffer.data(sgi + 386);
    const auto *sgi_387 = buffer.data(sgi + 387);
    const auto *sgi_388 = buffer.data(sgi + 388);
    const auto *sgi_389 = buffer.data(sgi + 389);
    const auto *sgi_390 = buffer.data(sgi + 390);
    const auto *sgi_391 = buffer.data(sgi + 391);
    const auto *sgi_392 = buffer.data(sgi + 392);
    const auto *sgi_395 = buffer.data(sgi + 395);
    const auto *sgi_397 = buffer.data(sgi + 397);
    const auto *sgi_398 = buffer.data(sgi + 398);
    const auto *sgi_401 = buffer.data(sgi + 401);
    const auto *sgi_402 = buffer.data(sgi + 402);
    const auto *sgi_404 = buffer.data(sgi + 404);
    const auto *sgi_406 = buffer.data(sgi + 406);
    const auto *sgi_407 = buffer.data(sgi + 407);
    const auto *sgi_409 = buffer.data(sgi + 409);
    const auto *sgi_410 = buffer.data(sgi + 410);
    const auto *sgi_412 = buffer.data(sgi + 412);
    const auto *sgi_413 = buffer.data(sgi + 413);
    const auto *sgi_414 = buffer.data(sgi + 414);
    const auto *sgi_415 = buffer.data(sgi + 415);
    const auto *sgi_416 = buffer.data(sgi + 416);
    const auto *sgi_417 = buffer.data(sgi + 417);
    const auto *sgi_418 = buffer.data(sgi + 418);
    const auto *sgi_419 = buffer.data(sgi + 419);

    const auto *sgk1_329 = buffer.data(sgk1 + 329);
    const auto *sgk1_333 = buffer.data(sgk1 + 333);
    const auto *sgk1_338 = buffer.data(sgk1 + 338);
    const auto *sgk1_344 = buffer.data(sgk1 + 344);
    const auto *sgk1_360 = buffer.data(sgk1 + 360);
    const auto *sgk1_363 = buffer.data(sgk1 + 363);
    const auto *sgk1_366 = buffer.data(sgk1 + 366);
    const auto *sgk1_370 = buffer.data(sgk1 + 370);
    const auto *sgk1_375 = buffer.data(sgk1 + 375);
    const auto *sgk1_471 = buffer.data(sgk1 + 471);
    const auto *sgk1_474 = buffer.data(sgk1 + 474);
    const auto *sgk1_478 = buffer.data(sgk1 + 478);
    const auto *sgk1_480 = buffer.data(sgk1 + 480);
    const auto *sgk1_483 = buffer.data(sgk1 + 483);
    const auto *sgk1_485 = buffer.data(sgk1 + 485);
    const auto *sgk1_486 = buffer.data(sgk1 + 486);
    const auto *sgk1_496 = buffer.data(sgk1 + 496);
    const auto *sgk1_498 = buffer.data(sgk1 + 498);
    const auto *sgk1_499 = buffer.data(sgk1 + 499);
    const auto *sgk1_500 = buffer.data(sgk1 + 500);
    const auto *sgk1_501 = buffer.data(sgk1 + 501);
    const auto *sgk1_503 = buffer.data(sgk1 + 503);
    const auto *sgk1_504 = buffer.data(sgk1 + 504);
    const auto *sgk1_507 = buffer.data(sgk1 + 507);
    const auto *sgk1_509 = buffer.data(sgk1 + 509);
    const auto *sgk1_510 = buffer.data(sgk1 + 510);
    const auto *sgk1_513 = buffer.data(sgk1 + 513);
    const auto *sgk1_514 = buffer.data(sgk1 + 514);
    const auto *sgk1_516 = buffer.data(sgk1 + 516);
    const auto *sgk1_518 = buffer.data(sgk1 + 518);
    const auto *sgk1_519 = buffer.data(sgk1 + 519);
    const auto *sgk1_521 = buffer.data(sgk1 + 521);
    const auto *sgk1_522 = buffer.data(sgk1 + 522);
    const auto *sgk1_524 = buffer.data(sgk1 + 524);
    const auto *sgk1_532 = buffer.data(sgk1 + 532);
    const auto *sgk1_534 = buffer.data(sgk1 + 534);
    const auto *sgk1_535 = buffer.data(sgk1 + 535);
    const auto *sgk1_536 = buffer.data(sgk1 + 536);
    const auto *sgk1_537 = buffer.data(sgk1 + 537);
    const auto *sgk1_539 = buffer.data(sgk1 + 539);

    const auto *shh0_315 = buffer.data(shh0 + 315);
    const auto *shh0_318 = buffer.data(shh0 + 318);
    const auto *shh0_320 = buffer.data(shh0 + 320);
    const auto *shh0_321 = buffer.data(shh0 + 321);
    const auto *shh0_324 = buffer.data(shh0 + 324);
    const auto *shh0_325 = buffer.data(shh0 + 325);
    const auto *shh0_327 = buffer.data(shh0 + 327);
    const auto *shh0_329 = buffer.data(shh0 + 329);
    const auto *shh0_330 = buffer.data(shh0 + 330);
    const auto *shh0_332 = buffer.data(shh0 + 332);
    const auto *shh0_333 = buffer.data(shh0 + 333);
    const auto *shh0_334 = buffer.data(shh0 + 334);
    const auto *shh0_335 = buffer.data(shh0 + 335);
    const auto *shh0_341 = buffer.data(shh0 + 341);
    const auto *shh0_345 = buffer.data(shh0 + 345);
    const auto *shh0_348 = buffer.data(shh0 + 348);
    const auto *shh0_350 = buffer.data(shh0 + 350);

    const auto *shh1_315 = buffer.data(shh1 + 315);
    const auto *shh1_318 = buffer.data(shh1 + 318);
    const auto *shh1_320 = buffer.data(shh1 + 320);
    const auto *shh1_321 = buffer.data(shh1 + 321);
    const auto *shh1_324 = buffer.data(shh1 + 324);
    const auto *shh1_325 = buffer.data(shh1 + 325);
    const auto *shh1_327 = buffer.data(shh1 + 327);
    const auto *shh1_329 = buffer.data(shh1 + 329);
    const auto *shh1_330 = buffer.data(shh1 + 330);
    const auto *shh1_332 = buffer.data(shh1 + 332);
    const auto *shh1_333 = buffer.data(shh1 + 333);
    const auto *shh1_334 = buffer.data(shh1 + 334);
    const auto *shh1_335 = buffer.data(shh1 + 335);
    const auto *shh1_341 = buffer.data(shh1 + 341);
    const auto *shh1_345 = buffer.data(shh1 + 345);
    const auto *shh1_348 = buffer.data(shh1 + 348);
    const auto *shh1_350 = buffer.data(shh1 + 350);

    const auto *shi_364 = buffer.data(shi + 364);
    const auto *shi_366 = buffer.data(shi + 366);
    const auto *shi_367 = buffer.data(shi + 367);
    const auto *shi_369 = buffer.data(shi + 369);
    const auto *shi_370 = buffer.data(shi + 370);
    const auto *shi_373 = buffer.data(shi + 373);
    const auto *shi_374 = buffer.data(shi + 374);
    const auto *shi_378 = buffer.data(shi + 378);
    const auto *shi_385 = buffer.data(shi + 385);
    const auto *shi_386 = buffer.data(shi + 386);
    const auto *shi_387 = buffer.data(shi + 387);
    const auto *shi_388 = buffer.data(shi + 388);
    const auto *shi_389 = buffer.data(shi + 389);
    const auto *shi_390 = buffer.data(shi + 390);
    const auto *shi_391 = buffer.data(shi + 391);
    const auto *shi_392 = buffer.data(shi + 392);
    const auto *shi_394 = buffer.data(shi + 394);
    const auto *shi_395 = buffer.data(shi + 395);
    const auto *shi_397 = buffer.data(shi + 397);
    const auto *shi_398 = buffer.data(shi + 398);
    const auto *shi_401 = buffer.data(shi + 401);
    const auto *shi_402 = buffer.data(shi + 402);
    const auto *shi_406 = buffer.data(shi + 406);
    const auto *shi_413 = buffer.data(shi + 413);
    const auto *shi_414 = buffer.data(shi + 414);
    const auto *shi_415 = buffer.data(shi + 415);
    const auto *shi_416 = buffer.data(shi + 416);
    const auto *shi_417 = buffer.data(shi + 417);
    const auto *shi_418 = buffer.data(shi + 418);
    const auto *shi_419 = buffer.data(shi + 419);
    const auto *shi_420 = buffer.data(shi + 420);
    const auto *shi_422 = buffer.data(shi + 422);
    const auto *shi_423 = buffer.data(shi + 423);
    const auto *shi_425 = buffer.data(shi + 425);
    const auto *shi_426 = buffer.data(shi + 426);
    const auto *shi_429 = buffer.data(shi + 429);
    const auto *shi_430 = buffer.data(shi + 430);
    const auto *shi_432 = buffer.data(shi + 432);
    const auto *shi_434 = buffer.data(shi + 434);
    const auto *shi_435 = buffer.data(shi + 435);
    const auto *shi_437 = buffer.data(shi + 437);
    const auto *shi_438 = buffer.data(shi + 438);
    const auto *shi_440 = buffer.data(shi + 440);
    const auto *shi_441 = buffer.data(shi + 441);
    const auto *shi_442 = buffer.data(shi + 442);
    const auto *shi_443 = buffer.data(shi + 443);
    const auto *shi_444 = buffer.data(shi + 444);
    const auto *shi_445 = buffer.data(shi + 445);
    const auto *shi_446 = buffer.data(shi + 446);
    const auto *shi_447 = buffer.data(shi + 447);
    const auto *shi_448 = buffer.data(shi + 448);
    const auto *shi_450 = buffer.data(shi + 450);
    const auto *shi_451 = buffer.data(shi + 451);
    const auto *shi_453 = buffer.data(shi + 453);
    const auto *shi_454 = buffer.data(shi + 454);
    const auto *shi_457 = buffer.data(shi + 457);
    const auto *shi_460 = buffer.data(shi + 460);
    const auto *shi_462 = buffer.data(shi + 462);

#pragma omp simd aligned(t_470, t_471, t_472, pb_x, pc_x, pc_y, pc_z, sgk0_471, sgi_224, \
                         sgi_254, sgi_367, sgk1_471, shi_364, shi_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * sgi_224[k]
                   + f_3 * pc_z[k] * shi_364[k];

        t_471[k] = pb_x[k] * sgk0_471[k]
                   + f_0 * sgi_367[k]
                   - f_12 * pc_x[k] * sgk1_471[k];

        t_472[k] = f_13 * sgi_254[k]
                   + f_3 * pc_y[k] * shi_366[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pb_x, pb_y, pc_x, pc_y, pc_z, sgk0_329, \
                         sgk0_474, sgi_227, sgi_370, sgk1_329, sgk1_474, \
                         shi_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = pb_y[k] * sgk0_329[k]
                   - f_12 * pc_y[k] * sgk1_329[k];

        t_474[k] = pb_x[k] * sgk0_474[k]
                   + f_16 * sgi_370[k]
                   - f_12 * pc_x[k] * sgk1_474[k];

        t_475[k] = f_15 * sgi_227[k]
                   + f_3 * pc_z[k] * shi_367[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pb_x, pb_y, pc_x, pc_y, sgk0_333, sgk0_478, \
                         sgi_257, sgi_374, sgk1_333, sgk1_478, \
                         shi_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_13 * sgi_257[k]
                   + f_3 * pc_y[k] * shi_369[k];

        t_477[k] = pb_y[k] * sgk0_333[k]
                   - f_12 * pc_y[k] * sgk1_333[k];

        t_478[k] = pb_x[k] * sgk0_478[k]
                   + f_15 * sgi_374[k]
                   - f_12 * pc_x[k] * sgk1_478[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pb_x, pc_x, pc_y, pc_z, sgk0_480, sgi_230, \
                         sgi_261, sgi_376, sgk1_480, shi_370, shi_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_15 * sgi_230[k]
                   + f_3 * pc_z[k] * shi_370[k];

        t_480[k] = pb_x[k] * sgk0_480[k]
                   + f_15 * sgi_376[k]
                   - f_12 * pc_x[k] * sgk1_480[k];

        t_481[k] = f_13 * sgi_261[k]
                   + f_3 * pc_y[k] * shi_373[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pb_x, pb_y, pc_x, pc_y, pc_z, sgk0_338, \
                         sgk0_483, sgi_234, sgi_379, sgk1_338, sgk1_483, \
                         shi_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = pb_y[k] * sgk0_338[k]
                   - f_12 * pc_y[k] * sgk1_338[k];

        t_483[k] = pb_x[k] * sgk0_483[k]
                   + f_14 * sgi_379[k]
                   - f_12 * pc_x[k] * sgk1_483[k];

        t_484[k] = f_15 * sgi_234[k]
                   + f_3 * pc_z[k] * shi_374[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pb_x, pc_x, pc_y, sgk0_485, sgk0_486, sgi_266, \
                         sgi_381, sgi_382, sgk1_485, sgk1_486, \
                         shi_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pb_x[k] * sgk0_485[k]
                   + f_14 * sgi_381[k]
                   - f_12 * pc_x[k] * sgk1_485[k];

        t_486[k] = pb_x[k] * sgk0_486[k]
                   + f_14 * sgi_382[k]
                   - f_12 * pc_x[k] * sgk1_486[k];

        t_487[k] = f_13 * sgi_266[k]
                   + f_3 * pc_y[k] * shi_378[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, pb_y, pc_x, pc_y, sgk0_344, sgi_385, \
                         sgi_386, sgi_387, sgk1_344, shi_385, shi_386, \
                         shi_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = pb_y[k] * sgk0_344[k]
                   - f_12 * pc_y[k] * sgk1_344[k];

        t_489[k] = f_13 * sgi_385[k]
                   + f_3 * pc_x[k] * shi_385[k];

        t_490[k] = f_13 * sgi_386[k]
                   + f_3 * pc_x[k] * shi_386[k];

        t_491[k] = f_13 * sgi_387[k]
                   + f_3 * pc_x[k] * shi_387[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pc_x, sgi_388, sgi_389, sgi_390, sgi_391, \
                         shi_388, shi_389, shi_390, shi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_13 * sgi_388[k]
                   + f_3 * pc_x[k] * shi_388[k];

        t_493[k] = f_13 * sgi_389[k]
                   + f_3 * pc_x[k] * shi_389[k];

        t_494[k] = f_13 * sgi_390[k]
                   + f_3 * pc_x[k] * shi_390[k];

        t_495[k] = f_13 * sgi_391[k]
                   + f_3 * pc_x[k] * shi_391[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pb_x, pc_x, pc_z, sgk0_496, sgk0_498, \
                         sgk0_499, sgi_245, sgk1_496, sgk1_498, sgk1_499, \
                         shi_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = pb_x[k] * sgk0_496[k]
                   - f_12 * pc_x[k] * sgk1_496[k];

        t_497[k] = f_15 * sgi_245[k]
                   + f_3 * pc_z[k] * shi_385[k];

        t_498[k] = pb_x[k] * sgk0_498[k]
                   - f_12 * pc_x[k] * sgk1_498[k];

        t_499[k] = pb_x[k] * sgk0_499[k]
                   - f_12 * pc_x[k] * sgk1_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pb_x, pc_x, pc_y, sgk0_500, sgk0_501, \
                         sgk0_503, sgi_279, sgk1_500, sgk1_501, sgk1_503, \
                         shi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = pb_x[k] * sgk0_500[k]
                   - f_12 * pc_x[k] * sgk1_500[k];

        t_501[k] = pb_x[k] * sgk0_501[k]
                   - f_12 * pc_x[k] * sgk1_501[k];

        t_502[k] = f_13 * sgi_279[k]
                   + f_3 * pc_y[k] * shi_391[k];

        t_503[k] = pb_x[k] * sgk0_503[k]
                   - f_12 * pc_x[k] * sgk1_503[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pb_x, pc_x, pc_y, pc_z, sgk0_504, \
                         sgk0_507, sgi_252, sgi_392, sgi_395, sgk1_504, sgk1_507, \
                         shi_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = pb_x[k] * sgk0_504[k]
                   + f_17 * sgi_392[k]
                   - f_12 * pc_x[k] * sgk1_504[k];

        t_505[k] = f_3 * pc_y[k] * shi_392[k];

        t_506[k] = f_16 * sgi_252[k]
                   + f_3 * pc_z[k] * shi_392[k];

        t_507[k] = pb_x[k] * sgk0_507[k]
                   + f_0 * sgi_395[k]
                   - f_12 * pc_x[k] * sgk1_507[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pb_x, pc_x, pc_y, sgk0_509, sgk0_510, sgi_397, \
                         sgi_398, sgk1_509, sgk1_510, shi_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_3 * pc_y[k] * shi_394[k];

        t_509[k] = pb_x[k] * sgk0_509[k]
                   + f_0 * sgi_397[k]
                   - f_12 * pc_x[k] * sgk1_509[k];

        t_510[k] = pb_x[k] * sgk0_510[k]
                   + f_16 * sgi_398[k]
                   - f_12 * pc_x[k] * sgk1_510[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pb_x, pc_x, pc_y, pc_z, sgk0_513, sgi_255, \
                         sgi_401, sgk1_513, shi_395, shi_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_16 * sgi_255[k]
                   + f_3 * pc_z[k] * shi_395[k];

        t_512[k] = f_3 * pc_y[k] * shi_397[k];

        t_513[k] = pb_x[k] * sgk0_513[k]
                   + f_16 * sgi_401[k]
                   - f_12 * pc_x[k] * sgk1_513[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, pb_x, pc_x, pc_z, sgk0_514, sgk0_516, sgi_258, \
                         sgi_402, sgi_404, sgk1_514, sgk1_516, \
                         shi_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = pb_x[k] * sgk0_514[k]
                   + f_15 * sgi_402[k]
                   - f_12 * pc_x[k] * sgk1_514[k];

        t_515[k] = f_16 * sgi_258[k]
                   + f_3 * pc_z[k] * shi_398[k];

        t_516[k] = pb_x[k] * sgk0_516[k]
                   + f_15 * sgi_404[k]
                   - f_12 * pc_x[k] * sgk1_516[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, pb_x, pc_x, pc_y, sgk0_518, sgk0_519, sgi_406, \
                         sgi_407, sgk1_518, sgk1_519, shi_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_3 * pc_y[k] * shi_401[k];

        t_518[k] = pb_x[k] * sgk0_518[k]
                   + f_15 * sgi_406[k]
                   - f_12 * pc_x[k] * sgk1_518[k];

        t_519[k] = pb_x[k] * sgk0_519[k]
                   + f_14 * sgi_407[k]
                   - f_12 * pc_x[k] * sgk1_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pb_x, pc_x, pc_z, sgk0_521, sgk0_522, sgi_262, \
                         sgi_409, sgi_410, sgk1_521, sgk1_522, \
                         shi_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_16 * sgi_262[k]
                   + f_3 * pc_z[k] * shi_402[k];

        t_521[k] = pb_x[k] * sgk0_521[k]
                   + f_14 * sgi_409[k]
                   - f_12 * pc_x[k] * sgk1_521[k];

        t_522[k] = pb_x[k] * sgk0_522[k]
                   + f_14 * sgi_410[k]
                   - f_12 * pc_x[k] * sgk1_522[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pb_x, pc_x, pc_y, sgk0_524, sgi_412, \
                         sgi_413, sgi_414, sgk1_524, shi_406, shi_413, \
                         shi_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_3 * pc_y[k] * shi_406[k];

        t_524[k] = pb_x[k] * sgk0_524[k]
                   + f_14 * sgi_412[k]
                   - f_12 * pc_x[k] * sgk1_524[k];

        t_525[k] = f_13 * sgi_413[k]
                   + f_3 * pc_x[k] * shi_413[k];

        t_526[k] = f_13 * sgi_414[k]
                   + f_3 * pc_x[k] * shi_414[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, pc_x, sgi_415, sgi_416, sgi_417, \
                         sgi_418, sgi_419, shi_415, shi_416, shi_417, shi_418, \
                         shi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_13 * sgi_415[k]
                   + f_3 * pc_x[k] * shi_415[k];

        t_528[k] = f_13 * sgi_416[k]
                   + f_3 * pc_x[k] * shi_416[k];

        t_529[k] = f_13 * sgi_417[k]
                   + f_3 * pc_x[k] * shi_417[k];

        t_530[k] = f_13 * sgi_418[k]
                   + f_3 * pc_x[k] * shi_418[k];

        t_531[k] = f_13 * sgi_419[k]
                   + f_3 * pc_x[k] * shi_419[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pb_x, pc_x, pc_z, sgk0_532, sgk0_534, \
                         sgk0_535, sgi_273, sgk1_532, sgk1_534, sgk1_535, \
                         shi_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = pb_x[k] * sgk0_532[k]
                   - f_12 * pc_x[k] * sgk1_532[k];

        t_533[k] = f_16 * sgi_273[k]
                   + f_3 * pc_z[k] * shi_413[k];

        t_534[k] = pb_x[k] * sgk0_534[k]
                   - f_12 * pc_x[k] * sgk1_534[k];

        t_535[k] = pb_x[k] * sgk0_535[k]
                   - f_12 * pc_x[k] * sgk1_535[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pb_x, pc_x, pc_y, sgk0_536, sgk0_537, \
                         sgk0_539, sgk1_536, sgk1_537, sgk1_539, \
                         shi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = pb_x[k] * sgk0_536[k]
                   - f_12 * pc_x[k] * sgk1_536[k];

        t_537[k] = pb_x[k] * sgk0_537[k]
                   - f_12 * pc_x[k] * sgk1_537[k];

        t_538[k] = f_3 * pc_y[k] * shi_419[k];

        t_539[k] = pb_x[k] * sgk0_539[k]
                   - f_12 * pc_x[k] * sgk1_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, pc_x, pc_y, pc_z, sgi_280, shh0_315, \
                         shh0_318, shh1_315, shh1_318, shi_420, \
                         shi_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_1 * shh0_315[k]
                   - f_2 * shh1_315[k]
                   + f_3 * pc_x[k] * shi_420[k];

        t_541[k] = f_0 * sgi_280[k]
                   + f_3 * pc_y[k] * shi_420[k];

        t_542[k] = f_3 * pc_z[k] * shi_420[k];

        t_543[k] = f_4 * shh0_318[k]
                   - f_5 * shh1_318[k]
                   + f_3 * pc_x[k] * shi_423[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pc_x, pc_y, pc_z, sgi_282, shh0_320, \
                         shh0_321, shh1_320, shh1_321, shi_422, shi_423, shi_425, \
                         shi_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_0 * sgi_282[k]
                   + f_3 * pc_y[k] * shi_422[k];

        t_545[k] = f_4 * shh0_320[k]
                   - f_5 * shh1_320[k]
                   + f_3 * pc_x[k] * shi_425[k];

        t_546[k] = f_6 * shh0_321[k]
                   - f_7 * shh1_321[k]
                   + f_3 * pc_x[k] * shi_426[k];

        t_547[k] = f_3 * pc_z[k] * shi_423[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pc_x, pc_y, pc_z, sgi_285, shh0_324, \
                         shh0_325, shh1_324, shh1_325, shi_425, shi_426, shi_429, \
                         shi_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_0 * sgi_285[k]
                   + f_3 * pc_y[k] * shi_425[k];

        t_549[k] = f_6 * shh0_324[k]
                   - f_7 * shh1_324[k]
                   + f_3 * pc_x[k] * shi_429[k];

        t_550[k] = f_8 * shh0_325[k]
                   - f_9 * shh1_325[k]
                   + f_3 * pc_x[k] * shi_430[k];

        t_551[k] = f_3 * pc_z[k] * shi_426[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pc_x, pc_y, sgi_289, shh0_327, shh0_329, \
                         shh1_327, shh1_329, shi_429, shi_432, \
                         shi_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_8 * shh0_327[k]
                   - f_9 * shh1_327[k]
                   + f_3 * pc_x[k] * shi_432[k];

        t_553[k] = f_0 * sgi_289[k]
                   + f_3 * pc_y[k] * shi_429[k];

        t_554[k] = f_8 * shh0_329[k]
                   - f_9 * shh1_329[k]
                   + f_3 * pc_x[k] * shi_434[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, pc_x, pc_z, shh0_330, shh0_332, shh0_333, \
                         shh1_330, shh1_332, shh1_333, shi_430, shi_435, shi_437, \
                         shi_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_10 * shh0_330[k]
                   - f_11 * shh1_330[k]
                   + f_3 * pc_x[k] * shi_435[k];

        t_556[k] = f_3 * pc_z[k] * shi_430[k];

        t_557[k] = f_10 * shh0_332[k]
                   - f_11 * shh1_332[k]
                   + f_3 * pc_x[k] * shi_437[k];

        t_558[k] = f_10 * shh0_333[k]
                   - f_11 * shh1_333[k]
                   + f_3 * pc_x[k] * shi_438[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, t_563, pc_x, pc_y, sgi_294, shh0_335, \
                         shh1_335, shi_434, shi_440, shi_441, shi_442, \
                         shi_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_0 * sgi_294[k]
                   + f_3 * pc_y[k] * shi_434[k];

        t_560[k] = f_10 * shh0_335[k]
                   - f_11 * shh1_335[k]
                   + f_3 * pc_x[k] * shi_440[k];

        t_561[k] = f_3 * pc_x[k] * shi_441[k];

        t_562[k] = f_3 * pc_x[k] * shi_442[k];

        t_563[k] = f_3 * pc_x[k] * shi_443[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, t_567, t_568, pc_x, pc_y, sgi_301, shh0_330, \
                         shh1_330, shi_441, shi_444, shi_445, shi_446, \
                         shi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_3 * pc_x[k] * shi_444[k];

        t_565[k] = f_3 * pc_x[k] * shi_445[k];

        t_566[k] = f_3 * pc_x[k] * shi_446[k];

        t_567[k] = f_3 * pc_x[k] * shi_447[k];

        t_568[k] = f_0 * sgi_301[k]
                   + f_1 * shh0_330[k]
                   - f_2 * shh1_330[k]
                   + f_3 * pc_y[k] * shi_441[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_y, pc_z, sgi_303, sgi_304, shh0_332, \
                         shh0_333, shh1_332, shh1_333, shi_441, shi_443, \
                         shi_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_3 * pc_z[k] * shi_441[k];

        t_570[k] = f_0 * sgi_303[k]
                   + f_4 * shh0_332[k]
                   - f_5 * shh1_332[k]
                   + f_3 * pc_y[k] * shi_443[k];

        t_571[k] = f_0 * sgi_304[k]
                   + f_6 * shh0_333[k]
                   - f_7 * shh1_333[k]
                   + f_3 * pc_y[k] * shi_444[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_y, pc_z, sgi_305, sgi_306, sgi_307, \
                         shh0_334, shh0_335, shh1_334, shh1_335, shi_445, shi_446, \
                         shi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_0 * sgi_305[k]
                   + f_8 * shh0_334[k]
                   - f_9 * shh1_334[k]
                   + f_3 * pc_y[k] * shi_445[k];

        t_573[k] = f_0 * sgi_306[k]
                   + f_10 * shh0_335[k]
                   - f_11 * shh1_335[k]
                   + f_3 * pc_y[k] * shi_446[k];

        t_574[k] = f_0 * sgi_307[k]
                   + f_3 * pc_y[k] * shi_447[k];

        t_575[k] = f_1 * shh0_335[k]
                   - f_2 * shh1_335[k]
                   + f_3 * pc_z[k] * shi_447[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, pb_z, pc_y, pc_z, sgk0_360, sgk0_363, \
                         sgi_280, sgi_308, sgk1_360, sgk1_363, \
                         shi_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = pb_z[k] * sgk0_360[k]
                   - f_12 * pc_z[k] * sgk1_360[k];

        t_577[k] = f_16 * sgi_308[k]
                   + f_3 * pc_y[k] * shi_448[k];

        t_578[k] = f_13 * sgi_280[k]
                   + f_3 * pc_z[k] * shi_448[k];

        t_579[k] = pb_z[k] * sgk0_363[k]
                   - f_12 * pc_z[k] * sgk1_363[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, pb_z, pc_x, pc_y, pc_z, sgk0_366, sgi_310, \
                         sgk1_366, shh0_341, shh1_341, shi_450, \
                         shi_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_16 * sgi_310[k]
                   + f_3 * pc_y[k] * shi_450[k];

        t_581[k] = f_4 * shh0_341[k]
                   - f_5 * shh1_341[k]
                   + f_3 * pc_x[k] * shi_453[k];

        t_582[k] = pb_z[k] * sgk0_366[k]
                   - f_12 * pc_z[k] * sgk1_366[k];
    }

#pragma omp simd aligned(t_583, t_584, t_585, pc_x, pc_y, pc_z, sgi_283, sgi_313, shh0_345, \
                         shh1_345, shi_451, shi_453, shi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_583[k] = f_13 * sgi_283[k]
                   + f_3 * pc_z[k] * shi_451[k];

        t_584[k] = f_16 * sgi_313[k]
                   + f_3 * pc_y[k] * shi_453[k];

        t_585[k] = f_6 * shh0_345[k]
                   - f_7 * shh1_345[k]
                   + f_3 * pc_x[k] * shi_457[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, pb_z, pc_x, pc_z, sgk0_370, sgi_286, sgk1_370, \
                         shh0_348, shh1_348, shi_454, shi_460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = pb_z[k] * sgk0_370[k]
                   - f_12 * pc_z[k] * sgk1_370[k];

        t_587[k] = f_13 * sgi_286[k]
                   + f_3 * pc_z[k] * shi_454[k];

        t_588[k] = f_8 * shh0_348[k]
                   - f_9 * shh1_348[k]
                   + f_3 * pc_x[k] * shi_460[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, pb_z, pc_x, pc_y, pc_z, sgk0_375, sgi_317, \
                         sgk1_375, shh0_350, shh1_350, shi_457, \
                         shi_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = f_16 * sgi_317[k]
                   + f_3 * pc_y[k] * shi_457[k];

        t_590[k] = f_8 * shh0_350[k]
                   - f_9 * shh1_350[k]
                   + f_3 * pc_x[k] * shi_462[k];

        t_591[k] = pb_z[k] * sgk0_375[k]
                   - f_12 * pc_z[k] * sgk1_375[k];
    }
}

static auto
compute_prim_shk_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgk0,
                                                          const size_t sgi, const size_t sgk1,
                                                          const size_t shh0, const size_t shh1,
                                                          const size_t shi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgk0_388 = buffer.data(sgk0 + 388);
    const auto *sgk0_390 = buffer.data(sgk0 + 390);
    const auto *sgk0_391 = buffer.data(sgk0 + 391);
    const auto *sgk0_392 = buffer.data(sgk0 + 392);
    const auto *sgk0_393 = buffer.data(sgk0 + 393);
    const auto *sgk0_504 = buffer.data(sgk0 + 504);
    const auto *sgk0_509 = buffer.data(sgk0 + 509);
    const auto *sgk0_513 = buffer.data(sgk0 + 513);
    const auto *sgk0_518 = buffer.data(sgk0 + 518);
    const auto *sgk0_524 = buffer.data(sgk0 + 524);
    const auto *sgk0_532 = buffer.data(sgk0 + 532);

    const auto *sgi_290 = buffer.data(sgi + 290);
    const auto *sgi_301 = buffer.data(sgi + 301);
    const auto *sgi_302 = buffer.data(sgi + 302);
    const auto *sgi_303 = buffer.data(sgi + 303);
    const auto *sgi_304 = buffer.data(sgi + 304);
    const auto *sgi_305 = buffer.data(sgi + 305);
    const auto *sgi_307 = buffer.data(sgi + 307);
    const auto *sgi_308 = buffer.data(sgi + 308);
    const auto *sgi_311 = buffer.data(sgi + 311);
    const auto *sgi_314 = buffer.data(sgi + 314);
    const auto *sgi_318 = buffer.data(sgi + 318);
    const auto *sgi_322 = buffer.data(sgi + 322);
    const auto *sgi_329 = buffer.data(sgi + 329);
    const auto *sgi_335 = buffer.data(sgi + 335);
    const auto *sgi_336 = buffer.data(sgi + 336);
    const auto *sgi_338 = buffer.data(sgi + 338);
    const auto *sgi_339 = buffer.data(sgi + 339);
    const auto *sgi_341 = buffer.data(sgi + 341);
    const auto *sgi_342 = buffer.data(sgi + 342);
    const auto *sgi_345 = buffer.data(sgi + 345);
    const auto *sgi_346 = buffer.data(sgi + 346);
    const auto *sgi_350 = buffer.data(sgi + 350);
    const auto *sgi_357 = buffer.data(sgi + 357);
    const auto *sgi_359 = buffer.data(sgi + 359);
    const auto *sgi_360 = buffer.data(sgi + 360);
    const auto *sgi_361 = buffer.data(sgi + 361);
    const auto *sgi_362 = buffer.data(sgi + 362);
    const auto *sgi_363 = buffer.data(sgi + 363);
    const auto *sgi_364 = buffer.data(sgi + 364);
    const auto *sgi_366 = buffer.data(sgi + 366);
    const auto *sgi_367 = buffer.data(sgi + 367);
    const auto *sgi_369 = buffer.data(sgi + 369);
    const auto *sgi_370 = buffer.data(sgi + 370);
    const auto *sgi_373 = buffer.data(sgi + 373);
    const auto *sgi_374 = buffer.data(sgi + 374);
    const auto *sgi_378 = buffer.data(sgi + 378);
    const auto *sgi_385 = buffer.data(sgi + 385);
    const auto *sgi_387 = buffer.data(sgi + 387);
    const auto *sgi_388 = buffer.data(sgi + 388);
    const auto *sgi_389 = buffer.data(sgi + 389);
    const auto *sgi_390 = buffer.data(sgi + 390);
    const auto *sgi_391 = buffer.data(sgi + 391);
    const auto *sgi_392 = buffer.data(sgi + 392);
    const auto *sgi_394 = buffer.data(sgi + 394);
    const auto *sgi_397 = buffer.data(sgi + 397);
    const auto *sgi_401 = buffer.data(sgi + 401);
    const auto *sgi_406 = buffer.data(sgi + 406);
    const auto *sgi_413 = buffer.data(sgi + 413);

    const auto *sgk1_388 = buffer.data(sgk1 + 388);
    const auto *sgk1_390 = buffer.data(sgk1 + 390);
    const auto *sgk1_391 = buffer.data(sgk1 + 391);
    const auto *sgk1_392 = buffer.data(sgk1 + 392);
    const auto *sgk1_393 = buffer.data(sgk1 + 393);
    const auto *sgk1_504 = buffer.data(sgk1 + 504);
    const auto *sgk1_509 = buffer.data(sgk1 + 509);
    const auto *sgk1_513 = buffer.data(sgk1 + 513);
    const auto *sgk1_518 = buffer.data(sgk1 + 518);
    const auto *sgk1_524 = buffer.data(sgk1 + 524);
    const auto *sgk1_532 = buffer.data(sgk1 + 532);

    const auto *shh0_353 = buffer.data(shh0 + 353);
    const auto *shh0_354 = buffer.data(shh0 + 354);
    const auto *shh0_356 = buffer.data(shh0 + 356);
    const auto *shh0_357 = buffer.data(shh0 + 357);
    const auto *shh0_360 = buffer.data(shh0 + 360);
    const auto *shh0_362 = buffer.data(shh0 + 362);
    const auto *shh0_363 = buffer.data(shh0 + 363);
    const auto *shh0_366 = buffer.data(shh0 + 366);
    const auto *shh0_367 = buffer.data(shh0 + 367);
    const auto *shh0_369 = buffer.data(shh0 + 369);
    const auto *shh0_371 = buffer.data(shh0 + 371);
    const auto *shh0_372 = buffer.data(shh0 + 372);
    const auto *shh0_374 = buffer.data(shh0 + 374);
    const auto *shh0_375 = buffer.data(shh0 + 375);
    const auto *shh0_376 = buffer.data(shh0 + 376);
    const auto *shh0_377 = buffer.data(shh0 + 377);
    const auto *shh0_378 = buffer.data(shh0 + 378);
    const auto *shh0_381 = buffer.data(shh0 + 381);
    const auto *shh0_383 = buffer.data(shh0 + 383);
    const auto *shh0_384 = buffer.data(shh0 + 384);
    const auto *shh0_387 = buffer.data(shh0 + 387);
    const auto *shh0_388 = buffer.data(shh0 + 388);
    const auto *shh0_390 = buffer.data(shh0 + 390);
    const auto *shh0_392 = buffer.data(shh0 + 392);
    const auto *shh0_393 = buffer.data(shh0 + 393);
    const auto *shh0_395 = buffer.data(shh0 + 395);
    const auto *shh0_396 = buffer.data(shh0 + 396);
    const auto *shh0_397 = buffer.data(shh0 + 397);
    const auto *shh0_398 = buffer.data(shh0 + 398);
    const auto *shh0_402 = buffer.data(shh0 + 402);
    const auto *shh0_405 = buffer.data(shh0 + 405);
    const auto *shh0_409 = buffer.data(shh0 + 409);
    const auto *shh0_411 = buffer.data(shh0 + 411);
    const auto *shh0_414 = buffer.data(shh0 + 414);
    const auto *shh0_416 = buffer.data(shh0 + 416);
    const auto *shh0_417 = buffer.data(shh0 + 417);

    const auto *shh1_353 = buffer.data(shh1 + 353);
    const auto *shh1_354 = buffer.data(shh1 + 354);
    const auto *shh1_356 = buffer.data(shh1 + 356);
    const auto *shh1_357 = buffer.data(shh1 + 357);
    const auto *shh1_360 = buffer.data(shh1 + 360);
    const auto *shh1_362 = buffer.data(shh1 + 362);
    const auto *shh1_363 = buffer.data(shh1 + 363);
    const auto *shh1_366 = buffer.data(shh1 + 366);
    const auto *shh1_367 = buffer.data(shh1 + 367);
    const auto *shh1_369 = buffer.data(shh1 + 369);
    const auto *shh1_371 = buffer.data(shh1 + 371);
    const auto *shh1_372 = buffer.data(shh1 + 372);
    const auto *shh1_374 = buffer.data(shh1 + 374);
    const auto *shh1_375 = buffer.data(shh1 + 375);
    const auto *shh1_376 = buffer.data(shh1 + 376);
    const auto *shh1_377 = buffer.data(shh1 + 377);
    const auto *shh1_378 = buffer.data(shh1 + 378);
    const auto *shh1_381 = buffer.data(shh1 + 381);
    const auto *shh1_383 = buffer.data(shh1 + 383);
    const auto *shh1_384 = buffer.data(shh1 + 384);
    const auto *shh1_387 = buffer.data(shh1 + 387);
    const auto *shh1_388 = buffer.data(shh1 + 388);
    const auto *shh1_390 = buffer.data(shh1 + 390);
    const auto *shh1_392 = buffer.data(shh1 + 392);
    const auto *shh1_393 = buffer.data(shh1 + 393);
    const auto *shh1_395 = buffer.data(shh1 + 395);
    const auto *shh1_396 = buffer.data(shh1 + 396);
    const auto *shh1_397 = buffer.data(shh1 + 397);
    const auto *shh1_398 = buffer.data(shh1 + 398);
    const auto *shh1_402 = buffer.data(shh1 + 402);
    const auto *shh1_405 = buffer.data(shh1 + 405);
    const auto *shh1_409 = buffer.data(shh1 + 409);
    const auto *shh1_411 = buffer.data(shh1 + 411);
    const auto *shh1_414 = buffer.data(shh1 + 414);
    const auto *shh1_416 = buffer.data(shh1 + 416);
    const auto *shh1_417 = buffer.data(shh1 + 417);

    const auto *shi_458 = buffer.data(shi + 458);
    const auto *shi_462 = buffer.data(shi + 462);
    const auto *shi_465 = buffer.data(shi + 465);
    const auto *shi_466 = buffer.data(shi + 466);
    const auto *shi_468 = buffer.data(shi + 468);
    const auto *shi_469 = buffer.data(shi + 469);
    const auto *shi_470 = buffer.data(shi + 470);
    const auto *shi_471 = buffer.data(shi + 471);
    const auto *shi_472 = buffer.data(shi + 472);
    const auto *shi_473 = buffer.data(shi + 473);
    const auto *shi_474 = buffer.data(shi + 474);
    const auto *shi_475 = buffer.data(shi + 475);
    const auto *shi_476 = buffer.data(shi + 476);
    const auto *shi_478 = buffer.data(shi + 478);
    const auto *shi_479 = buffer.data(shi + 479);
    const auto *shi_481 = buffer.data(shi + 481);
    const auto *shi_482 = buffer.data(shi + 482);
    const auto *shi_485 = buffer.data(shi + 485);
    const auto *shi_486 = buffer.data(shi + 486);
    const auto *shi_488 = buffer.data(shi + 488);
    const auto *shi_490 = buffer.data(shi + 490);
    const auto *shi_491 = buffer.data(shi + 491);
    const auto *shi_493 = buffer.data(shi + 493);
    const auto *shi_494 = buffer.data(shi + 494);
    const auto *shi_496 = buffer.data(shi + 496);
    const auto *shi_497 = buffer.data(shi + 497);
    const auto *shi_498 = buffer.data(shi + 498);
    const auto *shi_499 = buffer.data(shi + 499);
    const auto *shi_500 = buffer.data(shi + 500);
    const auto *shi_501 = buffer.data(shi + 501);
    const auto *shi_502 = buffer.data(shi + 502);
    const auto *shi_503 = buffer.data(shi + 503);
    const auto *shi_504 = buffer.data(shi + 504);
    const auto *shi_506 = buffer.data(shi + 506);
    const auto *shi_507 = buffer.data(shi + 507);
    const auto *shi_509 = buffer.data(shi + 509);
    const auto *shi_510 = buffer.data(shi + 510);
    const auto *shi_513 = buffer.data(shi + 513);
    const auto *shi_514 = buffer.data(shi + 514);
    const auto *shi_516 = buffer.data(shi + 516);
    const auto *shi_518 = buffer.data(shi + 518);
    const auto *shi_519 = buffer.data(shi + 519);
    const auto *shi_521 = buffer.data(shi + 521);
    const auto *shi_522 = buffer.data(shi + 522);
    const auto *shi_524 = buffer.data(shi + 524);
    const auto *shi_525 = buffer.data(shi + 525);
    const auto *shi_526 = buffer.data(shi + 526);
    const auto *shi_527 = buffer.data(shi + 527);
    const auto *shi_528 = buffer.data(shi + 528);
    const auto *shi_529 = buffer.data(shi + 529);
    const auto *shi_530 = buffer.data(shi + 530);
    const auto *shi_531 = buffer.data(shi + 531);
    const auto *shi_532 = buffer.data(shi + 532);
    const auto *shi_534 = buffer.data(shi + 534);
    const auto *shi_535 = buffer.data(shi + 535);
    const auto *shi_537 = buffer.data(shi + 537);
    const auto *shi_538 = buffer.data(shi + 538);
    const auto *shi_541 = buffer.data(shi + 541);
    const auto *shi_542 = buffer.data(shi + 542);
    const auto *shi_544 = buffer.data(shi + 544);
    const auto *shi_546 = buffer.data(shi + 546);
    const auto *shi_547 = buffer.data(shi + 547);
    const auto *shi_549 = buffer.data(shi + 549);
    const auto *shi_550 = buffer.data(shi + 550);
    const auto *shi_553 = buffer.data(shi + 553);
    const auto *shi_554 = buffer.data(shi + 554);
    const auto *shi_555 = buffer.data(shi + 555);
    const auto *shi_556 = buffer.data(shi + 556);
    const auto *shi_557 = buffer.data(shi + 557);
    const auto *shi_558 = buffer.data(shi + 558);
    const auto *shi_559 = buffer.data(shi + 559);

#pragma omp simd aligned(t_592, t_593, t_594, pc_x, pc_z, sgi_290, shh0_353, shh0_354, \
                         shh1_353, shh1_354, shi_458, shi_465, \
                         shi_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_13 * sgi_290[k]
                   + f_3 * pc_z[k] * shi_458[k];

        t_593[k] = f_10 * shh0_353[k]
                   - f_11 * shh1_353[k]
                   + f_3 * pc_x[k] * shi_465[k];

        t_594[k] = f_10 * shh0_354[k]
                   - f_11 * shh1_354[k]
                   + f_3 * pc_x[k] * shi_466[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, pc_x, pc_y, sgi_322, shh0_356, \
                         shh1_356, shi_462, shi_468, shi_469, shi_470, \
                         shi_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_16 * sgi_322[k]
                   + f_3 * pc_y[k] * shi_462[k];

        t_596[k] = f_10 * shh0_356[k]
                   - f_11 * shh1_356[k]
                   + f_3 * pc_x[k] * shi_468[k];

        t_597[k] = f_3 * pc_x[k] * shi_469[k];

        t_598[k] = f_3 * pc_x[k] * shi_470[k];

        t_599[k] = f_3 * pc_x[k] * shi_471[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, pb_z, pc_x, pc_z, sgk0_388, \
                         sgk1_388, shi_472, shi_473, shi_474, shi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_3 * pc_x[k] * shi_472[k];

        t_601[k] = f_3 * pc_x[k] * shi_473[k];

        t_602[k] = f_3 * pc_x[k] * shi_474[k];

        t_603[k] = f_3 * pc_x[k] * shi_475[k];

        t_604[k] = pb_z[k] * sgk0_388[k]
                   - f_12 * pc_z[k] * sgk1_388[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, pb_z, pc_z, sgk0_390, sgk0_391, sgi_301, \
                         sgi_302, sgi_303, sgk1_390, sgk1_391, \
                         shi_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = f_13 * sgi_301[k]
                   + f_3 * pc_z[k] * shi_469[k];

        t_606[k] = pb_z[k] * sgk0_390[k]
                   + f_14 * sgi_302[k]
                   - f_12 * pc_z[k] * sgk1_390[k];

        t_607[k] = pb_z[k] * sgk0_391[k]
                   + f_15 * sgi_303[k]
                   - f_12 * pc_z[k] * sgk1_391[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, pb_z, pc_y, pc_z, sgk0_392, sgk0_393, sgi_304, \
                         sgi_305, sgi_335, sgk1_392, sgk1_393, \
                         shi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = pb_z[k] * sgk0_392[k]
                   + f_16 * sgi_304[k]
                   - f_12 * pc_z[k] * sgk1_392[k];

        t_609[k] = pb_z[k] * sgk0_393[k]
                   + f_0 * sgi_305[k]
                   - f_12 * pc_z[k] * sgk1_393[k];

        t_610[k] = f_16 * sgi_335[k]
                   + f_3 * pc_y[k] * shi_475[k];
    }

#pragma omp simd aligned(t_611, t_612, t_613, t_614, pc_x, pc_y, pc_z, sgi_307, sgi_308, \
                         sgi_336, shh0_356, shh0_357, shh1_356, shh1_357, shi_475, \
                         shi_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = f_13 * sgi_307[k]
                   + f_1 * shh0_356[k]
                   - f_2 * shh1_356[k]
                   + f_3 * pc_z[k] * shi_475[k];

        t_612[k] = f_1 * shh0_357[k]
                   - f_2 * shh1_357[k]
                   + f_3 * pc_x[k] * shi_476[k];

        t_613[k] = f_15 * sgi_336[k]
                   + f_3 * pc_y[k] * shi_476[k];

        t_614[k] = f_14 * sgi_308[k]
                   + f_3 * pc_z[k] * shi_476[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pc_x, pc_y, sgi_338, shh0_360, shh0_362, \
                         shh1_360, shh1_362, shi_478, shi_479, \
                         shi_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_4 * shh0_360[k]
                   - f_5 * shh1_360[k]
                   + f_3 * pc_x[k] * shi_479[k];

        t_616[k] = f_15 * sgi_338[k]
                   + f_3 * pc_y[k] * shi_478[k];

        t_617[k] = f_4 * shh0_362[k]
                   - f_5 * shh1_362[k]
                   + f_3 * pc_x[k] * shi_481[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, pc_x, pc_y, pc_z, sgi_311, sgi_341, shh0_363, \
                         shh1_363, shi_479, shi_481, shi_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_6 * shh0_363[k]
                   - f_7 * shh1_363[k]
                   + f_3 * pc_x[k] * shi_482[k];

        t_619[k] = f_14 * sgi_311[k]
                   + f_3 * pc_z[k] * shi_479[k];

        t_620[k] = f_15 * sgi_341[k]
                   + f_3 * pc_y[k] * shi_481[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, pc_x, pc_z, sgi_314, shh0_366, shh0_367, \
                         shh1_366, shh1_367, shi_482, shi_485, \
                         shi_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_6 * shh0_366[k]
                   - f_7 * shh1_366[k]
                   + f_3 * pc_x[k] * shi_485[k];

        t_622[k] = f_8 * shh0_367[k]
                   - f_9 * shh1_367[k]
                   + f_3 * pc_x[k] * shi_486[k];

        t_623[k] = f_14 * sgi_314[k]
                   + f_3 * pc_z[k] * shi_482[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pc_x, pc_y, sgi_345, shh0_369, shh0_371, \
                         shh1_369, shh1_371, shi_485, shi_488, \
                         shi_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_8 * shh0_369[k]
                   - f_9 * shh1_369[k]
                   + f_3 * pc_x[k] * shi_488[k];

        t_625[k] = f_15 * sgi_345[k]
                   + f_3 * pc_y[k] * shi_485[k];

        t_626[k] = f_8 * shh0_371[k]
                   - f_9 * shh1_371[k]
                   + f_3 * pc_x[k] * shi_490[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pc_x, pc_z, sgi_318, shh0_372, shh0_374, \
                         shh1_372, shh1_374, shi_486, shi_491, \
                         shi_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_10 * shh0_372[k]
                   - f_11 * shh1_372[k]
                   + f_3 * pc_x[k] * shi_491[k];

        t_628[k] = f_14 * sgi_318[k]
                   + f_3 * pc_z[k] * shi_486[k];

        t_629[k] = f_10 * shh0_374[k]
                   - f_11 * shh1_374[k]
                   + f_3 * pc_x[k] * shi_493[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, pc_y, sgi_350, shh0_375, shh0_377, \
                         shh1_375, shh1_377, shi_490, shi_494, shi_496, \
                         shi_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_10 * shh0_375[k]
                   - f_11 * shh1_375[k]
                   + f_3 * pc_x[k] * shi_494[k];

        t_631[k] = f_15 * sgi_350[k]
                   + f_3 * pc_y[k] * shi_490[k];

        t_632[k] = f_10 * shh0_377[k]
                   - f_11 * shh1_377[k]
                   + f_3 * pc_x[k] * shi_496[k];

        t_633[k] = f_3 * pc_x[k] * shi_497[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, t_639, pc_x, shi_498, shi_499, \
                         shi_500, shi_501, shi_502, shi_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_3 * pc_x[k] * shi_498[k];

        t_635[k] = f_3 * pc_x[k] * shi_499[k];

        t_636[k] = f_3 * pc_x[k] * shi_500[k];

        t_637[k] = f_3 * pc_x[k] * shi_501[k];

        t_638[k] = f_3 * pc_x[k] * shi_502[k];

        t_639[k] = f_3 * pc_x[k] * shi_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pc_y, pc_z, sgi_329, sgi_357, sgi_359, shh0_372, \
                         shh0_374, shh1_372, shh1_374, shi_497, \
                         shi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_15 * sgi_357[k]
                   + f_1 * shh0_372[k]
                   - f_2 * shh1_372[k]
                   + f_3 * pc_y[k] * shi_497[k];

        t_641[k] = f_14 * sgi_329[k]
                   + f_3 * pc_z[k] * shi_497[k];

        t_642[k] = f_15 * sgi_359[k]
                   + f_4 * shh0_374[k]
                   - f_5 * shh1_374[k]
                   + f_3 * pc_y[k] * shi_499[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pc_y, sgi_360, sgi_361, sgi_362, shh0_375, \
                         shh0_376, shh0_377, shh1_375, shh1_376, shh1_377, shi_500, shi_501, \
                         shi_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_15 * sgi_360[k]
                   + f_6 * shh0_375[k]
                   - f_7 * shh1_375[k]
                   + f_3 * pc_y[k] * shi_500[k];

        t_644[k] = f_15 * sgi_361[k]
                   + f_8 * shh0_376[k]
                   - f_9 * shh1_376[k]
                   + f_3 * pc_y[k] * shi_501[k];

        t_645[k] = f_15 * sgi_362[k]
                   + f_10 * shh0_377[k]
                   - f_11 * shh1_377[k]
                   + f_3 * pc_y[k] * shi_502[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, pc_x, pc_y, pc_z, sgi_335, sgi_363, \
                         sgi_364, shh0_377, shh0_378, shh1_377, shh1_378, shi_503, \
                         shi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_15 * sgi_363[k]
                   + f_3 * pc_y[k] * shi_503[k];

        t_647[k] = f_14 * sgi_335[k]
                   + f_1 * shh0_377[k]
                   - f_2 * shh1_377[k]
                   + f_3 * pc_z[k] * shi_503[k];

        t_648[k] = f_1 * shh0_378[k]
                   - f_2 * shh1_378[k]
                   + f_3 * pc_x[k] * shi_504[k];

        t_649[k] = f_14 * sgi_364[k]
                   + f_3 * pc_y[k] * shi_504[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pc_x, pc_y, pc_z, sgi_336, sgi_366, shh0_381, \
                         shh1_381, shi_504, shi_506, shi_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_15 * sgi_336[k]
                   + f_3 * pc_z[k] * shi_504[k];

        t_651[k] = f_4 * shh0_381[k]
                   - f_5 * shh1_381[k]
                   + f_3 * pc_x[k] * shi_507[k];

        t_652[k] = f_14 * sgi_366[k]
                   + f_3 * pc_y[k] * shi_506[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, t_656, pc_x, pc_y, pc_z, sgi_339, sgi_369, \
                         shh0_383, shh0_384, shh1_383, shh1_384, shi_507, shi_509, \
                         shi_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_4 * shh0_383[k]
                   - f_5 * shh1_383[k]
                   + f_3 * pc_x[k] * shi_509[k];

        t_654[k] = f_6 * shh0_384[k]
                   - f_7 * shh1_384[k]
                   + f_3 * pc_x[k] * shi_510[k];

        t_655[k] = f_15 * sgi_339[k]
                   + f_3 * pc_z[k] * shi_507[k];

        t_656[k] = f_14 * sgi_369[k]
                   + f_3 * pc_y[k] * shi_509[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, pc_x, pc_z, sgi_342, shh0_387, shh0_388, \
                         shh1_387, shh1_388, shi_510, shi_513, \
                         shi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_6 * shh0_387[k]
                   - f_7 * shh1_387[k]
                   + f_3 * pc_x[k] * shi_513[k];

        t_658[k] = f_8 * shh0_388[k]
                   - f_9 * shh1_388[k]
                   + f_3 * pc_x[k] * shi_514[k];

        t_659[k] = f_15 * sgi_342[k]
                   + f_3 * pc_z[k] * shi_510[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, pc_x, pc_y, sgi_373, shh0_390, shh0_392, \
                         shh1_390, shh1_392, shi_513, shi_516, \
                         shi_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_8 * shh0_390[k]
                   - f_9 * shh1_390[k]
                   + f_3 * pc_x[k] * shi_516[k];

        t_661[k] = f_14 * sgi_373[k]
                   + f_3 * pc_y[k] * shi_513[k];

        t_662[k] = f_8 * shh0_392[k]
                   - f_9 * shh1_392[k]
                   + f_3 * pc_x[k] * shi_518[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, pc_x, pc_z, sgi_346, shh0_393, shh0_395, \
                         shh1_393, shh1_395, shi_514, shi_519, \
                         shi_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_10 * shh0_393[k]
                   - f_11 * shh1_393[k]
                   + f_3 * pc_x[k] * shi_519[k];

        t_664[k] = f_15 * sgi_346[k]
                   + f_3 * pc_z[k] * shi_514[k];

        t_665[k] = f_10 * shh0_395[k]
                   - f_11 * shh1_395[k]
                   + f_3 * pc_x[k] * shi_521[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pc_x, pc_y, sgi_378, shh0_396, shh0_398, \
                         shh1_396, shh1_398, shi_518, shi_522, shi_524, \
                         shi_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_10 * shh0_396[k]
                   - f_11 * shh1_396[k]
                   + f_3 * pc_x[k] * shi_522[k];

        t_667[k] = f_14 * sgi_378[k]
                   + f_3 * pc_y[k] * shi_518[k];

        t_668[k] = f_10 * shh0_398[k]
                   - f_11 * shh1_398[k]
                   + f_3 * pc_x[k] * shi_524[k];

        t_669[k] = f_3 * pc_x[k] * shi_525[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, t_675, pc_x, shi_526, shi_527, \
                         shi_528, shi_529, shi_530, shi_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_3 * pc_x[k] * shi_526[k];

        t_671[k] = f_3 * pc_x[k] * shi_527[k];

        t_672[k] = f_3 * pc_x[k] * shi_528[k];

        t_673[k] = f_3 * pc_x[k] * shi_529[k];

        t_674[k] = f_3 * pc_x[k] * shi_530[k];

        t_675[k] = f_3 * pc_x[k] * shi_531[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pc_y, pc_z, sgi_357, sgi_385, sgi_387, shh0_393, \
                         shh0_395, shh1_393, shh1_395, shi_525, \
                         shi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_14 * sgi_385[k]
                   + f_1 * shh0_393[k]
                   - f_2 * shh1_393[k]
                   + f_3 * pc_y[k] * shi_525[k];

        t_677[k] = f_15 * sgi_357[k]
                   + f_3 * pc_z[k] * shi_525[k];

        t_678[k] = f_14 * sgi_387[k]
                   + f_4 * shh0_395[k]
                   - f_5 * shh1_395[k]
                   + f_3 * pc_y[k] * shi_527[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, pc_y, sgi_388, sgi_389, sgi_390, shh0_396, \
                         shh0_397, shh0_398, shh1_396, shh1_397, shh1_398, shi_528, shi_529, \
                         shi_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_14 * sgi_388[k]
                   + f_6 * shh0_396[k]
                   - f_7 * shh1_396[k]
                   + f_3 * pc_y[k] * shi_528[k];

        t_680[k] = f_14 * sgi_389[k]
                   + f_8 * shh0_397[k]
                   - f_9 * shh1_397[k]
                   + f_3 * pc_y[k] * shi_529[k];

        t_681[k] = f_14 * sgi_390[k]
                   + f_10 * shh0_398[k]
                   - f_11 * shh1_398[k]
                   + f_3 * pc_y[k] * shi_530[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, pb_y, pc_y, pc_z, sgk0_504, sgi_363, \
                         sgi_391, sgi_392, sgk1_504, shh0_398, shh1_398, shi_531, \
                         shi_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_14 * sgi_391[k]
                   + f_3 * pc_y[k] * shi_531[k];

        t_683[k] = f_15 * sgi_363[k]
                   + f_1 * shh0_398[k]
                   - f_2 * shh1_398[k]
                   + f_3 * pc_z[k] * shi_531[k];

        t_684[k] = pb_y[k] * sgk0_504[k]
                   - f_12 * pc_y[k] * sgk1_504[k];

        t_685[k] = f_13 * sgi_392[k]
                   + f_3 * pc_y[k] * shi_532[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, pc_x, pc_y, pc_z, sgi_364, sgi_394, shh0_402, \
                         shh1_402, shi_532, shi_534, shi_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_16 * sgi_364[k]
                   + f_3 * pc_z[k] * shi_532[k];

        t_687[k] = f_4 * shh0_402[k]
                   - f_5 * shh1_402[k]
                   + f_3 * pc_x[k] * shi_535[k];

        t_688[k] = f_13 * sgi_394[k]
                   + f_3 * pc_y[k] * shi_534[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, pb_y, pc_x, pc_y, pc_z, sgk0_509, sgi_367, \
                         sgk1_509, shh0_405, shh1_405, shi_535, \
                         shi_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = pb_y[k] * sgk0_509[k]
                   - f_12 * pc_y[k] * sgk1_509[k];

        t_690[k] = f_6 * shh0_405[k]
                   - f_7 * shh1_405[k]
                   + f_3 * pc_x[k] * shi_538[k];

        t_691[k] = f_16 * sgi_367[k]
                   + f_3 * pc_z[k] * shi_535[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, pb_y, pc_x, pc_y, sgk0_513, sgi_397, sgk1_513, \
                         shh0_409, shh1_409, shi_537, shi_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = f_13 * sgi_397[k]
                   + f_3 * pc_y[k] * shi_537[k];

        t_693[k] = pb_y[k] * sgk0_513[k]
                   - f_12 * pc_y[k] * sgk1_513[k];

        t_694[k] = f_8 * shh0_409[k]
                   - f_9 * shh1_409[k]
                   + f_3 * pc_x[k] * shi_542[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, pc_x, pc_y, pc_z, sgi_370, sgi_401, shh0_411, \
                         shh1_411, shi_538, shi_541, shi_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_16 * sgi_370[k]
                   + f_3 * pc_z[k] * shi_538[k];

        t_696[k] = f_8 * shh0_411[k]
                   - f_9 * shh1_411[k]
                   + f_3 * pc_x[k] * shi_544[k];

        t_697[k] = f_13 * sgi_401[k]
                   + f_3 * pc_y[k] * shi_541[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pb_y, pc_x, pc_y, pc_z, sgk0_518, sgi_374, \
                         sgk1_518, shh0_414, shh1_414, shi_542, \
                         shi_547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = pb_y[k] * sgk0_518[k]
                   - f_12 * pc_y[k] * sgk1_518[k];

        t_699[k] = f_10 * shh0_414[k]
                   - f_11 * shh1_414[k]
                   + f_3 * pc_x[k] * shi_547[k];

        t_700[k] = f_16 * sgi_374[k]
                   + f_3 * pc_z[k] * shi_542[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, pc_x, pc_y, sgi_406, shh0_416, shh0_417, \
                         shh1_416, shh1_417, shi_546, shi_549, \
                         shi_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_10 * shh0_416[k]
                   - f_11 * shh1_416[k]
                   + f_3 * pc_x[k] * shi_549[k];

        t_702[k] = f_10 * shh0_417[k]
                   - f_11 * shh1_417[k]
                   + f_3 * pc_x[k] * shi_550[k];

        t_703[k] = f_13 * sgi_406[k]
                   + f_3 * pc_y[k] * shi_546[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, t_708, t_709, pb_y, pc_x, pc_y, sgk0_524, \
                         sgk1_524, shi_553, shi_554, shi_555, shi_556, \
                         shi_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = pb_y[k] * sgk0_524[k]
                   - f_12 * pc_y[k] * sgk1_524[k];

        t_705[k] = f_3 * pc_x[k] * shi_553[k];

        t_706[k] = f_3 * pc_x[k] * shi_554[k];

        t_707[k] = f_3 * pc_x[k] * shi_555[k];

        t_708[k] = f_3 * pc_x[k] * shi_556[k];

        t_709[k] = f_3 * pc_x[k] * shi_557[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pb_y, pc_x, pc_y, pc_z, sgk0_532, \
                         sgi_385, sgi_413, sgk1_532, shi_553, shi_558, \
                         shi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_3 * pc_x[k] * shi_558[k];

        t_711[k] = f_3 * pc_x[k] * shi_559[k];

        t_712[k] = pb_y[k] * sgk0_532[k]
                   + f_17 * sgi_413[k]
                   - f_12 * pc_y[k] * sgk1_532[k];

        t_713[k] = f_16 * sgi_385[k]
                   + f_3 * pc_z[k] * shi_553[k];
    }
}

static auto
compute_prim_shk_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgk0,
                                                          const size_t sgi, const size_t sgk1,
                                                          const size_t shh0, const size_t shh1,
                                                          const size_t shi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgk0_534 = buffer.data(sgk0 + 534);
    const auto *sgk0_535 = buffer.data(sgk0 + 535);
    const auto *sgk0_536 = buffer.data(sgk0 + 536);
    const auto *sgk0_537 = buffer.data(sgk0 + 537);
    const auto *sgk0_539 = buffer.data(sgk0 + 539);

    const auto *sgi_392 = buffer.data(sgi + 392);
    const auto *sgi_395 = buffer.data(sgi + 395);
    const auto *sgi_398 = buffer.data(sgi + 398);
    const auto *sgi_402 = buffer.data(sgi + 402);
    const auto *sgi_413 = buffer.data(sgi + 413);
    const auto *sgi_415 = buffer.data(sgi + 415);
    const auto *sgi_416 = buffer.data(sgi + 416);
    const auto *sgi_417 = buffer.data(sgi + 417);
    const auto *sgi_418 = buffer.data(sgi + 418);
    const auto *sgi_419 = buffer.data(sgi + 419);

    const auto *sgk1_534 = buffer.data(sgk1 + 534);
    const auto *sgk1_535 = buffer.data(sgk1 + 535);
    const auto *sgk1_536 = buffer.data(sgk1 + 536);
    const auto *sgk1_537 = buffer.data(sgk1 + 537);
    const auto *sgk1_539 = buffer.data(sgk1 + 539);

    const auto *shh0_420 = buffer.data(shh0 + 420);
    const auto *shh0_423 = buffer.data(shh0 + 423);
    const auto *shh0_425 = buffer.data(shh0 + 425);
    const auto *shh0_426 = buffer.data(shh0 + 426);
    const auto *shh0_429 = buffer.data(shh0 + 429);
    const auto *shh0_430 = buffer.data(shh0 + 430);
    const auto *shh0_432 = buffer.data(shh0 + 432);
    const auto *shh0_434 = buffer.data(shh0 + 434);
    const auto *shh0_435 = buffer.data(shh0 + 435);
    const auto *shh0_437 = buffer.data(shh0 + 437);
    const auto *shh0_438 = buffer.data(shh0 + 438);
    const auto *shh0_439 = buffer.data(shh0 + 439);
    const auto *shh0_440 = buffer.data(shh0 + 440);

    const auto *shh1_420 = buffer.data(shh1 + 420);
    const auto *shh1_423 = buffer.data(shh1 + 423);
    const auto *shh1_425 = buffer.data(shh1 + 425);
    const auto *shh1_426 = buffer.data(shh1 + 426);
    const auto *shh1_429 = buffer.data(shh1 + 429);
    const auto *shh1_430 = buffer.data(shh1 + 430);
    const auto *shh1_432 = buffer.data(shh1 + 432);
    const auto *shh1_434 = buffer.data(shh1 + 434);
    const auto *shh1_435 = buffer.data(shh1 + 435);
    const auto *shh1_437 = buffer.data(shh1 + 437);
    const auto *shh1_438 = buffer.data(shh1 + 438);
    const auto *shh1_439 = buffer.data(shh1 + 439);
    const auto *shh1_440 = buffer.data(shh1 + 440);

    const auto *shi_559 = buffer.data(shi + 559);
    const auto *shi_560 = buffer.data(shi + 560);
    const auto *shi_562 = buffer.data(shi + 562);
    const auto *shi_563 = buffer.data(shi + 563);
    const auto *shi_565 = buffer.data(shi + 565);
    const auto *shi_566 = buffer.data(shi + 566);
    const auto *shi_569 = buffer.data(shi + 569);
    const auto *shi_570 = buffer.data(shi + 570);
    const auto *shi_572 = buffer.data(shi + 572);
    const auto *shi_574 = buffer.data(shi + 574);
    const auto *shi_575 = buffer.data(shi + 575);
    const auto *shi_577 = buffer.data(shi + 577);
    const auto *shi_578 = buffer.data(shi + 578);
    const auto *shi_580 = buffer.data(shi + 580);
    const auto *shi_581 = buffer.data(shi + 581);
    const auto *shi_582 = buffer.data(shi + 582);
    const auto *shi_583 = buffer.data(shi + 583);
    const auto *shi_584 = buffer.data(shi + 584);
    const auto *shi_585 = buffer.data(shi + 585);
    const auto *shi_586 = buffer.data(shi + 586);
    const auto *shi_587 = buffer.data(shi + 587);

#pragma omp simd aligned(t_714, t_715, t_716, pb_y, pc_y, sgk0_534, sgk0_535, sgk0_536, \
                         sgi_415, sgi_416, sgi_417, sgk1_534, sgk1_535, \
                         sgk1_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = pb_y[k] * sgk0_534[k]
                   + f_0 * sgi_415[k]
                   - f_12 * pc_y[k] * sgk1_534[k];

        t_715[k] = pb_y[k] * sgk0_535[k]
                   + f_16 * sgi_416[k]
                   - f_12 * pc_y[k] * sgk1_535[k];

        t_716[k] = pb_y[k] * sgk0_536[k]
                   + f_15 * sgi_417[k]
                   - f_12 * pc_y[k] * sgk1_536[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pb_y, pc_y, sgk0_537, sgk0_539, sgi_418, \
                         sgi_419, sgk1_537, sgk1_539, shi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = pb_y[k] * sgk0_537[k]
                   + f_14 * sgi_418[k]
                   - f_12 * pc_y[k] * sgk1_537[k];

        t_718[k] = f_13 * sgi_419[k]
                   + f_3 * pc_y[k] * shi_559[k];

        t_719[k] = pb_y[k] * sgk0_539[k]
                   - f_12 * pc_y[k] * sgk1_539[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, pc_x, pc_y, pc_z, sgi_392, \
                         shh0_420, shh0_423, shh1_420, shh1_423, shi_560, shi_562, \
                         shi_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_1 * shh0_420[k]
                   - f_2 * shh1_420[k]
                   + f_3 * pc_x[k] * shi_560[k];

        t_721[k] = f_3 * pc_y[k] * shi_560[k];

        t_722[k] = f_0 * sgi_392[k]
                   + f_3 * pc_z[k] * shi_560[k];

        t_723[k] = f_4 * shh0_423[k]
                   - f_5 * shh1_423[k]
                   + f_3 * pc_x[k] * shi_563[k];

        t_724[k] = f_3 * pc_y[k] * shi_562[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pc_x, pc_y, pc_z, sgi_395, shh0_425, \
                         shh0_426, shh1_425, shh1_426, shi_563, shi_565, \
                         shi_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_4 * shh0_425[k]
                   - f_5 * shh1_425[k]
                   + f_3 * pc_x[k] * shi_565[k];

        t_726[k] = f_6 * shh0_426[k]
                   - f_7 * shh1_426[k]
                   + f_3 * pc_x[k] * shi_566[k];

        t_727[k] = f_0 * sgi_395[k]
                   + f_3 * pc_z[k] * shi_563[k];

        t_728[k] = f_3 * pc_y[k] * shi_565[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_x, pc_z, sgi_398, shh0_429, shh0_430, \
                         shh1_429, shh1_430, shi_566, shi_569, \
                         shi_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_6 * shh0_429[k]
                   - f_7 * shh1_429[k]
                   + f_3 * pc_x[k] * shi_569[k];

        t_730[k] = f_8 * shh0_430[k]
                   - f_9 * shh1_430[k]
                   + f_3 * pc_x[k] * shi_570[k];

        t_731[k] = f_0 * sgi_398[k]
                   + f_3 * pc_z[k] * shi_566[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, t_735, pc_x, pc_y, shh0_432, shh0_434, shh0_435, \
                         shh1_432, shh1_434, shh1_435, shi_569, shi_572, shi_574, \
                         shi_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_8 * shh0_432[k]
                   - f_9 * shh1_432[k]
                   + f_3 * pc_x[k] * shi_572[k];

        t_733[k] = f_3 * pc_y[k] * shi_569[k];

        t_734[k] = f_8 * shh0_434[k]
                   - f_9 * shh1_434[k]
                   + f_3 * pc_x[k] * shi_574[k];

        t_735[k] = f_10 * shh0_435[k]
                   - f_11 * shh1_435[k]
                   + f_3 * pc_x[k] * shi_575[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, t_739, pc_x, pc_y, pc_z, sgi_402, shh0_437, \
                         shh0_438, shh1_437, shh1_438, shi_570, shi_574, shi_577, \
                         shi_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = f_0 * sgi_402[k]
                   + f_3 * pc_z[k] * shi_570[k];

        t_737[k] = f_10 * shh0_437[k]
                   - f_11 * shh1_437[k]
                   + f_3 * pc_x[k] * shi_577[k];

        t_738[k] = f_10 * shh0_438[k]
                   - f_11 * shh1_438[k]
                   + f_3 * pc_x[k] * shi_578[k];

        t_739[k] = f_3 * pc_y[k] * shi_574[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, t_745, pc_x, shh0_440, shh1_440, \
                         shi_580, shi_581, shi_582, shi_583, shi_584, \
                         shi_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_10 * shh0_440[k]
                   - f_11 * shh1_440[k]
                   + f_3 * pc_x[k] * shi_580[k];

        t_741[k] = f_3 * pc_x[k] * shi_581[k];

        t_742[k] = f_3 * pc_x[k] * shi_582[k];

        t_743[k] = f_3 * pc_x[k] * shi_583[k];

        t_744[k] = f_3 * pc_x[k] * shi_584[k];

        t_745[k] = f_3 * pc_x[k] * shi_585[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pc_x, pc_y, pc_z, sgi_413, shh0_435, \
                         shh1_435, shi_581, shi_586, shi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_3 * pc_x[k] * shi_586[k];

        t_747[k] = f_3 * pc_x[k] * shi_587[k];

        t_748[k] = f_1 * shh0_435[k]
                   - f_2 * shh1_435[k]
                   + f_3 * pc_y[k] * shi_581[k];

        t_749[k] = f_0 * sgi_413[k]
                   + f_3 * pc_z[k] * shi_581[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, pc_y, shh0_437, shh0_438, shh0_439, shh1_437, \
                         shh1_438, shh1_439, shi_583, shi_584, \
                         shi_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_4 * shh0_437[k]
                   - f_5 * shh1_437[k]
                   + f_3 * pc_y[k] * shi_583[k];

        t_751[k] = f_6 * shh0_438[k]
                   - f_7 * shh1_438[k]
                   + f_3 * pc_y[k] * shi_584[k];

        t_752[k] = f_8 * shh0_439[k]
                   - f_9 * shh1_439[k]
                   + f_3 * pc_y[k] * shi_585[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, pc_y, pc_z, sgi_419, shh0_440, shh1_440, \
                         shi_586, shi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_10 * shh0_440[k]
                   - f_11 * shh1_440[k]
                   + f_3 * pc_y[k] * shi_586[k];

        t_754[k] = f_3 * pc_y[k] * shi_587[k];

        t_755[k] = f_0 * sgi_419[k]
                   + f_1 * shh0_440[k]
                   - f_2 * shh1_440[k]
                   + f_3 * pc_z[k] * shi_587[k];
    }
}

auto
compute_prim_shk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sgk0, const size_t sgi,
                                                   const size_t sgk1, const size_t shh0,
                                                   const size_t shh1, const size_t shi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_shk_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sgk0, sgi,
                                                              sgk1, shh0, shh1, shi, ncols,
                                                              gamma, p, q);

    compute_prim_shk_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sgk0, sgi,
                                                              sgk1, shh0, shh1, shi, ncols,
                                                              gamma, p, q);

    compute_prim_shk_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sgk0, sgi,
                                                              sgk1, shh0, shh1, shi, ncols,
                                                              gamma, p, q);

    compute_prim_shk_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sgk0, sgi,
                                                              sgk1, shh0, shh1, shi, ncols,
                                                              gamma, p, q);

    compute_prim_shk_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sgk0, sgi,
                                                              sgk1, shh0, shh1, shi, ncols,
                                                              gamma, p, q);

    compute_prim_shk_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, sgk0, sgi,
                                                              sgk1, shh0, shh1, shi, ncols,
                                                              gamma, p, q);

    compute_prim_shk_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, sgk0, sgi,
                                                              sgk1, shh0, shh1, shi, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
